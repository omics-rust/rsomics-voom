//! voom: log2-CPM transform plus mean-variance precision weights.
//!
//! Reimplements limma's `voom()` for an intercept-only design (one mean per
//! gene). Reference: Law, Chen, Shi & Smyth (2014), Genome Biology 15:R29,
//! "voom: precision weights unlock linear model analysis tools for RNA-seq read
//! counts" (DOI 10.1186/gb-2014-15-2-r29).

mod lowess;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// voom's adaptive LOWESS span (limma `chooseLowessSpan`, defaults
/// small.n=50, min.span=0.3, power=1/3). `n` is the full gene count.
fn lowess_span(n: usize) -> f64 {
    (0.3 + 0.7 * (50.0 / n as f64).powf(1.0 / 3.0)).min(1.0)
}

pub struct CountMatrix {
    pub header: String,
    pub genes: Vec<String>,
    pub samples: usize,
    /// row-major: counts[gene][sample]
    pub counts: Vec<Vec<f64>>,
}

pub fn read_matrix(path: &Path) -> Result<CountMatrix> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut lines = BufReader::new(file).lines();

    let header = lines
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?
        .map_err(RsomicsError::Io)?;
    let samples = header.split('\t').count() - 1;
    if samples < 1 {
        return Err(RsomicsError::InvalidInput(
            "count matrix needs at least one sample column".into(),
        ));
    }

    let mut genes = Vec::new();
    let mut counts = Vec::new();
    for line in lines {
        let line = line.map_err(RsomicsError::Io)?;
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let gene = fields
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("missing gene id".into()))?;
        let row: Vec<f64> = fields
            .map(|s| {
                s.parse::<f64>()
                    .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric count '{s}'")))
            })
            .collect::<Result<_>>()?;
        if row.len() != samples {
            return Err(RsomicsError::InvalidInput(format!(
                "gene '{gene}' has {} values, header declares {samples} samples",
                row.len()
            )));
        }
        if row.iter().any(|&c| c < 0.0) {
            return Err(RsomicsError::InvalidInput(format!(
                "negative count in gene '{gene}'"
            )));
        }
        genes.push(gene.to_string());
        counts.push(row);
    }

    if genes.len() < 2 {
        return Err(RsomicsError::InvalidInput(
            "need at least two genes to fit a mean-variance trend".into(),
        ));
    }

    Ok(CountMatrix {
        header,
        genes,
        samples,
        counts,
    })
}

pub struct Voom {
    /// log2-CPM, row-major [gene][sample]
    pub e: Vec<Vec<f64>>,
    /// precision weights, row-major [gene][sample]
    pub weights: Vec<Vec<f64>>,
}

pub fn voom(m: &CountMatrix) -> Result<Voom> {
    let n = m.samples;
    let ng = m.genes.len();
    if n < 2 {
        return Err(RsomicsError::InvalidInput(
            "intercept-only voom needs at least two samples for replication".into(),
        ));
    }

    let lib_size: Vec<f64> = (0..n)
        .map(|j| m.counts.iter().map(|row| row[j]).sum())
        .collect();
    let log_lib: Vec<f64> = lib_size.iter().map(|&l| (l + 1.0).log2()).collect();
    let mean_log_lib: f64 = log_lib.iter().sum::<f64>() / n as f64;

    // E = log2((count + 0.5) / (lib.size + 1) * 1e6), per column
    let mut e = vec![vec![0.0; n]; ng];
    for (gi, row) in m.counts.iter().enumerate() {
        for j in 0..n {
            e[gi][j] = ((row[j] + 0.5) / (lib_size[j] + 1.0) * 1e6).log2();
        }
    }

    // intercept-only lmFit: per-gene mean (Amean) and residual sd (sigma)
    let amean: Vec<f64> = e
        .iter()
        .map(|row| row.iter().sum::<f64>() / n as f64)
        .collect();
    let sigma: Vec<f64> = e
        .iter()
        .zip(&amean)
        .map(|(row, &mu)| {
            let sse: f64 = row.iter().map(|&v| (v - mu) * (v - mu)).sum();
            (sse / (n - 1) as f64).sqrt()
        })
        .collect();

    // mean-variance trend points; drop all-zero-count genes
    let row_zero: Vec<bool> = m
        .counts
        .iter()
        .map(|r| r.iter().all(|&c| c == 0.0))
        .collect();
    let mut sx: Vec<f64> = Vec::with_capacity(ng);
    let mut sy: Vec<f64> = Vec::with_capacity(ng);
    for gi in 0..ng {
        if row_zero[gi] {
            continue;
        }
        sx.push(amean[gi] + mean_log_lib - 1e6_f64.log2());
        sy.push(sigma[gi].sqrt());
    }
    if sx.len() < 2 {
        return Err(RsomicsError::InvalidInput(
            "fewer than two genes with nonzero counts; cannot fit trend".into(),
        ));
    }

    let (lx, ly) = lowess_sorted(&sx, &sy, lowess_span(ng));

    // weight = 1 / trend(fitted.logcount)^4, with
    // fitted.logcount = Amean[gene] + log2(1e-6 * (lib.size[sample] + 1))
    let trend = lowess::Trend::new(&lx, &ly);
    let sample_offset: Vec<f64> = lib_size
        .iter()
        .map(|&l| (1e-6 * (l + 1.0)).log2())
        .collect();
    let mut weights = vec![vec![0.0; n]; ng];
    for (gi, wrow) in weights.iter_mut().enumerate() {
        for (j, w) in wrow.iter_mut().enumerate() {
            let predicted = trend.eval(amean[gi] + sample_offset[j]);
            *w = 1.0 / predicted.powi(4);
        }
    }

    Ok(Voom { e, weights })
}

fn lowess_sorted(sx: &[f64], sy: &[f64], span: f64) -> (Vec<f64>, Vec<f64>) {
    let mut order: Vec<usize> = (0..sx.len()).collect();
    order.sort_by(|&a, &b| sx[a].partial_cmp(&sx[b]).unwrap());
    let x: Vec<f64> = order.iter().map(|&i| sx[i]).collect();
    let y: Vec<f64> = order.iter().map(|&i| sy[i]).collect();
    let delta = 0.01 * (x[x.len() - 1] - x[0]);
    let ys = lowess::lowess(&x, &y, span, 3, delta);
    (x, ys)
}

pub fn write_matrix(
    header: &str,
    genes: &[String],
    rows: &[Vec<f64>],
    out: &mut dyn Write,
) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    writeln!(w, "{header}").map_err(RsomicsError::Io)?;
    let mut fmt = ryu::Buffer::new();
    let mut line = String::with_capacity(rows.first().map_or(0, |r| r.len()) * 16);
    for (gene, row) in genes.iter().zip(rows) {
        line.clear();
        line.push_str(gene);
        for &v in row {
            line.push('\t');
            line.push_str(fmt.format(v));
        }
        line.push('\n');
        w.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log2cpm_offset() {
        // single gene, two samples: lib.size = counts here
        let m = CountMatrix {
            header: "gene\ts1\ts2".into(),
            genes: vec!["g1".into(), "g2".into()],
            samples: 2,
            counts: vec![vec![10.0, 20.0], vec![30.0, 40.0]],
        };
        let v = voom(&m).unwrap();
        let lib0 = 40.0;
        let expect = ((10.0 + 0.5) / (lib0 + 1.0) * 1e6_f64).log2();
        assert!((v.e[0][0] - expect).abs() < 1e-12);
    }
}
