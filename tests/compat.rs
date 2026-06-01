//! Differential compat against limma voom.
//!
//! - `golden_diff_*` always run: ours vs committed R-captured goldens.
//! - `live_r_*` run only when `conda run -n r-bioc Rscript` is available;
//!   they regenerate the oracle and diff against ours (loud-skip otherwise).

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

const EPS: f64 = 1e-6; // relative

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-voom"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn manifest(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

type Matrix = (Vec<String>, BTreeMap<String, Vec<f64>>);

fn parse(text: &str) -> Matrix {
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let mut rows = BTreeMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f.next().unwrap().to_string();
        let vals: Vec<f64> = f.map(|s| s.trim().parse().unwrap()).collect();
        rows.insert(gene, vals);
    }
    (header, rows)
}

fn assert_close(ours: &Matrix, theirs: &Matrix, label: &str) {
    assert_eq!(ours.0, theirs.0, "{label}: header mismatch");
    assert_eq!(ours.1.len(), theirs.1.len(), "{label}: row count mismatch");
    let mut max_rel = 0.0f64;
    for (gene, a) in &ours.1 {
        let b = theirs
            .1
            .get(gene)
            .unwrap_or_else(|| panic!("{label}: missing gene {gene}"));
        assert_eq!(a.len(), b.len(), "{label}: {gene} width mismatch");
        for (x, y) in a.iter().zip(b) {
            let rel = (x - y).abs() / y.abs().max(1e-12);
            max_rel = max_rel.max(rel);
            assert!(rel < EPS, "{label}: {gene} ours={x} theirs={y} rel={rel:e}");
        }
    }
    eprintln!("{label}: max relative deviation = {max_rel:e}");
}

fn run_ours(counts: &PathBuf) -> (String, String) {
    let scratch = std::env::temp_dir();
    let wpath = scratch.join(format!("voom_w_{}.tsv", std::process::id()));
    let out = Command::new(ours())
        .arg(counts)
        .args(["-w", wpath.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let e = String::from_utf8(out.stdout).unwrap();
    let w = std::fs::read_to_string(&wpath).unwrap();
    let _ = std::fs::remove_file(&wpath);
    (e, w)
}

#[test]
fn golden_diff_e_and_weights() {
    let (e, w) = run_ours(&golden("counts.tsv"));
    let e_exp = std::fs::read_to_string(golden("E.expected.tsv")).unwrap();
    let w_exp = std::fs::read_to_string(golden("weights.expected.tsv")).unwrap();
    assert_close(&parse(&e), &parse(&e_exp), "E (golden)");
    assert_close(&parse(&w), &parse(&w_exp), "weights (golden)");
}

fn rscript_available() -> bool {
    Command::new("conda")
        .args(["run", "-n", "r-bioc", "Rscript", "-e", "cat(1)"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[test]
fn live_r_matches_ours() {
    if !rscript_available() {
        eprintln!("SKIP live_r_matches_ours: `conda run -n r-bioc Rscript` unavailable");
        return;
    }
    let scratch = std::env::temp_dir();
    let e_r = scratch.join(format!("voom_er_{}.tsv", std::process::id()));
    let w_r = scratch.join(format!("voom_wr_{}.tsv", std::process::id()));
    let oracle = Command::new("conda")
        .args(["run", "-n", "r-bioc", "Rscript"])
        .arg(manifest("tests/voom_oracle.R"))
        .arg(golden("counts.tsv"))
        .arg(&e_r)
        .arg(&w_r)
        .output()
        .unwrap();
    assert!(
        oracle.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&oracle.stderr)
    );

    let (e, w) = run_ours(&golden("counts.tsv"));
    let e_ref = std::fs::read_to_string(&e_r).unwrap();
    let w_ref = std::fs::read_to_string(&w_r).unwrap();
    let _ = std::fs::remove_file(&e_r);
    let _ = std::fs::remove_file(&w_r);
    assert_close(&parse(&e), &parse(&e_ref), "E (live R)");
    assert_close(&parse(&w), &parse(&w_ref), "weights (live R)");
}
