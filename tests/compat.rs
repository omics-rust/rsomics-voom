//! Differential compat against limma voom.
//!
//! - `golden_diff_*` always run: ours vs committed R-captured goldens.
//! - `live_r_*` run only when `conda run -n r-bioc Rscript` is available;
//!   they regenerate the oracle and diff against ours (loud-skip otherwise).

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

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

fn run_ours(counts: &PathBuf, extra: &[&str]) -> (String, String) {
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let scratch = std::env::temp_dir();
    let wpath = scratch.join(format!(
        "voom_w_{}_{}.tsv",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    let out = Command::new(ours())
        .arg(counts)
        .args(["-w", wpath.to_str().unwrap()])
        .args(extra)
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

/// Default `voom(counts, design)` is `span = 0.5, adaptive.span = FALSE`; the
/// weights golden here was captured from that call, so a passing diff pins the
/// default to limma's default rather than the adaptive-span mode.
#[test]
fn golden_diff_e_and_weights() {
    let (e, w) = run_ours(&golden("counts.tsv"), &[]);
    let e_exp = std::fs::read_to_string(golden("E.expected.tsv")).unwrap();
    let w_exp = std::fs::read_to_string(golden("weights.expected.tsv")).unwrap();
    assert_close(&parse(&e), &parse(&e_exp), "E (golden)");
    assert_close(&parse(&w), &parse(&w_exp), "weights (golden)");
}

/// `--adaptive-span` must reproduce `voom(counts, adaptive.span = TRUE)`; on
/// this 120-gene matrix that span (0.82) diverges 13% from the default.
#[test]
fn golden_diff_adaptive_span() {
    let (_e, w) = run_ours(&golden("counts.tsv"), &["--adaptive-span"]);
    let w_exp = std::fs::read_to_string(golden("weights.adaptive.expected.tsv")).unwrap();
    assert_close(&parse(&w), &parse(&w_exp), "weights (adaptive golden)");
}

/// Degenerate fixture: zero-variance (constant) genes plus an all-zero gene.
/// With only 24 genes the adaptive span saturates to 1.0, so the default and
/// adaptive weights differ by 73% — the divergence that let the span-default
/// bug ship unnoticed. Both branches are pinned against limma here.
#[test]
fn golden_diff_constant_genes() {
    let (e, w) = run_ours(&golden("counts_const.tsv"), &[]);
    let e_exp = std::fs::read_to_string(golden("E_const.expected.tsv")).unwrap();
    let w_exp = std::fs::read_to_string(golden("weights_const.expected.tsv")).unwrap();
    assert_close(&parse(&e), &parse(&e_exp), "E const (golden)");
    assert_close(&parse(&w), &parse(&w_exp), "weights const default (golden)");

    let (_e2, w2) = run_ours(&golden("counts_const.tsv"), &["--adaptive-span"]);
    let w2_exp = std::fs::read_to_string(golden("weights_const.adaptive.expected.tsv")).unwrap();
    assert_close(
        &parse(&w2),
        &parse(&w2_exp),
        "weights const adaptive (golden)",
    );
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

    let (e, w) = run_ours(&golden("counts.tsv"), &[]);
    let e_ref = std::fs::read_to_string(&e_r).unwrap();
    let w_ref = std::fs::read_to_string(&w_r).unwrap();
    let _ = std::fs::remove_file(&e_r);
    let _ = std::fs::remove_file(&w_r);
    assert_close(&parse(&e), &parse(&e_ref), "E (live R)");
    assert_close(&parse(&w), &parse(&w_ref), "weights (live R)");
}

/// A non-finite (NaN/Inf) count is malformed input; limma's `DGEList` rejects it
/// ("NA counts not allowed"). The reader must fail loud, not let NaN reach the
/// lowess sort where a `partial_cmp` unwrap would panic.
#[test]
fn nonfinite_count_rejected() {
    static N: AtomicU64 = AtomicU64::new(0);
    let dir = std::env::temp_dir();
    for bad in ["nan", "inf"] {
        let p = dir.join(format!(
            "voom_bad_{}_{}_{bad}.tsv",
            std::process::id(),
            N.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(
            &p,
            format!("gene\tS1\tS2\tS3\tS4\ng1\t4\t7\t{bad}\t11\ng2\t9\t3\t6\t8\n"),
        )
        .unwrap();
        let out = Command::new(ours()).arg(&p).output().unwrap();
        let _ = std::fs::remove_file(&p);
        assert!(!out.status.success(), "{bad}: expected non-zero exit");
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(
            err.contains("finite"),
            "{bad}: stderr should mention finite, got: {err}"
        );
    }
}
