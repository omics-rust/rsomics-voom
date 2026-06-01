use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_voom(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-voom");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let counts = manifest.join("tests/golden/counts.tsv");
    c.bench_function("rsomics-voom golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([counts.to_str().unwrap(), "-o", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_voom);
criterion_main!(benches);
