use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_voom::{read_matrix, voom, write_matrix};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-voom", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    /// log2-CPM (E) matrix destination; "-" is stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    /// Precision-weights matrix destination. Omit to emit only E.
    #[arg(short = 'w', long)]
    weights: Option<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let m = read_matrix(&self.counts)?;
        let v = voom(&m)?;

        let mut e_out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_matrix(&m.header, &m.genes, &v.e, &mut e_out)?;
        drop(e_out);

        if let Some(wpath) = &self.weights {
            let mut w_out = std::fs::File::create(wpath).map_err(RsomicsError::Io)?;
            write_matrix(&m.header, &m.genes, &v.weights, &mut w_out)?;
        }

        if !self.common.quiet {
            eprintln!(
                "{} genes x {} samples voom-transformed",
                m.genes.len(),
                m.samples
            );
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "voom log2-CPM transform with mean-variance precision weights.",
    origin: Some(Origin {
        upstream: "limma voom",
        upstream_license: "GPL (>=2)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1186/gb-2014-15-2-r29"),
    }),
    usage_lines: &["<counts.tsv> [-o E.tsv] [-w weights.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "log2-CPM (E) matrix destination; \"-\" is stdout.",
                why_default: Some("E is the primary transformed expression matrix."),
            },
            FlagSpec {
                short: Some('w'),
                long: "weights",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Precision-weights matrix destination; omit to emit only E.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Write E to stdout and weights to a file",
            command: "rsomics-voom counts.tsv -w weights.tsv > E.tsv",
        },
        Example {
            description: "Write both matrices to files",
            command: "rsomics-voom counts.tsv -o E.tsv -w weights.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
