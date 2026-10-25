//! tarotx: a toolkit for processing long-read RNA-seq alignments.
//!
//! Every subcommand assumes stranded reads (oriented 5'->3', so the alignment strand is the
//! transcript strand) whose adapters and poly(A) tails have been trimmed (so the read ends where
//! the transcript does); see the README's "Input reads" section. Unstranded or untrimmed reads
//! aren't supported yet.
//!
//! Each subcommand lives in its own module, which holds its command-line arguments (`Args`) and
//! its entry point (`run`). Code shared by the subcommands is in `common`. To add a subcommand,
//! add a module with `Args` and `run`, and a variant to `Command`.

mod common;
mod label_by_gene;
mod label_internal_priming;

use std::process::ExitCode;

use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(
    name = "tarotx",
    version,
    about = "A toolkit for processing stranded, adapter- and poly(A)-trimmed long-read RNA-seq alignments.",
    long_about = None,
    max_term_width = 100,
    propagate_version = true
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    LabelByGene(label_by_gene::Args),
    LabelInternalPriming(label_internal_priming::Args),
}

impl Command {
    /// The subcommand's name, as typed on the command line.
    fn name(&self) -> &'static str {
        match self {
            Command::LabelByGene(_) => "label-by-gene",
            Command::LabelInternalPriming(_) => "label-internal-priming",
        }
    }
}

fn main() -> ExitCode {
    let command = Cli::parse().command;
    let name = command.name();
    let result = match command {
        Command::LabelByGene(args) => label_by_gene::run(args),
        Command::LabelInternalPriming(args) => label_internal_priming::run(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        // A closed pipe (e.g. `| head`) is a normal end, not an error.
        Err(e) if e.chain().any(|c| {
            c.downcast_ref::<std::io::Error>().is_some_and(|io| io.kind() == std::io::ErrorKind::BrokenPipe)
        }) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("tarotx {name}: error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
