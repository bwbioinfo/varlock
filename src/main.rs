use std::{
    fs::{File, OpenOptions},
    io::{self, Read, Write},
    num::NonZeroUsize,
    path::PathBuf,
    thread,
};

mod call_targets;
mod fasta_prep;

use anyhow::{Context, Result};
use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};

fn parse_fraction(value: &str) -> std::result::Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format!("invalid float: {value}"))?;
    if (0.0..=1.0).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!("expected a value in [0.0, 1.0], got {parsed}"))
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "varlock",
    about = "Variant calling from BAM alignments",
    version,
    propagate_version = true
)]
struct Cli {
    /// Increase logging verbosity (-v, -vv for more)
    #[arg(short, long, action = ArgAction::Count)]
    verbose: u8,

    /// Number of worker threads (default: available parallelism)
    #[arg(short = 't', long, value_name = "THREADS")]
    threads: Option<NonZeroUsize>,

    /// Mirror stderr logs to a file
    #[arg(long = "log-file", value_name = "FILE", global = true)]
    log_file: Option<PathBuf>,

    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Call simple SNVs from BAMs against target regions (bgzipped VCF + index)
    CallTargets(CallTargetsArgs),
}

#[derive(Args, Debug, Clone)]
pub struct CallTargetsArgs {
    /// BAM input paths (files or directories; directories are searched recursively)
    #[arg(
        short = 'i',
        long = "input",
        value_name = "PATH",
        required_unless_present = "bamlist",
        conflicts_with = "bamlist"
    )]
    pub inputs: Vec<PathBuf>,

    /// File with one BAM path per line
    #[arg(
        long = "bamlist",
        value_name = "FILE",
        required_unless_present = "inputs"
    )]
    pub bamlist: Option<PathBuf>,

    /// Reference FASTA path (plain or bgzipped; index created automatically if missing)
    #[arg(short = 'r', long = "reference", value_name = "FASTA")]
    pub reference: PathBuf,

    /// BED file of target intervals (0-based, half-open)
    #[arg(short = 'T', long = "targets", value_name = "BED")]
    pub targets: PathBuf,

    /// Output VCF.gz path
    #[arg(short = 'o', long = "output", value_name = "VCF_GZ")]
    pub output: PathBuf,

    /// Read-group to sample mapping file (RG<tab>SM)
    #[arg(long = "rg-map", value_name = "FILE")]
    pub rg_map: Option<PathBuf>,

    /// Index type for VCF.gz output
    #[arg(long = "index-type", value_name = "TYPE", value_enum, default_value_t = IndexType::Csi)]
    pub index_type: IndexType,

    /// Minimum mapping quality (0-60)
    #[arg(
        long = "min-mapq",
        value_name = "MAPQ",
        default_value_t = 20u8,
        value_parser = clap::value_parser!(u8).range(0..=60)
    )]
    pub min_mapq: u8,

    /// Minimum base quality (0-60)
    #[arg(
        long = "min-baseq",
        value_name = "BASEQ",
        default_value_t = 20u8,
        value_parser = clap::value_parser!(u8).range(0..=60)
    )]
    pub min_baseq: u8,

    /// Minimum total alt count across all samples
    #[arg(long = "min-alt-count", value_name = "COUNT", default_value_t = 1u32)]
    pub min_alt_count: u32,

    /// Minimum alt allele fraction across all samples (0.0-1.0)
    #[arg(
        long = "min-alt-fraction",
        value_name = "FRACTION",
        default_value_t = 0.0,
        value_parser = parse_fraction
    )]
    pub min_alt_fraction: f64,

    /// Maximum read depth per sample at a site
    #[arg(long = "max-depth", value_name = "DP", default_value_t = 100_000u32)]
    pub max_depth: u32,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum IndexType {
    Csi,
    Tbi,
}

pub struct ExecutionContext {
    pub verbose: u8,
    pub threads: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    #[cfg(unix)]
    let _stderr_tee = cli.log_file.as_deref().map(init_stderr_tee).transpose()?;
    let threads = cli
        .threads
        .map(NonZeroUsize::get)
        .or_else(|| thread::available_parallelism().ok().map(NonZeroUsize::get))
        .context("unable to determine available parallelism")?;

    let ctx = ExecutionContext {
        verbose: cli.verbose,
        threads,
    };

    match cli.command {
        Command::CallTargets(args) => call_targets::run(args, &ctx),
    }
}

pub(crate) fn log_verbose(ctx: &ExecutionContext, message: String) {
    if ctx.verbose > 0 {
        eprintln!("{}", message);
    }
}

#[cfg(unix)]
struct StderrTeeGuard {
    saved_stderr_fd: i32,
    join_handle: Option<thread::JoinHandle<()>>,
}

#[cfg(unix)]
impl Drop for StderrTeeGuard {
    fn drop(&mut self) {
        unsafe {
            let _ = libc::dup2(self.saved_stderr_fd, libc::STDERR_FILENO);
            let _ = libc::close(self.saved_stderr_fd);
        }
        if let Some(handle) = self.join_handle.take() {
            let _ = handle.join();
        }
    }
}

#[cfg(unix)]
fn init_stderr_tee(path: &std::path::Path) -> Result<StderrTeeGuard> {
    use std::os::fd::FromRawFd;

    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create log directory {}", parent.display()))?;
    }

    let log_file = OpenOptions::new()
        .create(true)
        .truncate(true)
        .write(true)
        .open(path)
        .with_context(|| format!("failed to open log file {}", path.display()))?;

    let mut fds = [0i32; 2];
    if unsafe { libc::pipe(fds.as_mut_ptr()) } != 0 {
        return Err(io::Error::last_os_error()).context("failed to create stderr pipe");
    }
    let read_fd = fds[0];
    let write_fd = fds[1];

    let saved_stderr_fd = unsafe { libc::dup(libc::STDERR_FILENO) };
    if saved_stderr_fd < 0 {
        unsafe {
            let _ = libc::close(read_fd);
            let _ = libc::close(write_fd);
        }
        return Err(io::Error::last_os_error()).context("failed to duplicate stderr fd");
    }

    let thread_stderr_fd = unsafe { libc::dup(saved_stderr_fd) };
    if thread_stderr_fd < 0 {
        unsafe {
            let _ = libc::close(read_fd);
            let _ = libc::close(write_fd);
            let _ = libc::close(saved_stderr_fd);
        }
        return Err(io::Error::last_os_error()).context("failed to duplicate stderr fd");
    }

    if unsafe { libc::dup2(write_fd, libc::STDERR_FILENO) } < 0 {
        unsafe {
            let _ = libc::close(read_fd);
            let _ = libc::close(write_fd);
            let _ = libc::close(saved_stderr_fd);
            let _ = libc::close(thread_stderr_fd);
        }
        return Err(io::Error::last_os_error()).context("failed to redirect stderr");
    }
    unsafe {
        let _ = libc::close(write_fd);
    }

    let join_handle = thread::spawn(move || {
        let mut reader = unsafe { File::from_raw_fd(read_fd) };
        let mut stderr_out = unsafe { File::from_raw_fd(thread_stderr_fd) };
        let mut log_out = log_file;
        let mut buf = [0u8; 8192];

        while let Ok(n) = reader.read(&mut buf) {
            if n == 0 {
                break;
            }
            let chunk = &buf[..n];
            let _ = stderr_out.write_all(chunk);
            let _ = stderr_out.flush();
            let _ = log_out.write_all(chunk);
            let _ = log_out.flush();
        }
    });

    Ok(StderrTeeGuard {
        saved_stderr_fd,
        join_handle: Some(join_handle),
    })
}
