//! Text input and output with transparent gzip/BGZF.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

const BUF: usize = 1 << 20;

/// Open a text file, or stdin for "-". gzip and BGZF input is detected from the magic bytes and
/// decoded in full (all members), so files from bgzip, pigz or `cat a.gz b.gz` all work.
pub fn open_input(path: &str) -> Result<Box<dyn BufRead>> {
    let raw: Box<dyn Read> = if path == "-" {
        Box::new(io::stdin().lock())
    } else {
        Box::new(File::open(path).with_context(|| format!("cannot open '{path}'"))?)
    };
    let mut r = BufReader::with_capacity(BUF, raw);
    let gz = r.fill_buf().with_context(|| format!("cannot read '{path}'"))?.starts_with(&[0x1f, 0x8b]);
    if gz {
        Ok(Box::new(BufReader::with_capacity(BUF, MultiGzDecoder::new(r))))
    } else {
        Ok(Box::new(r))
    }
}

/// A text output. `finish` must be called: it writes the gzip trailer and reports any error, which
/// dropping the writer would silently ignore.
pub enum Output {
    Plain(BufWriter<Box<dyn Write>>),
    Gz(Box<BufWriter<GzEncoder<File>>>),
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Output::Plain(w) => w.write(buf),
            Output::Gz(w) => w.write(buf),
        }
    }
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Output::Plain(w) => w.write_all(buf),
            Output::Gz(w) => w.write_all(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Output::Plain(w) => w.flush(),
            Output::Gz(w) => w.flush(),
        }
    }
}

impl Output {
    pub fn finish(self) -> io::Result<()> {
        match self {
            Output::Plain(mut w) => w.flush(),
            Output::Gz(w) => (*w).into_inner().map_err(|e| e.into_error())?.finish().map(drop),
        }
    }
}

/// Create a text output: stdout for "-", gzip for names ending in ".gz", else plain.
pub fn open_output(path: &str) -> Result<Output> {
    if path == "-" {
        return Ok(Output::Plain(BufWriter::with_capacity(BUF, Box::new(io::stdout().lock()))));
    }
    let f = File::create(path).with_context(|| format!("cannot create '{path}'"))?;
    if path.ends_with(".gz") {
        Ok(Output::Gz(Box::new(BufWriter::with_capacity(BUF, GzEncoder::new(f, Compression::default())))))
    } else {
        Ok(Output::Plain(BufWriter::with_capacity(BUF, Box::new(f))))
    }
}

/// Call `f` on each line (without the line terminator), with its 1-based line number.
pub fn for_each_line(path: &str, mut f: impl FnMut(usize, &str) -> Result<()>) -> Result<()> {
    let mut r = open_input(path)?;
    let mut line = String::new();
    let mut n = 0;
    loop {
        line.clear();
        if r.read_line(&mut line).with_context(|| format!("cannot read '{path}'"))? == 0 {
            return Ok(());
        }
        n += 1;
        f(n, line.trim_end_matches(['\n', '\r']))?;
    }
}
