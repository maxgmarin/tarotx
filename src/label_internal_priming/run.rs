//! Labelling every input alignment.

use anyhow::Result;

use super::label::{Label, Refs, Summary};
use crate::common::input::{for_each_alignment, table_columns, InputFormat, TableCols};
use crate::common::strand::Strand;
use crate::common::table::Sink;

/// Label every alignment of `path` and write the table to `output`. Returns the summary counts.
pub fn run(path: &str, format: InputFormat, cols: &TableCols, refs: &Refs, output: &str) -> Result<Summary> {
    let mut out = Sink::open(output, &table_columns(refs.fields()), "label-internal-priming")?;
    let mut summary = Summary::default();
    for_each_alignment(path, format, cols, |a| {
        let l = a.at.as_ref().map(|p| {
            let strand = if p.strand == Strand::Plus { "+" } else { "-" };
            refs.label(p.chrom, strand, p.start as i64, p.end as i64)
        });
        summary.add(refs, a.at.as_ref().map(|p| p.chrom).zip(l.as_ref()));
        let mut cells = a.lead_cells();
        cells.extend(refs.values(l.as_ref()));
        out.row(&cells)
    })?;
    out.finish()?;
    Ok(summary)
}

// ---------------------------------------------------------------- tagged PAF

/// The PAF line `f` (split on tabs) with this subcommand's tags appended, replacing any pd/td/ip
/// tags it already had: pd:i the poly(A) edge distance, td:i the signed distance to the annotated
/// end, ip:A the flag (Y/N). Not used by the CLI now; kept for a future tagged-PAF output.
#[allow(dead_code)]
pub fn tagged_paf_line(f: &[&str], l: &Label) -> String {
    let mut line: Vec<String> = f
        .iter()
        .enumerate()
        .filter(|(i, t)| *i < 12 || !(t.starts_with("pd:") || t.starts_with("td:") || t.starts_with("ip:")))
        .map(|(_, t)| t.to_string())
        .collect();
    if let Some((d, _)) = l.polya {
        line.push(format!("pd:i:{d}"));
    }
    if let Some((d, _)) = l.tx_end {
        line.push(format!("td:i:{d}"));
    }
    line.push(format!("ip:A:{}", if l.flag { 'Y' } else { 'N' }));
    let mut s = line.join("\t");
    s.push('\n');
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tagged_paf_line_replaces_old_tags() {
        let l = Label { end3: 100, tx_end: None, polya: None, near_polya: false, near_tx_end: false, flag: true };
        let f = ["r1", "55", "5", "55", "+", "chrT", "1000", "50", "100", "50", "50", "60", "tp:A:P", "ip:A:N", "pd:i:4"];
        assert_eq!(tagged_paf_line(&f, &l), "r1\t55\t5\t55\t+\tchrT\t1000\t50\t100\t50\t50\t60\ttp:A:P\tip:A:Y\n");
    }
}
