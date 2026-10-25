//! A strand on the reference.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    /// "+" or "-"; anything else (e.g. "*" for an unmapped PAF line) is None.
    pub fn parse(s: &str) -> Option<Strand> {
        match s {
            "+" => Some(Strand::Plus),
            "-" => Some(Strand::Minus),
            _ => None,
        }
    }
}
