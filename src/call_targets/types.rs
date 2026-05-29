use std::{collections::HashMap, path::PathBuf};

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub(crate) struct SiteKey {
    pub(crate) reference_sequence_id: usize,
    pub(crate) position: u32,
}

#[derive(Clone, Debug)]
pub(crate) struct SiteCounts {
    pub(crate) per_sample: Vec<[u32; 4]>,
}

#[derive(Clone, Debug)]
pub(crate) struct Interval {
    pub(crate) start: u64,
    pub(crate) end: u64,
}

#[derive(Clone, Debug)]
pub(crate) struct TargetIndex {
    pub(crate) by_ref: HashMap<usize, Vec<Interval>>,
}

pub(crate) struct PreparedCallTargets {
    pub(crate) prepared_reference: PathBuf,
    pub(crate) inputs: Vec<PathBuf>,
    pub(crate) ref_names: Vec<String>,
    pub(crate) targets: TargetIndex,
    pub(crate) sample_names: Vec<String>,
    pub(crate) rg_to_sm: HashMap<String, String>,
    pub(crate) sample_index: HashMap<String, usize>,
}

pub(crate) fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::base_index;

    #[test]
    fn base_index_acgt_maps_to_0123() {
        assert_eq!(base_index(b'A'), Some(0));
        assert_eq!(base_index(b'C'), Some(1));
        assert_eq!(base_index(b'G'), Some(2));
        assert_eq!(base_index(b'T'), Some(3));
    }

    #[test]
    fn base_index_is_case_insensitive() {
        assert_eq!(base_index(b'a'), Some(0));
        assert_eq!(base_index(b'c'), Some(1));
        assert_eq!(base_index(b'g'), Some(2));
        assert_eq!(base_index(b't'), Some(3));
    }

    #[test]
    fn base_index_rejects_non_acgt() {
        assert_eq!(base_index(b'N'), None);
        assert_eq!(base_index(b'n'), None);
        assert_eq!(base_index(b'X'), None);
        assert_eq!(base_index(b'-'), None);
    }
}
