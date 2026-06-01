#![allow(dead_code)]

use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result, bail};

use super::types::{SiteKey, TargetIndex};

const MAX_TARGET_SITES: usize = (u32::MAX >> 2) as usize;

#[repr(C)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "wgpu", derive(bytemuck::Pod, bytemuck::Zeroable))]
pub(crate) struct Observation {
    pub(crate) site_idx: u32,
    pub(crate) sample_base: u32,
}

impl Observation {
    pub(crate) fn new(site_idx: u32, sample_idx: u32, base_idx: u32) -> Result<Self> {
        if base_idx > 3 {
            bail!("base index {} exceeds packed observation range", base_idx);
        }

        Ok(Self {
            site_idx,
            sample_base: (sample_idx << 2) | base_idx,
        })
    }

    pub(crate) fn sample_idx(self) -> u32 {
        self.sample_base >> 2
    }

    pub(crate) fn base_idx(self) -> u32 {
        self.sample_base & 3
    }
}

#[derive(Clone, Debug)]
pub(crate) struct TargetSiteMap {
    lookup: HashMap<u64, u32>,
    pub(crate) site_keys: Vec<SiteKey>,
    target_refs: HashSet<usize>,
}

impl TargetSiteMap {
    pub(crate) fn len(&self) -> usize {
        self.site_keys.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.site_keys.is_empty()
    }

    pub(crate) fn contains_ref(&self, reference_sequence_id: usize) -> bool {
        self.target_refs.contains(&reference_sequence_id)
    }

    pub(crate) fn site_index(&self, reference_sequence_id: usize, position: u32) -> Option<u32> {
        self.lookup
            .get(&pack_site_lookup_key(reference_sequence_id, position))
            .copied()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct TargetFrontierIndex {
    ref_ranges: HashMap<usize, (usize, usize)>,
    site_count_u32: u32,
}

impl TargetFrontierIndex {
    pub(crate) fn frontier_for_record(
        &self,
        site_keys: &[SiteKey],
        ref_id: usize,
        pos1: u32,
    ) -> u32 {
        if let Some((start, end)) = self.ref_ranges.get(&ref_id) {
            let local = site_keys[*start..*end].partition_point(|key| key.position < pos1);
            return (start + local) as u32;
        }

        self.site_count_u32
    }
}

pub(crate) fn build_target_site_map(targets: &TargetIndex) -> Result<TargetSiteMap> {
    let estimated_sites = estimate_target_site_count(targets)?;
    if estimated_sites > MAX_TARGET_SITES {
        bail!(
            "target site count {} exceeds GPU-packed limit {}",
            estimated_sites,
            MAX_TARGET_SITES
        );
    }

    let mut ref_ids = targets.by_ref.keys().copied().collect::<Vec<_>>();
    ref_ids.sort_unstable();

    let mut lookup = HashMap::with_capacity(estimated_sites);
    let mut site_keys = Vec::with_capacity(estimated_sites);
    let mut target_refs = HashSet::new();

    for reference_sequence_id in ref_ids {
        if reference_sequence_id > u32::MAX as usize {
            bail!(
                "reference sequence id {} exceeds GPU-packed limit",
                reference_sequence_id
            );
        }

        target_refs.insert(reference_sequence_id);
        let intervals = targets
            .by_ref
            .get(&reference_sequence_id)
            .context("target reference id missing during GPU site map build")?;
        let mut intervals = intervals.iter().collect::<Vec<_>>();
        intervals.sort_unstable_by_key(|interval| (interval.start, interval.end));

        for interval in intervals {
            for pos0 in interval.start..interval.end {
                let pos1 = pos0 + 1;
                if pos1 > u32::MAX as u64 {
                    bail!(
                        "target position {} exceeds GPU-packed position limit on ref {}",
                        pos1,
                        reference_sequence_id
                    );
                }

                let site_idx = u32::try_from(site_keys.len())
                    .context("target site count exceeds u32 range")?;
                let key = SiteKey {
                    reference_sequence_id,
                    position: pos1 as u32,
                };
                lookup.insert(
                    pack_site_lookup_key(reference_sequence_id, key.position),
                    site_idx,
                );
                site_keys.push(key);
            }
        }
    }

    Ok(TargetSiteMap {
        lookup,
        site_keys,
        target_refs,
    })
}

pub(crate) fn build_target_frontier_index(
    site_keys: &[SiteKey],
    _ref_count: usize,
) -> Result<TargetFrontierIndex> {
    let site_count_u32 =
        u32::try_from(site_keys.len()).context("target site count exceeds u32 range")?;

    let mut ref_ranges = HashMap::new();
    let mut i = 0usize;
    while i < site_keys.len() {
        let reference_sequence_id = site_keys[i].reference_sequence_id;
        let start = i;
        while i < site_keys.len() && site_keys[i].reference_sequence_id == reference_sequence_id {
            i += 1;
        }
        ref_ranges.insert(reference_sequence_id, (start, i));
    }

    Ok(TargetFrontierIndex {
        ref_ranges,
        site_count_u32,
    })
}

fn estimate_target_site_count(targets: &TargetIndex) -> Result<usize> {
    let mut total = 0usize;
    for intervals in targets.by_ref.values() {
        for interval in intervals {
            let len_u64 = interval.end.saturating_sub(interval.start);
            let len = usize::try_from(len_u64)
                .context("target interval length exceeds addressable memory space")?;
            total = total
                .checked_add(len)
                .context("target site count overflow during GPU estimate")?;
            if total > MAX_TARGET_SITES {
                return Ok(total);
            }
        }
    }
    Ok(total)
}

fn pack_site_lookup_key(reference_sequence_id: usize, position: u32) -> u64 {
    ((reference_sequence_id as u64) << 32) | position as u64
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use anyhow::Result;

    use super::*;
    use crate::call_targets::types::Interval;

    fn target_index(entries: Vec<(usize, Vec<Interval>)>) -> TargetIndex {
        TargetIndex {
            by_ref: entries.into_iter().collect::<HashMap<_, _>>(),
        }
    }

    #[test]
    fn observation_packs_sample_and_base() -> Result<()> {
        let observation = Observation::new(17, 3, 2)?;

        assert_eq!(observation.site_idx, 17);
        assert_eq!(observation.sample_base, 14);
        assert_eq!(observation.sample_idx(), 3);
        assert_eq!(observation.base_idx(), 2);
        Ok(())
    }

    #[test]
    fn observation_rejects_invalid_base_index() {
        let err = Observation::new(0, 0, 4).unwrap_err();
        assert!(err.to_string().contains("base index"));
    }

    #[test]
    fn build_target_site_map_expands_sorted_one_based_sites() -> Result<()> {
        let targets = target_index(vec![
            (
                2,
                vec![Interval { start: 4, end: 6 }, Interval { start: 1, end: 3 }],
            ),
            (1, vec![Interval { start: 9, end: 11 }]),
        ]);

        let map = build_target_site_map(&targets)?;

        assert_eq!(
            map.site_keys,
            vec![
                SiteKey {
                    reference_sequence_id: 1,
                    position: 10,
                },
                SiteKey {
                    reference_sequence_id: 1,
                    position: 11,
                },
                SiteKey {
                    reference_sequence_id: 2,
                    position: 2,
                },
                SiteKey {
                    reference_sequence_id: 2,
                    position: 3,
                },
                SiteKey {
                    reference_sequence_id: 2,
                    position: 5,
                },
                SiteKey {
                    reference_sequence_id: 2,
                    position: 6,
                },
            ]
        );
        assert_eq!(map.site_index(1, 10), Some(0));
        assert_eq!(map.site_index(2, 6), Some(5));
        assert_eq!(map.site_index(2, 4), None);
        assert!(map.contains_ref(2));
        assert!(!map.contains_ref(3));
        assert_eq!(map.len(), 6);
        assert!(!map.is_empty());
        Ok(())
    }

    #[test]
    fn build_target_site_map_rejects_too_many_sites_without_expanding() {
        let targets = target_index(vec![(
            0,
            vec![Interval {
                start: 0,
                end: MAX_TARGET_SITES as u64 + 1,
            }],
        )]);

        let err = build_target_site_map(&targets).unwrap_err();
        assert!(err.to_string().contains("exceeds GPU-packed limit"));
    }

    #[test]
    fn target_frontier_finds_first_site_at_or_after_record_position() -> Result<()> {
        let targets = target_index(vec![
            (0, vec![Interval { start: 9, end: 12 }]),
            (2, vec![Interval { start: 4, end: 6 }]),
        ]);
        let map = build_target_site_map(&targets)?;
        let frontier = build_target_frontier_index(&map.site_keys, 3)?;

        assert_eq!(frontier.frontier_for_record(&map.site_keys, 0, 1), 0);
        assert_eq!(frontier.frontier_for_record(&map.site_keys, 0, 10), 0);
        assert_eq!(frontier.frontier_for_record(&map.site_keys, 0, 11), 1);
        assert_eq!(frontier.frontier_for_record(&map.site_keys, 0, 13), 3);
        assert_eq!(frontier.frontier_for_record(&map.site_keys, 1, 1), 5);
        assert_eq!(frontier.frontier_for_record(&map.site_keys, 2, 6), 4);
        Ok(())
    }
}
