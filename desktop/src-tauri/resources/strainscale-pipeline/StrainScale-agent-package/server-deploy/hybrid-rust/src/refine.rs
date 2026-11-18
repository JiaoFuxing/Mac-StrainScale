//! Specific-site finding, real-sample matching, recursive refinement and pruning.

use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::cluster::{cut_by_threshold, upgma};

pub const VALID: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Map a base to 0..3 (A/C/G/T), None if not a valid base.
fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Compute specific sites and match results for a candidate rep list.
/// `rep_idx[k]` = index into `seqs` for rep k; `poc[row]` = address of row.
/// Returns (specific_sites (row, rep_name, base), match rows (rep, total, matched, rate)).
pub fn compute_specific_match(
    rep_names: &[String],
    rep_idx: &[usize],
    poc: &[u32],
    seqs: &[Vec<u8>],
    filter: &HashMap<u32, String>,
) -> (Vec<(usize, String, u8)>, Vec<(String, usize, usize, f64)>) {
    let k = rep_names.len();
    let l = poc.len();
    let rows: Vec<(Vec<(usize, usize, u8)>, Vec<(usize, usize)>)> = (0..l)
        .into_par_iter()
        .map(|row| {
            let mut bases = Vec::with_capacity(k);
            for &ri in rep_idx {
                bases.push(seqs[ri].get(row).copied().unwrap_or(crate::fasta::IGNORE));
            }
            // count how many reps carry each base at this site (A/C/G/T)
            let mut cnt = [0usize; 4];
            for &b in &bases {
                if let Some(i) = base_idx(b) {
                    cnt[i] += 1;
                }
            }
            // per rep: (matched, total_specific)
            let mut counters = vec![(0usize, 0usize); k];
            let mut specific_here: Vec<(usize, usize, u8)> = Vec::new();
            for (i, &b) in bases.iter().enumerate() {
                let bi = match base_idx(b) {
                    Some(x) => x,
                    None => continue,
                };
                if cnt[bi] != 1 {
                    continue;
                }
                counters[i].1 += 1;
                let addr = poc[row];
                let matched = filter
                    .get(&addr)
                    .map(|s| base_in_counts(s, b))
                    .unwrap_or(false);
                if matched {
                    counters[i].0 += 1;
                }
                specific_here.push((row, i, b));
            }
            (specific_here, counters)
        })
        .collect();

    let mut specific = Vec::new();
    let mut totals = vec![0usize; k];
    let mut matched = vec![0usize; k];
    for (sh, counters) in rows {
        for (row, i, b) in sh {
            specific.push((row, rep_names[i].clone(), b));
        }
        for i in 0..k {
            totals[i] += counters[i].1;
            matched[i] += counters[i].0;
        }
    }

    let match_rows: Vec<(String, usize, usize, f64)> = (0..k)
        .map(|i| {
            let rate = if totals[i] > 0 {
                matched[i] as f64 / totals[i] as f64
            } else {
                0.0
            };
            (rep_names[i].clone(), totals[i], matched[i], rate)
        })
        .collect();
    (specific, match_rows)
}

/// Check whether `base` appears with count > 0 in "A:51,C:2" style string.
fn base_in_counts(counts: &str, base: u8) -> bool {
    for pair in counts.split(',') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        let mut sub = pair.split(':');
        let b = sub.next().unwrap_or("");
        let cnt = sub.next().unwrap_or("0");
        if b.as_bytes().len() == 1 && b.as_bytes()[0].to_ascii_uppercase() == base {
            if let Ok(n) = cnt.parse::<f64>() {
                if n > 0.0 {
                    return true;
                }
            }
        }
    }
    false
}

/// Cluster a set of member indices using their sub-distance matrix.
/// Prefers the finest threshold where every rep has >= min_sites specific sites
/// (checked against current_final_reps context). Returns (clusters_of_member_idx, threshold, rep_idx_list).
#[allow(clippy::too_many_arguments)]
pub fn cluster_sub_for_sites(
    members: &[usize], // indices into rep_idx space (global rep ordering)
    global_dist: &[Vec<f64>],
    rep_names: &[String],
    rep_tsv_dist: &[usize],
    rep_idx: &[usize],
    poc: &[u32],
    seqs: &[Vec<u8>],
    filter: &HashMap<u32, String>,
    current_final: &[String],
    min_sites: usize,
    max_clusters: usize,
) -> (Vec<Vec<usize>>, f64, Vec<usize>) {
    let m = members.len();
    if m <= 1 {
        return (vec![members.to_vec()], 0.0, members.to_vec());
    }
    // sub distance matrix
    let mut sub = vec![vec![0.0f64; m]; m];
    for i in 0..m {
        for j in (i + 1)..m {
            let d = global_dist[members[i]][members[j]];
            sub[i][j] = d;
            sub[j][i] = d;
        }
    }
    let sub_ref = &sub;
    let ut: Vec<f64> = (0..m)
        .flat_map(move |i| ((i + 1)..m).map(move |j| sub_ref[i][j]))
        .collect();
    if ut.is_empty() {
        return (vec![members.to_vec()], 0.0, members.to_vec());
    }
    let min_d = ut.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_d = ut.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let n_steps = 25usize;
    let mut thresholds: Vec<f64> = if max_d > min_d {
        (0..=n_steps)
            .map(|i| min_d + i as f64 * (max_d - min_d) / n_steps as f64)
            .collect()
    } else {
        vec![min_d]
    };
    thresholds = thresholds.into_iter().map(|x| (x * 10000.0).round() / 10000.0).collect();
    thresholds.sort_by(|a, b| b.partial_cmp(a).unwrap()); // fine -> coarse
    thresholds.dedup();

    let linkage = upgma(&sub);

    let mut best_clusters = vec![members.to_vec()];
    let mut best_threshold = max_d;
    let mut best_reps = members.to_vec();
    let mut best_min_count = -1i64;
    let mut first_no_zero: Option<(Vec<Vec<usize>>, f64, Vec<usize>, i64)> = None;

    for threshold in thresholds {
        let clusters = cut_by_threshold(m, &linkage, threshold);
        if clusters.len() > max_clusters {
            continue;
        }
        let mut reps_local: Vec<usize> = Vec::new();
        for c in &clusters {
            // pick member with smallest tsv dist
            let mut best = c[0];
            for &mm in c {
                if rep_tsv_dist[members[mm]] < rep_tsv_dist[members[best]] {
                    best = mm;
                }
            }
            reps_local.push(members[best]);
        }
        let mut candidate: Vec<String> = current_final.to_vec();
        let mut candidate_idx: Vec<usize> = current_final
            .iter()
            .map(|nm| {
                let p = rep_names.iter().position(|x| x == nm).unwrap_or(0);
                rep_idx[p]
            })
            .collect();
        let mut seen: HashSet<String> = current_final.iter().cloned().collect();
        for &r in &reps_local {
            if seen.insert(rep_names[r].clone()) {
                candidate.push(rep_names[r].clone());
                candidate_idx.push(rep_idx[r]);
            }
        }
        let (specific, _) = compute_specific_match(&candidate, &candidate_idx, poc, seqs, filter);
        let mut counts: HashMap<String, usize> = HashMap::new();
        for (_, name, _) in &specific {
            *counts.entry(name.clone()).or_insert(0) += 1;
        }
        let mut min_count = i64::MAX;
        let mut zero_reps = 0;
        for &r in &reps_local {
            let c = counts.get(&rep_names[r]).copied().unwrap_or(0) as i64;
            if c == 0 {
                zero_reps += 1;
            }
            min_count = min_count.min(c);
        }
        let this = (clusters.clone(), threshold, reps_local.clone(), min_count);
        if zero_reps == 0 && first_no_zero.is_none() {
            first_no_zero = Some(this.clone());
        }
        if zero_reps == 0 && min_count >= min_sites as i64 {
            best_clusters = clusters;
            best_threshold = threshold;
            best_reps = reps_local;
            best_min_count = min_count;
            break;
        }
    }
    if best_min_count < 0 {
        if let Some((c, t, r, mc)) = first_no_zero {
            best_clusters = c;
            best_threshold = t;
            best_reps = r;
            best_min_count = mc;
        }
    }
    // map back: clusters are over 0..m; convert to member indices
    let out_clusters: Vec<Vec<usize>> = best_clusters
        .iter()
        .map(|c| c.iter().map(|&x| members[x]).collect())
        .collect();
    (out_clusters, best_threshold, best_reps)
}
