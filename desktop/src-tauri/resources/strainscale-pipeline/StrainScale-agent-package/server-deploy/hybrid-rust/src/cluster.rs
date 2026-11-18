//! Clustering: greedy rep picking (pick-reps semantics) + UPGMA + threshold cut.

use rayon::prelude::*;

/// Greedy clustering over strains sorted by tsv distance.
/// Returns (representatives[(name, tsv_dist)], clusters: rep_idx -> member indices).
pub fn greedy_cluster(
    names: &[String],
    dists: &[usize],
    seqs: &[Vec<u8>],
    cutoff: usize,
) -> (Vec<(String, usize)>, Vec<Vec<usize>>) {
    let n = names.len();
    let mut remaining: Vec<usize> = (0..n).collect();
    let mut reps = Vec::new();
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    while !remaining.is_empty() {
        let rep_idx = remaining[0];
        let rep_seq = &seqs[rep_idx];
        reps.push((names[rep_idx].clone(), dists[rep_idx]));

        let dists_to_rep: Vec<(usize, usize)> = remaining
            .par_iter()
            .map(|&i| (i, super::fasta::snp_dist(rep_seq, &seqs[i])))
            .collect();

        let mut assigned = Vec::new();
        let mut keep = Vec::new();
        for (i, d) in dists_to_rep {
            if d <= cutoff {
                assigned.push(i);
            } else {
                keep.push(i);
            }
        }
        clusters.push(assigned);
        remaining = keep;
    }
    (reps, clusters)
}

fn avg_pairwise(a: &[usize], b: &[usize], dist: &[Vec<f64>]) -> f64 {
    let mut s = 0.0;
    let mut cnt = 0;
    for &m1 in a {
        for &m2 in b {
            s += if m1 < m2 { dist[m1][m2] } else { dist[m2][m1] };
            cnt += 1;
        }
    }
    if cnt == 0 {
        0.0
    } else {
        s / cnt as f64
    }
}

/// UPGMA (average-linkage) hierarchical clustering.
/// Returns linkage rows (leaf_a, leaf_b, merge_height) in merge order (n-1 rows).
/// Uses the Lance-Williams average-linkage update: O(n^3) total.
pub fn upgma(dist: &[Vec<f64>]) -> Vec<(usize, usize, f64)> {
    let n = dist.len();
    if n <= 1 {
        return Vec::new();
    }
    // clusters: members (leaf indices) + current size; cd: full m x m average distances
    let mut members: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let mut sizes: Vec<usize> = vec![1; n];
    // cd[i][j] = average distance between cluster i and j
    let mut cd = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dist[i][j];
            cd[i][j] = d;
            cd[j][i] = d;
        }
    }
    let mut linkage = Vec::with_capacity(n - 1);

    while members.len() > 1 {
        let m = members.len();
        let mut best = f64::INFINITY;
        let (mut bi, mut bj) = (0usize, 1usize);
        for i in 0..m {
            for j in (i + 1)..m {
                if cd[i][j] < best {
                    best = cd[i][j];
                    bi = i;
                    bj = j;
                }
            }
        }
        linkage.push((members[bi][0], members[bj][0], best));

        // merge bj into bi (remove bj afterwards so bi stays valid)
        let sz_a = sizes[bi] as f64;
        let sz_b = sizes[bj] as f64;
        let mut merged = std::mem::take(&mut members[bi]);
        merged.extend_from_slice(&members[bj]);
        members[bi] = merged;
        sizes[bi] = sz_a as usize + sz_b as usize;

        // Lance-Williams average linkage: d(new,k) = (n_a*d(a,k) + n_b*d(b,k)) / (n_a+n_b)
        // new_row indexed by OLD cluster indices (0..m, before removal)
        let mut new_row = vec![0.0f64; m];
        for k in 0..m {
            if k == bi || k == bj {
                continue;
            }
            new_row[k] = (sz_a * cd[bi][k] + sz_b * cd[bj][k]) / (sz_a + sz_b);
        }
        // merge bj into bi, then remove bj (bi < bj guaranteed by search loop)
        members.remove(bj);
        sizes.remove(bj);
        cd.remove(bj);
        for row in cd.iter_mut() {
            row.remove(bj);
        }
        // patch row/col bi with new distances (map new index -> old index)
        let m2 = members.len();
        for k2 in 0..m2 {
            if k2 == bi {
                continue;
            }
            let k_old = if k2 < bj { k2 } else { k2 + 1 };
            let d = new_row[k_old];
            cd[bi][k2] = d;
            cd[k2][bi] = d;
        }
    }
    linkage
}

/// Cut a UPGMA linkage at `threshold` (distance criterion): two leaves belong to
/// the same cluster iff the cophenetic distance at which they join is <= threshold
/// (scipy fcluster criterion="distance" semantics).
pub fn cut_by_threshold(
    n: usize,
    linkage: &[(usize, usize, f64)],
    threshold: f64,
) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = x;
        while parent[c] != r {
            let next = parent[c];
            parent[c] = r;
            c = next;
        }
        r
    }
    fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[rb] = ra;
        }
    }
    for &(a, b, h) in linkage {
        if h <= threshold {
            union(&mut parent, a, b);
        }
    }
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut map = std::collections::HashMap::new();
    for leaf in 0..n {
        let r = find(&mut parent, leaf);
        let entry = map.entry(r).or_insert_with(|| {
            clusters.push(Vec::new());
            clusters.len() - 1
        });
        clusters[*entry].push(leaf);
    }
    clusters
}
