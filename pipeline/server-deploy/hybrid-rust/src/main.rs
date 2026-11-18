//! hybrid: single binary that replaces
//!   pick-reps (greedy clustering) -> snp-dists (distance matrix)
//!   -> fasta2matrix (base matrix) + address merge -> strain-refine + final_refinement_hybrid.
//! Usage: hybrid <sample_out_dir> <cutoff> <max_dist> <filter_counts> [-t 0.9] [-o outdir]

mod cluster;
mod fasta;
mod refine;

use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use cluster::{cut_by_threshold, greedy_cluster, upgma};
use fasta::{build_name_index, read_fasta, snp_dist};
use refine::compute_specific_match;

const MIN_SITES: usize = 30;
const MAX_CLUSTERS: usize = 50;
const MAX_SUB_ATTEMPTS: usize = 25;
const MAX_PRUNE_ITERS: usize = 200;

struct Ctx {
    out_dir: PathBuf,          // sample output dir (rep_*.tsv, fa_*, poc2.txt)
    cutoff: usize,
    max_dist: usize,
    filter_path: PathBuf,
    threshold: f64,
    out_root: PathBuf,         // where final_* go
    cache: bool,               // use all.fa merged cache (faster on native ext4/NFS, slower on 9P)
}

fn main() {
    let t_all = Instant::now();
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!(
            "Usage: hybrid <sample_out_dir> <cutoff> <max_dist> <filter_counts> [-t 0.9] [-o outdir]"
        );
        std::process::exit(1);
    }
    let out_dir = PathBuf::from(&args[1]);
    let cutoff: usize = args[2].parse().expect("cutoff");
    let max_dist: usize = args[3].parse().expect("max_dist");
    let filter_path = PathBuf::from(&args[4]);
    let mut threshold = 0.9f64;
    let mut out_root = out_dir.clone();
    let mut cache = false;
    let mut i = 5;
    while i < args.len() {
        match args[i].as_str() {
            "-t" => {
                threshold = args[i + 1].parse().expect("threshold");
                i += 2;
            }
            "-o" => {
                out_root = PathBuf::from(&args[i + 1]);
                i += 2;
            }
            "-c" => {
                cache = true;
                i += 1;
            }
            _ => i += 1,
        }
    }
    let ctx = Ctx {
        out_dir,
        cutoff,
        max_dist,
        filter_path,
        threshold,
        out_root,
        cache,
    };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(20).build_global();

    // 1. collect strains from rep_*.tsv (dist <= max_dist), sort by dist
    let t = Instant::now();
    let mut seen: HashMap<String, usize> = HashMap::new();
    if let Ok(rd) = fs::read_dir(&ctx.out_dir) {
        for e in rd.flatten() {
            let p = e.path();
            let fname = p.file_name().unwrap().to_string_lossy();
            if fname.starts_with("rep_") && fname.ends_with(".tsv") {
                for line in BufReader::new(File::open(&p).unwrap()).lines().flatten() {
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() >= 2 {
                        let name = parts[0].to_string();
                        let d: usize = parts[1].parse().unwrap_or(usize::MAX);
                        if d <= ctx.max_dist {
                            seen.entry(name).or_insert(d);
                        }
                    }
                }
            }
        }
    }
    let mut all: Vec<(String, usize)> = seen.into_iter().collect();
    // deterministic: tsv dist asc, then name (pick-reps used HashMap random tie-break)
    all.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
    eprintln!(
        "[1] collected {} strains (dist <= {}) in {:?}",
        all.len(),
        ctx.max_dist,
        t.elapsed()
    );

    // fast path: nothing passed max_dist -> write empty finals and exit
    if all.is_empty() {
        fs::create_dir_all(&ctx.out_root).unwrap();
        let mut w = BufWriter::new(File::create(ctx.out_root.join("final_representatives.txt")).unwrap());
        drop(w);
        let mut w = BufWriter::new(File::create(ctx.out_root.join("final_match_results.tsv")).unwrap());
        writeln!(w, "Strain\tSpecificSites\tMatched\tMatchRate").unwrap();
        let mut w = BufWriter::new(File::create(ctx.out_root.join("final_specific_sites.tsv")).unwrap());
        writeln!(w, "strain\taddress\tbase").unwrap();
        let mut w = BufWriter::new(File::create(ctx.out_root.join("final_clusters.tsv")).unwrap());
        writeln!(w, "strain\tcluster_id").unwrap();
        eprintln!("[x] no strains passed max_dist; empty finals written in {:?}", t_all.elapsed());
        eprintln!("=== hybrid total {:?} ===", t_all.elapsed());
        eprintln!("final reps: 0");
        return;
    }

    // 2. load sequences (parallel per-file default; --cache uses merged all.fa, best on ext4/NFS)
    let t = Instant::now();
    let yc_dir = ctx.out_dir.parent().unwrap().parent().unwrap().to_path_buf();
    let cache_path = ctx.out_dir.join("all.fa");
    let loaded: Vec<(String, usize, Option<Vec<u8>>)> = if ctx.cache && cache_path.exists() {
        // sequential single-file read of merged cache
        let mut by_name: HashMap<String, Vec<u8>> = HashMap::new();
        for (nm, s) in fasta::read_fasta(&cache_path) {
            by_name.insert(nm, s);
        }
        all.par_iter()
            .map(|(name, d)| (name.clone(), *d, by_name.get(name).cloned()))
            .collect()
    } else {
        let idx = build_name_index(&ctx.out_dir, &yc_dir);
        let loaded = all
            .par_iter()
            .map(|(name, d)| {
                let seq = idx
                    .get(name)
                    .and_then(|fp| read_fasta(fp).into_iter().next())
                    .map(|(_, s)| s);
                (name.clone(), *d, seq)
            })
            .collect::<Vec<_>>();
        if ctx.cache {
            // write merged cache for next runs
            let f = File::create(&cache_path).unwrap();
            let mut w = BufWriter::with_capacity(1 << 20, f);
            for (name, _, seq) in &loaded {
                writeln!(w, ">{}", name).unwrap();
                if let Some(s) = seq {
                    for chunk in s.chunks(80) {
                        let line: String = chunk.iter().map(|&b| b as char).collect();
                        writeln!(w, "{}", line).unwrap();
                    }
                }
            }
        }
        loaded
    };
    let seqs: Vec<Vec<u8>> = loaded.iter().map(|x| x.2.clone().unwrap_or_default()).collect();
    let names_all: Vec<String> = loaded.iter().map(|x| x.0.clone()).collect();
    let dists_all: Vec<usize> = loaded.iter().map(|x| x.1).collect();
    eprintln!("[2] loaded {} FASTA sequences in {:?}", seqs.len(), t.elapsed());

    // 3. greedy clustering with cutoff
    let t = Instant::now();
    let (reps, clusters) = greedy_cluster(&names_all, &dists_all, &seqs, ctx.cutoff);
    let rep_seqs: Vec<Vec<u8>> = reps
        .iter()
        .map(|(name, _)| {
            names_all
                .iter()
                .position(|n| n == name)
                .map(|i| seqs[i].clone())
                .unwrap_or_default()
        })
        .collect();
    eprintln!(
        "[3] greedy clustering: {} representatives (cutoff={}) in {:?}",
        reps.len(),
        ctx.cutoff,
        t.elapsed()
    );

    // rep names S1_..Sk and seq indices
    let rep_names: Vec<String> = reps
        .iter()
        .enumerate()
        .map(|(i, (n, _))| format!("S{}_{}", i + 1, n))
        .collect();
    let rep_tsv_dist: Vec<usize> = reps.iter().map(|(_, d)| *d).collect();
    let rep_idx: Vec<usize> = reps
        .iter()
        .map(|(n, _)| names_all.iter().position(|x| x == n).unwrap())
        .collect();

    // 4. write representatives.txt / .fa / cluster_assignments.tsv
    write_reps(&ctx, &reps, &rep_names, &rep_seqs, &names_all, &clusters);

    // 5. distance matrix among reps -> representatives.dist.tsv (upper triangle only, mirrored)
    let t = Instant::now();
    let n_rep = rep_seqs.len();
    let tri: Vec<Vec<usize>> = (0..n_rep)
        .into_par_iter()
        .map(|a| (0..=a).map(|b| snp_dist(&rep_seqs[a], &rep_seqs[b])).collect())
        .collect();
    let mut dist_mat = vec![vec![0usize; n_rep]; n_rep];
    for a in 0..n_rep {
        for b in 0..=a {
            dist_mat[a][b] = tri[a][b];
            dist_mat[b][a] = tri[a][b];
        }
    }
    eprintln!("[5] rep distance matrix {}x{} in {:?}", n_rep, n_rep, t.elapsed());
    let dist_path = ctx.out_dir.join("representatives.dist.tsv");
    {
        let f = File::create(&dist_path).unwrap();
        let mut w = BufWriter::with_capacity(1 << 20, f);
        write!(w, "snp-dists 0.8.2").unwrap();
        for nm in &rep_names {
            write!(w, "\t{}", nm).unwrap();
        }
        writeln!(w).unwrap();
        for a in 0..n_rep {
            write!(w, "{}", rep_names[a]).unwrap();
            for b in 0..n_rep {
                write!(w, "\t{}", dist_mat[a][b]).unwrap();
            }
            writeln!(w).unwrap();
        }
    }
    eprintln!("    wrote {}", dist_path.display());

    // 6. base matrix (address x reps) from poc2.txt -> representatives.matrix.addr.tsv
    let t = Instant::now();
    let poc = read_poc2(&ctx.out_dir.join("poc2.txt"));
    eprintln!("[6] poc2 addresses: {}", poc.len());
    let matrix_path = ctx.out_dir.join("representatives.matrix.addr.tsv");
    {
        let f = File::create(&matrix_path).unwrap();
        let mut w = BufWriter::with_capacity(1 << 22, f);
        write!(w, "address").unwrap();
        for nm in &rep_names {
            write!(w, "\t{}", nm).unwrap();
        }
        writeln!(w).unwrap();
        // parallel row construction
        let rows: Vec<Vec<u8>> = (0..poc.len())
            .into_par_iter()
            .map(|row| {
                let mut v = Vec::with_capacity(n_rep * 2);
                for &ri in &rep_idx {
                    if v.len() > 0 {
                        v.push(b'\t');
                    }
                    v.push(seqs[ri].get(row).copied().unwrap_or(b'-'));
                }
                v
            })
            .collect();
        for (row, v) in rows.iter().enumerate() {
            write!(w, "{}", poc[row]).unwrap();
            w.write_all(v).unwrap();
            writeln!(w).unwrap();
        }
    }
    eprintln!(
        "[6] wrote {} rows matrix in {:?}",
        poc.len(),
        t.elapsed()
    );

    // 7. load filter counts
    let t = Instant::now();
    let filter = read_filter_counts(&ctx.filter_path);
    eprintln!("[7] loaded {} filter addresses in {:?}", filter.len(), t.elapsed());

    // 8. refinement (replaces strain-refine + final_refinement_hybrid)
    let t = Instant::now();
    let final_reps = run_refinement(&ctx, &rep_names, &rep_tsv_dist, &rep_idx, &poc, &seqs, &filter, &dist_mat);
    eprintln!("[8] refinement done in {:?}", t.elapsed());

    eprintln!(
        "=== hybrid total {:?} ===",
        t_all.elapsed()
    );
    eprintln!("final reps: {}", final_reps.len());
}

fn write_reps(
    ctx: &Ctx,
    reps: &[(String, usize)],
    rep_names: &[String],
    rep_seqs: &[Vec<u8>],
    names_all: &[String],
    clusters: &[Vec<usize>],
) {
    let rep_file = ctx.out_dir.join("representatives.txt");
    {
        let f = File::create(&rep_file).unwrap();
        let mut w = BufWriter::new(f);
        for (i, (n, d)) in reps.iter().enumerate() {
            writeln!(w, "S{}\t{}\ttsv_dist={}", i + 1, n, d).unwrap();
        }
    }
    let cat_file = ctx.out_dir.join("representatives.fa");
    {
        let f = File::create(&cat_file).unwrap();
        let mut w = BufWriter::with_capacity(1 << 20, f);
        for (i, nm) in rep_names.iter().enumerate() {
            writeln!(w, ">{}", nm).unwrap();
            for chunk in rep_seqs[i].chunks(80) {
                let s: String = chunk.iter().map(|&b| b as char).collect();
                writeln!(w, "{}", s).unwrap();
            }
        }
    }
    let cls_file = ctx.out_dir.join("cluster_assignments.tsv");
    {
        let f = File::create(&cls_file).unwrap();
        let mut w = BufWriter::new(f);
        writeln!(w, "representative\tmember").unwrap();
        for (ci, members) in clusters.iter().enumerate() {
            let rep_name = &reps[ci].0;
            for &m in members {
                writeln!(w, "{}\t{}", rep_name, names_all[m]).unwrap();
            }
        }
    }
}

fn read_poc2(path: &Path) -> Vec<u32> {
    let mut out = Vec::new();
    if let Ok(f) = File::open(path) {
        for line in BufReader::new(f).lines().flatten() {
            let t = line.trim();
            if !t.is_empty() {
                if let Ok(v) = t.parse::<u32>() {
                    out.push(v);
                }
            }
        }
    }
    out
}

fn read_filter_counts(path: &Path) -> HashMap<u32, String> {
    let mut out = HashMap::new();
    if let Ok(f) = File::open(path) {
        for line in BufReader::new(f).lines().flatten() {
            let mut parts = line.split('\t');
            let addr: u32 = match parts.next().and_then(|s| s.trim().parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            let counts = parts.next().unwrap_or("").to_string();
            out.insert(addr, counts);
        }
    }
    out
}

/// Map rep names (S{i}_name) back to their global seq indices.
fn idx_for(names: &[String], rep_names: &[String], rep_idx: &[usize]) -> Vec<usize> {
    names
        .iter()
        .map(|nm| {
            let p = rep_names.iter().position(|x| x == nm).unwrap_or(0);
            rep_idx[p]
        })
        .collect()
}

/// Full refinement: base clustering -> recursive refinement -> global pruning.
fn run_refinement(
    ctx: &Ctx,
    rep_names: &[String],
    rep_tsv_dist: &[usize],
    rep_idx: &[usize],
    poc: &[u32],
    seqs: &[Vec<u8>],
    filter: &HashMap<u32, String>,
    dist_mat: &[Vec<usize>],
) -> Vec<String> {
    let n = rep_names.len();
    let gdist: Vec<Vec<f64>> = dist_mat
        .iter()
        .map(|r| r.iter().map(|&x| x as f64).collect())
        .collect();
    let mut final_reps: Vec<String> = Vec::new();
    let mut current_final: Vec<String> = Vec::new();

    if n <= 1 {
        // single rep: run specific+match directly
        let (specific, match_rows) =
            compute_specific_match(rep_names, rep_idx, poc, seqs, filter);
        write_final_outputs(ctx, rep_names.to_vec(), &specific, &match_rows);
        return rep_names.to_vec();
    }

    // ---- base clustering: candidate thresholds from min..median..max ----
    let gdist_ref = &gdist;
    let ut: Vec<f64> = (0..n)
        .flat_map(move |i| ((i + 1)..n).map(move |j| gdist_ref[i][j]))
        .collect();
    let min_d = ut.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_d = ut.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let median = {
        let mut v = ut.clone();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if v.len() % 2 == 0 {
            (v[v.len() / 2 - 1] + v[v.len() / 2]) / 2.0
        } else {
            v[v.len() / 2]
        }
    };
    eprintln!(
        "base dist stats: min={:.0} median={:.0} max={:.0}",
        min_d, median, max_d
    );

    let mut thresholds: Vec<f64> = Vec::new();
    let n_lower = 30usize;
    let n_upper = 10usize;
    for i in 0..=n_lower {
        thresholds.push(median - i as f64 * (median - min_d) / n_lower as f64);
    }
    for i in 1..=n_upper {
        thresholds.push(median + i as f64 * (max_d - median) / n_upper as f64);
    }
    let mut th: Vec<f64> = thresholds
        .into_iter()
        .map(|x| (x * 100.0).round() / 100.0)
        .collect();
    th.sort_by(|a, b| a.partial_cmp(b).unwrap());
    th.dedup();
    let thresholds = th;

    let linkage = upgma(&gdist);
    let mut clusters: Vec<Vec<usize>> = vec![(0..n).collect()];
    let mut reps: Vec<usize> = vec![0];
    for tv in &thresholds {
        let c = cut_by_threshold(n, &linkage, *tv);
        if c.len() > MAX_CLUSTERS {
            continue;
        }
        let r: Vec<usize> = c.iter().map(|cl| pick_min_dist_member(cl, rep_tsv_dist)).collect();
        let r_names: Vec<String> = r.iter().map(|&x| rep_names[x].clone()).collect();
        let r_idx: Vec<usize> = r.iter().map(|&x| rep_idx[x]).collect();
        let (specific, _) = compute_specific_match(&r_names, &r_idx, poc, seqs, filter);
        let mut counts: HashMap<String, usize> = HashMap::new();
        for (_, nm, _) in &specific {
            *counts.entry(nm.clone()).or_insert(0) += 1;
        }
        let zero = r.iter().any(|&x| counts.get(&rep_names[x]).copied().unwrap_or(0) == 0);
        eprintln!("  threshold={:.1}: {} clusters, zero={}", tv, c.len(), zero);
        if zero {
            eprintln!("    zero-site reps -> coarser");
            continue;
        }
        clusters = c;
        reps = r;
        eprintln!(
            "  selected base threshold={:.1}: {} clusters, {} reps",
            tv,
            clusters.len(),
            reps.len()
        );
        break;
    }
    // fallback: coarsest with no zero-site reps
    if reps.len() == 1 && clusters[0].len() == n && n > 1 {
        for tv in thresholds.iter().rev() {
            let c = cut_by_threshold(n, &linkage, *tv);
            let r: Vec<usize> = c.iter().map(|cl| pick_min_dist_member(cl, rep_tsv_dist)).collect();
            let r_names: Vec<String> = r.iter().map(|&x| rep_names[x].clone()).collect();
            let r_idx: Vec<usize> = r.iter().map(|&x| rep_idx[x]).collect();
            let (specific, _) = compute_specific_match(&r_names, &r_idx, poc, seqs, filter);
            let mut counts: HashMap<String, usize> = HashMap::new();
            for (_, nm, _) in &specific {
                *counts.entry(nm.clone()).or_insert(0) += 1;
            }
            let zero = r.iter().any(|&x| counts.get(&rep_names[x]).copied().unwrap_or(0) == 0);
            if !zero {
                clusters = c;
                reps = r;
                eprintln!("  fallback threshold={:.1}: {} clusters", tv, clusters.len());
                break;
            }
        }
    }

    // root match rates
    let root_names: Vec<String> = reps.iter().map(|&x| rep_names[x].clone()).collect();
    let root_idx: Vec<usize> = reps.iter().map(|&x| rep_idx[x]).collect();
    let (_, root_match) = compute_specific_match(&root_names, &root_idx, poc, seqs, filter);
    let mut match_rates: HashMap<String, f64> = HashMap::new();
    let mut site_counts: HashMap<String, usize> = HashMap::new();
    for (nm, tot, _, rate) in &root_match {
        match_rates.insert(nm.clone(), *rate);
        site_counts.insert(nm.clone(), *tot);
    }

    // ---- recursive refinement over each root cluster ----
    for &rep in &reps {
        let members: Vec<usize> = clusters
            .iter()
            .find(|c| c.contains(&rep))
            .cloned()
            .unwrap_or_else(|| vec![rep]);
        recursive_refine(
            ctx,
            rep_names,
            rep_tsv_dist,
            rep_idx,
            poc,
            seqs,
            filter,
            &gdist,
            &members,
            rep,
            1,
            &mut current_final,
            &mut final_reps,
        );
    }

    // ---- global pruning ----
    let mut candidate: Vec<String> = final_reps.clone();
    let mut unique: Vec<String> = Vec::new();
    let mut seen_s = std::collections::HashSet::new();
    for r in candidate.drain(..) {
        if seen_s.insert(r.clone()) {
            unique.push(r);
        }
    }
    candidate = unique;
    for it in 0..MAX_PRUNE_ITERS {
        let cidx = idx_for(&candidate, rep_names, rep_idx);
        let (specific, mrows) = compute_specific_match(&candidate, &cidx, poc, seqs, filter);
        let mut counts: HashMap<String, usize> = HashMap::new();
        for (_, nm, _) in &specific {
            *counts.entry(nm.clone()).or_insert(0) += 1;
        }
        let rates: HashMap<String, f64> =
            mrows.iter().map(|(nm, _, _, r)| (nm.clone(), *r)).collect();
        let bad: Vec<String> = candidate
            .iter()
            .filter(|r| {
                counts.get(*r).copied().unwrap_or(0) < MIN_SITES
                    || rates.get(*r).copied().unwrap_or(0.0) <= ctx.threshold
            })
            .cloned()
            .collect();
        if bad.is_empty() {
            eprintln!("  pruning converged after {} iterations: {} reps", it, candidate.len());
            break;
        }
        // remove worst (largest tsv dist)
        let worst = bad
            .iter()
            .max_by_key(|r| {
                let idx = rep_names.iter().position(|x| *x == **r).unwrap_or(0);
                rep_tsv_dist[idx]
            })
            .unwrap()
            .clone();
        candidate.retain(|r| r != &worst);
        eprintln!(
            "  prune iter {}: removed {} (sites={}, rate={:.4}); {} bad remain",
            it + 1,
            worst,
            counts.get(&worst).copied().unwrap_or(0),
            rates.get(&worst).copied().unwrap_or(0.0),
            bad.len() - 1
        );
    }

    let cidx = idx_for(&candidate, rep_names, rep_idx);
    let (final_specific, final_match) = compute_specific_match(&candidate, &cidx, poc, seqs, filter);
    write_final_outputs(ctx, candidate.clone(), &final_specific, &final_match);
    candidate
}

fn pick_min_dist_member(cluster: &[usize], tsv_dist: &[usize]) -> usize {
    let mut best = cluster[0];
    for &m in cluster {
        if tsv_dist[m] < tsv_dist[best] {
            best = m;
        }
    }
    best
}

#[allow(clippy::too_many_arguments)]
fn recursive_refine(
    ctx: &Ctx,
    rep_names: &[String],
    rep_tsv_dist: &[usize],
    rep_idx: &[usize],
    poc: &[u32],
    seqs: &[Vec<u8>],
    filter: &HashMap<u32, String>,
    gdist: &[Vec<f64>],
    members: &[usize],
    parent_rep: usize,
    level: usize,
    current_final: &mut Vec<String>,
    final_reps: &mut Vec<String>,
) {
    let _ = level;
    if members.len() <= 1 {
        push_rep(final_reps, current_final, rep_names[members[0]].clone());
        return;
    }
    // sub-cluster
    let (sub_clusters, _sub_threshold, sub_reps) = refine::cluster_sub_for_sites(
        members,
        gdist,
        rep_names,
        rep_tsv_dist,
        rep_idx,
        poc,
        seqs,
        filter,
        current_final,
        MIN_SITES,
        MAX_CLUSTERS,
    );
    if sub_clusters.len() <= 1 {
        // cannot split further: keep min-dist member
        let best = pick_min_dist_member(members, rep_tsv_dist);
        push_rep(final_reps, current_final, rep_names[best].clone());
        return;
    }
    // match rates for sub reps
    let sub_names: Vec<String> = sub_reps.iter().map(|&x| rep_names[x].clone()).collect();
    let sub_idx: Vec<usize> = sub_reps.iter().map(|&x| rep_idx[x]).collect();
    let (_, sub_match) = compute_specific_match(&sub_names, &sub_idx, poc, seqs, filter);
    let mut rates: HashMap<String, f64> = HashMap::new();
    let mut sites: HashMap<String, usize> = HashMap::new();
    for (nm, tot, _, rate) in &sub_match {
        rates.insert(nm.clone(), *rate);
        sites.insert(nm.clone(), *tot);
    }
    let mut good: Vec<String> = Vec::new();
    let mut refine_list: Vec<usize> = Vec::new();
    let mut keep: Vec<String> = Vec::new();
    for &r in &sub_reps {
        let nm = &rep_names[r];
        let s = sites.get(nm).copied().unwrap_or(0);
        let rate = rates.get(nm).copied().unwrap_or(0.0);
        let cluster_members = sub_clusters.iter().find(|c| c.contains(&r)).map(|c| c.len()).unwrap_or(1);
        if s >= MIN_SITES && rate > ctx.threshold {
            good.push(nm.clone());
        } else if cluster_members > 1 {
            refine_list.push(r);
        } else {
            keep.push(nm.clone());
        }
    }
    if good.is_empty() && refine_list.is_empty() {
        // backtrack to parent min-dist rep
        let best = pick_min_dist_member(members, rep_tsv_dist);
        push_rep(final_reps, current_final, rep_names[best].clone());
        return;
    }
    for nm in good.iter().chain(keep.iter()) {
        push_rep(final_reps, current_final, nm.clone());
    }
    for &r in &refine_list {
        let sub_members: Vec<usize> = sub_clusters
            .iter()
            .find(|c| c.contains(&r))
            .cloned()
            .unwrap_or_else(|| vec![r]);
        recursive_refine(
            ctx,
            rep_names,
            rep_tsv_dist,
            rep_idx,
            poc,
            seqs,
            filter,
            gdist,
            &sub_members,
            r,
            level + 1,
            current_final,
            final_reps,
        );
    }
}

fn push_rep(final_reps: &mut Vec<String>, current_final: &mut Vec<String>, nm: String) {
    if !current_final.contains(&nm) {
        current_final.push(nm.clone());
    }
    if !final_reps.contains(&nm) {
        final_reps.push(nm);
    }
}

fn write_final_outputs(
    ctx: &Ctx,
    reps: Vec<String>,
    specific: &[(usize, String, u8)],
    match_rows: &[(String, usize, usize, f64)],
) {
    let poc = read_poc2(&ctx.out_dir.join("poc2.txt"));
    fs::create_dir_all(&ctx.out_root).unwrap();
    // final_representatives.txt
    {
        let f = File::create(ctx.out_root.join("final_representatives.txt")).unwrap();
        let mut w = BufWriter::new(f);
        for r in &reps {
            writeln!(w, "{}", r).unwrap();
        }
    }
    // final_specific_sites.tsv
    {
        let f = File::create(ctx.out_root.join("final_specific_sites.tsv")).unwrap();
        let mut w = BufWriter::new(f);
        writeln!(w, "strain\taddress\tbase").unwrap();
        for (row, nm, b) in specific {
            let addr = poc.get(*row).copied().unwrap_or(0);
            writeln!(w, "{}\t{}\t{}", nm, addr, *b as char).unwrap();
        }
    }
    // final_match_results.tsv
    {
        let f = File::create(ctx.out_root.join("final_match_results.tsv")).unwrap();
        let mut w = BufWriter::new(f);
        writeln!(w, "Strain\tSpecificSites\tMatched\tMatchRate").unwrap();
        for (nm, tot, mat, rate) in match_rows {
            writeln!(w, "{}\t{}\t{}\t{:.10}", nm, tot, mat, rate).unwrap();
        }
    }
    // final_clusters.tsv
    {
        let f = File::create(ctx.out_root.join("final_clusters.tsv")).unwrap();
        let mut w = BufWriter::new(f);
        writeln!(w, "strain\tcluster_id").unwrap();
        for (i, r) in reps.iter().enumerate() {
            writeln!(w, "{}\t{}", r, i).unwrap();
        }
    }
}
