//! FASTA loading and small sequence utilities.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

pub const IGNORE: u8 = b'.';

/// Normalize a base: keep A/C/G/T, map everything else to IGNORE ('.').
pub fn norm(b: u8) -> u8 {
    let c = b.to_ascii_uppercase();
    if c == b'A' || c == b'C' || c == b'G' || c == b'T' {
        c
    } else {
        IGNORE
    }
}

/// Read one FASTA file (assumed to contain a single record).
pub fn read_fasta(path: &Path) -> Vec<(String, Vec<u8>)> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => return vec![],
    };
    let reader = BufReader::new(file);
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        let t = line.trim();
        if t.starts_with('>') {
            if !name.is_empty() {
                out.push((name, seq));
            }
            name = t[1..].to_string();
            seq = Vec::new();
        } else if !t.is_empty() {
            seq.extend(t.as_bytes().iter().map(|&b| norm(b)));
        }
    }
    if !name.is_empty() {
        out.push((name, seq));
    }
    out
}

/// Scan fa_*/ directories under `out_dir` for a name -> path index.
pub fn build_name_index(out_dir: &Path, yc_dir: &Path) -> HashMap<String, PathBuf> {
    let mut idx = HashMap::new();
    if let Ok(rd) = fs::read_dir(out_dir) {
        for e in rd.flatten() {
            let p = e.path();
            let fname = p.file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_default();
            if p.is_dir() && fname.starts_with("fa_") {
                if let Ok(frd) = fs::read_dir(&p) {
                    for fe in frd.flatten() {
                        let fp = fe.path();
                        if fp.is_file() {
                            if let Some(n) = fp.file_name() {
                                idx.insert(n.to_string_lossy().to_string(), fp);
                            }
                        }
                    }
                }
            }
        }
    }
    if idx.is_empty() {
        if let Ok(rd) = fs::read_dir(yc_dir) {
            for e in rd.flatten() {
                let p = e.path();
                if p.is_dir() && p.file_name().map(|f| f.to_string_lossy().starts_with("_shared_")).unwrap_or(false) {
                    if let Ok(frd) = fs::read_dir(&p) {
                        for fe in frd.flatten() {
                            let fp = fe.path();
                            if fp.is_file() {
                                if let Some(n) = fp.file_name() {
                                    idx.insert(n.to_string_lossy().to_string(), fp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    idx
}

/// SNP distance: positions where both are called and differ.
pub fn snp_dist(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b.iter())
        .filter(|(&x, &y)| x != y && x != IGNORE && y != IGNORE)
        .count()
}
