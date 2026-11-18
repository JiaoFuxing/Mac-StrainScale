#!/bin/bash
# Parallel preStrain: one p process per DB, batched by -j; then per-sample
# representative clustering (pick-reps OR hybrid-rust single binary).
# Usage: bash run_parallel.sh <DBdir> <samples.list> <pickX> [mbgc_bin] [max_jobs] [rep_cutoff] [pickreps_bin] [max_dist] [hybrid_bin] [hyb_outdir]
#   rep_cutoff:    if set (>0), run representative clustering per sample after all DBs finish (cutoff in SNP distance)
#   pickreps_bin:  pick-reps binary path (default: pick-reps on PATH); used only when hybrid_bin is empty
#   max_dist:      strains with really-ref tsv distance > max_dist are excluded (default 10000)
#   hybrid_bin:    if set, replaces pick-reps+snp-dists+fasta2matrix+refinement with the single hybrid-rust binary
#                  (requires YC-result/<sample>/filter-result-with-counts.txt)
#   hyb_outdir:    output dir name under YC-result/<sample>/ for hybrid final artifacts (default: Strain-5k)
#   SKIP_MBGC=1    skip the mbgc stage (reuse existing rep_*.tsv / fa_* dirs)
set -e

DBDIR=${1:?usage: run_parallel.sh <DBdir> <samples.list> <pickX> [mbgc_bin] [max_jobs] [rep_cutoff] [pickreps_bin] [max_dist] [hybrid_bin] [hyb_outdir]}
SAMPLES=${2:?missing samples list}
PICKX=${3:?missing pickX}
MBGC=${4:-mbgc.zstd.batch}
MAXJOBS=${5:-12}
REPCUT=${6:-0}
PICKREPS=${7:-pick-reps}
MAXDIST=${8:-10000}
HYBRID=${9:-}
HYBROOT=${10:-Strain-5k}

# 任一实测样本深度 < 50x 时，放宽候选距离上限（默认 10000 -> 20000）
LOW_DEPTH=0
while IFS=$'\t' read -r sample _rest; do
    [ -n "$sample" ] || continue
    depth_file="YC-result/$sample/estimated_depth.txt"
    if [ -r "$depth_file" ]; then
        d=$(tr -dc '0-9' < "$depth_file")
        if [ -n "$d" ] && [ "$d" -lt 50 ]; then
            LOW_DEPTH=1
        fi
    fi
done < "$SAMPLES"
if [ "$LOW_DEPTH" -eq 1 ] && [ "$MAXDIST" -le 10000 ]; then
    MAXDIST=20000
    echo "检测到测序深度 < 50x，max_dist 调整为 20000"
fi

# 转绝对路径，避免软链接用相对路径导致找不到
DBDIR="$(cd "$DBDIR" && pwd)"

# ---------------- stage 1: mbgc p per DB ----------------
if [ "${SKIP_MBGC:-0}" != "1" ]; then
    WORKDIR=$(mktemp -d /tmp/mbgc_parallel.XXXXXX)
    trap 'rm -rf "$WORKDIR"' EXIT

    # DB file location: prefer $DBDIR/DB/*.zst (nested), else flat $DBDIR/*.zst
    if [ -d "$DBDIR/DB" ] && ls "$DBDIR"/DB/*.zst >/dev/null 2>&1; then
        DBSUB="$DBDIR/DB"
    else
        DBSUB="$DBDIR"
    fi
    # db list: DB.list preferred (order), else glob scan of DBSUB
    # bash 3.2 兼容（macOS 默认）：用 while read 替代 mapfile
    ALLDBS=()
    if [ -f "$DBDIR/DB.list" ]; then
        while IFS= read -r line; do
            line=$(echo "$line" | tr -d '\r')
            case "$line" in ''|\#*) continue;; esac
            ALLDBS+=("$line")
        done < "$DBDIR/DB.list"
    else
        for f in "$DBSUB"/*.zst; do
            [ -f "$f" ] || continue
            ALLDBS+=("$(basename "$f" .zst)")
        done
    fi
    echo "DB root: $DBSUB (${#ALLDBS[@]} dbs)"
    JOBS=0
    echo "== launching ${#ALLDBS[@]} db jobs, max $MAXJOBS at a time =="
    T0=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
    for db in "${ALLDBS[@]}"; do
        if [ ! -f "$DBSUB/$db.zst" ] || [ ! -f "$DBSUB/$db.clu" ] || [ ! -f "$DBSUB/$db.mbgc" ]; then
            echo "SKIP $db: incomplete"
            continue
        fi
        W=$WORKDIR/$db
        mkdir -p "$W/DB"
        ln -sf "$DBDIR/G0.fa" "$W/G0.fa"
        ln -sf "$DBSUB/$db.zst" "$W/DB/$db.zst"
        ln -sf "$DBSUB/$db.clu" "$W/DB/$db.clu"
        ln -sf "$DBSUB/$db.mbgc" "$W/DB/$db.mbgc"
        echo "$db" > "$W/DB.list"
        echo "[$db] start"
        "$MBGC" p -f -d "$W" "$SAMPLES" "$PICKX" "$MAXDIST" > "stdout_$db.log" 2> "stderr_$db.log" &
        JOBS=$((JOBS + 1))
        # throttle: wait if we hit MAXJOBS
        while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do
            sleep 0.5
        done
    done
    wait
    T1=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
    echo "== all $JOBS jobs done, wall: $(python3 -c 'print(f"{'"$T1"' - '"$T0"':.2f}")') s = $(python3 -c 'print(f"{(('"$T1"' - '"$T0"')/60):.1f}")') min =="

    echo
    echo "== per-db timing =="
    for db in "${ALLDBS[@]}"; do
        [ -f "stderr_$db.log" ] || continue
        s=$(grep -Eo 'scoring [0-9.]+s' "stderr_$db.log" | tail -1 | grep -Eo '[0-9.]+')
        e=$(grep -Eo 'extract [0-9.]+s' "stderr_$db.log" | tail -1 | grep -Eo '[0-9.]+')
        [ -z "$e" ] && e=$(grep -Eo 'sliced-extract [0-9.]+s' "stderr_$db.log" | tail -1 | grep -Eo '[0-9.]+')
        echo "$db: scoring=${s}s  extraction=${e}s"
    done

    echo
    echo "== outputs check =="
    fail=0
    while IFS= read -r line; do
        s=$(echo "$line" | awk '{print $1}')
        [ -n "$s" ] || continue
        case "$s" in \#*) continue;; esac
        for db in "${ALLDBS[@]}"; do
            f="YC-result/$s/output/rep_$db.tsv"
            [ -f "$f" ] || { echo "MISSING: $f"; fail=1; }
        done
    done < <(grep -v '^#' "$SAMPLES")
    [ "$fail" = 0 ] && echo "all rep tables present"
else
    echo "== SKIP_MBGC=1: reusing existing rep_*.tsv / fa_* dirs =="
fi

# ---------------- stage 2: per-sample representatives ----------------
if [ "$REPCUT" -gt 0 ] 2>/dev/null; then
    if [ -n "$HYBRID" ]; then
        if ! command -v "$HYBRID" >/dev/null 2>&1 && [ ! -x "$HYBRID" ]; then
            echo "ERROR: hybrid binary not found: $HYBRID" >&2
            exit 1
        fi
        echo
        echo "== hybrid-rust per sample: cutoff=$REPCUT max_dist=$MAXDIST outdir=$HYBROOT, max $MAXJOBS at a time =="
        T0=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
        while IFS=$'\t' read -r s rest; do
            case "$s" in ''|\#*) continue;; esac
            OUT="YC-result/$s/output"
            [ -d "$OUT" ] || { echo "SKIP $s: no $OUT"; continue; }
            FILTER="YC-result/$s/filter-result-with-counts.txt"
            if [ ! -f "$FILTER" ]; then
                echo "SKIP $s: no $FILTER (needed for refinement)"
                continue
            fi
            echo "[$s] start"
            "$HYBRID" "$OUT" "$REPCUT" "$MAXDIST" "$FILTER" -o "YC-result/$s/$HYBROOT" \
                > "stdout_hybrid_$s.log" 2> "stderr_hybrid_$s.log" &
            while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 0.5; done
        done < <(grep -v '^#' "$SAMPLES")
        wait
        T1=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
        echo "== all hybrid jobs done, wall: $(python3 -c 'print(f"{'"$T1"' - '"$T0"':.2f}")') s =="
        echo
        echo "== hybrid per-sample result =="
        while IFS=$'\t' read -r s rest; do
            case "$s" in ''|\#*) continue;; esac
            OUT="YC-result/$s/$HYBROOT"
            if [ -f "$OUT/final_representatives.txt" ]; then
                n=$(wc -l < "$OUT/final_representatives.txt")
                t=$(grep -Eo 'hybrid total [0-9.]+(ms|s)' "stderr_hybrid_$s.log" | head -1 | grep -Eo '[0-9.]+(ms|s)')
                echo "$s: $n final reps (hybrid ${t})"
            else
                echo "$s: MISSING final_representatives.txt"
            fi
        done < <(grep -v '^#' "$SAMPLES")
    else
        echo
        echo "== picking representatives per sample, cutoff=$REPCUT, max $MAXJOBS at a time =="
        if ! command -v "$PICKREPS" >/dev/null 2>&1 && [ ! -x "$PICKREPS" ]; then
            echo "ERROR: pick-reps binary not found: $PICKREPS" >&2
            exit 1
        fi
        T0=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
        while IFS=$'\t' read -r s rest; do
            case "$s" in ''|\#*) continue;; esac
            OUT="YC-result/$s/output"
            [ -d "$OUT" ] || { echo "SKIP $s: no $OUT"; continue; }
            echo "[$s] start"
            "$PICKREPS" "$OUT" "$REPCUT" "$MAXDIST" > "stdout_pickreps_$s.log" 2> "stderr_pickreps_$s.log" &
            while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 0.5; done
        done < <(grep -v '^#' "$SAMPLES")
        wait
        T1=$(date +%s.%N 2>/dev/null || python3 -c 'import time; print(time.time())')
        echo "== all pick-reps jobs done, wall: $(python3 -c 'print(f"{'"$T1"' - '"$T0"':.2f}")') s =="
        echo
        echo "== pick-reps per-sample result =="
        while IFS=$'\t' read -r s rest; do
            case "$s" in ''|\#*) continue;; esac
            OUT="YC-result/$s/output"
            if [ -f "$OUT/representatives.txt" ]; then
                n=$(wc -l < "$OUT/representatives.txt")
                t=$(grep -oP 'Clustering done in \K[0-9.]+' "stderr_pickreps_$s.log" | head -1)
                echo "$s: $n reps (cluster ${t}s)"
            else
                echo "$s: MISSING representatives.txt"
            fi
        done < <(grep -v '^#' "$SAMPLES")
    fi
fi
