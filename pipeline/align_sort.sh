#!/usr/bin/env bash
# 按样本列表将 FASTQ 比对到参考序列，并用 sambamba 排序、建索引，输出 BAM。
# 默认按 CPU 核数自动并行处理多个样本；-j 1 可恢复串行。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIST=""
REF="$script_dir/minimap2-rust-rammap/test-data/ref.fa"
OUTDIR=""
THREADS=8
JOBS="auto"
PRESET="auto"
RAMMAP=""
SAMBAMBA="$script_dir/sambamba-master/bin/sambamba-1.0.1"

cpu_count() {
    local cores=""
    cores=$(sysctl -n hw.logicalcpu 2>/dev/null || true)
    if [[ -z "$cores" ]]; then
        cores=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
    fi
    if [[ -z "$cores" ]]; then
        cores=4
    fi
    echo "$cores"
}
CORES="$(cpu_count)"

usage() {
    cat <<EOF
用法: $(basename "$0") [选项] <list> <outdir>

  list       样本列表路径
             每行: 样本名<TAB>fastq1[<TAB>fastq2]，双端时第三列为 read2
  outdir     输出目录
  -r REF     参考 FASTA（默认: ${REF}）
  -t N       总线程预算，按 -j 均分给每个样本（默认: ${THREADS}）
  -j N       并行样本数（默认 auto：最多 4，且不超过 CPU 核数/2；-j 1 为串行）
  -x PRESET  比对预设（默认 auto：双端=sr，单端=map-ont）
  -R PATH    rammap 可执行文件（默认: 自动探测）
  -S PATH    sambamba 可执行文件（默认: ${SAMBAMBA}）
  -h         显示帮助

输出: outdir/<样本名>.sorted.bam 及对应 .bai
      outdir/sample.list（每行: 样本名<TAB>BAM 绝对路径）
EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

while getopts "r:t:x:R:S:j:h" opt; do
    case "$opt" in
        r) REF="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        x) PRESET="$OPTARG" ;;
        R) RAMMAP="$OPTARG" ;;
        S) SAMBAMBA="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [[ -z "$RAMMAP" ]]; then
    for candidate in \
        "$script_dir/minimap2-rust-rammap/target/release/rammap" \
        "$script_dir/server-deploy/rammap-macos-arm64" \
        "$script_dir/bin/rammap"; do
        if [[ -x "$candidate" ]]; then
            RAMMAP="$candidate"
            break
        fi
    done
fi
RAMMAP="${RAMMAP:-$script_dir/minimap2-rust-rammap/target/release/rammap}"

[[ $# -ge 2 ]] || { usage; exit 1; }
LIST="$1"
OUTDIR="$2"

[[ "$THREADS" =~ ^[0-9]+$ ]] && [[ "$THREADS" -ge 1 ]] || die "线程数必须是正整数: $THREADS"
if [[ "$JOBS" != auto ]]; then
    [[ "$JOBS" =~ ^[0-9]+$ ]] && [[ "$JOBS" -ge 1 ]] || die "并行样本数必须是正整数或 auto: $JOBS"
fi

[[ -s "$LIST" ]] || die "样本列表不存在或为空: $LIST"
[[ -s "$REF" ]] || die "参考序列不存在或为空: $REF"
[[ -x "$RAMMAP" ]] || die "rammap 不可执行: $RAMMAP"
[[ -x "$SAMBAMBA" ]] || die "sambamba 不可执行: $SAMBAMBA"
mkdir -p "$OUTDIR/tmp"
OUTDIR="$(cd "$OUTDIR" && pwd)"
list_out="$OUTDIR/sample.list"
: > "$list_out"

jobs_file="$OUTDIR/tmp/.jobs.tsv"
order_file="$OUTDIR/tmp/.order.tsv"
seen_file="$OUTDIR/tmp/.seen.txt"
: > "$jobs_file"
: > "$order_file"
: > "$seen_file"

total=0
while IFS=$'\t' read -r name fq1 fq2 _rest || [[ -n "${name:-}" ]]; do
    name="${name%$'\r'}"
    fq1="${fq1%$'\r'}"
    fq2="${fq2%$'\r'}"
    [[ -z "$name" || "$name" == \#* ]] && continue
    total=$((total + 1))
    printf '%s\t%s\n' "$total" "$name" >> "$order_file"

    reason=""
    if [[ "$name" == */* || "$name" == *..* ]]; then
        reason="非法样本名: $name"
    elif grep -Fxq -- "$name" "$seen_file"; then
        reason="样本名重复: $name"
    elif [[ -z "$fq1" || ! -s "$fq1" ]]; then
        reason="找不到 ${fq1:-<空>}"
    elif [[ -n "$fq2" && ! -s "$fq2" ]]; then
        reason="找不到 $fq2"
    fi
    printf '%s\n' "$name" >> "$seen_file"

    if [[ -n "$reason" ]]; then
        echo "[SKIP] $name: $reason" >&2
        printf 'FAIL\t%s\t%s\n' "$name" "$reason" > "$OUTDIR/tmp/${total}.status"
        continue
    fi

    if [[ -n "$fq2" ]]; then
        if [[ "$PRESET" == auto ]]; then preset="sr"; else preset="$PRESET"; fi
    else
        if [[ "$PRESET" == auto ]]; then preset="map-ont"; else preset="$PRESET"; fi
    fi

    # bash 的 read 在 IFS 为制表符时会吞掉连续空列，用 "-" 占位保留单端样本的空 fq2。
    job_fq2="$fq2"
    if [[ -z "$job_fq2" ]]; then
        job_fq2="-"
    fi

    printf '%s\t%s\t%s\t%s\t%s\n' "$total" "$name" "$fq1" "$job_fq2" "$preset" >> "$jobs_file"
done < "$LIST"

valid_count=$(wc -l < "$jobs_file" | tr -d ' ')
if [[ "$JOBS" == auto ]]; then
    JOBS=4
    if [[ "$JOBS" -gt "$valid_count" ]]; then JOBS=$valid_count; fi
    half=$((CORES / 2))
    if [[ "$half" -lt 1 ]]; then half=1; fi
    if [[ "$JOBS" -gt "$half" ]]; then JOBS=$half; fi
fi
if [[ "$JOBS" -lt 1 ]]; then JOBS=1; fi
job_threads=$((THREADS / JOBS))
if [[ "$job_threads" -lt 1 ]]; then job_threads=1; fi

echo "并行样本数: ${JOBS}，每样本线程: ${job_threads}（总线程预算: ${THREADS}，CPU: ${CORES}）"

process_sample() {
    local idx="$1" name="$2" fq1="$3" fq2="$4" preset="$5"
    if [[ "$fq2" == "-" ]]; then
        fq2=""
    fi
    local inputs=("$fq1")
    if [[ -n "$fq2" ]]; then
        inputs+=("$fq2")
    fi

    local out="$OUTDIR/${name}.sorted.bam"
    local tmp="$OUTDIR/tmp/${idx}.unsorted.bam"
    local fifo="$OUTDIR/tmp/${idx}.sam.fifo"
    local log="$OUTDIR/tmp/${idx}.sambamba.log"
    local status="$OUTDIR/tmp/${idx}.status"
    local start=$SECONDS
    local reason=""

    rm -f "$tmp" "$out" "$out.bai" "$fifo" "$log" "$status"
    echo "[RUN ] $name  preset=$preset  fastq=${#inputs[@]}  threads=$job_threads"

    # sambamba view 只接受文件名参数，不读 stdin；用命名管道传 SAM 流，
    # 这样可以分别拿到 rammap 和 sambamba 的退出码。
    mkfifo "$fifo"
    "$RAMMAP" -t "$job_threads" -x "$preset" -a --cs=long "$REF" "${inputs[@]}" > "$fifo" &
    local ram_pid=$!

    if ! "$SAMBAMBA" view -S -f bam -o "$tmp" "$fifo" 2>"$log"; then
        reason="sambamba 转换 BAM 失败"
        kill "$ram_pid" 2>/dev/null || true
    fi
    local ram_status=0
    wait "$ram_pid" 2>/dev/null || ram_status=$?
    rm -f "$fifo"
    if [[ -z "$reason" && "$ram_status" -ne 0 ]]; then
        reason="rammap 比对失败（exit ${ram_status}）"
    fi

    if [[ -z "$reason" && ! -s "$tmp" ]]; then
        reason="生成的 BAM 为空"
    fi
    if [[ -z "$reason" ]] && ! "$SAMBAMBA" sort -t "$job_threads" -m 1GB --tmpdir="$OUTDIR/tmp" -o "$out" "$tmp" 2>>"$log"; then
        reason="排序失败"
    fi
    if [[ -z "$reason" ]] && ! "$SAMBAMBA" index -t "$job_threads" "$out" 2>>"$log"; then
        reason="建索引失败"
    fi
    if [[ -z "$reason" ]] && ! "$SAMBAMBA" view -H "$out" 2>/dev/null | grep -q 'SO:coordinate'; then
        reason="BAM 头缺少 SO:coordinate"
    fi

    if [[ -n "$reason" ]]; then
        echo "[FAIL] $name: $reason" >&2
        if [[ -s "$log" ]]; then
            cat "$log" >&2
        fi
        rm -f "$tmp" "$out" "$out.bai" "$log"
        printf 'FAIL\t%s\t%s\n' "$name" "$reason" > "$status"
        return 1
    fi

    local records
    records=$("$SAMBAMBA" view -c "$out" 2>/dev/null || echo "?")
    local elapsed=$((SECONDS - start))
    echo "[OK  ] $name: $records 条记录, ${elapsed}s -> $out"
    rm -f "$tmp" "$log"
    printf 'OK\t%s\t%s\t%s\n' "$name" "$out" "$records" > "$status"
}

# 滚动执行：完成一个样本立刻补位，避免长短样本造成空转。
while IFS=$'\t' read -r idx name fq1 fq2 preset; do
    while [[ $(jobs -pr | wc -l | tr -d ' ') -ge "$JOBS" ]]; do
        sleep 0.2
    done
    process_sample "$idx" "$name" "$fq1" "$fq2" "$preset" &
done < "$jobs_file"
wait || true

ok=0
fail=0
while IFS=$'\t' read -r idx name; do
    status="$OUTDIR/tmp/${idx}.status"
    if [[ ! -s "$status" ]]; then
        echo "[FAIL] $name: 未产生结果（进程异常退出）" >&2
        fail=$((fail + 1))
        continue
    fi
    IFS=$'\t' read -r flag field1 field2 _rest2 < "$status"
    if [[ "$flag" == "OK" ]]; then
        printf '%s\t%s\n' "$field1" "$field2" >> "$list_out"
        ok=$((ok + 1))
    else
        fail=$((fail + 1))
    fi
done < "$order_file"

rm -f "$jobs_file" "$order_file" "$seen_file" "$OUTDIR/tmp"/*.status
rmdir "$OUTDIR/tmp" 2>/dev/null || true

echo "完成: $ok/$total 成功, $fail 失败"
echo "BAM 列表: $list_out"
[[ $fail -eq 0 ]]
