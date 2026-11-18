#!/bin/bash

if [ "$#" -ne 2 ]; then  
    echo "Usage: $0 <sample.list><DB>"  
    exit 1  
fi  
file=$1
DB=$2
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# 转绝对路径，避免 sambamba-slice 相对路径崩
DB="$(cd "$DB" && pwd)"

if [[ "${DB: -1}" == "/" ]]; then
    DB="${DB%/}"  
fi

estimate_depth() {
 local bam=$1 db=$2
 local ref_name total end depth
 ref_name=$(awk 'NR==1 {gsub(/^>/, ""); print $1}' "$db/ref.fa")
 total=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "$db/ref.fa")
 end=10000
 if [ "${total:-0}" -lt "$end" ]; then
   end=$total
 fi
 depth=$("$script_dir/mpileup-standalone/mpileup" -q 0 -Q 0 -B -d 0 -A -a \
   -r "${ref_name}:1-${end}" "$bam" 2>/dev/null \
   | awk '$5 > 0 {count[$5]++; if (count[$5] > best) {best = count[$5]; depth = $5}} END {print depth + 0}')
 echo "${depth:-0}"
}

# 单样本处理：估计深度 -> 分段并行 mpileup -> get_site.pl
process_sample() {
 local SampleName=$1 BamFile=$2 SegBase=$3
 local depth segments min_mq min_bq error_base
 mkdir -p "YC-result/$SampleName"
 depth=$(estimate_depth "$BamFile" "$DB")
 echo "$depth" > "YC-result/$SampleName/estimated_depth.txt"

 # 高深度样本翻倍分段（上限 20）
 segments=$SegBase
 if [ "$depth" -gt 500 ]; then
   segments=$(( SegBase * 2 ))
   [ "$segments" -gt 20 ] && segments=20
 fi
 [ "$segments" -lt 2 ] && segments=2

 min_mq=20
 min_bq=20
 if [ "$depth" -lt 20 ]; then min_mq=10; min_bq=10; fi
 error_base=1
 if [ "$depth" -lt 30 ]; then error_base=0; fi

 echo "[$SampleName] 估计深度: ${depth}x, 分段: ${segments}, MQ/BQ: ${min_mq}/${min_bq}, error_base: ${error_base}"
 cd "YC-result/$SampleName" || return 1
 bash "$script_dir/slice2.sh" "$BamFile" "$DB" "$segments" "$min_mq" "$min_bq"
 perl "$script_dir/get_site.pl" ./output/2lie-mp-samtool "$error_base" 4
 cd - >/dev/null
}

deal_bam() {
 num=$(wc -l < "$file")
 mkdir -p YC-result

 # 按 CPU 核数规划并发：目标总进程数 = 1.5 × 核数（实测超配 1.5x 最优）
 local ncpu budget max_par seg_base tmpd t
 ncpu=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)
 budget=$(( ncpu * 3 / 2 ))
 max_par=$(( budget / 2 ))          # 每样本至少 2 段
 [ "$max_par" -lt 1 ] && max_par=1
 [ "$max_par" -gt "$num" ] && max_par=$num
 seg_base=$(( budget / max_par ))
 [ "$seg_base" -lt 2 ] && seg_base=2
 [ "$seg_base" -gt 10 ] && seg_base=10

 echo "并行规划: 最多 ${max_par} 个样本 × ~${seg_base} 段 (CPU=${ncpu})"

 tmpd=$(mktemp -d)
 trap 'rm -rf "$tmpd"' EXIT

 for ((t=1; t<=$num; t++)); do
   read SampleName BamFile <<< $(awk -v t="$t" -F'\t' 'NR==t {print $1, $2}' "$file")
   process_sample "$SampleName" "$BamFile" "$seg_base" > "$tmpd/$t.log" 2>&1 &
   while [ "$(jobs -rp | wc -l)" -ge "$max_par" ]; do sleep 0.1; done
 done
 wait

 # 按样本顺序输出日志，保持可读性
 for ((t=1; t<=$num; t++)); do cat "$tmpd/$t.log"; done
}






deal_bam
