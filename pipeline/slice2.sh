#!/bin/bash
if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then  
    echo "Usage: $0 <bam> <DB> <num_segments> [min_mq] [min_bq]"  
    exit 1  
fi  

input_bam=$1
DB=$2
num_segments=$3
min_mq=${4:-20}
min_bq=${5:-20}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# 去除路径末尾的斜杠
if [[ "${DB: -1}" == "/" ]]; then
    DB="${DB%/}"
fi

output_prefix="./output/"
mkdir -p "$output_prefix"

# 提取参考序列名（第一行fasta头，去掉>和后续注释）
ref_name=$(awk 'NR==1 {gsub(/^>/, ""); print $1}' "$DB/ref.fa")

# ========== 修复1：正确计算参考序列总长度（兼容多行FASTA） ==========
# 原脚本只取第2行长度，标准FASTA是多行序列，会导致总长度计算偏小
total_length=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "$DB/ref.fa")

# 计算基础分段长度（整数除法，保留商）
segment_length=$((total_length / num_segments))
echo "基础分段长度: $segment_length, 总长度: $total_length"

for ((i=1; i<=num_segments; i++))
do
    (
        start=$(( (i - 1) * segment_length + 1 ))
        
        # ========== 修复2：最后一段直接取总长度，避免余数丢失 ==========
        if [ $i -eq $num_segments ]; then
            end=$total_length
        else
            end=$(( i * segment_length ))
        fi

        # ========== 切片边界左右各外扩 1000 bp，供 mpileup 参考读段跨越边界 ==========
        slice_start=$(( start - 1000 ))
        if [ $slice_start -lt 1 ]; then
            slice_start=1
        fi
        slice_end=$(( end + 1000 ))
        if [ $slice_end -gt $total_length ]; then
            slice_end=$total_length
        fi

        output_file="${output_prefix}${i}.bam"
        add_file="${output_prefix}${i}.add"

        # 位点仍按原区间（start-end）输出，避免分段间重复
        awk -v start="$start" -v end="$end" -v ref_name="$ref_name" \
            '{if ($1 >= start && $1 <= end) print ref_name"\t"$1}' \
            "$DB/add" > "$add_file"

        # 执行sambamba切片（外扩边界） + mpileup
        "$script_dir/sambamba-master/bin/sambamba-1.0.1" -q slice "$input_bam" "$ref_name:$slice_start-$slice_end" \
            | "$script_dir/mpileup-standalone/mpileup" \
            -q "$min_mq" -Q "$min_bq" --no-output-ends --no-output-del --no-output-ins -A -a -l "$add_file" - \
            > "${output_prefix}${i}-out"
    ) &
done

wait
echo "BAM analysis is OK"

# 合并结果并清理临时文件
seq 1 "$num_segments" | xargs -I{} cat "${output_prefix}{}-out" > "${output_prefix}2lie-mp-samtool"
rm -f ${output_prefix}*.add
rm -f ${output_prefix}*-out
