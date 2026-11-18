# MBGC + zstd 集成版（mbgc-zstd）

> 2026-09-21 创建。将 mbgc2 壳的 zstd 自动压缩/解压能力**直接集成进 mbgc 内核**，
> 不再需要 mbgc2 包装器。产物：`mbgc.zstd`（完全静态，默认 `-march=native`）。

## 产物

| 文件 | 说明 |
|---|---|
| `mbgc.zstd` | 集成 zstd 的内核（静态链接，4.87 MB，v2.1.6 + g0-compare + zstd） |
| `src/` | 官方 v2.1.6 源码（`main.cpp` 已集成 zstd 逻辑） |
| `vendor/zstd-1.5.6/` | zstd 1.5.6 源码与静态库（`lib/libzstd.a`） |
| `build-mbgc-zstd.sh` | 构建脚本（默认静态 + `-march=native`） |
| `verify.sh` | 一键验证脚本（7 项测试） |

## 用法（与 mbgc2 完全一致，但少一层壳）

```bash
# g -D：diff 文件名以 .zst 结尾时，比较完成后自动 zstd 压缩
./mbgc.zstd g -D diff.zst sites.txt vp_fast.mbgc out.tsv

# g -q：自动检测 zstd 魔数并解压后查询（原始 diff 也兼容）
./mbgc.zstd g -q diff.zst sites.txt out.tsv

# g -q -B：批量查询，diff 只加载一次，多个样本共享
#   batch.list 每行: <sitesFile> <outputFile>（空行和 # 注释行跳过）
#   多样本用"并集单遍"算法：所有样本位点合并成一张标记表，diff 流只解析
#   一遍，每个位点同时累加所有包含它的样本；耗时和内存不随样本数增长
./mbgc.zstd g -q -B batch.list diff.zst

# 其余命令（c/d/i/a/r）与官方 mbgc 完全一致
./mbgc.zstd c seqlist.txt archive.mbgc
```

环境变量：`MBGC_ZSTD_LEVEL`（默认 19）、`MBGC_ZSTD_THREADS`（默认 4）。

## 兼容性

- `g -D` 输出 `.zst` 与 mbgc2 生成的压缩 diff **逐字节一致**（已验证）
- `g -D` 不带 `.zst` 后缀 → 保持原始 diff，向后兼容
- `g -q` 对 zstd 与原始 diff 均自动处理
- 压缩/解压档案（c/d 命令）行为不变，已验证逐字节还原

## 构建

```bash
bash build-mbgc-zstd.sh            # 静态（默认）
MARCH= bash build-mbgc-zstd.sh     # 关闭 -march=native（通用 x86-64）
BUILD_DYNAMIC=1 bash build-mbgc-zstd.sh   # 动态链接
```

需要 gcc/g++；静态版需要静态 glibc（libc6-dev/build-essential）。
