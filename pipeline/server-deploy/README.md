# preStrain 服务器部署包

## 内容

| 文件 | 说明 |
|---|---|
| `hybrid` | **静态编译**（musl + mimalloc 分配器, statically linked, 797KB）单二进制：聚类→距离矩阵→碱基矩阵→hybrid 精炼五步合一，任何 Linux x86_64 可直接运行，无任何动态依赖 |
| `mbgc.zstd` | 批处理版（preStrainCommand：really-ref 打分 + poc 截取写 fa + max_dist 前移过滤），静态链接，5.0MB |
| `run_parallel.sh` | 10 参数编排脚本（阶段1 mbgc 批处理 + 阶段2 hybrid），LF 行尾 |
| `build-hybrid.sh` | 服务器端源码重编译（可选，需服务器有 Rust；`bash build-hybrid.sh`=glibc 动态版，`bash build-hybrid.sh musl`=静态版需 zig） |
| `hybrid-rust/` | hybrid 源码（Cargo.toml + src/，依赖 rayon + mimalloc） |

## 上传（本地执行）

```bash
# 假设服务器账号 student05@mu01，目标目录与之前一致
scp -r C:\Users\Mayn\Desktop\preStrain\server-deploy\* \
    student05@mu01:/Stor/home/danielgroup/student05/StrainScale/10.10/
```

或在服务器上 `mkdir -p /Stor/home/danielgroup/student05/StrainScale/10.10/hybrid-deploy` 后上传。

## 服务器运行（全流程：阶段1 + 阶段2）

```bash
# ⚠️ 必做 1：zip 解压会丢失执行权限，先 chmod +x，否则 36 个进程 0.3s 全部秒退（Permission denied）
cd /Stor/home/danielgroup/student05/StrainScale/10.10/hybrid-deploy
chmod +x hybrid mbgc.zstd run_parallel.sh build-hybrid.sh

# ⚠️ 必做 2：务必在 YC-result 所在目录运行（服务器上是 test10.10），否则输出会写到当前目录下的 YC-result/
cd /Stor/home/danielgroup/student05/test10.10

bash /Stor/home/danielgroup/student05/StrainScale/10.10/hybrid-deploy/run_parallel.sh \
    /Stor/home/danielgroup/student05/DB/Ec \
    ./Africa-list \
    50 \
    /Stor/home/danielgroup/student05/StrainScale/10.10/hybrid-deploy/mbgc.zstd \
    36 \
    1000 \
    "" \
    10000 \
    /Stor/home/danielgroup/student05/StrainScale/10.10/hybrid-deploy/hybrid \
    Strain-5k
```

## 参数说明

| # | 参数 | 说明 |
|---|---|---|
| 1 | DBDIR | DB 目录（含 G0.fa + DB/*.zst,*.clu,*.mbgc） |
| 2 | SAMPLES | 样本列表，第一列样本名（tab 分隔） |
| 3 | PICKX | 每 db 取前 N 名代表株（50） |
| 4 | MBGC | mbgc 批处理二进制路径 |
| 5 | MAXJOBS | 阶段1/阶段2 最大并发数（服务器 36，本地 ≤12） |
| 6 | REPCUT | 聚类 cutoff（SNP 距离，如 1000） |
| 7 | PICKREPS | 备用 pick-reps 路径（用 hybrid 时传 `""`） |
| 8 | MAXDIST | 与 really-ref 距离 > 此值的菌株排除（10000） |
| 9 | HYBRID | hybrid 二进制路径；非空则阶段2 走单程序 |
| 10 | HYBROOT | 阶段2 输出目录名（YC-result/<样本>/ 下） |

## 常用变体

- **只跑阶段2（复用已有 rep 表）**：`SKIP_MBGC=1 bash run_parallel.sh ...`
- **只跑阶段1**：第 6 参传 `0`（REPCUT=0 跳过阶段2）
- 阶段2 前自动检查 `YC-result/<样本>/filter-result-with-counts.txt`，缺失样本 SKIP（精炼必需）

## 输出

- 阶段1：`YC-result/<样本>/output/`（rep_<db>.tsv、fa_<db>/ 截取序列、poc2.txt）
- 阶段2（hybrid）：`YC-result/<样本>/Strain-5k/`
  - `final_representatives.txt`、`final_match_results.tsv`、`final_specific_sites.tsv`、`final_clusters.tsv`
  - 中间产物 `representatives.txt/.fa/.dist.tsv/.matrix.addr.tsv`、`all.fa`（--cache 时）
- 日志：`stdout_<db>.log`、`stderr_<db>.log`、`stdout_hybrid_<样本>.log`、`stderr_hybrid_<样本>.log`

## 性能预期（服务器，36 并发）

- 阶段1（36 db × 10 样本）：~28s，峰值内存 ~34GB（36 × ~0.95GB/进程）
- 阶段2（hybrid 10 样本并行）：~20s；静态版单样本 5.9s（本地 9P，含 2s FASTA 加载），服务器原生盘更快
- 内存不够时可降 `MAXJOBS`（每进程约 1GB）

> 注意：静态版使用 mimalloc 分配器，性能与 glibc 动态版相当（本地实测 total 5.94s vs 5.68s）。若用无 mimalloc 的 musl 纯静态版，多线程精炼会慢 ~25 倍（musl 默认分配器无 per-thread arena），勿用。

## 已实测一致性（本地，10 样本）

final reps：SRS3638590=0、ERS8391334=28、SRS3638592=32、SRS3638676=1、ERS8391422=5、
ERS8391615=33、SRS3638605=1、SRS3638588=26、SRS3638609=30、ERS8391456=4

## 备注

- 贪心聚类 / 距离矩阵 / 精炼均在 **poc 截取序列**（每样本 1–23 万位点）上计算，非全长（really-ref 171 万）
- 无 rep 样本（全部 dist>MAXDIST）走快路径，~24ms 输出空 final
- 结果确定性与 Python 基线差 ≤1 株（ERS8391334：Rust 28 vs Python 29），源于贪心 tie-break 排序（Rust 确定性可复现）
