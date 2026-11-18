# 流程说明（Stage A / B / C）

> 本文档写的是**这个 App 怎么调用流程**：输入约定、真实命令行、产物位置、校验规则、常见报错。
> 算法细节（切片、位点判定、UPGMA、递归细化、阈值含义）见随包的上游文档：
> `src-tauri/resources/strainscale-pipeline/StrainScale-agent-package/PIPELINE.md`

---

## 0. 三个阶段的衔接

```
sample.list ──┬─► Stage A (SS.sh) ─► YC-result/<样本>/output/really-ref.fa + poc 文件
              │                          │
              │                          ▼
              └─► Stage B (mbgc p) ─► output/rep_<db>.tsv + output/fa_<db>/
                                         │
                                         ▼
                                    Stage C (hybrid) ─► Strain-result/final_*.{txt,tsv}
```

- Stage B/C **不读** `sample.list` 的第二列，只取第一列的样本名；第二列（BAM 路径）只有 Stage A 用。
- Stage B 要求每个样本的 Stage A 产物已存在（`really-ref.fa`、`ref-coords.txt`）。
- Stage C 要求 `filter-result-with-counts.txt` 存在；批量模式会跳过不满足的样本并在 `skipped` 里说明缺哪一项。

---

## 1. 输入约定

### 1.1 `sample.list`

制表符分隔、两列、无表头：

```
bwa-50-1	/abs/path/bwa-50-1.bam
bwa-10-9	/abs/path/bwa-10-9.bam
```

- 第一列是样本名，同时作为 `YC-result/<样本名>/` 目录名；校验禁止包含 `/` 或 `..`。
- 第二列必须是**存在的 BAM 文件绝对路径**，否则 Stage A 校验失败。
- 空行和 `#` 开头的行会被跳过。

### 1.2 数据库目录

Stage A 需要：

| 文件 | 说明 |
|---|---|
| `ref.fa` | BAM 比对参考；`slice2.sh` 从第一行取参考名并统计总长度 |
| `add` | 可变位点列表；mpileup 的 `-l` 白名单 |

Stage B 额外需要：

| 文件 | 说明 |
|---|---|
| `G0.fa` | G0 参考株在 `add` 位点上的序列 |
| `DB/<db>.zst` | 库序列压缩数据 |
| `DB/<db>.clu` | 菌株名与聚类编号 |
| `DB/<db>.mbgc` | 随机访问索引 |
| `DB.list` | 可选；一行一个库名。没有时扫描 `DB/*.zst` 并按字典序排序 |

库目录判定：存在 `<DBdir>/DB/` 就用它，否则把 `<DBdir>` 本身当库目录。

### 1.3 工作目录

所有阶段都以「工作目录」为 `current_dir` 执行，产物固定落在：

```
<工作目录>/YC-result/<样本名>/
```

工作目录不存在时后端会自动 `create_dir_all`。**不要**把工作目录设成数据库目录本身：Stage B 会在其中建暂存目录（见 3.2），虽然有包含判断兜底，但没有必要冒这个风险。

---

## 2. Stage A：BAM → 样本参考

### 2.1 校验（不通过则不会启动进程）

| 错误文案 | 触发条件 |
|---|---|
| `sample.list does not exist: <path>` | 文件不存在 |
| `database directory does not exist: <path>` | 目录不存在 |
| `database directory is missing ref.fa` | 缺参考 |
| `database directory is missing add` | 缺位点表 |
| `pipeline entry SS.sh does not exist` | 打包资源缺失 |
| `line N needs two tab-separated columns` | 该行不是两列（用空格分隔会命中） |
| `invalid sample name: <name>` | 样本名含 `/` 或 `..` |
| `BAM file does not exist: <path>` | 第二列路径不存在 |

### 2.2 执行

```bash
/bin/bash <pkg>/SS.sh <sample.list> <databaseDir>
# cwd = 工作目录
```

`SS.sh` 内部对每个样本：

```bash
mkdir -p YC-result/<样本>; cd YC-result/<样本>
bash <pkg>/slice2.sh <BAM> <DB> 1
perl <pkg>/get_site.pl ./output/2lie-mp-samtool 1 4
```

`slice2.sh` 在本机固定 `num_segments=1`，即整段一次处理；脚本本身支持多段并行（`num_segments > 1` 时后台 fork），但 App 不暴露这个参数。

### 2.3 产物

相对 `<工作目录>/YC-result/<样本>/`：

| 产物 | 说明 |
|---|---|
| `filter-result-with-counts.txt` | poc 位点的原始五列，Stage C 必需 |
| `output/2lie-mp-samtool` | 合并后的 pileup 五列 |
| `output/really-ref.fa` | 样本自身一致序列（含 `-`） |
| `output/ref-coords.txt` | 确定为单碱基的列号 |
| `output/poc.txt` | poc 位点的比对列号 |
| `output/poc2.txt` | poc 位点的原始坐标 |

后端会收集这 6 类文件；若一个 `really-ref.fa` 都没产出，追加错误 `really-ref.fa was not produced`。即使部分样本失败，其它样本的产物仍会列在 `outputs` 里。

### 2.4 已知依赖

- `SS.sh` 调用 `perl`、`awk`、`sambamba-1.0.1`、`mpileup-standalone`。后两者是随包二进制，不依赖系统 PATH；`perl` 和 `awk` 用系统自带。
- `slice2.sh` 会 `rm -f output/*.add`、`output/*-out`，只清自己目录下的临时文件。
- 包内的 `mpileup`、`sambamba-1.0.1`、`mbgc-macos-arm64`、`hybrid-macos-arm64`、`libomp.dylib` 均已确认为 **Mach-O arm64**，可直接在 Apple Silicon 上运行；换平台必须重新编译（见 `BUILD.md`）。

---

## 3. Stage B：库检索与前 N 候选

### 3.1 校验

| 错误文案 | 条件 |
|---|---|
| `native mbgc binary does not exist` | `server-deploy/mbgc-macos-arm64` 缺失 |
| `libomp.dylib is missing next to the native mbgc binary` | 同目录缺 `libomp.dylib` |
| `database directory does not exist: <path>` | 数据库目录不存在 |
| `database directory is missing G0.fa` | 缺 G0 |
| `samples list does not exist: <path>` | `sample.list` 不存在 |
| `missing really-ref.fa for sample <n>; run Stage A first` | Stage A 未完成 |
| `missing ref-coords.txt for sample <n>; run Stage A first` | Stage A 未完成 |
| `no database found: DB.list is empty and no *.zst under the database directory` | 没发现库 |
| `database file is missing: <db>.{zst,clu,mbgc}` | 库文件不全 |
| `samples list contains no sample` | 列表为空 |

### 3.2 暂存目录（重要）

每个库跑之前需要一份「只含该库」的目录结构，App 会在**工作目录**下建：

```
<工作目录>/.mbgc-run/<db>/
├── G0.fa            → 符号链接到 <DBdir>/G0.fa
├── DB.list          → 单行 "<db>"
├── DB/<db>.zst      → 符号链接到 <DBdir>/DB/<db>.zst
├── DB/<db>.clu      → …
├── DB/<db>.mbgc     → …
└── run.log          ← 该库的 stdout/stderr/退出码
```

每次运行前 `remove_dir_all(.mbgc-run)` 重建。为避免「数据库目录恰好也在 `.mbgc-run` 下」时把数据库本体删掉，代码里有兜底：

```rust
let mut stage_root = working_dir.join(".mbgc-run");
if database_dir.starts_with(&stage_root) {
    stage_root = std::env::temp_dir().join(format!("mbgc-run-{}", std::process::id()));
}
let _ = std::fs::remove_dir_all(&stage_root);
```

> 这个坑真实发生过：旧版本用 `<工作目录>/.mbgc-stage`，用户的 `Vp-db` 正好放在 `<工作目录>/.mbgc-stage/Vp`，清理时把 `G0.fa`/`DB/` 一起删了，随后重建的符号链接指向自己，报 `G0 not found` 且 `exitCode 1`。改独立目录名 + 包含判断后恢复。移植时**不要**把暂存根改回与输入可能重名的路径。

### 3.3 执行

参数默认：`pickX=30`（最小 1）、`maxJobs=4`（夹在 1–12）、`maxDist=10000`。

每个库一个进程，命令：

```bash
<pkg>/server-deploy/mbgc-macos-arm64 p -f -d <stage_dir> <samples_list> <pickX> <maxDist>
# cwd = 工作目录
```

- `p`：批处理模式（打分 + 截取 poc 序列）。
- `-f`：强制重算，不复用缓存。
- `-d`：指定暂存目录（内含 `G0.fa`、`DB.list`、`DB/`）。
- 其余为位置参数：样本表、每库保留数、最大 SNP 距离。

并发：按 `max_jobs` 对库列表分块，每块起 N 个线程，`join` 后汇总；单个库失败不会中断其它库。

### 3.4 产物

```
<工作目录>/YC-result/<样本>/output/rep_<db>.tsv   # 候选株、距离、附加统计（前两列被 Stage C 使用）
<工作目录>/YC-result/<样本>/output/fa_<db>/       # 每个候选株一条 FASTA（poc 位点）
```

日志：`<stage>/<db>/run.log`（完整 stdout/stderr）；后端把每个库的 `exit`、`ok`、stdout/stderr 汇总进 `PipelineResult.stdout`。

### 3.5 资源占用

- 每个 mbgc 进程约 **1 GB** 内存（载入库索引与差分数据）。
- 4 并发 ≈ 4 GB；12 并发 ≈ 12 GB。默认 4 是按 32 GB 机器的保守值，接入新机器时按内存调整。
- 单库 3 样本实测：loadDiff 0.10 s、打分 0.83 s、解码/截取数秒，整体几十秒级；库越大越接近 IO 瓶颈。

---

## 4. Stage C：代表株聚类与精炼

### 4.1 单样本入口（命令保留，界面不暴露）

请求字段：`sampleOutputDir`、`filterCounts`、`resultDir`，可选 `cutoff`、`maxDist`、`matchThreshold`。

校验（错误文案带路径和中文提示）：

| 错误 | 条件 |
|---|---|
| `native hybrid binary does not exist: <path>` | 二进制缺失 |
| `sample output directory does not exist: <path>（Stage A 未产出或样本名写错）` | `<样本>/output` 不存在 |
| `filter counts file does not exist: <path>（Stage A 的 get_site.pl 未写出）` | 缺 counts 文件 |
| `failed to create result directory: <err>` | 建结果目录失败 |

`create_dir_all(result_dir)` 只在没有错误时执行，避免失败时留下一个空的 `Strain-result/` 误导用户。

执行：

```bash
<pkg>/server-deploy/hybrid-macos-arm64 \
  <sampleOutputDir> <cutoff=1000> <maxDist=10000> <filterCounts> \
  -t <matchThreshold=0.9> -o <resultDir>
```

产物（`<样本>/Strain-result/`，扁平）：

```
final_representatives.txt
final_match_results.tsv
final_specific_sites.tsv
final_clusters.tsv
```

四个都缺失才算失败（`hybrid did not produce final result files`）。hybrid 还会写回 `representatives.txt/.fa/.dist.tsv/.matrix.addr.tsv`、`cluster_assignments.tsv` 等中间文件，并读取 `poc2.txt`、`all.fa`。

### 4.2 批量入口（界面按钮走这个）

请求字段：`workingDir`、`samplesList`，可选 `cutoff`、`maxDist`、`matchThreshold`。

逻辑：

1. `sample_names_from_list` 取每行第一列，去重（空行/`#` 跳过）。
2. 对每个样本：

   - 缺 `<样本>/output/` → `skipped: "<name>: 缺少 output/（Stage A 未完成）"`
   - 缺 `<样本>/filter-result-with-counts.txt` → `skipped: "<name>: 缺少 filter-result-with-counts.txt（Stage A 未完成）"`
   - 否则调用单样本逻辑，结果目录固定为 `<样本>/Strain-result/`。
3. `success = failed.is_empty() && !ran.is_empty()`。

返回 `HybridBatchResult`：`ran` / `skipped` / `failed` / `outputs`，界面分别用绿/黄/红三色渲染。`failed` 里带每个样本的具体错误（用 `;` 拼接）。

实测（3 样本）：`ran: 3`、`skipped: []`、`failed: []`，共 12 个 final 产物。

---

## 5. 常见报错对照表

| 现象 | 根因 | 处理 |
|---|---|---|
| `G0 not found: <...>/.mbgc-stage/Vp/G0.fa`，`exitCode 1` | 旧版暂存目录与数据库目录重名，清理时删掉了数据库 | 已修复；升级到 0.5.0，确认暂存根是 `.mbgc-run` |
| `missing really-ref.fa for sample X; run Stage A first` | Stage A 没跑或样本名不一致 | 先跑 Stage A；核对 `YC-result/` 下目录名 |
| `sample output directory does not exist` | 样本名写错 / Stage A 失败 | 看 Stage A 的 `errors`，不要只看 Stage C |
| `filter counts file does not exist` | Stage A 的 `get_site.pl` 未写出 | 检查 Stage A 是否真的产出 `filter-result-with-counts.txt` |
| `native mbgc binary does not exist` / `libomp.dylib is missing` | 资源未打包 / 被 Gatekeeper 隔离 | 重装 DMG；`xattr -dr com.apple.quarantine` |
| `database file is missing: ...` | `DB.list` 指向了不存在的库 | 修正 `DB.list` 或补文件 |
| `SS.sh exited with status` 但无 stderr | 样本目录已存在且被上一次运行占用 | 检查 `YC-result/<样本>` 权限与磁盘空间 |
| Stage B 很慢 / 内存爆 | 并发库数过高 | 降到 2–4；单进程约 1 GB |
| 输出目录不可写 | 输出目录填了只读位置或权限不足 | 换到用户可写目录；留空则回落到列表文件所在目录 |

---

## 6. 换平台时必须改的流程相关代码

| 位置 | 现状 | 移植方向 |
|---|---|---|
| `run_pipeline_blocking` | 写死 `/bin/bash` | Windows 用 `cmd /C` 或 WSL；更稳妥是把 `SS.sh`/`slice2.sh` 改成跨平台脚本或 Rust 实现 |
| `run_mbgc_blocking` | 二进制名 `mbgc-macos-arm64`，依赖同目录 `libomp.dylib` | 按 target 选 `mbgc-windows-x86_64.exe` + `libomp.dll`；或用 `cfg!(target_os)` 分支 |
| `run_hybrid_*` | 二进制名 `hybrid-macos-arm64` | 各平台各自编译（hybrid 是 Rust，交叉编译成本最低） |
| Stage B 暂存 | `#[cfg(unix)] symlink` | Windows 用硬链接/复制；`std::fs::symlink` 在 Windows 上不可用 |
| 产物路径 | `PathBuf::join`，已跨平台 | 一般无需改；注意 shell 脚本内部用的是 `/` |
| pipeline 资源 | `Contents/Resources/resources/...`（双层） | Windows 是 `<app>/resources/...`，`resource_dir()` 已处理，保留拼接即可 |

---
