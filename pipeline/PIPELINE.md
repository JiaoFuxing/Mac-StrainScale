# StrainScale 流程与算法说明

本文描述本仓库实际执行的完整流程：从 BAM 出发，构建样本自身的一致序列（`really-ref`），
到参考库中检索近缘菌株，再到代表株聚类与基于真实 reads 的精炼，最终输出样本的代表株集合。

流程共两个脚本、四个阶段（阶段 D 由桌面 App 在阶段 C 之后自动执行）：

```
sample.list (样本名 \t BAM 路径)
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│ 阶段 A  SS.sh / slice2.sh / get_site.pl                      │
│ BAM ──(分段)──> mpileup ──(位点判定)──> really-ref.fa        │
│                                      ├─> poc2.txt / poc.txt  │
│                                      └─> filter-result-...txt│
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│ 阶段 B  run_parallel.sh 阶段 1（mbgc p）                     │
│ 参考库菌株 × really-ref ──> SNP 距离 ──> 取前 pickX          │
│ ──> dist ≤ max_dist 过滤 ──> 截取 poc 位点序列 fa_<db>/      │
│ ──> output/rep_<db>.tsv                                      │
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│ 阶段 C  run_parallel.sh 阶段 2（hybrid）                     │
│ 贪心聚类 ──> 距离矩阵 ──> 碱基矩阵 ──> 递归细化 ──> 全局剪枝 │
│ ──> final_representatives.txt 等 4 个结果文件                │
└──────────────────────────────────────────────────────────────┘
        │
        ▼
┌──────────────────────────────────────────────────────────────┐
│ 阶段 D  App：cluster 成员扩展与距离                          │
│ final_representatives → *.clu 找 cluster → 全部成员          │
│ ──> mbgc d -E 流式提取 ──> 与 really-ref 计算 SNP 距离        │
│ ──> cluster_members.tsv / cluster_fa_<db>/                   │
└──────────────────────────────────────────────────────────────┘
```

---

## 1. 输入

### 1.1 样本列表 `sample.list`

制表符分隔两列，无表头：

| 列 | 含义 |
|---|---|
| 1 | 样本名，同时作为 `YC-result/<样本名>/` 目录名 |
| 2 | 该样本 BAM 的绝对路径 |

### 1.2 数据库目录（示例 `Vp-db/`）

| 文件 | 本库实际规模 | 作用 |
|---|---|---|
| `ref.fa` | 5,165,770 bp，头 `RIMD_2210633_NCBI` | BAM 比对的参考基因组；决定切片坐标与分段 |
| `add` | 头行 `address` + 1,801,821 个坐标 | 参与分析的可变位点列表；mpileup 只输出这些位点 |
| `G0.fa` | 1,801,821 bp，头 `Reference_RIMD_2210633.chr12` | **G0 参考株**在 `add` 位点上的序列（阶段 1 用） |
| `DB/<db>.zst` | 131 MB | 参考库全部菌株在 `add` 位点上的压缩序列 |
| `DB/<db>.clu` | 12,461 行 | 菌株名与聚类编号：`菌株名,cluster_id` |
| `DB/<db>.mbgc` | 118 MB | 库序列的随机访问索引（按地址取碱基） |
| `DB.list` | 1 行 | 参与运行的库名，一行一个 |

> `add` 的坐标是 `ref.fa` 上的碱基坐标；`really-ref`、`G0.fa`、库序列都是**按 `add` 位点
> 展开后的等长向量**（长度恒为 1,801,821），因此三者可以直接按“第几个位点”对齐比较。

---

## 2. 阶段 A：BAM → 样本参考（`SS.sh`）

`SS.sh <sample.list> <DB>` 逐样本执行，每个样本在 `YC-result/<样本>/` 下产生 `output/`。

### 2.1 切片与 pileup（`slice2.sh`）

参数：`slice2.sh <bam> <DB> <num_segments> [min_mq] [min_bq]`。

执行前 `SS.sh` 先估计测序深度：取 `ref.fa` 第一个 contig 的前 10000 bp，运行
`mpileup -q 0 -Q 0 -B -d 0 -A -a -r <ref>:1-<end>`，取第 5 列深度的众数，写入
`YC-result/<样本>/estimated_depth.txt`。据此动态选择参数：

| 实测深度 | `num_segments` | mpileup 质量阈值 |
|---|---|---|
| `<= 500` | 10 | `<20` 时 `-q 10 -Q 10`，否则 `-q 20 -Q 20` |
| `> 500` | 20 | `-q 20 -Q 20` |

1. 从 `ref.fa` 第一行取参考名，统计多行 FASTA 的总长度 `total_length`。
2. 按 `total_length / num_segments` 均分区间；最后一段取到 `total_length`，避免余数丢失。
   多段时每段在后台并行执行。
3. 对每段，从 `DB/add` 中筛出落在 `[start, end]` 的坐标，生成 `<i>.add`（作为 `mpileup -l` 的位点白名单）。
4. `sambamba -q slice <bam> <ref>:<start>-<end>` 取区间 reads，流式送给
   `mpileup-standalone/mpileup`：

   ```
   -q <min_mq> -Q <min_bq> --no-output-ends --no-output-del --no-output-ins -A -a
   -l <i>.add
   ```

   - `-q`：最小比对质量；`-Q`：最小碱基质量；深度 `< 20` 时两者取 10，否则取 20
   - `--no-output-ends/del/ins`：不输出 read 末端、删除、插入信息
   - `-A`：不过滤孤儿 read；`-a`：输出全部位点
   - `-l`：只保留 `add` 白名单中的坐标
5. 合并各段结果为 `output/2lie-mp-samtool`（无表头），每行 5 列：

   | 列 | 内容 | 示例 |
   |---|---|---|
   | 1 | 参考坐标（address） | `61` |
   | 2 | 各碱基支持数 | `G:57,T:17` |
   | 3 | 支持度最低碱基的计数（“错误/噪声”计数） | `17` |
   | 4 | 通过的候选碱基（按字母序拼接） | `GT` |
   | 5 | 深度（depth） | `74` |

   单碱基位点形如 `18  T:79  79  T  79`；多碱基位点形如 `61  G:57,T:17  17  GT  74`。

### 2.2 位点判定与参考构建（`get_site.pl`）

根据实测深度选择：

```
深度 < 100  -> perl get_site.pl ./output/2lie-mp-samtool 0 4
深度 >= 100 -> perl get_site.pl ./output/2lie-mp-samtool 1 4
```

`min_depth = 4` 固定；`error_base_number` 为 0 或 1。若输入是现成的
`2lie-mp-samtool`（没有 BAM），则取其前 10000 行第 5 列深度的众数作为实测深度。
逐行处理，`col` 为比对列号（1 起始）：

```
若 第4列 或 第2列 以 '*' 开头（缺失/删除）
    really-ref 追加 '-'，结束

若 第3列（最低支持计数） ≤ error_base_number：
    filtered = { 第2列中计数 > error_base_number 的碱基 }   # 计数 ≤ 阈值视为测序噪声
否则（最低支持计数 > error_base_number，无法判定单碱基）：
    filtered = 第4列的前 16 个字符             # 全部候选碱基，通常 ≥ 2

若 depth > 4 且 |filtered| == 1：
    really-ref 追加该碱基
    记录 ref-coords.txt（列号）
否则：
    really-ref 追加 '-'

若 |filtered| ≥ 2：        # 该位点在样本中呈多碱基（poc 位点）
    poc.txt  ← 列号
    poc2.txt ← 第1列原始坐标
    filter-result-with-counts.txt ← 原始 5 列
```

要点：

- **`really-ref.fa` 是“样本自己的参考序列”**：每位点只有两种取值——一个确定的 A/C/G/T，
  或 `-`（缺失、深度不足、或该位点存在多个可信碱基）。长度 = 比对列数（本库 1,801,821）。
- **`poc`（point of conflict）位点** = 样本内部存在 ≥ 2 种可信碱基的位点。这些位点无法给出
  单一碱基，但恰好携带区分不同菌株的信息，因此单独记录，供阶段 B/C 使用。
- 三个 poc 文件严格一一对应（本样本 80,370 行）：
  - `poc.txt`：比对列号（1 起始，用于定位矩阵行）
  - `poc2.txt`：原始坐标（address，用于回库取碱基、回 BAM 取计数）
  - `filter-result-with-counts.txt`：原始 5 列（阶段 C 判断“真实样本是否观测到某碱基”）

### 2.3 阶段 A 产物

| 文件 | 含义 |
|---|---|
| `output/2lie-mp-samtool` | 合并后的 pileup 五列结果 |
| `output/really-ref.fa` | 样本自身的一致序列（含 `-` 占位） |
| `output/ref-coords.txt` | really-ref 中确定为单碱基的列号 |
| `output/poc.txt` / `poc2.txt` | poc 位点的列号 / 原始坐标 |
| `filter-result-with-counts.txt`（位于 `YC-result/<样本>/`） | poc 位点的原始五列（阶段 C 精炼必需） |
| `estimated_depth.txt`（位于 `YC-result/<样本>/`） | 前 10000 bp 的深度众数，供阶段 A/B 选择参数 |

---

## 3. 阶段 B：库检索与前 N 候选（`run_parallel.sh` 阶段 1，`mbgc p`）

调用（本机实际参数）：

```bash
run_parallel.sh <DBdir> <sample.list> 30 <mbgc> 36 1000 "" 10000 <hybrid> Strain-1k
                  │        │         │    │      │   │    │    │       │         │
                  │        │         │    │      │   │    │    │       │         └ 阶段2输出目录名
                  │        │         │    │      │   │    │    │       └ hybrid 二进制
                  │        │         │    │      │   │    │    └ 阶段2开关（>0 才跑）
                  │        │         │    │      │   └ 备用 pick-reps（用 hybrid 时传 ""）
                  │        │         │    │      └ cutoff（SNP 距离）
                  │        │         │    └ 最大并发
                  │        │         └ 每个库保留前 N 名
                  │        └ 样本列表
                  └ 数据库根目录（含 G0.fa 与 DB/*.zst,*.clu,*.mbgc）
```

> `max_dist` 默认 10000；脚本会读取 `YC-result/<样本>/estimated_depth.txt`，
> 若任一实测深度 `< 50` 且传入值 `≤ 10000`，自动放宽为 **20000**。

执行逻辑：

1. 若 `DBdir/DB/*.zst` 存在则库目录取 `DBdir/DB`，否则取 `DBdir`。
2. 库名列表优先读 `DBdir/DB.list`，否则扫描 `*.zst`。
3. 每个库建一个临时工作目录（`mktemp -d /tmp/mbgc_parallel.XXXXXX`），软链
   `G0.fa`、`<db>.zst`、`<db>.clu`、`<db>.mbgc`，再写单行 `DB.list`。
4. 每库一个 `mbgc p` 进程，最多 `max_jobs` 个并发；`-f` 强制重算。

单库内部算法（流式、向量化）：

```
对每个样本 s：
    以 really-ref(s) 为查询，流式扫描库中全部菌株（G0.fa 作为差分编码的基准载入）
    只在 really-ref 有确定碱基的位点上计 SNP（'-'/'.' 不计）
    一次遍历同时为所有样本打分（3 个样本只扫一遍库）
按 (距离升序) 取前 pickX = 30 名作为候选
候选再做一次精确距离打分，剔除 dist > max_dist = 10000 的菌株
按 poc2.txt 的地址，用 .mbgc 索引从库序列中逐株截取 poc 位点
    ──> output/fa_<db>/<菌株名>          （单记录 FASTA，长度 = poc 位点数）
    ──> output/rep_<db>.tsv              （菌株名 \t really-ref 距离 \t 附加统计）
```

> `rep_<db>.tsv` 的第三列不参与下游计算；阶段 C 只读取前两列。
>
> 阶段 1 的 `max_dist` 过滤在这里第一次生效：距离过远的库菌株不会进入 `fa_<db>/`。

实测（3 样本、1 库）：打分 1.47 s，解码 2.47 s，截取 0.30 s；保留数 4 / 5 / 25。

---

## 4. 阶段 C：代表株聚类与精炼（`run_parallel.sh` 阶段 2，`hybrid`）

前置条件：`YC-result/<样本>/filter-result-with-counts.txt` 必须存在，否则该样本被 SKIP。
本阶段单样本完全独立、可并行，输出到 `YC-result/<样本>/Strain-1k/`。

调用：`hybrid <output_dir> <cutoff> <max_dist> <filter_counts> -o <Strain-1k>`

### C1. 汇总候选菌株（`[1]`）

扫描 `output/` 下**所有** `rep_*.tsv`，取 `dist ≤ max_dist` 的菌株，按名字去重，
再按 `(距离升序, 名字字典序)` 排序，保证结果确定可复现。

### C2. 载入序列（`[2]`）

按菌株名在 `fa_*/` 中定位文件（若无 `fa_*` 目录则回退到 `YC-result/<样本>/_shared_*`），
读取并规范化：A/C/G/T 保留，其余（含 `-`）统一记为 `.`（忽略）。

### C3. 贪心聚类（`[3]`，pick-reps 语义）

```
remaining = 全部菌株（已按距离升序）
while remaining 非空:
    rep = remaining[0]                        # 距离最小者为种子
    对 remaining 中每个菌株计算 SNP 距离(rep, 菌株)
        距离 ≤ cutoff(1000) → 归入 rep 的 cluster
        否则              → 留到下一轮
    输出 rep（命名 S1_..., S2_..., ...）
```

SNP 距离定义（`fasta::snp_dist`）：

```
dist(a,b) = count{ i : a[i] != b[i] 且 a[i] != '.' 且 b[i] != '.' }
```

即**只统计双方都有确定碱基的位点**，缺失位点不计差异也不计分母。

产物：`representatives.txt`（`S<i> \t 菌株名 \t tsv_dist=...`）、`representatives.fa`、
`cluster_assignments.tsv`（代表株 → 成员）。

### C4. 距离矩阵与碱基矩阵（`[5]`、`[6]`）

- `[5]` 代表株两两 SNP 距离 → `representatives.dist.tsv`（snp-dists 格式，对称阵）。
- `[6]` 读取 `poc2.txt`（地址列表），构造 `地址 × 代表株` 的碱基矩阵
  → `representatives.matrix.addr.tsv`。

### C5. 读取真实 reads 计数（`[7]`）

把 `filter-result-with-counts.txt` 载入内存：`地址 → "A:57,G:4,..."`。
它表示**真实样本在 poc 位点上实际观测到的碱基及支持数**，是后续“匹配率”的唯一依据。

### C6. 精炼：基础聚类 → 递归细化 → 全局剪枝（`[8]`）

#### (a) 基础聚类（UPGMA + 细到粗阈值扫描）

1. 对代表株距离矩阵做 **UPGMA（平均连接）**，得到 linkage；
   Lance–Williams 更新：`d(new,k) = (n_a·d(a,k) + n_b·d(b,k)) / (n_a+n_b)`。
2. 候选阈值：`min → median` 等分 30 步，`median → max` 等分 10 步，
   四舍五入到 0.01，**去重后按从细（小）到粗（大）排序**。
3. 每个阈值用并查集 `cut_by_threshold` 切树：仅当合并高度 ≤ 阈值时合并两簇
   （等价于 scipy `fcluster(criterion="distance")`）。
4. 每个簇取 **`tsv_dist` 最小的成员**为代表。
5. 对候选代表集计算“特异位点”：在某 poc 位点，若某碱基在**全部代表中只出现 1 次**，
   则该位点对该代表是特异的。基础聚类的要求是：
   - 簇数 ≤ `MAX_CLUSTERS = 50`；
   - 从细到粗扫描阈值，选第一个**所有代表都至少有 1 个特异位点**的划分；
   - 若最初就退化为单簇（1 个代表覆盖全部候选），则从粗到细回退，取第一个无 0 特异位点的划分。
   `MIN_SITES = 30` 是本阶段（递归细化、全局剪枝）的门槛，不用于基础阈值选择。

#### (b) 递归细化（`recursive_refine`）

对每个基础簇（成员数 > 1）：

1. 在簇内子距离矩阵上再做 UPGMA；阈值取 `min..max` 等分 25 步（细→粗），
   扫描满足“子簇数 ≤ 50 且每个子代表特异位点数 ≥ 30”的最细划分；
   若无满足者，取第一个“无 0 特异位点”的划分。
2. 每个子代表按真实 reads 计算**匹配率**：

   ```
   match_rate = #(特异位点上 filter 计数中该碱基 count > 0) / #(特异位点)
   ```

3. 分类处理：
   - `特异位点 ≥ 30` 且 `match_rate > 0.9` → **good**，直接保留；
   - 不达标且簇内成员 > 1 → 进入下一层递归；
   - 不达标且簇内仅 1 个成员 → **keep**，保留（无法再分）；
   - 若既无 good 也无可递归者 → **回退**到父簇中 `tsv_dist` 最小的成员。

#### (c) 全局剪枝

```
候选 = 去重后的全部 final_reps
重复（最多 MAX_PRUNE_ITERS = 200 轮）:
    对当前候选联合计算特异位点与匹配率
    bad = { rep : 特异位点数 < 30  或  match_rate ≤ 0.9 }
    若 bad 为空 → 收敛，结束
    否则移除 bad 中 tsv_dist 最大者（每次只删 1 个）
```

每轮只删最差的一个，避免一次删掉过多后失去上下文；`0.9` 可通过 `-t` 调整。

#### (d) 最终输出

用最终候选重新计算特异位点与匹配率，写入：

| 文件 | 内容 |
|---|---|
| `final_representatives.txt` | 最终代表株名，一行一个 |
| `final_match_results.tsv` | `Strain / SpecificSites / Matched / MatchRate` |
| `final_specific_sites.tsv` | `strain / address / base`，每个代表株的特异位点明细 |
| `final_clusters.tsv` | `strain / cluster_id`，最终代表株编号 |

实测（3 样本）：阶段 2 纯脚本 1.16 s，单样本 hybrid 0.19–1.10 s。

---

## 5. 阶段 D：cluster 成员扩展与距离（App 流程）

前置条件：`YC-result/<样本>/Strain-result/final_representatives.txt` 必须存在（阶段 C 产物）。

对每个样本：

1. 读取 `final_representatives.txt`，去掉 `S<编号>_` 前缀，得到库中的菌株名。
2. 在 `output/rep_<db>.tsv` 第三列查到该代表株的 `cluster_id`；查不到时回退扫描
   `<db>.clu`（`菌株名,cluster_id`）。多库时按库分别处理。
3. 从 `<db>.clu` 取出该 cluster 的全部成员（不只是一株代表）。
4. 将成员名写入 patterns 文件，执行：

   ```bash
   mbgc d -Q -E <patterns> <db>.mbgc -
   ```

   `-Q` 关闭并行解码以保证输出顺序稳定，`-` 表示流式输出到 stdout。`-E` 本身是
   子串匹配，程序再按精确菌株名过滤，因此不会混入名字相近的其他菌株。
5. 逐条序列与 `output/really-ref.fa` 比较，距离定义与阶段 B 相同：
   双方均为 A/C/G/T 且不同的位点计数（`N`、`.`、`-` 一律忽略）。
6. 按 `poc.txt` 的位点截取序列，写入 `Strain-result/cluster_fa_<db>/<菌株名>`。

输出：

| 文件 | 内容 |
|---|---|
| `YC-result/<样本>/Strain-result/cluster_members.tsv` | `strain / db / cluster / distance`，按距离升序 |
| `YC-result/<样本>/Strain-result/cluster_fa_<db>/` | 各成员的 poc 截取序列 |
| `Strain-result/cluster_members.tsv`（工作目录根） | `sample / strain / db / cluster / distance` 汇总 |

> 阶段 D 不使用 `max_dist` 过滤：代表株所在 cluster 的成员全部保留，
> 便于观察簇内距离分布和发现离群成员。

实测（3 样本、1 库 Vp）：每个样本 5 个代表株 → 5 个 cluster → 4,168 个成员，
流式提取 + 距离计算约 53 s/样本；代表株的距离与 `rep_Vp.tsv` 完全一致。

---

## 6. 关键概念与阈值汇总

| 概念 | 定义 |
|---|---|
| **address** | `ref.fa` 上的碱基坐标；`add` 列表给出参与分析的 1,801,821 个坐标 |
| **really-ref** | 样本在每个 `add` 位点的单碱基判定结果，`-` 表示不可判定 |
| **poc 位点** | 样本中出现 ≥ 2 种可信碱基的位点，携带菌株区分信息 |
| **SNP 距离** | 双方都非 `.`（缺失）的位点上，碱基不同的数量 |
| **特异位点** | 某位点上某碱基在当前代表集合中只出现 1 次 → 该位点属于该代表 |
| **匹配率** | 特异位点中，真实 reads 支持数 > 0 的比例 |

| 参数 | 本流程取值 | 作用 |
|---|---|---|
| `error_base_number` | 深度 <100 时 0，否则 1 | 碱基支持数 ≤ 该值视为噪声 |
| `min_depth` | 4 | 深度 ≤ 4 不给出单碱基站 |
| `num_segments` | 深度 >500 时 20，否则 10 | pileup 分段并行数 |
| `pickX` | 30 | 每库先取距离最近的前 30 株 |
| `max_dist` | 默认 10000；深度 <50 时 20000 | 与 really-ref 距离上限，超出即排除 |
| `cutoff` | 1000 | 贪心聚类的簇内 SNP 距离上限（精炼阶段 UPGMA 阈值由距离分布自适应生成） |
| `MIN_SITES` | 30 | 每个代表株至少需要的特异位点数 |
| `MAX_CLUSTERS` | 50 | 单次切树允许的最大簇数 |
| `match threshold` | 0.9 | 匹配率下限（`<=` 即视为 bad） |
| `MAX_PRUNE_ITERS` | 200 | 全局剪枝最大轮数 |

---

## 7. 运行方式

### 7.1 本机（macOS arm64）

```bash
cd /path/to/StrainScale-mac

# 阶段 A（本机原生：sambamba 1.0.1 + mpileup-standalone）
bash SS.sh sample.list Vp-db

# 阶段 B + C（Linux x86_64 静态二进制，经 Docker 运行）
docker run --rm --platform linux/amd64 \
  -v "$PWD:/work" -w /work ubuntu:24.04 \
  bash /work/server-deploy/run_parallel.sh \
    /work/Vp-db /work/sample.list 30 /work/server-deploy/mbgc.zstd 36 1000 "" 10000 \
    /work/server-deploy/hybrid Strain-1k
```

### 7.2 常用变体

- **只跑阶段 B**：第 6 参传 `0`（`REPCUT=0` 跳过阶段 C）。
- **只跑阶段 C**：`SKIP_MBGC=1`，复用已有 `rep_*.tsv` 与 `fa_*`。
- **调整并发**：第 5 参 `MAXJOBS`（本机建议 ≤ 12；每进程约 1 GB 量级）。

> 本机阶段 A 的 `sambamba` 是源码编译的 macOS arm64 原生版本（`slice2.sh` 中指向
> `sambamba-master/bin/sambamba-1.0.1`），`mpileup` 为本机编译的 `mpileup-standalone`；
> 无需 Rosetta。阶段 B/C 使用 Linux x86_64 静态二进制，通过上面的 Docker 命令运行。
>
> 阶段 D 由桌面 App 在阶段 C 之后自动执行；`run_parallel.sh` 只覆盖阶段 B/C。

### 7.3 目录结构

```
YC-result/
└── <样本名>/
    ├── estimated_depth.txt                # 阶段 A 深度估计
    ├── filter-result-with-counts.txt      # 阶段 A（阶段 C 必需）
    ├── output/
    │   ├── 2lie-mp-samtool
    │   ├── really-ref.fa
    │   ├── ref-coords.txt
    │   ├── poc.txt / poc2.txt
    │   ├── rep_<db>.tsv                   # 阶段 B
    │   ├── fa_<db>/<菌株名>               # 阶段 B
    │   ├── representatives.txt/.fa        # 阶段 C 中间产物
    │   ├── representatives.dist.tsv
    │   ├── representatives.matrix.addr.tsv
    │   └── cluster_assignments.tsv
    └── Strain-1k/                          # 阶段 C 最终结果 + 阶段 D
        ├── final_representatives.txt
        ├── final_match_results.tsv
        ├── final_specific_sites.tsv
        ├── final_clusters.tsv
        ├── cluster_members.tsv             # 阶段 D：成员与 really-ref 距离
        └── cluster_fa_<db>/<菌株名>        # 阶段 D：成员的 poc 截取序列

Strain-result/cluster_members.tsv           # 阶段 D：全部样本汇总
```

---

## 8. 算法要点与注意事项

1. **两套坐标系、一次对齐**：BAM 与 `add` 用 `ref.fa` 坐标；`really-ref`、`G0.fa`、
   库序列都按 `add` 位点展开成等长向量，因此只有“同一 `add` 位点”上的碱基会被比较。
2. **缺失即忽略**：`-`/`.` 在距离计算中既不计差异也不计分母，避免把测序缺口当成突变。
3. **poc 位点是信息源**：样本无法确定单碱基的位点，反而最能区分库中菌株；
   阶段 B 按这些位点截取库序列，阶段 C 按这些位点找代表株的特异碱基。
4. **两层聚类**：先按 SNP 距离贪心粗聚（cutoff），再用 UPGMA 细分并用真实 reads 校验；
   最后全局剪枝，删除特异位点不足或匹配率过低的代表。
5. **结果确定性**：候选排序使用 `(距离, 名字)` 双键；UPGMA 合并时取最小索引，
   因此同一输入重复运行结果一致（与 Python 基线最多相差 1 株，来自平局排序）。
6. **阶段 C 的输入是“目录下所有 `rep_*.tsv`”**：若 `output/` 中残留其他运行产生
   `rep_*.tsv`（例如手工构造的 `rep_Vp-target.tsv`），会被一并纳入候选集，
   从而改变最终代表株。复现主流程前应确保只保留本次运行的 `rep_<db>.tsv`。
7. **复杂度**：贪心聚类约 `O(k·n·L)`（`k` 为代表数、`n` 为候选株数、`L` 为 poc 位点数，
   已用 Rayon 并行）；UPGMA `O(k³)`；剪枝最多 200 轮，每轮重算特异位点与匹配率。
