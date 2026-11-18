# 架构说明

> 对应代码：`strainscale-desktop/src-tauri/src/lib.rs`（2655 行）、`strainscale-desktop/src/index.html`（1699 行）
> 行号会随改动漂移，本文档以函数名和调用关系为准，行号只作定位参考。

---

## 1. 分层

```
┌──────────────────────────────────────────────────────────────┐
│ WebView（src/index.html，单文件，无构建步骤）                  │
│  左栏：会话列表（演示数据）  中栏：对话  右栏：流程面板         │
│  IPC：window.__TAURI__.core.invoke(...)                      │
└───────────────▲───────────────────────────┬──────────────────┘
                │ JSON                       │ invoke
┌───────────────┴───────────────────────────▼──────────────────┐
│ Rust 后端（lib.rs）                                           │
│  ┌──────────────────┐  ┌──────────────────────────────────┐  │
│  │ 规则路由/路径解析 │  │ 流程执行器 A/B/C                  │  │
│  │ agent_route      │  │ Command + 校验 + 产物收集         │  │
│  │ plan_pipeline    │  │ 进度事件 pipeline-progress        │  │
│  └──────────────────┘  └──────────────────┬───────────────┘  │
└───────────────────────────────────────────┼──────────────────┘
                                            │ std::process::Command
                              ┌─────────────▼──────────────┐
                              │ /bin/bash SS.sh            │
                              │ mbgc-macos-arm64 p         │
                              │ hybrid-macos-arm64         │
                              └────────────────────────────┘
```

- 没有数据库、没有网络服务、没有账号体系；所有状态在进程内或工作目录里。
- 前端不做任何计算，只负责收集表单、调用命令、渲染 `PipelineResult`。
- 后端是唯一的执行入口：所有路径校验、超时、日志、产物收集都在 Rust 侧。

---

## 2. 进程模型

### 2.1 应用进程

Tauri 主进程只持有一个进度状态：

```rust
struct PipelineState(Mutex<Option<PipelineProgress>>);
```

- 流程执行过程中通过 `emit_pipeline_progress` 同时更新内存状态并向 WebView 发 `pipeline-progress` 事件，前端据此绘制进度条与阶段文案。
- `pipeline_status` 供界面重建时读取最近一次进度。
- 应用**不再启动任何 sidecar 进程**，没有模型端口、没有 `/health` 轮询、没有退出时 kill 子进程的逻辑。

### 2.2 流程子进程

- 全部用同步 `std::process::Command::output()`，外层包在 `tauri::async_runtime::spawn_blocking`，不会阻塞 WebView 的 IPC 线程。
- Stage A/B/C 都是**一次调用跑完**，没有流式进度；界面在请求期间只显示「运行中」，结束后一次性渲染结果。
- Stage B 会为每个库开一个线程（`std::thread::spawn`），按 `max_jobs` 分块并行，最后 `join` 汇总日志。

---

## 3. IPC 命令（11 个，全部注册在 `invoke_handler`）

| 命令 | 入参 | 返回 | 说明 |
|---|---|---|---|
| `pipeline_status` | 无 | `Option<PipelineProgress>` | 最近一次流程进度 |
| `system_metrics` | 无 | `SystemMetrics` | CPU / 内存 / 磁盘占用 |
| `agent_route` | `{ request: AgentCallRequest }` | `RouteResult` | 纯规则分类，判定 `chat` / `skill` / `clarify` |
| `plan_pipeline` | `{ request: AgentCallRequest }` | `PipelinePlan` | 从消息中解析模式、数据库、列表、输出目录 |
| `run_agent_tool` | `{ request: AgentToolRequest }` | `AgentToolOutcome` | 按解析结果直接执行对应流程 |
| `run_strainscale_pipeline` | `{ request: PipelineRequest }` | `PipelineResult` | Stage A |
| `run_mbgc_stage` | `{ request: MbgcStageRequest }` | `PipelineResult` | Stage B |
| `run_hybrid_stage` | `{ request: HybridStageRequest }` | `PipelineResult` | Stage C 单样本（界面已不暴露，命令保留） |
| `run_hybrid_batch` | `{ request: HybridBatchRequest }` | `HybridBatchResult` | Stage C 批量（界面按钮走这个） |
| `run_align_stage` | `{ request: AlignStageRequest }` | `PipelineResult` | 比对排序（align_sort.sh） |
| `run_pipeline_mode` | `{ request: PipelineModeRequest }` | `PipelineResult` | 按 `2lie-mp-samtool2Result` / `reads2Result` 模式一次跑完，右栏「开始运行」走这个 |

字段命名统一 `serde(rename_all = "camelCase")`：Rust 侧 `sample_list` ↔ 前端 `sampleList`。

---

## 4. 数据结构

```rust
AgentCallRequest   { message: String }
RouteResult        { route: "chat"|"skill"|"clarify", reason: String, files: Vec<String>, clarification: Option<String> }
PipelinePlan       { mode, databaseDir, list, outputDir, reason }

PipelineRequest    { sampleList, databaseDir, workingDir }                       // Stage A
MbgcStageRequest   { databaseDir, samplesList, workingDir,
                     pickX: Option<u32>, maxJobs: Option<u32>, maxDist: Option<u32> } // Stage B
HybridBatchRequest { workingDir, samplesList, cutoff, maxDist, matchThreshold }   // Stage C 批量
HybridStageRequest { sampleOutputDir, filterCounts, resultDir, cutoff, maxDist, matchThreshold }
AlignStageRequest  { fastqList, reference, outputDir, threads, jobs }
PipelineModeRequest{ mode, databaseDir, list, outputDir, threads, jobs }          // outputDir 可留空

PipelineResult     { success, exitCode: Option<i32>, workingDir,
                     stdout, stderr, outputs: Vec<String>, errors: Vec<String> }
HybridBatchResult  { success, ran: Vec<String>, skipped: Vec<String>,
                     failed: Vec<String>, outputs: Vec<String> }
```

`PipelineResult` 是统一结果外壳：**即使失败也会带回已经产生的 stdout/stderr 和 outputs**，前端据此渲染「部分成功 + 错误列表」，而不是只弹一句失败。

---

## 5. 规则路由与路径解析（无模型）

### 5.1 `classify` / `agent_route`

纯字符串规则，目的是避免把非分析请求误送进流程执行器：

1. 命中「分析动作词」（`质控/运行/分析/识别/比对/丰度/导出/qc/run/analyze/...`）**且**能提取到具体文件（带扩展名或含 `/` 的文件名）→ `skill`。
2. 命中动作词但没有具体文件 → `clarify`，返回固定提示「需要你先选择或提供真实输入文件…」，前端直接显示。
3. 命中「问句词」（`是什么/为什么/解释/介绍/what is/why/...`）→ `chat`，前端提示当前版本不支持自由对话。
4. 只提取到文件、没有动作词 → `skill`。
5. 都不命中 → `chat`。

文件提取 `extract_files` 按空白和中文标点切 token，`looks_like_file_reference` 判断扩展名白名单（fastq/fasta/bam/vcf/csv/tsv/json/xlsx/pdf/sh/py 等）或路径末段含 `.`。

### 5.2 `plan_pipeline`

`plan_pipeline_from_message` 把消息解析成 `PipelinePlan`：

1. 先读显式标签：`list/sample_list/样本列表/列表`、`database/database_dir/db/数据库`、`output/out/输出目录/输出`。
2. 再按路径特征归类：以 `.list` 结尾或能推断模式的文件 → 列表；同时包含 `ref.fa` 与 `add` 的目录 → 数据库；其余目录/路径 → 输出目录。
3. 消息里出现「桌面 / desktop」且没有显式输出目录时，默认输出到 `~/Desktop`。
4. `infer_pipeline_mode` 根据列表文件内容判断起点：`fastq` → `reads2Result`；BAM/`2lie-mp-samtool` → `2lie-mp-samtool2Result`。
5. 缺少列表、数据库或模式时返回明确的中文错误，前端逐条展示。

> 已知弱点：规则是关键词匹配，中文「帮我看看这个」+ 文件路径会走 `skill`；反过来「分析一下我的数据」没有路径会走 `clarify`。这是刻意保守的选择，宁可多问一句，也不猜路径。

---

## 6. 资源解析（跨平台移植的关键点）

`pipeline_resource_dir(app)`：

```
1) app.path().resource_dir()? + "resources/strainscale-pipeline/StrainScale-agent-package"
2) 若不存在 → env!("CARGO_MANIFEST_DIR") + "resources/strainscale-pipeline/StrainScale-agent-package"
3) 都不存在 → 返回错误 "bundled StrainScale pipeline directory does not exist"
```

打包后实际落盘是：

```
StrainScale Agent.app/Contents/Resources/resources/strainscale-pipeline/StrainScale-agent-package/
```

即 **`Contents/Resources` 下还有一层 `resources/`** —— 因为 `tauri.conf.json` 的 `resources` 写的是 `resources/...` 相对路径，Tauri 会保留目录层级。代码里的拼接正好对上；如果在新平台改了资源路径，这一层要同步改，否则运行时找不到脚本。

---

## 7. 错误模型

分三层：

1. **参数/前置校验失败**：返回 `PipelineResult { success:false, exitCode:null, errors:[...] }`，**不启动子进程**。错误文案带完整路径和原因（例如 `missing really-ref.fa for sample bwa-10-9; run Stage A first`）。
2. **子进程失败**：仍返回 `exitCode` 与 `stdout/stderr`，`errors` 里加一条 `SS.sh exited with status: ...` / `mbgc failed for database(s): ...`。
3. **后端 panic / 任务失败**：命令返回 `Err(String)`，前端 catch 后 toast。

前端 `renderStageResult` 统一渲染：成功显示 `exitCode` + 产物列表；失败显示错误条目，仍保留已生成的产物，便于判断「跑到哪一步断的」。

---

## 8. 安全边界（现状，不是承诺）

- `capabilities/default.json` 只有 `core:default`，没有 fs/shell 的前端权限；所有文件与进程操作都在 Rust 侧，前端只能调上面 11 个命令。
- 执行器只调用**固定路径**的二进制（`SS.sh`、`mbgc-macos-arm64`、`hybrid-macos-arm64`），不接收用户传入的可执行文件路径。
- 没有真正的沙箱：子进程继承当前用户权限，工作目录由用户指定；Stage A 的 `SS.sh` 会 `cd` 进 `YC-result/<样本>`，路径来自 `sample.list` 第一列（已禁止 `/` 和 `..`）。
- 没有网络隔离机制可依赖：应用本身不联网，也没有任何 HTTP 客户端依赖（`ureq` 已随模型一起移除）。
- 没有资源配额：Stage B 每个 mbgc 进程约 1 GB 内存，默认并发 4、上限 12，需要按机器内存调。

---
