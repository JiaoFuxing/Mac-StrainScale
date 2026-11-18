# IPC API 参考

前端通过 `window.__TAURI__.core.invoke(<command>, { request } | {})` 调用。所有结构体字段在 JS 侧用 **camelCase**；Rust 侧用 snake_case 并通过 `serde(rename_all = "camelCase")` 映射。

> `withGlobalTauri: true` 是必需的（`tauri.conf.json` 已开启），否则 `window.__TAURI__` 不存在。

---

## 1. 状态与系统指标

### `pipeline_status`

```js
const progress = await invoke('pipeline_status');
```

```json
{
  "mode": "2lie-mp-samtool2Result",
  "step": 2,
  "totalSteps": 3,
  "percent": 66.7,
  "stage": "库检索",
  "message": "Stage B 运行中…",
  "done": false
}
```

未运行过流程时返回 `null`。

### `system_metrics`

```js
const metrics = await invoke('system_metrics');
```

```json
{
  "cpuPercent": 21.4,
  "memoryPercent": 52.1,
  "memoryUsedGb": 8.4,
  "memoryTotalGb": 16,
  "diskPercent": 36.2,
  "diskUsedGb": 412.8,
  "diskTotalGb": 1140.2
}
```

---

## 2. 路由 / 解析

### `agent_route`

```js
const decision = await invoke('agent_route', {
  request: { message: '帮我分析 /data/s1.bam 的菌株组成' }
});
```

```json
{
  "route": "skill",
  "reason": "检测到明确分析动作和具体输入文件，进入流程执行。",
  "files": ["/data/s1.bam"],
  "clarification": null
}
```

`route` 取值：

| 值 | 前端行为 |
|---|---|
| `chat` | 不执行任何流程，提示当前版本不支持自由对话，引导填写流程面板 |
| `skill` | 调 `plan_pipeline` 解析参数，再调 `run_agent_tool` 执行 |
| `clarify` | 直接显示 `clarification` |

### `plan_pipeline`

```js
const plan = await invoke('plan_pipeline', {
  request: { message: '输出到桌面 /Users/me/db /Users/me/ec/fastq.list' }
});
```

```json
{
  "mode": "2lie-mp-samtool2Result",
  "databaseDir": "/Users/me/db",
  "list": "/Users/me/ec/fastq.list",
  "outputDir": "/Users/me/Desktop",
  "reason": "未使用 LLM：优先读取 DB/list/output 标签，再按路径规则识别…"
}
```

缺少列表 / 数据库 / 无法判断模式时返回 `Err(String)`，前端直接展示中文原因。

### `run_agent_tool`

```js
const outcome = await invoke('run_agent_tool', {
  request: {
    name: 'strain_pipeline',
    arguments: {
      databaseDir: '/Users/me/db',
      list: '/Users/me/ec/fastq.list',
      outputDir: '/Users/me/Desktop'
    },
    message: '…'
  }
});
```

```json
{
  "mode": "2lie-mp-samtool2Result",
  "result": { "success": true, "exitCode": 0, "workingDir": "…", "outputs": ["…"], "errors": [] }
}
```

`name` 目前只支持 `strain_pipeline`，其它值返回 `Err("Skill 尚未接入执行器: …")`。

---

## 3. 流程

### `run_pipeline_mode`（右栏「开始运行」）

```js
await invoke('run_pipeline_mode', {
  request: {
    mode: '',                                  // 留空则根据列表内容自动判断
    databaseDir: '/Users/me/Ecoli-DB',
    list: '/Users/me/Ec/fastq.list',
    outputDir: '/Users/me/Desktop'             // 留空则回退到列表文件所在目录
  }
});
```

- `mode` 支持 `reads2Result`（从 fastq 起步）与 `2lie-mp-samtool2Result`（从 BAM / 2lie-mp-samtool 结果起步）。
- 返回 `PipelineResult`；`workingDir` 是最终使用的输出目录，目录不存在会自动创建。
- 执行期间通过 `pipeline-progress` 事件推送 `step/totalSteps/percent/stage/message`。

### `run_align_stage`（比对排序）

```js
await invoke('run_align_stage', {
  request: {
    fastqList: '/Users/me/Ec/fastq.list',
    reference: '/Users/me/Ec/ref.fa',
    outputDir: '/Users/me/Ec/hi',
    threads: 8,
    jobs: 4
  }
});
```

调用 `align_sort.sh`（sambamba + samtools），产出 `sample.list` 与 BAM。

### `run_strainscale_pipeline`（Stage A）

```js
await invoke('run_strainscale_pipeline', {
  request: {
    sampleList: '/tmp/ss-A.list',
    databaseDir: '/path/to/Vp-db',
    workingDir: '/path/to/workdir'
  }
});
```

```json
{
  "success": true,
  "exitCode": 0,
  "workingDir": "…/hi",
  "stdout": "基础分段长度: …\nBAM analysis is OK\n…",
  "stderr": "",
  "outputs": [
    "…/hi/YC-result/bwa-50-1/filter-result-with-counts.txt",
    "…/hi/YC-result/bwa-50-1/output/really-ref.fa",
    "…"
  ],
  "errors": []
}
```

### `run_mbgc_stage`（Stage B）

```js
await invoke('run_mbgc_stage', {
  request: {
    databaseDir: '…/Vp-db',
    samplesList: '/tmp/ss-A.list',
    workingDir: '…/hi',
    pickX: 30,
    maxJobs: 4,
    maxDist: 10000
  }
});
```

返回 `PipelineResult`。`stdout` 是各库日志的拼接（形如 `[Vp] exit=Some(0) ok=true` + 完整 stdout/stderr）。`exitCode` 在全部库成功时为 0，否则 1（即使某个库失败，其它库的产物仍会出现在 `outputs`）。

### `run_hybrid_stage`（Stage C 单样本；界面未暴露）

```js
await invoke('run_hybrid_stage', {
  request: {
    sampleOutputDir: '…/hi/YC-result/bwa-50-1/output',
    filterCounts: '…/hi/YC-result/bwa-50-1/filter-result-with-counts.txt',
    resultDir: '…/hi/YC-result/bwa-50-1/Strain-result',
    cutoff: 1000,
    maxDist: 10000,
    matchThreshold: 0.9
  }
});
```

### `run_hybrid_batch`（Stage C 批量；界面按钮）

```js
await invoke('run_hybrid_batch', {
  request: {
    workingDir: '…/hi',
    samplesList: '/tmp/ss-A.list',
    cutoff: 1000,
    maxDist: 10000,
    matchThreshold: 0.9
  }
});
```

```json
{
  "success": true,
  "ran": ["bwa-50-1", "bwa-10-9", "bwa-10-10"],
  "skipped": [],
  "failed": [],
  "outputs": [
    "…/Strain-result/final_representatives.txt",
    "…/Strain-result/final_match_results.tsv",
    "…/Strain-result/final_specific_sites.tsv",
    "…/Strain-result/final_clusters.tsv",
    "…"
  ]
}
```

跳过样本的格式：`"<name>: 缺少 output/（Stage A 未完成）"`。失败样本：`"<name>: <错误1>; <错误2>"`。

---

## 4. 错误约定

| 层级 | 表现 |
|---|---|
| 命令返回 `Err(String)` | `invoke` 的 Promise reject，前端 toast 显示字符串 |
| `PipelineResult.success=false` | 业务失败，但 HTTP/IPC 成功；前端渲染 `errors` 列表，并保留 `stdout/stderr/outputs` |
| `RouteResult.route="clarify"` | 不是错误，是「需要用户补输入」的正常分支 |

常见顶层错误：

- `plan_pipeline`：`没有在消息中找到样本列表文件…`、`没有在消息中找到参考数据库目录…`、`无法从列表文件内容判断起点模式…`
- `run_agent_tool`：`Skill 尚未接入执行器: <name>`、`缺少 database_dir（参考数据库目录）`
- 流程命令：`pipeline task failed: …`（`spawn_blocking` join 失败，通常意味着 panic）

---

## 5. 前端最小调用骨架

```js
const tauri = window.__TAURI__;
const invoke = tauri && tauri.core && tauri.core.invoke;

async function runStageA(form) {
  const result = await invoke('run_strainscale_pipeline', {
    request: {
      sampleList: form.sampleList,
      databaseDir: form.databaseDir,
      workingDir: form.workingDir
    }
  });
  renderStageResult(resultEl, result, 'Stage A 执行完成', 'Stage A 未通过校验');
}
```

完整实现见 `src/index.html` 的 `submitTask`、`runPipelineFromForm`、`runMbgcFromForm`、`runHybridFromForm`、`renderStageResult`、`renderBatchResult`。

---
