# 交接说明：现状、坑、下一步

> 写给接手的人。目标不是复述功能，而是让你知道**哪些能信、哪些不能信、先改哪里**。

---

## 1. 当前状态（0.5.0）

**能跑的部分（已用真实数据验证）**

- `npx tauri build` 可产出可运行的 `.app`；ad-hoc 签名 + 手工 DMG 流程已跑通。
- Stage A：BAM 列表 → `really-ref.fa` + poc 文件。
- Stage B：原生 mbgc 打分 → `rep_<db>.tsv` + `fa_<db>`；3 个样本、1 个库（Vp）实测 `exitCode 0`。
- Stage C：批量跑全部样本 → 每个样本 4 个 `final_*` 文件；实测 `ran: 3 / skipped: 0 / failed: 0`，共 12 个产物。
- 本地模型已移除：0.5.0 不再启动 llama.cpp sidecar，任务解析改走 `agent_route` + `plan_pipeline` 规则，`run_agent_tool` / `run_pipeline_mode` 可直接执行。
- 右侧流程面板支持填写数据库目录、列表文件、输出目录后直接运行。
- 输入法回车误发送已修。

**还不能当成品用的部分**

- 规则解析是关键词 + 路径启发式，无法理解复杂自然语言意图；识别不到路径时只会要求补充输入。
- 左栏会话、进度、资源占用是静态演示数据。
- 没有任务队列、取消、超时、失败重试；关掉窗口不会杀正在跑的 mbgc/hybrid。
- 只有 macOS arm64 构建；Intel Mac、Windows、Linux 均未验证。
- 没有测试、没有 CI、没有 git 历史（交付目录不是仓库，第一次提交由你来做）。

---

## 2. 已验证的数据与路径（可复现）

| 项 | 值 |
|---|---|
| 工作目录 | `/path/to/workdir` |
| 数据库 | `/path/to/Vp-db`（`G0.fa` + `DB/{Vp.zst,Vp.clu,Vp.mbgc}`） |
| 样本 | `bwa-50-1`、`bwa-10-9`、`bwa-10-10` |
| Stage A 产物 | `hi/YC-result/<样本>/{filter-result-with-counts.txt,output/...}` |
| Stage B 产物 | `hi/YC-result/<样本>/output/{rep_Vp.tsv,fa_Vp/}` |
| Stage C 产物 | `hi/YC-result/<样本>/Strain-result/final_*.{txt,tsv}` |
| 样例表 | `/tmp/ss-A.list`（两列 TSV，可重建） |

实测 Stage B 日志关键行：

```
[Vp] exit=Some(0) ok=true
sample bwa-10-9: ref_num=12696, sites=1672817
sample bwa-50-1: ref_num=22, sites=1681107
sample bwa-10-10: ref_num=12803, sites=1673429
kept=5 / kept=4 / kept=25
```

> 这些路径是**原开发机**的，换机器后要重新准备数据；不要把它们硬编码进代码或测试。

---

## 3. 已经踩过、不要再踩的坑

### 3.1 暂存目录删掉数据库（最严重）

- 症状：Stage B `exitCode 1`，stderr `G0 not found: .../.mbgc-stage/Vp/G0.fa`；`ls` 却能看到文件。
- 原因：暂存根建在 `<工作目录>/.mbgc-stage`，数据库正好也在其下；每次 `remove_dir_all(stage_root)` 把 `G0.fa`/`DB/` 删掉，重建后的符号链接指向自己。
- 现状：已改成独立 `<工作目录>/.mbgc-run`，并加 `database_dir.starts_with(&stage_root)` 兜底到 `/tmp/mbgc-run-<pid>`。
- 教训：**凡是「先校验、再清理、再执行」的流程，都要预判清理会动到输入的父目录**。移植时不要为了「目录整齐」把暂存根改回通用名字。

### 3.2 目录改名导致 Stage C 假失败

用户把测试目录 `test` 改名成 `hi` 后，Stage C 报 `sample output directory does not exist`。这不是代码 bug，是路径变了。排查顺序永远是：**先看 Stage A 产物是否存在，再看 Stage C**。

### 3.3 Linux ELF 在 macOS 不能跑

`server-deploy/mbgc.zstd`、`server-deploy/hybrid` 是 Linux x86-64 静态 ELF，Rosetta 只转 x86 Mach-O，不转 ELF。macOS 上必须用 `mbgc-macos-arm64` / `hybrid-macos-arm64`；包内保留 Linux 版只是留档。

### 3.4 模型位置与 TCC（已随模型移除而消失）

0.5.0 之前模型放 `~/Desktop` 下 sidecar 读不到、返回 503。现在应用不再加载模型，这条不再适用；如需重新引入本地模型，记得放 `~/Library/Application Support/com.strainscale.agent/` 或直接由后端返回 `app_data_dir()`。

### 3.5 打包与签名

- `npm run tauri build` 报 `Missing script: "tauri"` → 用 `npx tauri build`。
- `tauri.conf.json` 只有 `targets: ["app"]` → DMG 必须手工 `hdiutil`。
- 构建后必须 `codesign --force --deep --sign -`，否则 `--verify` 失败。
- `.app` 内资源路径是 **`Contents/Resources/resources/...` 双层**，手动搬资源时别把这一层拍平。
- 只跑 `cargo build --release` 后要手动替换 `Contents/MacOS/strainscale-agent`，否则应用包里可能还是旧二进制（表现为界面与源码不一致）。

### 3.6 自动化操作 WebView（如果你要写 UI 测试）

- `set_value` 不触发前端 `input` 事件，必须 `click + type_text`。
- AX 元素 index 每次观测都会重排，只能用最新一次树里的 index。
- 控件在可视区外时，先点输入框再 `Page_Down`。
- 用户正在操作时会报 `MAC_GUI_ERROR_7_0`，等待后重试。
- Tauri 首启 WebView 需要 3–5 秒才出表单。

---

## 4. 技术债与风险（按优先级）

| 优先级 | 问题 | 影响 | 建议 |
|---|---|---|---|
| P0 | 规则解析覆盖有限 | 复杂表述识别不到路径 | 扩充 `labeled_path` 标签与 `extract_files` 规则，并补单元测试 |
| P0 | 无超时 / 取消 | Stage B 卡住只能杀应用 | `Command` 加超时线程 + `kill`；或改成 `tauri-plugin-shell` 的流式事件，支持 cancel |
| P1 | 仅 macOS arm64 | 无法交付给 Windows/Linux 用户 | 按 `docs/BUILD.md` 第 3 节做平台抽象；先做 Windows（用户群体最大）或 Linux（服务器侧） |
| P1 | 无测试 / CI | 每次改动靠手点 | 至少给 `validate_pipeline_request`、`classify`、`plan_pipeline` 写单元测试；CI 跑 `cargo test` + `npx tauri build` |
| P1 | 无 git 历史 | 无法回溯、无法协作 | 第一次提交前先加 `.gitignore`（已随包），不要把 `target/`、`node_modules/`、模型、DMG 提上去 |
| P2 | Stage B 并发/内存固定 | 32 GB 机器默认 4，16 GB 机器会吃紧 | 按 `sysctl hw.memsize` 或配置项自适应，上限仍保留 12 |
| P2 | `run_hybrid_stage` 命令未暴露 | 单样本调试不方便 | 要么在界面加「只跑当前样本」，要么从命令表删掉，减少维护面 |
| P2 | 演示数据混在真实界面 | 用户会误以为「已有历史任务」 | 把左栏会话、进度条、资源占用标记为 `demo`，或接真实状态 |
| P3 | 知识图谱/本体层只有 UI 示意 | 无数据层 | 明确从产品范围里拿掉，或单独立项 |

---

## 5. 建议的接手顺序

1. **先把仓库建起来**：`git init`，用随包 `.gitignore`，首次提交只含源码 + 资源 + 文档。
2. **跑通一次完整链路**：按 README 2.2 启动，用 `/tmp/ss-A.list` + `Vp-db` 跑 A→B→C，确认 `exitCode 0`。
3. **补最小测试**（P1）：先测校验、路由与 `plan_pipeline`，这三块是纯函数，最容易回归。
4. **做平台抽象**（P1）：把二进制名、bash、symlink、libomp 四处收口到一个 `platform` 模块，再开 Windows 分支。
5. **扩规则解析**（P0）：按真实用户表述补充标签与路径模式，避免「说了但没识别到」。
6. **接执行闭环**（P0 产品项）：解析 → 确认 → 执行 → 把 `PipelineResult` 回填对话。

---

## 6. 验收清单（改完代码照着跑）

```bash
# 1) 构建
cd strainscale-desktop
npx tauri build
codesign --force --deep --sign - "src-tauri/target/release/bundle/macos/StrainScale Agent.app"
codesign --verify --deep --strict "src-tauri/target/release/bundle/macos/StrainScale Agent.app"

# 2) 启动
open "src-tauri/target/release/bundle/macos/StrainScale Agent.app"

# 3) UI 检查
#    - 右栏只有「StrainScale 流程」：数据库目录、列表文件、输出目录 + 开始运行
#    - 侧栏底部只有「分析内核 / 网络访问」，没有模型运行时、启动模型
#    - 中文输入法回车只上屏、不发送

# 4) 流程（用你自己的数据替换）
#    数据库目录 + 列表文件 + 输出目录 → 开始运行
#    输出目录留空时回落到列表文件所在目录，目录不存在会自动创建

# 5) 规则解析
#    输入「分析 /path/x.bam」但路径不存在 → 应进入 clarify，而不是编造结果
```

判定标准：

- Stage A/B/C 的 `errors` 为空、`outputs` 非空；
- 失败时能看到**具体路径和原因**，不是一句「执行失败」；
- 关掉应用后 `ps aux | grep -i strainscale` 只剩 0 个进程；
- 重跑 Stage B 不会删掉数据库目录（检查 `Vp-db/G0.fa` 与 `DB/` 仍在）。

---

## 7. 版本历史

| 版本 | 大小 | 主要变化 |
|---|---|---|
| 0.1.0 | 3.3 MB | Tauri 壳 + 静态界面 |
| 0.2.0 | 12.1 MB | 接入 llama sidecar 与对话（此阶段知识问答曾被误路由成工具卡） |
| 0.3.0 | 12.1 MB | 三段式路由 + 参数守卫 |
| 0.4.0 | 18.2 MB | 接入 Stage A（SS.sh）与原生 mbgc/hybrid |
| 0.5.0 | 19.27 MB | Stage B/C 输入精简（共用 Stage A 的 sample.list）；修复 Stage B 暂存目录删除数据库的 bug；重新签名打包 |
| 0.5.0（去模型） | 约 26 MB | 移除 llama.cpp sidecar、GGUF 模型与推理命令；改为纯规则路由 + 路径解析；流程面板新增输出目录 |

---

## 8. 交接边界

- 本包**不含** DMG、`target/`、`node_modules/`，原因与获取方式见 README 第 6 节；GGUF 模型与 llama.cpp 已从项目中移除。
- 原开发机上的 `Vp-db`、`hi/YC-result` 是测试数据，**不要**打进仓库（体积大且含样本路径）。
- 第三方二进制（mbgc、hybrid、sambamba、mpileup）的许可需在对外分发前逐项确认；本包只做技术交接，不替代许可审查。
- 本包没有 git 历史，第一次提交即项目的初始 commit。

---
