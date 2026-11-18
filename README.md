# StrainScale Agent 桌面版 — 源码与交接文档

> 版本 **0.5.0** · 目标平台 **macOS arm64（Apple Silicon）** · 框架 **Tauri 2**
> 本文档面向接手开发的同学：读完这一份就能把工程跑起来、知道每个文件在干什么、知道哪里还没做完。

---

## 1. 这是一个什么程序

一个把 StrainScale 生信流程（BAM → 样本参考 → 库检索 → 代表株聚类）包成本地桌面工具的原型：

- **界面**：一个静态 HTML 页面（`src/index.html`），三栏布局：会话 / 对话 / 流程面板。没有 React、没有打包器，Tauri 直接把 `src/` 当 `frontendDist` 加载。
- **后端**：Rust（`src-tauri/src/lib.rs`），负责按规则做意图路由与路径解析、调用原生流程二进制、做输入校验并把结果以 JSON 返回界面。
- **无内嵌模型**：0.5.0 已移除 llama.cpp sidecar 与 GGUF 模型，应用纯规则解析（`agent_route` + `plan_pipeline`），不联网、不上传数据。
- **流程执行**：三段式，全部在本机执行，不依赖远程服务器、不依赖 Docker：

| 阶段 | 真正执行的程序 | 作用 |
|---|---|---|
| Stage A | `SS.sh` → `slice2.sh`（sambamba + mpileup）+ `get_site.pl` | BAM 切片、pileup、位点判定，产出 `really-ref.fa` 等 |
| Stage B | 原生 `mbgc-macos-arm64 p`（MBGC v2.1.6） | 用库对每个样本打分，取前 N 个候选株 |
| Stage C | 原生 `hybrid-macos-arm64` | 候选株聚类、距离矩阵、代表株精炼 |

界面右侧「StrainScale 流程」面板就是这三个按钮。Stage B/C 的样本列表统一沿用 Stage A 的 `sample.list`，Stage C 自动跑列表里的全部样本。

---

## 2. 五分钟跑起来（macOS）

### 2.1 环境

| 依赖 | 版本 / 说明 |
|---|---|
| macOS | 13+，Apple Silicon（arm64）。Intel Mac 需要自行重编原生二进制 |
| Xcode Command Line Tools | `xcode-select --install`（提供 clang 与 codesign） |
| Node.js | 18+（只用于跑 `@tauri-apps/cli`，前端无构建步骤） |
| Rust | stable（`rustup` 安装），需 `aarch64-apple-darwin` target |
| libomp | 仅重编 mbgc 时需要；`brew install libomp` 或复用 conda 环境 |

### 2.2 启动开发模式

```bash
cd desktop
npm install          # 只装 @tauri-apps/cli
npx tauri dev        # 首次会编译 Rust，约 2–5 分钟
```

> 注意：`package.json` 里的脚本名是 `dev` / `build`，所以 `npm run tauri build` 会报 `Missing script: "tauri"`。请用 `npx tauri dev` / `npx tauri build`。

### 2.3 打包

```bash
npx tauri build
# 产物：src-tauri/target/release/bundle/macos/StrainScale Agent.app
```

`tauri.conf.json` 的 `bundle.targets` 只有 `["app"]`，**不会自动出 DMG**。打包后必须重新 ad-hoc 签名，再手工做 DMG，完整命令见 `docs/BUILD.md`。

---

## 3. 目录结构

```
StrainScale-Agent/
├── README.md
├── .gitignore
├── desktop/                              # Tauri 桌面应用
│   ├── package.json                     # 只有 dev/build 两个脚本 + tauri CLI
│   ├── package-lock.json
│   ├── src/
│   │   └── index.html                   # 全部界面、样式、前端逻辑（单文件）
│   └── src-tauri/
│       ├── Cargo.toml                   # crate: strainscale-agent 0.5.0
│       ├── Cargo.lock
│       ├── build.rs                     # tauri_build::build()
│       ├── tauri.conf.json              # 窗口 / 资源 / 图标
│       ├── capabilities/default.json    # 仅 core:default（最小权限）
│       ├── icons/                       # 含 icon.icns / icon.ico，可跨平台复用
│       ├── resources/
│       │   └── strainscale-pipeline/    # 运行时完整包（脚本 + 二进制）
│       └── src/
│           ├── main.rs                  # 仅调 strainscale_agent_lib::run()
│           └── lib.rs                   # 全部后端逻辑
├── pipeline/                             # 流水线脚本
│   ├── SS.sh  slice2.sh  get_site.pl  align_sort.sh
│   ├── PIPELINE.md                      # 算法流程与参数说明
│   └── server-deploy/
│       ├── run_parallel.sh  build-hybrid.sh  README.md
│       └── hybrid-rust/                 # hybrid 源码（Rust）
└── mbgc-source/mbgc-zstd-prestrain/     # mbgc 源码 + vendor/zstd + 构建脚本
```

**包含** macOS arm64 运行时二进制（`sambamba`、`mpileup`、`mbgc`、`hybrid`、`libomp.dylib`），clone 后可直接 `cargo check` / `tauri build`。**不含** `target/`、`node_modules/`、`tools/`。

---

## 4. 已实现 / 未实现

**已实现并实测**

- Stage A/B/C 三个按钮，真实调用本机二进制，返回 `exitCode`、`stdout/stderr`、产物路径、逐条错误。
- 输入校验：文件存在性、`sample.list` 两列格式、样本名合法性、数据库结构、上游产物是否齐备。
- 纯规则任务解析：`agent_route` 判定 `chat` / `skill` / `clarify`，`plan_pipeline` 从自然语言里解析数据库目录、列表文件与输出目录，再交给流程执行器。
- 右侧流程面板可直接填写数据库目录、列表文件、输出目录后一键运行（输出目录留空则回退到列表文件所在目录）。
- 中文输入法回车误发送已修复（`compositionstart/end` + `isComposing || keyCode === 229`）。

**未实现 / 半成品（接手后优先处理）**

- 左栏「最近会话」「资源占用」等仍有部分静态演示数据，不是真实状态。
- 没有任务队列、取消、超时、并发上限；Stage B 的 1 GB/进程内存占用也没有按机器档位自适应。
- 没有单元测试、集成测试、CI；没有 git 历史（当前目录不是 git 仓库）。
- Windows / Linux 未构建、未验证；`/bin/bash`、`dylib`、`#[cfg(unix)] symlink` 都需要改（清单见 `docs/BUILD.md`）。
- 分发只做了 ad-hoc 签名；给别人用需要 Developer ID + notarization。

---

## 5. 文档索引

| 文件 | 内容 |
|---|---|
| `docs/ARCHITECTURE.md` | 分层结构、进程模型、11 个 IPC 命令、数据结构、路由与解析规则、资源解析、错误模型 |
| `docs/PIPELINE.md` | Stage A/B/C 的输入约定、真实命令、产物、校验规则、常见报错对照表 |
| `docs/API.md` | 每个命令的请求 / 响应 JSON 示例与字段说明 |
| `docs/BUILD.md` | macOS 打包、签名、DMG、原生二进制重编、Windows 移植清单 |
| `docs/HANDOFF.md` | 现状与验证证据、已知坑、技术债、建议的下一步、验收清单 |

算法层面的细节（切片长度、位点判定、UPGMA、递归细化、阈值）不在本文档重复，仓库内 `pipeline/PIPELINE.md` 是上游原始说明，可直接对照。

---

## 6. 源码包内容（本仓库）

仓库结构见第 3 节。含 macOS arm64 运行时二进制，clone 后可直接构建；只排除了 `target/`、`node_modules/`、`tools/`。

快速自检：

```bash
cd desktop && npm install
cd src-tauri && cargo check          # 只检查 Rust 代码，不产出二进制
cd ../../pipeline && bash -n SS.sh slice2.sh server-deploy/run_parallel.sh
```

> **从 GitHub 下载 zip 后**：macOS 会给文件加 `com.apple.quarantine` 属性，运行时可能提示"无法验证开发者"或二进制被拒绝执行。先清除属性并补执行权限：
> ```bash
> xattr -dr com.apple.quarantine StrainScale-Agent
> chmod +x desktop/src-tauri/resources/strainscale-pipeline/StrainScale-agent-package/server-deploy/{mbgc-macos-arm64,hybrid-macos-arm64,rammap-macos-arm64} \
>          desktop/src-tauri/resources/strainscale-pipeline/StrainScale-agent-package/sambamba-master/bin/sambamba-1.0.1 \
>          desktop/src-tauri/resources/strainscale-pipeline/StrainScale-agent-package/mpileup-standalone/mpileup
> ```

---

## 7. 交接时必须一起带走的东西

源码包里已包含：`src/`、`src-tauri/src`、配置、图标、`pipeline` 脚本与 `hybrid-rust` 源码、mbgc 源码与构建脚本、本文档。

**没有**包含（体积或授权原因）：

| 内容 | 原因 | 获取方式 |
|---|---|---|
| `target/`、`node_modules/` | 编译缓存 / 依赖，可重新生成 | `npm install` + `npx tauri build` |
| DMG / `.app` | 构建产物 | 自行 `npx tauri build` |

第三方组件（MBGC、hybrid、sambamba、mpileup、zstd/libdeflate 等）各自遵循其上游许可；对外分发前请逐一确认许可与署名要求。
