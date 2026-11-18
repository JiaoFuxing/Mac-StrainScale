# 构建与打包（含跨平台移植）

> 目标：拿到源码后能在新机器上从零构建出 `.app` / DMG，并知道换到 Windows/Linux 时需要改哪些地方。
> 当前已验证的构建环境：macOS（Apple Silicon）+ Node 18+ + Rust stable + Xcode CLT。

---

## 1. macOS 构建

### 1.1 首次准备

```bash
xcode-select --install                 # clang / codesign
rustup target add aarch64-apple-darwin
cd strainscale-desktop
npm install                            # 只装 @tauri-apps/cli
```

### 1.2 开发模式

```bash
npx tauri dev
```

首次会编译 Rust 依赖（约 2–5 分钟）。窗口尺寸 1440×900，最小 1080×700。WebView 首次渲染表单需要 3–5 秒，自动化测试时别在这一步误判「界面没出来」。

### 1.3 正式构建

```bash
npx tauri build
# => src-tauri/target/release/bundle/macos/StrainScale Agent.app
```

`tauri.conf.json` 的 `bundle.targets` 只有 `["app"]`，所以**不会**生成 DMG；需要手工做（1.5）。

### 1.4 必须重新签名

`tauri build` 产出的 `.app` 往往只有 linker-signed 的临时签名，`codesign --verify --deep --strict` 会失败，Gatekeeper 也会拦。构建后固定执行：

```bash
APP="src-tauri/target/release/bundle/macos/StrainScale Agent.app"
codesign --force --deep --sign - "$APP"
codesign --verify --deep --strict --verbose=2 "$APP"   # 必须返回 0
```

ad-hoc 签名足够本机运行和自用分发；**给别人用**需要 Apple Developer ID + `notarytool` 公证 + `stapler`，否则对方首次打开会被 Gatekeeper 拦。

### 1.5 手工打 DMG

```bash
ROOT="$(pwd)"
STG=$(mktemp -d)
cp -R "src-tauri/target/release/bundle/macos/StrainScale Agent.app" "$STG/"
ln -s /Applications "$STG/Applications"

hdiutil create -volname "StrainScale Agent" \
  -srcfolder "$STG" -ov -format UDZO \
  "$ROOT/../StrainScale-Agent-0.5.0.dmg"

hdiutil verify "$ROOT/../StrainScale-Agent-0.5.0.dmg"
```

本机实测（0.5.0 去模型版）：`.app` 约 26 MB，DMG 约 20 MB，`hdiutil verify` 返回 VALID。

### 1.6 验证打包结果

```bash
M=$(mktemp -d)
hdiutil attach -nobrowse -mountpoint "$M" StrainScale-Agent-0.5.0.dmg
codesign --verify --deep --strict "$M/StrainScale Agent.app" && echo SIG_OK

# 注意实际路径是双层 resources
PKG="$M/StrainScale Agent.app/Contents/Resources/resources/strainscale-pipeline/StrainScale-agent-package"
"$PKG/server-deploy/mbgc-macos-arm64" -v        # MBGC v2.1.6
"$PKG/server-deploy/hybrid-macos-arm64"         # 打印 Usage
hdiutil detach "$M"
```

### 1.7 打包后的目录长什么样

```
StrainScale Agent.app/Contents/
├── MacOS/
│   └── strainscale-agent              # Tauri 主程序
├── Resources/
│   ├── icon.icns
│   └── resources/                     # ← 多一层，来自 resources 相对路径
│       └── strainscale-pipeline/StrainScale-agent-package/…
└── Info.plist
```

0.5.0 已移除 llama.cpp sidecar 与 GGUF 模型，包内不再有 `llama-server`、`llama-sidecar/`、`models/`。

---

## 2. 重新编译原生流程二进制

### 2.1 mbgc（C++，需 libomp + zstd）

源码：`mbgc-source/mbgc-zstd-prestrain/`（`src/`、`vendor/zstd-1.5.6`、`build-mbgc-macos-arm64.sh`）。

```bash
cd mbgc-source/mbgc-zstd-prestrain

# 1) 先编 zstd 静态库（macOS 侧必须重编，不能复用 Linux 的 .a）
make -C vendor/zstd-1.5.6/lib libzstd.a CC=clang

# 2) 指定 OpenMP 前缀（Homebrew 或 conda 均可）
OMP_PREFIX="$(brew --prefix libomp)" bash build-mbgc-macos-arm64.sh

# 3) 产物：mbgc-macos-arm64 + 依赖 libomp.dylib
cp mbgc-macos-arm64 libomp.dylib \
   ../strainscale-desktop/src-tauri/resources/strainscale-pipeline/StrainScale-agent-package/server-deploy/
```

注意点：

- 脚本里 `OMP_PREFIX` 默认写死了一个 conda 路径，交接后必须用环境变量覆盖。
- 原 Linux 构建脚本的 `-static`、`-march=native`、x86 汇编库在 macOS 上全部不适用，已在本脚本里替换为 `clang` + `lib/arm/cpu_features.c` + `-Xpreprocessor -fopenmp`。
- 链接参数带 `-Wl,-rpath,@executable_path`，所以 `libomp.dylib` 必须和 `mbgc-macos-arm64` 放在同一目录；运行时还依赖 `/usr/lib`（libc++ 在 dyld 共享缓存里）。
- 验证：`./mbgc-macos-arm64 -v` 输出 `MBGC v2.1.6`；再跑一次 `p` 模式确认 `exitCode 0`。

### 2.2 hybrid（Rust）

源码已随包：`resources/strainscale-pipeline/StrainScale-agent-package/server-deploy/hybrid-rust/`。

```bash
cd .../server-deploy/hybrid-rust
cargo build --release
cp target/release/hybrid ../hybrid-macos-arm64
```

依赖只有 `rayon` + `mimalloc`，交叉编译成本低。Linux 静态版用 `build-hybrid.sh musl`（需要 `cargo-zigbuild` + `zig`），本机未验证。

### 2.3 mpileup / sambamba

- `sambamba-1.0.1`：随包是 arm64 Mach-O。上游是 D 语言项目，换平台需按 `sambamba-master` 重新编译（或直接用官方对应平台发布版）。
- `mpileup-standalone/mpileup`：随包是 arm64 Mach-O，来源是 htslib 的 `mpileup` 静态封装；换平台需重新编译。
- 二者都必须保留可执行权限；zip 解压会丢权限（见 3.4）。

---

## 3. 移植到其他平台

### 3.1 Windows：必须改的代码与配置

| 位置 | 现状 | 改法 |
|---|---|---|
| `run_pipeline_blocking` | `ProcessCommand::new("/bin/bash")` + `.sh` 脚本 | 三条路：① 要求用户装 WSL，调 `wsl bash ...`；② 把 `SS.sh/slice2.sh/get_site.pl` 用 Rust/Python 重写；③ 用 MSYS2 的 `bash.exe`。**推荐 ②**，长期最稳 |
| Stage B 暂存 | `#[cfg(unix)] std::fs::symlink` | Windows 没有普通用户可用的 symlink，改成复制文件或硬链接（`std::fs::hard_link`，同卷可用） |
| 原生二进制名 | `mbgc-macos-arm64`、`hybrid-macos-arm64` | 按 target 选择文件名，如 `mbgc-windows-x86_64.exe`；或统一改名 `mbgc`/`hybrid` 后由 `cfg!` 选目录 |
| OpenMP 依赖 | `libomp.dylib` 与二进制同目录 + rpath | Windows 用 `libomp.dll`（LLVM 发行版）放同目录，或静态链接 `/MT` |
| 路径分隔符 | Rust 侧 `PathBuf::join` 已跨平台 | shell 脚本内部大量 `/`，若保留 bash 方案需注意；Rust 侧无需改 |
| 签名 | `codesign` ad-hoc | Windows 用 Authenticode（可选），未签名会被 SmartScreen 提示 |

### 3.2 Linux：必须改的代码与配置

- `SS.sh`、`mbgc.zstd`、`hybrid` 本身就是 Linux x86_64，可以直接复用（`mbgc.zstd` 和 `hybrid` 已在包内留档）。
- `libomp` 换成 `libomp.so` 或用 GCC 的 `libgomp`。

### 3.3 跨平台资源命名建议

当前代码把 `mbgc-macos-arm64` / `hybrid-macos-arm64` 写死在两个函数里。建议第一步就抽成：

```rust
fn native_binary(pipeline_dir: &Path, base: &str) -> PathBuf {
    let name = if cfg!(target_os = "macos") && cfg!(target_arch = "aarch64") {
        format!("{base}-macos-arm64")
    } else if cfg!(target_os = "windows") {
        format!("{base}-windows-x86_64.exe")
    } else {
        format!("{base}-linux-x86_64")
    };
    pipeline_dir.join("server-deploy").join(name)
}
```

同时把 `libomp` 也做成按平台选择，避免每个 `run_*_blocking` 里重复判断。

### 3.4 打包类坑（真实踩过）

| 现象 | 原因 | 处理 |
|---|---|---|
| `Missing script: "tauri"` | `package.json` 的脚本叫 `dev`/`build`，不叫 `tauri` | 用 `npx tauri build` |
| 解压后二进制「秒退」/ `Permission denied` | zip 丢失可执行位 | `chmod +x` 所有原生二进制与 `.sh` |
| `Library not loaded: @rpath/libomp.dylib` | `libomp.dylib` 不在二进制同目录，或 rpath 被覆盖 | 保持同目录；必要时 `install_name_tool -add_rpath @executable_path` |
| `Library not loaded: @rpath/libc++.1.dylib` | 只设了 `@executable_path`，新系统 `/usr/lib` 无实体文件 | 追加 `-Wl,-rpath,/usr/lib`（mbgc 构建脚本已加） |
| 模型加载 503 | 模型在 Desktop/Documents，sidecar 受 TCC 限制读不到 | 已移除内嵌模型，此问题不再出现 |
| DMG 里签名校验失败 | 构建后没有重新 ad-hoc 签名 | 见 1.4 |
| `hdiutil attach -mountpoint` 警告已弃用 | macOS 新版本提示 | 仍可用；或改 `diskutil image attach --mountOptions nobrowse --mountPoint` |

---

## 4. 构建产物清单（0.5.0 实测）

| 产物 | 大小 | 校验 |
|---|---|---|
| `StrainScale Agent.app` | 约 26 MB | `codesign --verify --deep --strict` 通过（identifier `com.strainscale.agent`） |
| `StrainScale-Agent-0.5.0.dmg` | 约 20 MB | `hdiutil verify` VALID；挂载后 `mbgc -v` = v2.1.6，`hybrid` 打印 Usage |

---
