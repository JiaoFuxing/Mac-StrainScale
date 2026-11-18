use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command as ProcessCommand;
use std::sync::Mutex;
use tauri::{AppHandle, Emitter, Manager};

#[derive(Default)]
struct PipelineState(Mutex<Option<PipelineProgress>>);

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct AgentCallRequest {
    message: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct RouteResult {
    route: String,
    reason: String,
    files: Vec<String>,
    clarification: Option<String>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct PipelineRequest {
    sample_list: String,
    database_dir: String,
    working_dir: String,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct AlignStageRequest {
    fastq_list: String,
    reference: String,
    output_dir: String,
    threads: Option<u32>,
    jobs: Option<u32>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct PipelineModeRequest {
    mode: String,
    database_dir: String,
    list: String,
    #[serde(default)]
    output_dir: String,
    threads: Option<u32>,
    jobs: Option<u32>,
    cutoff: Option<u32>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct MbgcStageRequest {
    database_dir: String,
    samples_list: String,
    working_dir: String,
    pick_x: Option<u32>,
    max_jobs: Option<u32>,
    max_dist: Option<u32>,
    cutoff: Option<u32>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct HybridBatchRequest {
    working_dir: String,
    samples_list: String,
    cutoff: Option<u32>,
    max_dist: Option<u32>,
    match_threshold: Option<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct HybridBatchResult {
    success: bool,
    ran: Vec<String>,
    skipped: Vec<String>,
    failed: Vec<String>,
    outputs: Vec<String>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct HybridStageRequest {
    sample_output_dir: String,
    filter_counts: String,
    result_dir: String,
    cutoff: Option<u32>,
    max_dist: Option<u32>,
    match_threshold: Option<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PipelineResult {
    success: bool,
    exit_code: Option<i32>,
    working_dir: String,
    stdout: String,
    stderr: String,
    outputs: Vec<String>,
    errors: Vec<String>,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct PipelineProgress {
    mode: String,
    step: usize,
    total_steps: usize,
    percent: f64,
    stage: String,
    message: String,
    done: bool,
}

struct StageProgress<'a> {
    mode: &'a str,
    step: usize,
    total_steps: usize,
    band_start: f64,
    band_end: f64,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct SystemMetrics {
    cpu_percent: f64,
    memory_percent: f64,
    memory_used_gb: f64,
    memory_total_gb: f64,
    disk_percent: f64,
    disk_used_gb: f64,
    disk_total_gb: f64,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct AgentToolRequest {
    name: String,
    arguments: Value,
    #[serde(default)]
    message: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct AgentToolOutcome {
    mode: String,
    result: PipelineResult,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PipelinePlan {
    mode: String,
    database_dir: String,
    list: String,
    output_dir: String,
    cutoff: Option<u32>,
    reason: String,
}

fn emit_pipeline_progress(
    app: &AppHandle,
    mode: &str,
    step: usize,
    total_steps: usize,
    percent: f64,
    stage: &str,
    message: impl Into<String>,
    done: bool,
) {
    let progress = PipelineProgress {
        mode: mode.to_owned(),
        step,
        total_steps,
        percent: percent.clamp(0.0, 100.0),
        stage: stage.to_owned(),
        message: message.into(),
        done,
    };
    if let Some(state) = app.try_state::<PipelineState>() {
        if let Ok(mut current) = state.0.lock() {
            *current = Some(progress.clone());
        }
    }
    let _ = app.emit("pipeline-progress", progress);
}

#[tauri::command]
fn pipeline_status(app: AppHandle) -> Option<PipelineProgress> {
    app.state::<PipelineState>()
        .0
        .lock()
        .ok()
        .and_then(|current| current.clone())
}

fn expand_tilde(value: &str) -> String {
    if let Some(home) = std::env::var_os("HOME") {
        let home = PathBuf::from(home);
        if value == "~" {
            return home.display().to_string();
        }
        if let Some(rest) = value.strip_prefix("~/") {
            return home.join(rest).display().to_string();
        }
    }
    value.to_owned()
}

fn pipeline_log(message: &str) {
    use std::io::Write;
    if let Ok(mut file) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open("/tmp/strainscale-pipeline.log")
    {
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|duration| duration.as_secs())
            .unwrap_or(0);
        let _ = writeln!(file, "[{now}] {message}");
    }
}

fn stage_band(step: usize, total_steps: usize) -> (f64, f64) {
    let total = total_steps.max(1) as f64;
    let start = (step.saturating_sub(1) as f64 / total) * 100.0;
    let end = (step.min(total_steps) as f64 / total) * 100.0;
    (start, end)
}

fn finish_pipeline_progress(
    app: &AppHandle,
    mode: &str,
    total_steps: usize,
    percent: f64,
    stage: &str,
    result: PipelineResult,
) -> PipelineResult {
    let message = if result.success {
        "全部阶段执行完成".to_owned()
    } else {
        match result.errors.first() {
            Some(error) => format!("流程结束，存在错误：{error}"),
            None => "流程结束，存在错误".to_owned(),
        }
    };
    emit_pipeline_progress(
        app,
        mode,
        total_steps,
        total_steps,
        percent,
        stage,
        message,
        true,
    );
    result
}

fn vm_stat_pages(text: &str, key: &str) -> Option<f64> {
    text.lines()
        .find(|line| line.trim_start().starts_with(key))
        .and_then(|line| line.split(':').nth(1))
        .and_then(|value| value.trim().trim_end_matches('.').parse::<f64>().ok())
}

fn collect_system_metrics() -> Result<SystemMetrics, String> {
    let top_output = ProcessCommand::new("/usr/bin/top")
        .args(["-l", "1", "-n", "0"])
        .output()
        .map_err(|error| format!("failed to read cpu metrics: {error}"))?;
    let top_text = String::from_utf8_lossy(&top_output.stdout);
    let mut cpu_percent = 0.0;

    for line in top_text.lines() {
        if line.starts_with("CPU usage:") {
            for part in line.trim_start_matches("CPU usage:").split(',') {
                let part = part.trim();
                if let Some(value) = part.strip_suffix("% idle") {
                    if let Ok(idle) = value.trim().parse::<f64>() {
                        cpu_percent = (100.0 - idle).clamp(0.0, 100.0);
                    }
                }
            }
        }
    }

    let memory_total_bytes = ProcessCommand::new("/usr/sbin/sysctl")
        .args(["-n", "hw.memsize"])
        .output()
        .ok()
        .and_then(|output| {
            String::from_utf8_lossy(&output.stdout)
                .trim()
                .parse::<f64>()
                .ok()
        })
        .unwrap_or(0.0);
    let vm_output = ProcessCommand::new("/usr/bin/vm_stat")
        .output()
        .map_err(|error| format!("failed to read memory metrics: {error}"))?;
    let vm_text = String::from_utf8_lossy(&vm_output.stdout);
    let page_size = vm_text
        .lines()
        .next()
        .and_then(|line| line.split("page size of ").nth(1))
        .and_then(|rest| rest.split_whitespace().next())
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(4096.0);
    let active_pages = vm_stat_pages(&vm_text, "Pages active:").unwrap_or(0.0);
    let wired_pages = vm_stat_pages(&vm_text, "Pages wired down:").unwrap_or(0.0);
    let compressed_pages = vm_stat_pages(&vm_text, "Pages occupied by compressor:").unwrap_or(0.0);
    let memory_used_bytes = (active_pages + wired_pages + compressed_pages) * page_size;
    let memory_used_gb = memory_used_bytes / 1024.0 / 1024.0 / 1024.0;
    let memory_total_gb = memory_total_bytes / 1024.0 / 1024.0 / 1024.0;
    let memory_percent = if memory_total_bytes > 0.0 {
        (memory_used_bytes / memory_total_bytes * 100.0).clamp(0.0, 100.0)
    } else {
        0.0
    };

    let df_output = ProcessCommand::new("/bin/df")
        .args(["-k", "/"])
        .output()
        .map_err(|error| format!("failed to read disk metrics: {error}"))?;
    let df_text = String::from_utf8_lossy(&df_output.stdout);
    let last_line = df_text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .last()
        .unwrap_or("");
    let fields: Vec<&str> = last_line.split_whitespace().collect();
    let used_kb = fields
        .get(2)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0);
    let available_kb = fields
        .get(3)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0);
    let total_kb = used_kb + available_kb;
    let disk_percent = if total_kb > 0.0 {
        (used_kb / total_kb * 100.0).clamp(0.0, 100.0)
    } else {
        fields
            .get(4)
            .and_then(|value| value.trim_end_matches('%').parse::<f64>().ok())
            .unwrap_or(0.0)
    };

    Ok(SystemMetrics {
        cpu_percent,
        memory_percent,
        memory_used_gb: memory_used_gb.max(0.0),
        memory_total_gb: memory_total_gb.max(0.0),
        disk_percent,
        disk_used_gb: used_kb / 1024.0 / 1024.0,
        disk_total_gb: total_kb / 1024.0 / 1024.0,
    })
}

#[tauri::command]
async fn system_metrics() -> Result<SystemMetrics, String> {
    tauri::async_runtime::spawn_blocking(collect_system_metrics)
        .await
        .map_err(|error| format!("metrics task failed: {error}"))?
}

#[tauri::command]
fn open_external(url: String) -> Result<(), String> {
    let target = url.trim();
    if !(target.starts_with("https://") || target.starts_with("http://")) {
        return Err("only http(s) links are allowed".to_owned());
    }
    ProcessCommand::new("open")
        .arg(target)
        .status()
        .map_err(|error| format!("failed to open link: {error}"))?;
    Ok(())
}

fn contains_any(text: &str, patterns: &[&str]) -> bool {
    patterns.iter().any(|pattern| text.contains(pattern))
}

fn looks_like_file_reference(raw: &str) -> bool {
    let s = raw
        .trim()
        .trim_matches(|c: char| matches!(c, '"' | '\'' | '`' | '。' | '，' | '；'));
    if s.is_empty() {
        return false;
    }
    let lower = s.to_lowercase();

    if s.starts_with('/') || s.starts_with("~/") {
        return true;
    }

    if let Some(index) = s.rfind('/') {
        let last = &s[index + 1..];
        if last.contains('.') && !last.is_empty() {
            return true;
        }
    }

    let extensions = [
        ".fastq.gz",
        ".fq.gz",
        ".fasta.gz",
        ".fa.gz",
        ".vcf.gz",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
        ".fna",
        ".bam",
        ".sam",
        ".cram",
        ".vcf",
        ".bcf",
        ".csv",
        ".tsv",
        ".txt",
        ".json",
        ".xlsx",
        ".pdf",
        ".gguf",
        ".sh",
        ".py",
        ".list",
    ];
    extensions.iter().any(|ext| lower.ends_with(ext))
}

fn extract_files(message: &str) -> Vec<String> {
    let mut files = Vec::new();
    let delimiters = |c: char| {
        c.is_whitespace()
            || matches!(
                c,
                '，' | '。'
                    | '、'
                    | '；'
                    | '：'
                    | '“'
                    | '”'
                    | '"'
                    | '\''
                    | '('
                    | ')'
                    | '['
                    | ']'
                    | '{'
                    | '}'
                    | '<'
                    | '>'
                    | ','
                    | '?'
                    | '!'
                    | '？'
                    | '！'
            )
    };

    for token in message.split(delimiters) {
        let token = token.trim_matches('`');
        if looks_like_file_reference(token) && !files.iter().any(|item: &String| item == token) {
            files.push(token.to_owned());
        }
    }
    files
}

fn classify(message: &str) -> RouteResult {
    let text = message.trim();
    let lower = text.to_lowercase();
    let chat_patterns = [
        "是什么",
        "什么是",
        "是谁",
        "什么意思",
        "为什么",
        "解释",
        "介绍",
        "讲讲",
        "科普",
        "如何",
        "怎么",
        "怎样",
        "教程",
        "原理",
        "how to",
        "how do",
        "how does",
        "how can",
        "what is",
        "what's",
        "whats",
        "why",
        "explain",
        "define",
        "definition",
    ];
    let action_patterns = [
        "质控",
        "质量控制",
        "运行",
        "执行",
        "分析",
        "识别",
        "鉴定",
        "处理",
        "拼接",
        "比对",
        "丰度",
        "导出",
        "汇总",
        "生成报告",
        "检查",
        "qc",
        "quality control",
        "run ",
        "analyze",
        "identify",
        "align",
        "summarize",
        "export",
    ];

    let is_chat = contains_any(&lower, &chat_patterns);
    let is_action = contains_any(&lower, &action_patterns);
    let files = extract_files(text);

    if is_action && !files.is_empty() {
        return RouteResult {
            route: "skill".to_owned(),
            reason: "检测到明确分析动作和具体输入文件，进入受控 Skill 选择。".to_owned(),
            files,
            clarification: None,
        };
    }

    if is_chat {
        return RouteResult {
            route: "chat".to_owned(),
            reason: "知识解释或普通问答请求，不提供 tools，也不预填 tool_call。".to_owned(),
            files,
            clarification: None,
        };
    }

    if is_action {
        let clarification = "我可以执行这个分析，但需要你先选择或提供真实输入文件，例如 FASTQ、BAM、FASTA、VCF 或结果表。请避免只写“测序数据文件”这类占位描述。".to_owned();
        return RouteResult {
            route: "clarify".to_owned(),
            reason: "检测到分析意图，但没有发现具体路径或带扩展名的输入文件。".to_owned(),
            files,
            clarification: Some(clarification),
        };
    }

    if !files.is_empty() {
        return RouteResult {
            route: "skill".to_owned(),
            reason: "检测到具体输入文件，按分析任务进入 Skill 路由。".to_owned(),
            files,
            clarification: None,
        };
    }

    RouteResult {
        route: "chat".to_owned(),
        reason: "未检测到明确分析动作，按普通对话处理。".to_owned(),
        files,
        clarification: None,
    }
}

fn pipeline_resource_dir(app: &AppHandle) -> Result<PathBuf, String> {
    let bundled = app
        .path()
        .resource_dir()
        .map_err(|e| format!("failed to resolve resource directory: {e}"))?
        .join("resources/strainscale-pipeline/StrainScale-agent-package");
    if bundled.is_dir() {
        return Ok(bundled);
    }

    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("resources/strainscale-pipeline/StrainScale-agent-package");
    if source.is_dir() {
        Ok(source)
    } else {
        Err("bundled StrainScale pipeline directory does not exist".to_owned())
    }
}

fn validate_pipeline_request(
    request: &PipelineRequest,
    pipeline_dir: &Path,
) -> Result<Vec<(String, PathBuf)>, Vec<String>> {
    let mut errors = Vec::new();
    let sample_list = Path::new(&request.sample_list);
    let database_dir = Path::new(&request.database_dir);
    let working_dir = Path::new(&request.working_dir);

    if !sample_list.is_file() {
        errors.push(format!(
            "sample.list does not exist: {}",
            request.sample_list
        ));
    }
    if !database_dir.is_dir() {
        errors.push(format!(
            "database directory does not exist: {}",
            request.database_dir
        ));
    } else {
        if !database_dir.join("ref.fa").is_file() {
            errors.push("database directory is missing ref.fa".to_owned());
        }
        if !database_dir.join("add").is_file() {
            errors.push("database directory is missing add".to_owned());
        }
    }
    let script = pipeline_dir.join("SS.sh");
    if !script.is_file() {
        errors.push("pipeline entry SS.sh does not exist".to_owned());
    }

    if let Err(error) = std::fs::create_dir_all(working_dir) {
        errors.push(format!("failed to create working directory: {error}"));
    }

    let mut samples = Vec::new();
    if sample_list.is_file() {
        let content = std::fs::read_to_string(sample_list)
            .map_err(|e| vec![format!("failed to read sample.list: {e}")])?;
        for (line_number, raw_line) in content.lines().enumerate() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                errors.push(format!(
                    "line {} needs two tab-separated columns",
                    line_number + 1
                ));
                continue;
            }
            let sample_name = fields[0].trim();
            let bam_path = PathBuf::from(fields[1].trim().trim_end_matches('\r'));
            if sample_name.contains('/') || sample_name.contains("..") {
                errors.push(format!("invalid sample name: {sample_name}"));
            }
            if !bam_path.is_file() {
                errors.push(format!("BAM file does not exist: {}", bam_path.display()));
            }
            samples.push((sample_name.to_owned(), bam_path));
        }
    }

    if errors.is_empty() {
        Ok(samples)
    } else {
        Err(errors)
    }
}

fn collect_pipeline_outputs(working_dir: &Path, samples: &[(String, PathBuf)]) -> Vec<String> {
    let expected = [
        "filter-result-with-counts.txt",
        "output/2lie-mp-samtool",
        "output/really-ref.fa",
        "output/ref-coords.txt",
        "output/poc.txt",
        "output/poc2.txt",
    ];
    let mut outputs = Vec::new();
    for (sample_name, _) in samples {
        for relative in expected {
            let path = working_dir
                .join("YC-result")
                .join(sample_name)
                .join(relative);
            if path.is_file() {
                outputs.push(path.display().to_string());
            }
        }
    }
    outputs
}

fn run_pipeline_blocking(app: AppHandle, request: PipelineRequest) -> PipelineResult {
    let pipeline_dir = match pipeline_resource_dir(&app) {
        Ok(path) => path,
        Err(error) => {
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: request.working_dir,
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors: vec![error],
            };
        }
    };

    let samples = match validate_pipeline_request(&request, &pipeline_dir) {
        Ok(samples) => samples,
        Err(errors) => {
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: request.working_dir,
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors,
            };
        }
    };

    let script = pipeline_dir.join("SS.sh");
    let output = ProcessCommand::new("/bin/bash")
        .arg(&script)
        .arg(&request.sample_list)
        .arg(&request.database_dir)
        .current_dir(&request.working_dir)
        .output();

    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout).to_string();
            let stderr = String::from_utf8_lossy(&output.stderr).to_string();
            let outputs = collect_pipeline_outputs(Path::new(&request.working_dir), &samples);
            let mut errors = Vec::new();
            let exit_code = output.status.code();
            if !output.status.success() {
                errors.push(format!("SS.sh exited with status: {:?}", output.status));
            }
            if outputs.iter().any(|path| path.ends_with("really-ref.fa")) == false {
                errors.push("really-ref.fa was not produced".to_owned());
            }

            PipelineResult {
                success: errors.is_empty(),
                exit_code,
                working_dir: request.working_dir,
                stdout,
                stderr,
                outputs,
                errors,
            }
        }
        Err(error) => PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.working_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors: vec![error.to_string()],
        },
    }
}

fn is_executable_file(path: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::metadata(path)
            .map(|metadata| metadata.is_file() && metadata.permissions().mode() & 0o111 != 0)
            .unwrap_or(false)
    }
    #[cfg(not(unix))]
    {
        path.is_file()
    }
}

fn align_rammap_binary(pipeline_dir: &Path) -> PathBuf {
    let bundled = pipeline_dir.join("server-deploy/rammap-macos-arm64");
    if bundled.is_file() {
        bundled
    } else {
        pipeline_dir.join("bin/rammap")
    }
}

fn validate_align_request(
    request: &AlignStageRequest,
    pipeline_dir: &Path,
) -> Result<(), Vec<String>> {
    let mut errors = Vec::new();
    let fastq_list = Path::new(&request.fastq_list);
    let reference = Path::new(&request.reference);
    let output_dir = Path::new(&request.output_dir);

    if !fastq_list.is_file() {
        errors.push(format!("FASTQ list does not exist: {}", request.fastq_list));
    }
    if !reference.is_file() {
        errors.push(format!(
            "reference FASTA does not exist: {}",
            request.reference
        ));
    }
    if request.threads == Some(0) {
        errors.push("threads must be greater than 0".to_owned());
    }
    if request.jobs == Some(0) {
        errors.push("parallel sample count must be greater than 0".to_owned());
    }
    if let Err(error) = std::fs::create_dir_all(output_dir) {
        errors.push(format!("failed to create output directory: {error}"));
    }

    let script = pipeline_dir.join("align_sort.sh");
    if !script.is_file() {
        errors.push("align_sort.sh does not exist".to_owned());
    }
    let rammap = align_rammap_binary(pipeline_dir);
    if !is_executable_file(&rammap) {
        errors.push(format!(
            "rammap is missing or not executable: {}",
            rammap.display()
        ));
    }
    let sambamba = pipeline_dir.join("sambamba-master/bin/sambamba-1.0.1");
    if !is_executable_file(&sambamba) {
        errors.push(format!(
            "sambamba is missing or not executable: {}",
            sambamba.display()
        ));
    }

    if fastq_list.is_file() {
        match std::fs::read_to_string(fastq_list) {
            Ok(content) => {
                for (line_number, raw_line) in content.lines().enumerate() {
                    let line = raw_line.trim_end_matches('\r');
                    if line.trim().is_empty() || line.trim_start().starts_with('#') {
                        continue;
                    }
                    let fields: Vec<&str> = line.split('\t').collect();
                    if fields.len() < 2 {
                        errors.push(format!(
                            "line {} needs sample<TAB>fastq1[<TAB>fastq2] columns",
                            line_number + 1
                        ));
                        continue;
                    }
                    let sample_name = fields[0].trim();
                    if sample_name.is_empty()
                        || sample_name.contains('/')
                        || sample_name.contains("..")
                    {
                        errors.push(format!(
                            "line {} has an invalid sample name: {sample_name}",
                            line_number + 1
                        ));
                    }
                    let fastq1 = Path::new(fields[1].trim());
                    if !fastq1.is_file() {
                        errors.push(format!(
                            "line {} FASTQ does not exist: {}",
                            line_number + 1,
                            fastq1.display()
                        ));
                    }
                    if let Some(fastq2) = fields
                        .get(2)
                        .map(|value| value.trim())
                        .filter(|value| !value.is_empty())
                    {
                        let fastq2 = Path::new(fastq2);
                        if !fastq2.is_file() {
                            errors.push(format!(
                                "line {} FASTQ does not exist: {}",
                                line_number + 1,
                                fastq2.display()
                            ));
                        }
                    }
                }
            }
            Err(error) => errors.push(format!("failed to read FASTQ list: {error}")),
        }
    }

    if errors.is_empty() {
        Ok(())
    } else {
        Err(errors)
    }
}

fn collect_align_outputs(output_dir: &Path) -> Vec<String> {
    let mut outputs = Vec::new();
    let sample_list = output_dir.join("sample.list");
    if sample_list.is_file() {
        outputs.push(sample_list.display().to_string());
    }
    let mut bam_paths = Vec::new();
    if let Ok(entries) = std::fs::read_dir(output_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            let extension = path
                .extension()
                .and_then(|value| value.to_str())
                .unwrap_or_default();
            if extension == "bam" || extension == "bai" {
                bam_paths.push(path);
            }
        }
    }
    bam_paths.sort();
    outputs.extend(bam_paths.iter().map(|path| path.display().to_string()));
    outputs
}

fn run_align_blocking(app: AppHandle, request: AlignStageRequest) -> PipelineResult {
    let pipeline_dir = match pipeline_resource_dir(&app) {
        Ok(path) => path,
        Err(error) => {
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: request.output_dir,
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors: vec![error],
            };
        }
    };

    if let Err(errors) = validate_align_request(&request, &pipeline_dir) {
        return PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.output_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors,
        };
    }

    let threads = request.threads.unwrap_or(8).max(1);
    let output = ProcessCommand::new("/bin/bash")
        .arg(pipeline_dir.join("align_sort.sh"))
        .arg("-r")
        .arg(&request.reference)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-x")
        .arg("auto")
        .arg("-j")
        .arg(
            request
                .jobs
                .map(|jobs| jobs.to_string())
                .unwrap_or_else(|| "auto".to_owned()),
        )
        .arg("-R")
        .arg(align_rammap_binary(&pipeline_dir))
        .arg("-S")
        .arg(pipeline_dir.join("sambamba-master/bin/sambamba-1.0.1"))
        .arg(&request.fastq_list)
        .arg(&request.output_dir)
        .output();

    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout).to_string();
            let stderr = String::from_utf8_lossy(&output.stderr).to_string();
            let output_dir = Path::new(&request.output_dir);
            let outputs = collect_align_outputs(output_dir);
            let mut errors = Vec::new();
            let exit_code = output.status.code();
            if !output.status.success() {
                errors.push(format!(
                    "align_sort.sh exited with status: {:?}",
                    output.status
                ));
            }
            let sample_list = output_dir.join("sample.list");
            if !sample_list.is_file()
                || std::fs::metadata(&sample_list)
                    .map(|m| m.len())
                    .unwrap_or(0)
                    == 0
            {
                errors.push("sample.list was not produced or is empty".to_owned());
            }
            if !outputs.iter().any(|path| path.ends_with(".bam")) {
                errors.push("no sorted BAM was produced".to_owned());
            }

            PipelineResult {
                success: errors.is_empty(),
                exit_code,
                working_dir: request.output_dir,
                stdout,
                stderr,
                outputs,
                errors,
            }
        }
        Err(error) => PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.output_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors: vec![error.to_string()],
        },
    }
}

fn run_hybrid_blocking(app: AppHandle, request: HybridStageRequest) -> PipelineResult {
    let pipeline_dir = match pipeline_resource_dir(&app) {
        Ok(path) => path,
        Err(error) => {
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: request.result_dir,
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors: vec![error],
            };
        }
    };
    let binary = pipeline_dir.join("server-deploy/hybrid-macos-arm64");
    let sample_output_dir = Path::new(&request.sample_output_dir);
    let filter_counts = Path::new(&request.filter_counts);
    let result_dir = Path::new(&request.result_dir);
    let mut errors = Vec::new();

    if !binary.is_file() {
        errors.push(format!(
            "native hybrid binary does not exist: {}",
            binary.display()
        ));
    }
    if !sample_output_dir.is_dir() {
        errors.push(format!(
            "sample output directory does not exist: {}（Stage A 未产出或样本名写错）",
            sample_output_dir.display()
        ));
    }
    if !filter_counts.is_file() {
        errors.push(format!(
            "filter counts file does not exist: {}（Stage A 的 get_site.pl 未写出）",
            filter_counts.display()
        ));
    }
    if errors.is_empty() {
        if let Err(error) = std::fs::create_dir_all(result_dir) {
            errors.push(format!("failed to create result directory: {error}"));
        }
    }
    if !errors.is_empty() {
        return PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.result_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors,
        };
    }

    let cutoff = request.cutoff.unwrap_or(1000);
    let max_dist = request.max_dist.unwrap_or(10000);
    let match_threshold = request.match_threshold.unwrap_or(0.9);
    let output = ProcessCommand::new(&binary)
        .arg(sample_output_dir)
        .arg(cutoff.to_string())
        .arg(max_dist.to_string())
        .arg(filter_counts)
        .arg("-t")
        .arg(match_threshold.to_string())
        .arg("-o")
        .arg(result_dir)
        .output();

    match output {
        Ok(output) => {
            let stdout = String::from_utf8_lossy(&output.stdout).to_string();
            let stderr = String::from_utf8_lossy(&output.stderr).to_string();
            let expected = [
                "final_representatives.txt",
                "final_match_results.tsv",
                "final_specific_sites.tsv",
                "final_clusters.tsv",
            ];
            let outputs = expected
                .iter()
                .map(|name| result_dir.join(name))
                .filter(|path| path.is_file())
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>();
            let mut errors = Vec::new();
            if !output.status.success() {
                errors.push(format!("hybrid exited with status: {:?}", output.status));
            }
            if outputs.is_empty() {
                errors.push("hybrid did not produce final result files".to_owned());
            }

            PipelineResult {
                success: errors.is_empty(),
                exit_code: output.status.code(),
                working_dir: request.result_dir,
                stdout,
                stderr,
                outputs,
                errors,
            }
        }
        Err(error) => PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.result_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors: vec![error.to_string()],
        },
    }
}

fn collect_mbgc_databases(database_dir: &Path) -> Result<Vec<String>, Vec<String>> {
    let mut errors = Vec::new();
    let db_root = if database_dir.join("DB").is_dir() {
        database_dir.join("DB")
    } else {
        database_dir.to_path_buf()
    };

    let mut databases = Vec::new();
    let list_file = database_dir.join("DB.list");
    if list_file.is_file() {
        match std::fs::read_to_string(&list_file) {
            Ok(content) => {
                for raw_line in content.lines() {
                    let name = raw_line.trim().trim_end_matches('\r');
                    if name.is_empty() || name.starts_with('#') {
                        continue;
                    }
                    databases.push(name.to_owned());
                }
            }
            Err(error) => errors.push(format!("failed to read DB.list: {error}")),
        }
    } else if let Ok(entries) = std::fs::read_dir(&db_root) {
        for entry in entries.flatten() {
            let file_name = entry.file_name().to_string_lossy().to_string();
            if let Some(name) = file_name.strip_suffix(".zst") {
                databases.push(name.to_owned());
            }
        }
        databases.sort();
    }

    if databases.is_empty() {
        errors.push(
            "no database found: DB.list is empty and no *.zst under the database directory"
                .to_owned(),
        );
    }

    for db in &databases {
        for extension in ["zst", "clu", "mbgc"] {
            let path = db_root.join(format!("{db}.{extension}"));
            if !path.is_file() {
                errors.push(format!("database file is missing: {}", path.display()));
            }
        }
    }

    if errors.is_empty() {
        Ok(databases)
    } else {
        Err(errors)
    }
}

fn parse_cluster_header(content: &str) -> Result<Option<Vec<u32>>, String> {
    let first_line = content
        .lines()
        .next()
        .map(|line| line.trim_end_matches('\r'))
        .filter(|line| !line.trim().is_empty())
        .ok_or_else(|| ".clu file is empty".to_owned())?;
    let fields: Vec<&str> = first_line.split('\t').collect();
    if fields[0].trim() != "Strain-name" {
        if first_line.split(',').count() >= 2 {
            return Ok(None);
        }
        return Err(
            "invalid .clu header: expected Strain-name<TAB><cutoff...> or legacy name,cluster_id"
                .to_owned(),
        );
    }
    if fields.len() < 2 {
        return Err("invalid .clu header: no cutoff columns".to_owned());
    }

    let mut cutoffs = Vec::new();
    for field in &fields[1..] {
        let value = field.trim();
        let cutoff = value
            .parse::<u32>()
            .map_err(|_| format!("invalid .clu cutoff column: {value}"))?;
        if cutoffs.contains(&cutoff) {
            return Err(format!("duplicate .clu cutoff column: {cutoff}"));
        }
        cutoffs.push(cutoff);
    }
    Ok(Some(cutoffs))
}

fn common_cluster_cutoffs(groups: &[Vec<u32>]) -> Vec<u32> {
    let Some(first) = groups.first() else {
        return Vec::new();
    };
    first
        .iter()
        .copied()
        .filter(|cutoff| groups[1..].iter().all(|group| group.contains(cutoff)))
        .collect()
}

fn collect_cluster_cutoffs(database_dir: &Path) -> Result<Vec<u32>, String> {
    let databases = collect_mbgc_databases(database_dir).map_err(|errors| errors.join("\n"))?;
    let db_root = if database_dir.join("DB").is_dir() {
        database_dir.join("DB")
    } else {
        database_dir.to_path_buf()
    };
    let mut groups = Vec::new();
    for db in databases {
        let path = db_root.join(format!("{db}.clu"));
        let content = std::fs::read_to_string(&path)
            .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
        let cutoffs = parse_cluster_header(&content)
            .map_err(|error| format!("{}: {error}", path.display()))?
            .unwrap_or_else(|| vec![1000]);
        groups.push(cutoffs);
    }
    let common = common_cluster_cutoffs(&groups);
    if common.is_empty() {
        return Err("databases have no common cluster cutoff".to_owned());
    }
    Ok(common)
}

fn resolve_cluster_cutoff(database_dir: &Path, requested: Option<u32>) -> Result<u32, String> {
    let supported = collect_cluster_cutoffs(database_dir)?;
    let cutoff = requested.unwrap_or_else(|| {
        if supported.contains(&1000) {
            1000
        } else {
            supported[0]
        }
    });
    if supported.contains(&cutoff) {
        Ok(cutoff)
    } else {
        Err(format!(
            "cutoff {cutoff} is not supported by this database; available: {}",
            supported
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        ))
    }
}

fn selected_cluster_content(content: &str, cutoff: u32) -> Result<String, String> {
    let cutoffs = parse_cluster_header(content)?
        .ok_or_else(|| "legacy .clu does not need conversion".to_owned())?;
    let column = cutoffs
        .iter()
        .position(|value| *value == cutoff)
        .map(|index| index + 1)
        .ok_or_else(|| format!("cutoff {cutoff} is missing from .clu header"))?;
    let expected_fields = cutoffs.len() + 1;
    let mut output = String::new();
    for (index, raw_line) in content.lines().enumerate().skip(1) {
        let line = raw_line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != expected_fields {
            return Err(format!(
                "invalid .clu row {}: expected {expected_fields} columns, got {}",
                index + 1,
                fields.len()
            ));
        }
        let name = fields[0].trim();
        let cluster = fields[column].trim();
        if name.is_empty() || cluster.is_empty() {
            return Err(format!(
                "invalid .clu row {}: empty strain or cluster",
                index + 1
            ));
        }
        output.push_str(name);
        output.push(',');
        output.push_str(cluster);
        output.push('\n');
    }
    Ok(output)
}

fn collect_mbgc_outputs(
    working_dir: &Path,
    samples: &[String],
    databases: &[String],
) -> Vec<String> {
    let mut outputs = Vec::new();
    for sample in samples {
        for db in databases {
            let rep_table = working_dir
                .join("YC-result")
                .join(sample)
                .join("output")
                .join(format!("rep_{db}.tsv"));
            if rep_table.is_file() {
                outputs.push(rep_table.display().to_string());
            }
            let rep_dir = working_dir
                .join("YC-result")
                .join(sample)
                .join("output")
                .join(format!("fa_{db}"));
            if rep_dir.is_dir() {
                outputs.push(rep_dir.display().to_string());
            }
        }
    }
    outputs
}

fn run_mbgc_blocking(app: AppHandle, request: MbgcStageRequest) -> PipelineResult {
    let pipeline_dir = match pipeline_resource_dir(&app) {
        Ok(path) => path,
        Err(error) => {
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: request.working_dir,
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors: vec![error],
            };
        }
    };
    let binary = pipeline_dir.join("server-deploy/mbgc-macos-arm64");
    let database_dir = Path::new(&request.database_dir);
    let samples_list = Path::new(&request.samples_list);
    let working_dir = Path::new(&request.working_dir);
    let mut errors = Vec::new();

    if !binary.is_file() {
        errors.push("native mbgc binary does not exist".to_owned());
    }
    let omp_library = binary.parent().map(|dir| dir.join("libomp.dylib"));
    if let Some(path) = omp_library {
        if !path.is_file() {
            errors.push("libomp.dylib is missing next to the native mbgc binary".to_owned());
        }
    }
    if !database_dir.is_dir() {
        errors.push(format!(
            "database directory does not exist: {}",
            request.database_dir
        ));
    }
    if !database_dir.join("G0.fa").is_file() {
        errors.push("database directory is missing G0.fa".to_owned());
    }
    if !samples_list.is_file() {
        errors.push(format!(
            "samples list does not exist: {}",
            request.samples_list
        ));
    }
    if let Err(error) = std::fs::create_dir_all(working_dir) {
        errors.push(format!("failed to create working directory: {error}"));
    }

    let databases = if database_dir.is_dir() {
        match collect_mbgc_databases(database_dir) {
            Ok(list) => list,
            Err(mut list) => {
                errors.append(&mut list);
                Vec::new()
            }
        }
    } else {
        Vec::new()
    };
    let cutoff = if !databases.is_empty() {
        match resolve_cluster_cutoff(database_dir, request.cutoff) {
            Ok(value) => value,
            Err(error) => {
                errors.push(error);
                request.cutoff.unwrap_or(1000)
            }
        }
    } else {
        request.cutoff.unwrap_or(1000)
    };

    let mut samples = Vec::new();
    if samples_list.is_file() {
        if let Ok(content) = std::fs::read_to_string(samples_list) {
            for raw_line in content.lines() {
                let line = raw_line.trim().trim_end_matches('\r');
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let name = line.split_whitespace().next().unwrap_or("").to_owned();
                if name.is_empty() {
                    continue;
                }
                if name.contains('/') || name.contains("..") {
                    errors.push(format!("invalid sample name: {name}"));
                    continue;
                }
                let sample_output = working_dir.join("YC-result").join(&name).join("output");
                if !sample_output.join("really-ref.fa").is_file() {
                    errors.push(format!(
                        "missing really-ref.fa for sample {name}; run Stage A first"
                    ));
                }
                if !sample_output.join("ref-coords.txt").is_file() {
                    errors.push(format!(
                        "missing ref-coords.txt for sample {name}; run Stage A first"
                    ));
                }
                samples.push(name);
            }
        } else {
            errors.push("failed to read samples list".to_owned());
        }
    }
    if samples.is_empty() && errors.is_empty() {
        errors.push("samples list contains no sample".to_owned());
    }

    if !errors.is_empty() {
        return PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.working_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors,
        };
    }

    let pick_x = request.pick_x.unwrap_or(30).max(1);
    let max_jobs = request.max_jobs.unwrap_or(4).clamp(1, 12);
    let max_dist = request.max_dist.unwrap_or(10000);
    let db_root = if database_dir.join("DB").is_dir() {
        database_dir.join("DB")
    } else {
        database_dir.to_path_buf()
    };

    let mut stage_root = working_dir.join(".mbgc-run");
    if database_dir.starts_with(&stage_root) {
        stage_root = std::env::temp_dir().join(format!("mbgc-run-{}", std::process::id()));
    }
    let _ = std::fs::remove_dir_all(&stage_root);
    if let Err(error) = std::fs::create_dir_all(&stage_root) {
        return PipelineResult {
            success: false,
            exit_code: None,
            working_dir: request.working_dir,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors: vec![format!("failed to create mbgc staging directory: {error}")],
        };
    }

    let mut per_db_dirs: Vec<(String, PathBuf)> = Vec::new();
    for db in &databases {
        let stage_dir = stage_root.join(db);
        if let Err(error) = std::fs::create_dir_all(stage_dir.join("DB")) {
            errors.push(format!(
                "failed to prepare staging directory for {db}: {error}"
            ));
            continue;
        }
        let clu_path = db_root.join(format!("{db}.clu"));
        let converted_clu = match std::fs::read_to_string(&clu_path) {
            Ok(content) => match parse_cluster_header(&content) {
                Ok(Some(_)) => match selected_cluster_content(&content, cutoff) {
                    Ok(converted) => Some(converted),
                    Err(error) => {
                        errors.push(format!("{}: {error}", clu_path.display()));
                        continue;
                    }
                },
                Ok(None) => None,
                Err(error) => {
                    errors.push(format!("{}: {error}", clu_path.display()));
                    continue;
                }
            },
            Err(error) => {
                errors.push(format!("failed to read {}: {error}", clu_path.display()));
                continue;
            }
        };
        #[cfg(unix)]
        {
            use std::os::unix::fs::symlink;
            let _ = symlink(database_dir.join("G0.fa"), stage_dir.join("G0.fa"));
            let _ = symlink(
                db_root.join(format!("{db}.zst")),
                stage_dir.join("DB").join(format!("{db}.zst")),
            );
            if let Some(content) = converted_clu {
                if let Err(error) =
                    std::fs::write(stage_dir.join("DB").join(format!("{db}.clu")), content)
                {
                    errors.push(format!("failed to write staged .clu for {db}: {error}"));
                    continue;
                }
            } else {
                let _ = symlink(&clu_path, stage_dir.join("DB").join(format!("{db}.clu")));
            }
            let _ = symlink(
                db_root.join(format!("{db}.mbgc")),
                stage_dir.join("DB").join(format!("{db}.mbgc")),
            );
        }
        let _ = std::fs::write(stage_dir.join("DB.list"), format!("{db}\n"));
        per_db_dirs.push((db.clone(), stage_dir));
    }

    let mut logs: Vec<String> = Vec::new();
    let mut failed: Vec<String> = Vec::new();

    for chunk in per_db_dirs.chunks(max_jobs as usize) {
        let mut handles = Vec::new();
        for (db, stage_dir) in chunk {
            let db_name = db.clone();
            let stage = stage_dir.clone();
            let samples_file = samples_list.to_path_buf();
            let run_dir = working_dir.to_path_buf();
            let binary_path = binary.clone();
            handles.push(std::thread::spawn(move || {
                let result = ProcessCommand::new(&binary_path)
                    .arg("p")
                    .arg("-f")
                    .arg("-d")
                    .arg(&stage)
                    .arg(&samples_file)
                    .arg(pick_x.to_string())
                    .arg(max_dist.to_string())
                    .current_dir(&run_dir)
                    .output();
                let log_path = stage.join("run.log");
                match result {
                    Ok(output) => {
                        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
                        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
                        let _ = std::fs::write(
                            &log_path,
                            format!(
                                "exit: {:?}\n\n== stdout ==\n{}\n== stderr ==\n{}",
                                output.status, stdout, stderr
                            ),
                        );
                        (
                            db_name,
                            output.status.success(),
                            output.status.code(),
                            stdout,
                            stderr,
                        )
                    }
                    Err(error) => {
                        let _ = std::fs::write(&log_path, format!("failed to start mbgc: {error}"));
                        (db_name, false, None, String::new(), error.to_string())
                    }
                }
            }));
        }
        for handle in handles {
            if let Ok((db, ok, code, stdout, stderr)) = handle.join() {
                logs.push(format!("[{db}] exit={code:?} ok={ok}"));
                if !ok {
                    failed.push(db.clone());
                }
                if !stdout.trim().is_empty() {
                    logs.push(format!("[{db}] stdout:\n{}", stdout.trim()));
                }
                if !stderr.trim().is_empty() {
                    logs.push(format!("[{db}] stderr:\n{}", stderr.trim()));
                }
            }
        }
    }

    let outputs = collect_mbgc_outputs(working_dir, &samples, &databases);
    let mut final_errors = errors;
    if !failed.is_empty() {
        final_errors.push(format!(
            "mbgc failed for database(s): {}",
            failed.join(", ")
        ));
    }
    if outputs.is_empty() {
        final_errors.push("mbgc did not produce rep_<db>.tsv or fa_<db> outputs".to_owned());
    }

    PipelineResult {
        success: final_errors.is_empty(),
        exit_code: if failed.is_empty() { Some(0) } else { Some(1) },
        working_dir: request.working_dir,
        stdout: logs.join("\n"),
        stderr: String::new(),
        outputs,
        errors: final_errors,
    }
}

fn is_pipeline_mode(mode: &str) -> bool {
    matches!(
        mode,
        "reads2Result" | "bam2Result" | "2lie-mp-samtool2Result"
    )
}

fn resolve_list_input(base: &Path, value: &str) -> PathBuf {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        path
    } else {
        base.join(path)
    }
}

fn parse_pipeline_mode_list(
    list: &Path,
    mode: &str,
) -> Result<Vec<(String, Vec<PathBuf>)>, Vec<String>> {
    if !list.is_file() {
        return Err(vec![format!("list does not exist: {}", list.display())]);
    }
    let base = list
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let content = match std::fs::read_to_string(list) {
        Ok(content) => content,
        Err(error) => return Err(vec![format!("failed to read list: {error}")]),
    };

    let expected = if mode == "reads2Result" {
        "样本名<TAB>fastq1[<TAB>fastq2]"
    } else {
        "样本名<TAB>文件路径（两列）"
    };
    let mut errors = Vec::new();
    let mut entries: Vec<(String, Vec<PathBuf>)> = Vec::new();
    let mut seen_names: Vec<String> = Vec::new();

    for (line_number, raw_line) in content.lines().enumerate() {
        let line = raw_line.trim_end_matches('\r');
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let mut fields: Vec<&str> = if line.contains('\t') {
            line.split('\t').collect()
        } else {
            line.split_whitespace().collect()
        };
        while fields
            .last()
            .map(|value| value.trim().is_empty())
            .unwrap_or(false)
        {
            fields.pop();
        }

        let valid_columns = if mode == "reads2Result" {
            (2..=3).contains(&fields.len())
        } else {
            fields.len() == 2
        };
        if !valid_columns {
            errors.push(format!("line {} needs {expected}", line_number + 1));
            continue;
        }

        let name = fields[0].trim();
        if name.is_empty() || name.contains('/') || name.contains("..") {
            errors.push(format!(
                "line {} has an invalid sample name: {name}",
                line_number + 1
            ));
            continue;
        }
        if seen_names.iter().any(|existing| existing == name) {
            errors.push(format!(
                "line {} duplicates sample name: {name}",
                line_number + 1
            ));
            continue;
        }
        seen_names.push(name.to_owned());

        let paths: Vec<PathBuf> = fields[1..]
            .iter()
            .map(|value| resolve_list_input(base, value.trim()))
            .collect();
        for path in &paths {
            if !path.is_file() {
                errors.push(format!(
                    "line {} file does not exist: {}",
                    line_number + 1,
                    path.display()
                ));
            }
        }
        entries.push((name.to_owned(), paths));
    }

    if entries.is_empty() && errors.is_empty() {
        errors.push("list contains no sample".to_owned());
    }
    if errors.is_empty() {
        Ok(entries)
    } else {
        Err(errors)
    }
}

fn write_normalized_list(path: &Path, entries: &[(String, Vec<PathBuf>)]) -> Result<(), String> {
    let mut content = String::new();
    for (name, paths) in entries {
        content.push_str(name);
        for path in paths {
            content.push('\t');
            content.push_str(&path.display().to_string());
        }
        content.push('\n');
    }
    std::fs::write(path, content)
        .map_err(|error| format!("failed to write {}: {error}", path.display()))
}

fn merge_stage_result(combined: &mut PipelineResult, label: &str, stage: PipelineResult) -> bool {
    combined
        .stdout
        .push_str(&format!("\n===== {label} =====\n"));
    if !stage.stdout.trim().is_empty() {
        combined.stdout.push_str(stage.stdout.trim_end());
        combined.stdout.push('\n');
    }
    combined
        .stderr
        .push_str(&format!("\n===== {label} =====\n"));
    if !stage.stderr.trim().is_empty() {
        combined.stderr.push_str(stage.stderr.trim_end());
        combined.stderr.push('\n');
    }
    combined.outputs.extend(stage.outputs);
    for error in stage.errors {
        combined.errors.push(format!("[{label}] {error}"));
    }
    if !stage.success {
        combined.success = false;
        combined.exit_code = stage.exit_code.or(Some(1));
        return false;
    }
    true
}

fn estimate_depth_from_pileup(path: &Path) -> Option<u32> {
    use std::io::{BufRead, BufReader};
    let file = std::fs::File::open(path).ok()?;
    let reader = BufReader::new(file);
    let mut counts: HashMap<u32, u32> = HashMap::new();
    for line in reader.lines().take(10000).flatten() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        if let Ok(depth) = fields[4].parse::<u32>() {
            if depth > 0 {
                *counts.entry(depth).or_insert(0) += 1;
            }
        }
    }
    counts
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(depth, _)| depth)
}

fn import_pileup_samples(
    app: &AppHandle,
    pipeline_dir: &Path,
    entries: &[(String, Vec<PathBuf>)],
    working_dir: &Path,
    progress: &StageProgress<'_>,
) -> PipelineResult {
    let script = pipeline_dir.join("get_site.pl");
    let perl = if Path::new("/usr/bin/perl").is_file() {
        "/usr/bin/perl"
    } else {
        "perl"
    };
    let mut stdout = String::new();
    let mut stderr = String::new();
    let mut outputs = Vec::new();
    let mut errors = Vec::new();
    let total = entries.len().max(1);
    let stage = "Stage A' · get_site 一致序列";

    for (index, (name, paths)) in entries.iter().enumerate() {
        let start_percent = progress.band_start
            + (progress.band_end - progress.band_start) * (index as f64 / total as f64);
        emit_pipeline_progress(
            app,
            progress.mode,
            progress.step,
            progress.total_steps,
            start_percent,
            stage,
            format!("样本 {}/{}：{name} 处理中…", index + 1, entries.len()),
            false,
        );
        let input = &paths[0];
        let sample_dir = working_dir.join("YC-result").join(name);
        let output_dir = sample_dir.join("output");
        let sample_result: Result<Vec<String>, String> = (|| {
            std::fs::create_dir_all(&output_dir)
                .map_err(|error| format!("{name}: failed to create output directory: {error}"))?;
            let depth = estimate_depth_from_pileup(input).unwrap_or(0);
            std::fs::write(sample_dir.join("estimated_depth.txt"), depth.to_string())
                .map_err(|error| format!("{name}: failed to write estimated_depth.txt: {error}"))?;
            let error_base = if depth < 100 { "0" } else { "1" };
            let result = std::process::Command::new(perl)
                .arg(&script)
                .arg(input)
                .arg(error_base)
                .arg("4")
                .current_dir(&sample_dir)
                .output()
                .map_err(|error| format!("{name}: failed to run get_site.pl: {error}"))?;
            stdout.push_str(&format!("[{name}] get_site.pl\n"));
            stdout.push_str(&String::from_utf8_lossy(&result.stdout));
            stderr.push_str(&format!("[{name}] get_site.pl\n"));
            stderr.push_str(&String::from_utf8_lossy(&result.stderr));
            if !result.status.success() {
                return Err(format!(
                    "{name}: get_site.pl exited with status {:?}",
                    result.status
                ));
            }
            let really_ref = output_dir.join("really-ref.fa");
            let ref_coords = output_dir.join("ref-coords.txt");
            let filter_counts = sample_dir.join("filter-result-with-counts.txt");
            if !really_ref.is_file() {
                return Err(format!("{name}: really-ref.fa was not produced"));
            }
            if !ref_coords.is_file() {
                return Err(format!("{name}: ref-coords.txt was not produced"));
            }
            if !filter_counts.is_file() {
                return Err(format!(
                    "{name}: filter-result-with-counts.txt was not produced"
                ));
            }
            Ok(vec![
                really_ref.display().to_string(),
                ref_coords.display().to_string(),
                filter_counts.display().to_string(),
            ])
        })();

        let done_percent = progress.band_start
            + (progress.band_end - progress.band_start) * ((index + 1) as f64 / total as f64);
        match sample_result {
            Ok(sample_outputs) => {
                outputs.extend(sample_outputs);
                emit_pipeline_progress(
                    app,
                    progress.mode,
                    progress.step,
                    progress.total_steps,
                    done_percent,
                    stage,
                    format!("样本 {}/{}：{name} 完成", index + 1, entries.len()),
                    false,
                );
            }
            Err(error) => {
                emit_pipeline_progress(
                    app,
                    progress.mode,
                    progress.step,
                    progress.total_steps,
                    done_percent,
                    stage,
                    format!("样本 {}/{}：{name} 失败：{error}", index + 1, entries.len()),
                    false,
                );
                errors.push(error);
            }
        }
    }

    PipelineResult {
        success: errors.is_empty(),
        exit_code: Some(if errors.is_empty() { 0 } else { 1 }),
        working_dir: working_dir.display().to_string(),
        stdout,
        stderr,
        outputs,
        errors,
    }
}

fn run_stage_a_with_progress(
    app: &AppHandle,
    entries: &[(String, Vec<PathBuf>)],
    database_dir: &str,
    working_dir: &str,
    run_dir: &Path,
    progress: &StageProgress<'_>,
) -> PipelineResult {
    let total = entries.len().max(1);
    let stage = "Stage A · mpileup/一致序列";
    let mut combined = PipelineResult {
        success: true,
        exit_code: Some(0),
        working_dir: working_dir.to_owned(),
        stdout: String::new(),
        stderr: String::new(),
        outputs: Vec::new(),
        errors: Vec::new(),
    };

    for (index, (name, paths)) in entries.iter().enumerate() {
        let start_percent = progress.band_start
            + (progress.band_end - progress.band_start) * (index as f64 / total as f64);
        emit_pipeline_progress(
            app,
            progress.mode,
            progress.step,
            progress.total_steps,
            start_percent,
            stage,
            format!("样本 {}/{}：{name} 处理中…", index + 1, entries.len()),
            false,
        );

        let list_path = run_dir.join(format!("stage-a-{}.list", index + 1));
        let single = [(name.clone(), paths.clone())];
        let sample_result = match write_normalized_list(&list_path, &single) {
            Ok(()) => run_pipeline_blocking(
                app.clone(),
                PipelineRequest {
                    sample_list: list_path.display().to_string(),
                    database_dir: database_dir.to_owned(),
                    working_dir: working_dir.to_owned(),
                },
            ),
            Err(error) => PipelineResult {
                success: false,
                exit_code: Some(1),
                working_dir: working_dir.to_owned(),
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors: vec![format!("{name}: {error}")],
            },
        };

        let sample_success = sample_result.success;
        let sample_error = sample_result.errors.first().cloned();
        merge_stage_result(
            &mut combined,
            &format!("Stage A · 样本 {}/{} · {name}", index + 1, entries.len()),
            sample_result,
        );

        let done_percent = progress.band_start
            + (progress.band_end - progress.band_start) * ((index + 1) as f64 / total as f64);
        let message = if sample_success {
            format!("样本 {}/{}：{name} 完成", index + 1, entries.len())
        } else {
            format!(
                "样本 {}/{}：{name} 失败：{}",
                index + 1,
                entries.len(),
                sample_error.unwrap_or_else(|| "未知错误".to_owned())
            )
        };
        emit_pipeline_progress(
            app,
            progress.mode,
            progress.step,
            progress.total_steps,
            done_percent,
            stage,
            message,
            false,
        );
    }

    combined.success = combined.errors.is_empty();
    combined.exit_code = Some(if combined.success { 0 } else { 1 });
    combined
}

fn low_depth_from_estimated_files(
    working_dir: &Path,
    entries: &[(String, Vec<PathBuf>)],
) -> Option<bool> {
    let mut found = false;
    let mut low = false;
    for (name, _) in entries {
        let path = working_dir
            .join("YC-result")
            .join(name)
            .join("estimated_depth.txt");
        if let Ok(text) = std::fs::read_to_string(&path) {
            if let Ok(depth) = text.trim().parse::<u32>() {
                found = true;
                if depth < 50 {
                    low = true;
                }
            }
        }
    }
    if found {
        Some(low)
    } else {
        None
    }
}

fn run_pipeline_mode_blocking(app: AppHandle, request: PipelineModeRequest) -> PipelineResult {
    let raw_list = Path::new(&request.list);
    let list_abs = std::fs::canonicalize(raw_list).unwrap_or_else(|_| raw_list.to_path_buf());
    let mut resolved_mode = request.mode.trim().to_owned();
    if resolved_mode.is_empty() {
        resolved_mode = infer_pipeline_mode(&list_abs).unwrap_or_default();
    }
    let mode = resolved_mode.as_str();
    let list_base = list_abs
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    let output_value = expand_tilde(request.output_dir.trim());
    let output_path = if output_value.is_empty() {
        list_base.clone()
    } else {
        resolve_list_input(&list_base, &output_value)
    };
    let working_dir = std::fs::canonicalize(&output_path).unwrap_or_else(|_| output_path.clone());
    let working_dir_text = working_dir.display().to_string();
    pipeline_log(&format!(
        "resolved mode={mode} list={} workdir={working_dir_text}",
        list_abs.display()
    ));
    emit_pipeline_progress(&app, mode, 0, 3, 0.0, "准备", "正在校验输入…", false);

    let pipeline_dir = match pipeline_resource_dir(&app) {
        Ok(path) => path,
        Err(error) => {
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: working_dir_text,
                stdout: String::new(),
                stderr: String::new(),
                outputs: Vec::new(),
                errors: vec![error],
            };
        }
    };

    let mut errors = Vec::new();
    if !is_pipeline_mode(mode) {
        if mode.is_empty() {
            errors.push(
                "无法从列表内容判断起点模式：列表中的路径需为 fastq、BAM 或 2lie-mp-samtool"
                    .to_owned(),
            );
        } else {
            errors.push(format!("unknown pipeline mode: {mode}"));
        }
    }

    let database_dir = Path::new(&request.database_dir);
    let mut cutoff = request.cutoff.unwrap_or(1000);
    if !database_dir.is_dir() {
        errors.push(format!(
            "database directory does not exist: {}",
            request.database_dir
        ));
    } else {
        for name in ["ref.fa", "add", "G0.fa"] {
            if !database_dir.join(name).is_file() {
                errors.push(format!("database directory is missing {name}"));
            }
        }
        match resolve_cluster_cutoff(database_dir, request.cutoff) {
            Ok(value) => cutoff = value,
            Err(error) => errors.push(error),
        }
    }

    match mode {
        "reads2Result" => {
            if !pipeline_dir.join("align_sort.sh").is_file() {
                errors.push("align_sort.sh does not exist".to_owned());
            }
            let rammap = align_rammap_binary(&pipeline_dir);
            if !is_executable_file(&rammap) {
                errors.push(format!(
                    "rammap is missing or not executable: {}",
                    rammap.display()
                ));
            }
            let sambamba = pipeline_dir.join("sambamba-master/bin/sambamba-1.0.1");
            if !is_executable_file(&sambamba) {
                errors.push(format!(
                    "sambamba is missing or not executable: {}",
                    sambamba.display()
                ));
            }
        }
        "bam2Result" => {
            for script in ["SS.sh", "slice2.sh"] {
                if !pipeline_dir.join(script).is_file() {
                    errors.push(format!("{script} does not exist"));
                }
            }
        }
        "2lie-mp-samtool2Result" => {
            if !pipeline_dir.join("get_site.pl").is_file() {
                errors.push("get_site.pl does not exist".to_owned());
            }
        }
        _ => {}
    }

    if request.threads == Some(0) {
        errors.push("threads must be greater than 0".to_owned());
    }
    if request.jobs == Some(0) {
        errors.push("parallel sample count must be greater than 0".to_owned());
    }

    let entries = match parse_pipeline_mode_list(&list_abs, mode) {
        Ok(entries) => entries,
        Err(list_errors) => {
            errors.extend(list_errors);
            Vec::new()
        }
    };

    let run_dir = working_dir.join(".strainscale-run");
    if let Err(error) = std::fs::create_dir_all(&run_dir) {
        errors.push(format!("failed to create run directory: {error}"));
    }

    let total_steps = if mode == "reads2Result" { 5 } else { 4 };

    if !errors.is_empty() {
        let result = PipelineResult {
            success: false,
            exit_code: Some(1),
            working_dir: working_dir_text,
            stdout: String::new(),
            stderr: String::new(),
            outputs: Vec::new(),
            errors,
        };
        return finish_pipeline_progress(&app, mode, total_steps, 0.0, "输入校验", result);
    }

    pipeline_log(&format!("resolved cluster cutoff={cutoff}"));
    emit_pipeline_progress(
        &app,
        mode,
        0,
        total_steps,
        0.0,
        "准备",
        format!(
            "输入校验通过，共 {} 个样本，聚类 cutoff={cutoff}，开始运行",
            entries.len()
        ),
        false,
    );

    let mut combined = PipelineResult {
        success: true,
        exit_code: Some(0),
        working_dir: working_dir_text.clone(),
        stdout: String::new(),
        stderr: String::new(),
        outputs: Vec::new(),
        errors: Vec::new(),
    };

    let samples_list: PathBuf = match mode {
        "reads2Result" => {
            let reads_list = run_dir.join("reads.list");
            if let Err(error) = write_normalized_list(&reads_list, &entries) {
                combined.errors.push(error);
                combined.success = false;
                combined.exit_code = Some(1);
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    0.0,
                    "Stage 0 · 比对建 BAM",
                    combined,
                );
            }
            let stage0 = "Stage 0 · 比对建 BAM";
            let (band_start, band_end) = stage_band(1, total_steps);
            emit_pipeline_progress(
                &app,
                mode,
                1,
                total_steps,
                band_start,
                stage0,
                format!("FASTQ 比对建 BAM 开始，共 {} 个样本", entries.len()),
                false,
            );
            let bam_dir = working_dir.join("bam");
            let align = run_align_blocking(
                app.clone(),
                AlignStageRequest {
                    fastq_list: reads_list.display().to_string(),
                    reference: database_dir.join("ref.fa").display().to_string(),
                    output_dir: bam_dir.display().to_string(),
                    threads: Some(request.threads.unwrap_or(40)),
                    jobs: Some(request.jobs.unwrap_or(8)),
                },
            );
            if !merge_stage_result(&mut combined, "reads2Result · 比对建 BAM", align) {
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    band_end,
                    stage0,
                    combined,
                );
            }
            let bam_samples = bam_dir.join("sample.list");
            let bam_entries = match parse_pipeline_mode_list(&bam_samples, "bam2Result") {
                Ok(list) => list,
                Err(list_errors) => {
                    combined.errors.extend(list_errors);
                    combined.success = false;
                    combined.exit_code = Some(1);
                    return finish_pipeline_progress(
                        &app,
                        mode,
                        total_steps,
                        band_end,
                        stage0,
                        combined,
                    );
                }
            };
            emit_pipeline_progress(
                &app,
                mode,
                1,
                total_steps,
                band_end,
                stage0,
                format!("BAM 构建完成，共 {} 个样本", bam_entries.len()),
                false,
            );

            let (band_start, band_end) = stage_band(2, total_steps);
            let stage_a_progress = StageProgress {
                mode,
                step: 2,
                total_steps,
                band_start,
                band_end,
            };
            emit_pipeline_progress(
                &app,
                mode,
                2,
                total_steps,
                band_start,
                "Stage A · mpileup/一致序列",
                format!("开始生成一致序列，共 {} 个样本", bam_entries.len()),
                false,
            );
            let stage_a = run_stage_a_with_progress(
                &app,
                &bam_entries,
                &request.database_dir,
                &working_dir_text,
                &run_dir,
                &stage_a_progress,
            );
            if !merge_stage_result(&mut combined, "Stage A · mpileup/一致序列", stage_a) {
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    band_end,
                    "Stage A · mpileup/一致序列",
                    combined,
                );
            }
            emit_pipeline_progress(
                &app,
                mode,
                2,
                total_steps,
                band_end,
                "Stage A · mpileup/一致序列",
                "一致序列生成完成",
                false,
            );
            bam_samples
        }
        "bam2Result" => {
            let bam_list = run_dir.join("bam.list");
            if let Err(error) = write_normalized_list(&bam_list, &entries) {
                combined.errors.push(error);
                combined.success = false;
                combined.exit_code = Some(1);
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    0.0,
                    "Stage A · mpileup/一致序列",
                    combined,
                );
            }
            let (band_start, band_end) = stage_band(1, total_steps);
            let stage_a_progress = StageProgress {
                mode,
                step: 1,
                total_steps,
                band_start,
                band_end,
            };
            emit_pipeline_progress(
                &app,
                mode,
                1,
                total_steps,
                band_start,
                "Stage A · mpileup/一致序列",
                format!("开始生成一致序列，共 {} 个样本", entries.len()),
                false,
            );
            let stage_a = run_stage_a_with_progress(
                &app,
                &entries,
                &request.database_dir,
                &working_dir_text,
                &run_dir,
                &stage_a_progress,
            );
            if !merge_stage_result(&mut combined, "Stage A · mpileup/一致序列", stage_a) {
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    band_end,
                    "Stage A · mpileup/一致序列",
                    combined,
                );
            }
            emit_pipeline_progress(
                &app,
                mode,
                1,
                total_steps,
                band_end,
                "Stage A · mpileup/一致序列",
                "一致序列生成完成",
                false,
            );
            bam_list
        }
        "2lie-mp-samtool2Result" => {
            let pileup_list = run_dir.join("pileup.list");
            if let Err(error) = write_normalized_list(&pileup_list, &entries) {
                combined.errors.push(error);
                combined.success = false;
                combined.exit_code = Some(1);
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    0.0,
                    "Stage A' · get_site 一致序列",
                    combined,
                );
            }
            let (band_start, band_end) = stage_band(1, total_steps);
            let import_progress = StageProgress {
                mode,
                step: 1,
                total_steps,
                band_start,
                band_end,
            };
            emit_pipeline_progress(
                &app,
                mode,
                1,
                total_steps,
                band_start,
                "Stage A' · get_site 一致序列",
                format!("开始生成一致序列，共 {} 个样本", entries.len()),
                false,
            );
            let imports = import_pileup_samples(
                &app,
                &pipeline_dir,
                &entries,
                &working_dir,
                &import_progress,
            );
            if !merge_stage_result(&mut combined, "Stage A' · get_site 一致序列", imports) {
                return finish_pipeline_progress(
                    &app,
                    mode,
                    total_steps,
                    band_end,
                    "Stage A' · get_site 一致序列",
                    combined,
                );
            }
            emit_pipeline_progress(
                &app,
                mode,
                1,
                total_steps,
                band_end,
                "Stage A' · get_site 一致序列",
                "一致序列生成完成",
                false,
            );
            pileup_list
        }
        _ => unreachable!(),
    };

    let max_dist = match low_depth_from_estimated_files(&working_dir, &entries) {
        Some(true) => Some(20000),
        Some(false) => Some(10000),
        None => None,
    };

    let stage_b_step = if mode == "reads2Result" { 3 } else { 2 };
    let (band_start, band_end) = stage_band(stage_b_step, total_steps);
    emit_pipeline_progress(
        &app,
        mode,
        stage_b_step,
        total_steps,
        band_start,
        "Stage B · 原生 mbgc",
        format!("开始检索近缘菌株，共 {} 个样本", entries.len()),
        false,
    );
    let stage_b = run_mbgc_blocking(
        app.clone(),
        MbgcStageRequest {
            database_dir: request.database_dir.clone(),
            samples_list: samples_list.display().to_string(),
            working_dir: working_dir_text.clone(),
            pick_x: None,
            max_jobs: None,
            max_dist,
            cutoff: Some(cutoff),
        },
    );
    if !merge_stage_result(&mut combined, "Stage B · 原生 mbgc", stage_b) {
        return finish_pipeline_progress(
            &app,
            mode,
            total_steps,
            band_end,
            "Stage B · 原生 mbgc",
            combined,
        );
    }
    emit_pipeline_progress(
        &app,
        mode,
        stage_b_step,
        total_steps,
        band_end,
        "Stage B · 原生 mbgc",
        "近缘菌株检索完成",
        false,
    );

    let stage_c_step = total_steps - 1;
    let (band_start, band_end) = stage_band(stage_c_step, total_steps);
    let hybrid_progress = StageProgress {
        mode,
        step: stage_c_step,
        total_steps,
        band_start,
        band_end,
    };
    emit_pipeline_progress(
        &app,
        mode,
        stage_c_step,
        total_steps,
        band_start,
        "Stage C · 原生 hybrid",
        format!("开始代表株聚类与精炼，共 {} 个样本", entries.len()),
        false,
    );
    let hybrid = run_hybrid_batch_blocking(
        app.clone(),
        HybridBatchRequest {
            working_dir: working_dir_text.clone(),
            samples_list: samples_list.display().to_string(),
            cutoff: Some(cutoff),
            max_dist,
            match_threshold: None,
        },
        Some(&hybrid_progress),
    );
    let hybrid_stdout = format!(
        "ran: {}\nskipped: {}\nfailed: {}",
        hybrid.ran.join(", "),
        hybrid.skipped.join(", "),
        hybrid.failed.join(", ")
    );
    let mut hybrid_errors = hybrid.failed.clone();
    if !hybrid.success {
        hybrid_errors.extend(hybrid.skipped.clone());
        if hybrid_errors.is_empty() {
            hybrid_errors.push("hybrid produced no result".to_owned());
        }
    }
    let stage_c = PipelineResult {
        success: hybrid.success,
        exit_code: Some(if hybrid.success { 0 } else { 1 }),
        working_dir: working_dir_text.clone(),
        stdout: hybrid_stdout,
        stderr: String::new(),
        outputs: hybrid.outputs,
        errors: hybrid_errors,
    };
    let stage_c_success = stage_c.success;
    merge_stage_result(&mut combined, "Stage C · 原生 hybrid", stage_c);
    emit_pipeline_progress(
        &app,
        mode,
        stage_c_step,
        total_steps,
        band_end,
        "Stage C · 原生 hybrid",
        if stage_c_success {
            format!("代表株精炼完成，共 {} 个样本", entries.len())
        } else {
            format!(
                "代表株精炼结束，部分样本失败（共 {} 个样本）",
                entries.len()
            )
        },
        false,
    );

    let (band_start, band_end) = stage_band(total_steps, total_steps);
    let abundance_progress = StageProgress {
        mode,
        step: total_steps,
        total_steps,
        band_start,
        band_end,
    };
    emit_pipeline_progress(
        &app,
        mode,
        total_steps,
        total_steps,
        band_start,
        "abundance · 丰度筛选",
        format!("开始丰度筛选，共 {} 个样本", entries.len()),
        false,
    );
    let stage_d = run_abundance_stage_blocking(
        &app,
        &pipeline_dir,
        &working_dir,
        &samples_list,
        cutoff,
        &abundance_progress,
    );
    let stage_d_success = stage_d.success;
    merge_stage_result(&mut combined, "abundance · 丰度筛选", stage_d);
    emit_pipeline_progress(
        &app,
        mode,
        total_steps,
        total_steps,
        band_end,
        "abundance · 丰度筛选",
        if stage_d_success {
            format!("丰度筛选完成，共 {} 个样本", entries.len())
        } else {
            format!(
                "丰度筛选结束，部分样本失败（共 {} 个样本）",
                entries.len()
            )
        },
        false,
    );

    finish_pipeline_progress(&app, mode, total_steps, 100.0, "完成", combined)
}

fn pick_argument_string(arguments: &Value, keys: &[&str]) -> Option<String> {
    for key in keys {
        if let Some(value) = arguments.get(*key).and_then(Value::as_str) {
            let trimmed = value.trim();
            if !trimmed.is_empty() {
                return Some(trimmed.to_owned());
            }
        }
    }
    None
}

fn pick_argument_u32(arguments: &Value, keys: &[&str]) -> Result<Option<u32>, String> {
    for key in keys {
        let Some(value) = arguments.get(*key) else {
            continue;
        };
        if value.is_null() {
            continue;
        }
        if let Some(number) = value.as_u64() {
            return u32::try_from(number)
                .map(Some)
                .map_err(|_| format!("{key} 超出 u32 范围"));
        }
        if let Some(text) = value.as_str() {
            let trimmed = text.trim();
            if trimmed.is_empty() {
                continue;
            }
            return trimmed
                .parse::<u32>()
                .map(Some)
                .map_err(|_| format!("{key} 必须是非负整数"));
        }
        return Err(format!("{key} 必须是非负整数"));
    }
    Ok(None)
}

fn cutoff_from_message(message: &str) -> Option<u32> {
    let lower = message.to_lowercase();
    for key in ["cutoff", "聚类"] {
        let mut search_start = 0;
        while search_start < lower.len() {
            let Some(relative) = lower[search_start..].find(key) else {
                break;
            };
            let value_start = search_start + relative + key.len();
            search_start = value_start;
            let rest = lower[value_start..].trim_start_matches(|character: char| {
                character.is_whitespace() || matches!(character, '=' | ':' | '：')
            });
            let rest = rest.strip_prefix("阈值").unwrap_or(rest);
            let rest = rest.trim_start_matches(|character: char| {
                character.is_whitespace() || matches!(character, '=' | ':' | '：')
            });
            let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
            if let Ok(value) = digits.parse::<u32>() {
                return Some(value);
            }
        }
    }
    None
}

fn infer_pipeline_mode(list: &Path) -> Option<String> {
    let content = std::fs::read_to_string(list).ok()?;
    let line = content.lines().find(|raw| {
        let trimmed = raw.trim();
        !trimmed.is_empty() && !trimmed.starts_with('#')
    })?;
    let fields: Vec<&str> = if line.contains('\t') {
        line.split('\t').collect()
    } else {
        line.split_whitespace().collect()
    };
    if fields.len() < 2 {
        return None;
    }
    let paths: Vec<String> = fields[1..]
        .iter()
        .map(|value| value.trim().to_lowercase())
        .collect();
    if paths
        .iter()
        .any(|path| path.contains("2lie-mp-samtool") || path.contains("mp-samtool"))
    {
        return Some("2lie-mp-samtool2Result".to_owned());
    }
    if paths
        .iter()
        .any(|path| path.ends_with(".bam") || path.ends_with(".cram") || path.ends_with(".sam"))
    {
        return Some("bam2Result".to_owned());
    }
    if paths.iter().any(|path| {
        path.ends_with(".fastq")
            || path.ends_with(".fq")
            || path.ends_with(".fastq.gz")
            || path.ends_with(".fq.gz")
    }) {
        return Some("reads2Result".to_owned());
    }
    None
}

fn labeled_path(message: &str, keys: &[&str]) -> Option<PathBuf> {
    let lower = message.to_lowercase();
    for key in keys {
        for separator in ["=", ":", "："] {
            let pattern = format!("{key}{separator}");
            if let Some(index) = lower.find(&pattern) {
                if !message.is_char_boundary(index + pattern.len()) {
                    continue;
                }
                let rest = &message[index + pattern.len()..];
                let token = rest
                    .trim_start()
                    .split(|c: char| {
                        c.is_whitespace()
                            || matches!(c, '，' | '。' | '、' | '；' | ',' | ';' | '"' | '\'')
                    })
                    .next()
                    .unwrap_or("")
                    .trim();
                if !token.is_empty() {
                    return Some(PathBuf::from(token));
                }
            }
        }
    }
    None
}

fn plan_pipeline_from_message(message: &str) -> Result<PipelinePlan, String> {
    let mut list = labeled_path(message, &["sample_list", "list", "样本列表", "列表"]);
    let mut database = labeled_path(message, &["database_dir", "database", "db", "数据库"]);
    let mut output = labeled_path(
        message,
        &["output_dir", "output", "out", "输出目录", "输出"],
    );
    if output.is_none() && (message.contains("桌面") || message.to_lowercase().contains("desktop"))
    {
        if let Some(home) = std::env::var_os("HOME") {
            output = Some(PathBuf::from(home).join("Desktop"));
        }
    }

    for raw in extract_files(message) {
        let path = PathBuf::from(&raw);
        if path.is_file() {
            let is_list_extension = path
                .extension()
                .and_then(|value| value.to_str())
                .map(|value| value.eq_ignore_ascii_case("list"))
                .unwrap_or(false);
            if is_list_extension || infer_pipeline_mode(&path).is_some() {
                if list.is_none() {
                    list = Some(path);
                }
                continue;
            }
        }
        if path.is_dir() {
            if path.join("ref.fa").is_file() && path.join("add").is_file() {
                database = Some(path);
                continue;
            }
            if output.is_none() {
                output = Some(path);
            }
            continue;
        }
        if output.is_none() {
            output = Some(path);
        }
    }

    let list = list.ok_or_else(|| {
        "没有在消息中找到样本列表文件（.list 或包含 fastq/bam/pileup 路径的列表）".to_owned()
    })?;
    let database = database
        .ok_or_else(|| "没有在消息中找到参考数据库目录（需包含 ref.fa 和 add）".to_owned())?;
    let mode = infer_pipeline_mode(&list).ok_or_else(|| {
        "无法从列表文件内容判断起点模式，请使用 .list 并在其中写入 fastq / BAM / 2lie-mp-samtool 路径".to_owned()
    })?;
    let output_dir = output
        .map(|path| path.display().to_string())
        .unwrap_or_default();

    Ok(PipelinePlan {
        mode,
        database_dir: database.display().to_string(),
        list: list.display().to_string(),
        output_dir,
        cutoff: cutoff_from_message(message),
        reason: "未使用 LLM：优先读取 DB/list/output 标签及聚类/cutoff 参数，再按路径规则识别数据库、列表与输出目录，并根据列表内容选择起点模式。"
            .to_owned(),
    })
}

fn run_agent_tool_blocking(
    app: AppHandle,
    request: AgentToolRequest,
) -> Result<AgentToolOutcome, String> {
    if request.name != "strain_pipeline" {
        return Err(format!("Skill 尚未接入执行器: {}", request.name));
    }

    let arguments = &request.arguments;
    let database_dir = pick_argument_string(arguments, &["database_dir", "databaseDir"])
        .ok_or_else(|| "缺少 database_dir（参考数据库目录）".to_owned())?;
    let list = pick_argument_string(arguments, &["list", "sample_list", "sampleList"])
        .ok_or_else(|| "缺少 list（样本列表文件）".to_owned())?;
    let mut output_dir =
        pick_argument_string(arguments, &["output_dir", "outputDir"]).unwrap_or_default();
    output_dir = expand_tilde(&output_dir);
    let message_lower = request.message.to_lowercase();
    let mentions_desktop = request.message.contains("桌面") || message_lower.contains("desktop");
    if mentions_desktop
        && !output_dir.to_lowercase().contains("desktop")
        && !output_dir.contains("桌面")
    {
        if let Some(home) = std::env::var_os("HOME") {
            output_dir = PathBuf::from(home).join("Desktop").display().to_string();
        }
    }
    let provided_mode =
        pick_argument_string(arguments, &["mode"]).filter(|mode| is_pipeline_mode(mode));
    let inferred_mode = infer_pipeline_mode(Path::new(&list));
    let mode = inferred_mode.or(provided_mode).ok_or_else(|| {
        "无法从列表文件判断起点模式，请显式传入 mode（reads2Result / bam2Result / 2lie-mp-samtool2Result）".to_owned()
    })?;

    let cutoff = pick_argument_u32(arguments, &["cutoff", "cluster_cutoff", "clusterCutoff"])?
        .or_else(|| cutoff_from_message(&request.message));
    let result = run_pipeline_mode_blocking(
        app,
        PipelineModeRequest {
            mode: mode.clone(),
            database_dir,
            list,
            output_dir,
            threads: None,
            jobs: None,
            cutoff,
        },
    );

    Ok(AgentToolOutcome { mode, result })
}

#[tauri::command]
fn plan_pipeline(message: String) -> Result<PipelinePlan, String> {
    plan_pipeline_from_message(&message)
}

#[tauri::command]
async fn run_agent_tool(
    app: AppHandle,
    request: AgentToolRequest,
) -> Result<AgentToolOutcome, String> {
    tauri::async_runtime::spawn_blocking(move || run_agent_tool_blocking(app, request))
        .await
        .map_err(|error| format!("agent tool task failed: {error}"))?
}

#[tauri::command]
fn cluster_cutoffs(database_dir: String) -> Result<Vec<u32>, String> {
    collect_cluster_cutoffs(Path::new(&database_dir))
}

#[tauri::command]
fn agent_route(request: AgentCallRequest) -> RouteResult {
    classify(&request.message)
}

#[tauri::command]
async fn run_strainscale_pipeline(
    app: AppHandle,
    request: PipelineRequest,
) -> Result<PipelineResult, String> {
    tauri::async_runtime::spawn_blocking(move || run_pipeline_blocking(app, request))
        .await
        .map_err(|e| format!("pipeline task failed: {e}"))
}

#[tauri::command]
async fn run_mbgc_stage(
    app: AppHandle,
    request: MbgcStageRequest,
) -> Result<PipelineResult, String> {
    tauri::async_runtime::spawn_blocking(move || run_mbgc_blocking(app, request))
        .await
        .map_err(|e| format!("mbgc task failed: {e}"))
}

#[tauri::command]
async fn run_hybrid_stage(
    app: AppHandle,
    request: HybridStageRequest,
) -> Result<PipelineResult, String> {
    tauri::async_runtime::spawn_blocking(move || run_hybrid_blocking(app, request))
        .await
        .map_err(|e| format!("hybrid task failed: {e}"))
}

#[tauri::command]
async fn run_align_stage(
    app: AppHandle,
    request: AlignStageRequest,
) -> Result<PipelineResult, String> {
    tauri::async_runtime::spawn_blocking(move || run_align_blocking(app, request))
        .await
        .map_err(|e| format!("align task failed: {e}"))
}

#[tauri::command]
async fn run_pipeline_mode(
    app: AppHandle,
    request: PipelineModeRequest,
) -> Result<PipelineResult, String> {
    pipeline_log(&format!(
        "run start list={} db={} mode={}",
        request.list, request.database_dir, request.mode
    ));
    let result = tauri::async_runtime::spawn_blocking({
        let app = app.clone();
        move || run_pipeline_mode_blocking(app, request)
    })
    .await
    .map_err(|e| format!("pipeline task failed: {e}"));
    match &result {
        Ok(outcome) => pipeline_log(&format!(
            "run finished success={} exit={:?} errors={}",
            outcome.success,
            outcome.exit_code,
            outcome.errors.join(" | ")
        )),
        Err(error) => pipeline_log(&format!("run failed to start: {error}")),
    }
    result
}

fn sample_names_from_list(path: &Path) -> Result<Vec<String>, String> {
    let content =
        std::fs::read_to_string(path).map_err(|e| format!("failed to read samples list: {e}"))?;
    let mut names = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some(name) = trimmed.split_whitespace().next() {
            if !names.iter().any(|n: &String| n == name) {
                names.push(name.to_owned());
            }
        }
    }
    Ok(names)
}

fn strain_result_dir_name(cutoff: u32) -> String {
    format!("Strain-{cutoff}")
}

fn run_hybrid_batch_blocking(
    app: AppHandle,
    request: HybridBatchRequest,
    progress: Option<&StageProgress<'_>>,
) -> HybridBatchResult {
    let mut result = HybridBatchResult {
        success: false,
        ran: Vec::new(),
        skipped: Vec::new(),
        failed: Vec::new(),
        outputs: Vec::new(),
    };
    let samples_list = Path::new(&request.samples_list);
    if !samples_list.is_file() {
        result.skipped.push(format!(
            "samples list does not exist: {}",
            samples_list.display()
        ));
        return result;
    }
    let names = match sample_names_from_list(samples_list) {
        Ok(names) => names,
        Err(error) => {
            result.skipped.push(error);
            return result;
        }
    };
    if names.is_empty() {
        result.skipped.push("samples list is empty".to_owned());
        return result;
    }
    let base = Path::new(&request.working_dir).join("YC-result");
    let cutoff = request.cutoff.unwrap_or(1000);
    let max_dist = request.max_dist.unwrap_or(10000);
    let match_threshold = request.match_threshold.unwrap_or(0.9);

    let total = names.len().max(1);
    let stage = "Stage C · 原生 hybrid";
    for (index, name) in names.iter().enumerate() {
        if let Some(progress) = progress {
            let percent = progress.band_start
                + (progress.band_end - progress.band_start) * (index as f64 / total as f64);
            emit_pipeline_progress(
                &app,
                progress.mode,
                progress.step,
                progress.total_steps,
                percent,
                stage,
                format!("样本 {}/{}：{name} 精炼中…", index + 1, names.len()),
                false,
            );
        }
        let sample_dir = base.join(name);
        let sample_output_dir = sample_dir.join("output");
        let filter_counts = sample_dir.join("filter-result-with-counts.txt");
        if !sample_output_dir.is_dir() || !filter_counts.is_file() {
            result.skipped.push(format!(
                "{name}: 缺少 {}",
                if !sample_output_dir.is_dir() {
                    "output/（Stage A 未完成）"
                } else {
                    "filter-result-with-counts.txt（Stage A 未完成）"
                }
            ));
            continue;
        }
        let single = HybridStageRequest {
            sample_output_dir: sample_output_dir.to_string_lossy().to_string(),
            filter_counts: filter_counts.to_string_lossy().to_string(),
            result_dir: sample_dir
                .join(strain_result_dir_name(cutoff))
                .to_string_lossy()
                .to_string(),
            cutoff: Some(cutoff),
            max_dist: Some(max_dist),
            match_threshold: Some(match_threshold),
        };
        let outcome = run_hybrid_blocking(app.clone(), single);
        let success = outcome.success;
        if success {
            result.ran.push(name.clone());
            result.outputs.extend(outcome.outputs);
        } else {
            result
                .failed
                .push(format!("{name}: {}", outcome.errors.join("; ")));
        }
        if let Some(progress) = progress {
            let done_percent = progress.band_start
                + (progress.band_end - progress.band_start) * ((index + 1) as f64 / total as f64);
            let message = if success {
                format!("样本 {}/{}：{name} 完成", index + 1, names.len())
            } else {
                format!("样本 {}/{}：{name} 失败", index + 1, names.len())
            };
            emit_pipeline_progress(
                &app,
                progress.mode,
                progress.step,
                progress.total_steps,
                done_percent,
                stage,
                message,
                false,
            );
        }
    }

    result.success = result.failed.is_empty() && !result.ran.is_empty();
    result
}

const ABUNDANCE_PARAMS: [&str; 5] = ["0.05", "0.8", "20", "0.05", "200"];

fn run_abundance_stage_blocking(
    app: &AppHandle,
    pipeline_dir: &Path,
    working_dir: &Path,
    samples_list: &Path,
    cutoff: u32,
    progress: &StageProgress<'_>,
) -> PipelineResult {
    let stage = "abundance · 丰度筛选";
    let result_dir_name = strain_result_dir_name(cutoff);
    let working_dir_text = working_dir.display().to_string();
    let binary = pipeline_dir.join("server-deploy/abundance-macos-arm64");
    let mut stdout = String::new();
    let mut stderr = String::new();
    let mut outputs = Vec::new();
    let mut errors = Vec::new();
    let mut ran = Vec::new();
    let mut skipped = Vec::new();

    if !is_executable_file(&binary) {
        errors.push(format!(
            "native abundance binary is missing or not executable: {}",
            binary.display()
        ));
        return PipelineResult {
            success: false,
            exit_code: None,
            working_dir: working_dir_text,
            stdout,
            stderr,
            outputs,
            errors,
        };
    }

    let names = match sample_names_from_list(samples_list) {
        Ok(names) => names,
        Err(error) => {
            errors.push(error);
            return PipelineResult {
                success: false,
                exit_code: None,
                working_dir: working_dir_text,
                stdout,
                stderr,
                outputs,
                errors,
            };
        }
    };

    let base = working_dir.join("YC-result");
    let total = names.len().max(1);
    for (index, name) in names.iter().enumerate() {
        let percent = progress.band_start
            + (progress.band_end - progress.band_start) * (index as f64 / total as f64);
        emit_pipeline_progress(
            app,
            progress.mode,
            progress.step,
            progress.total_steps,
            percent,
            stage,
            format!("样本 {}/{}：{name} 丰度筛选开始", index + 1, names.len()),
            false,
        );

        let sample_dir = base.join(name);
        let filter_counts = sample_dir.join("filter-result-with-counts.txt");
        let sites = sample_dir.join(&result_dir_name).join("final_specific_sites.tsv");
        if !filter_counts.is_file() || !sites.is_file() {
            skipped.push(format!(
                "{name}: 缺少 {}",
                if !filter_counts.is_file() {
                    "filter-result-with-counts.txt（Stage A 未完成）"
                } else {
                    "final_specific_sites.tsv（Stage C 未完成）"
                }
            ));
            continue;
        }

        let out_dir = sample_dir.join(&result_dir_name).join("abundance");
        if let Err(error) = std::fs::create_dir_all(&out_dir) {
            errors.push(format!(
                "{name}: failed to create abundance output directory: {error}"
            ));
            continue;
        }
        let output_file = out_dir.join("abundance-result.tsv");
        let outcome = ProcessCommand::new(&binary)
            .current_dir(&out_dir)
            .arg(&filter_counts)
            .arg(&sites)
            .arg(&output_file)
            .args(ABUNDANCE_PARAMS)
            .output();
        match outcome {
            Ok(outcome) => {
                let out_text = String::from_utf8_lossy(&outcome.stdout);
                let err_text = String::from_utf8_lossy(&outcome.stderr);
                if !out_text.trim().is_empty() {
                    stdout.push_str(&format!("[{name}]\n{}", out_text.trim_end()));
                    stdout.push('\n');
                }
                if !err_text.trim().is_empty() {
                    stderr.push_str(&format!("[{name}]\n{}", err_text.trim_end()));
                    stderr.push('\n');
                }
                if !outcome.status.success() || !output_file.is_file() {
                    errors.push(format!(
                        "{name}: abundance 执行失败（{:?}）",
                        outcome.status
                    ));
                    continue;
                }
                ran.push(name.clone());
                outputs.push(output_file.display().to_string());
                let selected = out_dir.join("abundance-result.tsv.selected");
                if selected.is_file() {
                    outputs.push(selected.display().to_string());
                }
                let plots = out_dir.join("peak_plots");
                if let Ok(entries) = std::fs::read_dir(&plots) {
                    for entry in entries.flatten() {
                        let path = entry.path();
                        if path.is_file() {
                            outputs.push(path.display().to_string());
                        }
                    }
                }
            }
            Err(error) => {
                errors.push(format!("{name}: failed to run abundance: {error}"));
            }
        }
    }

    let success = errors.is_empty() && !ran.is_empty();
    if !success {
        errors.extend(skipped);
        if errors.is_empty() {
            errors.push("abundance produced no result".to_owned());
        }
    }

    PipelineResult {
        success,
        exit_code: Some(if success { 0 } else { 1 }),
        working_dir: working_dir_text,
        stdout,
        stderr,
        outputs,
        errors,
    }
}

#[tauri::command]
async fn run_hybrid_batch(
    app: AppHandle,
    request: HybridBatchRequest,
) -> Result<HybridBatchResult, String> {
    tauri::async_runtime::spawn_blocking(move || run_hybrid_batch_blocking(app, request, None))
        .await
        .map_err(|e| format!("hybrid batch task failed: {e}"))
}
#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .manage(PipelineState::default())
        .invoke_handler(tauri::generate_handler![
            agent_route,
            run_strainscale_pipeline,
            run_align_stage,
            run_pipeline_mode,
            plan_pipeline,
            cluster_cutoffs,
            run_agent_tool,
            run_mbgc_stage,
            run_hybrid_stage,
            run_hybrid_batch,
            system_metrics,
            pipeline_status,
            open_external
        ])
        .run(tauri::generate_context!())
        .expect("error while running StrainScale Agent");
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    const MULTI_CLU: &str = "Strain-name\t10000\t5000\t1000\t100\t0\r\nstrain-a\t1\t2\t3\t4\t5\r\nstrain-b\t6\t7\t8\t9\t10\r\n";

    fn create_database(name: &str, clusters: &[(&str, &str)]) -> PathBuf {
        let root =
            std::env::temp_dir().join(format!("strainscale-agent-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(root.join("DB")).unwrap();
        let mut db_list = String::new();
        for (db, clu) in clusters {
            db_list.push_str(db);
            db_list.push('\n');
            std::fs::write(root.join("DB").join(format!("{db}.clu")), clu).unwrap();
            std::fs::write(root.join("DB").join(format!("{db}.zst")), "test").unwrap();
            std::fs::write(root.join("DB").join(format!("{db}.mbgc")), "test").unwrap();
        }
        std::fs::write(root.join("DB.list"), db_list).unwrap();
        root
    }

    #[test]
    fn parses_dynamic_and_legacy_cluster_headers() {
        assert_eq!(
            parse_cluster_header(MULTI_CLU).unwrap(),
            Some(vec![10000, 5000, 1000, 100, 0])
        );
        assert_eq!(parse_cluster_header("strain-a,12\n").unwrap(), None);
        assert!(parse_cluster_header("Strain-name\t1000\tbad\n").is_err());
        assert!(parse_cluster_header("Strain-name\t1000\t1000\n").is_err());
    }

    #[test]
    fn intersects_cutoffs_in_first_database_order() {
        let groups = vec![vec![10000, 1000, 0], vec![0, 1000], vec![1000, 0]];
        assert_eq!(common_cluster_cutoffs(&groups), vec![1000, 0]);
        assert!(common_cluster_cutoffs(&[vec![100], vec![1000]]).is_empty());
    }

    #[test]
    fn loads_and_validates_cutoffs_from_database_directory() {
        let root = create_database("dynamic", &[("vp", MULTI_CLU)]);
        assert_eq!(
            collect_cluster_cutoffs(&root).unwrap(),
            vec![10000, 5000, 1000, 100, 0]
        );
        assert_eq!(resolve_cluster_cutoff(&root, None).unwrap(), 1000);
        assert_eq!(resolve_cluster_cutoff(&root, Some(10000)).unwrap(), 10000);
        assert_eq!(resolve_cluster_cutoff(&root, Some(0)).unwrap(), 0);
        assert!(resolve_cluster_cutoff(&root, Some(42)).is_err());
        let _ = std::fs::remove_dir_all(root);
    }

    #[test]
    fn supports_legacy_and_multi_database_intersection() {
        let legacy = create_database("legacy", &[("vp", "strain-a,12\n")]);
        assert_eq!(collect_cluster_cutoffs(&legacy).unwrap(), vec![1000]);
        let _ = std::fs::remove_dir_all(legacy);

        let multiple = create_database(
            "multiple",
            &[
                ("vp", MULTI_CLU),
                ("other", "Strain-name\t0\t1000\nstrain-a\t1\t2\n"),
            ],
        );
        assert_eq!(collect_cluster_cutoffs(&multiple).unwrap(), vec![1000, 0]);
        let _ = std::fs::remove_dir_all(multiple);
    }

    #[test]
    fn converts_only_the_selected_cluster_column() {
        assert_eq!(
            selected_cluster_content(MULTI_CLU, 1000).unwrap(),
            "strain-a,3\nstrain-b,8\n"
        );
        assert_eq!(
            selected_cluster_content(MULTI_CLU, 0).unwrap(),
            "strain-a,5\nstrain-b,10\n"
        );
        assert!(selected_cluster_content(MULTI_CLU, 42).is_err());
        assert!(selected_cluster_content("Strain-name\t1000\na\t1\textra\n", 1000).is_err());
    }

    #[test]
    fn reads_cutoff_from_natural_language() {
        for (message, expected) in [
            ("聚类1000", Some(1000)),
            ("聚类 1000", Some(1000)),
            ("聚类=1000", Some(1000)),
            ("聚类阈值：100", Some(100)),
            ("cutoff 10000", Some(10000)),
            ("CUTOFF=5000", Some(5000)),
            ("cutoff：0", Some(0)),
            ("没有阈值", None),
        ] {
            assert_eq!(cutoff_from_message(message), expected, "{message}");
        }
    }

    #[test]
    fn accepts_numeric_and_string_tool_cutoffs() {
        assert_eq!(
            pick_argument_u32(&json!({"cutoff": 0}), &["cutoff"]).unwrap(),
            Some(0)
        );
        assert_eq!(
            pick_argument_u32(&json!({"cluster_cutoff": "10000"}), &["cluster_cutoff"]).unwrap(),
            Some(10000)
        );
        assert!(pick_argument_u32(&json!({"cutoff": "bad"}), &["cutoff"]).is_err());
    }
}
