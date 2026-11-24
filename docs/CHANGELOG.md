# CHANGELOG

## 未发布（建议 1.5.0）— BAM 读写栈迁移：rust-htslib → noodles

> 破坏性变更：库的公开类型 `BamRecord` 由 `rust_htslib::bam::Record` 变为
> `noodles_sam::alignment::RecordBuf`；`BamWriter`/`BamReader` 别名移除。
> 最低 Rust 版本（MSRV）提升到 **1.91**（noodles 0.99/0.94 为 edition 2024）。

### 变更

- **彻底移除 `rust-htslib`（及 `gskits`）**，BAM/SAM 读写改用纯 Rust 的
  `noodles-bam 0.99` / `noodles-sam 0.94`（+ `noodles-core`）。
  - 构建不再依赖 bindgen / libclang / hts-sys / 系统 HTSlib 头文件，
    `docs/ci-htslib-opaque-structs.md` 记录的 CI 失败模式（opaque structs / E0609）
    从根上消除，CI 中的 `libhts-dev` / `libclang` 安装与诊断步骤可一并清理。
  - 随 gskits 一并移除的传递依赖：`bio`（及其整棵子树）、`chrono`、`uuid`。
- 原 `gskits` 中被 asts 使用的模块移植为 crate 内部实现：
  - `src/gsbam/`：`BamRecordExt`（tag getter / identity / coverage / 比对边界）、
    `aligned_pairs_full`（等价于 htslib 的 `BamRecordExtensions::aligned_pairs_full`）、
    `compute_max_ins_of_each_ref_position`、`fasta2bam` / `fastq2bam`；
  - `src/fastx_reader.rs`：hand-rolled FASTA/FASTQ 解析器（替代 bio crate），
    命名与 qual 约定（ASCII−33）与原实现一致；
  - `src/dna.rs`、`src/pbar.rs`、`src/sys_monitor.rs`（改用 sysinfo 0.38）、
    `src/samtools.rs`（`samtools_bai` / `sort_by_coordinates`）；
  - `src/utils.rs` 新增 `ScopedTimer` / `Timer` / `command_line_str` /
    `phreq_list_2_quality`。
- 记录类型：`BamRecord = noodles_sam::alignment::RecordBuf`（自有可变记录），
  读端经 `reader.record_bufs(&header)` 产出，写端经 `write_alignment_record` 落盘。

### 兼容性验证（用户确认的验收标准）

以迁移前（git HEAD 构建）与迁移后二进制对同一份输入跑 asts 主流程，并对比输出：

| 输入路径 | 记录数 | `samtools view` 逐行 diff | bam_diff.py |
| --- | --- | --- | --- |
| BAM（test_data/sbr.bam + smc.bam） | 682（507 条含 soft clip） | **0 行差异**（SEQ/QUAL/CIGAR/FLAG/MAPQ/POS 及全部 aux tag） | 0 不一致 |
| FASTA smc（fastx2bam → 排序 → 主流程） | 224 | **0 行差异** | — |
| FASTQ smc | 同上数据 | **0 行差异** | — |

- 说明：比较对象是解码后的记录内容（`samtools view` 文本），不承诺 BGZF 压缩流逐字节相同
  （两边默认压缩级别/实现不同）。
- header 内容逐字段一致（`@HD VN:1.5 SO:unknown`、`@PG ID:asts PN:asts CL:... VN:...`、
  `@SQ` 按 tid 顺序）；唯一已知差异：noodles 固定按 `@HD → @SQ... → @PG` 顺序输出，
  迁移前为 `@HD → @PG → @SQ...`。记录内容与顺序不变。
- qual 兼容性：noodles 校验 qual ≤ 93，BAM 的缺失标记 `0xFF` 不允许直接写入；
  全 `0xFF` 的 qual（如 fasta2bam 占位值）改为留空，由 encoder 自动填 `0xFF×len`，
  与 rust-htslib 输出逐字节一致。

### 已知取舍与待办

- noodles 无 `set_threads`（BGZF 单线程解压）。asts 本身按线程切分任务并用
  crossbeam 并行，实际影响限于单个 reader/writer 的解压/压缩吞吐。
- **读端性能待基准验证**：读 BAM 改用 `reader.record_bufs(&header)`，每条记录会
  一次性解码（含全部 aux tag）；而 rust-htslib 的 `Record` 是惰性解码。smc/sbr
  输入只需少量字段，大数据集（13k-5h 级）下可能有可感知的读端开销，建议用
  `benches/benchmark.rs` 在真实数据上回归；若明显，可后续把 `ReadInfo` 提取
  泛型化到 `dyn noodles_sam::alignment::Record` 以避免整记录解码。
- `sys_monitor` 从 sysinfo 0.34 升至 0.38（API 兼容，逐字段核对过）。
- CI 已移除 `libhts-dev` / `libclang` 安装步骤；`rust-version = "1.91"`。
