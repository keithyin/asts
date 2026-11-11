# 比对引擎迁移记录：minimap2（mm2）→ rammap

> 适用版本：asts **1.3.0** 起，比对引擎固定为 `rammap-core`，`minimap2`/`mm2` 后端与相关 feature 已移除。
> 本文是这次迁移的唯一权威说明：读它应当能完成「理解迁移 → 复现验证 → 升级 rammap」的全过程。

## 1. 摘要与结论

- 迁移目标：把 C 绑定（`minimap2-rs` + `mm2` 薄封装）换成纯 Rust 的 [rammap](https://github.com/jwanglab/rammap)（crate 名 `rammap-core`，**lib name `rammap`**），要求**比对结果保持一致**。
- 实测结论（真实 E. coli 数据，541,549 SMC / 8,303,236 subreads）：

| 指标 | 结果 |
| --- | --- |
| 输出记录数 | 8,089,046（两引擎相同） |
| qname 一致率 | **1.000000** |
| record 一致率（含 MAPQ） | **0.9999995**（8,089,042 / 8,089,046） |
| channel（SMC）命中集合 Jaccard | **1.000000**（541,538 channels） |
| 整份 SAM 逐行 diff | **仅 4 条记录**（3 条仅 MAPQ 不同，1 条选中了另一个等价比对） |
| Reporter 统计 | 两次运行逐项相同（含 fallback 数量） |
| 性能 | user CPU **-10.1%**，单 worker 吞吐 **+21.8%**，wall 持平（受 samtools sort 主导） |

- `test_data`（148 SMC / 702 subreads）上两引擎输出 `samtools view` **逐行完全一致**（682 条，含 `iy`/`ec` tag 与顺序）。

完整实测报告与产物：`/data1/eqv_run/REPORT.md`、`/data1/eqv_run/{full,logs,smoke}/`。

## 2. 背景与决策

迁移前的依赖链是 `minimap2-rs`（minimap2 C 源码绑定）→ `minimap2-sys`（FFI）→ `mm2`（作者自研薄封装，额外提供 `MappingExt`、`bam_writer`、`NoMemLeakAligner`）。问题：

1. **类型泄漏**：`minimap2::Mapping` 直接出现在 `mapping2record`、`cigar_adjust`、`sbr_and_ref_to_cs` 和各个 bin 里，换引擎等于全仓改签名。
2. **内存安全与可维护性**：C 绑定 + `NoMemLeakAligner`（一个空的 `Drop`，实际是绕坑写法）难以审计；`mm2` 的辅助代码与 asts 自身逻辑高度重复。
3. **可比性要求**：asts 的产物（BAM → MSA → identity/报告）对 CIGAR 与坐标极其敏感，迁移必须具备**可量化的一致性判据**，而不是"看起来差不多"。

关键决策：

1. **先抽 `Hit` 抽象，再换引擎**。`Hit`（`src/align/mod.rs`）是唯一的跨层契约，`AlignSetup` 描述引擎无关参数，`RefAligner` 是"单条参考序列的比对器"。这样换引擎只影响 `backend_rammap` 一个文件。
2. **参数显式全集**。所有原本由库默认值/隐式行为提供的参数，都逐条写进 `AlignSetup` 并显式下发，避免"默认值变了结果就变了"。
3. **迁移期双后端、迁移后删后端**。迁移期保留 `backend-minimap2` 用于对拍；1.3.0 起删除，改为**金标准快照测试**（`tests/alignment_golden.rs`）承担升级期回归职责。

## 3. 架构对比

迁移前：

```
bin/*  →  params::TOverrideAlignerParam  →  minimap2::Aligner<PresetSet> (C FFI)
                                                   ↓ minimap2::Mapping
                              mapping2record / cigar_adjust / MappingExt
```

迁移后：

```
bin/*  →  params::{AlignParams, IndexParams, OupParams}
             ↓ TOverrideAlignSetup
        align::AlignSetup  ──build_aligner()──▶  align::backend_rammap
                                                      ↓ rammap-core Aligner
        align::Hit  ◀─────────────────────────────────┘   （唯一跨层契约）
             ↓
        align::cigar (identity / coverage / record CIGAR)
             ↓
        mapping2record → bam_writer → 下游 MSA / 报告（零改动）
```

代码落点：

| 文件 | 职责 |
| --- | --- |
| `src/align/mod.rs` | `Hit` / `AlignSetup` / `Scoring` / `Preset` / `RefAligner` / `build_aligner` |
| `src/align/backend_rammap.rs` | 唯一的引擎适配层：建索引、参数翻译、`Mapping → Hit` |
| `src/align/cigar.rs` | record CIGAR（含反向坐标换算与 softclip）、`aligned_2_str`、identity/coverage |
| `src/align/cigar_adjust.rs` | 同聚物区 gap 左对齐（`--polyNGapLeftAlign`） |
| `src/bam_writer.rs` | BAM 写出（header/PG/`iy`/`ec`，顺序与历史一致） |
| `src/params.rs` | CLI 参数 → `AlignSetup` 的翻译 |

## 4. 参数与 API 映射表

| 迁移前（minimap2 语义） | rammap 对应 | 说明 |
| --- | --- | --- |
| `Aligner::builder().map_ont()` | `Aligner::builder(Preset::MapOnt)` | 预设同名 |
| `idxopt.set_hpc()` | `.hpc(true)` | 同聚物压缩索引；只在非短插入分支开启 |
| `idxopt.k / w` | `.k(..) / .w(..)` | 建索引参数只能走 builder，`MapOptions` 里没有 |
| `with_seq_and_id(seq, id)` | `builder(..).from_seqs(vec![(id, seq)])` | 内存索引，每次每 SMC 建一次（既有行为） |
| **隐式 `mapopt.mid_occ = 1000`** | **必须显式 `.mid_occ(1000)`** | 见"坑 1"，最易漏 |
| `mapopt.best_n` | `options_mut().filtering.best_n` | subread→SMC 固定 1 |
| `mapopt.min_cnt` | `chaining.min_cnt` | 短插入：2 |
| `mapopt.min_chain_score` | `chaining.min_chain_score` | 短插入：10 |
| `mapopt.min_dp_max` | `alignment.min_dp_max` | 短插入：10 |
| `mapopt.min_ksw_len` | `alignment.min_dp_len` | **改名** |
| `mapopt.a / b` | `scoring.match_score / mismatch_penalty` | |
| `mapopt.q / q2` | `scoring.gap_open / gap_open2` | `-o 4,48` 拆成两个字段 |
| `mapopt.e / e2` | `scoring.gap_extend / gap_extend2` | `-e 2,1` 拆成两个字段 |
| `with_cigar()` → `MM_F_CIGAR \| MM_F_OUT_CS` | `OutputConfig{do_cigar: true, do_cs: ..}` | |
| `with_sam_out()` → `MM_F_OUT_SAM` | `OutputConfig{output_sam: true}` | 影响 supplementary/sam_pri 语义，必须开 |
| `MM_F_EQX (0x4000000)` | `OutputConfig{eqx: true}` + `AlignFlags::EQX` | CIGAR 用 `=`/`X` |
| `MM_F_SAM_HIT_ONLY (0x40000000)` | `AlignFlags::SAM_HIT_ONLY` | 语义：只输出命中 |
| `MM_F_SECONDARY_SEQ (0x1000000000)` | `AlignFlags::SECONDARY_SEQ` | cs→ref 路径需要 |
| `map(seq, cs, md, ..)` 的 cs/md 开关 | `map_seq_with(name, seq, MapOpts{cs, md})` | per-call 覆盖 |
| `hit.alignment.cigar: Vec<(u32,u8)>` | `hit.cigar_ops: Vec<CigarOp{len, op}>` | op 编码一致（0..5, 7==, 8=X） |
| `hit.alignment.cs / md` | `hit.cs / hit.md` | |
| `hit.alignment.alignment_score` | `hit.score` | |
| `hit.query_len: Option<NonZeroI32>` | **不存在** | coverage 需外部传 `seq.len()` |
| `hit.strand: Strand` | `hit.strand: Strand` | 同名枚举 |
| `is_primary` | `is_primary && !is_supplementary` | 语义不同，见"坑 3" |
| `NoMemLeakAligner` | 无（纯 Rust 不需要） | |
| `mm2::convert_mapping_cigar_to_record_cigar` | `align::cigar::ops_to_record_cigar` | 逻辑内联 |
| `mm2::bam_writer::write_bam_worker` | `src/bam_writer.rs` | 逻辑内联 |

## 5. 必踩的坑（按踩坑概率排序）

### 坑 1：`mid_occ` 会被自动校准

迁移前的内存索引路径**手工把 `mid_occ` 设为 1000 且不调用 `mm_mapopt_update`**，即不做 occurrence 自动校准。rammap 的建库流程在 `finish()` 里会调 `calibrate_mid_occ`，它只在 `seeding.mid_occ == 0` 时生效。

> 不显式 `.mid_occ(1000)` → `mid_occ` 被按 `mid_occ_frac` 校准成别的值 → 高频 k-mer 过滤强度变化 → **命中集合整体漂移**，下游 MSA 大面积不同。

### 坑 2：CIGAR 必须是 `=`/`X`

`cigar::identity()` 依赖 op==7（`assert!(op != 0, "eqx required")`）。只设 `do_cigar` 而不设 `eqx` 时 CIGAR 会退化成 `M`，表现为 panic 或 identity 全 0 → 过滤全部失效。

### 坑 3：`is_primary` 语义不同

- rammap：`is_primary = !is_secondary`（**supplementary 也算 primary**），`is_supplementary = !is_secondary && !sam_pri`。
- 历史实现：`is_primary` 只对 SAM primary 为真。

适配层统一折叠为 SAM 语义（`is_primary = m.is_primary && !m.is_supplementary`），否则 `align_cs_to_ref` 里的 `if !hit.is_primary { continue }` 会多处理 supplementary 记录，改变"参考区间是否有错"的判定。

### 坑 4：`Mapping` 没有 `query_len`

`MapResult` 只有 `mappings` 一个字段，coverage 必须由调用方传原始 query 长度（asts 里恒等于整条 subread 长度）。

### 坑 5：输出配置是 aligner 级，不是 per-call

`OutputConfig`（`do_cigar/do_cs/cs_long/do_md/do_ds/eqx/output_sam/rg_id/split_mode`）与 `AlignFlags` 都挂在 aligner 上。asts 有两个用途、flag 组合不同：

- subread→SMC：`EQX + SAM_HIT_ONLY`，需要 `cs`/`MD`；
- cs→ref：`EQX + SECONDARY_SEQ`，不需要 `cs`/`MD`。

因此**必须各自建 aligner**（现状本来也是每个 SMC / 每个 ref 新建一次，不要图省事复用）。

### 坑 6：flag 位值两边不同，不能抄裸数字

`EQX` 在 minimap2 是 `0x4000000`，在 rammap 是 `0x8000`；`SECONDARY_SEQ` 分别是 `0x1000000000` / `0x800000000`。**一律用命名常量**。

### 坑 7：rayon 嵌套并行 + DP cache 峰值

`rammap-core` 默认开启 `parallel`(rayon)。asts 自身已按 `threads-4` 切 worker，嵌套并行会线程超订，因此依赖声明里 `default-features = false`。

高线程下每线程 DP scratch cache 会推高峰值 RSS，可用环境变量兜底：`ASTS_RAMMAP_DP_CACHE_CAP_MB`（0 表示不限制，默认取 rammap 的 `DEFAULT_DP_CACHE_CAP_MB`）。

### 坑 8：逐行 diff 之前先确认顺序确定性

asts 输出会 `samtools sort`，但同一坐标上的多条记录顺序取决于 worker 完成顺序。**判等应以"记录多重集"为准**（`scripts/compare_asts_runs.py` 按 qname 分组比较），需要"逐字节"结论时再单独确认 `qname` 序列 md5 是否相同。

## 6. 验证方法论（三层）

| 层 | 手段 | 耗时 | 触发时机 |
| --- | --- | --- | --- |
| 进程内 | `cargo test --test alignment_golden`：固定用例的 `Hit` 字段与 fixture 逐行比对 | 毫秒级 | 每次改引擎/参数 |
| 切片冒烟 | 从真实数据切 1–2 千个 channel，两版各跑一次后 `diff <(samtools view a) <(samtools view b)` | 秒级 | 换引擎/升大版本 |
| 全量 | 全量跑批 + `scripts/compare_asts_runs.py` + 整份 SAM diff | 分钟级 | 发布前 |

判据（沿用本次迁移验收标准）：

- **通过**：record 一致率 = 1.0 且 channel Jaccard = 1.0；
- **有条件通过**：一致率 ≥ 0.999，且全部差异可归因为 (a) MAPQ 数值、(b) 打平并列的 tie-breaking；
- **不通过**：低于 0.999 或差异无法归因。

差异归因三步法：先按 qname 求差集 → 再看是"CIGAR/坐标不同"还是"仅 MAPQ 不同" → 最后把该 channel 单独切出来复跑，确认是否落在低复杂度/重复区（tie-breaking 的高发区）。

## 7. 实测数据（本次迁移）

环境：256 core / 503 GB RAM / samtools 1.19.2 / rustc 1.92.0；两次运行参数完全一致（`--threads 64 --remove-sorted-bam`）。

### 7.1 test_data

输出 682 条记录，`samtools view` **逐行完全一致**（含 `iy`/`ec` tag）。

### 7.2 E. coli 全量（8,303,236 subreads / 541,549 SMC）

- 输出 8,089,046 条记录，541,538 个 channel 有结果（3 个 channel 无可用比对，两引擎一致）。
- Reporter 逐项相同：`sbr inp 8107310 / filter_by_length 38062 / filter_by_alignment 18264 / fallback 29293 / rescued 11029`。
- 差异 4 条：

| subread | flag | pos | MAPQ | CIGAR |
| --- | --- | --- | --- | --- |
| `read_452949/452949/subread/2` | 16 | 1118 / 839 | 4 / 1 | 不同（低复杂度区并列） |
| `read_536511/536511/subread/5` | 0 | 4 / 4 | 25 / 49 | 相同 |
| `read_559521/559521/subread/9` | 0 | 1 / 1 | 10 / 1 | 相同 |
| `read_879800/879800/subread/3` | 0 | 1 / 1 | 44 / 60 | 相同 |

  下游影响：MAPQ 不参与 MSA / identity / 报告 → 3 条零影响；1 条影响 1 个 SMC（1/541,538）。

### 7.3 性能

| 指标 | minimap2 | rammap |
| --- | --- | --- |
| wall | 326.66 s | 326.22 s |
| user CPU | 8,899.16 s | 8,002.21 s |
| sys CPU | 380.34 s | 339.27 s |
| 子进程峰值 RSS | 24.9 GB | 24.9 GB |
| 单 worker 比对吞吐 | 105.1 /s | 128.0 /s |

wall 持平是因为整条流程被 `samtools sort`（两份输入 `sort_by_tag` + 结果排序）与 I/O 主导；比对段本身更快。

### 7.4 收尾验证（1.3.0，删除 minimap2 后端之后）

删除旧后端后，用 1.3.0 在**同一份真实数据**上重跑，与 1.2.0（双后端版本）的 rammap 输出逐行对比：

```
samtools view rammap.bam | diff - <(samtools view rammap_v130.bam)
→ 0 行差异（8,089,046 条记录完全一致）
```

结论：删除 `minimap2` 依赖与 `backend_minimap2` 模块**没有改变 rammap 路径的任何行为**。
本次运行指标：wall 329.9 s、user 7,849.18 s、sys 284.56 s、峰值 RSS 24.9 GB（与基线同量级）。
产物与日志：`/data1/eqv_run/full/rammap_v130.bam`、`/data1/eqv_run/logs/full_rammap_v130.{log,json}`、
`/data1/eqv_run/logs/v130_vs_v120.diff`（空文件）。

## 8. 日常使用与运维

```bash
# 构建
cargo build --release

# 主流程
asts -q subreads.bam -t smc.bam -p out/prefix --threads 64

# 回归测试（改引擎/参数后必跑）
cargo test                       # 含金标准快照
cargo test --test alignment_golden -- --nocapture

# 比较两次运行（升级前后、参数调整前后）
scripts/compare_asts_runs.py --a old.bam --b new.bam --threshold 0.999
```

- `--threads` 必须 ≥ 10；实际 align worker = `threads - 4`。
- 每个 SMC 都会新建一次内存索引，线程数不是越高越好；建议 `threads` ≤ 核数的 1/4。
- 若峰值 RSS 偏高：`ASTS_RAMMAP_DP_CACHE_CAP_MB=64 asts ...`。

常见错误对照：

| 现象 | 原因 |
| --- | --- |
| `rammap backend: cigar_ops missing` | 建库时 `OutputConfig.do_cigar` 没开 |
| identity 全 0 / `eqx required` panic | 没开 `eqx` |
| 命中集合整体不同 | 漏了 `.mid_occ(1000)` |
| 内存暴涨 | `rammap-core` 的 `parallel` feature 没关，或线程数过高 |

## 9. rammap 升级 SOP

1. 改 `Cargo.toml` 里的精确版本（`rammap-core = "=X.Y.Z"`），`cargo build` 通过后先看 changelog/API 变动。
2. `cargo test`：金标准快照测试若失败，**逐条 review diff**（`UPDATE_GOLDEN=1` 会重写 fixture，先看 `git diff` 再决定是否接受）。重点关注 `mapq` 与 `cigar`。
3. 切片冒烟：从真实数据切 ~2000 channel，跑两版后 `diff <(samtools view a) <(samtools view b)`。
4. 全量：跑两版，`scripts/compare_asts_runs.py` + 整份 SAM diff，按第 6 节判据判定。
5. 若差异超出判据：优先检查第 5 节的坑（尤其 `mid_occ`、`eqx`、`is_primary`），再联系上游。
6. 结论写进本文件的"变更记录"小节，并把新的实测报告路径记录下来。

## 10. 历史与追溯

- 1.2.0：引入 `Hit` 抽象 + 双后端（`backend-rammap` / `backend-minimap2`），完成等价性验证。
- 1.3.0：删除 `minimap2` 依赖与 `backend-minimap2` feature，收敛为单引擎；双后端对拍测试改为金标准快照测试。
- 依赖变化：`minimap2`、`minimap2-sys`、`mm2` 及传递依赖（`hp_tr_finder`、`glob` 等）已从 `Cargo.lock` 移除。
- 追溯入口：git 历史（1.2.0 tag 保留旧后端）、`/data1/eqv_run/REPORT.md`（本次全量验证的完整记录与 BAM 产物）。

## 11. 如何再换一次引擎（如果将来有必要）

`Hit` / `AlignSetup` / `RefAligner` 这层是专门为此保留的：

1. 新增 `src/align/backend_xxx.rs`，实现 `build(setup, seq, name) -> Box<dyn RefAligner>`，把引擎结果翻译成 `Hit`；
2. `src/align/mod.rs::build_aligner` 改为调用它（或加一个临时 feature 做双实现对拍）；
3. 用第 6 节的判据验证；`mapping2record` 及以上**不需要任何改动**。
