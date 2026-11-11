```
# sbr2smc.aligned.bam and it's corresponding .bam file will be generated
asts -q subreads.bam -t smc.bam -p sbr2smc.aligned
```

# 比对引擎

asts 的比对引擎固定为纯 Rust 的 [rammap](https://github.com/jwanglab/rammap)
（crate `rammap-core = 1.2.0`，minimap2 的重实现）。minimap2 C 绑定与 `mm2` 封装已在
1.3.0 移除，无 feature 开关。

```bash
cargo build --release
cargo test                                  # 含金标准快照回归
cargo test --test alignment_golden -- --nocapture
```

参数统一收敛到与引擎无关的 `asts::align::AlignSetup`，再由 `backend_rammap` 翻译：
`k/w/hpc`、`best_n`、`mid_occ = 1000`（不做 occurrence 自动校准）、打分
`-A/-B/-O o1,o2/-E e1,e2`、短插入（SMC < 200bp）门槛
`min_cnt=2 / min_chain_score=10 / min_dp_max=10 / min_dp_len=0`，以及输出开关
`eqx`（CIGAR 必须用 `=`/`X`）、`cs`/`MD`、`SAM_HIT_ONLY` / `SECONDARY_SEQ`。

高线程下可用 `ASTS_RAMMAP_DP_CACHE_CAP_MB` 限制每线程 DP scratch cache（0 表示不限制）。

> 迁移背景、参数映射表、踩坑清单与 rammap 升级流程见
> [`docs/migration-mm2-to-rammap.md`](docs/migration-mm2-to-rammap.md)。

## 结果一致性验收

一致性结论已固化为两层护栏：

```bash
# 1) 进程内金标准快照（毫秒级，改引擎/参数后必跑）
cargo test --test alignment_golden

# 2) 整条流程对比：同一份输入跑两次，比较两份 BAM
asts -q test_data/sbr.bam -t test_data/smc.bam -p /tmp/run_a --threads 10
asts -q test_data/sbr.bam -t test_data/smc.bam -p /tmp/run_b --threads 10
scripts/compare_asts_runs.py --a /tmp/run_a.bam --b /tmp/run_b.bam
```

`compare_asts_runs.py` 对比 `flag / rname / pos / mapq / CIGAR / cs`，输出 qname 一致率、
record 一致率、每个 channel 的命中集合 Jaccard；任一项低于 `--threshold`（默认 0.999）
则退出码非 0。

历史实测（minimap2 → rammap 迁移验证）：`test_data` 上两引擎 `samtools view` 输出逐行完全
一致；真实 E. coli 全量（8,303,236 subreads / 541,549 SMC）输出 8,089,046 条记录中仅 4 条
不同（3 条仅 MAPQ、1 条并列 tie-breaking），user CPU -10.1%。详见迁移文档。

# chagelog

## 1.3.0

* 移除 `minimap2` 依赖与 `backend-minimap2` 特性，比对引擎收敛为单一 rammap
* 删除 `src/align/backend_minimap2.rs`，简化 `align/mod.rs` 的 cfg 门控
* `tests/backend_equivalence.rs` → `tests/alignment_golden.rs`（金标准快照，防 rammap 升级回归）
* `scripts/compare_backends.py` → `scripts/compare_asts_runs.py`（比较任意两次 asts 运行）
* 新增迁移文档 `docs/migration-mm2-to-rammap.md`

## 1.2.0

* 比对后端从 `minimap2-rs` + `mm2` 迁移到纯 Rust 的 `rammap-core`（默认特性 `backend-rammap`）
* 新增与引擎无关的 `asts::align::{Hit, AlignSetup, RefAligner}` 抽象层，
  移除对 `mm2` 的依赖（`bam_writer`、cigar 转换、identity/coverage 计算内联）
* 保留 `backend-minimap2` 特性用于回退与等价性对比
* 新增 `tests/backend_equivalence.rs` 与 `scripts/compare_backends.py`

## 0.9.2

* Don't execute fallback alignment when the first stage alignment get enough aligned result.

## 0.9.1

* reporter bug fix

## 0.9.0

* reporter
* gsmm2(0.20 -> 0.22)
* identity_gap_compressed instead of identity_without_long_indel(10)

## 0.8.0

* using hpc
* default kmer/wins and fallback kmer/wins


## 0.7.1

* length filter when add sbr to SubreadsAndSmc
* log file

## 0.4.0

* dw,cr,ar and qual from sbr.bam will dump the result bam
* worked for smc reads that length < 200 && >200

## 0.3.0

* add ch tag to bam record

## 0.2.1

* worked

## 0.1.1

## 0.1.0
