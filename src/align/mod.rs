//! 比对抽象层（单一引擎：rammap）。
//!
//! asts 早期直接依赖 minimap2 的 C 绑定（`minimap2-rs` + `mm2` 薄封装），引擎返回的
//! `Mapping` 类型泄漏到了 `mapping2record`、`cigar_adjust`、`sbr_and_ref_to_cs` 以及各个
//! bin 里。现在引擎固定为纯 Rust 的 [rammap](https://github.com/jwanglab/rammap)
//! （crate `rammap-core`，lib name `rammap`），并通过本模块的 [`Hit`] / [`AlignSetup`] /
//! [`RefAligner`] 把引擎隔离在一层之后：
//!
//! - [`AlignSetup`]：引擎无关的参数描述（k/w/hpc、打分、阈值、输出开关），由
//!   `params::AlignParams` / `IndexParams` 翻译而来；
//! - [`RefAligner`]：持有单条参考序列（SMC 或 ref）索引的比对器；
//! - [`Hit`]：唯一的跨层契约，`mapping2record` 及以上只认它。
//!
//! 关键契约：**下游结果只取决于 [`Hit`] 转换出来的 BAM record**，MSA / identity /
//! 报告都从 BAM record 派生。因此参数翻译的等价性只需在 `Hit` 这一层验证。
//!
//! 迁移背景、参数映射表与踩坑清单见 `docs/migration-mm2-to-rammap.md`。

pub mod backend_rammap;
pub mod cigar;
pub mod cigar_adjust;

/// 单条比对结果（引擎无关）。
///
/// 字段语义与 SAM/BAM 对齐：
/// - `query_start/query_end` 是**原始（正向）query** 坐标，反向比对的 query 也按正向给出；
/// - `is_primary` 表示 SAM primary（FLAG 0x0 的那条），`is_supplementary` 表示 SAM supplementary；
///   两者不会同时为 true，secondary 记录两者都为 false；
/// - `cigar_ops` 必须是 `=`/`X`（eqx）形式，op 编码与 SAM/BAM 一致：
///   `0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 7==, 8=X`。
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Hit {
    pub target_name: String,
    pub target_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub is_reverse: bool,
    pub is_primary: bool,
    pub is_supplementary: bool,
    pub mapq: u8,
    pub score: i32,
    pub cigar_ops: Vec<(u32, u8)>,
    pub cs: Option<String>,
    pub md: Option<String>,
}

/// 打分参数（对应 minimap2 的 `-A/-B/-O/-E/-O2/-E2`）。
///
/// `None` 表示保持 preset 默认值，不覆盖。一阶/二阶 gap 罚分拆成独立字段，
/// 是为了保留 CLI 语义：只给一个值（如 `-o 4`）时只改一阶，不改二阶。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: Option<i32>,
    pub mismatch_penalty: Option<i32>,
    /// `-O` 第一个值。
    pub gap_open: Option<i32>,
    /// `-O` 第二个值。
    pub gap_open2: Option<i32>,
    /// `-E` 第一个值。
    pub gap_extend: Option<i32>,
    /// `-E` 第二个值。
    pub gap_extend2: Option<i32>,
}

/// asts 用到的最小 preset 集合。目前所有入口都使用 `map-ont`。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Preset {
    #[default]
    MapOnt,
}

/// 与引擎无关的建库 + 比对参数。
///
/// 每个字段都对应一个必须显式设置的 minimap2 语义项，理由见
/// `docs/migration-mm2-to-rammap.md` 的「参数与 API 映射表」。
#[derive(Debug, Clone, Default)]
pub struct AlignSetup {
    pub preset: Preset,
    /// `None` 表示使用 preset 默认 k。
    pub k: Option<usize>,
    /// `None` 表示使用 preset 默认 w。
    pub w: Option<usize>,
    /// 同聚物压缩索引（minimap2 `-H`）。
    pub hpc: bool,
    pub scoring: Scoring,
    /// `None` 表示使用 preset 默认 `best_n`。
    pub best_n: Option<i32>,
    /// 短插入（<200bp SMC）专用门槛：min_cnt=2 / min_chain_score=10 /
    /// min_dp_max=10 / min_dp_len=0。
    pub short_insert: bool,
    /// 是否输出 `cs` tag。
    pub emit_cs: bool,
    /// 是否输出 `MD` tag。
    pub emit_md: bool,
    /// 只输出命中（minimap2 `MM_F_SAM_HIT_ONLY`）。
    pub sam_hit_only: bool,
    /// secondary 也输出硬剪切后的 SEQ（minimap2 `MM_F_SECONDARY_SEQ`）。
    pub secondary_seq: bool,
}

/// 一个持有单条参考序列（SMC 或 ref）索引的比对器。
///
/// 必须 `Send + Sync`：asts 把同一个比对器共享给多个 worker 线程。
pub trait RefAligner: Send + Sync {
    /// 比对一条 query。`name` 仅用于结果里的 query 名，不参与比对决策。
    fn map(&self, name: &str, seq: &[u8]) -> Vec<Hit>;

    /// 索引中唯一参考序列的名字（即 `Hit::target_name` 的取值）。
    fn target_name(&self) -> &str;

    /// 索引中唯一参考序列的长度。
    fn target_len(&self) -> usize;
}

/// 为一条参考序列建立索引并返回比对器。
pub fn build_aligner(setup: &AlignSetup, seq: &[u8], name: &str) -> Box<dyn RefAligner> {
    backend_rammap::build(setup, seq, name)
}
