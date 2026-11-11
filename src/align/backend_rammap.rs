//! rammap 后端（`rammap-core`，纯 Rust 的 minimap2 重实现）。
//!
//! 参数翻译刻意做成“显式全集”，不依赖库默认值 —— 历史实现（minimap2 系）在建内存
//! 索引时存在隐式行为，必须逐条复刻才能得到同样的结果（详见
//! `docs/migration-mm2-to-rammap.md`）：
//!
//! 1. 建内存索引时会把 `mid_occ` 硬设为 1000 且**不**做 occurrence 自动校准。
//!    rammap 的 `calibrate_mid_occ` 只在 `mid_occ == 0` 时生效，因此这里显式
//!    `.mid_occ(1000)`；漏掉这一步会让命中集合整体漂移。
//! 2. CIGAR 形态与 `is_supplementary` 语义由 `OutputConfig::{do_cigar, output_sam, eqx}`
//!    决定，缺 `eqx` 会退化成 `M`，而 identity 计算要求 `=`/`X`（op 7/8）。

use std::sync::Once;

use rammap::align::map::AlignFlags;
use rammap::align::pipeline::OutputConfig;
use rammap::api::{Aligner, MapOpts, Preset as RmPreset, Strand};
use rammap::{set_dp_cache_cap_mb, DEFAULT_DP_CACHE_CAP_MB};

use super::{AlignSetup, Hit, Preset, RefAligner};

/// 复刻历史实现的 `mapopt.mid_occ = 1000`（建内存索引时不做 occurrence 自动校准）。
const MID_OCC: usize = 1000;

/// DP scratch cache 上限：高线程数下限制峰值 RSS。
/// 可用环境变量 `ASTS_RAMMAP_DP_CACHE_CAP_MB` 覆盖，0 表示不限制。
const DP_CACHE_CAP_ENV: &str = "ASTS_RAMMAP_DP_CACHE_CAP_MB";

fn init_dp_cache_cap() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let cap = std::env::var(DP_CACHE_CAP_ENV)
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .unwrap_or(DEFAULT_DP_CACHE_CAP_MB);
        set_dp_cache_cap_mb(cap);
        tracing::info!(
            "rammap backend: dp_cache_cap_mb={}{}",
            cap,
            if cap == 0 { " (unbounded)" } else { "" }
        );
    });
}

fn preset_of(preset: Preset) -> RmPreset {
    match preset {
        Preset::MapOnt => RmPreset::MapOnt,
    }
}

struct RammapRefAligner {
    aligner: Aligner,
    target_name: String,
    target_len: usize,
    emit_cs: bool,
    emit_md: bool,
}

impl RefAligner for RammapRefAligner {
    fn map(&self, name: &str, seq: &[u8]) -> Vec<Hit> {
        let res = self.aligner.map_seq_with(
            name,
            seq,
            MapOpts {
                cs: Some(self.emit_cs),
                md: Some(self.emit_md),
            },
        );
        res.mappings.into_iter().map(hit_from_mapping).collect()
    }

    fn target_name(&self) -> &str {
        &self.target_name
    }

    fn target_len(&self) -> usize {
        self.target_len
    }
}

fn hit_from_mapping(m: rammap::Mapping) -> Hit {
    let cigar_ops = m
        .cigar_ops
        .as_ref()
        .expect("rammap backend: cigar_ops missing, OutputConfig::do_cigar must be true")
        .iter()
        .map(|op| (op.len, op.op))
        .collect::<Vec<_>>();

    Hit {
        target_name: m.target_name.to_string(),
        target_len: m.target_len,
        query_start: m.query_start,
        query_end: m.query_end,
        target_start: m.target_start,
        target_end: m.target_end,
        is_reverse: m.strand == Strand::Reverse,
        // rammap 的 is_primary = !is_secondary（附带的 supplementary 也算 primary），
        // 这里折叠成 SAM 语义：SAM primary = !secondary && sam_pri。
        is_primary: m.is_primary && !m.is_supplementary,
        is_supplementary: m.is_supplementary,
        mapq: m.mapq.clamp(0, 255) as u8,
        score: m.score,
        cigar_ops,
        cs: m.cs,
        md: m.md,
    }
}

pub fn build(setup: &AlignSetup, seq: &[u8], name: &str) -> Box<dyn RefAligner> {
    init_dp_cache_cap();

    let mut builder = Aligner::builder(preset_of(setup.preset))
        .hpc(setup.hpc)
        .cigar(true)
        // ★ 必须显式设置，否则 rammap 会按 mid_occ_frac 自动校准，结果与旧后端漂移
        .mid_occ(MID_OCC);
    if let Some(k) = setup.k {
        builder = builder.k(k);
    }
    if let Some(w) = setup.w {
        builder = builder.w(w);
    }

    let mut aligner = builder
        .from_seqs(vec![(name.to_string(), seq.to_vec())])
        .unwrap_or_else(|e| panic!("rammap backend: failed to build index for {}: {}", name, e));

    {
        let opt = aligner.options_mut();
        let scoring = &setup.scoring;
        if let Some(v) = scoring.match_score {
            opt.scoring.match_score = v;
        }
        if let Some(v) = scoring.mismatch_penalty {
            opt.scoring.mismatch_penalty = v;
        }
        if let Some(v) = scoring.gap_open {
            opt.scoring.gap_open = v;
        }
        if let Some(v) = scoring.gap_open2 {
            opt.scoring.gap_open2 = v;
        }
        if let Some(v) = scoring.gap_extend {
            opt.scoring.gap_extend = v;
        }
        if let Some(v) = scoring.gap_extend2 {
            opt.scoring.gap_extend2 = v;
        }
        if let Some(best_n) = setup.best_n {
            opt.filtering.best_n = best_n;
        }
        if setup.short_insert {
            opt.chaining.min_cnt = 2;
            opt.chaining.min_chain_score = 10;
            opt.alignment.min_dp_max = 10;
            // 语义等价于 minimap2 的 min_ksw_len，rammap 中改名为 min_dp_len
            opt.alignment.min_dp_len = 0;
        }
        opt.flags.insert(AlignFlags::EQX);
        if setup.sam_hit_only {
            opt.flags.insert(AlignFlags::SAM_HIT_ONLY);
        }
        if setup.secondary_seq {
            opt.flags.insert(AlignFlags::SECONDARY_SEQ);
        }
    }

    {
        let out = aligner.output_config_mut();
        *out = OutputConfig {
            do_cigar: true,
            do_cs: setup.emit_cs,
            cs_long: false,
            do_md: setup.emit_md,
            do_ds: false,
            // CIGAR 用 =/X 而非 M：identity() 依赖 op==7，缺了会直接 panic
            eqx: true,
            output_sam: true,
            rg_id: None,
            split_mode: false,
        };
    }

    Box::new(RammapRefAligner {
        aligner,
        target_name: name.to_string(),
        target_len: seq.len(),
        emit_cs: setup.emit_cs,
        emit_md: setup.emit_md,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::cigar;

    fn ont_setup() -> AlignSetup {
        AlignSetup {
            preset: Preset::MapOnt,
            k: Some(12),
            w: Some(5),
            hpc: true,
            best_n: Some(1),
            ..Default::default()
        }
    }

    #[test]
    fn test_build_and_map_basic() {
        let smc = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTGGGGGCCCCCAAAAA";
        let aln = build(&ont_setup(), smc, "smc");
        assert_eq!(aln.target_name(), "smc");
        assert_eq!(aln.target_len(), smc.len());

        let hits = aln.map("q1", smc);
        assert!(!hits.is_empty(), "perfect match should be aligned");
        let hit = &hits[0];
        assert_eq!(hit.target_name, "smc");
        assert!(hit.is_primary);
        assert!(!hit.is_supplementary);
        // eqx 生效：必须出现 = 或 X，而不是清一色 M
        assert!(
            hit.cigar_ops.iter().any(|&(_, op)| op == 7 || op == 8),
            "eqx required, got {:?}",
            hit.cigar_ops
        );
        assert!((cigar::identity(hit) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_cs_emitted_only_when_asked() {
        let smc = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTGGGGGCCCCCAAAAA";
        let mut setup = ont_setup();
        setup.emit_cs = true;
        let aln = build(&setup, smc, "smc");
        let hits = aln.map("q1", smc);
        assert!(hits[0].cs.is_some());

        let aln = build(&ont_setup(), smc, "smc");
        let hits = aln.map("q1", smc);
        assert!(hits[0].cs.is_none());
    }
}
