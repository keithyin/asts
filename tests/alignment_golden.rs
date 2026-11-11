//! 金标准快照回归测试。
//!
//! asts 只保留 rammap 一个比对引擎后，已经没有「另一个后端」可以对拍；这个测试把
//! **当前已验证过的比对结果**固化成 fixture，用于在未来出现下列变更时立刻报警：
//!
//! - 升级 `rammap-core` 版本（该 crate 月更，且 1.1 → 1.2 有 API/行为变动）
//! - 误改 `AlignSetup` 的任何参数（k/w/hpc、打分、短插入门槛、flag 组合）
//! - 改动 CIGAR 规范化、identity/coverage 或 `Hit` 转换逻辑
//!
//! 断言粒度是 `Hit` 的关键字段：`(qstart, qend, tstart, tend, strand, primary,
//! supplementary, mapq, cigar_ops, cs)`，即所有会影响下游 BAM / MSA 的字段。
//!
//! ```bash
//! # 常规回归
//! cargo test --test alignment_golden -- --nocapture
//!
//! # 人工确认差异可接受后，重新生成基线（会改写 tests/fixtures/*.tsv）
//! UPDATE_GOLDEN=1 cargo test --test alignment_golden
//! ```
//!
//! 注意：`mapq` 是升级后最可能发生偏移的字段（见
//! `docs/migration-mm2-to-rammap.md`）；每次重建基线都必须人工抽查 diff。

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use asts::align::{build_aligner, cigar, AlignSetup, Hit, Preset};
use asts::mapping2record::build_bam_record_from_mapping;
use asts::read_info::ReadInfo;
use gskits::fastx_reader::{fasta_reader::FastaFileReader, read_fastx};

const FIXTURE_DIR: &str = "tests/fixtures";

fn hit_key(hit: &Hit) -> String {
    format!(
        "{name}:{qs}-{qe}|{ts}-{te}|rev={rev}|pri={pri}|sup={sup}|mapq={mapq}|\
         score={score}|cigar={cigar:?}|cs={cs:?}",
        name = hit.target_name,
        qs = hit.query_start,
        qe = hit.query_end,
        ts = hit.target_start,
        te = hit.target_end,
        rev = hit.is_reverse,
        pri = hit.is_primary,
        sup = hit.is_supplementary,
        mapq = hit.mapq,
        score = hit.score,
        cigar = hit.cigar_ops,
        cs = hit.cs,
    )
}

fn fixture_path(case: &str) -> PathBuf {
    Path::new(FIXTURE_DIR).join(format!("align_golden_{case}.tsv"))
}

fn update_golden() -> bool {
    std::env::var("UPDATE_GOLDEN")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false)
}

/// 跑一个用例，并与 fixture 逐行比对（`UPDATE_GOLDEN=1` 时改写 fixture）。
fn check_case(
    case: &str,
    setup: &AlignSetup,
    target_name: &str,
    target_seq: &str,
    queries: &[(String, String)],
) {
    let aligner = build_aligner(setup, target_seq.as_bytes(), target_name);
    assert_eq!(aligner.target_name(), target_name);
    assert_eq!(aligner.target_len(), target_seq.len());

    let mut target_idx = HashMap::new();
    target_idx.insert(target_name.to_string(), (0, target_seq.len()));

    let mut actual: Vec<String> = Vec::new();
    let mut n_hit_records = 0usize;

    for (qname, qseq) in queries {
        if qseq.is_empty() {
            continue;
        }
        let hits = aligner.map(qname, qseq.as_bytes());
        if hits.is_empty() {
            actual.push(format!("{case}\t{qname}\tNO_HIT"));
            continue;
        }
        for (rank, hit) in hits.iter().enumerate() {
            // eqx 是 identity 计算的前提：出现纯 M 说明 OutputConfig.eqx 失效
            assert!(
                hit.cigar_ops.iter().any(|&(_, op)| op == 7 || op == 8),
                "{case}/{qname}: cigar 缺少 =/X（eqx 失效），got {:?}",
                hit.cigar_ops
            );
            if hit.is_primary && !hit.is_supplementary {
                let query_record = ReadInfo::new_fa_record(qname.clone(), qseq.clone());
                let record = build_bam_record_from_mapping(hit, &query_record, &target_idx);
                assert_eq!(record.pos(), hit.target_start as i64);
                n_hit_records += 1;
            }
            actual.push(format!("{case}\t{qname}\t{rank}\t{}", hit_key(hit)));
        }
    }

    let path = fixture_path(case);

    if update_golden() {
        fs::create_dir_all(FIXTURE_DIR).unwrap();
        fs::write(&path, actual.join("\n") + "\n").unwrap();
        println!(
            "[UPDATE_GOLDEN] wrote {} lines to {}",
            actual.len(),
            path.display()
        );
        return;
    }

    let expected = fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "fixture {} 不存在或不可读（{}）；先运行 `UPDATE_GOLDEN=1 cargo test --test alignment_golden` 生成并人工核对",
            path.display(),
            e
        )
    });
    let expected_lines = expected.lines().collect::<Vec<_>>();

    assert_eq!(
        actual.len(),
        expected_lines.len(),
        "{case}: hit 行数变化 {} -> {}（若为刻意变更，请 UPDATE_GOLDEN=1 重建并人工核对）",
        expected_lines.len(),
        actual.len()
    );
    for (idx, (got, want)) in actual.iter().zip(expected_lines.iter()).enumerate() {
        assert_eq!(
            got, want,
            "{case}: 第 {} 行与基线不一致（若为刻意变更，请 UPDATE_GOLDEN=1 重建并人工核对）",
            idx + 1
        );
    }

    println!(
        "[golden] {case}: {} 行一致（其中 {} 条 primary 记录成功装配为 BAM）",
        actual.len(),
        n_hit_records
    );
}

fn read_fa(path: &str) -> Option<Vec<(String, String)>> {
    if !Path::new(path).exists() {
        eprintln!("skip: {} not found", path);
        return None;
    }
    let reader = FastaFileReader::new(path.to_string());
    Some(
        read_fastx(reader)
            .into_iter()
            .map(|r| (r.name, r.seq))
            .collect(),
    )
}

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

fn random_seq(len: usize, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed;
    (0..len)
        .map(|_| bases[(lcg(&mut state) % 4) as usize] as char)
        .collect()
}

/// 按固定速率注入替换/插入/删除，构造与参考有真实分歧的 query（确定性，无随机依赖）。
fn mutate(seq: &str, rate: f64, seed: u64) -> String {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed;
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.as_bytes() {
        let r = (lcg(&mut state) % 10000) as f64 / 10000.0;
        if r < rate / 3.0 {
            out.push(bases[(lcg(&mut state) % 4) as usize]); // 替换
        } else if r < 2.0 * rate / 3.0 {
            out.push(bases[(lcg(&mut state) % 4) as usize]); // 插入
            out.push(b);
        } else if r < rate {
            // 删除
        } else {
            out.push(b);
        }
    }
    String::from_utf8(out).unwrap()
}

fn ont_setup(short_insert: bool) -> AlignSetup {
    AlignSetup {
        preset: Preset::MapOnt,
        k: Some(if short_insert { 4 } else { 12 }),
        w: Some(if short_insert { 1 } else { 5 }),
        hpc: !short_insert,
        best_n: Some(1),
        short_insert,
        emit_cs: true,
        emit_md: true,
        sam_hit_only: true,
        secondary_seq: false,
        ..Default::default()
    }
}

/// 短插入分支：<200bp 的 SMC，k=4/w=1、无 hpc、min_cnt=2 / min_chain_score=10 /
/// min_dp_max=10 / min_dp_len=0。
#[test]
fn golden_short_insert() {
    let ref_seq = "CCCCTTTTAAAAGGGGGAGAGACCCCTTTTGGGGGAAAAACCCCCGGGGGTTTTTAAA";
    let queries = vec![
        (
            "exact".to_string(),
            "CCCCTTTTAAAAGGGGGAGAGACCCCTTTTGGGGGAAAAACCCCCGGGGGTTTTTAAA".to_string(),
        ),
        (
            "ins".to_string(),
            "CCCCTTTTGAAAAGGGGGAGAGACCCCTTTTGGGGGAAAAACCCCCGGGGGTTTTTAAA".to_string(),
        ),
        (
            "del".to_string(),
            "CCCCTTTTAAAGGGGGAGAGACCCCTTTTGGGGGAAAAACCCCCGGGGGTTTTTAAA".to_string(),
        ),
        (
            "homopolymer".to_string(),
            "CCCCTTTTAAAAGGGGGGAGAGACCCCTTTTGGGGGAAAAACCCCCGGGGGTTTTTAAA".to_string(),
        ),
        (
            "truncated".to_string(),
            "CCCCTTTTAAAAGGGGGAGAGACCCCTTTTGGGGGAAAA".to_string(),
        ),
    ];
    check_case("short_insert", &ont_setup(true), "ref", ref_seq, &queries);
}

/// 反向链：query 是参考的反向互补，考察 strand / 坐标 / softclip 换算。
#[test]
fn golden_reverse_strand() {
    let ref_seq = "ACGTTGCAAGTCCTAGGCTAACCGGTTACGATCGTAGCTAGCATCGATCGATCGTAGC";
    let rc = String::from_utf8(gskits::dna::reverse_complement(ref_seq.as_bytes())).unwrap();
    let queries = vec![
        ("rc".to_string(), rc),
        (
            "fwd_partial".to_string(),
            "GTCCTAGGCTAACCGGTTACGATCGTAGCTAGCATCGATCG".to_string(),
        ),
    ];
    check_case("reverse", &ont_setup(false), "ref", ref_seq, &queries);
}

/// 合成 panel：3kb 参考 + 8 段各 5% 分歧的 query（正/反链/截断），考察多命中、
/// 部分覆盖与 fallback 之外的主流程路径。
#[test]
fn golden_synthetic_panel() {
    let ref_seq = random_seq(3000, 42);
    let mut queries: Vec<(String, String)> = Vec::new();
    for i in 0..8u64 {
        let start = (i as usize * 300).min(ref_seq.len() - 600);
        let sub = &ref_seq[start..start + 600];
        let q = mutate(sub, 0.05, 1000 + i);
        queries.push((format!("fwd{i}"), q.clone()));
        queries.push((
            format!("rev{i}"),
            String::from_utf8(gskits::dna::reverse_complement(q.as_bytes())).unwrap(),
        ));
        queries.push((format!("trunc{i}"), q[80..q.len() - 80].to_string()));
    }
    check_case("synthetic", &ont_setup(false), "ref", &ref_seq, &queries);
}

/// 仓库自带 test_data（缺失时跳过），覆盖真实长读长的完整流程。
#[test]
fn golden_test_data() {
    let ref_fa = match read_fa("test_data/ref.fa") {
        Some(v) => v,
        None => return,
    };
    let query_fa = match read_fa("test_data/query.fa") {
        Some(v) => v,
        None => return,
    };
    let (ref_name, ref_seq) = &ref_fa[0];
    check_case("test_data", &ont_setup(false), ref_name, ref_seq, &query_fa);
}

/// 与上一条用例互补：identity/coverage 辅助函数的数值不随引擎漂移。
#[test]
fn golden_identity_helpers() {
    let ref_seq = "ACGTTGCAAGTCCTAGGCTAACCGGTTACGATCGTAGCTAGCATCGATCGATCGTAGC";
    let aligner = build_aligner(&ont_setup(false), ref_seq.as_bytes(), "ref");
    let hits = aligner.map("self", ref_seq.as_bytes());
    assert!(!hits.is_empty());
    let hit = &hits[0];
    assert!((cigar::identity(hit) - 1.0).abs() < 1e-6);
    assert!((cigar::identity_gap_compressed(hit) - 1.0).abs() < 1e-6);
    assert!((cigar::query_coverage(hit, ref_seq.len()) - 1.0).abs() < 1e-6);
    assert!((cigar::target_coverage(hit) - 1.0).abs() < 1e-6);
}
