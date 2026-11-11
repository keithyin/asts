//! CIGAR / identity / coverage 工具。
//!
//! 这些函数产出 BAM 的 CIGAR 与下游过滤用的 identity/coverage。计算式与迁移前逐行一致
//! （包括 f32 的累加顺序），因为 identity/coverage 会直接决定 `OupParams` 的过滤结果与
//! `identity_gap_compressed` 的阈值判断，任何取值差异都会放大成“命中集合漂移”。
//! 历史实现的对应关系见 `docs/migration-mm2-to-rammap.md`。

use gskits::dna::reverse_complement;
use rust_htslib::bam::record::{Cigar, CigarString};

use super::Hit;

/// op 编码（0=M,1=I,2=D,3=N,4=S,5=H,7==,8=X）→ htslib Cigar。
pub fn ops_to_htslib_cigar(cigar: &[(u32, u8)]) -> CigarString {
    let ops = cigar
        .iter()
        .map(|&(cnt, op)| match op {
            0 => Cigar::Match(cnt),
            1 => Cigar::Ins(cnt),
            2 => Cigar::Del(cnt),
            3 => Cigar::RefSkip(cnt),
            4 => Cigar::SoftClip(cnt),
            5 => Cigar::HardClip(cnt),
            6 => Cigar::Pad(cnt),
            7 => Cigar::Equal(cnt),
            8 => Cigar::Diff(cnt),
            v => panic!("invalid cigar op :{}", v),
        })
        .collect::<Vec<_>>();
    CigarString(ops)
}

/// 给比对区间的两侧补 soft clip，得到写 BAM 用的完整 CIGAR。
///
/// `hit.query_start/query_end` 是正向 query 坐标；反向比对时先换算成
/// 反向补链坐标系下的区间再补 clip（与 BAM 输出约定一致）。
pub fn ops_to_record_cigar(
    cigar: &[(u32, u8)],
    query_start: usize,
    query_end: usize,
    query_len: usize,
    is_rev: bool,
) -> CigarString {
    let (query_start, query_end) = if is_rev {
        (query_len - query_end, query_len - query_start)
    } else {
        (query_start, query_end)
    };

    let mut cigar_str = CigarString(vec![]);
    if query_start > 0 {
        cigar_str.push(Cigar::SoftClip(query_start as u32));
    }
    cigar_str.extend(ops_to_htslib_cigar(cigar).0);
    if query_end != query_len {
        cigar_str.push(Cigar::SoftClip((query_len - query_end) as u32));
    }
    cigar_str
}

/// 把 hit 在 target/query 上的比对区间展开成两条等长字符串（gap 用 `-` 表示）。
pub fn aligned_2_str(hit: &Hit, target: &[u8], query: &[u8]) -> (String, String) {
    let sbr_target = &target[hit.target_start..hit.target_end];
    let sbr_target = String::from_utf8(sbr_target.to_vec()).unwrap();

    let sbr_query = &query[hit.query_start..hit.query_end];
    let sbr_query = if !hit.is_reverse {
        String::from_utf8(sbr_query.to_vec()).unwrap()
    } else {
        String::from_utf8(reverse_complement(sbr_query)).unwrap()
    };

    let mut aligned_target = String::new();
    let mut aligned_query = String::new();

    let mut target_start = 0;
    let mut query_start = 0;

    hit.cigar_ops.iter().for_each(|&(cnt, code)| {
        let cnt = cnt as usize;
        match code {
            0 | 7 | 8 => {
                aligned_target.push_str(&sbr_target[target_start..target_start + cnt]);
                aligned_query.push_str(&sbr_query[query_start..query_start + cnt]);
                target_start += cnt;
                query_start += cnt;
            }
            1 => {
                let ins = vec!["-"; cnt].join("");
                aligned_target.push_str(&ins);
                aligned_query.push_str(&sbr_query[query_start..query_start + cnt]);
                query_start += cnt;
            }
            2 => {
                let del = vec!["-"; cnt].join("");
                aligned_target.push_str(&sbr_target[target_start..target_start + cnt]);
                aligned_query.push_str(&del);
                target_start += cnt;
            }
            _ => panic!("invalid code. {}", code),
        }
    });

    (aligned_target, aligned_query)
}

fn num_matched(hit: &Hit) -> u32 {
    hit.cigar_ops
        .iter()
        .filter(|&&(_, op)| {
            assert!(op != 0, "eqx required");
            op == 7
        })
        .map(|&(cnt, _)| cnt)
        .sum::<u32>()
}

/// matching bases / 全部 CIGAR 长度（含 insertion/deletion 碱基）。
pub fn identity(hit: &Hit) -> f32 {
    let matched = num_matched(hit);
    let tot = hit.cigar_ops.iter().map(|&(cnt, _)| cnt).sum::<u32>();
    matched as f32 / tot as f32
}

/// 把每个 indel 压缩成 1 列后的 identity。
pub fn identity_gap_compressed(hit: &Hit) -> f32 {
    let matched = num_matched(hit);
    if matched == 0 {
        0.0
    } else {
        let tot = hit
            .cigar_ops
            .iter()
            .map(|&(cnt, op)| if op == 1 || op == 2 { 1 } else { cnt })
            .sum::<u32>();
        matched as f32 / tot as f32
    }
}

/// `(query_end - query_start) / query_len`。
///
/// 注：rammap 的 `Mapping` 不带 `query_len`，因此由调用方传入原始 query 长度；
/// 两者在 asts 中恒等（都是整条 subread 序列长度）。
pub fn query_coverage(hit: &Hit, query_len: usize) -> f32 {
    (hit.query_end - hit.query_start) as f32 / query_len as f32
}

/// `(target_end - target_start) / target_len`。
pub fn target_coverage(hit: &Hit) -> f32 {
    (hit.target_end - hit.target_start) as f32 / hit.target_len as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(cigar_ops: Vec<(u32, u8)>, target_len: usize) -> Hit {
        Hit {
            target_len,
            query_start: 0,
            query_end: 7,
            target_start: 0,
            target_end: 7,
            cigar_ops,
            ..Default::default()
        }
    }

    #[test]
    fn test_identity() {
        // 3= 2I 2= => matched 5 / total 7
        let h = hit(vec![(3, 7), (2, 1), (2, 7)], 10);
        assert!((identity(&h) - 5.0 / 7.0).abs() < 1e-6);
        // 压缩后：3 + 1 + 2 = 6
        assert!((identity_gap_compressed(&h) - 5.0 / 6.0).abs() < 1e-6);
        assert!((query_coverage(&h, 10) - 0.7).abs() < 1e-6);
        assert!((target_coverage(&h) - 0.7).abs() < 1e-6);
    }

    #[test]
    fn test_record_cigar_softclip() {
        let cigar = vec![(5, 7), (2, 2)];
        let c = ops_to_record_cigar(&cigar, 2, 9, 12, false);
        assert_eq!(c.to_string(), "2S5=2D3S");

        // reverse: forward q [2,9) -> rc 坐标 [3,10)
        let c = ops_to_record_cigar(&cigar, 2, 9, 12, true);
        assert_eq!(c.to_string(), "3S5=2D2S");
    }
}
