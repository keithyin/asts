//! `fasta2bam` / `fastq2bam` 移植（原 `gskits::fastx_reader::fastx2bam`）。
//!
//! 行为逐条对齐原实现：
//! - header：`@HD VN:1.5 SO:unknown` + 一个假 `@SQ SN:chr1 LN:1234`（无 PG 行）
//! - 输出路径：`{stem}.bam` 再经 tmp-name 处理（`{stem}-{pid}-{nanos}.bam`，
//!   原实现用 uuid5，这里等价替换，路径仍以 `.bam` 结尾且由调用方负责删除）
//! - 每 read：flags=0、tid=0、无比对位置；fasta 用 `vec![255; len]` 占位 qual，
//!   fastq 用真实 phred（无偏移）；仅写 `ch`（U32）tag
//! - `ch` 从 qname 按 delim 切分后取第 `channel_idx` 个字段，解析失败 panic（与原实现一致）

use std::{fs::File, num::NonZeroUsize};

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::io::Write as _,
    alignment::record::Flags,
    alignment::record_buf::{data::field::Value, QualityScores, RecordBuf, Sequence},
    header::record::value::{
        map::{self, header::Version},
        Map,
    },
};

use crate::fastx_reader::{read_fastx, FastaFileReader, FastqFileReader};
use crate::read_info::ReadInfo;

fn generate_tmp_filename(fname: &str) -> String {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    match fname.rsplit_once('.') {
        Some((prefix, suffix)) => format!("{}-{}-{}.{}", prefix, std::process::id(), nanos, suffix),
        None => format!("{}-{}-{}", fname, std::process::id(), nanos),
    }
}

fn build_header() -> sam::Header {
    let mut hd = Map::<map::Header>::new(Version::new(1, 5));
    hd.other_fields_mut()
        .insert(map::header::tag::SORT_ORDER, "unknown".into());

    sam::Header::builder()
        .set_header(hd)
        .add_reference_sequence(
            "chr1",
            Map::<map::ReferenceSequence>::new(NonZeroUsize::new(1234).unwrap()),
        )
        .build()
}

fn record2bam(read_info: &ReadInfo, delim: &str, channel_idx: usize) -> RecordBuf {
    let items = read_info.name.split(delim).skip(channel_idx).collect::<Vec<_>>();
    let ch: usize = items[0].parse().unwrap();

    // fasta（无 qual）时留空：noodles encoder 对空 qual 自动填 0xFF×len，
    // 与 rust-htslib 写 `vec![255; len]` 的输出逐字节一致（noodles 会校验
    // qual ≤ 93，255 占位值不能直接写入）
    let qualities = match &read_info.qual {
        Some(qual) => QualityScores::from(qual.clone()),
        None => QualityScores::default(),
    };

    RecordBuf::builder()
        .set_name(read_info.name.as_bytes())
        .set_flags(Flags::empty())
        .set_reference_sequence_id(0)
        .set_sequence(Sequence::from(read_info.seq.as_bytes().to_vec()))
        .set_quality_scores(qualities)
        .set_data(
            [(
                sam::alignment::record::data::field::Tag::new(b'c', b'h'),
                Value::UInt32(ch as u32),
            )]
            .into_iter()
            .collect(),
        )
        .build()
}

fn write_records2bam(reads: Vec<ReadInfo>, delim: &str, channel_idx: usize, o_filepath: &str) {
    let header = build_header();
    let mut writer = bam::io::Writer::new(File::create(o_filepath).unwrap());
    writer.write_header(&header).unwrap();

    for read_info in reads.iter() {
        let record = record2bam(read_info, delim, channel_idx);
        writer.write_alignment_record(&header, &record).unwrap();
    }

    writer.try_finish().unwrap();
}

pub fn fasta2bam(fa_filename: &str, delim: &str, channel_idx: usize) -> String {
    let reads = read_fastx(FastaFileReader::new(fa_filename.to_string()));
    let o_filepath = generate_tmp_filename(&format!(
        "{}.bam",
        fa_filename.rsplit_once('.').unwrap().0
    ));
    write_records2bam(reads, delim, channel_idx, &o_filepath);
    o_filepath
}

pub fn fastq2bam(fq_filename: &str, delim: &str, channel_idx: usize) -> String {
    let reads = read_fastx(FastqFileReader::new(fq_filename.to_string()));
    let o_filepath = generate_tmp_filename(&format!(
        "{}.bam",
        fq_filename.rsplit_once('.').unwrap().0
    ));
    write_records2bam(reads, delim, channel_idx, &o_filepath);
    o_filepath
}
