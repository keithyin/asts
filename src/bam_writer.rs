//! BAM 结果写出线程。
//!
//! 迁移说明：写库已从 rust-htslib 迁至 noodles。写出顺序、header 内容、
//! PG 记录、`iy`/`ec` 两个 tag 的追加顺序都与迁移前完全一致；唯一的已知差异是
//! noodles 固定按 `@HD → @SQ... → @PG` 的顺序输出 header（迁移前为 `@HD → @PG → @SQ...`），
//! 记录内容不变。相关背景见 `docs/migration-mm2-to-rammap.md` 与
//! `docs/ci-htslib-opaque-structs.md`。

use std::collections::HashMap;
use std::fs::File;
use std::num::NonZeroUsize;

use crossbeam::channel::Receiver;
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::io::Write as _,
    alignment::record::data::field::Tag,
    alignment::record_buf::data::field::Value,
    header::record::value::{map, Map},
};
use sam::header::record::value::map::header::Version;

use crate::gsbam::bam_record_ext::BamRecordExt;
use crate::pbar::{self, DEFAULT_INTERVAL};
use crate::utils::command_line_str;
use crate::{params::OupParams, AlignResult};

/// target_idx: &HashMap<target_name, (idx, target_len)>
pub fn write_bam_worker(
    recv: Receiver<AlignResult>,
    target_idx: &HashMap<String, (usize, usize)>,
    o_path: &str,
    oup_args: &OupParams,
    pg_name: &str,
    version: &str,
    enable_pb: bool,
) {
    let mut hd = Map::<map::Header>::new(Version::new(1, 5));
    hd.other_fields_mut()
        .insert(map::header::tag::SORT_ORDER, "unknown".into());

    let mut pg = Map::<map::Program>::default();
    pg.other_fields_mut()
        .insert(map::program::tag::NAME, pg_name.into());
    pg.other_fields_mut()
        .insert(map::program::tag::COMMAND_LINE, command_line_str().into());
    pg.other_fields_mut()
        .insert(map::program::tag::VERSION, version.into());

    let mut targets = target_idx.iter().collect::<Vec<_>>();
    targets.sort_by_key(|tup| tup.1 .0);

    let mut header_builder = sam::Header::builder().set_header(hd).add_program(pg_name, pg);
    for target in targets {
        header_builder = header_builder.add_reference_sequence(
            target.0.as_str(),
            Map::<map::ReferenceSequence>::new(NonZeroUsize::new(target.1 .1).unwrap()),
        );
    }
    let header = header_builder.build();

    let pb = if enable_pb {
        Some(pbar::get_spin_pb(
            format!("{}: writing alignment result to {}", pg_name, o_path),
            DEFAULT_INTERVAL,
        ))
    } else {
        None
    };

    let mut writer = bam::io::Writer::new(File::create(o_path).unwrap());
    writer.write_header(&header).unwrap();

    for align_res in recv {
        if let Some(pb) = pb.as_ref() {
            pb.inc(1);
        }
        for mut record in align_res.records {
            if oup_args.valid(&record) {
                let record_ext = BamRecordExt::new(&record);

                let iy = record_ext.compute_identity();
                let ec = record_ext.compute_query_coverage();

                let data = record.data_mut();
                data.insert(Tag::new(b'i', b'y'), Value::Float(iy));
                data.insert(Tag::new(b'e', b'c'), Value::Float(ec));

                writer.write_alignment_record(&header, &record).unwrap();
            }
        }
    }
    writer.try_finish().unwrap();

    if let Some(pb) = pb.as_ref() {
        pb.finish();
    }
}
