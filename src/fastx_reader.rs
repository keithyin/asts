//! FASTA/FASTQ 读取（原 `gskits::fastx_reader`，hand-rolled 解析器，替代 bio crate）。
//!
//! 命名约定与原实现一致：header 去掉 `>`/`@` 后取第一个空白分隔字段作为 name。
//! fastq 的 qual 存储 phred 值（ASCII 减 33），与 `crate::read_info::ReadInfo` 的
//! "phreq, no offset" 约定一致。

use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::read_info::ReadInfo;

fn open_file(fname: &str, msg: &str) -> BufReader<File> {
    let file = File::open(fname).expect(&format!("{}: {}", msg, fname));
    BufReader::new(file)
}

/// header 行（已去掉 `>`/`@`）→ name：第一个空白分隔字段。
fn header_line_to_name(header: &str) -> String {
    header
        .trim_end_matches(['\r', '\n'])
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

pub struct FastaFileReader {
    fname: String,
    reader: BufReader<File>,
    pending_header: Option<String>,
    eof: bool,
}

impl FastaFileReader {
    pub fn new(fname: String) -> Self {
        let reader = open_file(&fname, "open file error");
        Self {
            fname,
            reader,
            pending_header: None,
            eof: false,
        }
    }

    #[allow(dead_code)]
    pub fn get_fname(&self) -> &str {
        &self.fname
    }
}

impl Iterator for FastaFileReader {
    type Item = ReadInfo;

    fn next(&mut self) -> Option<ReadInfo> {
        // 定位到下一个 header（跳过文件头部的杂散内容）
        while self.pending_header.is_none() {
            if self.eof {
                return None;
            }
            let mut line = String::new();
            let n = self
                .reader
                .read_line(&mut line)
                .expect(&format!("read fasta error {}", self.fname));
            if n == 0 {
                self.eof = true;
                return None;
            }
            if line.starts_with('>') {
                self.pending_header = Some(line);
            }
        }

        let header = self.pending_header.take().unwrap();
        let name = header_line_to_name(&header[1..]);

        let mut seq: Vec<u8> = vec![];
        loop {
            let mut line = String::new();
            let n = self
                .reader
                .read_line(&mut line)
                .expect(&format!("read fasta error {}", self.fname));
            if n == 0 {
                self.eof = true;
                break;
            }
            if line.starts_with('>') {
                // 下一个 record 的 header，暂存给下一次 next()
                self.pending_header = Some(line);
                break;
            }
            seq.extend(line.trim_end_matches(['\r', '\n']).as_bytes());
        }

        Some(ReadInfo::new_fa_record(
            name,
            unsafe { String::from_utf8_unchecked(seq) },
        ))
    }
}

pub struct FastqFileReader {
    fname: String,
    reader: BufReader<File>,
    eof: bool,
}

impl FastqFileReader {
    pub fn new(fname: String) -> Self {
        let reader = open_file(&fname, "open file error");
        Self {
            fname,
            reader,
            eof: false,
        }
    }
}

impl Iterator for FastqFileReader {
    type Item = ReadInfo;

    fn next(&mut self) -> Option<ReadInfo> {
        if self.eof {
            return None;
        }

        let mut header = String::new();
        let n = self
            .reader
            .read_line(&mut header)
            .expect(&format!("read fastq record error. {}", self.fname));
        if n == 0 {
            self.eof = true;
            return None;
        }

        let mut seq = String::new();
        let n = self
            .reader
            .read_line(&mut seq)
            .expect(&format!("read fastq record error. {}", self.fname));
        if n == 0 {
            panic!("read fastq record error. {}", self.fname);
        }

        let mut plus = String::new();
        let n = self
            .reader
            .read_line(&mut plus)
            .expect(&format!("read fastq record error. {}", self.fname));
        if n == 0 {
            panic!("read fastq record error. {}", self.fname);
        }

        let mut qual_line = String::new();
        let n = self
            .reader
            .read_line(&mut qual_line)
            .expect(&format!("read fastq record error. {}", self.fname));
        if n == 0 {
            panic!("read fastq record error. {}", self.fname);
        }

        let seq = seq.trim_end_matches(['\r', '\n']);
        let qual_str = qual_line.trim_end_matches(['\r', '\n']);

        // 与 bio 的 fastq 解析保持一致：序列与质量长度必须相等
        if seq.len() != qual_str.len() {
            panic!(
                "read fastq record error. {}: sequence and quality length mismatch ({} != {})",
                self.fname,
                seq.len(),
                qual_str.len()
            );
        }

        let qual: Vec<u8> = qual_str
            .as_bytes()
            .iter()
            .map(|&v| {
                if v < 33 {
                    panic!(
                        "read fastq record error. {}: invalid quality character {:?}",
                        self.fname, v as char
                    );
                }
                v - 33
            })
            .collect();

        Some(ReadInfo::new_fq_record(
            header_line_to_name(&header[1..]),
            seq.to_string(),
            qual,
        ))
    }
}

/// 把 fastx reader 的迭代收集为 Vec。
pub fn read_fastx<T: IntoIterator<Item = ReadInfo>>(reader: T) -> Vec<ReadInfo> {
    reader.into_iter().collect()
}
