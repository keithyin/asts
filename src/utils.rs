use std::{fmt::Debug, process::Command, str::FromStr, time};

use serde::Serialize;

pub fn command_line_str() -> String {
    let args: Vec<String> = std::env::args().collect();

    args.join(" ")
}

/// 质量值列表（phred，无偏移）→ 平均准确率（原 `gskits::phreq::phreq_list_2_quality`）。
pub fn phreq_list_2_quality(baseq: &[u8]) -> Option<f32> {
    if baseq.len() == 0 {
        return None;
    }
    let length = baseq.len() as f64;
    let err_rate = baseq
        .iter()
        .map(|v| *v as f64)
        .map(|v| 10.0_f64.powf(v / -10.0_f64))
        .reduce(|acc, v| acc + v)
        .unwrap()
        / length;
    Some((1.0_f64 - err_rate) as f32)
}

/// 原地计时累加器（原 `gskits::utils::ScopedTimer`，自身不打印，调用方格式化）。
#[derive(Debug)]
pub struct ScopedTimer {
    iter: usize,
    elapsed: u128, // nano
}

impl ScopedTimer {
    pub fn new() -> Self {
        ScopedTimer {
            iter: 0,
            elapsed: 0,
        }
    }

    pub fn perform_timing(&mut self) -> Timer<'_> {
        Timer::new(self)
    }

    pub fn reset(&mut self) {
        self.iter = 0;
        self.elapsed = 0;
    }

    /// None: iter/nano_sec; 1000: iter/micro; 1000_000: iter/milli; 1000_000_000: iter/sec ...
    pub fn speed(&self, multiplier: Option<u128>) -> f64 {
        let multiplier = multiplier.unwrap_or(1);
        assert!(multiplier > 0);
        self.iter as f64 / (self.elapsed as f64 / multiplier as f64)
    }
}

impl Default for ScopedTimer {
    fn default() -> Self {
        Self::new()
    }
}

impl ToString for ScopedTimer {
    fn to_string(&self) -> String {
        format!("{:.10} iter/nano_secs", (self.iter as f64 / self.elapsed as f64))
    }
}

pub struct Timer<'a> {
    scoped_timer: &'a mut ScopedTimer,
    instant: time::Instant,
}

impl<'a> Timer<'a> {
    fn new(scoped_timer: &'a mut ScopedTimer) -> Self {
        let instant = time::Instant::now();
        Timer {
            scoped_timer,
            instant,
        }
    }

    pub fn done_with_cnt(self, cnt: usize) {
        let elapsed = self.instant.elapsed().as_nanos();
        self.scoped_timer.iter += cnt;
        self.scoped_timer.elapsed += elapsed;
    }

    pub fn speed(&self, multiplier: Option<u128>) -> f64 {
        self.scoped_timer.speed(multiplier)
    }
}

/// parse the range str to Vector of tuple
/// range fmt: b1:e1,b2:e2,be3,b4:e4. corresponds to [[b1, e1], [b2, e2], [be3, be3], [b4, e4]]. left inclusive, right inclusive
#[allow(unused)]
pub fn range_parser<T>(range_str: &str) -> Vec<(T, T)>
where
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    range_str
        .trim()
        .split(",")
        .map(|single_range| {
            let items = single_range
                .trim()
                .split(":")
                .filter(|item| item.trim().len() != 0)
                .collect::<Vec<_>>();
            return if items.len() == 2 {
                (
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                    items[1].parse::<T>().expect(&format!("{}", items[1])),
                )
            } else {
                (
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                )
            };
        })
        .collect::<Vec<_>>()
}

#[derive(Debug, Clone, Serialize)]
pub struct Range<T> {
    range: Vec<(T, T)>,
}

impl<T> Range<T>
where
    T: PartialOrd,
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    pub fn new(range_str: &str) -> Self {
        Self {
            range: range_parser::<T>(range_str),
        }
    }

    pub fn new_range(range: Vec<(T, T)>) -> Self {
        Self { range }
    }

    pub fn within_range(&self, v: T) -> bool {
        let mut result = false;
        for (b, e) in &self.range {
            if *b <= v && v <= *e {
                result = true;
                break;
            }
        }
        result
    }
}

impl Range<usize> {
    pub fn shift(&mut self, shift: i64) {
        self.range.iter_mut().for_each(|(start, end)| {
            let new_start = (*start as i64 + shift).max(0) as usize;
            let new_end = (*end as i64 + shift).max(0) as usize;
            *start = new_start;
            *end = new_end;
        });
    }
}

impl Range<i64> {
    pub fn shift(&mut self, shift: i64) {
        self.range.iter_mut().for_each(|(start, end)| {
            let new_start = *start as i64 + shift;
            let new_end = *end as i64 + shift;
            *start = new_start;
            *end = new_end;
        });
    }
}

// pub fn sort_by_tag(bam_file: &str, tag: &str, threads: Option<usize>) -> String {
//     let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
//     let out_bam = format!(
//         "{}.sort_by_{}.bam",
//         bam_file.rsplit_once(".").unwrap().0,
//         tag
//     );

//     let mut cmd = Command::new("samtools");
//     cmd.args([
//         "sort",
//         "-n",
//         "-t",
//         tag,
//         "-o",
//         out_bam.as_str(),
//         "-@",
//         threads.to_string().as_str(),
//         bam_file,
//     ]);

//     let oup = cmd.output().expect(&format!("sort {} error", bam_file));
//     if !oup.status.success() {
//         panic!(
//             "sort {} error. {}",
//             bam_file,
//             String::from_utf8(oup.stderr).unwrap()
//         );
//     }

//     out_bam
// }

pub fn sort_by_tag(bam_file: &str, tag: &str, threads: Option<usize>) -> String {
    let out_bam = format!(
        "{}.sort_by_{}.bam",
        bam_file.rsplit_once(".").unwrap().0,
        tag
    );

    // 如果输出文件存在且比输入文件新，则直接返回
    if let (Ok(input_meta), Ok(output_meta)) =
        (std::fs::metadata(bam_file), std::fs::metadata(&out_bam))
    {
        if let (Ok(input_time), Ok(output_time)) = (input_meta.modified(), output_meta.modified()) {
            if output_time > input_time {
                return out_bam;
            }
        }
    }

    // 否则执行排序
    let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
    let mut cmd = Command::new("samtools");
    cmd.args([
        "sort",
        "-n",
        "-t",
        tag,
        "-o",
        out_bam.as_str(),
        "-@",
        threads.to_string().as_str(),
        bam_file,
    ]);

    let oup = cmd.output().expect(&format!("sort {} error", bam_file));
    if !oup.status.success() {
        panic!(
            "sort {} error. {}",
            bam_file,
            String::from_utf8(oup.stderr).unwrap()
        );
    }

    out_bam
}
