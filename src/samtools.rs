//! samtools 命令行封装（原 `gskits::samtools`，行为逐条一致）。
//!
//! `sort_by_tag` 已在 `crate::utils` 中有带缓存的版本，此处只保留
//! `samtools_bai` 与 `sort_by_coordinates`。

use std::process::Command;

/// `samtools index -@ {threads} {bam_file} {bam_file}.bai`
pub fn samtools_bai(bam_file: &str, force: bool, threads: Option<usize>) -> anyhow::Result<()> {
    let index_file = format!("{}.bai", bam_file);

    let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
    let threads = threads.max(1);

    if force && std::path::Path::new(&index_file).exists() {
        std::fs::remove_file(&index_file).expect("remove index file error");
    }

    let mut cmd = Command::new("samtools");
    cmd.args([
        "index",
        "-@",
        threads.to_string().as_str(),
        bam_file,
        index_file.as_str(),
    ]);

    let status = cmd
        .status()
        .map_err(|_| anyhow::anyhow!("Run cmd error: {:?}", cmd))?;
    if !status.success() {
        return Err(anyhow::anyhow!(
            "Run cmd error, exit status not succ: {:?}",
            cmd
        ));
    }

    Ok(())
}

/// `samtools sort -o {bam_file} -@ {threads} {bam_file}`（原地排序）。
pub fn sort_by_coordinates(bam_file: &str, threads: Option<usize>) {
    let threads = threads.unwrap_or(num_cpus::get_physical() / 2);
    let threads = threads.max(1);

    let mut cmd = Command::new("samtools");
    cmd.args([
        "sort",
        "-o",
        bam_file,
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
}
