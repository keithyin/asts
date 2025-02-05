import pysam
from tqdm import tqdm


def read_bam(bam_file: str):
    with pysam.AlignmentFile(bam_file, mode="rb", threads=40) as f:
        for record in tqdm(f.fetch()):
            iy = float(record.get_tag("iy"))
            if iy < 0.7:
                raise ValueError(f"{record.query_name}. identity:{iy}")


if __name__ == "__main__":
    fn = "/data/ccs_data/ccs_eval2024q3/jinpu/smc501/20240711_Sync_Y0006_02_H01_Run0001_called.subreadsTOnn-polish-asts090-icing061.smc_all_reads.0.aligned.bam"
    read_bam(fn)
