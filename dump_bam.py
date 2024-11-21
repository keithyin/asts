import pysam

def dump_bam(bam_file: str):
    o_bam = "{}.part.bam".format(bam_file.rsplit(".", maxsplit=1)[0])
    with pysam.AlignmentFile(bam_file, mode="rb", threads=10, check_sq=False) as bam_h:
        with pysam.AlignmentFile(o_bam, mode="wb", threads=10, check_sq=False, header=bam_h.header) as out_bam_h:

            for record in bam_h.fetch(until_eof=True):
                ch = int(record.get_tag("ch"))
                if ch % 2 == 0:
                    out_bam_h.write(record)
    
def eq(ch_sbr_cnts_real, ch_sbr_cnt_n_real):
    for k in ch_sbr_cnts_real:
        if ch_sbr_cnts_real[k] != ch_sbr_cnt_n_real[k]:
            print(f"{k}, {ch_sbr_cnts_real[k]} != {ch_sbr_cnt_n_real[k]}")

if __name__ == "__main__":
    ch_sbr_cnts1 = dump_bam("/root/projects/asts/test_data/smc.bam")
