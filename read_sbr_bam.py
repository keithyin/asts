import pysam

def read_bam(bam_file: str):
    with pysam.AlignmentFile(bam_file, mode="rb", threads=40, check_sq=False) as bam_h:
        ch_sbr_cnts = {}
        for record in bam_h.fetch(until_eof=True):
            ch = int(record.get_tag("ch"))
            ch_sbr_cnts.setdefault(ch, 0)
            ch_sbr_cnts[ch] += 1
    
    # res_str = "\n".join([f"{ch}-{cnt}" for ch, cnt in list(ch_sbr_cnts.items())])
    # print(res_str)
    return ch_sbr_cnts

def parse_log(log_file: str):
    ch_sbr_cnts = {}

    with open(log_file) as fh:
        for line in fh:
            if "sbr_cnt:" in line:
                infos = line.strip().split("sbr_cnt:")[1].split("/")
                ch = int(infos[1])
                cnt = int(infos[2].split("-")[1])
                ch_sbr_cnts[ch] = cnt
    return ch_sbr_cnts

def eq(ch_sbr_cnts_real, ch_sbr_cnt_n_real):
    for k in ch_sbr_cnt_n_real:
        if ch_sbr_cnts_real[k] != ch_sbr_cnt_n_real[k]:
            print(f"{k}, {ch_sbr_cnts_real[k]} != {ch_sbr_cnt_n_real[k]}")

if __name__ == "__main__":
    ch_sbr_cnts1 = read_bam("/root/projects/asts/test_data/sbr.bam")
    ch_sbr_cnts2 = parse_log("/root/projects/asts/asts.log")
    eq(ch_sbr_cnts1, ch_sbr_cnts2)
