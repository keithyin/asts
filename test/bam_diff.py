import pysam
import argparse
from tqdm import tqdm
import os
import sys

sys.path.append(os.path.abspath(__file__).rsplit("/", maxsplit=1)[0])
import cigar_identity


def compare_bam_files(bam1_path, bam2_path):
    # 打开 BAM 文件
    bam1 = pysam.AlignmentFile(bam1_path, "rb", threads=40)
    bam2 = pysam.AlignmentFile(bam2_path, "rb", threads=40)

    # 构建第一个 BAM 文件的 QNAME 到记录的映射
    bam1_records = {}
    for read in tqdm(bam1.fetch(until_eof=True), desc=f"reading {bam1_path}"):
        bam1_records[read.query_name] = read

    # 初始化统计信息
    total_records = 0
    inconsistent_records = 0

    print("Inconsistent Records:")
    print("QNAME\tSEQ1\tSEQ2\tCIGAR1\tCIGAR2")

    # 遍历第二个 BAM 文件，比较记录
    for read in tqdm(bam2.fetch(until_eof=True), desc=f"reading {bam2_path}"):
        total_records += 1
        qname = read.query_name
        if qname in bam1_records:
            bam1_read = bam1_records[qname]
            if read.is_reverse != bam1_read.is_reverse:
                continue

            # 比较 SEQ 和 CIGAR
            if (
                bam1_read.query_sequence != read.query_sequence
                or bam1_read.cigarstring != read.cigarstring
            ):

                if cigar_identity.calculate_identity(
                    bam1_read.cigarstring
                ) >= cigar_identity.calculate_identity(read.cigarstring):
                    continue
                inconsistent_records += 1
                # print(
                #     f"{qname}\n{bam1_read.query_sequence}\n{read.query_sequence}\n"
                #     f"{bam1_read.cigarstring}\n{read.cigarstring}"
                # )
                # break
        else:
            # 如果 BAM1 中没有找到对应的 QNAME，视为不一致
            inconsistent_records += 1
            # print(
            #     f"{qname}\tMISSING\t{read.query_sequence}\tMISSING\t{read.cigarstring}"
            # )

    # 计算不一致比例
    if total_records > 0:
        inconsistent_ratio = inconsistent_records / total_records
    else:
        inconsistent_ratio = 0

    # 打印统计结果
    print("\nSummary:")
    print(f"Total records in BAM2: {total_records}")
    print(f"Inconsistent records: {inconsistent_records}")
    print(f"Inconsistent ratio: {inconsistent_ratio:.2%}")

    # 关闭 BAM 文件
    bam1.close()
    bam2.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="")
    parser.add_argument("bam1")
    parser.add_argument("bam2")

    args = parser.parse_args()
    compare_bam_files(args.bam1, args.bam2)
