#!/usr/bin/env python3
"""比较两次 asts 运行的输出是否等价。

典型用途：升级 `rammap-core` 版本、调整比对参数、或替换比对引擎前后，
用同一份输入各跑一次 asts，再逐条比对两份 BAM。

比较字段：flag / rname / pos / mapq / CIGAR / cs(md)，输出
「qname 一致率 / record 一致率 / 每个 channel 的命中集合 Jaccard」，并以退出码表达结论。

用法::

    # 1) 基线（例如升级前的版本）
    asts -q test_data/sbr.bam -t test_data/smc.bam -p /tmp/run_a --threads 10
    # 2) 新版本 / 新参数（同一份输入）
    asts -q test_data/sbr.bam -t test_data/smc.bam -p /tmp/run_b --threads 10
    # 3) 对比
    scripts/compare_asts_runs.py --a /tmp/run_a.bam --b /tmp/run_b.bam

退出码：0 = 达到门槛；1 = 未达到；2 = 输入/环境错误。

依赖：python3 标准库 + `samtools`（asts 运行时本来就依赖它）。
"""

from __future__ import annotations

import argparse
import collections
import shutil
import subprocess
import sys
from typing import Dict, List, Tuple

# 仅这些字段参与比较：它们共同决定下游 MSA / identity / 统计
TAG_KEYS = ("cs:Z:", "md:Z:")


def samtools_view(path: str) -> List[str]:
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools not found in PATH")
    out = subprocess.run(
        ["samtools", "view", path],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return out.stdout.decode("utf-8", "replace").splitlines()


def parse_cs(tags: List[str]) -> str:
    for tag in tags:
        for key in TAG_KEYS:
            if tag.startswith(key):
                return tag
    return ""


def normalize(path: str) -> Tuple[Dict[str, collections.Counter], Dict[str, set], int]:
    """返回 (qname -> 记录多重集, channel -> qname 集合, 记录总数)。"""
    by_qname: Dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    by_channel: Dict[str, set] = collections.defaultdict(set)
    n_records = 0

    for line in samtools_view(path):
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        qname, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
        tags = f[11:]
        key = (flag, rname, pos, mapq, cigar, parse_cs(tags))
        by_qname[qname][key] += 1
        n_records += 1

        ch = None
        for tag in tags:
            if tag.startswith("ch:i:"):
                ch = tag[5:]
                break
        by_channel[ch or qname].add(qname)

    return by_qname, by_channel, n_records


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def main() -> int:
    ap = argparse.ArgumentParser(description="compare the BAM outputs of two asts runs")
    ap.add_argument("--a", required=True, help="BAM of run A (baseline)")
    ap.add_argument("--b", required=True, help="BAM of run B (candidate)")
    ap.add_argument(
        "--threshold",
        type=float,
        default=0.999,
        help="minimum acceptable agreement ratio (default: 0.999)",
    )
    ap.add_argument("--max-examples", type=int, default=10)
    args = ap.parse_args()

    try:
        a_q, a_ch, a_n = normalize(args.a)
        b_q, b_ch, b_n = normalize(args.b)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 2

    qnames = set(a_q) | set(b_q)
    if not qnames:
        print("[error] no records found in either BAM", file=sys.stderr)
        return 2

    diff_queries = [q for q in qnames if a_q.get(q) != b_q.get(q)]
    q_agree = 1.0 - len(diff_queries) / len(qnames)

    a_recs = sum(sum(c.values()) for c in a_q.values())
    b_recs = sum(sum(c.values()) for c in b_q.values())
    matched_recs = sum(
        sum((a_q[q] & b_q[q]).values()) if (q in a_q and q in b_q) else 0 for q in qnames
    )
    rec_agree = matched_recs / max(a_recs, b_recs)

    channels = set(a_ch) | set(b_ch)
    jac = sum(jaccard(a_ch.get(c, set()), b_ch.get(c, set())) for c in channels) / max(len(channels), 1)

    print(f"records:        A={a_recs}  B={b_recs}  matched={matched_recs}")
    print(f"qname 一致率:   {q_agree:.6f}  ({len(qnames) - len(diff_queries)}/{len(qnames)})")
    print(f"record 一致率:  {rec_agree:.6f}  ({matched_recs}/{max(a_recs, b_recs)})")
    print(f"channel Jaccard: {jac:.6f}  ({len(channels)} channels)")

    for q in diff_queries[: args.max_examples]:
        only_a = a_q.get(q, collections.Counter()) - b_q.get(q, collections.Counter())
        only_b = b_q.get(q, collections.Counter()) - a_q.get(q, collections.Counter())
        print(f"  diff qname={q}")
        for k in list(only_a)[:3]:
            print(f"    A only: {k}")
        for k in list(only_b)[:3]:
            print(f"    B only: {k}")

    ok = q_agree >= args.threshold and rec_agree >= args.threshold and jac >= args.threshold
    print("RESULT:", "PASS" if ok else "FAIL", f"(threshold={args.threshold})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
