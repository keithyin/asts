build:
	cargo build --release

# 金标准快照回归：固定用例的比对结果与检入的期望值逐字段比对。
# 改动 rammap 版本或比对参数后必须跑一次。
test-golden:
	cargo test --test alignment_golden -- --nocapture

# 比较两次 asts 运行的输出（例如升级 rammap 前后、参数调整前后）。
eqv:
	@echo "usage: python3 scripts/compare_asts_runs.py --a A.bam --b B.bam [--threshold 0.999]"

install:
	cp target/release/asts /usr/bin/
	cp target/release/asrtc /usr/bin/
	cp target/release/asfstc /usr/bin/

bai:
	cargo build --release
	cp target/release/asts /usr/bin/
	cp target/release/asrtc /usr/bin/
	cp target/release/asrtc_poly /usr/bin/
	cp target/release/asfstc /usr/bin/

clean:
	rm -rf target
