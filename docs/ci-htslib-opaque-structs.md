# CI 排障：`hts_sys::BGZF` / `sam_hdr_t` 变成不透明结构体

> 症状记忆点：编译 `rust-htslib` 时报一批 `E0609`，提示 `available field is: _address`。
> 关键词：`no field block_address on type hts_sys::BGZF`、`no field text on type sam_hdr_t`。

## 1. 症状

```
error[E0609]: no field `block_address` on type `hts_sys::BGZF`
   --> rust-htslib-0.49.0/src/bam/mod.rs:198:15
    = note: available field is: `_address`

error[E0609]: no field `text` on type `sam_hdr_t`
   --> rust-htslib-0.49.0/src/bam/mod.rs:1129:20
    = note: available field is: `_address`
...
error: could not compile `rust-htslib` (lib) due to 11 previous errors
```

`_address: u8` 是 bindgen 对**只有前向声明**的结构体输出的占位字段。也就是说：`hts-sys` 生成的
Rust 绑定里，`BGZF` / `sam_hdr_t` 是空壳，而 `rust-htslib 0.49` 是直接按字段访问这两个结构体的
（`block_address`、`block_offset`、`text`、`l_text`、`n_targets`、`target_name`、`target_len`），
于是必然编译失败。

## 2. 根因分类

HTSlib 从 **1.21** 起把 `BGZF` / `sam_hdr_t` 等结构体改为**不透明**（访问要用 accessor 函数）。
`rust-htslib 0.49` 是为 HTSlib 1.19 写的，一旦绑定来自 1.21+ 的头文件就会挂。
所以问题不在于 asts 的代码，而在于**CI 上实际被 bindgen 解析的头文件 / 依赖版本**与本地不同。

本仓库锁定的依赖事实（`Cargo.lock`）：

| 组件 | 版本 | 说明 |
| --- | --- | --- |
| `rust-htslib` | `0.49.0` | 对 `hts-sys` 的约束是 `2.2.0`（caret，允许 2.2.x），且**无条件启用 `hts-sys/bindgen`** |
| `hts-sys` | `2.2.0` | 捆绑 HTSlib **1.19.1**（`HTS_VERSION 101901`）；`bgzf.h` 有 `struct BGZF { ... }`、`sam.h` 有 `typedef struct sam_hdr_t { ... }`，**不是不透明** |
| `bindgen` | `0.69.5` | 由 `hts-sys/bindgen` 触发，用 libclang 解析 `wrapper.h` |

补充事实（排查时确认过，避免重复怀疑）：

- `hts-sys` 的 `build.rs` **没有** pkg-config / 系统库探测逻辑，只编译捆绑的 `htslib/*.c`；
- `wrapper.h` 用**引号**包含 `"htslib/htslib/hts.h"`，按相对路径解析，系统头文件不会覆盖它；
- 关掉 bindgen 时的回退文件 `linux_prebuilt_bindings.rs` **同样包含** `block_address` / `text` /
  `l_text` / `n_targets` / `target_name` / `target_len`（已逐字段确认）。

也就是说：**正常的 2.2.x + bindgen（或预构建）都不会产生 `_address`**，CI 上一定发生了环境/版本漂移。

## 3. 排查步骤（复现失败时的第一手命令）

```bash
# 1) 干净 target 目录复现，避免增量缓存造成假象
CARGO_TARGET_DIR=/tmp/asts_ci_repro cargo build --locked

# 2) 确认实际依赖版本（CI 失败时看日志里这段）
cargo tree -i hts-sys --locked
cargo tree -i rust-htslib --locked

# 3) 直接看 hts-sys 生成的绑定是不是空壳
grep -n -A 8 "pub struct BGZF" target/debug/build/hts-sys-*/out/bindings.rs
```

判据：

| 观察 | 含义 |
| --- | --- |
| `BGZF` 有 `block_address` 等字段 | 绑定正常，问题在别处（找 `rust-htslib` 的真实报错行） |
| `BGZF` 只有 `_address` | bindgen 解析到的不是捆绑的 1.19 头文件，或用了 1.21+ 的系统头 |
| `cargo tree -i hts-sys` 不是 2.2.0/2.2.1 | `Cargo.lock` 没有被使用（改成 `--locked` 即可暴露） |

## 4. 本仓库采用的修复（方案 A：固定 CI 环境）

`.github/workflows/rust.yml` 现在做四件事：

1. `sudo apt-get install -y libhts-dev`：把系统 HTSlib 固定为发行版自带的 **1.19**（非不透明），
   即使 bindgen 意外回落到系统头也不会拿到 1.21+ 的定义；
2. 安装并固定 `libclang`（优先 `libclang-18-dev`，回退 `libclang-dev`），导出 `LIBCLANG_PATH`；
3. `cargo build/test --locked`：强制使用仓库里的 `Cargo.lock`，杜绝 runner 侧重新解析依赖；
4. 增加诊断 step（打印 `clang --version`、`LIBCLANG_PATH`、`pkg-config --modversion htslib`、
   `cargo tree -i hts-sys`）以及失败时的 bindings dump step。

## 5. 未采用的方案（记录取舍，避免重复讨论）

**方案 B：vendored `hts-sys` / `rust-htslib`，强制使用预构建 bindings。**
做法是把 `hts-sys` 复制进 `vendor/` 并改 `build.rs` 永远 copy `linux_prebuilt_bindings.rs`
（或在 vendored `rust-htslib` 的 `Cargo.toml` 里去掉 `hts-sys` 的 `bindgen` feature）。
优点：彻底不依赖 libclang 与头文件搜索顺序，确定性最强。
缺点：仓库里多出一个第三方 crate 的分叉（`hts-sys` 约 17 MB，`rust-htslib` 约 900 KB），
每次升级都要手工同步。**本项目选择不引入此类分叉**，仅在方案 A 无法稳定时再考虑。

**方案 C（根治，待上游）：升级 `rust-htslib` 到 0.50+/1.x。**
新版本适配了 HTSlib 1.21 的不透明结构体。但 `gskits 0.15.15`（当前最新）依赖
`rust-htslib = "0.49"`，且 asts 与 gskits 之间**直接传递 `rust_htslib::bam::Record`**
（`BamRecordExt::new(&record)`、`BamWriter` 等），两个 crate 版本会导致类型不兼容。
因此升级必须 `gskits` 与 `asts` 同步进行，属独立任务。

## 6. 升级 `hts-sys` / `rust-htslib` 时的检查清单

- [ ] 升级后立刻看 `cargo tree -i hts-sys`，确认 bundled HTSlib 版本（`grep HTS_VERSION
      ~/.cargo/registry/src/*/hts-sys-*/htslib/htslib/hts.h`），**≥1.21 就必须同时升级 rust-htslib**
- [ ] 若升级 `rust-htslib`，先确认 `gskits` 的 `rust-htslib` 约束同步放宽，否则会出现两份
      `rust_htslib`（类型不兼容，编译期就会报 `expected Record, found Record`）
- [ ] 本地用干净 target 目录跑 `cargo build --locked && cargo test --locked`
- [ ] 跑 `cargo test --test alignment_golden`（比对语义不受 htslib 版本影响，但能确认 BAM 读写正常）
- [ ] 推送后确认 CI 的 `Diagnose toolchain and deps` 与本机一致

## 7. 相关文档

- 比对引擎迁移（rammap）与等价性验证：[`migration-mm2-to-rammap.md`](migration-mm2-to-rammap.md)
