//! 系统资源监控（原 `gskits::sys_monitor`，改用 sysinfo 0.38 API）。

use std::{
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
    thread,
    time::{Duration, Instant},
};

use sysinfo::Pid;

#[derive(Debug, PartialEq, Eq)]
enum SysMonState {
    Init,
    Running,
    Stopped,
}

/// 系统资源监控，通过 tracing 打印资源信息。
pub struct SysMon {
    state: Arc<Mutex<SysMonState>>,
    need_stop: Arc<AtomicBool>,
    monitor_interval: Duration,
    prog_name: Arc<String>,
}

impl SysMon {
    pub fn new(monitor_interval: Duration, prog_name: String) -> Self {
        Self {
            state: Arc::new(Mutex::new(SysMonState::Init)),
            need_stop: Arc::new(AtomicBool::new(false)),
            monitor_interval,
            prog_name: Arc::new(prog_name),
        }
    }

    /// pid: 需要监控的进程，None 时使用当前进程。
    /// worker_threads: 用于计算每线程 CPU 负载，None 时读取 /proc 线程数。
    pub fn start_monitor(&self, pid: Option<u32>, worker_threads: Option<usize>) {
        let state = self.state.clone();
        let need_stop = self.need_stop.clone();
        let monitor_interval = self.monitor_interval;
        let prog_name = self.prog_name.clone();
        thread::spawn(move || {
            monitor(
                pid,
                worker_threads,
                state,
                need_stop,
                monitor_interval,
                prog_name,
            );
        });
    }
}

impl Drop for SysMon {
    fn drop(&mut self) {
        self.need_stop.store(true, Ordering::Relaxed);

        let start = Instant::now();
        while *self.state.lock().unwrap() != SysMonState::Stopped {
            thread::sleep(Duration::from_millis(500));

            if start.elapsed() > Duration::from_secs(60) {
                // avoid monitor thread unexpected stop
                break;
            }
        }
    }
}

fn bytes2_gb(num_bytes: usize) -> usize {
    num_bytes / 1024 / 1024 / 1024
}

struct _SysMonStateGuard {
    state: Arc<Mutex<SysMonState>>,
}

impl From<Arc<Mutex<SysMonState>>> for _SysMonStateGuard {
    fn from(value: Arc<Mutex<SysMonState>>) -> Self {
        assert!(*value.lock().unwrap() == SysMonState::Init);
        *value.lock().unwrap() = SysMonState::Running;

        Self { state: value }
    }
}

impl Drop for _SysMonStateGuard {
    fn drop(&mut self) {
        *self.state.lock().unwrap() = SysMonState::Stopped;
    }
}

fn monitor(
    pid: Option<u32>,
    worker_threads: Option<usize>,
    state: Arc<Mutex<SysMonState>>,
    need_stop_flag: Arc<AtomicBool>,
    monitor_interval: Duration,
    prog_name: Arc<String>,
) {
    let _sys_mon_state_guard: _SysMonStateGuard = state.into();
    let pid = pid.unwrap_or(std::process::id());

    let mut peak_tot_used_memory = 0;
    let mut peak_process_used_memory = 0;

    let mut now = Instant::now();
    let mut sys = sysinfo::System::new_all();

    loop {
        thread::sleep(Duration::from_secs(1));

        if now.elapsed() < monitor_interval {
            if need_stop_flag.load(Ordering::Relaxed) {
                break;
            }
            continue;
        }

        sys.refresh_all();
        std::thread::sleep(sysinfo::MINIMUM_CPU_UPDATE_INTERVAL);
        sys.refresh_all();

        now = Instant::now();
        let tot_mem = bytes2_gb(sys.total_memory() as usize);
        let tot_used_mem = bytes2_gb(sys.used_memory() as usize);

        peak_tot_used_memory = peak_tot_used_memory.max(tot_used_mem);
        tracing::info!(
            "System Memory Info ::> TotMem: {}GB, Used: {}GB, Free: {}GB. PeakUsed: {}GB",
            tot_mem,
            tot_used_mem,
            tot_mem - tot_used_mem,
            peak_tot_used_memory
        );
        let num_cpus = sys.cpus().len();
        let tot_cpu_usage: f32 = sys.cpus().iter().map(|cpu| cpu.cpu_usage()).sum();

        tracing::info!(
            "System Cpu Load Info ::> NumCores: {}, AllCoreUsage: {:.1}%, PerCoreUsage: {:.1}%",
            num_cpus,
            tot_cpu_usage,
            tot_cpu_usage / num_cpus as f32
        );

        if let Some(p) = sys.process(Pid::from_u32(pid)) {
            let worker_threads_ = worker_threads.unwrap_or(get_thread_count(pid).unwrap_or(0));

            let cpu_usage = p.cpu_usage();
            let process_mem = bytes2_gb(p.memory() as usize);
            peak_process_used_memory = peak_process_used_memory.max(process_mem);
            tracing::info!(
                "{} Memory Info ::> Used: {}GB, PeakUsed: {}GB",
                prog_name,
                process_mem,
                peak_process_used_memory
            );

            tracing::info!(
                "{} Cpu Load Info ::> CpuPercent: Tot:{:.1}%, NumThreads: {}, PerThread: {:.1}%",
                prog_name,
                cpu_usage,
                worker_threads_,
                if worker_threads_ > 0 {
                    cpu_usage / worker_threads_ as f32
                } else {
                    0.0
                }
            );
        }
    }
}

pub fn get_thread_count(pid: u32) -> std::io::Result<usize> {
    let procfs_path = format!("/proc/{}/task", pid);
    let mut count = 0;
    for entry in std::fs::read_dir(procfs_path)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            count += 1;
        }
    }
    Ok(count)
}
