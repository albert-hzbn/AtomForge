//! Cancellation and progress for long computations (the C++ util/TaskControl
//! and util/BackgroundTask): a worker runs with a control in scope; library
//! code calls `checkpoint()`, which fails once the task is cancelled.

use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::Arc;
use std::thread::JoinHandle;

pub const CANCELLED: &str = "Calculation cancelled";

#[derive(Debug)]
pub struct TaskControl {
    cancelled: AtomicBool,
    /// f32 bits; negative means no estimate.
    progress: AtomicU32,
}

impl Default for TaskControl {
    fn default() -> Self {
        Self { cancelled: AtomicBool::new(false), progress: AtomicU32::new((-1.0f32).to_bits()) }
    }
}

impl TaskControl {
    pub fn cancel(&self) {
        self.cancelled.store(true, Ordering::Relaxed);
    }
    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
    pub fn progress(&self) -> f32 {
        f32::from_bits(self.progress.load(Ordering::Relaxed))
    }
}

thread_local! {
    static ACTIVE: RefCell<Option<Arc<TaskControl>>> = const { RefCell::new(None) };
}

/// Makes `control` the active control of this thread until dropped.
pub struct Scope {
    previous: Option<Arc<TaskControl>>,
}

impl Scope {
    pub fn new(control: Arc<TaskControl>) -> Self {
        let previous = ACTIVE.with(|a| a.replace(Some(control)));
        Self { previous }
    }
}

impl Drop for Scope {
    fn drop(&mut self) {
        let previous = self.previous.take();
        ACTIVE.with(|a| *a.borrow_mut() = previous);
    }
}

/// Fails with "Calculation cancelled" once the active task is cancelled.
pub fn checkpoint() -> Result<(), String> {
    let cancelled = ACTIVE.with(|a| a.borrow().as_ref().is_some_and(|c| c.is_cancelled()));
    if cancelled {
        Err(CANCELLED.into())
    } else {
        Ok(())
    }
}

/// Reports progress (0-1) of the active task, after a checkpoint.
pub fn progress(fraction: f64) -> Result<(), String> {
    checkpoint()?;
    ACTIVE.with(|a| {
        if let Some(c) = a.borrow().as_ref() {
            c.progress.store((fraction.clamp(0.0, 1.0) as f32).to_bits(), Ordering::Relaxed);
        }
    });
    Ok(())
}

/// A computation on a worker thread, owned and polled by the UI thread.
pub struct BackgroundTask<R: Send + 'static> {
    handle: Option<JoinHandle<Result<R, String>>>,
    control: Option<Arc<TaskControl>>,
    result: Option<R>,
    error: String,
}

impl<R: Send + 'static> Default for BackgroundTask<R> {
    fn default() -> Self {
        Self { handle: None, control: None, result: None, error: String::new() }
    }
}

impl<R: Send + 'static> BackgroundTask<R> {
    pub fn running(&self) -> bool {
        self.handle.is_some()
    }
    pub fn result(&self) -> Option<&R> {
        self.result.as_ref()
    }
    pub fn take_result(&mut self) -> Option<R> {
        self.result.take()
    }
    pub fn error(&self) -> &str {
        &self.error
    }
    pub fn cancel(&self) {
        if let Some(c) = &self.control {
            c.cancel();
        }
    }
    pub fn progress(&self) -> f32 {
        self.control.as_ref().map_or(-1.0, |c| c.progress())
    }
    pub fn clear_result(&mut self) {
        self.result = None;
        self.error.clear();
    }

    /// Starts the work unless a task is running.
    pub fn start(&mut self, work: impl FnOnce() -> Result<R, String> + Send + 'static) -> bool {
        if self.running() {
            return false;
        }
        self.clear_result();
        let control = Arc::new(TaskControl::default());
        self.control = Some(control.clone());
        self.handle = Some(std::thread::spawn(move || {
            let _scope = Scope::new(control);
            checkpoint()?;
            let result = work()?;
            checkpoint()?;
            progress(1.0)?;
            Ok(result)
        }));
        true
    }

    /// True exactly once, when the result or error has been collected.
    pub fn poll(&mut self) -> bool {
        if !self.handle.as_ref().is_some_and(|h| h.is_finished()) {
            return false;
        }
        match self.handle.take().unwrap().join() {
            Ok(Ok(result)) => self.result = Some(result),
            Ok(Err(error)) => self.error = error,
            Err(_) => self.error = "Analysis failed with an unknown exception.".into(),
        }
        true
    }

    /// Waits for the task (for tests and the command line).
    pub fn wait(&mut self) {
        while self.running() && !self.poll() {
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
    }
}

impl<R: Send + 'static> Drop for BackgroundTask<R> {
    fn drop(&mut self) {
        self.cancel();
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cancellation_reaches_checkpoints() {
        let mut task: BackgroundTask<i32> = BackgroundTask::default();
        assert!(task.start(|| {
            for _ in 0..10_000 {
                checkpoint()?;
                std::thread::sleep(std::time::Duration::from_micros(100));
            }
            Ok(1)
        }));
        task.cancel();
        task.wait();
        assert_eq!(task.error(), CANCELLED);
        let mut done: BackgroundTask<i32> = BackgroundTask::default();
        done.start(|| Ok(7));
        done.wait();
        assert_eq!(done.result(), Some(&7));
        assert_eq!(done.progress(), 1.0);
        assert!(checkpoint().is_ok(), "no task on this thread");
    }
}
