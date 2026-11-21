//! Scripted interface actions for UI tests and documentation captures
//! (ATOMFORGE_UI_SCRIPT=file; the C++ ui/UiScript and ui/TutorialCapture
//! drive the same kind of tour). One command per line, run one per frame:
//!
//!   open PATH            open a structure in a new tab
//!   step ID              open the operation dialog of a pipeline step
//!   set NAME TEXT        set one of its parameters (text as typed)
//!   target all|view|where [CONDITION]
//!   apply                press Apply (or Select)
//!   tool ID              open a scientific tool
//!   toolset NAME TEXT    set one of its parameters
//!   run                  run it and wait for the result
//!   pick X Y             click the viewport at these window pixels
//!   select-all | clear-selection | delete | undo | redo | reset-view
//!   orbit DX DY          drag the camera
//!   save PATH            save the structure
//!   tab N                show tab N (from 1)
//!   expect-atoms N       fail (exit code 3) unless the structure has N atoms
//!   expect-selected N    fail unless N atoms are selected
//!   expect-tabs N        fail unless N tabs are open
//!   screenshot PATH      save the window as PNG
//!   wait N               let N frames pass
//!   quit
//!
//! Lines starting with # are comments.

pub struct UiScript {
    lines: Vec<String>,
    next: usize,
    pub wait_frames: u32,
    pub failed: Option<String>,
}

impl UiScript {
    pub fn from_env() -> Option<Self> {
        let path = std::env::var("ATOMFORGE_UI_SCRIPT").ok()?;
        let text = std::fs::read_to_string(&path).unwrap_or_else(|e| format!("fail cannot read {path}: {e}"));
        let lines = text.lines().map(str::trim).filter(|l| !l.is_empty() && !l.starts_with('#')).map(str::to_string).collect();
        Some(Self { lines, next: 0, wait_frames: 3, failed: None })
    }

    /// The next command, if no wait is pending.
    pub fn next_command(&mut self) -> Option<String> {
        if self.wait_frames > 0 {
            self.wait_frames -= 1;
            return None;
        }
        let line = self.lines.get(self.next).cloned()?;
        self.next += 1;
        Some(line)
    }

    /// Runs the last command again (after any wait): commands that wait for
    /// background work, such as wait-cloud.
    pub fn repeat_last(&mut self) {
        self.next = self.next.saturating_sub(1);
    }

    pub fn finished(&self) -> bool {
        self.next >= self.lines.len()
    }
}

/// The command word and the rest of a script line.
pub fn split(line: &str) -> (&str, &str) {
    match line.split_once(char::is_whitespace) {
        Some((head, rest)) => (head, rest.trim()),
        None => (line, ""),
    }
}
