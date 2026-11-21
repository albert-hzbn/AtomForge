//! Settings > Renderer (the C++ graphics/GraphicsBackend preference): the wgpu
//! backend to draw with, saved as `renderer=` in ~/.atomforge_preferences.txt
//! and applied at the next start. The environment variable ATOMFORGE_RENDERER
//! overrides it, and --renderer on the command line overrides both.

use std::path::PathBuf;

/// The backends offered (name in the file, menu label).
pub const CHOICES: [(&str, &str); 5] = [("auto", "Automatic (default)"), ("vulkan", "Vulkan"), ("dx12", "Direct3D 12"), ("metal", "Metal"), ("opengl", "OpenGL")];

/// Whether a backend can exist on this platform (menu entries).
pub fn available(name: &str) -> bool {
    match name {
        "dx12" => cfg!(windows),
        "metal" => cfg!(target_os = "macos"),
        "vulkan" => !cfg!(target_os = "macos"),
        _ => true,
    }
}

pub fn preferences_file() -> PathBuf {
    let home = if cfg!(windows) { std::env::var_os("USERPROFILE") } else { std::env::var_os("HOME") };
    home.filter(|h| !h.is_empty()).map(PathBuf::from).unwrap_or_else(|| PathBuf::from(".")).join(".atomforge_preferences.txt")
}

/// A renderer name as the C++ accepted it ("vk", "gl" too), normalised.
pub fn parse(text: &str) -> Option<&'static str> {
    match text.trim().to_ascii_lowercase().as_str() {
        "auto" | "automatic" | "default" => Some("auto"),
        "vulkan" | "vk" => Some("vulkan"),
        "opengl" | "gl" | "gles" => Some("opengl"),
        "dx12" | "d3d12" | "direct3d" => Some("dx12"),
        "metal" => Some("metal"),
        _ => None,
    }
}

/// The saved choice ("auto" when none).
pub fn saved() -> &'static str {
    std::fs::read_to_string(preferences_file())
        .ok()
        .and_then(|text| text.lines().filter_map(|l| l.strip_prefix("renderer=")).filter_map(parse).last())
        .unwrap_or("auto")
}

/// Rewrites the preferences file keeping its other settings.
pub fn save(name: &str) -> std::io::Result<()> {
    let path = preferences_file();
    let mut lines: Vec<String> = std::fs::read_to_string(&path).map(|t| t.lines().filter(|l| !l.starts_with("renderer=")).map(str::to_string).collect()).unwrap_or_default();
    lines.push(format!("renderer={name}"));
    std::fs::write(path, lines.join("\n") + "\n")
}

/// The choice for this run: the command line, else ATOMFORGE_RENDERER, else the saved one.
pub fn requested(command_line: Option<&str>) -> &'static str {
    let mut chosen = saved();
    if let Some(e) = std::env::var("ATOMFORGE_RENDERER").ok().as_deref().and_then(parse) {
        chosen = e;
    }
    if let Some(c) = command_line.and_then(parse) {
        chosen = c;
    }
    chosen
}

pub fn backends(name: &str) -> Option<eframe::wgpu::Backends> {
    match name {
        "vulkan" => Some(eframe::wgpu::Backends::VULKAN),
        "opengl" => Some(eframe::wgpu::Backends::GL),
        "dx12" => Some(eframe::wgpu::Backends::DX12),
        "metal" => Some(eframe::wgpu::Backends::METAL),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names() {
        assert_eq!(parse("VK"), Some("vulkan"));
        assert_eq!(parse("gl"), Some("opengl"));
        assert_eq!(parse("d3d12"), Some("dx12"));
        assert_eq!(parse("software"), None);
        assert!(available("opengl") && available("auto"));
    }
}
