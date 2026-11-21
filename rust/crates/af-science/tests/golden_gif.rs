//! af-science's GIF writer against the C++ GifWriter: the frames of
//! tools/golden/golden_gif.cpp, compared byte for byte, and read back.

use std::path::PathBuf;

use af_science::gif::{read_gif, GifWriter};

/// The frames of golden_gif.cpp.
fn frame(width: usize, height: usize, index: u32, channels: usize) -> Vec<u8> {
    let mut pixels = vec![0u8; width * height * channels];
    let mut state: u32 = 12345u32.wrapping_add(index.wrapping_mul(977));
    for y in 0..height {
        for x in 0..width {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            let p = &mut pixels[(y * width + x) * channels..];
            p[0] = ((x * 255) / (width - 1)) as u8;
            p[1] = ((y * 255) / (height - 1)) as u8;
            p[2] = ((state >> 24) & if index == 2 { 0x0F } else { 0xFF }) as u8;
            if channels == 4 {
                p[3] = 255;
            }
        }
    }
    pixels
}

#[test]
fn gifs_match_cpp() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden");
    let dir = std::env::temp_dir().join("af_golden_gif");
    std::fs::create_dir_all(&dir).unwrap();
    let mut gif = GifWriter::new(&dir.join("many_colours.gif"), 64, 48, 8).unwrap();
    for i in 0..3 {
        gif.add_frame(&frame(64, 48, i, 4), 4).unwrap();
    }
    gif.close().unwrap();
    let mut rgb = GifWriter::new(&dir.join("rgb.gif"), 17, 5, 3).unwrap();
    rgb.add_frame(&frame(17, 5, 2, 3), 3).unwrap();
    rgb.close().unwrap();
    for name in ["many_colours.gif", "rgb.gif"] {
        let mine = std::fs::read(dir.join(name)).unwrap();
        let theirs = std::fs::read(root.join(name)).unwrap();
        assert!(mine == theirs, "{name}: {} vs {} bytes, first difference at {:?}", mine.len(), theirs.len(), mine.iter().zip(&theirs).position(|(a, b)| a != b));
        let image = read_gif(&root.join(name)).unwrap();
        assert!(image.looping && !image.frames.is_empty());
    }
}
