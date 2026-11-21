//! Animated GIF writing and reading (the C++ science/GifWriter): a popularity
//! palette per frame on a 15-bit colour grid, LZW coding, looping.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use af_core::cxx::sort_by;

struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    accumulator: u32,
    bits: i32,
    block: Vec<u8>,
}

impl BitWriter<'_> {
    fn write(&mut self, code: i32, bits: i32) {
        self.accumulator |= (code as u32) << self.bits;
        self.bits += bits;
        while self.bits >= 8 {
            self.block.push((self.accumulator & 0xFF) as u8);
            self.accumulator >>= 8;
            self.bits -= 8;
            if self.block.len() == 255 {
                self.flush_block();
            }
        }
    }

    fn finish(&mut self) {
        if self.bits > 0 {
            self.block.push((self.accumulator & 0xFF) as u8);
        }
        self.bits = 0;
        self.accumulator = 0;
        self.flush_block();
        self.out.push(0);
    }

    fn flush_block(&mut self) {
        if self.block.is_empty() {
            return;
        }
        self.out.push(self.block.len() as u8);
        self.out.extend_from_slice(&self.block);
        self.block.clear();
    }
}

fn encode(out: &mut Vec<u8>, indices: &[u8]) {
    const CLEAR: i32 = 256;
    const END: i32 = 257;
    out.push(8); // minimum code size
    let mut writer = BitWriter { out, accumulator: 0, bits: 0, block: Vec::new() };
    let mut dictionary: HashMap<u32, i32> = HashMap::with_capacity(8192);
    let (mut code_size, mut next) = (9, 258);
    writer.write(CLEAR, code_size);
    if indices.is_empty() {
        writer.write(END, code_size);
        writer.finish();
        return;
    }
    let mut prefix = i32::from(indices[0]);
    for &c in &indices[1..] {
        let key = ((prefix as u32) << 8) | u32::from(c);
        if let Some(&found) = dictionary.get(&key) {
            prefix = found;
            continue;
        }
        writer.write(prefix, code_size);
        if next < 4096 {
            dictionary.insert(key, next);
            next += 1;
            // The decoder widens codes one entry later than it is added here.
            if next > (1 << code_size) && code_size < 12 {
                code_size += 1;
            }
        } else {
            writer.write(CLEAR, code_size);
            dictionary.clear();
            code_size = 9;
            next = 258;
        }
        prefix = i32::from(c);
    }
    writer.write(prefix, code_size);
    writer.write(END, code_size);
    writer.finish();
}

fn put16(out: &mut Vec<u8>, value: i32) {
    out.push((value & 0xFF) as u8);
    out.push(((value >> 8) & 0xFF) as u8);
}

/// Writes frames to an animated, looping GIF.
pub struct GifWriter {
    file: std::fs::File,
    width: i32,
    height: i32,
    delay: i32,
    frames: i32,
    closed: bool,
}

impl GifWriter {
    pub fn new(path: &Path, width: i32, height: i32, delay_centiseconds: i32) -> Result<Self, String> {
        let mut file = std::fs::File::create(path).map_err(|_| format!("Cannot write {}", path.display()))?;
        if !(1..=65535).contains(&width) || !(1..=65535).contains(&height) {
            return Err("GIF size must be 1-65535 pixels".into());
        }
        let mut head = b"GIF89a".to_vec();
        put16(&mut head, width);
        put16(&mut head, height);
        head.extend_from_slice(&[0, 0, 0]); // no global colour table
        // Loop forever (NETSCAPE2.0 application extension).
        head.extend_from_slice(&[0x21, 0xFF, 0x0B, b'N', b'E', b'T', b'S', b'C', b'A', b'P', b'E', b'2', b'.', b'0', 0x03, 0x01, 0x00, 0x00, 0x00]);
        file.write_all(&head).map_err(|e| e.to_string())?;
        Ok(Self { file, width, height, delay: delay_centiseconds, frames: 0, closed: false })
    }

    pub fn frames(&self) -> i32 {
        self.frames
    }

    /// RGBA (or RGB with channels = 3), rows top to bottom.
    pub fn add_frame(&mut self, pixels: &[u8], channels: usize) -> Result<(), String> {
        if self.closed {
            return Err("GIF already closed".into());
        }
        let count = self.width as usize * self.height as usize;
        if (channels != 3 && channels != 4) || pixels.len() < count * channels {
            return Err("GIF frame has the wrong size".into());
        }
        // Popularity palette on a 15-bit colour grid; each entry is its bin's mean colour.
        let mut histogram = vec![0u32; 32768];
        let mut sums = vec![[0u64; 3]; 32768];
        let mut keys = vec![0u16; count];
        for i in 0..count {
            let p = &pixels[i * channels..i * channels + 3];
            let key = ((u16::from(p[0] >> 3)) << 10) | ((u16::from(p[1] >> 3)) << 5) | u16::from(p[2] >> 3);
            keys[i] = key;
            histogram[key as usize] += 1;
            for k in 0..3 {
                sums[key as usize][k] += u64::from(p[k]);
            }
        }
        let mut bins: Vec<usize> = (0..32768).filter(|&k| histogram[k] > 0).collect();
        sort_by(&mut bins, |&a, &b| histogram[a] > histogram[b]);
        bins.truncate(256);
        let mut palette = [[0u8; 3]; 256];
        for (i, &bin) in bins.iter().enumerate() {
            for k in 0..3 {
                palette[i][k] = (sums[bin][k] / u64::from(histogram[bin])) as u8;
            }
        }
        let mut mapping = vec![-1i32; 32768];
        let mut indices = vec![0u8; count];
        for i in 0..count {
            let slot = &mut mapping[keys[i] as usize];
            if *slot < 0 {
                let p = &pixels[i * channels..i * channels + 3];
                let (mut best, mut best_distance) = (0usize, -1i64);
                for (c, colour) in palette.iter().enumerate().take(bins.len()) {
                    let distance: i64 = (0..3).map(|k| (i64::from(p[k]) - i64::from(colour[k])).pow(2)).sum();
                    if best_distance < 0 || distance < best_distance {
                        best_distance = distance;
                        best = c;
                    }
                }
                *slot = best as i32;
            }
            indices[i] = *slot as u8;
        }
        let mut out = Vec::new();
        // Graphic control extension: delay, no transparency.
        out.extend_from_slice(&[0x21, 0xF9, 0x04, 0x04]);
        put16(&mut out, self.delay);
        out.extend_from_slice(&[0, 0]);
        // Image descriptor with a 256-entry local colour table.
        out.push(0x2C);
        put16(&mut out, 0);
        put16(&mut out, 0);
        put16(&mut out, self.width);
        put16(&mut out, self.height);
        out.push(0x87);
        for colour in &palette {
            out.extend_from_slice(colour);
        }
        encode(&mut out, &indices);
        self.file.write_all(&out).map_err(|_| "Writing the GIF failed".to_string())?;
        self.frames += 1;
        Ok(())
    }

    pub fn close(&mut self) -> Result<(), String> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.file.write_all(&[0x3B]).map_err(|e| e.to_string())
    }
}

impl Drop for GifWriter {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

/// A decoded GIF: RGB frames, delays (centiseconds) and whether it loops.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GifImage {
    pub width: i32,
    pub height: i32,
    pub frames: Vec<Vec<u8>>,
    pub delays: Vec<i32>,
    pub looping: bool,
}

/// Reads a (non-interlaced) GIF.
pub fn read_gif(path: &Path) -> Result<GifImage, String> {
    let data = std::fs::read(path).unwrap_or_default();
    if data.len() < 13 || &data[..4] != b"GIF8" {
        return Err("Not a GIF file".into());
    }
    let at = |i: usize| -> Result<u8, String> { data.get(i).copied().ok_or_else(|| "Truncated GIF frame".to_string()) };
    let mut image = GifImage { width: i32::from(data[6]) | (i32::from(data[7]) << 8), height: i32::from(data[8]) | (i32::from(data[9]) << 8), ..Default::default() };
    let mut pos = 13usize;
    let mut global = Vec::new();
    if data[10] & 0x80 != 0 {
        let size = 2usize << (data[10] & 7);
        for i in 0..size {
            global.push([at(pos + 3 * i)?, at(pos + 3 * i + 1)?, at(pos + 3 * i + 2)?]);
        }
        pos += 3 * size;
    }
    let mut delay = 0;
    let sub_blocks = |pos: &mut usize, out: Option<&mut Vec<u8>>| {
        let mut out = out;
        while *pos < data.len() {
            let length = data[*pos] as usize;
            *pos += 1;
            if length == 0 {
                return;
            }
            if let Some(o) = out.as_deref_mut() {
                o.extend_from_slice(&data[*pos..(*pos + length).min(data.len())]);
            }
            *pos += length;
        }
    };
    while pos < data.len() {
        let tag = data[pos];
        pos += 1;
        if tag == 0x3B {
            break;
        }
        if tag == 0x21 {
            let label = at(pos)?;
            pos += 1;
            if label == 0xF9 {
                delay = i32::from(at(pos + 2)?) | (i32::from(at(pos + 3)?) << 8);
            }
            if label == 0xFF && pos + 12 <= data.len() && &data[pos + 1..pos + 12] == b"NETSCAPE2.0" {
                image.looping = true;
            }
            sub_blocks(&mut pos, None);
            continue;
        }
        if tag != 0x2C {
            return Err("Unexpected GIF block".into());
        }
        let w = i32::from(at(pos + 4)?) | (i32::from(at(pos + 5)?) << 8);
        let h = i32::from(at(pos + 6)?) | (i32::from(at(pos + 7)?) << 8);
        let flags = at(pos + 8)?;
        pos += 9;
        if flags & 0x40 != 0 {
            return Err("Interlaced GIF frames are not supported".into());
        }
        let mut palette = global.clone();
        if flags & 0x80 != 0 {
            palette.clear();
            let size = 2usize << (flags & 7);
            for i in 0..size {
                palette.push([at(pos + 3 * i)?, at(pos + 3 * i + 1)?, at(pos + 3 * i + 2)?]);
            }
            pos += 3 * size;
        }
        let minimum = i32::from(at(pos)?);
        pos += 1;
        let mut stream = Vec::new();
        sub_blocks(&mut pos, Some(&mut stream));
        // LZW decode.
        let clear = 1usize << minimum;
        let end = clear + 1;
        let reset = |table: &mut Vec<Vec<u8>>| {
            table.clear();
            table.resize(end + 1, Vec::new());
            for (i, entry) in table.iter_mut().enumerate().take(clear) {
                *entry = vec![i as u8];
            }
        };
        let mut table = Vec::new();
        reset(&mut table);
        let mut code_size = minimum + 1;
        let mut indices: Vec<u8> = Vec::new();
        let mut bit = 0usize;
        let mut previous: Option<usize> = None;
        while bit + code_size as usize <= stream.len() * 8 {
            let mut code = 0usize;
            for b in 0..code_size {
                code |= usize::from((stream[bit / 8] >> (bit % 8)) & 1) << b;
                bit += 1;
            }
            if code == clear {
                reset(&mut table);
                code_size = minimum + 1;
                previous = None;
                continue;
            }
            if code == end {
                break;
            }
            let entry = if code < table.len() && !table[code].is_empty() {
                table[code].clone()
            } else if let Some(p) = previous {
                let mut e = table[p].clone();
                e.push(e[0]);
                e
            } else {
                return Err("Corrupt GIF LZW stream".into());
            };
            indices.extend_from_slice(&entry);
            if let Some(p) = previous {
                if table.len() < 4096 {
                    let mut added = table[p].clone();
                    added.push(entry[0]);
                    table.push(added);
                    if table.len() >= (1usize << code_size) && code_size < 12 {
                        code_size += 1;
                    }
                }
            }
            previous = Some(code);
        }
        let n = w as usize * h as usize;
        if indices.len() < n {
            return Err("Truncated GIF frame".into());
        }
        let mut rgb = vec![0u8; n * 3];
        for i in 0..n {
            let colour = palette.get(indices[i] as usize).copied().unwrap_or([0; 3]);
            rgb[3 * i..3 * i + 3].copy_from_slice(&colour);
        }
        image.frames.push(rgb);
        image.delays.push(delay);
    }
    Ok(image)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        let dir = std::env::temp_dir().join("af_gif_round_trip");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.gif");
        let (w, h) = (37, 23);
        let mut writer = GifWriter::new(&path, w, h, 8).unwrap();
        for f in 0..3 {
            let pixels: Vec<u8> = (0..w * h).flat_map(|i| [((i * 7 + f * 40) % 256) as u8, ((i / 3) % 256) as u8, (f * 80) as u8, 255]).collect();
            writer.add_frame(&pixels, 4).unwrap();
        }
        writer.close().unwrap();
        let image = read_gif(&path).unwrap();
        assert_eq!((image.width, image.height, image.frames.len(), image.looping), (w, h, 3, true));
        assert_eq!(image.delays, vec![8, 8, 8]);
    }
}
