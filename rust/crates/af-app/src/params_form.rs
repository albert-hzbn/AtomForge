//! A pipeline step's parameters as a form (the C++ ui/ModifierEditor
//! drawModifierParameters): a two-column table, labels on the left (42 %) and
//! fields filling the right: a check box for booleans, a number field with
//! - and + for integers, three fields for vectors, a list for choices, a
//! text field (parsed as the step's text syntax) for everything else, and for
//! Build steps one row per command-line flag of the builder's usage.

use af_core::json::Json;
use af_pipeline::pipeline::{text_for, value_for, ModifierType};
use eframe::egui;

/// The texts of the fields being edited (kept between frames, so a value can
/// be typed through states that do not parse yet).
#[derive(Default)]
pub struct FormTexts {
    texts: Vec<(String, String)>,
}

impl FormTexts {
    pub fn clear(&mut self) {
        self.texts.clear();
    }

    /// Sets a field's text (UI scripts); the parsed value, if it parses.
    pub fn set(&mut self, t: &ModifierType, name: &str, text: &str) -> Result<Json, String> {
        let p = t.find_parameter(name).ok_or_else(|| format!("{} has no parameter {name}", t.id))?;
        let value = value_for(p, text)?;
        match self.texts.iter_mut().find(|(n, _)| n == name) {
            Some(slot) => slot.1 = text.to_string(),
            None => self.texts.push((name.to_string(), text.to_string())),
        }
        Ok(value)
    }

    fn slot(&mut self, name: &str, initial: impl FnOnce() -> String) -> usize {
        match self.texts.iter().position(|(n, _)| n == name) {
            Some(i) => i,
            None => {
                self.texts.push((name.to_string(), initial()));
                self.texts.len() - 1
            }
        }
    }
}

/// The command-line flags of a builder, from its help text (modifierFlags).
pub mod flags {
    #[derive(Clone, Debug, PartialEq)]
    pub struct Flag {
        /// Without "--".
        pub name: String,
        /// Placeholder or default ("" for a switch).
        pub value: String,
        pub optional: bool,
        pub repeatable: bool,
    }

    #[derive(Clone, Debug, Default, PartialEq)]
    pub struct Usage {
        pub flags: Vec<Flag>,
        pub description: String,
    }

    /// Splits a command line at spaces, keeping "quoted words" together (unquoted).
    pub fn tokens(text: &str) -> Vec<String> {
        let mut out = Vec::new();
        let (mut word, mut quoted, mut any) = (String::new(), false, false);
        for c in text.chars() {
            if c == '"' {
                quoted = !quoted;
                any = true;
            } else if !quoted && (c == ' ' || c == '\t' || c == '\n') {
                if any {
                    out.push(std::mem::take(&mut word));
                }
                word.clear();
                any = false;
            } else {
                word.push(c);
                any = true;
            }
        }
        if any {
            out.push(word);
        }
        out
    }

    fn quote(value: &str) -> String {
        if value.contains([' ', '\t']) {
            format!("\"{value}\"")
        } else {
            value.to_string()
        }
    }

    fn trim(text: &str) -> &str {
        text.trim_matches([' ', '\t'])
    }

    /// [preamble], "TITLE (--build x)", usage lines (starting with "--" or
    /// "["), description. Without a title line there is no usage.
    pub fn parse(help: &str) -> Usage {
        let mut usage = Usage::default();
        let Some(title) = help.find("(--") else { return usage };
        let rest = match help[title..].find('\n') {
            Some(n) => &help[title + n + 1..],
            None => "",
        };
        let mut usage_text = String::new();
        let mut in_usage = true;
        for line in rest.lines() {
            let text = trim(line);
            if text.is_empty() && usage.description.is_empty() {
                continue;
            }
            if in_usage && (text.starts_with("--") || text.starts_with('[')) {
                usage_text.push_str(text);
                usage_text.push(' ');
                continue;
            }
            in_usage = false;
            if !usage.description.is_empty() {
                usage.description.push(' ');
            }
            usage.description.push_str(text);
        }
        let b = usage_text.as_bytes();
        let mut depth = 0i32;
        let mut i = 0;
        while i < b.len() {
            let c = b[i];
            if c == b'[' {
                depth += 1;
                i += 1;
                continue;
            }
            if c == b']' {
                depth -= 1;
                i += 1;
                continue;
            }
            if !usage_text[i..].starts_with("--") {
                i += 1;
                continue;
            }
            let mut j = i + 2;
            while j < b.len() && (b[j].is_ascii_alphanumeric() || b[j] == b'-') {
                j += 1;
            }
            let mut flag = Flag { name: usage_text[i + 2..j].to_string(), value: String::new(), optional: depth > 0, repeatable: false };
            let mut k = j;
            while k < b.len() && b[k] == b' ' {
                k += 1;
            }
            if k < b.len() && b[k] == b'"' {
                let end = usage_text[k + 1..].find('"').map(|e| e + k + 1).unwrap_or(b.len());
                flag.value = usage_text[k + 1..end.min(b.len())].to_string();
                k = (end + 1).min(b.len());
            } else if k < b.len() && !usage_text[k..].starts_with("--") && !b"[]|(".contains(&b[k]) {
                let mut end = k;
                while end < b.len() && !b" ]|".contains(&b[end]) {
                    end += 1;
                }
                flag.value = usage_text[k..end].to_string();
                k = end;
            }
            let next = usage_text[k..].find("--").map(|n| n + k).unwrap_or(b.len());
            flag.repeatable = usage_text[k..next].contains("repeatable");
            if flag.name != "input" && flag.name != "output" {
                usage.flags.push(flag);
            }
            i = k;
        }
        usage
    }

    /// Each flag's value in `options` ("1" for a given switch; repeated
    /// values joined by "; "), and the words that are no flag of the usage.
    pub fn values(usage: &Usage, options: &str) -> (Vec<String>, String) {
        let mut out = vec![String::new(); usage.flags.len()];
        let mut other = String::new();
        let words = tokens(options);
        let mut i = 0;
        while i < words.len() {
            let f = words[i].strip_prefix("--").and_then(|name| usage.flags.iter().position(|f| f.name == name));
            let Some(f) = f else {
                if !other.is_empty() {
                    other.push(' ');
                }
                other.push_str(&quote(&words[i]));
                i += 1;
                continue;
            };
            if usage.flags[f].value.is_empty() {
                out[f] = "1".into();
                i += 1;
                continue;
            }
            let value = if i + 1 < words.len() {
                i += 1;
                words[i].clone()
            } else {
                String::new()
            };
            out[f] = if usage.flags[f].repeatable && !out[f].is_empty() { format!("{}; {value}", out[f]) } else { value };
            i += 1;
        }
        (out, other)
    }

    pub fn compose(usage: &Usage, values: &[String], other: &str) -> String {
        let mut out = String::new();
        let mut add = |word: &str| {
            if !out.is_empty() {
                out.push(' ');
            }
            out.push_str(word);
        };
        for (flag, value) in usage.flags.iter().zip(values) {
            if flag.value.is_empty() {
                if value == "1" {
                    add(&format!("--{}", flag.name));
                }
                continue;
            }
            let entries: Vec<&str> = if flag.repeatable { value.split(';').map(trim).collect() } else { vec![trim(value)] };
            for entry in entries {
                if !entry.is_empty() {
                    add(&format!("--{} {}", flag.name, quote(entry)));
                }
            }
        }
        if !trim(other).is_empty() {
            add(trim(other));
        }
        out
    }

    /// "min-separation" -> "Min separation"; single letters stay as they are.
    pub fn label(name: &str) -> String {
        match name {
            "symprec" => return "Symmetry tolerance".into(),
            "nmax" => return "Max Miller index".into(),
            "no-primitive" => return "Keep conventional cell".into(),
            _ => {}
        }
        let label = name.replace('-', " ");
        if label.chars().count() > 1 {
            let mut c = label.chars();
            let first = c.next().unwrap().to_uppercase().collect::<String>();
            return first + c.as_str();
        }
        label
    }

    /// The hint of a flag's field: "default 4", "number", the placeholder.
    pub fn hint(flag: &Flag) -> String {
        let numeric = flag.value.parse::<f64>().is_ok();
        let mut hint = if numeric {
            format!("default {}", flag.value)
        } else if flag.value.chars().count() == 1 {
            "number".into()
        } else {
            flag.value.clone()
        };
        hint = hint.to_lowercase();
        if flag.repeatable {
            hint.push_str("; ...");
        }
        if !flag.optional {
            hint.push_str(" (required)");
        }
        hint
    }
}

/// The description of a step: the builder's description for Build steps
/// whose usage could be read, else the step's help (modifierDescription).
pub fn description(t: &ModifierType) -> String {
    if t.parameters.iter().any(|p| p.kind == "options") {
        let usage = flags::parse(&t.help);
        if !usage.flags.is_empty() {
            return usage.description;
        }
    }
    t.help.clone()
}

/// A row: the label in the left column (wrapped), the field on the right.
fn row<R>(ui: &mut egui::Ui, label: &str, tip: &str, field: impl FnOnce(&mut egui::Ui, f32) -> R) -> R {
    let width = ui.available_width();
    let label_width = (width * 0.42).floor();
    let spacing = ui.spacing().item_spacing.x;
    let field_width = (width - label_width - spacing).max(40.0);
    ui.horizontal(|ui| {
        ui.allocate_ui_with_layout(egui::vec2(label_width, crate::theme::FRAME_HEIGHT), egui::Layout::left_to_right(egui::Align::Center).with_main_wrap(true), |ui| {
            ui.set_min_width(label_width);
            let response = ui.add(egui::Label::new(crate::plot_widget::unit_display(label)).wrap());
            if !tip.is_empty() {
                response.on_hover_text(tip);
            }
        });
        field(ui, field_width)
    })
    .inner
}

/// Draws the form; returns the parameters changed this frame.
pub fn parameter_form(ui: &mut egui::Ui, id: &str, t: &ModifierType, parameters: &Json, texts: &mut FormTexts) -> Vec<(String, Json)> {
    let mut changed = Vec::new();
    ui.push_id(id, |ui| {
        ui.spacing_mut().item_spacing.y = 6.0;
        for p in &t.parameters {
            let current = parameters.find(&p.name).cloned().unwrap_or(Json::Null);
            let kind = p.kind.as_str();
            if kind == "options" {
                let usage = flags::parse(&t.help);
                if !usage.flags.is_empty() {
                    let options = current.string().unwrap_or("").to_string();
                    let (mut values, mut other) = flags::values(&usage, &options);
                    let mut edited = false;
                    for (k, flag) in usage.flags.iter().enumerate() {
                        ui.push_id(("flag", k), |ui| {
                            row(ui, &flags::label(&flag.name), &format!("--{} {}", flag.name, flag.value), |ui, w| {
                                if flag.value.is_empty() {
                                    let mut on = values[k] == "1";
                                    if crate::dialog_widgets::checkbox(ui, &mut on, "").changed() {
                                        values[k] = if on { "1".into() } else { String::new() };
                                        edited = true;
                                    }
                                } else if ui.add(egui::TextEdit::singleline(&mut values[k]).hint_text(flags::hint(flag)).desired_width(w)).changed() {
                                    edited = true;
                                }
                            });
                        });
                    }
                    row(ui, "", "", |_, _| {});
                    let muted = crate::theme::palette(ui.visuals().dark_mode).text_muted;
                    ui.horizontal(|ui| {
                        let width = ui.available_width();
                        let lw = (width * 0.42).floor();
                        ui.allocate_ui_with_layout(egui::vec2(lw, crate::theme::FRAME_HEIGHT), egui::Layout::left_to_right(egui::Align::Center), |ui| {
                            ui.set_min_width(lw);
                            ui.label(egui::RichText::new("Other flags").color(muted));
                        });
                        let w = (width - lw - ui.spacing().item_spacing.x).max(40.0);
                        if ui.add(egui::TextEdit::singleline(&mut other).hint_text("optional").desired_width(w)).changed() {
                            edited = true;
                        }
                    });
                    if edited {
                        changed.push((p.name.clone(), Json::from(flags::compose(&usage, &values, &other))));
                    }
                    continue;
                }
            }
            ui.push_id(&p.name, |ui| {
                row(ui, &p.label, &p.name, |ui, width| match kind {
                    "bool" => {
                        let mut on = current.boolean().unwrap_or(false);
                        if crate::dialog_widgets::checkbox(ui, &mut on, "").changed() {
                            changed.push((p.name.clone(), Json::from(on)));
                        }
                    }
                    "choice" => {
                        let mut value = current.string().unwrap_or("").to_string();
                        egui::ComboBox::from_id_salt("value").selected_text(value.clone()).width(width - 8.0).show_ui(ui, |ui| {
                            for option in p.options.split('|') {
                                ui.selectable_value(&mut value, option.to_string(), option);
                            }
                        });
                        if Some(value.as_str()) != current.string().ok() {
                            changed.push((p.name.clone(), Json::from(value)));
                        }
                    }
                    "vector" | "ivec3" => {
                        // Three fields, as ImGui::InputScalarN.
                        let slot = texts.slot(&p.name, || text_for(p, &current).unwrap_or_default());
                        let mut parts: Vec<String> = texts.texts[slot].1.replace(',', " ").split_whitespace().map(str::to_string).collect();
                        parts.resize(3, String::new());
                        let each = ((width - 2.0 * ui.spacing().item_spacing.x) / 3.0).max(30.0);
                        let mut edited = false;
                        let mut last = None;
                        for (k, part) in parts.iter_mut().enumerate() {
                            let r = ui.add(egui::TextEdit::singleline(part).id_salt(("component", k)).desired_width(each - 8.0));
                            edited |= r.changed();
                            last = Some(r);
                        }
                        let text = parts.join(" ");
                        texts.texts[slot].1 = text.clone();
                        let parsed = value_for(p, &text);
                        if edited {
                            if let Ok(value) = &parsed {
                                changed.push((p.name.clone(), value.clone()));
                            }
                        }
                        if let (Err(e), Some(r)) = (parsed, last) {
                            r.on_hover_text(e);
                        }
                    }
                    _ => {
                        let slot = texts.slot(&p.name, || text_for(p, &current).unwrap_or_default());
                        let int = kind == "int";
                        let field_width = if int { width - 2.0 * (crate::theme::FRAME_HEIGHT + ui.spacing().item_spacing.x) } else { width };
                        let hint = match kind {
                            "expression" => "e.g. fz > 0.5 && element == O",
                            "options" => "command-line flags, e.g. --count 4",
                            "elements" => "e.g. Cu Ni",
                            "list" => "words separated by spaces",
                            _ => "",
                        };
                        let response = ui.add(egui::TextEdit::singleline(&mut texts.texts[slot].1).hint_text(hint).desired_width(field_width - 8.0));
                        let mut edited = response.changed();
                        // ImGui::InputInt's step buttons.
                        if int {
                            let step = |ui: &mut egui::Ui, label: &str| ui.add_sized([crate::theme::FRAME_HEIGHT, crate::theme::FRAME_HEIGHT], egui::Button::new(label)).clicked();
                            let value = texts.texts[slot].1.trim().parse::<f64>().unwrap_or(0.0);
                            if step(ui, "-") {
                                texts.texts[slot].1 = format!("{}", value - 1.0);
                                edited = true;
                            }
                            if step(ui, "+") {
                                texts.texts[slot].1 = format!("{}", value + 1.0);
                                edited = true;
                            }
                        }
                        let parsed = value_for(p, &texts.texts[slot].1);
                        if edited {
                            if let Ok(value) = &parsed {
                                changed.push((p.name.clone(), value.clone()));
                            }
                        }
                        if let Err(e) = parsed {
                            response.on_hover_text(e);
                        }
                    }
                });
            });
        }
    });
    changed
}

#[cfg(test)]
mod tests {
    use super::flags::*;

    #[test]
    fn builder_usage_round_trip() {
        let help = "Builds a slab.\nSurface Slab (--build surface)\n  --hkl \"1 1 1\" --layers 4 [--vacuum 10] [--center]\n  [--atom \"Cu 0 0 0\" (repeatable)]\nCuts a slab from a crystal.";
        let usage = parse(help);
        let names: Vec<&str> = usage.flags.iter().map(|f| f.name.as_str()).collect();
        assert_eq!(names, ["hkl", "layers", "vacuum", "center", "atom"]);
        assert_eq!(usage.flags[0].value, "1 1 1");
        assert!(!usage.flags[0].optional && usage.flags[2].optional);
        assert!(usage.flags[4].repeatable);
        assert_eq!(usage.description, "Cuts a slab from a crystal.");
        let (values, other) = values(&usage, "--hkl \"1 0 0\" --center --atom \"Cu 0 0 0\" --atom \"Ni 1 1 1\" --extra 2");
        assert_eq!(values[0], "1 0 0");
        assert_eq!(values[3], "1");
        assert_eq!(values[4], "Cu 0 0 0; Ni 1 1 1");
        assert_eq!(other, "--extra 2");
        assert_eq!(compose(&usage, &values, &other), "--hkl \"1 0 0\" --center --atom \"Cu 0 0 0\" --atom \"Ni 1 1 1\" --extra 2");
        assert_eq!(label("min-separation"), "Min separation");
        assert_eq!(hint(&usage.flags[1]), "default 4 (required)");
    }
}
