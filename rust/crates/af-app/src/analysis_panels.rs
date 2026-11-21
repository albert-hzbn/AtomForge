//! Layout pieces shared by the self-contained analysis dialogs (LOBSTER,
//! Wannier, electronic post-processing, dislocation post-processing): the
//! C++ ui/ResponsiveLayout controls/results split, ui/DialogLayout sections
//! and empty states, and ui/LinePlot's multi-series plot.

use af_science::plots::{PlotSeries, PlotSpec};
use eframe::egui;

/// A controls sidebar beside the results (stacked below 1000 points, as
/// responsive::stacked(1000)); the sidebar is 42% of the width, 420–540.
pub fn two_panels<S>(ui: &mut egui::Ui, id: &str, state: &mut S, controls: impl FnOnce(&mut S, &mut egui::Ui), results: impl FnOnce(&mut S, &mut egui::Ui)) {
    let width = ui.available_width();
    let stacked = width < 1000.0;
    if stacked {
        egui::ScrollArea::vertical().id_salt((id, "controls")).max_height(420.0).show(ui, |ui| {
            ui.set_width(ui.available_width());
            controls(state, ui);
        });
        ui.separator();
        egui::ScrollArea::vertical().id_salt((id, "results")).show(ui, |ui| {
            ui.set_width(ui.available_width());
            results(state, ui);
        });
        return;
    }
    let sidebar = (width * 0.42).clamp(420.0, 540.0);
    let height = ui.available_height().max(300.0);
    let top_down = egui::Layout::top_down(egui::Align::Min);
    ui.horizontal_top(|ui| {
        let p = crate::theme::palette(ui.visuals().dark_mode);
        ui.allocate_ui_with_layout(egui::vec2(sidebar, height), top_down, |ui| {
            // ImGui child windows: the surface colour, rounded, the controls bordered.
            egui::Frame::new().fill(p.surface).stroke(egui::Stroke::new(1.0, p.border)).corner_radius(8).inner_margin(egui::Margin::same(7)).show(ui, |ui| {
                ui.set_width(sidebar - 14.0);
                ui.set_height(height - 14.0);
                egui::ScrollArea::vertical().id_salt((id, "controls")).auto_shrink([false, false]).show(ui, |ui| {
                    ui.set_width(sidebar - 30.0);
                    controls(state, ui);
                });
            });
        });
        let rest = ui.available_width();
        ui.allocate_ui_with_layout(egui::vec2(rest, height), top_down, |ui| {
            egui::Frame::new().fill(p.surface).corner_radius(8).inner_margin(egui::Margin::same(7)).show(ui, |ui| {
                ui.set_height(height - 14.0);
                egui::ScrollArea::vertical().id_salt((id, "results")).auto_shrink([false, false]).show(ui, |ui| {
                    ui.set_width(rest - 26.0);
                    results(state, ui);
                });
            });
        });
    });
}

/// dialogLayout::section: a semibold heading over a thin rule.
pub fn section(ui: &mut egui::Ui, title: &str) {
    crate::dialog_widgets::section(ui, title);
}

/// dialogLayout::emptyState: a muted hint centred in the results.
pub fn empty_state(ui: &mut egui::Ui, text: &str) {
    let height = ui.available_height().max(120.0);
    crate::dialog_widgets::empty_state(ui, text, height);
}

/// The ImGui item width of a field in a panel: a fraction of the width left.
pub fn item_width(ui: &egui::Ui, fraction: f32) -> f32 {
    (ui.available_width() * fraction).max(60.0)
}

/// A label left of its control (the C++ field() rows: 40% label, 60% value).
pub fn field<R>(ui: &mut egui::Ui, label: &str, control: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let width = ui.available_width();
    ui.horizontal(|ui| {
        ui.add_sized([width * 0.40 - 6.0, 18.0], egui::Label::new(label).wrap());
        ui.spacing_mut().interact_size.x = width * 0.60 - 6.0;
        ui.scope(|ui| {
            ui.set_width(width * 0.60 - 6.0);
            control(ui)
        })
        .inner
    })
    .inner
}

/// The uiPlot::drawColoredLines plot: several series over a shared x axis.
pub fn line_plot(ui: &mut egui::Ui, title: &str, x_label: &str, y_label: &str, x: &[f64], series: Vec<(String, Vec<f64>)>, height: f32) {
    let spec = PlotSpec {
        title: title.to_string(),
        x_label: x_label.to_string(),
        y_label: y_label.to_string(),
        series: series.into_iter().map(|(name, y)| PlotSeries { name, x: x.to_vec(), y, points: false }).collect(),
        ..Default::default()
    };
    crate::plot_widget::science_plot(ui, &spec, height);
}

/// The file name of a path, for "Loaded: NAME" lines.
pub fn file_name(path: &str) -> String {
    std::path::Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| path.to_string())
}

/// C++ ostream default formatting (%g, 6 significant digits).
pub fn g(v: f64) -> String {
    af_core::format::format_g(v, 6)
}
