//! Result plots against the C++ science/ResultPlots: the plots, CSV text and
//! overlays that tools/golden/golden_plots.cpp made from the C++ tool results
//! (tests/golden/tools.json), and histogram edge cases, compared exactly.

use std::path::PathBuf;

use af_core::json::Json;
use af_science::plots::{histogram, overlay_plots, plot_csv, result_plots, PlotSeries, PlotSpec};

fn golden(name: &str) -> Json {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden").join(name);
    Json::parse(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn numbers(j: &Json) -> Vec<f64> {
    j.items()
        .unwrap()
        .iter()
        .map(|v| match v {
            Json::Number(n) => *n,
            Json::String(s) if s == "nan" => f64::NAN,
            Json::String(s) if s == "inf" => f64::INFINITY,
            Json::String(s) if s == "-inf" => f64::NEG_INFINITY,
            other => panic!("not a number: {other:?}"),
        })
        .collect()
}

fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()))
}

fn check_series(got: &PlotSeries, want: &Json, what: &str) {
    assert_eq!(got.name, want.at("name").unwrap().string().unwrap(), "{what}: name");
    assert!(same(&got.x, &numbers(want.at("x").unwrap())), "{what} {}: x\n{:?}\n{:?}", got.name, got.x, numbers(want.at("x").unwrap()));
    assert!(same(&got.y, &numbers(want.at("y").unwrap())), "{what} {}: y", got.name);
    assert_eq!(got.points, want.at("points").unwrap().boolean().unwrap(), "{what}: points");
}

fn check_plot(got: &PlotSpec, want: &Json, what: &str) {
    assert_eq!(got.title, want.at("title").unwrap().string().unwrap(), "{what}: title");
    assert_eq!(got.x_label, want.at("x_label").unwrap().string().unwrap(), "{what}");
    assert_eq!(got.y_label, want.at("y_label").unwrap().string().unwrap(), "{what}");
    let series = want.at("series").unwrap().items().unwrap();
    assert_eq!(got.series.len(), series.len(), "{what} {}: series", got.title);
    for (s, w) in got.series.iter().zip(series) {
        check_series(s, w, &format!("{what} {}", got.title));
    }
    let markers = want.at("markers").unwrap().items().unwrap();
    assert_eq!(got.markers.len(), markers.len(), "{what}: markers");
    for (m, w) in got.markers.iter().zip(markers) {
        let w = w.items().unwrap();
        assert!(m.x.to_bits() == w[0].number().unwrap().to_bits() && m.label == w[1].string().unwrap(), "{what}: marker {m:?}");
    }
    assert!(same(&got.horizontal, &numbers(want.at("horizontal").unwrap())), "{what}: horizontal");
    assert_eq!(got.log_y, want.at("log_y").unwrap().boolean().unwrap(), "{what}: log y");
    assert_eq!(plot_csv(got), want.at("csv").unwrap().string().unwrap(), "{what} {}: csv", got.title);
}

#[test]
fn plots_match_cpp() {
    let tools = golden("tools.json");
    let plots = golden("plots.json");
    let wanted = plots.at("tools").unwrap().items().unwrap();
    let results = tools.items().unwrap();
    assert_eq!(results.len(), wanted.len());
    let mut earlier: Vec<(String, Vec<PlotSpec>)> = Vec::new();
    let mut compared = 0;
    for (k, (item, want)) in results.iter().zip(wanted).enumerate() {
        let tool = item.at("tool").unwrap().string().unwrap();
        let Some(result) = item.find("result") else { continue };
        let what = format!("case {k} {tool}");
        let got = result_plots(tool, result);
        if let Some(e) = want.find("error") {
            assert!(got.is_err(), "{what}: the C++ failed ({})", e.string().unwrap());
            continue;
        }
        let got = got.unwrap_or_else(|e| panic!("{what}: {e}"));
        let list = want.at("plots").unwrap().items().unwrap();
        assert_eq!(got.len(), list.len(), "{what}: plot count");
        for (p, w) in got.iter().zip(list) {
            check_plot(p, w, &what);
            compared += 1;
        }
        if let Some(w) = want.find("overlay") {
            let same_tool: Vec<(String, Vec<PlotSpec>)> = earlier.iter().filter(|(label, _)| label.starts_with(&format!("{tool}#"))).cloned().collect();
            check_plot(&overlay_plots(&got[0], "now", &same_tool), w, &format!("{what} overlay"));
        }
        earlier.push((format!("{tool}#{}", earlier.len()), got));
    }
    assert!(compared > 150, "{compared} plots compared");
    let mut h = plots.at("histograms").unwrap().items().unwrap().iter();
    let nan = f64::NAN;
    let inf = f64::INFINITY;
    for values in [vec![], vec![nan, inf], vec![2.5], vec![-1.0, -1.0, -1.0], (0..=10).map(f64::from).collect(), vec![1e-9, 3e-9, nan, 2e-9], vec![-5.0, 0.25, 7.75, 1e6]] {
        for bins in [1, 3, 40] {
            check_series(&histogram("h", &values, bins), h.next().unwrap(), &format!("histogram {values:?} {bins}"));
        }
    }
}
