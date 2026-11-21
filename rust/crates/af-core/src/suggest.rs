//! "Did you mean ...?" for a misspelt name (an option, tool, step or
//! parameter): the closest known name, if close enough (the C++ util/Suggest).

/// Edit distance: insertions, deletions, substitutions and swapped
/// neighbours ("widht") each count one (bytes, as the C++).
pub fn distance(a: &str, b: &str) -> usize {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut d = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i;
    }
    for j in 0..=b.len() {
        d[0][j] = j;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            d[i][j] = (d[i - 1][j] + 1).min(d[i][j - 1] + 1).min(d[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                d[i][j] = d[i][j].min(d[i - 2][j - 2] + 1);
            }
        }
    }
    d[a.len()][b.len()]
}

/// The closest candidate when at most about a third of the name's letters
/// differ (ignoring a leading "--"), else "".
pub fn closest<S: AsRef<str>>(name: &str, candidates: impl IntoIterator<Item = S>) -> String {
    let letters = if name.starts_with("--") { name.len() - 2 } else { name.len() };
    let mut best = (letters / 3).max(1) + 1;
    let mut found = String::new();
    for candidate in candidates {
        let candidate = candidate.as_ref();
        let d = distance(name, candidate);
        if d < best {
            best = d;
            found = candidate.to_string();
        }
    }
    found
}

/// " (did you mean X?)" or "".
pub fn did_you_mean<S: AsRef<str>>(name: &str, candidates: impl IntoIterator<Item = S>) -> String {
    let found = closest(name, candidates);
    if found.is_empty() {
        String::new()
    } else {
        format!(" (did you mean {found}?)")
    }
}

/// "a, b, c".
pub fn list<S: AsRef<str>>(names: impl IntoIterator<Item = S>) -> String {
    names.into_iter().map(|n| n.as_ref().to_string()).collect::<Vec<_>>().join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn suggestions() {
        assert_eq!(distance("widht", "width"), 1);
        assert_eq!(did_you_mean("replicat", ["replicate", "wrap"]), " (did you mean replicate?)");
        assert_eq!(did_you_mean("xyz", ["replicate"]), "");
        assert_eq!(closest("--widht", ["--width"]), "--width");
        assert_eq!(list(["a", "b"]), "a, b");
    }
}
