//! Undo and redo (the C++ app/UndoRedo): whole snapshots of the structure and
//! the element tables, described by what changed ("delete 3 atoms").

use af_core::structure::Structure;

#[derive(Clone, Debug, Default)]
pub struct Snapshot {
    pub structure: Structure,
    pub element_radii: Vec<f32>,
    pub element_colours: Vec<[f32; 3]>,
    pub element_shininess: Vec<f32>,
}

impl PartialEq for Snapshot {
    fn eq(&self, other: &Self) -> bool {
        let (a, b) = (&self.structure, &other.structure);
        if a.has_unit_cell != b.has_unit_cell || a.cell_offset != b.cell_offset || a.cell_vectors != b.cell_vectors || a.atoms.len() != b.atoms.len() {
            return false;
        }
        let same_atoms = a.atoms.iter().zip(&b.atoms).all(|(p, q)| {
            p.symbol == q.symbol && p.atomic_number == q.atomic_number && p.x == q.x && p.y == q.y && p.z == q.z && p.r == q.r && p.g == q.g && p.b == q.b
        });
        same_atoms && self.element_radii == other.element_radii && self.element_colours == other.element_colours && self.element_shininess == other.element_shininess
    }
}

/// What changed from `before` to `after`, as a short phrase for the Undo and
/// Redo commands ("delete 3 atoms").
pub fn describe_change(before: &Snapshot, after: &Snapshot) -> String {
    let (a, b) = (before.structure.atoms.len(), after.structure.atoms.len());
    let atoms = |n: usize| format!("{n} {}", if n == 1 { "atom" } else { "atoms" });
    if a == 0 && b > 0 {
        return "load the structure".into();
    }
    if a != b {
        // Adding or deleting keeps the other atoms, in order; anything else
        // (another frame, a builder's result) replaces the structure.
        let (small, large) = if a < b { (&before.structure.atoms, &after.structure.atoms) } else { (&after.structure.atoms, &before.structure.atoms) };
        let mut matched = 0;
        for q in large {
            if matched >= small.len() {
                break;
            }
            let p = &small[matched];
            if p.symbol == q.symbol && p.x == q.x && p.y == q.y && p.z == q.z {
                matched += 1;
            }
        }
        if matched < small.len() {
            return format!("replace the structure ({})", atoms(b));
        }
        return if b > a { format!("add {}", atoms(b - a)) } else { format!("delete {}", atoms(a - b)) };
    }
    let (mut moved, mut changed) = (0, 0);
    for (p, q) in before.structure.atoms.iter().zip(&after.structure.atoms) {
        if p.x != q.x || p.y != q.y || p.z != q.z {
            moved += 1;
        }
        if p.symbol != q.symbol || p.atomic_number != q.atomic_number {
            changed += 1;
        }
    }
    let cell = before.structure.has_unit_cell != after.structure.has_unit_cell || before.structure.cell_vectors != after.structure.cell_vectors;
    if cell {
        "change the cell".into()
    } else if changed > 0 {
        format!("change the element of {}", atoms(changed))
    } else if moved > 0 {
        format!("move {}", atoms(moved))
    } else {
        "change display settings".into()
    }
}

/// The history of snapshots with the current one at `index`.
pub struct UndoRedo {
    history: Vec<Snapshot>,
    ids: Vec<u64>,
    descriptions: Vec<String>,
    index: usize,
    next_id: u64,
    max_bytes: usize,
}

const MAX_ENTRIES: usize = 200;

impl Default for UndoRedo {
    fn default() -> Self {
        Self { history: Vec::new(), ids: Vec::new(), descriptions: Vec::new(), index: 0, next_id: 0, max_bytes: 1536 * 1024 * 1024 }
    }
}

impl UndoRedo {
    pub fn reset(&mut self, initial: Snapshot) {
        self.next_id += 1;
        self.history = vec![initial];
        self.ids = vec![self.next_id];
        self.descriptions = vec![String::new()];
        self.index = 0;
    }

    pub fn commit(&mut self, snapshot: Snapshot) {
        if self.history.is_empty() || self.index >= self.history.len() {
            self.reset(snapshot);
            return;
        }
        // Nothing changed (a view update): no step to undo, and the document stays as saved.
        if self.history[self.index] == snapshot {
            return;
        }
        self.history.truncate(self.index + 1);
        self.ids.truncate(self.index + 1);
        self.descriptions.truncate(self.index + 1);
        // Display-only changes are undoable but keep the id: they are not saved in structure files.
        let description = describe_change(self.history.last().expect("history"), &snapshot);
        let display_only = description == "change display settings";
        self.descriptions.push(description);
        self.history.push(snapshot);
        let id = if display_only {
            *self.ids.last().expect("ids")
        } else {
            self.next_id += 1;
            self.next_id
        };
        self.ids.push(id);
        self.index = self.history.len() - 1;
        let mut bytes: usize = self.history.iter().map(snapshot_bytes).sum();
        while self.history.len() > 2 && (self.history.len() > MAX_ENTRIES || bytes > self.max_bytes) {
            bytes -= snapshot_bytes(&self.history[0]);
            self.history.remove(0);
            self.ids.remove(0);
            self.descriptions.remove(0);
            self.index -= 1;
        }
    }

    pub fn can_undo(&self) -> bool {
        !self.history.is_empty() && self.index > 0
    }

    pub fn can_redo(&self) -> bool {
        !self.history.is_empty() && self.index + 1 < self.history.len()
    }

    pub fn undo(&mut self) -> Option<&Snapshot> {
        if !self.can_undo() {
            return None;
        }
        self.index -= 1;
        self.history.get(self.index)
    }

    pub fn redo(&mut self) -> Option<&Snapshot> {
        if !self.can_redo() {
            return None;
        }
        self.index += 1;
        self.history.get(self.index)
    }

    pub fn jump_to(&mut self, i: usize) -> Option<&Snapshot> {
        if self.history.is_empty() {
            return None;
        }
        self.index = i.min(self.history.len() - 1);
        self.history.get(self.index)
    }

    pub fn len(&self) -> usize {
        self.history.len()
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn description(&self, i: usize) -> &str {
        self.descriptions.get(i).map(String::as_str).unwrap_or("")
    }

    pub fn undo_description(&self) -> String {
        if self.can_undo() { self.descriptions[self.index].clone() } else { String::new() }
    }

    pub fn redo_description(&self) -> String {
        if self.can_redo() { self.descriptions[self.index + 1].clone() } else { String::new() }
    }

    /// Identifies the current state (unchanged by undo/redo back to it).
    pub fn current_id(&self) -> u64 {
        self.ids.get(self.index.min(self.ids.len().saturating_sub(1))).copied().unwrap_or(0)
    }
}

fn snapshot_bytes(s: &Snapshot) -> usize {
    let st = &s.structure;
    std::mem::size_of::<Snapshot>()
        + st.atoms.len() * std::mem::size_of::<af_core::structure::AtomSite>()
        + st.grain_region_ids.len() * 4
        + st.grain_colors.len() * 12
        + (s.element_radii.len() + s.element_shininess.len()) * 4
        + s.element_colours.len() * 12
}

#[cfg(test)]
mod tests {
    use super::*;
    use af_core::structure::AtomSite;

    fn snap(xs: &[f64]) -> Snapshot {
        let mut s = Snapshot::default();
        for &x in xs {
            s.structure.atoms.push(AtomSite { symbol: "Cu".into(), atomic_number: 29, x, ..Default::default() });
        }
        s
    }

    #[test]
    fn history_descriptions_and_ids() {
        let mut u = UndoRedo::default();
        u.reset(snap(&[]));
        u.commit(snap(&[0.0, 1.0]));
        assert_eq!(u.undo_description(), "load the structure");
        let saved = u.current_id();
        u.commit(snap(&[0.0]));
        assert_eq!(u.undo_description(), "delete 1 atom");
        u.commit(snap(&[0.5]));
        assert_eq!(u.undo_description(), "move 1 atom");
        u.undo();
        u.undo();
        assert_eq!(u.current_id(), saved);
        assert_eq!(u.redo_description(), "delete 1 atom");
        u.commit(snap(&[0.0, 1.0, 2.0]));
        assert!(!u.can_redo());
        assert_eq!(u.undo_description(), "add 1 atom");
    }
}
