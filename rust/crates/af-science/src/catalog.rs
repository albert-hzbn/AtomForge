//! Desktop and CLI metadata for the native science tools (the C++
//! science/ScienceCatalog). Parameter names match the JSON request keys of the
//! atomforge.science Python API.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ToolView {
    Atoms,
    Trajectory,
    Properties,
    Plot,
    Simulation,
    Generator,
}

#[derive(Debug)]
pub struct ParameterDef {
    pub name: &'static str,
    pub label: &'static str,
    pub value: &'static str,
    pub kind: &'static str,
    pub required: bool,
}

#[derive(Debug)]
pub struct ToolDef {
    pub id: &'static str,
    pub title: &'static str,
    /// First line: one-line description; the rest: method notes.
    pub help: &'static str,
    /// Analysis menu section.
    pub category: &'static str,
    pub parameters: &'static [ParameterDef],
    pub view: ToolView,
}

pub fn catalog() -> &'static [ToolDef] {
    super::catalog_table::TOOLS
}

/// The catalog entry with this id.
pub fn find_tool(id: &str) -> Option<&'static ToolDef> {
    catalog().iter().find(|t| t.id == id)
}
