"""Non-destructive structure pipelines from Python.

A pipeline is an ordered list of modifiers applied to a copy of a structure,
the same pipelines as the desktop panel (Edit > Structure pipeline) and
``AtomForge --pipe``. The input structure is never changed.

>>> import atomforge as af
>>> from atomforge.pipeline import Pipeline
>>> s = af.load("cu.cif")
>>> p = (Pipeline()
...      .replicate(4, 4, 4)
...      .select_expression("fz > 0.5")
...      .assign_element("Ni")
...      .build_vacancy("--count 3 --seed 1"))
>>> result = p.run(s)
>>> result.structure, result.selected_count, [st["atoms"] for st in result.stages]
>>> p.move(3, 0); p.disable(1); p.remove(2)      # edit, reorder, bypass, delete
>>> Pipeline.parse("replicate 2 2 2 | select-random 0.1 | delete-selected").run(s)
>>> s.pipe("replicate 2 2 1 | wrap")              # the same as a one-off call

Steps may be written as methods (``select_expression``), with ``add("select-expression", ...)``,
or as text; positional values follow the modifier's parameter order, keywords
name them, and ``|`` joins pipelines and text. ``modifiers()`` lists every step
with its parameters.

Core modifiers (selection, deletion, replication, slicing, strain, analysis,
...) run in the native library ``atomforge_science_native``. Build and Edit
steps (``build-*``: bulk, surface, vacancy, solid solution, SQS, dislocation,
...) use the AtomForge executable, found like :func:`atomforge.view` (or set
``ATOMFORGE_PATH``); their option text is the command-line flags of
``AtomForge --build MODE``.
"""

import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

from ._structure import Atom, Structure
from . import native

__all__ = ["Pipeline", "PipelineResult", "PipelineError", "Step", "modifiers"]


class PipelineError(ValueError):
    """A pipeline step failed; ``stage`` is its index and ``stages`` the per-step results."""

    def __init__(self, message: str, stage: Optional[int] = None, stages: Optional[List[Dict[str, Any]]] = None):
        super().__init__(message)
        self.stage = stage
        self.stages = stages or []


# ---------------------------------------------------------------- engine access
def _call(function: str, *arguments: bytes) -> str:
    lib = native._load()
    status = native.ct.c_int(1)
    f = getattr(lib, function)
    f.restype = native.ct.c_void_p
    f.argtypes = [native.ct.c_char_p] * len(arguments) + [native.ct.POINTER(native.ct.c_int)]
    text = native._take(lib, f(*arguments, native.ct.byref(status)))
    if status.value != 0:
        raise PipelineError(text)
    return text


_native_types: Optional[List[Dict[str, Any]]] = None
_all_types: Optional[List[Dict[str, Any]]] = None


def _native_modifiers() -> List[Dict[str, Any]]:
    global _native_types
    if _native_types is None:
        lib = native._load()
        lib.afs_pipeline_modifiers.restype = native.ct.c_void_p
        _native_types = json.loads(native._take(lib, lib.afs_pipeline_modifiers()))
    return _native_types


def _executable() -> str:
    explicit = os.environ.get("ATOMFORGE_PATH")
    if explicit:
        return str(Path(explicit).resolve())
    from ._viewer import _find_atomforge
    found = _find_atomforge()
    if not found:
        raise FileNotFoundError("Build and Edit steps need the AtomForge executable; set ATOMFORGE_PATH to it")
    return str(found)


def modifiers() -> List[Dict[str, Any]]:
    """Every pipeline step: ``[{"id", "title", "category", "help", "parameters": [...]}]``.

    Build and Edit steps are included when the AtomForge executable is available."""
    global _all_types
    if _all_types is None:
        types = list(_native_modifiers())
        try:
            listed = subprocess.run([_executable(), "--pipe", "--list-json"], capture_output=True, text=True, timeout=60)
            if listed.returncode == 0:
                known = {t["id"] for t in types}
                types += [t for t in json.loads(listed.stdout) if t["id"] not in known]
        except (FileNotFoundError, OSError, ValueError):
            pass
        _all_types = types
    return _all_types


def _is_native(step_type: str) -> bool:
    return any(t["id"] == step_type for t in _native_modifiers())


# ---------------------------------------------------------------- steps and pipelines
class Step:
    """One pipeline step: modifier ``type``, ``parameters`` (dict), ``enabled`` and an optional ``label``."""

    def __init__(self, type: str, parameters: Optional[Dict[str, Any]] = None, enabled: bool = True, label: str = ""):
        self.type = type
        self.parameters = dict(parameters or {})
        self.enabled = enabled
        self.label = label

    def to_json(self) -> Dict[str, Any]:
        """The step as stored in pipeline JSON."""
        item: Dict[str, Any] = {"type": self.type, "enabled": self.enabled, "parameters": self.parameters}
        if self.label:
            item["label"] = self.label
        return item

    @classmethod
    def from_json(cls, item: Dict[str, Any]) -> "Step":
        """A step from pipeline JSON."""
        return cls(item["type"], item.get("parameters", {}), item.get("enabled", True), item.get("label", ""))

    def __repr__(self) -> str:
        flag = "" if self.enabled else ", disabled"
        return f"Step({self.type!r}, {self.parameters!r}{flag})"


def _values(value: Any) -> Any:
    if hasattr(value, "tolist"):
        return value.tolist()
    if isinstance(value, tuple):
        return list(value)
    return value


def _quote(word: str) -> str:
    if word and not any(c.isspace() for c in word) and '"' not in word and "'" not in word:
        return word
    mark = "'" if '"' in word else '"'
    return mark + word + mark


def _make_step(step_type: str, args: Sequence[Any], params: Dict[str, Any]) -> Step:
    """A step from positional and keyword values, with the same rules as the text syntax."""
    params = {k: _values(v) for k, v in params.items()}
    if step_type.startswith("build-"):
        # Build and Edit steps take their command-line flags as one option text.
        if args:
            params["options"] = " ".join(str(a) for a in args)
        return Step(step_type, params)
    words = [step_type]
    for a in args:
        a = _values(a)
        if isinstance(a, list):
            words.append(",".join(str(x) for x in a))
        elif isinstance(a, bool):
            words.append("true" if a else "false")
        else:
            words.append(_quote(str(a)))
    parsed = json.loads(_call("afs_pipeline_parse", " ".join(words).encode("utf-8")))["modifiers"][0]
    step = Step.from_json(parsed)
    for name, value in params.items():
        if not any(p["name"] == name for t in _native_modifiers() if t["id"] == step_type for p in t["parameters"]):
            raise PipelineError(f"{step_type} has no parameter {name}")
        step.parameters[name] = value
    return step


class Pipeline:
    """An ordered, editable list of steps; ``run(structure)`` applies them to a copy."""

    def __init__(self, steps: Optional[Iterable[Union[Step, Dict[str, Any]]]] = None):
        self.steps: List[Step] = [s if isinstance(s, Step) else Step.from_json(s) for s in (steps or [])]

    # ---- construction
    @classmethod
    def parse(cls, text: str) -> "Pipeline":
        """A pipeline from the shell-pipe text syntax, e.g. ``"replicate 2 2 2 | wrap"``."""
        pipeline = cls()
        for part in _split(text):
            words = part.strip()
            if not words:
                continue
            step_type, _, rest = words.partition(" ")
            if step_type.startswith("build-"):
                pipeline.steps.append(Step(step_type, {"options": rest.strip()} if rest.strip() else {}))
            else:
                pipeline.steps.append(Step.from_json(json.loads(_call("afs_pipeline_parse", words.encode("utf-8")))["modifiers"][0]))
        return pipeline

    @classmethod
    def from_json(cls, data: Union[str, Dict[str, Any]]) -> "Pipeline":
        """A pipeline from its JSON form (a dict or JSON text), as saved by the desktop and the CLI."""
        if isinstance(data, str):
            data = json.loads(data)
        return cls(data.get("modifiers", []))

    @classmethod
    def load(cls, path: Union[str, os.PathLike]) -> "Pipeline":
        """Load a pipeline file (JSON, or pipe text)."""
        text = Path(path).read_text(encoding="utf-8")
        return cls.from_json(text) if text.lstrip().startswith("{") else cls.parse(text)

    def add(self, step_type: str, *args: Any, **params: Any) -> "Pipeline":
        """Append a step; returns the pipeline, so calls chain."""
        self.steps.append(_make_step(step_type, args, params))
        return self

    def insert(self, index: int, step_type: str, *args: Any, **params: Any) -> "Pipeline":
        """Insert a step at ``index``."""
        self.steps.insert(index, _make_step(step_type, args, params))
        return self

    # ---- editing
    def remove(self, index: int) -> "Pipeline":
        """Delete the step at ``index``."""
        del self.steps[index]
        return self

    def move(self, source: int, destination: int) -> "Pipeline":
        """Move the step at ``source`` so that it ends up at ``destination``."""
        step = self.steps.pop(source)
        self.steps.insert(destination, step)
        return self

    def enable(self, index: int, enabled: bool = True) -> "Pipeline":
        """Apply (or with ``enabled=False`` bypass) the step at ``index``."""
        self.steps[index].enabled = enabled
        return self

    def disable(self, index: int) -> "Pipeline":
        """Bypass the step at ``index`` without removing it."""
        return self.enable(index, False)

    def set(self, index: int, **params: Any) -> "Pipeline":
        """Change parameters of the step at ``index``."""
        self.steps[index].parameters.update({k: _values(v) for k, v in params.items()})
        return self

    def copy(self) -> "Pipeline":
        """An independent copy."""
        return Pipeline.from_json(self.to_json())

    def __len__(self) -> int:
        return len(self.steps)

    def __getitem__(self, index: int) -> Step:
        return self.steps[index]

    def __iter__(self):
        return iter(self.steps)

    def __or__(self, other: Union[str, "Pipeline", Step]) -> "Pipeline":
        joined = self.copy()
        if isinstance(other, str):
            other = Pipeline.parse(other)
        joined.steps += [other] if isinstance(other, Step) else other.copy().steps
        return joined

    def __ror__(self, other: str) -> "Pipeline":
        return Pipeline.parse(other) | self

    def __getattr__(self, name: str):
        # p.select_expression("fz > 0.5") == p.add("select-expression", "fz > 0.5")
        if name.startswith("_"):
            raise AttributeError(name)
        step_type = name.replace("_", "-")
        if not step_type.startswith("build-") and not _is_native(step_type):
            raise AttributeError(f"No pipeline modifier {step_type!r}")

        def add(*args: Any, **params: Any) -> "Pipeline":
            return self.add(step_type, *args, **params)

        add.__name__ = name
        return add

    # ---- serialization
    def to_json(self) -> Dict[str, Any]:
        """The JSON form used by the desktop panel, project files and ``AtomForge --pipeline``."""
        return {"format": "atomforge-pipeline", "version": 1, "modifiers": [s.to_json() for s in self.steps]}

    def to_text(self) -> str:
        """The shell-pipe text of the enabled steps (``AtomForge --pipe`` accepts it)."""
        parts = []
        for step in self.steps:
            if not step.enabled:
                continue
            if step.type.startswith("build-"):
                parts.append((step.type + " " + str(step.parameters.get("options", ""))).strip())
            else:
                parts.append(_call("afs_pipeline_text", json.dumps({"modifiers": [step.to_json()]}).encode("utf-8")))
        return " | ".join(parts)

    def save(self, path: Union[str, os.PathLike]) -> None:
        """Save as JSON (open it in the desktop panel with Load pipeline...)."""
        Path(path).write_text(json.dumps(self.to_json(), indent=2) + "\n", encoding="utf-8")

    def __repr__(self) -> str:
        return f"Pipeline({self.to_text()!r})" if self.steps else "Pipeline()"

    # ---- running
    def run(self, structure: Structure, property: Optional[Sequence[float]] = None, check: bool = True) -> "PipelineResult":
        """Apply the enabled steps to a copy of ``structure`` (which is not changed).

        ``property`` gives per-atom values for ``select-property``. A failing
        step raises PipelineError, or with ``check=False`` the result holds the
        data before that step and the step's error."""
        if all(_is_native(s.type) for s in self.steps if s.enabled):
            request: Dict[str, Any] = {"pipeline": self.to_json(), "structure": native._structure_json(structure)}
            if property is not None:
                request["property"] = [float(v) for v in _values(property)]
            document = json.loads(_call("afs_pipeline_run", json.dumps(request).encode("utf-8")))
            result = PipelineResult(native._from_json(document["structure"]), document["selected"],
                                    document.get("property") or [], document.get("property_name", ""), document["stages"])
        else:
            result = self._run_executable(structure, property)
        for index, stage in enumerate(result.stages):
            if stage.get("error") and check:
                raise PipelineError(f"step {index + 1} ({stage['type']}): {stage['error']}", index, result.stages)
        return result

    def _run_executable(self, structure: Structure, property: Optional[Sequence[float]]) -> "PipelineResult":
        with tempfile.TemporaryDirectory(prefix="atomforge_pipe_") as folder:
            steps = Path(folder) / "pipeline.json"
            report = Path(folder) / "report.json"
            steps.write_text(json.dumps(self.to_json()), encoding="utf-8")
            done = subprocess.run([_executable(), "--pipeline", str(steps), "--input", "-", "--output", "-", "--quiet", "--report", str(report)],
                                  input=_extxyz(structure, property), capture_output=True, text=True, timeout=3600)
            stages = json.loads(report.read_text(encoding="utf-8"))["stages"] if report.exists() else []
            if done.returncode != 0:
                if stages and any(s.get("error") for s in stages):
                    out = _parse_extxyz("0\n\n") if not done.stdout.strip() else _parse_extxyz(done.stdout)
                    return PipelineResult(out[0], out[1], out[2], "", stages)
                message = done.stderr.strip()
                if message.startswith("Error: "):
                    message = message[len("Error: "):]
                raise PipelineError(message or "AtomForge --pipe failed", None, stages)
            output, selected, values = _parse_extxyz(done.stdout)
            name = json.loads(report.read_text(encoding="utf-8")).get("property_name", "") if report.exists() else ""
            return PipelineResult(output, selected, values, name, stages)


class PipelineResult:
    """Output of :meth:`Pipeline.run`: ``structure``, ``selected`` (one bool per atom),
    ``property`` (per-atom values or empty), ``property_name`` and ``stages``
    (per step: type, atoms, selected, milliseconds, skipped, error, notes)."""

    def __init__(self, structure: Structure, selected: Sequence[Any], property: Sequence[Any], property_name: str, stages: List[Dict[str, Any]]):
        self.structure = structure
        self.selected = [bool(s) for s in selected]
        self.property = [float("nan") if v is None else float(v) for v in property]
        self.property_name = property_name
        self.stages = stages

    @property
    def selected_count(self) -> int:
        """Number of selected atoms in the output."""
        return sum(self.selected)

    @property
    def selected_indices(self) -> List[int]:
        """Indices of the selected atoms in the output."""
        return [i for i, s in enumerate(self.selected) if s]

    def __repr__(self) -> str:
        return f"PipelineResult({self.structure!r}, {self.selected_count} selected, {len(self.stages)} steps)"


# ---------------------------------------------------------------- helpers
def _split(text: str) -> List[str]:
    """Split pipe text on | outside quotes (|| is the logical operator)."""
    parts, current, quote, i = [""], "", "", 0
    while i < len(text):
        c = text[i]
        if quote:
            if c == quote:
                quote = ""
        elif c in "\"'":
            quote = c
        elif c == "|":
            if i + 1 < len(text) and text[i + 1] == "|":
                parts[-1] += "||"
                i += 2
                continue
            parts.append("")
            i += 1
            continue
        parts[-1] += c
        i += 1
    return parts


def _extxyz(structure: Structure, property: Optional[Sequence[float]]) -> str:
    lines = [str(len(structure.atoms))]
    header = ""
    if structure.cell:
        header += 'Lattice="' + " ".join(repr(float(v)) for row in structure.cell for v in row) + '" pbc="T T T" '
    columns = "species:S:1:pos:R:3" + (":property:R:1" if property is not None else "")
    lines.append(header + "Properties=" + columns)
    values = list(_values(property)) if property is not None else None
    for i, a in enumerate(structure.atoms):
        row = f"{a.symbol} {a.x!r} {a.y!r} {a.z!r}"
        if values is not None:
            row += f" {float(values[i])!r}"
        lines.append(row)
    return "\n".join(lines) + "\n"


def _parse_extxyz(text: str):
    lines = text.splitlines()
    count = int(lines[0].split()[0]) if lines and lines[0].strip() else 0
    comment = lines[1] if len(lines) > 1 else ""
    structure = Structure()
    if 'Lattice="' in comment:
        numbers = [float(v) for v in comment.split('Lattice="')[1].split('"')[0].split()]
        structure.cell = [numbers[0:3], numbers[3:6], numbers[6:9]]
    has_property = ":property:R:1" in comment
    selected, values = [], []
    for row in lines[2:2 + count]:
        parts = row.split()
        structure.atoms.append(Atom(parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        selected.append(parts[4] == "1" if len(parts) > 4 else False)
        if has_property and len(parts) > 5:
            values.append(float(parts[5]))
    return structure, selected, values
