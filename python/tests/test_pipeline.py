"""atomforge.pipeline: non-destructive structure pipelines from Python."""

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import atomforge as af
from atomforge.pipeline import Pipeline, PipelineError, Step, modifiers


def copper(repeat=1, a=3.615):
    s = af.Structure()
    for i in range(repeat):
        for j in range(repeat):
            for k in range(repeat):
                for f in ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)):
                    s.add_atom("Cu", (i + f[0]) * a, (j + f[1]) * a, (k + f[2]) * a)
    s.cell = [[a * repeat, 0, 0], [0, a * repeat, 0], [0, 0, a * repeat]]
    return s


class PipelineTests(unittest.TestCase):
    def test_build_edit_and_run(self):
        s = copper()
        p = Pipeline().replicate(2, 2, 2).select_expression("fz > 0.5").assign_element("Ni")
        result = p.run(s)
        self.assertEqual(len(result.structure), 32)
        self.assertEqual(result.selected_count, 16)
        self.assertEqual(sum(a.symbol == "Ni" for a in result.structure.atoms), 16)
        self.assertEqual([st["atoms"] for st in result.stages], [32, 32, 32])
        self.assertEqual(len(s), 4, "the input is unchanged")
        # Edit the steps: move, bypass, change, remove, insert.
        p.add("delete-selected")
        self.assertEqual(len(p.run(s).structure), 16)
        p.disable(3)
        self.assertEqual(len(p.run(s).structure), 32)
        p.enable(3).set(1, expression="fz > 0.4")
        self.assertEqual(len(p.run(s).structure), 16)
        p.move(0, 3)  # replicate last: the selection is made on the single cell
        self.assertEqual(p[3].type, "replicate")
        self.assertEqual(len(p.run(s).structure), 16)
        p.remove(3).insert(0, "replicate", counts=(1, 1, 2))
        self.assertEqual(len(p.run(s).structure), 4)

    def test_text_json_and_operators(self):
        text = "replicate 2 2 2 | select-random 0.25 seed=3 | delete-selected"
        p = Pipeline.parse(text)
        self.assertEqual([st.type for st in p], ["replicate", "select-random", "delete-selected"])
        self.assertEqual(Pipeline.parse(p.to_text()).to_json(), p.to_json())
        self.assertEqual(Pipeline.from_json(json.dumps(p.to_json())).to_json(), p.to_json())
        joined = Pipeline.parse("replicate 2 1 1") | "wrap" | Pipeline().center()
        self.assertEqual([st.type for st in joined], ["replicate", "wrap", "center"])
        self.assertEqual(len(copper().pipe(text)), 24)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "steps.json"
            p.save(path)
            self.assertEqual(Pipeline.load(path).to_json(), p.to_json())
        with self.assertRaises(PipelineError):
            Pipeline.parse("replicate 2 2")
        with self.assertRaises(AttributeError):
            Pipeline().not_a_modifier
        with self.assertRaises(PipelineError):
            Pipeline().add("replicate", bogus=1)

    def test_errors_and_property(self):
        s = copper(2)
        with self.assertRaises(PipelineError) as failure:
            Pipeline().select_property(0, 1).run(s)
        self.assertEqual(failure.exception.stage, 0)
        partial = Pipeline().replicate(2, 1, 1).select_property(0, 1).run(s, check=False)
        self.assertEqual(len(partial.structure), 64)
        self.assertTrue(partial.stages[1]["error"])
        marked = Pipeline().select_property(0.5, 2).delete_selected().run(s, property=[i % 2 for i in range(32)])
        self.assertEqual(len(marked.structure), 16)
        computed = Pipeline().compute_property("coordination", cutoff=3).run(s)
        self.assertEqual(set(computed.property), {12.0})

    @unittest.skipUnless(os.environ.get("ATOMFORGE_PATH"), "needs the AtomForge executable")
    def test_build_and_edit_steps(self):
        ids = {t["id"] for t in modifiers()}
        for step in ("build-bulk", "build-vacancy", "build-surface", "build-sss", "replicate"):
            self.assertIn(step, ids)
        p = Pipeline().build_bulk().replicate(3, 3, 3).build_vacancy("--count 5 --seed 2").select_element("Cu")
        result = p.run(af.Structure())
        self.assertEqual(len(result.structure), 103)
        self.assertEqual(result.selected_count, 103)
        self.assertEqual([st["type"] for st in result.stages][:3], ["build-bulk", "replicate", "build-vacancy"])
        self.assertIn("build-vacancy --count 5 --seed 2", p.to_text())
        # Mixed with core steps and text, both ways.
        again = Pipeline.parse(p.to_text()).run(af.Structure())
        self.assertEqual(len(again.structure), 103)
        with self.assertRaises(PipelineError) as failure:
            Pipeline().build_vacancy("--count 999").run(copper())
        self.assertIn("vacanc", str(failure.exception))


if __name__ == "__main__":
    unittest.main()
