"""Exercise build serialization, argument forwarding and failure recovery."""

from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


CMAKE = sys.argv.pop(1)
LAUNCHER = Path(__file__).resolve().parents[1] / "cmake" / "RunBuildCommand.cmake"


class BuildResourceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="atomforge build test ")
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.prefix = [
            CMAKE, f"-DATOMFORGE_BUILD_LOCK={self.directory / 'compiler.lock'}",
            "-P", str(LAUNCHER), "--", sys.executable,
        ]

    def run_child(self, *args):
        return subprocess.run(
            self.prefix + list(args), capture_output=True, text=True, timeout=30
        )

    def test_parallel_commands_are_serialized(self):
        child = self.directory / "timed child.py"
        child.write_text(
            "import json, pathlib, sys, time\n"
            "start = time.monotonic_ns()\n"
            "time.sleep(0.1)\n"
            "pathlib.Path(sys.argv[1]).write_text("
            "json.dumps([start, time.monotonic_ns()]))\n"
        )
        records = [self.directory / f"job {i}.json" for i in range(8)]
        with ThreadPoolExecutor(max_workers=8) as executor:
            results = list(executor.map(
                lambda record: self.run_child(str(child), str(record)), records
            ))
        for result in results:
            self.assertEqual(result.returncode, 0, result.stderr)
        intervals = sorted(json.loads(record.read_text()) for record in records)
        for previous, following in zip(intervals, intervals[1:]):
            self.assertLessEqual(previous[1], following[0])

    def test_arguments_and_output_are_preserved(self):
        arguments = ["path with spaces", "semi;colon", 'quote"value', "-DVALUE=1"]
        child = self.directory / "arguments.py"
        child.write_text("import json, sys\nprint(json.dumps(sys.argv[1:]))\n")
        result = self.run_child(str(child), *arguments)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), arguments)

    def test_failure_releases_lock_and_reaches_build_tool(self):
        result = self.run_child("-c", "raise SystemExit(7)")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("build command failed: 7", result.stderr)
        result = self.run_child("-c", "print('recovered')")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("recovered", result.stdout)


if __name__ == "__main__":
    unittest.main()
