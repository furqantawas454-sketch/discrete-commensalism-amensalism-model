"""report.py -- tiny helper that prints results and saves them to output/<name>.txt."""
import os


class Report:
    def __init__(self, name):
        self.name = name
        self.lines = []
        self.n_pass = 0
        self.n_fail = 0

    def line(self, text=""):
        print(text)
        self.lines.append(text)

    def check(self, label, ok):
        tag = "PASS" if ok else "FAIL"
        self.n_pass += bool(ok)
        self.n_fail += (not ok)
        self.line(f"    [{tag}] {label}")

    def finish(self):
        self.line(f"\nSummary: {self.n_pass} checks passed, {self.n_fail} failed.")
        here = os.path.dirname(os.path.abspath(__file__))
        os.makedirs(os.path.join(here, "output"), exist_ok=True)
        with open(os.path.join(here, "output", self.name + ".txt"), "w") as fh:
            fh.write("\n".join(self.lines) + "\n")
        if self.n_fail:
            raise SystemExit(1)
