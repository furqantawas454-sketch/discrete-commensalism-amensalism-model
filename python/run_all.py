"""run_all.py -- runs the three verification scripts and regenerates the figures, in order."""
import subprocess, sys, os

here = os.path.dirname(os.path.abspath(__file__))
scripts = ["01_fixed_points_and_stability.py", "02_bifurcation_checks.py",
           "03_chaos_control_checks.py", "04_make_figures.py"]
for s in scripts:
    print(f"\n{'=' * 70}\nRunning {s}\n{'=' * 70}")
    r = subprocess.run([sys.executable, os.path.join(here, s)], cwd=here)
    if r.returncode != 0:
        sys.exit(f"{s} reported a failure")
print("\nAll scripts finished successfully.")
