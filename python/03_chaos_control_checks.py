"""03_chaos_control_checks.py  --  Section 7 and the 'Chaos control' paragraphs of Section 8.

For each example at its chaotic step size h:
  * the maximal Lyapunov exponent of the uncontrolled system (positive = chaos);
  * B* = h_E4 / h (Theorem 7.1);
  * eigenvalue moduli of J_B(E4) = (1-B) I + B J(E4) at B* -/+ 0.02;
  * direct iteration of the controlled map (7.1) at B* -/+ 0.02;
  * the number of iterations needed to reach E4 (within 1e-4) for the B used in the figures;
  * Lemma 7.1: the controlled map equals the Euler map with step B h.

Run:  python 03_chaos_control_checks.py   (writes output/03_chaos_control_checks.txt)
"""
import numpy as np
import model as M
from report import Report

R = Report("03_chaos_control_checks")

CASES = [("Example 8.1", M.EX1, 1.40, 0.70, 0.7613),
         ("Example 8.2", M.EX2, 0.74, 0.70, 0.7627),
         ("Example 8.3", M.EX3, 0.91, 0.80, 0.8578)]

rng = np.random.default_rng(0)
for name, p, h, Bfig, Bpaper in CASES:
    f = M.fixed_points(p)
    E = np.array([f["xs"], f["ys"], f["zs"]])
    R.line(f"\n{name}  (h = {h})")
    # Lemma 7.1 at random points
    X = E * (1 + 0.3 * rng.random(3)); B = 0.77
    lhs = B * np.array(M.step(*X, h, p)) + (1 - B) * X
    rhs = np.array(M.step(*X, B * h, p))
    R.check("Lemma 7.1: B G(X) + (1-B) X == G_{Bh}(X)", np.allclose(lhs, rhs, atol=1e-13))
    lam = float(M.mle(p, [h], E + 0.01, 3000, 6000)[0])
    R.line(f"  MLE of the uncontrolled system = {lam:.3f}")
    R.check("uncontrolled system is chaotic (MLE > 0)", lam > 0)
    Bs = f["hE4"] / h
    R.line(f"  B* = h_E4 / h = {f['hE4']:.5f} / {h} = {Bs:.4f}")
    R.check(f"B* matches the paper ({Bpaper})", abs(Bs - Bpaper) < 6e-5)
    for dB in (-0.02, 0.02):
        Bc = Bs + dB
        JB = (1 - Bc) * np.eye(3) + Bc * M.jac(*E, h, p)
        mods = np.sort(np.abs(np.linalg.eigvals(JB)))
        t = M.orbit(p, h, E + 0.01, 20000, B=Bc)[-100:]
        d = np.linalg.norm(t - E, axis=1).max()
        R.line(f"  B = B* {dB:+.2f}: |eig J_B(E4)| = {np.round(mods, 4)}, distance from E4 after 20000 steps = {d:.2e}")
        if dB < 0:
            R.check("E4 stabilized for B < B*", mods.max() < 1 and d < 1e-8)
        else:
            R.check("E4 not stabilized for B > B*", mods.max() > 1 and d > 1e-3)
    o = M.orbit(p, h, E + 0.01, 400, B=Bfig)
    k = int(np.argmax(np.linalg.norm(o - E, axis=1) < 1e-4))
    R.line(f"  with B = {Bfig} (used in the figure) the orbit reaches E4 within 1e-4 after {k} iterations")

R.finish()
