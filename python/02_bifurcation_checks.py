"""02_bifurcation_checks.py  --  Section 6 and the 'Verification' and
'Dynamics beyond the bifurcation' paragraphs of Section 8.

  * normal form coefficients l_F (6.7), l_B (6.8) and l_NS (6.10), computed
    independently of the Mathematica script 03_bifurcation.wl;
  * the period-two orbits just after the two flips, and the invariant closed
    curve just after the Neimark-Sacker point;
  * the period-doubling cascades and the period-5 locking reported in Section 8;
  * the invariance of the line x = x*, z = z* along the route of Example 8.1.

Run:  python 02_bifurcation_checks.py   (writes output/02_bifurcation_checks.txt)
"""
import numpy as np
import model as M
from report import Report

R = Report("02_bifurcation_checks")


def planar_forms(p, h):
    """M, B and C of the planar (x,z)-map at (x*, z*), Eqs. (5.3), (6.4), (6.5)."""
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    f = M.fixed_points(p); xs, zs = f["xs"], f["zs"]
    Mm = np.eye(2) + h * np.array([[-(b1 * xs + be * zs / xs), be], [-a3 * zs / xs**2, -a3 / xs]])
    B = lambda u, v: np.array([-2*h*b1*u[0]*v[0],
                               2*h*a3*zs/xs**3*u[0]*v[0] - h*a3/xs**2*(u[0]*v[1] + u[1]*v[0]) - 2*h*b3*u[1]*v[1]])
    C = lambda u, v, w: np.array([0.0 * u[0],
                                  -6*h*a3*zs/xs**4*u[0]*v[0]*w[0] + 2*h*a3/xs**3*(u[0]*v[0]*w[1] + u[0]*v[1]*w[0] + u[1]*v[0]*w[0])])
    return Mm, B, C


def l_flip_y(p):
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    f = M.fixed_points(p); xs, ys = f["xs"], f["ys"]
    K, S = c * xs, m * xs + n * ys + 1
    hF = f["hF"]
    a = hF * (-b2 - K * n / S**2 + ys * K * n**2 / S**3)
    b = hF * (K * n**2 / S**3 - ys * K * n**3 / S**4)
    return hF, a, b, a**2 + b


def l_flip_block(p):
    f = M.fixed_points(p)
    h1 = f["hstar"]
    Mm, B, C = planar_forms(p, h1)
    w, V = np.linalg.eig(Mm); q = V[:, np.argmin(abs(w + 1))].real
    w2, U = np.linalg.eig(Mm.T); pp = U[:, np.argmin(abs(w2 + 1))].real
    pp = pp / (pp @ q)
    lB = pp @ C(q, q, q) / 6 - 0.5 * pp @ B(q, np.linalg.solve(Mm - np.eye(2), B(q, q)))
    return h1, np.sort(w.real), lB


def l_ns(Mm, B, C):
    w, V = np.linalg.eig(Mm.astype(complex))
    i = np.argmax(w.imag); lam = w[i]; q = V[:, i]
    w2, U = np.linalg.eig(Mm.T.astype(complex))
    pp = U[:, np.argmin(abs(w2 - np.conj(lam)))]
    pp = pp / np.conj(np.vdot(pp, q))          # <p,q> = conj(p).q = 1
    th = np.angle(lam)
    ip = lambda u, v: np.vdot(u, v)
    I = np.eye(2)
    val = (ip(pp, C(q, q, np.conj(q)))
           + 2 * ip(pp, B(q, np.linalg.solve(I - Mm, B(q, np.conj(q)))))
           + ip(pp, B(np.conj(q), np.linalg.solve(np.exp(2j * th) * I - Mm, B(q, q)))))
    return 0.5 * np.real(np.exp(-1j * th) * val), lam, th


def start(p):
    f = M.fixed_points(p)
    return f, np.array([f["xs"], f["ys"], f["zs"]])


# ---------------- validation of the NS formula ---------------------------------
R.line("Validation: delayed logistic map (x,y) -> (y, 2y(1-x)) at (1/2,1/2)")
A0 = np.array([[0.0, 1.0], [-1.0, 1.0]])
B0 = lambda u, v: np.array([0 * u[0], -2 * (u[0] * v[1] + u[1] * v[0])])
C0 = lambda u, v, w: np.array([0 * u[0], 0 * u[0]])
d0 = l_ns(A0, B0, C0)[0]
R.line(f"  l_NS = {d0:.4f}")
R.check("known supercritical case gives l_NS < 0", d0 < 0)

# ---------------- Example 8.1 -----------------------------------------------------
f, E = start(M.EX1)
hF, a, b, lF = l_flip_y(M.EX1)
R.line("\nExample 8.1 (Theorem 6.1)")
R.line(f"  hF = {hF:.5f}, h* = {f['hstar']:.5f}, a = {a:.4f}, b = {b:.4f}, l_F = {lF:.4f}")
R.check("hF < h*", hF < f["hstar"])
R.check("b > 0 and l_F > 0 (supercritical)", b > 0 and lF > 0)
R.check("l_F matches the paper (0.2793)", abs(lF - 0.2793) < 1e-4)
o = M.orbit(M.EX1, hF + 0.01, E + 0.01, 40000)[-2:]
half = abs(o[0, 1] - o[1, 1]) / 2
R.line(f"  period-two orbit at hF + 0.01: {np.round(o, 4).tolist()}")
R.line(f"  half-width {half:.4f}, normal-form prediction {np.sqrt(f['rho4'] * 0.01 / lF):.4f}")
R.check("only y changes on the period-two orbit", np.allclose(o[:, [0, 2]], E[[0, 2]], atol=1e-9))
R.check("half-width agrees with the prediction to 0.005", abs(half - np.sqrt(f['rho4'] * 0.01 / lF)) < 5e-3)
X, Y, Z, alive = M.bifurcation(M.EX1, np.linspace(1.07, 1.50, 200), E + 0.01, 3000, 200)
dev = max(np.nanmax(abs(X - E[0])), np.nanmax(abs(Z - E[2])))
R.line(f"  max |x - x*|, |z - z*| on the attractor for h in [1.07, 1.50]: {dev:.1e}")
R.check("cascade stays on the invariant line x = x*, z = z*", dev < 1e-9)
per = {h: M.period_of(M.EX1, h, E + 0.01) for h in (1.25, 1.31, 1.355, 1.3925)}
R.line(f"  periods: {per}")
R.check("periods 2, 4, 8 and the period-6 window", [per[1.25], per[1.31], per[1.355], per[1.3925]] == [2, 4, 8, 6])

# ---------------- Example 8.2 -----------------------------------------------------
f, E = start(M.EX2)
h1, eigM, lB = l_flip_block(M.EX2)
R.line("\nExample 8.2 (Theorem 6.2)")
R.line(f"  h1 = {h1:.5f}, hF = {f['hF']:.5f}, eig M(h1) = {np.round(eigM, 4)}, "
       f"lambda_y = {1 - h1 * f['rho4']:.4f}, l_B = {lB:.4f}")
R.check("h1 < hF", h1 < f["hF"])
R.check("l_B > 0 (supercritical) and matches the paper (0.5805)", lB > 0 and abs(lB - 0.5805) < 1e-4)
o = M.orbit(M.EX2, h1 + 0.01, E + 0.01, 40000)[-3:]
R.line(f"  period-two orbit at h1 + 0.01: {np.round(o[-2:], 4).tolist()}")
R.check("orbit has period two", np.allclose(o[0], o[2], atol=1e-8) and not np.allclose(o[0], o[1], atol=1e-3))
per = {h: M.period_of(M.EX2, h, E + 0.01) for h in (0.62, 0.69, 0.705, 0.7085)}
R.line(f"  periods: {per}")
R.check("cascade 2, 4, 8, 16", [per[0.62], per[0.69], per[0.705], per[0.7085]] == [2, 4, 8, 16])

# ---------------- Example 8.3 -----------------------------------------------------
f, E = start(M.EX3)
hNS = f["hstar"]
Mm, B, C = planar_forms(M.EX3, hNS)
lNS, lam, th = l_ns(Mm, B, C)
R.line("\nExample 8.3 (Theorem 6.3)")
R.line(f"  hNS = {hNS:.5f}, hF = {f['hF']:.5f}, lambda = {lam:.4f}, |lambda| = {abs(lam):.6f}, "
       f"theta0 = {th:.4f}, T^2/D = {f['T']**2 / f['D']:.4f}, l_NS = {lNS:.3f}")
R.check("hNS < hF", hNS < f["hF"])
R.check("|lambda| = 1", abs(abs(lam) - 1) < 1e-10)
R.check("non-resonance: T^2/D differs from 2 and 3", abs(f['T']**2 / f['D'] - 2) > 0.1 and abs(f['T']**2 / f['D'] - 3) > 0.1)
R.check("l_NS < 0 (supercritical) and matches the paper (-10.326)", lNS < 0 and abs(lNS + 10.326) < 1e-3)
for dh in (-0.01, 0.01, 0.03):
    t = M.orbit(M.EX3, hNS + dh, E + 0.01, 80000)[70000:]
    d = np.linalg.norm(t - E, axis=1)
    R.line(f"  h = hNS {dh:+.2f}: distance of the attractor from E4 in [{d.min():.4f}, {d.max():.4f}]")
    if dh < 0:
        R.check("converges to E4 before hNS", d.max() < 1e-8)
    else:
        R.check("invariant closed curve around E4 after hNS", d.min() > 0.03)
per = {h: M.period_of(M.EX3, h, E + 0.01) for h in (0.82, 0.86, 0.865)}
R.line(f"  periods after locking: {per}")
R.check("period 5, 10, 20", [per[0.82], per[0.86], per[0.865]] == [5, 10, 20])

# ---------------- Remark 9 ------------------------------------------------------
pS = (0.2, 3., 2., 3., 3., 3., 1., 0.2, 1.5, 1., 3.)
fS = M.fixed_points(pS)
lS = l_ns(*planar_forms(pS, fS["hstar"]))[0]
R.line(f"\nRemark 9: hNS = {fS['hstar']:.4f} < hF = {fS['hF']:.4f}, l_NS = {lS:.2e}")
R.check("subcritical case: l_NS > 0", lS > 0)

R.finish()
