"""01_fixed_points_and_stability.py  --  Sections 4 and 5 of the paper.

For Examples 8.1-8.3:
  * fixed points E1-E4 and the existence quantity Delta (Theorem 4.1),
    checked against the residual of the fixed point equations;
  * the stability quantities T, D, rho4, h*, 2/rho4, h_E4 (Section 5);
  * a direct test of Theorem 5.4: eigenvalue moduli of J(E4) just below and
    just above h_E4.
Every value is compared with the number printed in the paper.

Run:  python 01_fixed_points_and_stability.py   (writes output/01_fixed_points_and_stability.txt)
"""
import numpy as np
import model as M
from report import Report

R = Report("01_fixed_points_and_stability")

# values printed in the paper (Section 8), 4 decimals
PAPER = {
    "Example 8.1": dict(E4=(2.1403, 3.4505, 0.6007), Delta=2.2074, T=-1.7889, D=0.6108, hstar=1.5044, rho4=1.8765, hF=1.0658),
    "Example 8.2": dict(E4=(2.2056, 2.8362, 3.4005), Delta=1.8504, T=-4.2968, D=2.6694, hstar=0.5644, rho4=1.7283, hF=1.1572),
    "Example 8.3": dict(E4=(0.7976, 2.6288, 0.6268), Delta=1.5177, T=-4.6664, D=5.9779, hstar=0.7806, rho4=0.7785, hF=2.5690),
}


def field(X, p):
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    x, y, z = X
    return np.array([x * (a1 - b1 * x) + be * z,
                     y * (a2 - b2 * y + c * x / (m * x + n * y + 1) - al * z),
                     z * (a3 / x - b3 * z)])


for name, p in M.EXAMPLES.items():
    f = M.fixed_points(p)
    E1, E2 = M.boundary_points(p)
    E3 = (f["xs"], 0.0, f["zs"])
    E4 = (f["xs"], f["ys"], f["zs"])
    R.line(f"\n{name}")
    for lab, E in [("E1", E1), ("E2", E2), ("E3", E3), ("E4", E4)]:
        res = np.abs(field(E, p)).max()
        R.line(f"  {lab} = ({E[0]:.4f}, {E[1]:.4f}, {E[2]:.4f})   residual |F(E)| = {res:.1e}")
        R.check(f"{lab} is a fixed point", res < 1e-10)
    R.line(f"  Delta = {f['Delta']:.4f}  (> 0, so E4 exists)")
    R.line(f"  T = {f['T']:.4f}, D = {f['D']:.4f}, T^2-4D = {f['T']**2 - 4*f['D']:.4f}, "
           f"h* = {f['hstar']:.5f}, rho4 = {f['rho4']:.4f}, 2/rho4 = {f['hF']:.5f}, h_E4 = {f['hE4']:.5f}")
    R.line(f"  E4 loses stability through: {f['route']}")
    pv = PAPER[name]
    R.check("E4 matches the paper", np.allclose(E4, pv["E4"], atol=6e-5))
    for key in ("Delta", "T", "D", "hstar", "rho4", "hF"):
        R.check(f"{key} matches the paper", abs(f[key] - pv[key]) < 6e-5)
    for fac in (0.98, 1.02):
        mods = np.sort(np.abs(np.linalg.eigvals(M.jac(*E4, fac * f["hE4"], p))))
        R.line(f"  |eig J(E4)| at {fac} h_E4: {np.round(mods, 5)}")
    R.check("E4 is a sink just below h_E4",
            np.abs(np.linalg.eigvals(M.jac(*E4, 0.98 * f["hE4"], p))).max() < 1)
    R.check("E4 is unstable just above h_E4",
            np.abs(np.linalg.eigvals(M.jac(*E4, 1.02 * f["hE4"], p))).max() > 1)

R.finish()
