"""model.py -- shared definitions for all Python scripts.

Euler-discretized three-species commensalism-amensalism model.

x_{k+1} = x + h[x(a1 - b1 x) + beta z]
y_{k+1} = y + h y[a2 - b2 y + c x/(m x + n y + 1) - alpha z]
z_{k+1} = z + h z[a3/x - b3 z]

Parameter vector order: (a1, b1, a2, b2, a3, b3, beta, alpha, c, m, n)
"""
import numpy as np

# the three parameter sets of Section 8 (order: a1, b1, a2, b2, a3, b3, beta, alpha, c, m, n)
EX1 = (1.2, 0.6, 0.5, 0.4, 0.9, 0.7, 0.3, 0.6, 2.0, 0.5, 0.4)   # Example 8.1: flip in y-direction
EX2 = (3.0, 1.5, 1.0, 0.5, 1.5, 0.2, 0.2, 0.2, 1.0, 0.2, 0.2)   # Example 8.2: flip in (x,z)-block
EX3 = (2.0, 3.0, 0.5, 0.2, 1.5, 3.0, 0.5, 0.5, 3.0, 1.0, 2.0)   # Example 8.3: Neimark-Sacker
EXAMPLES = {"Example 8.1": EX1, "Example 8.2": EX2, "Example 8.3": EX3}

NAMES = ("a1", "b1", "a2", "b2", "a3", "b3", "beta", "alpha", "c", "m", "n")


def step(x, y, z, h, p):
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    xn = x + h * (x * (a1 - b1 * x) + be * z)
    yn = y + h * y * (a2 - b2 * y + c * x / (m * x + n * y + 1) - al * z)
    zn = z + h * z * (a3 / x - b3 * z)
    return xn, yn, zn


def jac(x, y, z, h, p):
    """Jacobian of the map, vectorized: returns array (..., 3, 3)."""
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    D = m * x + n * y + 1
    J = np.zeros(np.broadcast(x, y, z, h).shape + (3, 3))
    J[..., 0, 0] = 1 + h * (a1 - 2 * b1 * x)
    J[..., 0, 2] = h * be
    J[..., 1, 0] = h * c * y * (n * y + 1) / D**2
    J[..., 1, 1] = 1 + h * (a2 - 2 * b2 * y - al * z + c * x * (m * x + 1) / D**2)
    J[..., 1, 2] = -h * al * y
    J[..., 2, 0] = -h * a3 * z / x**2
    J[..., 2, 2] = 1 + h * (a3 / x - 2 * b3 * z)
    return J


def fixed_points(p):
    """Theorem 4.1 and Section 5: returns x*, y*, z*, Delta, T, D, rho4, h*, hF, hE4 and the route."""
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    roots = np.roots([1, -a1 / b1, 0, -a3 * be / (b1 * b3)])
    xs = max(r.real for r in roots if abs(r.imag) < 1e-9 and r.real > 0)
    zs = a3 / (b3 * xs)
    L = a2 - al * zs
    Delta = L + c * xs / (m * xs + 1)
    T = -(b1 * xs + be * zs / xs + a3 / xs)
    Dd = a3 * b1 + 2 * a3 * be * zs / xs**2
    out = dict(xs=xs, zs=zs, Delta=Delta, T=T, D=Dd, ys=None, rho4=None)
    if T**2 >= 4 * Dd:
        out["hstar"] = (-T - np.sqrt(T**2 - 4 * Dd)) / Dd
        out["route"] = "flip in (x,z)-block"
    else:
        out["hstar"] = -T / Dd
        out["route"] = "Neimark-Sacker"
    if Delta > 0:
        B = b2 * (m * xs + 1) - n * L
        C = L * (m * xs + 1) + c * xs
        ys = (-B + np.sqrt(B**2 + 4 * b2 * n * C)) / (2 * b2 * n)
        S = m * xs + n * ys + 1
        rho4 = ys * (b2 + c * n * xs / S**2)
        out.update(ys=ys, rho4=rho4, hF=2 / rho4, hE4=min(out["hstar"], 2 / rho4))
        if 2 / rho4 < out["hstar"]:
            out["route"] = "flip in y-direction"
    return out


def bifurcation(p, hs, X0, n_trans=3000, n_keep=300):
    """Iterate all step sizes in hs at once. Returns arrays (n_keep, len(hs)) for x, y, z.
    Orbits that leave the positive octant are marked with NaN."""
    hs = np.asarray(hs, float)
    x = np.full_like(hs, X0[0]); y = np.full_like(hs, X0[1]); z = np.full_like(hs, X0[2])
    alive = np.ones_like(hs, bool)
    with np.errstate(all="ignore"):
        for _ in range(n_trans):
            x, y, z = step(x, y, z, hs, p)
            alive &= (x > 0) & (y > 0) & (z > 0) & np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
            x = np.where(alive, x, 1.0); y = np.where(alive, y, 1.0); z = np.where(alive, z, 1.0)
        X = np.empty((n_keep, hs.size)); Y = np.empty_like(X); Z = np.empty_like(X)
        for k in range(n_keep):
            x, y, z = step(x, y, z, hs, p)
            alive &= (x > 0) & (y > 0) & (z > 0) & np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
            x = np.where(alive, x, 1.0); y = np.where(alive, y, 1.0); z = np.where(alive, z, 1.0)
            X[k], Y[k], Z[k] = x, y, z
    X[:, ~alive] = np.nan; Y[:, ~alive] = np.nan; Z[:, ~alive] = np.nan
    return X, Y, Z, alive


def mle(p, hs, X0, n_trans=3000, n_avg=6000):
    """Maximal Lyapunov exponent for every h in hs (tangent-vector method)."""
    hs = np.asarray(hs, float)
    x = np.full_like(hs, X0[0]); y = np.full_like(hs, X0[1]); z = np.full_like(hs, X0[2])
    alive = np.ones_like(hs, bool)
    v = np.ones((hs.size, 3)) / np.sqrt(3)
    s = np.zeros_like(hs)
    with np.errstate(all="ignore"):
        for _ in range(n_trans):
            x, y, z = step(x, y, z, hs, p)
            alive &= (x > 0) & (y > 0) & (z > 0) & np.isfinite(x)
            x = np.where(alive, x, 1.0); y = np.where(alive, y, 1.0); z = np.where(alive, z, 1.0)
        for _ in range(n_avg):
            J = jac(x, y, z, hs, p)
            v = np.einsum("hij,hj->hi", J, v)
            nv = np.linalg.norm(v, axis=1)
            s += np.log(nv)
            v /= nv[:, None]
            x, y, z = step(x, y, z, hs, p)
            alive &= (x > 0) & (y > 0) & (z > 0) & np.isfinite(x)
            x = np.where(alive, x, 1.0); y = np.where(alive, y, 1.0); z = np.where(alive, z, 1.0)
    lam = s / n_avg
    lam[~alive] = np.nan
    return lam


def orbit(p, h, X0, n, B=1.0):
    """Single orbit of the (optionally controlled) map; returns array (n+1, 3)."""
    out = np.empty((n + 1, 3)); out[0] = X0
    x, y, z = X0
    for k in range(n):
        xn, yn, zn = step(x, y, z, h, p)
        x, y, z = B * xn + (1 - B) * x, B * yn + (1 - B) * y, B * zn + (1 - B) * z
        out[k + 1] = (x, y, z)
    return out


def boundary_points(p):
    """E1 and E2 of Theorem 4.1."""
    a1, b1, a2, b2, a3, b3, be, al, c, m, n = p
    x1 = a1 / b1
    A2 = a1 * m / (b1 * n) + 1 / n - a2 / b2
    A3 = -(a1 * a2 * m / (b1 * b2 * n) + a2 / (b2 * n) + a1 * c / (b1 * b2 * n))
    y2 = (-A2 + np.sqrt(A2**2 - 4 * A3)) / 2
    return (x1, 0.0, 0.0), (x1, y2, 0.0)


def period_of(p, h, X0, n=60000, tol=1e-7, maxp=128):
    """Period of the attractor reached from X0 (0 = not periodic up to maxp, -1 = left the positive octant)."""
    o = orbit(p, h, X0, n)
    t = o[-300:]
    if not np.isfinite(t).all() or (t <= 0).any():
        return -1
    for k in range(1, maxp + 1):
        if np.allclose(t[-1], t[-1 - k], atol=tol):
            return k
    return 0
