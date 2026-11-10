"""04_make_figures.py  --  regenerates all eight figures of Section 8.

Output: ../figures/ex1_bif.pdf, ex1_control.pdf, ex2_bif.pdf, ex2_phase.pdf,
        ex2_control.pdf, ex3_bif.pdf, ex3_phase.pdf, ex3_control.pdf (and .png previews)
Run:  python 04_make_figures.py        (about 30 seconds)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import model as M

# ---- style: thin marks, recessive axes, validated categorical hues -------------
CX, CY, CZ = "#2a78d6", "#eb6834", "#1baf7a"   # x blue, y orange, z aqua
INK, MUTED, GRID = "#1f1f1e", "#76756d", "#d9d8d2"
plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 9.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK, "axes.labelcolor": INK,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.4, "savefig.dpi": 450, "mathtext.fontset": "cm",
})
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures") + os.sep
os.makedirs(OUT, exist_ok=True)

P1, P2, P3 = M.EX1, M.EX2, M.EX3   # Examples 8.1, 8.2, 8.3


def start(p):
    f = M.fixed_points(p)
    return f, np.array([f["xs"], f["ys"], f["zs"]]) + 0.01


def bif_panel(ax, hs, V, color, ylabel, hcrit, label_crit, xlabel=r"$h$"):
    H = np.broadcast_to(hs, V.shape)
    ax.scatter(H.ravel(), V.ravel(), s=0.25, color=color, alpha=0.8, lw=0, rasterized=True)
    ax.axvline(hcrit, color=MUTED, lw=0.8, ls="--")
    ax.annotate(label_crit, xy=(hcrit, 1), xycoords=("data", "axes fraction"), xytext=(3, -11),
                textcoords="offset points", fontsize=8, color=INK)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)


def mle_panel(ax, hs, lam, hcrit, label_crit, xlabel=r"$h$"):
    ax.axhline(0, color=MUTED, lw=0.7)
    ax.plot(hs, lam, color=INK, lw=1.0)
    ax.fill_between(hs, 0, lam, where=lam > 0, color="#2a78d6", alpha=0.25, lw=0)
    ax.axvline(hcrit, color=MUTED, lw=0.8, ls="--")
    ax.annotate(label_crit, xy=(hcrit, 1), xycoords=("data", "axes fraction"), xytext=(3, -11),
                textcoords="offset points", fontsize=8, color=INK)
    ax.set_xlabel(xlabel); ax.set_ylabel("MLE")


def save(fig, name):
    fig.savefig(OUT + name + ".pdf", bbox_inches="tight")
    fig.savefig(OUT + name + ".png", bbox_inches="tight", dpi=170)
    plt.close(fig)


summary = {}

# ---------------- Figure 1: Example 8.1, flip in the y-direction ------------------
f, X0 = start(P1)
hs = np.linspace(0.95, 1.52, 1400)
X, Y, Z, alive = M.bifurcation(P1, hs, X0, n_trans=3000, n_keep=250)
hl = np.linspace(0.95, 1.52, 500); lam = M.mle(P1, hl, X0, 3000, 5000)
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7))
bif_panel(axs[0], hs, Y, CY, r"$y$", f["hF"], r"$h_F$")
mle_panel(axs[1], hl, lam, f["hF"], r"$h_F$")
axs[0].set_title("(a) bifurcation diagram of $y$"); axs[1].set_title("(b) maximal Lyapunov exponent")
fig.tight_layout(); save(fig, "ex1_bif")
summary["fig1"] = dict(hF=f["hF"], xrange=(np.nanmin(X), np.nanmax(X)), zrange=(np.nanmin(Z), np.nanmax(Z)),
                       first_pos_mle=hl[np.where(lam > 0.01)[0][0]] if (lam > 0.01).any() else None,
                       max_mle=(hl[np.nanargmax(lam)], np.nanmax(lam)), alive_to=hs[alive][-1])

# ---------------- Figure 3: Example 8.2, flip in the (x,z)-block ------------------
f, X0 = start(P2)
hs = np.linspace(0.50, 0.758, 1400)
X, Y, Z, alive = M.bifurcation(P2, hs, X0, n_trans=3000, n_keep=250)
hl = np.linspace(0.50, 0.758, 500); lam = M.mle(P2, hl, X0, 3000, 5000)
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.2))
bif_panel(axs[0, 0], hs, X, CX, r"$x$", f["hstar"], r"$h_1$")
bif_panel(axs[0, 1], hs, Y, CY, r"$y$", f["hstar"], r"$h_1$")
bif_panel(axs[1, 0], hs, Z, CZ, r"$z$", f["hstar"], r"$h_1$")
mle_panel(axs[1, 1], hl, lam, f["hstar"], r"$h_1$")
for a, t in zip(axs.ravel(), ["(a) $x$", "(b) $y$", "(c) $z$", "(d) maximal Lyapunov exponent"]):
    a.set_title(t)
fig.tight_layout(); save(fig, "ex2_bif")
summary["fig2"] = dict(h1=f["hstar"], first_pos_mle=hl[np.where(lam > 0.01)[0][0]],
                       max_mle=(hl[np.nanargmax(lam)], np.nanmax(lam)), alive_to=hs[alive][-1])

# ---------------- Figure 6: Example 8.3, Neimark-Sacker ---------------------------
f, X0 = start(P3)
hs = np.linspace(0.74, 0.916, 1400)
X, Y, Z, alive = M.bifurcation(P3, hs, X0, n_trans=4000, n_keep=300)
hl = np.linspace(0.74, 0.916, 500); lam = M.mle(P3, hl, X0, 4000, 6000)
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.2))
bif_panel(axs[0, 0], hs, X, CX, r"$x$", f["hstar"], r"$h_{NS}$")
bif_panel(axs[0, 1], hs, Y, CY, r"$y$", f["hstar"], r"$h_{NS}$")
bif_panel(axs[1, 0], hs, Z, CZ, r"$z$", f["hstar"], r"$h_{NS}$")
mle_panel(axs[1, 1], hl, lam, f["hstar"], r"$h_{NS}$")
for a, t in zip(axs.ravel(), ["(a) $x$", "(b) $y$", "(c) $z$", "(d) maximal Lyapunov exponent"]):
    a.set_title(t)
fig.tight_layout(); save(fig, "ex3_bif")
summary["fig3"] = dict(hNS=f["hstar"], first_pos_mle=hl[np.where(lam > 0.01)[0][0]],
                       max_mle=(hl[np.nanargmax(lam)], np.nanmax(lam)), alive_to=hs[alive][-1])

# ---------------- Figure 7: Example 8.3, phase portraits --------------------------
E = np.array([f["xs"], f["ys"], f["zs"]])
hvals = [0.77, 0.79, 0.82, 0.86, 0.89, 0.91]
fig = plt.figure(figsize=(7.4, 5.0))
mles4 = {}
for i, h in enumerate(hvals):
    o = M.orbit(P3, h, X0, 30000)
    ax = fig.add_subplot(2, 3, i + 1, projection="3d")
    if h < f["hstar"]:
        tr = o[:400]
        ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], color=CX, lw=0.35, alpha=0.8)
    elif h in (0.82, 0.86):
        tr = o[-40:]
        ax.scatter(tr[:, 0], tr[:, 1], tr[:, 2], s=9, color=CX, lw=0, depthshade=False)
    else:
        tr = o[10000:]
        ax.scatter(tr[:, 0], tr[:, 1], tr[:, 2], s=0.15, color=CX, alpha=0.6, lw=0, rasterized=True)
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_locator(MaxNLocator(3))
        fmt = ScalarFormatter(useOffset=False); axis.set_major_formatter(fmt)
    ax.scatter(*E, s=14, color=CY, depthshade=False, zorder=5)
    ax.set_xlabel("$x$", labelpad=-6); ax.set_ylabel("$y$", labelpad=-6); ax.set_zlabel("$z$", labelpad=-6)
    ax.tick_params(labelsize=6, pad=-2)
    ax.set_title(f"({'abcdef'[i]}) $h={h}$", pad=-2)
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    mles4[h] = float(M.mle(P3, [h], X0, 3000, 6000)[0])
fig.tight_layout(); save(fig, "ex3_phase")
summary["fig4"] = mles4

# ---------------- Figure 4: Example 8.2, phase portraits --------------------------
f2, X02 = start(P2)
E2 = np.array([f2["xs"], f2["ys"], f2["zs"]])
hvals = [0.55, 0.62, 0.69, 0.705, 0.73, 0.75]
fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6))
mles5 = {}; periods5 = {}
for i, h in enumerate(hvals):
    o = M.orbit(P2, h, X02, 30000)
    a = axs.ravel()[i]
    if h < f2["hstar"]:
        tr = o[:150]
        a.plot(tr[:, 0], tr[:, 2], "-o", color=CX, lw=0.4, ms=2.2, mew=0)
    else:
        tr = o[20000:]
        a.scatter(tr[:, 0], tr[:, 2], s=0.4 if h > 0.72 else 14, color=CX, lw=0, rasterized=True)
    a.scatter(E2[0], E2[2], s=16, color=CY, marker="x", lw=1.0, zorder=5)
    a.set_xlabel("$x$"); a.set_ylabel("$z$"); a.set_title(f"({'abcdef'[i]}) $h={h}$")
    a.ticklabel_format(useOffset=False)
    mles5[h] = float(M.mle(P2, [h], X02, 3000, 6000)[0])
    uniq = np.unique(np.round(tr[-64:], 6), axis=0); periods5[h] = len(uniq)
fig.tight_layout(); save(fig, "ex2_phase")
summary["fig5"] = dict(mle=mles5, period=periods5)

# ---------------- Figures 2, 5, 8: chaos control, one per example ---------------
def control_fig(p, h, Bmin, Bdemo, var, name):
    f, X0 = start(p)
    E = np.array([f["xs"], f["ys"], f["zs"]])
    Bstar = f["hE4"] / h
    Bs = np.linspace(Bmin, 1.0, 1200)
    Xc, Yc, Zc, _ = M.bifurcation(p, Bs * h, X0, n_trans=4000, n_keep=250)
    Bl = np.linspace(Bmin, 1.0, 400)
    lamB = M.mle(p, Bl * h, X0, 3000, 5000)
    o_un = M.orbit(p, h, X0, 400)
    o_ct = M.orbit(p, h, X0, 400, B=Bdemo)
    V, col, lab = [(Xc, CX, "x"), (Yc, CY, "y"), (Zc, CZ, "z")][var]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.2))
    bif_panel(axs[0, 0], Bs, V, col, f"${lab}$", Bstar, r"$B^*$", xlabel=r"$B$")
    mle_panel(axs[0, 1], Bl, lamB, Bstar, r"$B^*$", xlabel=r"$B$")
    k = np.arange(o_un.shape[0])
    for a, o, t in [(axs[1, 0], o_un, f"(c) uncontrolled, $h={h}$"), (axs[1, 1], o_ct, f"(d) controlled, $h={h}$, $B={Bdemo}$")]:
        a.plot(k, o[:, 0], color=CX, lw=0.8, label="$x$")
        a.plot(k, o[:, 1], color=CY, lw=0.8, label="$y$")
        a.plot(k, o[:, 2], color=CZ, lw=0.8, label="$z$")
        a.set_xlabel("iteration $k$"); a.set_ylabel("population"); a.set_title(t)
        a.legend(loc="upper right", fontsize=7.5, frameon=False, ncol=3)
    axs[0, 0].set_title(f"(a) controlled system: ${lab}$ versus $B$")
    axs[0, 1].set_title("(b) maximal Lyapunov exponent")
    ymax = max(o_un.max(), o_ct.max()) * 1.18
    axs[1, 0].set_ylim(0, ymax); axs[1, 1].set_ylim(0, ymax)
    fig.tight_layout(); save(fig, name)
    dist = np.linalg.norm(o_ct - E, axis=1)
    return dict(Bstar=Bstar, conv_iter=int(np.argmax(dist < 1e-4)) if (dist < 1e-4).any() else None,
                mle_uncontrolled=float(M.mle(p, [h], X0, 3000, 6000)[0]))

summary["ctrl1"] = control_fig(P1, 1.4, 0.60, 0.70, 1, "ex1_control")
summary["ctrl2"] = control_fig(P2, 0.74, 0.60, 0.70, 0, "ex2_control")
summary["ctrl3"] = control_fig(P3, 0.91, 0.70, 0.80, 0, "ex3_control")

print("All figures written to", os.path.abspath(OUT))
