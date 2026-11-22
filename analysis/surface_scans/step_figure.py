"""What Gaussian's VPT2 step sees: one X-H stretch per energy model, against B3LYP (appendix figure).

Per panel: the model's curvature along its own normal mode as a local wavenumber with the
straight-line part removed (as in surface_zoom), the B3LYP curve of the same mode, the three points
Gaussian's default step uses (Q = 0, +/-0.010 amu^1/2 A) on the model curve, and the parabola in
k(Q) through these three points (the cubic and quartic constant Gaussian forms from them), converted
to the same local wavenumber and with the same straight-line part removed. The label gives Gaussian's
printed k(i,i,i,i) of the model run and of the B3LYP run.

Grids: the 0.0002 A zoom (zoom.py) within +/-0.012 A where the panel needs it; MACE-OFF23 methanol
O-H on the 0.002 A fine grid over +/-0.05 A (smooth, no wobble at 0.0002 A, so no aliasing; its
stronger bend only shows over the wider range). B3LYP always on the fine grid (smooth).

python step_figure.py <thesis figures dir> [2x2|row]   -> surface_step.{pdf,png} (2 x 2; the row of four,
where the labels crowd the y ticks, as surface_step_row); data/step_figure.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402
from analyse import LABEL, load_scan, style  # noqa: E402

# (target, model, grid of the model curve, half range in dr / A, mode label)
PANELS = [
    ("glycine_nh", "mace_anicc", "zoom", 0.012, "glycine N–H"),
    ("glycine_nh", "mace_omol", "zoom", 0.012, "glycine N–H"),
    ("methanol_oh", "mace_off", "fine", 0.050, "methanol O–H"),
    ("methanol_ch", "mace_polar", "zoom", 0.012, "methanol C–H"),
]
Q_DEFAULT = 0.010  # amu^1/2 A


def curve(key, m, grid, half):
    """q, dr, k on the requested grid within |dr| <= half."""
    if grid == "zoom":
        d = c.load(f"zoom_{key}_{m}.json")
    else:
        d = load_scan(key, m)["grids"]["fine"]
    q, dr, k = (np.array(d[x]) for x in ("q", "dr", "k"))
    sel = np.abs(dr) <= half + 1e-9
    return q[sel], dr[sel], k[sel]


def gaussian_points(key, m):
    """Q, dr, k at Q = -0.010, 0, +0.010 (steps_q.py)."""
    s = c.load(f"stepsq_{m}_{key}.json")
    j = s["steps_q"].index(Q_DEFAULT)
    q = np.array([s["q"][2 * j], 0.0, s["q"][2 * j + 1]])
    dr = np.array([s["dr"][2 * j], 0.0, s["dr"][2 * j + 1]])
    k = np.array([s["k"][2 * j], s["k0"], s["k"][2 * j + 1]])
    return q, dr, k


def nu(k, k0):
    return c.k_to_wavenumber(k) - c.k_to_wavenumber(k0)


def main(out: Path, layout: str = "2x2") -> None:
    res = pd.read_csv(c.DATA / "results.csv")
    style.apply_style()
    W = style.figwidth()
    # B3LYP's own three-point parabola is not drawn: its curvature does not depend on the step (quartic
    # 421.3 / 421.0 / 420.2 at Q = 0.005 / 0.010 / 0.020), so the parabola lies on its curve
    if layout == "row":
        fig, axes = plt.subplots(1, 4, figsize=(W, 2.35))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(W, 4.1))
    rows = []
    for ax, (key, m, grid, half, mode_lab) in zip(axes.flat, PANELS):
        col = style.model_color(m)
        q, dr, k = curve(key, m, grid, half)
        gq, gdr, gk = gaussian_points(key, m)
        k0 = gk[1]
        y = nu(k, k0)
        lin = np.polynomial.polynomial.polyfit(dr, y, 1)

        def flat(x, v):
            return v - np.polynomial.polynomial.polyval(x, lin)  # noqa: B023

        # the parabola k(Q) = k0 + a Q + b Q^2 through Gaussian's three points, on the model's Q grid
        a = (gk[2] - gk[0]) / (gq[2] - gq[0])
        b = (gk[2] + gk[0] - 2 * k0) / (2 * gq[2] ** 2)  # gq symmetric: gq[0] = -gq[2]
        kpar = k0 + a * q + b * q**2
        # B3LYP, the same mode on its own minimum, own straight-line part over the same range
        _, bdr, bk = curve(key, "b3lyp", "fine", half)
        bk0 = load_scan(key, "b3lyp")["grids"]["zero"]["k"][0]
        by = nu(bk, bk0)
        blin = np.polynomial.polynomial.polyfit(bdr, by, 1)
        ax.plot(bdr, by - np.polynomial.polynomial.polyval(bdr, blin), color=style.INK, lw=1.1, zorder=2)
        ax.plot(dr, flat(dr, y), color=col, lw=1.0, zorder=3)
        ax.plot(dr, flat(dr, nu(kpar, k0)), color=col, lw=0.8, ls=(0, (3, 2)), zorder=4)
        ax.plot(gdr, flat(gdr, nu(gk, k0)), "o", color=col, ms=4.2, mec="white", mew=0.6, zorder=5)
        style.clean_axes(ax, grid=False)
        r = res[res.key == key].set_index("method")
        kg, kb = r.loc[m, "quartic_gauss"], r.loc["b3lyp", "quartic_gauss"]
        ax.text(0.5, 0.97, f"$k_{{iiii}}$ {kg:.0f} / B3LYP {kb:.0f}", transform=ax.transAxes, ha="center",
                va="top", fontsize=style.ANNOT_FS - 0.5, color=style.INK, zorder=6,
                bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.85))  # fmt: skip
        ax.set_title(f"{style.tex(LABEL[m])}\n{mode_lab}", color=col, fontsize=style.ANNOT_FS)
        ax.set_xlim(-half * 1.04, half * 1.04)
        ax.set_xticks([-0.01, 0, 0.01] if half < 0.02 else [-0.04, 0, 0.04])
        ax.tick_params(labelsize=style.ANNOT_FS - 1)
        # y range: all plotted data plus headroom for the label (nothing off the scale)
        ys = np.concatenate([flat(dr, y), flat(dr, nu(kpar, k0)), by - np.polynomial.polynomial.polyval(bdr, blin)])
        lo, hi = ys.min(), ys.max()
        ax.set_ylim(lo - 0.06 * (hi - lo), hi + 0.30 * (hi - lo))
        rows.append(dict(key=key, model=m, grid=grid, half_range_A=half, n_points=len(dr),
                         dr_minus=gdr[0], dr_plus=gdr[2], q_minus=gq[0], q_plus=gq[2],
                         quartic_gauss=kg, quartic_gauss_b3lyp=kb, quartic_3pt=2 * b, cubic_3pt=a))  # fmt: skip
    if layout == "row":
        axes[0].set_ylabel(r"$\tilde\nu_{\mathrm{loc}}$, lin. part removed / cm$^{-1}$", fontsize=style.ANNOT_FS)
    else:
        for r_ in range(2):
            axes[r_, 0].set_ylabel(r"$\tilde\nu_{\mathrm{loc}}$, lin. part" + "\n" + r"removed / cm$^{-1}$", fontsize=style.ANNOT_FS)
    fig.supxlabel(r"$\Delta r$(X–H) / Å", fontsize=style.BASE_FS)
    g = style.INK_MUTED
    handles = [Line2D([], [], color=style.INK, lw=1.1, label="B3LYP"),
               Line2D([], [], color=g, lw=1.0, label="model"),
               Line2D([], [], color=g, ls="none", marker="o", ms=4.2, mec="white", label="Gaussian's points"),
               Line2D([], [], color=g, lw=0.8, ls=(0, (3, 2)), label="parabola through them")]  # fmt: skip
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0),
               fontsize=style.ANNOT_FS, handlelength=1.8, columnspacing=1.2, borderpad=0)  # fmt: skip
    # legend directly above the panels; small outer pad: the crop removes the margin, and the pad also
    # sets the gap between the panels and the shared x label
    style.tight(fig, pad=0.4, w_pad=0.5, h_pad=0.8)
    style.save(fig, out, "surface_step" if layout == "2x2" else "surface_step_row")
    pd.DataFrame(rows).to_csv(c.DATA / "step_figure.csv", index=False)
    print(pd.DataFrame(rows).round(5).to_string(index=False))


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else c.HERE / "out", sys.argv[2] if len(sys.argv) > 2 else "2x2")
