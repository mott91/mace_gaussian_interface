"""Analysis of the surface scans: constants, wobble, energy shape; tables and figures.

Inputs: data/b3lyp_<target>.json, data/mace_<target>_<model>.json (scan_b3lyp.py, scan_mace.py),
the campaign VPT2 logs (vpt2_constants.py) and the band / chi tables of the campaign.

Constants along the mode, in Gaussian's "k" units, with +Q = X-H stretching:
  - smooth fit: polynomial of degree FIT_DEG in Q through the curvature k(Q) on the fine grid
    (|dr| <= 0.05 A); cubic = dk/dQ(0), quartic = d2k/dQ2(0);
  - finite differences at Gaussian's steps h (Cartesian A along the unit mode vector), exactly as
    VPT2 forms them: cubic = [k(+h) - k(-h)] / 2Q_h, quartic = [k(+h) + k(-h) - 2 k(0)] / Q_h^2;
  - Gaussian's printed k(i,i,i), k(i,i,i,i) of the campaign run (cubic times the phase of
    Gaussian's mode vector, so that the sign refers to stretching).
Wobble: residual of k(Q) from a smooth polynomial of degree WOBBLE_DEG over the fine grid, as a
wavenumber; its amplitude (max |residual|) and length scale (mean distance between sign changes,
in dr).

python analyse.py <thesis figures dir>   -> data/results.json, data/tables.md, figures surface_*
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402
import vpt2_constants as v  # noqa: E402

spec = importlib.util.spec_from_file_location("thesis_style", c.ROOT / "mace_gaussian" / "analysis" / "thesis_style.py")
style = importlib.util.module_from_spec(spec)
sys.modules["thesis_style"] = style
spec.loader.exec_module(style)

STEPS_Q = (0.005, 0.01, 0.02, 0.04)  # Gaussian's VPT2 step in the normal coordinate, amu^1/2 A (default 0.010)
YMAX = {"methanol_oh": 16, "methanol_ch": 16, "glycine_nh": 16}  # energy axis, 10^3 cm-1
FIT_DEG = 8  # B3LYP methanol O-H: degree 8 reproduces the h = 0.005 A differences, degree 6 is 0.2 % off
WOBBLE_DEG = 6  # smooth reference for the wobble (B3LYP residual about 0.02 cm-1)
LABEL = {"b3lyp": "B3LYP", "mace_anicc": "MACE-ANICC", "mace_omol": "MACE-OMOL-0", "mace_off": "MACE-OFF23",
         "mace_polar": "MACE-POLAR-1"}  # fmt: skip
TITLE = {"methanol_oh": "methanol, O–H stretch", "methanol_ch": "methanol, C–H stretch",
         "glycine_nh": "glycine, N–H stretch"}  # fmt: skip
HBAR = 1.054571817e-34


def color(m):
    return style.INK if m == "b3lyp" else style.model_color(m)


def load_scan(key, m):
    return c.load(f"b3lyp_{key}.json" if m == "b3lyp" else f"mace_{key}_{m}.json")


def poly(q, y, deg):
    p = np.polynomial.polynomial.polyfit(q, y, deg)
    return p, np.polynomial.polynomial.polyval(q, p)


def zpe_turning(freq, dr_dq):
    """Classical turning point of the harmonic v = 0 level, in dr (Angstrom)."""
    w = 2 * np.pi * c.C_CM * freq
    q = np.sqrt(HBAR / w) / np.sqrt(c.AMU_KG) * 1e10  # amu^1/2 Angstrom
    return q * dr_dq


def wobble(dr, resid_cm):
    s = np.sign(resid_cm)
    cross = dr[1:][s[1:] != s[:-1]]
    scale = float(np.mean(np.diff(cross))) if len(cross) > 2 else float("nan")
    return float(np.abs(resid_cm).max()), scale


def analyse():
    chi = pd.read_csv(c.CAMP / "figures" / "supervisor" / "data" / "chi.csv")
    chi = chi[chi.kind == "diagonal"]
    rows, curves = [], {}
    for key, mol, i in c.TARGETS:
        for m in c.METHODS:
            d = load_scan(key, m)
            mode = c.load_mode(key, m)
            g = d["grids"]
            q, dr, k = (np.array(g["fine"][x]) for x in ("q", "dr", "k"))
            k0 = g["zero"]["k"][0]
            p, kfit = poly(q, k, FIT_DEG)
            # wobble: on the 0.0002 A zoom grid where it exists (the 0.002 A grid aliases it)
            zp = c.DATA / f"zoom_{key}_{m}.json"
            if zp.exists():
                zz = c.load(zp.name)
                zq, zdr, zk = (np.array(zz[x]) for x in ("q", "dr", "k"))
            else:
                zq, zdr, zk = q, dr, k
            ps, ks = poly(zq, zk, WOBBLE_DEG)
            resid = c.k_to_wavenumber(zk) - c.k_to_wavenumber(ks)
            amp, scale = wobble(zdr, resid)
            # finite differences at Gaussian's points: steps in the normal coordinate Q (steps_q.py)
            fd = {}
            sq = c.load(f"stepsq_{m}_{key}.json")
            qs, kst, k0s = np.array(sq["q"]), np.array(sq["k"]), sq["k0"]
            for j, h in enumerate(sq["steps_q"]):
                qm, qp, km, kp = qs[2 * j], qs[2 * j + 1], kst[2 * j], kst[2 * j + 1]
                fd[h] = ((kp - km) / (qp - qm), (kp + km - 2 * k0s) / qp**2, float(mode.bond(qp)[0]))
            gv = v.constants(c.run_log(mol, m), mode.freq)
            cx = chi[(chi.molecule == mol) & (chi.energy_model == ("mace_anicc" if m == "b3lyp" else m)) & (chi.i == i)]
            band = c.b3lyp_band(mol, i, "mace_anicc" if m == "b3lyp" else m)
            # energy shape on the wide grid: E - E0 and its anharmonic part
            w = g["wide"]
            wq, wdr, we = np.array(w["q"]), np.array(w["dr"]), np.array(w["E"])
            e0 = g["zero"]["E"][0]
            to_cm = 1e-18 / (6.62607015e-34 * c.C_CM)  # attoJ -> cm-1
            ecm = (we - e0) * to_cm
            eharm = 0.5 * k0 * wq**2 * to_cm  # the method's own harmonic parabola
            at = {a: (float(np.interp(-a, wdr, ecm)), float(np.interp(a, wdr, ecm))) for a in (0.1, 0.2)}
            rows.append({
                "key": key, "method": m, "freq": mode.freq, "freq_scan": d["freq_scan"],
                "chi_ii": float(cx.dft_xmatrix.iloc[0] if m == "b3lyp" else cx.ml_xmatrix.iloc[0]) if len(cx) else np.nan,
                "vpt2_shift": float(band.dft_anharmonic - band.dft_harmonic if m == "b3lyp" else band.ml_anharmonic - band.ml_harmonic),
                "cubic_fit": float(p[1]), "quartic_fit": float(2 * p[2]),
                **{f"cubic_h{h}": fd[h][0] for h in c.STEPS}, **{f"quartic_h{h}": fd[h][1] for h in c.STEPS},
                **{f"dr_h{h}": fd[h][2] for h in c.STEPS},
                "cubic_gauss": gv["cubic_k"] * mode.gauss_phase, "quartic_gauss": gv["quartic_k"],
                "cubic_FI_gauss": gv["cubic_FI"] * mode.gauss_phase, "quartic_FI_gauss": gv["quartic_FI"],
                "wobble_amp_cm": amp, "wobble_scale_A": scale, "dr_dq": mode.dr_dq, "dnorm": mode.dnorm,
                **{f"purity_{a}": b for a, b in mode.purity.items()},
                "zpe_turning_A": zpe_turning(mode.freq, mode.dr_dq),
                **{f"E_m{a}": at[a][0] for a in at}, **{f"E_p{a}": at[a][1] for a in at},
                **{f"Eanh_m{a}": float(np.interp(-a, wdr, ecm - eharm)) for a in (0.1, 0.2)},
                **{f"Eanh_p{a}": float(np.interp(a, wdr, ecm - eharm)) for a in (0.1, 0.2)},
                "seconds": d.get("seconds", np.nan),
            })  # fmt: skip
            curves[(key, m)] = dict(dr=dr, w=c.k_to_wavenumber(k), w0=float(c.k_to_wavenumber(k0)), wdr=wdr, ecm=ecm,
                                    eharm=eharm, zdr=zdr, zw=c.k_to_wavenumber(zk), zoom=zp.exists())  # fmt: skip
    return pd.DataFrame(rows), curves


def figure(curves, res, out: Path, which: str = "energy"):
    """3 columns (modes) x 2 rows. Top: E - E_min (which='energy') or E minus the method's own
    harmonic parabola (which='anharmonic'); bottom: local wavenumber change from the curvature."""
    style.apply_style()
    fig, axes = plt.subplots(2, 3, figsize=(style.figwidth(), 4.4))
    for col, (key, _mol, _i) in enumerate(c.TARGETS):
        ax, bx = axes[0, col], axes[1, col]
        for m in c.METHODS:
            cv = curves[(key, m)]
            lw = 1.3 if m == "b3lyp" else 1.0
            y = cv["ecm"] if which == "energy" else cv["ecm"] - cv["eharm"]
            ax.plot(cv["wdr"], y / 1e3, color=color(m), lw=lw, zorder=3 if m == "b3lyp" else 2, label=LABEL[m])
            # curvature as local wavenumber, minus the method's own straight-line part (the cubic term,
            # within 3 % equal for all methods), so that the quartic bend and any wobble show
            yv = cv["w"] - cv["w0"]
            lin = np.polynomial.polynomial.polyfit(cv["dr"], yv, 1)
            bx.plot(cv["dr"], yv - np.polynomial.polynomial.polyval(cv["dr"], lin), color=color(m), lw=lw,
                    zorder=3 if m == "b3lyp" else 2)  # fmt: skip
        if which == "energy":
            b = curves[(key, "b3lyp")]
            ax.plot(b["wdr"], b["eharm"] / 1e3, color=style.INK_MUTED, lw=0.7, ls="--", zorder=1)
            ax.set_ylim(-0.5, YMAX[key])
        else:
            ax.axhline(0, color=style.INK_MUTED, lw=0.7, ls="--", zorder=1)
        # the Gaussian steps (B3LYP's conversion to dr; the models differ by < 2 %)
        rb = res[(res.key == key) & (res.method == "b3lyp")].iloc[0]
        for h in STEPS_Q:
            x = rb[f"dr_h{h}"]
            for s in (-1, 1):
                bx.axvline(s * x, color=style.INK if h == 0.01 else style.GRID,
                           lw=0.9 if h == 0.01 else 0.7, ls=":" if h == 0.01 else "-", zorder=0)  # fmt: skip
        tp = rb["zpe_turning_A"]
        for s in (-1, 1):
            ax.axvline(s * tp, color=style.GRID, lw=0.7, zorder=0)
        for a in (ax, bx):
            style.clean_axes(a, grid=False)
        ax.set_title(TITLE[key], fontsize=style.BASE_FS)
        ax.set_xlim(-0.26, 0.41)
        bx.set_xlim(-0.052, 0.052)
        bx.set_xticks([-0.04, 0, 0.04])
        bx.set_xlabel(r"$\Delta r$(X–H) / Å")
        ax.set_xlabel(r"$\Delta r$(X–H) / Å")
        ax.xaxis.set_tick_params(labelbottom=True)
    axes[0, 0].set_ylabel(("$E - E_{\\min}$" if which == "energy" else "$E - E_{\\mathrm{harm}}$") + r" / $10^3$ cm$^{-1}$")
    axes[1, 0].set_ylabel(r"curvature, $\tilde\nu_{\mathrm{loc}}$" + "\n" + r"minus linear part / cm$^{-1}$")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 1.01), handlelength=1.2,
               columnspacing=0.8, fontsize=style.ANNOT_FS)  # fmt: skip
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=1.2, w_pad=0.8)
    style.save(fig, out, "surface_scans" if which == "energy" else "surface_scans_anharmonic")


ZOOMS = [("glycine_nh", "mace_omol"), ("glycine_nh", "mace_polar"), ("glycine_nh", "mace_anicc"),
         ("methanol_ch", "mace_omol"), ("methanol_ch", "mace_polar"), ("methanol_oh", "mace_omol"),
         ("methanol_oh", "mace_polar"), ("methanol_oh", "mace_off")]  # fmt: skip


def zoom_figure(curves, res, out: Path):
    """Curvature on the 0.0002 A grid near the minimum, one panel per zoomed case, with the B3LYP
    curve of the same mode (0.002 A grid; smooth) and Gaussian's step points."""
    cases = [z for z in ZOOMS if curves[z]["zoom"]]
    style.apply_style()
    ncol = 4
    nrow = int(np.ceil(len(cases) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(style.figwidth(), 1.75 * nrow + 0.5), sharex=True)
    axes = np.atleast_2d(axes)
    for ax, (key, m) in zip(axes.flat, cases):
        cv, b = curves[(key, m)], curves[(key, "b3lyp")]
        sel = np.abs(b["dr"]) <= 0.0121
        # both relative to each method's own value at the minimum, straight-line part removed
        for x, y, col, lw in ((b["dr"][sel], b["w"][sel] - b["w0"], style.INK, 1.0),
                              (cv["zdr"], cv["zw"] - cv["w0"], color(m), 0.9)):  # fmt: skip
            lin = np.polynomial.polynomial.polyfit(x, y, 1)
            ax.plot(x, y - np.polynomial.polynomial.polyval(x, lin), color=col, lw=lw)
        rb = res[(res.key == key) & (res.method == m)].iloc[0]
        for h in (0.005, 0.01):
            for sgn in (-1, 1):
                ax.axvline(sgn * rb[f"dr_h{h}"], color=style.INK if h == 0.01 else style.GRID,
                           lw=0.8, ls=":" if h == 0.01 else "-", zorder=0)  # fmt: skip
        style.clean_axes(ax, grid=False)
        mol_, bond_ = TITLE[key].split(", ")
        ax.set_title(f"{style.tex(LABEL[m])}\n{mol_} {bond_.split()[0]}", fontsize=style.ANNOT_FS - 0.5, color=color(m))
        ax.set_xlim(-0.0125, 0.0125)
        ax.set_xticks([-0.01, 0, 0.01])
        ax.tick_params(labelsize=style.ANNOT_FS - 1)
    for ax in axes.flat[len(cases):]:
        ax.set_visible(False)
    for r in range(nrow):
        axes[r, 0].set_ylabel(r"$\tilde\nu_{\mathrm{loc}}$, lin. part" + "\n" + r"removed / cm$^{-1}$", fontsize=style.ANNOT_FS)
    fig.supxlabel(r"$\Delta r$(X–H) / Å", fontsize=style.BASE_FS)
    fig.tight_layout(w_pad=0.5, h_pad=0.8)
    style.save(fig, out, "surface_zoom")


def map2d(out: Path) -> dict:
    """Methanol O-H stretch x C-O-H bend: energy contours (top) and coupling part (bottom)."""
    from matplotlib.colors import TwoSlopeNorm

    to_cm = 1e-18 / (6.62607015e-34 * c.C_CM)
    data = {}
    for m in c.METHODS:
        p = c.DATA / f"map2d_{m}.json"
        if not p.exists():
            return {}
        d = c.load(p.name)
        e = (np.array(d["E"]) - np.array(d["E"])[10, 10]) * to_cm
        coup = e - e[:, [10]] - e[[10], :] + e[10, 10]
        data[m] = dict(dr=np.array(d["dr1"]), q2=np.array(d["q2"]), e=e, coup=coup, seconds=d["seconds"])
    lim = 600.0  # cm-1; the compressed, strongly bent corner reaches -9400 and is clipped (extend arrows)
    style.apply_style()
    fig, axes = plt.subplots(2, 5, figsize=(style.figwidth(), 3.6), sharex=True, sharey=True)
    norm = TwoSlopeNorm(0, -lim, lim)
    for col, m in enumerate(c.METHODS):
        x = data[m]
        X, Y = np.meshgrid(x["dr"], x["q2"], indexing="ij")
        ax, bx = axes[0, col], axes[1, col]
        ax.contour(X, Y, x["e"], levels=np.arange(1000, 20001, 1000), colors=color(m), linewidths=0.7)
        im = bx.contourf(X, Y, np.clip(x["coup"], -lim, lim), levels=np.linspace(-lim, lim, 21), cmap="RdBu_r",
                         norm=norm, extend="both")  # fmt: skip
        bx.contour(X, Y, x["coup"], levels=[0], colors=style.INK_MUTED, linewidths=0.5)
        ax.set_title(style.tex(LABEL[m]), color=color(m), fontsize=style.ANNOT_FS)
        for a in (ax, bx):
            a.set_xticks([-0.2, 0, 0.2])
            a.tick_params(labelsize=style.ANNOT_FS - 1)
    axes[0, 0].set_ylabel(r"$Q$(C–O–H bend)" + "\n" + r"/ amu$^{1/2}$ Å", fontsize=style.ANNOT_FS)
    axes[1, 0].set_ylabel(r"$Q$(C–O–H bend)" + "\n" + r"/ amu$^{1/2}$ Å", fontsize=style.ANNOT_FS)
    fig.supxlabel(r"$\Delta r$(O–H) / Å", fontsize=style.BASE_FS)
    fig.tight_layout(rect=(0, 0, 0.86, 1), w_pad=0.3, h_pad=0.6)
    cax = fig.add_axes([0.875, 0.13, 0.015, 0.36])
    cb = fig.colorbar(im, cax=cax, ticks=[-600, -300, 0, 300, 600], extend="both")
    cb.set_label(r"coupling / cm$^{-1}$", fontsize=style.ANNOT_FS)
    cb.ax.tick_params(labelsize=style.ANNOT_FS - 1)
    style.save(fig, out, "surface_2d_methanol")
    # numbers for the report: coupling at the corners of the grid
    summ = {}
    for m, x in data.items():
        cp = x["coup"]
        summ[m] = {"coup_min": float(cp.min()), "coup_max": float(cp.max()), "seconds": x["seconds"],
                   "corners (dr-,q2-) (dr-,q2+) (dr+,q2-) (dr+,q2+)": [round(float(cp[i, j]), 1) for i in (0, -1) for j in (0, -1)]}  # fmt: skip
    return summ


def tables(res: pd.DataFrame) -> str:
    lines = []
    for key, *_ in c.TARGETS:
        r = res[res.key == key].set_index("method").reindex(c.METHODS)
        lines += [f"\n### {TITLE[key]}\n",
                  "| method | ω̃ / cm-1 | χ_ii | VPT2 shift | cubic fit | cubic Q=0.010 | cubic Gaussian | quartic fit | quartic Q=0.005 | quartic Q=0.010 | quartic Q=0.020 | quartic Q=0.040 | quartic Gaussian | wobble / cm-1 | wobble scale / Å |",
                  "|---|" + "---|" * 14]  # fmt: skip
        for m, x in r.iterrows():
            lines.append(f"| {LABEL[m]} | {x.freq:.1f} | {x.chi_ii:.1f} | {x.vpt2_shift:.1f} | {x.cubic_fit:.2f} | "
                         f"{x['cubic_h0.01']:.2f} | {x.cubic_gauss:.2f} | {x.quartic_fit:.1f} | {x['quartic_h0.005']:.1f} | "
                         f"{x['quartic_h0.01']:.1f} | {x['quartic_h0.02']:.1f} | {x['quartic_h0.04']:.1f} | {x.quartic_gauss:.1f} | "
                         f"{x.wobble_amp_cm:.3f} | {x.wobble_scale_A:.3f} |")  # fmt: skip
    return "\n".join(lines)


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else c.HERE / "out"
    res, curves = analyse()
    res.to_csv(c.DATA / "results.csv", index=False)
    # estimate of the diagonal VPT2 term of chi_ii, (phi_iiii - 5 phi_iii^2 / 3 omega) / 16 in reduced
    # units (Gaussian's FI/k ratio of the run converts), at each step and from the smooth fit; "rest" =
    # benchmark chi_ii minus the term at the default step (couplings, held fixed: an approximation)
    est = []
    for _, x in res.iterrows():
        r4 = x.quartic_FI_gauss / x.quartic_gauss
        r3 = x.cubic_FI_gauss / x.cubic_gauss if abs(x.cubic_gauss) > 0.05 else 0.0
        diag = lambda c3, c4: (c4 * r4 - 5 / 3 * (c3 * r3) ** 2 / x.freq) / 16  # noqa: B023, E731
        d = {f"Q{h}": diag(x[f"cubic_h{h}"], x[f"quartic_h{h}"]) for h in STEPS_Q} | {"fit": diag(x.cubic_fit, x.quartic_fit)}
        est.append({"key": x.key, "method": x.method, "chi_ii": x.chi_ii, **d, "rest": x.chi_ii - d["Q0.01"]})
    pd.DataFrame(est).round(1).to_csv(c.DATA / "chi_diag_estimate.csv", index=False)
    (c.DATA / "tables.md").write_text(tables(res))
    print(tables(res))
    figure(curves, res, out, "energy")
    figure(curves, res, out, "anharmonic")
    zoom_figure(curves, res, out)
    summ = map2d(out)
    if summ:
        c.save("map2d_summary.json", summ)
        print(summ)
