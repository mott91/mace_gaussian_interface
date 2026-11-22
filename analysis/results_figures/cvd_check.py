"""Item 1a of figures-ch4-final: CAM02-UCS distances (colorspacious) of candidate reds to the colours in use,
normal vision and simulated deuteranopia/protanopia (colorspacious severity 100 / daltonlens Machado 2009).
Needs colorspacious and daltonlens: run in a throwaway venv, never in the mace4ir envs.
"""

import numpy as np
from colorspacious import cspace_convert, deltaE
from daltonlens import simulate

others = {
    "OFF23 green": "#1baf7a",
    "ANICC pink": "#e87ba4",
    "OMOL blue": "#2a78d6",
    "B3LYP ink": "#1a1a1a",
    "muted ink": "#6b6b6b",
}
cands = ["#c0392b", "#b2182b", "#c62828", "#d62728", "#a93226", "#e03a2f"]


def rgb(h):
    return np.array([int(h[i : i + 2], 16) for i in (1, 3, 5)], float) / 255


def sim_cs(c, kind):
    return np.clip(
        cspace_convert(c, {"name": "sRGB1+CVD", "cvd_type": kind, "severity": 100}, "sRGB1"), 0, 1
    )


mach = simulate.Simulator_Machado2009()


def sim_ml(c, kind):
    im = (np.array(c)[None, None, :] * 255).astype(np.uint8)
    d = {"deuteranomaly": simulate.Deficiency.DEUTAN, "protanomaly": simulate.Deficiency.PROTAN}[
        kind
    ]
    return mach.simulate_cvd(im, d, severity=1.0)[0, 0] / 255


def J(c):
    return cspace_convert(c, "sRGB1", "CAM02-UCS")[0]


for cand in cands:
    print(f"\n{cand}  J'={J(rgb(cand)):.1f}")
    for name, o in others.items():
        a, b = rgb(cand), rgb(o)
        row = [deltaE(a, b, input_space="sRGB1", uniform_space="CAM02-UCS")]
        for k in ("deuteranomaly", "protanomaly"):
            row.append(
                deltaE(sim_cs(a, k), sim_cs(b, k), input_space="sRGB1", uniform_space="CAM02-UCS")
            )
            row.append(
                deltaE(sim_ml(a, k), sim_ml(b, k), input_space="sRGB1", uniform_space="CAM02-UCS")
            )
        print(
            f"  vs {name:12s} J'={J(b):5.1f}  normal {row[0]:5.1f}  deut {row[1]:5.1f}/{row[2]:5.1f}  prot {row[3]:5.1f}/{row[4]:5.1f}"
        )
