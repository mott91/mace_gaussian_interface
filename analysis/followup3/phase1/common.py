import numpy as np, pandas as pd
D = "/home/mot/mace_gaussian/campaigns/2026/figures/supervisor/data/"
B = pd.read_csv(D + "bands.csv"); C = pd.read_csv(D + "chi.csv")
MODELS = ["mace_anicc", "mace_omol", "mace_off", "mace_polar"]
SHORT = dict(mace_anicc="ANI-cc", mace_omol="OMOL", mace_off="OFF", mace_polar="POLAR")
SYMTOP = {"ammonia", "phosphine"}
def fund(with_nh3=False):
    f = B[(B.band_type == "fundamental") & (B.dipole_model == "mace_ml")].copy()
    f = f[f.same_conformer.astype(str) != "False"]
    if not with_nh3: f = f[~f.molecule.isin(SYMTOP)]
    f["err"] = f.ml_anharmonic - f.dft_anharmonic
    f["herr"] = f.ml_harmonic - f.dft_harmonic
    return f
def st(e):
    e = np.asarray(e, float); a = np.abs(e)
    return dict(n=len(e), mae=a.mean(), med=np.median(a), mse=e.mean(), rmse=np.sqrt((e**2).mean()), mx=a.max(), gt100=100*(a > 100).mean())
