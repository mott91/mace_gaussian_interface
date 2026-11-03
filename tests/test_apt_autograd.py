"""Autograd dipole derivatives of every campaign dipole model match finite differences.

Needs a CUDA GPU and the model files, so it only runs with MACE_RUN_GPU_TESTS=1:

    MACE_RUN_GPU_TESTS=1 pytest tests/test_apt_autograd.py

Full table (three molecules, three step sizes): scripts/apt_autograd_vs_fd.py.
"""

import os

import numpy as np
import pytest
from ase.build import molecule

pytestmark = [
    pytest.mark.gpu,
    pytest.mark.slow,
    pytest.mark.skipif(
        os.getenv("MACE_RUN_GPU_TESTS") != "1", reason="set MACE_RUN_GPU_TESTS=1 to run"
    ),
]


@pytest.mark.parametrize("model", ["mace_ml", "mace_mdp", "mace_polar1"])
def test_autograd_apt_matches_finite_differences(model):
    from mace_gaussian.calculators import dipole_factory
    from mace_gaussian.calculators.base import DipoleCalculatorBase

    calc = dipole_factory.get_calculator(model)
    atoms = molecule("CO2")
    atoms.info["charge"], atoms.info["spin"] = 0.0, 1.0
    auto = np.asarray(calc.calculate_dipole_derivatives(atoms.copy()))
    fd = np.asarray(
        DipoleCalculatorBase.calculate_dipole_derivatives(calc, atoms.copy(), displacement=1e-3)
    )
    assert np.abs(fd - auto).max() < 1e-3  # e; measured 1e-4 or less (2026-11-03)
    # neutral molecule: the atomic polar tensors sum to zero
    assert np.abs(auto.reshape(-1, 3, 3).sum(axis=0)).max() < 1e-8
