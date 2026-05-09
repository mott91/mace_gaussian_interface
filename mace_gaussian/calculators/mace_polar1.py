"""MACE-POLAR-1 dipole calculator.

Uses the MACE-POLAR-1 foundation model (OMol25 / ωB97M-V) and pulls the
molecular dipole as a byproduct of the multipole charge-density expansion in
its forward pass — no separate dipole-only model needed.

Independent of the energy-side ``mace_polar`` calculator: this loads its own
copy of the model. Sharing a single instance across energy + dipole calls is
a deferred optimization (see project memory ``project_model_sharing_idea``).
"""

from __future__ import annotations

import logging

import numpy as np
import torch

from .base import DipoleCalculatorBase

logger = logging.getLogger(__name__)


class MACEPolar1DipoleCalculator(DipoleCalculatorBase):
    """Molecular dipole from MACE-POLAR-1.

    Notes
    -----
    PolarMACE's forward pass returns ``"dipole"`` as a side output of its
    multipole expansion. The dedicated MACE ``get_dielectric_derivatives()``
    API refuses to handle PolarMACE (only ``DipoleMACE`` /
    ``DipolePolarizabilityMACE``), so for derivatives we either run our own
    autograd loop or fall back to base-class central-difference.
    """

    use_autograd: bool = True

    def __init__(self, model: str = "polar-1-l", device: str | None = None):
        self.model_name = model
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.calc = None
        super().__init__("mace_polar1")

    def _check_availability(self):
        try:
            from mace.calculators import mace_polar  # noqa: F401

            self.available = True
            logger.info("✓ MACE-POLAR-1 dipole calculator available")
        except ImportError as e:
            self.available = False
            logger.warning(f"✗ MACE-POLAR-1 dipole calculator failed: {e}")

    def _ensure_calculator(self):
        if self.calc is None:
            from mace.calculators import mace_polar

            self.calc = mace_polar(
                model=self.model_name,
                device=self.device,
                default_dtype="float64",
            )
            logger.info("Loaded MACE-POLAR-1 (%s) on %s", self.model_name, self.device)

    def calculate_dipole(self, atoms, **kwargs):
        if not self.available:
            raise RuntimeError("MACE-POLAR-1 dipole calculator not available")
        self._ensure_calculator()

        atoms_copy = atoms.copy()
        atoms_copy.calc = self.calc
        atoms_copy.get_potential_energy()  # populates calc.results

        dipole = self.calc.results.get("dipole")
        if dipole is None:
            raise RuntimeError(
                "PolarMACE did not populate results['dipole']. "
                f"Available keys: {sorted(self.calc.results.keys())}"
            )
        if isinstance(dipole, torch.Tensor):
            dipole = dipole.detach().cpu().numpy()
        dipole = np.asarray(dipole).reshape(-1)[:3]
        return dipole, None

    def calculate_dipole_derivatives(self, atoms, **kwargs) -> np.ndarray:
        """Dipole derivatives ∂μ/∂r in shape (3*N_atoms, 3), units e/Å.

        Tries an autograd Jacobian (3 backward passes, ~20× faster for typical
        molecules). Falls back to base-class finite differences if autograd
        fails — this is the safety net since the MACE library does not
        officially support derivative calculation for PolarMACE.
        """
        if not self.use_autograd:
            return DipoleCalculatorBase.calculate_dipole_derivatives(self, atoms, **kwargs)
        try:
            return self._autograd_dipole_derivatives(atoms)
        except Exception as e:  # broad on purpose: any autograd issue → fallback
            logger.warning(
                "Autograd dipole derivatives failed (%s); falling back to finite differences",
                e,
            )
            return DipoleCalculatorBase.calculate_dipole_derivatives(self, atoms, **kwargs)

    def _autograd_dipole_derivatives(self, atoms) -> np.ndarray:
        self._ensure_calculator()
        n_atoms = len(atoms)

        batch = self.calc._atoms_to_batch(atoms)
        batch_dict = self.calc._clone_batch(batch).to_dict()

        positions = batch_dict["positions"]
        positions.requires_grad_(True)

        model = self.calc.models[0]
        output = model(batch_dict, training=True)
        dipole = output["dipole"].reshape(-1)[:3]

        dmu_dr = torch.zeros(
            (3, n_atoms, 3), device=positions.device, dtype=positions.dtype
        )
        for i in range(3):
            grad = torch.autograd.grad(
                dipole[i], positions, retain_graph=(i < 2), create_graph=False
            )[0]
            dmu_dr[i] = grad

        # Base-class contract: shape (3*N, 3), atom-major Cartesian-minor; dipole is the column
        return dmu_dr.detach().cpu().numpy().transpose(1, 2, 0).reshape(3 * n_atoms, 3)
