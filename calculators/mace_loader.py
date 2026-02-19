"""Safe MACE dipole model loading without sys.modules monkey-patching.

The dipole model file (dipole_model/model_1.model) was saved with the class path
``mace.modules.models.AtomicDielectricMACE``, but the correct implementation lives
in the local dipole fork at ``mace_dipole_core.modules.models``.  The standard
mace-torch package defines its own incompatible ``AtomicDielectricMACE`` class, so
simply importing at that path would load the wrong ``forward()`` logic.

Previous approach: swap ``sys.modules["mace.modules.models"]`` with a shallow copy
of the dipole fork module before every ``torch.load``, then clean up afterwards.

New approach: use ``torch.load(pickle_module=...)`` to intercept class resolution
*only* during deserialization.  A custom ``Unpickler.find_class`` redirects lookups
for ``mace.modules.models.*`` to ``mace_dipole_core.modules.models.*``.  No global
state is mutated, no cleanup is needed, and standard MACE imports are unaffected.
"""

from __future__ import annotations

import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import torch

logger = logging.getLogger(__name__)

# Path to mace_dipole_pkg/ for lazy sys.path setup
_DIPOLE_PKG_DIR = str(Path(__file__).resolve().parent.parent / "mace_dipole_pkg")


def _ensure_dipole_importable() -> None:
    """Add ``mace_dipole_pkg/`` to ``sys.path`` if not already present.

    The dipole fork package (``mace_dipole_core``) lives under
    ``mace_dipole_pkg/`` which is not on ``sys.path`` by default.
    This function makes it importable lazily on first use.
    """
    if _DIPOLE_PKG_DIR not in sys.path:
        sys.path.insert(0, _DIPOLE_PKG_DIR)


class _DipoleModelUnpickler(pickle.Unpickler):
    """Unpickler that remaps ``mace.modules.models`` to ``mace_dipole_core`` equivalents.

    When ``torch.load`` deserializes the dipole model, pickle calls
    ``find_class("mace.modules.models", "AtomicDielectricMACE")``.  This
    override intercepts that call and returns the class from the dipole fork
    instead, without touching ``sys.modules``.
    """

    def find_class(self, module: str, name: str):
        if module == "mace.modules.models":
            _ensure_dipole_importable()
            import mace_dipole_core.modules.models as dipole_models

            cls = getattr(dipole_models, name, None)
            if cls is not None:
                return cls
        return super().find_class(module, name)


class _DipolePickleModule:
    """Drop-in ``pickle_module`` replacement for ``torch.load``.

    Provides ``Unpickler`` (our remapping version) and delegates every other
    attribute to the stdlib ``pickle`` module so that ``torch.load`` can use
    the rest of the pickle API transparently.
    """

    Unpickler = _DipoleModelUnpickler

    def __getattr__(self, name: str):
        return getattr(pickle, name)


def load_dipole_calculator(
    model_path: str,
    device: str = "cuda",
    model_type: str = "DipolePolarizabilityMACE",
    default_dtype: str = "float64",
):
    """Load a MACE dipole calculator with safe class remapping.

    Temporarily patches ``torch.load`` on the dipole calculator module to
    inject ``pickle_module=_DipolePickleModule()`` and ``weights_only=False``,
    then restores the original immediately after construction.

    Parameters
    ----------
    model_path : str
        Path to the saved dipole model file.
    device : str
        Torch device string (``"cuda"`` or ``"cpu"``).
    model_type : str
        MACE model type identifier.
    default_dtype : str
        Default floating-point dtype.

    Returns
    -------
    MACECalculator_dipole
        ASE calculator instance ready for dipole calculations.
    """
    _ensure_dipole_importable()
    import mace_dipole_core.calculators.mace as dipole_calc_module
    from mace_dipole_core.calculators.mace import MACECalculator_dipole

    # Temporarily redirect torch.load in the dipole calculator module
    # to use our class-remapping pickle module
    original_load = dipole_calc_module.torch.load

    def _safe_load(*args, **kwargs):
        kwargs["pickle_module"] = _DipolePickleModule()
        kwargs["weights_only"] = False
        return original_load(*args, **kwargs)

    try:
        dipole_calc_module.torch.load = _safe_load
        calc = MACECalculator_dipole(
            model_paths=model_path,
            device=device,
            model_type=model_type,
            default_dtype=default_dtype,
        )
    finally:
        dipole_calc_module.torch.load = original_load

    logger.info("Loaded MACE dipole calculator from %s on %s", model_path, device)
    return calc


class MACEDipoleCalculator:
    """Wrapper for MACE dipole calculator with safe loading.

    Provides lazy initialization -- the underlying ASE calculator is not
    created until the first ``calculate_dipole`` call.  No
    ``cleanup_mace_modules()`` is needed because the safe loader does not
    mutate ``sys.modules``.
    """

    def __init__(self, model_path: str, device: str = "cuda") -> None:
        self.model_path = model_path
        self.device = device
        self.calc = None

    def _ensure_calculator(self) -> None:
        if self.calc is None:
            self.calc = load_dipole_calculator(
                model_path=self.model_path,
                device=self.device,
            )

    def calculate_dipole(self, atoms, **kwargs):
        """Calculate dipole moment using MACE dipole model.

        Parameters
        ----------
        atoms : ase.Atoms
            Molecular structure to compute the dipole for.

        Returns
        -------
        tuple[np.ndarray, None]
            ``(dipole_vector, None)`` where dipole_vector has shape ``(3,)``
            in units of e*Bohr.  Second element reserved for partial charges
            (not yet implemented).
        """
        self._ensure_calculator()

        atoms_copy = atoms.copy()
        atoms_copy.calc = self.calc

        try:
            dipole_moment = atoms_copy.get_dipole_moment()

            logger.debug("dipole_moment type: %s", type(dipole_moment))
            logger.debug("dipole_moment value: %s", dipole_moment)

            # Handle different return formats
            if isinstance(dipole_moment, tuple):
                logger.debug("Dipole returned as tuple, extracting first element")
                dipole_moment = dipole_moment[0]

            # Convert torch tensor to numpy if needed
            if isinstance(dipole_moment, torch.Tensor):
                dipole_moment = dipole_moment.detach().cpu().numpy()

            # Ensure it's a 1D numpy array
            dipole_moment = np.atleast_1d(np.array(dipole_moment))

            # Normalize shape to (3,)
            if dipole_moment.shape == (1, 3):
                dipole_moment = dipole_moment.flatten()
            elif dipole_moment.shape != (3,):
                logger.warning(
                    "Unexpected dipole shape: %s, attempting to reshape",
                    dipole_moment.shape,
                )
                dipole_moment = dipole_moment.reshape(-1)[:3]

            logger.debug("MACE dipole (final): %s, shape: %s", dipole_moment, dipole_moment.shape)

            return dipole_moment, None

        except Exception as e:
            logger.error("Error in MACE dipole calculation: %s", e)
            logger.error("Full traceback:", exc_info=True)
            raise
