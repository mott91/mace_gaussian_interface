"""PubChem PUG REST API client for fetching 3D molecular structures."""

from io import StringIO
from pathlib import Path

import requests
from ase.io import read, write

PUBCHEM_3D_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/record/SDF"
    "?record_type=3d"
)


def fetch_3d_structure(molecule_name: str, output_dir: Path, force: bool = False) -> Path:
    """Fetch 3D structure from PubChem and save as XYZ.

    Parameters
    ----------
    molecule_name : str
        Common name (e.g., "aspirin", "water")
    output_dir : Path
        Directory to save XYZ file (created if needed)
    force : bool
        If True, overwrite existing file

    Returns
    -------
    Path
        Path to the saved XYZ file

    Raises
    ------
    FileExistsError
        If file exists and force=False
    ValueError
        If molecule not found or no 3D conformer available
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{molecule_name}.xyz"

    if output_path.exists() and not force:
        raise FileExistsError(
            f"{molecule_name}.xyz already exists, skipping. Use --force to overwrite."
        )

    url = PUBCHEM_3D_URL.format(name=molecule_name)
    response = requests.get(url, timeout=30)

    if response.status_code == 404:
        body = response.text.strip()
        if "No CID found" in body:
            raise ValueError(
                f"Molecule '{molecule_name}' not found on PubChem. "
                f"Check the spelling or try an alternative name."
            )
        else:
            raise ValueError(
                f"No 3D structure available for '{molecule_name}' on PubChem. "
                f"The molecule exists but has no pre-computed 3D conformer."
            )
    response.raise_for_status()

    atoms = read(StringIO(response.text), format="sdf")
    write(str(output_path), atoms, format="xyz")
    return output_path
