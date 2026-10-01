"""The DFT baseline route line carries the frozen protocol keywords."""

from ase.build import molecule

from mace_gaussian.dft_baseline import create_gaussian_dft_input


def _route(tmp_path, **kw):
    create_gaussian_dft_input(molecule("H2O"), "dft.gjf", "b3lyp", "6-31G(d,p)",
                              output_dir=tmp_path, **kw)  # fmt: skip
    return next(
        line for line in (tmp_path / "dft.gjf").read_text().splitlines() if line.startswith("#")
    )


def test_baseline_route_has_protocol_keywords(tmp_path):
    route = _route(tmp_path)
    assert route.startswith("# opt=verytight freq(anharm) b3lyp/6-31G(d,p)")
    assert "nosymm" in route.split()
    assert "int=ultrafine" in route.split()


def test_no_opt_keeps_nosymm(tmp_path):
    route = _route(tmp_path, optimize=False)
    assert "opt" not in route.split()[1]
    assert "nosymm" in route.split()
