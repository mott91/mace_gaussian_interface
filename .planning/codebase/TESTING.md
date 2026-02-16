# Testing Patterns

**Analysis Date:** 2026-02-16

## Test Framework

**Runner:**
- `pytest` (version >= 7.0.0)
- Config: `mace_ML_pkg/setup.cfg`, `mace_dipole_pkg/setup.cfg` (for subpackage tests)

**Assertion Library:**
- `pytest` built-in assertions
- `numpy.testing` for numerical comparisons (`.numpy_testing.assert_allclose()`, `.numpy_testing.assert_array_equal()`)

**Run Commands:**
```bash
# Run all tests in codebase
pytest

# Run specific test file
pytest testingStuff/test_refactoring.py

# Run tests with verbose output
pytest -v

# Run with coverage
pytest --cov=. --cov-report=html

# Run specific test function
pytest testingStuff/test_refactoring.py::test_parse_gaussian_input

# Watch mode (if pytest-watch installed)
pytest-watch
```

## Test File Organization

**Location:**
- Main project tests: `testingStuff/` directory (co-located with source)
- Subpackage tests: `mace_ML_pkg/tests/`, `mace_dipole_pkg/tests/`
- Test files for external modules use pytest fixtures from those modules

**Naming:**
- Test files: `test_*.py` (e.g., `test_refactoring.py`)
- Test functions: `test_<functionality>()` (e.g., `test_parse_gaussian_input()`)
- Test classes: `Test<Functionality>` (if organizing tests into classes)

**Structure:**
```
mace_gaussian/
├── testingStuff/
│   ├── test_refactoring.py          # Unit tests for refactored functions
│   ├── test_*.py                    # Other test modules
│   └── log_analyzer.py              # Utility scripts (not tests)
├── mace_ML_pkg/tests/               # MACE package tests
│   ├── test_calculator.py
│   ├── test_models.py
│   └── ...
└── mace_dipole_pkg/tests/           # Dipole package tests
    └── ...
```

## Test Structure

**Suite Organization:**
```python
def test_parse_gaussian_input():
    """Test parsing Gaussian input file"""

    # Setup: Create test data
    test_input = """2 1 0 1
    1   0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
    1   1.400000000000E+00   0.000000000000E+00   0.000000000000E+00
"""
    with open('test_input.txt', 'w') as f:
        f.write(test_input)

    # Execution: Call function under test
    natoms, deriv, charge, spin, coords, names = parse_gaussian_input('test_input.txt')

    # Assertions: Verify results
    assert natoms == 2, f"Expected 2 atoms, got {natoms}"
    assert deriv == 1, f"Expected deriv=1, got {deriv}"
    assert len(names) == 2, f"Expected 2 atom names, got {len(names)}"
    assert all(n == 'H' for n in names), f"Expected all H, got {names}"

    # Cleanup: Remove test files
    import os
    os.remove('test_input.txt')
```

**Patterns:**

- **Setup**: Test data creation before calling function
- **Execution**: Single function call per test (focused)
- **Assertions**: Multiple assertions allowed if testing single aspect
- **Teardown**: Cleanup in reverse order (files removed, state reset)

## Mocking

**Framework:**
- `unittest.mock` (standard library) - not heavily used in current tests
- Tests prefer actual object instantiation and direct calls over mocks
- Monkeypatching used for module swapping in `mace_calculators.py`

**Patterns (if needed):**
```python
# Pattern: Monkeypatch system modules
import sys
original_module = sys.modules.get("mace.modules.models")

def setup():
    sys.modules["mace.modules.models"] = fake_module

def teardown():
    sys.modules["mace.modules.models"] = original_module

# Pattern: Direct instantiation preferred over mocking
atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])
calc = EspalomaDipoleCalculator()
dipole, charges = calc.calculate_dipole(atoms)
```

**What to Mock:**
- External calculators that require GPU (use CPU versions for testing)
- File system operations (use `tmp_path` pytest fixture)
- Network calls (if any added in future)

**What NOT to Mock:**
- MACE/ML calculators (test with actual models if available)
- Core parsing logic (test with real Gaussian output samples)
- Data transformations (test with actual coordinate arrays)

## Fixtures and Factories

**Test Data:**
```python
# Pattern 1: Inline fixtures for simple cases
@pytest.fixture
def simple_atoms():
    """Create a simple H2 molecule"""
    return Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])

# Pattern 2: Parametrized fixtures for multiple scenarios
@pytest.fixture(params=['H2', 'H2O', 'CH4'])
def molecules(request):
    """Test multiple molecule types"""
    if request.param == 'H2':
        return Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])
    # ... etc

# Pattern 3: Module-scoped fixtures for expensive setup (from MACE tests)
@pytest.fixture(scope="module", name="trained_model")
def trained_model_fixture(tmp_path_factory, fitting_configs):
    """Train a model once per test module"""
    # Expensive training code
    return trained_model
```

**Location:**
- Fixtures used by multiple tests: defined in `conftest.py` (if it exists)
- Single-test fixtures: defined in same file as test
- Shared test data: `mace_ML_pkg/tests/conftest.py` provides MACE-specific fixtures

## Coverage

**Requirements:** Not enforced (no minimum set)

**View Coverage:**
```bash
# Generate coverage report
pytest --cov=. --cov-report=html

# View in browser
open htmlcov/index.html
```

**Current Coverage:**
- Unit tests for core functions: `test_refactoring.py` covers parsing, geometry updates, output writing
- Integration tests: implicit through full-workflow scripts
- MACE package tests: extensive in `mace_ML_pkg/tests/`

## Test Types

**Unit Tests:**
- Scope: Individual functions (parsing, conversions, geometry updates)
- Approach: Direct function calls with simple inputs
- Examples: `test_parse_gaussian_input()`, `test_unit_conversions()`, `test_output_format()`
- Location: `testingStuff/test_refactoring.py`

**Integration Tests:**
- Scope: Multi-step workflows (optimization + frequency + analysis)
- Approach: Run through multiple stages, verify intermediate outputs
- Implicit through: `run_analysis.py`, `comparison_workflow.py`
- Can be run manually: `python run_analysis.py <molecule>` with assertions in logs

**E2E Tests:**
- Framework: Not explicitly used
- Could be implemented as: Full workflow scripts with saved output verification
- Currently tested through: Manual runs of `python cli.py run water.xyz`

## Common Patterns

**Async Testing:**
- Not applicable (no async code in codebase)

**Error Testing:**
```python
# Pattern 1: Test exception is raised
def test_missing_log_file():
    """Test that FileNotFoundError is raised for missing file"""
    with pytest.raises(FileNotFoundError):
        GaussianLogParser("nonexistent_file.log")

# Pattern 2: Test fallback behavior
def test_dipole_calculation_fallback():
    """Test graceful degradation when dipole calc fails"""
    # Mock calculator failure
    bad_calc = BrokenDipoleCalculator()

    # Should return zero dipole, not raise
    dipole, dipole_deriv, charges = calculate_dipole_properties(
        atoms, bad_calc, deriv=1, calculate_derivatives=True
    )

    assert np.allclose(dipole, [0, 0, 0]), "Should return zero dipole on failure"
    assert dipole_deriv.shape == (3*natoms, 3), "Should return zero derivatives"
```

**Numerical Testing:**
```python
# Pattern: Use numpy/pytest for floating point comparisons
def test_unit_conversions():
    """Test that unit conversions are correct"""

    BOHR_TO_ANGSTROM = 0.52917721092
    EV_TO_HARTREE = 27.211386246

    # Test with high precision
    energy_ev = 10.0
    energy_hartree = energy_ev / EV_TO_HARTREE
    expected_energy = 10.0 / 27.211386246

    # Use isclose for floating point comparison
    assert np.isclose(energy_hartree, expected_energy, rtol=1e-10), \
        f"Energy conversion wrong: {energy_hartree} vs {expected_energy}"

# Pattern: Array comparison
def test_update_molecule_geometry():
    """Test updating molecule geometry"""
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])
    new_coords = np.array([[0, 0, 0], [1.5, 0, 0]])

    update_molecule_geometry(atoms, new_coords, charge=1, spin=2)

    # Use allclose for array comparison
    assert np.allclose(atoms.get_positions(), new_coords), "Positions not updated"
    assert atoms.info['charge'] == 1.0, f"Charge not set"
    assert atoms.info['spin'] == 2.0, f"Spin not set"
```

**Output Validation:**
```python
# Pattern: Test file format and content
def test_output_format():
    """Test that output file has correct format"""

    natoms = 2
    energy = -1.0  # eV
    gradient = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    dipole = np.array([0.1, 0.2, 0.3])
    dipole_deriv = np.zeros((6, 3))
    hessian = None
    deriv = 1

    # Write output
    write_gaussian_output(
        'test_output.txt', natoms, energy, gradient,
        dipole, dipole_deriv, hessian, deriv
    )

    # Verify format
    with open('test_output.txt', 'r') as f:
        lines = f.readlines()

    assert 'D' in lines[0], "Should use Fortran D exponent, not E"
    assert len(lines) >= 10, f"Expected at least 10 lines, got {len(lines)}"

    first_line_values = lines[0].split()
    assert len(first_line_values) == 4, f"First line should have 4 values"

    # Cleanup
    import os
    os.remove('test_output.txt')
```

## Test Execution

**Running Tests:**
```bash
# Quick test run
python testingStuff/test_refactoring.py

# With pytest
pytest testingStuff/test_refactoring.py -v

# Specific test
pytest testingStuff/test_refactoring.py::test_parse_gaussian_input -v
```

**Test Entry Point:**
Tests in `testingStuff/test_refactoring.py` have a `if __name__ == '__main__':` block that allows standalone execution:

```python
if __name__ == '__main__':
    print("\n" + "="*60)
    print("RUNNING UNIT TESTS")
    print("="*60)

    test_parse_gaussian_input()
    test_unit_conversions()
    test_update_molecule_geometry()
    test_output_format()

    print("\n" + "="*60)
    print("ALL TESTS PASSED ✓")
    print("="*60)
```

## Test Coverage Areas

**Well-Tested:**
- Core parsing functions: `parse_gaussian_input()`, `GaussianLogParser`
- Unit conversions: Bohr ↔ Angstrom, eV ↔ Hartree
- Geometry updates: `update_molecule_geometry()`
- Output formatting: `write_gaussian_output()` produces valid Gaussian-compatible output

**Partially Tested:**
- MACE calculator integration (tested in `mace_ML_pkg/tests/`)
- Dipole calculations (availability checks, but not full workflows)
- File I/O (some integration tests in workflow)

**Under-Tested:**
- Mode matching algorithm (`mode_matching.py`)
- Spectral analysis statistics (`analyze_spectra.py`)
- HTML report generation (`html_report_generator.py`)
- Charge analysis visualization (`charge_analysis.py`)

**Gap:** Integration tests for full workflows (optimization → frequency → analysis) rely on manual execution of `python cli.py run` and examining outputs.

---

*Testing analysis: 2026-02-16*
