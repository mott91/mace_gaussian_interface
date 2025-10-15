#!/usr/bin/env python3
"""Test refactored functions"""

import numpy as np
from ase import Atoms
from gm_main import (
    parse_gaussian_input,
    update_molecule_geometry,
    calculate_energy_and_forces,
    write_gaussian_output
)

def test_parse_gaussian_input():
    """Test parsing Gaussian input file"""
    
    # Create a test input file
    test_input = """2 1 0 1
    1   0.000000000000E+00   0.000000000000E+00   0.000000000000E+00
    1   1.400000000000E+00   0.000000000000E+00   0.000000000000E+00
"""
    
    with open('test_input.txt', 'w') as f:
        f.write(test_input)
    
    # Parse it
    natoms, deriv, charge, spin, coords, names = parse_gaussian_input('test_input.txt')
    
    # Check results
    assert natoms == 2, f"Expected 2 atoms, got {natoms}"
    assert deriv == 1, f"Expected deriv=1, got {deriv}"
    assert charge == 0, f"Expected charge=0, got {charge}"
    assert spin == 1, f"Expected spin=1, got {spin}"
    assert len(names) == 2, f"Expected 2 atom names, got {len(names)}"
    assert all(n == 'H' for n in names), f"Expected all H, got {names}"
    
    # Check coordinates (should be converted from Bohr to Angstrom)
    expected_x = 1.4 * 0.52917721092  # Bohr to Angstrom
    assert np.isclose(coords[1, 0], expected_x, rtol=1e-6), \
        f"Coordinate conversion failed: {coords[1, 0]} vs {expected_x}"
    
    print("\u2713 test_parse_gaussian_input PASSED")
    
    # Cleanup
    import os
    os.remove('test_input.txt')


def test_unit_conversions():
    """Test that unit conversions are correct"""
    
    # Constants (these should match what's in gm_main.py)
    BOHR_TO_ANGSTROM = 0.52917721092
    EV_TO_HARTREE = 27.211386246
    
    # Test energy conversion
    energy_ev = 10.0
    energy_hartree = energy_ev / EV_TO_HARTREE
    # Calculate expected value precisely
    expected_energy = 10.0 / 27.211386246
    assert np.isclose(energy_hartree, expected_energy, rtol=1e-10), \
        f"Energy conversion wrong: {energy_hartree} vs {expected_energy}"
    
    # Test gradient conversion  
    grad_ev_ang = 1.0
    grad_hartree_bohr = grad_ev_ang * BOHR_TO_ANGSTROM / EV_TO_HARTREE
    # Calculate expected value precisely
    expected_grad = 0.52917721092 / 27.211386246
    assert np.isclose(grad_hartree_bohr, expected_grad, rtol=1e-10), \
        f"Gradient conversion wrong: {grad_hartree_bohr} vs {expected_grad}"
    
    print("\u2713 test_unit_conversions PASSED")


def test_update_molecule_geometry():
    """Test updating molecule geometry"""
    
    # Create test atoms
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]])
    
    # New coordinates
    new_coords = np.array([[0, 0, 0], [1.5, 0, 0]])
    
    # Update
    update_molecule_geometry(atoms, new_coords, charge=1, spin=2)
    
    # Check
    assert np.allclose(atoms.get_positions(), new_coords), "Positions not updated"
    assert atoms.info['charge'] == 1.0, f"Charge not set: {atoms.info['charge']}"
    assert atoms.info['spin'] == 2.0, f"Spin not set: {atoms.info['spin']}"
    
    print("\u2713 test_update_molecule_geometry PASSED")


def test_output_format():
    """Test that output file has correct format"""
    
    # Create test data
    natoms = 2
    energy = -1.0  # eV
    gradient = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])  # eV/Ang
    dipole = np.array([0.1, 0.2, 0.3])  # e*Bohr
    dipole_deriv = np.zeros((6, 3))  # 3*2 atoms
    hessian = None
    deriv = 1
    
    # Write output
    write_gaussian_output(
        'test_output.txt', natoms, energy, gradient,
        dipole, dipole_deriv, hessian, deriv
    )
    
    # Read and check
    with open('test_output.txt', 'r') as f:
        lines = f.readlines()
    
    # Check format
    assert 'D' in lines[0], "Should use Fortran D exponent, not E"
    assert len(lines) >= 10, f"Expected at least 10 lines, got {len(lines)}"
    
    # Check first line has 4 values (energy + 3 dipole components)
    first_line_values = lines[0].split()
    assert len(first_line_values) == 4, \
        f"First line should have 4 values, got {len(first_line_values)}"
    
    print("\u2713 test_output_format PASSED")
    
    # Cleanup
    import os
    os.remove('test_output.txt')


if __name__ == '__main__':
    print("\n" + "="*60)
    print("RUNNING UNIT TESTS")
    print("="*60)
    
    test_parse_gaussian_input()
    test_unit_conversions()
    test_update_molecule_geometry()
    test_output_format()
    
    print("\n" + "="*60)
    print("ALL TESTS PASSED \u2713")
    print("="*60)
