"""
Diagnostic utilities for MACE-Gaussian interface.

This module provides functions to diagnose the Python environment
and test the availability of optional dependencies like espaloma_charge.
"""

import site
import sys


def diagnose_python_environment():
    """Diagnose Python environment and package availability"""

    print("=" * 60)
    print("PYTHON ENVIRONMENT DIAGNOSTICS")
    print("=" * 60)

    print(f"Python executable: {sys.executable}")
    print(f"Python version: {sys.version}")
    print(f"Python path: {sys.path[:3]}...")  # First few entries

    print("\nSite packages:")
    for path in site.getsitepackages():
        print(f"  {path}")

    print(f"\nUser site packages: {site.getusersitepackages()}")

    print("\nTesting imports:")

    # Test individual imports
    test_imports = [
        ('numpy', 'np'),
        ('ase', None),
        ('rdkit', None),
        ('rdkit.Chem', 'Chem'),
        ('espaloma_charge', None),
    ]

    for module_name, alias in test_imports:
        try:
            if alias:
                exec(f"import {module_name} as {alias}")
            else:
                exec(f"import {module_name}")
            print(f"  \u2713 {module_name}")

            # Get version if available
            try:
                if module_name == 'rdkit':
                    from rdkit import rdBase
                    print(f"    Version: {rdBase.rdkitVersion}")
                elif module_name == 'espaloma_charge':
                    import espaloma_charge
                    if hasattr(espaloma_charge, '__version__'):
                        print(f"    Version: {espaloma_charge.__version__}")
                    else:
                        print(f"    Location: {espaloma_charge.__file__}")
            except Exception:
                pass

        except ImportError as e:
            print(f"  \u2717 {module_name}: {e}")
        except Exception as e:
            print(f"  ? {module_name}: {e}")

    print("=" * 60)
    print()


def test_espaloma_functionality():
    """Test espaloma_charge functionality step by step"""
    print("TESTING ESPALOMA_CHARGE FUNCTIONALITY")
    print("=" * 60)

    try:
        print("1. Testing espaloma_charge import...")
        import espaloma_charge
        print("   \u2713 espaloma_charge imported successfully")

        print("2. Testing charge function import...")
        from espaloma_charge import charge
        print("   \u2713 charge function imported successfully")

        print("3. Testing RDKit import...")
        from rdkit import Chem
        print("   \u2713 RDKit imported successfully")

        print("4. Testing simple molecule creation...")
        mol = Chem.MolFromSmiles("N#N")
        if mol is not None:
            print("   \u2713 RDKit molecule created successfully")
        else:
            print("   \u2717 RDKit molecule creation failed")
            return False

        print("5. Testing espaloma charge calculation...")
        charges = charge(mol)
        print(f"   \u2713 Charges calculated: {charges}")
        print(f"   Shape: {charges.shape}, Type: {type(charges)}")

        print("6. Testing with a larger molecule...")
        mol2 = Chem.MolFromSmiles("CCO")  # Ethanol
        charges2 = charge(mol2)
        print(f"   \u2713 Ethanol charges: {charges2}")

        return True

    except ImportError as e:
        print(f"   \u2717 Import error: {e}")
        return False
    except Exception as e:
        print(f"   \u2717 Runtime error: {e}")
        return False

    finally:
        print("=" * 60)
        print()
