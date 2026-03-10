---
created: 2026-03-10T15:17:30.446Z
title: Make dipole calculator factory lazy instantiation
area: general
files:
  - mace_gaussian/calculators/factory.py
---

## Problem

`dipole_factory = DipoleCalculatorFactory()` is a module-level global in `factory.py`. When this module is imported, `_register_calculators()` immediately instantiates all three calculators — `EspalomaDipoleCalculator()`, `XTBDipoleCalculator()`, `MACEMLDipoleCalculator()`. Each calls `_check_availability()`, which imports rdkit + espaloma_charge (runs a test molecule), imports xtb, and stat-checks the MACE model file.

This happens on every `mace-gaussian` invocation, for every calculator, even if only one is used. It slows startup and couples imports unnecessarily.

## Solution

Make the factory lazy: only instantiate a calculator when it's first requested by name.

```python
def get_calculator(self, method: str) -> DipoleCalculatorBase:
    if method not in self._instances:
        self._instances[method] = self._constructors[method]()
    return self._instances[method]
```

Store constructors (not instances) at init time. Instantiate on first `get_calculator()` call. The `auto` path instantiates in preferred order, stopping at the first available one.

Keep the `list_available()` method working (it can lazily instantiate all to check, or just report based on what's been instantiated).
