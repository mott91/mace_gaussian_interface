from ase import Atoms
from ase.build import molecule
from ase.io import write


def build_n_alkane(n_carbons: int):
    """Build a linear all-trans n-alkane with n carbons."""
    if n_carbons == 1:
        return molecule("CH4")

    # typical bond lengths
    cc = 1.54  #  A
    ch = 1.09  #  A

    atoms = []
    positions = []

    # build carbon backbone along x-axis
    for i in range(n_carbons):
        atoms.append("C")
        positions.append([i * cc, 0.0, 0.0])

    # add hydrogens
    for i in range(n_carbons):
        x, y, z = positions[i]
        if i == 0:  # CH3 group
            positions += [
                [x - 0.63,  ch,  0.63],
                [x - 0.63, -ch,  0.63],
                [x - 0.63,  0.0, -ch],
            ]
            atoms += ["H", "H", "H"]
        elif i == n_carbons - 1:  # CH3 group
            positions += [
                [x + 0.63,  ch,  0.63],
                [x + 0.63, -ch,  0.63],
                [x + 0.63,  0.0, -ch],
            ]
            atoms += ["H", "H", "H"]
        else:  # CH2 group
            # alternate up/down orientation for zig-zag shape
            sign = 1 if i % 2 == 0 else -1
            positions += [
                [x,  sign * ch,  ch * 0.4],
                [x, -sign * ch, -ch * 0.4],
            ]
            atoms += ["H", "H"]

    return Atoms(symbols=atoms, positions=positions)

# Generate methane \u2192 octane
for n in range(1, 9):
    mol = build_n_alkane(n)
    fname = f"n_alkane_{n}C.xyz"
    write(fname, mol)
    print(f"Saved: {fname}")
