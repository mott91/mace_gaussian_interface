# Taylor Expansion, Dipole Moment, and Why Overtones Exist

_Saved mid-session 2026-02-22. Resume here next time._

---

## Where we are in the story

We are trying to understand **why overtones (v=0→v=2 transitions) are forbidden in the
harmonic model but allowed in real molecules.** The answer requires understanding Taylor
expansions and what "linear dipole" actually means. We built all the foundations and were
about to reconnect them.

---

## What physically happens in v=0→v=1

When a molecule absorbs a photon with exactly the right energy (matching the gap between
rungs v=0 and v=1):
- The vibration **amplitude increases** — atoms swing further from equilibrium
- The **frequency stays the same** (harmonic approximation)
- The photon is absorbed, the molecule is now in state v=1

The photon doesn't start the vibration. It adds energy to one that's already happening (ZPE).

---

## Smooth functions

A function with no sharp corners, no jumps, no kinks. You can zoom in anywhere and it
always looks like a nice arc. sin(x), x², e^x are smooth. |x| is NOT — it has a kink at
x=0. Smooth functions are the ones where calculus works cleanly everywhere.

---

## Taylor expansion

**Core idea:** If you know everything about a function at one point — its height, its
slope, its curvature, how the curvature is changing — you can reconstruct the whole
function just from that local information.

**The formula** (near x=0):
```
f(x) = f(0) + f'(0)·x + f''(0)·x²/2 + f'''(0)·x³/6 + ...
         ↑       ↑            ↑               ↑
       height  slope       curvature     how curvature changes
```

**Example:** e^x near x=0
- Just constant:  f ≈ 1
- + linear term:  f ≈ 1 + x
- + quadratic:    f ≈ 1 + x + x²/2
- + cubic:        f ≈ 1 + x + x²/2 + x³/6   (getting very accurate)

**Key insight for physics:** For small x, higher-power terms are negligibly small.
x=0.1 → x²=0.01 → x³=0.001. So for small displacements you can stop after a few terms
and get a good approximation.

---

## Cubic wiggle

The x³ term is what creates **asymmetry**. The function x³ is positive on the right,
negative on the left. Adding a small x³ to a symmetric parabola (x²) tilts it: one wall
of the potential well gets steeper, the other shallower.

This is physically real for molecular bonds: you can compress a bond only so far (atoms
repel hard), but you can stretch it quite far before it breaks. The real potential well IS
asymmetric. The Morse potential is the classic example.

---

## Harmonic approximation from Taylor expansion

The potential energy of a vibrating bond is some complicated function of displacement x.
Taylor expand it near equilibrium (x=0):

```
V(x) = V(0) + V'(0)·x + V''(0)·x²/2 + V'''(0)·x³/6 + ...
```

**Why V'(0) = 0 (flat bottom):**
The equilibrium position is the MINIMUM of the potential energy — the atom naturally
sits there. At any minimum of any smooth function, the slope is zero. If the slope
weren't zero, it would be a hillside, and the atom would slide further down. So V'(0) = 0
and the linear term drops out completely.

**Result:**
```
V(x) ≈ V(0) + 0 + ½kx² + (higher terms)
            ↑
         set to 0 (just choose equilibrium as energy reference)
```

**Harmonic approximation = stop here:** V(x) ≈ ½kx²

This is NOT a physical law. It's "we kept only the lowest-order surviving term in the
Taylor expansion." Adding the cubic and quartic terms back gives you the anharmonic
potential.

---

## Dipole moment of a charge cloud

For point charges: dipole = charge × distance, as a vector.

For a fuzzy electron cloud:
- Find the "center of gravity" of all negative charge (weighted average position)
- Find the "center of gravity" of all positive charge (the nuclei)
- The dipole is the **vector** pointing from the negative center to the positive center,
  scaled by the total charge

It is a **single vector** — one arrow with magnitude and direction. Not a point in space,
not a field. It describes how separated the charges are and which way.

For water: oxygen is electronegative → net negative near O, net positive toward H's.
The dipole vector points from O toward the midpoint between the H's.

---

## THE KEY CONNECTION (where we stopped — start here next session)

Now we have all the pieces. Here's the argument to walk through:

**1. Taylor expand the dipole moment as a function of displacement:**
```
μ(x) = μ₀ + μ'·x + μ''·x²/2 + ...
```
- μ₀ = equilibrium dipole
- μ'·x = linear term (how dipole changes with small displacement)
- μ''·x²/2 = quadratic correction (small, needs anharmonicity to matter)

**2. Harmonic approximation throws away the x² term:**
So in the harmonic world, μ(x) ≈ μ₀ + μ'·x. Dipole is LINEAR in displacement.

**3. If the molecule vibrates as x(t) = A·sin(ωt), what does the dipole do?**
- Harmonic (linear dipole): μ(t) = μ₀ + μ'·A·sin(ωt) → oscillates at ω only
- Can only absorb photons at frequency ω → only fundamentals

**4. With the quadratic term:**
- μ(t) = μ₀ + μ'·A·sin(ωt) + μ''·A²·sin²(ωt)/2
- sin²(ωt) = ½(1 - cos(2ωt))  ← THIS IS THE KEY MATH STEP
- So μ(t) now contains a component oscillating at **2ω**
- That component can absorb a photon at 2ω → **overtone**

**5. Why VPT2:**
VPT2 treats the cubic and quartic terms in V(x) as small perturbations on the harmonic
solution. This lets Gaussian compute:
- How much the rung spacings shift (corrected fundamental frequencies)
- The intensities of overtones and combination bands
- The coupling between modes (which is why it needs 3rd and 4th derivatives)

---

## Questions to ask at start of next session

1. "Can you explain Taylor expansion back to me in your own words?"
2. "Why does V'(0) = 0 at equilibrium?"
3. "If we have the answer to both of those, walk me through why the quadratic dipole
   term is what allows overtone absorption."
