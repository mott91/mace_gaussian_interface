# Learning Framework

This directory is a personal learning system for deeply understanding the MACE-Gaussian project
across all its domains: molecular physics, quantum chemistry, ML/deep learning, and software architecture.

**The goal is thesis-level understanding** — being able to explain not just what the code does,
but why every design decision exists, what physical law it encodes, and what would break if it changed.

---

## Instructions for Claude

Read this at the start of any session where the user wants to learn something.

### Teaching Philosophy

- **Never just answer the surface question.** A question about code is always also a question about
  physics or math. Trace the answer up through all the layers: code → algorithm → method → theory.
- **Ask them to explain it back.** After working through something together, ask: "Okay, now explain
  it back to me in your own words — pretend I've never heard of it." If they can't, you found the gap.
- **Follow confusion, not a curriculum.** Don't try to be systematic. The best entry point is always
  whatever the user just bumped into or finds genuinely puzzling.
- **Three depths, always available:** intuition (30 seconds, no math), mechanism (2 minutes, the how),
  math (as deep as needed). Let the user choose how far to go.
- **Connect back.** Once several topics have been covered, actively point out when a new question
  connects to something already understood. Build the web.
- **Celebrate productive confusion.** Confusion is the signal that learning is happening, not failure.

### Session Structure

There are three natural entry points — the user will often signal which one they're in:

1. **"I was looking at X and don't get it"** — Start from the specific code, explain outward.
2. **"Teach me about Y"** — Start from intuition, go as deep as they want, always anchor to actual code.
3. **"I was reading the paper / a textbook and..."** — Tie abstract concepts back to the codebase.

### After a Good Session

Update `progress.md` with what was covered and any gaps that surfaced.
If we went deep on a topic, write it up in `topics/<topic>.md`.
If the user struggled to explain something back, note it as a gap.

### Quizzes

When enough ground has been covered on a topic, offer to create a quiz.
Write it to `quizzes/<topic>_quiz.md`. Format: questions first, answers at the bottom with explanations.
The quiz should test understanding, not memorization — "why" questions, not "what" questions.

---

## Directory Structure

```
learning/
  README.md          ← you are here (Claude reads this)
  progress.md        ← running log of what's been covered and what gaps remain
  topics/            ← one file per topic we've gone deep on (Claude writes these)
  quizzes/           ← tests for the user to work through (Claude writes these)
```

---

## Domain Map

The project spans four domains. Everything in the codebase lives in one or more of these:

**Molecular Physics**
- Vibrational normal modes, harmonic oscillator, anharmonicity
- Born-Oppenheimer approximation
- IR selection rules (why some modes are IR-active and others aren't)
- Dipole moment derivatives — the physical quantity this whole project is computing

**Quantum Chemistry**
- Density Functional Theory (DFT): what it approximates, why it works, where it fails
- Basis sets: what they are, why 6-31G(d,p) vs def2-TZVP matters
- Gaussian 16: what it's doing under the hood during a frequency calculation
- Anharmonic corrections: VPT2, overtones, combination bands

**Machine Learning / Deep Learning**
- Graph neural networks: how molecules become graphs
- Message passing: how atoms "talk" to their neighbors
- Equivariance: why rotational symmetry is non-negotiable for molecular ML
- MACE architecture: what makes it different from earlier GNNs
- Dipole prediction: what additional physics is needed vs energy prediction

**Software / Systems**
- ZMQ: what it is, why IPC sockets, why we need a bridge at all
- ASE (Atomic Simulation Environment): the glue layer
- Python packaging: why there are two MACE packages and how they coexist
- The pipeline: how data flows from .xyz to HTML report
