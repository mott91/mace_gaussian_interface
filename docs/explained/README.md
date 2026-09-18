# docs/explained

Written 2026-09-18 from a full read of `mace_gaussian/` on branch `spike/24-vpt2-psience`
(commit f929741). These are for understanding and for the defense, not for the thesis text;
the thesis chapters should be written from the code and from these, in your own words.

| file | what it answers | use it for |
|---|---|---|
| `01_physics.md` | What is being computed and why it is physically meaningful | explaining to the supervisor; Chapter 2 and 4 scaffolding; every section ends with a one-sentence version to say aloud |
| `02_software.md` | How the code makes Gaussian believe MACE is a quantum chemistry program; `mace-gaussian run water.xyz` traced with real files and numbers | Chapter 3; whiteboard explanations; "show me where that happens" questions |
| `03_hard_questions.md` | The questions a committee will ask, with spoken answers and the reasoning behind them | defense prep; the discussion chapter |
| `04_defense_outline.md` | 15 content slides for 25 minutes plus 14 backups, with the message and speaker notes per slide | building the slide deck once the campaign data exists |
| `REVIEW_FINDINGS.md` | Ranked code review: 4 high, 6 medium, 8 low, plus what was checked and found correct | to walk through together; nothing is fixed yet |

Reading order for a first pass: 01, then 02 §1-§6, then REVIEW_FINDINGS H1-H4, then 03.

Facts in these files were checked against the actual outputs under `comparison_results/water/`
and `comparison_results/methane/` from the July 2026 runs. Where a number in the thesis draft
or `docs/methods.md` disagrees with these files, the code is right and the older text is out of
date (see `REVIEW_FINDINGS.md` L8 for the list).
