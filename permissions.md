# Claude Permissions

Rules for what Claude may do autonomously vs. what requires user confirmation.

---

## Auto-approve

Low-risk, local, reversible actions — proceed without asking.

**File operations**
- Read, search (grep/glob), or write any file
- Edit source code, configs, docs, notebooks

**Code quality**
- Run linting, formatting, type checking (ruff, mypy, ty, eslint, etc.)
- Run tests (pytest, jest, etc.)
- Run diagnostic/status commands (diagnose, list, status, info)

**Git**
- Stage files (git add)
- Create commits (with appropriate message)
- Create branches

**Computation — known inputs**
- Re-run calculations or scripts on previously processed inputs/molecules

---

## Require confirmation

Irreversible, resource-intensive, or externally visible — always ask first.

**Git / remote**
- Push to remote (git push)
- Create, update, or merge pull requests
- Delete branches

**Destructive operations**
- Delete files or directories
- Overwrite uncommitted changes (git reset --hard, git checkout --)
- Any rm -rf

**Dependencies / environment**
- Modify pyproject.toml, package.json, uv.lock, requirements.txt
- Install or remove packages (pip, uv, npm, micromamba, conda)

**Long-running computation — new inputs**
- Run full calculation pipelines on new molecules, datasets, or inputs not previously processed

**External / shared systems**
- Post comments, messages, or notifications
- Modify CI/CD pipelines or shared infrastructure
