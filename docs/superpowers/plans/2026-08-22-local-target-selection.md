# Local Target Selection Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an option to import a local beta-lactamase PDB file and make it the only selected target for downstream docking and genetic-algorithm workflows.

**Architecture:** Extend `scripts/select_targets.py`, because it owns the transition from available structures to `config.yaml:selected_targets`. Keep the current downloaded-target selection intact, and add a local-target import path that validates, copies, indexes, and selects one target.

**Tech Stack:** Python 3, argparse, pathlib, shutil, yaml, Bio.PDB.

---

## Chunk 1: Local PDB Import

### Task 1: Extend Target Selection Script

**Files:**
- Modify: `scripts/select_targets.py`

- [x] **Step 1: Add imports and constants**

Add `argparse`, `re`, `shutil`, `Path`, `PDBParser`, and `PDBConstructionWarning`. Define `STRUCTURES_DIR = Path("data") / "structures"` while preserving `INDEX_FILE` and `CONFIG_FILE`.

- [x] **Step 2: Add local PDB helpers**

Add focused helpers:

```python
def sanitize_target_id(path: Path) -> str:
    target_id = re.sub(r"[^A-Za-z0-9_-]+", "_", path.stem).strip("_")
    return target_id.upper() or "LOCAL_TARGET"
```

Also add helpers to validate a PDB with Biopython, copy it into `data/structures/`, load/save the index with missing-index tolerance, and upsert a local metadata entry.

- [x] **Step 3: Update config writing**

Change `update_config(selected)` so it can optionally clear stale `binding_site.targets` values when a local target is imported. For local import, write only the imported target to `selected_targets`.

- [x] **Step 4: Add CLI mode**

Add `--local-pdb PATH` and `--force`. When `--local-pdb` is present, import the file, update the index, update config with only that target, print the next command, and exit without showing the downloaded-target table.

- [x] **Step 5: Add interactive source choice**

At startup, ask whether the user wants to select from downloaded targets or import a local PDB. Keep downloaded selection as the default so existing usage remains familiar.

- [x] **Step 6: Validate syntax**

Run: `python -m py_compile scripts/select_targets.py`

Expected: no output and exit code 0.

## Chunk 2: Documentation

### Task 2: Document Both Selection Paths

**Files:**
- Modify: `README.md`

- [x] **Step 1: Update Etapa 2**

Document:

```bash
python scripts/select_targets.py
python scripts/select_targets.py --local-pdb /caminho/para/minha_beta_lactamase.pdb --force
```

Explain that local import copies the file to `data/structures/` and replaces the current target selection with that one target.

- [x] **Step 2: Review diff**

Run: `git diff -- scripts/select_targets.py README.md docs/superpowers/specs/2026-08-22-local-target-selection-design.md docs/superpowers/plans/2026-08-22-local-target-selection.md`

Expected: diff contains only local-target selection changes.

- [ ] **Step 3: Commit**

Because the repository already has unrelated modified files, stage only the files changed for this feature:

```bash
git add scripts/select_targets.py README.md docs/superpowers/specs/2026-08-22-local-target-selection-design.md docs/superpowers/plans/2026-08-22-local-target-selection.md
git commit -m "feat: support local beta-lactamase target import"
```
