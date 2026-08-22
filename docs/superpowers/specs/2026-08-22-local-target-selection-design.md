# Local Target Selection Design

## Goal

Allow the user to choose between using beta-lactamase structures downloaded from RCSB or importing a desired local PDB file from their computer.

## Approved Behavior

The local PDB option copies the chosen file into `data/structures/` and makes it the only selected target for downstream preparation, docking, and genetic-algorithm validation workflows that rely on the selected target.

## Architecture

The existing pipeline already separates downloading (`scripts/fetch_structures.py`) from target selection (`scripts/select_targets.py`). The change belongs in `scripts/select_targets.py` because it is the point where `config.yaml:selected_targets` is written.

`select_targets.py` will support two target sources:

- Downloaded targets from `data/structures_index.json`, preserving the current table-based selection.
- A local `.pdb` file supplied interactively or through `--local-pdb`.

For local files, the script will:

1. Validate that the source file exists and has a `.pdb` extension.
2. Validate that Biopython can parse it and that it contains atoms.
3. Derive a safe target identifier from the filename.
4. Copy the file to `data/structures/<target_id>.pdb`.
5. Upsert a metadata entry into `data/structures_index.json`.
6. Replace `config.yaml:selected_targets` with a one-item list containing only `<target_id>`.
7. Remove stale `binding_site.targets` values so `prepare_protein.py` recalculates the active site for the imported target.

## Error Handling

Invalid paths, non-PDB extensions, empty files, parse failures, and conflicting destination names should produce clear terminal errors. Re-importing a local file with the same target id can overwrite the copied structure only when the user confirms it in interactive mode or passes `--force` in CLI mode.

## Documentation

The README target-selection step should show both supported paths:

- Download structures, then select from the downloaded list.
- Import a local PDB and use only that target.

## Testing

Manual validation is sufficient for this small terminal workflow:

- Run `python -m py_compile scripts/select_targets.py`.
- Create or reuse a valid local PDB, run `python scripts/select_targets.py --local-pdb <path> --force`, and verify:
  - copied PDB exists in `data/structures/`;
  - `data/structures_index.json` contains the imported target;
  - `config.yaml:selected_targets` contains only the imported target.
