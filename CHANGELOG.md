# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.2.0] - 2026-07-10

### Added
- `src/pks_dock/` installable package (`pip install -e .`) with a real
  `pks-dock` console-script entry point wrapping `run_pipeline.sh`.
- `pks_dock.net.HTTPClient`: shared retry/backoff/disk-cache HTTP layer,
  plus a `fallback_chain()` helper, now used by Phase 1 (NCBI genome
  download), Phase 4 (PubChem -> NCI/CIR -> ChEBI compound resolution),
  Phase 6 (RCSB PDB candidate selection/download), and Phase 6b
  (UniProt sequence + AlphaFold DB).
- `pks_dock.alphafold`: real AlphaFold DB API integration
  (`alphafold.ebi.ac.uk/api/prediction/{accession}`). Phase 6b now checks
  AlphaFold DB for an existing prediction before running a local ColabFold
  job, instead of always predicting from scratch.
- `pks_dock.reproducibility`: `DecisionLog` (writes
  `results/reports/decision_log.json` incrementally, crash-safe),
  `render_workflow_report()` (renders `workflow_report.md`), and
  `write_pipeline_metadata()` (`pipeline_metadata.json` with software
  version, git commit, environment, and run parameters). Wired into
  Phases 1, 4, 6, and 6b so far.
- Repo scaffolding: `pyproject.toml`, `LICENSE` (MIT), `CITATION.cff`,
  `CONTRIBUTING.md`, this `CHANGELOG.md`, `.github/workflows/ci.yml`.
- Real unit tests (`tests/test_net.py`, `tests/test_alphafold.py`,
  `tests/test_reproducibility.py`) covering retry/backoff, caching,
  fallback-chain exhaustion, AlphaFold DB hit/miss/error handling, and
  decision-log/report rendering. All run offline (mocked HTTP).

### Changed
- Phase 1, 4, 6, and 6b scripts now import the shared `pks_dock` package
  instead of each rolling its own `requests.get(...)` calls.

## [0.1.0] - earlier draft

- Initial fully-scripted (bash + standalone Python scripts) pipeline:
  genome acquisition, antiSMASH mining, ligand/receptor retrieval,
  fpocket-derived docking grids, AutoDock Vina, PLIP, ADMET, figures, and
  manuscript packaging. See `README.md`'s "What changed from the earlier
  draft" table for the bug fixes made at this stage.
