# Contributing to PKS_DOCK

Thanks for considering a contribution. This is a research pipeline first and
a piece of software second: correctness and honest reporting of automated
decisions matter more than feature count.

## Ground rules

1. **No fabricated science.** If a phase cannot automatically determine
   something (a receptor, a compound structure, a target panel), it must say
   so explicitly and stop or fall back — never guess and present the guess
   as fact. See `config/pathogen_targets.yaml` for the one intentionally
   human-curated part of the pipeline, and the README's "Honest limitations"
   section for the reasoning.
2. **No new hardcoded accessions/IDs/paths.** Anything that varies by
   organism, target, or compound belongs in `config/` or as a CLI argument.
3. **Every external API call goes through `pks_dock.net.HTTPClient`.** This
   gives you retries, backoff, and disk caching for free and keeps behavior
   consistent across phases. Don't add a bare `requests.get(...)`.
4. **Every automated choice gets logged.** Use
   `pks_dock.reproducibility.DecisionLog.record(...)` at the point a phase
   picks between alternatives (a PDB structure, a pocket, a fallback
   database), so it shows up in `workflow_report.md`.

## Development setup

```bash
git clone https://github.com/kizitlabs/pks-dock.git
cd pks-dock
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install
```

## Running the test suite

```bash
pytest tests/ -v
ruff check src/ tests/
black --check src/ tests/
mypy src/pks_dock
```

`scripts/` (the phase-by-phase pipeline entry points) is intentionally not
yet in the linted/formatted set — it is being migrated into `src/pks_dock/`
incrementally, phase by phase, alongside real test coverage for each one.

Tests must not require network access or any of the external bioinformatics
tools (antiSMASH, ColabFold, AutoDock Vina, fpocket, PLIP) to be installed —
mock the HTTP layer (`pks_dock.net.HTTPClient`) and subprocess calls, as the
existing tests in `tests/` do.

## Adding a new pathogen / target panel

Add a block to `config/pathogen_targets.yaml` with a literature citation for
each target (see the file's header comment for why this stays human-curated).
Do not add targets without a citation.

## Pull requests

- Keep PRs scoped to one phase or one concern where possible.
- Update `CHANGELOG.md` under "Unreleased".
- If you touch a phase script's CLI arguments, update the corresponding line
  in `run_pipeline.sh` and `README.md`.
- New external database integrations should follow the pattern in
  `src/pks_dock/alphafold.py`: a small, focused module with a
  `resolve_*`-style function that returns `(result_or_None, raw_hit_or_None,
  human_readable_note)` so callers can log the decision honestly.
