# PKS_DOCK

Fully automated, open-source *in silico* pipeline for PKS-I biosynthetic gene
cluster mining, compound retrieval, pocket-aware molecular docking, ADMET
prediction, and publication-ready reporting.

Built for antimicrobial-producing bacterial isolates whose genomes are not
yet sequenced (reference-genome-based workflow), screened against a
configurable panel of pathogen drug targets.

## What changed from the earlier draft (honesty log)

This version was rebuilt specifically to fix bugs and gaps identified in a
prior review. If you are comparing against an older version of this
pipeline, here is exactly what changed and why:

| Issue in the old draft | Fix in this version |
|---|---|
| Vina score parsing grabbed the wrong log line (`grep -A1 "mode |"` returned the literal string `"(kcal/mol)"`, not a number) | Rewritten with `awk` logic that finds the `-----` separator and reads the first real data row. **Verified against a synthetic Vina log before shipping** (see `tests/test_vina_parsing.sh`). |
| Grid box centered on a co-crystallized ligand that had already been stripped out (`--box_center_on_ligand` after `grep -v "^HETATM"`) | Replaced entirely. Grid boxes are now derived from **fpocket** pocket detection (Phase 8/8b) — the standard, defensible approach, and the one you explicitly asked for. |
| fpocket was completely absent | Added as Phase 8/8b: runs fpocket on every prepared receptor, ranks pockets by druggability score, and derives grid center/size from the top pocket's real coordinates. |
| ChEBI fallback wrote a raw API search response as if it were a structure file | Replaced with ChEBI's actual structure-download endpoint, NCI CIR as a second real fallback, and PubChem PUG-REST (direct HTTP, not an unverified library function) as the primary source. |
| `pcp.get_sdf()` — an unverified PubChemPy function | Replaced with direct PubChem PUG-REST HTTP calls, which are documented and stable. |
| `nhea_sequence.fasta` was required but never generated (undisclosed manual input) | Phase 6b now **auto-fetches the sequence from UniProt**; falls back to a manual file only if that genuinely fails, and says so out loud. |
| Hardcoded receptor panel presented as "works on any genome" | Still config-driven (`config/pathogen_targets.yaml`) — see below for why that is the *honest* choice, not a shortcut. Organism→genome resolution IS now automatic; organism→target-panel is not, and cannot honestly be, without a curated source. |
| PLIP parser only captured 3 of 8 interaction types | Now captures all 8: hydrogen bonds, hydrophobic contacts, salt bridges, water bridges, pi-stacking, pi-cation interactions, halogen bonds, metal complexes. |
| No `environment.yml` | Added, with a note on how to freeze it properly from YOUR working install before publishing. |
| 15 figures/tables, no repo structure documented | 20 outputs (15 figures + 5 tables), full directory tree below, and a manuscript/supplementary export script. |
| "Interpretation" was just terminal print statements | Phase 14 generates an actual prose Markdown interpretation report from your real result numbers. |
| Every phase rolled its own `requests.get(...)` with inconsistent (or no) retry logic | Replaced with one shared, tested `pks_dock.net.HTTPClient` (retry + exponential backoff + disk cache + `fallback_chain()`), now used by Phases 1, 4, 6, and 6b. |
| Phase 6b always ran a fresh local ColabFold prediction, even for receptors AlphaFold DB already had solved | Phase 6b now checks the real AlphaFold DB API by UniProt accession first and only falls back to local ColabFold if no prediction exists. |
| No reproducibility artifacts — decisions (which PDB, which pocket, which fallback DB) only ever existed as terminal output | `results/reports/decision_log.json` (machine-readable, crash-safe), `workflow_report.md` (human-readable), and `pipeline_metadata.json` (software/environment provenance) are now written automatically. |
| Not a real installable package — had to `cd` into the repo and remember `./run_pipeline.sh` | `pyproject.toml` + `src/pks_dock/` package; `pip install -e .` gives you a real `pks-dock` CLI. |
| No LICENSE, CITATION.cff, CONTRIBUTING.md, CHANGELOG.md, or CI | Added, plus `tests/test_net.py`, `test_alphafold.py`, `test_reproducibility.py` (20 offline unit tests, run in CI on every push). |

## Honest limitations (please read before you present this to your supervisor)

1. **Pathogen target panels are config-driven, not auto-discovered.** There
   is no reliable API that returns "the validated druggable targets for
   organism X" — every published docking study picks targets based on
   literature justification, not automated mining. `config/pathogen_targets.yaml`
   is intentionally the one human-curated part of this pipeline. To add a
   new pathogen, add a block and cite your source for each target.
2. **Genome accession resolution from an organism name is automatic**
   (Phase 1, via NCBI Assembly search) and picks the best available
   reference/representative genome — but you should still sanity-check that
   the assembly it picked is the one you'd cite in your thesis.
3. **fpocket's default pocket ranking is used as-is.** For receptors with
   multiple plausible pockets, always visually inspect the top 2–3 pockets
   (`results/fpocket/<receptor>/pockets/`) rather than trusting the ranking
   blindly for a thesis defense.
4. **External bioinformatics tools (antiSMASH, ColabFold, AutoDock Vina,
   fpocket, PLIP, admet-ai) are not bundled or testable in a sandboxed code
   environment** — this repository was built and its Python parsing logic
   was unit-tested against synthetic sample data (Vina logs, fpocket output),
   but the full end-to-end run has not been executed against real
   bioinformatics tool output. **Run it yourself on your WSL Ubuntu machine
   and fix anything that doesn't match your exact installed tool versions**
   before trusting the results for a thesis chapter.

## Directory structure

```
PKS_DOCK/
├── README.md
├── pyproject.toml                    <- installable `pks-dock` package (pip install -e .)
├── environment.yml
├── run_pipeline.sh                   <- single entry point (also wrapped by `pks-dock` CLI)
├── LICENSE  CITATION.cff  CONTRIBUTING.md  CHANGELOG.md
├── .github/workflows/ci.yml          <- lint + offline unit tests on every push/PR
├── src/pks_dock/                     <- installable package, shared across all phases
│   ├── net.py                        <- retry/backoff/caching HTTPClient + fallback_chain()
│   ├── alphafold.py                  <- real AlphaFold DB API integration
│   ├── reproducibility.py            <- DecisionLog, workflow_report.md, pipeline_metadata.json
│   └── cli.py                        <- `pks-dock` console-script entry point
├── config/
│   ├── pathogen_targets.yaml        <- target panels (human-curated, cited)
│   └── sequences/                   <- auto-fetched or manually supplied FASTA files
├── scripts/
│   ├── 01_fetch_genome.py
│   ├── 02_run_antismash.sh
│   ├── 03_parse_antismash.py
│   ├── 04_get_ligands.py
│   ├── 05_prep_ligands.sh
│   ├── 06_get_receptors.py
│   ├── 06b_predict_receptors.py
│   ├── 07_prep_receptors.sh
│   ├── 08_run_fpocket.sh
│   ├── 08b_generate_grid_configs.py
│   ├── 09_run_docking.sh
│   ├── 10_run_plip.sh
│   ├── 11_parse_plip.py
│   ├── 12_run_admet.py
│   ├── 13_generate_figures.py
│   ├── 14_generate_interpretation.py
│   ├── 15_generate_manuscript_package.py
│   └── 16_generate_reproducibility_report.py
├── results/                          <- all generated at runtime
│   ├── genomes/
│   ├── antismash/
│   ├── ligands/  ligands_pdbqt/
│   ├── receptors_raw/  receptors_pdbqt/
│   ├── fpocket/
│   ├── docking/
│   ├── plip/
│   ├── admet/
│   ├── figures/                      <- F1-F15, 300 DPI
│   ├── tables/                       <- Table01-Table05
│   └── reports/
│       ├── interpretation.md
│       ├── decision_log.json
│       ├── workflow_report.md
│       ├── pipeline_metadata.json
│       └── manuscript_package/
├── logs/                             <- timestamped full pipeline logs
├── cache/http/                       <- on-disk cache for retried API calls (HTTPClient)
├── tests/
│   ├── test_net.py                   <- retry/backoff/cache/fallback_chain unit tests
│   ├── test_alphafold.py             <- AlphaFold DB hit/miss/error unit tests
│   ├── test_reproducibility.py       <- decision log + report rendering unit tests
│   └── test_vina_parsing.sh          <- verifies the Vina log parsing fix
└── docs/
    └── DIRECTORY_STRUCTURE.md
```

## Setup

```bash
conda env create -f environment.yml
conda activate pksdock
```

This also installs the `pks-dock` package itself (`pip install -e .`,
declared in `environment.yml`), giving you the `pks-dock` command below.
If you're not using conda, `pip install -e ".[dev]"` from the repo root
installs the package and dev/test tooling on its own.

Before relying on this for publication, regenerate a real lock file from
your own tested install (see the note at the bottom of `environment.yml`).

## Usage

Two equivalent ways to run the pipeline:

```bash
# 1) Installed CLI (recommended -- works from anywhere once `pip install -e .`
#    has been run; also finalizes results/reports/workflow_report.md and
#    pipeline_metadata.json even if a phase fails partway through)
pks-dock --organism "Bacillus velezensis" \
    --pathogens Staphylococcus_aureus Bacillus_cereus Pseudomonas_aeruginosa \
    --threads 8

# 2) Direct shell script (must be run from inside the repo)
./run_pipeline.sh \
    --accession GCF_000063585.1 \
    --pathogens Staphylococcus_aureus \
    --threads 8
```

Every phase prints, in real time, what it's doing, why, and what it found.
A full combined log is saved to `logs/pipeline_log_<timestamp>.txt`. Every
automated decision (which reference genome, which PDB structure, which
fallback database resolved a compound, whether AlphaFold DB already had a
receptor prediction) is additionally written to:

- `results/reports/decision_log.json` -- machine-readable, appended to
  incrementally as each phase runs (survives a crash partway through).
- `results/reports/workflow_report.md` -- the same decisions rendered as a
  human-readable, phase-by-phase narrative.
- `results/reports/pipeline_metadata.json` -- software version, git commit,
  Python/platform info, and the exact CLI parameters used for the run.

## AlphaFold integration

Phase 6b (structure prediction for receptors with no experimental PDB) now
checks the [AlphaFold DB API](https://alphafold.ebi.ac.uk/api-docs) for an
existing prediction for the resolved UniProt accession *before* running a
new local ColabFold job. If AlphaFold DB has one, it's downloaded directly
(with retries/caching via `pks_dock.net.HTTPClient`) and ColabFold is
skipped entirely; if not, the pipeline says so explicitly and falls back to
local prediction, exactly as before. Either way the choice and the reason
are written to `decision_log.json`.

## Adding a new pathogen or target

Edit `config/pathogen_targets.yaml`. Add a block with `candidate_pdb_ids`
(the pipeline will auto-select the best one by resolution) or set
`predict_if_empty: true` if no experimental structure exists (the pipeline
will auto-fetch the sequence from UniProt and predict with ColabFold).
Always cite the source that validated the target as druggable.

## Development

```bash
pip install -e ".[dev]"
pytest tests/ -v          # offline unit tests: retries, caching, fallback
                            # chains, AlphaFold DB integration, decision log
ruff check src/ tests/
black --check src/ tests/
```

See `CONTRIBUTING.md` for the ground rules (no fabricated science, no new
hardcoded IDs, every external call goes through `pks_dock.net.HTTPClient`,
every automated choice gets logged). `CHANGELOG.md` tracks what's changed
release to release.

## Citation

If you use this pipeline, please cite the underlying tools it wraps:
antiSMASH, AutoDock Vina, fpocket, PLIP, admet-ai, ColabFold/AlphaFold,
Open Babel, and Biopython — and consider citing this repository itself
once you have assigned it a DOI (e.g. via Zenodo).
