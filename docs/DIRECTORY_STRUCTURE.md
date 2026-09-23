# Directory Structure

This document is the single source of truth for where everything lives.
If you extend the pipeline, keep new outputs inside `results/<new_subdir>/`
and new scripts inside `scripts/`, numbered to match their phase order.

```
PKS_DOCK/
├── README.md                        Project overview, honesty log, usage
├── environment.yml                  Pinned conda/pip environment
├── run_pipeline.sh                  Master entry point (parameter-driven)
│
├── config/
│   ├── pathogen_targets.yaml        Target panels per pathogen (human-curated, cited)
│   └── sequences/                   FASTA files for ColabFold prediction targets
│                                    (auto-fetched from UniProt when possible;
│                                     manually supplied only when auto-fetch fails)
│
├── scripts/                         One script per pipeline phase, numbered in run order
│   ├── 01_fetch_genome.py           Organism name or accession -> reference genome
│   ├── 02_run_antismash.sh          antiSMASH cluster mining
│   ├── 03_parse_antismash.py        Extract PKS-I clusters + predicted compounds
│   ├── 04_get_ligands.py            Auto-fetch compound structures (PubChem/CIR/ChEBI)
│   ├── 05_prep_ligands.sh           SDF -> PDBQT ligand preparation
│   ├── 06_get_receptors.py          Auto-select best experimental PDB per target
│   ├── 06b_predict_receptors.py     ColabFold prediction for targets with no PDB
│   ├── 07_prep_receptors.sh         Receptor cleaning, repair, hydrogens, PDBQT
│   ├── 08_run_fpocket.sh            Binding pocket detection
│   ├── 08b_generate_grid_configs.py Grid box derivation from detected pockets
│   ├── 09_run_docking.sh            Batch AutoDock Vina docking
│   ├── 10_run_plip.sh               Non-covalent interaction profiling
│   ├── 11_parse_plip.py             Parse all 8 PLIP interaction types
│   ├── 12_run_admet.py              ADMET property prediction
│   ├── 13_generate_figures.py       15 figures + 4 tables (Table05 comes from 08b)
│   ├── 14_generate_interpretation.py Auto-written prose interpretation report
│   └── 15_generate_manuscript_package.py  Supplementary Excel workbook + methods stub
│
├── results/                          Generated at runtime, one subdir per phase
│   ├── genomes/                      Downloaded reference genome + manifest
│   ├── antismash/                    Raw antiSMASH output + parsed compound list
│   ├── ligands/                      Resolved SDF structures + resolution manifest
│   ├── ligands_pdbqt/                Docking-ready ligand files
│   ├── receptors_raw/                Downloaded/predicted receptor PDBs + manifest
│   ├── receptors_pdbqt/              Cleaned, docking-ready receptor files
│   ├── fpocket/                      Pocket detection output + grid box configs
│   ├── docking/                      Vina output poses, logs, summary.csv
│   ├── plip/                         Interaction XML reports + interaction_summary.csv
│   ├── admet/                        ADMET predictions
│   ├── figures/                      F1-F15, 300 DPI PNGs
│   ├── tables/                       Table01-Table05 CSVs
│   └── reports/
│       ├── interpretation.md         Auto-generated prose interpretation
│       └── manuscript_package/       Supplementary workbook + figures + methods stub
│
├── logs/                             One timestamped full-pipeline log per run
├── tests/                            Standalone tests that don't require the
│                                     real bioinformatics tools to be installed
│   └── test_vina_parsing.sh          Proves the Vina log-parsing fix
└── docs/
    └── DIRECTORY_STRUCTURE.md        This file
```
