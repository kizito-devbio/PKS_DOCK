#!/usr/bin/env bash
# =============================================================================
# PKS_DOCK MASTER PIPELINE
# =============================================================================
#
# PURPOSE
# -------
# This pipeline discovers biosynthetic gene clusters (BGCs) from a
# PKS-producing organism and evaluates the predicted natural products against
# therapeutic protein targets from one or more user-selected pathogens.
#
# WORKFLOW
# --------
#
#   PKS-producing organism
#            |
#            v
#   Genome retrieval
#            |
#            v
#   antiSMASH
#            |
#            v
#   Compound identification
#            |
#            v
#   Ligand preparation
#            |
#            v
#   Retrieve pathogen receptors
#            |
#            v
#   Optional receptor prediction (ColabFold)
#            |
#            v
#   Receptor preparation
#            |
#            v
#   fpocket
#            |
#            v
#   AutoDock Vina
#            |
#            v
#   PLIP
#            |
#            v
#   ADMET
#            |
#            v
#   Figures, Tables and Manuscript Package
#
# NOTE
# ----
# The organism supplied via --organism or --accession is ONLY the source of
# biosynthetic gene clusters and candidate metabolites. Docking targets are
# obtained exclusively from the pathogen panel defined in
# config/pathogen_targets.yaml.
#
# USAGE:
#   ./run_pipeline.sh --organism "Bacillus velezensis" \
#                      --pathogens Staphylococcus_aureus Bacillus_cereus Pseudomonas_aeruginosa \
#                      --threads 8
#
#   ./run_pipeline.sh --accession GCF_000063585.1 \
#                      --pathogens Staphylococcus_aureus
#
# WHY --pathogens IS STILL A REQUIRED ARGUMENT, NOT SOMETHING SCRIPTED FOR YOU:
# There is no reliable API that returns "the validated druggable targets for
# organism X" -- target selection is a human, literature-justified decision in
# every published docking study. This is the one honestly-manual part of an
# otherwise fully automated pipeline (see config/pathogen_targets.yaml). If
# you want to add a new pathogen, add a config block and pass its key here.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/config/pathogen_targets.yaml"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Missing configuration file:"
    echo "$CONFIG_FILE"
    exit 1
fi

python3 -c "import yaml" 2>/dev/null || {
    echo "ERROR: Python package PyYAML is required (pip install pyyaml)"
    exit 1
}

ORGANISM=""
ACCESSION=""
PATHOGENS=()

# -------------------------------------------------------------------------
# Default runtime parameters.
# These may be overridden by command-line arguments.
# -------------------------------------------------------------------------
THREADS=4
EXHAUSTIVENESS=16
NUM_MODES=9
CPU=8

usage() {
    echo "Usage: $0 (--organism \"Genus species\" | --accession ACCESSION) --pathogens KEY [KEY...] [--threads N] [--exhaustiveness N] [--num-modes N] [--cpu N]"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --organism) ORGANISM="$2"; shift 2 ;;
        --accession) ACCESSION="$2"; shift 2 ;;
        --pathogens) shift; while [[ $# -gt 0 && "$1" != --* ]]; do PATHOGENS+=("$1"); shift; done ;;
        --threads) THREADS="$2"; shift 2 ;;
        --exhaustiveness) EXHAUSTIVENESS="$2"; shift 2 ;;
        --num-modes) NUM_MODES="$2"; shift 2 ;;
        --cpu) CPU="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

if [[ -z "$ORGANISM" && -z "$ACCESSION" ]]; then
    echo "ERROR: supply either --organism or --accession"; usage
fi
if [[ ${#PATHOGENS[@]} -eq 0 ]]; then
    echo "ERROR: supply at least one --pathogens KEY (see config/pathogen_targets.yaml for valid keys)"; usage
fi

mkdir -p logs
LOGFILE="logs/pipeline_log_$(date +%Y%m%d_%H%M%S).txt"

# -------------------------------------------------------------------------
# build_organism_lookup
#
# Builds the JSON mapping of {pathogen_key: scientific_organism_name} that
# Phase 6b passes to 06b_predict_receptors.py. This mapping must come from
# config/pathogen_targets.yaml -- NOT from the PKS-producing organism -- since
# the receptors being predicted belong to the pathogen panel, not to the
# organism supplied via --organism/--accession.
# -------------------------------------------------------------------------
build_organism_lookup() {
python3 - "$CONFIG_FILE" "${PATHOGENS[@]}" <<'PY'
import sys
import json
import yaml

config_file = sys.argv[1]
pathogens = sys.argv[2:]

with open(config_file) as f:
    config = yaml.safe_load(f)

lookup = {}

for pathogen in pathogens:
    if "pathogens" in config:
        entry = config["pathogens"].get(pathogen)
    else:
        entry = config.get(pathogen)

    if entry is None:
        raise SystemExit(f"ERROR: {pathogen} not found in pathogen_targets.yaml")

    organism = entry.get("organism")

    if organism is None:
        raise SystemExit(f"ERROR: {pathogen} has no organism field")

    lookup[pathogen] = organism

print(json.dumps(lookup))
PY
}

{
echo "==============================================================================="
echo " PKS_DOCK PIPELINE -- STARTED $(date)"
echo " PKS Source Organism: ${ORGANISM:-<using accession>}"
echo " Source Genome Accession: ${ACCESSION:-<resolved automatically from organism name>}"
echo " Target Pathogen Panel: ${PATHOGENS[*]}"
echo " Threads: $THREADS | Vina exhaustiveness: $EXHAUSTIVENESS | num_modes: $NUM_MODES | cpu: $CPU"
echo "==============================================================================="

echo ""; echo "########## PHASE 1: reference genome acquisition ##########"
if [[ -n "$ACCESSION" ]]; then
    python3 "$SCRIPT_DIR/scripts/01_fetch_genome.py" --accession "$ACCESSION" --outdir results/genomes
else
    python3 "$SCRIPT_DIR/scripts/01_fetch_genome.py" --organism "$ORGANISM" --outdir results/genomes
fi

echo ""; echo "########## PHASE 2: antiSMASH mining ##########"
bash "$SCRIPT_DIR/scripts/02_run_antismash.sh" results/genomes results/antismash "$THREADS"

echo ""; echo "########## PHASE 3: parse antiSMASH output ##########"
ANTISMASH_DIR=$(cat results/antismash/latest_antismash_dir.txt)
python3 "$SCRIPT_DIR/scripts/03_parse_antismash.py" --antismash-dir "$ANTISMASH_DIR" \
    --out results/antismash/phase3_compound_list.json

echo ""; echo "########## PHASE 4: automated ligand (compound) retrieval ##########"
python3 "$SCRIPT_DIR/scripts/04_get_ligands.py" --compound-list results/antismash/phase3_compound_list.json \
    --outdir results/ligands

echo ""; echo "########## PHASE 5: ligand preparation (PDBQT) ##########"
bash "$SCRIPT_DIR/scripts/05_prep_ligands.sh" results/ligands results/ligands_pdbqt

echo ""; echo "########## PHASE 6: retrieve validated pathogen receptors ##########"
python3 "$SCRIPT_DIR/scripts/06_get_receptors.py" --pathogens "${PATHOGENS[@]}" \
    --config "$SCRIPT_DIR/config/pathogen_targets.yaml" --outdir results/receptors_raw

echo ""; echo "########## PHASE 6b: predict missing pathogen receptors (ColabFold) ##########"
ORGANISM_LOOKUP=$(build_organism_lookup)
echo "Organism lookup passed to receptor prediction:"
echo "$ORGANISM_LOOKUP"

python3 "$SCRIPT_DIR/scripts/06b_predict_receptors.py" \
    --manifest results/receptors_raw/receptor_manifest.json \
    --organism-lookup "$ORGANISM_LOOKUP" \
    --outdir results/receptors_raw \
    || echo "[PHASE 6b] Nothing required prediction, or prediction step reported issues above -- see log."

echo ""; echo "########## PHASE 7: receptor cleaning + PDBQT preparation ##########"
bash "$SCRIPT_DIR/scripts/07_prep_receptors.sh" results/receptors_raw results/receptors_pdbqt

echo ""; echo "########## PHASE 8: fpocket pocket detection + grid generation ##########"
bash "$SCRIPT_DIR/scripts/08_run_fpocket.sh" \
    results/receptors_pdbqt \
    results/receptors_pdbqt \
    results/fpocket
    

echo ""; echo "########## PHASE 9: batch docking ##########"
bash "$SCRIPT_DIR/scripts/09_run_docking.sh" results/ligands_pdbqt results/receptors_pdbqt \
    results/fpocket results/docking "$EXHAUSTIVENESS" "$NUM_MODES" "$CPU"

echo ""; echo "########## PHASE 10: PLIP interaction profiling ##########"
python3 "$SCRIPT_DIR/scripts/10_run_plip.py" \
    results/docking \
    results/receptors_raw \
    results/ligands \
    results/plip


echo ""; echo "########## PHASE 11: parse PLIP output ##########"
python3 "$SCRIPT_DIR/scripts/11_parse_plip.py" --plip-dir results/plip --out results/plip/interaction_summary.csv

echo ""; echo "########## PHASE 12: ADMET prediction ##########"
python3 "$SCRIPT_DIR/scripts/12_run_admet.py" --resolved-compounds results/ligands/phase4_resolved_compounds.json \
    --out results/admet/admet_results.csv

echo ""; echo "########## PHASE 13: figure and table generation (15 figures + 5 tables) ##########"
python3 "$SCRIPT_DIR/scripts/13_generate_figures.py" --docking results/docking/summary.csv \
    --interactions results/plip/interaction_summary.csv --admet results/admet/admet_results.csv \
    --figures-dir results/figures --tables-dir results/tables

echo ""; echo "########## PHASE 14: auto-generated interpretation report ##########"
python3 "$SCRIPT_DIR/scripts/14_generate_interpretation.py" --tables-dir results/tables \
    --out results/reports/interpretation.md

echo ""; echo "########## PHASE 15: manuscript/supplementary package ##########"
python3 "$SCRIPT_DIR/scripts/15_generate_manuscript_package.py" --tables-dir results/tables \
    --figures-dir results/figures --out-dir results/reports/manuscript_package

echo ""; echo "########## PHASE 16: reproducibility report (decision_log.json, workflow_report.md, pipeline_metadata.json) ##########"
python3 "$SCRIPT_DIR/scripts/16_generate_reproducibility_report.py" \
    --organism "${ORGANISM:-unknown}" --accession "${ACCESSION:-unknown}" --pathogens "${PATHOGENS[@]}"

echo ""
echo "==============================================================================="
echo " PKS_DOCK PIPELINE COMPLETED SUCCESSFULLY $(date)"
echo " Outputs:"
echo "   results/tables/         (5 CSV tables)"
echo "   results/figures/        (15 PNG figures, 300 DPI)"
echo "   results/reports/interpretation.md         (auto-written prose interpretation)"
echo "   results/reports/manuscript_package/       (Excel workbook + figures + methods stub)"
echo "   results/reports/decision_log.json         (every automated decision, machine-readable)"
echo "   results/reports/workflow_report.md        (same decisions, human-readable narrative)"
echo "   results/reports/pipeline_metadata.json    (software version, git commit, run parameters)"
echo "==============================================================================="
} 2>&1 | tee "$LOGFILE"

