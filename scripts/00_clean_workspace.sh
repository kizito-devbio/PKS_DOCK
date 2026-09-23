#!/usr/bin/env bash
set -euo pipefail

echo "======================================"
echo " PKS_DOCK Workspace Cleanup Utility"
echo "======================================"
echo

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "Project:"
echo "$PROJECT_ROOT"
echo

##############################################
# Remove generated pipeline outputs
##############################################

echo "[1/8] Removing generated results..."

rm -rf results/antismash/*
rm -rf results/genomes/*
rm -rf results/ligands/*
rm -rf results/ligands_pdbqt/*
rm -rf results/receptors_raw/*
rm -rf results/receptors_pdbqt/*
rm -rf results/fpocket/*
rm -rf results/grid_configs/*
rm -rf results/docking/*
rm -rf results/plip/*
rm -rf results/admet/*
rm -rf results/figures/*
rm -rf results/manuscript/*
rm -rf results/reports/*
rm -f results/*.json
rm -f results/*.csv
rm -f results/*.txt
rm -f results/*.xlsx
rm -f results/*.zip

##############################################
# Remove temporary cache
##############################################

echo "[2/8] Removing caches..."

rm -rf cache/http/*
rm -rf .pytest_cache
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -name "*.pyc" -delete
find . -name "*.pyo" -delete

##############################################
# Remove temporary log files
##############################################

echo "[3/8] Removing logs..."

rm -rf logs/*
rm -f *.log
rm -f log.txt

##############################################
# Remove temporary receptor files
##############################################

echo "[4/8] Removing temporary PDB files..."

find results -name "*_fixed.pdb" -delete
find results -name "*_altloc_resolved.pdb" -delete
find results -name "*_repair_report.json" -delete

##############################################
# Remove temporary downloads
##############################################

echo "[5/8] Removing downloaded archives..."

find results -name "*.gz" -delete
find results -name "*.zip" -delete

##############################################
# Remove fpocket outputs if present
##############################################

echo "[6/8] Removing fpocket outputs..."

find results -type d -name "*_out" -exec rm -rf {} + 2>/dev/null || true

##############################################
# Remove miscellaneous temporary files
##############################################

echo "[7/8] Removing miscellaneous temporary files..."

find . -name ".DS_Store" -delete
find . -name "*~" -delete
find . -name "*.bak" -delete
find . -name "*.tmp" -delete

##############################################
# Recreate expected directory structure
##############################################

echo "[8/8] Recreating clean output directories..."

mkdir -p results/{antismash,genomes,ligands,ligands_pdbqt,receptors_raw,receptors_pdbqt,fpocket,grid_configs,docking,plip,admet,figures,manuscript,reports}

echo
echo "======================================"
echo " Workspace cleaned successfully."
echo "======================================"
echo
echo "Source code was NOT deleted."
echo "Configuration files were NOT deleted."
echo "Scripts were NOT deleted."
echo "Environment files were NOT deleted."
echo "Ready for a fresh pipeline run."
