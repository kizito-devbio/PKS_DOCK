#!/usr/bin/env bash
# =============================================================================
# Phase 5 - Prepare ligands for downstream docking and interaction analysis.
#
# Inputs:
#   results/ligands/*.sdf
#
# Outputs:
#   results/ligands_pdbqt/*.pdbqt   (AutoDock Vina)
#   results/ligands_pdb/*.pdb       (PLIP template)
#
# The original SDF is converted twice:
#   1. SDF -> PDBQT using Meeko (for docking)
#   2. SDF -> PDB using Open Babel (for PLIP)
#
# The PDB is NOT used for docking. It is preserved so that Phase 10 can build
# chemically correct protein-ligand complexes without converting PDBQT back
# into PDB (which fails for some Meeko atom types such as CG0/G0).
# =============================================================================

set -euo pipefail

LIGAND_DIR="${1:-results/ligands}"
PDBQT_OUTDIR="${2:-results/ligands_pdbqt}"
PDB_OUTDIR="${3:-results/ligands_pdb}"

mkdir -p "$PDBQT_OUTDIR"
mkdir -p "$PDB_OUTDIR"

echo "[PHASE 5] Preparing ligands for docking and PLIP..."

shopt -s nullglob

count=0

for sdf in "$LIGAND_DIR"/*.sdf; do

    name=$(basename "$sdf" .sdf)

    pdbqt="$PDBQT_OUTDIR/${name}.pdbqt"
    pdb="$PDB_OUTDIR/${name}.pdb"

    echo "[PHASE 5] ------------------------------------------------------------"
    echo "[PHASE 5] Ligand: $name"

    #
    # Generate docking-ready PDBQT
    #
    mk_prepare_ligand.py \
        -i "$sdf" \
        -o "$pdbqt" \
        2>&1 | sed 's/^/[PHASE 5]   Meeko: /'

    if [[ ! -s "$pdbqt" ]]; then
        echo "[PHASE 5] ERROR: Failed to generate:"
        echo "[PHASE 5]        $pdbqt"
        exit 1
    fi

    #
    # Generate canonical PDB directly from the original SDF.
    #
    obabel \
        "$sdf" \
        -O "$pdb" \
        2>&1 | sed 's/^/[PHASE 5]   OpenBabel: /'

    if [[ ! -s "$pdb" ]]; then
        echo "[PHASE 5] ERROR: Failed to generate:"
        echo "[PHASE 5]        $pdb"
        exit 1
    fi

    echo "[PHASE 5]   PDBQT : $pdbqt"
    echo "[PHASE 5]   PDB   : $pdb"

    count=$((count + 1))

done

echo "[PHASE 5] ------------------------------------------------------------"
echo "[PHASE 5] Successfully prepared $count ligand(s)."
echo "[PHASE 5] Docking files : $PDBQT_OUTDIR"
echo "[PHASE 5] PLIP templates: $PDB_OUTDIR"
