#!/usr/bin/env bash
# Phase 8 - Binding pocket detection with fpocket.
#
# Standard, defensible workflow: protein -> fpocket -> rank pockets by
# druggability -> pick the best one -> derive grid box from ITS coordinates.
# No fixed/generic 40x40x40 boxes, no centering on a ligand that was
# already stripped out.
#
# CORRECTED LOGIC:
# results/receptors_pdbqt/*.pdbqt is the SOURCE OF TRUTH for which
# receptors successfully passed Phase 7. We only run fpocket on the
# *_fixed.pdb structure that corresponds to a receptor that actually
# has a .pdbqt file. This prevents fpocket from silently processing
# nothing (previous bug: the loop glob'd *_fixed.pdb directly instead
# of using the .pdbqt files as the authority list) and prevents it
# from running on stale/failed/orphaned _fixed.pdb intermediates.

set -euo pipefail

RECEPTOR_DIR="${1:-results/receptors_pdbqt}"   # source of truth: *.pdbqt from Phase 7
FIXED_DIR="${2:-results/receptors_pdbqt}"      # where *_fixed.pdb structures live
OUTDIR="${3:-results/fpocket}"

mkdir -p "$OUTDIR"

echo "[PHASE 8] Starting fpocket pocket detection"
echo "[PHASE 8] Receptor authority dir : $RECEPTOR_DIR"
echo "[PHASE 8] Fixed-PDB source dir   : $FIXED_DIR"
echo "[PHASE 8] Output dir             : $OUTDIR"

shopt -s nullglob

PDBQT_FILES=("$RECEPTOR_DIR"/*.pdbqt)

if [[ ${#PDBQT_FILES[@]} -eq 0 ]]; then
    echo "[PHASE 8] ERROR: No prepared receptors (*.pdbqt) found in $RECEPTOR_DIR"
    echo "[PHASE 8] ERROR: Check that Phase 7 completed successfully before running Phase 8."
    exit 1
fi

echo "[PHASE 8] Found ${#PDBQT_FILES[@]} receptor(s) that passed Phase 7."

n=0
failed=0

for pdbqt in "${PDBQT_FILES[@]}"; do

    name=$(basename "$pdbqt" .pdbqt)

    echo "[PHASE 8] ---"
    echo "[PHASE 8] Processing receptor: $name"

    FIXED_PDB="$FIXED_DIR/${name}_fixed.pdb"

    if [[ ! -f "$FIXED_PDB" ]]; then
        echo "[PHASE 8] WARNING: Missing fixed PDB for $name (expected: $FIXED_PDB)"
        echo "[PHASE 8] WARNING: Skipping -- this receptor has a .pdbqt but no matching Phase 7 structure file."
        failed=$((failed+1))
        continue
    fi

    echo "[PHASE 8] Detecting pockets: $name (input: $FIXED_PDB)"

    fpocket -f "$FIXED_PDB" 2>&1 | sed "s/^/[PHASE 8]   fpocket: /"

    fpocket_out="${FIXED_PDB%.pdb}_out"

    if [[ -d "$fpocket_out" ]]; then
        target_dir="$OUTDIR/${name}"
        rm -rf "$target_dir"
        mv "$fpocket_out" "$target_dir"
        echo "[PHASE 8] SUCCESS: pocket detection complete -> $target_dir"
        n=$((n+1))
    else
        echo "[PHASE 8] WARNING: fpocket produced no output directory for $name -- check fpocket install/PATH."
        failed=$((failed+1))
    fi

done

echo ""
echo "===================================="
echo "Phase 8 Summary"
echo "===================================="
echo "Receptors passing Phase 7 : ${#PDBQT_FILES[@]}"
echo "Successful fpocket runs   : $n"
echo "Failed/skipped runs       : $failed"
echo "Output directory          : $OUTDIR"
echo "===================================="

if [[ "$n" -eq 0 ]]; then
    echo "[PHASE 8] ERROR: No pockets detected for any receptor. Aborting before grid generation."
    exit 1
fi

echo "[PHASE 8] fpocket run complete for $n receptor(s). Now generating grid configs..."

python3 "$(dirname "$0")/08b_generate_grid_configs.py" \
    --fpocket-dir "$OUTDIR" \
    --outdir "$OUTDIR" \
    --receptors-pdbqt-dir "$RECEPTOR_DIR"

echo "[PHASE 8] Complete."


