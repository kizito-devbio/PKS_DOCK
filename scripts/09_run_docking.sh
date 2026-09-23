#!/usr/bin/env bash
# ==============================================================================
# Phase 9 - Batch docking: every prepared ligand against every prepared
# receptor, using the grid box fpocket detected for that receptor (Phase 8b).
#
# Inputs   : results/ligands_pdbqt/*.pdbqt, results/receptors_pdbqt/*.pdbqt,
#            results/fpocket/<receptor>.box.txt
# Outputs  : results/docking/<ligand>__<receptor>.pdbqt
#            results/docking/<ligand>__<receptor>.log
#            results/docking/summary.csv
#            results/docking/failed_runs.csv
#
# This script's ONLY responsibility is running AutoDock Vina and recording
# results. PDB generation / PLIP merging now happens in the rewritten Phase 10
# (Python), which reads the original SDF + docked PDBQT directly. Nothing in
# this script's inputs, outputs, CLI signature, or CSV schema has changed, so
# run_pipeline.sh and all downstream phases require no modification.
#
# AFFINITY PARSING NOTE (verified fix, kept from previous draft):
#   Vina's log has a "-----+------------+----------+----------" separator
#   line followed immediately by the data rows. Grabbing "-A1" after the
#   header line that contains "mode |" lands on the units row
#   ("(kcal/mol)"), not a data row. The correct approach is to find the
#   "-----" separator and take the *next* line's second field:
#     awk '/^-----/{found=1; next} found && NF>=2 {print $2; exit}'
#   This has been verified against a synthetic Vina log and is preserved
#   unchanged below -- only wrapped with stricter validation.
# ==============================================================================
set -euo pipefail
IFS=$'\n\t'

# ------------------------------------------------------------------------
# CLI arguments (unchanged signature / defaults -- do not reorder or rename)
# ------------------------------------------------------------------------
LIGAND_DIR="${1:-results/ligands_pdbqt}"
RECEPTOR_DIR="${2:-results/receptors_pdbqt}"
GRID_DIR="${3:-results/fpocket}"
OUTDIR="${4:-results/docking}"
EXHAUSTIVENESS="${5:-16}"
NUM_MODES="${6:-9}"
CPU="${7:-8}"

readonly LIGAND_DIR RECEPTOR_DIR GRID_DIR OUTDIR EXHAUSTIVENESS NUM_MODES CPU

SUMMARY_CSV="$OUTDIR/summary.csv"
FAILED_CSV="$OUTDIR/failed_runs.csv"
readonly SUMMARY_CSV FAILED_CSV

SCRIPT_START_EPOCH=$(date +%s)

# ------------------------------------------------------------------------
# Logging helpers
# ------------------------------------------------------------------------
log() {
    # log <message...>  -- prefixes with phase tag + timestamp
    printf '[PHASE 9] [%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

warn() {
    log "WARNING: $*"
}

die() {
    log "FATAL: $*"
    exit 1
}

elapsed_since_start() {
    local now
    now=$(date +%s)
    printf '%ds' "$(( now - SCRIPT_START_EPOCH ))"
}

hms() {
    # hms <seconds> -> HH:MM:SS
    local total="$1"
    printf '%02d:%02d:%02d' $(( total/3600 )) $(( (total%3600)/60 )) $(( total%60 ))
}

# ------------------------------------------------------------------------
# CSV field escaping (defensive -- filenames should not contain commas, but
# quoting a bare reason field costs nothing and prevents corrupt rows).
# ------------------------------------------------------------------------
csv_field() {
    local field="$1"
    if [[ "$field" == *,* || "$field" == *\"* ]]; then
        field="${field//\"/\"\"}"
        printf '"%s"' "$field"
    else
        printf '%s' "$field"
    fi
}

record_failure() {
    # record_failure <ligand> <receptor> <reason>
    local ligand="$1" receptor="$2" reason="$3"
    printf '%s,%s,%s\n' \
        "$(csv_field "$ligand")" \
        "$(csv_field "$receptor")" \
        "$(csv_field "$reason")" >> "$FAILED_CSV"
}

# ------------------------------------------------------------------------
# Affinity comparison without depending on bc (awk is already a hard
# dependency of this script via the parser below).
# ------------------------------------------------------------------------
affinity_is_better() {
    # affinity_is_better <candidate> <current_best> -> exit 0 if candidate < current_best
    awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'
}

is_number() {
    [[ "$1" =~ ^-?[0-9]+([.][0-9]+)?$ ]]
}

# ------------------------------------------------------------------------
# Parse the best (first, i.e. top-ranked) binding affinity from a Vina log.
# Kept logic identical to the verified fix; only wrapped for validation.
# ------------------------------------------------------------------------
parse_best_affinity() {
    local logfile="$1"
    awk '/^-----/{found=1; next} found && NF>=2 {print $2; exit}' "$logfile"
}

# ------------------------------------------------------------------------
# Extract grid center x,y,z from an fpocket box config file in a single
# awk pass (no grep|awk|paste pipeline, no extra subprocesses).
# ------------------------------------------------------------------------
parse_grid_center() {
    local box_config="$1"
    awk -F'= ' '
        /^center_x/ { x=$2 }
        /^center_y/ { y=$2 }
        /^center_z/ { z=$2 }
        END {
            if (x=="" || y=="" || z=="")
                exit 1
            print x","y","z
        }
    ' "$box_config"
}

# ------------------------------------------------------------------------
# Validate a completed docking run's outputs before trusting them.
# Prints the parsed affinity on stdout on success; returns non-zero on any
# validation failure (caller is responsible for logging/recording reason).
# ------------------------------------------------------------------------
validate_and_extract() {
    local outfile="$1" logfile="$2"

    if [[ ! -f "$outfile" ]]; then
        echo "output PDBQT missing: $outfile"
        return 1
    fi
    if [[ ! -s "$outfile" ]]; then
        echo "output PDBQT is empty: $outfile"
        return 1
    fi
    if [[ ! -f "$logfile" ]]; then
        echo "vina log missing: $logfile"
        return 1
    fi

    local best
    best=$(parse_best_affinity "$logfile") || true
    if [[ -z "$best" ]] || ! is_number "$best"; then
        echo "could not parse a numeric binding affinity from $logfile"
        return 1
    fi

    printf '%s' "$best"
    return 0
}

# ------------------------------------------------------------------------
# Preflight checks
# ------------------------------------------------------------------------
command -v vina  >/dev/null 2>&1 || die "vina executable not found on PATH"
command -v awk   >/dev/null 2>&1 || die "awk not found on PATH"
command -v grep  >/dev/null 2>&1 || die "grep not found on PATH"

[[ -d "$LIGAND_DIR"   ]] || die "ligand directory not found: $LIGAND_DIR"
[[ -d "$RECEPTOR_DIR" ]] || die "receptor directory not found: $RECEPTOR_DIR"
[[ -d "$GRID_DIR"     ]] || die "grid/fpocket directory not found: $GRID_DIR"

mkdir -p "$OUTDIR"

echo "ligand,receptor,pathogen,best_affinity_kcal_mol,grid_center_x,grid_center_y,grid_center_z" > "$SUMMARY_CSV"
echo "ligand,receptor,reason" > "$FAILED_CSV"

log "Batch docking: every ligand against every prepared receptor"
log "Vina parameters: exhaustiveness=$EXHAUSTIVENESS, num_modes=$NUM_MODES, cpu=$CPU"

# ------------------------------------------------------------------------
# Build file lists via globs (never parse ls output).
# ------------------------------------------------------------------------
shopt -s nullglob
ligands=("$LIGAND_DIR"/*.pdbqt)
receptors=("$RECEPTOR_DIR"/*.pdbqt)
shopt -u nullglob

n_lig=${#ligands[@]}
n_rec=${#receptors[@]}
n_total=$(( n_lig * n_rec ))

[[ "$n_lig" -gt 0 ]] || die "no ligand *.pdbqt files found in $LIGAND_DIR"
[[ "$n_rec" -gt 0 ]] || die "no receptor *.pdbqt files found in $RECEPTOR_DIR"

log "$n_lig ligand(s) x $n_rec receptor(s) = $n_total docking run(s) planned"

# ------------------------------------------------------------------------
# Main docking loop
# ------------------------------------------------------------------------
best_overall=""
best_pair=""
n_success=0
n_failed=0
n_skipped=0
job_idx=0

for lig in "${ligands[@]}"; do
    ligname=$(basename "$lig" .pdbqt)

    for rec in "${receptors[@]}"; do
        recname=$(basename "$rec" .pdbqt)
        job_idx=$(( job_idx + 1 ))

        # Pathogen is whatever follows the final underscore in the receptor
        # name (e.g. MurA_Saureus -> Saureus, LasB_Paeruginosa -> Paeruginosa).
        # Future-proof: no fixed organism list to maintain or fall out of date.
        pathogen="${recname##*_}"

        log "(${job_idx}/${n_total}) Docking ligand ${ligname} against receptor ${recname} (pathogen: ${pathogen})"

        box_config="$GRID_DIR/${recname}.box.txt"
        if [[ ! -f "$box_config" ]]; then
            warn "no fpocket grid config for $recname (expected $box_config) -- skipping this pair."
            record_failure "$ligname" "$recname" "missing_grid_config:${box_config}"
            n_skipped=$(( n_skipped + 1 ))
            continue
        fi

        outfile="$OUTDIR/${ligname}__${recname}.pdbqt"
        logfile="$OUTDIR/${ligname}__${recname}.log"

        job_start=$(date +%s)

        if ! vina --receptor "$rec" --ligand "$lig" \
                  --config "$box_config" \
                  --exhaustiveness "$EXHAUSTIVENESS" --num_modes "$NUM_MODES" --cpu "$CPU" \
                  --out "$outfile" > "$logfile" 2>&1; then
            warn "Vina failed for $ligname vs $recname -- see $logfile"
            record_failure "$ligname" "$recname" "vina_nonzero_exit"
            n_failed=$(( n_failed + 1 ))
            continue
        fi

        job_elapsed=$(( $(date +%s) - job_start ))

        best=""
        validation_error=""
        if ! best=$(validate_and_extract "$outfile" "$logfile"); then
            validation_error="$best"
            best=""
        fi

        if [[ -z "$best" ]]; then
            warn "output validation failed for $ligname vs $recname: $validation_error"
            record_failure "$ligname" "$recname" "validation_failed:${validation_error}"
            n_failed=$(( n_failed + 1 ))
            continue
        fi

        center=$(parse_grid_center "$box_config") || true
        if [[ -z "$center" ]]; then
            warn "could not parse grid center from $box_config for $recname -- excluding from summary."
            record_failure "$ligname" "$recname" "missing_grid_center:${box_config}"
            n_failed=$(( n_failed + 1 ))
            continue
        fi
        IFS=',' read -r cx cy cz <<< "$center"

        printf '%s,%s,%s,%s,%s,%s,%s\n' \
            "$(csv_field "$ligname")" "$(csv_field "$recname")" "$(csv_field "$pathogen")" \
            "$best" "$cx" "$cy" "$cz" >> "$SUMMARY_CSV"

        n_success=$(( n_success + 1 ))
        log "Result: ${best} kcal/mol (job took ${job_elapsed}s)"

        if [[ -z "$best_overall" ]] || affinity_is_better "$best" "$best_overall"; then
            best_overall="$best"
            best_pair="$ligname vs $recname"
            log ">>> New best affinity so far: ${best_overall} kcal/mol (${best_pair})"
        fi
    done
done

# ------------------------------------------------------------------------
# Final statistics
# ------------------------------------------------------------------------
total_elapsed=$(( $(date +%s) - SCRIPT_START_EPOCH ))

log "===================================================================="
log "Phase 9 complete."
log "Total jobs planned : $n_total"
log "Successful         : $n_success"
log "Failed             : $n_failed"
log "Skipped (no grid)  : $n_skipped"
if [[ -n "$best_pair" ]]; then
    log "Best docking pair  : $best_pair"
    log "Best affinity      : ${best_overall} kcal/mol"
else
    log "Best docking pair  : none (no successful runs)"
fi
log "Runtime            : $(hms "$total_elapsed") ($(elapsed_since_start))"
log "Summary table       : $SUMMARY_CSV"
log "Failed runs table    : $FAILED_CSV"
log "===================================================================="

if [[ "$n_success" -eq 0 ]]; then
    die "no docking runs succeeded -- aborting pipeline"
fi

exit 0

