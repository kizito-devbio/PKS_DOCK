#!/usr/bin/env bash
# Phase 2 - antiSMASH mining for PKS-I biosynthetic gene clusters.
# Reads the manifest written by 01_fetch_genome.py -- no manual re-typing of
# the accession or genome path is needed here.
set -euo pipefail

GENOME_DIR="${1:-results/genomes}"
OUTDIR="${2:-results/antismash}"
THREADS="${3:-4}"

MANIFEST="$GENOME_DIR/genome_manifest.txt"
if [ ! -f "$MANIFEST" ]; then
    echo "[PHASE 2] ERROR: $MANIFEST not found. Run 01_fetch_genome.py first." >&2
    exit 1
fi

# shellcheck disable=SC1090
source <(sed 's/^/export /' "$MANIFEST")   # exports accession=... gbk_path=...

echo "[PHASE 2] antiSMASH mining for PKS-I clusters"
echo "[PHASE 2] Genome: $gbk_path (accession $accession)"
echo "[PHASE 2] Threads: $THREADS"

mkdir -p "$OUTDIR"

antismash "$gbk_path" \
    --output-dir "$OUTDIR/${accession}_antismash" \
    --genefinding-tool prodigal \
    --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go \
    --cpus "$THREADS"

echo "[PHASE 2] antiSMASH run finished. Results in $OUTDIR/${accession}_antismash"
echo "$OUTDIR/${accession}_antismash" > "$OUTDIR/latest_antismash_dir.txt"
echo "[PHASE 2] Wrote $OUTDIR/latest_antismash_dir.txt for Phase 3 to consume automatically."
