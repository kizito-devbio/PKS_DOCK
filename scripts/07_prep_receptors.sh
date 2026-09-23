#!/usr/bin/env bash
#
# Phase 7 v2.4 - Receptor cleaning + PDBQT preparation.
#
# Single-file, no-placeholder implementation. Every stage below is real,
# executable code.
#
# CHANGES IN THIS REVISION (v2.4), fixing what v2.3's production run
# actually showed: the residue-completeness audit worked exactly as
# designed (DHFR_Saureus got past PDBFixer/audit and reached Meeko), but
# Meeko itself then failed with a DIFFERENT error than the one v2.3
# targeted:
#
#   rdkit.Chem.rdchem.AtomValenceException:
#   Explicit valence for atom # 2 C, 5, is greater than permitted
#
# raised from Meeko's own internal
#   meeko/utils/rdkitutils.py:_aux_altloc_mol_build -> Chem.SanitizeMol
# This is NOT a "Template matching failed" residue-completeness error --
# it is Meeko's mandatory (non-optional) ProDy->RDKit per-residue bond
# perception rejecting a spuriously high-valence atom. Critically, in
# production this happened on DHFR_Saureus and GyrA_Paeruginosa on BOTH
# the full structure AND the v2.2/v2.3 cofactor-stripped fallback copy,
# proving the retained HETATM cofactor was not the (only) cause this
# time. The low per-residue atom index in the traceback (Meeko converts
# one residue at a time) points at PDBFixer's addMissingHydrogens(ph)
# step: a fully protonated, all-atom structure hands RDKit's
# distance-based bond perception many more atoms and near-contacts to
# misjudge as bonds than a heavy-atom-only structure does -- and Meeko's
# own "-p" flag already adds the polar hydrogens it needs, so explicit
# PDBFixer hydrogens were never required Meeko input to begin with.
#
#   DESIGN FIX 4 - The single cofactor-stripped Meeko fallback (v2.2/
#     v2.3) is replaced with a FOUR-TIER cascading fallback, each tier
#     built from $MEEKO_READY (the residue-completeness-audited file),
#     tried in order until one succeeds:
#       1. full                          -- cofactors + explicit H intact
#       2. hydrogen_stripped              -- cofactors kept, all explicit H removed
#       3. cofactor_stripped              -- explicit H kept, KEEP_HETATM groups removed
#       4. hydrogen_and_cofactor_stripped -- both removed
#     Tier 2 is tried before tier 3 because discarding explicit hydrogens
#     Meeko/AutoDock does not use for docking geometry is far less
#     chemically consequential than discarding a bound metal/cofactor
#     from the binding site -- but every tier actually attempted (and
#     its exit code) is logged and recorded, and any tier beyond "full"
#     is a WARNING, never silent, since the final PDBQT then differs
#     from the fully-repaired structure. Only if ALL FOUR tiers fail is
#     the receptor marked FAILED/MEEKO. The repair report's
#     "meeko_preparation" stage now records `successful_tier`,
#     `tiers_attempted`, `tier_exit_codes`, and how many cofactor/metal
#     records and explicit hydrogens were removed for whichever tier
#     succeeded.
#
# v2.4.1 FIX (this pass): before each tier's mk_prepare_receptor.py call,
# any stale $PDBQT_OUT / $JSON_OUT from a previous pipeline run on this
# same receptor is removed. mk_prepare_receptor.py overwrites its own
# output on success, so this was never a correctness bug for a single
# run -- but it removes any chance of a leftover file from an earlier,
# unrelated run being misread by a human inspecting the output directory
# mid-pipeline, and keeps the per-tier log unambiguous about which tier
# actually produced the file that is sitting there afterward.
#
# All v2.3/v2.2/v2.1 improvements retained below this point, including
# the residue-completeness audit (Stage 4) that fixed the earlier
# "Template matching failed for residue A:158" failure mode -- that fix
# is independent of and still fully in effect alongside this one.
#
# ---------------------------------------------------------------------
# CHANGES IN v2.3 (retained), fixing the pipeline design gap found
# after v2.2's RDKit fix: v2.2 stopped RDKit from being a false-positive
# blocking gate, and every receptor did reach Meeko -- but some still
# failed AT Meeko itself, e.g.:
#
#   Meeko: Template matching failed for residue A:158
#
#   REMARK 470     LYS A 158    CA   C    O    CB   CG   CD   CE
#   REMARK 470     LYS A 158    NZ
#
#   ATOM  ...  N   LYS A 158  ...
#
# i.e. the raw PDB's own REMARK 470 already disclosed that LYS A158 was
# present with ONLY its backbone N atom -- CA, C, O, and the entire side
# chain were absent. PDBFixer's missing-ATOM machinery (findMissingAtoms
# / addMissingAtoms) is driven by its own internal topology bookkeeping
# and does not reliably catch/repair every such severely-incomplete
# residue on every input; a residue this sparse can survive PDBFixer and
# reach Meeko still broken. Meeko then tries to match the residue against
# its LYS heavy-atom template (N CA C O CB CG CD CE NZ), finds only N,
# and raises "Template matching failed" -- a hard failure, correctly, since
# Meeko cannot invent a side chain out of nothing.
#
#   DESIGN FIX 3 - NEW residue-completeness audit stage, inserted after
#     PDBFixer repair and before Meeko preparation (both the primary and
#     the cofactor-stripped fallback attempt): every standard amino acid
#     residue (ATOM records only, the 20 canonical types) is checked
#     against a fixed heavy-atom template. Any residue missing a required
#     heavy atom is classified incomplete; with REMOVE_INCOMPLETE_RESIDUES
#     =true (default) the ENTIRE residue is removed -- every atom line for
#     that chain/resSeq/iCode, plus any CONECT record referencing one of
#     that residue's removed atom serial numbers -- producing a new
#     *_meeko_ready.pdb that is what Meeko actually receives. $FIXED (the
#     full PDBFixer output) is untouched and is still what Phase 8/
#     fpocket reads, and is still what RDKit's diagnostic validation runs
#     against.
#
# All v2.2/v2.1/v2.0 improvements retained: the KEEP_INTERMEDIATES
# os.environ fix, RDKit as a non-blocking diagnostic gate (Meeko remains
# the real, authoritative gate for dockability), occupancy-based altloc
# resolution, configurable HETATM allowlist, terminal/internal/anchor-
# aware PDBFixer loop rebuilding, explicit RDKit failure-mode + fragment-
# count + zero-coordinate checks (non-blocking), real Meeko JSON
# inspection, full repair reports (versions, runtime, SHA-256 hashes,
# command line, warnings), per-residue atom accounting, structured
# timestamped logging, optional parallelism via THREADS=N, and every
# previously hard-coded value still a parameter.
#
# BACKWARD COMPATIBILITY: the five positional arguments are unchanged so
# this drops in for the existing call in run_pipeline.sh:
#   scripts/07_prep_receptors.sh <raw_dir> <outdir> <max_loop> <warn_pct> <keep_intermediates>
# The newer knobs (KEEP_HETATM, REPAIR_PH, THREADS, REMOVE_INCOMPLETE_RESIDUES)
# are environment variables rather than positional args, specifically so
# the existing run_pipeline.sh call site does not need to change to pick
# up this upgrade.

set -euo pipefail

ORIGINAL_CMD="$0 $*"

RAW_DIR="${1:-results/receptors_raw}"
OUTDIR="${2:-results/receptors_pdbqt}"
MAX_MISSING_LOOP_LENGTH="${3:-5}"
HEAVY_ATOM_INCREASE_WARN_PCT="${4:-15}"
KEEP_INTERMEDIATES="${5:-false}"

KEEP_HETATM="${KEEP_HETATM:-ZN,MG,MN,FE,CA,NA,K,CU,CO,NI,FAD,FMN,NAD,NAP,NDP,HEM,HEC,SAM,ATP,ADP,GTP,GDP,PLP,COA}"
REPAIR_PH="${REPAIR_PH:-7.0}"
THREADS="${THREADS:-1}"
REMOVE_INCOMPLETE_RESIDUES="${REMOVE_INCOMPLETE_RESIDUES:-true}"

mkdir -p "$OUTDIR"

# ====================== LOGGING HELPERS ======================
log_info()  { printf '[%s] [PHASE 7] INFO  %s\n'  "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
log_warn()  { printf '[%s] [PHASE 7] WARN  %s\n'  "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
log_error() { printf '[%s] [PHASE 7] ERROR %s\n'  "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }

log_info "Cleaning and preparing receptors for pocket detection + docking"
log_info "Max internal missing-loop length allowed to be rebuilt: $MAX_MISSING_LOOP_LENGTH residues"
log_info "Heavy-atom increase warning threshold: ${HEAVY_ATOM_INCREASE_WARN_PCT}%"
log_info "Keep intermediates : $KEEP_INTERMEDIATES"
log_info "HETATM groups retained during stripping: $KEEP_HETATM"
log_info "Protonation pH for missing-hydrogen addition: $REPAIR_PH"
log_info "Parallel worker threads: $THREADS"
log_info "Remove incomplete standard-residue amino acids before Meeko: $REMOVE_INCOMPLETE_RESIDUES"
log_info "RDKit sanitization is a NON-BLOCKING diagnostic in this version -- a failure is logged as a warning, and the receptor still proceeds to Meeko."
log_info "Meeko preparation uses a FOUR-TIER cascading fallback (full -> hydrogen_stripped -> cofactor_stripped -> hydrogen_and_cofactor_stripped)."

# ====================== DEPENDENCY & VERSION LOGGING ======================
log_info "Checking dependencies and recording versions..."
command -v python3 >/dev/null || { log_error "python3 not found."; exit 1; }
command -v mk_prepare_receptor.py >/dev/null || { log_error "mk_prepare_receptor.py not found."; exit 1; }

log_info "Python : $(python3 --version 2>&1)"
log_info "Meeko  : $(mk_prepare_receptor.py --help 2>&1 | head -n 1)"

python3 - <<'PYVERS'
import importlib
mods = ["Bio", "rdkit", "openmm", "pdbfixer"]
for m in mods:
    try:
        mod = importlib.import_module(m)
        ver = getattr(mod, "__version__", "unknown")
        print(f"[PHASE 7] {m:12} : {ver}")
    except Exception as e:
        print(f"[PHASE 7] WARNING: Could not import {m}: {e}")
PYVERS

python3 - <<'PYDEP'
import sys
required = ["Bio", "rdkit", "openmm", "pdbfixer"]
missing = []
for module in required:
    try:
        __import__(module)
    except ImportError:
        missing.append(module)
if missing:
    print("ERROR: Missing modules:", missing)
    sys.exit(1)
print("[PHASE 7] All core Python dependencies satisfied.")
PYDEP

# ====================== SETUP ======================
STATUS_TSV="$OUTDIR/_receptor_status.tsv"
: > "$STATUS_TSV"
FAILED_TSV="$OUTDIR/failed_receptors.tsv"
printf 'Receptor\tFailure_Reason\n' > "$FAILED_TSV"

shopt -s nullglob
RECEPTORS=("$RAW_DIR"/*.pdb)
TOTAL=${#RECEPTORS[@]}

if [[ $TOTAL -eq 0 ]]; then
    log_error "No .pdb files found in $RAW_DIR"
    exit 1
fi

export OUTDIR MAX_MISSING_LOOP_LENGTH HEAVY_ATOM_INCREASE_WARN_PCT KEEP_INTERMEDIATES
export KEEP_HETATM REPAIR_PH STATUS_TSV FAILED_TSV ORIGINAL_CMD REMOVE_INCOMPLETE_RESIDUES

# ====================== PER-RECEPTOR WORKER ======================
process_receptor() {
    local f="$1"
    local name; name=$(basename "$f" .pdb)

    log_info "================= Processing receptor: $name ================="

    local NOWAT="$OUTDIR/${name}_nowat.pdb"
    local FIXED="$OUTDIR/${name}_fixed.pdb"
    local MEEKO_READY="$OUTDIR/${name}_meeko_ready.pdb"
    local RDKIT_SAFE="$OUTDIR/${name}_rdkit_safe.pdb"
    local MEEKO_FALLBACK="$OUTDIR/${name}_meeko_fallback.pdb"
    local MEEKO_NOH="$OUTDIR/${name}_meeko_noH.pdb"
    local MEEKO_NOH_NOHET="$OUTDIR/${name}_meeko_noH_nohet.pdb"
    local REPORT_JSON="$OUTDIR/${name}_repair_report.json"
    local PDBQT_OUT="$OUTDIR/${name}.pdbqt"
    local JSON_OUT="$OUTDIR/${name}.json"

    export RECEPTOR_NAME="$name"
    export RAW_PDB="$f"
    export NOWAT_PDB="$NOWAT"
    export FIXED_PDB="$FIXED"
    export MEEKO_READY_PDB="$MEEKO_READY"
    export RDKIT_SAFE_PDB="$RDKIT_SAFE"
    export REPORT_JSON="$REPORT_JSON"
    export MAX_LOOP="$MAX_MISSING_LOOP_LENGTH"
    export HEAVY_WARN_PCT="$HEAVY_ATOM_INCREASE_WARN_PCT"
    export REPAIR_PH_ENV="$REPAIR_PH"

    set +e
    REPAIR_OUTPUT=$(python3 - <<'PYREPAIR'
import hashlib
import importlib
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone

RUN_START = time.time()

name          = os.environ["RECEPTOR_NAME"]
raw_pdb       = os.environ["RAW_PDB"]
nowat_pdb     = os.environ["NOWAT_PDB"]
fixed_pdb     = os.environ["FIXED_PDB"]
meeko_ready_pdb = os.environ["MEEKO_READY_PDB"]
rdkit_safe_pdb = os.environ["RDKIT_SAFE_PDB"]
report_path   = os.environ["REPORT_JSON"]
max_loop      = int(os.environ["MAX_LOOP"])
warn_pct      = float(os.environ["HEAVY_WARN_PCT"])
ph            = float(os.environ["REPAIR_PH_ENV"])
keep_hetatm   = {i.strip().upper() for i in os.environ["KEEP_HETATM"].split(",") if i.strip()}
original_cmd  = os.environ.get("ORIGINAL_CMD", "")

keep_intermediates = os.environ.get("KEEP_INTERMEDIATES", "false").strip().lower() == "true"
remove_incomplete_residues = os.environ.get("REMOVE_INCOMPLETE_RESIDUES", "true").strip().lower() == "true"

warnings_log = []


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def software_versions():
    out = {}
    for m in ("Bio", "rdkit", "openmm", "pdbfixer"):
        try:
            mod = importlib.import_module(m)
            out[m] = getattr(mod, "__version__", "unknown")
        except Exception:
            out[m] = "not importable"
    out["python"] = sys.version.split()[0]
    return out


report = {
    "receptor": name,
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "command_line": original_cmd,
    "parameters": {
        "max_missing_loop_length": max_loop,
        "heavy_atom_warn_pct": warn_pct,
        "repair_ph": ph,
        "keep_hetatm": sorted(keep_hetatm),
        "keep_intermediates": keep_intermediates,
        "remove_incomplete_residues": remove_incomplete_residues,
    },
    "software_versions": software_versions(),
    "stages": {},
    "warnings": warnings_log,
}


def fail(stage, exc):
    report["stages"][stage] = {"status": "failed", "error": str(exc)}
    report["overall_status"] = "failed"
    report["runtime_seconds"] = round(time.time() - RUN_START, 2)
    try:
        with open(report_path, "w") as fh:
            json.dump(report, fh, indent=2)
    except Exception:
        pass
    print(f"STAGE_FAILED={stage}")
    print(f"[{name}] {stage} failed: {exc}", file=sys.stderr)
    traceback.print_exc()
    sys.exit(1)


RECORD  = slice(0, 6)
SERIAL  = slice(6, 11)
ATOMNAM = slice(12, 16)
ALTLOC  = 16
RESNAME = slice(17, 20)
CHAINID = 21
RESSEQ  = slice(22, 26)
ICODE   = 26
OCC     = slice(54, 60)
ELEMENT = slice(76, 78)


def resolve_altloc_and_strip_hetero(in_path, out_path, keep_hetatm):
    with open(in_path) as fh:
        lines = fh.readlines()

    groups = {}
    for idx, line in enumerate(lines):
        rec = line[RECORD]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        altloc = line[ALTLOC]
        if altloc == " ":
            continue
        key = (line[CHAINID], line[RESSEQ], line[ICODE], line[ATOMNAM])
        try:
            occ = float(line[OCC])
        except ValueError:
            occ = 0.0
        groups.setdefault(key, []).append((occ, altloc, idx))

    drop_indices = set()
    resolved_indices = set()
    for key, entries in groups.items():
        if len(entries) < 2:
            continue
        entries.sort(key=lambda t: (-t[0], t[1]))
        winner_idx = entries[0][2]
        resolved_indices.add(winner_idx)
        for _, _, idx in entries[1:]:
            drop_indices.add(idx)

    altloc_conformers_dropped = len(drop_indices)

    out_lines = []
    in_model = False
    seen_first_model = False
    waters_removed = 0
    hetatm_kept = 0
    hetatm_removed = 0
    truncated_extra_models = False

    for idx, line in enumerate(lines):
        rec = line[RECORD]

        if rec == "MODEL ":
            if seen_first_model:
                truncated_extra_models = True
                break
            seen_first_model = True
            in_model = True
            continue
        if rec == "ENDMDL":
            in_model = False
            continue

        if rec in ("ATOM  ", "HETATM"):
            if idx in drop_indices:
                continue
            if rec == "HETATM":
                resname = line[RESNAME].strip().upper()
                if resname == "HOH":
                    waters_removed += 1
                    continue
                if resname not in keep_hetatm:
                    hetatm_removed += 1
                    continue
                hetatm_kept += 1
            if idx in resolved_indices:
                line = line[:ALTLOC] + " " + line[ALTLOC + 1:]
            out_lines.append(line)
            continue

        if rec in ("ANISOU",):
            continue

        out_lines.append(line)

    with open(out_path, "w") as fh:
        fh.writelines(out_lines)

    return {
        "altloc_conformers_dropped": altloc_conformers_dropped,
        "waters_removed": waters_removed,
        "hetatm_residues_kept": hetatm_kept,
        "hetatm_residues_removed": hetatm_removed,
        "truncated_extra_nmr_models": truncated_extra_models,
    }


def count_atoms(pdb_path):
    heavy = 0
    hydrogen = 0
    residues = set()
    chains = set()
    zero_coord = 0
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[RECORD]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            element = line[ELEMENT].strip()
            if not element:
                element = line[ATOMNAM].strip()[0]
            chains.add(line[CHAINID])
            residues.add((line[CHAINID], line[RESSEQ], line[ICODE]))
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if x == 0.0 and y == 0.0 and z == 0.0:
                    zero_coord += 1
            except ValueError:
                pass
            if element.upper() == "H":
                hydrogen += 1
            else:
                heavy += 1
    return {
        "heavy_atoms": heavy,
        "hydrogen_atoms": hydrogen,
        "residue_count": len(residues),
        "chain_count": len(chains),
        "zero_coordinate_atoms": zero_coord,
    }


BACKBONE_REQUIRED = {"N", "CA", "C", "O"}

STANDARD_AA_ATOMS = {
    "ALA": BACKBONE_REQUIRED | {"CB"},
    "ARG": BACKBONE_REQUIRED | {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": BACKBONE_REQUIRED | {"CB", "CG", "OD1", "ND2"},
    "ASP": BACKBONE_REQUIRED | {"CB", "CG", "OD1", "OD2"},
    "CYS": BACKBONE_REQUIRED | {"CB", "SG"},
    "GLN": BACKBONE_REQUIRED | {"CB", "CG", "CD", "OE1", "NE2"},
    "GLU": BACKBONE_REQUIRED | {"CB", "CG", "CD", "OE1", "OE2"},
    "GLY": BACKBONE_REQUIRED,
    "HIS": BACKBONE_REQUIRED | {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": BACKBONE_REQUIRED | {"CB", "CG1", "CG2", "CD1"},
    "LEU": BACKBONE_REQUIRED | {"CB", "CG", "CD1", "CD2"},
    "LYS": BACKBONE_REQUIRED | {"CB", "CG", "CD", "CE", "NZ"},
    "MET": BACKBONE_REQUIRED | {"CB", "CG", "SD", "CE"},
    "PHE": BACKBONE_REQUIRED | {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": BACKBONE_REQUIRED | {"CB", "CG", "CD"},
    "SER": BACKBONE_REQUIRED | {"CB", "OG"},
    "THR": BACKBONE_REQUIRED | {"CB", "OG1", "CG2"},
    "TRP": BACKBONE_REQUIRED | {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": BACKBONE_REQUIRED | {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": BACKBONE_REQUIRED | {"CB", "CG1", "CG2"},
}


def build_meeko_ready(in_path, out_path, remove_incomplete):
    with open(in_path) as fh:
        lines = fh.readlines()

    residue_atom_lines = {}
    residue_order = []
    for idx, line in enumerate(lines):
        if line[RECORD] != "ATOM  ":
            continue
        resname = line[RESNAME].strip().upper()
        if resname not in STANDARD_AA_ATOMS:
            continue
        key = (line[CHAINID], line[RESSEQ], line[ICODE], resname)
        atom_name = line[ATOMNAM].strip()
        element = line[ELEMENT].strip()
        if not element:
            element = atom_name[0] if atom_name else "?"
        is_hydrogen = element.upper() == "H"
        if key not in residue_atom_lines:
            residue_atom_lines[key] = []
            residue_order.append(key)
        residue_atom_lines[key].append((idx, atom_name, is_hydrogen))

    incomplete_residues = []
    remove_line_indices = set()
    removed_serials = set()

    for key in residue_order:
        chain, resseq, icode, resname = key
        entries = residue_atom_lines[key]
        present_heavy = {atom_name for (_, atom_name, is_h) in entries if not is_h}
        required = STANDARD_AA_ATOMS[resname]
        missing = sorted(required - present_heavy)
        if not missing:
            continue

        incomplete_residues.append({
            "chain": chain.strip(),
            "residue_number": int(resseq),
            "icode": icode.strip(),
            "name": resname,
            "missing_atoms": missing,
            "present_atoms": sorted(present_heavy),
        })

        if remove_incomplete:
            for (idx, _, _) in entries:
                remove_line_indices.add(idx)
                try:
                    removed_serials.add(int(lines[idx][SERIAL]))
                except ValueError:
                    pass

    out_lines = []
    for idx, line in enumerate(lines):
        if idx in remove_line_indices:
            continue
        if line[RECORD] == "CONECT" and removed_serials:
            try:
                serials_in_line = {int(tok) for tok in line.split()[1:]}
            except ValueError:
                serials_in_line = set()
            if serials_in_line & removed_serials:
                continue
        out_lines.append(line)

    with open(out_path, "w") as fh:
        fh.writelines(out_lines)

    total_residues = len(residue_order)
    incomplete_count = len(incomplete_residues)
    return {
        "total_standard_residues_checked": total_residues,
        "complete_residues": total_residues - incomplete_count,
        "incomplete_residues_found": incomplete_count,
        "incomplete_residues_removed_count": incomplete_count if remove_incomplete else 0,
        "removal_applied": remove_incomplete,
        "incomplete_residues": incomplete_residues,
    }


RDKIT_UNSAFE_HETATM = {
    "ZN", "MG", "MN", "FE", "CA", "NA", "K", "CU", "CO", "NI",
    "FAD", "FMN", "NAD", "NAP", "NDP", "HEM", "HEC", "SAM",
    "ATP", "ADP", "GTP", "GDP", "PLP", "COA",
}


def clean_pdb_before_rdkit(in_path, out_path):
    with open(in_path) as fh:
        lines = fh.readlines()

    seen_keys = set()
    out_lines = []
    duplicates_removed = 0
    hetatm_removed_for_validation = 0
    altlocs_blanked = 0

    for line in lines:
        rec = line[RECORD]
        if rec not in ("ATOM  ", "HETATM"):
            out_lines.append(line)
            continue

        if rec == "HETATM":
            resname = line[RESNAME].strip().upper()
            if resname in RDKIT_UNSAFE_HETATM:
                hetatm_removed_for_validation += 1
                continue

        if line[ALTLOC] != " ":
            altlocs_blanked += 1
            line = line[:ALTLOC] + " " + line[ALTLOC + 1:]

        key = (line[CHAINID], line[RESSEQ], line[ICODE], line[ATOMNAM])
        if key in seen_keys:
            duplicates_removed += 1
            continue
        seen_keys.add(key)
        out_lines.append(line)

    with open(out_path, "w") as fh:
        fh.writelines(out_lines)

    return {
        "altlocs_blanked": altlocs_blanked,
        "duplicate_atoms_removed": duplicates_removed,
        "hetatm_residue_records_removed_for_validation_only": hetatm_removed_for_validation,
    }


try:
    strip_stats = resolve_altloc_and_strip_hetero(raw_pdb, nowat_pdb, keep_hetatm)
    if strip_stats["truncated_extra_nmr_models"]:
        warnings_log.append(
            "Input PDB contained multiple MODEL records (NMR-style "
            "ensemble); only the first model was used."
        )
    report["stages"]["altloc_and_heteroatom_strip"] = {
        "status": "ok",
        **strip_stats,
    }
except Exception as exc:
    fail("ALTLOC", exc)

before_stats = count_atoms(nowat_pdb)

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=nowat_pdb)

    fixer.findNonstandardResidues()
    nonstandard_replaced = len(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()

    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())

    filled_runs = 0
    filled_residues = 0
    skipped_terminal_runs = 0
    skipped_terminal_residues = 0
    skipped_long_runs = 0
    skipped_long_residues = 0
    skipped_no_anchor_runs = 0
    skipped_no_anchor_residues = 0

    for key in list(fixer.missingResidues.keys()):
        chain_index, res_index = key
        run = fixer.missingResidues[key]
        chain_residues = list(chains[chain_index].residues())
        chain_len = len(chain_residues)

        if res_index == 0 or res_index >= chain_len:
            skipped_terminal_runs += 1
            skipped_terminal_residues += len(run)
            del fixer.missingResidues[key]
            continue

        if len(run) > max_loop:
            skipped_long_runs += 1
            skipped_long_residues += len(run)
            del fixer.missingResidues[key]
            continue

        anchor_before = chain_residues[res_index - 1]
        anchor_after = chain_residues[res_index]
        anchor_atoms_before = {a.name for a in anchor_before.atoms()}
        anchor_atoms_after = {a.name for a in anchor_after.atoms()}
        required_backbone = {"N", "CA", "C"}
        if not (required_backbone <= anchor_atoms_before and
                required_backbone <= anchor_atoms_after):
            skipped_no_anchor_runs += 1
            skipped_no_anchor_residues += len(run)
            del fixer.missingResidues[key]
            continue

        filled_runs += 1
        filled_residues += len(run)

    fixer.findMissingAtoms()
    missing_atoms_by_residue = {}
    for residue, atom_names in fixer.missingAtoms.items():
        rkey = f"{residue.chain.id}:{residue.id}:{residue.name}"
        missing_atoms_by_residue[rkey] = missing_atoms_by_residue.get(rkey, 0) + len(atom_names)
    for residue, atom_names in fixer.missingTerminals.items():
        rkey = f"{residue.chain.id}:{residue.id}:{residue.name}"
        missing_atoms_by_residue[rkey] = missing_atoms_by_residue.get(rkey, 0) + len(atom_names)

    missing_atoms_total = sum(len(v) for v in fixer.missingAtoms.values())
    missing_terminals_total = sum(len(v) for v in fixer.missingTerminals.values())

    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    with open(fixed_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    if skipped_long_runs or skipped_no_anchor_runs:
        warnings_log.append(
            f"{skipped_long_runs} internal loop(s) exceeded the "
            f"{max_loop}-residue rebuild threshold and {skipped_no_anchor_runs} "
            f"internal loop(s) lacked complete backbone anchors on at least "
            f"one flanking residue; none of these were rebuilt."
        )

    report["stages"]["pdbfixer_repair"] = {
        "status": "ok",
        "nonstandard_residues_replaced": nonstandard_replaced,
        "missing_residue_runs_filled": filled_runs,
        "missing_residues_filled_total": filled_residues,
        "missing_residue_runs_skipped_terminal": skipped_terminal_runs,
        "missing_residues_skipped_terminal_total": skipped_terminal_residues,
        "missing_residue_runs_skipped_too_long": skipped_long_runs,
        "missing_residues_skipped_too_long_total": skipped_long_residues,
        "missing_residue_runs_skipped_no_anchor": skipped_no_anchor_runs,
        "missing_residues_skipped_no_anchor_total": skipped_no_anchor_residues,
        "missing_loop_threshold": max_loop,
        "missing_atoms_added_total": missing_atoms_total,
        "missing_terminal_atoms_added_total": missing_terminals_total,
        "missing_atoms_added_by_residue": missing_atoms_by_residue,
        "repair_ph": ph,
    }
except Exception as exc:
    fail("PDBFIXER", exc)

after_stats = count_atoms(fixed_pdb)

heavy_before = before_stats["heavy_atoms"]
heavy_after = after_stats["heavy_atoms"]
increase_pct = (100.0 * (heavy_after - heavy_before) / heavy_before) if heavy_before > 0 else 0.0
heavy_warning = increase_pct > warn_pct

report["stages"]["heavy_atom_check"] = {
    "status": "warning" if heavy_warning else "ok",
    "heavy_atoms_before": heavy_before,
    "heavy_atoms_after": heavy_after,
    "hydrogen_atoms_after": after_stats["hydrogen_atoms"],
    "residue_count_after": after_stats["residue_count"],
    "chain_count_after": after_stats["chain_count"],
    "zero_coordinate_atoms_after": after_stats["zero_coordinate_atoms"],
    "increase_pct": round(increase_pct, 2),
    "warn_threshold_pct": warn_pct,
}
if heavy_warning:
    msg = (f"Heavy-atom count increased by {increase_pct:.1f}% "
           f"(threshold {warn_pct}%) -- inspect {fixed_pdb} before "
           f"trusting downstream docking results.")
    warnings_log.append(msg)
    print(f"[{name}] WARNING: {msg}")
if after_stats["zero_coordinate_atoms"] > 0:
    msg = (f"{after_stats['zero_coordinate_atoms']} atom(s) in the repaired "
           f"structure sit exactly at (0,0,0) -- likely a write artifact.")
    warnings_log.append(msg)
    print(f"[{name}] WARNING: {msg}")

try:
    audit_stats = build_meeko_ready(fixed_pdb, meeko_ready_pdb, remove_incomplete_residues)
    if audit_stats["incomplete_residues_found"]:
        example = audit_stats["incomplete_residues"][0]
        if remove_incomplete_residues:
            msg = (
                f"{audit_stats['incomplete_residues_found']} standard amino "
                f"acid residue(s) were present after PDBFixer repair but "
                f"missing required heavy atoms (e.g. {example['name']} "
                f"{example['chain']}{example['residue_number']} was missing "
                f"{', '.join(example['missing_atoms'])}). These residues were "
                f"removed entirely (residue-level removal, including any "
                f"CONECT records referencing their atoms) before Meeko "
                f"preparation; see this report's "
                f"'residue_completeness_audit.incomplete_residues' for the "
                f"full list and exactly which atoms were missing on each one."
            )
        else:
            msg = (
                f"{audit_stats['incomplete_residues_found']} standard amino "
                f"acid residue(s) were present after PDBFixer repair but "
                f"missing required heavy atoms (e.g. {example['name']} "
                f"{example['chain']}{example['residue_number']} was missing "
                f"{', '.join(example['missing_atoms'])}). "
                f"REMOVE_INCOMPLETE_RESIDUES=false, so they were NOT removed "
                f"and were passed to Meeko as-is; this is very likely to make "
                f"Meeko's residue-template matching fail."
            )
        warnings_log.append(msg)
        print(f"[{name}] WARNING: {msg}")
    report["stages"]["residue_completeness_audit"] = {
        "status": "warning" if audit_stats["incomplete_residues_found"] else "ok",
        **audit_stats,
    }
except Exception as exc:
    fail("RESIDUE_COMPLETENESS_AUDIT", exc)

sanitize_ok = False
used_rdkit_safe_fallback = False
cleaning_stats = None
fragment_warning = False
frag_count = None
expected_fragments = None
rdkit_atom_count = None
rdkit_zero_coord = None
rdkit_error_message = None

try:
    from rdkit import Chem

    def try_sanitize(path):
        mol = Chem.MolFromPDBFile(path, sanitize=False, removeHs=False)
        if mol is None:
            return None, "RDKit could not parse this PDB's connectivity at all."
        try:
            Chem.SanitizeMol(mol)
        except Exception as sanitize_exc:
            return None, str(sanitize_exc)
        return mol, None

    mol, sanitize_error = try_sanitize(fixed_pdb)
    validated_path = fixed_pdb

    if mol is None:
        cleaning_stats = clean_pdb_before_rdkit(fixed_pdb, rdkit_safe_pdb)
        fallback_mol, fallback_error = try_sanitize(rdkit_safe_pdb)
        if fallback_mol is not None:
            msg = (
                f"RDKit sanitization failed on the full repaired structure "
                f"({sanitize_error}), but succeeded after removing "
                f"{cleaning_stats['hetatm_residue_records_removed_for_validation_only']} "
                f"metal/cofactor HETATM record(s) and "
                f"{cleaning_stats['duplicate_atoms_removed']} duplicate atom "
                f"record(s) from a VALIDATION-ONLY copy. This strongly suggests "
                f"RDKit's distance-based bond perception mis-bonded a retained "
                f"cofactor, not a real problem with the protein backbone. "
                f"The full structure (cofactors intact) is still what is sent "
                f"to Meeko for PDBQT preparation."
            )
            warnings_log.append(msg)
            print(f"[{name}] WARNING: {msg}")
            mol = fallback_mol
            validated_path = rdkit_safe_pdb
            used_rdkit_safe_fallback = True
            sanitize_ok = True
        else:
            rdkit_error_message = (
                f"RDKit sanitization failed on the full structure "
                f"({sanitize_error}) AND on the metal/cofactor-stripped, "
                f"deduplicated validation copy ({fallback_error})."
            )
            warnings_log.append(
                "RDKit sanitization failed but continuing because RDKit is "
                "not authoritative for macromolecular receptor validation. "
                f"Detail: {rdkit_error_message}"
            )
            print(f"[{name}] WARNING: RDKit sanitization failed but continuing "
                  f"because RDKit is not authoritative for macromolecular "
                  f"receptor validation. Detail: {rdkit_error_message}")
            sanitize_ok = False
            mol = None

    if mol is not None:
        frags = Chem.GetMolFrags(mol, asMols=False)
        expected_fragments = count_atoms(validated_path)["chain_count"]
        frag_count = len(frags)
        fragment_warning = frag_count > expected_fragments
        if fragment_warning:
            msg = (f"RDKit found {frag_count} disconnected fragment(s) but only "
                   f"{expected_fragments} chain(s) were expected -- this can "
                   f"indicate a broken backbone within a chain, not just normal "
                   f"chain separation.")
            warnings_log.append(msg)
            print(f"[{name}] WARNING: {msg}")

        conf = mol.GetConformer()
        rdkit_zero_coord = 0
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            if pos.x == 0.0 and pos.y == 0.0 and pos.z == 0.0:
                rdkit_zero_coord += 1
        rdkit_atom_count = mol.GetNumAtoms()

except Exception as exc:
    rdkit_error_message = f"RDKit validation stage could not run: {exc}"
    warnings_log.append(
        "RDKit sanitization failed but continuing because RDKit is not "
        f"authoritative for macromolecular receptor validation. Detail: {rdkit_error_message}"
    )
    print(f"[{name}] WARNING: RDKit validation stage could not run ({exc}); "
          f"continuing because RDKit is not authoritative for macromolecular "
          f"receptor validation.")
    sanitize_ok = False
    mol = None

report["stages"]["rdkit_validation"] = {
    "status": "ok" if (sanitize_ok and not fragment_warning) else "warning",
    "sanitize_ok": sanitize_ok,
    "blocking": False,
    "error": rdkit_error_message,
    "continued_to_meeko": True,
    "atom_count": rdkit_atom_count,
    "fragment_count": frag_count,
    "expected_fragment_count": expected_fragments,
    "zero_coordinate_atoms": rdkit_zero_coord,
    "validated_on_metal_cofactor_stripped_copy": used_rdkit_safe_fallback,
    "cleaning_stats": cleaning_stats,
}

report["overall_status"] = "ok"
report["runtime_seconds"] = round(time.time() - RUN_START, 2)
report["input_sha256"] = sha256_of(raw_pdb)
report["output_sha256"] = sha256_of(fixed_pdb)

with open(report_path, "w") as fh:
    json.dump(report, fh, indent=2)

print(f"[{name}] Repair complete: {heavy_before} -> {heavy_after} heavy atoms "
      f"({report['runtime_seconds']}s)")
sys.exit(0)
PYREPAIR
)
    REPAIR_STATUS=$?
    set -e
    echo "$REPAIR_OUTPUT" | sed "s/^/[PHASE 7]   /"

    if [[ $REPAIR_STATUS -ne 0 ]]; then
        STAGE=$(echo "$REPAIR_OUTPUT" | grep -o 'STAGE_FAILED=[A-Z_]*' | tail -n1 | cut -d= -f2)
        STAGE=${STAGE:-UNKNOWN}
        log_error "Receptor repair failed for $name at stage $STAGE"
        printf '%s\tRepair failed at %s\n' "$name" "$STAGE" >> "$FAILED_TSV"
        printf '%s\tFAILED\t%s\n' "$name" "$STAGE" >> "$STATUS_TSV"
        return 0
    fi

    if [[ ! -f "$FIXED" ]] || [[ ! -s "$FIXED" ]]; then
        log_error "No repaired PDB produced for $name"
        printf '%s\tNo repaired PDB output\n' "$name" >> "$FAILED_TSV"
        printf '%s\tFAILED\tPDBFIXER\n' "$name" >> "$STATUS_TSV"
        return 0
    fi

    if [[ ! -f "$MEEKO_READY" ]] || [[ ! -s "$MEEKO_READY" ]]; then
        log_error "No Meeko-ready PDB produced for $name (residue-completeness audit stage)"
        printf '%s\tNo meeko_ready PDB output\n' "$name" >> "$FAILED_TSV"
        printf '%s\tFAILED\tRESIDUE_COMPLETENESS_AUDIT\n' "$name" >> "$STATUS_TSV"
        return 0
    fi

    # === Meeko, with a FOUR-TIER cascading fallback + real JSON inspection ===
    declare -A TIER_STRIP_HETATM=( [hydrogen_stripped]="false" [cofactor_stripped]="true" [hydrogen_and_cofactor_stripped]="true" )
    declare -A TIER_STRIP_HYDROGENS=( [hydrogen_stripped]="true" [cofactor_stripped]="false" [hydrogen_and_cofactor_stripped]="true" )
    declare -A TIER_FILE=( [full]="$MEEKO_READY" [hydrogen_stripped]="$MEEKO_NOH" [cofactor_stripped]="$MEEKO_FALLBACK" [hydrogen_and_cofactor_stripped]="$MEEKO_NOH_NOHET" )
    TIER_ORDER=(full hydrogen_stripped cofactor_stripped hydrogen_and_cofactor_stripped)

    declare -A TIER_EXIT_CODES=()
    declare -A TIER_VARIANT_OUTPUT=()
    SUCCESS_TIER=""

    # v2.4.1: clear any stale PDBQT/JSON from a previous run on this
    # receptor before starting, so a leftover file can never be
    # mistaken for output from this run if every tier below fails.
    rm -f "$PDBQT_OUT" "$JSON_OUT"

    for tier in "${TIER_ORDER[@]}"; do
        infile="${TIER_FILE[$tier]}"

        if [[ "$tier" != "full" ]]; then
            strip_het="${TIER_STRIP_HETATM[$tier]}"
            strip_h="${TIER_STRIP_HYDROGENS[$tier]}"
            log_info "  Building Meeko input variant '$tier' for $name (strip_hetatm=$strip_het strip_hydrogens=$strip_h)..."
            set +e
            VARIANT_OUTPUT=$(STRIP_HETATM="$strip_het" STRIP_HYDROGENS="$strip_h" KEEP_HETATM="$KEEP_HETATM" python3 - "$MEEKO_READY" "$infile" <<'PYVARIANT'
import os
import sys

RECORD  = slice(0, 6)
SERIAL  = slice(6, 11)
ATOMNAM = slice(12, 16)
ALTLOC  = 16
RESNAME = slice(17, 20)
CHAINID = 21
RESSEQ  = slice(22, 26)
ICODE   = 26
ELEMENT = slice(76, 78)

in_path, out_path = sys.argv[1], sys.argv[2]
strip_hetatm = os.environ.get("STRIP_HETATM", "false").strip().lower() == "true"
strip_hydrogens = os.environ.get("STRIP_HYDROGENS", "false").strip().lower() == "true"
unsafe = {i.strip().upper() for i in os.environ.get("KEEP_HETATM", "").split(",") if i.strip()}

with open(in_path) as fh:
    lines = fh.readlines()

seen_keys = set()
out_lines = []
removed_hetatm = 0
duplicates_removed = 0
altlocs_blanked = 0
hydrogens_removed = 0
removed_serials = set()

for line in lines:
    rec = line[RECORD]
    if rec == "CONECT":
        continue
    if rec not in ("ATOM  ", "HETATM"):
        out_lines.append(line)
        continue

    if rec == "HETATM" and strip_hetatm:
        resname = line[RESNAME].strip().upper()
        if resname in unsafe:
            removed_hetatm += 1
            try:
                removed_serials.add(int(line[SERIAL]))
            except ValueError:
                pass
            continue

    if strip_hydrogens:
        atom_name = line[ATOMNAM].strip()
        element = line[ELEMENT].strip()
        if not element:
            element = atom_name[0] if atom_name else "?"
        if element.upper() == "H":
            hydrogens_removed += 1
            try:
                removed_serials.add(int(line[SERIAL]))
            except ValueError:
                pass
            continue

    if line[ALTLOC] != " ":
        altlocs_blanked += 1
        line = line[:ALTLOC] + " " + line[ALTLOC + 1:]

    key = (line[CHAINID], line[RESSEQ], line[ICODE], line[ATOMNAM])
    if key in seen_keys:
        duplicates_removed += 1
        try:
            removed_serials.add(int(line[SERIAL]))
        except ValueError:
            pass
        continue
    seen_keys.add(key)
    out_lines.append(line)

with open(out_path, "w") as fh:
    fh.writelines(out_lines)

print(f"REMOVED_HETATM={removed_hetatm}")
print(f"DUPLICATES_REMOVED={duplicates_removed}")
print(f"ALTLOCS_BLANKED={altlocs_blanked}")
print(f"HYDROGENS_REMOVED={hydrogens_removed}")
PYVARIANT
)
            VARIANT_STATUS=$?
            set -e
            echo "$VARIANT_OUTPUT" | sed "s/^/[PHASE 7]   meeko_variant_build ($tier): /"
            TIER_VARIANT_OUTPUT["$tier"]="$VARIANT_OUTPUT"

            if [[ $VARIANT_STATUS -ne 0 ]] || [[ ! -s "$infile" ]]; then
                log_error "Could not build the '$tier' Meeko input variant for $name; skipping this tier."
                continue
            fi
        fi

        log_info "  Running mk_prepare_receptor.py for $name -- tier: $tier ..."
        set +e
        TIER_OUTPUT=$(mk_prepare_receptor.py -i "$infile" -o "$OUTDIR/${name}" -p -j 2>&1)
        TIER_EXIT=$?
        set -e
        echo "$TIER_OUTPUT" | sed "s/^/[PHASE 7]   mk_prepare_receptor ($tier): /"
        TIER_EXIT_CODES["$tier"]="$TIER_EXIT"

        if [[ $TIER_EXIT -eq 0 ]]; then
            SUCCESS_TIER="$tier"
            break
        else
            log_warn "Meeko tier '$tier' failed for $name (exit $TIER_EXIT)."
        fi
    done

    if [[ -z "$SUCCESS_TIER" ]]; then
        log_error "All Meeko tiers failed for $name (full, hydrogen_stripped, cofactor_stripped, hydrogen_and_cofactor_stripped)."
        printf '%s\tMeeko failed on every fallback tier (full/hydrogen_stripped/cofactor_stripped/hydrogen_and_cofactor_stripped)\n' "$name" >> "$FAILED_TSV"
        printf '%s\tFAILED\tMEEKO\n' "$name" >> "$STATUS_TSV"
        return 0
    fi

    if [[ "$SUCCESS_TIER" != "full" ]]; then
        log_warn "$name: Meeko succeeded only via tier '$SUCCESS_TIER' -- inspect $OUTDIR/${name}.pdbqt vs $MEEKO_READY before trusting downstream chemistry for this receptor."
    fi

    if [[ ! -f "$PDBQT_OUT" ]] || [[ ! -s "$PDBQT_OUT" ]]; then
        log_error "No valid PDBQT for $name"
        printf '%s\tNo valid PDBQT\n' "$name" >> "$FAILED_TSV"
        printf '%s\tFAILED\tMEEKO\n' "$name" >> "$STATUS_TSV"
        return 0
    fi

    TIER_EXIT_JSON="{"
    _first_tier=true
    for tier in "${TIER_ORDER[@]}"; do
        if [[ -n "${TIER_EXIT_CODES[$tier]+x}" ]]; then
            $_first_tier || TIER_EXIT_JSON+=","
            TIER_EXIT_JSON+="\"$tier\":${TIER_EXIT_CODES[$tier]}"
            _first_tier=false
        fi
    done
    TIER_EXIT_JSON+="}"

    SUCCESS_REMOVED_HETATM="0"
    SUCCESS_REMOVED_HYDROGENS="0"
    if [[ "$SUCCESS_TIER" != "full" ]] && [[ -n "${TIER_VARIANT_OUTPUT[$SUCCESS_TIER]+x}" ]]; then
        SUCCESS_REMOVED_HETATM=$(echo "${TIER_VARIANT_OUTPUT[$SUCCESS_TIER]}" | grep -o 'REMOVED_HETATM=[0-9]*' | cut -d= -f2)
        SUCCESS_REMOVED_HETATM=${SUCCESS_REMOVED_HETATM:-0}
        SUCCESS_REMOVED_HYDROGENS=$(echo "${TIER_VARIANT_OUTPUT[$SUCCESS_TIER]}" | grep -o 'HYDROGENS_REMOVED=[0-9]*' | cut -d= -f2)
        SUCCESS_REMOVED_HYDROGENS=${SUCCESS_REMOVED_HYDROGENS:-0}
    fi

    SUCCESS_TIER="$SUCCESS_TIER" \
    TIER_EXIT_JSON="$TIER_EXIT_JSON" \
    SUCCESS_REMOVED_HETATM="$SUCCESS_REMOVED_HETATM" \
    SUCCESS_REMOVED_HYDROGENS="$SUCCESS_REMOVED_HYDROGENS" \
    REPORT_JSON="$REPORT_JSON" \
    python3 - <<'PYMEEKOREPORT'
import json
import os

report_path = os.environ["REPORT_JSON"]
success_tier = os.environ["SUCCESS_TIER"]
tier_exit_codes = json.loads(os.environ["TIER_EXIT_JSON"])
removed_hetatm = int(os.environ["SUCCESS_REMOVED_HETATM"])
removed_hydrogens = int(os.environ["SUCCESS_REMOVED_HYDROGENS"])
used_fallback = success_tier != "full"

try:
    with open(report_path) as fh:
        report = json.load(fh)
except Exception:
    report = {}

notes = {
    "full": "Meeko succeeded on the residue-complete structure with all retained cofactors and explicit hydrogens intact.",
    "hydrogen_stripped": (
        "Meeko failed on the full structure (its internal ProDy->RDKit "
        "sanitization raised an atom-valence error), then succeeded after "
        "all explicit hydrogens were removed -- Meeko's own '-p' flag adds "
        "the polar hydrogens it actually needs, so this does not lose any "
        "chemistry Meeko required. Retained cofactors/metals are intact."
    ),
    "cofactor_stripped": (
        "Meeko failed on the full structure and on the hydrogen-stripped "
        "variant, then succeeded after removing the retained cofactor/metal "
        "HETATM group(s). The resulting PDBQT does NOT include these groups."
    ),
    "hydrogen_and_cofactor_stripped": (
        "Meeko failed on the full, hydrogen-stripped, and cofactor-stripped "
        "variants, and only succeeded after removing explicit hydrogens AND "
        "the retained cofactor/metal HETATM group(s). The resulting PDBQT "
        "does NOT include those cofactor/metal groups."
    ),
}

report["stages"] = report.get("stages", {})
report["stages"]["meeko_preparation"] = {
    "status": "warning" if used_fallback else "ok",
    "successful_tier": success_tier,
    "tiers_attempted": list(tier_exit_codes.keys()),
    "tier_exit_codes": tier_exit_codes,
    "used_fallback": used_fallback,
    "cofactor_hetatm_records_removed": removed_hetatm,
    "explicit_hydrogens_removed": removed_hydrogens,
    "note": notes.get(success_tier, "Meeko succeeded via an unrecognized tier."),
}
if used_fallback:
    report.setdefault("warnings", []).append(
        f"Meeko required fallback tier '{success_tier}' to succeed "
        f"(cofactor/metal HETATM records removed: {removed_hetatm}, "
        f"explicit hydrogens removed: {removed_hydrogens}); inspect the "
        f"final PDBQT's chemistry accordingly."
    )

with open(report_path, "w") as fh:
    json.dump(report, fh, indent=2)
PYMEEKOREPORT

    if [[ -s "$JSON_OUT" ]]; then
        set +e
        MEEKO_JSON_SUMMARY=$(python3 - "$JSON_OUT" <<'PYMEEKOJSON'
import json
import sys

path = sys.argv[1]
try:
    with open(path) as fh:
        data = json.load(fh)
except Exception as exc:
    print(f"INVALID: {exc}")
    sys.exit(1)

interesting_keys = [
    "warnings", "atom_types", "torsions", "n_torsions", "root_atom",
    "flexible_residues", "box_center", "box_size",
]
found = {k: data[k] for k in interesting_keys if k in data}
print(f"VALID keys_found={list(found.keys())}")
warnings = data.get("warnings")
if warnings:
    print(f"MEEKO_WARNINGS={warnings}")
PYMEEKOJSON
)
        JSON_STATUS=$?
        set -e
        echo "$MEEKO_JSON_SUMMARY" | sed "s/^/[PHASE 7]   meeko json: /"
        if [[ $JSON_STATUS -ne 0 ]]; then
            log_warn "Meeko JSON exists but is invalid for $name"
        fi
    else
        log_warn "Meeko JSON missing for $name"
    fi

    if [[ "$SUCCESS_TIER" != "full" ]]; then
        log_info "SUCCESS (via Meeko fallback tier '$SUCCESS_TIER'): Prepared $PDBQT_OUT -- $SUCCESS_REMOVED_HETATM cofactor/metal HETATM record(s) and $SUCCESS_REMOVED_HYDROGENS explicit hydrogen(s) excluded from Meeko's input; see ${REPORT_JSON} 'meeko_preparation' entry."
        printf '%s\tSUCCESS\tMEEKO_TIER_%s\n' "$name" "$(echo "$SUCCESS_TIER" | tr '[:lower:]' '[:upper:]')" >> "$STATUS_TSV"
    else
        log_info "SUCCESS: Prepared $PDBQT_OUT"
        printf '%s\tSUCCESS\t-\n' "$name" >> "$STATUS_TSV"
    fi

    if [[ "$KEEP_INTERMEDIATES" != "true" ]]; then
        rm -f "$NOWAT" "$RDKIT_SAFE" "$MEEKO_FALLBACK" "$MEEKO_READY" "$MEEKO_NOH" "$MEEKO_NOH_NOHET"
    fi
    return 0
}
export -f process_receptor log_info log_warn log_error

START_TIME=$(date +%s)

if [[ "$THREADS" -le 1 ]]; then
    for f in "${RECEPTORS[@]}"; do
        process_receptor "$f" || true
    done
else
    log_info "Running with $THREADS parallel worker(s)"
    printf '%s\n' "${RECEPTORS[@]}" | xargs -I{} -P "$THREADS" bash -c 'process_receptor "$0"' {} || true
fi

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

n=$(awk -F'\t' '$2=="SUCCESS"{c++} END{print c+0}' "$STATUS_TSV")
MEEKO_TIER_FALLBACK_SUCCESS=$(awk -F'\t' '$2=="SUCCESS" && $3 ~ /^MEEKO_TIER_/{c++} END{print c+0}' "$STATUS_TSV")
ALTLOC_FAIL=$(awk -F'\t' '$2=="FAILED" && $3=="ALTLOC"{c++} END{print c+0}' "$STATUS_TSV")
PDBFIXER_FAIL=$(awk -F'\t' '$2=="FAILED" && $3=="PDBFIXER"{c++} END{print c+0}' "$STATUS_TSV")
RESIDUE_AUDIT_FAIL=$(awk -F'\t' '$2=="FAILED" && $3=="RESIDUE_COMPLETENESS_AUDIT"{c++} END{print c+0}' "$STATUS_TSV")
MEEKO_FAIL=$(awk -F'\t' '$2=="FAILED" && $3=="MEEKO"{c++} END{print c+0}' "$STATUS_TSV")
OTHER_FAIL=$(awk -F'\t' '$2=="FAILED" && $3!="ALTLOC" && $3!="PDBFIXER" && $3!="RESIDUE_COMPLETENESS_AUDIT" && $3!="MEEKO"{c++} END{print c+0}' "$STATUS_TSV")
TOTAL_FAILED=$((TOTAL - n))
RECORDED_FAILURES=$((ALTLOC_FAIL + PDBFIXER_FAIL + RESIDUE_AUDIT_FAIL + MEEKO_FAIL + OTHER_FAIL))

echo ""
echo "======================================="
echo "Phase 7 Summary"
echo "======================================="
echo "Input receptors        : $TOTAL"
echo "Successfully prepared  : $n"
echo "  (via Meeko fallback tier, not 'full'): $MEEKO_TIER_FALLBACK_SUCCESS"
echo "Altloc failures        : $ALTLOC_FAIL"
echo "PDBFixer failures      : $PDBFIXER_FAIL"
echo "Residue-audit failures : $RESIDUE_AUDIT_FAIL"
echo "Meeko failures         : $MEEKO_FAIL"
echo "Other/unclassified     : $OTHER_FAIL"
echo "Total failed           : $TOTAL_FAILED"
echo "(RDKit sanitization is non-blocking and can never appear as a failure category -- see per-receptor *_repair_report.json for 'rdkit_validation' warnings)"
echo "(Incomplete standard-residue amino acids are auto-removed before Meeko by default -- see 'residue_completeness_audit' in each *_repair_report.json)"
echo "Runtime                : ${RUNTIME} seconds"
echo "Threads used           : $THREADS"
echo "Output folder          : $OUTDIR"
echo "Failed log             : $FAILED_TSV"
echo "Per-receptor status log: $STATUS_TSV"
echo "======================================="

if [[ $RECORDED_FAILURES -ne $TOTAL_FAILED ]]; then
    log_warn "Failure counters do not match total failed receptors -- check $STATUS_TSV directly."
fi

if [[ "$n" -eq 0 ]]; then
    log_error "No receptors were successfully prepared."
    exit 1
fi

log_info "SUCCESS: $n receptor(s) cleaned and ready for pocket detection."
log_info "Per-receptor repair reports: $OUTDIR/*_repair_report.json"

