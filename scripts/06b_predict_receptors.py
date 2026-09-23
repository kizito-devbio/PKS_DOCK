#!/usr/bin/env python3
"""
Phase 6b - Structure prediction for receptors with no experimental PDB.
Version 2.2

WHAT THIS PHASE DOES, IN ORDER, FOR EVERY RECEPTOR THAT NEEDS ONE:

    1. Local FASTA supplied?
       -> use its sequence as-is (never overwritten).
       -> if it has no stored UniProt accession, still search UniProt
          purely to try to recover an accession so AlphaFold DB can be
          queried. The candidate's sequence is compared against the
          local FASTA using a real global BLOSUM62 protein alignment
          (Bio.Align.PairwiseAligner), not text diffing, and only
          accepted as a probable match above an identity threshold.
    2. No local FASTA -> UniProt search (name + organism), up to N
       candidates per tier, ranked by: reviewed/Swiss-Prot status, exact
       gene-name match, protein-name/gene-name text similarity, and
       species- vs genus-level match. Full ranking rationale is logged
       and recorded in the decision log.
    3. AlphaFold DB lookup, tried across the primary accession AND the
       next-best alternate candidate accessions (configurable via
       --alphafold-max-candidates) before falling back to ColabFold.
       Each hit is (a) structurally validated and (b) checked against a
       minimum mean-pLDDT confidence threshold BEFORE being accepted; a
       hit that fails either check is discarded and the next candidate
       accession is tried.
    4. No usable AlphaFold DB entry for any candidate -> local ColabFold,
       with resume support, a configurable timeout, GPU/CPU
       auto-detection, optional --templates/--amber/--use-gpu-relax, and
       full stdout/stderr capture to a per-receptor log file.
    5. Every candidate structure is validated for existence, non-empty
       content, and a plausible number of ATOM/CA records. When the
       input sequence length is known, the CA-atom count is checked
       against it directly (>= 80% of expected residues) instead of a
       fixed magic number, so small real domains aren't penalised and
       larger truncated models are still caught.
    6. SHA-256 checksums are recorded for the FASTA and the final PDB of
       every receptor, for reproducibility / future integrity checks.
    7. The manifest's existing receptor records (already in
       "to_predict") are updated IN PLACE with the prediction outcome --
       nothing Phase 6 put there is discarded.
    8. A detailed summary (AlphaFold reused, ColabFold fresh/reused,
       validation failures, low-confidence rejections, elapsed time) is
       printed at the end.

Every automated choice is still written to the shared decision log
(pks_dock.reproducibility.DecisionLog) so it ends up in workflow_report.md.

v2.2 changes vs v2.1 (all additive / optional-flag-gated, so existing
invocations from run_pipeline.sh keep working unmodified):
    - Real protein sequence identity (Bio.Align.PairwiseAligner, BLOSUM62
      global alignment) replaces difflib for anything that is an actual
      amino-acid sequence comparison. difflib is kept only for short
      protein-name/gene-name text comparisons, where it's appropriate.
    - PDB validation now uses an expected-residue-count-relative
      threshold (>=80% of input sequence length) when that length is
      known, instead of a single hard-coded minimum for every receptor.
    - AlphaFold DB hits below --min-mean-plddt (default 50.0) are now
      rejected and the next candidate accession is tried, instead of
      being silently accepted regardless of confidence.
    - New optional ColabFold flags: --colabfold-num-models,
      --colabfold-templates, --colabfold-amber.
    - SHA-256 checksums recorded for FASTA and final PDB files.
"""
import argparse
import difflib
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Allow running this script directly (python3 scripts/06b_...py) as well as
# via the installed `pks-dock` package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pks_dock.net import HTTPClient  # noqa: E402
from pks_dock.alphafold import resolve_receptor_structure  # noqa: E402
from pks_dock.reproducibility import DecisionLog  # noqa: E402

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
PHASE = "PHASE 6b"


def log(msg):
    print(f"[{PHASE}] {msg}", flush=True)


# --------------------------------------------------------------------------
# Biological sequence identity (replaces difflib for real AA sequences)
# --------------------------------------------------------------------------

try:
    from Bio.Align import PairwiseAligner, substitution_matrices
    _aligner = PairwiseAligner()
    _aligner.mode = "global"
    _aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    _aligner.open_gap_score = -10
    _aligner.extend_gap_score = -0.5
    _BIOPYTHON_AVAILABLE = True
except Exception as _bio_import_err:  # pragma: no cover - environment dependent
    _BIOPYTHON_AVAILABLE = False
    _BIO_IMPORT_ERROR = _bio_import_err


def protein_sequence_identity(seq_a: str, seq_b: str) -> float:
    """
    Percent identity from a global BLOSUM62 alignment
    (Bio.Align.PairwiseAligner), used for any comparison between two
    actual amino-acid sequences. This replaces difflib.SequenceMatcher
    for that purpose: difflib does plain text diffing and has no notion
    of conservative substitutions, indels, or alignment.

    Reconstructs matches from `alignment.aligned` (the documented,
    coordinate-based block representation Bio.Align returns: a pair of
    arrays of (start, end) index ranges into the two original
    sequences) and normalises by `alignment.shape[1]` (the total number
    of alignment columns, gaps included), rather than by stringifying
    and indexing the Alignment object. The coordinate form doesn't
    depend on how a particular Biopython version renders or indexes an
    alignment as text, so it stays correct across Biopython releases
    without needing to be re-verified against each one.

    Falls back to difflib (logged once) only if Biopython's Align
    module isn't importable in this environment.
    """
    if not seq_a or not seq_b:
        return 0.0
    if not _BIOPYTHON_AVAILABLE:
        log(f"NOTE: Bio.Align not available ({_BIO_IMPORT_ERROR}); falling back to difflib "
            f"text-similarity for sequence comparison, which is weaker for biological sequences.")
        return difflib.SequenceMatcher(None, seq_a, seq_b).ratio()
    try:
        alignment = _aligner.align(seq_a, seq_b)[0]
        blocks_a, blocks_b = alignment.aligned  # documented coordinate-block API
        matches = 0
        for (a_start, a_end), (b_start, b_end) in zip(blocks_a, blocks_b):
            sub_a = seq_a[a_start:a_end]
            sub_b = seq_b[b_start:b_end]
            matches += sum(1 for x, y in zip(sub_a, sub_b) if x == y)
        total_columns = alignment.shape[1]  # alignment length, gaps included
        return matches / total_columns if total_columns else 0.0
    except Exception as e:
        log(f"Biopython alignment failed ({e}), falling back to difflib for this comparison.")
        return difflib.SequenceMatcher(None, seq_a, seq_b).ratio()


def _text_similarity(a: str, b: str) -> float:
    """Plain text similarity for short strings like protein/gene names
    (NOT for amino-acid sequences -- use protein_sequence_identity for
    those)."""
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a.lower().strip(), b.lower().strip()).ratio()


# --------------------------------------------------------------------------
# Receptor name / organism resolution
# --------------------------------------------------------------------------

def compute_species_tag(organism_name: str) -> str:
    """
    'Staphylococcus aureus' -> 'Saureus'
    'Bacillus cereus'       -> 'Bcereus'
    'Pseudomonas aeruginosa'-> 'Paeruginosa'

    Mirrors the tagging convention already visible in your Phase 6 logs
    (SrtA_Saureus, GyrA_Bcereus, LasB_Paeruginosa, ...).
    """
    parts = organism_name.replace("_", " ").split()
    if len(parts) >= 2:
        return parts[0][0].upper() + parts[1].lower()
    return organism_name.replace(" ", "").replace("_", "")


def build_tag_lookup(organism_lookup: dict) -> dict:
    """
    Maps every tag derivable from organism_lookup values (and its keys,
    in case a key is already a short tag) back to the full organism
    name. Longest tags are preferred on lookup so a short tag can't
    accidentally match inside a longer one.
    """
    tag_lookup = {}
    for key, organism_name in organism_lookup.items():
        tag_lookup[compute_species_tag(organism_name)] = organism_name
        tag_lookup[compute_species_tag(key)] = organism_name
        tag_lookup[key] = organism_name
    return tag_lookup


def split_receptor_name(name: str, tag_lookup: dict):
    """
    'SrtA_Saureus' -> ('SrtA', 'Saureus'). Tries known tags (longest
    first) as a suffix match before falling back to a naive rsplit on
    the last underscore.
    """
    for tag in sorted(tag_lookup.keys(), key=len, reverse=True):
        suffix = "_" + tag
        if name.endswith(suffix):
            return name[: -len(suffix)], tag
    if "_" in name:
        target, tag = name.rsplit("_", 1)
        return target, tag
    return name, ""


# --------------------------------------------------------------------------
# UniProt candidate search + ranking
# --------------------------------------------------------------------------

def _uniprot_json_search(client: HTTPClient, query: str, size: int):
    params = {
        "query": query,
        "format": "json",
        "size": size,
        "fields": "accession,id,protein_name,organism_name,reviewed,sequence,gene_names",
    }
    text = client.get_text(UNIPROT_SEARCH, params=params)
    data = json.loads(text)
    return data.get("results", [])


def _parse_candidate(entry: dict, match_level: str, target_short_name: str):
    accession = entry.get("primaryAccession")
    entry_name = entry.get("uniProtkbId")
    organism_name = entry.get("organism", {}).get("scientificName")
    seq = entry.get("sequence", {}).get("value")

    protein_desc = entry.get("proteinDescription", {}) or {}
    rec_name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value")
    if not rec_name:
        submitted = protein_desc.get("submissionNames") or []
        if submitted:
            rec_name = submitted[0].get("fullName", {}).get("value")
    protein_name = rec_name or "unknown protein name"

    gene_names = []
    for g in entry.get("genes", []) or []:
        gene_name = (g.get("geneName") or {}).get("value")
        if gene_name:
            gene_names.append(gene_name)

    entry_type = entry.get("entryType", "") or ""
    reviewed = entry_type.startswith("UniProtKB reviewed")

    gene_match = any(g.lower() == target_short_name.lower() for g in gene_names)
    # Text similarity here is comparing SHORT NAMES/LABELS (e.g. "SrtA" vs
    # "Sortase A" vs a gene symbol) -- difflib is appropriate for this; it
    # is protein_sequence_identity() that's used for actual AA sequences.
    name_similarity = max(
        [_text_similarity(target_short_name, protein_name)]
        + [_text_similarity(target_short_name, g) for g in gene_names]
        + [0.0]
    )

    score = 0.0
    if reviewed:
        score += 10
    if gene_match:
        score += 5
    score += name_similarity * 5
    if match_level == "species-level":
        score += 3
    if seq:
        score += 1

    return {
        "accession": accession,
        "entry_name": entry_name,
        "protein_name": protein_name,
        "organism_name": organism_name,
        "sequence": seq,
        "gene_names": gene_names,
        "reviewed": reviewed,
        "match_level": match_level,
        "name_similarity": round(name_similarity, 3),
        "score": round(score, 3),
    }


def find_uniprot_candidates(client: HTTPClient, target_short_name: str, organism: str, size: int):
    """
    Returns a ranked list of candidate dicts (best first). Searches
    species-level then genus-level, merges, de-duplicates by accession
    (keeping the higher-scoring instance), and ranks together.
    """
    candidates = []

    species_query = f'({target_short_name}) AND (organism_name:"{organism}")'
    log(f"Querying UniProt (species-level): {species_query}")
    try:
        for e in _uniprot_json_search(client, species_query, size):
            candidates.append(_parse_candidate(e, "species-level", target_short_name))
    except Exception as e:
        log(f"  UniProt species-level query failed: {e}")

    genus = organism.split()[0] if organism else ""
    if genus:
        genus_query = f"({target_short_name}) AND (organism_name:{genus})"
        log(f"  Also querying UniProt (genus-level, for comparison): {genus_query}")
        try:
            for e in _uniprot_json_search(client, genus_query, size):
                candidates.append(_parse_candidate(e, "genus-level", target_short_name))
        except Exception as e:
            log(f"  UniProt genus-level query failed: {e}")

    by_accession = {}
    for c in candidates:
        if not c["accession"] or not c["sequence"]:
            continue
        existing = by_accession.get(c["accession"])
        if existing is None or c["score"] > existing["score"]:
            by_accession[c["accession"]] = c

    return sorted(by_accession.values(), key=lambda c: c["score"], reverse=True)


def candidate_to_fasta(c: dict) -> str:
    header = f">{c['accession']}|{c['entry_name']}|{c['protein_name']} OS={c['organism_name']}"
    seq = c["sequence"]
    wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
    return f"{header}\n{wrapped}\n"


def read_fasta_sequence(fasta_path: Path) -> str:
    """Plain, dependency-free FASTA sequence reader (concatenates all
    non-header lines, strips whitespace). Good enough to get the
    sequence out for identity/length checks; not a general parser."""
    seq_lines = []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            continue
        seq_lines.append(line.strip())
    return "".join(seq_lines)


# --------------------------------------------------------------------------
# Checksums
# --------------------------------------------------------------------------

def sha256_of_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# PDB validation
# --------------------------------------------------------------------------

def validate_pdb(path: Path, expected_residues: int = None):
    """
    Dependency-free sanity check: file exists, is non-empty, contains a
    plausible number of ATOM/HETATM records AND a plausible number of
    alpha-carbon (CA) atoms.

    If expected_residues (the length of the input amino-acid sequence)
    is known, the CA count is required to be >= 80% of it -- this
    scales correctly for both small real domains (won't be unfairly
    rejected by a flat minimum) and large truncated/corrupted models
    (still caught, proportionally). If it isn't known, a conservative
    fixed floor is used instead.
    """
    if not path.exists():
        return False, "file does not exist"
    if path.stat().st_size == 0:
        return False, "file is empty"

    atom_lines = 0
    ca_lines = 0
    residues = set()
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_lines += 1
                residues.add(line[21:27])
                atom_name = line[12:16].strip()
                if line.startswith("ATOM") and atom_name == "CA":
                    ca_lines += 1

    if atom_lines < 50:
        return False, f"only {atom_lines} atom record(s) -- too few for a real receptor structure"

    if expected_residues:
        min_ca = max(5, int(0.8 * expected_residues))
        if ca_lines < min_ca:
            return False, (f"only {ca_lines} CA atom(s), expected at least {min_ca} "
                            f"(80% of the {expected_residues}-residue input sequence)")
    else:
        if len(residues) < 10:
            return False, f"only {len(residues)} residue(s) found -- suspiciously small for a receptor"
        if ca_lines < 10:
            return False, f"only {ca_lines} alpha-carbon (CA) atom(s) found -- doesn't look like a real protein backbone"

    return True, f"{atom_lines} atom records, {len(residues)} residues, {ca_lines} CA atoms"


# --------------------------------------------------------------------------
# Local ColabFold execution
# --------------------------------------------------------------------------

def find_existing_rank001(pred_outdir: Path):
    hits = sorted(pred_outdir.glob("*_relaxed_rank_001*.pdb")) or \
        sorted(pred_outdir.glob("*_unrelaxed_rank_001*.pdb"))
    return hits[0] if hits else None


def run_colabfold(name: str, fasta_path: Path, pred_outdir: Path, args):
    """Returns (model_path, error_message, reused_existing: bool)."""
    pred_outdir.mkdir(parents=True, exist_ok=True)
    colabfold_log = pred_outdir / "colabfold.log"

    if not args.force_colabfold_rerun:
        existing = find_existing_rank001(pred_outdir)
        if existing is not None:
            log(f"{name}: found existing ColabFold rank_001 model at {existing}, reusing it "
                f"(pass --force-colabfold-rerun to redo the prediction).")
            return existing, None, True

    if shutil.which("colabfold_batch") is None:
        return None, ("'colabfold_batch' is not installed/on PATH. Install ColabFold locally "
                       "(see environment.yml) or ensure this receptor has an AlphaFold DB entry."), False

    cmd = ["colabfold_batch", str(fasta_path), str(pred_outdir),
           "--num-recycle", str(args.colabfold_num_recycle),
           "--num-models", str(args.colabfold_num_models)]

    has_gpu = shutil.which("nvidia-smi") is not None
    if not has_gpu:
        log(f"{name}: no NVIDIA GPU detected on this machine, running ColabFold on CPU "
            f"(this will be significantly slower).")
        cmd.append("--cpu")

    if args.colabfold_templates:
        cmd.append("--templates")
    if args.colabfold_amber:
        cmd.append("--amber")
        if has_gpu:
            cmd.append("--use-gpu-relax")

    log(f"{name}: running local ColabFold: {' '.join(cmd)}")
    log(f"{name}: full stdout/stderr will be written to {colabfold_log}")

    start = time.time()
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=args.colabfold_timeout if args.colabfold_timeout else None,
        )
    except subprocess.TimeoutExpired as e:
        colabfold_log.write_text(
            f"TIMEOUT after {args.colabfold_timeout}s\n\nSTDOUT:\n{e.stdout or ''}\n\nSTDERR:\n{e.stderr or ''}\n"
        )
        return None, f"ColabFold timed out after {args.colabfold_timeout}s (see {colabfold_log})", False

    elapsed = time.time() - start
    colabfold_log.write_text(f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n")

    if result.returncode != 0:
        tail = "\n".join(result.stderr.strip().splitlines()[-15:])
        return None, (f"ColabFold exited with code {result.returncode} after {elapsed:.0f}s. "
                       f"Last lines of stderr:\n{tail}\nFull log: {colabfold_log}"), False

    ranked = find_existing_rank001(pred_outdir)
    if ranked is None:
        return None, (f"ColabFold finished ({elapsed:.0f}s) but no rank_001 model file was found "
                       f"in {pred_outdir} -- inspect {colabfold_log} manually."), False

    log(f"{name}: ColabFold finished in {elapsed:.0f}s -> {ranked}")
    return ranked, None, False


# --------------------------------------------------------------------------
# Manifest persistence
# --------------------------------------------------------------------------

def write_manifest_atomic(manifest: dict, manifest_path: Path):
    tmp_path = manifest_path.with_suffix(manifest_path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(manifest, indent=2))
    tmp_path.replace(manifest_path)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="results/receptors_raw/receptor_manifest.json")
    ap.add_argument("--organism-lookup", required=True,
                     help='JSON mapping pathogen_key -> organism name, e.g. \'{"Bacillus_cereus": "Bacillus cereus"}\'')
    ap.add_argument("--fasta-dir", default="config/sequences")
    ap.add_argument("--outdir", default="results/receptors_raw")
    ap.add_argument("--decision-log", default="results/reports/decision_log.json")
    ap.add_argument("--no-alphafold-db", action="store_true",
                     help="Skip the AlphaFold DB lookup and go straight to local ColabFold (for offline/testing use).")
    ap.add_argument("--uniprot-candidates", type=int, default=20,
                     help="How many UniProt hits to retrieve per query tier before ranking (default 20).")
    ap.add_argument("--colabfold-timeout", type=int, default=0,
                     help="Kill a local ColabFold run after this many seconds. 0 = no timeout (default).")
    ap.add_argument("--colabfold-num-recycle", type=int, default=3,
                     help="ColabFold --num-recycle value (default 3).")
    ap.add_argument("--colabfold-num-models", type=int, default=5,
                     help="ColabFold --num-models value (default 5, matching ColabFold's own default).")
    ap.add_argument("--colabfold-templates", action="store_true",
                     help="Pass --templates to colabfold_batch (use PDB templates during prediction).")
    ap.add_argument("--colabfold-amber", action="store_true",
                     help="Pass --amber to colabfold_batch for Amber relaxation (adds --use-gpu-relax "
                          "automatically if a GPU is detected).")
    ap.add_argument("--force-colabfold-rerun", action="store_true",
                     help="Re-run ColabFold even if a rank_001 model already exists from a previous run.")
    ap.add_argument("--alphafold-max-candidates", type=int, default=5,
                     help="Max number of alternate UniProt accessions to try against AlphaFold DB (default 5).")
    ap.add_argument("--fasta-accession-similarity", type=float, default=0.90,
                     help="Minimum protein sequence identity (global BLOSUM62 alignment) required before an "
                          "accession recovered via UniProt search is trusted for a receptor that already had a "
                          "local FASTA (default 0.90).")
    ap.add_argument("--min-mean-plddt", type=float, default=50.0,
                     help="Minimum acceptable AlphaFold DB mean pLDDT confidence. Hits below this are rejected "
                          "and the next candidate accession (or ColabFold) is tried instead. Set to 0 to disable "
                          "this check (default 50.0).")
    args = ap.parse_args()

    run_start = time.time()

    manifest_path = Path(args.manifest)
    manifest = json.loads(manifest_path.read_text())
    organism_lookup = json.loads(args.organism_lookup)
    tag_lookup = build_tag_lookup(organism_lookup)

    fasta_dir = Path(args.fasta_dir)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    client = HTTPClient(phase=PHASE)
    decisions = DecisionLog(args.decision_log)

    to_predict = manifest.get("to_predict", [])
    log(f"{len(to_predict)} receptor(s) require structure prediction.")
    if not _BIOPYTHON_AVAILABLE:
        log(f"WARNING: Biopython's Bio.Align is not importable ({_BIO_IMPORT_ERROR}). "
            f"All sequence-identity checks in this run will fall back to difflib text similarity, "
            f"which is weaker for biological sequences. Check your environment (Bio 1.81 is listed "
            f"in environment.yml).")

    counts = {
        "total": len(to_predict),
        "alphafold_db": 0,
        "colabfold_fresh": 0,
        "colabfold_reused": 0,
        "failed": 0,
        "validation_failures": 0,
        "low_confidence_rejections": 0,
    }

    for entry in to_predict:
        name = entry["name"]
        target_short_name, pathogen_tag = split_receptor_name(name, tag_lookup)
        organism = tag_lookup.get(pathogen_tag) or organism_lookup.get(pathogen_tag) or pathogen_tag

        fasta_path = Path(entry["fasta_path"]) if entry.get("fasta_path") else fasta_dir / f"{name}.fasta"
        accession = entry.get("uniprot_accession")
        alternate_accessions = []
        # Sequence text is read from disk at most once per receptor per
        # branch below and reused everywhere it's needed (identity check,
        # expected_len for validation) rather than re-reading the file
        # each time it's needed.
        local_seq = None

        if fasta_path.exists():
            log(f"{name}: using existing supplied FASTA at {fasta_path}")
            local_seq = read_fasta_sequence(fasta_path)

            if not accession:
                log(f"{name}: existing FASTA has no stored UniProt accession -- searching UniProt "
                    f"anyway, purely to try to recover one for an AlphaFold DB lookup. The supplied "
                    f"FASTA sequence itself will NOT be modified or replaced.")
                candidates = find_uniprot_candidates(client, target_short_name, organism, args.uniprot_candidates)
                accepted = None
                for c in candidates:
                    identity = protein_sequence_identity(local_seq, c["sequence"])
                    if identity >= args.fasta_accession_similarity:
                        accepted = c
                        accepted["sequence_identity_to_local_fasta"] = round(identity, 3)
                        break
                    else:
                        log(f"{name}: candidate {c['accession']} rejected -- global-alignment sequence "
                            f"identity to local FASTA is only {identity:.2f} "
                            f"(need >= {args.fasta_accession_similarity}).")
                if accepted:
                    accession = accepted["accession"]
                    remaining = [c for c in candidates if c is not accepted][:args.alphafold_max_candidates]
                    alternate_accessions = [c["accession"] for c in remaining]
                    log(f"{name}: recovered accession {accession} ({accepted['entry_name']}) for AlphaFold DB "
                        f"lookup, sequence identity {accepted['sequence_identity_to_local_fasta']:.2f} "
                        f"to the locally supplied FASTA (global BLOSUM62 alignment).")
                    decisions.record(
                        PHASE, f"Accession recovery for existing FASTA of {name}",
                        database="UniProt", query=f"{target_short_name} / {organism}",
                        selected=accession,
                        reason=(f"Local FASTA had no stored accession; recovered {accession} via UniProt "
                                f"search with global-alignment sequence identity "
                                f"{accepted['sequence_identity_to_local_fasta']:.2f}."),
                        confidence="medium",
                        alternatives_considered=alternate_accessions,
                    )
                else:
                    log(f"{name}: no UniProt candidate matched the local FASTA closely enough "
                        f"(threshold {args.fasta_accession_similarity}); proceeding to ColabFold "
                        f"(no accession available for an AlphaFold DB lookup).")
                    decisions.record(
                        PHASE, f"Accession recovery for existing FASTA of {name}",
                        database="UniProt", query=f"{target_short_name} / {organism}", selected=None,
                        reason=f"No UniProt candidate reached identity threshold {args.fasta_accession_similarity}.",
                        confidence="none",
                    )
        else:
            log(f"{name}: no local FASTA found, auto-fetching sequence from UniProt "
                f"(target='{target_short_name}', organism='{organism}')...")
            candidates = find_uniprot_candidates(client, target_short_name, organism, args.uniprot_candidates)

            if not candidates:
                log(f"{name}: COULD NOT auto-fetch a sequence from UniProt (no candidates at species "
                    f"or genus level). Supply one manually at {fasta_path} before this receptor can be predicted.")
                decisions.record(
                    PHASE, f"Sequence acquisition for {name}",
                    database="UniProt", query=f"{target_short_name} / {organism}",
                    selected=None, reason="No UniProt hit at species or genus level.",
                    confidence="none",
                )
                counts["failed"] += 1
                entry["prediction_status"] = "failed"
                entry["prediction_reason"] = "no UniProt sequence found"
                continue

            best = candidates[0]
            accession = best["accession"]
            alternate_accessions = [c["accession"] for c in candidates[1:1 + args.alphafold_max_candidates]]
            local_seq = best["sequence"]  # already have it in memory, no need to re-read from disk

            fasta_path.write_text(candidate_to_fasta(best))
            note = (f"selected {accession} ({best['entry_name']}, "
                    f"{'reviewed' if best['reviewed'] else 'unreviewed'}, {best['match_level']}, "
                    f"name similarity {best['name_similarity']}) out of {len(candidates)} candidate(s)")
            log(f"{name}: {note} -> {fasta_path}")
            decisions.record(
                PHASE, f"Sequence acquisition for {name}",
                database="UniProt", query=f"{target_short_name} / {organism}",
                selected=accession, reason=note,
                confidence="high" if best["reviewed"] and best["match_level"] == "species-level" else "medium",
                alternatives_considered=alternate_accessions,
            )

        expected_len = len(local_seq) if local_seq else None
        # Hashing genuinely needs to read the file's raw bytes (a different
        # concern from parsing out the sequence text), so this is a
        # necessary second file access, not a redundant one.
        fasta_sha256 = sha256_of_file(fasta_path) if fasta_path.exists() else None

        # --- AlphaFold DB check(s): validated + confidence-gated -------------
        final_path = outdir / f"{name}.pdb"
        source = None
        af_hit = None

        accessions_to_try = [a for a in ([accession] + alternate_accessions) if a]
        if accessions_to_try and not args.no_alphafold_db:
            for acc in accessions_to_try:
                log(f"{name}: checking AlphaFold DB for accession {acc}...")
                af_path, hit, note = resolve_receptor_structure(client, acc, outdir / "alphafold_db")
                if af_path is None:
                    log(f"{name}: {note}")
                    continue

                ok, detail = validate_pdb(af_path, expected_residues=expected_len)
                if not ok:
                    log(f"{name}: AlphaFold DB hit for {acc} FAILED structural validation ({detail}); "
                        f"trying next candidate accession instead of accepting it.")
                    counts["validation_failures"] += 1
                    decisions.record(
                        PHASE, f"Structure source for {name}", database="AlphaFold DB", query=acc,
                        selected=None, reason=f"AlphaFold DB hit failed validation: {detail}", confidence="n/a",
                    )
                    continue

                if hit and hit.mean_plddt is not None and args.min_mean_plddt > 0 and hit.mean_plddt < args.min_mean_plddt:
                    log(f"{name}: AlphaFold DB hit for {acc} has mean pLDDT {hit.mean_plddt:.1f}, below the "
                        f"--min-mean-plddt threshold of {args.min_mean_plddt} -- treating as low-confidence and "
                        f"trying the next candidate instead (use --min-mean-plddt 0 to disable this check).")
                    counts["low_confidence_rejections"] += 1
                    decisions.record(
                        PHASE, f"Structure source for {name}", database="AlphaFold DB", query=acc,
                        selected=None,
                        reason=f"Rejected: mean pLDDT {hit.mean_plddt:.1f} below threshold {args.min_mean_plddt}.",
                        confidence=f"mean pLDDT {hit.mean_plddt:.1f} (rejected)",
                    )
                    continue

                shutil.copy2(af_path, final_path)
                log(f"{name}: {note} -- validated OK ({detail})")
                decisions.record(
                    PHASE, f"Structure source for {name}",
                    database="AlphaFold DB", query=acc,
                    selected=hit.entry_id if hit else acc,
                    alternatives_considered=accessions_to_try,
                    reason=note,
                    confidence=f"mean pLDDT {hit.mean_plddt:.1f}" if hit and hit.mean_plddt else None,
                )
                source = "AlphaFold DB"
                af_hit = hit
                break
        elif not accessions_to_try:
            log(f"{name}: no UniProt accession available (stored or recovered), cannot query AlphaFold DB; "
                f"proceeding straight to local ColabFold.")

        if source is None:
            if accessions_to_try and not args.no_alphafold_db:
                decisions.record(
                    PHASE, f"Structure source for {name}",
                    database="AlphaFold DB", query=accessions_to_try[0], selected=None,
                    alternatives_considered=["local ColabFold prediction"],
                    reason=f"No valid/confident AlphaFold DB entry for any of {len(accessions_to_try)} "
                           f"candidate accession(s) tried.",
                    confidence="n/a",
                )

            pred_outdir = outdir / f"{name}_colabfold"
            model_path, err, reused = run_colabfold(name, fasta_path, pred_outdir, args)
            if model_path is None:
                log(f"{name}: {err}")
                decisions.record(PHASE, f"Structure source for {name}", selected=None,
                                  reason=err, confidence="n/a")
                counts["failed"] += 1
                entry["prediction_status"] = "failed"
                entry["prediction_reason"] = err
                continue

            ok, detail = validate_pdb(model_path, expected_residues=expected_len)
            if not ok:
                log(f"{name}: ColabFold model FAILED validation ({detail}).")
                counts["validation_failures"] += 1
                counts["failed"] += 1
                decisions.record(PHASE, f"Structure source for {name}", selected=None,
                                  reason=f"ColabFold model failed validation: {detail}", confidence="n/a")
                entry["prediction_status"] = "failed"
                entry["prediction_reason"] = f"ColabFold model failed validation: {detail}"
                continue

            shutil.copy2(model_path, final_path)
            log(f"{name}: best local ColabFold model copied -> {final_path} (validated OK: {detail})")
            decisions.record(
                PHASE, f"Structure source for {name}", database="local ColabFold",
                selected=str(model_path.name),
                alternatives_considered=["AlphaFold DB (no valid/confident entry found for any candidate accession)"],
                reason="No usable AlphaFold DB prediction for any UniProt candidate; ran local ColabFold instead.",
            )
            source = "local ColabFold"
            counts["colabfold_reused" if reused else "colabfold_fresh"] += 1

        if source == "AlphaFold DB":
            counts["alphafold_db"] += 1

        entry["prediction_status"] = "ok"
        entry["prediction_source"] = source
        entry["prediction_pdb"] = str(final_path)
        entry["prediction_uniprot_accession"] = accession
        entry["prediction_mean_plddt"] = af_hit.mean_plddt if (af_hit and af_hit.mean_plddt) else None
        entry["fasta_sha256"] = fasta_sha256
        entry["prediction_pdb_sha256"] = sha256_of_file(final_path)

    # --- Persist outcomes back into the manifest, in place ------------------
    write_manifest_atomic(manifest, manifest_path)
    log(f"Updated manifest written -> {manifest_path} (existing receptor records updated in place)")

    elapsed_total = time.time() - run_start

    log("Summary:")
    log(f"  Total receptors needing prediction : {counts['total']}")
    log(f"  AlphaFold DB reused                 : {counts['alphafold_db']}")
    log(f"  ColabFold predicted (fresh)          : {counts['colabfold_fresh']}")
    log(f"  ColabFold predicted (reused run)     : {counts['colabfold_reused']}")
    log(f"  Validation failures encountered      : {counts['validation_failures']}")
    log(f"  Low-confidence AlphaFold rejections  : {counts['low_confidence_rejections']}")
    log(f"  Failed                               : {counts['failed']}")
    log(f"  Total elapsed time                   : {elapsed_total:.0f}s")

    log("Phase 6b complete.")


if __name__ == "__main__":
    main()
    