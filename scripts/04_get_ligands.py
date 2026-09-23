#!/usr/bin/env python3
"""
Phase 4 - Fully automated compound structure retrieval.

Reads results/antismash/phase3_compound_list.json produced by Phase 3.
Nothing is hand-typed here -- every compound name antiSMASH predicted is
looked up automatically against a chain of public, key-free, documented
chemical databases:

  1. PubChem PUG-REST (direct HTTP, not pubchempy -- this avoids relying on
     an unverified pubchempy function; PUG-REST's URL scheme is documented
     and stable: https://pubchem.ncbi.nlm.nih.gov/rest/pug/...)
  2. NCI Chemical Identifier Resolver (CIR) -- a real, simple REST resolver
     for common/trivial chemical names, useful for natural products that
     don't have a PubChem name-indexed entry. Only used to resolve a name
     to a SMILES string (see "3D generation" below for why CIR's own SDF
     endpoint is no longer used directly).
  3. ChEBI -- EBI's *current* ChEBI 2.0 REST API
     (https://www.ebi.ac.uk/chebi/backend/api/public/), searched via
     /advanced_search/ and resolved via /compound/{id}/.
  4. ChEMBL -- EBI's ChEMBL web-services REST API (also key-free and
     documented: https://www.ebi.ac.uk/chembl/api/data/docs), searched by
     name and, on a hit, resolved to a canonical SMILES / MOL block from
     the full compound record.

If a compound is not found anywhere, it is reported by name and skipped --
never silently guessed or fabricated.

===========================================================================
REVISION HISTORY / REVIEW RESPONSE
===========================================================================
This revision responds point-by-point to a pre-publication review of the
previous version. Each point below is a real change in this file, not just
a comment.

1. ChEMBL endpoints re-verified. As of this revision, the base
   `https://www.ebi.ac.uk/chembl/api/data/` service, the
   `molecule/search.json?q=...` search endpoint, and the
   `molecule/{chembl_id}.json` record endpoint are all still the current,
   documented ChEMBL web-services API (see
   https://www.ebi.ac.uk/chembl/api/data/docs and the ChEMBL interface
   documentation). CHEMBL_SEARCH / CHEMBL_RECORD below are unchanged from
   the prior revision because they check out -- but treat this as "verified
   as of the date of this revision," not "verified forever": ChEMBL has
   changed its web-service layer more than once historically, so re-check
   these against https://www.ebi.ac.uk/chembl/api/data/docs (or a quick
   `curl`) before relying on this in a new environment or after a long gap.

2. ChEBI endpoint was actually wrong / stale, and has been fixed. The
   previous revision used the legacy `saveStructure.do` servlet
   (`https://www.ebi.ac.uk/chebi/saveStructure.do`), which is part of the
   old, pre-2025 ChEBI web application. ChEBI underwent a ground-up
   redevelopment ("ChEBI 2.0", relaunched ~October 2025): the legacy SOAP
   services were retired and the old JSP-era `.do` endpoints are gone with
   the old site. This file now uses the current ChEBI 2.0 REST API
   (`https://www.ebi.ac.uk/chebi/backend/api/public/`): `/advanced_search/`
   (POST) to find candidate ChEBI IDs for a name, and `/compound/{id}/`
   (GET) to fetch the SMILES and mass for a candidate. Because that JSON
   endpoint does not itself expose a MOL/SDF block, ChEBI candidates are
   resolved to SMILES and then handed to the shared RDKit 3D-generation
   step below, same as CIR and ChEMBL. As with point 1, re-verify this
   against https://www.ebi.ac.uk/chebi/backend/api/docs/ if it has been a
   while since this was last run against the live service.

3. Candidate-scoring line-count heuristic removed. The previous scoring
   function added points for a structure file simply having more lines,
   which is not a chemically meaningful signal (a longer file is not a
   more correct structure). That term is gone. The scoring function now
   only rewards (a) a genuine, source-provided 3D structure over an
   RDKit-generated one, (b) stereochemistry-aware SMILES, and (c) presence
   of basic required fields -- see score_candidate_structure().

4. Stereo-prefix stripping is documented as a broadening fallback only
   (unchanged behavior, strengthened bookkeeping). generate_name_variants()
   still tries the exact name first, and only adds the stereo-stripped
   variant as an additional, broader search string -- it never replaces
   the original. What's new: which variant actually produced a hit is now
   recorded end-to-end (search_variant_used field in the resolved-compound
   record, DecisionLog entries, and the new per-compound TSV report), so a
   reviewer can audit every case where the broadened variant -- and not the
   literal antiSMASH name -- is what matched, and check by hand whether
   that could have merged two distinct stereoisomers/compounds.

5. Persistent, cross-run, on-disk cache added: PersistentCache below writes
   a small JSON index to <cache-dir>/compound_resolution_cache.json
   (default cache-dir: results/cache/) keyed by normalized compound name,
   plus the actual resolved structure files under
   <cache-dir>/structures/<hash>.sdf. A second run over the same (or an
   overlapping) compound list reuses cached structures with zero network
   calls. Compounds that were confirmed unresolved are also cached (so
   repeated runs don't re-hit all four databases for a compound already
   known to be absent from all of them); pass --recheck-unresolved to
   retry those anyway.

6. Concurrency added. Compounds are now resolved in parallel with a
   ThreadPoolExecutor (--workers, default 4) instead of strictly
   sequentially. Each worker thread gets its own HTTPClient (thread-local),
   and the shared DecisionLog and on-disk cache are protected with locks.
   The old "time.sleep(0.3) between compounds, be polite" throttle is
   replaced by simply keeping the worker count modest by default -- these
   are shared public APIs (PubChem/EBI/NCI), not a private cluster.

7. TSV troubleshooting report added: results/reports/phase4_resolution_
   report.tsv (path configurable via --report), one row per compound, with
   columns Compound / Database / ID / Status / SearchVariantUsed / SMILES /
   MolecularWeight / StructureOrigin. This is in addition to, not instead
   of, the existing DecisionLog and JSON manifest.

8. RDKit-based chemical validation and 3D generation added (soft
   dependency). Previously, "validation" was purely syntactic (does the
   file look like a well-formed MOL/SDF block). This revision adds a real
   chemistry check on top: every candidate structure is parsed with
   RDKit's sanitizing MOL-block parser (validate_molecule_chemistry()), and
   is rejected if RDKit can't parse/sanitize it (bad valences, a malformed
   but syntactically-plausible graph, etc.) -- something the old syntactic
   check could not catch. RDKit is also now used, when a source doesn't
   itself provide a genuine 3D structure (or its native 3D download fails
   validation), to compute one deterministically from the resolved SMILES
   via ETKDGv3 embedding + MMFF94 (UFF fallback) optimization
   (generate_3d_sdf_from_smiles()). To be explicit about what "never
   invented" still means here: the compound's identity/connectivity always
   comes from a public database lookup, never from a guess -- only the 3D
   conformer *geometry* is computed locally when the source doesn't supply
   one, which is standard practice for any docking-prep pipeline (the same
   thing a tool like `obabel --gen3d` would do), not an invention of
   chemistry. If RDKit is not installed, this file still runs -- it falls
   back to syntactic-only validation and skips 3D generation, but logs a
   loud warning that this is not a publication-ready configuration.
   Install with: pip install rdkit

No organism, gene cluster, compound, accession, CID, or ChEBI/ChEMBL ID is
hardcoded anywhere in this file; the logic is entirely generic and applies
identically to any compound name antiSMASH predicts.
"""
import argparse
import hashlib
import json
import re
import sys
import threading
import time
import unicodedata
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from pks_dock.net import HTTPClient  # noqa: E402
from pks_dock.reproducibility import DecisionLog  # noqa: E402

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CIR_BASE = "https://cactus.nci.nih.gov/chemical/structure"
# ChEBI 2.0 REST API (replaces the retired legacy .do/SOAP endpoints -- see
# point 2 in the module docstring). Verify against
# https://www.ebi.ac.uk/chebi/backend/api/docs/ if this has been a while.
CHEBI_BASE = "https://www.ebi.ac.uk/chebi/backend/api/public"
CHEBI_ACCESSION_RE = re.compile(r"CHEBI:(\d+)", re.IGNORECASE)
# ChEMBL web-services API -- verify against
# https://www.ebi.ac.uk/chembl/api/data/docs if this has been a while.
CHEMBL_SEARCH = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"
CHEMBL_RECORD = "https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"

# PubChem has changed which JSON key actually carries the SMILES string for
# a given requested property tag at least once already (CanonicalSMILES ->
# ConnectivitySMILES). Rather than hard-coding today's key name, every
# plausible key is checked, in order of chemical preference (isomeric/
# stereo-aware first, since that's the more useful string for docking).
SMILES_RESPONSE_KEYS = ("IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES", "SMILES")

# HTTP statuses worth retrying: rate limiting and server-side transient
# failures. Anything else (400/401/403/404/...) is a permanent failure for
# a given request and is not retried.
TRANSIENT_STATUS_CODES = {429, 500, 502, 503, 504}
MAX_RETRIES = 3
RETRY_BACKOFF_BASE_SECONDS = 1.5

# How many ambiguous name-search candidates a single database query is
# allowed to consider before giving up on that name variant.
MAX_CANDIDATES_PER_QUERY = 5

# Deterministic seed for RDKit's ETKDG conformer embedding, so re-running
# this script on the same SMILES always produces the same 3D geometry.
RDKIT_EMBED_SEED = 0xC0FFEE

# Generic Unicode Greek letters -> spelled-out Latin names, since PubChem/
# ChEBI/ChEMBL name indices are Latin-script and natural-product names
# routinely use Greek letters (e.g. "beta-lactone" vs "\u03b2-lactone").
GREEK_TO_LATIN = {
    "\u03b1": "alpha", "\u03b2": "beta", "\u03b3": "gamma", "\u03b4": "delta",
    "\u03b5": "epsilon", "\u03b6": "zeta", "\u03b7": "eta", "\u03b8": "theta",
    "\u03b9": "iota", "\u03ba": "kappa", "\u03bb": "lambda", "\u03bc": "mu",
    "\u03bd": "nu", "\u03be": "xi", "\u03bf": "omicron", "\u03c0": "pi",
    "\u03c1": "rho", "\u03c3": "sigma", "\u03c2": "sigma", "\u03c4": "tau",
    "\u03c5": "upsilon", "\u03c6": "phi", "\u03c7": "chi", "\u03c8": "psi",
    "\u03c9": "omega",
}

# Conservative leading-stereo-descriptor patterns. Only descriptors that are
# unambiguously stereochemical notation (not part of the compound's actual
# name) are stripped, and only to build an *additional*, broader search
# variant -- the un-stripped name is always tried first. See point 4 in the
# module docstring for how this is tracked/audited downstream.
STEREO_PREFIX_RE = re.compile(
    r"^\(?\s*(?:\d+[RSrs](?:,\s*\d+[RSrs])*|[RSrs]|rel|rac|racemic|meso|"
    r"cis|trans|\u00b1|\+/-)\s*\)?[-,\u2010\u2011\u2012\u2013\u2014\s]+",
)
SIMPLE_DL_PREFIX_RE = re.compile(r"^[DdLl]-(?=[A-Za-z])")

DASH_CHARS = "\u2010\u2011\u2012\u2013\u2014\u2212"  # hyphen/dash variants
QUOTE_CHARS = "\u2018\u2019\u201c\u201d"


def log(msg: str) -> None:
    print(f"[PHASE 4] {msg}", flush=True)


# --------------------------------------------------------------------------
# Name normalization / variant generation
# --------------------------------------------------------------------------

def _normalize_unicode_and_symbols(name: str) -> str:
    """NFKC-normalize and transliterate Greek letters to their Latin names."""
    normalized = unicodedata.normalize("NFKC", name)
    for greek, latin in GREEK_TO_LATIN.items():
        normalized = normalized.replace(greek, latin)
    for dash in DASH_CHARS:
        normalized = normalized.replace(dash, "-")
    for quote in QUOTE_CHARS:
        normalized = normalized.replace(quote, "'")
    return normalized


def _normalize_whitespace_and_punctuation(name: str) -> str:
    collapsed = re.sub(r"\s+", " ", name).strip()
    collapsed = re.sub(r"\s*-\s*", "-", collapsed)  # tidy spacing around hyphens
    return collapsed.strip(" \t\r\n.,;:")


def _strip_stereo_descriptors(name: str) -> Optional[str]:
    """Return `name` with a leading stereochemical descriptor removed, or
    None if no such descriptor was found (i.e. stripping would be a no-op).
    This is used only to build a *broader* fallback search string; the
    original name is always tried first, so this never substitutes one
    compound for another -- it only widens the string used to search for
    the SAME compound under a database's indexing conventions. Every hit
    that comes from this broadened variant is recorded (search_variant_used)
    so it can be spot-checked for an accidental cross-compound merge.
    """
    stripped = STEREO_PREFIX_RE.sub("", name, count=1)
    stripped = SIMPLE_DL_PREFIX_RE.sub("", stripped, count=1)
    stripped = stripped.strip(" \t\r\n-,;:")
    if stripped and stripped.lower() != name.lower():
        return stripped
    return None


def generate_name_variants(raw_name: str) -> List[str]:
    """Build an ordered, de-duplicated list of search strings for one
    antiSMASH-predicted compound name, from most-specific to broadest.
    Every variant refers to the same compound; none of them are guesses.
    """
    variants: List[str] = []

    def add(candidate: Optional[str]) -> None:
        if candidate and candidate not in variants:
            variants.append(candidate)

    stripped_raw = raw_name.strip()
    add(stripped_raw)
    unicode_normalized = _normalize_unicode_and_symbols(stripped_raw)
    add(unicode_normalized)
    whitespace_normalized = _normalize_whitespace_and_punctuation(unicode_normalized)
    add(whitespace_normalized)
    add(_strip_stereo_descriptors(whitespace_normalized))
    return variants


# --------------------------------------------------------------------------
# Thread-local HTTP client + retry helper for transient network failures
# --------------------------------------------------------------------------

_thread_local = threading.local()


def get_client() -> HTTPClient:
    """Each worker thread gets its own HTTPClient instance. We don't know
    whether pks_dock.net.HTTPClient's internals (session, connection pool,
    rate limiter) are safe to share across threads, so the safe default
    under concurrency (see point 6) is one client per thread rather than
    one shared instance.
    """
    if not hasattr(_thread_local, "client"):
        _thread_local.client = HTTPClient(phase="PHASE 4")
    return _thread_local.client


def call_with_retries(func, *args, description: str, **kwargs):
    """Call `func(*args, **kwargs)`, retrying with exponential backoff on
    transient failures (429/500/502/503/504, timeouts, connection errors).
    Permanent failures (e.g. 404) propagate immediately without wasting
    retries.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return func(*args, **kwargs)
        except requests.HTTPError as e:
            status = e.response.status_code if e.response is not None else None
            last_exc = e
            if status not in TRANSIENT_STATUS_CODES or attempt == MAX_RETRIES:
                raise
            wait = RETRY_BACKOFF_BASE_SECONDS * (2 ** (attempt - 1))
            log(f"  Transient HTTP {status} on {description} (attempt {attempt}/{MAX_RETRIES}); "
                f"retrying in {wait:.1f}s...")
            time.sleep(wait)
        except (requests.ConnectionError, requests.Timeout) as e:
            last_exc = e
            if attempt == MAX_RETRIES:
                raise
            wait = RETRY_BACKOFF_BASE_SECONDS * (2 ** (attempt - 1))
            log(f"  Network error on {description} ({e.__class__.__name__}; "
                f"attempt {attempt}/{MAX_RETRIES}); retrying in {wait:.1f}s...")
            time.sleep(wait)
    if last_exc:
        raise last_exc
    raise RuntimeError(f"call_with_retries exhausted attempts for {description} without an exception")


def _post_json(url: str, payload: dict, *, description: str) -> dict:
    """POST helper for the one endpoint that needs it (ChEBI's
    /advanced_search/). Goes through `requests` directly rather than the
    shared HTTPClient, since the latter's interface (verified only for GET
    in the previous revision) isn't known to support POST -- wire this
    through client.post_json() instead if/when pks_dock.net grows one.
    """
    def _do():
        resp = requests.post(url, json=payload, headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        return resp.json()
    return call_with_retries(_do, description=description)


# --------------------------------------------------------------------------
# Structure validation (syntax + chemistry)
# --------------------------------------------------------------------------

def validate_sdf_text(text: Optional[str]) -> Optional[str]:
    """Return None if `text` is a plausible, non-corrupted SDF/MOL body;
    otherwise return a human-readable rejection reason. Never accepts HTML
    error pages, empty bodies, or truncated/malformed MOL blocks. This is a
    syntactic check only -- see validate_molecule_chemistry() for the
    chemistry-aware check layered on top of this one.
    """
    if not text or not text.strip():
        return "empty response body"
    stripped = text.strip()
    lowered = stripped.lower()
    if lowered.startswith("<!doctype") or lowered.startswith("<html") or "<body" in lowered[:2000]:
        return "response looks like an HTML page, not a chemical structure file"
    if "m  end" not in lowered:
        return "no 'M  END' MOL-block terminator found (truncated or malformed structure)"
    lines = stripped.splitlines()
    if len(lines) < 4:
        return "too few lines to be a valid MOL/SDF header + counts block"
    counts_line = lines[3]
    if not re.match(r"^\s*\d+\s+\d+", counts_line):
        return f"MOL counts line is not numeric atom/bond counts: {counts_line!r}"
    return None


_rdkit_warning_logged = False


def validate_molecule_chemistry(sdf_text: Optional[str]) -> Optional[str]:
    """RDKit-level chemical validation, layered on top of validate_sdf_text().
    Returns None if the structure parses to a chemically sane molecule,
    otherwise a rejection reason. Catches malformed-but-syntactically-valid
    files (e.g. broken valences, disconnected garbage) that the syntactic
    check alone cannot catch (point 8).

    If RDKit isn't installed, this can't verify anything and returns None
    (i.e. does not itself reject the structure) -- but logs a one-time
    warning, since silently skipping chemical validation is not a
    publication-ready configuration.
    """
    global _rdkit_warning_logged
    if not sdf_text:
        return "no structure to validate"
    if not RDKIT_AVAILABLE:
        if not _rdkit_warning_logged:
            log("  WARNING: rdkit not installed -- structures are only being checked "
                "syntactically, not chemically (pip install rdkit for a publication-ready run).")
            _rdkit_warning_logged = True
        return None
    try:
        mol = Chem.MolFromMolBlock(sdf_text, sanitize=True)
    except Exception as e:
        return f"RDKit raised an exception parsing the structure: {e}"
    if mol is None:
        return "RDKit could not parse/sanitize the structure (invalid valences or malformed graph)"
    if mol.GetNumAtoms() == 0:
        return "RDKit parsed the structure but it has zero atoms"
    return None


def generate_3d_sdf_from_smiles(smiles: str, mol_name: str) -> Optional[str]:
    """Build a validated 3D SDF block from a SMILES string using RDKit
    (ETKDGv3 embedding + MMFF94, falling back to UFF, geometry
    optimization). Used whenever a source database doesn't itself provide a
    genuine 3D structure, so every resolved compound still ends up with
    real 3D coordinates suitable for docking. The embedding uses a fixed
    random seed, so this is deterministic and reproducible across runs on
    the same SMILES. Returns None if RDKit is unavailable, the SMILES
    doesn't parse, or embedding fails outright.
    """
    if not RDKIT_AVAILABLE or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = RDKIT_EMBED_SEED
    conf_id = AllChem.EmbedMolecule(mol, params)
    if conf_id == -1:
        params.useRandomCoords = True
        conf_id = AllChem.EmbedMolecule(mol, params)
        if conf_id == -1:
            return None
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        pass  # keep the embedded (unoptimized) geometry rather than failing outright
    mol.SetProp("_Name", mol_name)
    return Chem.MolToMolBlock(mol) + "\n$$$$\n"


def build_sdf_from_molfile(molfile: Optional[str]) -> Optional[str]:
    """Wrap a bare MOL block (as returned by e.g. ChEMBL) into a minimal,
    valid single-record SDF by appending the '$$$$' record terminator if
    it is not already present. Returns None if the MOL block itself looks
    unusable. Note: ChEMBL's molfile is typically a 2D depiction -- see
    _finalize_candidate() for why that alone is not treated as
    docking-ready.
    """
    if not molfile:
        return None
    if "M  END" not in molfile and "m  end" not in molfile.lower():
        return None
    sdf_text = molfile if molfile.rstrip().endswith("$$$$") else molfile.rstrip() + "\n$$$$\n"
    return sdf_text


# --------------------------------------------------------------------------
# Candidate scoring (used when a name search returns multiple plausible hits)
# --------------------------------------------------------------------------

def score_candidate_structure(smiles: Optional[str], is_native_3d: bool,
                               has_required_fields: bool) -> int:
    """Objective, documented quality score for choosing among multiple
    already-VALID candidates for the same searched name. Higher is better.

    Point 3: this deliberately does NOT consider file length / line count
    -- a longer SDF/MOL block is not chemically more correct, and that
    signal has been removed entirely from this revision.
    """
    score = 0
    if is_native_3d:
        score += 100  # a source-provided 3D structure beats an RDKit-embedded one
    if smiles and ("@" in smiles or "/" in smiles or "\\" in smiles):
        score += 10  # stereo-aware SMILES is more useful for docking
    if has_required_fields:
        score += 1
    return score


# --------------------------------------------------------------------------
# Shared finalize step: validate (syntax + chemistry), generate 3D if
# needed, write the file, and score the candidate. Every source funnels
# through this one function so validation/3D-generation logic exists in
# exactly one place instead of being reimplemented per database.
# --------------------------------------------------------------------------

def _finalize_candidate(*, name: str, source: str, source_id, smiles: Optional[str],
                         molecular_weight, native_sdf_text: Optional[str],
                         native_is_3d: bool, variant_used: str, outdir: Path):
    """Returns (info_dict, sdf_text, score) for a candidate that ended up
    with a usable, validated structure, or None if nothing usable could be
    built from it at all (native structure invalid/absent AND no SMILES to
    fall back to, or RDKit couldn't embed the SMILES either).
    """
    sdf_text = None
    is_native_3d_final = False
    origin = None

    if native_sdf_text:
        reason = validate_sdf_text(native_sdf_text) or validate_molecule_chemistry(native_sdf_text)
        if reason:
            log(f"  {source}: rejected native structure for {source_id} -- {reason}")
        else:
            sdf_text = native_sdf_text
            is_native_3d_final = native_is_3d
            origin = f"{source} native {'3D' if native_is_3d else '2D'} structure"

    # If there's no valid native structure yet, or what we have is only 2D,
    # try to get real 3D coordinates from the SMILES via RDKit. A 2D native
    # hit is kept only as a last-resort fallback if RDKit generation fails.
    if sdf_text is None or not is_native_3d_final:
        if smiles:
            generated = generate_3d_sdf_from_smiles(smiles, f"{name} ({source})")
            if generated is not None:
                gen_reason = validate_sdf_text(generated) or validate_molecule_chemistry(generated)
                if gen_reason is None:
                    if sdf_text is None:
                        log(f"  {source}: no usable native structure for {source_id}; generated "
                            f"3D coordinates from SMILES with RDKit (ETKDGv3 + MMFF94).")
                    else:
                        log(f"  {source}: native structure for {source_id} was 2D only; using "
                            f"RDKit-generated 3D coordinates from the same SMILES instead.")
                    sdf_text = generated
                    is_native_3d_final = False
                    origin = f"RDKit-generated 3D (embedded from {source} SMILES)"
                else:
                    log(f"  {source}: RDKit-generated structure for {source_id} rejected -- "
                        f"{gen_reason}")
            elif not RDKIT_AVAILABLE:
                pass  # already warned once in validate_molecule_chemistry()
            else:
                log(f"  {source}: RDKit could not embed a 3D conformer for {source_id} "
                    f"(SMILES: {smiles}).")

    if sdf_text is None:
        return None

    sdf_path = outdir / f"{name.replace(' ', '_')}.sdf"
    sdf_path.write_text(sdf_text)
    info = {
        "source": source,
        "cid": source_id,
        "smiles": smiles,
        "molecular_weight": molecular_weight,
        "is_native_3d": is_native_3d_final,
        "structure_origin": origin,
        "search_variant_used": variant_used,
        "sdf_path": str(sdf_path),
    }
    score = score_candidate_structure(smiles, is_native_3d_final,
                                       has_required_fields=molecular_weight is not None)
    return info, sdf_text, score


# --------------------------------------------------------------------------
# PubChem
# --------------------------------------------------------------------------

def _pubchem_lookup_cids(client: HTTPClient, name: str) -> List[int]:
    url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    try:
        cid_json = call_with_retries(client.get_json, url, description=f"PubChem CID lookup for '{name}'")
    except requests.HTTPError as e:
        log(f"  PubChem: no CID found for '{name}' ({e}).")
        return []
    except Exception as e:
        log(f"  PubChem: CID lookup for '{name}' failed unexpectedly ({e}).")
        return []
    try:
        cids = cid_json["IdentifierList"]["CID"]
    except (KeyError, TypeError):
        log(f"  PubChem: unexpected CID response shape for '{name}': {cid_json!r}")
        return []
    return list(cids)[:MAX_CANDIDATES_PER_QUERY]


def _pubchem_fetch_cid_data(client: HTTPClient, cid: int) -> Optional[dict]:
    """Fetch SMILES/MW and the best-available *native* SDF (3D preferred,
    2D fallback) for one PubChem CID. Does not validate or write anything --
    that's centralized in _finalize_candidate().
    """
    prop_url = (
        f"{PUBCHEM_BASE}/compound/cid/{cid}/property/"
        f"IsomericSMILES,CanonicalSMILES,MolecularWeight/JSON"
    )
    try:
        props = call_with_retries(
            client.get_json, prop_url, description=f"PubChem properties for CID {cid}"
        )["PropertyTable"]["Properties"][0]
    except Exception as e:
        log(f"  PubChem: property lookup failed for CID {cid}: {e}")
        return None

    smiles = None
    for key in SMILES_RESPONSE_KEYS:
        value = props.get(key)
        if value:
            smiles = value
            break
    if smiles is None:
        log(f"  PubChem: none of the expected SMILES keys {SMILES_RESPONSE_KEYS} were present "
            f"for CID {cid}; keys actually returned: {sorted(props.keys())}")
    mw = props.get("MolecularWeight")

    native_sdf, native_is_3d = None, False
    try:
        native_sdf = call_with_retries(
            client.get_text, f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d",
            description=f"PubChem 3D SDF for CID {cid}",
        )
        native_is_3d = True
    except Exception:
        try:
            native_sdf = call_with_retries(
                client.get_text, f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF",
                description=f"PubChem 2D SDF for CID {cid}",
            )
            native_is_3d = False
        except Exception as e:
            log(f"  PubChem: no SDF (3D or 2D) retrievable for CID {cid}: {e}")
            native_sdf = None

    return {"smiles": smiles, "mw": mw, "sdf": native_sdf, "is_3d": native_is_3d}


def try_pubchem(client: HTTPClient, variants: List[str], outdir: Path, original_name: str):
    for variant in variants:
        log(f"Searching PubChem (PUG-REST) for '{variant}'...")
        cids = _pubchem_lookup_cids(client, variant)
        if not cids:
            continue

        best = None  # (info, sdf_text, score, cid)
        tried = []
        for cid in cids:
            data = _pubchem_fetch_cid_data(client, cid)
            tried.append(cid)
            if not data or not data["smiles"]:
                continue
            result = _finalize_candidate(
                name=original_name, source="PubChem", source_id=cid,
                smiles=data["smiles"], molecular_weight=data["mw"],
                native_sdf_text=data["sdf"], native_is_3d=data["is_3d"],
                variant_used=variant, outdir=outdir,
            )
            if result is None:
                continue
            info, sdf_text, score = result
            if best is None or score > best[2]:
                best = (info, sdf_text, score, cid)

        if best is None:
            log(f"  PubChem: {len(cids)} CID candidate(s) for '{variant}' but none yielded a "
                f"usable structure; broadening search...")
            continue

        info, sdf_text, score, chosen_cid = best
        alternatives = [c for c in tried if c != chosen_cid]
        if alternatives:
            log(f"  PubChem: {len(tried)} candidate CID(s) considered for '{variant}' "
                f"({tried}); selected CID {chosen_cid} (quality score {score}) over {alternatives}.")

        info["_sdf_text"] = sdf_text  # internal only, stripped before writing the manifest
        info["_alternatives"] = [f"PubChem CID {c}" for c in alternatives]
        log(f"  RESOLVED via PubChem: CID={chosen_cid}, MW={info['molecular_weight']}, "
            f"SMILES={info['smiles']} -> {info['sdf_path']}")
        return info
    return None


# --------------------------------------------------------------------------
# NCI Chemical Identifier Resolver (CIR)
# --------------------------------------------------------------------------

def try_cir(client: HTTPClient, variants: List[str], outdir: Path, original_name: str):
    for variant in variants:
        log(f"Falling back to NCI Chemical Identifier Resolver for '{variant}'...")
        try:
            smiles_url = f"{CIR_BASE}/{requests.utils.quote(variant)}/smiles"
            smiles_text = call_with_retries(
                client.get_text, smiles_url, description=f"CIR SMILES lookup for '{variant}'"
            )
        except Exception as e:
            log(f"  CIR: SMILES lookup failed for '{variant}': {e}")
            continue
        if not smiles_text.strip():
            log(f"  CIR: no SMILES found for '{variant}'.")
            continue
        smiles = smiles_text.strip().splitlines()[0]

        # CIR's own SDF-file endpoint is not used here: relying on four
        # different databases' own SDF/MOL generators made validation and
        # 3D-quality inconsistent. Instead every SMILES-only source (CIR,
        # ChEBI) is routed through the same RDKit 3D-generation step as the
        # 2D fallback path for PubChem/ChEMBL, so all resolved structures
        # get consistent, validated 3D geometry.
        result = _finalize_candidate(
            name=original_name, source="NCI_CIR", source_id=None,
            smiles=smiles, molecular_weight=None,
            native_sdf_text=None, native_is_3d=False,
            variant_used=variant, outdir=outdir,
        )
        if result is None:
            log(f"  CIR: SMILES resolved ({smiles}) for '{variant}' but no usable 3D structure "
                f"could be built from it.")
            continue
        info, sdf_text, score = result
        info["_sdf_text"] = sdf_text
        info["_alternatives"] = []
        log(f"  RESOLVED via NCI/CIR: SMILES={smiles} -> {info['sdf_path']}")
        return info
    return None


# --------------------------------------------------------------------------
# ChEBI (2.0 REST API)
# --------------------------------------------------------------------------

def _chebi_search(variant: str) -> List[str]:
    """Return up to MAX_CANDIDATES_PER_QUERY numeric ChEBI IDs (as strings,
    without the 'CHEBI:' prefix) matching `variant`, via the current
    /advanced_search/ endpoint.
    """
    url = f"{CHEBI_BASE}/advanced_search/"
    payload = {
        "text_search_specification": {"or_specification": [{"text": variant, "category": "all"}]},
        "stars": [2, 3],
    }
    try:
        raw = _post_json(url, payload, description=f"ChEBI advanced_search for '{variant}'")
    except Exception as e:
        log(f"  ChEBI: search failed for '{variant}': {e}")
        return []
    hits = raw.get("results") or []
    ids: List[str] = []
    for hit in hits[:MAX_CANDIDATES_PER_QUERY]:
        source = hit.get("_source", {}) if isinstance(hit, dict) else {}
        acc = source.get("chebi_accession", "")
        match = CHEBI_ACCESSION_RE.search(str(acc))
        if match:
            ids.append(match.group(1))
    return ids


def _chebi_fetch_compound(client: HTTPClient, numeric_id: str) -> Optional[dict]:
    url = f"{CHEBI_BASE}/compound/{numeric_id}/"
    try:
        raw = call_with_retries(
            client.get_json, url, description=f"ChEBI compound fetch for CHEBI:{numeric_id}"
        )
    except Exception as e:
        log(f"  ChEBI: compound fetch failed for CHEBI:{numeric_id}: {e}")
        return None
    struct = raw.get("default_structure") or {}
    chem = raw.get("chemical_data") or {}
    smiles = struct.get("smiles")
    mw = chem.get("mass")
    if isinstance(mw, str):
        try:
            mw = float(mw)
        except ValueError:
            mw = None
    return {"smiles": smiles, "mw": mw}


def try_chebi(client: HTTPClient, variants: List[str], outdir: Path, original_name: str):
    for variant in variants:
        log(f"Falling back to ChEBI (2.0 REST API) for '{variant}'...")
        ids = _chebi_search(variant)
        if not ids:
            log(f"  ChEBI: no entry found for '{variant}'.")
            continue

        best = None
        for numeric_id in ids:
            data = _chebi_fetch_compound(client, numeric_id)
            if not data or not data["smiles"]:
                continue
            result = _finalize_candidate(
                name=original_name, source="ChEBI", source_id=f"CHEBI:{numeric_id}",
                smiles=data["smiles"], molecular_weight=data["mw"],
                native_sdf_text=None, native_is_3d=False,
                variant_used=variant, outdir=outdir,
            )
            if result is None:
                continue
            info, sdf_text, score = result
            if best is None or score > best[2]:
                best = (info, sdf_text, score, numeric_id)

        if best is None:
            log(f"  ChEBI: {len(ids)} candidate(s) for '{variant}' but none yielded a usable "
                f"structure; broadening search...")
            continue

        info, sdf_text, score, chosen_id = best
        alternatives = [f"CHEBI:{i}" for i in ids if i != chosen_id]
        if alternatives:
            log(f"  ChEBI: {len(ids)} candidate(s) considered for '{variant}' ({ids}); "
                f"selected CHEBI:{chosen_id} (quality score {score}) over {alternatives}.")
        info["_sdf_text"] = sdf_text
        info["_alternatives"] = alternatives
        log(f"  RESOLVED via ChEBI: CHEBI:{chosen_id} -> {info['sdf_path']}")
        return info
    return None


# --------------------------------------------------------------------------
# ChEMBL
# --------------------------------------------------------------------------

def try_chembl(client: HTTPClient, variants: List[str], outdir: Path, original_name: str):
    for variant in variants:
        log(f"Falling back to ChEMBL for '{variant}'...")
        try:
            search_json = call_with_retries(
                client.get_json, CHEMBL_SEARCH,
                params={"q": variant, "limit": MAX_CANDIDATES_PER_QUERY, "format": "json"},
                description=f"ChEMBL search for '{variant}'",
            )
        except Exception as e:
            log(f"  ChEMBL: search failed for '{variant}': {e}")
            continue
        molecules = search_json.get("molecules") or []
        if not molecules:
            log(f"  ChEMBL: no entry found for '{variant}'.")
            continue

        best = None
        tried_ids = []
        for molecule in molecules:
            chembl_id = molecule.get("molecule_chembl_id")
            if not chembl_id:
                continue
            tried_ids.append(chembl_id)
            try:
                record = call_with_retries(
                    client.get_json, CHEMBL_RECORD.format(chembl_id=chembl_id),
                    description=f"ChEMBL record fetch for {chembl_id}",
                )
            except Exception as e:
                log(f"  ChEMBL: record fetch failed for {chembl_id}: {e}")
                continue
            structures = record.get("molecule_structures") or {}
            smiles = structures.get("canonical_smiles")
            molfile = structures.get("molfile")
            mw = (record.get("molecule_properties") or {}).get("full_mwt")
            # ChEMBL's molfile is a 2D depiction, not a docking-ready 3D
            # structure. It is passed through as a *native 2D* candidate,
            # but _finalize_candidate() will prefer generating real 3D
            # coordinates from the SMILES with RDKit rather than accepting
            # the 2D file as the final structure.
            native_sdf = build_sdf_from_molfile(molfile) if molfile else None
            if not smiles and not native_sdf:
                log(f"  ChEMBL: {chembl_id} has neither a usable SMILES nor MOL block; skipping.")
                continue
            result = _finalize_candidate(
                name=original_name, source="ChEMBL", source_id=chembl_id,
                smiles=smiles, molecular_weight=mw,
                native_sdf_text=native_sdf, native_is_3d=False,
                variant_used=variant, outdir=outdir,
            )
            if result is None:
                continue
            info, sdf_text, score = result
            if best is None or score > best[2]:
                best = (info, sdf_text, score, chembl_id)

        if best is None:
            log(f"  ChEMBL: {len(tried_ids)} candidate(s) for '{variant}' but none yielded a "
                f"usable structure; broadening search...")
            continue

        info, sdf_text, score, chosen_id = best
        alternatives = [f"ChEMBL {i}" for i in tried_ids if i != chosen_id]
        info["_sdf_text"] = sdf_text
        info["_alternatives"] = alternatives
        log(f"  RESOLVED via ChEMBL: {chosen_id}, SMILES={info['smiles']} -> {info['sdf_path']}")
        return info
    return None


# --------------------------------------------------------------------------
# Persistent, cross-run cache (point 5)
# --------------------------------------------------------------------------

class PersistentCache:
    """On-disk cache keyed by normalized compound name, persisting across
    runs (not just within one run, unlike the previous revision's in-memory
    dict). Two things are stored:

      - a small JSON index at <cache_dir>/compound_resolution_cache.json
        with the resolved metadata (source, id, SMILES, MW, whether the
        structure is a genuine source-native 3D structure, which search
        variant matched) -- or `null` for a compound confirmed unresolved;
      - the actual resolved SDF text, one small file per cached compound,
        under <cache_dir>/structures/, so a cache hit costs zero network
        calls and doesn't require regenerating a 3D conformer.

    A single instance is shared read-only-ish across worker threads; all
    mutating access goes through a lock.
    """

    def __init__(self, cache_dir: Path):
        self.cache_dir = cache_dir
        self.index_path = cache_dir / "compound_resolution_cache.json"
        self.structures_dir = cache_dir / "structures"
        self._lock = threading.Lock()
        self._data: Dict[str, Optional[dict]] = {}
        if self.index_path.exists():
            try:
                self._data = json.loads(self.index_path.read_text())
            except Exception as e:
                log(f"Cache: could not read {self.index_path} ({e}); starting with an empty cache.")
                self._data = {}

    def has(self, key: str) -> bool:
        with self._lock:
            return key in self._data

    def get(self, key: str) -> Optional[dict]:
        with self._lock:
            return self._data.get(key)

    def _structure_path(self, key: str) -> Path:
        digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]
        return self.structures_dir / f"{digest}.sdf"

    def store_resolved(self, key: str, info: dict, sdf_text: str) -> None:
        self.structures_dir.mkdir(parents=True, exist_ok=True)
        struct_path = self._structure_path(key)
        struct_path.write_text(sdf_text)
        entry = dict(info)
        entry.pop("sdf_path", None)  # per-run output path; not meaningful across runs
        entry["cached_structure_file"] = str(struct_path)
        entry["cached_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with self._lock:
            self._data[key] = entry
            self._flush_locked()

    def store_unresolved(self, key: str) -> None:
        with self._lock:
            self._data[key] = None
            self._flush_locked()

    def _flush_locked(self) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = self.index_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self._data, indent=2, sort_keys=True))
        tmp.replace(self.index_path)


def _cache_key(name: str) -> str:
    return _normalize_whitespace_and_punctuation(_normalize_unicode_and_symbols(name)).lower()


# --------------------------------------------------------------------------
# Per-compound resolution
# --------------------------------------------------------------------------

decisions_lock = threading.Lock()


def record_decision(decisions: DecisionLog, *args, **kwargs) -> None:
    """DecisionLog is shared across worker threads; serialize writes to it.
    We don't know its internals are thread-safe, so don't assume they are.
    """
    with decisions_lock:
        decisions.record(*args, **kwargs)


def resolve_one_compound(name: str, outdir: Path, cache: PersistentCache,
                          decisions: DecisionLog, recheck_unresolved: bool
                          ) -> Tuple[Optional[dict], dict]:
    """Resolve a single antiSMASH-predicted compound name to a validated
    structure. Returns (info_or_None, report_row_dict). Runs inside a
    worker thread (see main()); uses a thread-local HTTPClient.
    """
    key = _cache_key(name)

    if cache.has(key):
        entry = cache.get(key)
        if entry is None:
            if not recheck_unresolved:
                log(f"'{name}' -- cached as unresolved from a previous run; skipping "
                    f"(pass --recheck-unresolved to retry).")
                record_decision(
                    decisions, "PHASE 4", f"Structure resolution for compound '{name}'",
                    database=None, selected=None,
                    alternatives_considered=["PubChem", "NCI_CIR", "ChEBI", "ChEMBL"],
                    reason=f"Cached as unresolved from a previous run ({cache.index_path}).",
                    confidence="none",
                )
                return None, {"database": "", "id": "", "status": "unresolved (cached)",
                               "variant": "", "smiles": "", "mw": "", "origin": ""}
            log(f"'{name}' -- cached as unresolved previously; --recheck-unresolved set, retrying.")
        else:
            struct_file = Path(entry.get("cached_structure_file", ""))
            if struct_file.exists():
                sdf_path = outdir / f"{name.replace(' ', '_')}.sdf"
                sdf_path.write_text(struct_file.read_text())
                info = {k: v for k, v in entry.items()
                        if k not in ("cached_structure_file", "cached_at")}
                info["sdf_path"] = str(sdf_path)
                log(f"'{name}' -- resolved from cache (source: {info['source']}, "
                    f"no network lookup).")
                record_decision(
                    decisions, "PHASE 4", f"Structure resolution for compound '{name}'",
                    database=info["source"], selected=info.get("cid") or info.get("smiles"),
                    alternatives_considered=[],
                    reason=f"Reused cached structure from a previous run ({cache.index_path}).",
                    confidence="high" if info.get("is_native_3d") else "medium",
                )
                row = {"database": info["source"], "id": info.get("cid") or "",
                       "status": "resolved (cached)", "variant": info.get("search_variant_used", ""),
                       "smiles": info.get("smiles") or "", "mw": info.get("molecular_weight") or "",
                       "origin": info.get("structure_origin") or ""}
                return info, row
            log(f"'{name}' -- cache entry found but its structure file is missing on disk; "
                f"re-resolving.")

    client = get_client()
    variants = generate_name_variants(name)
    log(f"Resolving '{name}' -- {len(variants)} search variant(s) will be tried per "
        f"database if needed: {variants}")

    sources = [
        ("PubChem", lambda: try_pubchem(client, variants, outdir, name)),
        ("NCI_CIR", lambda: try_cir(client, variants, outdir, name)),
        ("ChEBI", lambda: try_chebi(client, variants, outdir, name)),
        ("ChEMBL", lambda: try_chembl(client, variants, outdir, name)),
    ]
    try:
        result = client.fallback_chain(name, sources)
        info = result.value
        sdf_text = info.pop("_sdf_text")
        alternatives_within_source = info.pop("_alternatives", [])
        record_decision(
            decisions, "PHASE 4", f"Structure resolution for compound '{name}'",
            database=result.source, selected=info.get("cid") or info.get("smiles"),
            alternatives_considered=[a.source for a in result.attempts if not a.success]
            + alternatives_within_source,
            reason=f"First successful source in PubChem -> NCI/CIR -> ChEBI -> ChEMBL chain: "
                   f"{result.source}. Structure: {info.get('structure_origin')}.",
            confidence="high" if info.get("is_native_3d") else "medium",
        )
        cache.store_resolved(key, info, sdf_text)
        row = {"database": info["source"], "id": info.get("cid") or "", "status": "resolved",
               "variant": info.get("search_variant_used", ""), "smiles": info.get("smiles") or "",
               "mw": info.get("molecular_weight") or "", "origin": info.get("structure_origin") or ""}
        return info, row
    except Exception:
        log(f"UNRESOLVED: '{name}' -- not found in PubChem, NCI/CIR, ChEBI, or ChEMBL "
            f"(across {len(variants)} name variant(s) each). This compound must be sourced "
            f"from the primary literature that first characterised it and added manually to "
            f"{outdir}/ as a real SDF file.")
        record_decision(
            decisions, "PHASE 4", f"Structure resolution for compound '{name}'",
            database=None, selected=None,
            alternatives_considered=["PubChem", "NCI_CIR", "ChEBI", "ChEMBL"],
            reason="All four sources exhausted across all name variants; compound requires "
                   "manual literature sourcing.",
            confidence="none",
        )
        cache.store_unresolved(key)
        return None, {"database": "", "id": "", "status": "unresolved", "variant": "",
                       "smiles": "", "mw": "", "origin": ""}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compound-list", default="results/antismash/phase3_compound_list.json")
    ap.add_argument("--outdir", default="results/ligands")
    ap.add_argument("--decision-log", default="results/reports/decision_log.json")
    ap.add_argument("--report", default="results/reports/phase4_resolution_report.tsv",
                     help="Per-compound troubleshooting TSV (point 7).")
    ap.add_argument("--cache-dir", default="results/cache",
                     help="Cross-run cache directory (point 5).")
    ap.add_argument("--workers", type=int, default=4,
                     help="Concurrent compound-resolution workers (point 6). Keep modest -- "
                          "these are shared public APIs, not a private cluster.")
    ap.add_argument("--recheck-unresolved", action="store_true",
                     help="Re-attempt compounds cached as unresolved in a previous run instead "
                          "of skipping them.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)

    if not RDKIT_AVAILABLE:
        log("WARNING: rdkit is not installed. 3D-coordinate generation and chemical-level "
            "structure validation (beyond basic SDF syntax) are both disabled -- this is NOT a "
            "publication-ready configuration. Install with: pip install rdkit")

    decisions = DecisionLog(args.decision_log)
    cache = PersistentCache(Path(args.cache_dir))

    clusters = json.loads(Path(args.compound_list).read_text())
    names = sorted({c["predicted_compound"] for c in clusters
                    if c["predicted_compound"] != "unannotated PKS-I cluster"})

    if not names:
        log("No named compounds to resolve (all clusters were unannotated). "
            "Nothing to dock -- check Phase 3 output / knownclusterblast results.")

    log(f"{len(names)} named compound(s) to resolve: {names}")
    log(f"Using up to {max(1, args.workers)} concurrent worker(s); cache dir: {args.cache_dir}")

    resolved: Dict[str, dict] = {}
    unresolved: List[str] = []
    report_rows: List[Tuple[str, str, str, str, str, str, str, str]] = []

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {
            pool.submit(resolve_one_compound, name, outdir, cache, decisions,
                        args.recheck_unresolved): name
            for name in names
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                info, row = future.result()
            except Exception as e:
                log(f"UNEXPECTED ERROR resolving '{name}': {e}")
                info, row = None, {"database": "", "id": "", "status": f"error: {e}",
                                    "variant": "", "smiles": "", "mw": "", "origin": ""}
            if info is not None:
                resolved[name] = info
            else:
                unresolved.append(name)
            report_rows.append((
                name, str(row["database"]), str(row["id"]), str(row["status"]),
                str(row["variant"]), str(row["smiles"] or ""), str(row["mw"] or ""),
                str(row["origin"] or ""),
            ))

    # FIX (carried forward): the resolved-compounds manifest is written
    # inside --outdir (results/ligands/), not its parent -- with the
    # default --outdir this is results/ligands/phase4_resolved_compounds.json,
    # exactly where Phase 12 (ADMET) looks for it by default.
    manifest = outdir / "phase4_resolved_compounds.json"
    manifest.write_text(json.dumps({"resolved": resolved, "unresolved": unresolved}, indent=2))
    log(f"Summary: {len(resolved)} resolved, {len(unresolved)} unresolved.")
    if unresolved:
        log(f"Manual literature follow-up needed for: {unresolved}")
    log(f"Wrote {manifest} -- read automatically downstream by ADMET (Phase 12). "
        f"No SMILES or structure is ever hand-typed again after this point.")

    # Point 7: human-readable troubleshooting report, ordered to match the
    # original `names` list regardless of the order concurrent workers
    # actually finished in.
    name_order = {n: i for i, n in enumerate(names)}
    report_rows.sort(key=lambda r: name_order.get(r[0], len(names)))
    report_path = Path(args.report)
    with report_path.open("w") as fh:
        fh.write("Compound\tDatabase\tID\tStatus\tSearchVariantUsed\tSMILES\tMolecularWeight\t"
                  "StructureOrigin\n")
        for row in report_rows:
            fh.write("\t".join(row) + "\n")
    log(f"Wrote troubleshooting report: {report_path}")


if __name__ == "__main__":
    main()

    