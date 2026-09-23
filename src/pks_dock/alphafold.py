"""
pks_dock.alphafold
====================

Real integration with the AlphaFold Protein Structure Database public API.

WHY THIS MODULE EXISTS
-----------------------
The previous version of Phase 6b fetched a sequence from UniProt and then
went straight to a *local* ColabFold prediction, even when a perfectly good
prediction for that exact UniProt accession already existed for free in the
AlphaFold DB. That wastes hours of GPU time and is not what the AlphaFold DB
API is for. This module closes that gap with the documented, public
endpoints:

    Prediction lookup:  https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}
    Model file:          https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb

Reference: EMBL-EBI AlphaFold DB API documentation,
https://alphafold.ebi.ac.uk/api-docs

Flow implemented here (matches the brief's "AlphaFold Integration" spec):

    UniProt accession
        -> query AlphaFold DB prediction API
        -> exists?
            yes -> download model + confidence (pLDDT) metadata
            no  -> report clearly; caller falls back to local ColabFold
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .net import HTTPClient

ALPHAFOLD_PREDICTION_API = "https://alphafold.ebi.ac.uk/api/prediction"


@dataclass
class AlphaFoldHit:
    uniprot_accession: str
    entry_id: str
    model_url: str
    cif_url: Optional[str]
    confidence_url: Optional[str]
    mean_plddt: Optional[float]
    model_version: Optional[str]
    latest_version_date: Optional[str]


def query_alphafold_db(client: HTTPClient, uniprot_accession: str) -> Optional[AlphaFoldHit]:
    """
    Queries the AlphaFold DB for an existing prediction for a UniProt
    accession. Returns None (not an exception) if the accession genuinely
    has no prediction -- that is a valid, expected outcome, not a failure.
    Raises on genuine transient/network failure so the retry logic in
    HTTPClient can do its job.
    """
    url = f"{ALPHAFOLD_PREDICTION_API}/{uniprot_accession}"
    try:
        data = client.get_json(url)
    except Exception as exc:
        # A 404 surfaces as an HTTPError via requests.raise_for_status();
        # distinguish "no prediction exists" from "the API is down".
        status = getattr(getattr(exc, "response", None), "status_code", None)
        if status == 404:
            return None
        raise

    if not data:
        return None
    entry = data[0] if isinstance(data, list) else data
    return AlphaFoldHit(
        uniprot_accession=uniprot_accession,
        entry_id=entry.get("entryId", f"AF-{uniprot_accession}-F1"),
        model_url=entry.get("pdbUrl") or entry.get("modelUrl", ""),
        cif_url=entry.get("cifUrl"),
        confidence_url=entry.get("paeImageUrl") or entry.get("paeDocUrl"),
        mean_plddt=entry.get("globalMetricValue") or entry.get("confidenceAvgLocalScore"),
        model_version=str(entry.get("latestVersion")) if entry.get("latestVersion") is not None else None,
        latest_version_date=entry.get("modelCreatedDate"),
    )


def download_alphafold_model(client: HTTPClient, hit: AlphaFoldHit, out_dir: Path) -> Path:
    """Downloads the PDB model file for a resolved AlphaFoldHit."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{hit.entry_id}.pdb"
    if not hit.model_url:
        raise ValueError(f"AlphaFold hit for {hit.uniprot_accession} has no model URL")
    content = client.get_binary(hit.model_url)
    out_path.write_bytes(content)
    return out_path


def resolve_receptor_structure(
    client: HTTPClient,
    uniprot_accession: str,
    out_dir: Path,
) -> tuple[Optional[Path], Optional[AlphaFoldHit], str]:
    """
    High-level convenience used by Phase 6b: given a UniProt accession,
    tries the AlphaFold DB first. Returns (path_or_None, hit_or_None, note).

    A None path with a non-empty note is the expected, honest outcome when
    no prediction exists -- the caller should fall through to local
    ColabFold prediction and log the note verbatim in the decision log.
    """
    try:
        hit = query_alphafold_db(client, uniprot_accession)
    except Exception as exc:
        return None, None, f"AlphaFold DB API request failed after retries: {exc}"

    if hit is None:
        return (
            None,
            None,
            (
                f"No existing AlphaFold DB prediction for {uniprot_accession}; "
                "falling back to local ColabFold prediction."
            ),
        )

    try:
        path = download_alphafold_model(client, hit, out_dir)
    except Exception as exc:
        return None, hit, f"Found AlphaFold DB entry {hit.entry_id} but download failed: {exc}"

    plddt_note = f", mean pLDDT {hit.mean_plddt:.1f}" if hit.mean_plddt else ""
    return (
        path,
        hit,
        (
            f"Used existing AlphaFold DB prediction {hit.entry_id}{plddt_note} "
            "instead of running a new local prediction."
        ),
    )
