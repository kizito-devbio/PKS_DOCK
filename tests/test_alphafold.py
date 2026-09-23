"""
Tests for pks_dock.alphafold: the real AlphaFold DB API integration used
by Phase 6b before it falls back to local ColabFold prediction.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pks_dock.alphafold import (  # noqa: E402
    AlphaFoldHit,
    query_alphafold_db,
    resolve_receptor_structure,
)

SAMPLE_ENTRY = {
    "entryId": "AF-P0A6F1-F1",
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P0A6F1-F1-model_v4.pdb",
    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P0A6F1-F1-model_v4.cif",
    "paeImageUrl": "https://alphafold.ebi.ac.uk/files/AF-P0A6F1-F1-predicted_aligned_error_v4.png",
    "globalMetricValue": 91.4,
    "latestVersion": 4,
    "modelCreatedDate": "2022-06-01",
}


def test_query_alphafold_db_found():
    client = MagicMock()
    client.get_json.return_value = [SAMPLE_ENTRY]
    hit = query_alphafold_db(client, "P0A6F1")
    assert isinstance(hit, AlphaFoldHit)
    assert hit.entry_id == "AF-P0A6F1-F1"
    assert hit.mean_plddt == 91.4
    client.get_json.assert_called_once()


def test_query_alphafold_db_not_found_returns_none():
    client = MagicMock()
    resp = MagicMock(status_code=404)
    client.get_json.side_effect = requests.HTTPError("404", response=resp)
    hit = query_alphafold_db(client, "NOSUCHACCESSION")
    assert hit is None


def test_query_alphafold_db_transient_error_propagates():
    client = MagicMock()
    resp = MagicMock(status_code=503)
    client.get_json.side_effect = requests.HTTPError("503", response=resp)
    with pytest.raises(requests.HTTPError):
        query_alphafold_db(client, "P0A6F1")


def test_resolve_receptor_structure_downloads_when_found(tmp_path):
    client = MagicMock()
    client.get_json.return_value = [SAMPLE_ENTRY]
    client.get_binary.return_value = b"HEADER    FAKE PDB CONTENT\n"

    path, hit, note = resolve_receptor_structure(client, "P0A6F1", tmp_path)
    assert path is not None
    assert path.exists()
    assert path.read_bytes().startswith(b"HEADER")
    assert hit.entry_id == "AF-P0A6F1-F1"
    assert "AF-P0A6F1-F1" in note
    assert "91.4" in note


def test_resolve_receptor_structure_falls_back_when_absent(tmp_path):
    client = MagicMock()
    resp = MagicMock(status_code=404)
    client.get_json.side_effect = requests.HTTPError("404", response=resp)

    path, hit, note = resolve_receptor_structure(client, "NOSUCHACCESSION", tmp_path)
    assert path is None
    assert hit is None
    assert "ColabFold" in note


def test_resolve_receptor_structure_reports_api_failure_honestly(tmp_path):
    client = MagicMock()
    resp = MagicMock(status_code=503)
    client.get_json.side_effect = requests.HTTPError("503", response=resp)

    path, hit, note = resolve_receptor_structure(client, "P0A6F1", tmp_path)
    assert path is None
    assert hit is None
    assert "failed" in note.lower()
