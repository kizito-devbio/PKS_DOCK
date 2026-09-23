"""
Tests for pks_dock.reproducibility: decision_log.json accumulation,
workflow_report.md rendering, and pipeline_metadata.json provenance.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from pks_dock.reproducibility import (  # noqa: E402
    DecisionLog,
    render_workflow_report,
    write_pipeline_metadata,
)


def test_decision_log_persists_across_instances(tmp_path):
    log_path = tmp_path / "decision_log.json"
    log1 = DecisionLog(log_path)
    log1.record("PHASE 1", "Selected genome", database="NCBI", selected="GCF_000001.1")

    log2 = DecisionLog(log_path)  # simulates the next phase script starting fresh
    assert len(log2.entries()) == 1
    assert log2.entries()[0]["selected"] == "GCF_000001.1"

    log2.record("PHASE 4", "Resolved compound", database="PubChem", selected="CID123")
    log3 = DecisionLog(log_path)
    assert len(log3.entries()) == 2


def test_decision_log_survives_corrupted_file(tmp_path):
    log_path = tmp_path / "decision_log.json"
    log_path.write_text("{not valid json")
    log = DecisionLog(log_path)
    assert log.entries() == []
    log.record("PHASE 1", "test")
    assert len(log.entries()) == 1


def test_render_workflow_report_groups_by_phase(tmp_path):
    log_path = tmp_path / "decision_log.json"
    log = DecisionLog(log_path)
    log.record("PHASE 1", "Selected genome", database="NCBI", selected="GCF_1", reason="best N50")
    log.record(
        "PHASE 6",
        "Selected receptor",
        database="RCSB",
        selected="4JVC",
        alternatives_considered=["1ABC", "2XYZ"],
        reason="highest resolution",
    )

    out_path = tmp_path / "workflow_report.md"
    render_workflow_report(log_path, out_path)
    text = out_path.read_text()

    assert "PHASE 1" in text
    assert "PHASE 6" in text
    assert "GCF_1" in text
    assert "4JVC" in text
    assert "1ABC" in text


def test_render_workflow_report_handles_empty_log(tmp_path):
    out_path = tmp_path / "workflow_report.md"
    render_workflow_report(tmp_path / "does_not_exist.json", out_path)
    assert "No automated decisions" in out_path.read_text()


def test_write_pipeline_metadata_captures_run_params(tmp_path):
    out_path = tmp_path / "pipeline_metadata.json"
    write_pipeline_metadata({"organism": "Bacillus velezensis", "threads": 8}, out_path)
    meta = json.loads(out_path.read_text())
    assert meta["run_parameters"]["organism"] == "Bacillus velezensis"
    assert "pks_dock_version" in meta
    assert "python_version" in meta
