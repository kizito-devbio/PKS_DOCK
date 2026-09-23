"""
pks_dock.reproducibility
=========================

Every phase script imports `DecisionLog` and calls `.record(...)` at each
point where the pipeline made an automated choice (which PDB structure,
which pocket, which fallback database, which pocket-derived grid box, ...).

At the end of a run, `scripts/16_generate_reproducibility_report.py` calls
`render_workflow_report()` and `write_pipeline_metadata()` to turn the
accumulated decisions into the three artifacts the pipeline promises:

- results/reports/decision_log.json      (machine-readable, append-only)
- results/reports/workflow_report.md     (human-readable narrative)
- results/reports/pipeline_metadata.json (run provenance: versions, times,
                                           inputs, environment)

Design choice: the log is a flat JSON-lines-style list under one file,
written with an atomic replace on every `.record()` call, so a crash in
phase 9 does not lose the decisions recorded in phases 1-8.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

PKS_DOCK_VERSION = "0.2.0"

DEFAULT_LOG_PATH = Path("results/reports/decision_log.json")


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class Decision:
    phase: str
    action: str
    database: Optional[str] = None
    query: Optional[str] = None
    selected: Optional[str] = None
    alternatives_considered: list[str] = field(default_factory=list)
    reason: str = ""
    confidence: Optional[str] = None
    api_response_summary: Optional[str] = None
    timestamp: str = field(default_factory=_utcnow_iso)


class DecisionLog:
    """Append-only, crash-safe decision log shared across all phase scripts."""

    def __init__(self, path: Path | str = DEFAULT_LOG_PATH):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._entries: list[dict] = self._load()

    def _load(self) -> list[dict]:
        if self.path.exists():
            try:
                return json.loads(self.path.read_text())
            except (json.JSONDecodeError, OSError):
                return []
        return []

    def record(
        self,
        phase: str,
        action: str,
        *,
        database: Optional[str] = None,
        query: Optional[str] = None,
        selected: Optional[str] = None,
        alternatives_considered: Optional[list[str]] = None,
        reason: str = "",
        confidence: Optional[str] = None,
        api_response_summary: Optional[str] = None,
    ) -> Decision:
        decision = Decision(
            phase=phase,
            action=action,
            database=database,
            query=query,
            selected=selected,
            alternatives_considered=alternatives_considered or [],
            reason=reason,
            confidence=confidence,
            api_response_summary=api_response_summary,
        )
        self._entries.append(asdict(decision))
        self._flush()
        return decision

    def _flush(self) -> None:
        tmp = self.path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self._entries, indent=2))
        os.replace(tmp, self.path)  # atomic on POSIX

    def entries(self) -> list[dict]:
        return list(self._entries)


# -- report rendering ------------------------------------------------------


def _git_commit() -> Optional[str]:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def write_pipeline_metadata(
    run_params: dict[str, Any],
    out_path: Path | str = "results/reports/pipeline_metadata.json",
) -> Path:
    """Captures software/environment provenance for a single pipeline run."""
    meta = {
        "pks_dock_version": PKS_DOCK_VERSION,
        "git_commit": _git_commit(),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "run_started_utc": run_params.get("run_started_utc", _utcnow_iso()),
        "run_parameters": run_params,
    }
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(meta, indent=2))
    return out


def render_workflow_report(
    log_path: Path | str = DEFAULT_LOG_PATH,
    out_path: Path | str = "results/reports/workflow_report.md",
) -> Path:
    """
    Turns the accumulated decision_log.json into a human-readable narrative,
    grouped by phase, in the same "what/why/alternatives rejected" style
    the pipeline prints to the terminal.
    """
    log_path = Path(log_path)
    entries = json.loads(log_path.read_text()) if log_path.exists() else []

    by_phase: dict[str, list[dict]] = {}
    for e in entries:
        by_phase.setdefault(e["phase"], []).append(e)

    lines = ["# PKS_DOCK Workflow Report", "", f"_Generated {_utcnow_iso()}_", ""]
    if not entries:
        lines.append("No automated decisions were recorded for this run.")
    for phase in sorted(by_phase):
        lines.append(f"## {phase}")
        lines.append("")
        for e in by_phase[phase]:
            lines.append(f"- **{e['action']}**")
            if e.get("database"):
                lines.append(f"  - Database: {e['database']}")
            if e.get("query"):
                lines.append(f"  - Query: `{e['query']}`")
            if e.get("selected"):
                lines.append(f"  - Selected: **{e['selected']}**")
            if e.get("alternatives_considered"):
                lines.append(f"  - Alternatives considered: {', '.join(e['alternatives_considered'])}")
            if e.get("reason"):
                lines.append(f"  - Reason: {e['reason']}")
            if e.get("confidence"):
                lines.append(f"  - Confidence: {e['confidence']}")
            lines.append(f"  - Timestamp (UTC): {e['timestamp']}")
            lines.append("")
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines))
    return out
