#!/usr/bin/env python3
"""
Phase 16 - Reproducibility finalization.

Renders results/reports/workflow_report.md from the decision_log.json that
phases 1, 4, 6, and 6b have been writing to incrementally, and writes
pipeline_metadata.json with software/environment provenance for this run.

This is also done automatically by the installed `pks-dock` CLI wrapper
(src/pks_dock/cli.py) after run_pipeline.sh exits, including on failure.
This script exists as a standalone step for people running run_pipeline.sh
directly rather than through `pks-dock`.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from pks_dock.reproducibility import render_workflow_report, write_pipeline_metadata  # noqa: E402

PHASE = "PHASE 16"


def log(msg):
    print(f"[{PHASE}] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organism", default=None)
    ap.add_argument("--accession", default=None)
    ap.add_argument("--pathogens", nargs="*", default=[])
    ap.add_argument("--decision-log", default="results/reports/decision_log.json")
    ap.add_argument("--metadata-out", default="results/reports/pipeline_metadata.json")
    ap.add_argument("--report-out", default="results/reports/workflow_report.md")
    args = ap.parse_args()

    log("Finalizing reproducibility artifacts...")
    write_pipeline_metadata(
        {"organism": args.organism, "accession": args.accession, "pathogens": args.pathogens},
        out_path=args.metadata_out,
    )
    log(f"Wrote {args.metadata_out}")

    report_path = render_workflow_report(log_path=args.decision_log, out_path=args.report_out)
    log(f"Wrote {report_path}")
    log("Phase 16 complete. decision_log.json, workflow_report.md, and "
        "pipeline_metadata.json are all in results/reports/.")


if __name__ == "__main__":
    main()
