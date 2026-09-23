"""
pks_dock.cli
============

Installable command-line entry point (`pks-dock`) for PKS_DOCK.

This is a thin wrapper -- not a reimplementation -- around
`run_pipeline_container.sh`, which remains the source of truth for
pipeline phase ordering and execution.

The repository-level `run_pipeline.sh` is the Docker launcher used by
end users. It mounts the PKS_DOCK project directory into the container
so that results and logs are written directly to the user's filesystem.

Inside the container, `pks-dock` executes `run_pipeline_container.sh`,
which contains the complete 16-phase pipeline.

The Python CLI adds:

  * A proper `pks-dock` executable when PKS_DOCK is installed as a package.
  * Validation that the pipeline script can be located before execution.
  * Forwarding of pipeline parameters without reimplementing the workflow.
  * Finalization of reproducibility artifacts after the pipeline exits,
    including when a phase fails partway through.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .reproducibility import render_workflow_report, write_pipeline_metadata


PIPELINE_SCRIPT = "run_pipeline_container.sh"


def _find_repo_root() -> Path:
    """
    Locate the PKS_DOCK repository root.

    The repository root is identified by the presence of the actual
    pipeline script, `run_pipeline_container.sh`.

    This supports execution from a PKS_DOCK repository checkout.
    """
    here = Path(__file__).resolve()

    for candidate in (here.parent.parent.parent, Path.cwd()):
        if (candidate / PIPELINE_SCRIPT).exists():
            return candidate

    sys.exit(
        "[pks-dock] ERROR: could not locate "
        f"{PIPELINE_SCRIPT}. Run `pks-dock` from inside a PKS_DOCK "
        "repository checkout."
    )


def main(argv: list[str] | None = None) -> int:
    """Run the PKS_DOCK pipeline and finalize reproducibility artifacts."""

    parser = argparse.ArgumentParser(
        prog="pks-dock",
        description=(
            "Fully automated PKS-I mining, compound retrieval, "
            "and molecular docking pipeline."
        ),
    )

    parser.add_argument(
        "--organism",
        help='PKS-producing organism name, e.g. "Bacillus velezensis"',
    )

    parser.add_argument(
        "--accession",
        help="Specific NCBI genome accession, if already known",
    )

    parser.add_argument(
        "--pathogens",
        nargs="+",
        help=(
            "Pathogen keys from config/pathogen_targets.yaml, "
            "e.g. Staphylococcus_aureus"
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads used by supported pipeline steps (default: 4).",
    )

    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=16,
        help="AutoDock Vina exhaustiveness (default: 16).",
    )

    parser.add_argument(
        "--num-modes",
        type=int,
        default=9,
        help="Number of docking poses requested from AutoDock Vina (default: 9).",
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=8,
        help="CPU count passed to the docking step (default: 8).",
    )

    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the PKS_DOCK version and exit.",
    )

    args = parser.parse_args(argv)

    if args.version:
        from . import __version__

        print(f"pks-dock {__version__}")
        return 0

    if not args.pathogens:
        parser.error("the following arguments are required: --pathogens")

    if not args.organism and not args.accession:
        parser.error("supply either --organism or --accession")

    repo_root = _find_repo_root()
    pipeline_script = repo_root / PIPELINE_SCRIPT

    run_started = datetime.now(timezone.utc).isoformat(timespec="seconds")

    cmd = [str(pipeline_script)]

    if args.organism:
        cmd += ["--organism", args.organism]

    if args.accession:
        cmd += ["--accession", args.accession]

    cmd += ["--pathogens", *args.pathogens]
    cmd += ["--threads", str(args.threads)]
    cmd += ["--exhaustiveness", str(args.exhaustiveness)]
    cmd += ["--num-modes", str(args.num_modes)]
    cmd += ["--cpu", str(args.cpu)]

    print(
        f"[pks-dock] Running: {' '.join(cmd)}",
        flush=True,
    )

    result = subprocess.run(
        cmd,
        cwd=repo_root,
    )

    # Finalize reproducibility metadata even when the pipeline exits
    # unsuccessfully, so partial runs retain an auditable record.
    write_pipeline_metadata(
        {
            "organism": args.organism,
            "accession": args.accession,
            "pathogens": args.pathogens,
            "threads": args.threads,
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "cpu": args.cpu,
            "run_started_utc": run_started,
            "exit_code": result.returncode,
        },
        out_path=repo_root / "results/reports/pipeline_metadata.json",
    )

    report_path = render_workflow_report(
        log_path=repo_root / "results/reports/decision_log.json",
        out_path=repo_root / "results/reports/workflow_report.md",
    )

    print(
        f"[pks-dock] Reproducibility report written -> {report_path}",
        flush=True,
    )

    if result.returncode != 0:
        print(
            f"[pks-dock] Pipeline exited with status {result.returncode}. "
            "See the reproducibility report for what completed before "
            "the failure.",
            file=sys.stderr,
        )

    return result.returncode


if __name__ == "__main__":
    sys.exit(main())



