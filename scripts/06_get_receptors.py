#!/usr/bin/env python3
"""
Phase 6 - Receptor retrieval, with automated quality-based selection among
candidate PDB IDs (not a single hardcoded ID picked blindly).

For every target listed in config/pathogen_targets.yaml:
  - If it has candidate_pdb_ids, query the RCSB Data API
    (https://data.rcsb.org/rest/v1/core/entry/{id}, a real, documented REST
    endpoint) for resolution and experimental method for each candidate, and
    automatically pick the best one (lowest resolution number = highest
    quality; prefer X-RAY DIFFRACTION over other methods when tied).
  - If predict_if_empty is true (no experimental structure exists), it
    prepares the FASTA for Phase 6b (ColabFold) instead of downloading.

USAGE:
    python3 06_get_receptors.py --pathogens Staphylococcus_aureus Bacillus_cereus Pseudomonas_aeruginosa \
        --config config/pathogen_targets.yaml --outdir results/receptors_raw
"""
import argparse
import json
import sys
from pathlib import Path

import requests
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from pks_dock.net import HTTPClient  # noqa: E402
from pks_dock.reproducibility import DecisionLog  # noqa: E402

RCSB_ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry"
RCSB_DOWNLOAD = "https://files.rcsb.org/download"


def log(msg):
    print(f"[PHASE 6] {msg}", flush=True)


def get_entry_quality(client: HTTPClient, pdb_id: str):
    """Queries RCSB Data API (with retry+cache) for resolution + experimental method."""
    try:
        data = client.get_json(f"{RCSB_ENTRY_API}/{pdb_id}")
        methods = data.get("exptl", [{}])
        method = methods[0].get("method", "UNKNOWN") if methods else "UNKNOWN"
        resolution = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])
        resolution = resolution[0] if resolution else None
        return {"pdb_id": pdb_id, "method": method, "resolution": resolution}
    except Exception as e:
        log(f"  Could not fetch RCSB metadata for {pdb_id}: {e}")
        return {"pdb_id": pdb_id, "method": "UNKNOWN", "resolution": None}


def pick_best_candidate(candidates: list):
    """Ranks candidate PDB IDs: prefer X-ray with lowest (best) resolution; None resolution ranks last."""
    def rank(c):
        method_score = 1 if c["method"] == "X-RAY DIFFRACTION" else 0
        res = c["resolution"] if c["resolution"] is not None else 999.0
        return (method_score, -res)
    return sorted(candidates, key=rank, reverse=True)[0]


def download_pdb(client: HTTPClient, pdb_id: str, outpath: Path):
    content = client.get_binary(f"{RCSB_DOWNLOAD}/{pdb_id}.pdb")
    outpath.write_bytes(content)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pathogens", nargs="+", required=True,
                    help="Keys from config/pathogen_targets.yaml, e.g. Staphylococcus_aureus")
    ap.add_argument("--config", default="config/pathogen_targets.yaml")
    ap.add_argument("--outdir", default="results/receptors_raw")
    ap.add_argument("--decision-log", default="results/reports/decision_log.json")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    client = HTTPClient(phase="PHASE 6")
    decisions = DecisionLog(args.decision_log)

    manifest = {}
    to_predict = []

    for pathogen_key in args.pathogens:
        pdata = cfg["pathogens"].get(pathogen_key)
        if pdata is None:
            log(f"WARNING: '{pathogen_key}' not found in {args.config} -- skipping. "
                f"Add a block for it to config/pathogen_targets.yaml if you want it included.")
            continue

        tag = pdata["organism_tag"]
        log(f"Processing target panel for {pathogen_key} (tag: {tag})")

        for target in pdata["targets"]:
            tname = f"{target['name']}_{tag}"

            if target.get("predict_if_empty") and not target.get("candidate_pdb_ids"):
                log(f"  {tname}: no experimental structure available -- queued for ColabFold prediction "
                    f"(Phase 6b) using {target.get('fasta_path', 'MISSING FASTA PATH -- add one to config')}")
                to_predict.append({"name": tname, "fasta_path": target.get("fasta_path")})
                continue

            candidates = target.get("candidate_pdb_ids", [])
            if not candidates:
                log(f"  {tname}: WARNING -- no candidate_pdb_ids and predict_if_empty is not set. "
                    f"Fix config/pathogen_targets.yaml for this target.")
                continue

            log(f"  {tname}: evaluating {len(candidates)} candidate PDB structure(s): {candidates}")
            quality_info = [get_entry_quality(client, pid) for pid in candidates]
            for qi in quality_info:
                log(f"    {qi['pdb_id']}: method={qi['method']}, resolution={qi['resolution']}")
            best = pick_best_candidate(quality_info)
            log(f"  {tname}: SELECTED {best['pdb_id']} "
                f"(method={best['method']}, resolution={best['resolution']}) as highest quality.")
            decisions.record(
                "PHASE 6", f"Receptor structure selection for {tname}",
                database="RCSB PDB", query=str(candidates), selected=best["pdb_id"],
                alternatives_considered=[c["pdb_id"] for c in quality_info if c["pdb_id"] != best["pdb_id"]],
                reason=f"method={best['method']}, resolution={best['resolution']} "
                       f"(lowest resolution / X-ray preferred among candidates)",
                confidence=best["method"],
            )

            pdb_path = outdir / f"{tname}.pdb"
            download_pdb(client, best["pdb_id"], pdb_path)
            log(f"  Downloaded -> {pdb_path}")

            manifest[tname] = {
                "pathogen": pathogen_key, "pdb_id": best["pdb_id"],
                "resolution": best["resolution"], "method": best["method"],
                "pdb_path": str(pdb_path), "source": "experimental",
            }

    manifest_path = outdir / "receptor_manifest.json"
    manifest_path.write_text(json.dumps({"resolved": manifest, "to_predict": to_predict}, indent=2))
    log(f"Wrote {manifest_path} -- {len(manifest)} experimental receptor(s) ready, "
        f"{len(to_predict)} queued for ColabFold prediction (run 06b_predict_receptors.sh next).")


if __name__ == "__main__":
    main()
