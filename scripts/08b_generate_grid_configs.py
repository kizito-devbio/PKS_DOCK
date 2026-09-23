#!/usr/bin/env python3
"""
Phase 8b - Parse fpocket's info.txt, rank pockets by druggability score, and
derive an AutoDock Vina grid box (center + size) from the SELECTED pocket's
actual 3D coordinates (from its pocketN_vert.pqr pseudo-atoms, i.e. the real
detected cavity -- not a fixed 40x40x40 box, and not a box centered on a
ligand that no longer exists in the file).

CORRECTED:
  1. Druggability sort key no longer treats a real score of 0.0 as missing
     (previously `p["druggability"] or -1` collapsed 0.0 -> -1).
  2. Only fpocket output directories whose receptor actually has a matching
     .pdbqt in --receptors-pdbqt-dir (i.e. passed Phase 7) are processed.
     Stale/orphaned fpocket output for a receptor that failed Phase 7 is
     logged and skipped instead of silently producing a grid.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import yaml


def log(msg):
    print(f"[PHASE 8b] {msg}", flush=True)


def parse_fpocket_info(info_path: Path):
    """
    Parses fpocket's <name>_info.txt. Format is blocks like:
        Pocket 1 :
            Score :                     0.5231
            Druggability Score :        0.789
            ...
    Returns a list of dicts: [{"pocket_num": 1, "score": .., "druggability": ..}, ...]
    """
    text = info_path.read_text()
    blocks = re.split(r"\n(?=Pocket\s+\d+\s*:)", text)
    pockets = []
    for block in blocks:
        m = re.match(r"Pocket\s+(\d+)\s*:", block.strip())
        if not m:
            continue
        pocket_num = int(m.group(1))
        score_m = re.search(r"^\s*Score\s*:\s*([\-\d.]+)", block, re.MULTILINE)
        drug_m = re.search(r"Druggability Score\s*:\s*([\-\d.]+)", block)
        volume_m = re.search(r"^\s*Volume\s*:\s*([\-\d.]+)", block, re.MULTILINE)
        pockets.append({
            "pocket_num": pocket_num,
            "score": float(score_m.group(1)) if score_m else None,
            "druggability": float(drug_m.group(1)) if drug_m else None,
            "volume": float(volume_m.group(1)) if volume_m else None,
        })
    return pockets


def pocket_coordinates(fpocket_dir: Path, pocket_num: int):
    """Reads pocketN_vert.pqr (the alpha-sphere pseudo-atoms marking the cavity) and
    returns their XYZ coordinates. Falls back to pocketN_atm.pdb (lining residues) if needed."""
    vert_path = fpocket_dir / "pockets" / f"pocket{pocket_num}_vert.pqr"
    atm_path = fpocket_dir / "pockets" / f"pocket{pocket_num}_atm.pdb"

    coords = []
    target = vert_path if vert_path.exists() else atm_path
    if not target.exists():
        return None

    is_pqr = target.suffix == ".pqr"
    for line in target.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if is_pqr:
            # PQR format: the last 5 whitespace-separated tokens are always
            # x, y, z, charge, radius -- robust regardless of how the earlier
            # serial/atomname/resname columns are padded (verified against
            # sample fpocket output; fixed-column PDB-style slicing is NOT
            # safe here because PQR column widths vary by writer).
            tokens = line.split()
            try:
                x, y, z = (float(t) for t in tokens[-5:-2])
                coords.append((x, y, z))
            except (ValueError, IndexError):
                continue
        else:
            # Standard fixed-column PDB format for the _atm.pdb fallback
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                continue
    return np.array(coords) if coords else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fpocket-dir", required=True, help="results/fpocket, containing one subdir per receptor")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default="config/pathogen_targets.yaml")
    ap.add_argument(
        "--receptors-pdbqt-dir",
        default="results/receptors_pdbqt",
        help="Directory of *.pdbqt files -- the authority list of receptors that passed Phase 7. "
             "fpocket output for any receptor NOT in this list is skipped.",
    )
    args = ap.parse_args()

    docking_params = yaml.safe_load(Path(args.config).read_text()).get("docking_parameters", {})
    min_drug_score = docking_params.get("fpocket_min_druggability_score", 0.0)
    padding = docking_params.get("grid_padding_angstrom", 8)

    fpocket_root = Path(args.fpocket_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdbqt_dir = Path(args.receptors_pdbqt_dir)
    allowed = {x.stem for x in pdbqt_dir.glob("*.pdbqt")}
    if not allowed:
        log(f"WARNING: no *.pdbqt files found in {pdbqt_dir} -- nothing is considered "
            f"to have passed Phase 7, so no grids will be generated.")

    grid_summary = []

    for receptor_dir in sorted(fpocket_root.iterdir()):
        if not receptor_dir.is_dir():
            continue
        name = receptor_dir.name

        if name not in allowed:
            log(f"{name}: fpocket output exists but receptor is not in {pdbqt_dir} "
                f"(did not pass Phase 7). Skipping.")
            continue

        info_files = list(receptor_dir.glob("*_info.txt"))
        if not info_files:
            log(f"{name}: no fpocket info.txt found -- skipping.")
            continue

        pockets = parse_fpocket_info(info_files[0])
        if not pockets:
            log(f"{name}: fpocket ran but no pockets were parsed from {info_files[0].name}.")
            continue

        log(f"{name}: {len(pockets)} pocket(s) detected. Ranking by druggability score...")
        pockets_ranked = sorted(
            pockets,
            key=lambda p: (p["druggability"] if p["druggability"] is not None else -1),
            reverse=True,
        )
        for p in pockets_ranked[:3]:
            log(f"  Pocket {p['pocket_num']}: druggability={p['druggability']}, "
                f"score={p['score']}, volume={p['volume']}")

        best = next(
            (p for p in pockets_ranked if (p["druggability"] if p["druggability"] is not None else -1) >= min_drug_score),
            pockets_ranked[0],
        )
        log(f"{name}: SELECTED pocket {best['pocket_num']} "
            f"(druggability={best['druggability']}) as the docking site.")

        coords = pocket_coordinates(receptor_dir, best["pocket_num"])
        if coords is None:
            log(f"{name}: could not read coordinates for pocket {best['pocket_num']} -- skipping grid generation.")
            continue

        center = coords.mean(axis=0)
        box_min, box_max = coords.min(axis=0), coords.max(axis=0)
        size = (box_max - box_min) + 2 * padding
        size = np.maximum(size, 20.0)  # sane minimum box size regardless of a very tight pocket

        log(f"{name}: grid center = ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}), "
            f"size = ({size[0]:.1f}, {size[1]:.1f}, {size[2]:.1f}) A "
            f"(derived from detected pocket + {padding} A padding)")

        config_path = outdir / f"{name}.box.txt"
        config_path.write_text(
            f"center_x = {center[0]:.3f}\n"
            f"center_y = {center[1]:.3f}\n"
            f"center_z = {center[2]:.3f}\n"
            f"size_x = {size[0]:.3f}\n"
            f"size_y = {size[1]:.3f}\n"
            f"size_z = {size[2]:.3f}\n"
        )
        log(f"{name}: wrote grid config -> {config_path}")

        grid_summary.append({
            "receptor": name, "pocket_num": best["pocket_num"], "druggability": best["druggability"],
            "score": best["score"], "volume": best["volume"],
            "center_x": center[0], "center_y": center[1], "center_z": center[2],
            "size_x": size[0], "size_y": size[1], "size_z": size[2],
        })

    import pandas as pd
    summary_path = Path("results/tables/Table05_pocket_grid_summary.csv")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(grid_summary).to_csv(summary_path, index=False)
    log(f"Wrote {summary_path} ({len(grid_summary)} receptor(s)) -- this is Table 5 in the figure/table set.")


if __name__ == "__main__":
    main()

    