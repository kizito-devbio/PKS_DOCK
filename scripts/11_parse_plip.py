#!/usr/bin/env python3
"""
Phase 11 - Parse PLIP XML reports.

Parses every PLIP XML report generated during Phase 10 and summarizes the
number of detected interactions for each docked complex.

Supports all interaction types currently produced by PLIP and automatically
discovers XML report filenames (report.xml, complex_report.xml, etc.) without
hardcoding a specific filename.
"""

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

import pandas as pd


INTERACTION_TAGS = {
    "hydrogen_bonds": "hydrogen_bond",
    "hydrophobic_contacts": "hydrophobic_interaction",
    "salt_bridges": "salt_bridge",
    "water_bridges": "water_bridge",
    "pi_stacks": "pi_stack",
    "pi_cation_interactions": "pi_cation_interaction",
    "halogen_bonds": "halogen_bond",
    "metal_complexes": "metal_complex",
}


def log(message):
    print(f"[PHASE 11] {message}", flush=True)


def find_xml_reports(plip_dir: Path):
    """
    Recursively locate every XML report produced by PLIP.

    Works regardless of filename:
        report.xml
        complex_report.xml
        xxxx.xml
    """
    reports = []

    for xml in plip_dir.rglob("*.xml"):
        try:
            root = ET.parse(xml).getroot()

            if root.tag.lower() == "report":
                reports.append(xml)

        except Exception:
            continue

    return sorted(reports)


def count_interactions(root):

    counts = {}

    for column, tag in INTERACTION_TAGS.items():
        counts[column] = len(root.findall(f".//{tag}"))

    return counts


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--plip-dir",
        default="results/plip",
        help="Directory containing PLIP results."
    )

    parser.add_argument(
        "--out",
        default="results/plip/interaction_summary.csv",
        help="Output CSV."
    )

    args = parser.parse_args()

    plip_dir = Path(args.plip_dir)

    log("Parsing PLIP XML reports...")

    xml_reports = find_xml_reports(plip_dir)

    if not xml_reports:
        log("WARNING: No PLIP XML reports were found.")

    rows = []

    for xml in xml_reports:

        complex_name = xml.parent.name

        try:
            root = ET.parse(xml).getroot()

        except ET.ParseError as e:

            log(f"WARNING: Could not parse {xml}: {e}")

            continue

        row = {

            "complex": complex_name

        }

        row.update(count_interactions(root))

        rows.append(row)

        total = sum(row[k] for k in INTERACTION_TAGS)

        detail = ", ".join(
            f"{row[k]} {k.replace('_', ' ')}"
            for k in INTERACTION_TAGS
            if row[k] > 0
        )

        log(
            f"{complex_name}: {total} interaction(s)"
            + (f" ({detail})" if detail else "")
        )

    columns = ["complex"] + list(INTERACTION_TAGS.keys())

    df = pd.DataFrame(rows, columns=columns)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(args.out, index=False)

    if not df.empty:

        interaction_cols = list(INTERACTION_TAGS.keys())

        df["total_interactions"] = df[interaction_cols].sum(axis=1)

        best = df.loc[df["total_interactions"].idxmax()]

        log("")

        log("Summary")

        log("--------------------------------")

        log(f"Complexes parsed : {len(df)}")

        log(f"Total interactions: {df['total_interactions'].sum()}")

        log(
            f"Most interactions : "
            f"{best['complex']} "
            f"({best['total_interactions']})"
        )

    else:

        log("No interaction records were extracted.")

    log(f"Interaction summary written to {args.out}")


if __name__ == "__main__":
    main()
    