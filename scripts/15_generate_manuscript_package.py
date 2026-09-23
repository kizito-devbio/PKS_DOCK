#!/usr/bin/env python3
"""
Phase 15 - Manuscript / supplementary package export.
Bundles all tables into one multi-sheet Excel workbook (the format journals'
supplementary-materials portals almost always want), copies the figures at
publication DPI, and writes a supplementary-methods Markdown stub you can
paste into (or attach alongside) the manuscript's Methods/Supplementary
Information section.
"""
import argparse
import shutil
from pathlib import Path

import pandas as pd


def log(msg):
    print(f"[PHASE 15] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-dir", default="results/tables")
    ap.add_argument("--figures-dir", default="results/figures")
    ap.add_argument("--out-dir", default="results/reports/manuscript_package")
    args = ap.parse_args()

    tabdir, figdir = Path(args.tables_dir), Path(args.figures_dir)
    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Consolidated supplementary table workbook
    workbook_path = outdir / "Supplementary_Tables.xlsx"
    table_files = sorted(tabdir.glob("Table*.csv"))
    if not table_files:
        log("No table CSVs found -- run Phase 13 (and Phase 8b for Table05) first.")
        return

    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        for tf in table_files:
            sheet_name = tf.stem[:31]  # Excel sheet name limit
            df = pd.read_csv(tf)
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            log(f"  Added sheet '{sheet_name}' ({df.shape[0]} rows x {df.shape[1]} cols)")
    log(f"Wrote consolidated supplementary workbook -> {workbook_path}")

    # 2. Copy all figures into the package folder at their generated 300 DPI
    fig_out = outdir / "figures"
    fig_out.mkdir(exist_ok=True)
    n_figs = 0
    for f in sorted(figdir.glob("F*.png")):
        shutil.copy(f, fig_out / f.name)
        n_figs += 1
    log(f"Copied {n_figs} figure(s) (300 DPI, as generated) -> {fig_out}")

    # 3. Supplementary methods stub
    methods_path = outdir / "Supplementary_Methods_stub.md"
    methods_path.write_text(
        "## Supplementary Methods (auto-generated stub -- review and edit before submission)\n\n"
        "In silico analysis was performed using the open-source PKS_DOCK pipeline "
        "(https://github.com/<your-username>/PKS_DOCK, version X.X). Briefly: reference genome "
        "acquisition and PKS-I biosynthetic gene cluster mining were performed with antiSMASH; "
        "predicted compound structures were retrieved automatically from PubChem, the NCI Chemical "
        "Identifier Resolver, and ChEBI; receptor structures were retrieved from the RCSB Protein "
        "Data Bank with automated resolution-based quality selection among literature-validated "
        "candidate structures, or predicted with ColabFold where no experimental structure existed; "
        "binding pockets were detected with fpocket and the top-ranked pocket by druggability score "
        "was used to define the docking grid; molecular docking was performed with AutoDock Vina; "
        "non-covalent interactions were profiled with PLIP; and ADMET properties were predicted with "
        "admet-ai. All software versions are pinned in environment.yml in the repository above. "
        "[EDIT: insert your specific parameter values, genome accession(s) used, and any deviations "
        "from the default pipeline configuration here.]"
    )
    log(f"Wrote supplementary methods stub -> {methods_path}")
    log(f"Manuscript package complete: {outdir}")


if __name__ == "__main__":
    main()
