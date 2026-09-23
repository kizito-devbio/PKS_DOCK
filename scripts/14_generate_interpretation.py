#!/usr/bin/env python3
"""
Phase 14 - Auto-generated interpretation report.

FIX FROM THE PREVIOUS DRAFT: "interpretation" was only scattered terminal
print statements. This script turns the actual merged results + ADMET +
pocket data into a real, readable prose Markdown report -- every sentence is
built from a real number in your results, not a template placeholder.
"""
import argparse
from pathlib import Path

import pandas as pd


def log(msg):
    print(f"[PHASE 14] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-dir", default="results/tables")
    ap.add_argument("--out", default="results/reports/interpretation.md")
    args = ap.parse_args()

    tabdir = Path(args.tables_dir)
    t1 = pd.read_csv(tabdir / "Table01_full_merged_results.csv")
    t2 = pd.read_csv(tabdir / "Table02_top10_docking_pairs.csv")
    t3_path = tabdir / "Table03_admet_summary.csv"
    t4_path = tabdir / "Table04_lipinski_rule_of_five.csv"
    t5_path = tabdir / "Table05_pocket_grid_summary.csv"

    lines = []
    lines.append("# Automated Interpretation Report\n")
    lines.append("_Generated directly from pipeline output tables. Every statement below is "
                  "computed from your actual results -- nothing here is templated filler._\n")

    lines.append("## Overall docking summary\n")
    n_pairs = len(t1)
    n_compounds = t1["ligand"].nunique()
    n_receptors = t1["receptor"].nunique()
    n_pathogens = t1["pathogen"].nunique() if "pathogen" in t1.columns else None
    lines.append(f"- {n_pairs} ligand-receptor docking pair(s) were evaluated, "
                 f"covering {n_compounds} predicted compound(s) against {n_receptors} receptor target(s)"
                 + (f" across {n_pathogens} pathogen(s)." if n_pathogens else "."))

    best_row = t1.loc[t1["best_affinity_kcal_mol"].idxmin()]
    lines.append(f"- The strongest predicted interaction overall was **{best_row['ligand']}** "
                 f"against **{best_row['receptor']}**, at {best_row['best_affinity_kcal_mol']} kcal/mol.")

    if "pathogen" in t1.columns:
        per_pathogen_best = t1.groupby("pathogen")["best_affinity_kcal_mol"].min().sort_values()
        lines.append("- Strongest binding per pathogen:")
        for pathogen, val in per_pathogen_best.items():
            lines.append(f"  - {pathogen}: {val:.2f} kcal/mol")

    interaction_cols = [c for c in ["hydrogen_bonds", "hydrophobic_contacts", "salt_bridges",
                                     "water_bridges", "pi_stacks", "pi_cation_interactions",
                                     "halogen_bonds", "metal_complexes"] if c in t1.columns]
    if interaction_cols:
        t1["_total_interactions"] = t1[interaction_cols].sum(axis=1)
        top_int = t1.loc[t1["_total_interactions"].idxmax()]
        lines.append(f"- The most extensively stabilised complex by non-covalent contact count was "
                     f"**{top_int['complex']}**, with {int(top_int['_total_interactions'])} total interactions "
                     f"across {len(interaction_cols)} PLIP-detected interaction types.")

    lines.append("\n## Top 10 ranked pairs\n")
    lines.append(t2[["ligand", "receptor", "pathogen", "best_affinity_kcal_mol"]].to_markdown(index=False))

    if t3_path.exists() and t3_path.stat().st_size > 0:
        admet = pd.read_csv(t3_path)
        if not admet.empty:
            lines.append("\n## ADMET summary\n")
            lines.append(f"- ADMET profiles were predicted for {len(admet)} compound(s).")
            if t4_path.exists():
                lipinski = pd.read_csv(t4_path)
                n_pass = int(lipinski["rule_of_five_pass"].sum()) if "rule_of_five_pass" in lipinski.columns else "n/a"
                lines.append(f"- {n_pass}/{len(lipinski)} compound(s) pass the Lipinski Rule-of-Five "
                             f"drug-likeness check (MW < 500, LogP < 5).")
                if isinstance(n_pass, int) and n_pass < len(lipinski):
                    failing = lipinski.loc[~lipinski["rule_of_five_pass"], "compound"].tolist() \
                        if "compound" in lipinski.columns else []
                    if failing:
                        lines.append(f"- Compound(s) failing the Rule-of-Five: {', '.join(failing)}. "
                                     f"State this explicitly alongside any binding-affinity claims for these "
                                     f"compounds in the discussion chapter, since drug-likeness concerns temper "
                                     f"their translational relevance even if binding affinity is strong.")

    if t5_path.exists() and t5_path.stat().st_size > 0:
        pocket_df = pd.read_csv(t5_path)
        if not pocket_df.empty:
            lines.append("\n## Pocket detection summary (fpocket)\n")
            lines.append(f"- Binding pockets were detected and grid boxes derived for "
                         f"{len(pocket_df)} receptor(s).")
            low_drug = pocket_df[pocket_df["druggability"] < 0.5] if "druggability" in pocket_df.columns else pd.DataFrame()
            if not low_drug.empty:
                names = ", ".join(low_drug["receptor"].tolist())
                lines.append(f"- **Caution:** the selected pocket for {names} had a druggability score "
                             f"below 0.5 -- interpret docking results for these receptor(s) with more caution, "
                             f"and consider inspecting alternative pockets manually.")

    lines.append("\n## Honesty statement for the thesis\n")
    lines.append("Because the isolates have not yet been sequenced, this entire computational chapter is "
                 "presumptive and reference-genome-based. State plainly wherever this analysis appears: "
                 "the biosynthetic gene cluster, predicted metabolite, and downstream docking/ADMET results "
                 "were derived from a representative reference genome of the presumptive genus, and remain "
                 "to be confirmed once whole-genome or targeted PKS-I sequencing of the isolates themselves "
                 "is completed.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines))
    log(f"Interpretation report written -> {out_path}")


if __name__ == "__main__":
    main()
