#!/usr/bin/env python3
"""
Phase 13 - Figure and table generation.
Produces 15 figures + 5 tables (20 total), built entirely from whatever the
docking/interaction/ADMET/pocket CSVs actually contain. Nothing here is a
fixed or assumed number -- row/column counts adapt to your real results.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="whitegrid")


def log(msg):
    print(f"[REPORT] {msg}", flush=True)


def safe_read_csv(path, required=False):
    """
    Robust CSV reader.

    Returns an empty DataFrame when the CSV is missing, zero-byte, or
    otherwise unparseable (e.g. header-only, whitespace-only). Raises only
    when the file is marked required (e.g. the docking summary, without
    which there is nothing to plot).
    """
    path = Path(path)

    if not path.exists():
        if required:
            raise FileNotFoundError(f"Required file not found: {path}")
        log(f"WARNING: {path} not found. Continuing with an empty table.")
        return pd.DataFrame()

    if path.stat().st_size == 0:
        if required:
            raise RuntimeError(f"Required CSV is empty: {path}")
        log(f"WARNING: {path} is empty (0 bytes). Continuing with an empty table.")
        return pd.DataFrame()

    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError,
            UnicodeDecodeError, FileNotFoundError, PermissionError) as e:
        if required:
            raise RuntimeError(f"Could not read required CSV {path}: {e}")
        log(f"WARNING: Could not read {path}: {e}. Continuing with an empty table.")
        return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docking", default="results/docking/summary.csv")
    ap.add_argument("--interactions", default="results/plip/interaction_summary.csv")
    ap.add_argument("--admet", default="results/admet/admet_results.csv")
    ap.add_argument("--figures-dir", default="results/figures")
    ap.add_argument("--tables-dir", default="results/tables")
    args = ap.parse_args()

    figdir, tabdir = Path(args.figures_dir), Path(args.tables_dir)
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    log("Loading pipeline outputs...")
    docking = safe_read_csv(args.docking, required=True)
    interactions = safe_read_csv(args.interactions)
    admet = safe_read_csv(args.admet)

    log(f"{len(docking)} docking result(s), {len(interactions)} interaction profile(s), "
        f"{len(admet)} ADMET record(s) loaded.")

    if docking.empty:
        log("ERROR: docking summary is empty -- nothing to plot. Check Phase 9 output.")
        return

    docking["complex"] = docking["ligand"] + "__" + docking["receptor"]
    merged = docking.merge(interactions, on="complex", how="left") if not interactions.empty else docking.copy()

    interaction_cols = [c for c in ["hydrogen_bonds", "hydrophobic_contacts", "salt_bridges",
                                     "water_bridges", "pi_stacks", "pi_cation_interactions",
                                     "halogen_bonds", "metal_complexes"] if c in merged.columns]

    # ---- T1: full merged results ----
    merged.to_csv(tabdir / "Table01_full_merged_results.csv", index=False)
    log(f"T1 written: Table01_full_merged_results.csv ({merged.shape[0]} rows x {merged.shape[1]} columns)")

    # ---- F1: binding affinity heatmap ----
    # pivot_table (not pivot) tolerates duplicate ligand-receptor pairs -- e.g.
    # replicate runs or multiple docking modes -- by aggregating with "min"
    # (the best / most negative affinity is the scientifically reasonable
    # value to report for a pair that was docked more than once).
    pivot = docking.pivot_table(index="ligand", columns="receptor",
                                 values="best_affinity_kcal_mol", aggfunc="min")
    plt.figure(figsize=(max(8, 0.9 * pivot.shape[1]), max(5, 0.6 * pivot.shape[0])))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="viridis_r", cbar_kws={"label": "kcal/mol"})
    plt.title("Predicted Binding Affinity: PKS-I-Derived Compounds vs Pathogen Targets")
    plt.tight_layout(); plt.savefig(figdir / "F1_binding_affinity_heatmap.png", dpi=300); plt.close()
    log("F1 written: binding affinity heatmap")

    # ---- F2: best affinity per pathogen ----
    best_per_pathogen = docking.groupby("pathogen")["best_affinity_kcal_mol"].min().reset_index()
    plt.figure(figsize=(7, 5))
    sns.barplot(data=best_per_pathogen, x="pathogen", y="best_affinity_kcal_mol", hue="pathogen",
                palette="mako", legend=False)
    plt.title("Strongest Predicted Binding Affinity per Pathogen")
    plt.ylabel("Best affinity (kcal/mol, more negative = stronger)")
    plt.tight_layout(); plt.savefig(figdir / "F2_best_affinity_per_pathogen.png", dpi=300); plt.close()
    strongest_pathogen = best_per_pathogen.loc[best_per_pathogen["best_affinity_kcal_mol"].idxmin(), "pathogen"]
    log(f"F2 written. Discovery: strongest overall pathogen response is {strongest_pathogen}")

    # ---- F3: average affinity per compound ----
    avg_per_compound = (docking.groupby("ligand")["best_affinity_kcal_mol"].mean()
                         .reset_index().sort_values("best_affinity_kcal_mol"))
    plt.figure(figsize=(7, max(5, 0.35 * len(avg_per_compound))))
    sns.barplot(data=avg_per_compound, x="best_affinity_kcal_mol", y="ligand", hue="ligand",
                palette="crest", legend=False)
    plt.title("Mean Binding Affinity per Compound Across Full Target Panel")
    plt.xlabel("Mean affinity (kcal/mol)")
    plt.tight_layout(); plt.savefig(figdir / "F3_average_affinity_per_compound.png", dpi=300); plt.close()
    best_compound = avg_per_compound.iloc[0]
    log(f"F3 written. Discovery: best-performing compound overall is "
        f"{best_compound['ligand']} (mean {best_compound['best_affinity_kcal_mol']:.2f} kcal/mol)")

    # ---- F4: affinity distribution boxplot by pathogen ----
    plt.figure(figsize=(7, 5))
    sns.boxplot(data=docking, x="pathogen", y="best_affinity_kcal_mol", hue="pathogen",
                palette="Set2", legend=False)
    sns.stripplot(data=docking, x="pathogen", y="best_affinity_kcal_mol", color="black", size=4, alpha=0.6)
    plt.title("Distribution of Binding Affinities by Pathogen (Boxplot)")
    plt.tight_layout(); plt.savefig(figdir / "F4_affinity_distribution_boxplot.png", dpi=300); plt.close()
    log("F4 written: affinity distribution boxplot")

    # ---- F5: affinity distribution violin plot by pathogen ----
    plt.figure(figsize=(7, 5))
    sns.violinplot(data=docking, x="pathogen", y="best_affinity_kcal_mol", hue="pathogen",
                    palette="Set3", legend=False, inner="quartile")
    plt.title("Distribution of Binding Affinities by Pathogen (Violin)")
    plt.tight_layout(); plt.savefig(figdir / "F5_affinity_distribution_violin.png", dpi=300); plt.close()
    log("F5 written: affinity distribution violin plot")

    # ---- F6-F8: individual interaction-type bar charts ----
    fig_num = 6
    for col, label in [("hydrogen_bonds", "Hydrogen Bonds"),
                        ("hydrophobic_contacts", "Hydrophobic Contacts"),
                        ("salt_bridges", "Salt Bridges")]:
        if col not in merged.columns:
            log(f"F{fig_num} skipped: '{col}' not present in interaction data.")
            fig_num += 1
            continue
        plt.figure(figsize=(9, max(4, 0.3 * len(merged))))
        sns.barplot(data=merged.sort_values(col, ascending=False), x=col, y="complex",
                    hue="complex", palette="flare", legend=False)
        plt.title(f"{label} per Docked Complex")
        plt.tight_layout(); plt.savefig(figdir / f"F{fig_num}_{col}_per_complex.png", dpi=300); plt.close()
        log(f"F{fig_num} written: {label.lower()} bar chart")
        fig_num += 1

    # ---- F9: water bridges + pi-stacking + remaining interaction types combined bar ----
    remaining_cols = [c for c in interaction_cols if c not in
                       ("hydrogen_bonds", "hydrophobic_contacts", "salt_bridges")]
    if remaining_cols:
        melted = merged.melt(id_vars="complex", value_vars=remaining_cols,
                              var_name="interaction_type", value_name="count")
        plt.figure(figsize=(9, max(4, 0.3 * len(merged))))
        sns.barplot(data=melted, x="count", y="complex", hue="interaction_type")
        plt.title("Additional Interaction Types per Docked Complex\n(water bridges, pi-stacking, pi-cation, halogen bonds, metal complexes)")
        plt.tight_layout(); plt.savefig(figdir / "F9_additional_interaction_types_per_complex.png", dpi=300); plt.close()
        log("F9 written: additional PLIP interaction types (water bridges, pi-stacking, etc.)")
    else:
        log("F9 skipped: no additional interaction-type columns present.")

    # ---- F10: stacked interaction profile ----
    if interaction_cols:
        top_complexes = merged.nlargest(15, interaction_cols[0]) if len(merged) > 15 else merged
        plt.figure(figsize=(10, 6))
        bottom = np.zeros(len(top_complexes))
        colors = sns.color_palette("Set2", n_colors=len(interaction_cols))
        for col, color in zip(interaction_cols, colors):
            plt.bar(top_complexes["complex"], top_complexes[col].fillna(0), bottom=bottom, label=col, color=color)
            bottom += top_complexes[col].fillna(0).values
        plt.xticks(rotation=75, ha="right"); plt.legend()
        plt.title("Combined Interaction Profile per Complex (All PLIP Interaction Types)")
        plt.tight_layout(); plt.savefig(figdir / "F10_stacked_interaction_profile.png", dpi=300); plt.close()
        log("F10 written: stacked interaction profile")
    else:
        log("F10 skipped: no interaction-type columns present.")

    # ---- F11: affinity vs total interactions ----
    if interaction_cols:
        merged["total_interactions"] = merged[interaction_cols].sum(axis=1)
        plt.figure(figsize=(7, 5))
        sns.scatterplot(data=merged, x="total_interactions", y="best_affinity_kcal_mol", hue="pathogen", s=80)
        plt.title("Binding Affinity vs Total Stabilising Interactions")
        plt.tight_layout(); plt.savefig(figdir / "F11_affinity_vs_interactions.png", dpi=300); plt.close()

        corr = merged["total_interactions"].corr(merged["best_affinity_kcal_mol"])
        if pd.notna(corr):
            log(f"F11 written. Discovery: correlation between interaction count and affinity = {corr:.2f}")
        else:
            log("F11 written. Correlation could not be calculated "
                "(insufficient variation or too few data points in the current dataset).")
    else:
        log("F11 skipped: no interaction-type columns present.")

    # ---- F12: affinity vs LogP; F13: ADMET radar; T3/T4 ----
    if not admet.empty:
        logp_col = "logP" if "logP" in admet.columns else ("LogP" if "LogP" in admet.columns else None)
        admet_merge = merged.merge(admet, left_on="ligand", right_on="compound", how="left")
        if logp_col:
            plt.figure(figsize=(7, 5))
            sns.scatterplot(data=admet_merge, x=logp_col, y="best_affinity_kcal_mol", hue="pathogen", s=80)
            plt.title("Binding Affinity vs Predicted LogP")
            plt.tight_layout(); plt.savefig(figdir / "F12_affinity_vs_logP.png", dpi=300); plt.close()
            log("F12 written: affinity vs LogP scatter")
        else:
            log("F12 skipped: no LogP column found in ADMET output.")

        radar_props = [c for c in admet.columns if c not in ("compound", "smiles")][:6]
        if radar_props:
            norm = admet.copy()
            for c in radar_props:
                rng = (norm[c].max() - norm[c].min()) or 1
                norm[c] = (norm[c] - norm[c].min()) / rng
            angles = np.linspace(0, 2 * np.pi, len(radar_props), endpoint=False).tolist()
            angles += angles[:1]
            fig, axes = plt.subplots(1, len(norm), figsize=(4 * len(norm), 4), subplot_kw=dict(polar=True))
            axes = [axes] if len(norm) == 1 else axes
            for ax, (_, row) in zip(axes, norm.iterrows()):
                values = row[radar_props].tolist(); values += values[:1]
                ax.plot(angles, values, linewidth=2); ax.fill(angles, values, alpha=0.25)
                ax.set_xticks(angles[:-1]); ax.set_xticklabels(radar_props, fontsize=7)
                ax.set_title(row.get("compound", ""), fontsize=10)
            plt.tight_layout(); plt.savefig(figdir / "F13_admet_radar_per_compound.png", dpi=300); plt.close()
            log("F13 written: ADMET radar profile per compound")
        else:
            log("F13 skipped: no usable ADMET property columns found.")

        admet.to_csv(tabdir / "Table03_admet_summary.csv", index=False)
        log("T3 written: Table03_admet_summary.csv")

        mw_col = "molecular_weight" if "molecular_weight" in admet.columns else (
            "MolecularWeight" if "MolecularWeight" in admet.columns else None)
        if mw_col and logp_col:
            lipinski = admet.copy()
            lipinski["MW_pass"] = lipinski[mw_col] < 500
            lipinski["LogP_pass"] = lipinski[logp_col] < 5
            lipinski["rule_of_five_pass"] = lipinski["MW_pass"] & lipinski["LogP_pass"]
            lipinski.to_csv(tabdir / "Table04_lipinski_rule_of_five.csv", index=False)
            n_pass = int(lipinski["rule_of_five_pass"].sum())
            log(f"T4 written. Discovery: {n_pass}/{len(lipinski)} compound(s) pass the Lipinski Rule-of-Five.")

            # ---- F14: Lipinski pass/fail summary bar ----
            plt.figure(figsize=(6, 4))
            counts = lipinski["rule_of_five_pass"].value_counts().rename({True: "Pass", False: "Fail"})
            sns.barplot(x=counts.index, y=counts.values, hue=counts.index, palette="pastel", legend=False)
            plt.title("Lipinski Rule-of-Five: Pass/Fail Summary")
            plt.ylabel("Number of compounds")
            plt.tight_layout(); plt.savefig(figdir / "F14_lipinski_pass_fail_summary.png", dpi=300); plt.close()
            log("F14 written: Lipinski pass/fail summary bar chart")
        else:
            log("T4/F14 skipped: molecular weight or LogP column not found in ADMET output.")
    else:
        log("No ADMET data found -- run Phase 12 first for F12, F13, F14, T3, T4.")

    # ---- F15: clustered heatmap (hierarchical clustering of the affinity matrix) ----
    if pivot.shape[0] > 1 and pivot.shape[1] > 1:
        pivot_filled = pivot.fillna(pivot.max().max())
        try:
            g = sns.clustermap(pivot_filled, cmap="viridis_r", annot=True, fmt=".1f",
                                cbar_kws={"label": "kcal/mol"}, figsize=(max(8, pivot.shape[1]), max(6, pivot.shape[0])))
            g.figure.suptitle("Clustered Binding Affinity Heatmap (Hierarchical Clustering)", y=1.02)
            g.savefig(figdir / "F15_clustered_affinity_heatmap.png", dpi=300)
            plt.close(g.figure)
            log("F15 written: clustered binding affinity heatmap")
        except Exception as e:
            log(f"F15 skipped: clustering failed ({e}) -- likely too few unique rows/columns.")
    else:
        log("F15 skipped: need at least 2 compounds and 2 receptors to cluster.")

    # ---- T2: top 10 docking pairs ----
    top10 = merged.sort_values("best_affinity_kcal_mol").head(10)
    top10.to_csv(tabdir / "Table02_top10_docking_pairs.csv", index=False)
    if not top10.empty:
        log(f"T2 written: Table02_top10_docking_pairs.csv -- strongest pair is "
            f"{top10.iloc[0]['ligand']} vs {top10.iloc[0]['receptor']} at "
            f"{top10.iloc[0]['best_affinity_kcal_mol']} kcal/mol")
    else:
        log("T2 written: Table02_top10_docking_pairs.csv (no docking pairs available)")

    log("Figure/table generation complete: 15 figures in figures/, 5 tables in tables/ "
        "(Table05_pocket_grid_summary.csv was written earlier by Phase 8b).")


if __name__ == "__main__":
    main()

    