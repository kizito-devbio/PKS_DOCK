#!/usr/bin/env python3
"""
Phase 12 - ADMET prediction.
Reads results/ligands/phase4_resolved_compounds.json (produced by Phase 4).
Every SMILES used here was fetched automatically from PubChem/NCI-CIR/ChEBI --
nothing is hand-typed.

NOTE: this script's default --resolved-compounds path was already correct
(results/ligands/phase4_resolved_compounds.json). The FileNotFoundError seen
in the run log was caused entirely by 04_get_ligands.py writing its manifest
to results/phase4_resolved_compounds.json instead -- one directory level up
from where this script (correctly) looks. That has been fixed in
04_get_ligands.py; no change was needed here. Kept identical for the record.
"""
import argparse
import json
from pathlib import Path

import pandas as pd


def log(msg):
    print(f"[PHASE 12] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolved-compounds", default="results/ligands/phase4_resolved_compounds.json")
    ap.add_argument("--out", default="results/admet/admet_results.csv")
    args = ap.parse_args()

    log("Loading compound structures resolved automatically in Phase 4...")
    data = json.loads(Path(args.resolved_compounds).read_text())
    resolved = data["resolved"]

    if not resolved:
        log("No resolved compounds to run ADMET on. Exiting.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame().to_csv(args.out, index=False)
        return

    from admet_ai import ADMETModel
    model = ADMETModel()

    results = []
    for name, info in resolved.items():
        smiles = info.get("smiles")
        if not smiles:
            log(f"Skipping '{name}' -- no SMILES available (source was {info['source']}, which "
                f"did not return a SMILES string directly in this pipeline's lookup chain). "
                f"If you need ADMET for this compound, retrieve its SMILES manually and add it "
                f"to {args.resolved_compounds}.")
            continue

        log(f"Predicting ADMET profile for '{name}' (SMILES: {smiles})...")
        preds = model.predict(smiles=smiles)
        preds["compound"] = name
        preds["smiles"] = smiles
        results.append(preds)

        logp = preds.get("logP", preds.get("LogP", "n/a"))
        hep = preds.get("hepatotoxicity", preds.get("DILI", "n/a"))
        ames = preds.get("AMES_toxicity", preds.get("AMES", "n/a"))
        log(f"  LogP={logp}, hepatotoxicity={hep}, AMES mutagenicity={ames}")

    df = pd.DataFrame(results)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    log(f"ADMET prediction complete for {len(df)} compound(s). Saved to {args.out}")


if __name__ == "__main__":
    main()