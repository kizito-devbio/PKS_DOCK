#!/usr/bin/env python3
"""
Phase 3 - Extract Type-I PKS clusters from antiSMASH output.

HONESTY NOTE ON THIS SCRIPT (please read before trusting it blindly):
antiSMASH's exact JSON schema has changed across major versions (v5/v6/v7),
and I cannot guarantee from here which exact key names your installed
version will produce. Rather than hard-assume one schema and silently
produce wrong results if it doesn't match (the mistake in the previous
draft), this script:
  1. Tries the GenBank output first (region/proto_cluster features with a
     /product qualifier) -- this format has been far more stable across
     antiSMASH versions than the JSON internals.
  2. Falls back to the KnownClusterBlast text tables (knownclusterblast/*.txt)
     for known-cluster similarity, which are plain, documented text tables.
  3. Prints exactly what it found and which method worked, so if your
     antiSMASH version's output differs, you will SEE that in the terminal
     rather than get silently wrong numbers.

If neither method finds anything but you know antiSMASH detected clusters,
run this once with --debug to dump the raw structure it saw, and adjust the
key names accordingly -- do not just trust a silent "0 clusters found".
"""
import argparse
import glob
import json
import sys
from pathlib import Path

from Bio import SeqIO

PKS1_TAGS = ("T1PKS", "transAT-PKS", "T1PKS-like")


def log(msg):
    print(f"[PARSE] {msg}", flush=True)


def parse_from_genbank(antismash_dir: Path, debug=False):
    """Primary method: read the antiSMASH-annotated .gbk/.gbff and pull region features."""
    gbk_files = list(antismash_dir.glob("*.gbk")) + list(antismash_dir.glob("*.gbff"))
    if not gbk_files:
        log("No antiSMASH-annotated GenBank file found in output dir.")
        return []

    hits = []
    for gbk in gbk_files:
        log(f"Scanning annotated GenBank: {gbk.name}")
        for record in SeqIO.parse(str(gbk), "genbank"):
            for feature in record.features:
                if feature.type not in ("region", "proto_cluster", "cand_cluster"):
                    continue
                products = feature.qualifiers.get("product", [])
                if not any(any(tag in p for tag in PKS1_TAGS) for p in products):
                    continue
                region_number = feature.qualifiers.get("region_number", ["?"])[0]
                contig_edge = feature.qualifiers.get("contig_edge", ["?"])[0]
                hits.append({
                    "contig": record.id,
                    "region_number": region_number,
                    "products_predicted": products,
                    "location": str(feature.location),
                    "contig_edge": contig_edge,
                    "predicted_compound": "unannotated PKS-I cluster",   # filled in by knownclusterblast lookup below
                    "most_similar_known_cluster": "none / unknown cluster",
                    "similarity_percent": "n/a",
                })
                if debug:
                    log(f"DEBUG feature qualifiers: {dict(feature.qualifiers)}")
    return hits


def enrich_with_knownclusterblast(antismash_dir: Path, hits: list):
    """Reads knownclusterblast/regionX_cY.txt tables (stable plain-text format) for
    the top known-cluster similarity hit per region, and attaches it to the matching hit."""
    kcb_dir = antismash_dir / "knownclusterblast"
    if not kcb_dir.exists():
        log("No knownclusterblast/ directory found -- similarity/compound-name annotation will stay 'unknown'.")
        return hits

    txt_files = sorted(kcb_dir.glob("*.txt"))
    log(f"Found {len(txt_files)} knownclusterblast table(s) to cross-reference.")

    for hit in hits:
        region_num = str(hit["region_number"])
        match_file = None
        for tf in txt_files:
            if f"_c{region_num}" in tf.stem or f"region{int(region_num):03d}" in tf.stem.lower() if region_num.isdigit() else False:
                match_file = tf
                break
        if match_file is None:
            # fall back to loose filename match
            for tf in txt_files:
                if region_num in tf.stem:
                    match_file = tf
                    break
        if match_file is None:
            continue

        lines = match_file.read_text(errors="ignore").splitlines()
        # KnownClusterBlast text tables list "1. BGC0000123_c1  compound_name  similarity%" style rows
        # under a header like "Significant hits:" -- take the first data row after that header.
        capture = False
        for line in lines:
            if line.strip().lower().startswith("significant hits"):
                capture = True
                continue
            if capture and line.strip() and line.strip()[0].isdigit():
                parts = line.strip().split("\t") if "\t" in line else line.strip().split("  ")
                parts = [p for p in parts if p.strip()]
                if len(parts) >= 2:
                    hit["most_similar_known_cluster"] = parts[1].strip() if len(parts) > 1 else parts[0]
                    hit["predicted_compound"] = parts[1].strip() if len(parts) > 1 else "unannotated PKS-I cluster"
                break

    return hits


def parse_from_json_fallback(antismash_dir: Path, debug=False):
    """Secondary/legacy fallback for older antiSMASH JSON layouts. Best-effort only."""
    json_files = glob.glob(str(antismash_dir / "*.json"))
    log(f"[fallback] Found {len(json_files)} antiSMASH summary JSON file(s), attempting best-effort parse...")
    hits = []
    for jf in json_files:
        try:
            data = json.load(open(jf))
        except Exception as e:
            log(f"[fallback] Could not parse {jf}: {e}")
            continue
        for record in data.get("records", []):
            for region in record.get("areas", record.get("regions", [])):
                products = region.get("products", [])
                if not any(any(tag in p for tag in PKS1_TAGS) for p in products):
                    continue
                knowncluster = region.get("knowncluster", {})
                hits.append({
                    "contig": record.get("id", "unknown"),
                    "region_number": region.get("region_number", "?"),
                    "products_predicted": products,
                    "location": region.get("location", "n/a"),
                    "contig_edge": "n/a",
                    "predicted_compound": knowncluster.get("description", "unannotated PKS-I cluster"),
                    "most_similar_known_cluster": knowncluster.get("mibig_id", "none / unknown cluster"),
                    "similarity_percent": knowncluster.get("similarity", "n/a"),
                })
        if debug:
            log(f"[fallback] DEBUG top-level keys in {jf}: {list(data.keys())}")
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--antismash-dir", required=True, help="Or pass results/antismash/latest_antismash_dir.txt content")
    ap.add_argument("--out", default="results/antismash/phase3_compound_list.json")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    antismash_dir = Path(args.antismash_dir)
    log(f"Scanning antiSMASH output in: {antismash_dir}")

    hits = parse_from_genbank(antismash_dir, debug=args.debug)
    method = "genbank"
    if hits:
        hits = enrich_with_knownclusterblast(antismash_dir, hits)
    else:
        log("GenBank-feature method found nothing -- trying JSON fallback...")
        hits = parse_from_json_fallback(antismash_dir, debug=args.debug)
        method = "json_fallback"

    log(f"Discovery: {len(hits)} Type-I PKS cluster(s) detected (method: {method})")

    for i, h in enumerate(hits, 1):
        sim = h["similarity_percent"]
        log(f"Cluster {i} (region {h['region_number']}, contig {h['contig']}): "
            f"predicted product = {h['predicted_compound']} "
            f"(most similar known cluster: {h['most_similar_known_cluster']}, similarity: {sim})")
        if isinstance(sim, (int, float)) and sim < 50:
            log(f"  NOTE: similarity below 50% -- treat compound identity as TENTATIVE ONLY. "
                f"State this explicitly in the thesis text for cluster {i}.")
        if h["predicted_compound"] == "unannotated PKS-I cluster":
            log(f"  NOTE: no KnownClusterBlast match -- this cluster's product is NOT confidently named. "
                f"Do not treat 'unannotated PKS-I cluster' as a real compound name downstream.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(hits, f, indent=2)
    log(f"Wrote {len(hits)} cluster record(s) -> {out_path} for Phase 4 to consume automatically.")

    if not hits:
        log("WARNING: zero PKS-I clusters detected. Either the genome genuinely has none, "
            "or the parser didn't match this antiSMASH version's output format -- "
            "re-run with --debug and inspect the printed structure before concluding 'no PKS-I genes'.")
        sys.exit(2)


if __name__ == "__main__":
    main()
