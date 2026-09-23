#!/usr/bin/env python3
"""
Phase 1 - Reference genome acquisition.

WHY THIS PHASE EXISTS:
Because the isolate's own genome is not yet sequenced, PKS-I mining is done
on a representative, fully-sequenced reference genome of the presumptive
genus/species instead. This script accepts EITHER a genome accession
(if you already know it) OR an organism name (e.g. "Bacillus velezensis"),
in which case it queries NCBI itself to find a suitable reference/representative
genome assembly -- so you do not have to manually look up an accession.

REAL API USED: NCBI Entrez esearch/esummary against the "assembly" database.
This is a standard, documented NCBI E-utilities pattern.

USAGE:
    python3 01_fetch_genome.py --organism "Bacillus velezensis" --outdir results/genomes
    python3 01_fetch_genome.py --accession GCF_000063585.1 --outdir results/genomes
"""
import argparse
import gzip
import shutil
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from pks_dock.reproducibility import DecisionLog  # noqa: E402

Entrez.email = "kizitlabs.research@gmail.com"  # required by NCBI E-utilities; replace with your own contact email
Entrez.tool = "pksdock_pipeline"

_DECISIONS = DecisionLog()


def log(msg):
    print(f"[PHASE 1] {msg}", flush=True)


def find_assembly_by_organism(organism: str, retries: int = 3):
    """
    Searches NCBI assembly DB for a reference or representative genome
    matching the given organism name. Falls back to 'complete genome'
    latest assemblies if no reference/representative flag is found.
    Returns the assembly accession (e.g. GCF_...) and its FTP directory.
    """
    log(f"Searching NCBI Assembly database for organism: '{organism}'")

    queries = [
        f'"{organism}"[Organism] AND "reference genome"[filter]',
        f'"{organism}"[Organism] AND "representative genome"[filter]',
        f'"{organism}"[Organism] AND "complete genome"[assembly_level] AND "latest"[filter]',
        f'"{organism}"[Organism]',
    ]

    for q in queries:
        log(f"Query attempt: {q}")
        for attempt in range(retries):
            try:
                handle = Entrez.esearch(db="assembly", term=q, retmax=5)
                record = Entrez.read(handle)
                handle.close()
                break
            except Exception as e:
                log(f"NCBI request failed ({e}), retrying ({attempt + 1}/{retries})...")
                time.sleep(2 * (attempt + 1))
        else:
            continue

        ids = record.get("IdList", [])
        if not ids:
            log("No hits for this query tier, trying a broader query...")
            continue

        log(f"Found {len(ids)} candidate assembly record(s). Fetching summaries to pick the best one...")
        handle = Entrez.esummary(db="assembly", id=",".join(ids), report="full")
        summary = Entrez.read(handle, validate=False)
        handle.close()

        docs = summary["DocumentSummarySet"]["DocumentSummary"]
        # Prefer RefSeq category "reference genome" > "representative genome" > anything else,
        # then prefer the assembly with the highest contig N50 (more complete/contiguous).
        def rank(doc):
            cat = doc.get("RefSeq_category", "na")
            n50 = float(doc.get("ContigN50", 0) or 0)
            cat_score = {"reference genome": 2, "representative genome": 1}.get(cat, 0)
            return (cat_score, n50)

        docs_sorted = sorted(docs, key=rank, reverse=True)
        best = docs_sorted[0]
        accession = best["AssemblyAccession"]
        ftp_path = best.get("FtpPath_RefSeq") or best.get("FtpPath_GenBank")
        organism_name = best.get("Organism", organism)

        log(f"Selected assembly: {accession} ({organism_name}) "
            f"[RefSeq_category={best.get('RefSeq_category', 'na')}, ContigN50={best.get('ContigN50', 'na')}]")
        if not ftp_path:
            log("WARNING: no FTP path found for this assembly, trying next candidate query tier.")
            continue
        _DECISIONS.record(
            "PHASE 1", f"Reference genome selection for '{organism}'",
            database="NCBI Assembly", query=q, selected=accession,
            alternatives_considered=[d["AssemblyAccession"] for d in docs_sorted[1:4]],
            reason=f"RefSeq_category={best.get('RefSeq_category', 'na')}, highest ContigN50 among candidates",
            confidence=best.get("RefSeq_category", "na"),
        )
        return accession, ftp_path

    raise RuntimeError(
        f"Could not find any assembly for organism '{organism}' via NCBI Assembly search. "
        f"Try supplying --accession directly if you already know the genome you want to use."
    )


def download_genome_from_ftp(ftp_path: str, accession: str, outdir: Path, retries: int = 3) -> Path:
    """Downloads the *_genomic.gbff.gz file from an NCBI assembly FTP directory and unzips it."""
    base_name = ftp_path.rstrip("/").split("/")[-1]
    https_path = ftp_path.replace("ftp://", "https://")
    gbff_url = f"{https_path}/{base_name}_genomic.gbff.gz"
    log(f"Downloading GenBank flat file: {gbff_url}")

    gz_path = outdir / f"{accession}_genomic.gbff.gz"
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(gbff_url, stream=True, timeout=120)
            r.raise_for_status()
            with open(gz_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
            break
        except requests.RequestException as e:
            last_exc = e
            if attempt < retries:
                sleep_s = 2 * attempt
                log(f"Download attempt {attempt}/{retries} failed ({e}); retrying in {sleep_s}s...")
                time.sleep(sleep_s)
            else:
                raise RuntimeError(f"Failed to download {gbff_url} after {retries} attempts: {last_exc}")
    log(f"Downloaded {gz_path.stat().st_size:,} bytes -> {gz_path}")

    gbk_path = outdir / f"{accession}.gbk"
    with gzip.open(gz_path, "rb") as f_in, open(gbk_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    log(f"Decompressed -> {gbk_path}")
    return gbk_path


def download_genome_by_accession(accession: str, outdir: Path) -> Path:
    """Direct nucleotide fetch when the user already supplied an accession (e.g. a single-record RefSeq accession)."""
    log(f"Fetching {accession} directly from NCBI nucleotide database...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    gbk_path = outdir / f"{accession}.gbk"
    SeqIO.write(record, str(gbk_path), "genbank")
    log(f"Download complete: {accession} is {len(record.seq):,} bp, "
        f"organism = {record.annotations.get('organism', 'unknown')}")
    return gbk_path


def main():
    ap = argparse.ArgumentParser(description="Phase 1: acquire a reference genome by organism name or accession")
    ap.add_argument("--organism", help='e.g. "Bacillus velezensis" (searched automatically via NCBI Assembly)')
    ap.add_argument("--accession", help="A specific GenBank/RefSeq accession, if already known")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if not args.organism and not args.accession:
        sys.exit("[PHASE 1] ERROR: supply either --organism or --accession")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log("Reason: the isolate's own genome is not yet sequenced, so a documented reference")
    log("genome of the presumptive genus/species is mined instead, as a stated, transparent substitute.")

    if args.accession:
        gbk_path = download_genome_by_accession(args.accession, outdir)
        accession_used = args.accession
    else:
        accession_used, ftp_path = find_assembly_by_organism(args.organism)
        gbk_path = download_genome_from_ftp(ftp_path, accession_used, outdir)

    manifest = outdir / "genome_manifest.txt"
    with open(manifest, "w") as f:
        f.write(f"accession={accession_used}\n")
        f.write(f"gbk_path={gbk_path}\n")
    log(f"Wrote manifest -> {manifest} (Phase 2 reads this automatically, no manual accession retyping needed)")
    log("Phase 1 complete.")


if __name__ == "__main__":
    main()
