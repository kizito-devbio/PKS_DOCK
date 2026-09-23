#!/usr/bin/env python3
"""
Phase 10 - PLIP interaction profiling

Purpose:
    Convert docking results into PLIP-compatible complexes and
    perform protein-ligand interaction profiling.

Input:
    results/docking/*.pdbqt          (docked ligand poses; best pose used)
    results/receptors_pdbqt/*_fixed.pdb   (preferred receptor, if present)
    results/receptors_raw/*.pdb      (fallback receptor)
    results/ligands/*.sdf            (accepted for CLI compatibility; not
                                       used to build the complex anymore)

Output:
    results/plip/
        ligand__receptor/
            complex.pdb
            ligand.pdb
            complex_report.xml
            PLIP logs

Usage (unchanged, backward compatible):

python scripts/10_run_plip.py \
    results/docking \
    results/receptors_raw \
    results/ligands \
    results/plip

Optional new flags:

    --prepared-receptor-dir DIR
        Explicit directory containing <receptor>_fixed.pdb files.
        If omitted, defaults to <docking_dir's parent>/receptors_pdbqt
        for backward compatibility with the existing project layout,
        but this default is only a guess and is logged as such.

    --strip-solvent-hetatm
        If set, removes common crystallographic solvent/buffer HETATM
        records (waters, sulfate, glycerol, PEG, etc.) from the receptor
        before building the complex. Off by default: by default ALL
        receptor HETATM records (including solvent/buffer molecules)
        are kept, exactly as in the original script. Metals and
        cofactors (e.g. ZN, HEM, NAG) are never stripped by this flag,
        even when it is enabled.
"""


import argparse
import csv
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="[PHASE 10] %(message)s"
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------
# Deliberate list of common non-cofactor solvent/buffer
# HETATM residue names. This list is intentionally narrow:
# it only ever removes crystallization additives, never
# metals or cofactors (ZN, HEM, NAG, FAD, NAD, etc. are
# always kept regardless of this flag).
# ---------------------------------------------------------

SOLVENT_HETATM_RESIDUES = {
    "HOH", "WAT", "DOD",           # water
    "SO4", "PO4",                  # sulfate / phosphate
    "GOL", "EDO", "PEG", "PG4",    # glycerol / PEG cryoprotectants
    "TRS", "ACT", "MPD", "DMS",    # buffer / cryo components
    "BME", "FMT", "IPA", "1PE"
    # Deliberately NOT included: CL, NA, K, IOD (and other simple
    # ions). These are frequently biologically or structurally
    # relevant (e.g. stabilizing binding sites, catalysis) rather
    # than pure crystallization contaminants, so they are always
    # kept even when --strip-solvent-hetatm is enabled.
}


# ---------------------------------------------------------
# Dependency check
# ---------------------------------------------------------

def check_program(program):

    if shutil.which(program) is None:

        raise RuntimeError(
            f"{program} is not available in PATH"
        )



# ---------------------------------------------------------
# Run command
# ---------------------------------------------------------

def run_command(cmd, logfile=None):

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )


    if logfile:

        with open(logfile, "w") as f:

            f.write("COMMAND:\n")
            f.write(" ".join(cmd))

            f.write("\n\nSTDOUT:\n")
            f.write(result.stdout)

            f.write("\n\nSTDERR:\n")
            f.write(result.stderr)



    if result.returncode != 0:

        raise RuntimeError(
            result.stderr.strip()
        )


    return result



# ---------------------------------------------------------
# Extract MODEL 1 (best Vina pose) directly from the PDBQT
# text, with NO dependence on obabel's internal molecule
# indexing. This guarantees the block we convert is exactly
# the one Vina labelled "MODEL 1", not whatever obabel
# happens to treat as its first entry.
# ---------------------------------------------------------

def extract_model_1_block(docking_file):

    lines = docking_file.read_text().splitlines()

    start = None
    end = None

    for i, line in enumerate(lines):

        stripped = line.strip()

        if start is None and stripped.startswith("MODEL"):

            # MODEL records look like "MODEL     1" or "MODEL 1"
            parts = stripped.split()

            if len(parts) >= 2 and parts[1] == "1":

                start = i

                continue


        if start is not None and stripped.startswith("ENDMDL"):

            end = i

            break


    if start is None:

        # No MODEL records at all: some single-pose pdbqt files
        # (e.g. rigid ligand, single conformer) omit MODEL/ENDMDL
        # entirely. In that case the whole file IS the best pose.
        logger.info(
            f"No MODEL records found in {docking_file.name}; "
            "treating entire file as the single pose."
        )

        return "\n".join(lines)


    if end is None:

        raise RuntimeError(
            f"Found MODEL 1 in {docking_file.name} but no matching "
            "ENDMDL. File may be truncated or malformed."
        )


    # include the MODEL/ENDMDL wrapper lines themselves; obabel
    # handles this fine and it keeps the block self-describing
    block = lines[start:end + 1]

    return "\n".join(block)



# ---------------------------------------------------------
# Extract best pose (Model 1) and convert to PDB
# ---------------------------------------------------------

def extract_best_pose(
        docking_file,
        pdb
):
    """
    Extract the top-ranked Vina pose (MODEL 1) from a docking
    output .pdbqt and convert it to a .pdb file, preserving the
    DOCKED coordinates (not the original ligand geometry).

    MODEL 1 is isolated in pure Python first (see
    extract_model_1_block), then that isolated single-model
    PDBQT is handed to obabel for conversion. This avoids any
    reliance on obabel's own multi-molecule indexing lining up
    with Vina's MODEL numbering.
    """

    logger.info(
        f"Extracting best pose (Model 1): {docking_file.name}"
    )


    model_1_text = extract_model_1_block(docking_file)


    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".pdbqt",
        delete=False
    ) as tmp:

        tmp.write(model_1_text)
        tmp.write("\n")

        tmp_path = Path(tmp.name)


    try:

        cmd = [

            "obabel",

            str(tmp_path),

            "-O",

            str(pdb)
        ]


        run_command(cmd)


    finally:

        if tmp_path.exists():

            tmp_path.unlink()



# ---------------------------------------------------------
# Resolve receptor file (prefer fixed pdbqt-derived PDB)
# ---------------------------------------------------------

def resolve_receptor(
        receptor_name,
        receptor_dir,
        prepared_receptor_dir
):
    """
    Prefer <prepared_receptor_dir>/<receptor>_fixed.pdb when present.
    Fall back to <receptor_dir>/<receptor>.pdb for backward
    compatibility.

    prepared_receptor_dir is explicit: either passed in via
    --prepared-receptor-dir, or a logged best-guess default. This
    function never silently assumes a project layout.
    """

    preferred = (
        prepared_receptor_dir /
        f"{receptor_name}_fixed.pdb"
    )

    if preferred.exists():

        logger.info(
            f"Using preferred receptor: {preferred}"
        )

        return preferred


    fallback = (
        receptor_dir /
        f"{receptor_name}.pdb"
    )

    logger.info(
        f"Preferred receptor not found ({preferred}); "
        f"falling back to {fallback}"
    )

    return fallback



# ---------------------------------------------------------
# Filter receptor HETATM lines (deliberate, documented)
# ---------------------------------------------------------

def filter_receptor_lines(
        lines,
        strip_solvent_hetatm
):
    """
    Deliberately decide what receptor HETATM records make it into
    the complex.

    strip_solvent_hetatm=False (default): keep everything, exactly
        as the original script did (ATOM + all HETATM).

    strip_solvent_hetatm=True: drop only residues in
        SOLVENT_HETATM_RESIDUES (waters, sulfate/phosphate,
        cryoprotectants, common crystallization ions). Metals and
        cofactors are always kept.

    Logs which residue names were kept vs removed either way, so
    the choice is always visible in the run log rather than silent.
    """

    kept = []

    removed_residues = set()
    kept_hetatm_residues = set()


    for line in lines:

        if line.startswith("ATOM"):

            kept.append(line)

            continue


        if line.startswith("HETATM"):

            resname = line[17:20].strip()

            if strip_solvent_hetatm and resname in SOLVENT_HETATM_RESIDUES:

                removed_residues.add(resname)

                continue


            kept_hetatm_residues.add(resname)

            kept.append(line)


    if kept_hetatm_residues:

        logger.info(
            f"Receptor HETATM residues kept: {sorted(kept_hetatm_residues)}"
        )


    if removed_residues:

        logger.info(
            f"Receptor HETATM residues removed (solvent/buffer): "
            f"{sorted(removed_residues)}"
        )


    return kept



# ---------------------------------------------------------
# Build complex PDB
# ---------------------------------------------------------

def build_complex(
        receptor,
        ligand,
        output,
        strip_solvent_hetatm
):

    logger.info(
        f"Building complex: {output.name}"
    )


    with open(receptor) as rec:

        receptor_lines = filter_receptor_lines(
            rec.readlines(),
            strip_solvent_hetatm
        )


    with open(output, "w") as out:


        # receptor

        for line in receptor_lines:

            out.write(
                line
                if line.endswith("\n")
                else line + "\n"
            )


        # TER before ligand: many structural tools (ChimeraX,
        # PyMOL) expect an explicit chain break between receptor
        # and ligand records. PLIP tolerates its absence, but this
        # keeps the file well-formed for downstream visualization.

        out.write(
            "TER\n"
        )



        # ligand

        with open(ligand) as lig:

            for line in lig:

                if line.startswith(
                    (
                        "ATOM",
                        "HETATM"
                    )
                ):

                    out.write(line)



        out.write(
            "END\n"
        )



# ---------------------------------------------------------
# Run PLIP
# ---------------------------------------------------------

def run_plip(
        complex_file,
        output_dir
):


    logger.info(
        f"Running PLIP: {complex_file.parent.name}"
    )


    logfile = (
        output_dir /
        "plip_execution.log"
    )


    cmd = [

        "plip",

        "-f",

        str(complex_file),

        "--xml",

        "-o",

        str(output_dir)

    ]


    run_command(
        cmd,
        logfile
    )



    # PLIP versions differ:
    # complex_report.xml
    # *_report.xml

    reports = []

    for pattern in [
        "*_report.xml",
        "*report.xml",
        "*.xml"
    ]:

        reports.extend(
            output_dir.glob(pattern)
        )


    reports = list(
        set(reports)
    )


    if not reports:

        raise RuntimeError(
            "PLIP completed but no XML report found. "
            f"Check {logfile}"
        )


    logger.info(
        f"Found PLIP report: {reports[0].name}"
    )


    return reports[0]



# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():


    parser = argparse.ArgumentParser(
        description="Phase 10 PLIP analysis"
    )


    parser.add_argument(
        "docking_dir"
    )

    parser.add_argument(
        "receptor_dir"
    )

    parser.add_argument(
        "ligand_dir"
    )

    parser.add_argument(
        "output_dir"
    )

    parser.add_argument(
        "--prepared-receptor-dir",
        dest="prepared_receptor_dir",
        default=None,
        help=(
            "Directory containing <receptor>_fixed.pdb files. "
            "Defaults to <docking_dir's parent>/receptors_pdbqt if "
            "not given, for backward compatibility (logged as a guess)."
        )
    )

    parser.add_argument(
        "--strip-solvent-hetatm",
        dest="strip_solvent_hetatm",
        action="store_true",
        default=False,
        help=(
            "Remove common solvent/buffer HETATM residues (water, "
            "sulfate, glycerol, PEG, common ions) from the receptor "
            "before building the complex. Off by default: all "
            "receptor HETATM records are kept, matching original "
            "behavior. Metals and cofactors are never removed."
        )
    )


    args = parser.parse_args()



    docking_dir = Path(
        args.docking_dir
    )

    receptor_dir = Path(
        args.receptor_dir
    )

    ligand_dir = Path(
        args.ligand_dir
    )

    output_dir = Path(
        args.output_dir
    )


    if args.prepared_receptor_dir:

        prepared_receptor_dir = Path(
            args.prepared_receptor_dir
        )

        logger.info(
            f"Using explicit prepared-receptor-dir: {prepared_receptor_dir}"
        )


    else:

        prepared_receptor_dir = (
            docking_dir.parent /
            "receptors_pdbqt"
        )

        logger.info(
            "--prepared-receptor-dir not given; defaulting (guess) to "
            f"{prepared_receptor_dir}. Pass --prepared-receptor-dir "
            "explicitly if your project layout differs."
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    logger.info(
        f"ligand_dir={ligand_dir} accepted for CLI compatibility; "
        "ligand coordinates are now sourced from docked .pdbqt poses, "
        "not from ligand_dir."
    )

    logger.info(
        f"strip_solvent_hetatm={args.strip_solvent_hetatm}"
    )



    check_program("obabel")
    check_program("plip")



    docking_files = sorted(
        docking_dir.glob(
            "*.pdbqt"
        )
    )


    if not docking_files:

        raise RuntimeError(
            "No docking PDBQT files found"
        )



    logger.info(
        f"Found {len(docking_files)} docking complexes"
    )



    summary = []

    failures = []



    for docking_file in docking_files:


        name = docking_file.stem


        logger.info(
            "----------------------------------------"
        )

        logger.info(
            f"Processing {name}"
        )


        try:


            if "__" not in name:

                raise RuntimeError(
                    "Expected ligand__receptor filename"
                )


            ligand_name, receptor_name = (
                name.split(
                    "__",
                    1
                )
            )


            receptor_pdb = resolve_receptor(
                receptor_name,
                receptor_dir,
                prepared_receptor_dir
            )



            if not docking_file.exists():

                raise FileNotFoundError(
                    docking_file
                )


            if not receptor_pdb.exists():

                raise FileNotFoundError(
                    receptor_pdb
                )



            plip_dir = (
                output_dir /
                name
            )


            plip_dir.mkdir(
                parents=True,
                exist_ok=True
            )



            ligand_pdb = (
                plip_dir /
                f"{ligand_name}.pdb"
            )


            complex_pdb = (
                plip_dir /
                "complex.pdb"
            )



            extract_best_pose(
                docking_file,
                ligand_pdb
            )


            build_complex(
                receptor_pdb,
                ligand_pdb,
                complex_pdb,
                args.strip_solvent_hetatm
            )


            run_plip(
                complex_pdb,
                plip_dir
            )


            summary.append(
                [
                    name,
                    "SUCCESS"
                ]
            )


            logger.info(
                f"Completed: {name}"
            )



        except Exception as e:


            logger.error(
                f"FAILED {name}: {e}"
            )


            failures.append(
                [
                    name,
                    str(e)
                ]
            )


            summary.append(
                [
                    name,
                    "FAILED"
                ]
            )



    # summary

    with open(
        output_dir /
        "plip_summary.csv",
        "w",
        newline=""
    ) as f:


        writer = csv.writer(f)

        writer.writerow(
            [
                "complex",
                "status"
            ]
        )

        writer.writerows(
            summary
        )



    if failures:

        with open(
            output_dir /
            "failed_plip.csv",
            "w",
            newline=""
        ) as f:

            writer = csv.writer(f)

            writer.writerow(
                [
                    "complex",
                    "reason"
                ]
            )

            writer.writerows(
                failures
            )



    success = sum(
        1
        for row in summary
        if row[1]=="SUCCESS"
    )


    logger.info(
        "========================================"
    )

    logger.info(
        f"PLIP completed: {success}/{len(summary)} successful"
    )

    logger.info(
        f"Results: {output_dir}"
    )



if __name__ == "__main__":

    main()

    