"""
TMV_pgn_writer.py
Writes TMV_all_out.pgn (49 subunits: TMV0.pdb … TMV48.pdb) to the current
working directory.  The GUI launches this script from eviVLP/GUI/.
"""
import os


def nterm_slot():
    """'; first NONE; mutate 1 NPC' when this run built an N-terminus residue.
    Mirrors FileProcessor::ProcessPgnFile so the monomer and capsid agree."""
    return "; first NONE; mutate 1 NPC" if os.path.exists("ntrm_clean.str") else ""


def main():
    slot = nterm_slot()
    print(f"N-terminus slot: '{slot}'")

    lines = [
        "# vmd -dispdev text -e TMV_all_out.pgn",
        "package require psfgen",
        "topology top_all36_prot.rtf",
        "topology top_all36_cgenff_CBD.rtf",
        "pdbalias residue HIS HSE",
        "pdbalias atom ILE CD1 CD",
        "",
    ]

    for i in range(49):
        pdb = f"TMV{i}.pdb"
        lines.append(f"segment U{i} {{pdb {pdb}{slot}}}")
        lines.append(f"coordpdb {pdb} U{i}")

    lines += [
        "",
        "guesscoord",
        "regenerate angles dihedrals",
        "guesscoord",
        "writepdb TMV_naa_capsid.pdb",
        "writepsf TMV_naa_capsid.psf",
        "package require solvate",
        "solvate TMV_naa_capsid.psf TMV_naa_capsid.pdb -t 10 -o TMV_all_water",
        "package require autoionize",
        "autoionize -psf TMV_all_water.psf -pdb TMV_all_water.pdb -sc 0.2 -o TMV_naa_capsid_wb",
        "exit",
    ]

    with open("TMV_all_out.pgn", "w") as f:
        f.write("\n".join(lines) + "\n")

    print("TMV_all_out.pgn written successfully.")


if __name__ == "__main__":
    main()
