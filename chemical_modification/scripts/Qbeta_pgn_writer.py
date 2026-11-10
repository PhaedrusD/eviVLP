"""
Qbeta_pgn_writer.py
Writes QB_capsid.pgn (180 subunits: 60 A + 60 B + 60 C) to the current
working directory.  The GUI launches this script from eviVLP/GUI/, so all
segment PDB files (QB_A_com*.pdb, QB_B_com*.pdb, QB_C_com*.pdb) produced by
Qbeta_capsid_maker.py are expected there too.
"""
import os


def nterm_slot():
    """'; first NONE; mutate 1 NPC' when this run built an N-terminus residue.
    Mirrors FileProcessor::ProcessPgnFile so the monomer and capsid agree."""
    return "; first NONE; mutate 1 NPC" if os.path.exists("ntrm_clean.str") else ""


def number_to_letters(n):
    """Convert 0 → 'A', 25 → 'Z', 26 → 'BA', … (reverse-alphabetic, like Excel cols)."""
    result = []
    while n >= 0:
        result.append(chr(ord('A') + n % 26))
        n = n // 26 - 1
        if n < 0:
            break
    return ''.join(reversed(result))


def main():
    slot = nterm_slot()
    print(f"N-terminus slot: '{slot}'")

    lines = [
        "# vmd -dispdev text -e QB_capsid.pgn",
        "package require psfgen",
        "topology top_all36_prot.rtf",
        "topology top_all36_cgenff_CBD.rtf",
        "pdbalias residue HIS HSE",
        "pdbalias atom ILE CD1 CD",
        "",
    ]

    for group_offset, prefix in [(0, "QB_A"), (60, "QB_B"), (120, "QB_C")]:
        for i in range(60):
            seg = number_to_letters(group_offset + i)
            pdb = f"{prefix}_com{i}.pdb"
            # Unmodified chains take psfgen's default NTER patch (NH3+ at
            # residue 1), matching the monomer template and the TMV stream.
            # When residue 1 is modified (NPC), the slot supplies
            # "; first NONE; mutate 1 NPC" so no terminal patch is applied.
            lines.append(f"segment {seg} {{pdb {pdb}{slot}}}")
            lines.append(f"coordpdb {pdb} {seg}")

    lines += [
        "",
        "guesscoord",
        "regenerate angles dihedrals",
        "guesscoord",
        "writepdb capsid_naa.pdb",
        "writepsf capsid_naa.psf",
        # Solvate and ionise in the same run, as the TMV stream does, so that
        # Minimize Capsid always gets a water box built from THIS capsid. The
        # capsid spans ~300 A, so "-t 10" gives the ~320 A cube that used to be
        # set by hand in Minimize_GUI_solvate.tcl; capsid_cell.py reads the
        # exact cell back from the CRYST1 record autoionize writes.
        "package require solvate",
        "solvate capsid_naa.psf capsid_naa.pdb -t 10 -o QB_all_water",
        "package require autoionize",
        "autoionize -psf QB_all_water.psf -pdb QB_all_water.pdb -sc 0.2 -o QB_naa_capsid_wb",
        "exit",
    ]

    with open("QB_capsid.pgn", "w") as f:
        f.write("\n".join(lines) + "\n")

    print("QB_capsid.pgn written successfully.")


if __name__ == "__main__":
    main()
