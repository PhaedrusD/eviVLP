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
            # The Qbeta capsid has always been built with "first NONE" (no NTER
            # patch); the NPC slot already carries its own "first NONE".
            seg_opts = slot if slot else "; first NONE"
            lines.append(f"segment {seg} {{pdb {pdb}{seg_opts}}}")
            lines.append(f"coordpdb {pdb} {seg}")

    lines += [
        "",
        "guesscoord",
        "regenerate angles dihedrals",
        "guesscoord",
        "writepdb capsid_naa.pdb",
        "writepsf capsid_naa.psf",
        "exit",
    ]

    with open("QB_capsid.pgn", "w") as f:
        f.write("\n".join(lines) + "\n")

    print("QB_capsid.pgn written successfully.")


if __name__ == "__main__":
    main()
