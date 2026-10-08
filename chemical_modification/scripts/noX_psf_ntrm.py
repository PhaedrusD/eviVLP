#!/usr/bin/env python3
"""
noX_psf_ntrm.py
Remove cap atoms (names containing 'X', e.g. HX) from the N-terminal NAA
stream file, together with every BOND line that references them.

Atoms are matched BY NAME, not by position. The earlier index-based version
silently deleted the wrong atom (drug C1 instead of HX) whenever the PDB and
the CGenFF .str had different atom counts (e.g. extra HN2/HN3 on the N).

Before writing anything, the non-cap atom names in the PDB and the .str are
compared; if they differ the script stops with exit code 1 so the GUI halts.
"""
import sys


def pdb_atom_names(pdb_file):
    names = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                names.append(line[12:16].strip())
    return names


def str_atom_names(str_file):
    names = []
    with open(str_file) as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 2:
                    names.append(parts[1])
    return names


def is_cap(name):
    return 'X' in name


def process_str_files(pdb_file='ntrm_fixed.pdb', input_file='ntrm.str',
                      output_file='ntrm_clean.str'):
    pdb_names = pdb_atom_names(pdb_file)
    str_names = str_atom_names(input_file)

    # --- consistency check -------------------------------------------------
    pdb_core = {n for n in pdb_names if not is_cap(n)}
    str_core = {n for n in str_names if not is_cap(n)}
    only_pdb = sorted(pdb_core - str_core)
    only_str = sorted(str_core - pdb_core)
    if only_pdb or only_str:
        print(f"ERROR: atom names in {pdb_file} and {input_file} do not match.")
        if only_pdb:
            print(f"  only in {pdb_file}: {only_pdb}")
        if only_str:
            print(f"  only in {input_file}: {only_str}")
        print("  The .str must come from the same structure as the PDB "
              "(re-run CGenFF on the current ntrm.mol2).")
        sys.exit(1)

    cap_names = {n for n in str_names if is_cap(n)}
    if not cap_names:
        print(f"Warning: no cap (X) atoms found in {input_file}")
    print(f"Removing cap atoms by name: {sorted(cap_names)}")

    # --- filter ------------------------------------------------------------
    removed = written = 0
    with open(input_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('RESI') and 'naa.pdb' in line:
                line = line.replace('naa.pdb', 'NAA')
            parts = line.split()
            if line.startswith('ATOM') and len(parts) >= 2 and parts[1] in cap_names:
                removed += 1
                continue
            if line.startswith('BOND') and any(p in cap_names for p in parts[1:]):
                removed += 1
                continue
            fout.write(line)
            written += 1

    print(f"Lines removed: {removed}; lines written to {output_file}: {written}")
    print(f"Atoms: {len(str_names)} -> {len(str_names) - len(cap_names)}")


if __name__ == "__main__":
    process_str_files()
