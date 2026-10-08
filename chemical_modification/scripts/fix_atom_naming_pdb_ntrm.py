#!/usr/bin/env python3
"""
fix_atom_naming_pdb_ntrm.py  (N-terminus / NPC path)

Post-noX step, mirror of fix_atom_naming_pdb_naa.py:
  * ntrm_clean.str : RESI header normalised to 'NPC'
  * ntrm_clean.pdb : HETATM->ATOM, resname NPC, resid 1, atom names copied
                     by position from ntrm_clean.str (CGenFF may rename atoms)
Both files are the ones IC_ntrm_NPC reads next, so they must agree.
Exits 1 on an atom-count mismatch so the GUI stops.
"""
import sys

STR_FILE = 'ntrm_clean.str'
PDB_FILE = 'ntrm_clean.pdb'
RESNAME = 'NPC'


def extract_atom_info(str_path):
    with open(str_path) as f:
        lines = f.readlines()
    out, names = [], []
    for line in lines:
        if line.startswith('RESI'):
            parts = line.split()
            if len(parts) >= 2 and parts[1] != RESNAME:
                line = line.replace(parts[1], RESNAME, 1)
        out.append(line)
        if line.startswith('ATOM'):
            p = line.split()
            if len(p) >= 2 and 'X' not in p[1]:
                names.append(p[1])
    with open(str_path, 'w') as f:
        f.writelines(out)
    return names


def update_pdb_file(pdb_path, names):
    with open(pdb_path) as f:
        lines = f.readlines()
    out, n = [], 0
    for line in lines:
        if line.startswith('COMPND') and 'UNNAMED' in line:
            line = line.replace('UNNAMED', RESNAME)
        if line.startswith('HETATM'):
            line = 'ATOM  ' + line[6:]
        if line.startswith('ATOM'):
            line = line[:17] + RESNAME + line[20:22] + '   1' + line[26:]
            if n < len(names):
                line = line[:12] + f"{names[n]:<4}" + line[16:]
            n += 1
        out.append(line)
    with open(pdb_path, 'w') as f:
        f.writelines(out)
    return n


def main():
    names = extract_atom_info(STR_FILE)
    n_pdb = update_pdb_file(PDB_FILE, names)
    if n_pdb != len(names):
        print(f"ERROR: {STR_FILE} has {len(names)} atoms but {PDB_FILE} has {n_pdb}.")
        sys.exit(1)
    print(f"Update successful. {n_pdb} atoms named from {STR_FILE} into {PDB_FILE}.")


if __name__ == '__main__':
    main()
