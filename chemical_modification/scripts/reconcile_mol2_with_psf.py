#!/usr/bin/env python3
"""reconcile_mol2_with_psf.py  in.mol2  patched.psf  out.mol2

PSF is the source of truth for the molecule: which ATOMS exist and which BONDS
exist. obabel only supplies coordinates + perceived bond ORDERS.
  - mol2 atoms beyond the PSF set (e.g. valence-fill H) are dropped.
  - mol2 bonds absent from the PSF are dropped (spurious).
  - PSF bonds obabel missed are added (order 1).
Indices line up: psfgen writes PSF and PDB in the same atom order and plain
obabel conversion preserves it (and never adds atoms). Extra atoms only ever
appear appended at the end, so truncating to NATOM keeps 1..N contiguous.

Bond orders: obabel's where it saw the bond. A bond added from the PSF gets
order 1, except a terminal oxygen on a carbon within CARBONYL_MAX A, which is
written as a double bond (carbonyl obabel missed because the patch geometry
stretched it). No molecule-specific indices are hard-coded.
"""
import sys, math

CARBONYL_MAX = 1.35   # A; C=O is ~1.23, C-O single ~1.43

def psf_meta(path):
    L = open(path).read().splitlines()
    ai = next(i for i, l in enumerate(L) if "!NATOM" in l)
    natom = int(L[ai].split()[0])
    names = [L[ai + 1 + k].split()[4] for k in range(natom)]
    bi = next(i for i, l in enumerate(L) if "!NBOND" in l)
    nb = int(L[bi].split()[0]); nums = []
    for l in L[bi + 1:]:
        if "!" in l:
            break
        nums += [int(x) for x in l.split()]
    nums = nums[:2 * nb]
    bonds = {frozenset((nums[2*k], nums[2*k+1])) for k in range(nb)}
    return natom, names, bonds

def section(L, tag):
    s = next(i for i, l in enumerate(L) if l.strip() == tag)
    e = s + 1
    while e < len(L) and not L[e].lstrip().startswith("@<TRIPOS>"):
        e += 1
    return s, e

def main(mol2_in, psf, mol2_out):
    natom, names, truth = psf_meta(psf)
    L = open(mol2_in).read().splitlines()

    # --- atoms: keep only the first NATOM, sanity-check name alignment ---
    a0, a1 = section(L, "@<TRIPOS>ATOM")
    atom_lines = [l for l in L[a0+1:a1] if l.strip()]
    mismatches = [(i+1, atom_lines[i].split()[1], names[i])
                  for i in range(natom) if atom_lines[i].split()[1] != names[i]]
    if mismatches:
        print("WARNING: atom-name mismatch vs PSF (index alignment broken!):", mismatches[:5])
    dropped_atoms = len(atom_lines) - natom
    atom_block = ["@<TRIPOS>ATOM"] + atom_lines[:natom]

    # element + coordinates per atom index (1-based) for the carbonyl rule
    elem, xyz = {}, {}
    for i, l in enumerate(atom_lines[:natom], 1):
        f = l.split()
        elem[i] = f[5].split(".")[0].upper()
        xyz[i] = tuple(float(v) for v in f[2:5])
    degree = {}
    for k in truth:
        for i in k:
            degree[i] = degree.get(i, 0) + 1

    def order_for_missing(k):
        i, j = sorted(k)
        pair = {elem.get(i), elem.get(j)}
        if pair == {"C", "O"}:
            o = i if elem[i] == "O" else j
            d = math.dist(xyz[i], xyz[j])
            if degree.get(o, 0) == 1 and d < CARBONYL_MAX:
                return "2"
        return "1"

    # --- bonds: PSF set is the truth, obabel supplies orders where it agreed ---
    b0, b1 = section(L, "@<TRIPOS>BOND")
    kept, present, raw = [], set(), 0
    for l in L[b0+1:b1]:
        if not l.strip():
            continue
        raw += 1
        f = l.split(); k = frozenset((int(f[1]), int(f[2])))
        if k not in truth:
            continue
        kept.append((int(f[1]), int(f[2]), f[3])); present.add(k)
    for k in truth - present:
        i, j = sorted(k); kept.append((i, j, order_for_missing(k)))

    # --- carbonyl rule, applied to EVERY C-O bond, not just the missing ones ---
    # A terminal O (one PSF bond) on a carbon is a carbonyl: force order 2 and
    # sp2 SYBYL types. obabel writes these as single/sp3 when a neighbouring
    # bond (e.g. the stretched junction N-C) confuses its hybridisation guess.
    retype = {}
    fixed = []
    for n, (i, j, o) in enumerate(kept):
        k = frozenset((i, j))
        if order_for_missing(k) == "2" and o != "2":
            kept[n] = (i, j, "2"); fixed.append(k)
        if order_for_missing(k) == "2":
            c, ox = (i, j) if elem[i] == "C" else (j, i)
            retype[c] = "C.2"; retype[ox] = "O.2"
    if fixed:
        print("carbonyl bonds forced to order 2:", [tuple(sorted(k)) for k in fixed])
    if retype:
        for idx, t in retype.items():
            f = atom_lines[idx-1].split()
            if f[5] != t:
                atom_lines[idx-1] = atom_lines[idx-1].replace(f" {f[5]} ", f" {t:<5} ", 1) if f" {f[5]} " in atom_lines[idx-1] else atom_lines[idx-1]
        atom_block = ["@<TRIPOS>ATOM"] + atom_lines[:natom]
    bond_block = ["@<TRIPOS>BOND"] + [f"{n:>6} {i:>6} {j:>6} {o:>2}"
                                      for n, (i, j, o) in enumerate(kept, 1)]

    # --- counts + reassemble ---
    out = L[:a0] + atom_block + bond_block + L[b1:]
    m0, _ = section(out, "@<TRIPOS>MOLECULE")
    c = out[m0+2].split(); c[0], c[1] = str(natom), str(len(kept))
    out[m0+2] = " " + " ".join(c)
    # Molecule name becomes CGenFF's RESI name: fix it here instead of downstream
    stem = mol2_out.rsplit("/", 1)[-1].split(".")[0].lower()
    out[m0+1] = {"naa": "NAA", "ntrm": "NPC"}.get(stem, stem.upper()[:4])
    open(mol2_out, "w").write("\n".join(out) + "\n")
    print(f"atoms {len(atom_lines)} -> {natom} (dropped {dropped_atoms}) | "
          f"bonds {raw} -> {len(kept)} (spurious dropped {raw-len(present)}, "
          f"missing added {len(truth-present)})")

if __name__ == "__main__":
    main(*sys.argv[1:4])