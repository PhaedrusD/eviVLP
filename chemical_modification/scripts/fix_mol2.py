#!/usr/bin/env python3
"""fix_mol2.py  in.mol2  patched.psf  out.mol2
Keep obabel's perceived bond ORDERS, use the PSF as truth for which bonds EXIST.
Bonds in the mol2 but not the PSF are spurious -> dropped. PSF bonds the mol2
missed -> added (order 1). Indices line up: psfgen writes PSF and PDB in the
same atom order, and obabel preserves PDB order in the mol2.
"""
import sys

FORCE_ORDER = {frozenset((65, 67)): "2"}   # optional safety pin for C13=O3; {} to disable

def psf_bonds(path):
    lines = open(path).read().splitlines()
    i = next(k for k, l in enumerate(lines) if "!NBOND" in l)
    nb = int(lines[i].split()[0]); nums = []
    for l in lines[i + 1:]:
        if "!" in l:
            break
        nums += [int(x) for x in l.split()]
    nums = nums[:2 * nb]
    return {frozenset((nums[2*k], nums[2*k+1])) for k in range(nb)}

def section(lines, tag):
    s = next(i for i, l in enumerate(lines) if l.strip() == tag)
    e = s + 1
    while e < len(lines) and not lines[e].lstrip().startswith("@<TRIPOS>"):
        e += 1
    return s, e

def main(mol2_in, psf, mol2_out):
    truth = psf_bonds(psf)
    lines = open(mol2_in).read().splitlines()
    b0, b1 = section(lines, "@<TRIPOS>BOND")
    kept, present = [], set()
    for l in lines[b0+1:b1]:
        if not l.strip():
            continue
        f = l.split(); k = frozenset((int(f[1]), int(f[2])))
        if k not in truth:                 # spurious -> drop
            continue
        kept.append((int(f[1]), int(f[2]), FORCE_ORDER.get(k, f[3]))); present.add(k)
    for k in truth - present:              # PSF bond obabel missed -> add
        i, j = sorted(k); kept.append((i, j, FORCE_ORDER.get(k, "1")))
    block = ["@<TRIPOS>BOND"] + [f"{n:>6} {i:>6} {j:>6} {o:>2}"
                                 for n, (i, j, o) in enumerate(kept, 1)]
    m0, _ = section(lines, "@<TRIPOS>MOLECULE")
    c = lines[m0+2].split(); c[1] = str(len(kept)); lines[m0+2] = " " + " ".join(c)
    open(mol2_out, "w").write("\n".join(lines[:b0] + block + lines[b1:]) + "\n")

if __name__ == "__main__":
    main(*sys.argv[1:4])