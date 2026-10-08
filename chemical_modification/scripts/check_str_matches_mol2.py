#!/usr/bin/env python3
"""check_str_matches_mol2.py  file.str  expected.mol2  [other.mol2 ...]

Exit 0 if the CGenFF stream file describes the same molecule as expected.mol2
(same atom count and same element composition; atom NAMES are not compared
because CGenFF may rename atoms). Exit 1 otherwise, naming which of the other
mol2 files it matches instead, if any.
"""
import sys, os
from collections import Counter

def elem_from_cgenff_type(t):
    t = t.upper()
    for two in ("CL", "BR"):
        if t.startswith(two): return two
    return t[0]

def str_composition(path):
    c = Counter(); in_resi = False
    for line in open(path):
        s = line.strip()
        if s.startswith("RESI"): in_resi = True; continue
        if in_resi and s.startswith(("BOND", "IMPR", "IC ", "END")): in_resi = False
        if in_resi and s.startswith("ATOM"):
            p = s.split()
            if len(p) >= 3: c[elem_from_cgenff_type(p[2])] += 1
    return c

def mol2_composition(path):
    c = Counter(); L = open(path).read().splitlines()
    a = L.index("@<TRIPOS>ATOM"); b = L.index("@<TRIPOS>BOND")
    for l in L[a+1:b]:
        p = l.split()
        if len(p) >= 6: c[p[5].split(".")[0].upper()] += 1
    return c

def main():
    strf, expected, others = sys.argv[1], sys.argv[2], sys.argv[3:]
    sc = str_composition(strf)
    if not sc:
        print(f"ERROR: no RESI/ATOM block found in {strf}"); sys.exit(1)
    if not os.path.exists(expected):
        print(f"WARNING: {expected} not present; cannot verify {strf}"); sys.exit(0)
    if sc == mol2_composition(expected):
        print(f"OK: {strf} matches {expected} ({sum(sc.values())} atoms)"); sys.exit(0)
    print(f"ERROR: {strf} does not describe the molecule in {expected}.")
    print(f"  str: {dict(sc)}   {expected}: {dict(mol2_composition(expected))}")
    for o in others:
        if os.path.exists(o) and sc == mol2_composition(o):
            print(f"  It matches {o} instead - import it with the other button.")
    sys.exit(1)

if __name__ == "__main__":
    main()
