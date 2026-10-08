#!/usr/bin/env python3
"""validate_rtf_ic.py RESNAME [RESNAME ...]

Sanity-check the RESI blocks a previous IC_* run appended to
top_all36_cgenff_CBD.rtf, before psfgen ever sees them.

Checks per residue:
  1. no IC entry with a repeated atom (a degenerate frame; psfgen's guesscoord
     takes the first satisfiable IC, so such an entry can beat the correct one
     and place the atom from a 0-degree frame -> atom lands anywhere)
  2. every atom name used in BOND/IC lines exists in the residue (catches
     CHARMM22-era names such as HN for the amide hydrogen, which is H in
     CHARMM36/CGenFF)
  3. every atom has at least one non-degenerate IC that could place it

Exit 1 if any check fails.
"""
import re, sys

RTF = 'top_all36_cgenff_CBD.rtf'
EXTERNAL = {'-C', '+N', '+CA', '-CA', '-O', '+C'}


def blocks(path):
    cur, out = None, {}
    for line in open(path):
        m = re.match(r'^\s*(RESI|PRES)\s+(\S+)', line, re.I)
        if m:
            cur = m.group(2).upper() if m.group(1).upper() == 'RESI' else None
            if cur:
                out[cur] = []
            continue
        if cur is not None:
            if re.match(r'^\s*END\b', line, re.I):
                cur = None
            else:
                out[cur].append(line.rstrip('\n'))
    return out


def check(name, lines):
    atoms, ics, bonds, bad = [], [], [], []
    for l in lines:
        p = l.split('!')[0].split()   # drop trailing comments
        if not p:
            continue
        if p[0] == 'ATOM' and len(p) >= 2:
            atoms.append(p[1])
        elif p[0] == 'IC' and len(p) >= 5:
            ics.append([a.lstrip('*') for a in p[1:5]])
        elif p[0] == 'BOND':
            bonds += p[1:]
    known = set(atoms) | EXTERNAL
    ok = True

    degen = [q for q in ics if len(set(q)) != 4]
    if degen:
        ok = False
        bad.append(f"  {len(degen)} IC entries with a repeated atom, e.g. {' '.join(degen[0])}")

    unknown = sorted({a for q in ics for a in q if a not in known} |
                     {a for a in bonds if a not in known})
    if unknown:
        ok = False
        bad.append(f"  atom names used but not defined in the residue: {unknown}")

    placeable = {q[3] for q in ics if len(set(q)) == 4} | {q[0] for q in ics if len(set(q)) == 4}
    orphans = [a for a in atoms if a not in placeable]
    if orphans:
        ok = False
        bad.append(f"  atoms with no usable IC entry (guesscoord will improvise): {orphans}")

    print(f"RESI {name}: {len(atoms)} atoms, {len(ics)} IC entries -> {'OK' if ok else 'PROBLEMS'}")
    for b in bad:
        print(b)
    return ok


def main():
    found = blocks(RTF)
    allok = True
    for want in [n.upper() for n in sys.argv[1:]]:
        if want not in found:
            print(f"RESI {want}: not present in {RTF}")
            continue
        allok = check(want, found[want]) and allok
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
