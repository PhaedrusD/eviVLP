#!/usr/bin/env python3
"""validate_prm.py <file.prm> [...]

Check a CHARMM parameter file the way NAMD's reader needs it: sections present,
each data line carrying the right number of atom types and numbers, and a
terminating END. CMAP blocks are recognised (8 atom types then a grid of bare
numbers) so their grid lines are not mistaken for malformed IMPROPERS - an
orphaned CMAP block with no `CMAP` header is exactly what makes NAMD die with
"FATAL ERROR: ABNORMAL EOF FOUND - buffer=*END*".
"""
import sys

SPEC = {'BONDS': (2, [2]), 'ANGLES': (3, [2, 4]), 'DIHEDRALS': (4, [3]),
        'IMPROPERS': (4, [3]), 'IMPROPER': (4, [3]),
        'NONBONDED': (1, [3, 6]), 'NBFIX': (2, [2, 4])}
ENDERS = {'END', 'RETURN', 'HBOND', 'ATOMS', 'THOLE'}


def isnum(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def check(path):
    sec, bad, counts, has_end = None, [], {}, False
    for n, raw in enumerate(open(path, errors='replace'), 1):
        line = raw.split('!')[0].rstrip()
        if not line.strip():
            continue
        head = line.split()[0].upper()
        if head == 'END':
            has_end = True
            sec = None
            continue
        if head == 'CMAP':
            sec = 'CMAP'
            continue
        if head in SPEC:
            sec = head
            continue
        if head in ENDERS:
            sec = None
            continue
        if sec is None or sec == 'CMAP':
            continue                      # CMAP grid data is free-form
        if line.rstrip().endswith('-'):
            continue                      # continuation
        toks = line.split()
        nt, oks = SPEC[sec]
        counts[sec] = counts.get(sec, 0) + 1
        nums = [t for t in toks[nt:] if isnum(t)]
        if (len(toks) < nt + min(oks) or any(isnum(t) for t in toks[:nt])
                or len(nums) not in oks):
            bad.append((n, sec, line[:72]))

    ok = not bad and has_end
    print(f"{path}: {counts}")
    print(f"   terminating END: {'yes' if has_end else 'NO - NAMD will hit ABNORMAL EOF if this file is read last'}")
    if bad:
        print(f"   {len(bad)} malformed line(s):")
        for n, s, l in bad[:8]:
            print(f"     {n} [{s}] {l}")
    if ok:
        print("   OK")
    return ok


if __name__ == '__main__':
    sys.exit(0 if all(check(p) for p in sys.argv[1:]) else 1)
