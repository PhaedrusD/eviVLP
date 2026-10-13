#!/usr/bin/env python3
"""filter_cgenff_nbfix.py <release.prm> <output.prm> [other.prm ...]

Write a copy of a released CGenFF parameter file with the NBFIX lines removed
whose atom types have no vdW (NONBONDED) entry in the files actually loaded.

Why this is needed: CGenFF's NBFIX section carries pair corrections against the
nucleic-acid, carbohydrate and lipid force fields (types NN1..NN4, ON1, ON1C,
OC2D1..OC2D4, OBL, OCL, NH3L). Its own header says so:

    If protein, nucleic acid and carbohydrate topology and parameter files are
    NOT read prior to the CGenFF files "WARNING: ATOMS IN NBFIX" errors will be
    encountered. These can be passed by setting BOMLEV ( -2)...

CHARMM can be told to continue; NAMD cannot, and dies with
"FATAL ERROR: DIDN'T FIND vdW PARAMETER FOR ATOM TYPE OC2D1". A VLP capsid in
water with ions contains no nucleic acid, carbohydrate or lipid, so those pair
corrections can never apply and dropping them changes no interaction present in
the system.

The released file is never edited: keep it as the provenance copy and regenerate
the filtered one from it, so a run stays reproducible.
"""
import sys

SECTIONS = ('ATOMS', 'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'IMPROPERS',
            'CMAP', 'NONBONDED', 'NBFIX', 'HBOND', 'END', 'THOLE')


def vdw_types(path):
    out, sec = set(), None
    for l in open(path, errors='replace'):
        t = l.split('!')[0].split()
        if not t:
            continue
        if t[0].upper() in SECTIONS:
            sec = t[0].upper()
            continue
        if sec == 'NONBONDED':
            out.add(t[0])
    return out


def main():
    if len(sys.argv) < 3:
        print(__doc__.strip().splitlines()[0])
        sys.exit(1)
    release, output, others = sys.argv[1], sys.argv[2], sys.argv[3:]

    available = vdw_types(release)
    for o in others:
        available |= vdw_types(o)

    lines = open(release, errors='replace').read().split('\n')
    out, sec, dropped = [], None, []
    for l in lines:
        t = l.split('!')[0].split()
        head = t[0].upper() if t else ''
        if head in SECTIONS:
            sec = head
            out.append(l)
            continue
        if sec == 'NBFIX' and len(t) >= 2:
            unknown = [x for x in t[:2] if x not in available]
            if unknown:
                dropped.append((l.strip(), unknown))
                continue
        out.append(l)

    banner = [
        '* FILTERED COPY - do not edit by hand. Regenerate with:',
        f'*   python3 filter_cgenff_nbfix.py {release} {output} ' + ' '.join(others),
        f'* {len(dropped)} NBFIX line(s) removed: their atom types have no vdW entry in the',
        '* loaded parameter set (CGenFF ships pair corrections against the nucleic-acid,',
        '* carbohydrate and lipid force fields). NAMD treats such a reference as fatal.',
        '*',
    ]
    open(output, 'w').write('\n'.join(banner + out))

    types = sorted({u for _, us in dropped for u in us})
    print(f"{output}: removed {len(dropped)} NBFIX line(s)")
    print(f"   unsatisfied types: {types}")
    print(f"   vdW types available: {len(available)}")


if __name__ == '__main__':
    main()
