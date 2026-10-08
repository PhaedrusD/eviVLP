#!/usr/bin/env python3
"""install_nterm_reference.py <reference.str> [RESNAME]

One-time setup: install an N-terminus model residue (e.g. RESI BBT, the TMV
serine reference) into top_all36_cgenff_CBD_template.rtf, so psfgen can build
segment A of patched_ntrm.pgn from ntrm_ref.pdb.

Only the RESI topology block is needed - psfgen reads topology, not parameters.
The conjugate's own parameters still come from the NPC stream file later in the
pipeline. An existing block of the same name is replaced.
"""
import re, sys, shutil

RTF = 'top_all36_cgenff_CBD_template.rtf'


def resi_block(path, name=None):
    lines = open(path).read().split('\n')
    start = None
    for i, l in enumerate(lines):
        m = re.match(r'^\s*RESI\s+(\S+)', l, re.I)
        if m and (name is None or m.group(1).upper() == name.upper()):
            start, found = i, m.group(1).upper()
            break
    if start is None:
        sys.exit(f"ERROR: no RESI block{' named '+name if name else ''} in {path}")
    out = []
    for l in lines[start:]:
        if re.match(r'^\s*(END|RETURN)\b', l, re.I):
            break
        if l.strip().lower().startswith('read param'):
            break
        if len(out) and re.match(r'^\s*(RESI|PRES)\s+\S+', l, re.I):
            break
        out.append(l)
    return found, out


def strip_existing(lines, name):
    out, skipping, n = [], False, 0
    for l in lines:
        m = re.match(r'^\s*(RESI|PRES)\s+(\S+)', l, re.I)
        if m or re.match(r'^\s*END\b', l, re.I):
            skipping = bool(m) and m.group(2).upper() == name.upper()
            if skipping:
                n += 1
        if not skipping:
            out.append(l)
    return out, n


def main():
    src = sys.argv[1]
    want = sys.argv[2] if len(sys.argv) > 2 else None
    name, block = resi_block(src, want)
    atoms = [l.split()[1] for l in block if l.strip().startswith('ATOM')]
    print(f"RESI {name}: {len(atoms)} atoms -> {atoms}")

    shutil.copyfile(RTF, RTF + '.bak')
    lines = open(RTF).read().split('\n')
    lines, removed = strip_existing(lines, name)
    if removed:
        print(f"  replaced {removed} existing RESI {name} block(s)")

    # insert before the final END
    ends = [i for i, l in enumerate(lines) if re.match(r'^\s*END\b', l, re.I)]
    at = ends[-1] if ends else len(lines)
    lines[at:at] = block + ['']
    open(RTF, 'w').write('\n'.join(lines))
    print(f"  inserted into {RTF} (backup: {RTF}.bak)")
    print(f"  now run 'Run Renaming' in the GUI to refresh the working rtf from the template")


if __name__ == '__main__':
    main()
