#!/usr/bin/env python3
"""
add_naa_params_to_prm.py  <str_file> <amino_acid> [label] [prm_file]

Collect every parameter the simulation needs that is NOT already in a released
force-field file, and insert it into par_GUI.prm.

A NAMD run loads three parameter files: par_all36m_prot.prm (CHARMM36 protein),
par_all36_cgenff.prm (released CGenFF), and par_GUI.prm. So par_GUI.prm must
carry exactly the remainder:

  * novel CGenFF terms for the drug - a CGenFF stream file's parameter section
    only ever contains terms absent from the released set, so every drug-only
    line qualifies;
  * cross-force-field terms at the junction - lines mixing CHARMM36 protein
    types with CGenFF drug types, which no single force field defines;
  * all-protein terms that CHARMM36 happens to lack - kept rather than assumed
    redundant, by checking par_all36m_prot.prm directly.

Protein-side types come from the rename map in rename_atom_types_naa.py, so the
two cannot drift apart. Blocks written by a previous run under the same label
are replaced, making the step re-runnable.
"""
import importlib.util
import os
import sys

PARAM_SECTIONS = {'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS'}
ORDER = ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']


def _rename_module():
    spec = importlib.util.spec_from_file_location('ren', 'rename_atom_types_naa.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def n_types(section):
    return {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}.get(section, 0)


def extract(str_file, protein_types, c36, has):
    """Group the lines par_GUI.prm must carry, by section."""
    out = {s: [] for s in PARAM_SECTIONS}
    kinds = {'cross': 0, 'drug': 0, 'gap': 0, 'redundant': 0}
    section, in_params = None, False

    for line in open(str_file):
        s = line.strip()
        if 'read param' in s.lower():
            in_params = True
            continue
        if not in_params:
            continue
        if s in PARAM_SECTIONS:
            section = s
            continue
        if s in ('END', 'RETURN'):
            section = None
            continue
        if not s or s.startswith('!') or section not in PARAM_SECTIONS:
            continue

        parts = s.split()
        n = n_types(section)
        if len(parts) < n + 2:
            continue
        types = parts[0:n]

        n_prot = sum(1 for t in types if t in protein_types)
        if n_prot == len(types):
            # every type is CHARMM36: redundant only if CHARMM36 really has it
            if has(c36, section, types):
                kinds['redundant'] += 1
                continue
            kinds['gap'] += 1
        elif n_prot == 0:
            kinds['drug'] += 1
        else:
            kinds['cross'] += 1

        out[section].append(line.rstrip('\n') + '\n')

    return out, kinds


def insert(prm_file, params, label):
    header = f"! === {label} cross-force-field "
    lines, skipping, removed = [], False, 0
    for l in open(prm_file):
        if l.startswith(header):
            skipping = True
            removed += 1
            continue
        if skipping:
            if l.strip() == '':
                skipping = False
            continue
        lines.append(l)
    if removed:
        print(f"  replaced {removed} previous {label} block(s)")

    # Anything from here on is not a bonded-parameter section, so a block must
    # never be inserted past it. par_GUI.prm has no CMAP, so keying only on CMAP
    # (as this did originally) would append an IMPROPERS block after NONBONDED.
    TAIL = ('CMAP', 'NONBONDED', 'NBFIX', 'HBOND', 'END')

    starts, tail = {}, None
    for i, l in enumerate(lines):
        s = l.strip()
        if s in PARAM_SECTIONS:
            starts[s] = i
        head = s.split()[0].upper() if s else ''
        if head in TAIL and tail is None:
            tail = i

    ends = {}
    for idx, sec in enumerate(ORDER):
        if sec not in starts:
            continue
        end = tail if tail is not None else len(lines)
        for nxt in ORDER[idx + 1:]:
            if nxt in starts:
                end = min(end, starts[nxt])
                break
        ends[sec] = end

    total = 0
    for sec in reversed(ORDER):
        block_params = params.get(sec, [])
        if not block_params or sec not in ends:
            continue
        block = [header + f"{sec} ===\n"] + block_params + ['\n']
        lines[ends[sec]:ends[sec]] = block
        total += len(block_params)

    open(prm_file, 'w').write(''.join(lines))
    return total


def main():
    if len(sys.argv) < 3:
        print(__doc__.strip().splitlines()[1])
        sys.exit(1)
    str_file = sys.argv[1]
    aa = sys.argv[2].lower()
    label = sys.argv[3] if len(sys.argv) > 3 else 'NAA'
    prm_file = sys.argv[4] if len(sys.argv) > 4 else 'par_GUI.prm'

    if not os.path.exists(str_file):
        print(f"ERROR: {str_file} not found")
        sys.exit(1)

    ren = _rename_module()
    if aa not in ren.AMINO_ACID_MAPS:
        print(f"ERROR: unknown amino acid '{aa}'; valid: {', '.join(ren.AMINO_ACID_MAPS)}")
        sys.exit(1)
    protein_types = set(ren.AMINO_ACID_MAPS[aa].values())
    c36 = ren.load_c36_params()

    print(f"Extracting parameters from {str_file} (map {aa}, protein types {sorted(protein_types)})")
    params, kinds = extract(str_file, protein_types, c36, ren.c36_has)
    print(f"  cross-force-field: {kinds['cross']}   drug-only (novel CGenFF): {kinds['drug']}"
          f"   CHARMM36 gaps kept: {kinds['gap']}   redundant dropped: {kinds['redundant']}")
    for sec in ORDER:
        if params[sec]:
            print(f"    {sec}: {len(params[sec])}")
    total = insert(prm_file, params, label)
    print(f"  inserted {total} parameters into {prm_file}")


if __name__ == '__main__':
    main()
