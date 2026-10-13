#!/usr/bin/env python3
"""complete_cross_ff.py <raw.str> <renamed.str> <label> [prm ...]

Fill the gap that renaming protein-side atom types opens up.

The novel amino acid is parameterised whole by CGenFF, then the protein-side
atoms are renamed back to CHARMM36 types. A term CGenFF satisfied from its own
RELEASE library never appears in the stream file's parameter section, so nothing
re-emits it under the new type - and after renaming it exists nowhere:

    CGenFF release has   CG321 CG2O1 NG2S1   80.00  116.50
    N is renamed NG2S1 -> NH1, so the structure now needs
                         CG321 CG2O1 NH1    <- undefined
    => FATAL ERROR: UNABLE TO FIND ANGLE PARAMETERS FOR CG321 CG2O1 NH1

This script enumerates every bonded term the residue actually contains (bonds,
angles, dihedrals, impropers from the BOND/IMPR lists), checks whether the
RENAMED type combination is available in the loaded parameter files, and where
it is not, copies the numbers from the equivalent ORIGINAL-type term, writing
them into par_GUI.prm under a labelled block. It only ever copies values that a
force field already defines - it never invents a parameter.
"""
import sys

SECTIONS = ('ATOMS', 'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'IMPROPERS',
            'CMAP', 'NONBONDED', 'NBFIX', 'HBOND', 'END', 'THOLE')
NFIELDS = {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}


def resi_block(path):
    """atom name -> type, plus BOND and IMPR tuples, from the first RESI."""
    types, bonds, imprs = {}, [], []
    in_resi, in_params = False, False
    for line in open(path, errors='replace'):
        s = line.split('!')[0].strip()
        if not s:
            continue
        if 'read param' in line.lower():
            in_params = True
        if s.upper().startswith('RESI'):
            in_resi, in_params = True, False
            continue
        if in_resi and not in_params:
            t = s.split()
            k = t[0].upper()
            if k == 'ATOM' and len(t) >= 3:
                types[t[1]] = t[2]
            elif k in ('BOND', 'DOUBLE', 'TRIPLE'):
                n = t[1:]
                for i in range(0, len(n) - 1, 2):
                    bonds.append((n[i], n[i + 1]))
            elif k == 'IMPR' and len(t) >= 5:
                for i in range(1, len(t) - 3, 4):
                    imprs.append(tuple(t[i:i + 4]))
            elif k in ('END', 'RETURN'):
                in_resi = False
    return types, bonds, imprs


def load_params(paths, exclude_label=None):
    """(section, type-tuple) -> list of full data lines, both directions.

    exclude_label skips lines inside '! === <label> ...' blocks, so a re-run does
    not count its OWN previous output as already-available and then replace it -
    which would silently drop every term it had filled the first time.
    """
    table = {}
    for p in paths:
        sec, skipping = None, False
        for line in open(p, errors='replace'):
            if line.startswith('! === '):
                skipping = bool(exclude_label) and exclude_label in line
                continue
            if skipping:
                if not line.strip():
                    skipping = False
                continue
            body = line.split('!')[0].rstrip()
            if not body.strip():
                continue
            head = body.split()[0].upper()
            if head in SECTIONS:
                sec = 'IMPROPERS' if head == 'IMPROPER' else head
                continue
            if sec not in NFIELDS:
                continue
            t = body.split()
            n = NFIELDS[sec]
            if len(t) < n + 1:
                continue
            key = tuple(t[:n])
            for k in (key, tuple(reversed(key))):
                table.setdefault((sec, k), []).append(body)
    return table


def lookup(table, sec, key):
    """Exact match, reversed, then dihedral/improper wildcards."""
    for k in (tuple(key), tuple(reversed(key))):
        if (sec, k) in table:
            return table[(sec, k)]
    if sec in ('DIHEDRALS', 'IMPROPERS') and len(key) == 4:
        a, b, c, d = key
        cands = [('X', b, c, 'X'), ('X', c, b, 'X')]
        if sec == 'IMPROPERS':
            cands += [(a, 'X', 'X', d), (d, 'X', 'X', a)]
        for k in cands:
            if (sec, k) in table:
                return table[(sec, k)]
    return None


def retype(line, n, newtypes):
    """Replace the leading n type fields, keeping the numeric tail unchanged."""
    tail = line.split()[n:]
    return ' '.join('{:<6}'.format(x) for x in newtypes) + '  ' + ' '.join(tail)


def main():
    raw, renamed, label = sys.argv[1], sys.argv[2], sys.argv[3]
    prms = sys.argv[4:] or ['par_all36m_prot.prm', 'par_all36_cgenff.prm', 'par_GUI.prm']

    old_types, _, _ = resi_block(raw)
    new_types, bonds, imprs = resi_block(renamed)
    if not new_types:
        sys.exit(f"ERROR: no RESI atoms found in {renamed}")

    adj = {}
    for a, b in bonds:
        if a in new_types and b in new_types:
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)

    # enumerate the terms the residue actually contains
    terms = {'BONDS': set(), 'ANGLES': set(), 'DIHEDRALS': set(), 'IMPROPERS': set()}
    for a in adj:
        for b in adj[a]:
            terms['BONDS'].add(tuple(sorted((a, b))))
            for c in adj[b]:
                if c == a:
                    continue
                terms['ANGLES'].add((a, b, c) if a < c else (c, b, a))
                for d in adj[c]:
                    if d in (a, b):
                        continue
                    terms['DIHEDRALS'].add((a, b, c, d) if a < d else (d, c, b, a))
    for q in imprs:
        if all(x in new_types for x in q):
            terms['IMPROPERS'].add(q)

    table = load_params(prms, exclude_label=label + ' junction')
    # Donor order matters: the stream file first, because it carries the values
    # the CGenFF SERVICE assigned, which may come from a newer library than the
    # locally installed release. The release is the fallback.
    str_donors = load_params([raw])
    release_donors = load_params([p for p in prms if 'cgenff' in p or 'prot' in p])

    def find_donor(sec, key):
        hit = lookup(str_donors, sec, key)
        if hit:
            return hit, 'stream file'
        hit = lookup(release_donors, sec, key)
        return (hit, 'CGenFF release') if hit else (None, None)

    filled, unresolved, stats, sources = {s: [] for s in NFIELDS}, [], {}, {}
    for sec in ('BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS'):
        n = NFIELDS[sec]
        for names in sorted(terms[sec]):
            new = tuple(new_types[x] for x in names)
            if lookup(table, sec, new):
                continue                                  # already available
            old = tuple(old_types.get(x, new_types[x]) for x in names)
            if old == new:
                unresolved.append((sec, names, new, 'no rename on this term'))
                continue
            src, origin = find_donor(sec, old)
            if not src:
                unresolved.append((sec, names, new, f'no donor for {" ".join(old)}'))
                continue
            sources[origin] = sources.get(origin, 0) + 1
            for line in src:
                filled[sec].append(retype(line, n, new) +
                                   f'  ! {label} cross-FF: from {" ".join(old)}\n')
            stats[sec] = stats.get(sec, 0) + 1

    # Was the stream file generated with CGenFF's "copy existing parameters" on?
    n_str_params = sum(len(v) for (sec, _), v in str_donors.items()) // 2
    if n_str_params < 20:
        print(f"   NOTE: {raw} carries only ~{n_str_params} parameter line(s), so it was")
        print(f"         generated with CGenFF's 'include parameters already in CGenFF'")
        print(f"         OPTION OFF. Every other term must then come from the locally")
        print(f"         installed par_all36_cgenff.prm, so a newer service library than")
        print(f"         your local release leaves terms with no donor. Turning that")
        print(f"         option ON makes the stream file self-contained.")

    total = sum(len(v) for v in filled.values())
    print(f"{label}: terms examined "
          f"{ {k: len(v) for k, v in terms.items()} }")
    print(f"   filled from original-type donors: { {k: v for k, v in stats.items()} } "
          f"-> {total} parameter line(s)")
    if sources:
        print(f"   donor provenance: {sources}")
    if unresolved:
        real = [u for u in unresolved if not u[3].startswith('no rename')]
        if real:
            print(f"   {len(real)} term(s) with no donor either:")
            for sec, names, new, why in real[:8]:
                print(f"     [{sec}] {' '.join(names)} = {' '.join(new)}  ({why})")

    import importlib.util
    spec = importlib.util.spec_from_file_location('adder', 'add_naa_params_to_prm.py')
    adder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adder)
    added = adder.insert('par_GUI.prm', filled, label + ' junction')
    print(f"   inserted {added} line(s) into par_GUI.prm")
    return 0 if total or not unresolved else 0


if __name__ == '__main__':
    sys.exit(main())
