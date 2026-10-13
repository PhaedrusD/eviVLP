#!/usr/bin/env python3
"""write_fixed_atoms.py <pdb_in> <pdb_out> [--only-residues NAA NPC]

Write a NAMD fixedAtoms flag file for a modified monomer: B = 0.00 for atoms
that may move, B = 1.00 for atoms that are held.

The criterion is psfgen's own: `writepdb` sets occupancy 0.00 for coordinates
it GUESSED from the internal-coordinate table and 1.00 for coordinates that
came from the input structure. So occupancy 0 atoms are exactly what this
pipeline invented - the drug adduct, any side-chain atoms the mutation added,
and the hydrogens psfgen placed - while occupancy 1 atoms are crystallographic.

Point NAMD at it with
    fixedAtoms     on
    fixedAtomsFile <pdb_out>
    fixedAtomsCol  B

Why fix the rest: the capsid is assembled by applying crystallographic symmetry
matrices to this monomer, and those matrices are valid only for the conformation
they were derived from. If crystallographic heavy atoms relax here, every
subunit interface in the assembled capsid is placed wrong - and the monomer
conformation is itself held by inter-subunit contacts that do not exist for an
isolated subunit, so it would drift away from the right answer, not toward it.

--only-residues additionally restricts movement to the named residues, so just
the adduct and built side chain relax and even added hydrogens elsewhere stay
put. The looser default also tidies protein hydrogens, which is harmless:
hydrogens play no part in the symmetry matrices.
"""
import sys

FREE, FIXED = 0.00, 1.00


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__.strip().splitlines()[0])
        sys.exit(1)
    src, dst = args[0], args[1]
    only = None
    if '--only-residues' in args:
        i = args.index('--only-residues')
        only = {r.upper() for r in args[i + 1:]} or {'NAA', 'NPC'}

    out, free_by_res, n_fixed, no_occ = [], {}, 0, True
    for line in open(src):
        if line.startswith(('ATOM', 'HETATM')):
            resname = line[17:20].strip().upper()
            try:
                occ = float(line[54:60])
            except ValueError:
                occ = 1.0
            if occ != 1.0:
                no_occ = False
            guessed = (occ == 0.0)
            movable = guessed and (only is None or resname in only)
            if movable:
                free_by_res[resname] = free_by_res.get(resname, 0) + 1
            else:
                n_fixed += 1
            line = line[:60] + f"{FREE if movable else FIXED:6.2f}" + line[66:]
        out.append(line)

    n_free = sum(free_by_res.values())
    if no_occ:
        print(f"ERROR: every atom in {src} has occupancy 1.00, so nothing is marked as\n"
              f"       guessed. This is not a psfgen writepdb output, or the occupancies\n"
              f"       were overwritten. Refusing to guess which atoms may move.")
        sys.exit(1)
    if not n_free:
        print(f"ERROR: no movable atoms selected in {src}"
              + (f" for residues {sorted(only)}" if only else "") + ".")
        sys.exit(1)

    open(dst, 'w').write(''.join(out))
    print(f"{dst}: {n_free} free (guessed) atoms, {n_fixed} fixed")
    for r, n in sorted(free_by_res.items(), key=lambda kv: -kv[1]):
        print(f"    {r}: {n}")


if __name__ == '__main__':
    main()
