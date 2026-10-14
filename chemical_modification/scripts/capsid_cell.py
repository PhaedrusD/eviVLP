"""Write the NAMD periodic cell for a solvated capsid from its own PDB.

usage: python3 capsid_cell.py <solvated_capsid.pdb>

VMD's solvate/autoionize record the water-box size in the CRYST1 line but
NAMD does not read CRYST1, and the box centre is wherever the capsid happened
to sit (Qbeta is centred at x=74, the TMV rod at z=36). Hand-copied numbers
drift - the Qbeta conf carried 320 A for a 318 A box - so the minimize confs
`source` the file this writes instead:

    <stem>_cell.namd      cellBasisVector1/2/3 from CRYST1, cellOrigin from
                          the geometric centre of all atoms.
"""
import os, sys


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    pdb = sys.argv[1]
    a = b = c = None
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    with open(pdb) as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                a, b, c = (float(line[6:15]), float(line[15:24]), float(line[24:33]))
            elif line.startswith(("ATOM", "HETATM")):
                for i, (s, e) in enumerate(((30, 38), (38, 46), (46, 54))):
                    v = float(line[s:e])
                    if v < lo[i]: lo[i] = v
                    if v > hi[i]: hi[i] = v
    if a is None:
        sys.exit(f"{pdb}: no CRYST1 record - was this written by solvate/autoionize?")
    if lo[0] == float("inf"):
        sys.exit(f"{pdb}: no ATOM records")
    centre = [(l + h) / 2 for l, h in zip(lo, hi)]
    extent = [h - l for l, h in zip(lo, hi)]
    for dim, ext, name in zip((a, b, c), extent, "xyz"):
        if abs(dim - ext) > 5.0:
            print(f"warning: CRYST1 {name}={dim:.1f} but atoms span {ext:.1f} A")
    out = os.path.splitext(os.path.basename(pdb))[0] + "_cell.namd"
    with open(out, "w") as fh:
        fh.write(f"# periodic cell written by capsid_cell.py from {os.path.basename(pdb)}\n")
        fh.write(f"cellBasisVector1  {a:9.3f}  0.000  0.000\n")
        fh.write(f"cellBasisVector2    0.000 {b:9.3f}  0.000\n")
        fh.write(f"cellBasisVector3    0.000  0.000 {c:9.3f}\n")
        fh.write(f"cellOrigin       {centre[0]:9.3f} {centre[1]:9.3f} {centre[2]:9.3f}\n")
    print(f"wrote {out}: box {a:.2f} x {b:.2f} x {c:.2f} A, origin "
          f"({centre[0]:.2f}, {centre[1]:.2f}, {centre[2]:.2f})")


if __name__ == "__main__":
    main()
