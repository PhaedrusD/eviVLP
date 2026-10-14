"""Point this case study's charge overrides at whatever atom names CGenFF just assigned.

WHY THIS IS HERE AND NOT IN THE GUI
  The GUI's apply_charge_overrides.py addresses atoms by name and knows no
  chemistry - deliberately, because it has to serve any small molecule. But
  CGenFF names atoms from the order they appear in the upload, so every re-run
  renumbers them and a name-based override file goes stale. Recognising WHICH
  atom is the radical nitrogen requires knowing what a nitroxide looks like,
  and that knowledge is specific to this one case study. So it lives here.

  A different molecule with chemistry CGenFF cannot type needs its own small
  mapper like this one, or its names written out by hand after each run. The
  base GUI stays general; the chemistry stays with the case that needs it.

HOW IT WORKS
  Each line of the @CHARGES block in nitroxide_overrides.prm carries a
  `role=` tag. This script resolves each role from the residue's own BOND
  records and rewrites only the NAME column, leaving the types, charges and
  target untouched - those are the published values and are not ours to edit.
  Atoms themselves are never renamed: they keep the names CGenFF gave them,
  which the .pdb, the patch and IC_naa all depend on.

ROLES, for a 2,2,6,6-tetrasubstituted piperidin-1-oxyl
  radical_N    the N bonded to a terminal (degree-1) O        1
  radical_O    that O, which must carry no hydrogen           1
  alpha_C      the carbon neighbours of radical_N             2
  methyl_C     carbons on alpha_C carrying three H            4
  ring_CH2     the other carbons on alpha_C                   2
  ring_C_sub   the carbon shared by both ring_CH2             1
  linker_C1    first substituent carbon past ring_C_sub       1

  linker_C1 steps PAST a heteroatom on ring_C_sub. In this conjugate that atom
  is the protonated linker amine, whose charges carry the residue's +1, so the
  balancing remainder must not land on it.

Resolution is strict: a role that resolves to the wrong number of atoms is an
error naming what was found, and nothing is written.

usage: python3 map_override_atoms.py <str_file> <override_file> [--dry-run]
       e.g. python3 map_override_atoms.py nad.str nitroxide_overrides.prm
"""
import sys, re, os, shutil

EXPECT = {"radical_N": 1, "radical_O": 1, "alpha_C": 2, "methyl_C": 4,
          "ring_CH2": 2, "ring_C_sub": 1, "linker_C1": 1}


def read_residue(path, resi):
    """Atom types and adjacency for one RESI block."""
    types, adj, on = {}, {}, False
    for raw in open(path, errors="replace"):
        if raw.startswith("RESI"):
            f = raw.split()
            on = len(f) > 1 and f[1] == resi
            continue
        if not on:
            continue
        if raw.startswith(("RESI", "PRES")):
            break
        if raw.startswith("ATOM"):
            f = raw.split()
            if len(f) >= 3:
                types[f[1]] = f[2]
        elif raw.startswith(("BOND", "DOUBLE", "TRIPLE")):
            t = raw.split("!")[0].split()[1:]
            for a, b in zip(t[0::2], t[1::2]):
                adj.setdefault(a, set()).add(b)
                adj.setdefault(b, set()).add(a)
    return types, adj


def resolve(types, adj):
    el = lambda n: "H" if types.get(n, "").startswith("H") else types.get(n, "")[:1]
    hyd = lambda n: [x for x in adj.get(n, ()) if el(x) == "H"]
    carbons = lambda n: [x for x in adj.get(n, ()) if el(x) == "C"]

    cand = [n for n in types if el(n) == "N"
            and any(el(o) == "O" and len(adj.get(o, ())) == 1 for o in adj.get(n, ()))]
    if len(cand) != 1:
        raise ValueError(f"expected exactly one nitrogen bonded to a terminal oxygen, "
                         f"found {len(cand)}: {sorted(cand)}")
    rN = cand[0]
    rO = [o for o in adj[rN] if el(o) == "O" and len(adj.get(o, ())) == 1][0]
    if hyd(rO):
        raise ValueError(f"{rO} carries a hydrogen - that is a hydroxylamine, "
                         f"not an aminoxyl radical")

    alpha = sorted(carbons(rN))
    methyl, ch2 = [], []
    for a in alpha:
        for c in carbons(a):
            if c == rN:
                continue
            (methyl if len(hyd(c)) == 3 else ch2).append(c)
    methyl, ch2 = sorted(set(methyl)), sorted(set(ch2))

    sub = sorted(set.intersection(*[set(carbons(c)) for c in ch2]) - set(alpha)) if len(ch2) == 2 else []

    linker = []
    if len(sub) == 1:
        for nb in sorted(adj.get(sub[0], ())):
            if el(nb) in ("H", "") or nb in ch2:
                continue
            if el(nb) == "C":
                linker = [nb]
                break
            onward = [c for c in carbons(nb) if c != sub[0]]   # step past the heteroatom
            if len(onward) == 1:
                linker = onward
                break

    roles = {"radical_N": [rN], "radical_O": [rO], "alpha_C": alpha, "methyl_C": methyl,
             "ring_CH2": ch2, "ring_C_sub": sub, "linker_C1": linker}
    for role, want in EXPECT.items():
        if len(roles[role]) != want:
            raise ValueError(f"role {role}: expected {want} atom(s), resolved "
                             f"{len(roles[role])} {sorted(roles[role])}. This is not the "
                             f"2,2,6,6-tetrasubstituted nitroxide the case study describes.")
    return roles


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    if len(args) != 2:
        sys.exit("usage: python3 map_override_atoms.py <str_file> <override_file> [--dry-run]")
    str_file, ovr_file = args

    lines = open(ovr_file, errors="replace").read().splitlines(keepends=True)
    resi = None
    for L in lines:
        m = re.search(r"@CHARGES.*?residue\s*=\s*(\S+)", L)
        if m:
            resi = m.group(1); break
    if not resi:
        sys.exit(f"{ovr_file}: no '@CHARGES residue=<RESI>' line found")

    types, adj = read_residue(str_file, resi)
    if not types:
        sys.exit(f"{str_file}: RESI {resi} not found")
    if not adj:
        sys.exit(f"{str_file}: RESI {resi} has no BOND records, so roles cannot be resolved")
    try:
        roles = resolve(types, adj)
    except ValueError as e:
        sys.exit(f"{str_file}: could not resolve the nitroxide - {e}")

    print(f"{str_file}: RESI {resi}, {len(types)} atoms\n")
    for role in ("radical_N", "radical_O", "alpha_C", "methyl_C",
                 "ring_CH2", "ring_C_sub", "linker_C1"):
        print(f"  {role:<12s} -> {', '.join(roles[role])}")

    # rewrite only the name column of each role-tagged line, in role order
    queue = {r: list(v) for r, v in roles.items()}
    changed, out = 0, []
    for L in lines:
        m = re.match(r"^(!\s*)(\S+)(\s+\S+\s+[+-]?[\d.]+\s*!.*?\brole\s*=\s*)(\w+)(.*)$", L)
        if m and m.group(4) in queue:
            role = m.group(4)
            if not queue[role]:
                sys.exit(f"{ovr_file}: more lines tagged role={role} than the "
                         f"{EXPECT[role]} atom(s) it resolves to")
            new = queue[role].pop(0)
            if new != m.group(2):
                changed += 1
            out.append(f"{m.group(1)}{new:<5s}{m.group(3)}{role}{m.group(5)}\n")
            continue
        m = re.match(r"^(!\s*@BALANCE\s+)(\S+)(.*?\brole\s*=\s*)(\w+)(.*)$", L)
        if m and m.group(4) in queue:
            role = m.group(4)
            new = queue[role].pop(0) if queue[role] else m.group(2)
            if new != m.group(2):
                changed += 1
            out.append(f"{m.group(1)}{new}{m.group(3)}{role}{m.group(5)}\n")
            continue
        out.append(L)

    leftover = {r: v for r, v in queue.items() if v}
    if leftover:
        sys.exit(f"{ovr_file}: fewer role-tagged lines than atoms resolved - "
                 f"unused: { {r: sorted(v) for r, v in leftover.items()} }")

    print(f"\n{changed} name(s) rewritten" + ("" if changed else " (already current)"))
    if dry:
        print("--dry-run: no file written"); return
    if changed:
        backup = ovr_file + ".premap"
        if not os.path.exists(backup):
            shutil.copyfile(ovr_file, backup)
        open(ovr_file, "w").write("".join(out))
        print(f"wrote {ovr_file}")
        print("Now run 'Create Novel Amino Acid' in the GUI, or apply_charge_overrides.py directly.")


main()
