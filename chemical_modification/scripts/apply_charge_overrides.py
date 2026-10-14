"""Apply hand-curated atom types and partial charges to a CHARMM stream file.

CGenFF has no chemistry for some groups - a persistent radical is the usual
case - so it types them as the nearest closed-shell analogue and assigns
charges accordingly. Those types and charges have to be replaced by hand. This
script does that replacement from a data file, so the GUI carries the
mechanism and the case carries the numbers.

The override file is a CHARMM parameter file with a tagged comment block. NAMD
ignores comment lines, so one file stays loadable as parameters AND readable as
the charge record - a single source of truth.

    ! @CHARGES  residue=NAA  target=+1
    ! N1  NG301  +0.167    ! radical nitrogen
    ! O1  OG311  -0.379    ! radical oxygen, no hydrogen
    ! @BALANCE  C3 C8
    ! @END

  residue=  the RESI to edit.
  target=   the net charge the residue must end on. Omit to keep whatever net
            charge it already has.
  one line per atom: <name> <type> <charge>. Type may be '-' to keep the
            existing type and change only the charge.
  @BALANCE  atoms that absorb the residual, split evenly, so the residue
            reaches `target`. Required whenever the overrides change the net
            charge: nothing here can infer which atom should take the
            difference, and putting it on the wrong atom is worse than failing.

Re-running is safe: the script reads the current file, applies the same
absolute values, and recomputes the balance from scratch, so two runs give an
identical result.

usage: python3 apply_charge_overrides.py <str_file> <override_file> [--dry-run]
"""
import sys, re, os, shutil

# A residual bigger than this per atom is a protonation-state error, not rounding.
MAX_SHARE = 0.25

def parse_overrides(path):
    resi, target, atoms, balance = None, None, {}, []
    inblock = False
    for raw in open(path, errors="replace"):
        s = raw.strip()
        if not s.startswith("!"):
            continue
        s = s.lstrip("!").strip()
        if s.startswith("@CHARGES"):
            inblock = True
            m = re.search(r"residue\s*=\s*(\S+)", s)
            if m: resi = m.group(1)
            m = re.search(r"target\s*=\s*([+-]?\d+(?:\.\d+)?)", s)
            if m: target = float(m.group(1))
            continue
        if not inblock:
            continue
        if s.startswith("@END"):
            break
        if s.startswith("@BALANCE"):
            # strip any trailing comment first: the atom lines below already do
            # this, and a comment here would otherwise be read as atom names
            balance += s.split("!")[0].split()[1:]
            continue
        f = s.split("!")[0].split()
        if len(f) >= 3:
            try: q = float(f[2])
            except ValueError: continue
            atoms[f[0]] = (f[1], q)
    if resi is None:
        sys.exit(f"{path}: no '! @CHARGES residue=<RESI>' block found")
    if not atoms:
        sys.exit(f"{path}: @CHARGES block for {resi} lists no atoms")
    return resi, target, atoms, balance

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry-run" in sys.argv
    if len(args) != 2:
        sys.exit(__doc__.strip().splitlines()[-1])
    str_file, ovr_file = args
    resi, target, overrides, balance = parse_overrides(ovr_file)
    print(f"override file : {ovr_file}")
    print(f"residue       : {resi}" + (f"   target net charge {target:+.3f}" if target is not None else ""))

    lines = open(str_file, errors="replace").read().splitlines(keepends=True)

    # locate the RESI block and index its ATOM lines by atom name
    start = None
    for i, L in enumerate(lines):
        if L.startswith("RESI") and len(L.split()) > 1 and L.split()[1] == resi:
            start = i; break
    if start is None:
        sys.exit(f"{str_file}: RESI {resi} not found")
    idx, end = {}, len(lines)
    for i in range(start + 1, len(lines)):
        L = lines[i]
        if L.startswith(("RESI", "PRES")): end = i; break
        if L.startswith("ATOM"):
            f = L.split()
            if len(f) >= 4: idx[f[1]] = i

    missing = [a for a in list(overrides) + balance if a not in idx]
    if missing:
        sys.exit(f"{str_file}: RESI {resi} has no atom(s) {', '.join(missing)}\n"
                 f"  the override file does not match this molecule - check the atom "
                 f"names CGenFF assigned")

    def charge_of(i): return float(lines[i].split()[3])
    def retype(i, newtype, q):
        f = lines[i].split()
        name, typ = f[1], (f[2] if newtype == "-" else newtype)
        tail = lines[i].split("!", 1)
        comment = ("!" + tail[1].rstrip("\n")) if len(tail) > 1 else ""
        lines[i] = f"ATOM {name:<6s} {typ:<6s} {q:7.3f} {comment}".rstrip() + "\n"
        return typ

    before = sum(charge_of(i) for i in idx.values())

    print("\noverrides applied:")
    for name, (typ, q) in overrides.items():
        i = idx[name]
        f = lines[i].split()
        old_t, old_q = f[2], float(f[3])
        new_t = retype(i, typ, q)
        flag = "" if (old_t == new_t and abs(old_q - q) < 5e-4) else "  <-- changed"
        print(f"  {name:<6s} {old_t:<7s} {old_q:+.3f}  ->  {new_t:<7s} {q:+.3f}{flag}")

    mid = sum(charge_of(i) for i in idx.values())
    print(f"\nnet charge    : {before:+.4f} (as imported)  ->  {mid:+.4f} (after overrides)")

    if target is not None:
        residual = target - mid
        if abs(residual) < 5e-4:
            print("              : already on target, nothing to balance")
        elif not balance:
            sys.exit(f"\nresidual {residual:+.4f} must be absorbed but the override file "
                     f"gives no @BALANCE atoms.\n  Add '! @BALANCE <atom> [atom ...]' "
                     f"naming the atom(s) that should take it.")
        else:
            share = residual / len(balance)
            if abs(share) > MAX_SHARE:
                sys.exit(
                    f"\nREFUSING to balance: residual {residual:+.4f} over "
                    f"{len(balance)} atom(s) is {share:+.4f} each, beyond the "
                    f"{MAX_SHARE:.2f} e limit.\n"
                    f"  A residual this large is not a rounding difference. It has two\n"
                    f"  likely causes, and neither is fixed by redistributing charge:\n"
                    f"\n"
                    f"  1. The override UNDER-SPECIFIES the charge group. If CGenFF gave the\n"
                    f"     group the wrong FORMAL charge - typing a neutral radical's oxygen\n"
                    f"     as an alkoxide, say - the correction is a whole electron. That has\n"
                    f"     to be absorbed by replacing every charge in the group, not by\n"
                    f"     dumping the difference on one or two atoms. List the full group.\n"
                    f"\n"
                    f"  2. The molecule's PROTONATION STATE does not match `target`. An amine\n"
                    f"     drawn neutral that should be protonated (or the reverse) shifts the\n"
                    f"     net charge by a unit, and no charge on carbon stands in for a\n"
                    f"     missing hydrogen. Fix the structure and re-run CGenFF.")
            print(f"\nbalance       : residual {residual:+.4f} split over "
                  f"{len(balance)} atom(s), {share:+.4f} each")
            for name in balance:
                i = idx[name]
                f = lines[i].split()
                old_q = float(f[3])
                retype(i, "-", old_q + share)
                print(f"  {name:<6s} {f[2]:<7s} {old_q:+.3f}  ->  {old_q + share:+.3f}")

    final = sum(charge_of(i) for i in idx.values())
    print(f"\nfinal net     : {final:+.4f} over {len(idx)} atoms")
    if target is not None and abs(final - target) > 1e-3:
        sys.exit(f"ERROR: net {final:+.4f} does not match target {target:+.3f}")
    if abs(final - round(final)) > 1e-3:
        print(f"WARNING: net charge {final:+.4f} is not an integer - psfgen will "
              f"accept it but the system will not neutralise cleanly")

    # update the RESI header's stated charge if it carries one
    h = lines[start].split()
    if len(h) > 2:
        try:
            float(h[2])
            lines[start] = re.sub(r"(RESI\s+\S+\s+)[+-]?\d+\.?\d*",
                                  lambda m: f"{m.group(1)}{final:.3f}", lines[start], count=1)
        except ValueError:
            pass

    if dry:
        print("\n--dry-run: no file written")
        return
    backup = str_file + ".preoverride"
    if not os.path.exists(backup):          # keep the ORIGINAL import, not the last run
        shutil.copyfile(str_file, backup)
        print(f"\nwrote {str_file}  (original kept as {backup})")
    else:
        open(str_file, "w").write("".join(lines))
        print(f"\nwrote {str_file}  ({backup} already holds the original import)")
        return
    open(str_file, "w").write("".join(lines))

main()
