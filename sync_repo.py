"""Keep the flat GUI working folder and the structured repo in step.

WHY THIS EXISTS
  The two trees are not duplicates, they are two layouts of the same code.
  The GUI shells out to `python3 rename_atom_types_naa.py` with no path, so
  every script and data file it touches has to sit beside the executable in one
  flat folder. The repo is organised (src/, scripts/, tcl/, data/) because that
  is what reads properly for publication. They cannot be the same directory
  without flattening the repo or teaching the GUI about subdirectories, so a
  file edited in one place has to reach the other - and when it silently does
  not, you rebuild and your change isn't there.

  Run it bare and it reports drift and says which side is newer. That check is
  the point; copying is secondary.

WHAT IT CONSIDERS
  The repo is the manifest. Every file under a mapped module directory is
  compared against the same basename in the flat folder. Anything in the flat
  folder that is NOT in the repo is ignored - that is how the working folder's
  scratch and retired sources stay out of the published tree. Pass --untracked
  to list them anyway.

GENERATED FILES
  A repo file is treated as generated, and its drift ignored, when a
  `<stem>_template<ext>` sibling exists beside it. That one rule covers the
  per-run outputs: par_GUI.prm, top_all36_cgenff_CBD.rtf, and the
  orient_place_drug_*/patch_maker_* scripts the GUI writes from templates.

SAFETY
  Nothing is ever deleted. A copy is skipped when the destination is NEWER than
  the source, because overwriting newer work is the failure this tool exists to
  prevent; --force overrides that, deliberately.

usage:
  python3 sync_repo.py                 report drift, exit 1 if any
  python3 sync_repo.py --push          flat working folder -> repo
  python3 sync_repo.py --pull          repo -> flat working folder
  python3 sync_repo.py --untracked     also list working-folder-only files
  python3 sync_repo.py --gui PATH      point at a different working folder
  --force                              copy even when the destination is newer
  --dry-run                            with --push/--pull, say what would happen
"""
import sys, os, shutil, filecmp

HERE = os.path.dirname(os.path.abspath(__file__))

# repo module -> flat working folder, relative to the repo root
MODULES = [
    ("chemical_modification", os.path.join("..", "eviVLP", "GUI")),
    ("genetic_modification",  os.path.join("..", "eviVLP", "evi_gen_GUI")),
]
SUBDIRS = ["src", "scripts", "tcl", "data"]
SKIP_EXT = {".o", ".pyc", ".dSYM", ".bak"}


def norm(path):
    """File bytes with trailing whitespace on each line removed, for telling a
    real change from a stray newline."""
    try:
        with open(path, "rb") as fh:
            return b"\n".join(l.rstrip() for l in fh.read().split(b"\n")).rstrip()
    except OSError:
        return None


def is_generated(repo_path):
    d, f = os.path.split(repo_path)
    stem, ext = os.path.splitext(f)
    return os.path.exists(os.path.join(d, stem + "_template" + ext))


def collect(repo_root, module, gui_rel):
    gui_dir = os.path.normpath(os.path.join(repo_root, gui_rel))
    if not os.path.isdir(gui_dir):
        return None, gui_dir, []
    rows, by_base = [], {}
    for sub in SUBDIRS:
        d = os.path.join(repo_root, module, sub)
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if not os.path.isfile(p) or os.path.splitext(f)[1] in SKIP_EXT or f.startswith("."):
                continue
            by_base.setdefault(f, []).append(os.path.join(module, sub, f))
            rows.append((os.path.join(module, sub, f), p, os.path.join(gui_dir, f), f))
    collisions = {b: v for b, v in by_base.items() if len(v) > 1}
    return rows, gui_dir, collisions


def classify(repo_p, gui_p):
    if not os.path.exists(gui_p):
        return "missing in working folder", None
    if filecmp.cmp(repo_p, gui_p, shallow=False):
        return "in sync", None
    if norm(repo_p) == norm(gui_p):
        return "whitespace only", None
    rm, gm = os.path.getmtime(repo_p), os.path.getmtime(gui_p)
    return ("working folder newer" if gm > rm else "repo newer"), (rm, gm)


def main():
    argv = sys.argv[1:]
    push, pull = "--push" in argv, "--pull" in argv
    force, dry, untracked = "--force" in argv, "--dry-run" in argv, "--untracked" in argv
    if push and pull:
        sys.exit("--push and --pull are opposites; pick one")
    gui_override = None
    if "--gui" in argv:
        i = argv.index("--gui")
        if i + 1 >= len(argv): sys.exit("--gui needs a path")
        gui_override = argv[i + 1]

    drift = 0
    for module, gui_rel in MODULES:
        if gui_override: gui_rel = gui_override
        rows, gui_dir, collisions = collect(HERE, module, gui_rel)
        if rows is None:
            print(f"\n{module}  ->  {gui_dir}\n  working folder not found, skipped")
            continue
        print(f"\n{module}  <->  {gui_dir}")
        if not rows:
            print("  no mapped files"); continue

        for base, paths in collisions.items():
            print(f"  !! {base}: {len(paths)} repo files share this name, so the flat")
            print(f"     folder cannot hold both - {', '.join(paths)}")
            drift += 1

        buckets = {}
        for rel, repo_p, gui_p, base in rows:
            if base in collisions:
                continue              # already reported above; do not classify twice
            if is_generated(repo_p):
                buckets.setdefault("generated (ignored)", []).append((rel, repo_p, gui_p))
                continue
            state, times = classify(repo_p, gui_p)
            buckets.setdefault(state, []).append((rel, repo_p, gui_p, times))

        for state in ("working folder newer", "repo newer", "missing in working folder",
                      "whitespace only", "in sync", "generated (ignored)"):
            items = buckets.get(state)
            if not items: continue
            if state in ("in sync", "generated (ignored)"):
                print(f"  {state:<28} {len(items)} file(s)")
                continue
            print(f"  {state}:")
            for item in items:
                rel, repo_p, gui_p = item[0], item[1], item[2]
                print(f"    {os.path.basename(rel):<38} {rel}")
                if state in ("working folder newer", "repo newer", "missing in working folder"):
                    drift += 1
                # act
                if state == "whitespace only":
                    continue
                src, dst = (None, None)
                if push and state in ("working folder newer",):
                    src, dst = gui_p, repo_p
                elif push and state == "missing in working folder":
                    print("      (nothing in the working folder to push)"); continue
                elif pull and state in ("repo newer", "missing in working folder"):
                    src, dst = repo_p, gui_p
                elif (push and state == "repo newer") or (pull and state == "working folder newer"):
                    if not force:
                        print("      SKIPPED: destination is newer; re-run with --force to overwrite")
                        continue
                    src, dst = (gui_p, repo_p) if push else (repo_p, gui_p)
                if src:
                    if dry:
                        print(f"      would copy -> {dst}")
                    else:
                        os.makedirs(os.path.dirname(dst), exist_ok=True)
                        shutil.copy2(src, dst)
                        print(f"      copied -> {dst}")

        if untracked:
            tracked = {base for _, _, _, base in rows}
            extra = sorted(f for f in os.listdir(gui_dir)
                           if os.path.isfile(os.path.join(gui_dir, f))
                           and f not in tracked and not f.startswith(".")
                           and os.path.splitext(f)[1] not in SKIP_EXT)
            print(f"  working-folder-only ({len(extra)} file(s), not in the repo):")
            for f in extra: print(f"    {f}")
        else:
            tracked = {base for _, _, _, base in rows}
            extra = [f for f in os.listdir(gui_dir)
                     if os.path.isfile(os.path.join(gui_dir, f)) and f not in tracked
                     and not f.startswith(".") and os.path.splitext(f)[1] not in SKIP_EXT]
            print(f"  working-folder-only            {len(extra)} file(s)  (--untracked to list)")

    print()
    if push or pull:
        print("done." if not dry else "dry run: nothing written.")
        return 0
    if drift:
        print(f"{drift} file(s) out of step. --pull to update the working folder, "
              f"--push to update the repo.")
        return 1
    print("repo and working folder are in step.")
    return 0


sys.exit(main())
