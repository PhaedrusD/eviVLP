#!/bin/bash
# eviVLP one-command setup for the chemical-modification GUI.
#
#   git clone https://github.com/PhaedrusD/eviVLP.git
#   cd eviVLP
#   ./setup.sh            # run folder defaults to ~/eviVLP_run
#   ./setup.sh ~/my_dir   # or choose one
#
# What it does: checks the build prerequisites (installing them with Homebrew
# on macOS when it can), lays out the flat run folder the GUI needs, builds the
# GUI and the IC generators there, and reports where VMD and NAMD were found.
# Re-running it is safe: files you edited in the run folder are never
# overwritten (sync_repo.py skips a newer destination).
set -e
RUN_DIR="${1:-$HOME/eviVLP_run}"
REPO="$(cd "$(dirname "$0")" && pwd)"
OS="$(uname -s)"
ok()   { printf '  [ok]   %s\n' "$*"; }
warn() { printf '  [!!]   %s\n' "$*"; }
need_brew() {
    if command -v brew >/dev/null 2>&1; then return 0; fi
    warn "Homebrew not found - install it from https://brew.sh then re-run, or install '$1' yourself."
    return 1
}

echo "eviVLP setup"
echo "  repo:       $REPO"
echo "  run folder: $RUN_DIR"
echo
echo "1. Build prerequisites"
if [ "$OS" = "Darwin" ] && ! xcode-select -p >/dev/null 2>&1; then
    warn "Xcode command-line tools missing; running 'xcode-select --install' (finish the dialog, then re-run)."
    xcode-select --install || true
    exit 1
fi
command -v g++ >/dev/null 2>&1 && ok "g++ $(g++ --version | head -1)" || { warn "g++ not found"; exit 1; }

for tool in wx-config:wxwidgets obabel:open-babel python3:python; do
    exe="${tool%%:*}"; pkg="${tool##*:}"
    if command -v "$exe" >/dev/null 2>&1; then
        ok "$exe ($(command -v "$exe"))"
    elif [ "$OS" = "Darwin" ] && need_brew "$pkg"; then
        echo "  installing $pkg with Homebrew..."
        brew install "$pkg"
    else
        warn "$exe not found - install package '$pkg' and re-run."; exit 1
    fi
done
if python3 -c "import numpy" >/dev/null 2>&1; then
    ok "numpy $(python3 -c 'import numpy; print(numpy.__version__)')"
else
    echo "  installing numpy..."
    python3 -m pip install --user numpy || python3 -m pip install --user --break-system-packages numpy
fi

echo
echo "2. Run folder (flat layout of chemical_modification/{src,scripts,tcl,data})"
python3 "$REPO/sync_repo.py" --pull --module chemical_modification --gui "$RUN_DIR" | tail -4

echo
echo "3. Build"
if ( cd "$RUN_DIR" && make -s all ); then
    ok "GUI, IC_naa, IC_ntrm_NPC, IC_table_one built in $RUN_DIR"
else
    warn "build failed - see the compiler output above. Nothing else to do until it builds."
    exit 1
fi

echo
echo "4. Simulation programs (the GUI looks these up at startup; override with EVIVLP_VMD / EVIVLP_NAMD)"
VMD="${EVIVLP_VMD:-$(command -v vmd || ls -d /Applications/VMD*.app/Contents/MacOS/startup.command.csh 2>/dev/null | tail -1)}"
NAMD="${EVIVLP_NAMD:-$(command -v namd3 || ls "$HOME"/Downloads/NAMD*/namd3 /Applications/NAMD*/namd3 2>/dev/null | tail -1)}"
[ -n "$VMD"  ] && ok "VMD:  $VMD"   || warn "VMD not found. Install VMD 1.9.4+ into /Applications (https://www.ks.uiuc.edu/Research/vmd/) or export EVIVLP_VMD=/path/to/vmd"
[ -n "$NAMD" ] && ok "NAMD: $NAMD"  || warn "NAMD not found. Unpack NAMD 3.x into ~/Downloads (https://www.ks.uiuc.edu/Research/namd/) or export EVIVLP_NAMD=/path/to/namd3"

echo
echo "Done. Start the GUI with:"
echo "  cd \"$RUN_DIR\" && ./GUI"
echo "Test input: $REPO/case_studies/case_study1_TMV_TEMPO/inputs/TEMPO.pdb (override file nitroxide_overrides.prm beside it)."
