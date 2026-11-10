#!/bin/sh
# check_tools.sh <vmd launcher> <namd binary>
#
# macOS only. Reports whether the VMD and NAMD that eviVLP found can actually
# run on this Mac, so a wrong download is caught at setup / GUI startup instead
# of as "Bad CPU type in executable" halfway through a run. Two things go wrong
# in practice:
#   * the build is for the other chip (Intel x86_64 vs Apple-silicon arm64);
#     an Intel build runs on Apple silicon only if Rosetta is installed, and an
#     arm64 build never runs on an Intel Mac;
#   * macOS has quarantined a program downloaded from the web and will refuse
#     to launch it from a script.
# Prints one line per finding; lines starting "WARNING:" need action.
# Called by setup.sh and by the GUI at startup. Exit status is always 0.

[ "$(uname -s)" = "Darwin" ] || exit 0

mach=$(uname -m)
# A process running under Rosetta reports x86_64; the hardware is still arm64.
[ "$(sysctl -n sysctl.proc_translated 2>/dev/null)" = "1" ] && mach=arm64

have_rosetta() { "${EVIVLP_ARCH_BIN:-/usr/bin/arch}" -x86_64 /usr/bin/true >/dev/null 2>&1; }

check_arch() {   # label, Mach-O binary
    label=$1; bin=$2
    desc=$(file -b "$bin" 2>/dev/null)
    case "$desc" in
        *Mach-O*) ;;
        *) return 0 ;;                      # a script or unreadable: nothing to say
    esac
    case "$desc" in
        *"$mach"*)
            echo "$label: built for this Mac ($mach) - OK" ;;
        *x86_64*)
            if have_rosetta; then
                echo "$label: Intel (x86_64) build, runs through Rosetta - OK (the arm64 build would be faster)"
            else
                echo "WARNING: $label is an Intel (x86_64) build and Rosetta is not installed, so it will fail with 'Bad CPU type in executable'. Either run:  softwareupdate --install-rosetta --agree-to-license   or download the Apple-silicon (arm64) build instead."
            fi ;;
        *arm64*)
            echo "WARNING: $label is an Apple-silicon (arm64) build but this Mac has an Intel processor, so it will fail with 'Bad CPU type in executable'. Download the Intel (x86_64) build instead." ;;
    esac
}

check_quarantine() {   # label, path to de-quarantine
    label=$1; target=$2
    if xattr -p com.apple.quarantine "$target" >/dev/null 2>&1; then
        echo "WARNING: macOS has quarantined $label (it was downloaded from the web) and may refuse to launch it. Run:  xattr -dr com.apple.quarantine \"$target\""
    fi
}

vmd=$1; namd=$2

if [ -n "$vmd" ] && [ -e "$vmd" ]; then
    case "$vmd" in
        */Contents/MacOS/*)                 # the launcher inside VMD.app: find the real binary
            app=${vmd%/Contents/MacOS/*}
            bin=$(find "$app/Contents" -maxdepth 3 -type f -name 'vmd_MACOSX*' 2>/dev/null | head -1)
            [ -n "$bin" ] && check_arch "VMD" "$bin"
            check_quarantine "VMD" "$app" ;;
        *)
            check_arch "VMD" "$vmd" ;;      # no-op if it is a wrapper script
    esac
fi

if [ -n "$namd" ] && [ -e "$namd" ]; then
    check_arch "NAMD" "$namd"
    check_quarantine "NAMD" "$(dirname "$namd")"
fi
exit 0
