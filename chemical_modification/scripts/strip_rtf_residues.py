#!/usr/bin/env python3
"""strip_rtf_residues.py NAME [NAME ...]
Remove existing RESI/PRES block(s) with these names from top_all36_cgenff_CBD.rtf
so a step that appends them can be re-run. A block runs from its RESI/PRES line
to the line before the next RESI / PRES / END. Case-insensitive.
"""
import re, sys
RTF = 'top_all36_cgenff_CBD.rtf'
names = {n.upper() for n in sys.argv[1:]}
out, skipping, removed = [], False, 0
for l in open(RTF).read().split('\n'):
    m = re.match(r'^\s*(RESI|PRES)\s+(\S+)', l, re.I)
    if m or re.match(r'^\s*END\b', l, re.I):
        skipping = bool(m) and m.group(2).upper() in names
        if skipping: removed += 1
    if not skipping:
        out.append(l)
open(RTF, 'w').write('\n'.join(out))
print(f"{RTF}: removed {removed} previous block(s) for {sorted(names)}")
