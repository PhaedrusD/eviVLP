#!/usr/bin/env python3
"""
rename_atom_types_ntrm.py

1. Renames CGenFF atom types to CHARMM36 protein types in ntrm_clean.str
   (both the RESI topology section and parameter sections).
2. After renaming, extracts cross-force-field parameters (lines containing
   at least one protein-side type AND at least one drug-side type) and
   appends them to par_GUI.prm.
3. Removes purely protein-side parameter lines (redundant with CHARMM36).

Uses in-place word replacement to preserve column alignment.
"""

import re

# CGenFF -> CHARMM36 mapping for protein-side atoms
TYPE_MAP = {
    'NG311':  'NH1',
    'CG311':  'CT1',
    'CG2O4':  'C',
    'CG331':  'CT3',
    'HGPAM1': 'H',
    'HGA1':   'HB1',
    'HGA3':   'HA3',
    'OG2D1':  'O',
}

# The set of CHARMM36 protein types (after renaming)
PROTEIN_TYPES = set(TYPE_MAP.values())


def rename_all_types(line):
    """
    Apply all type renames to a line using in-place regex substitution.
    Preserves surrounding whitespace by padding shorter replacements
    to match the original token width.
    Sorted by length descending to avoid partial matches.
    """
    new_line = line
    for old_type in sorted(TYPE_MAP.keys(), key=len, reverse=True):
        new_type = TYPE_MAP[old_type]
        # Match whole word only (at start of line or after whitespace, followed by whitespace)
        pattern = re.compile(r'(?:(?<=\s)|(?<=^))' + re.escape(old_type) + r'(?=\s)', re.MULTILINE)
        
        def make_replacer(ot, nt):
            def replacer(match):
                # Pad new type to same width as old type
                return nt.ljust(len(ot))
            return replacer
        
        new_line = pattern.sub(make_replacer(old_type, new_type), new_line)
    return new_line


def has_type(line, type_set):
    """Check if a line contains any atom type from the given set as a whole word."""
    words = line.split()
    return any(w in type_set for w in words)


def get_type_fields(section, line):
    """Return the atom type tokens from a parameter line."""
    parts = line.split()
    if not parts or parts[0].startswith('!'):
        return None
    
    n = {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}.get(section, 0)
    if n == 0 or len(parts) < n + 2:
        return None
    return parts[0:n]


def classify_param_line(section, line):
    """
    Classify a parameter line after renaming.
    Returns: 'cross', 'redundant', 'drug_only', or None.
    """
    type_fields = get_type_fields(section, line)
    if type_fields is None:
        return None
    
    has_protein = any(t in PROTEIN_TYPES for t in type_fields)
    all_protein = all(t in PROTEIN_TYPES for t in type_fields)
    
    if all_protein:
        return 'redundant'
    elif has_protein:
        return 'cross'
    else:
        return 'drug_only'


def process_str_file(str_file):
    """Process the .str file: rename types, identify cross-FF params."""
    with open(str_file, 'r') as f:
        lines = f.readlines()

    new_lines = []
    cross_ff_params = {'BONDS': [], 'ANGLES': [], 'DIHEDRALS': [], 'IMPROPERS': []}
    param_sections = {'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS'}
    current_section = None
    in_resi = False
    in_params = False

    renamed_count = 0
    removed_count = 0

    for line in lines:
        stripped = line.strip()

        # Track RESI block
        if stripped.startswith('RESI'):
            in_resi = True

        # Detect "read param" to know we're in the parameter half
        if 'read param' in stripped.lower():
            in_params = True

        # Track sections
        if stripped in param_sections:
            current_section = stripped
            new_lines.append(line)
            continue

        # Handle END
        if stripped == 'END':
            if in_resi and not in_params:
                in_resi = False
            current_section = None
            new_lines.append(line)
            continue

        # Handle RETURN
        if stripped == 'RETURN':
            new_lines.append(line)
            continue

        # In RESI topology block - rename atom types in ATOM lines
        if in_resi and not in_params:
            if stripped.startswith('ATOM'):
                new_line = rename_all_types(line)
                if new_line != line:
                    renamed_count += 1
                new_lines.append(new_line)
                continue
            # Handle IMPR lines that reference removed HX atom
            elif stripped.startswith('IMPR'):
                if 'HX' in stripped:
                    removed_count += 1
                    continue
            new_lines.append(line)
            continue

        # In parameter sections - rename and classify
        if current_section in param_sections:
            if not stripped or stripped.startswith('!'):
                new_lines.append(line)
                continue

            # Rename types in place
            new_line = rename_all_types(line)
            if new_line != line:
                renamed_count += 1

            # Classify the renamed line
            classification = classify_param_line(current_section, new_line)

            if classification == 'redundant':
                removed_count += 1
                continue
            elif classification == 'cross':
                cross_ff_params[current_section].append(new_line)
                new_lines.append(new_line)
            else:
                new_lines.append(new_line)
            continue

        new_lines.append(line)

    # Write updated .str file
    with open(str_file, 'w') as f:
        f.writelines(new_lines)

    print(f"Atom type renaming complete in {str_file}")
    print(f"  Lines with types renamed: {renamed_count}")
    print(f"  Redundant parameter lines removed: {removed_count}")

    return cross_ff_params


def main():
    str_file = 'ntrm_clean.str'
    prm_file = 'par_GUI.prm'

    print(f"Processing {str_file}...")
    print(f"Type mapping: {TYPE_MAP}\n")

    cross_ff_params = process_str_file(str_file)


    print("\nDone!")


if __name__ == "__main__":
    main()