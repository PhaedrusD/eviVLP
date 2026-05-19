#!/usr/bin/env python3
"""
add_naa_params_to_prm.py

Extracts non-redundant parameters from a renamed .str file and inserts
them into the correct BONDS, ANGLES, DIHEDRALS, and IMPROPERS sections
of par_GUI.prm.

Usage: python3 add_naa_params_to_prm.py <str_file> [prm_file]
  str_file: the renamed .str file (e.g. nad.str, ntrm_clean.str)
  prm_file: defaults to par_GUI.prm
"""

import sys

# CHARMM36 protein-side types (after renaming)
# Parameters where ALL types are in this set are redundant
PROTEIN_TYPES = {
    'NH1', 'NH2', 'CT1', 'CT2', 'CT3', 'C', 'O', 'H', 'HN',
    'HB1', 'HB2', 'HB3', 'HA3', 'HP', 'CA', 'OH1', 'S',
}

PARAM_SECTIONS = {'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS'}


def get_num_type_fields(section):
    return {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}.get(section, 0)


def extract_params(str_file):
    """Extract all non-redundant parameters from a .str file."""
    params = {s: [] for s in PARAM_SECTIONS}

    with open(str_file, 'r') as f:
        lines = f.readlines()

    current_section = None
    in_params = False

    for line in lines:
        stripped = line.strip()

        if 'read param' in stripped.lower():
            in_params = True
            continue

        if not in_params:
            continue

        if stripped in PARAM_SECTIONS:
            current_section = stripped
            continue

        if stripped == 'END' or stripped == 'RETURN':
            current_section = None
            continue

        if not stripped or stripped.startswith('!'):
            continue

        if current_section not in PARAM_SECTIONS:
            continue

        parts = stripped.split()
        n = get_num_type_fields(current_section)
        if len(parts) < n + 2:
            continue

        type_fields = parts[0:n]

        # Skip if ALL types are protein-side (redundant with CHARMM36)
        if all(t in PROTEIN_TYPES for t in type_fields):
            continue

        params[current_section].append(line)

    return params


def insert_into_prm(prm_file, params, label):
    """Insert parameters into the correct sections of par_GUI.prm."""
    with open(prm_file, 'r') as f:
        lines = f.readlines()

    section_starts = {}
    cmap_line = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped in PARAM_SECTIONS:
            section_starts[stripped] = i
        if stripped == 'CMAP':
            cmap_line = i

    ordered_sections = ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']
    section_ends = {}

    for idx, section in enumerate(ordered_sections):
        if section not in section_starts:
            continue
        end = len(lines)
        for next_section in ordered_sections[idx + 1:]:
            if next_section in section_starts:
                end = section_starts[next_section]
                break
        if cmap_line is not None and cmap_line < end:
            end = cmap_line
        section_ends[section] = end

    new_lines = list(lines)
    total_added = 0

    for section in reversed(ordered_sections):
        section_params = params.get(section, [])
        if not section_params or section not in section_ends:
            continue

        insert_at = section_ends[section]
        block = [f"\n! === {label} cross-force-field {section} ===\n"]
        for p in section_params:
            if not p.endswith('\n'):
                p += '\n'
            block.append(p)

        new_lines[insert_at:insert_at] = block
        total_added += len(section_params)

    with open(prm_file, 'w') as f:
        f.writelines(new_lines)

    return total_added


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 add_naa_params_to_prm.py <str_file> [prm_file] [label]")
        sys.exit(1)

    str_file = sys.argv[1]
    prm_file = sys.argv[2] if len(sys.argv) > 2 else 'par_GUI.prm'
    label = sys.argv[3] if len(sys.argv) > 3 else 'NAA'

    print(f"Extracting parameters from {str_file}...")
    params = extract_params(str_file)

    for section in ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']:
        count = len(params[section])
        if count:
            print(f"  {section}: {count}")

    print(f"\nInserting into {prm_file}...")
    total = insert_into_prm(prm_file, params, label)
    print(f"  Total parameters added: {total}")
    print("Done!")


if __name__ == "__main__":
    main()
