#!/usr/bin/env python3
"""
add_npc_params_to_prm.py

Extracts cross-force-field parameters from ntrm_clean.str (already renamed)
and inserts them into the correct BONDS, ANGLES, DIHEDRALS, and IMPROPERS
sections of par_GUI.prm.

Cross-FF parameters are lines containing at least one CHARMM36 protein-side
atom type AND at least one CGenFF drug-side atom type.
Lines where ALL types are protein-side are redundant and skipped.
Lines where ALL types are drug-side are CGenFF-internal and skipped.
"""

# CHARMM36 protein-side types (after renaming by rename_atom_types_ntrm.py)
PROTEIN_TYPES = {'NH1', 'CT1', 'C', 'CT2', 'H', 'HB1', 'HB3', 'O'}

PARAM_SECTIONS = {'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS'}


def get_num_type_fields(section):
    """Number of atom type fields at the start of a parameter line."""
    return {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}.get(section, 0)


def extract_cross_ff_params(str_file):
    """
    Read the parameter section of ntrm_clean.str and extract
    cross-force-field parameter lines, grouped by section.
    """
    cross_params = {s: [] for s in PARAM_SECTIONS}

    with open(str_file, 'r') as f:
        lines = f.readlines()

    current_section = None
    in_params = False

    for line in lines:
        stripped = line.strip()

        # Detect start of parameter block
        if 'read param' in stripped.lower():
            in_params = True
            continue

        if not in_params:
            continue

        # Track sections
        if stripped in PARAM_SECTIONS:
            current_section = stripped
            continue

        if stripped == 'END' or stripped == 'RETURN':
            current_section = None
            continue

        # Skip empty lines and comments
        if not stripped or stripped.startswith('!'):
            continue

        if current_section not in PARAM_SECTIONS:
            continue

        # Get atom type fields
        parts = stripped.split()
        n = get_num_type_fields(current_section)
        if len(parts) < n + 2:
            continue

        type_fields = parts[0:n]

        has_protein = any(t in PROTEIN_TYPES for t in type_fields)
        all_protein = all(t in PROTEIN_TYPES for t in type_fields)

        if not all_protein:
            # Cross-FF or drug-side parameter - both needed
            cross_params[current_section].append(line)
            
    return cross_params


def insert_into_prm(prm_file, cross_params):
    """
    Insert cross-FF parameters into the correct sections of par_GUI.prm.
    Each section's parameters are appended at the end of that section,
    just before the next section header or CMAP/NONBONDED/END.
    """
    with open(prm_file, 'r') as f:
        lines = f.readlines()

    # Find the line indices of each section header and the CMAP line
    section_starts = {}
    cmap_line = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped in PARAM_SECTIONS:
            section_starts[stripped] = i
        if stripped == 'CMAP':
            cmap_line = i

    # Determine insertion points: end of each section
    # (just before the next section or CMAP)
    ordered_sections = ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']
    section_ends = {}

    for idx, section in enumerate(ordered_sections):
        if section not in section_starts:
            continue
        start = section_starts[section]

        # Find where this section ends: next section header or CMAP
        end = len(lines)
        for next_section in ordered_sections[idx + 1:]:
            if next_section in section_starts:
                end = section_starts[next_section]
                break
        if cmap_line is not None and cmap_line < end:
            end = cmap_line

        # Walk backward from end to find last non-empty line in this section
        insert_at = end
        section_ends[section] = insert_at

    # Build the new file content, inserting at each section end
    # Work backwards so line indices don't shift
    new_lines = list(lines)
    total_added = 0

    for section in reversed(ordered_sections):
        params = cross_params.get(section, [])
        if not params or section not in section_ends:
            continue

        insert_at = section_ends[section]
        block = [f"! === NPC cross-force-field {section} (N-terminus) ===\n"]
        for p in params:
            if not p.endswith('\n'):
                p += '\n'
            block.append(p)
        block.append('\n')

        new_lines[insert_at:insert_at] = block
        total_added += len(params)

    with open(prm_file, 'w') as f:
        f.writelines(new_lines)

    return total_added, cross_params


def main():
    str_file = 'ntrm_clean.str'
    prm_file = 'par_GUI.prm'

    print(f"Extracting cross-FF parameters from {str_file}...")
    cross_params = extract_cross_ff_params(str_file)

    for section in ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']:
        count = len(cross_params[section])
        if count:
            print(f"  {section}: {count} cross-FF parameters found")

    print(f"\nInserting into {prm_file}...")
    total, _ = insert_into_prm(prm_file, cross_params)
    print(f"  Total parameters added: {total}")
    print("Done!")


if __name__ == "__main__":
    main()
