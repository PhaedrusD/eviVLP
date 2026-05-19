#!/usr/bin/env python3
"""
rename_atom_types_naa.py

Generic script to rename CGenFF atom types to CHARMM36 protein types
in the NAA .str file. Works for any amino acid and any incoming drug.

Usage: python3 rename_atom_types_naa.py <amino_acid> [str_file]
  amino_acid: lys, cys, tyr, ntrm
  str_file:   defaults to naa_clean.str (or ntrm_clean.str for ntrm)

Algorithm:
1. Rename atom types in the RESI topology section by atom NAME (not type),
   using a per-amino-acid mapping of protein atom names to CHARMM36 types.
2. Scan the renamed RESI to discover which CGenFF types appear on the drug
   side that ALSO appeared on the protein side (ambiguous types).
3. In the parameter sections, rename unambiguous protein types directly.
   For ambiguous types, generate all permutations so NAMD can find the
   right parameters regardless of which side the type came from.
4. Remove purely redundant lines (all protein-side types after renaming).
"""

import re
import sys
from itertools import product

# ============================================================
# Per-amino-acid mappings: atom_name -> CHARMM36_type
# These define which atoms are protein-side and their target types.
# ============================================================

AMINO_ACID_MAPS = {
    'lys': {
        # Backbone
        'N':   'NH1',
        'CA':  'CT1',
        'C':   'C',
        'O':   'O',
        'H':   'HN',
        'HA':  'HB1',
        # Side chain
        'CB':  'CT2',
        'CG':  'CT2',
        'CD':  'CT2',
        'CE':  'CT2',
        'NZ':  'NH2',
        'HZ':  'H',
        'HB1': 'HB2',
        'HB2': 'HB2',
        'HG1': 'HB2',
        'HG2': 'HB2',
        'HD1': 'HB2',
        'HD2': 'HB2',
        'HE1': 'HB2',
        'HE2': 'HB2',
        # Cap atoms (removed by noX, but included for completeness)
        'CX':  'C',
        'OX':  'O',
        'NX':  'NH2',
        'HX1': 'H',
        'HX2': 'H',
    },
    'cys': {
        # Backbone
        'N':   'NH1',
        'CA':  'CT1',
        'C':   'C',
        'O':   'O',
        'H':   'HN',
        'HA':  'HB1',
        # Side chain
        'CB':  'CT2',
        'SG':  'S',
        'HB1': 'HB2',
        'HB2': 'HB2',
        # Cap atoms
        'CX':  'C',
        'OX':  'O',
        'NX':  'NH2',
        'HX1': 'H',
        'HX2': 'H',
    },
    'tyr': {
        # Backbone
        'N':   'NH1',
        'CA':  'CT1',
        'C':   'C',
        'O':   'O',
        'H':   'HN',
        'HA':  'HB1',
        # Side chain
        'CB':  'CT2',
        'CG':  'CA',
        'CD1': 'CA',
        'CD2': 'CA',
        'CE1': 'CA',
        'CE2': 'CA',
        'CZ':  'CA',
        'OH':  'OH1',
        'HB1': 'HB2',
        'HB2': 'HB2',
        'HD1': 'HP',
        'HD2': 'HP',
        'HE1': 'HP',
        'HE2': 'HP',
        # Cap atoms
        'CX':  'C',
        'OX':  'O',
        'NX':  'NH2',
        'HX1': 'H',
        'HX2': 'H',
    },
    'ntrm': {
        'N':   'NH1',
        'CA':  'CT1',
        'C':   'C',
        'O':   'O',
        'HA':  'HB1',
        'CB':  'CT3',
        'HB1': 'HA3',
        'HB2': 'HA3',
        'HB3': 'HA3',
        'HN1': 'H',
    },
}


def parse_resi_atoms(lines):
    """
    Parse ATOM lines from the RESI section of a .str file.
    Returns list of (atom_name, atom_type) tuples.
    """
    atoms = []
    in_resi = False
    in_params = False
    for line in lines:
        stripped = line.strip()
        if 'read param' in stripped.lower():
            in_params = True
        if stripped.startswith('RESI'):
            in_resi = True
            continue
        if stripped == 'END' and in_resi and not in_params:
            in_resi = False
            continue
        if in_resi and not in_params and stripped.startswith('ATOM'):
            parts = stripped.split()
            if len(parts) >= 3:
                atoms.append((parts[1], parts[2]))  # (name, type)
    return atoms


def build_type_mapping(resi_atoms, protein_atom_map):
    """
    From the RESI atom list and protein atom name mapping:
    1. Build CGenFF->CHARMM36 type mapping for protein atoms
    2. Identify drug-side types
    3. Find ambiguous types (CGenFF types that appear on both sides)
    
    Returns: (cgenff_to_charmm, ambiguous_types, protein_types_after)
      cgenff_to_charmm: dict of CGenFF_type -> CHARMM36_type (from protein atoms)
      ambiguous_types: set of CGenFF types that also appear on drug atoms
      protein_types_after: set of CHARMM36 types assigned to protein atoms
    """
    # Build the CGenFF -> CHARMM36 mapping from protein atom names
    # Use a dict of sets since one CGenFF type can map to multiple CHARMM36 types
    # (e.g. HGP1 -> HN for backbone H, HGP1 -> H for cap HX1)
    cgenff_to_charmm = {}
    protein_cgenff_types = set()
    
    for atom_name, cgenff_type in resi_atoms:
        if atom_name in protein_atom_map:
            charmm_type = protein_atom_map[atom_name]
            if cgenff_type not in cgenff_to_charmm:
                cgenff_to_charmm[cgenff_type] = set()
            cgenff_to_charmm[cgenff_type].add(charmm_type)
            protein_cgenff_types.add(cgenff_type)
    
    # Find drug-side types
    drug_cgenff_types = set()
    for atom_name, cgenff_type in resi_atoms:
        if atom_name not in protein_atom_map:
            drug_cgenff_types.add(cgenff_type)
    
    # Ambiguous: types that appear on both protein and drug sides
    ambiguous_types = protein_cgenff_types & drug_cgenff_types
    
    protein_types_after = set()
    for charmm_types in cgenff_to_charmm.values():
        protein_types_after.update(charmm_types)
    
    return cgenff_to_charmm, ambiguous_types, protein_types_after


def rename_resi_by_name(lines, protein_atom_map):
    """
    Rename atom types in the RESI topology section by atom NAME.
    Only renames ATOM lines where the atom name is in protein_atom_map.
    Also removes IMPR lines referencing X atoms (HX, HX1, HX2, HX3).
    """
    new_lines = []
    in_resi = False
    in_params = False
    renamed_count = 0
    
    for line in lines:
        stripped = line.strip()
        
        if 'read param' in stripped.lower():
            in_params = True
        if stripped.startswith('RESI'):
            in_resi = True
        if stripped == 'END' and in_resi and not in_params:
            in_resi = False
        
        if in_resi and not in_params and stripped.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 3:
                atom_name = parts[1]
                if atom_name in protein_atom_map:
                    old_type = parts[2]
                    new_type = protein_atom_map[atom_name]
                    # Replace in-place preserving column width
                    pattern = re.compile(
                        r'(?:(?<=\s)|(?<=^))' + re.escape(old_type) + r'(?=\s)',
                        re.MULTILINE
                    )
                    new_line = pattern.sub(
                        lambda m: new_type.ljust(len(old_type)), line, count=1
                    )
                    if new_line != line:
                        renamed_count += 1
                    new_lines.append(new_line)
                    continue
        
        # Remove IMPR lines referencing X atoms
        if in_resi and not in_params and stripped.startswith('IMPR'):
            if any(x in stripped for x in ['HX', 'CX', 'OX', 'NX']):
                continue
        
        new_lines.append(line)
    
    print(f"  Topology: renamed {renamed_count} atom types by name")
    return new_lines


def get_type_fields_count(section):
    """Number of atom type fields at the start of a parameter line."""
    return {'BONDS': 2, 'ANGLES': 3, 'DIHEDRALS': 4, 'IMPROPERS': 4}.get(section, 0)


def generate_permutations(type_fields, cgenff_to_charmm, ambiguous_types, protein_types_after):
    """
    For a list of atom type fields, generate all valid permutations
    where ambiguous types can be either their CGenFF or CHARMM36 form(s).
    Non-ambiguous protein types are always renamed.
    Drug-only types are always kept.
    
    cgenff_to_charmm maps CGenFF_type -> set of CHARMM36_types
    
    Returns list of type field tuples.
    """
    options_per_field = []
    
    for t in type_fields:
        if t in cgenff_to_charmm:
            charmm_types = list(cgenff_to_charmm[t])
            if t in ambiguous_types:
                # Could be protein side (any of the CHARMM36 types) or drug side (keep original)
                options_per_field.append(charmm_types + [t])
            else:
                # Unambiguously protein side, use all possible CHARMM36 types
                options_per_field.append(charmm_types)
        else:
            # Drug-only type, keep as-is
            options_per_field.append([t])
    
    # Generate all combinations
    return list(product(*options_per_field))


def rename_types_in_param_line(line, type_fields, new_types):
    """
    Replace the atom type fields in a parameter line with new types,
    preserving column alignment.
    """
    new_line = line
    for old_t, new_t in zip(type_fields, new_types):
        if old_t != new_t:
            pattern = re.compile(
                r'(?:(?<=\s)|(?<=^))' + re.escape(old_t) + r'(?=\s)',
                re.MULTILINE
            )
            new_line = pattern.sub(
                lambda m, nt=new_t, ot=old_t: nt.ljust(len(ot)),
                new_line, count=1
            )
    return new_line


def process_parameters(lines, cgenff_to_charmm, ambiguous_types, protein_types_after):
    """
    Process parameter sections:
    - Rename unambiguous protein types
    - Generate permutations for ambiguous types
    - Remove purely redundant lines (all protein types after rename)
    - Keep drug-only lines as-is
    
    Returns (new_lines, all_extracted_params)
    """
    param_sections = {'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS'}
    current_section = None
    in_params = False
    new_lines = []
    extracted_params = {s: [] for s in param_sections}
    removed_count = 0
    added_count = 0
    
    for line in lines:
        stripped = line.strip()
        
        if 'read param' in stripped.lower():
            in_params = True
        
        if stripped in param_sections:
            current_section = stripped
            new_lines.append(line)
            continue
        
        if stripped == 'END' or stripped == 'RETURN':
            current_section = None
            new_lines.append(line)
            continue
        
        if not in_params or current_section not in param_sections:
            new_lines.append(line)
            continue
        
        if not stripped or stripped.startswith('!'):
            new_lines.append(line)
            continue
        
        parts = stripped.split()
        n = get_type_fields_count(current_section)
        if n == 0 or len(parts) < n + 2:
            new_lines.append(line)
            continue
        
        type_fields = parts[0:n]
        
        # Generate all valid permutations
        perms = generate_permutations(
            type_fields, cgenff_to_charmm, ambiguous_types, protein_types_after
        )
        
        # Filter out purely redundant permutations (all protein types)
        valid_perms = []
        for perm in perms:
            if all(t in protein_types_after for t in perm):
                continue  # Redundant with CHARMM36
            valid_perms.append(perm)
        
        if not valid_perms:
            removed_count += 1
            continue
        
        # First permutation replaces the original line
        first_line = rename_types_in_param_line(line, type_fields, valid_perms[0])
        new_lines.append(first_line)
        extracted_params[current_section].append(first_line)
        
        # Additional permutations are added as new lines
        for perm in valid_perms[1:]:
            extra_line = rename_types_in_param_line(line, type_fields, perm)
            new_lines.append(extra_line)
            extracted_params[current_section].append(extra_line)
            added_count += 1
        
    print(f"  Parameters: {removed_count} redundant lines removed, {added_count} permutation lines added")
    return new_lines, extracted_params


def process_str_file(str_file, protein_atom_map):
    """Full processing pipeline for a .str file."""
    with open(str_file, 'r') as f:
        lines = f.readlines()
    
    # Step 1: Parse RESI atoms to discover types
    resi_atoms = parse_resi_atoms(lines)
    print(f"  Found {len(resi_atoms)} atoms in RESI")
    
    # Step 2: Build type mapping and find ambiguous types
    cgenff_to_charmm, ambiguous_types, protein_types_after = build_type_mapping(
        resi_atoms, protein_atom_map
    )
    
    print(f"  CGenFF -> CHARMM36 mapping: {cgenff_to_charmm}")
    if ambiguous_types:
        print(f"  *** Ambiguous types (on both sides): {ambiguous_types}")
    else:
        print(f"  No ambiguous types found")
    
    # Step 3: Rename topology section by atom name
    lines = rename_resi_by_name(lines, protein_atom_map)
    
    # Step 4: Process parameter sections
    lines, extracted_params = process_parameters(
        lines, cgenff_to_charmm, ambiguous_types, protein_types_after
    )
    
    # Write updated .str file
    with open(str_file, 'w') as f:
        f.writelines(lines)
    
    print(f"  Updated {str_file}")
    return extracted_params


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 rename_atom_types_naa.py <amino_acid> [str_file]")
        print("  amino_acid: lys, cys, tyr, ntrm")
        sys.exit(1)
    
    aa = sys.argv[1].lower()
    
    if aa not in AMINO_ACID_MAPS:
        print(f"Error: Unknown amino acid '{aa}'")
        print(f"  Valid options: {', '.join(AMINO_ACID_MAPS.keys())}")
        sys.exit(1)
    
    # Default str file based on amino acid
    if aa == 'ntrm':
        default_str = 'ntrm_clean.str'
    else:
        default_str = 'naa_clean.str'
    
    str_file = sys.argv[2] if len(sys.argv) > 2 else default_str
    
    protein_atom_map = AMINO_ACID_MAPS[aa]
    
    print(f"Processing {str_file} for amino acid: {aa}")
    print(f"  Protein atom count: {len(protein_atom_map)}")
    
    extracted_params = process_str_file(str_file, protein_atom_map)
    
    total = sum(len(v) for v in extracted_params.values())
    print(f"\n  Total parameters to extract: {total}")
    for section in ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS']:
        count = len(extracted_params.get(section, []))
        if count:
            print(f"    {section}: {count}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()