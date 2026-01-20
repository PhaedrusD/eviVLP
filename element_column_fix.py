#!/usr/bin/env python3
"""
PDB Element Column Fixer
Automatically adds missing element symbols in PDB files
"""

import sys
import os

def fix_pdb_element_column(input_file, output_file=None):
    """
    Read a PDB file and fix missing element symbols in the last column.
    If no element symbol exists, use the first character of the atom name.
    """
    
    fixed_lines = []
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Only process ATOM and HETATM records
            if line.startswith(('ATOM', 'HETATM')):
                # Split the line into components
                parts = line.split()
                
                # Check if we have at least 11 columns (indicating possible missing element)
                if len(parts) >= 11:
                    # Reconstruct the fixed line
                    fixed_line = fix_atom_line(line)
                    fixed_lines.append(fixed_line)
                else:
                    # Keep line as is if it doesn't meet expected format
                    fixed_lines.append(line)
            else:
                # Keep non-ATOM/HETATM lines as is
                fixed_lines.append(line)
    
    # Write the fixed file
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"Fixed PDB file saved as: {output_file}")
    return output_file

def fix_atom_line(line):
    """
    Fix a single ATOM/HETATM line by ensuring it has an element symbol
    """
    # Standard PDB format positions
    record = line[0:6].strip()
    serial = line[6:11].strip()
    atom_name = line[12:16].strip()
    alt_loc = line[16:17].strip()
    res_name = line[17:20].strip()
    chain_id = line[21:22].strip()
    res_seq = line[22:26].strip()
    i_code = line[26:27].strip()
    x = line[30:38].strip()
    y = line[38:46].strip()
    z = line[46:54].strip()
    occupancy = line[54:60].strip()
    temp_factor = line[60:66].strip()
    segment = line[72:76].strip() if len(line) > 72 else ""
    element = line[76:78].strip() if len(line) > 76 else ""
    charge = line[78:80].strip() if len(line) > 78 else ""
    
    # If element is missing, derive it from atom name
    if not element:
        # Get first character of atom name, strip digits
        element_candidate = ''.join([c for c in atom_name if not c.isdigit()])
        if element_candidate:
            element = element_candidate[0].upper()
        else:
            element = "X"  # Fallback
    
    # Reconstruct the line with proper formatting
    fixed_line = f"{record:6s}{serial:>5s} {atom_name:4s}{alt_loc:1s}{res_name:3s} {chain_id:1s}{res_seq:>4s}{i_code:1s}   {x:>8s}{y:>8s}{z:>8s}{occupancy:>6s}{temp_factor:>6s}      {segment:4s}{element:>2s}{charge:2s}\n"
    
    return fixed_line

def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_pdb_elements.py input.pdb [output.pdb]")
        print("If output is not specified, creates naa.pdb")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    try:
        fix_pdb_element_column(input_file, output_file)
        print("PDB file fixed successfully!")
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
