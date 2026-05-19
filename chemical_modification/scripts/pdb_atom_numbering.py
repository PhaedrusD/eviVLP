import re

def rename_pdb_atoms(input_pdb, output_pdb):
    atom_counts = {}  # Dictionary to track occurrences of atom names
    new_lines = []

    with open(input_pdb, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                
                # Generate a unique index for each atom type
                key = (res_name, atom_name)
                atom_counts[key] = atom_counts.get(key, 0) + 1
                
                # Create new unique atom name
                new_atom_name = f"{atom_name}{atom_counts[key]}"
                new_atom_name = new_atom_name[:4]  # Ensure it fits PDB format
                
                # Replace the atom name in the line
                line = line[:12] + new_atom_name.ljust(4) + line[16:]
            
            new_lines.append(line)
    
    # Write updated PDB file
    with open(output_pdb, 'w') as f:
        f.writelines(new_lines)

# Example usage:
rename_pdb_atoms("zif_8.pdb", "zif_8_fix.pdb")
