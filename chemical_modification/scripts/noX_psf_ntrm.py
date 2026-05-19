def process_str_files(pdb_file='ntrm_fixed.pdb', input_file='ntrm.str', output_file='ntrm_clean.str'):
    try:
        # Step 1: Get atom indices (1-based) and names containing 'X' from the pdb file
        print("Step 1: Reading naa.pdb for atom indices containing 'X'...")
        x_atom_indices = set()
        x_atom_names_pdb = set()
        atom_index = 0
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_index += 1
                    atom_name = line[12:16].strip()
                    if 'X' in atom_name:
                        x_atom_indices.add(atom_index)
                        x_atom_names_pdb.add(atom_name)
        
        print(f"  Found {len(x_atom_indices)} atoms with 'X' at indices: {sorted(x_atom_indices)}")

        # Step 2: Get the corresponding str atom names at those indices
        print("\nStep 2: Reading naa.str for atom names at those indices...")
        x_atom_names_str = set()
        atom_index = 0
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    atom_index += 1
                    if atom_index in x_atom_indices:
                        parts = line.split()
                        if len(parts) >= 2:
                            x_atom_names_str.add(parts[1])

        print(f"  Corresponding str atom names to remove: {x_atom_names_str}")

        # Step 3: Process str file, skipping flagged ATOM lines and BOND lines referencing them
        print(f"\nStep 3: Processing {input_file} -> {output_file}...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines_read = 0
            lines_written = 0
            lines_removed = 0
            atom_index = 0
            
            for line in infile:
                lines_read += 1
                
                # Fix RESI lines
                if line.startswith('RESI') and 'naa.pdb' in line:
                    line = line.replace('naa.pdb', 'NAA')
                
                # Skip flagged ATOM lines by index
                if line.startswith('ATOM'):
                    atom_index += 1
                    if atom_index in x_atom_indices:
                        lines_removed += 1
                        continue
                
                # Skip BOND lines that reference any removed atom name
                if line.startswith('BOND'):
                    parts = line.split()
                    if any(name in x_atom_names_str for name in parts[1:]):
                        lines_removed += 1
                        continue
                
                outfile.write(line)
                lines_written += 1
        
        print(f"  Lines read: {lines_read}")
        print(f"  Lines removed: {lines_removed}")
        print(f"  Lines written to {output_file}: {lines_written}")
        print(f"\nProcessing complete! Output saved to '{output_file}'")
        
    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    process_str_files()