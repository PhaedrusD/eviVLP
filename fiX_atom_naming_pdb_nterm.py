import re

def extract_atom_info(str_file_path):
    """Extract atom names from str file after updating RESI lines"""
    # First pass: read and update RESI lines
    with open(str_file_path, 'r') as file:
        lines = file.readlines()
    
    updated_lines = []
    for line in lines:
        # Update RESI line if it contains *****
        if line.startswith("RESI") and "*****" in line:
            line = line.replace("*****", "BBS")
        updated_lines.append(line)

        updated_lines = []
    for line in lines:
        # Update RESI line if it contains *****
        if line.startswith("RESI") and "BB_stick.pdb" in line:
            line = line.replace("BB_stick.pdb", "BBS")
        updated_lines.append(line)
    
    # Write updated content back to str file
    with open(str_file_path, 'w') as file:
        file.writelines(updated_lines)
    
    # Second pass: extract atom names from ATOM lines
    atom_info = []
    for line in updated_lines:
        if line.startswith("ATOM"):
            parts = line.split()
            if len(parts) >= 2:
                atom_info.append(parts[1])  # Store the atom name (second column)
    
    return atom_info

def update_pdb_file(pdb_file_path, atom_info):
    """Update pdb file with multiple transformations in sequence"""
    with open(pdb_file_path, 'r') as file:
        lines = file.readlines()
    
    # Step 0: Replace UNNAMED with ADW in COMPND lines
    updated_lines = []
    for line in lines:
        if line.startswith("COMPND") and "UNNAMED" in line:
            line = line.replace("UNNAMED", "BBS")
        updated_lines.append(line)
    
    # Step 1: Replace HETATM with ATOM (both are 6 characters including spaces)
    new_lines = []
    for line in updated_lines:  # Changed from 'lines' to 'updated_lines'
        if line.startswith("HETATM"):
            # Replace "HETATM" with "ATOM  " (ATOM followed by 2 spaces)
            line = "ATOM  " + line[6:]
        new_lines.append(line)
    
    # Step 2: For all lines starting with ATOM, replace 4th column with BBS and 6th column with 1
    updated_lines2 = []
    for line in new_lines:
        if line.startswith("ATOM"):
            # PDB format: columns 18-20 (0-indexed: 17-20) contain the residue name (4th column)
            # PDB format: columns 23-26 (0-indexed: 22-26) contain the residue sequence number (6th column)
            new_line = line[:17] + "BBS" + line[20:22] + "   1" + line[26:]
            updated_lines2.append(new_line)
        else:
            updated_lines2.append(line)
    
    # Step 3: Replace atom names in ATOM lines with names from str file
    final_lines = []
    atom_count = 0
    for line in updated_lines2:
        if line.startswith("ATOM"):
            if atom_count < len(atom_info):
                # Replace atom name (columns 13-16 in PDB format, 0-indexed: 12-16)
                new_line = line[:12] + f"{atom_info[atom_count]:<4}" + line[16:]
                final_lines.append(new_line)
                atom_count += 1
            else:
                final_lines.append(line)
        else:
            final_lines.append(line)
    
    # Write updated content to pdb file
    with open(pdb_file_path, 'w') as file:
        file.writelines(final_lines)
    
    return atom_count

def main():
    str_file_path = 'BB_stick.str'
    pdb_file_path = 'BB_stick.pdb'
    
    # Extract atom info from str file
    atom_info = extract_atom_info(str_file_path)
    atom_count = len(atom_info)
    
    # Update pdb file
    processed_atom_count = update_pdb_file(pdb_file_path, atom_info)
    
    if atom_count == processed_atom_count:
        print("Update successful. Number of ATOM lines from str matches processed ATOM lines in pdb.")
        print(f"Total atoms processed: {atom_count}")
    else:
        print(f"Warning: Mismatch in counts. ATOM lines from str: {atom_count}, Processed ATOM lines: {processed_atom_count}")

if __name__ == "__main__":
    main()
