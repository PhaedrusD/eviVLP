import re

# First function: process_str_file
def process_str_file(input_file, output_file):
    # Step 1: Read input file and count occurrences of N values in BOND lines
    with open(input_file, 'r') as f:
        lines = f.readlines()

    bond_pattern = re.compile(r'(BOND\s+)(N\d+)\s+(N\d+)(\s*!.*)')
    atom_pattern = re.compile(r'(ATOM\s+)(N\d+)(\s+.*)')
    
    n_mapping = {}  # To store the N values that need renaming
    n_counts = {}   # To count occurrences of each N value

    # First pass: Identify common N values in BOND lines
    for line in lines:
        match = bond_pattern.match(line)
        if match:
            _, n1, n2, _ = match.groups()
            n_counts[n1] = n_counts.get(n1, 0) + 1
            n_counts[n2] = n_counts.get(n2, 0) + 1

    # Step 2: Map N values that appear more than once (common N values)
    for n, count in n_counts.items():
        if count > 1:
            n_mapping[n] = "N12"  # Rename common N values to N12

    # Step 3: Apply renaming rules in BOND and ATOM lines
    with open(output_file, 'w') as out:
        for line in lines:
            bond_match = bond_pattern.match(line)
            atom_match = atom_pattern.match(line)

            if bond_match:
                prefix, n1, n2, suffix = bond_match.groups()
                n1 = n_mapping.get(n1, n1)
                n2 = n_mapping.get(n2, n2)
                new_line = f"{prefix}{n1}  {n2}{suffix}\n"
                out.write(new_line)
            elif atom_match:
                atom_prefix, atom_n, atom_suffix = atom_match.groups()
                atom_n = n_mapping.get(atom_n, atom_n)
                new_line = f"{atom_prefix}{atom_n}{atom_suffix}\n"
                out.write(new_line)
            else:
                out.write(line)

# Second function: rename_bonds_and_atoms
def rename_bonds_and_atoms(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    bond_pattern = re.compile(r"^BOND.*(H\w{1,2}).*(N\w{1,2})|^BOND.*(N\w{1,2}).*(H\w{1,2})")
    rename_dict = {}

    # Identify BOND lines and match patterns
    for line in lines:
        bond_match = bond_pattern.search(line)
        if bond_match:
            n_value = bond_match.group(2) if bond_match.group(2) else bond_match.group(3)
            second_bond_pattern = re.compile(rf"^BOND.*(N{n_value[1:]})\b.*N\w{{1,2}}")
            for bond_line in lines:
                if second_bond_pattern.search(bond_line):
                    rename_dict[n_value] = 'N11'
                    break

    # Replace matched N values in all BOND and ATOM lines
    renamed_lines = []
    for line in lines:
        if line.startswith("BOND") or line.startswith("ATOM"):
            for old_n, new_n in rename_dict.items():
                line = re.sub(rf"\b{old_n}\b", new_n, line)
        renamed_lines.append(line)

    # Writing the modified lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(renamed_lines)

# Third function: modify_str_file
def modify_str_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    stored_n_value = None
    bond_pattern = re.compile(r'^BOND\s+.*(N12)\s+(N\w{1,2})\s|^(BOND\s+.*(N\w{1,2})\s+(N12)\s)')
    bond_atom_pattern = re.compile(r'^(BOND|ATOM)\s+(.*)')
    modified_lines = []
    
    # First pass: find the N value paired with N12 and store it, unless it's N11
    for line in lines:
        if line.startswith('BOND'):
            match = bond_pattern.search(line)
            if match:
                n12_first_pair = match.group(1) == 'N12' and match.group(2)
                n12_second_pair = match.group(4) == 'N12' and match.group(3)
                
                if n12_first_pair and match.group(2) != 'N11':
                    stored_n_value = match.group(2)
                    break
                elif n12_second_pair and match.group(3) != 'N11':
                    stored_n_value = match.group(3)
                    break

    if stored_n_value:
        # Second pass: Replace stored N value (with a space) in both BOND and ATOM lines
        for line in lines:
            match = bond_atom_pattern.search(line)
            if match and f'{stored_n_value} ' in line:
                line = line.replace(f'{stored_n_value} ', 'N13 ')
            modified_lines.append(line)
    else:
        modified_lines = lines

    # Save the modified content to a new file
    with open(file_path, 'w') as file:
        file.writelines(modified_lines)

# Main execution function
def process_files_in_sequence(input_file):
    output_file = input_file  # Overwrite same file after each step
    
    process_str_file(input_file, output_file)  # First function
    rename_bonds_and_atoms(output_file)  # Second function
    modify_str_file(output_file)  # Third function

# Example usage
process_files_in_sequence('input1.str')
