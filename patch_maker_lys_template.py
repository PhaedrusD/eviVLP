import re

atom_name = ['2*T1', '2*T2', '2*T3', '2*T4', '2*T5', '2*T6']
atom_info = ['2*T1', '2*T2', '2*T3']

def replace_2star(lines):
    processed_lines = []
    pattern = re.compile(r'2t\*(\d)')
    for line in lines:
        def replacer(match):
            num = int(match.group(1))
            if 1 <= num <= len(atom_name):
                return atom_name[num - 1]
            else:
                return match.group(0)
        new_line = pattern.sub(replacer, line)
        processed_lines.append(new_line)
    return processed_lines

def get_str_lines(str_filename):
    atom_lines = []
    for target in ['*T1', '*T2', '*T3']:  # Search for original names without '2'
        found = False
        with open(str_filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0] == 'ATOM' and parts[1] == target:
                    # Modify the line to add '2' prefix to atom name
                    modified_line = line.replace(f'ATOM {target}', f'ATOM 2{target}')
                    atom_lines.append(modified_line)
                    found = True
                    break
            if not found:
                raise ValueError(f"ATOM {target} not found in {str_filename}")
    return atom_lines

def replace_ss(lines, atom_lines):
    processed_lines = []
    group_count = 0
    for line in lines:
        stripped = line.strip()
        if stripped == 'GROUP':
            group_count += 1
            processed_lines.append(line)
        elif stripped == '2i**1' and group_count >= 2:
            processed_lines.append(atom_lines[0])
        elif stripped == '2i**2' and group_count >= 2:
            processed_lines.append(atom_lines[1])
        elif stripped == '2i**3' and group_count >= 2:
            processed_lines.append(atom_lines[2])
        else:
            processed_lines.append(line)
    return processed_lines

def process_files(template_file, str_file, rtf_file):
    # Read and process the template
    with open(template_file, 'r') as f:
        template_lines = f.readlines()
    
    modified_lines = replace_2star(template_lines)
    
    try:
        atom_lines = get_str_lines(str_file)
    except ValueError as e:
        print(f"Error: {e}")
        return
    
    final_lines = replace_ss(modified_lines, atom_lines)
    
    # Filter out lines containing '2i*' before writing to file
    final_lines = [line for line in final_lines if '2t*' not in line]
    
    # Read the RTF file and find the END line
    with open(rtf_file, 'r') as f:
        rtf_content = f.readlines()
    
    end_indices = []
    for idx, line in enumerate(rtf_content):
        if re.match(r'^\s*END\b', line, re.IGNORECASE):
            end_indices.append(idx)
    
    if not end_indices:
        raise ValueError(f"END statement not found in {rtf_file}")
    
    # Use the last occurrence of END
    end_index = end_indices[-1]
    
    # Ensure final_lines ends with a newline if it doesn't already
    if final_lines and not final_lines[-1].endswith('\n'):
        final_lines[-1] = final_lines[-1] + '\n'
    
    # Insert new content before END with an extra newline for separation
    updated_rtf = (
        rtf_content[:end_index] +
        final_lines +
        ['\n'] +  # Add an extra newline before END
        rtf_content[end_index:]
    )
    
    # Write updated content back to RTF file
    with open(rtf_file, 'w') as f:
        f.writelines(updated_rtf)

# Example usage
if __name__ == "__main__":
    process_files('patch_ZAP_lys_template.txt', 'one.str', 'top_all36_cgenff_CBD.rtf')
