def process_str_files(input_file='naa.str', output_file='nad.str'):
    """
    Combined processing:
    1. Fix RESI lines: replace 'naa.pdb' with 'NAA'
    2. Remove lines containing 'X'
    
    Args:
        input_file (str): Input filename (default: 'naa.str')
        output_file (str): Output filename (default: 'nad.str')
    """
    try:
        # Step 1: Fix RESI lines (modifies input file in-place)
        print("Step 1: Fixing RESI lines...")
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        with open(input_file, 'w') as f:
            resi_fixes = 0
            for line in lines:
                if line.startswith('RESI'):
                    if 'naa.pdb' in line:
                        line = line.replace('naa.pdb', 'NAA')
                        resi_fixes += 1
                f.write(line)
        
        print(f"  Fixed {resi_fixes} RESI line(s)")
        
        # Step 2: Filter out lines with 'X'
        print(f"\nStep 2: Removing lines with 'X'...")
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines_read = 0
            lines_written = 0
            lines_with_x = 0
            
            for line in infile:
                lines_read += 1
                
                # Check if line contains 'X' (case-sensitive)
                if 'X' in line:
                    lines_with_x += 1
                    continue  # Skip this line
                
                # Write the line to output file
                outfile.write(line)
                lines_written += 1
        
        # Print summary
        print(f"  Lines read: {lines_read}")
        print(f"  Lines with 'X' removed: {lines_with_x}")
        print(f"  Lines written to {output_file}: {lines_written}")
        print(f"\nProcessing complete! Output saved to '{output_file}'")
        
    except FileNotFoundError:
        print(f"Error: Could not find input file '{input_file}'")
        print("Make sure the file exists in the current directory.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    process_str_files()
