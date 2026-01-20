def filter_str_file(input_file='naa.pdb', output_file='nad.pdb'):
    """
    Read a .pdb file and remove lines that contain the letter 'X'.
    Focuses on lines starting with ATOM or BOND but processes all lines.
    
    Args:
        input_file (str): Input filename (default: 'naa.pdb')
        output_file (str): Output filename (default: 'nad.pdb')
    """
    try:
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
        print(f"Processing complete!")
        print(f"Lines read: {lines_read}")
        print(f"Lines with 'X' removed: {lines_with_x}")
        print(f"Lines written to {output_file}: {lines_written}")
        
    except FileNotFoundError:
        print(f"Error: Could not find input file '{input_file}'")
        print("Make sure the file exists in the current directory.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Alternative version if you want case-insensitive matching
def filter_str_file_case_insensitive(input_file='naa.pdb', output_file='nad.pdb'):
    """
    Same as above but removes lines containing 'X' or 'x' (case-insensitive)
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines_read = 0
            lines_written = 0
            lines_with_x = 0
            
            for line in infile:
                lines_read += 1
                
                # Check if line contains 'X' or 'x' (case-insensitive)
                if 'X' in line.upper():
                    lines_with_x += 1
                    continue  # Skip this line
                
                # Write the line to output file
                outfile.write(line)
                lines_written += 1
        
        print(f"Processing complete!")
        print(f"Lines read: {lines_read}")
        print(f"Lines with 'X' or 'x' removed: {lines_with_x}")
        print(f"Lines written to {output_file}: {lines_written}")
        
    except FileNotFoundError:
        print(f"Error: Could not find input file '{input_file}'")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Run the filter function
    filter_str_file()
    
    # Uncomment the line below if you want case-insensitive filtering instead
    # filter_str_file_case_insensitive()
