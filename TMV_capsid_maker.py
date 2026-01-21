from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

def transpose(matrix, mat_transpose):
    """Calculate the transpose of a 3x3 matrix."""
    for i in range(3):
        for j in range(3):
            mat_transpose[i][j] = matrix[j][i]

def process_tmv_structure(structure, matrix_file, output_prefix):
    """
    Process a TMV structure by applying transformations (rotation and translation)
    defined in a matrix file. Generates TMV0.pdb, TMV1.pdb, etc.
    """
    io = PDBIO()
    io.set_structure(structure)

    matrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # Rotation matrix
    mat_transpose = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # Transpose matrix
    translate = [0, 0, 0]  # Translation vector

    matrix_index = -1
    counter = 0

    with open(matrix_file, "r") as f:
        for line in f.readlines():
            counter += 1
            string_list = line.split()
            index = counter % 3

            if index == 1:
                matrix_index += 1

            n = index - 1
            if n == -1:
                n = 2

            # Populate rotation matrix and translation vector
            matrix[n][0] = float(string_list[1])
            matrix[n][1] = float(string_list[2])
            matrix[n][2] = float(string_list[3])
            translate[n] = float(string_list[4])

            if index == 0:  # Every 3 lines, perform transformations
                print(f"Processing transformation {matrix_index}...")
                
                # Apply rotation
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a.left_multiply(matrix)
                                atom.set_coord(a)

                # Apply translation
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a + translate
                                atom.set_coord(a)

                # Save the transformed PDB file with TMV naming convention
                output_file = f"{output_prefix}/TMV{matrix_index}.pdb"
                io.save(output_file)
                print(f"Saved: {output_file}")

                # Undo translation to restore original structure
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a - translate
                                atom.set_coord(a)

                # Undo rotation to restore original structure
                transpose(matrix, mat_transpose)
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a.left_multiply(mat_transpose)
                                atom.set_coord(a)

    print(f"Total transformations processed: {matrix_index + 1}")

def main():
    parser = PDBParser(PERMISSIVE=True, QUIET=True)

    # Process TMV structure
    input_pdb = "/Users/carydarwin/eviVLP/GUI/TMV_naa.pdb"
    matrix_file = "/Users/carydarwin/eviVLP/GUI/matrix_TMV.txt"
    output_dir = "/Users/carydarwin/eviVLP/GUI"  # Directory where TMV0.pdb, TMV1.pdb, etc. will be saved
    
    print(f"Loading structure from: {input_pdb}")
    print(f"Using transformation matrix: {matrix_file}")
    print(f"Output directory: {output_dir}")
    
    try:
        data_tmv = parser.get_structure("tmv_structure", input_pdb)
        process_tmv_structure(data_tmv, matrix_file, output_dir)
        print("TMV capsid generation completed successfully!")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
