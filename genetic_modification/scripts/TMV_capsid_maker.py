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

def process_protein(structure, matrix_file, output_prefix):
    """
    Process a protein structure by applying transformations (rotation and translation)
    defined in a matrix file. Uses QB-style naming (e.g., QB_A_com0.pdb).
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

                # Save the transformed PDB file
                output_file = f"{output_prefix}_com{matrix_index}.pdb"
                io.save(output_file)

                # Undo translation
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a - translate
                                atom.set_coord(a)

                # Undo rotation
                transpose(matrix, mat_transpose)
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a.left_multiply(mat_transpose)
                                atom.set_coord(a)

def main():
    parser = PDBParser(PERMISSIVE=True, QUIET=True)

    # Process TMV structure
    input_pdb = "/Users/carydarwin/eviVLP/evi_gen_GUI/TMV_template.pdb"
    matrix_file = "/Users/carydarwin/eviVLP/evi_gen_GUI/matrix_TMV.txt"
    output_dir = "/Users/carydarwin/eviVLP/evi_gen_GUI"
    
    print("="*60)
    print("PROCESSING TMV CAPSID")
    print("="*60)
    print(f"Loading structure from: {input_pdb}")
    print(f"Using transformation matrix: {matrix_file}")
    print(f"Output directory: {output_dir}")
    
    try:
        data_tmv = parser.get_structure("tmv_structure", input_pdb)
        process_tmv_structure(data_tmv, matrix_file, output_dir)
        print("TMV capsid generation completed successfully!")
    except Exception as e:
        print(f"Error processing TMV: {e}")
        import traceback
        traceback.print_exc()

    # Process aligned protein (P2_nu.pdb) with coat protein rotation
    print("\n" + "="*60)
    print("PROCESSING ALIGNED PROTEIN (P2_nu.pdb)")
    print("="*60)
    
    try:
        data_p2 = parser.get_structure("p2_structure", "/Users/carydarwin/eviVLP/evi_gen_GUI/P2_nu.pdb")
        process_protein(data_p2, matrix_file, "TMV_patch")
        print("P2_nu protein processing completed successfully!")
    except Exception as e:
        print(f"Error processing P2_nu: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "="*60)
    print("ALL PROCESSING COMPLETED")
    print("="*60)

if __name__ == "__main__":
    main()