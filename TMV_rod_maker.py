from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

def transpose(matrix, mat_transpose):
    """Calculate the transpose of a 3x3 matrix."""
    for i in range(3):
        for j in range(3):
            mat_transpose[i][j] = matrix[j][i]

def process_tmv_structure(pdb_file, matrix_file, output_prefix):
    """
    Process the TMV coat protein structure by applying transformations
    (rotation and translation) defined in a matrix file.
    """
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    data_a = parser.get_structure("test", pdb_file)
    io = PDBIO()
    io.set_structure(data_a)

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
                for model in data_a:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                # Rotate protein
                                a = atom.get_vector()
                                a = a.left_multiply(matrix)
                                atom.set_coord(a)

                for model in data_a:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                # Translate protein
                                a = atom.get_vector()
                                a = a + translate
                                atom.set_coord(a)

                # Save the transformed PDB file
                output_file = f"{output_prefix}{matrix_index}.pdb"
                io.save(output_file)

                # Undo translation
                for model in data_a:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a - translate
                                atom.set_coord(a)

                # Undo rotation
                transpose(matrix, mat_transpose)
                for model in data_a:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a.left_multiply(mat_transpose)
                                atom.set_coord(a)

def main():
    pdb_file = "/Users/carydarwin/eviVLP/tmv_resall.pdb"
    matrix_file = "/Users/carydarwin/eviVLP/matrix_TMV.txt"
    output_prefix = "TMV"
    process_tmv_structure(pdb_file, matrix_file, output_prefix)

if __name__ == "__main__":
    main()
