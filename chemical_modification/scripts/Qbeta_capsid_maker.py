import os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

# All paths are relative to the script's own directory (eviVLP/GUI/).
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def transpose(matrix, mat_transpose):
    """Calculate the transpose of a 3x3 matrix."""
    for i in range(3):
        for j in range(3):
            mat_transpose[i][j] = matrix[j][i]

def process_protein(structure, matrix_file, output_prefix):
    """
    Process a protein structure by applying transformations (rotation and translation)
    defined in a matrix file.
    """
    io = PDBIO()
    io.set_structure(structure)

    matrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    mat_transpose = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    translate = [0, 0, 0]

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

            matrix[n][0] = float(string_list[1])
            matrix[n][1] = float(string_list[2])
            matrix[n][2] = float(string_list[3])
            translate[n] = float(string_list[4])

            if index == 0:
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a.left_multiply(matrix)
                                atom.set_coord(a)

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a + translate
                                atom.set_coord(a)

                output_file = f"{output_prefix}_com{matrix_index}.pdb"
                io.save(output_file)

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                a = atom.get_vector()
                                a = a - translate
                                atom.set_coord(a)

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
    matrix_file = os.path.join(SCRIPT_DIR, "matrix.txt")

    data_a = parser.get_structure("test_a", os.path.join(SCRIPT_DIR, "fix_a.pdb"))
    process_protein(data_a, matrix_file, "QB_A")

    data_b = parser.get_structure("test_b", os.path.join(SCRIPT_DIR, "fix_b.pdb"))
    process_protein(data_b, matrix_file, "QB_B")

    data_c = parser.get_structure("test_c", os.path.join(SCRIPT_DIR, "fix_c.pdb"))
    process_protein(data_c, matrix_file, "QB_C")

if __name__ == "__main__":
    main()
