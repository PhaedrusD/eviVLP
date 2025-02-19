import numpy as np 

def load_pdb_coordinates(pdb_file):
    """Load atomic coordinates and other data from a PDB file."""
    atoms = []
    coordinates = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append(line)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coordinates.append([x, y, z])
    return atoms, np.array(coordinates)

def save_pdb_coordinates(pdb_file, atoms, new_coordinates):
    """Save transformed atomic coordinates back to a PDB file."""
    with open(pdb_file, 'w') as file:
        for atom, coord in zip(atoms, new_coordinates):
            x, y, z = coord
            file.write(f"{atom[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{atom[54:]}")
    print(f"Transformed PDB saved to {pdb_file}")

def load_rotation_matrix(rotation_matrix_file):
    """Load a 4x4 rotation matrix from a text file."""
    return np.loadtxt(rotation_matrix_file)

def transform_pdb(pdb_file, output_file, rotation_matrix_file):
    """Translate to origin, invert the matrix, rotate, and translate back."""
    atoms, coordinates = load_pdb_coordinates(pdb_file)

    # Load the rotation matrix from the file
    rotation_matrix = load_rotation_matrix(rotation_matrix_file)

    # Calculate translation to origin
    translation_vector = -coordinates[0]

    # Translate to origin
    translated_coordinates = coordinates + translation_vector

    # Apply the inverted rotation
    # Extend 3D coordinates to homogeneous coordinates for 4x4 matrix multiplication
    homogeneous_coordinates = np.hstack([translated_coordinates, np.ones((translated_coordinates.shape[0], 1))])
    rotated_homogeneous = np.dot(homogeneous_coordinates, rotation_matrix.T)

    # Convert back to 3D coordinates
    rotated_coordinates = rotated_homogeneous[:, :3]

    # Translate back to the original position
    final_coordinates = rotated_coordinates - translation_vector

    # Save the transformed coordinates
    save_pdb_coordinates(output_file, atoms, final_coordinates)

# Example usage
pdb_file = "P2.pdb"
output_file = "P2_nu.pdb"
rotation_matrix_file = "rotation_matrix.txt"  # Your rotation matrix file

# Apply transformation
transform_pdb(pdb_file, output_file, rotation_matrix_file)


