import numpy as np
import re

class PDBMolecule:
    def __init__(self, filename):
        self.atoms = []
        self.coordinates = []
        self.atom_names = []
        self.load_pdb(filename)
    
    def load_pdb(self, filename):
        """Load PDB file and extract atom information."""
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Parse PDB format
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    self.atom_names.append(atom_name)
                    self.coordinates.append([x, y, z])
                    self.atoms.append(line)
        
        self.coordinates = np.array(self.coordinates)
    
    def get_atom_coords(self, atom_name):
        """Get coordinates of a specific atom by name."""
        try:
            idx = self.atom_names.index(atom_name)
            return self.coordinates[idx]
        except ValueError:
            raise ValueError(f"Atom '{atom_name}' not found in molecule")
    
    def translate(self, translation_vector):
        """Translate all atoms by the given vector."""
        self.coordinates += translation_vector
    
    def rotate(self, rotation_matrix, center=None):
        """Rotate all atoms around a center point (default: origin)."""
        if center is None:
            center = np.array([0, 0, 0])
        
        # Translate to center, rotate, then translate back
        centered_coords = self.coordinates - center
        rotated_coords = np.dot(centered_coords, rotation_matrix.T)
        self.coordinates = rotated_coords + center
    
    def save_pdb(self, filename):
        """Save the molecule to a PDB file with updated coordinates."""
        with open(filename, 'w') as f:
            for i, atom_line in enumerate(self.atoms):
                # Update coordinates in the PDB line
                new_line = (atom_line[:30] + 
                           f"{self.coordinates[i][0]:8.3f}" +
                           f"{self.coordinates[i][1]:8.3f}" +
                           f"{self.coordinates[i][2]:8.3f}" +
                           atom_line[54:])
                f.write(new_line)

def calculate_rotation_matrix(vec1, vec2):
    """
    Calculate rotation matrix to align vec1 with vec2.
    Uses Rodrigues' rotation formula.
    """
    # Normalize vectors
    v1 = vec1 / np.linalg.norm(vec1)
    v2 = vec2 / np.linalg.norm(vec2)
    
    # If vectors are already aligned, return identity matrix
    if np.allclose(v1, v2):
        return np.eye(3)
    
    # If vectors are opposite, we need a 180-degree rotation
    if np.allclose(v1, -v2):
        # Find a perpendicular vector for 180-degree rotation
        if abs(v1[0]) < 0.9:
            perp = np.array([1, 0, 0])
        else:
            perp = np.array([0, 1, 0])
        
        # Make it truly perpendicular
        perp = perp - np.dot(perp, v1) * v1
        perp = perp / np.linalg.norm(perp)
        
        # 180-degree rotation around perpendicular axis
        return 2 * np.outer(perp, perp) - np.eye(3)
    
    # Calculate rotation axis and angle
    cross_product = np.cross(v1, v2)
    sin_angle = np.linalg.norm(cross_product)
    cos_angle = np.dot(v1, v2)
    
    # Rotation axis (normalized)
    k = cross_product / sin_angle
    
    # Rodrigues' rotation formula
    K = np.array([[0, -k[2], k[1]],
                  [k[2], 0, -k[0]],
                  [-k[1], k[0], 0]])
    
    R = np.eye(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)
    
    return R

def align_molecules(reference_pdb, one_pdb, output_pdb, 
                   ref_atoms=['CD', 'NZ'], one_atoms=['*T1', '*T3'], 
                   target_distance=-1.65):
    """
    Align test molecule to reference molecule based on specified atom pairs.
    
    Parameters:
    - reference_pdb: path to reference PDB file (ref.pdb)
    - one_pdb: path to test PDB file to be aligned
    - output_pdb: path for output aligned PDB file
    - ref_atoms: list of two atom names in reference molecule for vector
    - one_atoms: list of two atom names in test molecule for vector
    - target_distance: desired distance in Angstroms between *REF_ATOM2* and *T1
    """
    
    # Load molecules
    print(f"Loading reference molecule: {reference_pdb}")
    ref_mol = PDBMolecule(reference_pdb)
    
    print(f"Loading test molecule: {one_pdb}")
    one_mol = PDBMolecule(one_pdb)
    
    # Get reference vector (*REF_ATOM1* -> *REF_ATOM2*)
    ref_atom1 = ref_mol.get_atom_coords(ref_atoms[0])
    ref_atom2 = ref_mol.get_atom_coords(ref_atoms[1])
    ref_vector = ref_atom2 - ref_atom1
    
    print(f"Reference vector ({ref_atoms[0]} -> {ref_atoms[1]}): {ref_vector}")
    
    # Get test vector (*T1 -> *T3)
    one_atom1 = one_mol.get_atom_coords(one_atoms[0])
    one_atom2 = one_mol.get_atom_coords(one_atoms[1])
    one_vector = one_atom2 - one_atom1
    
    print(f"Test vector ({one_atoms[0]} -> {one_atoms[1]}): {one_vector}")
    
    # Step 1: Translate test molecule so *T1 is at origin
    print("Translating test molecule...")
    one_mol.translate(-one_atom1)
    
    # Update one_atom2 position after translation
    one_atom2_translated = one_atom2 - one_atom1
    
    # Step 2: Calculate rotation matrix to align vectors
    print("Calculating rotation matrix...")
    rotation_matrix = calculate_rotation_matrix(one_vector, ref_vector)
    
    # Step 3: Rotate test molecule
    print("Rotating test molecule...")
    one_mol.rotate(rotation_matrix)
    
    # Step 4: Calculate precise positioning for target distance
    print(f"Positioning for target distance of {target_distance} Å...")
    
    # Get the unit vector direction from *REF_ATOM1* to *REF_ATOM2*
    ref_unit_vector = ref_vector / np.linalg.norm(ref_vector)
    
    # Position *T1 at the target distance from *REF_ATOM2* along the vector direction
    # *T1 should be at: *REF_ATOM2*_position - target_distance * unit_vector
    target_position = ref_atom2 - target_distance * ref_unit_vector
    
    # Translate test molecule so *T1 is at the target position
    one_mol.translate(target_position)
    
    # Verify alignment and distance
    final_one_atom1 = one_mol.get_atom_coords(one_atoms[0])  # *T1
    final_one_atom2 = one_mol.get_atom_coords(one_atoms[1])  # *T3
    final_one_vector = final_one_atom2 - final_one_atom1
    
    # Calculate actual distance between *REF_ATOM2* and *T1
    actual_distance = np.linalg.norm(ref_atom2 - final_one_atom1)
    
    print(f"\nAlignment and distance verification:")
    print(f"Final test vector: {final_one_vector}")
    print(f"Reference vector:  {ref_vector}")
    print(f"Vector alignment error: {np.linalg.norm(final_one_vector - ref_vector):.6f}")
    print(f"Target distance (NZ-*T1): {target_distance:.3f} Å")
    print(f"Actual distance (NZ-*T1): {actual_distance:.6f} Å")
    print(f"Distance error: {abs(actual_distance - target_distance):.6f} Å")
    
    # Save aligned molecule
    print(f"Saving aligned molecule to: {output_pdb}")
    one_mol.save_pdb(output_pdb)
    
    return one_mol

# Example usage
if __name__ == "__main__":
    try:
        # Align one.pdb to lys.pdb with target distance
        aligned_mol = align_molecules(
            reference_pdb="lys.pdb",
            one_pdb="one.pdb", 
            output_pdb="one_aligned.pdb",
            ref_atoms=['CD', 'NZ'],
            one_atoms=['*T1', '*T3'],
            target_distance=-1.65  # Angstroms
        )
        
        print("\nAlignment completed successfully!")
        
    except FileNotFoundError as e:
        print(f"Error: Could not find PDB file - {e}")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
