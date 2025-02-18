import numpy as np
import MDAnalysis as mda

# Load the first and second protein structures
u1 = mda.Universe("QB_ready.pdb")  # First protein (reference)
u2 = mda.Universe("****.pdb")      # Second protein (target)

# Select the specific atoms for translation
c_atom = u1.select_atoms("resid 132 and name C")  # C-terminus carbon in the first protein
n_atom = u2.select_atoms("resid 1 and name N")    # N-terminus nitrogen in the second protein

# Ensure both selections are valid
if len(c_atom) != 1 or len(n_atom) != 1:
    raise ValueError("Atom selections failed. Check residue IDs or atom names.")

# Get the coordinates of the selected atoms
c_coords = c_atom.positions[0]  # Coordinates of the C atom
n_coords = n_atom.positions[0]  # Coordinates of the N atom

# Calculate the translation vector to align N to C
translation_vector = c_coords - n_coords

# Translate the entire second protein
u2.atoms.translate(translation_vector)

# Translate the second protein 6 Å along the x-axis
translation_6x = np.array([-6.0, 0.0, 0.0])  # 6 Å shift in the x-direction
u2.atoms.translate(translation_6x)

# Write the transformed coordinates of the second protein
u2.atoms.write("P2.pdb")
with mda.Writer("P2.xyz", n_atoms=u2.atoms.n_atoms) as xyz_writer:
    xyz_writer.write(u2)

# Write the coordinates of the first protein
u1.atoms.write("P1.pdb")
with mda.Writer("P1.xyz", n_atoms=u1.atoms.n_atoms) as xyz_writer:
    xyz_writer.write(u1)

print("Files written: P1.pdb, P1.xyz, P2.pdb, P2.xyz")
