#!/usr/bin/env vmd -e
# vmd -e load_rod.tcl

# Step 1: Run the PGN script to generate PSF/PDB
puts "Running TMV_capsid.pgn to generate rod.psf and rod.pdb..."

# Step 2: Load the generated PSF/PDB
mol new rod.psf
mol addfile rod.pdb

# Step 3: Delete default representation
mol delrep 0 top

# Step 4: Add representation for chains I, J, K (CPK, ColorID 1 = red)
mol representation CPK
mol color ColorID 1
mol selection "segname A "
mol material Opaque
mol addrep top

mol representation CPK
mol color ColorID 1
mol selection "segname AX "
mol material Opaque
mol addrep top



# Step 5: Add representation for patch segments (CPK, ColorID 1 = blue)
mol representation LINES
mol color name
mol selection "all"
mol material Transparent
mol addrep top

puts "Visualization setup complete!"