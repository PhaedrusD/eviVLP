#!/usr/bin/env vmd -e
# vmd -e load_capsid.tcl

# Step 1: Run the PGN script to generate PSF/PDB
puts "Running QB_capsid.pgn to generate capsid.psf and capsid.pdb..."

# Step 2: Load the generated PSF/PDB
mol new TMV_naa_capsid.psf
mol addfile TMV_naa_capsid.pdb

# Step 3: Delete default representation
mol delrep 0 top

# Step 4: Add representation for chains I, J, K (CPK, ColorID 1 = red)
mol representation Lines
mol color ColorID 1
mol selection "segid U1"
mol material Opaque
mol addrep top

mol representation CPK
mol color Name
mol selection "segid U1 and resname NAA"
mol material Opaque
mol addrep top

# Step 5: Add representation for patch segments (CPK, ColorID 1 = blue)
mol representation LINES
mol color name
mol selection "all"
mol material Transparent
mol addrep top

puts "Visualization setup complete!"