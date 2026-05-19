#!/usr/bin/env vmd -e
# vmd -e load_capsid.tcl

# Step 1: Run the PGN script to generate PSF/PDB
puts "Running QB_capsid.pgn to generate capsid.psf and capsid.pdb..."

# Step 2: Load the generated PSF/PDB
mol new capsid.psf
mol addfile capsid.pdb

# Step 3: Delete default representation
mol delrep 0 top

# Step 4: Add representation for chains I, J, K (CPK, ColorID 7 = red)
mol representation NewCartoon
mol color ColorID 7
mol selection "segname IG "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname IL "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname IQ "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname IV "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname JA "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname JF "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname JK "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname JP "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname JU "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname JZ "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname KE "
mol material Glossy
mol addrep top

mol representation NewCartoon
mol color ColorID 7
mol selection "segname KJ "
mol material Glossy
mol addrep top

# Step 5: Add representation for patch segments (CPK, ColorID 7 = blue)
mol representation LINES
mol color name
mol selection "all"
mol material Transparent
mol addrep top

puts "Visualization setup complete!"