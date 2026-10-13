# coor_to_pdb.tcl - write minimised coordinates back over the monomer PDB.
# NAMD writes a binary .coor; the capsid makers read PDB.
mol new {{STEM}}.psf type psf waitfor all
mol addfile {{STEM}}_min.coor type namdbin waitfor all
animate write pdb {{STEM}}.pdb beg 0 end 0
puts "wrote {{STEM}}.pdb from {{STEM}}_min.coor"
exit
