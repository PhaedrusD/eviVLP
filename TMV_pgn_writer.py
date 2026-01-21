def generate_vmd_script():
    base_dir = "/Users/carydarwin/eviVLP/GUI"
    output_file = f"{base_dir}/TMV_all_out.pgn"
    
    # Build the script content
    script_content = []
    
    # Header
    script_content.extend([
        '# vmd -dispdev text -e TMV_all_out.pgn',
        'package require psfgen',
        'topology top_all36_prot.rtf',
        'topology top_all36_cgenff_CBD.rtf',
        'pdbalias residue HIS HSE',
        'pdbalias atom ILE CD1 CD',
        ''
    ])
    
    # Segments and coordinates
    for i in range(0, 49):
        script_content.append(f'segment U{i} {{pdb TMV{i}.pdb}}')
        script_content.append(f'coordpdb TMV{i}.pdb U{i}')
    
    script_content.extend(['', ''])  # Blank lines
    
    # Footer
    script_content.extend([
        'guesscoord',
        'regenerate angles dihedrals',
        'guesscoord',
        'writepdb TMV_naa_capsid.pdb',
        'writepsf TMV_naa_capsid.psf',
        'package require solvate',
        'solvate TMV_naa_capsid.psf TMV_naa_capsid.pdb -t 10 -o TMV_all_water',
        'package require autoionize',
        'autoionize -psf TMV_all_water.psf -pdb TMV_all_water.pdb -sc 0.2 -o TMV_naa_capsid_wb',
        'exit'
    ])
    
    # Write everything at once
    with open(output_file, 'w') as f:
        f.write('\n'.join(script_content))
    
    print(f"Script completed successfully! {output_file} has been generated.")

generate_vmd_script()
