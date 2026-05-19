def number_to_letters(n):
    """Convert a number to a letter combination (e.g., 0 -> AA, 1 -> AB, ..., 25 -> AZ, 26 -> BA, etc.)."""
    result = []
    while n >= 0:
        remainder = n % 26
        result.append(chr(ord('A') + remainder))
        n = n // 26 - 1
        if n < 0:
            break
    return ''.join(reversed(result))

def generate_vmd_script():
    base_dir = "/Users/carydarwin/eviVLP/evi_gen_GUI"
    output_file = f"{base_dir}/TMV_capsid.pgn"
    
    # Build the script content
    script_content = []
    
    # Header
    script_content.extend([
        '# vmd -dispdev text -e TMV_capsid.pgn',
        'package require psfgen',
        'topology top_all27_prot_na_CBD.inp',
        'pdbalias residue HIS HSE',
        'pdbalias atom ILE CD1 CD',
        ''
    ])
    
    # Segments and coordinates for main TMV units
    for i in range(0, 49):
        letters = number_to_letters(i)
        script_content.append(f'segment {letters} {{pdb TMV{i}.pdb; first NONE}}')
        script_content.append(f'coordpdb TMV{i}.pdb {letters}')
    
    # Patches (every 5th unit)
    for i in range(0, 49):
        letters2 = number_to_letters(i+49)
        letters = number_to_letters(i)
        script_content.append(f'segment {letters2} {{pdb TMV_patch_com{i}.pdb; first NONE}}')
        script_content.append(f'coordpdb TMV_patch_com{i}.pdb {letters2}')
        script_content.append(f'patch CLNK {letters}:158 {letters2}:r*s')
    
    script_content.extend(['', ''])  # Blank lines
    
    # Footer
    script_content.extend([
        'guesscoord',
        'regenerate angles dihedrals',
        'guesscoord',
        'writepdb rod.pdb',
        'writepsf rod.psf',
        'package require solvate',
        'solvate rod.psf rod.pdb -t 10 -o TMV_all_water',
        'package require autoionize',
        'autoionize -psf TMV_all_water.psf -pdb TMV_all_water.pdb -sc 0.2 -o TMV_gen_capsid_wb',
        'exit'
    ])
    
    # Write everything at once
    with open(output_file, 'w') as f:
        f.write('\n'.join(script_content))
    
    print(f"Script completed successfully! {output_file} has been generated.")

generate_vmd_script()