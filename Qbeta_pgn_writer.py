# Define the initial lines to write to the file
lines = [
    "# vmd -dispdev text -e QB_capsid.pgn",
    "package require psfgen",
    "topology top_all36_prot.rtf",
    "topology top_all36_cgenff_CBD.rtf",
    "pdbalias residue HIS HSE",
    "pdbalias atom ILE CD1 CD"
]

# Write the initial lines to the file
file_path = '/Users/carydarwin/eviVLP/GUI/QB_capsid.pgn'
with open(file_path, 'w') as f:
    for line in lines:
        f.write(line + '\n')
    f.write('\n')

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

# Generate segments and coordinates for the fit_puzzle
fit_puzzle = []
for i in range(0, 60):
    letters = number_to_letters(i)
    seg = f'segment {letters} {{pdb QB_A_com{i}.pdb; first NONE}}'
    coord = f'coordpdb QB_A_com{i}.pdb {letters}'
    fit_puzzle.append(seg)
    fit_puzzle.append(coord)

for i in range(0, 60):
    letters2 = number_to_letters(i+60)
    seg = f'segment {letters2} {{pdb QB_B_com{i}.pdb; first NONE}}'
    coord = f'coordpdb QB_B_com{i}.pdb {letters2}'
    fit_puzzle.append(seg)
    fit_puzzle.append(coord)

for i in range(0, 60):
    letters3 = number_to_letters(i+120)
    seg = f'segment {letters3} {{pdb QB_C_com{i}.pdb; first NONE}}'
    coord = f'coordpdb QB_C_com{i}.pdb {letters3}'
    fit_puzzle.append(seg)
    fit_puzzle.append(coord)

# Write the fit_puzzle data to the file
with open(file_path, 'a') as f:
    f.write('\n'.join(fit_puzzle) + '\n\n')

# Define the final lines to append to the file
lines = [
    "guesscoord",
    "regenerate angles dihedrals",
    "guesscoord",
    "writepdb capsid_naa.pdb",
    "writepsf capsid_naa.psf",
    "exit"
]

# Append the final lines to the file
with open(file_path, 'a') as f:
    for line in lines:
        f.write(line + '\n')


#    "package require solvate",
#    "solvate capsid_naa.psf capsid_naa.pdb -t 10 -o QB_capsid_water",
#    "package require autoionize",
#    "autoionize -psf QB_capsid_water.psf -pdb QB_capsid_water.pdb -sc 0.2 -o QB_capsid_ion",
