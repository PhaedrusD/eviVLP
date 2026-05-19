#!/usr/bin/env python3
"""
charge_fix_ntrm.py
Reads the charge of the removed HX atom from ntrm.str (pre-noX),
then adds it to atom C's charge in ntrm_clean.str (post-noX).
"""

def get_atom_charge(filename, atom_name):
    """Extract the partial charge for a given atom name from a .str file."""
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                parts = line.split()
                if len(parts) >= 4 and parts[1] == atom_name:
                    return float(parts[3])
    return None

def update_atom_charge(filename, atom_name, new_charge):
    """Update the charge of a specific atom in a .str file."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    updated = False
    new_lines = []
    for line in lines:
        if line.startswith('ATOM'):
            parts = line.split()
            if len(parts) >= 4 and parts[1] == atom_name:
                old_charge = parts[3]
                line = line.replace(old_charge, f"{new_charge:.3f}", 1)
                updated = True
        new_lines.append(line)
    
    if updated:
        with open(filename, 'w') as f:
            f.writelines(new_lines)
    
    return updated

def main():
    source_atom = 'HX'
    target_atom = 'C'
    original_file = 'ntrm.str'
    clean_file = 'ntrm_clean.str'
    
    # Get charge of removed atom from original file
    hx_charge = get_atom_charge(original_file, source_atom)
    if hx_charge is None:
        print(f"Error: Could not find atom {source_atom} in {original_file}")
        return
    print(f"Charge on removed atom {source_atom}: {hx_charge:.3f}")
    
    # Get current charge of target atom in cleaned file
    c_charge = get_atom_charge(clean_file, target_atom)
    if c_charge is None:
        print(f"Error: Could not find atom {target_atom} in {clean_file}")
        return
    print(f"Current charge on {target_atom}: {c_charge:.3f}")
    
    # Add removed charge to target
    new_charge = c_charge + hx_charge
    print(f"New charge on {target_atom}: {new_charge:.3f}")
    
    # Update the cleaned file
    if update_atom_charge(clean_file, target_atom, new_charge):
        print(f"Successfully updated {target_atom} charge in {clean_file}")
    else:
        print(f"Error: Failed to update charge in {clean_file}")

if __name__ == "__main__":
    main()
