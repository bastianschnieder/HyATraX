# S R 
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

print("########################################################")
print("           Transition State Guess Generation            ")
print("        C-H Bonds in Polyaromatic Hydrocarbons          ")
print("########################################################")

# 
# Helper Functions 
#

def detect_aldehyde(mol,cidx,hidx): 
    carbon = mol.GetAtomWithIdx(cidx)
    if carbon.GetSymbol()!='C':
        return False
    vicinals = carbon.GetNeighbors()
    hydrogens = 0
    CO_double_bond = False

    for vicinal in vicinals:
        if vicinal.GetSymbol() == 'H':
            hydrogens += 1
        elif vicinal.GetSymbol() == 'O':
            bond = mol.GetBondBetweenAtoms(cidx,vicinal.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                CO_double_bond = True
    return hydrogens == 1 and CO_double_bond

def get_vector_length(vector):
    length = np.linalg.norm(vector)
    return length

def xyz_to_mol(coord_path):
    raw_mol = AllChem.MolFromXYZFile(coord_path)
    return raw_mol

def replace_h(path, idx, new_line, output):
    with open(path, 'r') as file:
        lines = file.readlines()
        lines[idx] = f'{new_line}\n'
    with open(output, 'w') as file:
        num_atoms = int(lines[0].strip())
        num_atoms += 2
        lines[0] = f'{num_atoms}\n' 
        file.writelines(lines)

def delete_h(path, idx, output):
    with open(path, 'r') as file:
        lines = file.readlines()
        del lines[idx]
    with open(output,'w') as file:
        num_atoms = int(lines[0].strip())
        num_atoms -= 1
        lines[0] = f'{num_atoms}\n' 
        file.writelines(lines)
        print(f'Product for H idx {idx-1} saved as XYZ File in {output}.')

def teta_v(vector1, vector2):
    vector1_norm = vector1 / np.linalg.norm(vector1)
    vector2_norm = vector2 / np.linalg.norm(vector2)
    angle = np.arccos(np.dot(vector1_norm, vector2_norm))
    return angle

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
        atom.SetAtomMapNum(atom.GetIdx()+1)
    return mol

def scale_ch(vector, ch_ts_length):
    if get_vector_length(vector) == 0:
        raise ValueError("Vector has a length l = 0, cannot be scaled!")
    scale = float(ch_ts_length) / float(get_vector_length(vector))
    new_vector = scale * vector
    return new_vector

def add_hydrogen(vector, h_distance):
    if get_vector_length(vector) == 0:
        raise ValueError("Vector has a length l = 0, cannot be scaled!")
    scale = float(h_distance) / float(get_vector_length(vector))
    new_vector = scale * vector
    return new_vector

def rotate_vector(vector, angle, axis):
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    rotated_vector = (np.cos(angle) * vector + np.sin(angle) * np.cross(axis, vector) + (1-np.cos(angle)) * np.dot(axis,vector) * axis)
    return rotated_vector

def add_hydrogen_to_oxygen(vector, h_distance):
    if get_vector_length(vector) == 0:
        raise ValueError("Vector has a length l = 0, cannot be scaled!")
    scale = float(h_distance) / float(get_vector_length(vector))
    new_vector = scale * vector
    return new_vector

def get_atomic_numbers(mol):
    for atom in mol.GetAtoms():
        Num = atom.GetAtomicNum()
        AtmNum.append(Num)

#
# Start of the Program 
#

if len(sys.argv) <= 1:
    print ("usage: python3 main.py <xyz> | e.g. python3 main.py coord.xyz ")
    sys.exit("Molecular coordinates are missing")
else:
    file_path = sys.argv[1] # Molecular Coordinates in xyz format (Angstrom)  

# Generates Products Directory for saving xyz files of products 

results= 'Products'
try: 
     os.makedirs(results)
except FileExistsError:
    pass
print("Products will be saved in Directory: Products")

#
# Define Lists 
#

# Aromatic C-H Bonds

H_vectors              = [] # Contains Vectors of hydrogen atoms (existing in the initial molecule)
C_vectors              = [] # Contains Vectors of carbon atoms (existing in the initial molecule)
Hydrogen_Positions     = []
new_coords             = [] # Contains Vectors (c to h) of the elongated C-H bond
add_hydrogen_positions = []
C_H_Vector             = []
ch_bonds               = [] # Contains C-H bonds detected by rdkit
new_H                  = [] # Contains Coordinates(Angstrom) of the existing hydrogen atom but elongated to the ts length
add_H                  = [] # Contains Coordinates(Angstrom) of the newly added oxygen atom ||||||| RENAME
third_H                = []
add_HtoO_r             = [] # Contains Coordinates(Angstrom) of the newly added hydrogen atom of the OH bond  
AtmNum                 = [] # Contains atomic numbers of the molecule 

# Aldehyde C-H Bonds 

aldehyde_bonds                  = []       
H_vectors_aldehyde              = []
C_vectors_aldehyde              = []
Hydrogen_Position_aldehyde      = []
new_coords_aldehyde             = []
add_hydrogen_positions_aldehyde = []
new_H_aldehyde                  = []
add_H_aldehyde                  = []
third_H_aldehyde                = []
add_HtoO_r_aldehyde             = []

# Aldehyde O-H Bonds 

oh_bonds                  = []       
H_vectors_oh              = []
C_vectors_oh              = []
Hydrogen_Position_oh      = []
new_coords_oh             = []
add_hydrogen_positions_oh = []
new_H_oh                  = []
add_H_oh                  = []
third_H_oh                = []
add_HtoO_r_oh             = []

# Aldehyde N-H Bonds 

nh_bonds                  = []       
H_vectors_nh              = []
C_vectors_nh              = []
Hydrogen_Position_nh      = []
new_coords_nh             = []
add_hydrogen_positions_nh = []
new_H_nh                  = []
add_H_nh                  = []
third_H_nh                = []
add_HtoO_r_nh             = []

# 
# Generate Mol Object
#

raw_mol = AllChem.MolFromXYZFile(file_path)
mol= AllChem.Mol(raw_mol)
rdDetermineBonds.DetermineBonds(mol,charge=0) # set charge correctly!
get_atomic_numbers(mol)

#
# TS Bond Distances
#

# if else statement allows to modify TS Bond Distances for different atom compositions
# e.g. if molecule contains nitrogen, adjust values accordingly 

if 7 in AtmNum:
    print("Using Nitrogen-Updated Parameters for Hydrocarbons")
    ts_value =  1.52                # Elongated C-H Distance (Angstrom)  
    ts_value_aldehyde = 1.23        # Elongated C-H Distance for Aldehyde Abstraction  
    h_distance = ts_value+0.85      # Additional H Atom Distance to C Atom (Angstrom)    
    h_distance_aldehyde = ts_value_aldehyde + 1.16
    ts_value_oh = 1.230
    ts_value_nh = 0.310
    h_distance_oh = ts_value_oh + 0.85
    h_distance_nh = ts_value_nh + 0.91
    oh_angle = 81
    oh_distance = 0.98

else:
    print("Using Default Parameters for Hydrocarbons")
    ts_value =  1.52 
    ts_value_aldehyde = 1.23 
    h_distance = ts_value+0.85 
    h_distance_aldehyde = ts_value_aldehyde + 1.16
    ts_value_oh = 2 
    ts_value_nh = 2
    h_distance_oh = 3
    h_distance_nh = 3
    oh_angle = 81
    oh_distance = 0.98


#################################################################################################
# 
# Main Routine 
#
#################################################################################################

for bond in mol.GetBonds():
    atom1_idx = bond.GetBeginAtomIdx()
    atom2_idx = bond.GetEndAtomIdx()
    # Check for CH Bonds in different functional groups 
    if mol.GetAtomWithIdx(atom1_idx).GetSymbol() == 'C' and mol.GetAtomWithIdx(atom2_idx).GetSymbol() == 'H':
        if detect_aldehyde(mol, atom1_idx,atom2_idx):
            aldehyde_bonds.append((atom1_idx, atom2_idx)) 
        else:
            ch_bonds.append((atom1_idx, atom2_idx)) 
    elif mol.GetAtomWithIdx(atom1_idx).GetSymbol() == 'H' and mol.GetAtomWithIdx(atom2_idx).GetSymbol() == 'C':
        if detect_aldehyde(mol,atom2_idx,atom1_idx):
            aldehyde_bonds.append((atom2_idx, atom1_idx)) 
        else:
            ch_bonds.append((atom2_idx, atom1_idx)) 
    # Check for N-H and O-H Bonds
    if mol.GetAtomWithIdx(atom1_idx).GetSymbol() == 'N' and mol.GetAtomWithIdx(atom2_idx).GetSymbol() == 'H':
        nh_bonds.append((atom1_idx,atom2_idx)) 
    elif mol.GetAtomWithIdx(atom1_idx).GetSymbol() == 'H' and mol.GetAtomWithIdx(atom2_idx).GetSymbol() == 'N':
        nh_bonds.append((atom1_idx,atom2_idx)) 
    if mol.GetAtomWithIdx(atom1_idx).GetSymbol() == 'O' and mol.GetAtomWithIdx(atom2_idx).GetSymbol() == 'H':
        oh_bonds.append((atom1_idx, atom2_idx))
    elif mol.GetAtomWithIdx(atom1_idx).GetSymbol() == 'H' and mol.GetAtomWithIdx(atom2_idx).GetSymbol() == 'O':
        oh_bonds.append((atom1_idx, atom2_idx)) 

configuration = mol.GetConformer()

#
# Generating X-H Bond Vectors
#

# Dealing with Aromatic C-H Bonds 

for index in ch_bonds:
    atom_position_C = configuration.GetAtomPosition(index[0])
    atom_position_H = configuration.GetAtomPosition(index[1])
    vector_Hydrogen = atom_position_H.x, atom_position_H.y, atom_position_H.z
    vector_Carbon = atom_position_C.x, atom_position_C.y, atom_position_C.z
    H_vectors.append(vector_Hydrogen)
    C_vectors.append(vector_Carbon)

# Dealing with Aldehyde C-H Bonds 

for index in aldehyde_bonds:
    atom_position_C_aldehyde = configuration.GetAtomPosition(index[0])
    atom_position_H_aldehyde = configuration.GetAtomPosition(index[1])
    vector_Hydrogen = atom_position_H_aldehyde.x, atom_position_H_aldehyde.y, atom_position_H_aldehyde.z
    vector_Carbon = atom_position_C_aldehyde.x, atom_position_C_aldehyde.y, atom_position_C_aldehyde.z
    H_vectors_aldehyde.append(vector_Hydrogen)
    C_vectors_aldehyde.append(vector_Carbon)

# Dealing with N-H Bonds

for index in nh_bonds:
    atom_position_C_nh = configuration.GetAtomPosition(index[0])
    atom_position_H_nh = configuration.GetAtomPosition(index[1])
    vector_Hydrogen = atom_position_H_nh.x, atom_position_H_nh.y, atom_position_H_nh.z
    vector_Carbon = atom_position_C_nh.x, atom_position_C_nh.y, atom_position_C_nh.z
    H_vectors_nh.append(vector_Hydrogen)
    C_vectors_nh.append(vector_Carbon)

# Dealing with O-H Bonds
    
for index in oh_bonds:
    atom_position_C_oh = configuration.GetAtomPosition(index[0])
    atom_position_H_oh = configuration.GetAtomPosition(index[1])
    vector_Hydrogen = atom_position_H_oh.x, atom_position_H_oh.y, atom_position_H_oh.z
    vector_Carbon = atom_position_C_oh.x, atom_position_C_oh.y, atom_position_C_oh.z
    H_vectors_oh.append(vector_Hydrogen)
    C_vectors_oh.append(vector_Carbon)

# ...
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') 

def main():
  
    radians = np.radians(oh_angle)
    axis = [-1,0,0] 

    # Aromatic C-H Bonds

    for c, h in zip(C_vectors, H_vectors):  # handle c,h as tuples
        v_c = np.array(c)
        v_h = np.array(h)
        C_H_Vector = v_h - v_c
        Hydrogen_Positions.append(C_H_Vector)
        ax.quiver(c[0], c[1], c[2], C_H_Vector[0], C_H_Vector[1], C_H_Vector[2], arrow_length_ratio=0.8, color='r')
        new_coord = scale_ch(C_H_Vector, ts_value)
        new_coords.append(new_coord)
        ax.quiver(c[0], c[1], c[2], new_coord[0], new_coord[1], new_coord[2], arrow_length_ratio=0.55, color='g')
        add_hydrogen_position = add_hydrogen(C_H_Vector, h_distance)
        add_hydrogen_positions.append(add_hydrogen_position)
        ax.quiver(c[0], c[1], c[2], add_hydrogen_position[0], add_hydrogen_position[1], add_hydrogen_position[2], arrow_length_ratio=0.1, color='b')
        new_H.append((c[0] + new_coord[0], c[1] + new_coord[1], c[2] + new_coord[2]))
        add_H.append((c[0] + add_hydrogen_position[0], c[1] + add_hydrogen_position[1], c[2] + add_hydrogen_position[2]))
        
        add_HtoO = add_hydrogen_to_oxygen(C_H_Vector, oh_distance)
        third_H.append(add_HtoO)
        rotated_vector = rotate_vector(add_HtoO, radians, axis)
        add_HtoO_r.append((c[0] + add_hydrogen_position[0] + rotated_vector[0],c[1] + add_hydrogen_position[1] + rotated_vector[1],c[2] + add_hydrogen_position[2] + rotated_vector[2]))
        ax.quiver(c[0] + add_hydrogen_position[0], c[1] + add_hydrogen_position[1], c[2] + add_hydrogen_position[2], rotated_vector[0], rotated_vector[1], rotated_vector[2], arrow_length_ratio=2, color='c')
    
    # Aldehyde C-H Bonds

    for c, h in zip(C_vectors_aldehyde, H_vectors_aldehyde): 
        v_c_aldehyde = np.array(c)
        v_h_aldehyde = np.array(h)
        C_H_Vector_aldehyde = v_h_aldehyde - v_c_aldehyde
        Hydrogen_Position_aldehyde.append(C_H_Vector_aldehyde)
        ax.quiver(c[0], c[1], c[2], C_H_Vector_aldehyde[0], C_H_Vector_aldehyde[1], C_H_Vector_aldehyde[2], arrow_length_ratio=0.8, color='r')
        new_coord_aldehyde = scale_ch(C_H_Vector_aldehyde,ts_value_aldehyde)
        new_coords_aldehyde.append(new_coord_aldehyde)
        ax.quiver(c[0], c[1], c[2], new_coord_aldehyde[0], new_coord_aldehyde[1], new_coord_aldehyde[2], arrow_length_ratio=0.55, color='g')
        add_hydrogen_position_aldehyde = add_hydrogen(C_H_Vector_aldehyde, h_distance_aldehyde)
        add_hydrogen_positions_aldehyde.append(add_hydrogen_position_aldehyde)
        ax.quiver(c[0], c[1], c[2], add_hydrogen_position_aldehyde[0], add_hydrogen_position_aldehyde[1], add_hydrogen_position_aldehyde[2], arrow_length_ratio=0.1, color='b')
        new_H_aldehyde.append((c[0] + new_coord_aldehyde[0], c[1] + new_coord_aldehyde[1], c[2] + new_coord_aldehyde[2]))
        add_H_aldehyde.append((c[0] + add_hydrogen_position_aldehyde[0], c[1] + add_hydrogen_position_aldehyde[1], c[2] + add_hydrogen_position_aldehyde[2]))
        
        add_HtoO = add_hydrogen_to_oxygen(C_H_Vector_aldehyde, oh_distance)
        third_H_aldehyde.append(add_HtoO)
        rotated_vector_aldehyde = rotate_vector(add_HtoO, radians, axis)
        add_HtoO_r_aldehyde.append((c[0] + add_hydrogen_position_aldehyde[0] + rotated_vector_aldehyde[0],c[1] + add_hydrogen_position_aldehyde[1] + rotated_vector_aldehyde[1],c[2] + add_hydrogen_position_aldehyde[2] + rotated_vector_aldehyde[2]))
        ax.quiver(c[0] + add_hydrogen_position_aldehyde[0], c[1] + add_hydrogen_position_aldehyde[1], c[2] + add_hydrogen_position_aldehyde[2], rotated_vector_aldehyde[0], rotated_vector_aldehyde[1], rotated_vector_aldehyde[2], arrow_length_ratio=2, color='c')

    for n, h in zip(C_vectors_nh, H_vectors_nh):
        v_n = np.array(n)
        v_h = np.array(h)
        N_H_Vector = v_n - v_h
        Hydrogen_Position_nh.append(N_H_Vector)
        ax.quiver(n[0],n[1],n[2], N_H_Vector[0],N_H_Vector[1],N_H_Vector[2],arrow_length_ratio=0.8, color = 'r')
        new_coord_nh = scale_ch(N_H_Vector, ts_value_nh)
        new_coords_nh.append(new_coord_nh)
        ax.quiver(n[0],n[1],n[2],new_coord_nh[0],new_coord_nh[1],new_coord_nh[2], arrow_length_ratio=0.55, color = 'g')
        add_hydrogen_nh = add_hydrogen(N_H_Vector,h_distance_nh)
        add_hydrogen_positions_nh.append(add_hydrogen_nh)
        ax.quiver(n[0],n[1],n[2],add_hydrogen_nh[0],add_hydrogen_nh[1],add_hydrogen_nh[2], arrow_length_ratio = 0.1, color = 'b')
        new_H_nh.append((n[0] + new_coord_nh[0], n[1] + new_coord_nh[1], n[2] + new_coord_nh[2]))
        add_H_nh.append((n[0] + add_hydrogen_nh[0], n[1] + add_hydrogen_nh[1], n[2] + add_hydrogen_nh[2]))
        
        add_HtoO = add_hydrogen_to_oxygen(N_H_Vector, oh_distance)
        third_H_nh.append(add_HtoO)
        rotated_vector = rotate_vector(add_HtoO, radians, axis)
        add_HtoO_r_nh.append((n[0] + add_hydrogen_nh[0] + rotated_vector[0],n[1] + add_hydrogen_nh[1] + rotated_vector[1],n[2] + add_hydrogen_nh[2] + rotated_vector[2]))
        ax.quiver(n[0] + add_hydrogen_nh[0], n[1] + add_hydrogen_nh[1], n[2] + add_hydrogen_nh[2], rotated_vector[0], rotated_vector[1], rotated_vector[2], arrow_length_ratio=2, color='c')
   
    for o, h in zip(C_vectors_oh, H_vectors_oh):
        v_o = np.array(o)
        v_h = np.array(h)
        O_H_Vector = v_o - v_h
        Hydrogen_Position_oh.append(O_H_Vector)
        ax.quiver(o[0],o[1],o[2], O_H_Vector[0],O_H_Vector[1],O_H_Vector[2],arrow_length_ratio=0.8, color = 'r')
        new_coord_oh = scale_ch(O_H_Vector, ts_value_oh)
        new_coords_oh.append(new_coord_oh)
        ax.quiver(o[0],o[1],o[2],new_coord_oh[0],new_coord_oh[1],new_coord_oh[2], arrow_length_ratio=0.55, color = 'g')
        add_hydrogen_oh = add_hydrogen(O_H_Vector,h_distance_oh)
        add_hydrogen_positions_oh.append(add_hydrogen_oh)
        ax.quiver(o[0],o[1],o[2],add_hydrogen_oh[0],add_hydrogen_oh[1],add_hydrogen_oh[2], arrow_length_ratio = 0.1, color = 'b')
        new_H_oh.append((o[0] + new_coord_oh[0], o[1] + new_coord_oh[1], o[2] + new_coord_oh[2]))
        add_H_oh.append((o[0] + add_hydrogen_oh[0], o[1] + add_hydrogen_oh[1], o[2] + add_hydrogen_oh[2]))

        add_HtoO = add_hydrogen_to_oxygen(O_H_Vector, oh_distance)
        third_H_oh.append(add_HtoO)
        rotated_vector = rotate_vector(add_HtoO, radians, axis)
        add_HtoO_r_oh.append((o[0] + add_hydrogen_oh[0] + rotated_vector[0],o[1] + add_hydrogen_oh[1] + rotated_vector[1],o[2] + add_hydrogen_oh[2] + rotated_vector[2]))
        ax.quiver(o[0] + add_hydrogen_oh[0], o[1] + add_hydrogen_oh[1], o[2] + add_hydrogen_oh[2], rotated_vector[0], rotated_vector[1], rotated_vector[2], arrow_length_ratio=2, color='c')

if __name__ == "__main__":
    main()

#################################################################################################
# 
# Generating Coordinate Files
#
#################################################################################################

#
# Saving Changed XYZ Files
#

for ch_bond,hydrogen,added_hydrogen,oh_h in zip(ch_bonds, new_H, add_H,add_HtoO_r):
    idx = ch_bond[1]+2 
    new_line = f"H {hydrogen[0]} {hydrogen[1]} {hydrogen[2]}"
    second_h = f"O {added_hydrogen[0]} {added_hydrogen[1]} {added_hydrogen[2]}"
    third_h  = f"H {oh_h[0]} {oh_h[1]} {oh_h[2]}"
    path = file_path
    output = f'coord_hidx_{idx-1}.xyz' 
    replace_h(path, idx, new_line, output)
    with open(f'coord_hidx_{idx-1}.xyz','a') as file:
        file.write(second_h + '\n')
        file.write(third_h + '\n') 
        print(f'Finished with OH Radical Positioning on Former H with Idx {idx-1}. Saved as XYZ File in {output}.')

for ch_bond,hydrogen,added_hydrogen,oh_h in zip(aldehyde_bonds, new_H_aldehyde, add_H_aldehyde,add_HtoO_r_aldehyde):
    idx = ch_bond[1]+2 
    new_line = f"H {hydrogen[0]} {hydrogen[1]} {hydrogen[2]}"
    second_h = f"O {added_hydrogen[0]} {added_hydrogen[1]} {added_hydrogen[2]}"
    third_h  = f"H {oh_h[0]} {oh_h[1]} {oh_h[2]}"
    path = file_path
    output = f'coord_hidx_{idx-1}.xyz' 
    replace_h(path, idx, new_line, output)
    with open(f'coord_hidx_{idx-1}.xyz','a') as file:
        file.write(second_h + '\n')
        file.write(third_h + '\n') 
        print(f'Finished with OH Radical Positioning on Former H with Idx {idx-1}. Saved as XYZ File in {output}.')

for ch_bond,hydrogen,added_hydrogen,oh_h in zip(nh_bonds, new_H_nh, add_H_nh, add_HtoO_r_nh): # BUG add_HtoO_r_nh!
    idx = ch_bond[0]+2 
    new_line = f"H {hydrogen[0]} {hydrogen[1]} {hydrogen[2]}"
    second_h = f"O {added_hydrogen[0]} {added_hydrogen[1]} {added_hydrogen[2]}"
    third_h  = f"H {oh_h[0]} {oh_h[1]} {oh_h[2]}"
    path = file_path
    output = f'coord_hidx_{idx-1}.xyz' 
    replace_h(path, idx, new_line, output)
    with open(f'coord_hidx_{idx-1}.xyz','a') as file:
        file.write(second_h + '\n')
        file.write(third_h + '\n') 
        print(f'Finished with OH Radical Positioning on Former H with Idx {idx-1}. Saved as XYZ File in {output}.')

for ch_bond,hydrogen,added_hydrogen in zip(oh_bonds, new_H_oh, add_H_oh):
    idx = ch_bond[0]+2 
    new_line = f"H {hydrogen[0]} {hydrogen[1]} {hydrogen[2]}"
    second_h = f"O {added_hydrogen[0]} {added_hydrogen[1]} {added_hydrogen[2]}"
    third_h  = f"H {oh_h[0]} {oh_h[1]} {oh_h[2]}"
    path = file_path
    output = f'coord_hidx_{idx-1}.xyz' 
    replace_h(path, idx, new_line, output)
    with open(f'coord_hidx_{idx-1}.xyz','a') as file:
        file.write(second_h + '\n')
        file.write(third_h + '\n') 
        print(f'Finished with OH Radical Positioning on Former H with Idx {idx-1}. Saved as XYZ File in {output}.')

#
# Generating Product Geometries 
# 

# Aromatic C-H Bonds

for ch_bond,hydrogen in zip(ch_bonds, new_H):
    idx = ch_bond[1]+2
    path = file_path
    output = f'Products/product_idx_{idx-1}.xyz'
    delete_h(path, idx, output)

# Aldehyde C-H Bonds

for ch_bond,hydrogen in zip(aldehyde_bonds, new_H_aldehyde):
    idx = ch_bond[1]+2
    path = file_path
    output = f'Products/product_idx_{idx-1}.xyz'
    delete_h(path, idx, output)

# N-H Bonds

for nh_bond,hydrogen in zip(nh_bonds, new_H_nh):
    idx = nh_bond[0]+2
    path = file_path
    output = f'Products/product_idx_{idx-1}.xyz'
    delete_h(path, idx, output) 

# O-H Bonds

for oh_bond,hydrogen in zip(oh_bonds, new_H_oh):
    idx = oh_bond[0]+2
    path = file_path
    output = f'Products/product_idx_{idx-1}.xyz'
    delete_h(path, idx, output)

#
# Vector Plotting
# 

print("########################################################")
print("                 Vector Visualisation                   ")
print("########################################################")

ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Coordinate')

ax.set_xlim([-5,5])
ax.set_ylim([-5,5])
ax.set_zlim([-5,5])

legend1 = mpatches.Patch(color='red', label='Initial C-H Bonds')
legend2 = mpatches.Patch(color='green', label='Elongated C-H Bonds')
legend3 = mpatches.Patch(color='blue', label='Additional Oxygen Atoms')
legend4 = mpatches.Patch(color='cyan', label='Additional H Atoms in OH Bonds')

visualize = input("Do you want to plot the C-H, H-O and O-H vectors? yes/no:").lower()
if visualize == "yes":
    plt.legend(handles=[legend1,legend2,legend3,legend4])
    plt.show()
    print("########################################################")
    print("                  END OF THE PROGRAM                    ")
    print("########################################################")
else:
    print("########################################################")
    print("                  END OF THE PROGRAM                    ")
    print("########################################################")
