import amino_acids as aa
import aa_ringdata as rings
import numpy as np
import networkx as nx
from scipy.linalg import expm3, norm


###############################################################################
# graph setup for amino acids and coordinate assignments for conserved atoms  #
###############################################################################
#Create graph for amino acid from library data 
def create_tree_for_aa(aa_name):
    aa_tree = nx.Graph()
    for atom in aa.atomdata[aa_name]:
        aa_tree.add_node(atom[0],name=atom[1],element=atom[2],hybrid=atom[3])
    aa_tree.add_edges_from(aa.bonddata[aa_name])
    return aa_tree

#Assign a given set of coordinates to an amino acid
def add_coordinates_to_residue(aa_tree , coords):
    for atom in aa_tree.nodes():
        aa_tree.node[atom]['XYZ'] = coords[atom-1]

#Find all atoms that are the side chain
def find_subgraph_sidechain(aa_tree):
    atoms_sidechain = list()
    atoms_backbones = ['N' , 'H' , 'H1' , 'H2' , 'H3' , 'CA' , 'HA' , 'C' ,
                       'O' , 'OXT']
    for atom in aa_tree.nodes():
        if aa_tree.node[atom]['name'] not in atoms_backbones:
            atoms_sidechain.append(atom)
    return nx.subgraph(aa_tree , atoms_sidechain)

#Find all atoms that two amino acids have in common
def find_identical_atoms(tree1 , tree2):
    identical_atoms = list()
    for atom1 in tree1.nodes():
        for atom2 in tree2.nodes():
            if ((tree1.node[atom1]['name'] == tree2.node[atom2]['name']) and 
                (tree1.node[atom1]['hybrid'] == tree2.node[atom2]['hybrid'])):
                    if tree1.node[atom1]['element'] != 'H':
                        identical_atoms.append((atom1 , atom2))
                    else:
                        host1 = find_all_bonded_neighbours(tree1 , atom1)
                        host2 = find_all_bonded_neighbours(tree2 , atom2)                       
                        if (len(host1) == 1) and (len(host2) == 1):
                            host_atom1 = host1[0]
                            host_atom2 = host2[0]
                            if ((tree1.node[host_atom1]['name'] == tree2.node[host_atom2]['name']) and 
                                (tree1.node[host_atom1]['hybrid'] == tree2.node[host_atom2]['hybrid'])):
                                   identical_atoms.append((atom1 , atom2))                                   
    return identical_atoms

#Assign the coordinates of the backbone atoms (this will always be unchanged!)
def assign_coordinates_backbone(aa_tree1 , aa_tree2):
    atoms_backbone = ['N' , 'H' , 'H1' , 'H2' , 'H3' , 'CA' , 'HA' , 'C' ,
                       'O' , 'OXT']   
    for atom1 in aa_tree1.nodes():
        if aa_tree1.node[atom1]['name'] in atoms_backbone:
            for atom2 in aa_tree2.nodes():
                if aa_tree1.node[atom1]['name'] == aa_tree2.node[atom2]['name']:
                    aa_tree2.node[atom2]['XYZ'] = aa_tree1.node[atom1]['XYZ']

#Assign the coordinates for all identical atoms (keep as much data as we can!)
def assign_coordinates_identical_atoms(aa_tree1 , aa_tree2):
    tree1 = find_subgraph_sidechain(aa_tree1)
    tree2 = find_subgraph_sidechain(aa_tree2)
    identical_atoms = find_identical_atoms(tree1 , tree2)
    for match in identical_atoms:
        aa_tree2.node[match[1]]['XYZ'] = aa_tree1.node[match[0]]['XYZ']

#Find all atoms without coordinates
def find_atoms_without_coords(aa_tree):
    atoms_unassigned = list()
    assigned_coordinates = nx.get_node_attributes(aa_tree ,'XYZ')
    for atom in aa_tree.nodes():
        try:
            assigned_coordinates[atom]
        except KeyError:
            atoms_unassigned.append(atom)
    return atoms_unassigned


def get_number_from_name(aa_tree , atomname):
    return filter(lambda (n, d): d['name'] == atomname, aa_tree.nodes(data=True))[0][0]

def get_bondlength(aa_tree , atom1 , atom2):
    element1 = aa_tree.node[atom1]['element']
    element2 = aa_tree.node[atom2]['element']
    hybrid1 = aa_tree.node[atom1]['hybrid']
    hybrid2 = aa_tree.node[atom2]['hybrid']
    if (element1 == "C") and (element2 == "C"):
        if (hybrid1 == "sp2") and (hybrid2 == "sp2"):
            #only in benzene in amino acids
            return 1.396
        elif (hybrid1 == "sp2") or (hybrid2 == "sp2"):
            return 1.500
        else:
            return 1.540
    elif (element1 == "H") or (element2 == "H"):
        if (element1 == "C") or (element2 == "C"):
            return 1.090
        elif (element1 == "O") or (element2 == "O"):
            return 0.960
        elif (element1 == "N") or (element2 == "N"):
            return 1.010
        elif (element1 == "S") or (element2 == "S"):
            return 1.340
    elif (element1 == "N") or (element2 == "N"):
        if (hybrid1 == "sp3") and (hybrid2 == "sp3"):
            return 1.470
        elif (hybrid1 == "sp2") and (hybrid2 == "sp2"):
            return 1.335
        else:
            return 1.449
    elif (element1 == "O") or (element2 == "O"):
        if (hybrid1 == "sp3") and (hybrid2 == "sp3"):
            return 1.440
        elif (hybrid1 == "sp2") and (hybrid2 == "sp2"):
            return 1.229
        else:
            return 1.440
    elif (element1 == "S") or (element2 == "S"):
        if (hybrid1 == "sp3") and (hybrid2 == "sp3"):
            return 1.819
        elif (hybrid1 == "sp2") and (hybrid2 == "sp2"):
            return 1.553
        else:
            return 1.800
    else:
        return 1.00

###############################################################################
# Vector, matrix and search functions to create coordinates for atoms         #
###############################################################################
#Find neighbours that are bonded to atoms   
def find_all_bonded_neighbours(aa_tree , atom):
    bonds = aa_tree.edges(nbunch = atom)
    neighbouring_atoms = list()
    for bond in bonds:
        if bond[0] != atom:
            neighbouring_atoms.append(bond[0])
        else:
            neighbouring_atoms.append(bond[1])
    return neighbouring_atoms

#Find bonded atoms that have assigned coordinates
def find_neighbours_with_coords(aa_tree , atom):
    neighbours = find_all_bonded_neighbours(aa_tree , atom)
    neighbours_with_coords = list()
    assigned_coordinates = nx.get_node_attributes(aa_tree , 'XYZ')
    for neighbour in neighbours:
        try:
            assigned_coordinates[neighbour]
            neighbours_with_coords.append(neighbour)
        except KeyError:
            pass
    return neighbours_with_coords

#Derive geometry from element and hybridisation
def vacancy_from_hybrid(aa_tree , atom):
    hybrid = aa_tree.node[atom]['hybrid']
    element = aa_tree.node[atom]['element']
    if element == 'C':
        if hybrid == 'sp2':
            return (False , True , False)
        elif hybrid == 'sp3':
            return (True , False , False)
        else:
            raise NameError('C should be sp2 or sp3')
    elif element == 'O':
        if hybrid == 'sp2':
            return (False , False , True)
        elif hybrid == 'sp3':
            return (False , True , False)
        else:
            raise NameError('O should be sp2 or sp3')
    elif element == 'N':
        if hybrid == 'sp2':
            return (False , True , False)
        elif hybrid == 'sp3':
            return (True , False , False)
        else:
            raise NameError('N should be sp2 or sp3')
    elif element == 'H':
        return (False , False , True)
    elif element == 'S':
        if hybrid == 'sp3':
            return (True , False , True)
        else:
            raise NameError('S should be sp3')
    else:
        raise NameError('Element not programmed')

#Normalisation of vectors
def normalise_vector(vector):
    return vector/np.linalg.norm(vector)

#Find any perpendicular vector
def find_normal_vector(vector):
    #define a perpendicular vector (at this stage any perpendicular vector is fine)
    if (vector[0] != 0.0 and vector[1] != 0.0): 
        perp_vector = np.asarray([vector[1],-vector[0],0.0])
    elif (vector[0] != 0.0 and vector[2] != 0.0):
        perp_vector = np.asarray([vector[2],0.0,-vector[0]])
    elif (vector[1] != 0.0 and vector[2] != 0.0):
        perp_vector = np.asarray([0.0,vector[2],-vector[1]])
    elif (vector == np.asarray([0.0,0.0,0.0])):
        raise NameError('Bonded atoms are on top of each other')
    else:
        if (vector[0] != 0.0):
            perp_vector = np.asarray([0.0,1.0,0.0])
        else:
            perp_vector = np.asarray([1.0,0.0,0.0])
    return normalise_vector(perp_vector)


#Rotate vector around axis
def rotate_vector(vector , axis , angle):
    theta = np.radians(angle)
    rotmat = expm3(np.cross(np.eye(3), axis/norm(axis)*theta))
    return np.dot(rotmat,vector)

def Gram_Schmidt_vector(vec1 , vec2):
    return vec2 - (np.dot(vec1,vec2)/np.dot(vec1,vec1)) * vec1

def angle_between_vectors(v0 , v1):
    if ((abs(v0[0] - v1[0]) <= 1.0e-6) and (abs(v0[1] - v1[1]) <= 1.0e-6) and 
        (abs(v0[2] - v1[2]) <= 1.0e-6)):
        return 0.000
    else:
        plane_normal = normalise_vector(np.cross(v0,v1))
        angle= np.math.atan2(np.dot(np.cross(v1,v0),plane_normal),np.dot(v0,v1))
        return np.degrees(angle)


###############################################################################
# Functions relating to cyclic structures                                     #
###############################################################################     
#Find rings in amino acid
def check_for_rings(aa_tree):
    length = len(nx.cycle_basis(aa_tree))
    if length == 0:
        return False
    else:
        return True

#Check if all atoms in the ring are sp2 - then we have a planar ring
def check_planarity_ring(aa_tree , ring_atoms):
    for atom in ring_atoms:
        if aa_tree.node[atom]['hybrid'] != 'sp2':
            return False
        else:
            return True

#Are there any assigned coordinates in the ring?
def get_existing_coords_ring(aa_tree , ring_atoms):
    existing_coords = list()
    ncoords = 0
    for atom in ring_atoms:
        try:
            aa_tree.node[atom]['XYZ']
            existing_coords.append(atom)
            ncoords += 1
        except KeyError:
            pass
    return ncoords , existing_coords

def get_atom_list_upto_ring(aa_tree , first_atom):
    atoms_unassigned = find_atoms_without_coords(aa_tree)
    chain_atoms = list()
    for atom in atoms_unassigned:
        if atom <= first_atom:
            chain_atoms.append(atom)
    return chain_atoms

def assign_coordinates_ring(aa_tree , resname):
    #we assume that only the first C atom is assigned
    CA = get_number_from_name(aa_tree , "CA")
    CB = get_number_from_name(aa_tree , "CB")
    CG = get_number_from_name(aa_tree , "CG")
    CA_coords = np.asarray(aa_tree.node[CA]['XYZ'])
    CB_coords = np.asarray(aa_tree.node[CB]['XYZ'])
    CG_coords = np.asarray(aa_tree.node[CG]['XYZ'])
    CB_lib = np.asarray(rings.axis[resname][0])
    CG_lib = np.asarray(rings.axis[resname][1])
    diff_CB = CB_coords - CB_lib
    CG_lib = CG_lib + diff_CB
    CA_CB = CB_coords - CA_coords
    CB_CG_lib = CG_lib - CB_coords
    CB_CG_mol = CG_coords - CB_coords
    lib_vec = Gram_Schmidt_vector(CA_CB , CB_CG_lib)
    mol_vec = Gram_Schmidt_vector(CA_CB , CB_CG_mol)
    angle = angle_between_vectors(lib_vec , mol_vec)
    for atom in rings.ring_atoms[resname].keys():
        number = get_number_from_name(aa_tree , atom)
        atom_lib = np.asarray(rings.ring_atoms[resname][atom])
        atom_lib = atom_lib + diff_CB
        position_lib = atom_lib - CB_coords
        new_position = rotate_vector(position_lib , CA_CB , angle)
        aa_tree.node[number]['XYZ'] = new_position + CB_coords
    return

def position_CG_for_TRP(old_tree , new_tree, oldname):
    if oldname == 'ILE' or oldname == 'CILE' or oldname == 'NILE':
        old_CG1 = get_number_from_name(old_tree , 'CG1')
        new_CG = get_number_from_name(new_tree , 'CG')
        new_tree.node[new_CG]['XYZ'] = old_tree.node[old_CG1]['XYZ']
    elif oldname == 'THR' or oldname == 'CTHR' or oldname == 'NTHR':
        old_CB = get_number_from_name(old_tree , 'CB')
        old_HB = get_number_from_name(old_tree , 'HB')
        new_CG = get_number_from_name(new_tree , 'CG')
        bond_old = np.asarray(old_tree.node[old_HB]['XYZ']) - np.asarray(old_tree.node[old_CB]['XYZ'])
        new_bond = 1.50 * normalise_vector(bond_old)
        new_tree.node[new_CG]['XYZ'] = np.asarray(old_tree.node[old_CB]['XYZ']) + new_bond
    elif oldname == 'VAL' or oldname == 'CVAL' or oldname == 'NVAL':
        old_CG2 = get_number_from_name(old_tree , 'CG2')
        new_CG = get_number_from_name(new_tree , 'CG')
        new_tree.node[new_CG]['XYZ'] = old_tree.node[old_CG2]['XYZ']        
    return

def check_special_assignment(resname1 , resname2):
    if len(resname1) == 4:
        resname1_ = resname1[1:]
        resname2_ = resname2[1:]
    else:
        resname1_ = resname1
        resname2_ = resname2
    if resname1_ in ['HIS' , 'HID' , 'HIE' , 'HIP']:
        resname1_ = 'HIS'
    if resname2_ in ['HIS' , 'HID' , 'HIE' , 'HIP']:
        resname2_ = 'HIS'   
    if (((resname1_ == 'PHE') and (resname2_ == 'TYR')) or 
        ((resname2_ == 'PHE') and (resname1_ == 'TYR')) or
        ((resname1_ == 'HIS') and (resname2_ == 'TRP')) or
        ((resname2_ == 'TRP') and (resname1_ == 'HIS'))):
        return True
    else:
        return False
     
def special_assignment(aa_trees , resnames):
    aa_tree_old , aa_tree_new = aa_trees
    oldname , newname = resnames
    if len(oldname) == 4:
        oldname_ = oldname[1:]
        newname_ = newname[1:]
    else:
        oldname_ = oldname
        newname_ = newname
    if (((oldname_ == 'PHE') and (newname_ == 'TYR')) or 
        ((oldname_ == 'TYR') and (newname_ == 'PHE'))):
        return #basically identical, so the assignment is done already
    #must be TRP <--> HIS in some form, want to assign the five memebered ring
    else:
        if oldname_ == 'TRP':
    #assign all atoms correctly that we have in both/that will be near identical
            old_CD2 = get_number_from_name(aa_tree_old , "CD2")
            old_CD1 = get_number_from_name(aa_tree_old , "CD1")
            old_CE2 = get_number_from_name(aa_tree_old , "CE2")
            old_NE1 = get_number_from_name(aa_tree_old , "NE1")
            old_HD1 = get_number_from_name(aa_tree_old , "HD1")
            old_HE1 = get_number_from_name(aa_tree_old , "HE1")           
            new_ND1 = get_number_from_name(aa_tree_new , "ND1")
            new_CD2 = get_number_from_name(aa_tree_new , "CD2")
            new_CE1 = get_number_from_name(aa_tree_new , "CE1")
            new_NE2 = get_number_from_name(aa_tree_new , "NE2")
            new_HD2 = get_number_from_name(aa_tree_new , "HD2")
            new_HE2 = get_number_from_name(aa_tree_new , "HE2")            
            aa_tree_new.node[new_ND1]['XYZ'] = aa_tree_old.node[old_CD2]['XYZ']
            aa_tree_new.node[new_CD2]['XYZ'] = aa_tree_old.node[old_CD1]['XYZ']
            aa_tree_new.node[new_CE1]['XYZ'] = aa_tree_old.node[old_CE2]['XYZ']
            aa_tree_new.node[new_NE2]['XYZ'] = aa_tree_old.node[old_NE1]['XYZ']
            aa_tree_new.node[new_HD2]['XYZ'] = aa_tree_old.node[old_HD1]['XYZ']
            aa_tree_new.node[new_HE2]['XYZ'] = aa_tree_old.node[old_HE1]['XYZ']                  
        else:
    #assign all atoms correctly that we have in both/that will be near identical           
            old_ND1 = get_number_from_name(aa_tree_old , "ND1")
            old_CD2 = get_number_from_name(aa_tree_old , "CD2")
            old_CE1 = get_number_from_name(aa_tree_old , "CE1")
            old_NE2 = get_number_from_name(aa_tree_old , "NE2")
            old_HD2 = get_number_from_name(aa_tree_old , "HD2")
            new_CD2 = get_number_from_name(aa_tree_new , "CD2")
            new_CD1 = get_number_from_name(aa_tree_new , "CD1")
            new_CE2 = get_number_from_name(aa_tree_new , "CE2")
            new_NE1 = get_number_from_name(aa_tree_new , "NE1")
            new_HD1 = get_number_from_name(aa_tree_new , "HD1")          
            aa_tree_new.node[new_CD2]['XYZ'] = aa_tree_old.node[old_ND1]['XYZ']
            aa_tree_new.node[new_CD1]['XYZ'] = aa_tree_old.node[old_CD2]['XYZ']
            aa_tree_new.node[new_CE2]['XYZ'] = aa_tree_old.node[old_CE1]['XYZ']
            aa_tree_new.node[new_NE1]['XYZ'] = aa_tree_old.node[old_NE2]['XYZ']
            aa_tree_new.node[new_HD1]['XYZ'] = aa_tree_old.node[old_HD2]['XYZ']
            try:
                old_HE2 = get_number_from_name(aa_tree_old , "HE2")
                new_HE1 = get_number_from_name(aa_tree_new , "HE1") 
                aa_tree_new.node[new_HE1]['XYZ'] = aa_tree_old.node[old_HE2]['XYZ']
            except IndexError: 
                pass
            #assign ring atoms for TRP
            old_HE1 = get_number_from_name(aa_tree_old , "HE1")
            H_coords = np.asarray(aa_tree_old.node[old_HE1]['XYZ'])
            #HE1 is essentially where the next C is
            host_coords = np.asarray(aa_tree_new.node[new_CE2]['XYZ'])
            bond_vector = H_coords - host_coords
            new_bond_vector = 1.39 * normalise_vector(bond_vector)
            new_CZ2 = get_number_from_name(aa_tree_new , "CZ2")
            aa_tree_new.node[new_CZ2]['XYZ'] = host_coords + new_bond_vector
            #now go for CE3 --> rotation of CZ2
            CE2_CD2 = np.asarray(aa_tree_new.node[new_CD2]['XYZ']) - np.asarray(aa_tree_new.node[new_CE2]['XYZ'])         
            CE2_CZ2 = np.asarray(aa_tree_new.node[new_CZ2]['XYZ']) - np.asarray(aa_tree_new.node[new_CE2]['XYZ'])
            normal = Gram_Schmidt_vector(CE2_CD2 , CE2_CZ2)
            midpoint = 0.5 * (np.asarray(aa_tree_new.node[new_CD2]['XYZ']) + np.asarray(aa_tree_new.node[new_CE2]['XYZ']))
            CZ2_mid = np.asarray(aa_tree_new.node[new_CZ2]['XYZ']) - midpoint
            CE3_mid = rotate_vector(CZ2_mid , normal , 180.0)
            new_CE3 = get_number_from_name(aa_tree_new , "CE3")
            aa_tree_new.node[new_CE3]['XYZ'] = midpoint + CE3_mid
            #now do the last two by rotations
            new_CZ3 = get_number_from_name(aa_tree_new , "CZ3")
            new_CH2 = get_number_from_name(aa_tree_new , "CH2")
            CE3_CZ2 = np.asarray(aa_tree_new.node[new_CZ2]['XYZ']) - np.asarray(aa_tree_new.node[new_CE3]['XYZ'])
            CE3_CD2 = np.asarray(aa_tree_new.node[new_CD2]['XYZ']) - np.asarray(aa_tree_new.node[new_CE3]['XYZ'])
            CZ2_CE2 = np.asarray(aa_tree_new.node[new_CE2]['XYZ']) - np.asarray(aa_tree_new.node[new_CZ2]['XYZ'])
            perp_CH2 = Gram_Schmidt_vector(CE3_CZ2 , -CZ2_CE2)
            par_CH2 = -CZ2_CE2 - perp_CH2
            perp_CZ3 = Gram_Schmidt_vector(CE3_CZ2 , -CE3_CD2)
            par_CZ3 = -CE3_CD2 - perp_CZ3
            aa_tree_new.node[new_CH2]['XYZ'] = aa_tree_new.node[new_CZ2]['XYZ'] + perp_CH2 - par_CH2
            aa_tree_new.node[new_CZ3]['XYZ'] = aa_tree_new.node[new_CE3]['XYZ'] + perp_CZ3 - par_CZ3
        return
 
def mutate_GLY_chiral(aa_tree_old , aa_tree_new):
    N = np.asarray(aa_tree_old.node[get_number_from_name(aa_tree_old,'N')]['XYZ'])
    CA = np.asarray(aa_tree_old.node[get_number_from_name(aa_tree_old,'CA')]['XYZ'])
    HA2 = np.asarray(aa_tree_old.node[get_number_from_name(aa_tree_old,'HA2')]['XYZ'])
    HA3 = np.asarray(aa_tree_old.node[get_number_from_name(aa_tree_old,'HA3')]['XYZ'])
    C = np.asarray(aa_tree_old.node[get_number_from_name(aa_tree_old,'C')]['XYZ'])
    normal = normalise_vector(np.cross((N - CA) , (N - C)))
    HA2_HA3 = normalise_vector(HA2 - HA3)
    dotp = np.dot(normal , HA2_HA3)
    CA = get_number_from_name(aa_tree_old , 'CA')
    HA2 = get_number_from_name(aa_tree_old , 'HA2')
    HA3 = get_number_from_name(aa_tree_old , 'HA3')
    HA = get_number_from_name(aa_tree_new , 'HA')
    CB = get_number_from_name(aa_tree_new , 'CB')
    if dotp > 0.0:       
        aa_tree_new.node[HA]['XYZ'] = aa_tree_old.node[HA2]['XYZ']
        bond = np.asarray(aa_tree_old.node[HA3]['XYZ']) - np.asarray(aa_tree_old.node[CA]['XYZ'])
        aa_tree_new.node[CB]['XYZ'] = np.asarray(aa_tree_old.node[CA]['XYZ']) + 1.540 * bond
    else:
        aa_tree_new.node[HA]['XYZ'] = aa_tree_old.node[HA3]['XYZ']
        bond = np.asarray(aa_tree_old.node[HA2]['XYZ']) - np.asarray(aa_tree_old.node[CA]['XYZ'])
        aa_tree_new.node[CB]['XYZ'] = np.asarray(aa_tree_old.node[CA]['XYZ']) + 1.540 * bond        
    return

def mutate_HYP_from_PRO(aa_tree_old , aa_tree_new):
    CG_old  = get_number_from_name(aa_tree_old , 'CG' )
    HG2_old = get_number_from_name(aa_tree_old , 'HG2')
    HG3_old = get_number_from_name(aa_tree_old , 'HG3')
    CG_new  = get_number_from_name(aa_tree_new , 'CG' )
    HG_new  = get_number_from_name(aa_tree_new , 'HG' )
    OD1_new = get_number_from_name(aa_tree_new , 'OD1')
    aa_tree_new.node[HG_new]['XYZ'] = aa_tree_old.node[HG2_old]['XYZ']
    CG_HG3 = np.asarray(aa_tree_old.node[HG3_old]['XYZ']) - np.asarray(aa_tree_old.node[CG_old]['XYZ'])
    CG_coords = aa_tree_new.node[CG_new]['XYZ']
    aa_tree_new.node[OD1_new]['XYZ'] = CG_coords + 1.440 * normalise_vector(CG_HG3)
    return
         

def check_chirality_ILE(aa_tree_new):
    CB  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CB')]['XYZ'])
    CA  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CA')]['XYZ'])
    HB  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'HB')]['XYZ'])
    CG1 = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CG1')]['XYZ'])
    CG2 = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CG2')]['XYZ'])
    normal = normalise_vector(np.cross((CG1 - CB) , (CA - CB)))
    CG2_HB = normalise_vector(CG2 - HB)
    if np.dot(normal , CG2_HB) > 0.0: #wronmg chirality
        #switch HB and CG2, and then delete H
        HB = get_number_from_name(aa_tree_new , 'HB')
        CG2 = get_number_from_name(aa_tree_new , 'CG2')
        CB = get_number_from_name(aa_tree_new , 'CB')
        HB_coords = np.asarray(aa_tree_new.node[HB]['XYZ'])
        CB_coords = np.asarray(aa_tree_new.node[CB]['XYZ'])
        CG2_coords = np.asarray(aa_tree_new.node[CG2]['XYZ'])        
        bond_CC = np.linalg.norm(CG2_coords - CB_coords)
        bond_CH = np.linalg.norm(HB_coords - CB_coords)
        aa_tree_new.node[CG2]['XYZ'] = CB_coords + bond_CC * normalise_vector(HB_coords - CB_coords)
        aa_tree_new.node[HB]['XYZ'] = CB_coords + bond_CH * normalise_vector(CG2_coords - CB_coords)
        try:
            del aa_tree_new.node[CG2 + 1]['XYZ']
            del aa_tree_new.node[CG2 + 2]['XYZ']
            del aa_tree_new.node[CG2 + 3]['XYZ']
        except KeyError:
            pass
    return

def check_chirality_THR(aa_tree_new):
    CB  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CB')]['XYZ'])
    CA  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CA')]['XYZ'])
    HB  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'HB')]['XYZ'])
    OG1 = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'OG1')]['XYZ'])
    CG2 = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CG2')]['XYZ'])
    normal = normalise_vector(np.cross((CG2 - CB) , (CA - CB)))
    OG1_HB = normalise_vector(OG1 - HB)
    if np.dot(normal , OG1_HB) < 0.0: #wronmg chirality
        #switch HB and OG1, and then delete H
        HB = get_number_from_name(aa_tree_new , 'HB')
        OG1 = get_number_from_name(aa_tree_new , 'OG1')
        CB = get_number_from_name(aa_tree_new , 'CB')
        HB_coords = np.asarray(aa_tree_new.node[HB]['XYZ'])
        CB_coords = np.asarray(aa_tree_new.node[CB]['XYZ'])
        OG1_coords = np.asarray(aa_tree_new.node[OG1]['XYZ'])        
        bond_CO = np.linalg.norm(OG1_coords - CB_coords)
        bond_CH = np.linalg.norm(HB_coords - CB_coords)
        aa_tree_new.node[OG1]['XYZ'] = CB_coords + bond_CO * normalise_vector(HB_coords - CB_coords)
        aa_tree_new.node[HB]['XYZ'] = CB_coords + bond_CH * normalise_vector(OG1_coords - CB_coords)
        try:
            del aa_tree_new.node[OG1 + 1]['XYZ']
        except KeyError:
            pass
    return

def check_chirality_HYP(aa_tree_new):   
    CG  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CG')]['XYZ'])
    CB  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CB')]['XYZ'])
    CD  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'CD')]['XYZ'])
    OD1 = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'OD1')]['XYZ'])
    HG  = np.asarray(aa_tree_new.node[get_number_from_name(aa_tree_new,'HG')]['XYZ'])
    normal = normalise_vector(np.cross((CD - CG) , (CB - CG)))
    OD1_HG = normalise_vector(OD1 - HG)
    if np.dot(normal , OD1_HG) > 0.0: #wrong chirality
        #switch HG and OD1, and then delete H
        HG = get_number_from_name(aa_tree_new , 'HG')
        OD1 = get_number_from_name(aa_tree_new , 'OD1')
        CG = get_number_from_name(aa_tree_new , 'CG')
        HG_coords = np.asarray(aa_tree_new.node[HG]['XYZ'])
        CG_coords = np.asarray(aa_tree_new.node[CG]['XYZ'])
        OD1_coords = np.asarray(aa_tree_new.node[OD1]['XYZ'])        
        bond_CO = np.linalg.norm(OD1_coords - CG_coords)
        bond_CH = np.linalg.norm(HG_coords - CG_coords)
        aa_tree_new.node[OD1]['XYZ'] = CG_coords + bond_CO * normalise_vector(HG_coords - CG_coords)
        aa_tree_new.node[HG]['XYZ'] = CG_coords + bond_CH * normalise_vector(OD1_coords - CG_coords)
        try:
            del aa_tree_new.node[OD1 + 1]['XYZ']
        except KeyError:
            pass
    return


def create_proline_coords(aa_tree):
    CA = get_number_from_name(aa_tree , "CA")
    N  = get_number_from_name(aa_tree , "N")
    CB = get_number_from_name(aa_tree , "CB")
    CG = get_number_from_name(aa_tree , "CG")
    CD = get_number_from_name(aa_tree , "CD")
    CA_coords = np.asarray(aa_tree.node[CA]['XYZ'])
    N_coords  = np.asarray(aa_tree.node[N]['XYZ'])
    CB_coords  = np.asarray(aa_tree.node[CB]['XYZ'])

    N_CA = CB_coords - N_coords
    N_CB = CA_coords - N_coords
    perp_1 = Gram_Schmidt_vector(N_CB , N_CA)
    midpoint1 = 0.5 * (N_coords + CA_coords)
    midpoint2 = midpoint1 + perp_1
    CD_coords = midpoint2 + (midpoint2 - CB_coords)
    perp_2 = np.cross(perp_1 , midpoint2 - CB_coords)
    CG_coords = midpoint2 + 0.5 * perp_1 + 0.33 * perp_2
    aa_tree.node[CG]['XYZ'] = CG_coords
    aa_tree.node[CD]['XYZ'] = CD_coords    
    return   


###############################################################################
# Create coordinates given we have one bonded neighbour with coordinates(host)#
###############################################################################
def create_new_coords_one_host(aa_tree , new_atom , host_atom):
    neighbours_host = find_neighbours_with_coords(aa_tree , host_atom)
    host_coords = np.asarray(aa_tree.node[host_atom]['XYZ'])
    neighbours_coords = [np.asarray(aa_tree.node[x]['XYZ']) for x in neighbours_host]
    bondlength = get_bondlength(aa_tree , new_atom , host_atom)
    tetrahedral , planar , terminal = vacancy_from_hybrid(aa_tree , host_atom)
    if len(neighbours_host) == 0:
        del aa_tree.node[host_atom]['XYZ']
    elif tetrahedral:
        if len(neighbours_host) == 1:
            #find vector along bond
            bond_vector = np.asarray(host_coords) - np.asarray(neighbours_coords[0])
            perp_unit = find_normal_vector(bond_vector)
            new_bond_vector = -(rotate_vector(bond_vector , perp_unit , 109.5))
            new_bond_vector = bondlength * normalise_vector(new_bond_vector)
            aa_tree.node[new_atom]['XYZ'] = host_coords + new_bond_vector
        elif len(neighbours_host) == 2:
            neighbour1 = np.asarray(neighbours_coords[0])
            neighbour2 = np.asarray(neighbours_coords[1])
            bond1 = np.asarray(neighbour1 - host_coords)
            bond2 = np.asarray(neighbour2 - host_coords)
            perp_1 = Gram_Schmidt_vector(bond2 , bond1)
            para_1 = bond1 - perp_1
            new_perp = rotate_vector(perp_1 , bond2 , -120.0)
            new_bond = new_perp + para_1
            new_bond = bondlength * normalise_vector(new_bond)
            aa_tree.node[new_atom]['XYZ'] = host_coords + new_bond
        elif len(neighbours_host) == 3:
            centroid = np.asarray([0.0 , 0.0 , 0.0])
            for atom in neighbours_host:
                centroid += aa_tree.node[atom]['XYZ']
            centroid = np.asarray([centroid[0]/3.0 , centroid[1]/3.0 , centroid[2]/3.0])
            new_bond_vector= host_coords - centroid
            new_bond_vector = bondlength * normalise_vector(new_bond_vector)
            aa_tree.node[new_atom]['XYZ'] = host_coords + new_bond_vector
        else:
            raise NameError('''Too few or many neighbours known for creating a 
                               tetrahedral centre''')
    elif planar:
        if len(neighbours_host) == 1:
            bond_vector = host_coords - neighbours_coords[0]
            perp_unit = find_normal_vector(bond_vector)
            new_bond_vector = -(rotate_vector(bond_vector , perp_unit , 120.0))
            new_bond_vector = bondlength * normalise_vector(new_bond_vector)
            aa_tree.node[new_atom]['XYZ'] = host_coords + new_bond_vector
        elif len(neighbours_host) == 2:
            bond1 = host_coords - neighbours_coords[0]
            bond2 = host_coords - neighbours_coords[1]
            perp_1 = Gram_Schmidt_vector(bond2 , bond1)
            para_1 = bond1 - perp_1
            new_bond_vector = perp_1 - para_1
            new_bond_vector = bondlength * normalise_vector(new_bond_vector)
            aa_tree.node[new_atom]['XYZ'] = host_coords + new_bond_vector
        else:
            raise NameError('''Too few or many neighbours known for creating a 
                               planar centre''')
    else:
        raise NameError('Cannot create new bond with terminal atom')

###############################################################################
# Create coordinates given we have more than one host                         #
###############################################################################
def remove_multiple_hosts(aa_tree , hosts):
    for host in hosts:
        if aa_tree.node[host]["element"] == "H":
            del aa_tree.node[host]['XYZ']
    if len(hosts) == 2:
        atom1 = hosts[0]
        atom2 = hosts[1]
        if atom1 > atom2:
            try:
                del aa_tree.node[atom1]['XYZ']
            except KeyError:
                pass
        else:
            try:
                del aa_tree.node[atom2]['XYZ']
            except KeyError:
                pass                
    return
