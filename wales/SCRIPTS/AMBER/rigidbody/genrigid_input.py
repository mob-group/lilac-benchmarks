# SCRIPT TO DEFINE THE RIGIDIFICATION FOR PROTEINS
# by Konstantin Roeder
# creates input files rbodyconfig and coordsinirigid and additionally a log file rigid.log

# ! /usr/bin/env python

import sys
import time
import math

if len(sys.argv) > 1:
    print 'Script creates input files for RIGIDINIT keyword, needs inpcrd and pdb' + \
          'for RNA and DNA only post AMBER12 atom names'
    print
    pdb_inp = sys.argv[1]
    coords_inp = sys.argv[2]
    # defines the tolerance for the checks whether rigid bodies are linear
    try:
        tolerance = float(sys.argv[3])
    except IndexError:
        tolerance = 0.01
    # Pymol use; change to False if pymol is not installed/used
    try:
        if sys.argv[4] == 'pymol':
            pymol_check = True
    except IndexError:
        pymol_check = False

else:
    tolerance = 0.01
    pymol_check = False
    pdb_inp = raw_input('PDB file: ')
    coords_inp = raw_input('Coords input: ')

inp_sys = raw_input('Use script for peptides (1), RNA (2), DNA (3), or mixed (4)? ')
try:
    if int(inp_sys) == 1:
        protein_t = True
        rna_t = False
        dna_t = False
    elif int(inp_sys) == 2:
        protein_t = False
        rna_t = True
        dna_t = False
    elif int(inp_sys) == 3:
        protein_t = False
        rna_t = False
        dna_t = True
    else:
        protein_t = True
        rna_t = True
        dna_t = True
except ValueError:
    print 'Assume RNA and protein in mixed system.'
    protein_t = True
    rna_t = True
    dna_t = True

if pymol_check:
    import __main__

    __main__.pymol_argv = ['pymol', '-qix']  # Pymol: quiet and no GUI(internal and external)
    import pymol


# class containing the main methods for the script, new functionality should be embedded here
class protein():
    def __init__(self, atom_dic, res_dic):
        self.atom_dic = atom_dic  # dictionary of all atoms with x,y,z coordinates and res id and name and atom name
        self.res_dic = res_dic  # dictionary containing a list of atom ids for all residues

    def num_atoms(self):  # number of atoms
        return len(self.atom_dic)

    def num_res(self):  # number of residues
        return len(self.res_dic)

    # returns the name of a given residue
    def get_residue_name(self, residue):
        return atom_dic[res_dic[residue][0]][1]

    # prints the atoms in a specified residue(need number of residue)
    # includes the atom number, atom name and coordinates
    def get_atoms_in_res(self, residue):
        all_atoms = self.atom_dic.keys()
        atoms_in_res = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atoms_in_res.append(all_atoms[i - 1])
                residuename = atom[1]
        if atoms_in_res == []:
            return 'Residue not in molecule'
        print 'Residue:', residue, residuename
        print
        print 'Atoms:'
        for i in range(0, len(atoms_in_res)):
            j = atoms_in_res[i]
            atom = self.atom_dic[j]
            print j, atom[0], atom[3], atom[4], atom[5]
        return

    # returns a list of atoms in a residue
    def get_atom_list(self, residue):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atomlist.append(all_atoms[i - 1])
        return atomlist

    # allows the input of a centre from the user, is complemented by the selection of the CoM as centre
    def choose_centre(self):
        check = raw_input('Display atoms in residue (y/n)?  ')
        if check == 'y':
            while 1:  # acts as input mask --> once entered any number of residues can be viewed at
                residue_number = raw_input('Residue number: (enter x to leave) ')
                if residue_number == 'x':
                    print
                    break
                else:
                    try:
                        self.get_atoms_in_res(int(residue_number))
                        print
                    except ValueError:
                        pass
        while 1:  # loop for choice of molecule, only accepts atoms in the molecule
            centre_input = int(raw_input('Choose atom:  '))
            if centre_input > self.num_atoms() or centre_input < 1:
                print 'Not in molecule.'
            else:
                break
        return centre_input

    # transforms a given list of atoms into a list of residues
    # can be used to find the residues in a sphere from the output of the atoms in sphere method
    def transform_atomlist_to_reslist(self, atomlist):
        reslist = []
        for i in atomlist:
            atom_info = self.atom_dic[i]
            k = 0
            if reslist != []:
                for j in range(0, len(reslist)):
                    if reslist[j] == atom_info[2]:
                        k = 1
            if k == 0:
                reslist.append(atom_info[2])
            else:
                pass
        return reslist

    def transform_reslist_to_atomlist(self, reslist):
        atomlist = []
        for i in reslist:
            atomlist += self.res_dic[i]
        return atomlist

    # get residue name from id
    def res_name_from_id(self, residue_id):
        atom_id = res_dic[residue_id][0]
        return atom_dic[atom_id][1]

    # checks for unknown residues/ligands, add standard residues here if needed
    def get_ligands_unknown_res(self):
        list_res_names = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY', 'GLX', 'HIS', 'HIE',
                          'HID', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                          'NALA', 'NARG', 'NASN', 'NASP', 'NASX', 'NCYS', 'NCYX', 'NGLN', 'NGLU', 'NGLY', 'NGLX',
                          'NHIS', 'NHIE', 'NHID', 'NHIP', 'NILE', 'NLEU', 'NLYS', 'NMET', 'NPHE', 'NPRO', 'NSER',
                          'NTHR', 'NTRP', 'NTYR', 'NVAL', 'CALA', 'CARG', 'CASN', 'CASP', 'CASX', 'CCYS', 'CCYX',
                          'CGLN', 'CGLU', 'CGLY', 'CGLX', 'CHIS', 'CHIE', 'CHID', 'CHIP', 'CILE', 'CLEU', 'CLYS',
                          'CMET', 'CPHE', 'CPRO', 'CSER', 'CTHR', 'CTRP', 'CTYR', 'CVAL', 'RAN', 'RCN', 'RGN', 'RUN',
                          'A', 'A3', 'A5', 'AN', 'C', 'C3', 'C5', 'CN', 'G', 'G3', 'G5', 'GN', 'U', 'U3', 'U5', 'UN',
                          'OHE', 'RA', 'RA3', 'RA5', 'RC', 'RC3', 'RC5', 'RG', 'RG3', 'RG5', 'RU', 'RU3', 'RU5',
                          'DA', 'DA5', 'DC', 'DC5', 'DG', 'DT', 'DT5', 'DA3', 'DC3', 'DG3', 'DT3',
                          'DAN', 'DCN', 'DGN', 'DTN']
        ligand_dic = {}
        for i in self.res_dic.keys():
            if self.res_name_from_id(i) not in list_res_names:
                ligand_dic[i] = self.res_dic[i]
        return ligand_dic

    # returns a list of atoms within a sphere of entered radius
    # needs an atom number as centre
    def get_atoms_in_sphere(self, centre_atom, radius):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            distance = (atom[3] - centre_atom[0]) ** 2 + (atom[4] - centre_atom[1]) ** 2 + (atom[5] - centre_atom[
                2]) ** 2
            if distance <= radius:
                atomlist.append(all_atoms[i - 1])
            else:
                pass
        return atomlist

    # returns the atoms within an ellipsoidal volume
    def get_atoms_in_ellipsoid(self, centre_atom, a, b, c):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            x1 = ((atom[3] - centre_atom[0]) / a) ** 2
            x2 = ((atom[4] - centre_atom[1]) / b) ** 2
            x3 = ((atom[5] - centre_atom[2]) / c) ** 2
            distance = x1 + x2 + x3
            if distance <= 1:
                atomlist.append(all_atoms[i - 1])
            else:
                pass
        return atomlist

    # returns a list of atoms in a residue
    def get_atom_list(self, residue):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atomlist.append(all_atoms[i - 1])
        return atomlist

    # gives the mass of any residue containing only C, N, O, H, D and S, other masses can only be used if they are added to the dictionary mass_dic below
    def res_mass(self, residue):
        atom_list = self.get_atom_list(residue)
        mass_dic = {'H': 1.007825, 'D': 2.014102, 'C': 12.0116, 'N': 14.00728, 'O': 15.99977, 'P': 30.973762,
                    'S': 32.076}
        mass = 0
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    mass += mass_dic[j]
        return mass

    # gives the coordinates for the mass_weighted centre of a residue
    def mass_weighted_centre(self, residue):
        atom_list = self.get_atom_list(residue)
        mass_dic = {'H': 1.007825, 'D': 2.014102, 'C': 12.0116, 'N': 14.00728, 'O': 15.99977, 'P': 30.973762,
                    'S': 32.076}
        X = 0
        Y = 0
        Z = 0
        mass_res = float(self.res_mass(residue))
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    X += mass_dic[j] * self.atom_dic[i][3]
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    Y += mass_dic[j] * self.atom_dic[i][4]
        for i in atom_list:
            for j in mass_dic.keys():
                if j == self.atom_dic[i][0][0]:
                    Z += mass_dic[j] * self.atom_dic[i][5]
        X = X / mass_res
        Y = Y / mass_res
        Z = Z / mass_res
        print 'Mass weighted centre: ', (X, Y, Z)
        return (X, Y, Z)

    # returns all residues of a list of type 'name'
    def get_all_res_name(self, res_input, name):
        res_list = []
        for i in res_input:
            if name == self.get_residue_name(i):
                res_list.append(i)
        return res_list

    # atom id from atom name given a residue id
    def get_atom_id_from_name(self, atom_name, res_num):
        for i in self.res_dic[res_num]:
            if self.atom_dic[i][0] == atom_name:
                return i

    # dot product of two vectors
    def dot_product(self, vector1, vector2):
        return (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2])

    # magnitude of a vector
    def magnitude(self, vector):
        return math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)

    # angle between three atoms
    def get_angle(self, (atom1, atom2, atom3)):
        vector1 = (self.atom_dic[atom2][3] - self.atom_dic[atom1][3], self.atom_dic[atom2][4] - self.atom_dic[atom1][4]
                   , self.atom_dic[atom2][5] - self.atom_dic[atom1][5])
        vector2 = (self.atom_dic[atom3][3] - self.atom_dic[atom1][3], self.atom_dic[atom3][4] - self.atom_dic[atom1][4]
                   , self.atom_dic[atom3][5] - self.atom_dic[atom1][5])
        dot_product = self.dot_product(vector1, vector2)
        magnitudes = self.magnitude(vector1) * self.magnitude(vector2)
        angle = math.acos(dot_product / magnitudes)
        angle = math.degrees(angle)
        if (angle > 0.00 and angle <= 180.00):
            return angle
        else:
            return angle - 180.00

    # checks a list of atoms whether it is linear, creates a triple list of all possible combinations and then computes angles
    def check_linear(self, atom_list, tolerance):
        triple_list = []  # list of atom triples to check for linearity
        for atom1 in atom_list:
            for atom2 in atom_list:
                for atom3 in atom_list:
                    if (atom1 != atom2) and (atom1 != atom3) and (atom2 != atom3):
                        triple_list.append((atom1, atom2, atom3))
                    else:
                        pass
        for triple in triple_list:
            if (self.get_angle(triple) > (0.00 + tolerance)) and (self.get_angle(triple) < (180.00 - tolerance)):
                return False
        return True

    # returns a dictionary including all the information needed for the localised rigidification scheme,
    # if an extension is made to the default options for rigidification, this must change the local_scheme list!
    # the input scheme must be a list of 0,1,... for the possible options within the rigidification for a given residue
    # any additional scheme can be added below that will be available on default if chosen
    def get_rigid_for_local(self, res_list, local_scheme, atoms_used):
        # local_scheme=['PRO','ARG',('HIS','HIE','HID','HIP'),'LYS','ASP','ASN','GLU','GLN','PHE','TYR','TRP']
        groups_res = {}
        i = 1
        # all following cases are for a single residue name with given patterns
        # when adding additional ones check, that atoms can only be chosen once as no atom can be part of more than one rigid body
        if local_scheme[0] == 1:
            l_res = self.get_all_res_name(res_list, 'PRO')
            for j in l_res:
                atom_list = []
                for k in ['C', 'O', 'N', 'CD', 'CG', 'CB', 'CA']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'PRO'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[0] == 2:
            l_res = self.get_all_res_name(res_list, 'PRO')
            for j in l_res:
                atom_list = []
                for k in ['N', 'CD', 'CG', 'CB', 'CA']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'PRO'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[1] == 1:
            l_res = self.get_all_res_name(res_list, 'ARG')
            for j in l_res:
                atom_list = []
                for k in ['NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'ARG'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[2] == 1:
            l_res = self.get_all_res_name(res_list, 'HIS')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HIS'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
            l_res = self.get_all_res_name(res_list, 'HIE')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HIE'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
            l_res = self.get_all_res_name(res_list, 'HID')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HID'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
            l_res = self.get_all_res_name(res_list, 'HIP')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'ND1', 'CE1', 'NE2', 'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'HIP'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[3] == 1:
            l_res = self.get_all_res_name(res_list, 'LYS')
            for j in l_res:
                atom_list = []
                for k in ['NZ', 'HZ1', 'HZ2', 'HZ3']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'LYS'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[4] == 1:
            l_res = self.get_all_res_name(res_list, 'ASP')
            for j in l_res:
                atom_list = []
                for k in ['CB', 'CG', 'OD1', 'OD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'ASP'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[5] == 1:
            l_res = self.get_all_res_name(res_list, 'ASN')
            for j in l_res:
                atom_list = []
                for k in ['CB', 'CG', 'OD1', 'ND2', 'HD21', 'HD22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'ASN'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[6] == 1:
            l_res = self.get_all_res_name(res_list, 'GLU')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD', 'OE1', 'OE2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'GLU'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[7] == 1:
            l_res = self.get_all_res_name(res_list, 'GLN')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD', 'OE1', 'NE2', 'HE21', 'HE22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'GLN'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[8] == 1:
            l_res = self.get_all_res_name(res_list, 'PHE')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'HZ', 'CE2', 'HE2', 'CD2', 'HD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'PHE'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[9] == 1:
            l_res = self.get_all_res_name(res_list, 'TYR')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD1', 'HD1', 'CE1', 'HE1', 'CZ', 'OH', 'CE2', 'HE2', 'CD2', 'HD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'TYR'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[10] == 1:
            l_res = self.get_all_res_name(res_list, 'TRP')
            for j in l_res:
                atom_list = []
                for k in ['CG', 'CD1', 'HD1', 'NE1', 'HE1', 'CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'CE3', 'HE3',
                          'CD2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'TRP'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[10] == 2:
            l_res = self.get_all_res_name(res_list, 'TRP')
            for j in l_res:
                atom_list1 = []
                atom_list2 = []
                for k in ['CG', 'CD1', 'HD1', 'NE1', 'HE1']:
                    atom_list1.append(self.get_atom_id_from_name(k, int(j)))
                for k in ['CE2', 'CZ2', 'HZ2', 'CH2', 'HH2', 'CZ3', 'CE3', 'HE3', 'CD2']:
                    atom_list2.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'TRP'] + atom_list1
                groups_res[i + 1] = [j, 'TRP'] + atom_list2
                atoms_used += atom_list1
                atoms_used += atom_list2
                i += 2
            l_res = []
        return (groups_res, atoms_used)

    def loc_NA(self, res_list, local_scheme, atoms_used):
        # local_scheme=[G_large, G_small, A_large, A_small, C_large, C_small, T_large, T_small, U_large, U_small])
        groups_res = {}
        i = 1
        # all following cases are for a single residue name with given patterns
        # when adding additional ones check, that atoms can only be chosen once as no atom can be part of more than one rigid body
        if local_scheme[0]:
            l_res = []
            for RA_name in ['G', 'G3', 'G5', 'GN', 'DG', 'DG3', 'DG5', 'DGN', 'RG', 'RG3', 'RG5', 'RGN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'N9', 'C8', 'N7', 'C5', 'C6', 'H8', 'O6', 'C2', 'N3', 'H1', 'N2', 'H21', 'H22']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'G'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[1]:
            l_res = []
            for RA_name in ['G', 'G3', 'G5', 'GN', 'DG', 'DG3', 'DG5', 'DGN', 'RG', 'RG3', 'RG5', 'RGN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N9', 'C8', 'N7', 'C5', 'C6', 'H8', 'C2', 'N3', 'H1', 'N2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'G'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[2]:
            l_res = []
            for RA_name in ['A', 'A3', 'A5', 'AN', 'DA', 'DA3', 'DA5', 'DAN', 'RA', 'RA3', 'RA5', 'RAN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'H2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'H8', 'N6', 'H61', 'H62']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'A'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[3]:
            l_res = []
            for RA_name in ['A', 'A3', 'A5', 'AN', 'DA', 'DA3', 'DA5', 'DAN', 'RA', 'RA3', 'RA5', 'RAN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'H2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9', 'H8', 'N6']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'A'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[4]:
            l_res = []
            for RA_name in ['C', 'C3', 'C5', 'CN', 'DC', 'DC3', 'DC5', 'DCN', 'RC', 'RC3', 'RC5', 'RCN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'H5', 'H6', 'O2', 'N4', 'H41', 'H42']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'C'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[5]:
            l_res = []
            for RA_name in ['C', 'C3', 'C5', 'CN', 'DC', 'DC3', 'DC5', 'DCN', 'RC', 'RC3', 'RC5', 'RCN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'H5', 'H6', 'N4']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'C'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[6]:
            l_res = []
            for RA_name in ['DT', 'DT3', 'DT5', 'DTN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'H5', 'H3', 'H6', 'O4', 'O2', 'C7', 'H71', 'H72', 'H73']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'T'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[7]:
            l_res = []
            for RA_name in ['DT', 'DT3', 'DT5', 'DTN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'H5', 'H6', 'C7', 'H71', 'H72', 'H73']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'T'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[8]:
            l_res = []
            for RA_name in ['U', 'U3', 'U5', 'UN', 'RU', 'RU3', 'RU5', 'RUN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'H5', 'H3', 'H6', 'O4', 'O2']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'U'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        if local_scheme[9]:
            l_res = []
            for RA_name in ['U', 'U3', 'U5', 'UN', 'RU', 'RU3', 'RU5', 'RUN']:
                l_res += self.get_all_res_name(res_list, RA_name)
            for j in l_res:
                atom_list = []
                for k in ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'H5', 'H6']:
                    atom_list.append(self.get_atom_id_from_name(k, int(j)))
                groups_res[i] = [j, 'U'] + atom_list
                atoms_used += atom_list
                i += 1
            l_res = []
        return (groups_res, atoms_used)


# class to implement a live feed in pymol - not yet tested
class PymolView():
    # initialises pymol and starts it without GUI, loads a given file and sets it to a white cartoon representation
    def __init__(self, file_name):
        print 'Launching Pymol'
        pymol.finish_launching()
        pymol.cmd.load(file_name)
        print '%s loaded in Pymol' % file_name

        # show molecule as cartoon and hide lines
        pymol.cmd.show('cartoon')
        pymol.cmd.hide('lines')
        pymol.cmd.set('cartoon_color', 0)

    def set_colour_range(self, start, finish, colour):
        # setting residue range to any colour
        pymol.cmd.set('cartoon_color', colour, 'resi %s-%s' % (str(start), str(finish)))

    # reset all of the molecule to a white colour
    # def reset_colour(self,length):
    #    pymol.cmd.set('cartoon_color',0,'resi 1-%s'%str(length))

    def set_colour_res(self, res, colour):
        pymol.cmd.set('cartoon_color', colour, 'resi %s' % str(res))

    def __del__(self):
        print 'Quit pymol'
        pymol.cmd.quit()

    def apply_colour_scheme(self, atomistic_list, local_list, colour_atomistic, colour_local):
        # self.reset_colour()
        for i in atomistic_list:
            self.set_colour_res(i, colour_atomistic)
        for j in local_list:
            self.set_colour_res(j, colour_local)


###FUNCTION TO PARSE THE PDB AND COORDS.INPCRD FILE ###################################################################
def reading_pdb_inpcrd(pdb, inpcrd):
    dic = {}
    try:
        atom_list = []
        with open(pdb, 'r') as f_pdb:
            for line in f_pdb:
                if line.split()[0] == 'ATOM':
                    atom_list.append((line[12:16].strip(), line[17:20].strip(), int(line[22:26])))
        with open(inpcrd, 'r') as f_crd:
            for _ in xrange(2):
                next(f_crd)
            i = 1
            for line in f_crd:
                atom_line = line.split()
                try:
                    atom_list[2 * i - 2] = atom_list[2 * i - 2] + (
                        float(atom_line[0]), float(atom_line[1]), float(atom_line[2]))
                except IndexError:
                    pass
                try:
                    atom_list[2 * i - 1] = atom_list[2 * i - 1] + (
                        float(atom_line[3]), float(atom_line[4]), float(atom_line[5]))
                    i = i + 1
                except IndexError:
                    pass
        for i in range(1, len(atom_list) + 1):
            dic[i] = atom_list[i - 1]
        return dic
    except IOError:
        return dic


def create_res_dic(atom_dic):
    res_dic = {}
    res_list = []
    for i in atom_dic.keys():
        res_list.append(atom_dic[i][2])
    res_list = list(set(res_list))
    for j in res_list:
        res_dic[j] = []
    for k in atom_dic.keys():
        res_dic[atom_dic[k][2]].append(k)
    return res_dic


###FUNCTIONS HELPING TO MANAGE THE LISTS AND DICTIONARY ###############################################################
def merge_lists(list_a, list_b):
    return sorted(list(set(list_a + list_b)))


def remove_selection_from_list(list_res, selection):
    new_res_list = []
    for i in list_res:
        if i not in selection:
            new_res_list.append(i)
    return sorted(new_res_list)


def check_res_exists(res_list, selection):
    new_selection = []
    for i in selection:
        if i in res_list:
            new_selection.append(i)
    return new_selection


def merge_new_and_old_res_lists(frozen, local, atomistic, sel_local, sel_atomistic):
    new_sel_atomistic = check_res_exists(frozen, sel_atomistic)
    new_sel_atomistic_from_local = check_res_exists(local, sel_atomistic)
    frozen = remove_selection_from_list(frozen, new_sel_atomistic)
    atomistic = merge_lists(atomistic, new_sel_atomistic)
    atomistic = merge_lists(atomistic, new_sel_atomistic_from_local)
    local = remove_selection_from_list(local, new_sel_atomistic_from_local)
    new_sel_local = check_res_exists(frozen, sel_local)
    local = merge_lists(local, new_sel_local)
    frozen = remove_selection_from_list(frozen, new_sel_local)
    return (sorted(atomistic), sorted(local), sorted(frozen))


def removal_selection(frozen, local, atomistic, sel_local, sel_atomistic):
    frozen = merge_lists(frozen, sel_local)
    frozen = merge_lists(frozen, sel_atomistic)
    local = remove_selection_from_list(local, sel_local)
    atomistic = remove_selection_from_list(atomistic, sel_atomistic)
    return (sorted(atomistic), sorted(local), sorted(frozen))


# FUNCTION TO SELECT RESIDUES ########################################################################################
def choose_res(mol):
    print '1 - choose sphere/ellipsoid of unfrozen atoms and add localised rigidification'
    print '2 - choose a range of residues or single residues to unfreeze or locally rigidify'
    print
    log = []
    try:
        method = int(raw_input('Method:  '))
    except IOError:
        method = 0
    unfrozen_res = []
    local_res = []

    if method == 0:
        print 'Invalid input'

    # ellipsoidal or spherical volume
    elif method == 1:
        while 1:  # choose the geometry: currently ellipsoidal 'e' or spherical 's'
            log.append('Method: sphere / ellipsoid')
            while 1:
                centre_method = raw_input('Mass weighted centre of residue (1) or atom in a selected residue (2)?  ')
                try:
                    if int(centre_method) == 1:
                        residue = raw_input('Residue:  ')
                        centre_xyz = mol.mass_weighted_centre(int(residue))
                        log.append('Mass-weighted centre: Residue' + residue + ' ' + mol.get_residue_name(
                            int(residue)) + '  ' + str(centre_xyz[0]) + '  ' + str(centre_xyz[1]) + '  ' + str(
                            centre_xyz[2]))
                        break
                    elif int(centre_method) == 2:
                        centre = mol.choose_centre()
                        centre_xyz = (mol.atom_dic[centre][3], mol.atom_dic[centre][4], mol.atom_dic[centre][5])
                        log.append('Centre: Atom  ' + str(centre) + '  ' + str(centre_xyz[0]) + '  ' + str(
                            centre_xyz[1]) + '  ' + str(centre_xyz[2]))
                        break
                    else:
                        print 'Invalid choice'
                except ValueError:
                    print 'Invalid input'

            radius_input = raw_input(
                'Radius of atomistic region (enter three values for an ellipsoid and one for a sphere):  ')
            l_radii = radius_input.split()
            try:
                if len(l_radii) == 1:
                    log.append('Radius of atomistic region: ' + l_radii[0])
                    l_atom = mol.get_atoms_in_sphere(centre_xyz, float(l_radii[0]))
                    l_res = mol.transform_atomlist_to_reslist(l_atom)
                    unfrozen_res += l_res
                    break
                elif len(l_radii) == 3:
                    log.append('Radii for the atomistic region: ' + l_radii[0] + '  ' + l_radii[1] + '  ' + l_radii[2])
                    l_atom = mol.get_atoms_in_ellipsoid(centre_xyz, float(l_radii[0]), float(l_radii[1]),
                                                        float(l_radii[2]))
                    l_res = mol.transform_atomlist_to_reslist(l_atom)
                    unfrozen_res += l_res
                    break
                else:
                    print 'Wrong input'
            except ValueError:
                print 'Invalid input'
        print
        while 1:  # geometry is same as for atomistic region, if a change is wanted employ the above
            check = raw_input('Add localised rigidification (y/n)?  ')
            if check == 'yes' or check == 'y':
                log.append('Local rigidification added')
                radius_input = raw_input(
                    'Radius of atomistic region (enter three values for an ellipsoid and one for a sphere):  ')
                l_radii_outer = radius_input.split()
                try:
                    if len(l_radii_outer) == 1:
                        log.append('Outer radius: ' + l_radii_outer[0])
                        l_atom_out = mol.get_atoms_in_sphere(centre_xyz, float(l_radii_outer[0]))
                        l_res_out = mol.transform_atomlist_to_reslist(l_atom_out)
                        for i in l_res_out:
                            if i not in l_res:
                                local_res.append(i)
                        break
                    elif len(l_radii_outer) == 3:
                        log.append(
                            'Radii for locally rigid region: ' + l_radii_outer[0] + '  ' + l_radii_outer[1] + '  ' +
                            l_radii_outer[2])
                        l_atom_out = mol.get_atoms_in_ellipsoid(centre_xyz, float(l_radii_outer[0]),
                                                                float(l_radii_outer[1]),
                                                                float(l_radii_outer[2]))
                        l_res_out = mol.transform_atomlist_to_reslist(l_atom_out)
                        for i in l_res_out:
                            if i not in l_res:
                                local_res.append(i)
                        break
                    else:
                        print 'Wrong input'
                except ValueError:
                    print 'Invalid input'
            else:
                break

    # rigidification using a range of residues, sanity check in function that chooses
    elif method == 2:
        log.append('Method:  range selection')
        check = raw_input('Define atomistic residues (y/n)?  ')
        if check == 'yes' or check == 'y':
            res = raw_input(
                'Enter residues separated by spaces, use r followed by two residues to define a range: ')  # input allows to enter 'r num1 num2' being the range between num1 and num2
            res = res.split()
            i = 0
            log_res = ''
            while i < len(res):
                if res[i] == 'r':  # if a range is indicated
                    i += 1
                    start_range = int(res[i])  # finds starting residue
                    i += 1
                    end_range = int(res[i])  # finds final residue
                    for j in range(start_range, end_range + 1):
                        unfrozen_res.append(j)
                        log_res += str(j) + ' '
                    i += 1
                else:
                    try:
                        unfrozen_res.append(int(res[i]))
                        log_res += res[i] + ' '
                        i += 1
                    except ValueError:
                        i += 1
            log.append('Residues entered as atomistically: ' + log_res)
        check = raw_input('Locally rigidify (y/n)?  ')
        if check == 'yes' or check == 'y':
            res = raw_input(
                'Enter residues separated by spaces, use r followed by two residues to define a range: ')  # input allows to enter 'r num1 num2' being the range between num1 and num2
            res = res.split()
            i = 0
            log_res = ''
            while i < len(res):
                if res[i] == 'r':  # if a range is indicated
                    i += 1
                    start_range = int(res[i])  # finds starting residue
                    i += 1
                    end_range = int(res[i])  # finds final residue
                    i += 1
                    for j in range(start_range, end_range + 1):
                        if j <= mol.num_res():
                            local_res.append(j)
                            log_res += str(j) + ' '
                        else:
                            break
                else:
                    try:
                        if int(res[i]) <= mol.num_res():
                            local_res.append(int(res[i]))
                            log_res += res[i] + ' '
                            i += 1
                    except ValueError:
                        i += 1
            log.append('Residues entered for local rigidification: ' + log_res)

    else:
        return 'Wrong input'
    return (unfrozen_res, local_res, log)


def remove_res(mol):
    unfrozen_res = []
    local_res = []
    log = []
    log.append('Removal of selected residues')
    check = raw_input('Remove from the atomistic residues (y/n)?  ')
    if check == 'yes' or check == 'y':
        res = raw_input(
            'Enter residues separated by spaces, use r followed by two residues to define a range: ')  # input allows to enter 'r num1 num2' being the range between num1 and num2
        res = res.split()
        i = 0
        log_res = ''
        while i < len(res):
            if res[i] == 'r':  # if a range is indicated
                i += 1
                start_range = int(res[i])  # finds starting residue
                i += 1
                end_range = int(res[i])  # finds final residue
                for j in range(start_range, end_range + 1):
                    unfrozen_res.append(j)
                    log_res += str(j) + ' '
                i += 1
            else:
                try:
                    unfrozen_res.append(int(res[i]))
                    log_res += res[i] + ' '
                    i += 1
                except ValueError:
                    i += 1
        log.append('Residues entered to be removed: ' + log_res + "\n")
    check = raw_input('Remove from the local selection (y/n)?  ')
    if check == 'yes' or check == 'y':
        res = raw_input(
            'Enter residues separated by spaces, use r followed by two residues to define a range: ')  # input allows to enter 'r num1 num2' being the range between num1 and num2
        res = res.split()
        i = 0
        log_res = ''
        while i < len(res):
            if res[i] == 'r':  # if a range is indicated
                i += 1
                start_range = int(res[i])  # finds starting residue
                i += 1
                end_range = int(res[i])  # finds final residue
                i += 1
                for j in range(start_range, end_range + 1):
                    if j <= mol.num_res():
                        local_res.append(j)
                        log_res += str(j) + ' '
                    else:
                        break
            else:
                try:
                    if int(res[i]) <= mol.num_res():
                        local_res.append(int(res[i]))
                        log_res += res[i] + ' '
                        i += 1
                except ValueError:
                    i += 1
        log.append('Residues entered to be removed from local rigidification: ' + log_res + "\n")
    return (unfrozen_res, local_res, log)


# FUNCTIONS FOR GROUPROTATION FILE ###################################################################################
def get_prob_amb():
    while 1:
        inp_p = 'Probability: '
        inp_a = 'Amplitude: '
        try:
            probability = float(inp_p)
            amplitude = float(inp_a)
            if (probability > 1.0) or (probability < 0.0):
                print 'Invalid probability'
            else:
                break
        except ValueError:
            print 'Both values need to be numbers'
    return (probability, amplitude)


# VARIABLE DEFINITION ################################################################################################
atom_dic = {}  # dictionary containing all information sorted by atom number:
# the information is in a tuple: (atomname, residue name, residue number, x, y, z)
res_dic = {}  # dictionary containing all atoms per residue sorted by residue number
res_atomistic = []  # residue list that are treated fully atomistic
res_local = []  # residue list that are treated locally rigid
res_frozen = []  # residue list that are treated fully rigid
atoms_used = []  # atoms used in rigid bodies, prevents double usage!
groups = {}  # dictionary of rigid groups
groups_local = {}
num_rbody = 0  # number of rigid bodies defined

# PARSING OF INPUT FILES #############################################################################################
while 1:  # starts input selection
    atom_dic = reading_pdb_inpcrd(pdb_inp, coords_inp)  # opens files and loads them into a dictionary
    res_dic = create_res_dic(atom_dic)
    if atom_dic != {}:  # if the dictionary contains data the data is read into the molecule class
        loaded_mol = protein(atom_dic, res_dic)
        print
        print 'Molecule was loaded'
        print 'Number of residues: ', loaded_mol.num_res()
        print 'Number of atoms: ', loaded_mol.num_atoms()
        print
        break
    else:
        print 'Files could not be found or they are empty - please try again'
        print
        pdb_inp = raw_input('PDB file: ')
        coords_inp = raw_input('Coords input: ')

check_file = open('rigid.log', 'w')
check_file.write('Input file for structure: ' + pdb_inp + "\n")
check_file.write('Coordinates file used: ' + coords_inp + "\n" + "\n")
check_file.write('Number of residues:  ' + str(loaded_mol.num_res()) + "\n")
check_file.write('Number of atoms:  ' + str(loaded_mol.num_atoms()) + "\n" + "\n")
check_file.write('Time:  ' + time.strftime("%a, %d %b %Y %H:%M:%S") + "\n" + "\n")

# CHECK FOR LIGANDS AND UNKNOWN RESIDUES #############################################################################
check_file.write('Ligands and unknown residues:' + "\n")
ligand_dic = loaded_mol.get_ligands_unknown_res()
if ligand_dic == {}:
    print 'There are no ligands or unknown residues in the molecule'
    check_file.write('None' + "\n" + "\n")
else:
    print 'There are ligands or unknown residues.'
    ligand_list = ligand_dic.keys()
    for res_id in ligand_dic.keys():
        print str(res_id) + ' - ' + loaded_mol.res_name_from_id(res_id)
        check_file.write(str(res_id) + ' - ' + loaded_mol.res_name_from_id(res_id) + "\n")
    check_ligand = raw_input('Shall some or all of the above be treated as normal residues, all other residues will'
                             ' be treated as fully atomistic? (y/n) ')
    if check_ligand == ('y' or 'yes'):
        while 1:
            ligand_inp = raw_input('Enter the residue numbers separated by spaces (enter a for all residues): ')
            if ligand_inp == 'a':
                print 'All residues/ligands are treated normally.'
                check_file.write('All are treated normally.' + "\n" + "\n")
                break
            else:
                try:
                    res_selected = [int(s) for s in ligand_inp.split()]
                    new_ligand_list = []
                    string_ligands = ''
                    for i in ligand_list:
                        if i not in res_selected:
                            new_ligand_list.append(i)
                            string_ligands += str(i) + ' '
                    res_atomistic += new_ligand_list
                    print 'The following residues will be treated atomistically: ' + string_ligands
                    check_file.write('The following residues will be treated atomistically: ' +
                                     string_ligands + "\n" + "\n")
                    break
                except ValueError:
                    print 'Please enter integer numbers.'
    else:
        res_atomistic += ligand_list
        print 'The above listed residues/ligands will be treated fully atomistically.'
        check_file.write('All are treated fully atomistically.' + "\n" + "\n")

# SELECTION OF RESIDUES FOR FULLY AND LOCALLY RIGID REGIONS ##########################################################
print
print
print 'Selection of atomistic and locally rigidified residues'
print
check_file.write('Selection of atomistic and locally rigidified residues' + "\n")
res_frozen = remove_selection_from_list(res_dic.keys(), res_atomistic)
x = choose_res(loaded_mol)  # initiates the first round of selections
y = merge_new_and_old_res_lists(res_frozen, res_local, res_atomistic, x[1], x[0])
res_frozen = y[2]
res_local = y[1]
res_atomistic = y[0]
log_string = ''
for strings in x[2]:
    log_string += strings + "\n"
check_file.write(log_string + "\n" + "\n")

while 1:  # allows to choose more residues with different methods
    print
    print 'The following residues are not frozen (atomistic):'
    print res_atomistic
    print
    print 'The following residues are selected to be locally rigidified:'
    print res_local
    print
    print 'The following residues are selected to be fully rigidified:'
    print res_frozen
    print
    check = raw_input('Enter more residues (1) or remove residues from previous selection (2)? Any other key to exit. ')
    try:
        value = int(
            check)  # to enter more residues the choose_res function is called again with the same options as before
        if int(check) == 1:
            x = choose_res(loaded_mol)
            y = merge_new_and_old_res_lists(res_frozen, res_local, res_atomistic, x[1], x[0])
            res_frozen = y[2]
            res_local = y[1]
            res_atomistic = y[0]
            log_string = ''
            for strings in x[2]:
                log_string += strings + "\n"
            check_file.write(log_string + "\n" + "\n")
        elif int(check) == 2:
            x = remove_res(loaded_mol)
            y = removal_selection(res_frozen, res_local, res_atomistic, x[1],
                                  x[0])  # the user enters a selection for removal
            res_frozen = y[2]
            res_local = y[1]
            res_atomistic = y[0]
            check_file.write('Remove residues' + "\n" + 'Atomistic residues removed:' + "\n")
            string = ''
            for i in x[0]:
                string = string + str(i) + '  '
            check_file.write(string + "\n" + "\n")
            check_file.write('Locally rigidified residues removed:' + "\n")
            string = ''
            for i in x[1]:
                string = string + str(i) + '  '
            check_file.write(string + "\n" + "\n")
        else:
            break

    except ValueError:
        break
check_file.write('Selection completed' + "\n" + "\n")
string_atomistic = ''
for res in res_atomistic:
    string_atomistic += str(res) + ', '
string_local = ''
for res in res_local:
    string_local += str(res) + ', '
check_file.write('Atomistic residues: ' + string_atomistic + "\n")
check_file.write('Locally rigidified residues: ' + string_local + "\n")

# GROUPING OF THE RIGID BODIES ######################################################################################
if res_frozen != []:
    for i in res_frozen:
        atoms_used += res_dic[i]
    print
    print 'Grouping of the rigid body'
    print
    print '1-Group all parts of the rigid system as one body'
    print '2-Define groups'
    print
    while 1:
        check = raw_input('Grouping method:  ')
        try:
            value = int(check)
            if int(check) == 1:  # all residues are grouped as one rigid body
                num_rbody += 1
                check_file.write('All frozen residues are rigidified as one body' + "\n" + "\n" + "\n")
                groups = {1: res_frozen}
                break
            elif int(check) == 2:  # residues are frozen in a number of rigid bodies
                check_file.write('The frozen residues are rigidified as more than one body' + "\n" + "\n" + "\n")
                print 'The following residues are frozen:'
                print res_frozen
                print
                res_frozen_left = res_frozen
                num_groups = int(raw_input('How many groups will be defined?  '))  # number of rigid bodies is defined
                counter = 0
                while counter < num_groups:
                    counter += 1
                    num_rbody += 1
                    group = raw_input('Residues for rigid body: ')
                    # input allows to enter 'r num1 num2' being the range between num1 and num2
                    group = group.split()
                    i = 0
                    group_res = []
                    while i < len(group):
                        if group[i] == 'r':  # if a range is indicated
                            i += 1
                            start_range = int(group[i])  # finds starting residue
                            i += 1
                            end_range = int(group[i])  # finds final residue
                            for j in range(start_range, end_range):
                                if j in res_frozen_left:  # appends list if residues have not be assigned earlier
                                    group_res.append(j)
                        else:
                            try:
                                if int(group[i]) in res_frozen_left:
                                    group_res.append(int(group[i]))  # appends single residues entered
                                i += 1
                            except ValueError:
                                i += 1
                    res_frozen_left = remove_selection_from_list(res_frozen_left,
                                                                 group_res)
                    # removes newly assigned residues from rest of residues
                    groups[num_rbody] = group_res  # defines a new entry in the dictionary that saves the groups
                    if res_frozen_left == []:
                        # automatically terminates when no residues are left to be assigned even if the number of rigid
                        # bodies entered is not reached
                        print 'All residues are assigned.'
                        break
                    if counter == num_groups and res_frozen_left != []:
                        # if all rigid bodies are assigned and there are still residues not assigned
                        # it deletes all groups and restarts the process
                        print 'Not all residues were assigned. Please start again.'
                        counter = 0
                        groups = {}
                        res_frozen_left = res_frozen
                break
            else:
                pass
        except ValueError:
            print 'Wrong input'
    print
    print 'The following rigid bodies have been specified:'  # prints all the entered rigid bodies from above
    for i in groups.keys():
        print 'Rigid body: ', i
        print groups[i]
        print

else:
    groups = {}

# LOCALLY RIGIDIFICATION, STANDARD GROUPS ###########################################################################
NA_res = ['RAN', 'RCN', 'RGN', 'RUN', 'A', 'A3', 'A5', 'AN', 'C', 'C3', 'C5', 'CN', 'G', 'G3', 'G5', 'GN', 'U',
          'U3', 'U5', 'UN', 'OHE', 'RA', 'RA3', 'RA5', 'RC', 'RC3', 'RC5', 'RG', 'RG3', 'RG5', 'RU', 'RU3',
          'RU5', 'DA', 'DA5', 'DC', 'DC5', 'DG', 'DG5', 'DT', 'DT5', 'DA3', 'DC3', 'DG3', 'DT3', 'DAN', 'DCN',
          'DGN', 'DTN']
protein_res = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'CYX', 'GLN', 'GLU', 'GLY', 'GLX', 'HIS', 'HIE',
               'HID', 'HIP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
               'NALA', 'NARG', 'NASN', 'NASP', 'NASX', 'NCYS', 'NCYX', 'NGLN', 'NGLU', 'NGLY', 'NGLX',
               'NHIS', 'NHIE', 'NHID', 'NHIP', 'NILE', 'NLEU', 'NLYS', 'NMET', 'NPHE', 'NPRO', 'NSER',
               'NTHR', 'NTRP', 'NTYR', 'NVAL', 'CALA', 'CARG', 'CASN', 'CASP', 'CASX', 'CCYS', 'CCYX',
               'CGLN', 'CGLU', 'CGLY', 'CGLX', 'CHIS', 'CHIE', 'CHID', 'CHIP', 'CILE', 'CLEU', 'CLYS',
               'CMET', 'CPHE', 'CPRO', 'CSER', 'CTHR', 'CTRP', 'CTYR', 'CVAL']

if res_local != []:
    print 'Grouping of the local rigid bodies'
    print
    local_scheme = []
    if protein_t:
        proline_check = raw_input('Group proline as rigid body C-O-N-CD-CG-CB-CA(1) or as  N-CD-CG-CB-CA(2)? ')
        if proline_check == '1':
            check_file.write('Grouped proline as rigid body C-O-N-CD-CG-CB-CA' + "\n")
            local_scheme.append(1)
        elif proline_check == '2':
            check_file.write('Grouped proline as rigid body  N-CD-CG-CB-CA' + "\n")
            local_scheme.append(2)
        else:
            local_scheme.append(0)

        arginine_check = raw_input('Group arginine as rigid body NE-HE-CZ-NH1-HH11-HH12-NH2-HH21-HH22 (y/n)? ')
        if arginine_check == 'y':
            check_file.write('Grouped arginine as rigid body NE-HE-CZ-NH1-HH11-HH12-NH2-HH21-HH22' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        histidine_check = raw_input('Group histidine (HIS, HIE, HID, HIP) as rigid body CG-ND1-CE1-NE2-CD2 (y/n)? ')
        if histidine_check == 'y':
            check_file.write('Grouped histidine (HIS, HIE, HID, HIP) as rigid body CG-ND1-CE1-NE2-CD2' + "\n")
            local_scheme.append(1)
           # local_scheme.append(1)
           # local_scheme.append(1)
           # local_scheme.append(1)
        else:
            local_scheme.append(0)
           # local_scheme.append(0)
           # local_scheme.append(0)
           # local_scheme.append(0)

        lysine_check = raw_input('Group lysine as rigid body NZ-HZ1-HZ2-HZ3 (y/n)? ')
        if lysine_check == 'y':
            check_file.write('Grouped lysine as rigid body NZ-HZ1-HZ2-HZ3' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        aspartic_check = raw_input('Group aspartic acid as rigid body CB-CG-OD1-OD2 (y/n)? ')
        if aspartic_check == 'y':
            check_file.write('Grouped aspartic acid as rigid body CB-CG-OD1-OD2' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        asparagine_check = raw_input('Group asparagine as rigid body CB-CG-OD1-ND2-HD21-HD22 (y/n)? ')
        if asparagine_check == 'y':
            check_file.write('Grouped asparagine as rigid body GB-CG-OD1-ND2-HD21-HD22' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        glutamic_check = raw_input('Group glutamic acid as rigid body CG-CD-OE1-OE2 (y/n)? ')
        if glutamic_check == 'y':
            check_file.write('Grouped glutamic acid as rigid body CG-CD-OE1-OE2' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        glutamine_check = raw_input('Group glutamine as rigid body CG-CD-OE1-NE2-HE21-HE22 (y/n)? ')
        if glutamine_check == 'y':
            check_file.write('Grouped glutamine as rigid body CG-CD-OE1-NE2-HE21-HE22' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        phenylalanine_check = raw_input(
            'Group phenylalanine as rigid body CG-CD1-HD1-CE1-HE1-CZ-HZ-CE2-HE2-CD2-HD2 (y/n)? ')
        if phenylalanine_check == 'y':
            check_file.write('Grouped phenylalanine as rigid body CG-CD1-HD1-CE1-HE1-CZ-HZ-CE2-HE2-CD2-HD2' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        tyrosine_check = raw_input('Group tyrosine as rigid body CG-CD1-HD1-CE1-HE1-CZ-OH-CE2-HE2-CD2-HD2 (y/n)? ')
        if tyrosine_check == 'y':
            check_file.write('Grouped tyrosine as rigid body CG-CD1-HD1-CE1-HE1-CZ-OH-CE2-HE2-CD2-HD2' + "\n")
            local_scheme.append(1)
        else:
            local_scheme.append(0)

        tryptophan_check = raw_input(
            'Group tryptophan as either: 1 rigid body:(CG-CD1-HD1-NE1-HE1-CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-CE3-HE3-CD2) '
            'or 2 rigid bodies:(CG-CD1-HD1-NE1-HE1)-(CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-CE3-HE3-CD2)? ')
        if tryptophan_check == '1':
            check_file.write(
                'Grouped tryptophan as rigid body CG-CD1-HD1-NE1-HE1-CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-HZ3-CE3-HE3-CD2' +
                "\n" + "\n")
            local_scheme.append(1)
        elif tryptophan_check == '2':
            check_file.write(
                'Grouped tryptophan as two rigid bodies CG-CD1-HD1-NE1-HE1 and' +
                ' CE2-CZ2-HZ2-CH2-HH2-CZ3-HZ3-HZ3-CE3-HE3-CD2' + "\n" + "\n")
            local_scheme.append(2)
        else:
            local_scheme.append(0)

        loc_protein = loaded_mol.get_rigid_for_local(res_local, local_scheme, atoms_used)
        groups_local = loc_protein[0]
        atoms_used = loc_protein[1]
    print

    if protein_t:
        check_peptide_bonds = raw_input('Shall the peptide bonds be treated as rigid bodies (excl. PRO) (y/n)? ')
        if check_peptide_bonds == 'y' or 'yes':
            pairlist = []  # all peptide bonds are part of two residues
            for i in res_local:
                if (i - 1, i) not in pairlist and (i - 1 > 0):
                    if loaded_mol.get_residue_name(i) != 'PRO':
                        if (loaded_mol.get_residue_name(i) in protein_res) and \
                                (loaded_mol.get_residue_name(i - 1) in protein_res):
                            pairlist.append((i - 1, i))
                if (i, i + 1) not in pairlist and (i + 1 <= loaded_mol.num_res()):
                    if loaded_mol.get_residue_name(i + 1) != 'PRO':
                        if (loaded_mol.get_residue_name(i) in protein_res) and \
                                (loaded_mol.get_residue_name(i + 1) in protein_res):
                            pairlist.append((i, i + 1))
            for j in pairlist:
                atomlist = [loaded_mol.get_atom_id_from_name('C', j[0]), loaded_mol.get_atom_id_from_name('O', j[0]),
                            loaded_mol.get_atom_id_from_name('N', j[1]), loaded_mol.get_atom_id_from_name('H', j[1])]
                check_double = 0
                for k in atomlist:
                    if k in atoms_used:
                        check_double = 1
                if check_double == 0:
                    atoms_used += atomlist
                    check_file.write(
                        'Peptide bond rigidified for residues ' + str(j[0]) + ' and ' + str(j[1]) + "\n")
                    groups_local[len(groups_local) + 1] = [j, (loaded_mol.get_residue_name(j[0]),
                                                               loaded_mol.get_residue_name(j[1]))] + atomlist

    if rna_t or dna_t:
        phosphate_check = raw_input("Group phosphate as (O3'-POO-O5') (1) or POO (2)? ")
        if phosphate_check == '1':
            NArigid_PO4_large = True
            NArigid_PO4_small = False
            check_file.write("Grouped phosphate group between nucleic acids as (O3'-POO-O5')" + "\n")
        elif phosphate_check == '2':
            NArigid_PO4_large = False
            NArigid_PO4_small = True
            check_file.write("Grouped phosphate group between nucleic acids as (-POO-)" + "\n")
        else:
            NArigid_PO4_large = False
            NArigid_PO4_small = False

        sugar_check = raw_input("Rigidify sugar (y/n)? ")
        if sugar_check == 'y':
            NArigid_sugar = True
            check_file.write('Grouped sugar group in nucleic acids as rigid body' + "\n")
        else:
            NArigid_sugar = False

        guanine_check = raw_input('Group guanine base as one rigid body (1) or leave contacts free (2)? ')
        if guanine_check == '1':
            NArigid_G_1 = True
            NArigid_G_2 = False
        elif guanine_check == '2':
            NArigid_G_1 = False
            NArigid_G_2 = True
        else:
            NArigid_G_1 = False
            NArigid_G_2 = False

        adenine_check = raw_input('Group adenine base as one rigid body (1) or leave contacts free (2)? ')
        if adenine_check == '1':
            NArigid_A_1 = True
            NArigid_A_2 = False
        elif adenine_check == '2':
            NArigid_A_1 = False
            NArigid_A_2 = True
        else:
            NArigid_A_1 = False
            NArigid_A_2 = False

        cystosine_check = raw_input('Group cystosine base as one rigid body (1) or leave contacts free (2)? ')
        if cystosine_check == '1':
            NArigid_C_1 = True
            NArigid_C_2 = False
        elif cystosine_check == '2':
            NArigid_C_1 = False
            NArigid_C_2 = True
        else:
            NArigid_C_1 = False
            NArigid_C_2 = False

        NArigid_T_1 = False
        NArigid_T_2 = False
        NArigid_U_1 = False
        NArigid_U_2 = False
        if dna_t:
            thymine_check = raw_input('Group thymine base as one rigid body (1) or leave contacts free (2)? ')
            if thymine_check == '1':
                NArigid_T_1 = True
                NArigid_T_2 = False
            elif thymine_check == '2':
                NArigid_T_1 = False
                NArigid_T_2 = True
            else:
                NArigid_T_1 = False
                NArigid_T_2 = False

        if rna_t:
            uracil_check = raw_input('Group uracil base as one rigid body (1) or leave contacts free (2)? ')
            if uracil_check == '1':
                NArigid_U_1 = True
                NArigid_U_2 = False
            elif uracil_check == '2':
                NArigid_U_1 = False
                NArigid_U_2 = True
            else:
                NArigid_U_1 = False
                NArigid_U_2 = False

        gloc_NA, atoms_used = loaded_mol.loc_NA(res_local,
                                                [NArigid_G_1, NArigid_G_2, NArigid_A_1, NArigid_A_2, NArigid_C_1,
                                                 NArigid_C_2, NArigid_T_1, NArigid_T_2, NArigid_U_1, NArigid_U_2],
                                                atoms_used)

        print gloc_NA

        for key in gloc_NA.keys():
            groups_local[len(groups_local) + 1] = gloc_NA[key]  # add NA groups to overall groups

        # now do all sugars and phosphate if necessary



        if NArigid_PO4_large:
            pairlist = []  # all phosphates that are part of two NAs
            for i in res_local:
                if loaded_mol.get_residue_name(i) in NA_res:
                    if (i - 1, i) not in pairlist and (i - 1 > 0):
                        pairlist.append((i - 1, i))
                    if (i, i + 1) not in pairlist and (i + 1 <= loaded_mol.num_res()):
                        pairlist.append((i, i + 1))
            for j in pairlist:
                atomlist = [loaded_mol.get_atom_id_from_name("O3'", j[0]), loaded_mol.get_atom_id_from_name('P', j[1]),
                            loaded_mol.get_atom_id_from_name('OP1', j[1]),
                            loaded_mol.get_atom_id_from_name('OP2', j[1]),
                            loaded_mol.get_atom_id_from_name("O5'", j[1])]

                atoms_used += atomlist
                check_file.write(
                    'Phosphate (all) rigidified for residues ' + str(j[0]) + ' and ' + str(j[1]) + "\n")
                groups_local[len(groups_local) + 1] = [j, (loaded_mol.get_residue_name(j[0]),
                                                           loaded_mol.get_residue_name(j[1]))] + atomlist
        if NArigid_PO4_small:
            for i in res_local:
                if loaded_mol.get_residue_name(i) in NA_res:
                    atomlist = [loaded_mol.get_atom_id_from_name('P', i), loaded_mol.get_atom_id_from_name('OP1', i),
                                loaded_mol.get_atom_id_from_name('OP2', i)]
                    atoms_used += atomlist
                    check_file.write(
                        'Phosphate (only POO) rigidified for residue ' + str(i) + "\n")
                    groups_local[len(groups_local) + 1] = [i, loaded_mol.get_residue_name(i)] + atomlist

        if NArigid_sugar:
            for i in res_local:
                if loaded_mol.get_residue_name(i) in NA_res:
                    atomlist = [loaded_mol.get_atom_id_from_name("C5'", i), loaded_mol.get_atom_id_from_name("H5'", i),
                                loaded_mol.get_atom_id_from_name("H5''", i), loaded_mol.get_atom_id_from_name("C4'", i),
                                loaded_mol.get_atom_id_from_name("H4'", i), loaded_mol.get_atom_id_from_name("O4'", i),
                                loaded_mol.get_atom_id_from_name("C3'", i), loaded_mol.get_atom_id_from_name("H3'", i),
                                loaded_mol.get_atom_id_from_name("C2'", i), loaded_mol.get_atom_id_from_name("H2'", i),
                                loaded_mol.get_atom_id_from_name("C1'", i), loaded_mol.get_atom_id_from_name("H1'", i)]
                    if loaded_mol.get_residue_name(i) in ['DA', 'DA5', 'DC', 'DC5', 'DG', 'DG5', 'DT', 'DT5', 'DA3',
                                                          'DC3', 'DG3', 'DT3', 'DAN', 'DCN','DGN', 'DTN']:
                        atomlist.append(loaded_mol.get_atom_id_from_name("H2''", i))
                    else:
                        atomlist.append(loaded_mol.get_atom_id_from_name("O2'", i))
                        atomlist.append(loaded_mol.get_atom_id_from_name("HO2'", i))
                    
                    atoms_used += atomlist
                    check_file.write(
                        'Sugar rigidified for residue ' + str(i) + "\n")
                    groups_local[len(groups_local) + 1] = [i, loaded_mol.get_residue_name(i)] + atomlist

    check_file.write("\n")
    check_new_groups = raw_input('Enter additional user-defined rigid bodies (y/n)? ')
    if check_new_groups in ['yes', 'y']:
        print '1 - select a residues by name'
        print '2 - select one or more residues by id'
        while 1:
            while 1:
                local_method = raw_input('Method (1/2): ')
                try:
                    if int(local_method) in [1, 2]:
                        break
                    else:
                        print 'Invalid choice'
                except ValueError:
                    print 'Invalid input'

            if int(local_method) == 1:
                while 1:
                    exit_1 = 0
                    res_name = raw_input('Residue name (three letter code, n to exit): ')
                    if res_name == 'n':
                        exit_1 = 1
                        break
                    try:
                        res_list = loaded_mol.get_all_res_name(res_local, res_name)
                        if res_list != []:
                            check_file.write('Additional rigid bodies for residues: ' + res_name + "\n")
                            break
                        else:
                            print 'None of the residues for local rigidification is of that type.'
                    except ValueError:
                        print 'Invalid input'
                if exit_1 == 0:
                    all_atoms_res_name = []
                    atoms_used_name = []
                    for i in res_list:
                        all_atoms_res_name += res_dic[i]
                    for j in all_atoms_res_name:
                        if j in atoms_used:
                            atoms_used_name.append(j)
                    if atoms_used_name != []:
                        print 'The following atoms are already in local rigid bodies: '
                        for k in atoms_used_name:
                            print str(k) + ' ' + atom_dic[k][0] + ' ' + atom_dic[k][1] + ' ' + str(atom_dic[k][2])
                    print
                    print 'The residues %s contain the following atoms:' % res_name
                    atoms_in_res = loaded_mol.get_atom_list(res_list[0])
                    for l in atoms_in_res:
                        print atom_dic[l][0]
                    print
                    atom_seq = raw_input('Enter a sequence of atom names separated by spaces: ')
                    for m in res_list:
                        check_res = 0
                        atomlist = res_dic[m]
                        atoms_rbody = []
                        for n in atomlist:
                            if atom_dic[n][0] in atom_seq.split():
                                atoms_rbody.append(n)
                        if len(atoms_rbody) != len(atom_seq.split()):
                            print 'Not all atoms entered are in the residue.'
                            check_file.write('For residue: ' + str(m) + ' - not all entered atoms in residue' + "\n")
                            check_res = 1
                        if len(atoms_rbody) < 3:
                            print 'Rigid bodies need to have at least 3 atoms.'
                            check_file.write('For residue: ' + str(m) + ' - not 3 or more atoms in residue' + "\n")
                            check_res = 1
                        for o in atoms_rbody:
                            if o in atoms_used_name:
                                print 'Atoms cannot belong to two rigid bodies.'
                                check_file.write('For residue: ' + str(m) +
                                                 ' - atoms cannot belong to two rigid bodies' + "\n" + "\n")
                                check_res = 1
                                break
                        if loaded_mol.check_linear(atoms_rbody, tolerance):
                            print 'Linear arrangement of atoms cannot be a rigid body.'
                            check_file.write('For residue: ' + str(m) + ' - linear group entered' + "\n")
                            check_res = 1
                        if check_res == 0:
                            groups_local[len(groups_local) + 1] = [m, loaded_mol.get_residue_name(m)] + atoms_rbody
                            atoms_used += atoms_rbody
                            check_file.write('For residue: ' + str(m) + ' new group: ' + atom_seq + "\n")
                    check_file.write("\n")
            elif int(local_method) == 2:
                while 1:
                    exit_2 = 0
                    res_inp = raw_input('Residue numbers (n to exit): ')
                    if res_inp == 'n':
                        exit_2 = 1
                        break
                    try:
                        res_list_str = res_inp.split()
                        res_list_int = []
                        for res in res_list_str:
                            res_list_int.append(int(res))
                        if res_list_int != []:
                            check_file.write('Additional rigid bodies for residues: ' + res_inp + "\n")
                            break
                        else:
                            print 'No residues entered'
                    except ValueError:
                        print 'Invalid input'
                if exit_2 == 0:
                    all_atoms_res = []
                    atoms_used_res = []
                    for i in res_list_int:
                        all_atoms_res += res_dic[i]
                    for j in all_atoms_res:
                        if j in atoms_used:
                            atoms_used_res.append(j)
                    atoms_allowed = []
                    for atom in all_atoms_res:
                        if atom not in atoms_used_res:
                            atoms_allowed.append(atom)
                    if atoms_used_res != []:
                        print 'The following atoms are already in local rigid bodies: '
                        for k in atoms_used_res:
                            print str(k) + ' ' + atom_dic[k][0] + ' ' + atom_dic[k][1] + ' ' + str(atom_dic[k][2])
                    print
                    print 'The residues %s contain the following atoms that are not assigned to rigid bodies:' \
                          % res_list_str
                    string_atoms = ''
                    for atom in atoms_allowed:
                        string_atoms += str(atom) + ' - ' + atom_dic[atom][0] + ', '
                    print string_atoms
                    print
                    atom_seq = raw_input('Enter a sequence of atom numbers separated by spaces: ')
                    check_res = 0
                    atoms_rbody = []
                    check_rbody = 0
                    for atom_str in atom_seq.split():
                        try:
                            if int(atom_str) in atoms_allowed:
                                atoms_rbody.append(int(atom_str))
                            elif int(atom_str) in atoms_used_res:
                                print 'Atoms cannot belong to two rigid bodies.'
                                check_file.write('Atoms cannot belong to two rigid bodies' + "\n")
                                check_rbody = 1
                                break
                        except ValueError:
                            pass
                    print len(atoms_rbody)
                    print len(atom_seq.split())
                    if len(atoms_rbody) != len(atom_seq.split()):
                        print 'Not all atoms entered are in the residue.'
                        check_file.write('Not all entered atoms in residue' + "\n")
                        check_rbody = 1
                    if len(atoms_rbody) < 3:
                        print 'Rigid bodies need to have at least 3 atoms.'
                        check_file.write('Not 3 or more atoms in residue' + "\n")
                        check_rbody = 1
                    if loaded_mol.check_linear(atoms_rbody, tolerance):
                        print 'Linear arrangement of atoms cannot be a rigid body.'
                        check_file.write('Linear group entered' + "\n")
                        check_rbody = 1
                    if check_rbody == 0:
                        groups_local[len(groups_local) + 1] = ['', ''] + atoms_rbody
                        atoms_used += atoms_rbody
                        check_file.write('For residues: ' + res_inp + ' new group: ' + atom_seq + "\n")
                    check_file.write("\n")
            check_more = raw_input('Enter more rigid bodies (y/n)? ')
            if check_more not in ['y', 'yes']:
                break

# WRITING OUTPUT FILES ##############################################################################################

print 'Output files will be written now.'

# COORDSINIRIGID #
coords_out = open('coordsinirigid', 'w')

for i in range(1, len(atom_dic.keys()) + 1):
    atom = str(atom_dic[i][3]) + '   ' + str(atom_dic[i][4]) + '   ' + str(atom_dic[i][5])
    coords_out.write(atom)
    if i < len(atom_dic.keys()):
        coords_out.write("\n")
coords_out.close()

# RBODYCONFIG FILE #
groups_file = open('rbodyconfig', 'w')
natom_rbody = 0
for group in range(1, len(groups.keys()) + 1):
    string = ''
    counter = 0
    for res in groups[group]:
        for atom in res_dic[res]:
            string += str(atom) + "\n"
            counter += 1
            natom_rbody += 1
    groups_file.write('GROUP ' + str(counter) + "\n" + string)
for group_local in range(1, len(groups_local.keys()) + 1):
    string = ''
    counter = 0
    for atom in groups_local[group_local][2:]:
        string += str(atom) + "\n"
        counter += 1
        natom_rbody += 1
    groups_file.write('GROUP ' + str(counter) + "\n" + string)
groups_file.close()

# FINISH LOG FILE FOR RIGID BODIES #
check_file.write("\n" + "\n" + 'Number of rigid bodies: ' + str(len(groups.keys()) + len(groups_local.keys())))
check_file.write(
    "\n" + 'Degrees of freedom: ' + str(3 *(len(atom_dic) - natom_rbody) + 6 * (len(groups.keys()) + len(groups_local.keys()))))

# CREATING GROUP ROTATION FILES #####################################################################################
print ''
group_rot_t = raw_input('Create group rotation file (y/n)? ')
if group_rot_t in ['y', 'yes']:
    check_file.write("\n" + "\n" + 'Create group rotation file' + "\n")
    while 1:
        prob_uni_inp = raw_input('Uniform selection probability (suggested: 0.025): ')
        try:
            prob_uni = float(prob_uni_inp)
            if prob_uni > 0.1:
                print 'This is a large probability (values recommended: 0.01 to 0.1).'
                prob_uni_check = raw_input(' Continue with the selected value (y/n)? ')
                if prob_uni_check in ['y', 'yes']:
                    break
            elif prob_uni < 0.01:
                print 'This is a small probability (values recommended: 0.01 to 0.1).'
                prob_uni_check = raw_input(' Continue with the selected value (y/n)? ')
                if prob_uni_check in ['y', 'yes']:
                    break
            else:
                break
        except ValueError:
            print 'Input needs to be a float'
    check_file.write('Selection probability: ' + str(prob_uni) + "\n")
    while 1:
        grouprot_reg_inp = raw_input(
            'Create the group rotation file for the all atom region (1), the locally rigidified region (2),' +
            ' or both regions (3)? ')
        try:
            group_reg = int(grouprot_reg_inp)
            if group_reg == 1:
                grot_res = res_atomistic
                check_file.write('Only all atom region' + "\n")
                break
            elif group_reg == 2:
                grot_res = res_local
                check_file.write('Only local rigid region' + "\n")
                break
            elif group_reg == 3:
                grot_res = res_atomistic + res_local
                check_file.write('Local rigid and all atom regions' + "\n")
                break
            else:
                print 'Please choose one of the options.'
        except ValueError:
            print 'Please choose one of the options.'

    group_rotations = [False, False, False, False, False]
    def_values = [True, True, True, True, True]
    # default amplitudes (copied from original Fortran script)
    amp_dic = {'ALA': (1.0,),
               'ASN': (0.5, 1.0),
               'ARG': (0.2, 0.3, 0.5),
               'ASP': (0.5, 1.0, 0.0),
               'CYS': (1.0,),
               'GLN': (0.3, 0.5, 1.0),
               'GLU': (0.3, 0.5, 1.0),
               'HIS': (0.3, 0.5),
               'HIE': (0.3, 0.5),
               'HIP': (0.3, 0.5),
               'HID': (0.3, 0.5),
               'ILE': (0.5, 1.0),
               'LEU': (0.5, 1.0),
               'LYS': (0.2, 0.3, 0.5),
               'MET': (0.5, 0.7),
               'PHE': (0.3, 0.5),
               'SER': (1.0,),
               'THR': (1.0,),
               'TRP': (0.3, 0.4),
               'TYR': (0.3, 0.5),
               'VAL': (1.0,),
               'AMB0': (0.1,)  # default for all non residue specific rotations
               }
    group_rotations_NA = [False, False, False, False]
    def_values_NA = [True, True, True, True]
    # default amplitudes (copied from original Fortran script)
    amp_dic_NA = {'A': (0.2, 0.3, 0.3, 0.1),
                  'C': (0.2, 0.3, 0.3, 0.1),
                  'G': (0.2, 0.3, 0.3, 0.1),
                  'U': (0.2, 0.3, 0.3, 0.1),
                  'T': (0.2, 0.3, 0.3, 0.1)
                  }

    # number: [name, ax_atom1, ax_atom2, natom, amp, prob, [atom1, atom2, ...]]
    grot_dic = {}

    if protein_t:
        N_CA = raw_input('Group rotations of side chains about N-CA axis (not PRO) (y/n)? ')
        if N_CA in ['y', 'yes']:
            group_rotations[0] = True
            check_file.write('Group rotation for N - CA' + "\n")
            N_CA_def = raw_input('Use default values (y/n)? ')
            if N_CA_def not in ['y', 'yes']:
                def_values[0] = False
        C_CA = raw_input('Group rotations of side chains about C-CA axis (not PRO) (y/n)? ')
        if C_CA in ['y', 'yes']:
            check_file.write('Group rotation for C - CA' + "\n")
            group_rotations[1] = True
            C_CA_def = raw_input('Use default values (y/n)? ')
            if C_CA_def not in ['y', 'yes']:
                def_values[1] = False
        CA_CB = raw_input('Group rotations of side chains about CA-CB axis (not ALA, GLY, PRO) (y/n)? ')
        if CA_CB in ['y', 'yes']:
            check_file.write('Group rotation for CA - CB' + "\n")
            group_rotations[2] = True
            CA_CB_def = raw_input('Use default values (y/n)? ')
            if CA_CB_def not in ['y', 'yes']:
                def_values[2] = False
        CB_CG = raw_input(
            'Group rotations of side chains about CB-CG axis (not ALA, GLY, PRO, SER, CYS, THR, VAL) (y/n)? ')
        if CB_CG in ['y', 'yes']:
            check_file.write('Group rotation for CB - CG' + "\n")
            group_rotations[3] = True
            CB_CG_def = raw_input('Use default values (y/n)? ')
            if CB_CG_def not in ['y', 'yes']:
                def_values[3] = False
        CG_CD = raw_input('Group rotations of side chains about CG-CD axis (only ARG, LYS, GLU, GLN) (y/n)? ')
        if CG_CD in ['y', 'yes']:
            check_file.write('Group rotation for CG - CD' + "\n")
            group_rotations[4] = True
            CG_CD_def = raw_input('Use default values (y/n)? ')
            if CG_CD_def not in ['y', 'yes']:
                def_values[4] = False

        for res_gr in grot_res:
            gr_res_name = loaded_mol.get_residue_name(res_gr)
            if gr_res_name in protein_res:
                if group_rotations[0] and gr_res_name not in ['PRO','CYX']:
                    # N-CA rotation
                    if def_values[0]:
                        prob = prob_uni
                        amp = 0.1  # use default amplitude
                    else:
                        prob, amp = get_prob_amb()
                    ax1 = loaded_mol.get_atom_id_from_name('N', res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name('CA', res_gr)
                    gr_name = 'NCA' + str(res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ['N',  'C', 'CA', 'HA', 'O', 'OXT', 'H1', 'H2', 'H3']:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations[1] and gr_res_name not in ['PRO','CYX']:
                    # C-CA rotation
                    if def_values[1]:
                        prob = prob_uni
                        amp = 0.1  # use default amplitude
                    else:
                        prob, amp = get_prob_amb()
                    ax1 = loaded_mol.get_atom_id_from_name('C', res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name('CA', res_gr)
                    gr_name = 'CCA' + str(res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ['N',  'C', 'CA', 'HA', 'O', 'OXT', 'H1', 'H2', 'H3']:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations[2] and gr_res_name not in ['PRO', 'ALA', 'GLY','CYX']:
                    # CA-CB rotation
                    if def_values[2]:
                        prob = prob_uni
                        try:
                            amp = amp_dic[gr_res_name][0]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    ax1 = loaded_mol.get_atom_id_from_name('CA', res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name('CB', res_gr)
                    gr_name = 'CACB' + str(res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ['N', 'H', 'C', 'CA', 'HA', 'CB', 'O', 'OXT', 'H1', 'H2', 'H3']:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations[3] and gr_res_name not in ['PRO', 'ILE', 'ALA', 'GLY', 'SER', 'CYS', 'THR', 'VAL','CYX']:
                    # CB-CG rotation
                    if def_values[3]:
                        prob = prob_uni
                        try:
                            amp = amp_dic[gr_res_name][1]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    ax1 = loaded_mol.get_atom_id_from_name('CB', res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name('CG', res_gr)
                    gr_name = 'CBCG' + str(res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ['N', 'H', 'C', 'CA', 'HA', 'CB', 'CG', 'HB', 'HB1', 'HB2',
                                                        'HB3',
                                                        'O', 'OXT', 'H1', 'H2', 'H3']:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations[4] and gr_res_name in ['ARG', 'LYS', 'GLU', 'GLN']:
                    # CG-CD rotation
                    if def_values[4]:
                        prob = prob_uni
                        try:
                            amp = amp_dic[gr_res_name][2]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    ax1 = loaded_mol.get_atom_id_from_name('CG', res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name('CD', res_gr)
                    gr_name = 'CGCD' + str(res_gr)
                    gr_atoms = []
                    atom_exc = ['N', 'H', 'C', 'CA', 'HA', 'CB', 'CG', 'HB', 'HB1', 'HB2', 'HB3', 'O', 'OXT', 'H1',
                                'H2', 'H3', 'CD',
                                'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23']
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in atom_exc:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

    if rna_t or dna_t:

        SugarBaseCheck = raw_input("Group rotations of base about C1'-N axis (y/n)? ")
        if SugarBaseCheck in ['y', 'yes']:
            check_file.write("Group rotation for C1' - N" + "\n")
            group_rotations_NA[0] = True
            CG_CD_def = raw_input('Use default values (y/n)? ')
            if CG_CD_def not in ['y', 'yes']:
                def_values_NA[0] = False

        RotBaseCheck = raw_input("Group rotations of base about intracyclic axis (y/n)? ")
        if RotBaseCheck in ['y', 'yes']:
            check_file.write("Group rotation for intracyclic axis in base" + "\n")
            group_rotations_NA[1] = True
            group_rotations_NA[2] = True
            CG_CD_def = raw_input('Use default values (y/n)? ')
            if CG_CD_def not in ['y', 'yes']:
                def_values_NA[1] = False
                def_values_NA[2] = False

        SugarBackCheck = raw_input("Group rotations about C4'-C3' axis (y/n)? ")
        if SugarBackCheck in ['y', 'yes']:
            check_file.write("Group rotation for C4'-C3' " + "\n")
            group_rotations_NA[3] = True
            CG_CD_def = raw_input('Use default values (y/n)? ')
            if CG_CD_def not in ['y', 'yes']:
                def_values_NA[3] = False

        for res_gr in grot_res:
            gr_res_name = loaded_mol.get_residue_name(res_gr)
            if gr_res_name in ['RAN', 'A', 'A3', 'A5', 'AN', 'RA', 'RA3', 'RA5', 'DA', 'DA5', 'DA3', 'DAN']:
                gr_res_name_dummy = 'A'
            elif gr_res_name in ['RCN', 'C', 'C3', 'C5', 'CN', 'RC', 'RC3', 'RC5', 'DC', 'DC5', 'DC3', 'DCN']:
                gr_res_name_dummy = 'C'
            elif gr_res_name in ['RGN', 'G', 'G3', 'G5', 'GN', 'RG', 'RG3', 'RG5', 'DG', 'DG5', 'DG3', 'DGN']:
                gr_res_name_dummy = 'G'
            elif gr_res_name in ['RCN', 'U', 'U3', 'U5', 'UN', 'RU', 'RU3', 'RU5']:
                gr_res_name_dummy = 'U'
            elif gr_res_name in ['DT', 'DT5', 'DT3', 'DTN']:
                gr_res_name_dummy = 'T'

            if gr_res_name in NA_res:

                if group_rotations_NA[0]:
                    # C1'- N rotation
                    if def_values[1]:
                        try:
                            prob = prob_uni
                            amp = amp_dic_NA[gr_res_name_dummy][0]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    if gr_res_name_dummy in ['A', 'G']:
                        ax1 = loaded_mol.get_atom_id_from_name('N9', res_gr)
                        gr_name = 'N9C1p' + str(res_gr)
                    elif gr_res_name_dummy in ['C', 'U', 'T']:
                        ax1 = loaded_mol.get_atom_id_from_name('N1', res_gr)
                        gr_name = 'N1C1p' + str(res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name("C1'", res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ["P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5''", "C4'", "H4'",
                                                        "O4'", "C1'", "H1'", "C3'", "H3'", "O3'", "C2'", "H2'", "O2'",
                                                        "HO2'", "HO3'", "HO5'", "H2''"]:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations_NA[1]:
                    # intracyclic rotation
                    if def_values[1]:
                        try:
                            prob = prob_uni
                            amp = amp_dic_NA[gr_res_name_dummy][1]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    if gr_res_name_dummy in ['A', 'G']:
                        ax1 = loaded_mol.get_atom_id_from_name('N9', res_gr)
                        ax2 = loaded_mol.get_atom_id_from_name('C6', res_gr)
                        gr_name = 'N9C6' + str(res_gr)
                    elif gr_res_name_dummy in ['C', 'U', 'T']:
                        ax1 = loaded_mol.get_atom_id_from_name('N1', res_gr)
                        ax2 = loaded_mol.get_atom_id_from_name('C4', res_gr)
                        gr_name = 'N1C4' + str(res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ["P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5''", "C4'", "H4'",
                                                        "O4'", "C1'", "H1'", "C3'", "H3'", "O3'", "C2'", "H2'", "O2'",
                                                        "HO2'", "HO3'", "HO5'", "H2''"]:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations_NA[2]:
                    # intracyclic rotation 2
                    if def_values[2]:
                        try:
                            prob = prob_uni
                            amp = amp_dic_NA[gr_res_name_dummy][2]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    if gr_res_name_dummy in ['A', 'G']:
                        ax1 = loaded_mol.get_atom_id_from_name('N9', res_gr)
                        ax2 = loaded_mol.get_atom_id_from_name('N1', res_gr)
                        gr_name = 'N9N1' + str(res_gr)
                    elif gr_res_name_dummy in ['C', 'U', 'T']:
                        ax1 = loaded_mol.get_atom_id_from_name('N1', res_gr)
                        ax2 = loaded_mol.get_atom_id_from_name('N3', res_gr)
                        gr_name = 'N1N3' + str(res_gr)
                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ["P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5''", "C4'", "H4'",
                                                        "O4'", "C1'", "H1'", "C3'", "H3'", "O3'", "C2'", "H2'", "O2'",
                                                        "HO2'", "HO3'", "HO5'", "H2''"]:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

                if group_rotations_NA[3]:
                    # sugarback bone rotation
                    if def_values[3]:
                        try:
                            prob = prob_uni
                            amp = amp_dic_NA[gr_res_name_dummy][3]  # use default amplitude
                        except KeyError:
                            print 'No default data available for ' + gr_res_name
                            prob, amp = get_prob_amb()
                    else:
                        prob, amp = get_prob_amb()
                    ax1 = loaded_mol.get_atom_id_from_name("C4'", res_gr)
                    ax2 = loaded_mol.get_atom_id_from_name("C3'", res_gr)
                    gr_name = 'N9N1' + str(res_gr)

                    gr_atoms = []
                    for gr_atom in loaded_mol.get_atom_list(res_gr):
                        if atom_dic[gr_atom][0] not in ["P", "OP1", "OP2", "O5'", "C5'", "H5'", "H5''", "C4'",
                                                        "C3'", "H3'", "O3'", "HO3'", "HO5'"]:
                            gr_atoms.append(gr_atom)
                    grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

    # user defined groups
    gr_user_inp = raw_input('Add user defined rotational groups? ')
    if gr_user_inp in ['y', 'yes']:
        gr_user_t = True
    else:
        gr_user_t = False
    useri = 1
    while gr_user_t:
        print 'No tests for valid groups!'
        gr_name = 'USERD' + str(useri)
        try:
            ax1 = int(raw_input('Atom number for axis1: '))
            ax2 = int(raw_input('Atom number for axis2: '))
            prob, amp = get_prob_amb()
            gr_atoms = [int(x) for x in (raw_input('List of atoms in group separated by spaces: ')).split()]
        except:
            print 'Invalid input'
        else:
            grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]
        finally:
            gr_user_inp = raw_input('Add user defined rotational group? ')
            if gr_user_inp in ['y', 'yes']:
                gr_user_t = True
            else:
                gr_user_t = False

    # write output for group rotation files
    groups_file = open('atomgroups', 'w')
    for group in range(1, len(grot_dic.keys()) + 1):
        string = 'GROUP ' + grot_dic[group][0] + ' ' + str(grot_dic[group][1]) + ' ' + str(grot_dic[group][2])
        string += ' ' + str(grot_dic[group][3]) + ' ' + str(grot_dic[group][4]) + ' ' + str(grot_dic[group][5]) + "\n"
        for atom in grot_dic[group][6]:
            string += str(atom) + "\n"
        groups_file.write(string)

check_file.close()

# Writing in file parmed.py to exclude interactions ####
parmed_in = open('parmed_in', "w")
for group in range(1, len(groups.keys()) + 1):
    string = 'addExclusions '
    string2 = '@'
    counter = 1
    for res in groups[group]:
        for atom in res_dic[res]:
            if counter == 1:
                string2 += str(atom)
                counter = 2
            else:
                string2 += ',' + str(atom)
    parmed_in.write(string + string2 + ' ' + string2 + "\n")
for group_local in range(1, len(groups_local.keys()) + 1):
    string = 'addExclusions '
    string2 = '@'
    counter = 1
    for atom in groups_local[group_local][2:]:
        if counter == 1:
            string2 += str(atom)
            counter = 2
        else:
            string2 += ',' + str(atom)
    parmed_in.write(string + string2 + ' ' + string2 + "\n")
parmed_in.write('parmout coords.prmtop.ni' + "\n" + 'go' + "\n")

print 'Selection process completed - change topology file to exclude interactions within the rigid bodies.'

if pymol_check:
    mol_view = PymolView(pdb_inp)
    mol_view.apply_colour_scheme(res_atomistic, res_local, 2, 3)
    pymol.cmd.save('rigid_%s' % pdb_inp, format='png')
    time.sleep(10)
    mol_view.__del__()
