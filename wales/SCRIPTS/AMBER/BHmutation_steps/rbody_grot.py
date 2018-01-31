import sys

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

    # atom id from atom name given a residue id
    def get_atom_id_from_name(self, atom_name, res_num):
        for i in self.res_dic[res_num]:
            if self.atom_dic[i][0] == atom_name:
                return i

    # returns a list of atoms in a residue
    def get_atom_list(self, residue):
        all_atoms = self.atom_dic.keys()
        atomlist = []
        for i in range(1, len(self.atom_dic) + 1):
            atom = self.atom_dic[i]
            if residue == atom[2]:
                atomlist.append(all_atoms[i - 1])
        return atomlist

#################################################
def reading_pdb(pdb):
    dic = {}
    atom_list = []
    with open(pdb, 'r') as f_pdb:
        for line in f_pdb:
            if line.split()[0] == 'ATOM':
                atom_list.append((line[12:16].strip(), line[17:20].strip(), int(line[22:26])))
    for i in range(1, len(atom_list) + 1):
        dic[i] = atom_list[i - 1]
    return dic

def reading_rbodyconfig(fname):
    dic = dict()
    with open(fname , "r") as f:
        for line in f:
            if line.split()[0] == 'GROUP':
                dic[len(dic) + 1] = list()
            else:
                dic[len(dic)].append(int(line.split()[0]))
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

#################################################       

prob_uni = 0.025

rbody_old = sys.argv[1]
shift = int(sys.argv[2])
start = int(sys.argv[3])
pdb_file = sys.argv[4]
atom_dic = reading_pdb(pdb_file)
res_dic = create_res_dic(atom_dic)
loaded_mol = protein(atom_dic, res_dic)

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

#first read in old rbodyconfig and adjust it to give new output
rb_dict = reading_rbodyconfig(rbody_old)
used_atoms = list()
for key in rb_dict.keys():
    if rb_dict[key][0] > start:
        rb_dict[key] = [(atom + shift) for atom in rb_dict[key]]
    used_atoms = used_atoms + rb_dict[key]

rboutput = open("rbodyconfig.new" , "w")
for key in rb_dict.keys():
    rboutput.write("GROUP " + str(len(rb_dict[key])) + "\n")
    for atom in rb_dict[key]:
        rboutput.write(str(atom) + "\n")
rboutput.close()

# number: [name, ax_atom1, ax_atom2, natom, amp, prob, [atom1, atom2, ...]]
grot_dic = {}
residues = res_dic.keys()
for res in residues:
    gr_res_name = loaded_mol.get_residue_name(res)
    if gr_res_name not in ['PRO' , 'HYP' , 'CYX']:
        #N-CA rotation
        prob = prob_uni
        amp = 0.1  # use default amplitude

        ax1 = loaded_mol.get_atom_id_from_name('N', res)
        ax2 = loaded_mol.get_atom_id_from_name('CA', res)
        gr_name = 'NCA' + str(res)
        gr_atoms = []
        for gr_atom in loaded_mol.get_atom_list(res):
            if atom_dic[gr_atom][0] not in ['N',  'C', 'CA', 'HA', 'O', 'OXT', 'H1', 'H2', 'H3']:
                gr_atoms.append(gr_atom)
        inrbody = False
        for atom in gr_atoms:
            if atom in used_atoms:
                inrbody = True
        if not(inrbody):
            grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

    if gr_res_name not in ['PRO' , 'HYP' , 'CYX']:
        #C-CA rotation
        prob = prob_uni
        amp = 0.1  # use default amplitude
        ax1 = loaded_mol.get_atom_id_from_name('C', res)
        ax2 = loaded_mol.get_atom_id_from_name('CA', res)
        gr_name = 'CCA' + str(res)
        gr_atoms = []
        for gr_atom in loaded_mol.get_atom_list(res):
            if atom_dic[gr_atom][0] not in ['N',  'C', 'CA', 'HA', 'O', 'OXT', 'H1', 'H2', 'H3']:
                gr_atoms.append(gr_atom)
        inrbody = False
        for atom in gr_atoms:
            if atom in used_atoms:
                inrbody = True
        if not(inrbody):
            grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

    if gr_res_name not in ['PRO', 'HYP' , 'ALA', 'GLY' , 'CYX']:
        # CA-CB rotation
        prob = prob_uni
        try:
            amp = amp_dic[gr_res_name][0]  # use default amplitude
        except KeyError:
            amp = 0.1
        ax1 = loaded_mol.get_atom_id_from_name('CA', res)
        ax2 = loaded_mol.get_atom_id_from_name('CB', res)
        gr_name = 'CACB' + str(res)
        gr_atoms = []
        for gr_atom in loaded_mol.get_atom_list(res):
            if atom_dic[gr_atom][0] not in ['N', 'H', 'C', 'CA', 'HA', 'CB', 'O', 'OXT', 'H1', 'H2', 'H3']:
                gr_atoms.append(gr_atom)
        inrbody = False
        for atom in gr_atoms:
            if atom in used_atoms:
                inrbody = True
        if not(inrbody):
            grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

    if gr_res_name not in ['PRO', 'HYP' , 'ILE', 'ALA', 'GLY', 'SER', 'CYS', 'THR', 'VAL' , 'CYX']:
        # CB-CG rotation
        prob = prob_uni
        try:
            amp = amp_dic[gr_res_name][1]  # use default amplitude
        except KeyError:
            amp = 0.1
        ax1 = loaded_mol.get_atom_id_from_name('CB', res)
        ax2 = loaded_mol.get_atom_id_from_name('CG', res)
        gr_name = 'CBCG' + str(res)
        gr_atoms = []
        for gr_atom in loaded_mol.get_atom_list(res):
            if atom_dic[gr_atom][0] not in ['N', 'H', 'C', 'CA', 'HA', 'CB', 'CG', 'HB', 'HB1', 'HB2',
                                            'HB3', 'O', 'OXT', 'H1', 'H2', 'H3']:
                gr_atoms.append(gr_atom)
        inrbody = False
        for atom in gr_atoms:
            if atom in used_atoms:
                inrbody = True
        if not(inrbody):
            grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

    if gr_res_name in ['ARG', 'LYS', 'GLU', 'GLN']:
        # CG-CD rotation
        prob = prob_uni
        try:
            amp = amp_dic[gr_res_name][2]  # use default amplitude
        except KeyError:
            amp = 0.1
        ax1 = loaded_mol.get_atom_id_from_name('CG', res)
        ax2 = loaded_mol.get_atom_id_from_name('CD', res)
        gr_name = 'CGCD' + str(res)
        gr_atoms = []
        atom_exc = ['N', 'H', 'C', 'CA', 'HA', 'CB', 'CG', 'HB', 'HB1', 'HB2', 'HB3', 'O', 'OXT', 'H1',
                    'H2', 'H3', 'CD', 'HG1', 'HG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23']
        for gr_atom in loaded_mol.get_atom_list(res):
            if atom_dic[gr_atom][0] not in atom_exc:
                gr_atoms.append(gr_atom)
                inrbody = False
        for atom in gr_atoms:
            if atom in used_atoms:
                inrbody = True
        if not(inrbody):
            grot_dic[len(grot_dic) + 1] = [gr_name, ax1, ax2, len(gr_atoms), amp, prob, gr_atoms]

  
# write output for group rotation files
groups_file = open('atomgroups', 'w')
for group in range(1, len(grot_dic.keys()) + 1):
    string = 'GROUP ' + grot_dic[group][0] + ' ' + str(grot_dic[group][1]) + ' ' + str(grot_dic[group][2])
    string += ' ' + str(grot_dic[group][3]) + ' ' + str(grot_dic[group][4]) + ' ' + str(grot_dic[group][5]) + "\n"
    for atom in grot_dic[group][6]:
        string += str(atom) + "\n"
    groups_file.write(string)
groups_file.close()



