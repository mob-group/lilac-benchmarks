import sys

inp_f = sys.argv[1]

peptide_bond_atom = {}
ter_list = []
i = 0

print 'Only for single chains at the moment, check if the system has multiple chains'

with open(inp_f , "r") as f:
    for line in f:
        if line.split()[0] == 'ATOM':
            i = i + 1
            if line[12:15].strip() in ['H','CA','O']:
                peptide_bond_atom[i] = line[12:15].strip()
        elif line.split()[0] == 'TER':
            ter_list.append((i, i + 1))
        else:
            continue

peptide_bonds = {}
for key in sorted(peptide_bond_atom.keys()):
    if peptide_bond_atom[key] == 'CA':
        if len(peptide_bonds) != 0:
            peptide_bonds[len(peptide_bonds)].append(key)
        peptide_bonds[len(peptide_bonds) + 1] = [key]
    else:
        if len(peptide_bonds.keys()) > 0:
            peptide_bonds[len(peptide_bonds)].append(key)
        else:
            continue 
out_f = open("constraintfile" , "w")
for key in peptide_bonds.keys():
    if len(peptide_bonds[key]) == 4:
        l = peptide_bonds[key]
        out_f.write('%6d %6d \n' %(l[0],l[2]))
        out_f.write('%6d %6d \n' %(l[0],l[3]))
        out_f.write('%6d %6d \n' %(l[1],l[2]))
        out_f.write('%6d %6d \n' %(l[1],l[3]))
out_f.close()

