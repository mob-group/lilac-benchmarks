import perm_data_amino_acids as perm
import sys

oldname = sys.argv[1]
newname = sys.argv[2]
start_atom = int(sys.argv[3])

perm_old = dict()
perm_new = dict()
perm_shift = dict()
perm_add = dict()
perm_add_new = dict()
perm_add_shift = dict()
with open("perm.allow" , "r") as f:
    next(f)
    for line in f:
        perm_add[len(perm_add) + 1] = line
        line = f.next()
        perm_old[len(perm_old) + 1] = [int(x) for x in line.split()]

n_old = perm.atomdata[oldname][0]
groups_old = perm.atomdata[oldname][1:]
n_new = perm.atomdata[newname][0]
groups_new = perm.atomdata[newname][1:]

shift = n_new - n_old
final_atom = start_atom + n_old - 1

groupt = False
shiftt = False

for key in perm_old.keys():
    for atom in perm_old[key]:
        if atom > start_atom:
            groupt = True
        if atom > final_atom:
            shiftt = True
    if not groupt:
        perm_new[key] = perm_old[key]
        perm_add_new[key] = perm_add[key]
    if groupt and shiftt:
        perm_shift[key] = perm_old[key]
        perm_add_shift[key] = perm_add[key]

for group in groups_new:
    if isinstance(group[-1] , tuple):
        perm_add_new[len(perm_add_new) + 1] = str(group[-1][0]) + ' ' + str(group[-1][1]) + "\n"
        perm_new[len(perm_new) + 1] = [(start_atom - 1 + x) for x in group[0]]
    else:
        perm_add_new[len(perm_add_new) + 1] = str(len(group)) + ' 0' +"\n"       
        perm_new[len(perm_new) + 1] = [(start_atom - 1 + x) for x in group]
    
for key in perm_shift.keys():
    perm_new[len(perm_new) + 1] = [(x + shift) for x in perm_shift[key]]        
    perm_add_new[len(perm_add_new) + 1] = perm_add_shift[key]

output = open("perm.allow.new" , "w")    
output.write(str(len(perm_new.keys())) + "\n")
for key in perm_new.keys():
    output.write(perm_add_new[key]) 
    atoms = perm_new[key]
    string = str()
    for atom in atoms:
        string += ' ' + str(atom)
    string += "\n"        
    output.write(string)

output.close()
        
    
