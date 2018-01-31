import coordinates_mut as cmut
import numpy as np
import sys

def assignment_from_list(aa_tree , atomlist):
    for atom in atomlist:
        hosts = cmut.find_neighbours_with_coords(aa_tree , atom)
        if len(hosts) == 1:
            cmut.create_new_coords_one_host(aa_tree, atom, hosts[0])
        else:
            cmut.remove_multiple_hosts(aa_tree, hosts)

oldname = sys.argv[1]
newname = sys.argv[2]
 
coords = np.genfromtxt("coords.oldres") 
#coords = coord.coords_aa[oldname]

#Create graphs for old and new residue and reassign the conserved coordinates
old_residue = cmut.create_tree_for_aa(oldname)
cmut.add_coordinates_to_residue(old_residue , coords)
new_residue = cmut.create_tree_for_aa(newname)
cmut.assign_coordinates_backbone(old_residue , new_residue)
cmut.assign_coordinates_identical_atoms(old_residue , new_residue)
#set correct chirality for GLY
if ((oldname == 'GLY' and newname != 'GLY') or 
    (oldname == 'CGLY' and newname != 'CGLY') or 
    (oldname == 'NGLY' and newname != 'NGLY')):
    cmut.mutate_GLY_chiral(old_residue , new_residue)
#remove crowding for PRO and HYP
if (((oldname in ['PRO','HYP']) and (newname not in ['HYP' , 'PRO'])) or 
    ((oldname in ['CPRO','CHYP']) and (newname not in ['CHYP' , 'CPRO']))):
    try:
        CG = cmut.get_number_from_name(new_residue , 'CG')
        del new_residue.node[CG]['XYZ']
        CD = cmut.get_number_from_name(new_residue , 'CD')
        del new_residue.node[CD]['XYZ']
    except (IndexError , KeyError):
        pass

if (((newname in ['PRO' , 'HYP']) and (oldname not in ['HYP' , 'PRO'])) or
    ((newname in ['CPRO' , 'CHYP']) and (oldname not in ['CHYP' , 'CPRO']))):
    for atom in ['CG' , 'CD' , 'HG' ,'HG2' , 'HG3' , 'HD2' , 'HD3' , 'OD1' , 'HD1']:
        try:
            id = cmut.get_number_from_name(new_residue , atom)
            del new_residue.node[id]['XYZ']
        except (IndexError , KeyError):
            pass

#remove clashes due to position of CG (this may not be necessary ...)
if ((newname == 'TRP') or (newname == 'TYR') or
    (newname == 'CTRP') or (newname == 'CTYR') or
    (newname == 'NTRP') or (newname == 'NTYR')):
    cmut.position_CG_for_TRP(old_residue , new_residue , oldname)

#Now check for rings in the amino acid
ringt = cmut.check_for_rings(new_residue)
atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
if not ringt:
    while len(atoms_unassigned) > 0:
        assignment_from_list(new_residue , atoms_unassigned)
        atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
        if (newname == 'ILE') and (len(atoms_unassigned) == 0):            
            cmut.check_chirality_ILE(new_residue)
            atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
        if (newname == 'THR') and (len(atoms_unassigned) == 0):            
            cmut.check_chirality_THR(new_residue)
            atoms_unassigned = cmut.find_atoms_without_coords(new_residue)

elif len(atoms_unassigned) > 0:
    #PHE <--> TYR and TRP <--> HIS
    specialmutt = cmut.check_special_assignment(oldname,newname)
    if specialmutt:
        cmut.special_assignment((old_residue,new_residue) , (oldname,newname))
        atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
        while len(atoms_unassigned) > 0:
            assignment_from_list(new_residue , atoms_unassigned)
            atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
    elif newname not in ["PRO" , "HYP" , "CPRO" , "CHYP" , "NPRO"]:
        #locate first atom belonging to ring ("CG" for everything but proline) 
        for atom in new_residue.nodes():
            if new_residue.node[atom]['name'] == "CG":
                first_ring_atom = atom
        atoms_to_assign = cmut.get_atom_list_upto_ring(new_residue , first_ring_atom)
        while len(atoms_to_assign) > 0:
            assignment_from_list(new_residue , atoms_to_assign)
            atoms_to_assign = cmut.get_atom_list_upto_ring(new_residue , first_ring_atom)
        cmut.assign_coordinates_ring(new_residue , newname)
        atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
        while len(atoms_unassigned) > 0:
            assignment_from_list(new_residue , atoms_unassigned)
            atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
    else:
        if (((oldname == 'PRO') and (newname == 'HYP')) or
            ((oldname == 'CPRO') and (newname == 'CHYP'))):
            cmut.mutate_HYP_from_PRO(old_residue , new_residue)  
        elif (((oldname == 'HYP') and (newname == 'PRO')) or
              ((oldname == 'CHYP') and (newname == 'CPRO'))):
            pass    
        else:
            cmut.create_proline_coords(new_residue)
        atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
        while len(atoms_unassigned) > 0:
            assignment_from_list(new_residue , atoms_unassigned)
            atoms_unassigned = cmut.find_atoms_without_coords(new_residue)
            if ((newname == 'HYP') or (newname == 'CHYP')) and (len(atoms_unassigned) == 0):            
                cmut.check_chirality_HYP(new_residue)
                atoms_unassigned = cmut.find_atoms_without_coords(new_residue)            

output_file = open("coords.newres" , "w")
for atom in new_residue.nodes():
    coords = new_residue.node[atom]['XYZ']
    output_file.write("%14.7f %14.7f %14.7f \n" %(coords[0] , coords[1] , coords[2]))
output_file.close()    
    
