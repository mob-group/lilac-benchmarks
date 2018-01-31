#! /usr/bin/env python
import sys
#written by by James Klimavicz jk669 March 2017 to make group rotation rotations along sugar phosphate backbone.

#reads in pdb file (hopefully)
def read_pdb(pdb,exit):
	pdb_atom_list=[]
	pdbopen=False
	while (not pdbopen):
		try:
			with open(pdb,'r') as f:
				for line in f:
					s=line.split()
					if (s[0]=='ATOM'):
						pdb_atom_list.append((s[1],s[2],s[3],s[4]))
			pdbopen=True
		except IOError:
			pdb=raw_input("{} could not be opened. Please type the AMBER-produced pdb filepath, or 'exit' to end the program: ".format(pdb))
			if pdb.lower() in exit:
				quit()
	return pdb_atom_list

#gets user input value for probability, and makes sure it is a reasonable value.
def get_prob():
	global accept
	prob_in= raw_input('Input a Probability of group rotation (suggested: 0.025): ')
	try: 
		prob=float(prob_in)
		if (prob>0.01 and prob<=0.10):
			prob_allow=True
		else:
			prob_allow=False
		while not prob_allow:
			if prob > 0.1:
				print('This is a large probability (values recommended: 0.01 to 0.1).')
				prob_check = raw_input(' Continue with the selected value (y/n)? ')
				if prob_check in accept:
					prob_allow=True
			elif prob < 0.01:
				print('This is a small probability (values recommended: 0.01 to 0.1).')
				prob_check = raw_input(' Continue with the selected value (y/n)? ')
				if prob_check in accept:
					prob_allow=True
			else:
				prob_allow=True
			if not prob_allow:
				prob_in= raw_input('Input a Probability of group rotation (suggested: 0.025): ')
				try: 
					prob=float(prob_in)
				except ValueError:
					print('Probability needs to be a float')
	except ValueError:
		print('Input needs to be a float')
	return prob

#gets user input value for amplitude, and makes sure it is a reasonable value.
def get_amp():
	global accept
	global reject
	amp_in= raw_input('Input the amplitude of group rotation (suggested: 0.05): ')
	try: 
		amp=float(amp_in)
		if (amp<=1.0 and amp>=0.0):
			amp_allow=True
		else:
			amp_allow=False
		while not amp_allow:
			if amp > 1.0:
				print('The amplitude should be in the half-open interval (0,1].')
				amp_check = raw_input('Set the amplitude to 1 (y/n)? ')
				if amp_check in accept:
					amp=1.0
					amp_allow=True
			elif amp <= 0.00:
				print('The amplitude should be in the half-open interval (0,1].')
				amp_check = raw_input('Set to default value of 0.025 (y/n)? ')
				if amp_check not in reject:
					amp=0.025
					amp_allow=True
			else:
				amp_allow=True
			if not amp_allow:
				amp_in= raw_input('Input the amplitude of group rotation (suggested: 0.05): ')
				try: 
					amp=float(amp_in)
				except ValueError:
					print('Amplitude needs to be a float')
	except ValueError:
	    print('Input needs to be a float')
	return amp

#####################################################################################

backbone_rotation = [False, False, False]
prob_def = [0.025, 0.025, 0.025]
amp_def = [0.05, 0.05, 0.05]
def_change_check=False
accept = ['y','yes', 'affirmative','yeah', 'accept']
exit = ['exit', 'break', 'stop', 'die', 'quit', 'kill']
reject = ['n','no','negative','nope', 'reject']

pdb=raw_input("Please type the AMBER-produced pdb filepath, or 'exit' to end the program: ")
pdb_atom_list=read_pdb(pdb,exit)
nres=int(pdb_atom_list[-1][3])
natoms=len(pdb_atom_list)

def_change = (raw_input("Do you wish to change any default probabilities (0.025) or amplitudes (0.025) (y/n)? ")).lower()
if def_change in accept:
	def_change_check=True
elif def_change in exit:
	quit()

O3Pcheck = (raw_input("Group rotations of base about O3'-P axis (y/n)? ")).lower()
if O3Pcheck in accept:
	backbone_rotation[0] = True
	if def_change_check:
		def_prob = (raw_input('Use default values (P=0.025, amplitude =0.05) (y/n)? ')).lower()
		if def_prob in reject:
			prob_def[0] = get_prob()
			amp_def[0]=get_amp()
		elif def_prob in exit:
			quit()
elif O3Pcheck in exit:
	quit()

PO5check = (raw_input("Group rotations of base about P-O5' axis (y/n)? ")).lower()
if PO5check in accept:
	backbone_rotation[1] = True
	if def_change_check:
		def_prob = (raw_input('Use default values (P=0.025, amplitude =0.05) (y/n)? ')).lower()
		if def_prob in reject:
			prob_def[1] = get_prob()
			amp_def[1]=get_amp()
		elif def_prob in exit:
			quit()
elif PO5check in exit:
	quit()

O5C5check = (raw_input("Group rotations of base about O5'-C5' axis (y/n)? ")).lower()
if O5C5check in accept:
	backbone_rotation[2] = True
	if def_change_check:
		def_prob = (raw_input('Use default values (P=0.025, amplitude =0.05) (y/n)? ')).lower()
		if def_prob in reject:
			prob_def[2] = get_prob()
			amp_def[2]=get_amp()
		elif def_prob in exit:
			quit()
elif O5C5check in exit:
	quit()

if (sum(backbone_rotation)==0):
	print('No group rotations were selected. Exiting gracefully without file creation.')
	quit() #no point in doing the rest if we don't actually want any rotations....

O3_ind=[None]*nres
P_ind=[None]*nres
O5_ind=[None]*nres
C5_ind=[None]*nres
for i in range(0,natoms):
	ind=i+1
	atom=pdb_atom_list[i][1]
	res=int(pdb_atom_list[i][3])-1
	if (atom=="O3'"):
		O3_ind[res]=ind
	elif (atom=='P'):
		P_ind[res]=ind
	elif (atom=="O5'"):
		O5_ind[res]=ind
	elif (atom=="C5'"):
		C5_ind[res]=ind

#write new atomgroups file. This is named different so as not to overwrite any other atomgroups file in the 
#folder; if desired, the atromgroups_phos can be concatenated onto an existing atomgroups file with bash commands
print('The default atomgroups name on output is atomgroups_phos so that this will not overwrite any existing atomgroups file.')
print('This can be concatenated onto an existing atomgroups file with the command "cat atomgroups_phos >> atomgroups".')
changename=raw_input('Type "yes" or "y" if you would like to change the output file name. ')
if changename in accept:
	filename=raw_input('Type the filename you would like: ')
elif changename in exit:
	quit()
else: 
	filename='atomgroups_phos'
ag=open(filename,'w')

for i in range(0,nres-1):
	if (backbone_rotation[0]):
		name = "O3pP{}".format(i+1)
		ax1 = O3_ind[i]
		ax2 = P_ind[i+1]
		if (ax2<natoms/2):
			atomlist=range(1,O3_ind[i])
		else:
			atomlist=range(O5_ind[i+1]+1,natoms+1)
		string="GROUP {} {} {} {} {} {}\n".format(name, ax1, ax2, len(atomlist), amp_def[0], prob_def[0])
		ag.write(string)
		for k in atomlist:
			ag.write("{}\n".format(k))

	if (backbone_rotation[1]):
		name = "PO5p{}".format(i+1)
		ax1 = P_ind[i+1]
		ax2 = O5_ind[i+1]
		if (ax2<natoms/2):
			atomlist=list(range(1,P_ind[i+1])) + list(range(P_ind[i+1]+1,O5_ind[i+1]))
		else:
			atomlist=range(O5_ind[i+1]+1,natoms+1)
		string="GROUP {} {} {} {} {} {}\n".format(name, ax1, ax2, len(atomlist), amp_def[1], prob_def[1])
		ag.write(string)
		for k in atomlist:
			ag.write("{}\n".format(k))

	if (backbone_rotation[2]):
		name = "O5pC5p{}".format(i+1)
		ax1 = O5_ind[i+1]
		ax2 = C5_ind[i+1]
		if (ax2<natoms/2):
			atomlist=range(1,O5_ind[i+1])
		else:
			atomlist=range(C5_ind[i+1]+1,natoms+1)
		string="GROUP {} {} {} {} {} {}\n".format(name, ax1, ax2, len(atomlist), amp_def[2], prob_def[2])
		ag.write(string)
		for k in atomlist:
			ag.write("{}\n".format(k))

ag.close()
print('File {} was created successfully.'.format(filename))