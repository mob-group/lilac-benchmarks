import numpy as np

mass = {'H':1.008,'C':12.01,'N':14.007,'O':15.999,'S':32.065}

top_f = "coords.prmtop"
start_f = "start"
names = False
atomnames = ''
#get coordinates
coords = np.genfromtxt(start_f , dtype=float)
#get atom type list
with open(top_f , "r") as f:
   for line in f:
      if line.split()[0] == '%FLAG':
         if line.split()[1] =='ATOM_NAME':
            names = True
         elif line.split()[1] == 'CHARGE':
            break
         else:
            continue
      elif names:
         if line[:7] == '%FORMAT':
            continue
         else:
            atomnames += line[::4].rstrip()

atomnames = list(atomnames)

#calculate centre of mass

atommass = []
totmass = 0.0
for atom in atomnames:
   atommass.append(mass[atom])
   totmass += mass[atom]

X = np.sum(coords[:,0]*atommass)/totmass
Y = np.sum(coords[:,1]*atommass)/totmass
Z = np.sum(coords[:,2]*atommass)/totmass
coords[:] -= [X,Y,Z]

#output files
out_f = open('start.new' , "w")
for xyz in coords[:]:
    out_f.write('%20.10f %20.10f %20.10f \n' %(xyz[0],xyz[1],xyz[2]))
out_f.close()
print 'Centre of mass shifted to origin'
