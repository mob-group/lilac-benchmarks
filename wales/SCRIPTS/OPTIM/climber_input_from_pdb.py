# Script to get a list of atoms between which constraints are allowed for CLIMBERINT

import sys

inp_pdb = sys.argv[1]
bbconstr = []
scconstr = []
with open(inp_pdb , "r") as f:
    for line in f:
        if line.split()[0] == 'ATOM':
           atomname=line[12:16].strip()
           if atomname=='CA':
              bbconstr.append(line[6:12].strip())
           else:
             if (atomname[0]!='H') and (atomname not in ['C','O','N']):
                scconstr.append(line[6:12].strip())
out_file=open("climber.constraints" , "w")
string = ''
string = string + '%8s' % str(len(bbconstr))
string = string + "\n"
for i in xrange(len(bbconstr)):
    string = string +'%8s' %  bbconstr[i]
    if (i+1)%10==0:
       string = string + "\n"

if len(bbconstr)%10==0:
   string = string + '%8s' % str(len(scconstr))
else:
   string = string + "\n"
   string = string + '%8s' % str(len(scconstr))

string = string + "\n"
for i in xrange(len(scconstr)):
    string = string + '%8s' % scconstr[i] 
    if (i+1)%10==0:
       string = string + "\n"
out_file.write(string)
out_file.close()
