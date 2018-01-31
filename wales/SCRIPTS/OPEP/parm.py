# !/usr/bin/env python
# -*- coding: UTF8 -*-

import argparse
import math
import os
import os.path
import shutil
import string
import sys

parser = argparse.ArgumentParser(description='Parse parameter generator arguments')
parser.add_argument("pdbfilename")
parser.add_argument('-p', '--parmfile', type=str, default=None)
parser.add_argument('-o', '--outbase', type=str, default="Output")
parser.add_argument('--local', action='store_true')

args = parser.parse_args()


# atom name to weight type
atn2ty = {"P":"P", "O5*":"O", "C5'":"C", "C5*":"C", "CA":"R4", "C4'":"R4", "C4*":"R4", "C1'":"R1", "C1*":"R1", "CY":"R1", "A1":"A1", "A2":"A2", "G1":"G1", "G2":"G2", "C1":"C1", "T1":"U1", "D":"D", "MG":"MG"}
# Main, sugar, end, or heteroatom ?
aty2mse = {"P":"M", "O":"M", "C":"M", "R4":"M", "R1":"S", "A1":"S", "A2":"E", "G1":"S", "G2":"E", "C1":"E", "U1":"E", "D":"M", "MG":"H"}
# residue letter to number
resl2i_ = {'G':1, 'A':2, 'C':3, 'U':4, 'D':5, 'DG':1, 'DA':2, 'DC':3, 'DU':4, 'DT':4}
# residue letter to name
resl2n_ = {'G':"GUA", 'A':"ADE", 'C':"CYT", 'U':"URA", 'D':"DUM", 'DG':"GUA", 'DA':"ADE", 'DC':"CYT", 'DU':"URA", 'DT': "THY"}

def res_cut(rn, rdict):
    # if we have a residue name C5,G5,C3,G3, etc... remove the number
    if rn[-1] in "53":
      rn = rn[:-1]
    return rdict[rn]

def resl2i(rn):
  return res_cut(rn, resl2i_)
def resl2n(rn):
  return res_cut(rn, resl2n_)

class Atom():
  def __init__(self, pdbline):
    self.idx = int(pdbline[6:11])
    self.name = pdbline[12:16]
    self.aty = atn2ty[self.name.strip()]
    self.resname = pdbline[17:20]
    self.chainid = pdbline[21:22]
    self.residx = int(pdbline[22:26])
    self.x = float(pdbline[30:38])
    self.y = float(pdbline[38:46])
    self.z = float(pdbline[46:54])

"""
Represent one chain in a PDB file
"""
class Chain():
  def __init__(self):
    self.atlist = list()
    self.reslist = list()

"""
Degrees to radians conversion
"""
def deg2rad(f):
  return f/180*math.pi

"""
Wrap a long text line on 80 characters
"""
def wraptxt(txt, width=80, ls=False):
  i = 0
  txt2 = ""
  while i+width < len(txt):
    i0 = i
    i += width
    # Look for either beginning or end of word, i.e space+letter or letter+space
    while (txt[i] in string.whitespace) == (txt[i-1] in string.whitespace):
      i -= 1
      if i <= 0:
        raise "Error"
    if ls:
      txt2 += txt[i0:i].lstrip()+'\n'
    else:
      txt2 += txt[i0:i]+'\n'
  
  if ls:
    txt2 += txt[i:].lstrip()
  else:
    txt2 += txt[i:]
  
  return txt2


def fortranfloat(f):
  s = "{:12.8E}".format(f)
  return "% 9.8fE%+03d" % (float(s[:-4])/10, int(s[-3:])+1)


######### read parameter file
if args.parmfile is None:
  parmfile =  os.path.join(os.path.dirname(sys.argv[0]),"parm-OPEP3RNA.dat")
else:
  parmfile = args.parmfile

pf = open(parmfile, 'r')

atl = list()
atw = dict()
pf = pf.readlines()

i=1
for l in pf[1:]:
  #if l[0] == '\n':
  if not l.strip():
    break
  atl.append(l.split()[0])
  atw[l.split()[0]] = float(l.split()[1])
  i += 1

i +=3
bonds = []
for l in pf[i:]:
  #if l[0] == '\n':
  if not l.strip():
    break
  b = l.split("-")
  b2 = b[1].split()
  rk = float(b2[1])
  req = float(b2[2])
  b1 = b[0].strip()
  b2 = b2[0].strip()
  bonds.append((b1, b2, rk, req))
  i += 1

i += 1
angles = []
for l in pf[i:]:
  #if l[0] == '\n':
  if not l.strip():
    break
  b = l.split("-")
  b3 = b[2].split()
  rk = float(b3[1])
  req = deg2rad(float(b3[2]))
  b1 = b[0].strip()
  b2 = b[1].strip()
  b3 = b3[0].strip()
  angles.append((b1, b2, b3, rk, req))
  i += 1

i += 1
diheds = []
for l in pf[i:]:
  #if l[0] == '\n':
  if not l.strip():
    break
  b = l[:13].split("-")
  b2 = l[13:].split()
  pk = float(b2[1])
  phase = deg2rad(float(b2[2]))
  pn = float(b2[3])
  b1 = b[0].strip()
  b2 = b[1].strip()
  b3 = b[2].strip()
  b4 = b[3].strip()
  diheds.append((b1, b2, b3, b4, pk,  pn, phase))
  i += 1




#########




pdbbase = os.path.basename(args.pdbfilename).split('.')[0]

basefilename = "baselist.dat"
chainfilename = "ichain.dat"
if args.local:
  topfilename = "parametres.top"
else:
  topfilename = pdbbase+".top"

atoms = []
hetatms = []
with open(args.pdbfilename, "r") as pdbfile:
  for line in pdbfile:
    if line[0:4] == "ATOM":
      atoms.append(Atom(line))
    elif line[0:6] == "HETATM":
      hetatms.append(Atom(line))

# # ordered dict, because everything is index-based, and ordinary dict do not preserve the order of insertion
chains = [Chain()]
for ai, at in enumerate(atoms):
  # cut chains when chainid changes, or when the next residue's first atom is O5*, not P (i.e, it should be two chains in the first place
  if ai > 1 and (atoms[ai-1].chainid != at.chainid or (at.aty == 'O' and atoms[ai-1].aty != 'P')):
    chains.append(Chain())
  chains[-1].atlist.append(at)

for c in chains:
  for at in c.atlist:
    while len(c.reslist) < at.residx-c.atlist[0].residx+1:
      c.reslist.append(list())
    c.reslist[at.residx-c.atlist[0].residx].append(at)

residues = []
for c in chains:
  residues += c.reslist


#print [[len(r) for r in c.reslist] for c in chains]







if args.local:
  outdir=""
else:
  if not os.path.exists(args.outbase):
    os.mkdir(args.outbase)
  outdir = os.path.join(args.outbase, pdbbase)
  print "Creating dir %s" % (outdir)
  if os.path.exists(outdir):
    print "Output directory exists! Backing it up"
    bakdir = outdir+".bak"
    if os.path.exists(bakdir):
      shutil.rmtree(bakdir)
    shutil.move(outdir, bakdir)
  os.mkdir(outdir)
  shutil.copy(args.pdbfilename, os.path.join(outdir, pdbbase+".pdb"))




with open(os.path.join(outdir,basefilename), 'w') as bfn:
  for r in residues:
    rn = r[0].resname.strip()
    print >>bfn,  r[-1].idx, resl2i(rn)

with open(os.path.join(outdir,chainfilename), 'w') as cfn:
  print >>cfn, len(chains)
  for ci, c in enumerate(chains):
    print >>cfn, ci+1, len(c.atlist), '0'



aty = aty2mse.keys()
#for at in atoms:
#  if at.aty in atw and at.aty not in aty:
#    aty.append(at.aty)

#if len(aty) == 12:
#  print "Dummy residues were included (representing simplified double helices)."
#elif len(aty) != 11:
#  print "There are not 11 atom types, as should be!"

atl2 = list()
for l in atl:
  if l in aty:
    atl2.append(l)
#print aty
#print atl2






############ Find bonds, angles, etc...

bondl = list()
b1 = [(b[0], b[1]) for b in bonds]
#import pdb; pdb.set_trace()
for ri, r in enumerate(residues):
  r4i = None
  for ati, at in enumerate(r[:-1]):
    if at.aty == "R4":
      r4i = at.idx
    for at2 in r[ati+1:]:
      if (at.aty, at2.aty) in b1:
        bondl.append(((at.idx-1)*3, (at2.idx-1)*3, b1.index((at.aty, at2.aty))+1))
  if ri+1 >= len(residues):
    continue
  if r4i is not None and residues[ri+1][0].chainid == residues[ri][0].chainid:
    if residues[ri+1][0].aty == 'P':
      bondl.append(((r4i-1)*3, (residues[ri+1][0].idx-1)*3, b1.index(('R4', 'P'))+1))
    elif residues[ri+1][0].aty == 'D':
      bondl.append(((r4i-1)*3, (residues[ri+1][0].idx-1)*3, b1.index(('R4', 'D'))+1))
  elif r4i is None:
    bondl.append(((r[0].idx-1)*3, (residues[ri+1][0].idx-1)*3, b1.index(('D', residues[ri+1][0].aty))+1))



anglel = list()
b1 = [(b[0], b[1], b[2]) for b in angles]
for ri, r in enumerate(residues):
  r4i = None
  for ati, at in enumerate(r[:-2]):
    if at.aty == "R4":
      r4i = at.idx
    for at2i, at2 in enumerate(r[ati+1:-1]):
      for at3 in r[at2i+1:]:
        if (at.aty, at2.aty, at3.aty) in b1:
          anglel.append(((at.idx-1)*3, (at2.idx-1)*3, (at3.idx-1)*3, b1.index((at.aty, at2.aty, at3.aty))+1))
  if ri+1 >= len(residues):
    continue
  if r4i is not None and residues[ri+1][0].aty == 'P' and residues[ri+1][0].chainid == residues[ri][0].chainid:
    anglel.append(((r4i-2)*3, (r4i-1)*3, (residues[ri+1][0].idx-1)*3, b1.index(('C', 'R4', 'P'))+1))
    anglel.append(((r4i)*3, (r4i-1)*3, (residues[ri+1][0].idx-1)*3, b1.index(('R1', 'R4', 'P'))+1))
    anglel.append(((r4i-1)*3, (residues[ri+1][0].idx-1)*3, (residues[ri+1][1].idx-1)*3, b1.index(('R4', 'P', 'O'))+1))
  elif r4i is None:
    if residues[ri+1][0].aty == 'D' and ri+2 < len(residues) and residues[ri+2][0].aty == 'D':
      anglel.append(((r[0].idx-1)*3, (residues[ri+1][0].idx-1)*3, (residues[ri+2][0].idx-1)*3, b1.index(('D', 'D', 'D'))+1))
#      anglel.append(((r[0].idx-1)*3, (residues[ri+1][0].idx-1)*3, (residues[ri+2][0].idx-1)*3, b1.index(('D', 'D', residues[ri+2][0].aty))+1))



dihedl = list()
b0 = [(b[0], b[1], b[2], b[3]) for b in diheds]
b1 = [(b[0], b[1]) for b in anglel]
b2 = [(b[0], b[1]) for b in anglel]
for ani, an in enumerate(anglel):
  for an2 in anglel[ani+1:]:
    if an[1] == an2[0] and an[2] == an2[1]:
      if (atoms[an[0]/3].aty, atoms[an[1]/3].aty, atoms[an[2]/3].aty, atoms[an2[2]/3].aty) in b0:
        dihedl.append((an[0], an[1], an[2], an2[2], b0.index((atoms[an[0]/3].aty, atoms[an[1]/3].aty, atoms[an[2]/3].aty, atoms[an2[2]/3].aty))+1))
    if an2[1] == an[0] and an2[2] == an[1]:
      if (atoms[an2[0]/3].aty, atoms[an2[1]/3].aty, atoms[an2[2]/3].aty, atoms[an[2]/3].aty) in b0:
        dihedl.append((an2[0], an2[1], an2[2], an[2], b0.index((atoms[an2[0]/3].aty, atoms[an2[1]/3].aty, atoms[an2[2]/3].aty, atoms[an[2]/3].aty))+1))
    if an[0] == an2[1] and an[1] == an2[0]:
      if (atoms[an2[2]/3].aty, atoms[an[0]/3].aty, atoms[an[1]/3].aty, atoms[an[2]/3].aty) in b0:
        dihedl.append((an2[2], an[0], an[1], an[2], b0.index((atoms[an2[2]/3].aty, atoms[an[0]/3].aty, atoms[an[1]/3].aty, atoms[an[2]/3].aty))+1))

#for r in residues:
#  for at in r:
#    if at.aty == "R4":
#      r4i = at.idx
#  dihedl.append(((r4i)*3, (r4i+1)*3, (r4i+2)*3, (r4i+3)*3, b1.index(('R4', 'R1', 'G1', 'G2'))+1))



############







atoms += hetatms


with open(os.path.join(outdir, topfilename), 'w') as tfn:
  print >>tfn, "RNA molecule"
  
  #      read(nf,9118) NATOM,NTYPES,NBONH,MBONA,NTHETH,MTHETA,NPHIH,MPHIA,
  #     $              NHPARM,NPARM,NNB,NRES,NBONA,NTHETA,NPHIA,
  #     $          NUMBND,NUMANG,NPTRA,NATYP,NPHB,IDUM,IDUM,IDUM,IDUM,IDUM,
  #     $          IDUM,IDUM,IFBOX,NMXRS,IFCAP
  
  txt =  ("{:6d}"*8).format(len(atoms), len(aty), 0, len(bondl), 0, len(anglel), 0, len(dihedl))
  txt += ("{:6d}"*7).format(0, 0, 0, len(residues), len(bondl), len(anglel), len(dihedl))
  txt += ("{:6d}"*10).format(len(bonds), len(angles), len(diheds), len(atl2), 0, 0, 0, 0, 0, 0)
  txt += ("{:6d}"*5).format(0, 0, 0, 0, 0)
  print >>tfn, wraptxt(txt, 72)
  
  igraph = charg = mass = iac = numex = ""
  for a in atoms:
     igraph += "{:4s}".format(a.aty)
     charg += " {: 12.8E}".format(0)
     mass += " {:12s}".format(fortranfloat(atw[a.aty]))
     iac += "{:6d}".format(atl2.index(a.aty)+1)
     numex += "{:6d}".format(0)
  
  print >>tfn, wraptxt(igraph, ls=True)
  print >>tfn, wraptxt(charg)
  print >>tfn, wraptxt(mass)
  print >>tfn, wraptxt(iac, 72)
  print >>tfn, wraptxt(numex, 72)
  
  
  print >>tfn, wraptxt(("{:6d}"*len(aty)**2).format(*[0,]*len(aty)**2), 72)
  
  
  labres = ipres = ""
  for r in residues:
    labres += "{:4s}".format(resl2n(r[0].resname.strip()))
    ipres += "{:6d}".format(r[0].idx)
  
  labres = labres.replace(' ', 'i', 1)
  
  print >>tfn, wraptxt(labres, ls=True)
  print >>tfn, wraptxt(ipres, 72)
  
  
  rkl = reql = ""
  for b1, b2, rk, req in bonds:
    rkl += " {:12s}".format(fortranfloat(rk))
    reql += " {:12s}".format(fortranfloat(req))
  
  print >>tfn, wraptxt(rkl)
  print >>tfn, wraptxt(reql)
  
  
  rkl = reql = ""
  for b1, b2, b3, rk, req in angles:
    rkl += " {:12s}".format(fortranfloat(rk))
    reql += " {:12s}".format(fortranfloat(req))
  
  print >>tfn, wraptxt(rkl)
  print >>tfn, wraptxt(reql)
  
  
  pkl = pnl = phasel = ""
  for b1, b2, b3, b4, pk, pn, phase in diheds:
    pkl += " {:12s}".format(fortranfloat(pk))
    pnl += " {:12s}".format(fortranfloat(pn))
    phasel += " {:12s}".format(fortranfloat(phase))
  
  print >>tfn, wraptxt(pkl)
  print >>tfn, wraptxt(pnl)
  print >>tfn, wraptxt(phasel)
  
  
  print >>tfn, wraptxt((" {: 12.8E}"*len(aty)).format(*[0,]*len(aty)))
  
  nttyp = len(aty)*(len(aty)+1)/2
  
  print >>tfn, wraptxt((" {: 12.8E}"*nttyp).format(*[0,]*nttyp))
  print >>tfn, wraptxt((" {: 12.8E}"*nttyp).format(*[0,]*nttyp))
  
  
  print >>tfn, ""
  bondt = ""
  for b in bondl:
    bondt += ("{:6d}"*3).format(*b)
  print >>tfn, wraptxt(bondt, 72)
  
  
  print >>tfn, ""
  anglet = ""
  for a in anglel:
    anglet += ("{:6d}"*4).format(*a)
  print >>tfn, wraptxt(anglet, 72)
  
  
  print >>tfn, ""
  dihedt = ""
  for a in dihedl:
    dihedt += ("{:6d}"*5).format(*a)
  print >>tfn, wraptxt(dihedt, 72)
  
  for i in range(20):  
    print >>tfn, wraptxt(("{:6d}"*len(aty)**2).format(*[0,]*len(aty)**2), 72)
  
  print >>tfn, ""
  print >>tfn, ""
  print >>tfn, ""
  igraph = mse = acount = ""
  for ai, a in enumerate(atoms):
     igraph += "{:4s}".format(a.aty)
     mse += "{:4s}".format(aty2mse[a.aty])
     acount += "{:6d}".format(ai)
  
  print >>tfn, wraptxt(igraph, 80, ls=True)
  print >>tfn, wraptxt(mse, 80, ls=True)
  print >>tfn, wraptxt(acount, 72, ls=True)
  



