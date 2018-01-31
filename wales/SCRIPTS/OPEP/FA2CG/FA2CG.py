# !/usr/bin/env python
# -*- coding: UTF8 -*-

import sys

import numpy as np

import pdbparser
import RNA

def filtbase(base, type):
	B = base.filter(RNA.heavy[type])
	# mean along first coordinate
	B.coords[0] = B.coords.mean(0)
	B.atomnames[0] = " "+type+" "
	return B.pdbline(0)

def Base_(FA_pdb, resid):
#	import pdb; pdb.set_trace()
	base = FA_pdb.filter_res([resid])
	lines = []
	for b in RNA.bases:
		if b in base.resnames[0]:
			if b in RNA.purines:
				lines.append(filtbase(base, b+"2"))
			lines.append(filtbase(base, b+"1"))
			break
	return lines

def Base(FA_pdb):
	for i in xrange(Natoms):
		if "O5'" in FA_pdb.atomnames[i]:
			Base_(FA_pdb, i)

def Coarse(FA_pdb):
	CG_backbone = FA_pdb.filter(RNA.FA2CG_backbone)
	for i in xrange(CG_backbone.Natoms):
		CG_backbone.atomnames[i] = CG_backbone.atomnames[i].replace("'", "*")
	CG_full = CG_backbone._make_copy(np.repeat(True, CG_backbone.Natoms))
	for i in range(CG_backbone.Natoms-1, -1, -1):
		if "C1*" in CG_backbone.atomnames[i]:
			for l in Base_(FA_pdb, CG_backbone.resnum[i]):
				CG_full.insertatom(l, i+1)
	
	return CG_full

def main(filename, out):
	#--------------------- Main Routine ------------------------#
	
	famodels = pdbparser.readpdb(filename)
	CG = Coarse(famodels[0])
	CG.writetofile(out)


if __name__ == "__main__":
	if len(sys.argv) < 2:
		fname = "6TNA_fit.pdb"
	else:
		fname = sys.argv[1]
	if len(sys.argv) < 3:
		outname = None
	else:
		outname = sys.argv[2]
	
	main(fname, outname)
