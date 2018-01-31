#! /usr/bin/env python

import numpy
import sys

__author__ = [
        "Tristan Cranolini <tristan.cranolini@gmail.com>",
        "Jonathan Barnoud <jonathan@barnoud.net>",
]


class PDB():
    def __init__(self, coords_=None, resize=2000):
        self.atomnames = []
        self.resnames = []
        self.resnum = []
        self.chainname = ""

        self.resize = resize
        # coords_ allows us to resize the coordinates array less often
        if coords_ is None:
            self.coords_ = numpy.zeros((resize, 3))
            self.Natoms = 0
        else:
            self.coords_ = numpy.array(coords_, dtype=float).reshape(-1, 3)
            self.Natoms = self.coords_.shape[0]
        self.Bfactors = []

    def get_coords(self):
        return self.coords_[0:self.Natoms, ]

    coords = property(get_coords)

    def addatom(self, line):
        self.atomnames.append(line[12:16])
        self.resnames.append(line[17:20])
        self.resnum.append(int(line[22:26]))

        if self.Natoms >= self.coords_.shape[0]:
            self.coords_ = numpy.row_stack([self.coords_, numpy.zeros((self.resize,3))])

        self.coords_[self.Natoms,:] = [line[30:38], line[38:46], line[46:54]]
        self.Natoms += 1
        try:
            self.Bfactors.append(float(line[60:66]))
        except (ValueError):
            self.Bfactors.append(0.0)
    
    def insertatom(self, line, idx=None):
        if idx is None:
            idx = self.Natoms
        
        self.atomnames.insert(idx, line[12:16])
        self.resnames.insert(idx, line[17:20])
        self.resnum.insert(idx, int(line[22:26]))

        if self.Natoms >= self.coords_.shape[0]:
            self.coords_ = numpy.row_stack([self.coords_, numpy.zeros((self.resize,3))])

        self.coords_[idx+1:self.Natoms+1,:] = self.coords_[idx:self.Natoms,:]
        self.coords_[idx,:] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        self.Natoms += 1
        try:
            self.Bfactors.insert(idx, float(line[60:66]))
        except (ValueError):
            self.Bfactors.insert(idx, 0.0)


    def center(self):
        return sum(self.coords) / self.coords.shape[0]

    def movetocenter(self):
        self.coords -= self.center()
    
    def pdbline(self, idx):
        #return "%-6s%5d %4s%1s%3s %1s%4d%1s        %8.3f%8.3f%8.3f%6.2f%6.2f                     %2s%2s\n" \
	return	"{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(\
            "ATOM", idx+1, self.atomnames[idx], " ", self.resnames[idx],\
            self.chainname, self.resnum[idx], " ",\
            self.coords[idx][0], self.coords[idx][1], self.coords[idx][2],\
            0, self.Bfactors[idx], " ", " ")
            
    
    def writetofile(self, filename=None, filemode="w"):
        if filename is None:
            pdbfile = sys.stdout
        else:
            pdbfile = open(filename, filemode)
        
        #pdbfile.write("MODEL\n")
        for idx in range(len(self.atomnames)):
            pdbfile.write(self.pdbline(idx))
        
        #pdbfile.write("TER\nENDMDL\n")
        
        pdbfile.close()

    def _make_copy(self, filter) :
        """
        Create a filtered copy of the structure.

        Parameters :
        ------------
        - filter : an 1D array of boolean which will be used as mask. The lenght of
        this array have to fit the number of atoms of the current structure.

        Returns :
        ---------
        - copy : a filtered copy as a PDB object.
        """
        structure = PDB(numpy.array(self.coords[filter, 0:3]))
        structure.atomnames = list(numpy.array(self.atomnames)[filter])
        structure.resnames = list(numpy.array(self.resnames)[filter])
        structure.resnum = list(numpy.array(self.resnum)[filter])
        structure.chainname = self.chainname
        
        structure.resize = self.resize
        structure.Natoms = len(structure.atomnames)
        structure.Bfactors = self.Bfactors
        return structure

    def filter(self, values, exclude=False, filtprop=lambda x: x.atomnames, filtcrit=lambda x, y: x in y) :
        """
        Return a copy filtered by atom names.

        Parameters :
        ------------
        - values : values used to do the filtering
        - filtcrit : criterion for the filter
        - exclude : if False, only the atoms with the specified atom names will be
        kept. If True they will be exclude and only the other atoms will appier in
        the copy. (False by default)
        - 

        Returns :
        ---------
        - copy : a filtered copy of the structure.
        """
        filter = [filtcrit(i, values) for i in filtprop(self)]
        filter = numpy.array(filter, dtype=bool) != exclude
        return self._make_copy(filter)

    def section(self, begin, end, exclude=False) :
        """
        Return a copy filtered with by resid between two values.

        Parameters :
        ------------
        - begin, end : the first and last index of the selection.
        - exclude : if True the selection will be exclude from the copy, else only
        the selection will be kept in the copy. (False by delfault)

        Returns :
        ---------
        - copy : a filtered copy of the structure.
        """
        filter_min = (numpy.array(self.resnum) >= begin)
        filter_max = (numpy.array(self.resnum) <= end)
        filter = (filter_min == filter_max)
        filter = (filter != exclude)
        return self._make_copy(filter)

    def filter_res(self, residues, exclude=False) :
        """
        Return a copy filtered by a resid list.

        Parameters :
        ------------
        - residues : a list of resid to put in the selection
        - exclude : if True, the selection will be excluded from the copy. (False
        by default)

        Returns :
        ---------
        -copy : a filtered copy of the structure
        """
        filter = numpy.zeros(len(self.resnum))
        for index, resid in enumerate(self.resnum) :
            if resid in residues :
                filter[index] = 1
        filter = numpy.array(filter, dtype=bool)
        filter = (filter != exclude)
        return self._make_copy(filter)

    def check_integrity(self) :
        """
        Check if there is no residue missing and if the residues are well sorted.

        Returns :
        ---------
        - answer : False if something is wrong.
        """
        last_id = self.resnum[0]
        for id in self.resnum :
            if not last_id <= id <= (last_id + 1) :
                return False
            last_id = id
        return True

def readpdb(filename):
    pdbfile = open(filename, "r")
    models = []
    models.append(readmodel(pdbfile))
    for line in pdbfile:
        if line.startswith("MODEL") or line.startswith("ATOM"):
            #~ models.append(readmodel(pdbfile))
            readmodel2(pdbfile, models)
    
    return models

def isatom(line):
        return (line.startswith("ATOM") or line.startswith("HETATM"))# and line[17:20] in PDB.standard_aa_names:

def readmodel(pdbfile):
    model = PDB()
    for line in pdbfile:
        if isatom(line):# and line[17:20] in PDB.standard_aa_names:
            model.addatom(line)
        elif line.startswith("ENDMDL"):
            break
    
    model.coords = numpy.array(model.coords)
    return model

def readmodel2(pdbfile, models):
    model = PDB(numpy.empty_like(models[0].coords))
    i = 0
    for line in pdbfile:
        if isatom(line):
            model.insertatom(line, i)
            i += 1
        elif line.startswith("ENDMDL"):
            break
    
    models.append(model)


if __name__ == "__main__":
    models = readpdb(sys.argv[1])
    print models[0].coords.shape
    print len(models[0].coords)
