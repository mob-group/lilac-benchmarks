# !/usr/bin/env python
# -*- coding: UTF8 -*-

purines = ["A", "G"]
pyrimidines = ["C", "U", "T"]
bases = purines+pyrimidines
lenbases_AA = {'A': 33, 'G': 34, 'C': 31, 'U': 30}


phosphate = ["P", "OP1", "OP2"]
sugar = ["C5'", "C4'", "C3'", "C2'", "C1'", "O4'", "O3'", "O2'"]

heavy = dict()
heavy["A1"] = [" N7 ", " N9 ", " C4 ", " C5 ", " C8 "]
heavy["A2"] = [" N1 ", " N3 ", " N6 ", " C2 ", " C4 ", " C5 ", " C6 "]
heavy["A"] = heavy["A1"] + heavy["A2"]
heavy["C1"] = [" N1 ", " N3 ", " N4 ", " C2 ", " C4 ", " C5 ", " C6 ", " O2 "]
heavy["G1"] = heavy["A1"]
heavy["G2"] = [" N1 ", " N2 ", " N3 ", " C2 ", " C4 ", " C5 ", " C6 ", " O6"]
heavy["G"] = heavy["G1"] + heavy["G2"]
heavy["T1"] = [" N1 ", " N3 ", " C2 ", " C4 ", " C5 ", " C6 ", " O2 ", " O4 "]
heavy["U1"] = [" N1 ", " N3 ", " C2 ", " C4 ", " C5 ", " C6 ", " O2 ", " O4 "]


FA2CG_backbone = [" P  ", " O5'", " C5'", " C4'", " C1'"]
FA_backbone = FA2CG_backbone + ["O1P", "O2P"]
CG_backbone = [" P  ", " O5* ", " C5*", " CA", " CY"]
