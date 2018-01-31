#! /usr/bin/env python

import sys
sys.path.append('/usr/local/lib64/python2.4/site-packages')
from numpy import *

def aaxis_to_mat(p):
    """Converts an angle-axis rotation into a rotation matrix"""
    # cf wikipedia page for Rodrigues's rotation formula
    theta = linalg.norm(p)
    if theta == 0:
        return identity(3)
    else:
        k = p / theta
        kx = matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        return identity(3) * cos(theta) + kx * sin(theta) + (1 - cos(theta)) * outer(k, k)

def getFilename():
    if len(sys.argv) != 2:  # the program name and one arguments
      # stop the program and print an error message
      sys.exit("Usage: (please enter full filename):  pysitesToXYZ.py pysitesFilename")
    return sys.argv[1]

def readPysites(filename):
    #open file
    f = open(filename, 'r')
    inp=f.readlines()
    f.close()
    return inp[2:]

def writeXYZ(filename,pysites):
    pi=3.1415926535897
    coords=['0.0 0.0 0.0','0.0 0.0 '+str(2*pi)]
    f=open(filename,'w')
    numBB=len(coords)/2
    numAtoms=len(pysites)*numBB
    f.write(str(numAtoms)+'\n')
    f.write('Energy of minimum     0=   999.9 first found at step     0'+'\n')
    for (BBPos,BBP) in zip(coords[0:numBB],coords[numBB:2*numBB]):
       BBPosSp=BBPos.split()
       COG=array([float(BBPosSp[0]), float(BBPosSp[1]),float(BBPosSp[2])])
       BBPSp=BBP.split()
       BBP=array([float(BBPSp[0]),float(BBPSp[1]),float(BBPSp[2])])
       Rmat=aaxis_to_mat(BBP)

       for site in pysites:
           siteData=site.split()
           POS=COG+dot(Rmat,array([float(siteData[1]),float(siteData[2]),float(siteData[3])]))
           POS=array(POS.tolist()[0])
           ABC=[float(siteData[7]),float(siteData[8]),float(siteData[9])]
           P=array(dot(Rmat,array([float(siteData[14]),float(siteData[15]),float(siteData[16])])))
           P=array(P.tolist()[0])
           R=aaxis_to_mat(P)
           x=str(POS[0])
           y=str(POS[1])
           z=str(POS[2])
           a=str(ABC[0])
           b=str(ABC[1])
           c=str(ABC[2])
           a11=str(R[0,0])
           a12=str(R[0,1])
           a13=str(R[0,2])
           a21=str(R[1,0])
           a22=str(R[1,1])
           a23=str(R[1,2])
           a31=str(R[2,0])
           a32=str(R[2,1])
           a33=str(R[2,2])
           PS1=str(P[0])
           PS2=str(P[1])
           PS3=str(P[2])
           outStr='O '+x+' '+y+' '+z+' ellipse '+a+' '+b+' '+c+' '+a11+' '+a12+' '+a13+' '+a21+' '+a22+' '+a23+' '+a31+' '+a32+' '+a33+' atom_vector '+PS1+' '+PS2+' '+PS3+'\n'
           f.write(outStr)
    return

if __name__ == '__main__':
    filename=getFilename()
    pysites = readPysites(filename)
    writeXYZ(filename+'.xyz',pysites)
