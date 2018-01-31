import sys
import numpy as np

def read_pathinfo(fname, nconfigs, n_coord_lines):
    fin = open(fname,'r')

    energies = []
    symmetries = []
    thisconfig = []

    sline = fin.readline().split() # Read the first energy line
    for i in xrange(nconfigs):
        if len(sline)!=1:
            raise ValueError("Energy line appears to have multiple entries. Check nconfigs = "+str(nconfigs))
        energies.append(float(sline[0]))   # Save the energy of this configuration
        symmetries.append(fin.readline())  # Read point group information for this configuration
        sline = fin.readline().split()     # Read the first coordinate/frequency line
        for i in xrange(n_coord_lines):
            thisconfig.append(map(float, sline))
            sline = fin.readline().split() # Read the next coordinate/frequency line.
                                           # When this readline() encounters the start of the next configuration,
                                           # we exit the loop with sline containing the next energy line, ready for
                                           # the next cycle of the loop over i.
        thisconfig=[]                      # and empty ready to receive the next

    fin.close()

    return energies, symmetries

def check_ordering(E1, E2, tol):
    nconfigs = len(E1)
    mapping = range(nconfigs)
    for i in xrange(nconfigs):
        if mapping[i]==i and abs((E1[i]-E2[i])/E1[i]) > tol:
            for j in xrange(i+1, nconfigs):
                if abs((E1[i]-E2[j])/E1[i]) < tol:
                    mapping[i] = j
                    mapping[j] = i
                    break
            if mapping[i]==i:
                print "Failed to find a stationary point with the same energy as ", i
                print "Mapping points together is not possible."
                return None  # Failure: there is no configuration which matches configuration i.
    return mapping

nconfigs = 3*int(sys.argv[1])
ncoords = int(sys.argv[2])

tol = 1.0e-4
fname1 = 'path.info'
fname2 = 'expected_output/path.info'

n_coord_lines = int(np.floor(ncoords/3.0)) + ncoords%3
print "Expecting ", n_coord_lines, "lines of coordinates per configuration"

print "Comparing ", fname1, "with ", fname2
success = True

# We are assuming only one triple per path.info file. This is a pretty safe bet, since the PATH run should only ever produce one triple.
E1, sym1 = read_pathinfo(fname1, nconfigs, n_coord_lines)
E2, sym2 = read_pathinfo(fname2, nconfigs, n_coord_lines)
    
if success:
    Efailures = 0
    symfailures = 0

    mapping = check_ordering(E1, E2, tol)

    if mapping is not None and mapping!=range(nconfigs):  # if it is None, then we will report the failure in the normal way
        print "Mapping stationary points between the files according to: ", mapping
        E2 = [ E2[i] for i in mapping]
        sym2 = [ sym2[i] for i in mapping]     

    for i in xrange(nconfigs):
        if abs((E1[i]-E2[i])/E1[i]) > tol:
            Efailures +=1
            print "Different energies for structure "+str(i)+':', E1[i], E2[i]
            success = False

    for i in xrange(nconfigs):
        if sym1[i]!=sym2[i]:
            symfailures+=1
            print "Symmetries differ for structure "+str(i)+':', sym1[i][:-2], sym2[i][:-2]
            success = False

else:  # If the first safety checks failed
    print "path.info files differed dramatically: comparison not possible."
    sys.exit(1)


if success:
    print "All tests completed successfully."
    sys.exit(0)
else:
    print "Number of stationary point energies which differed between the files:", Efailures
    print "Number of stationary point symmetries which differed between the files:", symfailures
    sys.exit(1)
