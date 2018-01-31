import sys
import numpy as np
import rb_to_atomistic as rb

def get_Is(config):
    Q = np.array(config).reshape(-1,3)
    Q = rb.to_atomistic(Q)

    CoM = np.average(Q,axis=0)
    Q-=CoM
    
    I = np.zeros((3,3))  # non-mass-weighted Inertia tensor
    for atom in Q:
        I[0][0]+=atom[1]*atom[1]+atom[2]*atom[2]
        I[1][1]+=atom[0]*atom[0]+atom[2]*atom[2]
        I[2][2]+=atom[1]*atom[1]+atom[0]*atom[0]
        I[0][1]-=atom[0]*atom[1]
        I[0][2]-=atom[0]*atom[2]
        I[1][2]-=atom[1]*atom[2]
    I[1][0] = I[0][1]
    I[2][0] = I[0][2]
    I[2][1] = I[1][2]

    return np.linalg.eigvalsh(I)  # The eigenvalues of I


def read_pathinfo(fname, nconfigs):
    fin = open(fname,'r')

    energies = []
    symmetries = []
    inertias = []
    thisconfig = []

    sline = fin.readline().split() # Read the first energy line
    for i in xrange(nconfigs):
        if len(sline)!=1:
            raise ValueError("Energy line appears to have multiple entries. Check nconfigs = "+str(nconfigs))
        energies.append(float(sline[0]))   # Save the energy of this configuration
        symmetries.append(fin.readline())  # Read point group information for this configuration
        sline = fin.readline().split()     # Read the first coordinate/frequency line
        while len(sline)==3:
            thisconfig.append(map(float, sline))
            sline = fin.readline().split() # Read the next coordinate/frequency line.
                                           # When this readline() encounters the start of the next configuration,
                                           # we exit the loop with sline containing the next energy line, ready for
                                           # the next cycle of the loop over i.
        Is = get_Is(thisconfig)
        inertias.append(Is)                # Save principle moments of inertia for the current configuration
        thisconfig=[]                      # and empty ready to receive the next

    fin.close()

    return energies, symmetries, inertias

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
tol = 1.0e-4
fname1 = 'path.info'
fname2 = 'expected_output/path.info'
print "Comparing ", fname1, "with ", fname2
success = True

# We are assuming only one triple per path.info file. This is a pretty safe bet, since the PATH run should only ever produce one triple.
E1, sym1, I1 = read_pathinfo(fname1, nconfigs)
E2, sym2, I2 = read_pathinfo(fname2, nconfigs)
    
if success:
    Efailures = 0
    symfailures = 0
    Ifailures = 0

    mapping = check_ordering(E1, E2, tol)

    if mapping is not None and mapping!=range(nconfigs):  # if it is None, then we will report the failure in the normal way
        print "Mapping stationary points between the files according to: ", mapping
        E2 = [ E2[i] for i in mapping]
        sym2 = [ sym2[i] for i in mapping]
        I2 = [ I2[i] for i in mapping]       

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

    for i in xrange(nconfigs):
        if (abs((I1[i][0]-I2[i][0])/I1[i][0])>tol) or (abs((I1[i][1]-I2[i][1])/I1[i][1])>tol) or (abs((I1[i][2]-I2[i][2])/I1[i][2])>tol):
            Ifailures+=1
            print "Different moments of inertia for structure "+str(i)+':'
            print I1[i]
            print I2[i]
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
    print "Number of stationary points where the moments of inertia differed between the files:", Ifailures
    sys.exit(1)
