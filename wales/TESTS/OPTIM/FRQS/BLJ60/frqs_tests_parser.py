import sys

def compare_frqs_mindatainfo(f1, f2, tol=1.0e-3, fractional=False):
    fin1 = open(f1,'r')
    fin2 = open(f2,'r')

    logp1 = float(fin1.readline().split()[1])
    logp2 = float(fin2.readline().split()[1])

    if fractional:
        tol*=abs(logp1)  # If logp1 is negative, we are in trouble anyway...

    success = abs(logp1-logp2)<tol
    if not success:
        print "Log product frequencies differ:", logp1, logp2

    return success

def read_pathinfo(fname, nconfigs, n_frq_lines, n_coord_lines):
    fin = open(fname,'r')

    energies = []
    symmetries = []
    frqs = []
    thisconfig = []

    for i in xrange(nconfigs):
        energies.append(float(fin.readline().split()[0])) # Read the energy line
        symmetries.append(fin.readline()) # Read point group information for this configuration
        for j in xrange(n_frq_lines):
            thisconfig.append(map(float, fin.readline().split()))
        for j in xrange(n_coord_lines):
            fin.readline()
        frqs.append(thisconfig)
        thisconfig=[]

    fin.close()

    return energies, symmetries, frqs

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


def compare_frqs_pathinfo(f1, f2, nconfigs, n_frqs_lines, n_coords_lines, tol=1.0e-3):  # tol is fractional.

    print "Comparing ", f1, "with ", f2
    success = True

    E1, sym1, frqs1 = read_pathinfo(f1, nconfigs, n_frqs_lines, n_coords_lines)
    E2, sym2, frqs2 = read_pathinfo(f2, nconfigs, n_frqs_lines, n_coords_lines)
    

    Efailures = 0
    symfailures = 0
    frqfailures = 0
    nzeros = 6

    mapping = check_ordering(E1, E2, tol)

    if mapping is None:
        success=False
    elif mapping!=range(nconfigs):  # if it is None, then we will report the failure in the normal way
        print "Mapping stationary points between the files according to: ", mapping
        E2 = [ E2[i] for i in mapping]
        sym2 = [ sym2[i] for i in mapping]
        frqs2 = [ frqs2[i] for i in mapping]       

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

    if success:
        for i in xrange(nconfigs):

            # Set up a list of which frequencies are zeros - these need to be tested with a different tolerance, since
            # there is no guarantee they will always have the same value
            zeros = [[False, False, False] for x in xrange(n_frqs_lines)]
            zeros[-1][:]=[True, True, True]
            if i%3==1:  # Transition state
                zeros[-2][2]=True
                zeros[-1][2]=False
            smallest = frqs2[i][-2][1]  # Save the smallest positive eigenvalue to test the zeros

            for j in xrange(n_frqs_lines):
                nfails = 0
                for k in xrange(3):
                    
                    if (zeros[j][k]):
                        # Zero frequencies must be no larger than 10% of the smallest non-zero frequency (this is a very lenient criterion!)   
                        if abs(frqs2[i][j][k]/smallest) > 0.1:
                            print "zero frequency ", frqs2[i][j][k], "is too large compared with smallest positive frequency ", smallest
                            nfails += 1
                    else:
                        # Non-zero frequency is tested with a fractional tolerance against expected output.
                        if (abs((frqs1[i][j][k]-frqs2[i][j][k])/frqs1[i][j][k]) > tol):  
                            nfails += 1
                if nfails:
                    print "Config ", i, "Line ", j+1, "failed ", nfails, "times"
                    frqfailures += nfails
                    success = False

        if frqfailures:
            print frqfailures, "frequencies differed between the two files"
        return success
    else:
        print "Stationary points in path.info file differ. Test failed."
        return False


def run_tests(nconfigs, file_len, natoms):
    nfreq_lines = file_len/nconfigs - 2 - natoms  # Divide the file into nconfigs chunks. Each has 2 header lines and natoms coordinate lines. All remaining lines are for frequencies.

    print "*"*30
    print "GEOPT tests"
    print "*"*30
    geopt_success = compare_frqs_mindatainfo('GEOPT/expected_output/min.data.info', 'GEOPT/min.data.info')
    # NOTE: DUMPVECTORS is not tested currently. This might be a sensible thing to add in future.
    if geopt_success:
        print "Completed successfully"

    print "\n"
    print "*"*30
    print "PATH tests"
    print "*"*30
    path_success = compare_frqs_pathinfo('PATH/expected_output/path.info', 'PATH/path.info', 3, nfreq_lines, natoms)
    if path_success:
        print "Completed successfully"

    print "\n"
    print "*"*30
    print "DOUBLE_ENDED tests"
    print "*"*30
    connect_success = compare_frqs_pathinfo('DOUBLE_ENDED/expected_output/path.info', 'DOUBLE_ENDED/path.info', nconfigs, nfreq_lines, natoms)
    if connect_success:
        print "Completed successfully"

    print "\n"
    print "*"*30
    print "FINAL SUMMARY"
    print "*"*30

    print "GEOPT frequency calculations. Success =", geopt_success
    print "PATH frequency calculations. Success =", path_success
    print "DOUBLE_ENDED frequency calculations. Success =", connect_success

    if (geopt_success and path_success and connect_success):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__=='__main__':

    nTS = int(sys.argv[1])
    nconfigs = nTS*3

    file_len = int(sys.argv[2])
    n_coord_lines = int(sys.argv[3])

    run_tests(nconfigs, file_len, n_coord_lines)
