import sys

fractol = 1.0e-4

fin1 = open('expected_output/min.data.info','r')
fin2 = open('min.data.info','r')

print "Comparing expected_output/min.data.info and min.data.info"

head1 = fin1.readline()
head2 = fin2.readline()
    
fin1.close()
fin2.close()

E1 = float(head1.split()[0])
E2 = float(head2.split()[0])
S1 = int(head1.split()[2])
S2 = int(head2.split()[2])
Inertia1 = map(float,head1.split()[3:])
Inertia2 = map(float,head2.split()[3:])

success = True

if abs(E1-E2)>fractol*abs(E1):
    print "Energies of the two minima are different: ", E1, E2
    success = False

if S1!=S2:
    print "The minima have different point group orders:", S1, S2
    success = False

for i in xrange(3):
    if abs(Inertia1[i]-Inertia2[i])>fractol*Inertia1[i]:
        print "Principle moment of inertia ", i, "is different for the two minima"
        print Inertia1[i], Inertia2[i]
        success = False

if success:
    print "All tests completed successfully: minima are the same"
    sys.exit(0)
else:
    print "One or more tests failed. Minima are different"
    sys.exit(1)
