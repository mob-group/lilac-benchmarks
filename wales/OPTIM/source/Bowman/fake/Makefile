../../libjbdummy.a : bowmanwater.o
	ar sr $@ $^

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

clean :
	rm -f *.o *.mod 
