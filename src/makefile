FC=/usr/bin/gfortran
FLAGS= -fopenmp -O2 -ffree-form
WARN= -Wall -pedantic
DBG= -g -fcheck=all -fbacktrace
COMPILE=$(FC) $(FLAGS) $(WARN) $(DBG)

BAND=grid_bandpass
OBJS=integrate_bandpass.o extinction.o blackbody.o const.o

$(BAND): $(OBJS) grid_bandpass.o
	$(COMPILE) -o ../$(BAND) $^

%.o: %.f90
	$(COMPILE) -c $<

grid_bandpass.o: integrate_bandpass.o
integrate_bandpass.o: extinction.o blackbody.o const.o
blackbody.o: const.o

.PHONY: clean check

check: $(TEST)
	./$(TEST)

clean:
	rm -f *.o *.mod $(BAND)
