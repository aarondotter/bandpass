FC=gfortran
FLAGS=-O2 -ffree-form -fopenmp
WARN= -Wall -pedantic
DBG= -g -fcheck=all -fbacktrace
COMPILE=$(FC) $(WARN) $(DBG) $(FLAGS)

BAND=grid_bandpass
OBJS=integrate_bandpass.o extinction.o
CCM=test_CCM

$(BAND): $(OBJS) grid_bandpass.o
	$(COMPILE) -o ../$(BAND) $^

$(CCM): test_CCM.o extinction.o
	$(COMPILE) -o $@ $^

test_CCM.o: test_CCM.f
	$(COMPILE) -c $<

%.o: %.f
	$(COMPILE) -c $<

grid_bandpass.o: integrate_bandpass.o
integrate_bandpass.o: extinction.o

.PHONY: clean check

check: $(TEST)
	./$(TEST)

clean:
	rm -f *.o *.mod $(BAND) $(CCM)
