FC=gfortran
FLAGS=-O3 -ffree-form -fopenmp
WARN= -Wall -pedantic
DBG= #-g -fcheck=all -fbacktrace
COMPILE=$(FC) $(WARN) $(DBG) $(FLAGS)

BAND=grid_bandpass
OBJS=integrate_bandpass.o extinction.o
CCM=test_CCM
TEST=test_bandpass

$(BAND): $(OBJS) grid_bandpass.o
	$(COMPILE) -o ../$(BAND) $^

$(CCM): test_CCM.o extinction.o
	$(COMPILE) -o $@ $^

$(TEST): $(OBJS) test_bandpass.o
	$(COMPILE) -o ../$(TEST) $^

test_CCM.o: test_CCM.f
	$(COMPILE) -c $<

%.o: %.f
	$(COMPILE) -c $<

test_bandpass.o: integrate_bandpass.o
grid_bandpass.o: integrate_bandpass.o
integrate_bandpass.o: extinction.o

.PHONY: clean check

check: $(TEST)
	./$(TEST)

clean:
	rm -f *.o *.mod $(BAND) $(TEST) $(CCM)
