FC= /home/squar/MPI/OpenMPI/bin/mpif90 OMPI_FC=flang
#FC=flang
#FC= mpif90
#FCFLAGS = -fbounds-check -Werror -fopenmp
FCFLAGS = -fbounds-check -fopenmp
#FCFLAGS = -fopenmp -S -emit-llvm

PROGRAMS = mainSMP.x mainFence.x mainLock.x mainShared.x mainSharedNoRMAInnerLocks.x mainSharedNoRMAOuterLocks.x
DEPENDENCIES = dataTypesMod.o printsMod.o computationMod.o benchmarkUtilityMod.o 

SRC = src/

# "make" builds all
all: $(PROGRAMS)

debug:  FCFLAGS += -g
debug: $(PROGRAMS)

mainSMP.o: $(DEPENDENCIES) 
mainSMP.x: $(DEPENDENCIES)

mainFence.o: $(DEPENDENCIES)
mainFence.x: $(DEPENDENCIES)

mainLock.o: $(DEPENDENCIES)
mainLock.x: $(DEPENDENCIES)  

mainShared.o: $(DEPENDENCIES)
mainShared.x: $(DEPENDENCIES)  

mainSharedNoRMAInnerLocks.o: $(DEPENDENCIES)
mainSharedNoRMAInnerLocks.x: $(DEPENDENCIES)  

mainSharedNoRMAOuterLocks.o: $(DEPENDENCIES)
mainSharedNoRMAOuterLocks.x: $(DEPENDENCIES)  


%.x: %.o
	$(FC) $(FCFLAGS) $(OMFLAGS) -o $@ $^ $(LDFLAGS)

%.o: $(SRC)%.f90
	 $(FC) $(FCFLAGS) $(OMFLAGS) -c $<


.PHONY: clean cleanall

clean:
	 rm -f *.o *.mod *.MOD *.ll

cleanall: clean
	 rm -f $(PROGRAMS)

print-%  : ; @echo $* = $($*)
