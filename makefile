COMPILER= gfortran

FLAGS = 

LAPACK = -llapack #-lblas

EXEC = woods_saxon

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 types.o qm_solver.o
	$(COMPILER) $(FLAGS) -c $<

hamiltonian.o: hamiltonian.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

qm_solver.o: qm_solver.f90 types.o hamiltonian.o 
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o qm_solver.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
