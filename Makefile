# makefile for bem program
SRC = ./src
OBJ = ./obj
BIN = ./bin
EXT = ./ext
EXC = tryHDF_fortran
#BLIBS = -L/usr/lib/x86_64-linux-gnu/ -llapack -lblas
BLIBS = -L//usr/local/opt/lapack/lib -llapack -lblas

clean: mdl_main
	rm *.o

# -----------------------------------------------
mdl_main : mdl_main.o solver.o
	h5fc -o $(BIN)/$(EXC) -I$(OBJ) mdl_main.o gaussquad.o model.o mesh.o bem_chunk.o \
	element.o material.o h5interface.o BEM_elasticity.o discretization.o \
	funda_elasticity.o element_computations.o element_calc.o \
	comp_model.o matOps.o solver.o freefield.o $(BLIBS) -fcheck=all
	
mdl_main.o : model.o gaussquad.o comp_model.o
	h5fc  -J$(OBJ) -c $(SRC)/mdl_main.f08 -I$(OBJ) 

model.o : bem_chunk.o mesh.o material.o
	h5fc -J$(OBJ) -c $(SRC)/model.f08 -I$(OBJ) 

bem_chunk.o : h5interface.o 
	h5fc -J$(OBJ) -c $(SRC)/bem_chunk.f08

mesh.o : h5interface.o element.o
	h5fc -J$(OBJ) -c $(SRC)/mesh.f08 

element.o : h5interface.o 
	gfortran -J$(OBJ) -c $(SRC)/element.f08 

material.o :
	h5fc -J$(OBJ) -c  $(SRC)/material.f08 

h5interface.o : 
	h5fc -J$(OBJ) -c $(SRC)/h5interface.f08 -I$(OBJ) 

# ---------------------------------------------------------
#
comp_model.o : BEM_elasticity.o model.o bem_chunk.o matOps.o solver.o freefield.o
	h5fc -J$(OBJ) -c $(SRC)/comp_model.f08 

matOps.o : 
	h5fc -J$(OBJ) -c $(SRC)/matOps.f08 

# BEM element_computations
BEM_elasticity.o: model.o element_computations.o funda_elasticity.o element_calc.o gaussquad.o discretization.o
	gfortran -J$(OBJ) -c  $(SRC)/BEM_elasticity.f08 

discretization.o: element_computations.o
	gfortran -J$(OBJ) -c  $(SRC)/discretization.f08 $(BLIBS) -fcheck=all

element_calc.o: element_computations.o
	gfortran -J$(OBJ) -c  $(SRC)/element_calc.f08 $(BLIBS)

funda_elasticity.o: 
	gfortran -J$(OBJ) -c $(SRC)/funda_elasticity.f08

freefield.o: element_computations.o
	gfortran -J$(OBJ) -c $(SRC)/freefield.f08 $(BLIBS)

solver.o:
	gfortran -J$(OBJ) -c  $(SRC)/solver.f08 $(BLIBS)

# ----------------------------------------------------------
# external subroutines
element_computations.o:
	gfortran -J$(OBJ) -c  $(EXT)/element_computations.f08

gaussquad.o:
	gfortran -J$(OBJ) -c $(EXT)/gaussquad.f90 $(BLIBS)
	
# --------------------------------------------------------
# run the code
	
run : py_run.o
	$(BIN)/$(EXC)

py_run.o:
	python3 ./py/pyMain.py

clear:
	rm ./*.o
	rm $(OBJ)/*.mod
	rm $(BIN)/$(EXC)
