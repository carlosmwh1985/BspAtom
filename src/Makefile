# makefile to build the excutable for Bsp_Atom
# by Carlos Granados
# USAGE:
# make all	      # compile all binary
# make clean      # remove ALL binaries and objects

.PHONY = all clean

#.SUFFIXES: .f90 .f .F .F90 .o
#%.o: %.f90
#	$(FC) $(FFLAGS) -c $<
#.f.o:
#	$(FC77) $(FFLAGS) -c $<

#----------------------------------------------------
#     COMPILERS AND LIBRARIES
#----------------------------------------------------
FC=ifort
FC77=ifort
CC=icc
FFLAGS = -O3
IFLAGS = 
LIBS = -lpthread -lm -mkl

#----------------------------------------------------
#     SOURCE FILES
#----------------------------------------------------
SRC1 = Modules.f90
SRC2 = Bsp_Atom.f90 ReadInputs.f90 matrices.f90 PhotoIon.f90 WriteWF.f90
SRC3 = grid.f90 CubicSpline.f90 bsplvb.f90 interv.f90
SRC4 = Ang_Ints.f90 Ang_Ints_Aux.f90 TorusFuns.f90 TorusFunsInts.f90
SRC5 = Funs_AssLegendre.f90 Funs_AssLaguerre.f90 Funs_SphHarms.f90 Funs_Bessel.f90 Funs_WignerSymbols.for
SRC = SRC1 + SRC2 + SCR3 + SCR4 + SRC5

#----------------------------------------------------
#     OBJECT FILES
#----------------------------------------------------
#OBJ = Bsp_Atom.o
OBJ = Modules.o\
	Bsp_Atom.o ReadInputs.o matrices.o PhotoIon.o WriteWF.o\
	grid.o CubicSpline.o bsplvb.o interv.o\
	Ang_Ints.o Ang_Ints_Aux.o TorusFuns.o TorusFunsInts.o\
	Funs_AssLegendre.o Funs_AssLaguerre.o Funs_SphHarms.o Funs_Bessel.o Funs_WignerSymbols.o

#----------------------------------------------------
#     COMPILATION GROUPING
#----------------------------------------------------
OMP = Bsp_Atom_omp.x

#----------------------------------------------------
#     COMPILE PROGRAM
#----------------------------------------------------
all:	exec

exec: $(OBJ) 
	$(FC) -o $(OMP) $(OBJ) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.for
	$(FC77) $(FFLAGS) -c $<

clean:$
	rm -rf *.o *.mod *.x

