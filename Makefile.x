# Makefile for Pandora2.0
# Created: 30th October 2017
# Author:  Thorsten Wittemeier

EXE = pandora2.0

SRC_FILES_C =	xfc_arrays.c

SRC_C=$(addprefix src/, $(SRC_FILES_C))

SRC_FILES =     pg_petsc.f90          \
               	xf_arrays.f90        
 
SRC=$(addprefix src/, $(SRC_FILES))

OBJ_FILES =	pg_petsc.o            #\
#		ptf_checkarrays.o     \
#                ptf_rk3.o

OBJ=$(addprefix src/, $(OBJ_FILES))

CC = mpicc -O3
FC = mpifort -O3
 
include $(PETSC_DIR)/lib/petsc/conf/petscvariables
FFTW_INCLUDES=-I$(FFTW_INC_DIR)
PANDORA_INCLUDES=-I ./include

# Explicit dependencies

src/pg_petsc.o:         src/pg_petsc.f90

# Explicit dependencies for C interfaces
src/xf_arrays.o:	src/xf_arrays.f90            \
			src/pg_petsc.o               \
			src/xfc_arrays.o             \

$(OBJ):
	$(CC) $(SRC_C) $(PETSC_CC_INCLUDES) -c -cpp
	$(FC) $(SRC) $(PANDORA_INCLUDES) $(PETSC_FC_INCLUDES) $(FFTW_INCLUDES) -c -cpp
	mv *.o ./src

pandora:
	mkdir -p ./bin
	$(FC) $(OBJ) $(PANDORA_INCLUDES) $(PETSC_FC_INCLUDES) $(FFTW_INCLUDES) -o $(EXE) ${PETSC_LIB} -L$(FFTW_LIB_DIR) -lfftw3_mpi -lfftw3
	cp -v $(EXE) ./bin/.

configure:
	mkdir -p ~/.pandora
	cp -v ./conf/pandora.conf ~/.pandora/

clean:
	rm -fv src/*.o
	rm -fv $(EXE)
	rm -fv *~
	rm -fv *.mod

doxygen:
	doxygen doc/code.dxg
