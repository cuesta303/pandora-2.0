# Makefile for Pandora2.0
# Created: 30th October 2017
# Author:  Thorsten Wittemeier

EXE = pandora2.0

SRC_FILES =	pg_petsc.f90          \
		pg_parameters.f90     \
		pg_constants.f90      \
		pg_headers.f90        \
		pg_files.f90          \
		pg_domain.f90         \
		pg_random.f90         \
		pf_arrays.f90         \
		pp_arrays.f90         \
		pp_mpi.f90            \
		pp_interp.f90         \
		pp_rk3.f90            \
		pp_control.f90        \
		pg_restart.f90        \
		pf_fftw.f90           \
		pf_force.f90          \
		pf_fluidstats.f90     \
		pp_stats.f90          \
		pf_initialise.f90     \
		pf_rk3.f90            \
		pg_rk3.f90            \
		ptg_checkvalues.f90   \
		ptg_checkrandom.f90   \
		ptg_checkdomain.f90   \
		ptf_checkarrays.f90   \
		ptf_rk3.f90           \
		pg_control.f90        \
		pa_pandora.f90        

SRC=$(addprefix src/, $(SRC_FILES))

OBJ_FILES =	pa_pandora.o          \
		pg_petsc.o            \
		pg_control.o          \
		pg_headers.o          \
		pg_parameters.o       \
		pg_constants.o        \
		pg_files.o            \
		pg_domain.o           \
		pg_random.o           \
		pf_arrays.o           \
		pp_arrays.o           \
		pp_mpi.o              \
		pp_interp.o           \
		pp_rk3.o              \
		pp_control.o          \
		pg_restart.o          \
		pf_fftw.o             \
		pf_force.o            \
		pf_fluidstats.o       \
		pp_stats.o            \
		pf_initialise.o       \
		pf_rk3.o              \
		pg_rk3.o              \
		ptg_checkvalues.o     \
		ptg_checkrandom.o     \
		ptg_checkdomain.o     \
		ptf_checkarrays.o     \
                ptf_rk3.o

OBJ=$(addprefix src/, $(OBJ_FILES))

FC = mpifort -O3
 
include $(PETSC_DIR)/lib/petsc/conf/petscvariables
FFTW_INCLUDES=-I$(FFTW_INC_DIR)
PANDORA_INCLUDES=-I ./include

# Explicit dependencies

src/pa_pandora.o:	src/pa_pandora.f90           \
			src/pg_control.o

src/pg_control.o:	src/pg_control.f90           \
			src/pg_headers.o             \
			src/pg_parameters.o          \
			src/pg_constants.o           \
			src/pg_files.o               \
			src/pg_domain.o              \
			src/pf_arrays.o              \
			src/pf_initialise.o          \
			src/pp_control.o             \
			src/pg_rk3.o                 \
			src/pg_restart.o             \
			src/ptg_checkvalues.o        \
			src/ptg_checkrandom.o        \
		        src/ptg_checkdomain.o        \
			src/ptf_checkarrays.o        \
			src/ptf_rk3.o 

src/pg_headers.o:	src/pg_headers.f90           

src/pg_parameters.o:	src/pg_parameters.f90

src/pg_constants.o:	src/pg_constants.f90         \
			src/pg_parameters.o

src/pg_files.o:		src/pg_files.f90             \
			src/pg_parameters.o

src/pg_domain.o:	src/pg_domain.f90            \
			src/pg_parameters.o          \
			src/pg_constants.o 

src/pg_random.o:	src/pg_random.f90            \
			src/pg_constants.o 

src/pf_arrays.o:	src/pf_arrays.f90            \
			src/pg_domain.o              \
			src/pg_parameters.o          \
			src/pg_constants.o 

src/pg_restart.o:	src/pg_restart.f90           \
			src/pf_arrays.o              \
			src/pp_arrays.o

src/pp_mpi.o:		src/pp_mpi.f90               \
			src/pf_arrays.o              \
			src/pg_parameters.o 

src/pp_interp.o:	src/pp_interp.f90            \
			src/pg_parameters.o 

src/pp_rk3.o:		src/pp_rk3.f90               \
			src/pg_constants.o 

src/pp_control.o:	src/pp_control.f90           \
			src/pg_parameters.o          \
			src/pp_arrays.o              \
			src/pp_mpi.o                 \
			src/pp_interp.o              \
			src/pp_rk3.o   

src/pf_force.o:		src/pf_force.f90             \
			src/pg_parameters.o 

src/pf_initialise.o:	src/pf_initialise.f90        \
			src/pg_parameters.o          \
			src/pg_domain.o              \
			src/pg_constants.o           \
			src/pf_arrays.o             

src/pf_rk3.o:		src/pf_rk3.f90               \
			src/pg_parameters.o          \
			src/pg_domain.o              \
			src/pf_arrays.o              \
			src/pf_fftw.o                \
			src/pf_fluidstats.o          \
			src/pf_force.o               \
			src/pg_constants.o 

src/pf_fluidstats.o:	src/pf_fluidstats.f90        \
			src/pg_parameters.o          \
			src/pf_arrays.o              \
			src/pg_constants.o 

src/pg_rk3.o:		src/pg_rk3.f90               \
			src/pg_parameters.o          \
			src/pg_constants.o 

src/pf_fftw.o:		src/pf_fftw.f90              \
			src/pg_domain.o              \
			src/pg_parameters.o          \
			src/pg_constants.o 

src/ptg_checkvalues.o:	src/ptg_checkvalues.f90      \
			src/pg_parameters.o

src/ptg_checkdomain.o:	src/ptg_checkdomain.f90      \
			src/pg_parameters.o          \
			src/pg_domain.o

src/ptg_checkrandom.o:	src/ptg_checkrandom.f90      \
			src/pg_parameters.o          \
			src/pg_constants.o 

src/ptf_checkarrays.o:	src/ptf_checkarrays.f90      \
			src/pf_arrays.o

src/ptf_rk3.o:		src/ptf_rk3.f90              \
			src/pf_rk3.o                 \
			src/pg_parameters.o          \
			src/pg_domain.o              \
			src/pf_arrays.o              \
			src/pg_constants.o 

$(OBJ):
	$(FC) $(SRC) $(PANDORA_INCLUDES) $(PETSC_FC_INCLUDES) $(FFTW_INCLUDES) -c -cpp
	mv *.o ./src

pandora:
	mkdir -p ./bin
	$(FC) $(OBJ) $(PANDORA_INCLUDES) $(PETSC_FC_INCLUDES) $(FFTW_INCLUDES) -o $(EXE) ${PETSC_LIB} -lfftw3_mpi -lfftw3
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
