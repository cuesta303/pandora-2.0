!---------------------------------------------------------------------------------
! module tg_checkvalues
!> Contains tests to ensure that the program control parameters are known on
!! all processes.
module tg_checkvalues
 
!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

!------- data section begins ----------------------
implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

public :: tg_CheckValuesParameters

!---------- data section ends ----------------------
contains


!---------------------------------------------------------------------------------
! subroutine tg_CheckValuesParameters 
!> Write the configuration parameters on all processes to standard out to make sure
!! they are identical.
!> @param ierr should return 0
subroutine tg_CheckValuesParameters ( ierr )

  use g_parameters, only : MASTER, MYRANK, NPROCS, &
                           V_ALLOCS, V_DEALLOCS, V_MPI, V_DFT, V_GENERAL, V_DEBUG, V_LOG, V_TEST, V_64BIT, &
                           RUNMODE, TSTEPS, TS_REPORT, TS_FLUIDVTK, TS_CORR, TS_SPEC3D, TS_SPEC1D, TS_PLANE, F_MOMENTS, &
                           NL_BUFFER, NODES_X, NODES_Y, NODES_Z, B110, B220, B330, &
                           GRAVITY, &
                           RANDOM_SEED_SIZE, RANDOM_SEED_FLUID, RANDOM_SEED_PARTICLES, &
                           F_RHO, F_NU, &
                           F_TYPE, F_DEFORM_A, F_DEFORM_B, F_DEFORM_C, F_DEFORM_S, F_REMESH_T, &
                           F_INIT_TYPE, F_INIT_U0, F_INIT_KP, F_INIT_KV, F_INIT_KSTART, F_INIT_KEND, & 
                           F_INIT_A, F_INIT_B, F_INIT_C, F_INIT_I, F_INIT_J, F_INIT_K, &
                           FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP, &
                           F_RESTART_WRITE, F_RESTART_FILE_READ, F_RESTART_FILE_WRITE, F_RESTART_SCALE, &
                           F_FORCE_FILE_READ, F_FORCE_FILE_WRITE, &
                           P_RESTART_WRITE, P_RESTART_FILE_READ, P_RESTART_FILE_WRITE, &
                           RK3_TIMESTEP, RK3_CFL, RK3_DT, RK3_SCHEME, NS_DEFORMED, DEALIASING, &
                           LES_CUTOFF, LES_CEPSILON, &
                           P_TRACK_PART, P_INIT_TYPE, P_GAUSS_SIGMAX, P_GAUSS_SIGMAY, P_GAUSS_SIGMAZ, &
                           P_NP_GLOBAL, P_NP_GLOBAL, P_GRAVITY, P_TS_RELEASE, P_TWO_WAY, &
                           P_CHARGE_PART, P_TS_CHARGE_RELEASE, P_TS_FREEZE, P_INTERP_ORDER, P_TS_D2VALS, &
                           P_D_NCELLS, P_FOLLOW, P_KEEP_INITIAL_POSITION, P_RHO, P_DIADIST_TYPE, &
                           P_DIA_LVLS, P_DIA_VALS, P_DIA_NPS, P_SURF_TENSION, &
                           P_REF_CHARGE, P_NP_PACKET_CH, MASS_LOADING, PHI_PERCENT

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  integer               :: iproc
 
  namelist /nl_verbose/         MYRANK, V_ALLOCS, V_DEALLOCS, V_MPI, V_DFT, V_GENERAL, V_DEBUG, V_LOG, V_TEST, V_64BIT 
  namelist /nl_control/         MYRANK, RUNMODE, TSTEPS, TS_REPORT, TS_FLUIDVTK, TS_CORR, TS_SPEC3D, TS_SPEC1D, TS_PLANE, &
                                F_MOMENTS, NL_BUFFER
  namelist /nl_domain/          MYRANK, NODES_X, NODES_Y, NODES_Z, B110, B220, B330
  namelist /nl_gravity/         MYRANK, GRAVITY 
  namelist /nl_random_setup/    MYRANK, RANDOM_SEED_SIZE, RANDOM_SEED_FLUID, RANDOM_SEED_PARTICLES
  namelist /nl_fluid/           MYRANK, F_RHO, F_NU
  namelist /nl_flow/            MYRANK, F_TYPE, F_DEFORM_A, F_DEFORM_B, F_DEFORM_C, F_DEFORM_S, F_REMESH_T 
  namelist /nl_f_init/          MYRANK, F_INIT_TYPE, F_INIT_U0, F_INIT_KP, F_INIT_KV, F_INIT_KSTART, F_INIT_KEND, &
                                F_INIT_A, F_INIT_B, F_INIT_C, F_INIT_I, F_INIT_J, F_INIT_K
  namelist /nl_force/           MYRANK, FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP
  namelist /nl_restart/         F_RESTART_WRITE, F_RESTART_FILE_READ, F_RESTART_FILE_WRITE, F_RESTART_SCALE, &
                                F_FORCE_FILE_READ, F_FORCE_FILE_WRITE, &
                                P_RESTART_WRITE, P_RESTART_FILE_READ, P_RESTART_FILE_WRITE
  namelist /nl_numerics/        MYRANK, RK3_TIMESTEP, RK3_CFL, RK3_DT, RK3_SCHEME, NS_DEFORMED, DEALIASING, &
                                LES_CUTOFF, LES_CEPSILON
  namelist /nl_particle/        MYRANK, P_TRACK_PART, P_INIT_TYPE, P_GAUSS_SIGMAX, P_GAUSS_SIGMAY, P_GAUSS_SIGMAZ, &
                                P_NP_GLOBAL, P_NP_GLOBAL, P_GRAVITY, P_TS_RELEASE, &
                                P_TWO_WAY, P_CHARGE_PART, P_TS_CHARGE_RELEASE, P_TS_FREEZE, P_INTERP_ORDER, &
                                P_TS_D2VALS, P_D_NCELLS, P_FOLLOW, P_KEEP_INITIAL_POSITION
  namelist /nl_part_properties/ MYRANK, P_RHO, P_DIADIST_TYPE, P_DIA_LVLS, P_DIA_VALS, P_DIA_NPS, P_SURF_TENSION, &
                                P_REF_CHARGE, P_NP_PACKET_CH 
  namelist /nl_coupling/        MYRANK, MASS_LOADING, PHI_PERCENT

 
 
!> Using MPI_Barrier to make the output more likely to be organised by increasing process rank.
!! There is no guarantee for that though. 
  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking parameters on all processes. \n', ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking verbose parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOVERBOSE: do iproc=0, NPROCS-1, 1
    IFRANKVERBOSE: if ( MYRANK == iproc ) then
      write(*,nml=nl_verbose) 
    end if IFRANKVERBOSE
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOVERBOSE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking control parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOCONTROL: do iproc=0, NPROCS-1, 1
    IFRANKCONTROL: if ( MYRANK == iproc ) then
      write(*,nml=nl_control) 
    end if IFRANKCONTROL
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOCONTROL

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking domain parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DODOMAIN: do iproc=0, NPROCS-1, 1
    IFRANKDOMAIN: if ( MYRANK == iproc ) then
      write(*,nml=nl_domain) 
    end if IFRANKDOMAIN
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DODOMAIN

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking gravity parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOGRAVITY: do iproc=0, NPROCS-1, 1
    IFRANKGRAVITY: if ( MYRANK == iproc ) then
      write(*,nml=nl_gravity) 
    end if IFRANKGRAVITY
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOGRAVITY

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking random parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DORANDOM: do iproc=0, NPROCS-1, 1
    IFRANKRANDOM: if ( MYRANK == iproc ) then
      write(*,nml=nl_random_setup) 
    end if IFRANKRANDOM
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DORANDOM

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking fluid parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOFLUID: do iproc=0, NPROCS-1, 1
    IFRANKFLUID: if ( MYRANK == iproc ) then
      write(*,nml=nl_fluid) 
    end if IFRANKFLUID
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOFLUID

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking flow parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOFLOW: do iproc=0, NPROCS-1, 1
    IFRANKFLOW: if ( MYRANK == iproc ) then
      write(*,nml=nl_flow) 
    end if IFRANKFLOW
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOFLOW

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking flow initialisation parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOFINIT: do iproc=0, NPROCS-1, 1
    IFRANKFINIT: if ( MYRANK == iproc ) then
      write(*,nml=nl_f_init) 
    end if IFRANKFINIT
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOFINIT

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking forcing parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOFORCE: do iproc=0, NPROCS-1, 1
    IFRANKFORCE: if ( MYRANK == iproc ) then
      write(*,nml=nl_force) 
    end if IFRANKFORCE
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOFORCE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking restart file parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DORESTART: do iproc=0, NPROCS-1, 1
    IFRANKRESTART: if ( MYRANK == iproc ) then
      write(*,nml=nl_restart) 
    end if IFRANKRESTART
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DORESTART

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking numerics parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DONUMERICS: do iproc=0, NPROCS-1, 1
    IFRANKNUMERICS: if ( MYRANK == iproc ) then
      write(*,nml=nl_numerics) 
    end if IFRANKNUMERICS
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DONUMERICS

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking particle parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOPARTICLE: do iproc=0, NPROCS-1, 1
    IFRANKPARTICLE: if ( MYRANK == iproc ) then
      write(*,nml=nl_particle) 
    end if IFRANKPARTICLE
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOPARTICLE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking particle properties... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOPARTPROP: do iproc=0, NPROCS-1, 1
    IFRANKPARTPROP: if ( MYRANK == iproc ) then
      write(*,nml=nl_part_properties) 
    end if IFRANKPARTPROP
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOPARTPROP

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking two-way coupling parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOCOUPLING: do iproc=0, NPROCS-1, 1
    IFRANKCOUPLING: if ( MYRANK == iproc ) then
      write(*,nml=nl_coupling) 
    end if IFRANKCOUPLING
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOCOUPLING

end subroutine tg_CheckValuesParameters
!---------------------------------------------------------------------------------

end module tg_checkvalues
