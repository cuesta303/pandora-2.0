!---------------------------------------------------------------------------------
! module g_control
!> This is the central module that is responsible for initialising the program, 
!! time-stepping and finalising the program.
module g_control

use iso_c_binding
 
!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscdm.h>
!#include <petsc/finclude/petscdmda.h>
use g_petsc
!use petscsys
!use petscdm
!use petscdmda

!------- data section begins ----------------------
implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

!> @var date_start 
!> start date of the calculation
!> @var date_finish 
!> finish date of the calculation
!> @var time_start 
!> start time of the calculation
!> @var time_finish 
!> end time of the calculation
character(len=8)          :: date_start 
character(len=8)          :: date_finish
character(len=10)         :: time_start 
character(len=10)         :: time_finish

public :: ProgramInit
public :: Calculate
public :: ProgramFinalise

!---------- data section ends ----------------------
contains
!---------- procedures -----------------------------

!---------------------------------------------------------------------------------
! subroutine ProgramInit
!> Either prepare for simulation or run testing routines.
!> @param ierr should return 0
subroutine ProgramInit ( ierr )

  use g_headers,       only : g_HeadersInit, g_HeadersSimulation, g_HeadersTest
  use g_parameters,    only : g_ParametersInit, P_TRACK_PART, YES, NO, V_TEST
                                              
  use g_constants,     only : g_ConstantsInit
  use g_files,         only : g_FilesInput, g_FilesOutput
  use g_domain,        only : g_DomainInit

  use f_arrays,        only : f_ArraysInit
  use f_fftw,          only : f_FftwInit
  use f_force,         only : f_ForceInit
  use f_fluidstats,    only : f_FluidstatsInit
  use f_initialise,    only : f_InitialiseInit

  use p_control,       only : p_ControlInit
  use g_rk3,           only : g_RK3Init

  use tg_checkvalues,  only : tg_CheckValuesParameters
  use tg_checkdomain,  only : tg_CheckDomainParameters
  use tg_checkrandom,  only : tg_CheckRandomSeed, tg_CheckRandomNormal
  use tf_checkarrays,  only : tf_CheckArraysInit, tf_CheckArraysDimension , &
                              tf_CheckVecGetArray
  use tf_rk3,          only : tf_RK3Gauss, tf_RK3TaylorGreen, &
                              tf_RK3TaylorGreenReverse, tf_RK3CheckNonlinear

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  PetscErrorCode        :: perr

!!> Initialise PETSc and MPI. 
!  call PetscInitialize ( PETSC_NULL_CHARACTER, ierr )

!> Record memory usage.
  call PetscMemorySetGetMaximumUsage ( ierr )

!> Write program header (first few lines in standard output file).
  call g_HeadersInit ( ierr )                         

  ! INITIALISE GENERAL MODULES                                        
!> Find out the rank of each process and find MASTER 
  call g_ParametersInit ( ierr )                      
!> Initialise constants for reuse
  call g_ConstantsInit ( ierr )                       
!> Open input files
  call g_FilesInput ( ierr )                          

!> Initialise program control
  call g_ControlInit ( ierr )                         
!> Find domain parameters in both real space and wave space
  call g_DomainInit ( ierr )                          
!> Allocate wave-space arrays (kept through the whole simulation)
  call f_ArraysInit ( ierr )                          
!> Allocate FFTW buffers and plans
  call f_FftwInit ( ierr )
!> Initialise fluid statistics
  call f_FluidstatsInit ( ierr )

!> If there are particles, call the initialisation routines.
  IFPARTICLES: if ( P_TRACK_PART == YES ) then
    call p_ControlInit ( ierr )
  end if IFPARTICLES

!> Initialise Runge-Kutta routines
  call g_RK3Init ( ierr )

!> Initialise forcing
!> Get random seed for flow
  call f_ForceInit ( ierr )

  IFVTEST: if ( V_TEST == NO ) then
    call g_HeadersSimulation ( ierr )
    !> Open output files as necessary
    call g_FilesOutput ( ierr )
    call f_InitialiseInit ( ierr )
  else if ( V_TEST == YES ) then
    call PetscPrintf( PETSC_COMM_WORLD, '==> Running test routines. \n', ierr )
    call g_HeadersTest ( ierr )
!    call tg_CheckValuesParameters ( ierr )
!    call tg_CheckDomainParameters ( ierr )
    call tg_CheckRandomSeed ( ierr )
!    call tg_CheckRandomNormal ( ierr )
!    call tf_CheckArraysInit ( ierr )
!    call tf_CheckArraysDimension ( ierr )
!    call tf_CheckVecGetArray ( ierr )
!    call tf_RK3Gauss ( ierr )
!    call tf_RK3TaylorGreen ( ierr )
!    call tf_RK3TaylorGreenReverse ( ierr )
!    call tf_RK3CheckNonlinear ( ierr )
  end if IFVTEST


end subroutine ProgramInit

!---------------------------------------------------------------------------------
! subroutine Calculate
!> This subroutine does the time-stepping for the simulation depending on the configuration.
!! Within each time-step, typically calls to the Runge-Kutta routines are followed by flow
!! and particle statistics and output routines.
!> @param ierr should return 0
subroutine Calculate ( ierr )

  use g_headers,       only : g_HeadersCalc
  use g_parameters,    only : YES, NO, V_TEST, TSTEPS

  implicit none
  
  PetscErrorCode,intent(inout)       :: ierr
  integer,save                       :: tstep
  real(kind=C_DOUBLE)                :: simtime

  !> Start the clock
  call date_and_time ( date=date_start, time=time_start )

  !> Set simulation time to zero
  simtime = 0.0

  !> Do not enter the time loop if this is only a test run.
  IFTEST: if ( V_TEST == NO ) then
    !> Write header for the simulation
    call g_HeadersCalc ( ierr )

    !> Enter the time loop for the simulation
    TIMELOOP: DO tstep=1, TSTEPS

      !> Decide what to report during this time step and write header if applicable
      call g_ControlCalcReport ( tstep, simtime, ierr )

!      !> Remesh for shear flow if at 1/2
!      call g_ControlCalcRemesh ( tstep, ierr )
 
      !> Runge-Kutta routines for this time step
      call g_ControlCalcRK3 ( tstep, simtime, ierr )

!      !> Compute statistics and write to output files as required.
!      !! Real-space statistics best in first RK3 stage rather than here

!> @ todo Fluid statistics are now in f_rk3 - decide if to move back here
!      call g_ControlCalcOutput ( tstep, simtime, ierr )

    end do TIMELOOP
  end if IFTEST

  !> Stop the clock
  call date_and_time ( date=date_finish, time=time_finish )

end subroutine Calculate

!---------------------------------------------------------------------------------
! subroutine ProgramFinalise
!> This subroutine is responsible for deallocating all arrays, closing all files
!! and finalising MPI.
!> @param ierr should return 0
subroutine ProgramFinalise ( ierr )

  use g_parameters,  only : P_TRACK_PART, F_RESTART_WRITE, P_RESTART_WRITE, &
                            F_RESTART_SCALE, YES, WRITE_FILE

  use g_headers,     only : g_HeadersFinal
  use g_restart,     only : g_RestartPetscFluid
  use g_files,       only : g_FilesFinalise

  use f_arrays,      only : f_ArraysFinalise
  use f_force,       only : f_ForceFinalise
  use f_fluidstats,  only : f_FluidstatsFinalise
  use f_fftw,        only : f_FftwFinalise

  use p_control,     only : p_ControlFinalise

  implicit none

  PetscErrorCode,intent(inout)    :: ierr

  PetscErrorCode        :: perr

!> write headers for finalise phase
  call g_HeadersFinal ( ierr )
!> write timing information
  call g_ControlFinalTiming ( ierr )
!> write memory information
  call PetscMemoryView ( PETSC_VIEWER_STDOUT_WORLD, ' Memory use during this run: ', ierr ) 
!> write restart files
  IFFLUIDRESTART: if ( F_RESTART_WRITE  == YES ) then
    call g_RestartPetscFluid ( WRITE_FILE, ierr )
  end if IFFLUIDRESTART
!> finalise forcing
  call f_ForceFinalise ( ierr )
!> finalise fluid statistics
  call f_FluidstatsFinalise ( ierr )

  IFPARTICLES: if ( P_TRACK_PART == YES ) then
  !> if particles in simulation, finalise particle routines
    call p_ControlFinalise ( ierr )
  end if IFPARTICLES

!> close all files
  call g_FilesFinalise ( ierr )
!> deallocate FFTW buffers and plans
  call f_FftwFinalise ( ierr )
!> deallocate fluid arrays
  call f_ArraysFinalise ( ierr )
!!> finalise PETSc and MPI
!  call PetscFinalize ( ierr )


end subroutine ProgramFinalise

!---------------------------------------------------------------------------------
! subroutine g_ControlInit
!>  Read control parameters from the configuration and input files and distribute
!! them to all processes.
!> @param ierr should return 0
subroutine g_ControlInit ( ierr )

  use g_constants, only : ZERO
  use g_parameters, only : MYRANK, MASTER, READ_FILE, WRITE_FILE, YES, NO, &
                           V_ALLOCS, V_DEALLOCS, V_MPI, V_DFT, V_GENERAL, V_DEBUG, V_LOG, V_TEST, &  
                           F_FORCING, F_TYPE, F_ISOTROPIC, F_ISOTROPIC_DEFORMED, FORCE_STOP, F_INIT_TYPE, &
                           F_INIT_MINUSFIVETHIRD, F_INIT_EPSILON, F_NU, F_INIT_KSTART, F_INIT_KEND

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,10)

!> read control parameters from input file
    write(*,20)
    call g_ControlRW ( READ_FILE, ierr )

    IFVTEST: if ( V_TEST == NO ) then
      if ( V_ALLOCS == YES )   write(*,30)
      if ( V_DEALLOCS == YES ) write(*,40)
      if ( V_MPI == YES )      write(*,50)
      if ( V_DFT == YES )      write(*,60)
      if ( V_GENERAL == YES )  write(*,70)
      if ( V_DEBUG == YES )    write(*,80)
      if ( V_LOG == YES )      write(*,90)
    else if ( V_TEST == YES ) then 
      write(*,100)
    end if IFVTEST

!> write control parameters to output file
    write(*,110)
    call g_ControlRW ( WRITE_FILE, ierr )


!> broadcast control parameters to all processes
    write(*,120)
  end if IFMASTER

  call g_ControlBcast ( ierr )

!> Derived parameters
  !> Auxiliary variable to see if forcing is used
  F_FORCING = NO

  IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
    IFFORCED: if ( FORCE_STOP > 0 ) then 
      F_FORCING = YES
    end if IFFORCED
  else if ( F_TYPE == F_ISOTROPIC_DEFORMED ) then
    IFFORCED2: if ( FORCE_STOP > 0 ) then 
      F_FORCING = YES
    end if IFFORCED2
  end if IFISOTROPIC

  IFFIVETHIRD: if ( F_INIT_TYPE == F_INIT_MINUSFIVETHIRD ) then
    F_INIT_EPSILON = ( (3.0/2.0) * F_NU * (3.0/2.0) * &
                       ( F_INIT_KSTART**(4.0/3.0) ) &
                     * ( ( F_INIT_KEND /  F_INIT_KSTART )**(4.0/3.0) &
                       - (11.0/15.0) ) )**3
  end if IFFIVETHIRD


10  format('==> Initialising g_control module')
20  format('    | Reading control parameters from files')

30  format('    | Verbose array allocs')
40  format('    | Verbose array deallocs')
50  format('    | Verbose mpi calls')
60  format('    | Verbose dft calls')
70  format('    | Verbose general information')
80  format('    | Verbose debug information')
90  format('    | Verbose profiling information')
100 format('    | Perform built-in tests')

110 format('    | Writing control parameters to file')
120 format('    | Broadcasting control parameters')
end subroutine g_ControlInit

!---------------------------------------------------------------------------------
! subroutine g_ControlRW
!> Either read control parameters from the configuration and input files
!! or write them to the log file.
!> @param status is either READ_FILE or WRITE_FILE, where READ_FILE means reading from the configuration
!! and input files and WRITE_FILE means writing to the log file.
!> @param ierr should return 0
subroutine g_ControlRW ( status, ierr )
  use g_parameters, only : MYRANK, MASTER, READ_FILE, WRITE_FILE, &
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
                           F_FORCE_FILE_READ, F_FORCE_FILE_WRITE, P_RESTART_WRITE, P_RESTART_FILE_READ, & 
                           P_RESTART_FILE_WRITE, RK3_TIMESTEP, RK3_CFL, RK3_DT, RK3_SCHEME, NS_DEFORMED, &
                           DEALIASING, LES_CUTOFF, LES_CEPSILON, &
                           P_TRACK_PART, P_INIT_TYPE, P_GAUSS_SIGMAX, P_GAUSS_SIGMAY, P_GAUSS_SIGMAZ, &
                           P_NP_GLOBAL, P_NP_GLOBAL, P_GRAVITY, P_TS_RELEASE, P_TWO_WAY, &
                           P_CHARGE_PART, P_TS_CHARGE_RELEASE, P_TS_FREEZE, P_INTERP_ORDER, P_TS_D2VALS, &
                           P_D_NCELLS, P_FOLLOW, P_KEEP_INITIAL_POSITION, P_RHO, P_DIADIST_TYPE, &
                           P_DIA_LVLS, P_DIA_VALS, P_DIA_NPS, P_SURF_TENSION, &
                           P_REF_CHARGE, P_NP_PACKET_CH, MASS_LOADING, PHI_PERCENT

  use g_files,      only : FUNIT_INPUT, FUNIT_CONFIG, FUNIT_OUTPUT

  implicit none

  integer,intent(in)              :: status
  PetscErrorCode,intent(inout)           :: ierr

  character(len=1)                :: env64bit

  namelist /nl_verbose/         V_ALLOCS, V_DEALLOCS, V_MPI, V_DFT, V_GENERAL, V_DEBUG, V_LOG, V_TEST, V_64BIT 
  namelist /nl_control/         RUNMODE, TSTEPS, TS_REPORT, TS_FLUIDVTK, TS_CORR, TS_SPEC3D, TS_SPEC1D, TS_PLANE, F_MOMENTS, &
                                NL_BUFFER
  namelist /nl_domain/          NODES_X, NODES_Y, NODES_Z, B110, B220, B330
  namelist /nl_gravity/         GRAVITY 
  namelist /nl_random_setup/    RANDOM_SEED_SIZE, RANDOM_SEED_FLUID, RANDOM_SEED_PARTICLES
  namelist /nl_fluid/           F_RHO, F_NU
  namelist /nl_flow/            F_TYPE, F_DEFORM_A, F_DEFORM_B, F_DEFORM_C, F_DEFORM_S, F_REMESH_T 
  namelist /nl_f_init/          F_INIT_TYPE, F_INIT_U0, F_INIT_KP, F_INIT_KV, F_INIT_KSTART, F_INIT_KEND, &
                                F_INIT_A, F_INIT_B, F_INIT_C, F_INIT_I, F_INIT_J, F_INIT_K
  namelist /nl_force/           FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP
  namelist /nl_restart/         F_RESTART_WRITE, F_RESTART_FILE_READ, F_RESTART_FILE_WRITE, F_RESTART_SCALE, &
                                F_FORCE_FILE_READ, F_FORCE_FILE_WRITE, &
                                P_RESTART_WRITE, P_RESTART_FILE_READ, P_RESTART_FILE_WRITE
  namelist /nl_numerics/        RK3_TIMESTEP, RK3_CFL, RK3_DT, RK3_SCHEME, NS_DEFORMED, DEALIASING, &
                                LES_CUTOFF, LES_CEPSILON
  namelist /nl_particle/        P_TRACK_PART, P_INIT_TYPE, P_GAUSS_SIGMAX, P_GAUSS_SIGMAY, P_GAUSS_SIGMAZ, &
                                P_NP_GLOBAL, P_NP_GLOBAL, P_GRAVITY, P_TS_RELEASE, &
                                P_TWO_WAY, P_CHARGE_PART, P_TS_CHARGE_RELEASE, P_TS_FREEZE, &
                                P_INTERP_ORDER, P_TS_D2VALS, P_D_NCELLS, P_FOLLOW, P_KEEP_INITIAL_POSITION
  namelist /nl_part_properties/ P_RHO, P_DIADIST_TYPE, P_DIA_LVLS, P_DIA_VALS, P_DIA_NPS, P_SURF_TENSION, &
                                P_REF_CHARGE, P_NP_PACKET_CH 
  namelist /nl_coupling/        MASS_LOADING, PHI_PERCENT


  SELECTRW: select case ( status )

    case ( READ_FILE )
    ! read control parameters from configuration file
      rewind ( unit=FUNIT_CONFIG )
      read ( unit=FUNIT_CONFIG, nml=nl_verbose )
  
      !> Get 64 bit variable from corresponding environment variable
      call get_environment_variable ( 'PANDORA_64BIT', env64bit ) 
      read(env64bit,20) V_64BIT
      write(*,*) '64bit: ', V_64BIT

    ! read control parameters from input file
      rewind ( unit=FUNIT_INPUT )
      read ( unit=FUNIT_INPUT, nml=nl_control )
      read ( unit=FUNIT_INPUT, nml=nl_domain )
      read ( unit=FUNIT_INPUT, nml=nl_gravity )
      read ( unit=FUNIT_INPUT, nml=nl_random_setup )
      read ( unit=FUNIT_INPUT, nml=nl_fluid )
      read ( unit=FUNIT_INPUT, nml=nl_flow )
      read ( unit=FUNIT_INPUT, nml=nl_f_init )
      read ( unit=FUNIT_INPUT, nml=nl_force )
      read ( unit=FUNIT_INPUT, nml=nl_restart )
      read ( unit=FUNIT_INPUT, nml=nl_numerics )
      read ( unit=FUNIT_INPUT, nml=nl_particle )
      read ( unit=FUNIT_INPUT, nml=nl_part_properties )
      read ( unit=FUNIT_INPUT, nml=nl_coupling )

    ! write control parameters to output file
    case ( WRITE_FILE )
      write ( unit=FUNIT_OUTPUT, nml=nl_control )
      write ( unit=FUNIT_OUTPUT, nml=nl_domain )
      write ( unit=FUNIT_OUTPUT, nml=nl_gravity )
      write ( unit=FUNIT_OUTPUT, nml=nl_random_setup )
      write ( unit=FUNIT_OUTPUT, nml=nl_fluid )
      write ( unit=FUNIT_OUTPUT, nml=nl_flow )
      write ( unit=FUNIT_OUTPUT, nml=nl_f_init )
      write ( unit=FUNIT_OUTPUT, nml=nl_force )
      write ( unit=FUNIT_OUTPUT, nml=nl_restart )
      write ( unit=FUNIT_OUTPUT, nml=nl_numerics )
      write ( unit=FUNIT_OUTPUT, nml=nl_particle )
      write ( unit=FUNIT_OUTPUT, nml=nl_part_properties )
      write ( unit=FUNIT_OUTPUT, nml=nl_coupling )

  end select SELECTRW

  20 format(i1)

end subroutine g_ControlRW

!---------------------------------------------------------------------------------
! subroutine g_ControlBcast
!> Broadcast the parameters read from the configuration and input files to all processes
!> @param ierr should return 0
subroutine g_ControlBcast ( ierr )
  use g_parameters, only : MASTER, YES, &
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
 
  !> broadcast control parameters to all processes
  call MPI_Bcast ( V_ALLOCS                , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_DEALLOCS              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_MPI                   , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_DFT                   , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_GENERAL               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_DEBUG                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_LOG                   , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_TEST                  , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( V_64BIT                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( RUNMODE                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TSTEPS                  , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TS_REPORT               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TS_FLUIDVTK             , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TS_CORR                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TS_SPEC3D               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TS_SPEC1D               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( TS_PLANE                , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_MOMENTS               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( NL_BUFFER               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( NODES_X                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( NODES_Y                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( NODES_Z                 , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( B110                    , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( B220                    , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( B330                    , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  
  call MPI_Bcast ( GRAVITY                 , 3, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( RANDOM_SEED_SIZE        , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( RANDOM_SEED_FLUID       , 20, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( RANDOM_SEED_PARTICLES   , 20, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( F_RHO                   , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_NU                    , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( F_TYPE                  , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_DEFORM_A              , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_DEFORM_B              , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_DEFORM_C              , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_DEFORM_S              , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_REMESH_T              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( F_INIT_TYPE             , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_U0               , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_KP               , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_KV               , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_KSTART           , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_KEND             , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_A                , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_B                , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_C                , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_I                , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_J                , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_INIT_K                , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  
  call MPI_Bcast ( FORCE_TF                , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( FORCE_SIGMA             , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( FORCE_KF                , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( FORCE_STOP              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( F_RESTART_WRITE         , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_RESTART_FILE_READ     , 100, MPI_Character, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_RESTART_FILE_WRITE    , 100, MPI_Character, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_RESTART_SCALE         , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_FORCE_FILE_READ       , 100, MPI_Character, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( F_FORCE_FILE_WRITE      , 100, MPI_Character, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_RESTART_WRITE         , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_RESTART_FILE_READ     , 100, MPI_Character, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_RESTART_FILE_WRITE    , 100, MPI_Character, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( RK3_TIMESTEP            , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( RK3_CFL                 , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( RK3_DT                  , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( RK3_SCHEME              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( NS_DEFORMED             , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( DEALIASING              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( LES_CUTOFF              , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( LES_CEPSILON            , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( P_TRACK_PART            , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_INIT_TYPE             , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_GAUSS_SIGMAX          , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_GAUSS_SIGMAY          , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_GAUSS_SIGMAZ          , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )

  IF64BIT: if ( V_64BIT == YES ) then
    call MPI_Bcast ( P_NP_GLOBAL             , 1, MPI_Integer8, MASTER, PETSC_COMM_WORLD, ierr )
  else
    call MPI_Bcast ( P_NP_GLOBAL             , 1, MPI_Integer , MASTER, PETSC_COMM_WORLD, ierr )
  end if IF64BIT

  call MPI_Bcast ( P_GRAVITY               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_TS_RELEASE            , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_TWO_WAY               , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_CHARGE_PART           , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_TS_CHARGE_RELEASE     , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_TS_FREEZE             , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_INTERP_ORDER          , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_TS_D2VALS             , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_D_NCELLS              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_FOLLOW                , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_KEEP_INITIAL_POSITION , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( P_RHO                   , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_DIADIST_TYPE          , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_DIA_LVLS              , 1, MPI_Integer, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_DIA_VALS              , 9, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_DIA_NPS               , 9, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_SURF_TENSION          , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_REF_CHARGE            , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( P_NP_PACKET_CH          , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )

  call MPI_Bcast ( MASS_LOADING            , 1, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )
  call MPI_Bcast ( PHI_PERCENT             , 9, MPI_Double, MASTER, PETSC_COMM_WORLD, ierr )

end subroutine g_ControlBcast
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_ControlCalcReport
!> Decide what to report during this time step and write header if applicable
!> @param tstep number of the current time step 
!> @param simtime time in the simulation time frame
!> @param ierr should return 0
subroutine g_ControlCalcReport ( tstep, simtime, ierr )

  use g_headers,       only : g_HeadersTstep
  use g_parameters,    only : YES, NO, TS_REPORT, REPORT, TS_CORR, TS_SPEC3D, TS_SPEC1D, WRITE_CORR, MYRANK, MASTER
  use g_domain,        only : KMAX
  use f_arrays,        only : u1_3w, u2_3w, u3_3w
  use f_fluidstats,    only : ln_kol

  implicit none
  
  integer,intent(in)                 :: tstep
  real(kind=C_DOUBLE),intent(in)     :: simtime
  PetscErrorCode,intent(inout)              :: ierr
  PetscErrorCode                     :: perr

  real(kind=C_DOUBLE)   :: u1squared, u2squared, u3squared



      !> Report this time step if requested
      IFREPORTTS: if ( mod ( tstep, TS_REPORT ) == 0 ) then
        REPORT = YES
        call g_HeadersTstep ( tstep, simtime, ierr )
        !> Determine norm of velocity components
!        call VecNorm (u1_3w, NORM_2, u1squared, ierr)
!        call VecNorm (u2_3w, NORM_2, u2squared, ierr)
!        call VecNorm (u3_3w, NORM_2, u3squared, ierr)

      IFMASTER: if ( MYRANK == MASTER ) then
        write(*,10) ln_kol * KMAX
!        write(*,*) 'Mean square velocity: ', u1squared, u2squared, u3squared
      end if IFMASTER

      else 
        REPORT = NO
      end if IFREPORTTS

      !> Write correlations for this time step if requested
      IFCORRTS: if ( mod ( tstep, TS_CORR ) == 0 ) then
        WRITE_CORR = YES
      else 
        WRITE_CORR = NO
      end if IFCORRTS

10  format('    | -  Resolution = ', f22.19 )

end subroutine g_ControlCalcReport
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_ControlCalcRK3
!> Call the RK3 routines.
!> @param tstep number of the current time step 
!> @param simtime time in the simulation time frame
!> @param ierr should return 0
subroutine g_ControlCalcRK3 ( tstep, simtime, ierr )

  use g_parameters,    only : 
  use g_rk3,           only : g_RK3Calc

  implicit none
  
  integer,intent(in)                 :: tstep
  real(kind=C_DOUBLE),intent(inout)  :: simtime
  PetscErrorCode,intent(inout)              :: ierr

!> Implement timing code here.

!> Call g_RK3Calc subroutine for the Runge-Kutta routines.
  call g_RK3Calc ( tstep, simtime, ierr )

!> Implement timing code here.

end subroutine g_ControlCalcRK3
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_ControlCalcOutput
!> Compute fluid and particle statistics and write to output files as required.
!> @param tstep number of the current time step 
!> @param simtime time in the simulation time frame
!> @param ierr should return 0
subroutine g_ControlCalcOutput ( tstep, simtime, ierr )

  use g_parameters,    only : 
  use g_rk3,           only : g_RK3Calc

  implicit none
  
  integer,intent(in)                 :: tstep
  real(kind=C_DOUBLE),intent(inout)  :: simtime
  PetscErrorCode,intent(inout)              :: ierr

!> Compute fluid statistics.

!> If there are particles compute particle statistics.

!> Write on xy plane if applicable

!> Write to restart file if applicable

end subroutine g_ControlCalcOutput
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_ControlCalcRemesh
!> If flow type is shear flow, see if remesh is required 
!! and if so call remesh subroutine.
!> @param tstep number of the current time step 
!> @param ierr should return 0
!subroutine g_ControlCalcRemesh ( tstep, ierr )

!  use g_parameters,    only : F_REMESH_T, F_TYPE, F_ISOTROPIC, &
!                              RK3_TIMESTEP, RK3_TIMESTEP_DECIDE
!  use g_rk3,           only : g_RK3Remesh
!  use f_rk3,           only : f_RK3RemeshFluidRealSpace

!  implicit none
  
!  integer,intent(in)                 :: tstep
!  integer,intent(inout)              :: ierr

!  integer                            :: a, b

  !> if shear flow, check if remeshing at current time step required
!  IFSHEAR: if ( F_TYPE /= F_ISOTROPIC ) then

    !> Only remesh if timestepping method is chosen by PANDORA
!    IFTSTEPDECIDE: if ( RK3_TIMESTEP == RK3_TIMESTEP_DECIDE ) then
  
!      a = tstep - F_REMESH_T - 1
!      b = 2 * F_REMESH_T
    
      !> call g_RK3Remesh if remesh required at this time step
!      IFREMESH: if ( mod ( a, b ) == 0 ) then
!        call g_RK3Remesh ( ierr )
!        call f_RK3RemeshFluidRealSpace ( ierr )
!      end if IFREMESH

!    end if IFTSTEPDECIDE

!  end if IFSHEAR
  
!end subroutine g_ControlCalcRemesh
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_ControlFinalTiming
!> Write start and finish date of the simulation to standard out. 
!> @param tstep number of the current time step 
!> @param ierr should return 0
subroutine g_ControlFinalTiming ( ierr )

  use g_parameters,    only : MYRANK, MASTER 

  implicit none
  
  PetscErrorCode,intent(inout)              :: ierr

  !> if shear flow, check if remeshing at current time step required
  IFMASTER: if ( MYRANK == MASTER ) then
        write(*,10) date_start, time_start
        write(*,20) date_finish, time_finish
        write(*,*)
        write(*,*)
  end if IFMASTER
  
10  format('Run start date : ',A,' and time: ',A)
20  format('Run finish date: ',A,' and time: ',A)

end subroutine g_ControlFinalTiming
!---------------------------------------------------------------------------------

end module g_control
