!--------------------------------------------------
! module g_parameters
!> Defines all control parameters that are provided by either the
!! configuration file or the input file.
module g_parameters

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

use iso_c_binding

!------- data section begins ----------------------

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

!--------------------------------------------------
!general parameters 
!--------------------------------------------------
integer,parameter,public                       :: NO         = 0  
integer,parameter,public                       :: YES        = 1  
!> report \f$ \Delta t \f$ and CFL at current time step?
integer,save,public                            :: REPORT     = NO 
!> write correlations at current time step?
integer,save,public                            :: WRITE_CORR = NO 

!--------------------------------------------------
!read/write parameters 
!--------------------------------------------------
!> read from file = 0
integer,parameter,public                       :: READ_FILE  = 0  
!> write to file  = 1
integer,parameter,public                       :: WRITE_FILE = 1  
!> open file = 0
integer,parameter,public                       :: OPEN_FILE  = 0  
!> close file  = 1
integer,parameter,public                       :: CLOSE_FILE = 1  

!--------------------------------------------------
!MPI parameters 
!--------------------------------------------------
!> master always = 0
integer,parameter,public                       :: MASTER     = 0  
!> rank of current process (to be determined)
!integer,save,public                            :: MYRANK     = 0  
PetscMPIInt,save,public                        :: MYRANK     = 0  
!> number of processes (to be determined)
!integer,save,public                            :: NPROCS     = 0  
PetscMPIInt,save,public                        :: NPROCS     = 0  

!--------------------------------------------------
!control parameters 
!--------------------------------------------------
!> verbose array allocations?
integer,save,public                            :: V_ALLOCS   = 0 
!> verbose array deallocations?
integer,save,public                            :: V_DEALLOCS = 0 
!> verbose MPI operations?
integer,save,public                            :: V_MPI      = 0 
!> verbose dft operations?
integer,save,public                            :: V_DFT      = 0 
!> verbose general information?
integer,save,public                            :: V_GENERAL  = 0 
!> write debug information to file?
integer,save,public                            :: V_DEBUG    = 0 
!> write logging information to file?
integer,save,public                            :: V_LOG      = 0 
!> run test routines?
integer,save,public                            :: V_TEST     = 0 
!> 64 bit integers?
integer,save,public                            :: V_64BIT    = 0 

integer,parameter,public                       :: SIMULATIONMODE = 0
integer,parameter,public                       :: TIMINGMODE  = 1
integer,save,public                            :: RUNMODE     = SIMULATIONMODE 
integer,save,public                            :: TSTEPS      = 0
integer,save,public                            :: TS_REPORT   = 0
integer,save,public                            :: TS_FLUIDVTK = 1 ! Initialise to 1 so modulo won't crash the simulation 
integer,save,public                            :: TS_CORR     = 1 ! in case these parameters are not included in the input
integer,save,public                            :: TS_SPEC3D   = 1 ! file
integer,save,public                            :: TS_SPEC1D   = 1
integer,save,public                            :: TS_PLANE    = 1
integer,save,public                            :: F_MOMENTS   = 0
integer,parameter,public                       :: BUFFER_DNS1 = 1
integer,parameter,public                       :: BUFFER_DNS3 = 1
integer,parameter,public                       :: BUFFER_RDT1 = 1
integer,parameter,public                       :: BUFFER_RDT3 = 1
integer,save,public                            :: NL_BUFFER   = BUFFER_DNS1
!--------------------------------------------------
!domain parameters 
!--------------------------------------------------
integer,save,public                            :: NODES_X     = 0
integer,save,public                            :: NODES_Y     = 0
integer,save,public                            :: NODES_Z     = 0
integer,save,public                            :: NHALO       = 0
!> \f$ 1 / B110 \f$ = length in x direction
integer,save,public                            :: B110        = 1
!> \f$ 1 / B220 \f$ = length in y direction
integer,save,public                            :: B220        = 1 
!> \f$ 1 / B330 \f$ = length in z direction
integer,save,public                            :: B330        = 1 

!-----------------------------------------------------------------
!physical quantities 
!-----------------------------------------------------------------
real(kind=C_DOUBLE),save,public,dimension(3)   :: GRAVITY = [ 0.0, 0.0, -9.81 ]

!-----------------------------------------------------------------
!random parameters 
!-----------------------------------------------------------------
integer,save,public                            :: RANDOM_SEED_SIZE = 0
integer,save,public,dimension(20)              :: RANDOM_SEED_FLUID = 0
integer,save,public,dimension(20)              :: RANDOM_SEED_PARTICLES = 0

!-----------------------------------------------------------------
!fluid properties 
!-----------------------------------------------------------------
real(kind=C_DOUBLE),save,public                :: F_RHO = 0.0
real(kind=C_DOUBLE),save,public                :: F_NU  = 0.0

!-----------------------------------------------------------------
!flow parameters 
!-----------------------------------------------------------------
integer,parameter,public                       :: F_ISOTROPIC          = 1
integer,parameter,public                       :: F_SHEAR_DNS          = 2
integer,parameter,public                       :: F_SHEAR_RDT_VISCOUS  = 3
integer,parameter,public                       :: F_SHEAR_RDT_INVISCID = 4
integer,parameter,public                       :: F_ISOTROPIC_DEFORMED = 5
integer,save,public                            :: F_TYPE = F_ISOTROPIC
!> initial B11 deformation
real(kind=C_DOUBLE),save,public                :: F_DEFORM_A = 0.0
!> initial B22 deformation
real(kind=C_DOUBLE),save,public                :: F_DEFORM_B = 0.0
!> initial B33 deformation
real(kind=C_DOUBLE),save,public                :: F_DEFORM_C = 0.0
real(kind=C_DOUBLE),save,public                :: F_DEFORM_S = 0.0
integer,save,public                            :: F_REMESH_T = 0

!-----------------------------------------------------------------
!initialisation parameters 
!-----------------------------------------------------------------
integer,save,public                            :: F_INIT_TYPE             = 0   !> Fluid initialisation type
integer,parameter,public                       :: F_INIT_NULL             = 0
integer,parameter,public                       :: F_INIT_RESTART          = 1
integer,parameter,public                       :: F_INIT_RESTART_DEFORMED = 2
integer,parameter,public                       :: F_INIT_GAUSS            = 3
integer,parameter,public                       :: F_INIT_KERR             = 4
integer,parameter,public                       :: F_INIT_PULSE            = 5
integer,parameter,public                       :: F_INIT_SQUARE           = 6
integer,parameter,public                       :: F_INIT_TAYLOR_GREEN     = 7
integer,parameter,public                       :: F_INIT_MINUSFIVETHIRD   = 8 
integer,parameter,public                       :: F_INIT_MATSUMOTO        = 9
integer,parameter,public                       :: F_INIT_AHMED_ELGHOBASHI = 10

real(kind=C_DOUBLE),save,public                :: F_INIT_U0      = 0.0 !> Distribution magnitude
real(kind=C_DOUBLE),save,public                :: F_INIT_KP      = 0.0 !> Distribution peak position
real(kind=C_DOUBLE),save,public                :: F_INIT_KV      = 0.0 !> Distribution variance
real(kind=C_DOUBLE),save,public                :: F_INIT_EPSILON = 0.0 !> Prescribed dissipation for -5/3 spectrum
real(kind=C_DOUBLE),save,public                :: F_INIT_KSTART  = 0
real(kind=C_DOUBLE),save,public                :: F_INIT_KEND    = 0
!> Amplitudes of Taylor-Green vortex
real(kind=C_DOUBLE),save,public                :: F_INIT_A = 1.0
real(kind=C_DOUBLE),save,public                :: F_INIT_B = 1.0
real(kind=C_DOUBLE),save,public                :: F_INIT_C = 1.0
!> Wavenumbers of Taylor-Green vortex
integer,save,public                            :: F_INIT_I = 1
integer,save,public                            :: F_INIT_J = 1
integer,save,public                            :: F_INIT_K = 1
!---------------------------------------------------------------------
!forcing parameters 
!---------------------------------------------------------------------
integer,save,public                            :: F_FORCING = 0 ! 0 = no   1 = yes
real(kind=C_DOUBLE),save,public                :: FORCE_TF = 0.0
real(kind=C_DOUBLE),save,public                :: FORCE_SIGMA = 0.0
real(kind=C_DOUBLE),public                     :: FORCE_KF = 2.829
integer,save ,public                           :: FORCE_STOP = 99999999

!---------------------------------------------------------------------
!restart file parameters 
!---------------------------------------------------------------------
!> Write fluid restart file yes or no
integer,save,public                            :: F_RESTART_WRITE  = NO
!> Filename for reading the fluid restart file.
character(len=100),save,public                 :: F_RESTART_FILE_READ
!> Filename for writing the fluid restart file.
character(len=100),save,public                 :: F_RESTART_FILE_WRITE
!> Scaling factor for restart file
integer,save,public                            :: F_RESTART_SCALE = 1
!> Filename for reading the fluid forcing file.
character(len=100),save,public                 :: F_FORCE_FILE_READ
!> Filename for writing the fluid forcing file.
character(len=100),save,public                 :: F_FORCE_FILE_WRITE
!> Write particle restart file yes or no
integer,save,public                            :: P_RESTART_WRITE = NO
!> Filename for reading the particle restart file.
character(len=100),save,public                 :: P_RESTART_FILE_READ
!> Filename for writing the particle restart file.
character(len=100),save,public                 :: P_RESTART_FILE_WRITE

!---------------------------------------------------------------------
!numerics parameters 
!---------------------------------------------------------------------
!> Options for how to determine the RK3 timestep
integer,parameter,public                       :: RK3_TIMESTEP_DECIDE = 1
integer,parameter,public                       :: RK3_TIMESTEP_FIXED  = 2
integer,parameter,public                       :: RK3_TIMESTEP_CFL    = 3
integer,save,public                            :: RK3_TIMESTEP = RK3_TIMESTEP_DECIDE
real(kind=C_DOUBLE),save,public                :: RK3_CFL = 0.75
real(kind=C_DOUBLE),save,public                :: RK3_DT  = 0.01
integer,save,public                            :: RK3_SCHEME = 7
!> Formulation of the deformed Navier-Stokes equation
integer,parameter,public                       :: NS_BF_ROT = 1 ! Rotational form equivalent to Bardino/Ferziger
integer,parameter,public                       :: NS_R_CONV = 2 ! Convective form as in the original Rogallo
integer,parameter,public                       :: NS_R_ROT  = 3 ! Rotational form equivalent to Rogallo
integer,save,public                            :: NS_DEFORMED = NS_BF_ROT
!> Dealising scheme
integer,save,public                            :: DEALIASING_TWOTHIRD  = 1
integer,save,public                            :: DEALIASING_SPHERICAL = 2
integer,save,public                            :: DEALIASING_LES = 3
integer,save,public                            :: DEALIASING = 1
real(kind=C_DOUBLE),save,public                :: LES_CUTOFF = 1.0
real(kind=C_DOUBLE),save,public                :: LES_CEPSILON = 0.0
 
!---------------------------------------------------------------------
!particle parameters 
!---------------------------------------------------------------------
integer,save,public                            :: P_TRACK_PART = NO
integer,parameter,public                       :: P_INIT_RESTART = 1
integer,parameter,public                       :: P_INIT_RANDOM = 2 
integer,parameter,public                       :: P_INIT_UNIFORM = 3
integer,parameter,public                       :: P_INIT_GAUSS = 4
real(kind=C_DOUBLE),save,public                :: P_GAUSS_SIGMAX = 0.1
real(kind=C_DOUBLE),save,public                :: P_GAUSS_SIGMAY = 0.1
real(kind=C_DOUBLE),save,public                :: P_GAUSS_SIGMAZ = 0.1
integer,save,public                            :: P_INIT_TYPE = 0
PetscInt,save,public                           :: P_NP_GLOBAL = 0
integer,save,public                            :: P_GRAVITY = NO
integer,save,public                            :: P_TS_RELEASE = 0
integer,save,public                            :: P_TWO_WAY = NO
integer,save,public                            :: P_CHARGE_PART = NO
integer,save,public                            :: P_TS_CHARGE_RELEASE = 0
integer,save,public                            :: P_TS_FREEZE = 0
integer,save,public                            :: P_INTERP_ORDER = 4
integer,save,public                            :: P_TS_D2VALS = 0
integer,save,public                            :: P_D_NCELLS = 0
integer,save,public                            :: P_FOLLOW = 0
integer,save,public                            :: P_KEEP_INITIAL_POSITION = 0

!---------------------------------------------------------------------
!particle properties 
!---------------------------------------------------------------------
real(kind=C_DOUBLE),save,public                :: P_RHO = 0.0
integer,save,public                            :: P_DIADIST_TYPE = 0
integer,save,public                            :: P_DIA_LVLS = 0
real(kind=C_DOUBLE),save,public,dimension(9)  :: P_DIA_VALS = 0.0
real(kind=C_DOUBLE),save,public,dimension(9)  :: P_DIA_NPS  = 0.0
PetscInt,save,public,dimension(9)             :: P_DIA_CDF  = 0
real(kind=C_DOUBLE),save,public                :: P_SURF_TENSION = 0.0
real(kind=C_DOUBLE),save,public                :: P_REF_CHARGE = 0.0
real(kind=C_DOUBLE),save,public                :: P_NP_PACKET_CH = 0.0

!---------------------------------------------------------------------
!two-way coupling parameters 
!---------------------------------------------------------------------
real(kind=C_DOUBLE),save,public                :: MASS_LOADING = 0.0
real(kind=C_DOUBLE),save,public,dimension(9)  :: PHI_PERCENT

!---------------------------------------------------------------------


public :: g_ParametersInit
!------- data section ends ------------------------

contains

!---------------------------------------------------------------------------------
! subroutine  g_ParametersInit
!> Determines the rank of each process and the total number of processes and any
!! other parameters that need to be determined.
!> @param ierr should return 0
subroutine  g_ParametersInit ( ierr )

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising g_parameters module \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Setting myrank and nprocs \n', ierr )

  call MPI_Comm_Rank ( PETSC_COMM_WORLD, MYRANK, ierr )
  call MPI_Comm_Size ( PETSC_COMM_WORLD, NPROCS, ierr )

end subroutine g_ParametersInit
!---------------------------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!---------------------------------------------------------------------------------
end module g_parameters
