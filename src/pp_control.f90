!---------------------------------------------------------------------------------
! module p_control
!> Contains auxiliary routines for Lagrangian particles, in particular routines
!! for allocating and deallocating memory and for initialising particles.
module p_control

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>


private

!> Number of particles in the centre of the subdomain
integer, save, public :: np_centre
!> Total number of particles associated with a process
integer, save, public :: np_total

integer,dimension(1:8) :: np_send
integer,dimension(1:8) :: np_recv
integer,dimension(1:8) :: np_send_centre
integer,dimension(1:8) :: np_recv_centre

integer,parameter      :: INITIALISATION = -1

!> @todo The following arrays probably belong into the f_arrays module, being Eulerian arrays
!real(kind=C_DOUBLE),save,allocatable,public :: force_partr(:,:,:,:)
!complex(kind=C_DOUBLE),save,allocatable,public :: force_partk(:,:,:,:)
!real(kind=C_DOUBLE),save,allocatable,public :: part_chr(:,:,:,:)
!complex(kind=C_DOUBLE),save,allocatable,public :: part_chk(:,:,:,:)
!real(kind=C_DOUBLE),save,allocatable,public :: part_efieldr(:,:,:,:)
!complex(kind=C_DOUBLE),save,allocatable,public :: part_efieldk(:,:,:,:)
!real(kind=C_DOUBLE),save,allocatable,public :: part_er_verify(:,:,:,:)
!complex(kind=C_DOUBLE),save,allocatable,public :: part_ek_verify(:,:,:,:)
!real(kind=C_DOUBLE),save,allocatable,public :: part_potr(:,:,:,:)
!complex(kind=C_DOUBLE),save,allocatable,public :: part_potk(:,:,:,:)

!----------------------------------------------------------------------

public :: p_ControlInit
public :: p_ControlFinalise
public :: p_ControlFinishTimeStep
public :: p_ControlRK3
public :: p_ControlPrepareWriteRestart 
public :: p_ControlRemesh
public :: p_ControlStats

contains
!---------------------------------------------------------------------------------
! subroutine p_ControlInit 
!> Initialise particle routines.
!> @param ierr should return 0
subroutine p_ControlInit ( ierr )

  use g_parameters, only : MYRANK, MASTER, V_TEST, NO, READ_FILE, WRITE_FILE, P_TRACK_PART, YES, &
                           P_INIT_TYPE, P_INIT_RESTART
  use g_constants,  only : ZERO
  use p_mpi,        only : p_MPIFindNeighbours, p_MPIInitParticleDomain
  use p_arrays,     only : ARRAYHANDLE_XP, ARRAYHANDLE_STATUS, ARRAYHANDLE_DP, p_ArraysAllocate
  use p_interp,     only : p_InterpAllocate
  use p_rk3,        only : xpgrav, rk3xpgrav, upgrav, rk3upgrav

  implicit none

  PetscErrorCode, intent(inout) :: ierr

!  ierr = 0

  if ( V_TEST == NO ) then
    call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising p_control module \n', ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Checking release time is within simulation time \n', ierr )
    call p_ControlCheckReleaseTime ( ierr )

    IFPART: if ( P_TRACK_PART == YES ) then

      call PetscPrintf( PETSC_COMM_WORLD, '    | Finding neighbour nodes \n', ierr )
      call p_MPIFindNeighbours ( ierr )
!      CHKERRQ( ierr )

      call PetscPrintf( PETSC_COMM_WORLD, '    | Initialising particle domain \n', ierr )
      call p_MPIInitParticleDomain ( ierr )
!      CHKERRQ( ierr )

      IFRESTART: if ( P_INIT_TYPE == P_INIT_RESTART ) then 

        !> Allocate restart arrays

        !> Read particle restart file

        !> Read the file for the original position if it is kept

        !> Find out if all particles belong to this process
 
           !> If not find out if all remaining particles belong to the neighbours

              !> If particles are left broadcast them and let each process find
              !! if it is the owner

      else

        call PetscPrintf( PETSC_COMM_WORLD, '    | Estimating local array dimension \n', ierr )
        call p_ControlInitLocalDimensions ( ierr )
!        CHKERRQ( ierr )

        call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating particle position array \n', ierr )
        call p_ArraysAllocate ( ARRAYHANDLE_XP, np_centre, ierr )
!      CHKERRQ( ierr )

        call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating particle status array \n', ierr )
        call p_ArraysAllocate ( ARRAYHANDLE_STATUS, np_centre, ierr )
!      CHKERRQ( ierr )

        call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating particle diameter array \n', ierr )
        call p_ArraysAllocate ( ARRAYHANDLE_DP, np_centre, ierr )
!      CHKERRQ( ierr )

      end if IFRESTART

      call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating interpolation arrays \n', ierr )
      call p_InterpAllocate ( ierr )
!      CHKERRQ( ierr )

      call PetscPrintf( PETSC_COMM_WORLD, '    | Initialising particles \n', ierr )
      call p_ControlInitialiseParticles ( ierr )
!      CHKERRQ( ierr )
!    call p_quants(ierr)

      !> Initialise absolute particle position
      xpgrav = ZERO
      rk3xpgrav  = ZERO
      upgrav = ZERO
      rk3upgrav  = ZERO

    end if IFPART
  end if 

end subroutine p_ControlInit
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_ControlCheckReleasetime ( ierr )

  use g_parameters, only: TSTEPS, P_TS_RELEASE, P_INIT_TYPE, P_INIT_RESTART, &
                          P_TRACK_PART, NO

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  IFRESTART: if ( P_INIT_TYPE == P_INIT_RESTART ) then

    IFLATERELEASE: if ( P_TS_RELEASE /= 0 ) then
      call PetscPrintf( PETSC_COMM_WORLD, '    | Ignoring particle release time. Particle simulation starts immediately. \n', ierr )
      P_TS_RELEASE = 0 
    end if IFLATERELEASE

  else

    IFNOTWITHIN: if ( P_TS_RELEASE > TSTEPS ) then
      call PetscPrintf( PETSC_COMM_WORLD, '    | Particle release time after end of simulation. &
                                                 Switching to fluid simulation. \n', ierr )
      P_TRACK_PART = NO
    end if IFNOTWITHIN

  end if IFRESTART

end subroutine p_ControlCheckReleasetime
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_ControlInitLocalDimensions
!> Allocate initial arrays for Lagrangian particle.
!> @param ierr should return 0
subroutine p_ControlInitLocalDimensions ( ierr )

  use g_parameters, only : MYRANK, MASTER, NPROCS, YES, NO, V_TEST, &
                           NODES_X, NODES_Y, NODES_Z, &
                           P_TS_RELEASE, P_CHARGE_PART, P_TS_CHARGE_RELEASE, &
                           P_TRACK_PART, P_NP_GLOBAL, P_INIT_TYPE, P_INIT_UNIFORM

  implicit none

  PetscErrorCode, intent(inout) :: ierr
  !> Scale the initial position array by factor 3 to 5 depending on number of particles
  !! per process. This is very conservative, but the
  !! memory is definitely available, since we do not allocate velocity/acceleration yet.
  integer                :: scalefactor

  IFSCALE: if ( ( P_NP_GLOBAL / NPROCS ) <= 1 ) then
    scalefactor = 5
  else if ( ( P_NP_GLOBAL / NPROCS ) <= 2 ) then
    scalefactor = 4
  else
    scalefactor = 3
  end if IFSCALE

  !> Do not proceed if in testing mode
  IFNOTTESTING: if ( V_TEST == NO ) then    

    !> If no particles in simulation, just write this to the header
    IFPARTICLES: if ( P_TRACK_PART == 0 ) then

      call PetscPrintf( PETSC_COMM_WORLD, '    | Particle tracking: off \n', ierr )

    !> Proceed if particles in simulation
    else

      call PetscPrintf( PETSC_COMM_WORLD, '    | Particle tracking: on \n', ierr )

      !> Set global number of particles if uniformly distributed
      IFUNIFORM: if ( P_INIT_TYPE == P_INIT_UNIFORM ) then
        P_NP_GLOBAL = NODES_X * NODES_Y * NODES_Z
      end if IFUNIFORM

      !> Initial estimate of local number of particles
      np_centre = nint ( real ( P_NP_GLOBAL, kind=C_DOUBLE ) / &
                         real ( NPROCS, kind=C_DOUBLE ) )
      IFZERO: if ( np_centre == 0 ) then
        np_centre = 1
      end if IFZERO
      np_centre = np_centre * scalefactor
      
      !> We do not need more local particles than global particles
      IFLOCALGTGLOBAL: if ( np_centre > P_NP_GLOBAL ) then
        np_centre = P_NP_GLOBAL
      end if IFLOCALGTGLOBAL

      IFMASTERWRITEPART: if ( MYRANK == MASTER ) then
        write(*,10) P_NP_GLOBAL
        write(*,20) np_centre
        write(*,30) P_TS_RELEASE
      end if IFMASTERWRITEPART

      IFCHARGE: if ( P_CHARGE_PART == YES) then
               
        IFMASTERWRITECHARGE: if ( MYRANK == MASTER ) then
          write(*,50)
          write(*,60) P_TS_CHARGE_RELEASE
        end if IFMASTERWRITECHARGE

        IFEARLYCHARGE: if ( P_TS_CHARGE_RELEASE >= P_TS_RELEASE ) then
          call PetscPrintf( PETSC_COMM_WORLD, ' Terminating program: particles cannot be charged &
                                                before they are released. \n ', ierr )
          ierr = 1
!          CHKERRQ ( ierr )
        end if IFEARLYCHARGE

      end if IFCHARGE
    end if IFPARTICLES
  end if IFNOTTESTING

10 format('    | Global number of particles: ',i10)
20 format('    | Initial number of particles/process: ',i10)
30 format('    | Particles to be released after ',i10,' time steps')
50 format('    | Particle charging: on')
60 format('    | Particles to be charged after ',i10,' time steps')

end subroutine p_ControlInitLocalDimensions
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_ControlFinalise 
!> Finalise particle routines.
!> @param ierr should return 0
subroutine p_ControlFinalise ( ierr )

  use g_parameters, only : MYRANK, MASTER, V_TEST, NO, READ_FILE, WRITE_FILE, &
                           P_RESTART_WRITE, YES
  use g_restart,    only : g_RestartPetscParticle
  use p_arrays,     only : ARRAYHANDLE_XP, ARRAYHANDLE_UP, ARRAYHANDLE_RK3, &
                           ARRAYHANDLE_DP, ARRAYHANDLE_VOL, ARRAYHANDLE_CHARGE, &
                           ARRAYHANDLE_STATUS, p_ArraysDeallocate
  use p_interp,     only : p_InterpDeallocate

  implicit none

  PetscErrorCode, intent(inout) :: ierr

  if ( V_TEST == NO ) then
    call PetscPrintf( PETSC_COMM_WORLD, '==> Finalising p_control module \n', ierr )

    IFRESTART: if ( P_RESTART_WRITE == YES ) then
      call PetscPrintf( PETSC_COMM_WORLD, '    | Writing restart file. \n', ierr )
      call p_ControlPrepareWriteRestart ( ierr )
      call g_RestartPetscParticle ( P_RESTART_WRITE, ierr )
    end if IFRESTART

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating particle position arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_XP, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating particle velocity arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_UP, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating Runge-Kutta arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_RK3, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating particle diameter arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_DP, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating vol arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_VOL, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating charge arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_CHARGE, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating particle status arrays \n', ierr )
    call p_ArraysDeallocate ( ARRAYHANDLE_STATUS, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating interpolation arrays \n', ierr )
    call p_InterpDeallocate ( ierr )

  end if 

end subroutine p_ControlFinalise
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_ControlPrepareInterp 
!subroutine p_ControlPrepareInterp ( ierr )

!  use p_arrays,     only : ptr_xp, ptr_uf, ptr_status
!  use p_interp,     only : p_InterpFindMolecule1D, p_InterpSetToZero 

!  implicit none

!  integer,intent(inout) :: ierr
  !> Degrees of freedom of array to be interpolated
!  integer               :: dof

!  call PetscPrintf( PETSC_COMM_WORLD, '    | Exchanging ghost particles. \n', ierr )


!  IFCENTRE: if ( np_centre > 0 ) then
!    call p_InterpFindMolecule1D ( np_centre, ptr_xp, &
!                                  ptr_status, ierr )    
!    call p_InterpSetToZero ( np_centre, ptr_uf, ierr )    
!  end if IFCENTRE

!end subroutine p_ControlPrepareInterp
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_ControlVelocity2D 
!subroutine p_ControlVelocity2D ( zslice, ierr )

!  use g_parameters, only : MYRANK, NHALO 
!  use p_arrays,     only : ptr_xp_centre, ptr_uf_centre, ptr_status_centre
!  use p_interp,     only : p_InterpInterpolate2D 

!  implicit none

!  integer,intent(in)    :: zslice
!  integer,intent(inout) :: ierr
  !> Degrees of freedom of array to be interpolated
!  integer               :: dof

!  IFCENTRE: if ( np_centre > 0 ) then
!    call p_InterpInterpolate2D ( np_centre, ptr_xp_centre, ptr_uf_centre, &
!                                 ptr_status_centre, zslice, ierr )    
!  end if IFCENTRE

!end subroutine p_ControlVelocity2D
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_ControlInitPartVelocity
subroutine p_ControlInitPartVelocity ( ierr )

  use g_constants,  only : ZERO
  use g_parameters, only : MYRANK, NHALO, P_DIA_VALS, P_KEEP_INITIAL_POSITION, YES
  use p_arrays,     only : particle_xp, particle_up, particle_dp, particle_rk3, particle_up_rk3, particle_xp0
  use p_interp,     only : p_InterpGetLocalVelocity, &
                           p_InterpReturnLocalVelocity, &
                           p_InterpInterpolate3D 

  implicit none

  PetscErrorCode,intent(inout)               :: ierr

  integer                             :: np
  real(kind=C_DOUBLE),dimension(1:3)  :: xp
  real(kind=C_DOUBLE),dimension(1:3)  :: uf

!  write(*,*) ' Number of particles at initialisation: ', np_centre

  IFPARTEXIST: if ( np_centre > 0 ) then
    DOPARTICLES: do np = 1, np_centre, 1
      xp(:) = particle_xp(:,np)
      call p_InterpInterpolate3D ( xp, uf, ierr )
      particle_up(:,np)  = uf(:)
      particle_rk3(:,np) = ZERO
      particle_up_rk3(:,np) = ZERO

      IFKEEPX0: if ( P_KEEP_INITIAL_POSITION == YES ) then
        particle_xp0(:,np) = particle_xp(:,np)
      end if IFKEEPX0

    end do DOPARTICLES

  end if IFPARTEXIST

end subroutine p_ControlInitPartVelocity
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_ControlInitPartVelocityR
subroutine p_ControlInitPartVelocityR ( ierr )

  use g_constants,  only : ZERO
  use g_parameters, only : MYRANK, NHALO, P_DIA_VALS, P_KEEP_INITIAL_POSITION, YES
  use p_arrays,     only : particle_xp, particle_up, particle_dp, particle_rk3, particle_up_rk3, particle_xp0
  use p_interp,     only : p_InterpGetLocalVelocity, &
                           p_InterpReturnLocalVelocity, &
                           p_InterpInterpolate3D 
  use p_rk3,        only : p_RK3Mov2Lab, p_RK3UAbs

  implicit none

  PetscErrorCode,intent(inout)               :: ierr

  integer                             :: np
  real(kind=C_DOUBLE),dimension(1:3)  :: xpmov
  real(kind=C_DOUBLE),dimension(1:3)  :: xplab
  real(kind=C_DOUBLE),dimension(1:3)  :: uf
  real(kind=C_DOUBLE),dimension(1:3)  :: uabs

!  write(*,*) ' Number of particles at initialisation: ', np_centre

  IFPARTEXIST: if ( np_centre > 0 ) then
    DOPARTICLES: do np = 1, np_centre, 1
      xpmov(:) = particle_xp(:,np)
      call p_InterpInterpolate3D ( xpmov, uf, ierr )

      !> Velocity is already in laboratory system. Only need to obtain
      !! mean velocity at coordinate.
      !> Obtain mean velocity
      call p_RK3Mov2Lab ( xpmov, xplab, ierr )
!      call p_RK3UAbs ( xplab, uabs, ierr )

      particle_up(:,np)  = uf(:) !+ uabs(:)
      particle_rk3(:,np) = ZERO
      particle_up_rk3(:,np) = ZERO

      IFKEEPX0: if ( P_KEEP_INITIAL_POSITION == YES ) then
        particle_xp0(:,np) = particle_xp(:,np)
      end if IFKEEPX0

    end do DOPARTICLES

  end if IFPARTEXIST

end subroutine p_ControlInitPartVelocityR
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_ControlInitPartVelocityBF
subroutine p_ControlInitPartVelocityBF ( ierr )

  use g_constants,  only : ZERO
  use g_parameters, only : MYRANK, NHALO, P_DIA_VALS, P_KEEP_INITIAL_POSITION, YES
  use p_arrays,     only : particle_xp, particle_up, particle_dp, particle_rk3, particle_up_rk3, particle_xp0
  use p_interp,     only : p_InterpGetLocalVelocity, &
                           p_InterpReturnLocalVelocity, &
                           p_InterpInterpolate3D 
  use p_rk3,        only : p_RK3Mov2Lab, p_RK3UAbs

  implicit none

  PetscErrorCode,intent(inout)               :: ierr

  integer                             :: np
  real(kind=C_DOUBLE),dimension(1:3)  :: xpmov
  real(kind=C_DOUBLE),dimension(1:3)  :: xplab
  real(kind=C_DOUBLE),dimension(1:3)  :: ufmov
  real(kind=C_DOUBLE),dimension(1:3)  :: uflab
  real(kind=C_DOUBLE),dimension(1:3)  :: uabs

  IFPARTEXIST: if ( np_centre > 0 ) then
    DOPARTICLES: do np = 1, np_centre, 1
      xpmov(:) = particle_xp(:,np)
      call p_InterpInterpolate3D ( xpmov, ufmov, ierr )

      !> Transform velocity to laboratory system
      call p_RK3Mov2Lab ( ufmov, uflab, ierr )

      !> Obtain mean velocity
      call p_RK3Mov2Lab ( xpmov, xplab, ierr )
!      call p_RK3UAbs ( xplab, uabs, ierr )

      particle_up(:,np)  = uflab(:) !+ uabs(:)
      particle_rk3(:,np) = ZERO
      particle_up_rk3(:,np) = ZERO

      IFKEEPX0: if ( P_KEEP_INITIAL_POSITION == YES ) then
        particle_xp0(:,np) = particle_xp(:,np)
      end if IFKEEPX0

    end do DOPARTICLES

  end if IFPARTEXIST

end subroutine p_ControlInitPartVelocityBF
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_ControlRK3
subroutine p_ControlRK3 ( stage, tstep, deltat, ierr )

  use g_parameters, only : F_TYPE, F_ISOTROPIC, P_KEEP_INITIAL_POSITION, &
                           P_GRAVITY, YES, P_TS_RELEASE, NS_DEFORMED, &
                           NS_BF_ROT, NS_R_CONV, NS_R_ROT
  use p_arrays,     only : ptr_xp, ptr_up, ptr_rk3, ptr_dp
  use p_rk3,        only : p_RK3SolveParticle, p_RK3AbsoluteVelocity, &
                           p_RK3SolveParticleR, p_RK3AbsoluteVelocityR, &
                           p_RK3SolveParticleBF, p_RK3AbsoluteVelocityBF
  use p_interp,     only : p_InterpGetLocalVelocity, &
                           p_InterpReturnLocalVelocity

  implicit none

  integer,intent(in)              :: stage
  integer,intent(in)              :: tstep 
  real(kind=C_DOUBLE),intent(in)  :: deltat
  PetscErrorCode,intent(inout)           :: ierr
  !> Degrees of freedom of array to be interpolated
  integer                         :: dof

  !> Obtain real-space fluid arrays with halos (local PETSc vectors) as required for interpolation.
  call p_InterpGetLocalVelocity ( ierr )

  !> Particle velocity needs to equal the fluid velocity at release
  IFTSRELEASE: if ( tstep == P_TS_RELEASE ) then
    IFSTAGE1: if ( stage == 1 ) then
      IFINITISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
        call p_ControlInitPartVelocity ( ierr )
      else if ( NS_DEFORMED == NS_BF_ROT ) then
!        write(*,*) 'Initialising particle velocity.'
        call p_ControlInitPartVelocityBF ( ierr )
      else
        call p_ControlInitPartVelocityR ( ierr )
      end if IFINITISOTROPIC
    end if IFSTAGE1
  end if IFTSRELEASE
      
  IFCENTRE: if ( np_centre > 0 ) then

    IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
      !> Compute absolute velocity for statistics
      call p_RK3AbsoluteVelocity ( np_centre, stage, deltat, ierr )    
      call p_RK3SolveParticle ( np_centre, stage, deltat, ierr )    

    else if ( NS_DEFORMED == NS_BF_ROT ) then
      !> Compute absolute velocity for statistics
!      write(*,*) 'Determining absolute velocity.'
      call p_RK3AbsoluteVelocityBF ( np_centre, stage, deltat, ierr )    
!      write(*,*) 'Solving particle equation of motion.'
      call p_RK3SolveParticleBF ( np_centre, stage, deltat, ierr )    
    else
      !> Compute absolute velocity for statistics
      call p_RK3AbsoluteVelocityR ( np_centre, stage, deltat, ierr )    
      call p_RK3SolveParticleR ( np_centre, stage, deltat, ierr )    
    end if IFISOTROPIC

  end if IFCENTRE

  !> Return real-space fluid arrays
  call p_InterpReturnLocalVelocity ( ierr )

end subroutine p_ControlRK3
!---------------------------------------------------------------------------------


!----------------------------------------------------------------
! subroutine p_ControlInitialiseParticles
!> Initialise the particle position either according to a uniform
!! distribution or according to a random distribution.
!> @param ierr should return 0
subroutine p_ControlInitialiseParticles ( ierr )

  use g_parameters, only : P_INIT_TYPE, P_INIT_UNIFORM, P_INIT_RANDOM, P_INIT_RESTART, &
                           NODES_X, NODES_Y, NODES_Z, P_NP_GLOBAL, NHALO, &
                           MYRANK, P_DIA_LVLS, P_DIA_CDF, P_DIA_VALS
  !> @todo check if dx and lx need to be used in the moving or laboratory frame
  !! - might be best in the moving frame
  use g_domain,     only : DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME, &
                           LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME

  use g_constants,  only : HALF

  use g_random,     only : g_RandomInitParticle

  use p_mpi,        only : PART_OUTSIDE, PART_INSIDE, p_MPIFindLocalParticles, &
                           p_MPIFindLocalParticlesLong
  use p_arrays,     only : particle_xp, particle_status, particle_dp

  implicit none

  !> @todo go through all subroutines and check if ierr is actually returning anything
  PetscErrorCode,intent(inout)                    :: ierr

  integer                                  :: i, j, k, m
  PetscInt                                 :: np
  real(kind=C_DOUBLE),dimension(1:3)       :: rn
  real(kind=C_LONG_DOUBLE),dimension(1:3)  :: xp
  integer                                  :: particlelocation

  np_centre = 0

  !> @todo initialise diameter distribution

  call p_ControlInitialiseDiameter ( ierr )

  select case ( P_INIT_TYPE )
    !-------------------------------
    case ( P_INIT_UNIFORM )

      call PetscPrintf( PETSC_COMM_WORLD, '    | Uniform particle distribution \n', ierr )

      np = 0

      DOUNIFORMI: do i = 1, NODES_X, 1
        DOUNIFORMJ: do j = 1, NODES_Y, 1
          DOUNIFORMK: do k = 1, NODES_Z, 1

            np = np + 1

            xp(1) = ( real(i-1, kind=C_DOUBLE) * DX_MOVINGFRAME ) + ( DX_MOVINGFRAME * HALF )
            xp(2) = ( real(j-1, kind=C_DOUBLE) * DY_MOVINGFRAME ) + ( DY_MOVINGFRAME * HALF )
            xp(3) = ( real(k-1, kind=C_DOUBLE) * DZ_MOVINGFRAME ) + ( DZ_MOVINGFRAME * HALF )

            particlelocation = p_MPIFindLocalParticlesLong ( xp )

            IFINSIDEUNIFORM: if ( particlelocation == PART_INSIDE ) then
              np_centre = np_centre + 1
              particle_status(1,np_centre) = np
              particle_xp(:,np_centre) = xp(:)

              DOUNIDIAMETER: do m = 1, P_DIA_LVLS, 1
                IFUNIDIAMETER: if ( np <= P_DIA_CDF(m) ) then
                  particle_dp(1,np_centre) = P_DIA_VALS(m)
                end if IFUNIDIAMETER
              end do DOUNIDIAMETER

            end if IFINSIDEUNIFORM

          end do DOUNIFORMK
        end do DOUNIFORMJ
      end do DOUNIFORMI

    !-------------------------------
    ! generate random position in domain
    case ( P_INIT_RANDOM )

      !> Need to initialise random seed to ensure the same particles on each rank.
      call g_RandomInitParticle ( ierr )

      call PetscPrintf( PETSC_COMM_WORLD, '    | Random particle distribution \n', ierr )

      DORANDOM: do np = 1, P_NP_GLOBAL, 1

        call random_number ( rn )
        xp(1) = real(rn(1),kind=C_LONG_DOUBLE) * &
                real(LX_MOVINGFRAME,kind=C_LONG_DOUBLE) - &
                ( real(DX_MOVINGFRAME,kind=C_LONG_DOUBLE) * &
                real(HALF,kind=C_LONG_DOUBLE) )
        xp(2) = real(rn(2),kind=C_LONG_DOUBLE) * &
                real(LY_MOVINGFRAME,kind=C_LONG_DOUBLE) - &
                ( real(DY_MOVINGFRAME,kind=C_LONG_DOUBLE) * &
                real(HALF,kind=C_LONG_DOUBLE) )
        xp(3) = real(rn(3),kind=C_LONG_DOUBLE) * &
                real(LZ_MOVINGFRAME,kind=C_LONG_DOUBLE) - &
                ( real(DZ_MOVINGFRAME,kind=C_LONG_DOUBLE) * &
                real(HALF,kind=C_LONG_DOUBLE) )

        particlelocation = p_MPIFindLocalParticlesLong ( xp )

        IFINSIDERANDOM: if ( particlelocation == PART_INSIDE ) then
          np_centre = np_centre + 1
          particle_status(1,np_centre) = np
          particle_xp(:,np_centre) = xp(:)
          
!          write(*,*) ' Inside, MYRANK, xp: ', MYRANK, xp 

          DORANDIAMETER: do i = 1, P_DIA_LVLS, 1
            IFRANDIAMETER: if ( np <= P_DIA_CDF(i) ) then
              particle_dp(1,np_centre) = P_DIA_VALS(i)
            end if IFRANDIAMETER
          end do DORANDIAMETER

        else

!          write(*,*) ' Outside, MYRANK, xp: ', MYRANK, xp 

        end if IFINSIDERANDOM

      end do DORANDOM

  endselect

!  write(*,*) np_centre, ' particles on process ', MYRANK, '.'!, &
!             MYRANK, 'before redistribute: ', particle_status(1,:)

  call p_ControlCheckSum ( np_centre, ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Moving particles to appropriately sized arrays. \n', ierr )

  call p_ControlRedistribute ( INITIALISATION, ierr )

!  IFSTATUSEXIST: if ( allocated ( particle_status ) ) then
!    write(*,*) MYRANK, 'centre after first redistribute: ', particle_status(1,:)
!  end if IFSTATUSEXIST

  !> Check if particles are lost
  call p_ControlCheckSum ( np_centre, ierr )

end subroutine p_ControlInitialiseParticles
!----------------------------------------------------------------

!----------------------------------------------------------------
! Compute particle diameter cdf from pdf.
subroutine p_ControlInitialiseDiameter ( ierr )
 
  use g_constants,  only : ZERO
  use g_parameters, only : P_NP_GLOBAL, P_DIA_LVLS, P_DIA_NPS, P_DIA_CDF, &
                           P_DIA_VALS, MYRANK, MASTER

  implicit none

  PetscErrorCode,intent(inout)    :: ierr
  real(kind=C_DOUBLE)      :: pdfnorm
  
  integer                  :: np
  integer                  :: i
  real(kind=C_DOUBLE)      :: relativecdf

  pdfnorm = sum ( P_DIA_NPS )
  P_DIA_NPS = P_DIA_NPS / pdfnorm

  relativecdf = ZERO

  DOCDF: do i = 1, P_DIA_LVLS, 1
    relativecdf = relativecdf + P_DIA_NPS(i)
    P_DIA_CDF(i) = nint ( relativecdf * real(P_NP_GLOBAL,kind=C_DOUBLE) )
    IFMASTER: if ( MYRANK == MASTER ) then
      write(*,10) nint ( P_DIA_NPS(i) * real(P_NP_GLOBAL,kind=C_DOUBLE) ), P_DIA_VALS(i)
    end if IFMASTER
  end do DOCDF

  !> The following should be true, but just to be safe
  P_DIA_CDF(P_DIA_LVLS) = P_NP_GLOBAL

10 format('    | ', I10, ' particles with diameter ', F15.10)

end subroutine p_ControlInitialiseDiameter
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlCheckSum ( nplocal, ierr )

  use g_parameters, only : MYRANK, P_NP_GLOBAL, V_64BIT, YES

  implicit none

  integer,intent(in)    :: nplocal
  PetscErrorCode,intent(inout) :: ierr

  PetscInt              :: npsum

  IF64BIT: if ( V_64BIT == YES ) then
    call MPI_Allreduce ( nplocal, npsum, 1, MPI_Integer8, MPI_SUM,PETSC_COMM_WORLD, ierr )
  else
    call MPI_Allreduce ( nplocal, npsum, 1, MPI_Integer , MPI_SUM,PETSC_COMM_WORLD, ierr )
  end if IF64BIT

  IFLOSTANDFOUND: if ( npsum /= P_NP_GLOBAL ) then
    write(*,*) ' MYRANK: ', MYRANK, &
               ' Local particles: ', nplocal, ' Global particles: ', npsum, &
               ' should be ', P_NP_GLOBAL
    call PetscPrintf( PETSC_COMM_WORLD, 'WRONG NUMBER OF PARTICLES! \n', ierr )
    ierr = 60
!    CHKERRQ ( ierr )
  end if IFLOSTANDFOUND

end subroutine p_ControlCheckSum
!----------------------------------------------------------------


!----------------------------------------------------------------
subroutine p_ControlRedistribute ( timestep, ierr )

  use g_parameters, only: MYRANK, YES, P_TWO_WAY, P_CHARGE_PART, P_KEEP_INITIAL_POSITION, &
                          F_TYPE, F_ISOTROPIC

  use p_mpi,        only: p_MPIFindParticlesOutside, p_MPICountIncomingParticles, p_MPIMigrateParticles, &
                          p_MPIAllocateBufferMigrate, p_MPIDeallocateBufferMigrate, p_MPICopyToSendBuffer, & 
                          p_MPICopyFromRecvBuffer, &
                          PART_MIGRATE_REAL, PART_MIGRATE_STATUS, &
                          PART_IN_CENTRE, PARTICLE_NOT_FOUND

  use p_arrays,     only: ARRAYHANDLE_XP, ARRAYHANDLE_UP, ARRAYHANDLE_RK3, ARRAYHANDLE_DP, &
                          ARRAYHANDLE_XP0, ARRAYHANDLE_UP_RK3, &
                          ARRAYHANDLE_VOL, ARRAYHANDLE_CHARGE, ARRAYHANDLE_STATUS, DOF_STATUSARRAY, &
                          ptr_xp, ptr_xp0, ptr_up, ptr_rk3, ptr_up_rk3, ptr_dp, &
                          ptr_vol, ptr_charge, ptr_status, &
                          p_ArraysDeallocate, p_ArraysAllocate

  implicit none

  integer,intent(in)                  :: timestep
  PetscErrorCode,intent(inout)               :: ierr

  integer                             :: np
  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation
  !> Position of the particle
  real(kind=C_DOUBLE),dimension(1:3)  :: xp
  !> Degree of freedom of an array
  integer                             :: dof

  integer                             :: np_centre_new
  integer                             :: np_centre_old

  logical                             :: allocated_centre

  !> Update particle information for copying
  np_centre_old = np_centre
  np_total = np_centre
  np_centre = 0
  np_send = 0
  np_recv = 0
  np_send_centre = 0
  np_recv_centre = 0

!  call PetscPrintf( PETSC_COMM_WORLD, 'Starting redistribute. \n', ierr )

!if ( timestep < 2 ) then

!  write(*,*) 'before', MYRANK, np_centre_old

!> Count particles in centre, halos and neighbour nodes
  IFLOCALEXIST: if ( np_centre_old > 0 ) then

    !> When using the Rogallo transform, the particle velocities need to be
    !! updated for particles moving across the boundary.

!    IFNOTISOTROPIC: if ( F_TYPE /= F_ISOTROPIC ) then
!      IFUPNOTINITIAL: if( timestep > INITIALISATION ) then
!        call p_ControlUpdateVelocity ( ptr_xp, ptr_up, np_centre_old, ierr )
!      end if IFUPNOTINITIAL
!    end if IFNOTISOTROPIC

    !> If this is not the initialisation update the initial position array
    IFKEEPX0CORR: if ( P_KEEP_INITIAL_POSITION == YES ) then
      IFX0NOTINITIAL: if( timestep > INITIALISATION ) then
        call p_ControlUpdateInitial ( ptr_xp, ptr_xp0, np_centre_old, ierr )
      end if IFX0NOTINITIAL
    end if IFKEEPX0CORR

    call p_ControlUpdateCount ( ptr_xp, ptr_status, np_centre_old, ierr )
  end if IFLOCALEXIST

!  write(*,*) 'after centre', MYRANK, np_centre

!  write(*,*) 'after ', MYRANK, np_centre

!  write(*,*) MYRANK, timestep, ' Sending:     ', np_send

  !> Find np_recv
  call  p_MPICountIncomingParticles ( np_send_centre, np_recv_centre, ierr )

  np_recv = np_recv_centre
  np_centre_new = np_centre + sum ( np_recv_centre )

!  write(*,*) MYRANK, timestep, ' Receiving:   ', np_recv

  !> Update all arrays now.

!  call PetscPrintf( PETSC_COMM_WORLD, 'Real 3. \n', ierr )

 ! Allocate real 3 send and receive buffers
  dof = 3
  np_total = np_centre
!  write(*,*) 'number of particles: ', np_total, np_centre
  call p_MPIAllocateBufferMigrate ( PART_MIGRATE_REAL, dof, np_total, &
                                    np_send, np_recv, ierr )
!  CHKERRQ( ierr )

 ! --------------------------------------------------------------------------------------------
 !   Use subroutine to copy xp to send buffers
  call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                               np_centre, np_send_centre, &
                               ptr_xp, ptr_status, np_centre_old, ierr )
!  CHKERRQ( ierr )


 !   Use subroutine to send to recv buffers
  call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )
!  CHKERRQ( ierr )

 !   Deallocate xp arrays
  call p_ArraysDeallocate ( ARRAYHANDLE_XP, ierr )
!  CHKERRQ( ierr )

 !   Allocate xp arrays with new dimensions
  call p_ArraysAllocate ( ARRAYHANDLE_XP, np_centre_new, ierr )
!  CHKERRQ( ierr )

 !   Use subroutine to copy from centre and recv buffers to new arrays 
  call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                 np_centre, np_recv_centre, &
                                 ptr_xp, np_centre_new, ierr )
!  CHKERRQ( ierr )


 ! --------------------------------------------------------------------------------------------
  IFINITUP1: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy up to send buffers
    call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                 np_centre, np_send_centre, &
                                 ptr_up, ptr_status, np_centre_old, ierr )
!    CHKERRQ( ierr )

   !   Use subroutine to send to recv buffers
    call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )
!    CHKERRQ( ierr )

   !   Deallocate up arrays
    call p_ArraysDeallocate ( ARRAYHANDLE_UP, ierr )
!    CHKERRQ( ierr )
  end if IFINITUP1


 !   Allocate up arrays with new dimensions
  call p_ArraysAllocate ( ARRAYHANDLE_UP, np_centre_new, ierr )
!  CHKERRQ( ierr )


  IFINITUP2: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy from centre and recv buffers to new arrays 
    call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                   np_centre, np_recv_centre, &
                                   ptr_up, np_centre_new, ierr )
!    CHKERRQ( ierr )
    end if IFINITUP2

 ! --------------------------------------------------------------------------------------------
  IFINITRK31: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy rk3 to send buffers
    call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                 np_centre, np_send_centre, &
                                 ptr_rk3, ptr_status, np_centre_old, ierr )
!    CHKERRQ( ierr )

   !   Use subroutine to send to recv buffers
    call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

   !   Deallocate rk3 arrays
    call p_ArraysDeallocate ( ARRAYHANDLE_RK3, ierr )
  end if IFINITRK31


 !   Allocate rk3 arrays with new dimensions
  call p_ArraysAllocate ( ARRAYHANDLE_RK3, np_centre_new, ierr )


  IFINITRK32: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy from centre and recv buffers to new arrays 
    call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                   np_centre, np_recv_centre, &
                                   ptr_rk3, np_centre_new, ierr )

  end if IFINITRK32

 ! --------------------------------------------------------------------------------------------
  IFINITUPRK31: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy up to send buffers
    call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                 np_centre, np_send_centre, &
                                 ptr_up_rk3, ptr_status, np_centre_old, ierr )
!    CHKERRQ( ierr )

   !   Use subroutine to send to recv buffers
    call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )
!    CHKERRQ( ierr )

   !   Deallocate up arrays
    call p_ArraysDeallocate ( ARRAYHANDLE_UP_RK3, ierr )
!    CHKERRQ( ierr )
  end if IFINITUPRK31


 !   Allocate up arrays with new dimensions
  call p_ArraysAllocate ( ARRAYHANDLE_UP_RK3, np_centre_new, ierr )
!  CHKERRQ( ierr )


  IFINITUPRK32: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy from centre and recv buffers to new arrays 
    call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                   np_centre, np_recv_centre, &
                                   ptr_up_rk3, np_centre_new, ierr )
!    CHKERRQ( ierr )
  end if IFINITUPRK32

 ! --------------------------------------------------------------------------------------------
  !> Redistribute original particle position only if chosen for particle
  !! statistics
  IFKEEPX0: if ( P_KEEP_INITIAL_POSITION == YES ) then

    IFINITXP01: if ( timestep > INITIALISATION ) then
      !> Use subroutine to copy xp0 to send buffers
      call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                   np_centre, np_send_centre, &
                                   ptr_xp0, ptr_status, np_centre_old, ierr )
!      CHKERRQ( ierr )

      !> Use subroutine to send to recv buffers
      call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

      !> Deallocate rk3 arrays
      call p_ArraysDeallocate ( ARRAYHANDLE_XP0, ierr )
    end if IFINITXP01


   !   Allocate rk3 arrays with new dimensions
    call p_ArraysAllocate ( ARRAYHANDLE_XP0, np_centre_new, ierr )


    IFINITXP02: if ( timestep > INITIALISATION ) then
      !> Use subroutine to copy from centre and recv buffers to new arrays 
      call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                     np_centre, np_recv_centre, &
                                     ptr_xp0, np_centre_new, ierr )

    end if IFINITXP02

  end if IFKEEPX0 


 ! --------------------------------------------------------------------------------------------
   ! Deallocate send and receive buffers
!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  call p_MPIDeallocateBufferMigrate ( PART_MIGRATE_REAL, ierr )
 ! --------------------------------------------------------------------------------------------


!  call PetscPrintf( PETSC_COMM_WORLD, 'Real 1. \n', ierr )


 ! Allocate real 1 send and receive buffers
  dof = 1
  call p_MPIAllocateBufferMigrate ( PART_MIGRATE_REAL, dof, np_total, &
                                    np_send, np_recv, ierr )

 ! --------------------------------------------------------------------------------------------
 ! As above for dp
!  IFINITDP1: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy dp to send buffers
!    write(*,*) 'Calling copy function for dp: ', dof, &
!                                 np_centre, np_halo_min, np_halo_max, &
!                                 np_send_centre, np_send_halo_min, np_send_halo_max, &
!                                 np_centre_old, np_halo_min_old, np_halo_max_old

    call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                 np_centre, np_send_centre, &
                                 ptr_dp, ptr_status, np_centre_old, ierr )

 !  write(*,*) 'Migrate'

   !   Use subroutine to send to recv buffers
   call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

   !   Deallocate dp arrays
    call p_ArraysDeallocate ( ARRAYHANDLE_DP, ierr )
!    write(*,*) 'Allocate'
!  end if IFINITDP1

   !   Allocate dp arrays with new dimensions
  call p_ArraysAllocate ( ARRAYHANDLE_DP, np_centre_new, ierr )

!    write(*,*) 'Copy'

!  IFINITDP2: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy from centre and recv buffers to new arrays 
    call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                   np_centre, np_recv_centre, &
                                   ptr_dp, np_centre_new, ierr )
!  end if IFINITDP2

 ! If 2-way same for vol
  IFTWOWAY: if ( P_TWO_WAY == YES ) then
 ! --------------------------------------------------------------------------------------------
    IFINITVOL1: if ( timestep > INITIALISATION ) then
     !   Use subroutine to copy vol to send buffers
      call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                   np_centre, np_send_centre, &
                                   ptr_vol, ptr_status, np_centre_old, ierr )

     !   Use subroutine to send to recv buffers
      call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

     !   Deallocate vol arrays
      call p_ArraysDeallocate ( ARRAYHANDLE_VOL, ierr )
    end if IFINITVOL1


   !   Allocate vol arrays with new dimensions
    call p_ArraysAllocate ( ARRAYHANDLE_VOL, np_centre_new, ierr )


    IFINITVOL2: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy from centre and recv buffers to new arrays 
      call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                     np_centre, np_recv_centre, &
                                     ptr_vol, np_centre_new, ierr )
    end if IFINITVOL2
  end if IFTWOWAY

 ! If charge same for charge
  IFCHARGE: if ( P_CHARGE_PART == YES ) then
 ! --------------------------------------------------------------------------------------------
    IFINITCHARGE1: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy charge to send buffers
      call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                   np_centre, np_send_centre, &
                                   ptr_charge, ptr_status, np_centre_old, ierr )

   !   Use subroutine to send to recv buffers
      call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

   !   Deallocate charge arrays
      call p_ArraysDeallocate ( ARRAYHANDLE_CHARGE, ierr )
    end if IFINITCHARGE1


   !   Allocate charge arrays with new dimensions
      call p_ArraysAllocate ( ARRAYHANDLE_CHARGE, np_centre_new, ierr )


    IFINITCHARGE2: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy from centre and recv buffers to new arrays 
      call p_MPICopyFromRecvBuffer ( PART_MIGRATE_REAL, dof, &
                                     np_centre, np_recv_centre, &
                                     ptr_charge, np_centre_new, ierr )
    end if IFINITCHARGE2
  end if IFCHARGE

!end if
 ! Deallocate real 1 buffers
!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  call p_MPIDeallocateBufferMigrate ( PART_MIGRATE_REAL, ierr )



!  call PetscPrintf( PETSC_COMM_WORLD, 'Status \n', ierr )

 ! Similar for status arrays
 ! Allocate integer 3 send and receive buffers
  dof = DOF_STATUSARRAY
  call p_MPIAllocateBufferMigrate ( PART_MIGRATE_STATUS, dof, np_total, &
                                    np_send, np_recv, ierr )

 ! --------------------------------------------------------------------------------------------
 !   Use subroutine to copy status to send buffers
  call p_MPICopyToSendBuffer ( PART_MIGRATE_STATUS, dof, DOF_STATUSARRAY, &
                               np_centre, np_send_centre, &
                               ptr_status, ptr_status, np_centre_old, ierr )

 !   Use subroutine to send to recv buffers
  call p_MPIMigrateParticles ( PART_MIGRATE_STATUS, dof, np_send, np_recv, ierr )

 !   Deallocate status arrays
  call p_ArraysDeallocate ( ARRAYHANDLE_STATUS, ierr )

 !   Allocate status arrays with new dimensions
  call p_ArraysAllocate ( ARRAYHANDLE_STATUS, np_centre_new, ierr )

 !   Use subroutine to copy from centre and recv buffers to new arrays 
  call p_MPICopyFromRecvBuffer ( PART_MIGRATE_STATUS, dof, &
                                 np_centre, np_recv_centre, &
                                 ptr_status, np_centre_new, ierr )

 ! --------------------------------------------------------------------------------------------
   ! Deallocate send and receive buffers
!    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
    call p_MPIDeallocateBufferMigrate ( PART_MIGRATE_STATUS, ierr )

 ! --------------------------------------------------------------------------------------------

!end if

! Before handling halos update np
  np_centre = np_centre_new


  !> Check if any particles got lost
!  call PetscPrintf( PETSC_COMM_WORLD, 'Checking if any particles lost. \n', ierr )
  np_total = np_centre 
!  write(*,*) 'timestep: ', timestep, '. ', np_recv, ' particles on ', MYRANK  
  call p_ControlCheckSum ( np_total, ierr )

end subroutine p_ControlRedistribute
!----------------------------------------------------------------


!----------------------------------------------------------------
subroutine p_ControlPrepareWriteRestart (  ierr )

  use g_parameters, only: MYRANK, YES, P_TWO_WAY, P_CHARGE_PART, P_KEEP_INITIAL_POSITION

  use p_mpi,        only: p_MPIFindParticlesOutside, p_MPICountIncomingParticles, p_MPIMigrateParticles, &
                          p_MPIAllocateBufferMigrate, p_MPIDeallocateBufferMigrate, p_MPICopyToSendBuffer, & 
                          p_MPICopyToRestart, PART_RESTART_XP, PART_RESTART_XP0, PART_RESTART_UP, &
                          PART_RESTART_DP, PART_RESTART_STATUS, &
                          HALO_UPDATE_REAL, HALO_UPDATE_STATUS, PART_MIGRATE_REAL, PART_MIGRATE_STATUS, &
                          PART_IN_CENTRE, PARTICLE_NOT_FOUND

  use p_arrays,     only: ARRAYHANDLE_XP, ARRAYHANDLE_XP0, ARRAYHANDLE_UP, ARRAYHANDLE_RK3, ARRAYHANDLE_DP, & 
                          ARRAYHANDLE_VOL, ARRAYHANDLE_CHARGE, ARRAYHANDLE_STATUS, DOF_STATUSARRAY, &
                          ptr_xp, ptr_xp0, ptr_up, ptr_rk3, ptr_dp, ptr_vol, ptr_charge, ptr_status, &
                          p_ArraysDeallocate, p_ArraysAllocateRestartArrays

  implicit none

  PetscErrorCode,intent(inout)               :: ierr

  integer                             :: np
  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation
  !> Position of the particle
  real(kind=C_DOUBLE),dimension(1:3)  :: xp
  !> Degree of freedom of an array
  integer                             :: dof

  integer                             :: np_centre_new
  integer                             :: np_centre_old
  integer                             :: np_total_new

  logical                             :: allocated_centre

  !> Update particle information for copying
  np_centre_old = np_centre
  np_total = np_centre
  np_centre = 0
  np_send = 0
  np_recv = 0
  np_send_centre = 0
  np_recv_centre = 0

  call PetscPrintf( PETSC_COMM_WORLD, '    | Preparing particle arrays for writing to restart file. \n', ierr )

!  write(*,*) 'before', MYRANK, np_halo_min_old, np_centre_old, np_halo_max_old

!> Count particles in centre, halos and neighbour nodes
  IFLOCALEXIST: if ( np_centre_old > 0 ) then
    call p_ControlUpdateCount ( ptr_xp, ptr_status, np_centre_old, ierr )
  end if IFLOCALEXIST

!  write(*,*) 'after centre', MYRANK, np_centre, np_halo_min, np_halo_max

!  write(*,*) 'after ', MYRANK, np_halo_min, np_centre, np_halo_max

!  write(*,*) 'Neighbours: ', np_send

  !> Find np_recv
  call  p_MPICountIncomingParticles ( np_send_centre, np_recv_centre, ierr )

  np_recv = np_recv_centre
  np_centre_new = np_centre + sum ( np_recv_centre )
  np_total_new = np_centre_new
!  write(*,*) MYRANK, 'np_total, np_total_new, np_centre_new, :   ', &
!             np_total, np_total_new, np_centre_new

  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating restart arrays. \n', ierr )
  call p_ArraysAllocateRestartArrays ( np_total_new, ierr )

  !> Copy position, velocity, diameter and status arrays to restart vectors.

!  call PetscPrintf( PETSC_COMM_WORLD, 'Real 3. \n', ierr )

 ! Allocate real 3 send and receive buffers
  dof = 3
  np_total = np_centre
!  write(*,*) 'number of particles: ', np_total, np_centre, np_halo_min, np_halo_max
  call p_MPIAllocateBufferMigrate ( PART_MIGRATE_REAL, dof, np_total, &
                                    np_send, np_recv, ierr )

 ! --------------------------------------------------------------------------------------------
 !   Use subroutine to copy xp to send buffers
  call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                               np_centre, np_send_centre, &
                               ptr_xp, ptr_status, np_centre_old, ierr )

 !   Use subroutine to send to recv buffers
  call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

 !   Deallocate xp arrays
  call p_ArraysDeallocate ( ARRAYHANDLE_XP, ierr )

 !   Use subroutine to copy from centre and recv buffers to new arrays 
  call PetscPrintf( PETSC_COMM_WORLD, '    | Copying xp to restart. \n', ierr )
  call p_MPICopyToRestart ( PART_RESTART_XP, np_centre, np_recv_centre, &
                            np_total_new, ierr )

 ! --------------------------------------------------------------------------------------------
  IFKEEPINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then
    ! --------------------------------------------------------------------------------------------
    !   Use subroutine to copy xp0 to send buffers
    call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                                 np_centre, np_send_centre, &
                                 ptr_xp0, ptr_status, np_centre_old, ierr )

    !   Use subroutine to send to recv buffers
    call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

    !   Deallocate xp0 arrays
    call p_ArraysDeallocate ( ARRAYHANDLE_XP0, ierr )

    !   Use subroutine to copy from centre and recv buffers to new arrays 
    call PetscPrintf( PETSC_COMM_WORLD, '    | Copying xp0 to restart. \n', ierr )
    call p_MPICopyToRestart ( PART_RESTART_XP0, np_centre, np_recv_centre, &
                              np_total_new, ierr )

  end if IFKEEPINITIAL

 ! --------------------------------------------------------------------------------------------
   !   Use subroutine to copy up to send buffers
  call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                               np_centre, np_send_centre, &
                               ptr_up, ptr_status, np_centre_old, ierr )

 !   Use subroutine to send to recv buffers
  call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

 !   Deallocate up arrays
  call p_ArraysDeallocate ( ARRAYHANDLE_UP, ierr )

 !   Use subroutine to copy from centre and recv buffers to new arrays 
  call PetscPrintf( PETSC_COMM_WORLD, 'Copying up to restart. \n', ierr )
  call p_MPICopyToRestart ( PART_RESTART_UP, np_centre, np_recv_centre, &
                            np_total_new, ierr )

 ! --------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------
   ! Deallocate send and receive buffers
!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  call p_MPIDeallocateBufferMigrate ( PART_MIGRATE_REAL, ierr )
 ! --------------------------------------------------------------------------------------------


!  call PetscPrintf( PETSC_COMM_WORLD, 'Real 1. \n', ierr )

 ! Allocate real 1 send and receive buffers
  dof = 1
  call p_MPIAllocateBufferMigrate ( PART_MIGRATE_REAL, dof, np_total, &
                                    np_send, np_recv, ierr )

 ! --------------------------------------------------------------------------------------------
 ! As above for dp
!  IFINITDP1: if ( timestep > INITIALISATION ) then
   !   Use subroutine to copy dp to send buffers
!    write(*,*) 'Calling copy function for dp: ', dof, &
!                                 np_centre, np_halo_min, np_halo_max, &
!                                 np_send_centre, np_send_halo_min, np_send_halo_max, &
!                                 np_centre_old, np_halo_min_old, np_halo_max_old

  call p_MPICopyToSendBuffer ( PART_MIGRATE_REAL, dof, DOF_STATUSARRAY, &
                               np_centre, np_send_centre, &
                               ptr_dp, ptr_status, np_centre_old, ierr )

 !   Use subroutine to send to recv buffers
  call p_MPIMigrateParticles ( PART_MIGRATE_REAL, dof, np_send, np_recv, ierr )

 !   Deallocate dp arrays
  call p_ArraysDeallocate ( ARRAYHANDLE_DP, ierr )

   !   Use subroutine to copy from centre and recv buffers to new arrays 
  call PetscPrintf( PETSC_COMM_WORLD, 'Copying dp to restart. \n', ierr )
  call p_MPICopyToRestart ( PART_RESTART_DP, np_centre, np_recv_centre, &
                            np_total_new, ierr )

 ! Deallocate real 1 buffers
!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  call p_MPIDeallocateBufferMigrate ( PART_MIGRATE_REAL, ierr )

 ! Similar for status arrays
 ! Allocate integer 3 send and receive buffers
  dof = DOF_STATUSARRAY
  call p_MPIAllocateBufferMigrate ( PART_MIGRATE_STATUS, dof, np_total, &
                                    np_send, np_recv, ierr )

 ! --------------------------------------------------------------------------------------------
 !   Use subroutine to copy status to send buffers
  call p_MPICopyToSendBuffer ( PART_MIGRATE_STATUS, dof, DOF_STATUSARRAY, &
                               np_centre, np_send_centre, &
                               ptr_status, ptr_status, np_centre_old, ierr )

 !   Use subroutine to send to recv buffers
  call p_MPIMigrateParticles ( PART_MIGRATE_STATUS, dof, np_send, np_recv, ierr )

 !   Deallocate status arrays
  call p_ArraysDeallocate ( ARRAYHANDLE_STATUS, ierr )

 !   Use subroutine to copy from centre and recv buffers to new arrays 
  call PetscPrintf( PETSC_COMM_WORLD, 'Copying status to restart. \n', ierr )
  call p_MPICopyToRestart ( PART_RESTART_STATUS, np_centre, np_recv_centre, &
                            np_total_new, ierr )


 ! --------------------------------------------------------------------------------------------
 ! Deallocate send and receive buffers
!    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  call p_MPIDeallocateBufferMigrate ( PART_MIGRATE_STATUS, ierr )

 ! --------------------------------------------------------------------------------------------

! Before handling halos update np
  np_centre = np_centre_new

end subroutine p_ControlPrepareWriteRestart
!----------------------------------------------------------------


!----------------------------------------------------------------
subroutine p_ControlFinishTimeStep ( tstep, ierr )

  implicit none

  integer,intent(in)    :: tstep
  PetscErrorCode,intent(inout) :: ierr

  !> Put particles where they belong 
  call p_ControlRedistribute ( tstep, ierr ) 

end subroutine p_ControlFinishTimeStep
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlUpdateCount ( positionpointer, statuspointer, arraynp, ierr )

  use g_parameters, only : MYRANK, P_KEEP_INITIAL_POSITION, YES
  use p_mpi,        only : p_MPIFindLocalParticles, & 
                           p_MPIFindParticlesOutside, p_MPICorrectPosition, & 
                           HALO_UPDATE_REAL, HALO_UPDATE_STATUS, PART_MIGRATE_REAL, PART_MIGRATE_STATUS, &
                           PART_INSIDE, PART_OUTSIDE, PARTICLE_NOT_FOUND

  implicit none

  type(C_PTR),intent(in)                     :: positionpointer
  type(C_PTR),intent(in)                     :: statuspointer
  integer,intent(in)                         :: arraynp
  PetscErrorCode,intent(inout)                      :: ierr
 
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
  integer,pointer,dimension(:,:)             :: particlestatus
  real(kind=C_DOUBLE),dimension(1:3)         :: xp
  integer                                    :: np

  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation

  integer                             :: status = 0

  !> Associate position pointer array with target array.
  call c_f_pointer ( positionpointer, particleposition, [3, arraynp])
  !> Associate status pointer array with target array.
  call c_f_pointer ( statuspointer, particlestatus, [3, arraynp])

!  write(*,*) 'Updating count: ', particleposition

  DOLOCAL: do np = 1, arraynp, 1

    xp = particleposition(:,np)

    !> Correct the particle position if needed.
!    write(*,*) 'Position before correction: ', xp
    call p_MPICorrectPosition ( xp )
    particleposition(:,np) = xp
!    write(*,*) 'Position after correction:  ', particleposition(:,np) 

    !> Find if the particle is inside
    particlelocation = p_MPIFindLocalParticles ( xp )

    !> Update the number of particles in each part.
    !! Use the status array to mark the location of each particle.
    !! This will be used for redistributing the particles.
    IFCOUNT: if ( particlelocation == PART_INSIDE ) then
      np_centre = np_centre + 1
      particlestatus(2,np) = 0 
      particlestatus(3,np) = PART_INSIDE
    else

      particleowner = p_MPIFindParticlesOutside ( xp )
!      write(*,*) 'Particle outside. Rank: ', MYRANK, ' owner: ', particleowner

      IFLOST: if ( particleowner == PARTICLE_NOT_FOUND ) then
        write(*,*) 'Particle ', np, ' not found: ', xp
        call PetscPrintf( PETSC_COMM_WORLD, 'Particle lost. Abort. \n', ierr )
!        ierr = np
!        CHKERRQ(ierr)
        call MPI_Abort ( PETSC_COMM_WORLD, status, ierr )
      else
        np_send(particleowner) = np_send(particleowner) + 1
        particlestatus(2,np) = particleowner
        np_send_centre(particleowner) = np_send_centre(particleowner) + 1
        particlestatus(3,np) = PART_OUTSIDE
!        write(*,*) ' New send count: ', np_send(particleowner), &
!                   ' Rank: ', MYRANK, ' owner: ', particleowner
      end if IFLOST

    end if IFCOUNT 
  end do DOLOCAL

!  write(*,*) 'Finished updating particle count on rank ', MYRANK

!  call PetscPrintf( PETSC_COMM_WORLD, 'Finished updating &
!                    particle count on all processes. \n', ierr )

end subroutine p_ControlUpdateCount
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlUpdateInitial ( positionpointer, initialpositionpointer, arraynp, ierr )

  use g_parameters, only : MYRANK, P_KEEP_INITIAL_POSITION, YES
  use p_mpi,        only : p_MPIFindLocalParticles, & 
                           p_MPIFindParticlesOutside, p_MPICorrectPosition, & 
                           HALO_UPDATE_REAL, HALO_UPDATE_STATUS, PART_MIGRATE_REAL, PART_MIGRATE_STATUS, &
                           PART_INSIDE, PART_OUTSIDE, PARTICLE_NOT_FOUND

  implicit none

  type(C_PTR),intent(in)                     :: positionpointer
  type(C_PTR),intent(in)                     :: initialpositionpointer
  integer,intent(in)                         :: arraynp
  PetscErrorCode,intent(inout)                      :: ierr
 
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: initialposition
  real(kind=C_DOUBLE),dimension(1:3)         :: xp
  integer                                    :: np

  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation

  integer                             :: status = 0

  !> Associate position pointer array with target array.
  call c_f_pointer ( positionpointer, particleposition, [3, arraynp])
  !> Associate initial position pointer array with target array.
  call c_f_pointer ( initialpositionpointer, initialposition, [3, arraynp])

!  write(*,*) 'Updating count: ', particleposition

  DOLOCAL: do np = 1, arraynp, 1

    xp = particleposition(:,np)

    !> Correct the initial particle position if the particle position will be
    !! corrected.
    call p_MPICorrectPosition ( xp )
    
!    write(*,*) np, initialposition(:,np), xp
    initialposition(:,np) = initialposition(:,np) + xp - particleposition(:,np)

  end do DOLOCAL


end subroutine p_ControlUpdateInitial
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlUpdateVelocity ( positionpointer, particlevelocitypointer, arraynp, ierr )

  use g_parameters, only : MYRANK, P_KEEP_INITIAL_POSITION, YES
  use p_mpi,        only : p_MPIFindLocalParticles, & 
                           p_MPIFindParticlesOutside, p_MPICorrectPosition, & 
                           HALO_UPDATE_REAL, HALO_UPDATE_STATUS, PART_MIGRATE_REAL, PART_MIGRATE_STATUS, &
                           PART_INSIDE, PART_OUTSIDE, PARTICLE_NOT_FOUND
  use p_rk3,        only : p_RK3Mov2Lab, p_RK3UAbs

  implicit none

  type(C_PTR),intent(in)                     :: positionpointer
  type(C_PTR),intent(in)                     :: particlevelocitypointer
  integer,intent(in)                         :: arraynp
  PetscErrorCode,intent(inout)                      :: ierr
 
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particlevelocity
  real(kind=C_DOUBLE),dimension(1:3)         :: xpmov
  real(kind=C_DOUBLE),dimension(1:3)         :: xplab
  real(kind=C_DOUBLE),dimension(1:3)         :: upabsold
  real(kind=C_DOUBLE),dimension(1:3)         :: upabsnew
  integer                                    :: np

  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation

  integer                             :: status = 0

  !> Associate position pointer array with target array.
  call c_f_pointer ( positionpointer, particleposition, [3, arraynp])
  !> Associate initial position pointer array with target array.
  call c_f_pointer ( particlevelocitypointer, particlevelocity, [3, arraynp])


  DOLOCAL: do np = 1, arraynp, 1

    xpmov = particleposition(:,np)

    !> Obtain position in laboratory system
    call p_RK3Mov2Lab ( xpmov, xplab, ierr )

    !> Obtain absolute velocity
    call p_RK3UAbs ( xplab, upabsold, ierr )

    !> Correct the initial particle position if the particle position will be
    !! corrected.
    call p_MPICorrectPosition ( xpmov )
    
    !> Obtain position in laboratory system
    call p_RK3Mov2Lab ( xpmov, xplab, ierr )

    !> Obtain absolute velocity
    call p_RK3UAbs ( xplab, upabsnew, ierr )

    particlevelocity(:,np) = particlevelocity(:,np) + upabsnew - upabsold

  end do DOLOCAL


end subroutine p_ControlUpdateVelocity
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlRemesh ( ierr )

  use g_parameters, only : MYRANK, P_KEEP_INITIAL_POSITION, YES
  use p_mpi,        only : p_MPIFindLocalParticles, & 
                           p_MPIFindParticlesOutside, p_MPICorrectPosition, & 
                           HALO_UPDATE_REAL, HALO_UPDATE_STATUS, PART_MIGRATE_REAL, PART_MIGRATE_STATUS, &
                           PART_INSIDE, PART_OUTSIDE, PARTICLE_NOT_FOUND
  use p_arrays,     only : ptr_xp, ptr_xp0, ptr_up

  implicit none

  PetscErrorCode,intent(inout)                      :: ierr
 
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
  real(kind=C_DOUBLE),dimension(1:3)         :: xp
  integer                                    :: np

  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation

  integer                             :: status = 0

  IFINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then

    call p_ControlRemeshInitial ( ptr_xp, ptr_xp0, np_total, ierr )

  end if IFINITIAL

  call p_ControlRemeshParticles ( ptr_xp, ptr_up, np_total, ierr )

!  call PetscPrintf( PETSC_COMM_WORLD, 'Finished remeshing particles. \n', ierr )

end subroutine p_ControlRemesh
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlRemeshParticles ( positionpointer, velocitypointer, arraynp, ierr )

  use g_parameters, only : MYRANK
  use p_mpi,        only : p_MPIRemeshPosition

  implicit none

  type(C_PTR),intent(in)                     :: positionpointer
  type(C_PTR),intent(in)                     :: velocitypointer 
  integer,intent(in)                         :: arraynp
  PetscErrorCode,intent(inout)                      :: ierr
 
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particlevelocity
  real(kind=C_DOUBLE),dimension(1:3)         :: xp
  real(kind=C_DOUBLE),dimension(1:3)         :: up
  integer                                    :: np

  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation

  !> Associate position pointer array with target array.
  call c_f_pointer ( positionpointer, particleposition, [3, arraynp])
  call c_f_pointer ( velocitypointer, particlevelocity, [3, arraynp])

  DOLOCAL: do np = 1, arraynp, 1

    xp = particleposition(:,np)
    up = particlevelocity(:,np)

    !> Correct the particle position if needed.
!    write(*,*) 'Position before correction: ', xp, ' velocity: ', up
    call p_MPIRemeshPosition ( xp )
    particleposition(:,np) = xp

!  Particle velocity stays the same - it is in the laboratory system.
!  (for use in moving system: particlevelocity(1,np) = up(1) + up(3))

!    write(*,*) 'Position after correction:  ', particleposition(:,np), &
!               ' velocity: ', particlevelocity(:,np)

  end do DOLOCAL

!  call PetscPrintf( PETSC_COMM_WORLD, 'Finished remeshing particles. \n', ierr )

end subroutine p_ControlRemeshParticles
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ControlRemeshInitial ( positionpointer, initialpositionpointer, arraynp, ierr )

  use g_parameters, only : MYRANK
  use p_mpi,        only : p_MPIRemeshInitial

  implicit none

  type(C_PTR),intent(in)                     :: positionpointer
  type(C_PTR),intent(in)                     :: initialpositionpointer
  integer,intent(in)                         :: arraynp
  PetscErrorCode,intent(inout)                      :: ierr
 
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
  real(kind=C_DOUBLE),pointer,dimension(:,:) :: initialposition
  real(kind=C_DOUBLE),dimension(1:3)         :: xp
  real(kind=C_DOUBLE),dimension(1:3)         :: xp0
  integer                                    :: np

  !> Neighbour that owns a particle, given by its local index from 1 to 8
  integer                             :: particleowner
  !> Location of the particle either in the centre or one of the halos
  integer                             :: particlelocation

  !> Associate position pointer array with target array.
  call c_f_pointer ( positionpointer, particleposition, [3, arraynp])
  !> Associate initial position pointer array with target array.
  call c_f_pointer ( initialpositionpointer, initialposition, [3, arraynp])

  DOLOCAL: do np = 1, arraynp, 1

    xp0 = initialposition(:,np)
    xp  = particleposition(:,np)

    !> Remesh particle position.
!    write(*,*) 'Position before remesh: ', xp
    call p_MPIRemeshInitial ( xp0, xp )
    initialposition(:,np) = xp0
!    write(*,*) 'Position after remesh:  ', particleposition(:,np) 

  end do DOLOCAL

!  call PetscPrintf( PETSC_COMM_WORLD, 'Finished remeshing initial particle positions. \n', ierr )

end subroutine p_ControlRemeshInitial
!----------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_ControlStats 
!> Compute real-space statistics locally on 2-D slices.
!> @param ierr should return 0
subroutine p_ControlStats ( tstep, ierr )

  use g_files,      only : FUNIT_PSTATS
  use p_arrays,     only : particle_dp, particle_up
  use g_parameters, only : MYRANK, MASTER, P_NP_GLOBAL, P_RHO
  use g_constants,  only : HALF

  implicit none

  integer,intent(in)               :: tstep
  PetscErrorCode,intent(inout)            :: ierr
  integer                          :: np
  real(kind=C_DOUBLE)              :: u1, u2, u3
  real(kind=C_DOUBLE)              :: dp
  real(kind=C_DOUBLE)              :: u1u1, u1u2, u1u3, &
                                      u2u2, u2u3, u3u3
  real(kind=C_DOUBLE)              :: rs_u1u1, rs_u1u2, rs_u1u3, &
                                      rs_u2u2, rs_u2u3, rs_u3u3
  real(kind=C_DOUBLE)              :: tke

  

  DOPART: do np = 1, np_centre, 1

    dp = particle_dp(1,np)

    !> Read local velocity
    u1 = particle_up(1,np)
    u2 = particle_up(2,np)
    u3 = particle_up(3,np)

    !> Compute Reynolds stresses
    u1u1 = u1u1 + u1 * u1
    u2u2 = u2u2 + u2 * u2
    u3u3 = u3u3 + u3 * u3
    u1u2 = u1u2 + u1 * u2
    u1u3 = u1u3 + u1 * u3
    u2u3 = u2u3 + u2 * u3

  end do DOPART

  !> Determine Reynolds stresses
  call MPI_Allreduce ( u1u1, rs_u1u1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2, rs_u2u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3, rs_u3u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2, rs_u1u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u3, rs_u1u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u3, rs_u2u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

  rs_u1u1 = rs_u1u1 / real ( P_NP_GLOBAL, kind=C_DOUBLE)
  rs_u2u2 = rs_u2u2 / real ( P_NP_GLOBAL, kind=C_DOUBLE)
  rs_u3u3 = rs_u3u3 / real ( P_NP_GLOBAL, kind=C_DOUBLE)
  rs_u1u2 = rs_u1u2 / real ( P_NP_GLOBAL, kind=C_DOUBLE)
  rs_u1u3 = rs_u1u3 / real ( P_NP_GLOBAL, kind=C_DOUBLE)
  rs_u2u3 = rs_u2u3 / real ( P_NP_GLOBAL, kind=C_DOUBLE)

  !> Compute turbulent kinetic energy
  tke = HALF * ( rs_u1u1 + rs_u2u2 + rs_u3u3 )

  IFMASTER: if ( MYRANK == MASTER ) then
    write(FUNIT_PSTATS,50) tstep, rs_u1u1, rs_u2u2, rs_u3u3, rs_u1u2, rs_u1u3, rs_u2u3, tke
  end if IFMASTER

50  FORMAT('ts: ',i6,' u1u1: ',e14.6,' u2u2: ',e14.6,' u3u3: ',e14.6, &
           ' u1u2: ',e14.6,' u1u3: ',e14.6,' u2u3: ',e14.6, ' tke: ',e14.6)

end subroutine p_ControlStats
!---------------------------------------------------------------------------------



!----------------------------------------------------------------
!subroutine p_quants(ierr)

!use d_arrays, only : length_x,length_y,length_z

!implicit none

!PetscErrorCode, intent(inout) :: ierr

!    eul_nd = real(np_global,kind=C_DOUBLE) /  (length_x*length_y*length_z)

!end subroutine
!---------------------------------------------------------------------------------

end module p_control
