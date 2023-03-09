!---------------------------------------------------------------------------------
! module f_force
!> Compute forcing for isotropic turbulence
module f_force

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscdm.h>
!#include <petsc/finclude/petscdmda.h>
use g_petsc
!use petscsys
!use petscdm
!use petscdmda

use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

!---------------------------------------------------------------------
! force array 
!---------------------------------------------------------------------
!> i component of the forcing array
real(kind=C_DOUBLE),allocatable,save,public,dimension(:,:,:,:)     :: fluid_f1
!> j component of the forcing array
real(kind=C_DOUBLE),allocatable,save,public,dimension(:,:,:,:)     :: fluid_f2
!> k component of the forcing array
real(kind=C_DOUBLE),allocatable,save,public,dimension(:,:,:,:)     :: fluid_f3
!> number of forced wavemodes
integer,save,public                                                :: FORCE_NK

!---------------------------------------------------------------------
! forcing parameters 
!---------------------------------------------------------------------


public :: f_ForceInit
public :: f_ForceCalc
public :: f_ForceFinalise

contains
!---------------------------------------------------------------------------------
! subroutine f_ForceInit
!> Initialise forcing.
!> @param ierr should return 0
subroutine f_ForceInit ( ierr )

  use g_parameters, only : F_INIT_TYPE, F_INIT_RESTART, READ_FILE

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising f_force module \n', ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking forcing parameters \n', ierr )
  call f_ForceCheck(ierr)

  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating forcing arrays \n', ierr )
  call f_ForceAlloc(ierr)

  IFRESTART: if ( F_INIT_TYPE == F_INIT_RESTART ) then
    call PetscPrintf( PETSC_COMM_WORLD, '    | Initialising forcing arrays from file \n', ierr )
    call f_ForceRestart ( READ_FILE, ierr )

  end if IFRESTART

  !> @todo Initialise random seed.

end subroutine f_ForceInit
!---------------------------------------------------------------------------------
! subroutine f_ForceFinalise
!> Finalise forcing.
!> @param ierr should return 0
subroutine f_ForceFinalise ( ierr )

  use g_parameters, only : WRITE_FILE

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Finalising f_force module \n', ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Writing forcing arrays to file \n', ierr )
  call f_ForceRestart ( WRITE_FILE, ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating forcing arrays \n', ierr )
  call f_ForceDealloc(ierr)

end subroutine f_ForceFinalise
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------
! subroutine f_ForceCalc
!> Compute forcing on all forced wave modes and ensure conjugate 
!! symmetry on the i=0 plane.
!> @param tstep current time step @todo is the time step needed here?
!> @param dt step size $\Delta t$
!> @param ierr should return 0
subroutine f_ForceCalc ( tstep, dt, ierr)

  use g_parameters, only : MASTER, MYRANK, NODES_X, NODES_Y, NODES_Z, YES, &
                           FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP, REPORT
  use g_constants,  only : ZERO, ONE, TWO
  use g_random,     only : g_RandomNormal

  implicit none


  PetscErrorCode,intent(inout)                 :: ierr
  integer,intent(in)                    :: tstep
  real(kind=C_DOUBLE),intent(in)        :: dt
  real(kind=C_DOUBLE)                   :: real1, real2, real3
  real(kind=C_DOUBLE)                   :: imag1, imag2, imag3
  real(kind=C_DOUBLE)                   :: k1,k2,k3
  integer                               :: i,j,k,counter
  integer                               :: ii,jj,kk
  integer                               :: finished
  real(kind=C_DOUBLE),dimension(0:1)    :: comp_ran1
  real(kind=C_DOUBLE),dimension(0:1)    :: comp_ran2
  real(kind=C_DOUBLE),dimension(0:1)    :: comp_ran3
  real(kind=C_DOUBLE)                   :: dt_over_TF

  IFDECAY: if ( tstep >= FORCE_STOP ) then
    fluid_f1  = ZERO
    fluid_f2  = ZERO
    fluid_f3  = ZERO
  else 

    !> Make sure dt/TF is well within the (0,1) interval
    dt_over_TF = dt / FORCE_TF
  !  IFRANDOM: if ( dt_over_TF > 0.95 ) then
  !    dt_over_TF = 0.95
  !  end if IFRANDOM
  !  IFCORR: if ( dt_over_TF < 0.05 ) then
  !    dt_over_TF = 0.05
  !  end if IFCORR

    !> ensure conjugate symmetric forcing on i=0

    DOJZERO: do j=0,FORCE_NK,1
      DOKZERO: do k=0,FORCE_NK,1

        !> get normally distributed random variables
        call g_RandomNormal ( ZERO, ONE, real1, imag1, ierr )
        call g_RandomNormal ( ZERO, ONE, real2, imag2, ierr )
        call g_RandomNormal ( ZERO, ONE, real3, imag3, ierr )

        !> calculate new force according to uhlenbeck-ornstein process
        fluid_f1(0,k,j,0) = fluid_f1(0,k,j,0) * ( one - dt_over_TF ) +      &
                            real1 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f1(1,k,j,0) = fluid_f1(1,k,j,0) * ( one - dt_over_TF ) +      &
                            imag1 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f2(0,k,j,0) = fluid_f2(0,k,j,0) * ( one - dt_over_TF ) +           &
                            real2 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f2(1,k,j,0) = fluid_f2(1,k,j,0) * ( one - dt_over_TF ) +           &
                            imag2 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f3(0,k,j,0) = fluid_f3(0,k,j,0) * ( one - dt_over_TF ) +           &
                            real3 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f3(1,k,j,0) = fluid_f3(1,k,j,0) * ( one - dt_over_TF ) +           &
                            imag3 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
    
        !> match conjugate symmetry on i=0 plane
        fluid_f1(0,-k,-j,0) = fluid_f1(0,k,j,0)
        fluid_f1(1,-k,-j,0) = -fluid_f1(1,k,j,0)
        fluid_f2(0,-k,-j,0) = fluid_f2(0,k,j,0)
        fluid_f2(1,-k,-j,0) = -fluid_f2(1,k,j,0)
        fluid_f3(0,-k,-j,0) = fluid_f3(0,k,j,0)
        fluid_f3(1,-k,-j,0) = -fluid_f3(1,k,j,0)

!        write(*,*) 0, j, k, real1, imag1, real2, imag2, real3, imag3

        IFJGT0: if ( j > 0 ) then

          !> get normally distributed random variables
          call g_RandomNormal ( ZERO, ONE, real1, imag1, ierr )
          call g_RandomNormal ( ZERO, ONE, real2, imag2, ierr )
          call g_RandomNormal ( ZERO, ONE, real3, imag3, ierr )

          !> calculate new force on (0, -j, k) mode
          fluid_f1(0,k,-j,0) = fluid_f1(0,k,-j,0) * ( one - dt_over_TF ) +      &
                               real1 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
          fluid_f1(1,k,-j,0) = fluid_f1(1,k,-j,0) * ( one - dt_over_TF ) +      &
                               imag1 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
          fluid_f2(0,k,-j,0) = fluid_f2(0,k,-j,0) * ( one - dt_over_TF ) +           &
                               real2 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
          fluid_f2(1,k,-j,0) = fluid_f2(1,k,-j,0) * ( one - dt_over_TF ) +           &
                               imag2 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
          fluid_f3(0,k,-j,0) = fluid_f3(0,k,-j,0) * ( one - dt_over_TF ) +           &
                               real3 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
          fluid_f3(1,k,-j,0) = fluid_f3(1,k,-j,0) * ( one - dt_over_TF ) +           &
                               imag3 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )

          !> Match conjugate symmetry
          fluid_f1(0,-k,j,0) = fluid_f1(0,k,-j,0)
          fluid_f1(1,-k,j,0) = -fluid_f1(1,k,-j,0)
          fluid_f2(0,-k,j,0) = fluid_f2(0,k,-j,0)
          fluid_f2(1,-k,j,0) = -fluid_f2(1,k,-j,0)
          fluid_f3(0,-k,j,0) = fluid_f3(0,k,-j,0)
          fluid_f3(1,-k,j,0) = -fluid_f3(1,k,-j,0)

        end if IFJGT0

!        write(*,*) 0, j, k, real1, imag1, real2, imag2, real3, imag3
!        write(*,*) 0, j, k, fluid_f1(:,-k,j,0), fluid_f1(:,k,-j,0), fluid_f1(:,k,j,0), fluid_f1(:,-k,-j,0)
!        write(*,*) 0, j, k, fluid_f2(:,-k,j,0), fluid_f2(:,k,-j,0), fluid_f2(:,k,j,0), fluid_f2(:,-k,-j,0)
!        write(*,*) 0, j, k, fluid_f3(:,-k,j,0), fluid_f3(:,k,-j,0), fluid_f3(:,k,j,0), fluid_f3(:,-k,-j,0)

!          write(*,20) 0,j,k, fluid_f(0,j,k,:)
!          write(*,20) 0,-j,-k, fluid_f(0,-j,-k,:)

      end do DOKZERO
    end do DOJZERO

    !> conjugate symmetry is guaranteed on all other modes     
    DOIOTHER: do i=1,FORCE_NK,1
      DOJOTHER: do j=(-FORCE_NK),FORCE_NK,1
        DOKOTHER: do k=(-FORCE_NK),FORCE_NK,1


        ! get normally distributed random variables
        call g_RandomNormal ( ZERO, ONE, real1, imag1, ierr )
        call g_RandomNormal ( ZERO, ONE, real2, imag2, ierr )
        call g_RandomNormal ( ZERO, ONE, real3, imag3, ierr )
        
        ! form complex random variables
!            comp_ran1 = cmplx ( r1, c1, kind=C_DOUBLE_COMPLEX )
!            comp_ran2 = cmplx ( r2, c2, kind=C_DOUBLE_COMPLEX )
!            comp_ran3 = cmplx ( r3, c3, kind=C_DOUBLE_COMPLEX )
  
        ! calculate new force according to uhlenbeck-ornstein process
        fluid_f1(0,k,j,i) = fluid_f1(0,k,j,i) * ( one - dt_over_TF ) +      &
                            real1 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f1(1,k,j,i) = fluid_f1(1,k,j,i) * ( one - dt_over_TF ) +      &
                            imag1 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f2(0,k,j,i) = fluid_f2(0,k,j,i) * ( one - dt_over_TF ) +           &
                            real2 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f2(1,k,j,i) = fluid_f2(1,k,j,i) * ( one - dt_over_TF ) +           &
                            imag2 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f3(0,k,j,i) = fluid_f3(0,k,j,i) * ( one - dt_over_TF ) +           &
                            real3 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
        fluid_f3(1,k,j,i) = fluid_f3(1,k,j,i) * ( one - dt_over_TF ) +           &
                            imag3 * sqrt( two * FORCE_SIGMA**2 * dt_over_TF )
    
!          write(*,20) i,j,k, fluid_f(i,j,k,:)

!        write(*,*) i, j, k, real1, imag1, real2, imag2, real3, imag3

        end do DOKOTHER
      end do DOJOTHER
    end do DOIOTHER

!  write(*,*) fluid_f1
!  write(*,*) fluid_f2
!  write(*,*) fluid_f3

  end if IFDECAY

20  format('    | ', I3, I3, I3, 6E14.6)

200 format(6e14.6)
end subroutine f_ForceCalc
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ForceAlloc
!> Allocate fluid forcing array depending on number of forced wave modes.
!> @param ierr should return 0
subroutine f_ForceAlloc ( ierr )

  use g_parameters, only : MYRANK, MASTER, NODES_X, NODES_Y, NODES_Z, &
                           FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP
  use g_constants,  only : ZERO
  use g_domain,     only : KI_MIN, KJ_MIN, KK_MIN, &
                           KI_MAX, KJ_MAX, KK_MAX

  implicit none

  PetscErrorCode,intent(inout)                 :: ierr

  integer                               :: alloc_stat

  FORCE_NK = aint ( FORCE_KF )
  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,10) FORCE_NK
  end if IFMASTER

  ! allocate force array
  allocate ( fluid_f1(0:1,(-FORCE_NK):FORCE_NK,(-FORCE_NK):FORCE_NK,0:FORCE_NK),stat=alloc_stat )
  allocate ( fluid_f2(0:1,(-FORCE_NK):FORCE_NK,(-FORCE_NK):FORCE_NK,0:FORCE_NK),stat=alloc_stat )
  allocate ( fluid_f3(0:1,(-FORCE_NK):FORCE_NK,(-FORCE_NK):FORCE_NK,0:FORCE_NK),stat=alloc_stat )
  fluid_f1 = ZERO
  fluid_f2 = ZERO
  fluid_f3 = ZERO

10  format('    | number of modes in forcing array: ', I5)
end subroutine f_ForceAlloc
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------
! subroutine f_ForceDealloc
!> Deallocate fluid forcing arrays.
!> @param ierr should return 0
subroutine f_ForceDealloc ( ierr )

  use g_parameters, only : MYRANK, MASTER, NODES_X, NODES_Y, NODES_Z, &
                           FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP

  implicit none

  PetscErrorCode,intent(inout)                 :: ierr

  integer                               :: alloc_stat

  ! deallocate force array
  deallocate ( fluid_f1, stat=alloc_stat )
  deallocate ( fluid_f2, stat=alloc_stat )
  deallocate ( fluid_f3, stat=alloc_stat )

end subroutine f_ForceDealloc
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------
subroutine f_ForceRestart ( readorwrite, ierr )

  use g_parameters, only : MYRANK, MASTER, READ_FILE, WRITE_FILE, &
                           F_FORCE_FILE_READ, F_FORCE_FILE_WRITE, FORCE_KF
  use g_files,      only : g_FilesForcing, FUNIT_FORCEW, FUNIT_FORCER

  implicit none

  integer,intent(in)              :: readorwrite
  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: i, j, k

    !> select read/write case
    select case(readorwrite)
      !> read force coefficients
      case(READ_FILE)
        
        !> Read in forced wavenumber range. 
!        read ( FUNIT_FORCER, 10 ) FORCE_KF 

        !> @todo Check if array size is correct.

        !> If correct array size, read forcing.
        DOIREAD: do i = 0, FORCE_NK
          DOJREAD: do j = (-FORCE_NK), FORCE_NK
            DOKREAD: do k = (-FORCE_NK), FORCE_NK
!              read ( FUNIT_FORCER, 20 ) fluid_f1(:,k,j,i), fluid_f2(:,k,j,i), fluid_f3(:,k,j,i)
            end do DOKREAD
          end do DOJREAD
        end do DOIREAD

      ! write force coefficients
      case(WRITE_FILE)

        !> Open forcing file.
        call g_FilesForcing ( ierr )

        !> @todo Write random seed to file.

        !> Write forced wavenumber range.
!        write ( FUNIT_FORCEW, 10 ) FORCE_KF 

        !> Write forcing array.
        DOIWRITE: do i = 0, FORCE_NK
          DOJWRITE: do j = (-FORCE_NK), FORCE_NK
            DOKWRITE: do k = (-FORCE_NK), FORCE_NK
!              write ( FUNIT_FORCEW, 20 ) fluid_f1(:,k,j,i), fluid_f2(:,k,j,i), fluid_f3(:,k,j,i)
            end do DOKWRITE
          end do DOJWRITE
        end do DOIWRITE

    end select

10  format(f15.8)
20  format(6f15.8)

end subroutine f_ForceRestart
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------
! subroutine f_ForceCheck
!> Compute forcing dissipation, Reynolds number and time scale based
!! on forcing parameters.
!> @param ierr should return 0
subroutine f_ForceCheck(ierr)

  use g_parameters, only : MYRANK, MASTER, F_NU, &
                           FORCE_TF, FORCE_SIGMA, FORCE_KF, FORCE_STOP
  use g_constants,  only : ONE, TWO, THREE, FOUR

  implicit none


  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: status = 0
  real(kind=C_DOUBLE)             :: re
  real(kind=C_DOUBLE)             :: t
  real(kind=C_DOUBLE)             :: epsilon
  real(kind=C_DOUBLE)             :: k0 = 1.0

  epsilon = FORCE_SIGMA * FORCE_SIGMA * FORCE_TF
  Re      = epsilon**(one/three) * k0**(-four/three) / F_NU
  t       = FORCE_TF * epsilon**(one/three) * k0**(two/three)

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,10) Re
    write(*,20) t
    write(*,30) epsilon
  end if IFMASTER

10  format('    | Re* = 'f10.5)
20  format('    | t*  = 'f10.5)
30  format('    | e*  = 'f10.5)

end subroutine f_ForceCheck
!---------------------------------------------------------------------------------

end module f_force
