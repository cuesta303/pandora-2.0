!---------------------------------------------------------------------------------
! module p_arrays
!> Contains auxiliary routines for Lagrangian particles, in particular routines
!! for allocating and deallocating memory and for initialising particles.
module p_arrays

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
use g_petsc
!use petscsys
!use petscvec

use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>


private

real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_xp
type(C_PTR),public                                                  :: ptr_xp
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_xp0
type(C_PTR),public                                                  :: ptr_xp0
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_up
type(C_PTR),public                                                  :: ptr_up
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_rk3
type(C_PTR),public                                                  :: ptr_rk3
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_up_rk3
type(C_PTR),public                                                  :: ptr_up_rk3
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_dp
type(C_PTR),public                                                  :: ptr_dp
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_vol
type(C_PTR),public                                                  :: ptr_vol
real(kind=C_DOUBLE),save,allocatable,dimension(:,:),target,public   :: particle_charge
type(C_PTR),public                                                  :: ptr_charge
PetscInt,save,allocatable,dimension(:,:),target,public              :: particle_status
type(C_PTR),public                                                  :: ptr_status

integer,parameter,public                                            :: DOF_STATUSARRAY    = 3 

integer,parameter,public                                            :: ARRAYHANDLE_XP     = 1
integer,parameter,public                                            :: ARRAYHANDLE_XP0    = 2
integer,parameter,public                                            :: ARRAYHANDLE_UP     = 3
integer,parameter,public                                            :: ARRAYHANDLE_RK3    = 4
integer,parameter,public                                            :: ARRAYHANDLE_UP_RK3 = 5
integer,parameter,public                                            :: ARRAYHANDLE_DP     = 6
integer,parameter,public                                            :: ARRAYHANDLE_VOL    = 7
integer,parameter,public                                            :: ARRAYHANDLE_CHARGE = 8
integer,parameter,public                                            :: ARRAYHANDLE_STATUS = 9

Vec,public                                                          :: p_restart_xp1
Vec,public                                                          :: p_restart_xp2
Vec,public                                                          :: p_restart_xp3
Vec,public                                                          :: p_restart_up1
Vec,public                                                          :: p_restart_up2
Vec,public                                                          :: p_restart_up3
Vec,public                                                          :: p_restart_dp
Vec,public                                                          :: p_restart_status
Vec,public                                                          :: p_restart_xp01
Vec,public                                                          :: p_restart_xp02
Vec,public                                                          :: p_restart_xp03

real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_xp1
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_xp2
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_xp3
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_up1
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_up2
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_up3
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_dp
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_status
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_xp01
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_xp02
real(kind=C_DOUBLE),pointer,dimension(:),public                     :: arr_p_restart_xp03

public :: p_ArraysAllocate
public :: p_ArraysDeallocate
public :: p_ArraysAllocateRestartArrays

contains
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_ArraysAllocate
!> Allocate the array corresponding to arrayhandle.
!> @param arrayhandle identifier of the array to allocate
!> @param ierr should return 0
subroutine p_ArraysAllocate ( arrayhandle, sizelocal, ierr )

  use g_parameters, only : MYRANK, MASTER, V_ALLOCS, YES, NPROCS

  implicit none

  integer, intent(in)    :: arrayhandle
  integer, intent(in)    :: sizelocal
  PetscErrorCode, intent(inout) :: ierr

  select case ( arrayhandle )
    !-------------------------------
    case ( ARRAYHANDLE_XP )

      IFVERBOSEXP: if ( V_ALLOCS == 1 ) then
!        IFMASTERXP: if ( MYRANK == MASTER ) then
          write(*,10) MYRANK, sizelocal
!        end if IFMASTERXP
      end if IFVERBOSEXP

      ierr = 0
      IFXPNOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_xp(1:3,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_xp   = C_LOC(particle_xp)
      end if IFXPNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_XP0 )

      IFVERBOSEXP0: if ( V_ALLOCS == 1 ) then
!        IFMASTERXP0: if ( MYRANK == MASTER ) then
          write(*,10) MYRANK, sizelocal
!        end if IFMASTERXP0
      end if IFVERBOSEXP0

      ierr = 0
      IFXP0NOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_xp0(1:3,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_xp0   = C_LOC(particle_xp0)
      end if IFXP0NOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_UP )

      IFVERBOSEUP: if ( V_ALLOCS == 1 ) then
!        IFMASTERUP: if ( MYRANK == MASTER ) then
          write(*,30) MYRANK, sizelocal
!        end if IFMASTERUP
      end if IFVERBOSEUP

      ierr = 0
      IFUPNOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_up(1:3,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_up = C_LOC(particle_up)
      end if IFUPNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_RK3 )

      IFVERBOSERK3: if ( V_ALLOCS == 1 ) then
!        IFMASTERRK3: if ( MYRANK == MASTER ) then
          write(*,40) MYRANK, sizelocal
!        end if IFMASTERRK3
      end if IFVERBOSERK3

      ierr = 0
      IFRK3NOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_rk3(1:3,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_rk3 = C_LOC(particle_rk3)
      end if IFRK3NOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_UP_RK3 )

      IFVERBOSEUPRK3: if ( V_ALLOCS == 1 ) then
!        IFMASTERUP: if ( MYRANK == MASTER ) then
          write(*,30) MYRANK, sizelocal
!        end if IFMASTERUP
      end if IFVERBOSEUPRK3

      ierr = 0
      IFUPRK3NOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_up_rk3(1:3,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_up_rk3 = C_LOC(particle_up_rk3)
      end if IFUPRK3NOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_DP )

      IFVERBOSEDP: if ( V_ALLOCS == 1 ) then
!        IFMASTERDP: if ( MYRANK == MASTER ) then
          write(*,50) MYRANK, sizelocal
!        end if IFMASTERDP
      end if IFVERBOSEDP

      ierr = 0
      IFDPNOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_dp(1:1,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_dp   = C_LOC(particle_dp)
      end if IFDPNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_STATUS )

      IFVERBOSESTATUS: if ( V_ALLOCS == 1 ) then
!        IFMASTERSTATUS: if ( MYRANK == MASTER ) then
          write(*,60) MYRANK, sizelocal
!        end if IFMASTERSTATUS
      end if IFVERBOSESTATUS

      ierr = 0
      IFSTATUSNOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_status(1:DOF_STATUSARRAY,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_status = C_LOC(particle_status)
      end if IFSTATUSNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_VOL )

      IFVERBOSEVOL: if ( V_ALLOCS == 1 ) then
!        IFMASTERVOL: if ( MYRANK == MASTER ) then
          write(*,70) MYRANK, sizelocal
!        end if IFMASTERVOL
      end if IFVERBOSEVOL

      ierr = 0
      IFVOLNOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_vol(1:1,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_vol = C_LOC(particle_vol)
      end if IFVOLNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_CHARGE )

      IFVERBOSECHARGE: if ( V_ALLOCS == 1 ) then
!        IFMASTERCHARGE: if ( MYRANK == MASTER ) then
          write(*,80) MYRANK, sizelocal
!        end if IFMASTERCHARGE
      end if IFVERBOSECHARGE

      ierr = 0
      IFCHARGENOTEMPTY: if ( sizelocal > 0 ) then
        allocate ( particle_charge(1:1,1:sizelocal),stat=ierr)
!        CHKERRQ ( ierr )
        ptr_charge   = C_LOC(particle_charge)
      end if IFCHARGENOTEMPTY
  end select


10  format('Allocating particle position array on rank:   ', i5, 'Local number of particles: ', i16)
20  format('Allocating local fluid velocity array on rank:   ', i5, 'Local number of particles: ', i16)
30  format('Allocating particle velocity array on rank:   ', i5, 'Local number of particles: ', i16)
40  format('Allocating particle acceleration / Runge Kutta buffer array on rank:   ', i5, 'Local number of particles: ', i16)
50  format('Allocating particle diameter array on rank:   ', i5, 'Local number of particles: ', i16)
60  format('Allocating particle status array on rank:   ', i5, 'Local number of particles: ', i16)
70  format('Allocating vol_particle array on rank:   ', i5, 'Local number of particles: ', i16)
80  format('Allocating particle charge array on rank:   ', i5, 'Local number of particles: ', i16)

end subroutine p_ArraysAllocate
!----------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_ArraysDeallocate
!> Deallocate the array corresponding to arrayhandle.
!> @param arrayhandle identifier of the array to be deallocated.
!> @param ierr should return 0
subroutine p_ArraysDeallocate ( arrayhandle, ierr )

  use g_parameters, only : MYRANK, MASTER, V_ALLOCS, YES, NPROCS

  implicit none

  integer, intent(in)    :: arrayhandle
  PetscErrorCode, intent(inout) :: ierr

  select case ( arrayhandle )
    !-------------------------------
    case ( ARRAYHANDLE_XP )

      ierr = 0
      IFXPNOTEMPTY: if ( allocated ( particle_xp ) ) then
        deallocate ( particle_xp,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFXPNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_XP0 )

      ierr = 0
      IFXP0NOTEMPTY: if ( allocated ( particle_xp0 ) ) then
        deallocate ( particle_xp0,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFXP0NOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_UP )

      ierr = 0
      IFUPNOTEMPTY: if ( allocated ( particle_up ) ) then
        deallocate ( particle_up,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFUPNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_RK3 )

      ierr = 0
      IFRK3NOTEMPTY: if ( allocated ( particle_rk3 ) ) then
        deallocate ( particle_rk3,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFRK3NOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_UP_RK3 )

      ierr = 0
      IFUPRK3NOTEMPTY: if ( allocated ( particle_up_rk3 ) ) then
        deallocate ( particle_up_rk3,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFUPRK3NOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_DP )

      ierr = 0
      IFDPNOTEMPTY: if ( allocated ( particle_dp ) ) then
        deallocate ( particle_dp,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFDPNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_STATUS )

      ierr = 0
      IFSTATUSNOTEMPTY: if ( allocated ( particle_status ) ) then
        deallocate ( particle_status,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFSTATUSNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_VOL )

      ierr = 0
      IFVOLNOTEMPTY: if ( allocated ( particle_vol ) ) then
        deallocate ( particle_vol,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFVOLNOTEMPTY
    !-------------------------------
    case ( ARRAYHANDLE_CHARGE )

      ierr = 0
      IFCHARGENOTEMPTY: if ( allocated ( particle_charge ) ) then
        deallocate ( particle_charge,stat=ierr)
!        CHKERRQ ( ierr )
      end if IFCHARGENOTEMPTY
  end select

end subroutine p_ArraysDeallocate
!----------------------------------------------------------------

!----------------------------------------------------------------
subroutine p_ArraysAllocateRestartArrays ( localdimension, ierr )

  use g_parameters, only : P_NP_GLOBAL, P_KEEP_INITIAL_POSITION, YES

  implicit none

  integer,intent(in)    :: localdimension
  PetscErrorCode,intent(inout) :: ierr

!  integer              :: globaldimension
  integer               :: status = 0

  PetscInt              :: ld
  PetscInt              :: gd
  PetscErrorCode        :: perr

  ld = localdimension

  !> Create a parallel vector for one particle component
  call VecCreateMPI ( PETSC_COMM_WORLD, ld, PETSC_DETERMINE, p_restart_xp1, perr )

  !> Check if the global dimension is equal to the global number of particles
  call VecGetSize ( p_restart_xp1, gd, perr )
  IFPARTICLESLOST: if ( gd /= P_NP_GLOBAL ) then
    call PetscPrintf( PETSC_COMM_WORLD, 'Particles lost. Abort. \n', perr )
    write(*,*) 'globaldimension ', gd, 'P_NP_GLOBAL ', P_NP_GLOBAL
    call MPI_Abort ( PETSC_COMM_WORLD, status, ierr )
  end if IFPARTICLESLOST

  !> Duplicate the first vector for all other components
  call VecDuplicate ( p_restart_xp1, p_restart_xp2, perr )
  call VecDuplicate ( p_restart_xp1, p_restart_xp3, perr )
  call VecDuplicate ( p_restart_xp1, p_restart_up1, perr )
  call VecDuplicate ( p_restart_xp1, p_restart_up2, perr )
  call VecDuplicate ( p_restart_xp1, p_restart_up3, perr )
  call VecDuplicate ( p_restart_xp1, p_restart_dp, perr )
  call VecDuplicate ( p_restart_xp1, p_restart_status, perr )

  IFINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then 
    call VecDuplicate ( p_restart_xp1, p_restart_xp01, perr )
    call VecDuplicate ( p_restart_xp1, p_restart_xp02, perr )
    call VecDuplicate ( p_restart_xp1, p_restart_xp03, perr )
  end if IFINITIAL

end subroutine p_ArraysAllocateRestartArrays
!----------------------------------------------------------------


!----------------------------------------------------------------
subroutine p_ArraysDeallocateRestartArrays ( localdimension, ierr )

  use g_parameters, only : P_NP_GLOBAL, P_KEEP_INITIAL_POSITION, YES

  implicit none

  integer,intent(in)    :: localdimension
  PetscErrorCode,intent(inout) :: ierr

  integer               :: status = 0

  PetscErrorCode        :: perr

  !> Destroy all parallel vectors
  call VecDestroy ( p_restart_xp1, perr )
  call VecDestroy ( p_restart_xp2, perr )
  call VecDestroy ( p_restart_xp3, perr )
  call VecDestroy ( p_restart_up1, perr )
  call VecDestroy ( p_restart_up2, perr )
  call VecDestroy ( p_restart_up3, perr )
  call VecDestroy ( p_restart_dp, perr )
  call VecDestroy ( p_restart_status, perr )

  IFINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then 
    call VecDestroy ( p_restart_xp01, perr )
    call VecDestroy ( p_restart_xp02, perr )
    call VecDestroy ( p_restart_xp03, perr )
  end if IFINITIAL

end subroutine p_ArraysDeallocateRestartArrays
!----------------------------------------------------------------

!---------------------------------------------------------------------------------

end module p_arrays
