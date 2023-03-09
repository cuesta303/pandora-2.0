!---------------------------------------------------------------------------------
! module tf_checkarrays
!> Contains tests to ensure that all Eulerian arrays work properly.
module tf_checkarrays
 
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

!------- data section begins ----------------------
implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

public :: tf_CheckArraysInit
public :: tf_CheckArraysDimension
public :: tf_CheckVecGetArray

!---------- data section ends ----------------------
contains


!---------------------------------------------------------------------------------
! subroutine tf_CheckArraysInit 
!> Simple test to see if all arrays are initialised to zero. 
!> @param ierr should return 0
subroutine tf_CheckArraysInit ( ierr )

  use g_constants, only : PETSCONE
  use f_arrays,    only : u1_3w, u2_3w, u3_3w, u1_1r2w, u2_1r2w, u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          u1_2w, u2_2w, u3_2w, u4_2w, u5_2w, u6_2w, &
                          u1_2r, u2_2r, u3_2r, u4_2r, u5_2r, u6_2r
  use g_parameters, only: P_TWO_WAY, YES

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  PetscScalar           :: arraysum
 
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking if all arrays are initialised to zero. \n', ierr )
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_3w... \n', ierr )
  call VecSum ( u1_3w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_3w... \n', ierr )
  call VecSum ( u2_3w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_3w... \n', ierr )
  call VecSum ( u3_3w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_1r2w... \n', ierr )
  call VecSum ( u1_1r2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_1r2w... \n', ierr )
  call VecSum ( u2_1r2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_1r2w... \n', ierr )
  call VecSum ( u3_1r2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_1r2w... \n', ierr )
  call VecSum ( u4_1r2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_1r2w... \n', ierr )
  call VecSum ( u5_1r2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_1r2w... \n', ierr )
  call VecSum ( u6_1r2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_1... \n', ierr )
  call VecSum ( rk3_3w_1, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_2... \n', ierr )
  call VecSum ( rk3_3w_2, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_3... \n', ierr )
  call VecSum ( rk3_3w_3, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  IFTWOWAY: if ( P_TWO_WAY == YES ) then
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_2w... \n', ierr )
    call VecSum ( u1_2w, arraysum, ierr )
!    call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
    call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_2w... \n', ierr )
    call VecSum ( u2_2w, arraysum, ierr )
!    call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
    call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_2w... \n', ierr )
    call VecSum ( u3_2w, arraysum, ierr )
!    call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
    call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  end if IFTWOWAY

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_2w... \n', ierr )
  call VecSum ( u4_2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_2w... \n', ierr )
  call VecSum ( u5_2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_2w... \n', ierr )
  call VecSum ( u6_2w, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_2r... \n', ierr )
  call VecSum ( u1_2r, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_2r... \n', ierr )
  call VecSum ( u2_2r, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_2r... \n', ierr )
  call VecSum ( u3_2r, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_2r... \n', ierr )
  call VecSum ( u4_2r, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_2r... \n', ierr )
  call VecSum ( u5_2r, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_2r... \n', ierr )
  call VecSum ( u6_2r, arraysum, ierr )
!  call PetscScalarView ( 1, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, arraysum, PETSC_VIEWER_STDOUT_WORLD, ierr )

end subroutine tf_CheckArraysInit
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_CheckArraysDimension 
!> Simple test to see if all arrays have the correct dimensions. 
!> @param ierr should return 0
subroutine tf_CheckArraysDimension ( ierr )

  use g_parameters, only : NODES_X, NODES_Y, NODES_Z, P_TWO_WAY, YES
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK

  use f_arrays, only :    u1_3w, u2_3w, u3_3w, u1_1r2w, u2_1r2w, u3_1r2w, u4_1r2w, u5_1r2w, u6_1r2w, &
                          rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          u1_2w, u2_2w, u3_2w, u4_2w, u5_2w, u6_2w, &
                          u1_2r, u2_2r, u3_2r, u4_2r, u5_2r, u6_2r

  use g_constants, only : ONE, TWO, PETSCONE

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  PetscScalar           :: arraysum, dimensionerror
 
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking if all arrays have the correct dimensions. \n', ierr )
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Result should be zero if correct dimensions. \n', ierr )
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_3w... \n', ierr )
  call VecSet ( u1_3w, ONE, ierr )
  call VecSum ( u1_3w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_3w... \n', ierr )
  call VecSet ( u2_3w, ONE, ierr )
  call VecSum ( u2_3w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_3w... \n', ierr )
  call VecSet ( u3_3w, ONE, ierr )
  call VecSum ( u3_3w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_1r2w... \n', ierr )
  call VecSet ( u1_1r2w, ONE, ierr )
  call VecSum ( u1_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_1r2w... \n', ierr )
  call VecSet ( u2_1r2w, ONE, ierr )
  call VecSum ( u2_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_1r2w... \n', ierr )
  call VecSet ( u3_1r2w, ONE, ierr )
  call VecSum ( u3_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_1r2w... \n', ierr )
  call VecSet ( u4_1r2w, ONE, ierr )
  call VecSum ( u4_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_1r2w... \n', ierr )
  call VecSet ( u5_1r2w, ONE, ierr )
  call VecSum ( u5_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_1r2w... \n', ierr )
  call VecSet ( u6_1r2w, ONE, ierr )
  call VecSum ( u6_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_1... \n', ierr )
  call VecSet ( rk3_3w_1, ONE, ierr )
  call VecSum ( rk3_3w_1, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_2... \n', ierr )
  call VecSet ( rk3_3w_2, ONE, ierr )
  call VecSum ( rk3_3w_2, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_3... \n', ierr )
  call VecSet ( rk3_3w_3, ONE, ierr )
  call VecSum ( rk3_3w_3, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  IFTWOWAY: if ( P_TWO_WAY == YES ) then
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_2w... \n', ierr )
    call VecSet ( u1_2w, ONE, ierr )
    call VecSum ( u1_2w, arraysum, ierr )
    dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
!    call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_2w... \n', ierr )
    call VecSet ( u2_2w, ONE, ierr )
    call VecSum ( u2_2w, arraysum, ierr )
    dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
!    call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
    call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_2w... \n', ierr )
    call VecSet ( u3_2w, ONE, ierr )
    call VecSum ( u3_2w, arraysum, ierr )
    dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
!    call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
    call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  end if IFTWOWAY

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_2w... \n', ierr )
  call VecSet ( u4_2w, ONE, ierr )
  call VecSum ( u4_2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_2w... \n', ierr )
  call VecSet ( u5_2w, ONE, ierr )
  call VecSum ( u5_2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_2w... \n', ierr )
  call VecSet ( u6_2w, ONE, ierr )
  call VecSum ( u6_2w, arraysum, ierr )
  dimensionerror = arraysum - real ( TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_2r... \n', ierr )
  call VecSet ( u1_2r, ONE, ierr )
  call VecSum ( u1_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( NODES_X * NODES_Y, kind=C_DOUBLE )
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_2r... \n', ierr )
  call VecSet ( u2_2r, ONE, ierr )
  call VecSum ( u2_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( NODES_X * NODES_Y, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_2r... \n', ierr )
  call VecSet ( u3_2r, ONE, ierr )
  call VecSum ( u3_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( NODES_X * NODES_Y, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_2r... \n', ierr )
  call VecSet ( u4_2r, ONE, ierr )
  call VecSum ( u4_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( NODES_X * NODES_Y, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_2r... \n', ierr )
  call VecSet ( u5_2r, ONE, ierr )
  call VecSum ( u5_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( NODES_X * NODES_Y, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_2r... \n', ierr )
  call VecSet ( u6_2r, ONE, ierr )
  call VecSum ( u6_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( NODES_X * NODES_Y, kind=C_DOUBLE ) 
!  call PetscScalarView ( 1, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

end subroutine tf_CheckArraysDimension
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_CheckVecGetArray
!> Simple test to see if all arrays have the correct dimensions. 
!> @param ierr should return 0
subroutine tf_CheckVecGetArray ( ierr )

  use g_parameters, only : NODES_X, NODES_Y, NODES_Z, P_TWO_WAY, YES
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK

  use f_arrays, only :    u1_3w, u2_3w, u3_3w, u1_1r2w, u2_1r2w, u3_1r2w, u4_1r2w, u5_1r2w, u6_1r2w, &
                          rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          u1_2w, u2_2w, u3_2w, u4_2w, u5_2w, u6_2w, &
                          u1_2r, u2_2r, u3_2r, u4_2r, u5_2r, u6_2r, &
                          da3w, da1r2w, da2w, da2r, &
                          i_min_3w, j_min_3w, k_min_3w, &
                          i_max_3w, j_max_3w, k_max_3w, &
                          i_min_1r2w, j_min_1r2w, &
                          i_max_1r2w, j_max_1r2w, &
                          z_min, z_max, y_min_2r, &
                          y_max_2r, x_min_2r, x_max_2r, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, & 
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, & 
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r

  use g_constants, only : HALF, ONE, TWO, PETSCONE

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  PetscScalar           :: arraysum, dimensionerror 
  integer               :: i, j, k, x, y, z

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking if pointer arrays are correctly obtained from & 
                     global vector. \n', ierr )
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Result should be zero if correct. \n', ierr )


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_3w... \n', ierr )
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )

  DOU1I3W: do i = i_min_3w, i_max_3w, 1
    DOU1J3W: do j = j_min_3w, j_max_3w, 1
      DOU1K3W: do k = k_min_3w, k_max_3w, 1

        arr_u1_3w(:,k,j,i) = HALF

      end do DOU1K3W
    end do DOU1J3W
  end do DOU1I3W

  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )

  call VecSum ( u1_3w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF * TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_3w... \n', ierr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )

  DOU2I3W: do i = i_min_3w, i_max_3w, 1
    DOU2J3W: do j = j_min_3w, j_max_3w, 1
      DOU2K3W: do k = k_min_3w, k_max_3w, 1

        arr_u2_3w(:,k,j,i) = HALF

      end do DOU2K3W
    end do DOU2J3W
  end do DOU2I3W

  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )

  call VecSum ( u2_3w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF * TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_3w... \n', ierr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

  DOU3I3W: do i = i_min_3w, i_max_3w, 1
    DOU3J3W: do j = j_min_3w, j_max_3w, 1
      DOU3K3W: do k = k_min_3w, k_max_3w, 1

        arr_u3_3w(:,k,j,i) = HALF

      end do DOU3K3W
    end do DOU3J3W
  end do DOU3I3W

  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

  call VecSum ( u3_3w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF * TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_1r2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, ierr )
  DOU1I1R2W: do i = i_min_1r2w, i_max_1r2w, 1
    DOU1J1R2W: do z = z_min, z_max, 1
      DOU1K1R2W: do j = j_min_1r2w, j_max_1r2w, 1

        arr_u1_1r2w(:,j,z,i) = HALF

      end do DOU1K1R2W
    end do DOU1J1R2W
  end do DOU1I1R2W

  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, ierr )
  call VecSum ( u1_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_1r2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, ierr )
  DOU2I1R2W: do i = i_min_1r2w, i_max_1r2w, 1
    DOU2J1R2W: do z = z_min, z_max, 1
      DOU2K1R2W: do j = j_min_1r2w, j_max_1r2w, 1

        arr_u2_1r2w(:,j,z,i) = HALF

      end do DOU2K1R2W
    end do DOU2J1R2W
  end do DOU2I1R2W
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, ierr )
  call VecSum ( u2_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_1r2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, ierr )
  DOU3I1R2W: do i = i_min_1r2w, i_max_1r2w, 1
    DOU3J1R2W: do z = z_min, z_max, 1
      DOU3K1R2W: do j = j_min_1r2w, j_max_1r2w, 1

        arr_u3_1r2w(:,j,z,i) = HALF

      end do DOU3K1R2W
    end do DOU3J1R2W
  end do DOU3I1R2W
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, ierr )
  call VecSum ( u3_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_1r2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, ierr )
  DOU4I1R2W: do i = i_min_1r2w, i_max_1r2w, 1
    DOU4J1R2W: do z = z_min, z_max, 1
      DOU4K1R2W: do j = j_min_1r2w, j_max_1r2w, 1

        arr_u4_1r2w(:,j,z,i) = HALF

      end do DOU4K1R2W
    end do DOU4J1R2W
  end do DOU4I1R2W
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, ierr )
  call VecSum ( u4_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_1r2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, ierr )
  DOU5I1R2W: do i = i_min_1r2w, i_max_1r2w, 1
    DOU5J1R2W: do z = z_min, z_max, 1
      DOU5K1R2W: do j = j_min_1r2w, j_max_1r2w, 1

        arr_u5_1r2w(:,j,z,i) = HALF

      end do DOU5K1R2W
    end do DOU5J1R2W
  end do DOU5I1R2W
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, ierr )
  call VecSum ( u5_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_1r2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, ierr )
  DOU6I1R2W: do i = i_min_1r2w, i_max_1r2w, 1
    DOU6J1R2W: do z = z_min, z_max, 1
      DOU6K1R2W: do j = j_min_1r2w, j_max_1r2w, 1

        arr_u6_1r2w(:,j,z,i) = HALF

      end do DOU6K1R2W
    end do DOU6J1R2W
  end do DOU6I1R2W
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, ierr )
  call VecSum ( u6_1r2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KJ * NODES_Z * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_1... \n', ierr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, ierr )
  DORK31I3W: do i = i_min_3w, i_max_3w, 1
    DORK31J3W: do j = j_min_3w, j_max_3w, 1
      DORK31K3W: do k = k_min_3w, k_max_3w, 1

        arr_rk3_3w_1(:,k,j,i) = HALF

      end do DORK31K3W
    end do DORK31J3W
  end do DORK31I3W
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, ierr )
  call VecSum ( rk3_3w_1, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_2... \n', ierr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, ierr )
  DORK32I3W: do i = i_min_3w, i_max_3w, 1
    DORK32J3W: do j = j_min_3w, j_max_3w, 1
      DORK32K3W: do k = k_min_3w, k_max_3w, 1

        arr_rk3_3w_2(:,k,j,i) = HALF

      end do DORK32K3W
    end do DORK32J3W
  end do DORK32I3W
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, ierr )
  call VecSum ( rk3_3w_2, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking rk3_3w_3... \n', ierr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, ierr )
  DORK33I3W: do i = i_min_3w, i_max_3w, 1
    DORK33J3W: do j = j_min_3w, j_max_3w, 1
      DORK33K3W: do k = k_min_3w, k_max_3w, 1

        arr_rk3_3w_3(:,k,j,i) = HALF

      end do DORK33K3W
    end do DORK33J3W
  end do DORK33I3W
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, ierr )
  call VecSum ( rk3_3w_3, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ * NODES_KI, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  IFTWOWAY: if ( P_TWO_WAY == YES ) then
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_2w... \n', ierr )
    call DMDAVecGetArrayF90 ( da2w, u1_2w, arr_u1_2w, ierr )
    DOU1J2W: do j = j_min_3w, j_max_3w, 1
      DOU1K2W: do k = k_min_3w, k_max_3w, 1

        arr_u1_2w(:,k,j) = HALF

      end do DOU1K2W
    end do DOU1J2W
    call DMDAVecRestoreArrayF90 ( da2w, u1_2w, arr_u1_2w, ierr )
    call VecSum ( u1_2w, arraysum, ierr )
    dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
    call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_2w... \n', ierr )
    call DMDAVecGetArrayF90 ( da2w, u2_2w, arr_u2_2w, ierr )
    DOU2J2W: do j = j_min_3w, j_max_3w, 1
      DOU2K2W: do k = k_min_3w, k_max_3w, 1

        arr_u2_2w(:,k,j) = HALF

      end do DOU2K2W
    end do DOU2J2W
    call DMDAVecRestoreArrayF90 ( da2w, u2_2w, arr_u2_2w, ierr )
    call VecSum ( u2_2w, arraysum, ierr )
    dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
    call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_2w... \n', ierr )
    call DMDAVecGetArrayF90 ( da2w, u3_2w, arr_u3_2w, ierr )
    DOU3J2W: do j = j_min_3w, j_max_3w, 1
      DOU3K2W: do k = k_min_3w, k_max_3w, 1

        arr_u3_2w(:,k,j) = HALF

      end do DOU3K2W
    end do DOU3J2W
    call DMDAVecRestoreArrayF90 ( da2w, u3_2w, arr_u3_2w, ierr )
    call VecSum ( u3_2w, arraysum, ierr )
    dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
    call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )
  end if IFTWOWAY

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, ierr )
  DOU4J2W: do j = j_min_3w, j_max_3w, 1
    DOU4K2W: do k = k_min_3w, k_max_3w, 1

      arr_u4_2w(:,k,j) = HALF

    end do DOU4K2W
  end do DOU4J2W
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, ierr )
  call VecSum ( u4_2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, ierr )
  DOU5J2W: do j = j_min_3w, j_max_3w, 1
    DOU5K2W: do k = k_min_3w, k_max_3w, 1

      arr_u5_2w(:,k,j) = HALF

    end do DOU5K2W
  end do DOU5J2W
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, ierr )
  call VecSum ( u5_2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_2w... \n', ierr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, ierr )
  DOU6J2W: do j = j_min_3w, j_max_3w, 1
    DOU6K2W: do k = k_min_3w, k_max_3w, 1

      arr_u6_2w(:,k,j) = HALF

    end do DOU6K2W
  end do DOU6J2W
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, ierr )
  call VecSum ( u6_2w, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  TWO * NODES_KK * NODES_KJ, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u1_2r... \n', ierr )
  call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, ierr )
  DOU1Y2R: do y = y_min_2r, y_max_2r, 1
    DOU1X2R: do x = x_min_2r, x_max_2r, 1

      arr_u1_2r(x,y) = HALF

    end do DOU1X2R
  end do DOU1Y2R
  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, ierr )
  call VecSum ( u1_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  NODES_X * NODES_Y, kind=C_DOUBLE )
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u2_2r... \n', ierr )
  call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, ierr )
  DOU2Y2R: do y = y_min_2r, y_max_2r, 1
    DOU2X2R: do x = x_min_2r, x_max_2r, 1

      arr_u2_2r(x,y) = HALF

    end do DOU2X2R
  end do DOU2Y2R
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, ierr )
  call VecSum ( u2_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  NODES_X * NODES_Y, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u3_2r... \n', ierr )
  call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, ierr )
  DOU3Y2R: do y = y_min_2r, y_max_2r, 1
    DOU3X2R: do x = x_min_2r, x_max_2r, 1

      arr_u3_2r(x,y) = HALF

    end do DOU3X2R
  end do DOU3Y2R
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, ierr )
  call VecSum ( u3_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  NODES_X * NODES_Y, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u4_2r... \n', ierr )
  call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, ierr )
  DOU4Y2R: do y = y_min_2r, y_max_2r, 1
    DOU4X2R: do x = x_min_2r, x_max_2r, 1

      arr_u4_2r(x,y) = HALF

    end do DOU4X2R
  end do DOU4Y2R
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, ierr )
  call VecSum ( u4_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  NODES_X * NODES_Y, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u5_2r... \n', ierr )
  call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, ierr )
  DOU5Y2R: do y = y_min_2r, y_max_2r, 1
    DOU5X2R: do x = x_min_2r, x_max_2r, 1

      arr_u5_2r(x,y) = HALF

    end do DOU5X2R
  end do DOU5Y2R
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, ierr )
  call VecSum ( u5_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  NODES_X * NODES_Y, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )

  call PetscPrintf ( PETSC_COMM_WORLD, '    | Checking u6_2r... \n', ierr )
  call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, ierr )
  DOU6Y2R: do y = y_min_2r, y_max_2r, 1
    DOU6X2R: do x = x_min_2r, x_max_2r, 1

      arr_u6_2r(x,y) = HALF

    end do DOU6X2R
  end do DOU6Y2R
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, ierr )
  call VecSum ( u6_2r, arraysum, ierr )
  dimensionerror = arraysum - real ( HALF *  NODES_X * NODES_Y, kind=C_DOUBLE ) 
  call PetscScalarView ( PETSCONE, dimensionerror, PETSC_VIEWER_STDOUT_WORLD, ierr )


end subroutine tf_CheckVecGetArray
!---------------------------------------------------------------------------------

end module tf_checkarrays
