!---------------------------------------------------------------------------------
! module tf_rk3
!> Contains all subroutines necessary to solve the fluid equations of motion.
module tf_rk3

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

use f_rk3, only: f_RK3PrepareRealSpace, f_RK3RealSpace, f_RK3TransformToRealSpace, &
                 f_RK3TransformAllArraysjiRs2Fs, f_RK3TransformAllArrayskRs2Fs

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

public :: tf_RK3Gauss
public :: tf_RK3TaylorGreen
public :: tf_RK3TaylorGreenReverse
public :: tf_RK3CheckNonlinear

contains
!---------------------------------------------------------------------------------
! subroutine tf_RK3Gauss 
!> Initialise velocity with a 3D Gaussian and check if it matches the analytical
!! solution in real space. Transform back to wave space and check if the original
!! field is reproduced. Also check if the imaginary part is zero as required.
!! Perform equivalent tests on the vorticity.
!> @param ierr should return 0
subroutine tf_RK3Gauss ( ierr )

  use f_arrays,     only: u1_3w, u1_1r2w, u1_2r, u1_2w, &
                          u2_3w, u2_1r2w, u2_2r, u2_2w, &
                          u3_3w, u3_1r2w, u3_2r, u3_2w

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  !> Initialise with a Gaussian field in real space
  call tf_RK3GaussRealSpace ( ierr )
  !> Transform to wave space
  call tf_RK3GaussTransformBack ( ierr )
  !> Compute vorticity and do a 1D transform in k direction
  call f_RK3PrepareRealSpace ( ierr )
  !> Transform back to real space. Check error in non-linear term, vorticity and velocity on each slice.
  call tf_RK3GaussRealSpace2 ( ierr )

end subroutine tf_RK3Gauss
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine tf_RK3GaussRealSpace 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3GaussRealSpace ( ierr )

  use f_arrays, only: da1r2w, da2r, z_min, z_max, &
                      u1_1r2w, u2_1r2w, u3_1r2w, &
                      u4_1r2w, u5_1r2w, u6_1r2w, &
                      arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                      arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                      u1_2r, u2_2r, u3_2r, &
                      u4_2r, u5_2r, u6_2r, &
                      arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                      arr_u4_2r, arr_u5_2r, arr_u6_2r
  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: zslice                 

  PetscErrorCode        :: perr

  !> Get write access to all velocity and vorticity components.
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )


  !> Loop through slices in z direction.
  DOSLICES: do zslice = z_min, z_max, 1
               
    call tf_RK3GaussInit ( zslice, ierr )

    !> Transform them back.
    call f_RK3TransformAllArraysjiRs2Fs ( zslice, ierr )

  end do DOSLICES

  !> Return velocity and vorticity components to PETSc.
  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

end subroutine tf_RK3GaussRealSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3GaussInit 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3GaussInit ( zslice, ierr )

  use g_constants,  only: ZERO, ONE, TWO, FOUR
  use g_parameters, only: NODES_X, NODES_Y, NODES_Z
  use f_arrays,     only: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          y_min_2r, y_max_2r, x_min_2r, x_max_2r

  implicit none
 
  integer,intent(in)    :: zslice                 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: i, j

  write(*,*) 'Initialising slice ', zslice

  !> Loop through slices in z direction.
  DOY: do j = y_min_2r, y_max_2r, 1
    DOX: do i = x_min_2r, x_max_2r, 1
      
      arr_u1_2r(i,j) = tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, ONE )
      arr_u2_2r(i,j) = tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, TWO )
      arr_u3_2r(i,j) = tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, FOUR )

    end do DOX
  end do DOY

end subroutine tf_RK3GaussInit
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine tf_RK3GaussTransformBack
!> @param ierr should return 0
subroutine tf_RK3GaussTransformBack ( ierr )

  use f_arrays,     only: da3w, da1r2w, da2w, i_min_3w, i_max_3w, &
                          u1_3w, u2_3w, u3_3w, &
                          rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          u1_1r2w, u2_1r2w, u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          u1_2w, u2_2w, u3_2w, &
                          u4_2w, u5_2w, u6_2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          k_min_3w, k_max_3w, j_min_3w, j_max_3w

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: islice, j, k                 
 
  PetscErrorCode        :: perr

!> Get read access to all non-linear components
  call DMDAVecGetArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> Get read access to all coupling components
  call DMDAVecGetArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Get write access to all RK3 components
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )


  !> Loop through slices in i direction
  DOSLICES: do islice = i_min_3w, i_max_3w, 1

    !> Transform non-linear components.
    call f_RK3TransformAllArrayskRs2Fs ( islice, ierr )

    arr_u1_3w(:,:,:,islice) = arr_u1_2w(:,:,:)
    arr_u2_3w(:,:,:,islice) = arr_u2_2w(:,:,:)
    arr_u3_3w(:,:,:,islice) = arr_u3_2w(:,:,:)
 
  end do DOSLICES
  !end loop


!> Return all non-linear components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> Return all coupling components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

!> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Return all RK3 components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )

end subroutine tf_RK3GaussTransformBack
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3GaussRealSpace2 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3GaussRealSpace2 ( ierr )

  use f_arrays, only: da1r2w, da2r, z_min, z_max, &
                      u1_1r2w, u2_1r2w, u3_1r2w, &
                      u4_1r2w, u5_1r2w, u6_1r2w, &
                      arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                      arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                      u1_2r, u2_2r, u3_2r, &
                      u4_2r, u5_2r, u6_2r, &
                      arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                      arr_u4_2r, arr_u5_2r, arr_u6_2r
  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: zslice                 

  real(kind=C_DOUBLE)   :: minimumerror, maximumerror
  integer               :: minerrloc, maxerrloc

  PetscErrorCode        :: perr

  !> Get write access to all velocity and vorticity components.
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )


  !> Loop through slices in z direction.
  DOSLICES: do zslice = z_min, z_max, 1
               
    call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

    !> Transform velocity and vorticity components to real space.
    call f_RK3TransformToRealSpace ( zslice, 2, 2, ierr )

    !> Check absolute error in velocity and vorticity.
    call tf_RK3GaussCheckAbsoluteError ( zslice, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Gauss test - viewing absolute error in u1 component. \n', perr )
    call VecMin ( u1_2r, minerrloc, minimumerror, perr )
    call VecMax ( u1_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u1 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u1 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Gauss test - viewing absolute error in u2 component. \n', perr )
    call VecMin ( u2_2r, minerrloc, minimumerror, perr )
    call VecMax ( u2_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u2 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u2 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Gauss test - viewing absolute error in u3 component. \n', perr )
    call VecMin ( u3_2r, minerrloc, minimumerror, perr )
    call VecMax ( u3_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u3 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u3 component: ', maximumerror, ' at index ', maxerrloc

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Gauss test - viewing absolute error in u4 component. \n', perr )
    call VecMin ( u4_2r, minerrloc, minimumerror, perr )
    call VecMax ( u4_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u4 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u4 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Gauss test - viewing absolute error in u5 component. \n', perr )
    call VecMin ( u5_2r, minerrloc, minimumerror, perr )
    call VecMax ( u5_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u5 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u5 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Gauss test - viewing absolute error in u6 component. \n', perr )
    call VecMin ( u6_2r, minerrloc, minimumerror, perr )
    call VecMax ( u6_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u6 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u6 component: ', maximumerror, ' at index ', maxerrloc

    call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

  end do DOSLICES

  !> Return velocity and vorticity components to PETSc.
  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

end subroutine tf_RK3GaussRealSpace2
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3GaussCheckAbsoluteError 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3GaussCheckAbsoluteError ( zslice, ierr )

  use g_constants,  only: ZERO, ONE, TWO, FOUR
  use g_parameters, only: NODES_X, NODES_Y, NODES_Z
  use f_arrays,     only: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          y_min_2r, y_max_2r, x_min_2r, x_max_2r

  implicit none
 
  integer,intent(in)    :: zslice                 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: i, j

  real(kind=C_DOUBLE)   :: abserr1, abserr2, abserr3, abserr4, abserr5, abserr6

  write(*,*) 'Checking slice ', zslice

  !> Loop through slices in z direction.
  DOY: do j = y_min_2r, y_max_2r, 1
    DOX: do i = x_min_2r, x_max_2r, 1

!      write(*,*) 'Velocity at (,', i, j, '): ', &
!                 arr_u1_2r(i,j), arr_u2_2r(i,j), arr_u3_2r(i,j), &
!                 ' analytical: ', &
!                 tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, ONE ), &
!                 tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, TWO ), &
!                 tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, FOUR )

      !> Compute absolute error      
      abserr1 = arr_u1_2r(i,j) - tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, ONE )
      abserr2 = arr_u2_2r(i,j) - tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, TWO )
      abserr3 = arr_u3_2r(i,j) - tf_RK3GaussianRealSpace ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, FOUR )

!      call tf_RK3GaussianVorticity ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, &
!                                       ZERO, ZERO, FOUR, abserr4, abserr5, abserr6 )

      call tf_RK3GaussianVorticity ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, &
                                       ONE, TWO, FOUR, abserr4, abserr5, abserr6 )

!      write(*,*) 'Vorticity analytical vs computation ', &
!                 i, j, zslice, abserr4, arr_u4_2r(i,j), abserr5, arr_u5_2r(i,j), abserr6, arr_u6_2r(i,j)

      abserr4 = - abserr4 + arr_u4_2r(i,j)
      abserr5 = - abserr5 + arr_u5_2r(i,j)
      abserr6 = - abserr6 + arr_u6_2r(i,j)

      !> Copy to array
      arr_u1_2r(i,j) = abserr1
      arr_u2_2r(i,j) = abserr2
      arr_u3_2r(i,j) = abserr3
      arr_u4_2r(i,j) = abserr4
      arr_u5_2r(i,j) = abserr5
      arr_u6_2r(i,j) = abserr6

    end do DOX
  end do DOY

end subroutine tf_RK3GaussCheckAbsoluteError
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3GaussError
!> Initialise velocity with a 3D Gaussian!> @param ierr should return 0
subroutine tf_RK3GaussError ( ierr )

  use f_arrays, only: rk3_3w_1, rk3_3w_2, rk3_3w_3

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  real(kind=C_DOUBLE)   :: minimumerror, maximumerror
  integer               :: minerrloc, maxerrloc

  PetscErrorCode        :: perr

  call VecMin ( rk3_3w_1, minerrloc, minimumerror, perr )
  call VecMax ( rk3_3w_1, maxerrloc, maximumerror, perr )

  write(*,*) 'Minimum error in u1 component: ', minimumerror, ' at index ', minerrloc
  write(*,*) 'Maximum error in u1 component: ', maximumerror, ' at index ', maxerrloc

  call VecMin ( rk3_3w_2, minerrloc, minimumerror, perr )
  call VecMax ( rk3_3w_2, maxerrloc, maximumerror, perr )
  write(*,*) 'Minimum error in u2 component: ', minimumerror, ' at index ', minerrloc
  write(*,*) 'Maximum error in u2 component: ', maximumerror, ' at index ', maxerrloc

  call VecMin ( rk3_3w_3, minerrloc, minimumerror, perr )
  call VecMax ( rk3_3w_3, maxerrloc, maximumerror, perr )
  write(*,*) 'Minimum error in u3 component: ', minimumerror, ' at index ', minerrloc
  write(*,*) 'Maximum error in u3 component: ', maximumerror, ' at index ', maxerrloc

end subroutine tf_RK3GaussError
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreen 
!> Initialise velocity with a 3D Gaussian and check if it matches the analytical
!! solution in real space. Transform back to wave space and check if the original
!! field is reproduced. Also check if the imaginary part is zero as required.
!! Perform equivalent tests on the vorticity.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreen ( ierr )

  use f_arrays,     only: u1_3w, u1_1r2w, u1_2r, u1_2w, &
                          u2_3w, u2_1r2w, u2_2r, u2_2w, &
                          u3_3w, u3_1r2w, u3_2r, u3_2w

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  !> Initialise with a Taylor-Gree velocity field in real space
  call tf_RK3TaylorGreenRealSpace ( ierr )
  !> Transform to wave space
  call tf_RK3TaylorGreenTransformBack ( ierr )
  !> Compute vorticity and do a 1D transform in k direction
  call f_RK3PrepareRealSpace ( ierr )
  !> Transform back to real space. Check error in non-linear term, vorticity and velocity on each slice.
  call tf_RK3TaylorGreenRealSpace2 ( ierr )
  !> Transform to wave space
  call tf_RK3TaylorGreenTransformBack ( ierr )
  !> Check error in wavespace array
  call tf_RK3TaylorGreenCheckWavespaceError ( ierr )

end subroutine tf_RK3TaylorGreen
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenRealSpace 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenRealSpace ( ierr )

  use f_arrays, only: da1r2w, da2r, z_min, z_max, &
                      u1_1r2w, u2_1r2w, u3_1r2w, &
                      u4_1r2w, u5_1r2w, u6_1r2w, &
                      arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                      arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                      u1_2r, u2_2r, u3_2r, &
                      u4_2r, u5_2r, u6_2r, &
                      arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                      arr_u4_2r, arr_u5_2r, arr_u6_2r
  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: zslice                 

  PetscErrorCode        :: perr

  !> Get write access to all velocity and vorticity components.
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )


  !> Loop through slices in z direction.
  DOSLICES: do zslice = z_min, z_max, 1
               
    call tf_RK3TaylorGreenInit ( zslice, ierr )

    !> Transform them back.
    call f_RK3TransformAllArraysjiRs2Fs ( zslice, ierr )

  end do DOSLICES

  !> Return velocity and vorticity components to PETSc.
  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

end subroutine tf_RK3TaylorGreenRealSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenInit 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenInit ( zslice, ierr )

  use g_constants,  only: ZERO, ONE, HALF, THIRD, SIXTH, TWO, THREE, SIX
  use g_parameters, only: NODES_X, NODES_Y, NODES_Z, &
                          F_INIT_A, F_INIT_B, F_INIT_C, &
                          F_INIT_I, F_INIT_J, F_INIT_K
  use f_arrays,     only: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          y_min_2r, y_max_2r, x_min_2r, x_max_2r

  implicit none
 
  integer,intent(in)    :: zslice                 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: i, j
  real(kind=C_DOUBLE)   :: u, v, w

  write(*,*) 'Initialising slice ', zslice

  !> Loop through slices in z direction.
  DOY: do j = y_min_2r, y_max_2r, 1
    DOX: do i = x_min_2r, x_max_2r, 1

      call tf_RK3ComputeTaylorGreen  ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, &
                                       F_INIT_A, F_INIT_B, F_INIT_C, &
                                       F_INIT_I, F_INIT_J, F_INIT_K, u, v, w )

      arr_u1_2r(i,j) = u 
      arr_u2_2r(i,j) = v
      arr_u3_2r(i,j) = w

    end do DOX
  end do DOY

end subroutine tf_RK3TaylorGreenInit
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenTransformBack
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenTransformBack ( ierr )

  use f_arrays,     only: da3w, da1r2w, da2w, i_min_3w, i_max_3w, &
                          u1_3w, u2_3w, u3_3w, &
                          rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          u1_1r2w, u2_1r2w, u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          u1_2w, u2_2w, u3_2w, &
                          u4_2w, u5_2w, u6_2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          k_min_3w, k_max_3w, j_min_3w, j_max_3w

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: islice, j, k                 
 
  PetscErrorCode        :: perr

!> Get read access to all non-linear components
  call DMDAVecGetArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> Get read access to all coupling components
  call DMDAVecGetArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Get write access to all RK3 components
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )


  !> Loop through slices in i direction
  DOSLICES: do islice = i_min_3w, i_max_3w, 1

    !> Transform non-linear components.
    call f_RK3TransformAllArrayskRs2Fs ( islice, ierr )


    arr_u1_3w(:,:,:,islice) = arr_u1_2w(:,:,:)
    arr_u2_3w(:,:,:,islice) = arr_u2_2w(:,:,:)
    arr_u3_3w(:,:,:,islice) = arr_u3_2w(:,:,:)
    arr_rk3_3w_1(:,:,:,islice) = arr_u4_2w(:,:,:)
    arr_rk3_3w_2(:,:,:,islice) = arr_u5_2w(:,:,:)
    arr_rk3_3w_3(:,:,:,islice) = arr_u6_2w(:,:,:)
 
!    do j = j_min_3w, j_max_3w, 1 
!      do k = k_min_3w, k_max_3w, 1
!        write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)
!        write(*,*) islice, j, k, arr_u1_3w(:,k,j,islice), arr_u2_3w(:,k,j,islice), arr_u3_3w(:,k,j,islice)
!      end do
!    end do

  end do DOSLICES
  !end loop


!> Return all non-linear components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> Return all coupling components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

!> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Return all RK3 components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )

end subroutine tf_RK3TaylorGreenTransformBack
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenRealSpace2 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenRealSpace2 ( ierr )

  use f_arrays, only: da1r2w, da2r, z_min, z_max, &
                      u1_1r2w, u2_1r2w, u3_1r2w, &
                      u4_1r2w, u5_1r2w, u6_1r2w, &
                      arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                      arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                      u1_2r, u2_2r, u3_2r, &
                      u4_2r, u5_2r, u6_2r, &
                      arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                      arr_u4_2r, arr_u5_2r, arr_u6_2r
  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: zslice                 

  real(kind=C_DOUBLE)   :: minimumerror, maximumerror
  integer               :: minerrloc, maxerrloc

  PetscErrorCode        :: perr

  !> Get write access to all velocity and vorticity components.
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )


  !> Loop through slices in z direction.
  DOSLICES: do zslice = z_min, z_max, 1
               
    call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

    !> Transform velocity and vorticity components to real space.
    call f_RK3TransformToRealSpace ( zslice, 2, 2, ierr )
    !> Transform them back for later use.
    call f_RK3TransformAllArraysjiRs2Fs ( zslice, ierr )

!     write(*,*) arr_u5_2r
!    call VecView ( u4_2r, PETSC_VIEWER_STDOUT_WORLD, ierr )
!    call VecView ( u5_2r, PETSC_VIEWER_STDOUT_WORLD, ierr )
!    call VecView ( u6_2r, PETSC_VIEWER_STDOUT_WORLD, ierr )

    !> Check absolute error in velocity and vorticity.
    call tf_RK3TaylorGreenCheckAbsoluteError ( zslice, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing absolute error in u1 component. \n', perr )
    call VecMin ( u1_2r, minerrloc, minimumerror, perr )
    call VecMax ( u1_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u1 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u1 component: ', maximumerror, ' at index ', maxerrloc

!    call VecView ( u1_2r, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing absolute error in u2 component. \n', perr )
    call VecMin ( u2_2r, minerrloc, minimumerror, perr )
    call VecMax ( u2_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u2 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u2 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing absolute error in u3 component. \n', perr )
    call VecMin ( u3_2r, minerrloc, minimumerror, perr )
    call VecMax ( u3_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u3 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u3 component: ', maximumerror, ' at index ', maxerrloc

    call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing absolute error in u4 component. \n', perr )
    call VecMin ( u4_2r, minerrloc, minimumerror, perr )
    call VecMax ( u4_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u4 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u4 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing absolute error in u5 component. \n', perr )
    call VecMin ( u5_2r, minerrloc, minimumerror, perr )
    call VecMax ( u5_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u5 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u5 component: ', maximumerror, ' at index ', maxerrloc
    call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing absolute error in u6 component. \n', perr )
    call VecMin ( u6_2r, minerrloc, minimumerror, perr )
    call VecMax ( u6_2r, maxerrloc, maximumerror, perr )
    write(*,*) 'Minimum error in u6 component: ', minimumerror, ' at index ', minerrloc
    write(*,*) 'Maximum error in u6 component: ', maximumerror, ' at index ', maxerrloc

    call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

  end do DOSLICES

  !> Return velocity and vorticity components to PETSc.
  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

end subroutine tf_RK3TaylorGreenRealSpace2
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenCheckAbsoluteError 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenCheckAbsoluteError ( zslice, ierr )

  use g_constants,  only: ZERO, ONE, HALF, THIRD, SIXTH, TWO, THREE, SIX
  use g_parameters, only: NODES_X, NODES_Y, NODES_Z, &
                          F_INIT_A, F_INIT_B, F_INIT_C, &
                          F_INIT_I, F_INIT_J, F_INIT_K
  use f_arrays,     only: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          y_min_2r, y_max_2r, x_min_2r, x_max_2r

  implicit none
 
  integer,intent(in)    :: zslice                 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: i, j

  real(kind=C_DOUBLE)   :: u, v, w
  real(kind=C_DOUBLE)   :: abserr1, abserr2, abserr3, abserr4, abserr5, abserr6

  write(*,*) 'Checking slice ', zslice

  !> Loop through slices in z direction.
  DOY: do j = y_min_2r, y_max_2r, 1
    DOX: do i = x_min_2r, x_max_2r, 1

      !> Compute absolute error      
      call tf_RK3ComputeTaylorGreen  ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, & 
                                       F_INIT_A, F_INIT_B, F_INIT_C, &
                                       F_INIT_I, F_INIT_J, F_INIT_K, u, v, w )

!      write(*,*) 'Velocity at (,', i, j, '): ', &
!                 arr_u1_2r(i,j), arr_u2_2r(i,j), arr_u3_2r(i,j), &
!                 ' analytical: ', u, v, w

      abserr1 = arr_u1_2r(i,j) - u 
      abserr2 = arr_u2_2r(i,j) - v
      abserr3 = arr_u3_2r(i,j) - w

      call tf_RK3TaylorGreenVorticity ( i, j, zslice, NODES_X, NODES_Y, NODES_Z, & 
                                       F_INIT_A, F_INIT_B, F_INIT_C, &
                                       F_INIT_I, F_INIT_J, F_INIT_K, u, v, w )

!      write(*,*) 'Vorticity at (,', i, j, '): ', &
!                 arr_u4_2r(i,j), arr_u5_2r(i,j), arr_u6_2r(i,j), &
!                 ' analytical: ', u, v, w

      abserr4 = arr_u4_2r(i,j) - u 
      abserr5 = arr_u5_2r(i,j) - v
      abserr6 = arr_u6_2r(i,j) - w

      !> Copy to array
      arr_u1_2r(i,j) = abserr1
      arr_u2_2r(i,j) = abserr2
      arr_u3_2r(i,j) = abserr3
      arr_u4_2r(i,j) = abserr4
      arr_u5_2r(i,j) = abserr5
      arr_u6_2r(i,j) = abserr6

    end do DOX
  end do DOY

end subroutine tf_RK3TaylorGreenCheckAbsoluteError
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenCheckWavespaceError 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenCheckWavespaceError ( ierr )

  use g_constants,  only: HALF, THIRD, ONE, FOUR
  use g_parameters, only: F_INIT_A, F_INIT_B, F_INIT_C, &
                          F_INIT_I, F_INIT_J, F_INIT_K
  use f_arrays,     only: u1_3w, u2_3w, u3_3w, rk3_3w_1, rk3_3w_2, rk3_3w_3

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  real(kind=C_DOUBLE)   :: error1, error2, error3, error4, error5, error6

  PetscErrorCode        :: perr

  call VecNorm (u1_3w, NORM_1, error1, perr)
  call VecNorm (u2_3w, NORM_1, error2, perr)
  call VecNorm (u3_3w, NORM_1, error3, perr)
  call VecNorm (rk3_3w_1, NORM_1, error4, perr)
  call VecNorm (rk3_3w_2, NORM_1, error5, perr)
  call VecNorm (rk3_3w_3, NORM_1, error6, perr)

!  call VecView ( u1_3w, PETSC_VIEWER_STDOUT_WORLD, ierr )

  !> In each sector of wavespace the trigonometric functions should be represented on on wave mode, in total 8.
  !! Due to conjugate symmetry there are only four. The u3 component is initialised with a factor 1/3.

  error1 = error1 !- HALF
  error2 = error2 !- HALF
  error3 = error3 !- THIRD * HALF
  error4 = error4 !- (real(4*11, kind=C_DOUBLE) / real(24, kind=C_DOUBLE))  ! |(11/3) / 2|
  error5 = error5 !- (FOUR * ONE * THIRD)                                   ! | (8/3) / 2|
  error6 = error6 !- (real(4*3, kind=C_DOUBLE) / real(8, kind=C_DOUBLE))    ! |(-1-2) / 2|


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing error in u1 wavespace array. \n', perr )
  write(*,*) 'Sum of error in u1: ', error1
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing error in u2 wavespace array. \n', perr )
  write(*,*) 'Sum of error in u2: ', error2
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing error in u3 wavespace array. \n', perr )
  write(*,*) 'Sum of error in u3: ', error3
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing error in u4 wavespace array. \n', perr )
  write(*,*) 'Sum of error in u4: ', error4
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing error in u5 wavespace array. \n', perr )
  write(*,*) 'Sum of error in u5: ', error5
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Taylor Green - viewing error in u6 wavespace array. \n', perr )
  write(*,*) 'Sum of error in u6: ', error6

end subroutine tf_RK3TaylorGreenCheckWavespaceError
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3TaylorGreenReverse 
!> Initialise velocity with a 3D Gaussian and check if it matches the analytical
!! solution in real space. Transform back to wave space and check if the original
!! field is reproduced. Also check if the imaginary part is zero as required.
!! Perform equivalent tests on the vorticity.
!> @param ierr should return 0
subroutine tf_RK3TaylorGreenReverse ( ierr )

  use f_arrays,     only: u1_3w, u1_1r2w, u1_2r, u1_2w, &
                          u2_3w, u2_1r2w, u2_2r, u2_2w, &
                          u3_3w, u3_1r2w, u3_2r, u3_2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w
  use f_initialise, only: f_InitialiseInitTaylorGreen

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  !> Initialise with a Taylor-Green velocity field in wave space
  call f_InitialiseInitTaylorGreen ( ierr )
  !> Compute vorticity and do a 1D transform in k direction
  call f_RK3PrepareRealSpace ( ierr )
  !> Transform back to real space. Check error in non-linear term, vorticity and velocity on each slice.
  call tf_RK3TaylorGreenRealSpace2 ( ierr )
  !> Transform to wave space
  call tf_RK3TaylorGreenTransformBack ( ierr )
  !> Check error in wavespace array
  call tf_RK3TaylorGreenCheckWavespaceError ( ierr )

!  call PetscPrintf ( PETSC_COMM_WORLD, 'Viewing u4_1r2w array. \n', ierr )
!  call VecView ( u4_1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )
!  call PetscPrintf ( PETSC_COMM_WORLD, 'Viewing u5_1r2w array. \n', ierr )
!  call VecView ( u5_1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )
!  call PetscPrintf ( PETSC_COMM_WORLD, 'Viewing u6_1r2w array. \n', ierr )
!  call VecView ( u6_1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )

end subroutine tf_RK3TaylorGreenReverse
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3InitTaylorGreenReverse 
!> Initialise velocity with a Taylor Green function in wave space
!> @param ierr should return 0
subroutine tf_RK3InitTaylorGreenReverse ( ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: NODES_KJ, NODES_KK 
  use f_arrays,     only: da3w, u1_3w, u2_3w, u3_3w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          i_min_3w, j_min_3w, k_min_3w, &
                          i_max_3w, j_max_3w, k_max_3w

  implicit none
 
  PetscErrorCode,intent(inout)     :: ierr

  integer                   :: i, j, k

  real(kind=C_DOUBLE)       :: oneeighth, onetwentyfourth

  PetscErrorCode        :: perr

  oneeighth = real(1, kind=C_DOUBLE) / real(8, kind=C_DOUBLE)
  onetwentyfourth = real(1, kind=C_DOUBLE) / real(24, kind=C_DOUBLE)

  !> Initialise all arrays to zero
  call VecSet ( u1_3w, ZERO, perr )
  call VecSet ( u2_3w, ZERO, perr )
  call VecSet ( u3_3w, ZERO, perr )

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

  DOIINIT: do i = i_min_3w, i_max_3w, 1  
    DOJINIT: do j = j_min_3w, j_max_3w, 1  
      DOKINIT: do k = k_min_3w, k_max_3w, 1  

        IFI: if (i == 1) then

          IFJPOS: if (j == 2) then
            IFJPKPOS: if (k == 3) then
              arr_u1_3w(0,k,j,i) = - oneeighth              
              arr_u2_3w(0,k,j,i) = oneeighth              
              arr_u3_3w(0,k,j,i) = - onetwentyfourth              
            end if IFJPKPOS
            IFJPKNEG: if (k == NODES_KK - 3) then
              arr_u1_3w(0,k,j,i) = oneeighth              
              arr_u2_3w(0,k,j,i) = - oneeighth              
              arr_u3_3w(0,k,j,i) = - onetwentyfourth              
            end if IFJPKNEG
          end if IFJPOS

          IFJNEG: if (j == NODES_KJ - 2) then
            IFJNKPOS: if (k == 3) then
              arr_u1_3w(0,k,j,i) = oneeighth              
              arr_u2_3w(0,k,j,i) = oneeighth              
              arr_u3_3w(0,k,j,i) = onetwentyfourth              
            end if IFJNKPOS
            IFJNKNEG: if (k == NODES_KK - 3) then
              arr_u1_3w(0,k,j,i) = - oneeighth              
              arr_u2_3w(0,k,j,i) = - oneeighth              
              arr_u3_3w(0,k,j,i) = onetwentyfourth              
            end if IFJNKNEG
          end if IFJNEG

        end if IFI

      end do DOKINIT
    end do DOJINIT
  end do DOIINIT

!> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

end subroutine tf_RK3InitTaylorGreenReverse
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3CheckNonlinear 
!> Initialise velocity with a 3D Gaussian and check if it matches the analytical
!! solution in real space. Transform back to wave space and check if the original
!! field is reproduced. Also check if the imaginary part is zero as required.
!! Perform equivalent tests on the vorticity.
!> @param ierr should return 0
subroutine tf_RK3CheckNonlinear ( ierr )

  use g_constants,  only: ZERO
  use f_arrays,     only: u1_3w, u1_1r2w, u1_2r, u1_2w, &
                          u2_3w, u2_1r2w, u2_2r, u2_2w, &
                          u3_3w, u3_1r2w, u3_2r, u3_2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  real(kind=C_DOUBLE)   :: dt

  dt = ZERO
  !> Initialise with a Gaussian field in real space
  call tf_RK3InitTaylorGreenReverse ( ierr )
  !> Compute vorticity and do a 1D transform in k direction
  call f_RK3PrepareRealSpace ( ierr )
  !> Transform back to real space. Check error in non-linear term, vorticity and velocity on each slice.
  call f_RK3RealSpace ( 0, 0, dt, ierr )
  !> Transform to wave space
  call tf_RK3TaylorGreenTransformBack ( ierr )
  !> Check error in wavespace array
  call tf_RK3CheckNonlinearError ( ierr )

!  call PetscPrintf ( PETSC_COMM_WORLD, 'Viewing u4_1r2w array. \n', ierr )
!  call VecView ( u4_1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )
!  call PetscPrintf ( PETSC_COMM_WORLD, 'Viewing u5_1r2w array. \n', ierr )
!  call VecView ( u5_1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )
!  call PetscPrintf ( PETSC_COMM_WORLD, 'Viewing u6_1r2w array. \n', ierr )
!  call VecView ( u6_1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )

end subroutine tf_RK3CheckNonlinear
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3CheckNonlinearError 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3CheckNonlinearError ( ierr )

  use g_constants,  only: THREE, FOUR, EIGHT, NINE, ELEVEN
  use f_arrays,     only: rk3_3w_1, rk3_3w_2, rk3_3w_3

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  real(kind=C_DOUBLE)   :: error4, error5, error6

  PetscErrorCode        :: perr

  call VecNorm (rk3_3w_1, NORM_1, error4, perr)
  call VecNorm (rk3_3w_2, NORM_1, error5, perr)
  call VecNorm (rk3_3w_3, NORM_1, error6, perr)

!  call VecView ( u1_3w, PETSC_VIEWER_STDOUT_WORLD, ierr )

  !> Due to conjugate symmetry not all wavemodes are represented. The derivation of the analytical
  !! solution will be documented in the project report.

  error4 = error4 - THREE / FOUR
  error5 = error5 - NINE / EIGHT
  error6 = error6 - ELEVEN / EIGHT


  call PetscPrintf ( PETSC_COMM_WORLD, '    | Viewing error in u4 nonlinear component. \n', perr )
  write(*,*) 'Sum of error in u4: ', error4
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Viewing error in u5 nonlinear component. \n', perr )
  write(*,*) 'Sum of error in u5: ', error5
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Viewing error in u6 nonlinear component. \n', perr )
  write(*,*) 'Sum of error in u6: ', error6

end subroutine tf_RK3CheckNonlinearError
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3RealSpace 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param ierr should return 0
subroutine tf_RK3RealSpace ( ierr )

  use f_arrays, only: da1r2w, da2r, z_min, z_max, &
                      u1_1r2w, u2_1r2w, u3_1r2w, &
                      u4_1r2w, u5_1r2w, u6_1r2w, &
                      arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                      arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                      u1_2r, u2_2r, u3_2r, &
                      u4_2r, u5_2r, u6_2r, &
                      arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                      arr_u4_2r, arr_u5_2r, arr_u6_2r
  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: zslice                 

  PetscErrorCode        :: perr

  !> Get write access to all velocity and vorticity components.
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )


  !> Loop through slices in z direction.
  DOSLICES: do zslice = z_min, z_max, 1
               
    !> Transform velocity and vorticity components to real space.
    call f_RK3TransformToRealSpace ( zslice, 2, 2, ierr )

    !> Transform them back.
    call f_RK3TransformAllArraysjiRs2Fs ( zslice, ierr )

  end do DOSLICES

  !> Return velocity and vorticity components to PETSc.
  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

end subroutine tf_RK3RealSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine tf_RK3SolveNavierStokes
!> @param ierr should return 0
subroutine tf_RK3SolveNavierStokes ( ierr )

  use f_arrays,     only: da3w, da1r2w, da2w, i_min_3w, i_max_3w, &
                          u1_3w, u2_3w, u3_3w, &
                          rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          u1_1r2w, u2_1r2w, u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          u1_2w, u2_2w, u3_2w, &
                          u4_2w, u5_2w, u6_2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w
  implicit none
 
  PetscErrorCode,intent(inout) :: ierr

  integer               :: islice                 

  PetscErrorCode        :: perr

 
!> Get read access to all non-linear components
  call DMDAVecGetArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> Get read access to all coupling components
  call DMDAVecGetArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Get write access to all RK3 components
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )


  !> Loop through slices in i direction
  DOSLICES: do islice = i_min_3w, i_max_3w, 1

    !> Transform non-linear components.
    call f_RK3TransformAllArrayskRs2Fs ( islice, ierr )
 
  end do DOSLICES
  !end loop


!> Return all non-linear components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> Return all coupling components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

!> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Return all RK3 components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )

end subroutine tf_RK3SolveNavierStokes
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
!pure function tf_RK3ComputeGaussian ( a, k ) result ( gaussian )
!  real(kind=C_DOUBLE),intent(in)   :: a, k 
!  real(kind=C_DOUBLE)              :: pi, gaussian

!  pi = asin(1.0)
!  gaussian = exp(- ((pi**2 * k**2)/(4.0*a))) 
!  gaussian = ((pi/a)**0.5) * exp(- ((k**2)/(4.0*a))) 
!  gaussian = (sqrt(pi)/sigma**2)**(1.0/3.0) * exp(-((pi**2) * ((k**2)/(4.0*sigma**2))))  
!  gaussian = (sqrt(pi)/sigma) * exp(-((k**2)/(4.0*sigma**2)))  

!end function tf_RK3ComputeGaussian
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
pure function tf_RK3GaussianRealSpace ( i, j, k, nx, ny, nz, a ) result ( gaussian )

  use g_constants,       only : LONGTWOPI

  integer, intent(in)              :: i, j, k, nx, ny, nz
  real(kind=C_DOUBLE), intent(in)  :: a
  real(kind=C_DOUBLE)              :: gaussian

  real(kind=C_LONG_DOUBLE)              :: x, y, z, argument, exactresult
!  real(kind=C_LONG_DOUBLE)              :: twopi

!  twopi = 4.0 * asin(1.0)

  x = LONGTWOPI * real((i - (nx/2)), kind=C_LONG_DOUBLE) / real (nx, kind=C_LONG_DOUBLE)
  y = LONGTWOPI * real((j - (ny/2)), kind=C_LONG_DOUBLE) / real (ny, kind=C_LONG_DOUBLE)
  z = LONGTWOPI * real((k - (nz/2)), kind=C_LONG_DOUBLE) / real (nz, kind=C_LONG_DOUBLE)

  argument = real ( -a, kind=C_LONG_DOUBLE) * ( x**2 + y**2 + z**2 )
  exactresult = exp ( argument )
!  exp( -a * ( x**2 + y**2 + z**2 ) ) 
  gaussian = real ( exactresult, kind = C_DOUBLE )

end function tf_RK3GaussianRealSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
pure subroutine tf_RK3GaussianVorticity ( i, j, k, nx, ny, nz, ax, ay, az, vx, vy, vz )

  use g_constants,       only : LONGTWOPI

  integer, intent(in)              :: i, j, k, nx, ny, nz
  real(kind=C_DOUBLE), intent(in)  :: ax, ay, az
  real(kind=C_DOUBLE), intent(out) :: vx, vy, vz

  real(kind=C_LONG_DOUBLE)              :: dudy, dudz, dvdx, dvdz, dwdx, dwdy

  real(kind=C_LONG_DOUBLE)              :: x, y, z, amplitude, argument, exactresult
!  real(kind=C_LONG_DOUBLE)              :: twopi

!  twopi = 4.0 * asin(1.0)

  x = LONGTWOPI * real((i - (nx/2)), kind=C_LONG_DOUBLE) / real (nx, kind=C_LONG_DOUBLE)
  y = LONGTWOPI * real((j - (ny/2)), kind=C_LONG_DOUBLE) / real (ny, kind=C_LONG_DOUBLE)
  z = LONGTWOPI * real((k - (nz/2)), kind=C_LONG_DOUBLE) / real (nz, kind=C_LONG_DOUBLE)

  amplitude = 2.0
  amplitude = - amplitude * real(ax, kind=C_LONG_DOUBLE) * y
  argument = real ( -ax, kind=C_LONG_DOUBLE) * ( x**2 + y**2 + z**2 )
  dudy = amplitude * exp( argument ) 
!  dudy = -2.0 * ax * y * exp( -ax * ( x**2 + y**2 + z**2 ) ) 
  amplitude = 2.0
  amplitude = - amplitude * real(ax, kind=C_LONG_DOUBLE) * z
  dudz = amplitude * exp( argument )  
!  dudz = -2.0 * ax * z * exp( -ax * ( x**2 + y**2 + z**2 ) ) 
  amplitude = 2.0
  amplitude = - amplitude * real(ay, kind=C_LONG_DOUBLE) * x
  argument = real ( -ay, kind=C_LONG_DOUBLE) * ( x**2 + y**2 + z**2 )
!  dvdx = amplitude * exp( argument )  
  dvdx = -2.0 * ay * x * exp( -ay * ( x**2 + y**2 + z**2 ) ) 
  amplitude = 2.0
  amplitude = - amplitude * real(ay, kind=C_LONG_DOUBLE) * z
  dvdz = amplitude * exp( argument )  
!  dvdz = -2.0 * ay * z * exp( -ay * ( x**2 + y**2 + z**2 ) ) 
  amplitude = 2.0
  amplitude = - amplitude * real(az, kind=C_LONG_DOUBLE) * x
  argument = real ( -az, kind=C_LONG_DOUBLE) * ( x**2 + y**2 + z**2 )
  dwdx = amplitude * exp( argument )  
!  dwdx = -2.0 * az * x * exp( -az * ( x**2 + y**2 + z**2 ) ) 
  amplitude = 2.0
  amplitude = - amplitude * real(az, kind=C_LONG_DOUBLE) * y
  dwdy = amplitude * exp( argument )  
!  dwdy = -2.0 * az * y * exp( -az * ( x**2 + y**2 + z**2 ) ) 

  exactresult = dwdy - dvdz
  vx = real(exactresult, kind=C_DOUBLE)

  exactresult = dudz - dwdx
  vy = real(exactresult, kind=C_DOUBLE)

  exactresult = dvdx - dudy
  vz = real(exactresult, kind=C_DOUBLE)

!  vx = dwdy - dvdz
!  vy = dudz - dwdx
!  vz = dvdx - dudy

end subroutine tf_RK3GaussianVorticity
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
pure subroutine tf_RK3ComputeTaylorGreen ( i, j, k, nx, ny, nz, a1, b1, c1, a2, b2, c2, u, v, w )

  use g_constants,       only : LONGTWOPI

  integer, intent(in)              :: i, j, k, nx, ny, nz
  real(kind=C_DOUBLE), intent(in)  :: a1, b1, c1
  integer, intent(in)              :: a2, b2, c2
  real(kind=C_DOUBLE), intent(out) :: u, v, w

  real(kind=C_LONG_DOUBLE)              :: x, y, z, amplitude, &
                                           argument1, argument2, argument3, exactresult
!  real(kind=C_LONG_DOUBLE)              :: twopi

!  twopi = 4.0 * asin(1.0)

  x = LONGTWOPI * real((i - (nx/2)), kind=C_LONG_DOUBLE) / real (nx, kind=C_LONG_DOUBLE)
  y = LONGTWOPI * real((j - (ny/2)), kind=C_LONG_DOUBLE) / real (ny, kind=C_LONG_DOUBLE)
  z = LONGTWOPI * real((k - (nz/2)), kind=C_LONG_DOUBLE) / real (nz, kind=C_LONG_DOUBLE)

  !> Same arguments for all trigonometric functions
  argument1 = real(a2, kind=C_LONG_DOUBLE) * x
  argument2 = real(b2, kind=C_LONG_DOUBLE) * y
  argument3 = real(c2, kind=C_LONG_DOUBLE) * z

  amplitude = real(a1, kind=C_LONG_DOUBLE)
  exactresult = amplitude * cos ( argument1 ) * sin ( argument2 ) * sin ( argument3 )
  u = real(exactresult, kind=C_DOUBLE)

  amplitude = real(b1, kind=C_LONG_DOUBLE)
  exactresult = amplitude * sin ( argument1 ) * cos ( argument2 ) * sin ( argument3 )
  v = real(exactresult, kind=C_DOUBLE)

  amplitude = real(c1, kind=C_LONG_DOUBLE)
  exactresult = amplitude * sin ( argument1 ) * sin ( argument2 ) * cos ( argument3 )
  w = real(exactresult, kind=C_DOUBLE)

end subroutine tf_RK3ComputeTaylorGreen
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine tf_RK3TaylorGreenVorticity ( i, j, k, nx, ny, nz, a1, b1, c1, a2, b2, c2, vx, vy, vz )

  use g_constants,       only : LONGTWOPI

  integer, intent(in)              :: i, j, k, nx, ny, nz
  real(kind=C_DOUBLE), intent(in)  :: a1, b1, c1 
  integer, intent(in)              :: a2, b2, c2
  real(kind=C_DOUBLE), intent(out) :: vx, vy, vz

  real(kind=C_LONG_DOUBLE)              :: x, y, z, amplitude, exactresult, &
                                           argument1, argument2, argument3

!  real(kind=C_LONG_DOUBLE)              :: twopi

  real(kind=C_LONG_DOUBLE)              :: dudy, dudz, dvdx, dvdz, dwdx, dwdy

  x = LONGTWOPI * real((i - (nx/2)), kind=C_LONG_DOUBLE) / real (nx, kind=C_LONG_DOUBLE)
  y = LONGTWOPI * real((j - (ny/2)), kind=C_LONG_DOUBLE) / real (ny, kind=C_LONG_DOUBLE)
  z = LONGTWOPI * real((k - (nz/2)), kind=C_LONG_DOUBLE) / real (nz, kind=C_LONG_DOUBLE)

!  write(*,*) ' i, j, k, nx, ny, nz, x, y, z: ', &
!             i, j, k, nx, ny, nz, x, y, z

  !> Same arguments for all trigonometric functions
  argument1 = real(a2, kind=C_LONG_DOUBLE) * x
  argument2 = real(b2, kind=C_LONG_DOUBLE) * y
  argument3 = real(c2, kind=C_LONG_DOUBLE) * z

  exactresult = amplitude * cos ( argument1 ) * sin ( argument2 ) * sin ( argument3 )

  amplitude = real(a1, kind=C_LONG_DOUBLE) * real(b2, kind=C_LONG_DOUBLE)
  dudy = amplitude * cos ( argument1 ) * cos ( argument2 ) * sin ( argument3 )
  amplitude = real(a1, kind=C_LONG_DOUBLE) * real(c2, kind=C_LONG_DOUBLE)
  dudz = amplitude * cos ( argument1 ) * sin ( argument2 ) * cos ( argument3 )
  amplitude = real(b1, kind=C_LONG_DOUBLE) * real(a2, kind=C_LONG_DOUBLE)
  dvdx = amplitude * cos ( argument1 ) * cos ( argument2 ) * sin ( argument3 )
  amplitude = real(b1, kind=C_LONG_DOUBLE) * real(c2, kind=C_LONG_DOUBLE)
  dvdz = amplitude * sin ( argument1 ) * cos ( argument2 ) * cos ( argument3 )
  amplitude = real(c1, kind=C_LONG_DOUBLE) * real(a2, kind=C_LONG_DOUBLE)
  dwdx = amplitude * cos ( argument1 ) * sin ( argument2 ) * cos ( argument3 )
  amplitude = real(c1, kind=C_LONG_DOUBLE) * real(b2, kind=C_LONG_DOUBLE)
  dwdy = amplitude * sin ( argument1 ) * cos ( argument2 ) * cos ( argument3 )

  exactresult = dwdy - dvdz
  vx = real(exactresult, kind=C_DOUBLE)

!  write(*,*) ' exactresult, dwdy, dvdz, vx: ', &
!             exactresult, dwdy, dvdz, vx

  exactresult = dudz - dwdx
  vy = real(exactresult, kind=C_DOUBLE)

!  write(*,*) ' exactresult, dudz, dwdx, vy: ', &
!             exactresult, dudz, dwdx, vy

  exactresult = dvdx - dudy
  vz = real(exactresult, kind=C_DOUBLE)

!  write(*,*) ' exactresult, dvdx, dudy, vz: ', &
!             exactresult, dvdx, dudy, vz

!  write(*,*) ' vx, vy, vz, dwdy, dvdz, dudz, dwdx, dvdx, dudy: ', &
!             vx, vy, vz, dwdy, dvdz, dudz, dwdx, dvdx, dudy, &
!             ' a1, b1, c1, a2, b2, c2: ', &
!             a1, b1, c1, a2, b2, c2, &
!             ' x, y, z: ', &
!             x, y, z, &
!             ' argument1, argument2, argument3: ', &
!             argument1, argument2, argument3, &
!             ' cos / sin: ', &
!             cos(argument1), cos(argument2), cos(argument3), &
!             sin(argument1), sin(argument2), sin(argument3), &
!             ' amplitudes: ', &
!             real(a1, kind=C_LONG_DOUBLE) * real(b2, kind=C_LONG_DOUBLE), &
!             real(a1, kind=C_LONG_DOUBLE) * real(c2, kind=C_LONG_DOUBLE), &
!             real(b1, kind=C_LONG_DOUBLE) * real(a2, kind=C_LONG_DOUBLE), &
!             real(b1, kind=C_LONG_DOUBLE) * real(c2, kind=C_LONG_DOUBLE), &
!             real(c1, kind=C_LONG_DOUBLE) * real(a2, kind=C_LONG_DOUBLE), &
!             real(c1, kind=C_LONG_DOUBLE) * real(b2, kind=C_LONG_DOUBLE)


end subroutine tf_RK3TaylorGreenVorticity
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
end module tf_rk3
