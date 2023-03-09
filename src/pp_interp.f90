module p_interp

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

!----------------------------------------------------------------------

public :: p_InterpAllocate
public :: p_InterpDeallocate
!public :: p_InterpFindMolecule1D
!public :: p_InterpInterpolate2D
public :: p_InterpSetToZero 
public :: p_InterpGetLocalVelocity    
public :: p_InterpReturnLocalVelocity    
public :: p_InterpInterpolate3D 

!----------------------------------------------------------------------
integer,parameter,public                          :: FIELD_VEL = 0
integer,parameter,public                          :: FIELD_EL  = 2
integer,parameter,public                          :: FIELD_EP  = 3

!> Array containing all values in the 3D interpolation molecule
real(kind=C_DOUBLE),allocatable,dimension(:,:,:,:) :: interpolationvalues
!> Array containing the x coordinates of the 3D interpolation molecule
real(kind=C_DOUBLE),allocatable,dimension(:)       :: interpolationx
!> Array containing the y coordinates of the 3D interpolation molecule
real(kind=C_DOUBLE),allocatable,dimension(:)       :: interpolationy
!> Array containing the z coordinates of the 3D interpolation molecule
real(kind=C_DOUBLE),allocatable,dimension(:)       :: interpolationz
!> Array containing intermediate velocity value in the 2D molecule
real(kind=C_DOUBLE),dimension(1:3)                 :: interpolationuf

contains
!---------------------------------------------------------------------------------
subroutine p_InterpAllocate ( ierr )

  use g_parameters,   only : P_INTERP_ORDER

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  allocate ( interpolationvalues(1:3,1:P_INTERP_ORDER,&
             1:P_INTERP_ORDER,1:P_INTERP_ORDER), stat=ierr )
!  CHKERRQ ( ierr )
  allocate ( interpolationx(1:P_INTERP_ORDER), stat=ierr )
!  CHKERRQ ( ierr )
  allocate ( interpolationy(1:P_INTERP_ORDER), stat=ierr )
!  CHKERRQ ( ierr )
  allocate ( interpolationz(1:P_INTERP_ORDER), stat=ierr )
!  CHKERRQ ( ierr )

end subroutine p_InterpAllocate 
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_InterpDeallocate ( ierr )

  use g_parameters,   only : P_INTERP_ORDER

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  deallocate ( interpolationvalues, stat=ierr )
!  CHKERRQ ( ierr )
  deallocate ( interpolationx, stat=ierr )
!  CHKERRQ ( ierr )
  deallocate ( interpolationy, stat=ierr )
!  CHKERRQ ( ierr )
  deallocate ( interpolationz, stat=ierr )
!  CHKERRQ ( ierr )

end subroutine p_InterpDeallocate 
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_InterpSetToZero ( npin, ptr_values, ierr )    

  use g_constants, only : ZERO

  integer,intent(in)     :: npin
  type(C_PTR),intent(in) :: ptr_values
  PetscErrorCode,intent(inout)  :: ierr

  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particlevalues

  !> Associate interpolation values pointer array with target array.
  call c_f_pointer ( ptr_values, particlevalues, [3, npin])

  particlevalues = ZERO 

end subroutine p_InterpSetToZero
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_InterpGetLocalVelocity 
subroutine p_InterpGetLocalVelocity ( ierr )    

  use f_arrays,     only : da3r, arr_u1_3r, u1_3r, u1_3r_local, &
                           arr_u2_3r, u2_3r, u2_3r_local, &
                           arr_u3_3r, u3_3r, u3_3r_local, z_min, z_max, &
                           x_min_2r, x_max_2r, y_min_2r, y_max_2r, &
                           x_width_2r, y_width_2r, z_width_1r2w, &
                           x_min_ghost, y_min_ghost, z_min_ghost, &
                           x_max_ghost, y_max_ghost, z_max_ghost
  use g_constants,  only : ONE
  use g_parameters, only : MYRANK, NHALO

  implicit none

  PetscErrorCode,intent(inout)  :: ierr
  PetscErrorCode         :: perr

  PetscScalar            :: localsum, globalsum

  integer                :: i, j, k

  !> Get local vector with halos
  call  DMGetLocalVector( da3r, u1_3r_local, perr )
!  call  DMCreateLocalVector( da3r, u1_3r_local, perr )
!  CHKERRQ(perr)
  call  DMGetLocalVector( da3r, u2_3r_local, perr )
!  call  DMCreateLocalVector( da3r, u2_3r_local, perr )
!  CHKERRQ(perr)
  call  DMGetLocalVector( da3r, u3_3r_local, perr )
!   call  DMCreateLocalVector( da3r, u3_3r_local, perr )
!  CHKERRQ(perr)

  !> Update local vectors from global vector
  call DMGlobalToLocalBegin ( da3r, u1_3r, INSERT_VALUES, u1_3r_local, perr )
!  CHKERRQ(perr)
  call DMGlobalToLocalEnd   ( da3r, u1_3r, INSERT_VALUES, u1_3r_local, perr )

!  call VecSum ( u1_3r, globalsum, ierr )
!  call VecSum ( u1_3r_local, localsum, ierr )

!  write(*,*) MYRANK, 'global sum: ', globalsum, ' local sum: ', localsum

!  CHKERRQ(perr)
!  call VecView ( u1_3r_local, PETSC_VIEWER_STDOUT_SELF )
  !> Return global velocity vector
  call  DMRestoreGlobalVector( da3r, u1_3r, perr )
!  CHKERRQ(perr)

  !> Update local vectors from global vector
  call DMGlobalToLocalBegin ( da3r, u2_3r, INSERT_VALUES, u2_3r_local, perr )
!  CHKERRQ(perr)
  call DMGlobalToLocalEnd   ( da3r, u2_3r, INSERT_VALUES, u2_3r_local, perr )
!  CHKERRQ(perr)
!  call VecSum ( u2_3r, globalsum, ierr )
!  call VecSum ( u2_3r_local, localsum, ierr )

!  write(*,*) MYRANK, 'global sum: ', globalsum, ' local sum: ', localsum

!  call VecView ( u2_3r_local, PETSC_VIEWER_STDOUT_SELF )
  !> Return global velocity vector
  call  DMRestoreGlobalVector( da3r, u2_3r, perr )
!  CHKERRQ(perr)

  !> Update local vectors from global vector
  call DMGlobalToLocalBegin ( da3r, u3_3r, INSERT_VALUES, u3_3r_local, perr )
!  CHKERRQ(perr)
  call DMGlobalToLocalEnd   ( da3r, u3_3r, INSERT_VALUES, u3_3r_local, perr )
!  CHKERRQ(perr)
  call VecSum ( u3_3r, globalsum, ierr )
  call VecSum ( u3_3r_local, localsum, ierr )

!  write(*,*) MYRANK, 'global sum: ', globalsum, ' local sum: ', localsum

!  if ( MYRANK == 4 ) then
!    write(*,*) ' Rank: ', MYRANK, ' global boundaries: ', &
!               x_min_2r, x_max_2r, z_min, z_max, y_min_2r, y_max_2r, &
!               ' local boundaries: ', &
!               x_min_ghost, x_max_ghost, z_min_ghost, z_max_ghost, &
!               y_min_ghost, y_max_ghost
!    call VecView ( u3_3r_local, PETSC_VIEWER_STDOUT_SELF, ierr )
!  end if

!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  !> Return global velocity vector
  call  DMRestoreGlobalVector( da3r, u3_3r, perr )
!  CHKERRQ(perr)

!  call VecSet ( u1_3r_local, ONE, ierr )
!  CHKERRQ ( ierr )
!  call VecSet ( u2_3r_local, ONE, ierr )
!  CHKERRQ ( ierr )
!  call VecSet ( u3_3r_local, ONE, ierr )
!  CHKERRQ ( ierr )
  !> Get local velocity arrays with halos.
  call DMDAVecGetArrayReadF90 ( da3r, u1_3r_local, arr_u1_3r, perr )
!  CHKERRQ(perr)
  call DMDAVecGetArrayReadF90 ( da3r, u2_3r_local, arr_u2_3r, perr )
!  CHKERRQ(perr)
  call DMDAVecGetArrayReadF90 ( da3r, u3_3r_local, arr_u3_3r, perr )
!  CHKERRQ(perr)
!  call VecGetArrayReadF90 ( u1_3r_local, arr_u1_3r, perr )
!  CHKERRQ(perr)
!  call VecGetArrayReadF90 ( u2_3r_local, arr_u2_3r, perr )
!  CHKERRQ(perr)
!  call VecGetArrayReadF90 ( u3_3r_local, arr_u3_3r, perr )
!  CHKERRQ(perr)

!  write(*,*) ' 3D u1 local array: ', shape(arr_u1_3r_local), &
!             ' 3D u2 local array: ', shape(arr_u2_3r_local), &
!             ' 3D u3 local array: ', shape(arr_u3_3r_local)

!  write(*,*) ' Rank: ', MYRANK, ' global boundaries: ', &
!             x_min_2r, x_max_2r, z_min, z_max, y_min_2r, y_max_2r, &
!             ' local boundaries: ', &
!             x_min_ghost, x_max_ghost, z_min_ghost, z_max_ghost, &
!             y_min_ghost, y_max_ghost
!  write(*,*) ' 3D u1 local array: ', arr_u1_3r(:,:,:), &
!             ' 3D u2 local array: ', arr_u2_3r(:,:,:), &
!             ' 3D u3 local array: ', arr_u3_3r(:,:,:)

!  write(*,*) MYRANK, 'Test 1 to 16', MYRANK, arr_u3_3r(1,1,1), &
!             arr_u3_3r(2,2,2), arr_u3_3r(3,3,3), &
!             arr_u3_3r(4,4,4), arr_u3_3r(5,5,5), &
!             arr_u3_3r(6,6,6), arr_u3_3r(7,7,7), &
!             arr_u3_3r(8,8,8), arr_u3_3r(9,9,9), &
!             arr_u3_3r(10,10,10), arr_u3_3r(11,11,11), &
!             arr_u3_3r(12,12,12), arr_u3_3r(13,13,13), &
!             arr_u3_3r(14,14,14), arr_u3_3r(15,15,15), &
!             arr_u3_3r(16,16,16)

!  if ( MYRANK == 4 ) then
!    write(*,*) ' u3: ', arr_u3_3r_local(:,:,:)
!  end if

!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

!  do j = y_min_ghost, y_max_ghost, 1
!    do k = z_min_ghost, z_max_ghost, 1
!      do i = x_min_ghost, x_max_ghost, 1 
!        write(*,*) MYRANK, i, k, j, int(arr_u3_3r(i,k,j)), &
!                   int(arr_u1_3r(i,k,j))!, &
!                   int(arr_u2_3r(x_min_2r+i,z_min+k,y_min_2r+j))
!      end do
!    end do
!  end do

  ierr = perr

end subroutine p_InterpGetLocalVelocity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! p_InterpReturnLocalVelocity 
subroutine p_InterpReturnLocalVelocity ( ierr )    
  use g_constants,  only : ZERO
  use f_arrays,     only : da3r, arr_u1_3r, u1_3r, u1_3r_local, &
                           arr_u2_3r, u2_3r, u2_3r_local, &
                           arr_u3_3r, u3_3r, u3_3r_local

  implicit none

  PetscErrorCode,intent(inout)  :: ierr

  PetscErrorCode         :: perr

  !> Return local velocity arrays.
  call DMDAVecRestoreArrayReadF90 ( da3r, u1_3r_local, arr_u1_3r, perr )
  call DMDAVecRestoreArrayReadF90 ( da3r, u2_3r_local, arr_u2_3r, perr )
  call DMDAVecRestoreArrayReadF90 ( da3r, u3_3r_local, arr_u3_3r, perr )
!  call VecRestoreArrayReadF90 ( u1_3r_local, arr_u1_3r, perr )
!  call VecRestoreArrayReadF90 ( u2_3r_local, arr_u2_3r, perr )
!  call VecRestoreArrayReadF90 ( u3_3r_local, arr_u3_3r, perr )

  !> Return local velocity vectors
  call  DMRestoreLocalVector( da3r, u1_3r_local, perr )
  call  DMRestoreLocalVector( da3r, u2_3r_local, perr )
  call  DMRestoreLocalVector( da3r, u3_3r_local, perr )
  call  VecDestroy( u1_3r_local, perr )
  call  VecDestroy( u2_3r_local, perr )
  call  VecDestroy( u3_3r_local, perr )

  !> Get global vectors again and initialise them to zero
  call DMGetGlobalVector ( da3r, u1_3r, perr )
  call DMGetGlobalVector ( da3r, u2_3r, perr )
  call DMGetGlobalVector ( da3r, u3_3r, perr )

  call VecSet ( u1_3r, ZERO, perr )
  call VecSet ( u2_3r, ZERO, perr )
  call VecSet ( u3_3r, ZERO, perr )

end subroutine p_InterpReturnLocalVelocity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! p_InterpInterpolate3D 
subroutine p_InterpInterpolate3D ( xp, uf, ierr )    
  use g_constants,  only : ZERO
  use g_parameters, only : P_INTERP_ORDER, NHALO 
  use f_arrays,     only : arr_u1_3r, arr_u2_3r, arr_u3_3r, z_min, z_max, &
                           x_min_2r, y_min_2r, x_max_2r, y_max_2r

  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(in)  :: xp
  real(kind=C_DOUBLE),dimension(1:3),intent(out) :: uf
  PetscErrorCode,intent(inout)                          :: ierr

  integer                                    :: i0
  integer                                    :: j0
  integer                                    :: k0
  integer                                    :: kp


  !> Get particle coordinates
  !> Find interpolation molecule
  call p_InterpFindMolecule ( xp, i0, j0, k0, ierr )
  !> Copy fluid velocity to interpolation molecule 
  interpolationvalues(1,:,:,:) = arr_u1_3r(i0:(i0+P_INTERP_ORDER-1),&
                                           k0:(k0+P_INTERP_ORDER-1),&
                                           j0:(j0+P_INTERP_ORDER-1))
  interpolationvalues(2,:,:,:) = arr_u2_3r(i0:(i0+P_INTERP_ORDER-1),&
                                           k0:(k0+P_INTERP_ORDER-1),&
                                           j0:(j0+P_INTERP_ORDER-1))
  interpolationvalues(3,:,:,:) = arr_u3_3r(i0:(i0+P_INTERP_ORDER-1),&
                                           k0:(k0+P_INTERP_ORDER-1),&
                                           j0:(j0+P_INTERP_ORDER-1))

!  write(*,*) 'Interpolation molecule: ', i0, i0+P_INTERP_ORDER-1, &
!             j0, j0+P_INTERP_ORDER-1, k0, k0+P_INTERP_ORDER-1, &
!             ' position: ', xp, &
!             ' x: ', x_min_2r, x_max_2r, ' y: ', y_min_2r, y_max_2r, &
!             ' z: ', z_min, z_max, ' values: ', interpolationvalues, &
!             ' u1 at corners: ', arr_u1_3r(x_min_2r,z_min,y_min_2r), &
!             arr_u1_3r(x_max_2r,z_max,y_max_2r), ' NHALO: ', NHALO, &
!             ' P_INTERP_ORDER: ', P_INTERP_ORDER, & 
!             ' u1 at halo: ', arr_u1_3r(x_min_2r-NHALO,z_min-NHALO,y_min_2r-NHALO), &
!             arr_u1_3r(x_max_2r+NHALO,z_max+NHALO,y_max_2r+NHALO)
!  write(*,*) ' imin: ', i0, ' u1 at lower i border: ', arr_u1_3r(i0,:,:), & 
!             ' jmin: ', j0, ' u1 at lower j border: ', arr_u1_3r(:,:,j0), & 
!             ' kmin: ', k0, 'u1 at lower k border: ', arr_u1_3r(:,k0,:), & 
!             ' imax: ', i0+P_INTERP_ORDER-1, &
!             'u1 at upper i border: ', arr_u1_3r(i0+P_INTERP_ORDER-1,:,:), & 
!             ' jmax: ', j0+P_INTERP_ORDER-1, &
!             'u1 at upper j border: ', arr_u1_3r(:,:,j0+P_INTERP_ORDER-1), & 
!             ' kmax: ', k0+P_INTERP_ORDER-1, &
!             'u1 at upper k border: ', arr_u1_3r(:,k0+P_INTERP_ORDER-1,:)
!  write(*,*) ' xmin: ', x_min_2r-NHALO, ' u1 at lower x border: ', arr_u1_3r(x_min_2r-NHALO,:,:), & 
!             ' ymin: ', y_min_2r-NHALO, ' u1 at lower y border: ', arr_u1_3r(:,:,y_min_2r-NHALO), & 
!             ' zmin: ', z_min-NHALO, 'u1 at lower z border: ', arr_u1_3r(:,z_min-NHALO,:), & 
!             ' xmax: ', x_max_2r+NHALO, &
!             'u1 at upper i border: ', arr_u1_3r(x_max_2r+NHALO,:,:), & 
!             ' ymax: ', y_max_2r+NHALO, &
!             'u1 at upper y border: ', arr_u1_3r(:,:,y_max_2r+NHALO), & 
!             ' zmax: ', z_max+NHALO, &
!             'u1 at upper z border: ', arr_u1_3r(:,z_max+NHALO,:)
!  write(*,*) ' xmin: ', x_min_2r, ' u1 at lower x border: ', arr_u1_3r(x_min_2r,:,:), & 
!             ' ymin: ', y_min_2r, ' u1 at lower y border: ', arr_u1_3r(:,:,y_min_2r), & 
!             ' zmin: ', z_min, 'u1 at lower z border: ', arr_u1_3r(:,z_min,:), & 
!             ' xmax: ', x_max_2r, &
!             'u1 at upper i border: ', arr_u1_3r(x_max_2r,:,:), & 
!             ' ymax: ', y_max_2r, &
!             'u1 at upper y border: ', arr_u1_3r(:,:,y_max_2r), & 
!             ' zmax: ', z_max, &
!             'u1 at upper z border: ', arr_u1_3r(:,z_max,:)
  if ( i0 < x_min_2r-NHALO ) write(*,*) 'i0 ', i0,  x_min_2r-NHALO
  if ( j0 < y_min_2r-NHALO ) write(*,*) 'j0 ', j0,  y_min_2r-NHALO
  if ( k0 < z_min-NHALO ) write(*,*) 'k0 ', k0,  z_min-NHALO
  if ( i0+P_INTERP_ORDER-1 > x_max_2r+NHALO ) write(*,*) 'in ', &
       i0+P_INTERP_ORDER-1,  x_max_2r+NHALO
  if ( j0+P_INTERP_ORDER-1 > y_max_2r+NHALO ) write(*,*) 'jn ', &
       j0+P_INTERP_ORDER-1,  y_max_2r+NHALO
  if ( k0+P_INTERP_ORDER-1 > z_max+NHALO ) write(*,*) 'kn ', &
       k0+P_INTERP_ORDER-1,  z_max+NHALO

  !> Perform Neville algorithm to interpolate
  call p_InterpNeville3D ( xp, ierr )

  !> Copy result to output array
  uf(1) = interpolationvalues(1,1,1,1)
  uf(2) = interpolationvalues(2,1,1,1)
  uf(3) = interpolationvalues(3,1,1,1)

! write(*,*) 'Position: ', xp, ' Velocity: ', uf

end subroutine p_InterpInterpolate3D
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_InterpFindMolecule ( xp, i0, j0, k0, ierr )

  use g_constants,    only : HALF
  use g_parameters,   only : P_INTERP_ORDER, NHALO
  use g_domain,       only : XMIN_MOVINGFRAME, YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                             XMAX_MOVINGFRAME, YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME, &
                             DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME, &
                             LX_HALO, LY_HALO, LZ_HALO
  use f_arrays,       only : x_min_2r, y_min_2r, z_min, x_max_2r, y_max_2r, z_max

  implicit none
  real(kind=C_DOUBLE),intent(in)     :: xp(3)
  integer,intent(out)                :: i0, j0, k0
  integer                            :: i, j, k
  integer                            :: m
  PetscErrorCode,intent(inout)              :: ierr
  real(kind=C_DOUBLE)                :: xdiff1, xdiff2, xdiff3


  IFODDEVEN: if ( modulo(P_INTERP_ORDER,2) == 0 ) then
!    xdiff1 = xp(1) + DX_MOVINGFRAME
    xdiff1 = xp(1) + LX_HALO + DX_MOVINGFRAME + ( HALF * DX_MOVINGFRAME )
!    xdiff2 = xp(2) + DY_MOVINGFRAME
    xdiff2 = xp(2) + LY_HALO + DY_MOVINGFRAME + ( HALF * DY_MOVINGFRAME )
    xdiff3 = xp(3) + LZ_HALO + DZ_MOVINGFRAME + ( HALF * DZ_MOVINGFRAME )
  else
!    xdiff1 = xp(1) + ( HALF * DX_MOVINGFRAME )
    xdiff1 = xp(1) + LX_HALO + ( HALF * DX_MOVINGFRAME )
!    xdiff2 = xp(2) + ( HALF * DY_MOVINGFRAME )
    xdiff2 = xp(2) + LY_HALO + ( HALF * DY_MOVINGFRAME )
    xdiff3 = xp(3) + LZ_HALO + ( HALF * DZ_MOVINGFRAME )
  end if IFODDEVEN

  !> determine lower corner of interpolation molecule
  i = int(xdiff1/DX_MOVINGFRAME) - NHALO - NHALO
  j = int(xdiff2/DY_MOVINGFRAME) - NHALO - NHALO
  k = int(xdiff3/DZ_MOVINGFRAME) - NHALO - NHALO

  i0 = i
  j0 = j
  k0 = k

  DOX: do m = 1, P_INTERP_ORDER, 1
    interpolationx(m) = real(i,kind=C_DOUBLE) * DX_MOVINGFRAME
    i = i + 1
  end do DOX

  DOY: do m = 1, P_INTERP_ORDER, 1
    interpolationy(m) = real(j,kind=C_DOUBLE) * DY_MOVINGFRAME
    j = j + 1
  end do DOY

  DOZ: do m = 1, P_INTERP_ORDER, 1
    interpolationz(m) = real(k,kind=C_DOUBLE) * DZ_MOVINGFRAME
    k = k + 1
  end do DOZ

!  write(*,*) 'xp: ', xp(1), 'x: ', interpolationx, 'xmin, xmax: ', XMIN_MOVINGFRAME, XMAX_MOVINGFRAME
!  write(*,*) 'yp: ', xp(2), 'y: ', interpolationy, 'ymin, ymax: ', YMIN_MOVINGFRAME, YMAX_MOVINGFRAME
!  write(*,*) 'zp: ', xp(3), 'z: ', interpolationz, 'zmin, zmax: ', ZMIN_MOVINGFRAME, ZMAX_MOVINGFRAME

end subroutine p_InterpFindMolecule
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_InterpNeville3D
!> Neville's algorithm as described in Numerical Recipes,
!! but implemented in three directions consecutively
subroutine p_InterpNeville3D ( xp, ierr )

  use g_parameters,   only : P_INTERP_ORDER

  implicit none

  !> Coordinates of the particle
  real(kind=C_DOUBLE),dimension(1:3)     :: xp
  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_DOUBLE)                    :: xdiff1
  real(kind=C_DOUBLE)                    :: xdiff2
  real(kind=C_DOUBLE)                    :: xdiff3
  real(kind=C_DOUBLE)                    :: zdiff1
  real(kind=C_DOUBLE)                    :: zdiff2
  real(kind=C_DOUBLE)                    :: zdiff3
  real(kind=C_DOUBLE)                    :: ydiff1
  real(kind=C_DOUBLE)                    :: ydiff2
  real(kind=C_DOUBLE)                    :: ydiff3
  real(kind=C_DOUBLE),dimension(1:3)     :: P
  integer                                :: column
  integer                                :: nvals
  integer                                :: xindex
  integer                                :: zindex

!  write(*,*) 'Interpolation molecule (', xp, '): ', interpolationvalues

  !> Go through tableau to compute Lagrange polynomials recursively
  DOYLEFTTORIGHT: do column = 1, P_INTERP_ORDER - 1, 1
    !> Go through each column
    DOYCOLUMN: do nvals = 1, (P_INTERP_ORDER - column), 1
      !> Compute distances between particle position and Eulerian velocity as needed for interpolation
      ydiff1 = xp(2) - interpolationy(nvals + column)
      ydiff2 = interpolationy(nvals) - xp(2)
      ydiff3 = interpolationy(nvals) - interpolationy(nvals + column)
      ! P = ( (xp - x_(i+m)) * value(j) + (x_i - xp) * value(j+1) ) / ( x_i - x_(i+m) )
      !> Go through rows of the interpolation molecule and update each column value
      DOYZ: do zindex = 1, P_INTERP_ORDER, 1
        DOYX: do xindex = 1, P_INTERP_ORDER, 1
          P = ( ( ydiff1 * interpolationvalues(:,xindex,zindex,nvals) ) &
              + ( ydiff2 * interpolationvalues(:,xindex,zindex,nvals+1) ) ) &
              / ydiff3
          interpolationvalues(:,xindex,zindex,nvals) = P
        end do DOYX
      end do DOYZ
    end do DOYCOLUMN
  end do DOYLEFTTORIGHT

  !> Go through tableau to compute Lagrange polynomials recursively
  DOZLEFTTORIGHT: do column = 1, P_INTERP_ORDER - 1, 1
    !> Go through each column
    DOZCOLUMN: do nvals = 1, (P_INTERP_ORDER - column), 1
      !> Compute distances between particle position and Eulerian velocity as needed for interpolation
      zdiff1 = xp(3) - interpolationz(nvals + column)
      zdiff2 = interpolationz(nvals) - xp(3)
      zdiff3 = interpolationz(nvals) - interpolationz(nvals + column)
      ! P = ( (xp - x_(i+m)) * value(j) + (x_i - xp) * value(j+1) ) / ( x_i - x_(i+m) )
      !> Go through rows of the interpolation molecule and update each column value
      DOZX: do xindex = 1, P_INTERP_ORDER, 1
        P = ( ( zdiff1 * interpolationvalues(:,xindex,nvals,1) ) &
            + ( zdiff2 * interpolationvalues(:,xindex,nvals+1,1) ) ) &
            / zdiff3
        interpolationvalues(:,xindex,nvals,1) = P
      end do DOZX
    end do DOZCOLUMN
  end do DOZLEFTTORIGHT

  !> Go through tableau again for the x direction
  DOXLEFTTORIGHT: do column = 1, P_INTERP_ORDER - 1, 1
    !> Go through each column
    DOXCOLUMN: do nvals = 1, (P_INTERP_ORDER - column), 1
      !> Compute distances between particle position and Eulerian velocity as needed for interpolation
      xdiff1 = xp(1) - interpolationx(nvals + column)
      xdiff2 = interpolationx(nvals) - xp(1)
      xdiff3 = interpolationx(nvals) - interpolationx(nvals + column)
      P = ( ( xdiff1 * interpolationvalues(:,nvals,1,1) ) &
          + ( xdiff2 * interpolationvalues(:,nvals+1,1,1) ) ) &
          / xdiff3
      interpolationvalues(:,nvals,1,1) = P
    end do DOXCOLUMN
  end do DOXLEFTTORIGHT

  !> Now interpolationvalues(1,1) is the value at the position on the x-y plane.

!  write(*,*) 'Interpolated value (', xp, '): ', interpolationvalues(:,1,1,1)

end subroutine p_InterpNeville3D
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! p_InterpInterpolate2D 
!subroutine p_InterpInterpolate2D ( npin, ptr_positions, ptr_values, ptr_status, zslice, ierr )    
!  use g_parameters, only : MYRANK, NHALO, P_INTERP_ORDER 
!  use g_domain,     only : DX_MOVINGFRAME, XMIN_MOVINGFRAME, XMAX_MOVINGFRAME, &
!                           DY_MOVINGFRAME, YMIN_MOVINGFRAME, YMAX_MOVINGFRAME
!  use f_arrays,     only : arr_u1_2r, arr_u2_2r, arr_u3_2r, &
!                           x_min_2r, x_max_2r, y_min_2r, y_max_2r
!  use p_arrays,     only : particle_xp, particle_xp_halo_min, particle_xp_halo_max, &
!                           particle_xp_ghost_min, particle_xp_ghost_max

!  implicit none

!  integer,intent(in)     :: npin
!  type(C_PTR),intent(in) :: ptr_positions
!  type(C_PTR),intent(in) :: ptr_values
!  type(C_PTR),intent(in) :: ptr_status
!  integer,intent(in)     :: zslice
!  PetscErrorCode,intent(inout)  :: ierr


!  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
!  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particlevalues
!  integer,pointer,dimension(:,:)             :: particlestatus
!  real(kind=C_DOUBLE),dimension(1:3)         :: xp
!  real(kind=C_DOUBLE),dimension(1:3)         :: uf
!  integer                                    :: np
!  integer                                    :: i0
!  integer                                    :: j0
!  integer                                    :: k0
!  integer                                    :: kp

  !> Associate position pointer array with target array.
!  call c_f_pointer ( ptr_positions, particleposition, [3, npin])
  !> Associate interpolation values pointer array with target array.
!  call c_f_pointer ( ptr_values, particlevalues, [3, npin])
  !> Associate status pointer array with target array.
!  call c_f_pointer ( ptr_status, particlestatus, [3, npin])

!  DOINTERP: do np = 1, npin, 1
!    IFONSLICE: if ( particlestatus(2,np) == zslice ) then
      !> Get particle coordinates
!      xp = particleposition(:,np)
      !> Find interpolation molecule
!      call p_InterpFindMolecule ( xp, i0, j0, k0, ierr )
!      kp = zslice - k0 + 1
!      write(*,*) 'Interpolation molecule: ', i0, j0, k0, 'Current slice: ', zslice, 'z position within molecule: ', kp
      !> Interpolate velocity components one by one
!      interpolationvalues = arr_u1_2r(i0:(i0+P_INTERP_ORDER-1),j0:(j0+P_INTERP_ORDER-1))
!      call p_InterpNeville2D ( xp, ierr )
!      interpolationuf(1) = interpolationvalues(1,1)
!      interpolationvalues = arr_u2_2r(i0:(i0+P_INTERP_ORDER-1),j0:(j0+P_INTERP_ORDER-1))
!      call p_InterpNeville2D ( xp, ierr )
!      interpolationuf(2) = interpolationvalues(1,1)
!      interpolationvalues = arr_u3_2r(i0:(i0+P_INTERP_ORDER-1),j0:(j0+P_INTERP_ORDER-1))
!      call p_InterpNeville2D ( xp, ierr )
!      interpolationuf(3) = interpolationvalues(1,1)
!      write(*,*) 'Position: ', xp, '2D Velocity: ', interpolationuf
      !> Update interpolation in third direction
!      call p_InterpUpdate3D ( xp, kp, uf, ierr )
!      particlevalues(:,np) = uf
!      write(*,*) 'Position: ', xp, '3D Velocity: ', particlevalues(:,np)
      !> Update particle status if still need next slice for interpolation
!      IFNEXTSLICE: if ( ( particlestatus(2,np) ) < ( k0 + P_INTERP_ORDER - 1 ) ) then
!        particlestatus(2,np) = particlestatus(2,np) + 1
!      end if IFNEXTSLICE
!    end if IFONSLICE
!  end do DOINTERP

!end subroutine p_InterpInterpolate2D
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! p_InterpAddHaloGhost 
!subroutine p_InterpAddHaloGhost ( nphalomin, nphalomax, ierr )    
!  use g_parameters, only : MYRANK, NHALO, P_INTERP_ORDER 
!  use g_domain,     only : DX_MOVINGFRAME, XMIN_MOVINGFRAME, XMAX_MOVINGFRAME, &
!                           DY_MOVINGFRAME, YMIN_MOVINGFRAME, YMAX_MOVINGFRAME
!  use f_arrays,     only : arr_u1_2r, arr_u2_2r, arr_u3_2r, &
!                           x_min_2r, x_max_2r, y_min_2r, y_max_2r
!  use p_arrays,     only : particle_xp_halo_min, particle_xp_halo_max, &
!                           particle_uf_halo_min, particle_uf_halo_max, &
!                           particle_xp_ghost_min, particle_xp_ghost_max

!  implicit none

!  integer,intent(in)     :: nphalomin
!  integer,intent(in)     :: nphalomax
!  PetscErrorCode,intent(inout)  :: ierr

!  integer                :: np

!  write(*,*) 'halo: ', size(particle_uf_halo_min), 'ghost: ',  size(particle_xp_ghost_min)
!  DOMIN: do np = 1, nphalomin, 1
!    write(*,*) 'Position: ', particle_xp_halo_min(:,np), 'Halo: ', particle_uf_halo_min(:,np), &
!               'Ghost: ', particle_xp_ghost_min(:,np)
!    particle_uf_halo_min(:,np) = particle_uf_halo_min(:,np) + particle_xp_ghost_min(:,np)
!  end do DOMIN

!  write(*,*) 'halo: ', size(particle_uf_halo_max), 'ghost: ',  size(particle_xp_ghost_max)
!  DOMAX: do np = 1, nphalomax, 1
!    write(*,*) 'Position: ', particle_xp_halo_max(:,np), 'Halo: ', particle_uf_halo_max(:,np), &
!               'Ghost: ', particle_xp_ghost_max(:,np)
!    particle_uf_halo_max(:,np) = particle_uf_halo_max(:,np) + particle_xp_ghost_max(:,np)
!  end do DOMAX

!end subroutine p_InterpAddHaloGhost
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_InterpNeville2D
!> Neville's algorithm as described in Numerical Recipes,
!! but implemented in two directions consecutively
!subroutine p_InterpNeville2D ( xp, ierr )

!  use g_parameters,   only : P_INTERP_ORDER

!  implicit none

  !> Coordinates of the particle
!  real(kind=C_DOUBLE),dimension(1:3)     :: xp
!  PetscErrorCode,intent(inout)                  :: ierr

!  real(kind=C_DOUBLE)                    :: xdiff1
!  real(kind=C_DOUBLE)                    :: xdiff2
!  real(kind=C_DOUBLE)                    :: xdiff3
!  real(kind=C_DOUBLE)                    :: ydiff1
!  real(kind=C_DOUBLE)                    :: ydiff2
!  real(kind=C_DOUBLE)                    :: ydiff3
!  real(kind=C_DOUBLE)                    :: P
!  integer                                :: column
!  integer                                :: nvals
!  integer                                :: xindex

!  write(*,*) 'Interpolation molecule (', xp, '): ', interpolationvalues

  !> Go through tableau to compute Lagrange polynomials recursively
!  DOYLEFTTORIGHT: do column = 1, P_INTERP_ORDER - 1, 1
    !> Go through each column
!    DOYCOLUMN: do nvals = 1, (P_INTERP_ORDER - column), 1
      !> Compute distances between particle position and Eulerian velocity as needed for interpolation
!      ydiff1 = xp(2) - interpolationy(nvals + column)
!      ydiff2 = interpolationy(nvals) - xp(2)
!      ydiff3 = interpolationy(nvals) - interpolationy(nvals + column)
      ! P = ( (xp - x_(i+m)) * value(j) + (x_i - xp) * value(j+1) ) / ( x_i - x_(i+m) )
      !> Go through rows of the interpolation molecule and update each column value
!      DOX: do xindex = 1, P_INTERP_ORDER, 1
!        P = ( ( ydiff1 * interpolationvalues(xindex,nvals) ) &
!            + ( ydiff2 * interpolationvalues(xindex,nvals+1) ) ) &
!            / ydiff3
!        interpolationvalues(xindex,nvals) = P
!      end do DOX
!    end do DOYCOLUMN
!  end do DOYLEFTTORIGHT

  !> Go through tableau again for the x direction
!  DOXLEFTTORIGHT: do column = 1, P_INTERP_ORDER - 1, 1
    !> Go through each column
!    DOXCOLUMN: do nvals = 1, (P_INTERP_ORDER - column), 1
      !> Compute distances between particle position and Eulerian velocity as needed for interpolation
!      xdiff1 = xp(1) - interpolationx(nvals + column)
!      xdiff2 = interpolationx(nvals) - xp(1)
!      xdiff3 = interpolationx(nvals) - interpolationx(nvals + column)
!      P = ( ( xdiff1 * interpolationvalues(nvals,1) ) &
!          + ( xdiff2 * interpolationvalues(nvals+1,1) ) ) &
!          / xdiff3
!      interpolationvalues(nvals,1) = P
!    end do DOXCOLUMN
!  end do DOXLEFTTORIGHT

  !> Now interpolationvalues(1,1) is the value at the position on the x-y plane.

!  write(*,*) 'Interpolated value (', xp, '): ', interpolationvalues(1,1)

!end subroutine p_InterpNeville2D
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_InterpUpdate3D
!> Update the interpolation polynomial in the 3rd (z) direction.
!subroutine p_InterpUpdate3D ( xp, kp, uf, ierr )

!  use g_parameters,   only : P_INTERP_ORDER

!  implicit none

!  real(kind=C_DOUBLE),dimension(1:3),intent(in)    :: xp
!  integer,intent(in)                               :: kp            ! Interpolation node of current z slice
!  real(kind=C_DOUBLE),dimension(1:3),intent(inout) :: uf
!  PetscErrorCode,intent(inout)                            :: ierr

!  real(kind=C_DOUBLE)                              :: poly_lagrange

  !> Compute the Lagrange polynomial at the z position.
!  poly_lagrange = p_InterpPolyLagrange ( xp(3), interpolationz, kp, P_INTERP_ORDER ) 
!  write(*,*) 'uf: ', uf, 'polynomial: ', poly_lagrange, 'interpolationuf: ', interpolationuf
!  uf = uf + poly_lagrange * interpolationuf

!end subroutine p_InterpUpdate3D
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
function p_InterpPolyLagrange (zpin, znodein, jin, orderin ) result ( poly_lagrange )
    
    !
    ! works our the coefficient products for each of the node values
    !

  use g_constants, only : ONE

  implicit none

  real(kind=C_DOUBLE)                             :: poly_lagrange  
  real(kind=C_DOUBLE),intent(in)                  :: zpin          ! particle position, relative to xnode(1) 
  real(kind=C_DOUBLE),intent(in),dimension(:)     :: znodein       ! 1d list of node positions (may be irregular)
  integer,intent(in)                              :: jin           ! coeft for the jth node on the lagrange formula
  integer,intent(in)                              :: orderin       ! number of nodes on the lagrange polynomial
  integer                                         :: k
    
  poly_lagrange = ONE

!  write(*,*) 'znodein: ', znodein

  DOPOLY: do k = 1, orderin, 1
    IFSELF: if ( k == jin ) then
      cycle
    end if IFSELF
    poly_lagrange = poly_lagrange * ( (zpin - znodein(k) ) / ( znodein(jin) - znodein(k) ) )
!    write(*,*) 'k, jin: ', k, jin, 'poly_lagrange: ', poly_lagrange, 'zpin: ', zpin, &
!               'znodein(k): ', znodein(k), 'znodein(jin): ', znodein(jin)
  end do DOPOLY

end function p_InterpPolyLagrange
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
!subroutine p_InterpFindMolecule1D ( npin, ptr_positions, ptr_status, ierr )    

!  use g_constants,    only : HALF, ZERO
!  use g_parameters,   only : P_INTERP_ORDER, NHALO
!  use g_domain,       only : ZMIN_MOVINGFRAME, ZMAX_MOVINGFRAME, DZ_MOVINGFRAME, LZ_HALO
!  use f_arrays,       only : z_min, z_max

!  implicit none
!
!  integer,intent(in)     :: npin
!  type(C_PTR),intent(in) :: ptr_positions
!  type(C_PTR),intent(in) :: ptr_status
!  PetscErrorCode,intent(inout)  :: ierr

!  real(kind=C_DOUBLE)                :: xdiff3

!  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particleposition
!  integer,pointer,dimension(:,:)             :: particlestatus
!  real(kind=C_DOUBLE),dimension(1:3)         :: xp
!  integer                                    :: np

  !> Associate position pointer array with target array.
!  call c_f_pointer ( ptr_positions, particleposition, [3, npin])
  !> Associate status pointer array with target array.
!  call c_f_pointer ( ptr_status, particlestatus, [3, npin])

!  DOINTERP: do np = 1, npin, 1
    !> Get particle coordinates
!    xp = particleposition(:,np)
    !> Find interpolation molecule
    !> For particles the fluid grid points are in the cell centre. 
    !! Therefore a shift is needed for even interpolation molecules,
    !! whereas for odd molecules everything is in the right place.
!    IFODDEVEN: if ( modulo(P_INTERP_ORDER,2) == 0 ) then
!      xdiff3 = xp(3) + LZ_HALO + DZ_MOVINGFRAME + ( HALF * DZ_MOVINGFRAME )
!    else
!      xdiff3 = xp(3) + LZ_HALO + ( HALF * DZ_MOVINGFRAME )
!    end if IFODDEVEN

    !> determine lowest slice of interpolation molecule
!    particlestatus(2,np) = int(xdiff3/DZ_MOVINGFRAME) - NHALO - NHALO
!    particlestatus(3,np) = particlestatus(2,np) + P_INTERP_ORDER - 1 

!    IFLOW: if ( particlestatus(2,np) < z_min ) then
!      IFLOWINSIDE: if ( particlestatus(3,np) >= z_min ) then
!        particlestatus(2,np) = z_min
!      end if IFLOWINSIDE
!   end if IFLOW

!    IFHIGH: if ( particlestatus(3,np) > z_max ) then
!      IFHIGHINSIDE: if ( particlestatus(2,np) <= z_max ) then
!        particlestatus(3,np) = z_max
!      end if IFHIGHINSIDE
!    end if IFHIGH

!    write(*,*) 'Position: ', xp(3), 'halo: ', LZ_HALO, 'Status: ', particlestatus(:,np), &
!               'zmin, zmax: ', z_min, z_max, 'xdiff: ', xdiff3
!  end do DOINTERP

!end subroutine p_InterpFindMolecule1D
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
!subroutine interp_uf(order,xparticle,phi,method,field_type)

!    use d_arrays,    only : fluid_x
!    use f_arrays,    only : fluid_u
!    use p_control,    only : part_efieldr, part_potr   ! real-space electric field & potential
!    use g_constants, only : two

!    implicit none
!    integer                ,intent(in)                      :: order       ! == 1,2,3,4    
!    real(kind=C_DOUBLE),intent(in) ,dimension(:)            :: xparticle   ! particle position
!    real(kind=C_DOUBLE),intent(out),dimension(:)            :: phi         ! fluid u at particle pos
!    integer                ,intent(in)                      :: method      ! type of interpolation
!    integer                ,intent(in)                      :: field_type
!    integer                                                 :: error
!    integer                                                 :: dir
!    integer                                                 :: oshift
!    integer                ,save,dimension(1:3)             :: ijk
!    integer                ,save,dimension(1:3)             :: irxc
!    integer                ,save,dimension(1:3)             :: ilxc
!    integer                ,save,dimension(1:3)             :: im
!    integer                ,save,dimension(1:3)             :: ip
!    real(kind=C_DOUBLE),save,dimension(1:3)                 :: xc
!    real(kind=C_DOUBLE),save,dimension(1:4,1:4,1:4,1:3)     :: lvertex  ! saved between calls
!    real(kind=C_DOUBLE),save,dimension(1:4,1:4,1:4,1:3)     :: lfluidu  ! saved between calls
!    real(kind=C_DOUBLE)     ,dimension(1:3)                 :: xp       ! relative particle position vector
!    real(kind=C_DOUBLE)     ,dimension(1:3)                 :: dx       ! node spacing


!    call interp_find(xparticle,ijk(1),ijk(2),ijk(3))

    ! stage 1
    ! construct the 3d interpolation molecule - depends on order and
    ! for even nunber of nodes which side of the cell centre the particle is (ie left or right)

!    oshift=0  ! base size of the interpolation block, +/- either side of the particle cv

!    if(order==3.or.order==4)oshift=1

!    if(order==2.or.order==4)then

        ! flag for 2,4 block interpolations to find out of the space cell is to 
        ! the ent or the wsb side of the block

!        irxc(:)=1 ! default is to th ent side
!        ilxc(:)=0

!        xc(1) = fluid_x(ijk(1)-1,ijk(2)  ,ijk(3)  ,1)+((fluid_x(ijk(1),ijk(2),ijk(3),1)-fluid_x(ijk(1)-1,ijk(2)  ,ijk(3)  ,1))/two)
!        xc(2) = fluid_x(ijk(1)  ,ijk(2)-1,ijk(3)  ,2)+((fluid_x(ijk(1),ijk(2),ijk(3),2)-fluid_x(ijk(1)-1,ijk(2)-1,ijk(3)  ,2))/two)
!        xc(3) = fluid_x(ijk(1)  ,ijk(2)  ,ijk(3)-1,3)+((fluid_x(ijk(1),ijk(2),ijk(3),3)-fluid_x(ijk(1)-1,ijk(2)  ,ijk(3)-1,3))/two)

        ! note: choosing the interpolation side based on position not velocity

!        do dir=1,3
!            if(xparticle(dir).lt.xc(dir))then
!                irxc(dir) =  0
!                ilxc(dir) = -1
!            end if
!        end do
!    end if

!    do dir=1,3
!        im(dir) = ijk(dir) - oshift + ilxc(dir)
!        ip(dir) = ijk(dir) + oshift + irxc(dir)

!    end do

    ! define the ent position of the cell corner and the information block

!    lvertex(1:order,1:order,1:order,1:3) = fluid_x(im(1):ip(1),im(2):ip(2),im(3):ip(3),1:3)
!    if(field_type.eq.field_vel)then
!        lfluidu(1:order,1:order,1:order,1:3) = fluid_u(im(1):ip(1),im(2):ip(2),im(3):ip(3),1:3)
!    else if(field_type.eq.field_el)then
!        lfluidu(1:order,1:order,1:order,1:3) = part_efieldr(im(1):ip(1),im(2):ip(2),im(3):ip(3),1:3)
!    else if(field_type.eq.field_ep)then
!        lfluidu(1:order,1:order,1:order,1) = part_potr(im(1):ip(1),im(2):ip(2),im(3):ip(3),1)
!    else
!        write(*,*)'error: invalid interpolation field flag'
!    end if

!    dx(1) = lvertex(2,1,1,1)-lvertex(1,1,1,1) ! define the node spacing. here the mesh width is the node spacing
!    dx(2) = lvertex(1,2,1,2)-lvertex(1,1,1,2) ! varies between coord directions, but uniform within
!    dx(3) = lvertex(1,1,2,3)-lvertex(1,1,1,3) ! assumes coord axis and mesh axis are identical

    ! stage 2
    ! find the fluid velocity

!    if(order==1)then

!        phi(1:3)=lfluidu(1,1,1,1:3)  ! not so hard this one...

!    else if(order/=1.and.method==lagrange)then

        ! define the particle position, relative to the first node, then
        ! ensure lvertex is refereced to a local x(1,1,1)=(0,0,0) system, consistent
        ! with the relatative particle position

!        do dir=1,3
!            xp(dir) =     xparticle(dir) !- (lvertex(1,1,1,dir)-dx(dir)/two)
!            lvertex(:,:,:,dir) = lvertex(:,:,:,dir) -dx(dir)/two
!        end do

!        call interp_lagrange(lfluidu,lvertex,xp,order,phi)

!    end if

!end subroutine
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!subroutine interp_lagrange(fluidu,fluidx,xparticle,orderin,phi)

    !
    ! preforms a lagrange interpolation of a velocity to a particle position.
    !

!    use g_constants, only : zero

!    implicit none
    
!    real(kind=C_DOUBLE),intent(in),dimension(:,:,:,:)        :: fluidu    ! cell centre velocity (i,j,k,1:3)
!    real(kind=C_DOUBLE),intent(in),dimension(:,:,:,:)        :: fluidx    ! cell centre position (i,j,k,1:3) 
!    real(kind=C_DOUBLE),intent(in),dimension(:)              :: xparticle ! particle position
!    integer                ,intent(in)                       :: orderin   ! size of interpolation molecule
!    real(kind=C_DOUBLE),intent(out),dimension(:)             :: phi       ! fluid velocity at xparticle
!    real(kind=C_DOUBLE),save,dimension(1:4,1:4,1:3)          :: phi_x2x3  ! values interpolated to the x2x3 plane at x1=xp(1)
!    real(kind=C_DOUBLE),save,dimension(1:4,1:3)              :: phi_x2    ! values interpolated to vector on x2x3 at x2=xp(2)
!    integer,save :: j,k,dir,cdir,n
    
!    cdir=1

    ! stage 1
    ! perform 12*3 interpolations in the 1(x) direction, maps 12 values of each fluid
    ! velocity component onto the x2-x3 plane at x=xp(1)
    
!    phi_x2x3(:,:,:) = zero

!    do dir=1,3
!        do j=1,orderin
!            do k=1,orderin
!                do n=1,orderin ! note n must be the inner loop
!                    phi_x2x3(j,k,dir) = phi_x2x3(j,k,dir) + fluidu(n,j,k,dir) * poly_lagrange(xparticle(1),fluidx(:,cdir,cdir,1),n,orderin)
!                end do
!            end do
!        end do
!    end do
    
    ! stage 2
    ! perform 4*3 interpolations in the 2(y) direction. maps 4 values of each fluid
    ! velocity component onto a vector lying on the x2-x3 plane at x=xp(1) and passing
    ! through the point y=xp(2), pointing in the z direction
    
!    phi_x2(:,:) = zero
!    do dir=1,3
!        do k=1,orderin
!            do n=1,orderin
!                phi_x2(k,dir) = phi_x2(k,dir) + phi_x2x3(n,k,dir) * poly_lagrange(xparticle(2),fluidx(cdir,:,cdir,2),n,orderin)
!            end do
!        end do
!    end do
    
    ! stage 3
    ! perform the interpolation in the 3(z) direction to get phi
    
!    phi(:)=zero
!    do dir=1,3
!        do n=1,orderin
!            phi(dir) = phi(dir) + phi_x2(n,dir) * poly_lagrange(xparticle(3),fluidx(cdir,cdir,:,3),n,orderin)
!        end do
!    end do

!end subroutine interp_lagrange
!---------------------------------------------------------------------------------


end module p_interp
