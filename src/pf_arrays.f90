!----------------------------------------------------------------------
! module f_arrays
!> Contains all subroutines to allocate and deallocate the fluid arrays.
!----------------------------------------------------------------------
! Definition of the PETSc DMDAs (Distributed Arrays).
! These define the data structure, but do not contain the data. They 
! serve as templates for global and local vectors for the actual data.
module f_arrays

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscdm.h>
!#include <petsc/finclude/petscdmda.h>
use g_petsc
!use petscsys
!use petscdm
!use petscdmda
!use petscvec
!use petscdmlabel     

use iso_c_binding
implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

! var da3w 
!  array with three wave-space dimensions, da1rw2 one real-
! space and two wave-space dimensions, and da1r1w one real-space and
! one wave-space dimension.
! @var kjcomm
! @var ijmat 
! buffer after FFT in i direction
! @var jimat 
! transpose of ijmat, so FFT can be performed in j direction
! @var jkmat 
! buffer after FFT in j direction
! @var kjmat 
! transpose of ijmat, so FFT can be performed in k direction

!>  array with three wave-space dimensions
DM,save,public                       :: da3w 
!>  array with one real-space and two wave-space dimensions
DM,save,public                       :: da1r2w 
!>  array with three real-space dimensions
DM,save,public                       :: da3r 
!>  array with three wave-space dimensions used for restart files
DM,save,public                       :: da3w_rs
!>  array with two wave-space dimensions
DM,save,public                       :: da2w 
!>  array with two real-space dimensions
DM,save,public                       :: da2r 
!>  array with one real-space and one wave-space dimension
!! for evaluating 1D spectra
DM,save,public                       :: da1r1w 
MPI_Comm,save,public              :: kjcomm, jicomm 
MPI_Comm,save,public              :: restartcomm 
integer,save,public                  :: restartcolour
!> communicators for processes in the kj and ji planes

Mat,save,public                                  :: kjmat

!> Global vectors for the 3D arrays
Vec,save,public                                  :: u1_3w, u2_3w, u3_3w, & 
                                                    u1_1r2w, u2_1r2w, u3_1r2w, & 
                                                    u4_1r2w, u5_1r2w, u6_1r2w, &
                                                    u7_1r2w, u8_1r2w, u9_1r2w, &
                                                    u10_1r2w, u11_1r2w, u12_1r2w, &
                                                    rk3_3w_1, rk3_3w_2, rk3_3w_3

Vec,save,public                                  :: u1_3r 
Vec,save,public                                  :: u2_3r 
Vec,save,public                                  :: u3_3r 
Vec,save,public                                  :: u1_3r_local 
Vec,save,public                                  :: u2_3r_local 
Vec,save,public                                  :: u3_3r_local 

!> Global vectors for the 2D arrays
Vec,save,public                                  :: u1_2w, u2_2w, u3_2w, &
                                                    u4_2w, u5_2w, u6_2w, &
                                                    u7_2w, u8_2w, u9_2w, &
                                                    u10_2w, u11_2w, u12_2w, &
                                                    u1_2r, u2_2r, u3_2r, &
                                                    u4_2r, u5_2r, u6_2r, &
                                                    u7_2r, u8_2r, u9_2r, &
                                                    u10_2r, u11_2r, u12_2r

!> Global vectors for 2D arrays only used for 1D spectra
Vec,save,public                                  :: u1_1r1w, u2_1r1w, u3_1r1w

!> Global vectors for fluid restart file only 
Vec,save,public                                  :: u1_3w_rs, u2_3w_rs, u3_3w_rs

!> Arrays for direct access to vectors
PetscScalar,pointer,public,dimension(:,:,:,:)    :: arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                                                            arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                                                            arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                                                            arr_u7_1r2w, arr_u8_1r2w, arr_u9_1r2w, &
                                                            arr_u10_1r2w, arr_u11_1r2w, arr_u12_1r2w, &
                                                            arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u1_3r
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u2_3r
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u3_3r
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u1_3r_local
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u2_3r_local
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u3_3r_local
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                                                            arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                                                            arr_u7_2w, arr_u8_2w, arr_u9_2w, &
                                                            arr_u10_2w, arr_u11_2w, arr_u12_2w
PetscScalar,pointer,public,dimension(:,:)        :: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                                                            arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                                                            arr_u7_2r, arr_u8_2r, arr_u9_2r, &
                                                            arr_u10_2r, arr_u11_2r, arr_u12_2r
PetscScalar,pointer,public,dimension(:,:,:)      :: arr_u1_1r1w, arr_u2_1r1w, arr_u3_1r1w
PetscScalar,allocatable,public,dimension(:,:)    :: arr_u1_remesh, arr_u2_remesh, arr_u3_remesh
PetscScalar,allocatable,public,dimension(:,:,:)  :: arr_remesh

!> Arrays for restart file only
PetscScalar,pointer,public,dimension(:, :, :, :) :: arr_u1_3w_rs, arr_u2_3w_rs, arr_u3_3w_rs
                                                    
!> Corners of local portion of the domain
integer,public    :: i_min_3w, j_min_3w, k_min_3w, i_min_1r2w, j_min_1r2w, z_min
integer,public    :: i_max_3w, j_max_3w, k_max_3w, i_max_1r2w, j_max_1r2w, z_max
integer,public    :: i_width_3w, j_width_3w, k_width_3w, z_width_3w, &
                     i_width_1r2w, j_width_1r2w, k_width_1r2w, z_width_1r2w, &
                     x_min_2r, y_min_2r, x_max_2r, y_max_2r, x_width_2r, y_width_2r, &
                     i_min_1r1w, y_min_1r1w, i_max_1r1w, &
                     y_max_1r1w, i_width_1r1w, y_width_1r1w
integer,public    :: x_min_ghost, y_min_ghost, z_min_ghost, &
                     x_width_ghost, y_width_ghost, z_width_ghost, &
                     x_max_ghost, y_max_ghost, z_max_ghost
!> Local sizes on all processes
integer,public                          :: nprocsi_3w, nprocsj_3w, &
                                           nprocsi_1r2w, nprocsz_1r2w, nprocsy_2r
integer,dimension(:),allocatable,public :: j_allwidths_3w, i_allwidths_1r2w, z_allwidths_1r2w, &
                                           y_allwidths_2r
integer,public                          :: j_maxwidth_3w, i_maxwidth_1r2w, z_maxwidth_1r2w, &
                                           y_maxwidth_2r
!> Global size of transpose buffer
integer,public                          :: nodes_kj_transpose, nodes_z_transpose, &
                                           nodes_ki_transpose, nodes_y_transpose 

!---------------------------------------------------------------------------------

public :: f_ArraysInit
public :: f_ArraysFinalise

public :: f_ArraysWavenumberPositive
public :: f_ArraysWavenumberPositiveLong
public :: f_ArraysWavenumberNegative
public :: f_ArraysWavenumber
public :: f_ArraysWavenumberLong
public :: f_ArraysVectorCreateRestart
public :: f_ArraysDealiasSpherical
public :: f_ArraysLESCutoff

contains
!---------------------------------------------------------------------------------
! subroutine f_ArraysInit
!> @param ierr should return 0
subroutine f_ArraysInit ( ierr )

  use g_parameters, only : MYRANK, MASTER, F_TYPE, F_ISOTROPIC
 
  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  PetscErrorCode        :: perr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising f_arrays module \n', ierr )

!> Call subroutine to create distributed PETSc array structures (DMDA) in 3D
  call PetscPrintf( PETSC_COMM_WORLD, '    | Creating distributed 3D arrays \n', ierr )
  call f_ArraysDMDA3DCreate ( ierr )
!> Determine local part of the domain
  call PetscPrintf( PETSC_COMM_WORLD, '    | Determining local dimensions of 3D arrays \n', ierr )
  call f_ArraysDMDA3DLocalIndices ( ierr )
!> Call subroutine to create distributed PETSc array structures (DMDA) in 2D
  call PetscPrintf( PETSC_COMM_WORLD, '    | Creating distributed 2D arrays \n', ierr )
  call f_ArraysDMDA2DCreate ( ierr )
!> Determine local part of the domain
  call PetscPrintf( PETSC_COMM_WORLD, '    | Determining local dimensions of 2D arrays \n', ierr )
  call f_ArraysDMDA2DLocalIndices ( ierr )
!> Gather information about local dimensions on all processes
  call f_ArraysGatherLocalDimensions ( ierr )
!> Call subroutine to create global and local vectors from the DMDA
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating fluid arrays \n', ierr )
  call f_ArraysVectorCreate ( ierr )
!> Call subroutine to create buffer array for remeshing
  IFSHEAR: if ( F_TYPE /= F_ISOTROPIC ) then
    call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating buffer array for remeshing \n', ierr )
    call f_ArraysRemeshBufferAllocate ( ierr )
  end if IFSHEAR

end subroutine f_ArraysInit
!---------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysFinalise
!> @param ierr should return 0
subroutine f_ArraysFinalise ( ierr )

  use g_parameters, only : MYRANK, MASTER, F_TYPE, F_ISOTROPIC
 
  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  PetscErrorCode        :: perr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Finalising f_arrays module \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating fluid wave-space arrays \n', ierr )

!> Call subroutine to deallocate PETSc arrays for the fluid components
!! in wave space.
  call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating distributed fluid arrays \n', ierr )
  call f_ArraysDMDADealloc ( ierr )
!> Call subroutine to deallocate remesh buffer
  IFSHEAR: if ( F_TYPE /= F_ISOTROPIC ) then
    call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating buffer array for remeshing \n', ierr )
    call f_ArraysRemeshBufferDeallocate ( ierr )
  end if IFSHEAR


end subroutine f_ArraysFinalise
!---------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysDMDA3DCreate
!> Allocate distributed array structures for all Eulerian fields.
!> @param ierr should return 0
subroutine f_ArraysDMDA3DCreate ( ierr )

  use g_parameters, only : V_ALLOCS, YES, MYRANK, MASTER, &
                           NODES_X, NODES_Y, NODES_Z, NHALO, &
                           P_TRACK_PART, YES, &
                           F_RESTART_SCALE
  use g_constants,  only : ZERO , PETSCONE, PETSCTWO, PETSCTHREE 
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK

  implicit none

  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: i

  PetscInt                        :: nkk, nkj, nki, nx, ny, nz
  PetscInt                        :: npi3w
  PetscInt                        :: nh

  PetscErrorCode        :: perr

  nkk = NODES_KK
  nkj = NODES_KJ
  nki = NODES_KI
  nx  = NODES_X
  ny  = NODES_Y
  nz  = NODES_Z

!  call DMDACreate ( PETSC_COMM_WORLD, da3w, ierr )

  !> Create distributed array with three wave-space dimensions
  call DMDACreate3d ( PETSC_COMM_WORLD, &                ! MPI communicator
                      DM_BOUNDARY_PERIODIC, &            ! periodic boundary conditions 
                      DM_BOUNDARY_PERIODIC, &
                      DM_BOUNDARY_PERIODIC, &
                      DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                      nkk, nkj, nki,  &   ! global dimension in each direction
                      PETSCONE, PETSC_DECIDE, PETSC_DECIDE, &   ! number of processors in each direction - let PETSc do this
                      PETSCTWO,                             &   ! degrees of freedom = array components, we use two for real and imaginary components
                      PETSC_NULL_INTEGER,                             &   ! stencil width = 0
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                      da3w,                          &   ! DMDA to be created     
                      ierr )                             ! PETSc error code

  call DMSetUp ( da3w, ierr )

  !> Need to make sure parallelisation in i direction is the same for second array type
  call DMDAGetInfo ( da3w, PETSC_NULL_INTEGER, &                                        ! DMDA, dimension of DMDA
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! global dimensions M, N, P
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, npi3w, &              ! number of procs m, n, p
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! number of dof per node, stencil width, ghost nodes bx
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr ) ! ghost nodes by, bz, stencil type, error code

  nprocsi_3w = int(npi3w)

  !> Create distributed array with one real-space and two wave-space dimensions
  call DMDACreate3d ( PETSC_COMM_WORLD, &                ! MPI communicator
                      DM_BOUNDARY_PERIODIC, &            ! periodic boundary conditions
                      DM_BOUNDARY_PERIODIC, &
                      DM_BOUNDARY_PERIODIC, &
                      DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                      nkj, nz, nki,   &   ! global dimension in each direction
                      PETSCONE, PETSC_DECIDE, npi3w,   &   ! number of processors in each direction - let PETSc do this
                      PETSCTWO,                             &   ! degrees of freedom = array components, we use two for real and imaginary components
                      PETSC_NULL_INTEGER,                             &   ! stencil width = 0
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                      da1r2w,                        &   ! DMDA to be created     
                      ierr )                             ! PETSc error code

  call DMSetUp ( da1r2w, ierr )

  IFPART: if ( P_TRACK_PART == YES ) then

    !> Create distributed array with three real-space dimensions

    !> @todo need to make sure this has the same parallelisation as 2d

    nh = NHALO

    call DMDACreate3d ( PETSC_COMM_WORLD, &                ! MPI communicator
                        DM_BOUNDARY_PERIODIC, &            ! periodic boundary conditions 
                        DM_BOUNDARY_PERIODIC, &
                        DM_BOUNDARY_PERIODIC, &
                        DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                        nx, nz, ny,   &   ! global dimension in each direction
                        PETSCONE, PETSC_DECIDE, npi3w,   &   ! number of processors in each direction - let PETSc do this
                        PETSCONE,                             &   ! degrees of freedom = three velocity components 
                        nh,                         &   ! stencil width = number of halo nodes for particle interpolation
                        PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                        da3r,                          &   ! DMDA to be created     
                        ierr )                             ! PETSc error code
  
    call DMSetUp ( da3r, ierr )

  end if IFPART

  IFVERBOSE: if ( V_ALLOCS == YES ) then
    call PetscPrintf( PETSC_COMM_WORLD, '     | Viewing wave-space distributed array: \n', ierr )
    call DMView ( da3w, PETSC_VIEWER_STDOUT_WORLD, ierr )

    call PetscPrintf( PETSC_COMM_WORLD, '     | Viewing distributed array after one wave-space to real transform: \n', ierr )
    call DMView ( da1r2w, PETSC_VIEWER_STDOUT_WORLD, ierr )

    IFVERBPART: if ( P_TRACK_PART == YES ) then
      call PetscPrintf( PETSC_COMM_WORLD, '     | Viewing real-space distributed array: \n', ierr )
      call DMView ( da3r, PETSC_VIEWER_STDOUT_WORLD, ierr )
    end if IFVERBPART
  end if IFVERBOSE

end subroutine f_ArraysDMDA3DCreate
!----------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_ArraysDMDA3DLocalIndices
!> Call PETSc to find the local start and end points of the domain
!> @param ierr should return 0
subroutine f_ArraysDMDA3DLocalIndices ( ierr )

  use g_parameters, only: MYRANK, V_ALLOCS, P_TRACK_PART, YES
  use g_domain,     only: DZ_MOVINGFRAME, ZMIN_MOVINGFRAME, ZMAX_MOVINGFRAME
  use g_constants,  only: HALF

  implicit none

  PetscErrorCode,intent(inout)   :: ierr

  PetscInt                :: km3w, jm3w, im3w, kw3w, jw3w, iw3w, jm1r2w, zm, &
                             im1r2w, jw1r2w, zw1r2w, iw1r2w
  PetscInt                :: xm3r, ym3r, zm3r, xw3r, yw3r, zw3r

  PetscErrorCode        :: perr

!  call DMDAGetCorners ( da3w, k_min_3w, j_min_3w, i_min_3w, k_width_3w, j_width_3w, i_width_3w, ierr )
  call DMDAGetCorners ( da3w, km3w, jm3w, im3w, kw3w, jw3w, iw3w, ierr )
  k_min_3w = int(km3w) 
  j_min_3w = int(jm3w) 
  i_min_3w = int(im3w) 
  k_width_3w = int(kw3w) 
  j_width_3w = int(jw3w)
  i_width_3w = int(iw3w) 
  k_max_3w = k_min_3w + k_width_3w - 1
  i_max_3w = i_min_3w + i_width_3w - 1
  j_max_3w = j_min_3w + j_width_3w - 1

!  call DMDAGetCorners ( da1r2w, j_min_1r2w, z_min, i_min_1r2w, j_width_1r2w, z_width_1r2w, i_width_1r2w, ierr )
  call DMDAGetCorners ( da1r2w, jm1r2w, zm, im1r2w, jw1r2w, zw1r2w, iw1r2w, ierr )
  j_min_1r2w = int(jm1r2w)
  z_min = int(zm)
  i_min_1r2w = int(im1r2w) 
  j_width_1r2w = int(jw1r2w)
  z_width_1r2w = int(zw1r2w)
  i_width_1r2w = int(iw1r2w)
  z_max = z_min + z_width_1r2w - 1
  ZMIN_MOVINGFRAME = ( real(z_min,kind=C_DOUBLE) * DZ_MOVINGFRAME ) - ( HALF * DZ_MOVINGFRAME )
  ZMAX_MOVINGFRAME = ( real(z_max+1,kind=C_DOUBLE) * DZ_MOVINGFRAME ) - ( HALF * DZ_MOVINGFRAME )
  j_max_1r2w = j_min_1r2w + j_width_1r2w - 1
  i_max_1r2w = i_min_1r2w + i_width_1r2w - 1

  IFPART: if ( P_TRACK_PART == YES ) then

    call DMDAGetGhostCorners ( da3r, xm3r, zm3r, ym3r, xw3r, zw3r, yw3r, ierr )
    x_min_ghost = int(xm3r)
    y_min_ghost = int(ym3r)
    z_min_ghost = int(zm3r)
    x_width_ghost = int(xw3r)
    y_width_ghost = int(yw3r)
    z_width_ghost = int(zw3r)
    x_max_ghost = x_min_ghost + x_width_ghost - 1
    y_max_ghost = y_min_ghost + y_width_ghost - 1
    z_max_ghost = z_min_ghost + z_width_ghost - 1

  end if IFPART

  IFVERBOSE: if ( V_ALLOCS == YES ) then
    write(*,*) MYRANK, k_min_3w, k_max_3w, j_min_3w, j_max_3w, i_min_3w, i_max_3w
    write(*,*) MYRANK, j_min_1r2w, j_max_1r2w, z_min, z_max, i_min_1r2w, i_max_1r2w
  end if IFVERBOSE

end subroutine f_ArraysDMDA3DLocalIndices
!----------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysDMDA2DCreate
!> Allocate fluid wave-space arrays.
!> @param ierr should return 0
subroutine f_ArraysDMDA2DCreate ( ierr )

  use g_parameters, only : V_ALLOCS, YES, MYRANK, MASTER, &
                           NODES_X, NODES_Y, NODES_Z, NHALO, &
                           F_TYPE, F_ISOTROPIC, TS_SPEC1D, TSTEPS
  use g_constants,  only : ZERO , PETSCONE, PETSCTWO
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK 

  implicit none

  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: i

  PetscInt                        :: nkk, nkj, nki, nx, ny

  PetscErrorCode        :: perr

  nkk = NODES_KK
  nkj = NODES_KJ
  nki = NODES_KI
  nx = NODES_X
  ny = NODES_Y

  call MPI_Comm_Split ( PETSC_COMM_WORLD, i_min_3w, 0, kjcomm, ierr )

  !> Create distributed array with two wave-space dimensions.
  !> Need this for computing and storing the vorticity components before the first FFT
  call DMDACreate2d ( kjcomm,                        &   ! MPI communicator
                      DM_BOUNDARY_PERIODIC,          &   ! periodic boundary conditions
                      DM_BOUNDARY_PERIODIC,          &
                      DMDA_STENCIL_BOX,              &   ! need to choose a stencil even if not used
!                      NODES_KK, NODES_KJ,            &   ! global dimension in each direction
                      nkk, nkj,            &   ! global dimension in each direction
!                      1, PETSC_DECIDE,               &   ! number of processors in each direction - let PETSc do this
                      PETSCONE, PETSC_DECIDE,               &   ! number of processors in each direction - let PETSc do this
!                      2,                             &   ! degrees of freedom = array components, two for complex number
                      PETSCTWO,                             &   ! degrees of freedom = array components, two for complex number
!                      0,                             &   ! stencil width = 0, no halo needed in wave space
                      PETSC_NULL_INTEGER,                             &   ! stencil width = 0, no halo needed in wave space
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                      da2w,                          &   ! DMDA to be created     
                      ierr )                             ! PETSc error code

  call DMSetUp ( da2w, ierr )

  call MPI_Comm_Split ( PETSC_COMM_WORLD, z_min, 0, jicomm, ierr )

  !> Create distributed array with two real-space dimensions
  !> Need this for computing the non-linear component and for the particle part
  call DMDACreate2d ( jicomm, &                          ! MPI communicator
                      DM_BOUNDARY_PERIODIC, &            ! boundary condition - we don't need it here
                      DM_BOUNDARY_PERIODIC, &
                      DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
!                      NODES_X, NODES_Y, &                ! global dimension in each direction
                      nx, ny, &                ! global dimension in each direction
!                      1, PETSC_DECIDE,  &                ! number of processors in each direction - let PETSc do this
                      PETSCONE, PETSC_DECIDE,  &                ! number of processors in each direction - let PETSc do this
!                      1,                             &   ! degrees of freedom = array components, one because real space
                      PETSCONE,                             &   ! degrees of freedom = array components, one because real space
                      PETSC_NULL_INTEGER,            &   ! stencil width = 0 
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                      da2r,                          &   ! DMDA to be created     
                      ierr )                             ! PETSc error code


  call DMSetUp ( da2r, ierr )

  !> Arrays for 1D spectra if needed
  IFISOTROPIC: if ( F_TYPE /= F_ISOTROPIC ) then

    !> Create distributed array with one wave-space and one real-space dimension
    call DMDACreate2d ( jicomm, &                          ! MPI communicator
                        DM_BOUNDARY_PERIODIC, &            ! boundary condition - we don't need it here
                        DM_BOUNDARY_PERIODIC, &
                        DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                        nki, ny, &                         ! global dimension in each direction
                        PETSCONE, PETSC_DECIDE,  &         ! number of processors in each direction - let PETSc do this
                        PETSCTWO,                      &   ! degrees of freedom = array components, one because real space
                        PETSC_NULL_INTEGER,            &   ! stencil width = 0 
                        PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                        da1r1w,                        &   ! DMDA to be created     
                        ierr )                             ! PETSc error code

    call DMSetUp ( da1r1w, ierr )

  else
    IFSPEC1D: if ( TS_SPEC1D <= TSTEPS ) then

      !> Create distributed array with one wave-space and one real-space dimension
      call DMDACreate2d ( jicomm, &                          ! MPI communicator
                          DM_BOUNDARY_PERIODIC, &            ! boundary condition - we don't need it here
                          DM_BOUNDARY_PERIODIC, &
                          DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                          nki, ny, &                         ! global dimension in each direction
                          PETSCONE, PETSC_DECIDE,  &         ! number of processors in each direction - let PETSc do this
                          PETSCTWO,                      &   ! degrees of freedom = array components, one because real space
                          PETSC_NULL_INTEGER,            &   ! stencil width = 0 
                          PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, & ! arrays containing number of nodes per cell - we don't use this
                          da1r1w,                        &   ! DMDA to be created     
                          ierr )                             ! PETSc error code

      call DMSetUp ( da1r1w, ierr )

    end if IFSPEC1D
  end if IFISOTROPIC

  IFVERBOSE: if ( V_ALLOCS == YES ) then
!> So far this is only an example, but will be used as buffer for FFTW
    write(*,*) '     | Corners: ', MYRANK, k_min_3w, j_min_3w, i_min_3w, k_width_3w, j_width_3w, i_width_3w 
    write(*,*) '     | Rank: ', MYRANK, 'Colour: ', i_min_3w

!    call PetscPrintf( PETSC_COMM_WORLD, '     | Viewing distributed real-space 2-D array: \n', ierr )
!    call DMView ( da2w, PETSC_VIEWER_STDOUT_(kjcomm), ierr )
    write(*,*) '     | Communicator: ', kjcomm

    write(*,*) '     | Corners: ', MYRANK, j_min_1r2w, z_min, i_min_1r2w, j_width_1r2w, z_width_1r2w, i_width_1r2w  
    write(*,*) '     | Rank: ', MYRANK, 'Colour: ', z_min

!    call PetscPrintf( PETSC_COMM_WORLD, '     | Viewing distributed real-space 2-D array: \n', ierr )
!    call DMView ( da2r, PETSC_VIEWER_STDOUT_(jicomm), ierr )
    write(*,*) '     | Communicator: ', jicomm
  end if IFVERBOSE


end subroutine f_ArraysDMDA2DCreate
!----------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysDMDA2DLocalIndices
!> Call PETSc to find the local start and end points of the domain
!> @param ierr should return 0
subroutine f_ArraysDMDA2DLocalIndices ( ierr )

  use g_parameters, only: MYRANK, V_ALLOCS, YES, F_TYPE, F_ISOTROPIC, &
                          TS_SPEC1D, TSTEPS
  use g_domain,     only: DX_MOVINGFRAME, XMIN_MOVINGFRAME, XMAX_MOVINGFRAME, &
                          DY_MOVINGFRAME, YMIN_MOVINGFRAME, YMAX_MOVINGFRAME
  use g_constants,  only: HALF

  implicit none

  PetscErrorCode,intent(inout)   :: ierr

  PetscInt              :: xm, ym, xw, yw
  PetscErrorCode        :: perr

!  call DMDAGetCorners ( da2r, x_min_2r, y_min_2r, PETSC_NULL_INTEGER, &
!                        x_width_2r, y_width_2r, PETSC_NULL_INTEGER, perr )
  call DMDAGetCorners ( da2r, xm, ym, PETSC_NULL_INTEGER, &
                        xw, yw, PETSC_NULL_INTEGER, ierr )
  x_min_2r   = int(xm)
  x_width_2r = int(xw)
  y_min_2r   = int(ym)
  y_width_2r = int(yw)
  x_max_2r = x_min_2r + x_width_2r - 1
  y_max_2r = y_min_2r + y_width_2r - 1

  XMIN_MOVINGFRAME = real(x_min_2r,kind=C_DOUBLE) * DX_MOVINGFRAME - ( HALF * DX_MOVINGFRAME )
  XMAX_MOVINGFRAME = real(x_max_2r+1,kind=C_DOUBLE) * DX_MOVINGFRAME - ( HALF * DX_MOVINGFRAME )
  YMIN_MOVINGFRAME = real(y_min_2r,kind=C_DOUBLE) * DY_MOVINGFRAME - ( HALF * DY_MOVINGFRAME )
  YMAX_MOVINGFRAME = real(y_max_2r+1,kind=C_DOUBLE) * DY_MOVINGFRAME - ( HALF * DY_MOVINGFRAME )

  IFISOTROPIC: if ( F_TYPE /= F_ISOTROPIC ) then

    call DMDAGetCorners ( da1r1w, xm, ym, PETSC_NULL_INTEGER, &
                          xw, yw, PETSC_NULL_INTEGER, ierr )
    i_min_1r1w   = int(xm)
    i_width_1r1w = int(xw)
    y_min_1r1w   = int(ym)
    y_width_1r1w = int(yw)
    i_max_1r1w = i_min_1r1w + i_width_1r1w - 1
    y_max_1r1w = y_min_1r1w + y_width_1r1w - 1

  else
    IFSPEC1D: if ( TS_SPEC1D <= TSTEPS ) then

      call DMDAGetCorners ( da1r1w, xm, ym, PETSC_NULL_INTEGER, &
                            xw, yw, PETSC_NULL_INTEGER, ierr )
      i_min_1r1w   = int(xm)
      i_width_1r1w = int(xw)
      y_min_1r1w   = int(ym)
      y_width_1r1w = int(yw)
      i_max_1r1w = i_min_1r1w + i_width_1r1w - 1
      y_max_1r1w = y_min_1r1w + y_width_1r1w - 1

    end if IFSPEC1D
  end if IFISOTROPIC

  IFVERBOSE: if ( V_ALLOCS == YES ) then
    write(*,*) MYRANK, x_min_2r, x_max_2r, y_min_2r, y_max_2r
  end if IFVERBOSE

end subroutine f_ArraysDMDA2DLocalIndices
!----------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_ArraysGatherLocalDimensions
!> Find out local dimensions on other processes for transpose
!> @param ierr should return 0
subroutine f_ArraysGatherLocalDimensions ( ierr )

  use g_parameters, only: MYRANK, V_ALLOCS, YES

  implicit none

  PetscErrorCode,intent(inout)   :: ierr
  integer                 :: n
  integer                 :: myrankkj, myrankji, nprocskj, nprocsji

  PetscInt              :: npj, npz, npi, npy
  PetscErrorCode        :: perr

  !> Find number of processes in relevant directions.

  call DMDAGetInfo ( da3w, PETSC_NULL_INTEGER, &                                       ! DMDA, dimension of DMDA
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &          ! global dimensions M, N, P
!                PETSC_NULL_INTEGER, nprocsj_3w, PETSC_NULL_INTEGER, &                  ! number of procs m, n, p
                PETSC_NULL_INTEGER, npj, PETSC_NULL_INTEGER, &                  ! number of procs m, n, p
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &          ! number of dof per node, stencil width, ghost nodes bx
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr )     ! ghost nodes by, bz, stencil type, error code

  nprocsj_3w = int(npj)

  call DMDAGetInfo ( da1r2w, PETSC_NULL_INTEGER, &                                     ! DMDA, dimension of DMDA
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &          ! global dimensions M, N, P
!                PETSC_NULL_INTEGER, nprocsz_1r2w, nprocsi_1r2w, &                      ! number of procs m, n, p
                PETSC_NULL_INTEGER, npz, npi, &                      ! number of procs m, n, p
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &          ! number of dof per node, stencil width, ghost nodes bx
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr )     ! ghost nodes by, bz, stencil type, error code

  nprocsz_1r2w = int(npz)
  nprocsi_1r2w = int(npi)

  call DMDAGetInfo ( da2r, PETSC_NULL_INTEGER, &                                       ! DMDA, dimension of DMDA
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &          ! global dimensions M, N, P
!                PETSC_NULL_INTEGER, nprocsy_2r, PETSC_NULL_INTEGER, &                  ! number of procs m, n, p
                PETSC_NULL_INTEGER, npy, PETSC_NULL_INTEGER, &                  ! number of procs m, n, p
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &          ! number of dof per node, stencil width, ghost nodes bx
                PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr )     ! ghost nodes by, bz, stencil type, error code

  nprocsy_2r = int(npy)

  !> Allocate arrays for dimensions on all processes.
  allocate ( j_allwidths_3w ( 0:nprocsj_3w-1 ), STAT=ierr )
  CHKERRQ ( ierr )
  allocate ( i_allwidths_1r2w ( 0:nprocsi_1r2w-1 ), STAT=ierr )
  CHKERRQ ( ierr )
  allocate ( z_allwidths_1r2w ( 0:nprocsz_1r2w-1 ), STAT=ierr )
  CHKERRQ ( ierr )
  allocate ( y_allwidths_2r ( 0:nprocsy_2r-1 ), STAT=ierr )
  CHKERRQ ( ierr )

  !> Get local ranks and sizes on kj communicator.
  call MPI_Comm_Rank ( kjcomm, myrankkj, ierr )
  call MPI_Comm_Size ( kjcomm, nprocskj, ierr )

  DOALLPROCSKJ: do n = 0, nprocskj - 1, 1   

    IFMEKJ: if ( n == myrankkj ) then
      j_allwidths_3w(n) = j_width_3w
      z_allwidths_1r2w(n) = z_width_1r2w
    end if IFMEKJ

    call MPI_Bcast ( j_allwidths_3w(n), 1, MPI_Integer, n, kjcomm, ierr )
    call MPI_Bcast ( z_allwidths_1r2w(n), 1, MPI_Integer, n, kjcomm, ierr )
 
  end do DOALLPROCSKJ

  !> Get local ranks and sizes on ji communicator.
  call MPI_Comm_Rank ( jicomm, myrankji, ierr )
  call MPI_Comm_Size ( jicomm, nprocsji, ierr )

  DOALLPROCSJI: do n = 0, nprocsji - 1, 1   

    IFMEJI: if ( n == myrankji ) then
      i_allwidths_1r2w(n) = i_width_1r2w
      y_allwidths_2r(n) = y_width_2r
    end if IFMEJI

    call MPI_Bcast ( i_allwidths_1r2w(n), 1, MPI_Integer, n, jicomm, ierr )
    call MPI_Bcast ( y_allwidths_2r(n), 1, MPI_Integer, n, jicomm, ierr )
 
  end do DOALLPROCSJI

!PetscInt,public                          :: j_maxwidth_3w, i_maxwidth_1r2w, y_maxwidth_2r
  j_maxwidth_3w = MAXVAL ( j_allwidths_3w ) 
  i_maxwidth_1r2w = MAXVAL ( i_allwidths_1r2w )
  z_maxwidth_1r2w = MAXVAL ( z_allwidths_1r2w )
  y_maxwidth_2r = MAXVAL ( y_allwidths_2r )

  nodes_kj_transpose = nprocsj_3w * j_maxwidth_3w
  nodes_ki_transpose = nprocsi_1r2w * i_maxwidth_1r2w
  nodes_z_transpose = nprocsz_1r2w * z_maxwidth_1r2w
  nodes_y_transpose = nprocsy_2r * y_maxwidth_2r

  IFVERBOSE: if ( V_ALLOCS == YES ) then
    write(*,*) 'This is process ', MYRANK, '. There are ', nprocsj_3w, ' processes in j direction.'
    write(*,*) 'This is process ', MYRANK, '. The local dimensions in j direction are ', j_allwidths_3w, &
               '(maximum: ', j_maxwidth_3w, ')'
    write(*,*) 'This is process ', MYRANK, '. There are ', nodes_kj_transpose, ' transpose nodes in j direction.'
    write(*,*) 'This is process ', MYRANK, '. There are ', nprocsi_1r2w, ' processes in i direction.'
    write(*,*) 'This is process ', MYRANK, '. The local dimensions in i direction are ', i_allwidths_1r2w, &
               '(maximum: ', i_maxwidth_1r2w, ')'
    write(*,*) 'This is process ', MYRANK, '. There are ', nodes_ki_transpose, ' transpose nodes in i direction.'
    write(*,*) 'This is process ', MYRANK, '. There are ', nprocsz_1r2w, ' processes in z direction.'
    write(*,*) 'This is process ', MYRANK, '. The local dimensions in z direction are ', z_allwidths_1r2w, &
               '(maximum: ', z_maxwidth_1r2w, ')'
    write(*,*) 'This is process ', MYRANK, '. There are ', nodes_z_transpose, ' transpose nodes in z direction.'
    write(*,*) 'This is process ', MYRANK, '. There are ', nprocsy_2r, ' processes in y direction.'
    write(*,*) 'This is process ', MYRANK, '. The local dimensions in y direction are ', y_allwidths_2r, &
               '(maximum: ', y_maxwidth_2r, ')'
    write(*,*) 'This is process ', MYRANK, '. There are ', nodes_y_transpose, ' transpose nodes in y direction.'
  end if IFVERBOSE

end subroutine f_ArraysGatherLocalDimensions
!----------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_ArraysVectorCreate
!> Allocate vectors from DMDAs.
!> @param ierr should return 0
subroutine f_ArraysVectorCreate ( ierr )

  use g_parameters, only : V_ALLOCS, V_TEST, P_TRACK_PART, P_TWO_WAY, YES, &
                           NS_DEFORMED, NS_R_CONV, &
                           F_TYPE, F_ISOTROPIC, TS_SPEC1D, TSTEPS
  use g_constants, only  : ZERO

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  PetscErrorCode        :: perr

  !> u1 to u3: components of the velocity in wave space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for velocity &
                                       in wave space. \n', ierr )
  call DMGetGlobalVector ( da3w, u1_3w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da3w, u2_3w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da3w, u3_3w, ierr )
  CHKERRQ ( ierr )

  !> u1 to u3: components of the velocity in 2d wave space 1d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for velocity &
                                       in real space (1d) / wave space (2d). \n', ierr )
  call DMGetGlobalVector ( da1r2w, u1_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u2_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u3_1r2w, ierr )
  CHKERRQ ( ierr )

  !> u4 to u6: components of the vorticity or non-linear term in 2d wave space 1d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for vorticity or non-linear &
                                       component in real space (1d) / wave space (2d). \n', ierr )
  call DMGetGlobalVector ( da1r2w, u4_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u5_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u6_1r2w, ierr )
  CHKERRQ ( ierr )

  IFCONV1R2W: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMGetGlobalVector ( da1r2w, u7_1r2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r2w, u8_1r2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r2w, u9_1r2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r2w, u10_1r2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r2w, u11_1r2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r2w, u12_1r2w, ierr )
    CHKERRQ ( ierr )
  end if IFCONV1R2W

  !> Runge-Kutta wave-space arrays 
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating Runge-Kutta arrays in wave space. \n', ierr )
  call DMGetGlobalVector ( da3w, rk3_3w_1, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da3w, rk3_3w_2, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da3w, rk3_3w_3, ierr )
  CHKERRQ ( ierr )

  !> Initialise all arrays to zero
  call VecSet ( u1_3w, ZERO, ierr )
  call VecSet ( u2_3w, ZERO, ierr )
  call VecSet ( u3_3w, ZERO, ierr )
  call VecSet ( u1_1r2w, ZERO, ierr )
  call VecSet ( u2_1r2w, ZERO, ierr )
  call VecSet ( u3_1r2w, ZERO, ierr )
  call VecSet ( u4_1r2w, ZERO, ierr )
  call VecSet ( u5_1r2w, ZERO, ierr )
  call VecSet ( u6_1r2w, ZERO, ierr )
  IFCONV1R2WZERO: if ( NS_DEFORMED == NS_R_CONV ) then
    call VecSet ( u7_1r2w, ZERO, ierr )
    call VecSet ( u8_1r2w, ZERO, ierr )
    call VecSet ( u9_1r2w, ZERO, ierr )
    call VecSet ( u10_1r2w, ZERO, ierr )
    call VecSet ( u11_1r2w, ZERO, ierr )
    call VecSet ( u12_1r2w, ZERO, ierr )
    call VecSet ( rk3_3w_1, ZERO, ierr )
    call VecSet ( rk3_3w_2, ZERO, ierr )
    call VecSet ( rk3_3w_3, ZERO, ierr )
 end if IFCONV1R2WZERO

 !> Continuing with the 2D arrays

 IFTWOWAYARRAYS: if ( P_TWO_WAY == YES ) then
    !> u1 to u3: components of the coupling term in 2d wave space
    call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for coupling term &
                                       in wave space (2d). \n', ierr )
    call DMGetGlobalVector ( da2w, u1_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u2_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u3_2w, ierr )
    CHKERRQ ( ierr )
  else
    IFTESTARRAYS: if ( V_TEST == YES ) then
      !> u1 to u3: components of the coupling term in 2d wave space
      call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for coupling term &
                                         in wave space (2d). \n', ierr )
      call DMGetGlobalVector ( da2w, u1_2w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da2w, u2_2w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da2w, u3_2w, ierr )
      CHKERRQ ( ierr )
    end if IFTESTARRAYS
  end if IFTWOWAYARRAYS


  !> u4 to u6: components of the vorticity or non-linear term in 2d wave space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for vorticity or non-linear &
                                       component in wave space (2d). \n', ierr )
  call DMGetGlobalVector ( da2w, u4_2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2w, u5_2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2w, u6_2w, ierr )
  CHKERRQ ( ierr )

  IFCONV2W: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMGetGlobalVector ( da2w, u7_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u8_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u9_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u10_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u11_2w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2w, u12_2w, ierr )
    CHKERRQ ( ierr )
  end if IFCONV2W

  !> u1 to u3: components of the velocity in 2d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for velocity &
                                       in real space (2d). \n', ierr )
  call DMGetGlobalVector ( da2r, u1_2r, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2r, u2_2r, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2r, u3_2r, ierr )
  CHKERRQ ( ierr )

  !> u4 to u6: components of the vorticity or non-linear term in 2d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for vorticity or non-linear &
                                       component in real space (2d). \n', ierr )
  call DMGetGlobalVector ( da2r, u4_2r, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2r, u5_2r, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2r, u6_2r, ierr )
  CHKERRQ ( ierr )
  IFCONV2R: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMGetGlobalVector ( da2r, u7_2r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2r, u8_2r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2r, u9_2r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2r, u10_2r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2r, u11_2r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da2r, u12_2r, ierr )
    CHKERRQ ( ierr )
  end if IFCONV2R

  !> Arrays for 1D spectra if needed
  IFISOTROPIC: if ( F_TYPE /= F_ISOTROPIC ) then

    call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for 1-D &
                                         spectra in streamwise direction. \n', ierr )
    call DMGetGlobalVector ( da1r1w, u1_1r1w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r1w, u2_1r1w, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da1r1w, u3_1r1w, ierr )
    CHKERRQ ( ierr )

  else
    IFSPEC1D: if ( TS_SPEC1D <= TSTEPS ) then

      call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for 1-D &
                                           spectra in streamwise direction. \n', ierr )
      call DMGetGlobalVector ( da1r1w, u1_1r1w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da1r1w, u2_1r1w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da1r1w, u3_1r1w, ierr )
      CHKERRQ ( ierr )

    end if IFSPEC1D
  end if IFISOTROPIC

  !> Initialise all arrays to zero
  IFTWOWAYINIT: if ( P_TWO_WAY == YES ) then
    call VecSet ( u1_2w, ZERO, ierr )
    call VecSet ( u2_2w, ZERO, ierr )
    call VecSet ( u3_2w, ZERO, ierr )
  end if IFTWOWAYINIT

  call VecSet ( u4_2w, ZERO, ierr )
  call VecSet ( u5_2w, ZERO, ierr )
  call VecSet ( u6_2w, ZERO, ierr )
  call VecSet ( u1_2r, ZERO, ierr )
  call VecSet ( u2_2r, ZERO, ierr )
  call VecSet ( u3_2r, ZERO, ierr )
  call VecSet ( u4_2r, ZERO, ierr )
  call VecSet ( u5_2r, ZERO, ierr )
  call VecSet ( u6_2r, ZERO, ierr )

  IFCONV02D: if ( NS_DEFORMED == NS_R_CONV ) then
    call VecSet ( u7_2w, ZERO, ierr )
    call VecSet ( u8_2w, ZERO, ierr )
    call VecSet ( u9_2w, ZERO, ierr )
    call VecSet ( u10_2w, ZERO, ierr )
    call VecSet ( u11_2w, ZERO, ierr )
    call VecSet ( u12_2w, ZERO, ierr )
    call VecSet ( u7_2r, ZERO, ierr )
    call VecSet ( u8_2r, ZERO, ierr )
    call VecSet ( u9_2r, ZERO, ierr )
    call VecSet ( u10_2r, ZERO, ierr )
    call VecSet ( u11_2r, ZERO, ierr )
    call VecSet ( u12_2r, ZERO, ierr )
  end if IFCONV02D

  !> Arrays for 1D spectra if needed
  IFISOTROPIC0: if ( F_TYPE /= F_ISOTROPIC ) then

    call VecSet ( u1_1r1w, ZERO, ierr )
    call VecSet ( u2_1r1w, ZERO, ierr )
    call VecSet ( u3_1r1w, ZERO, ierr )

  else
    IFSPEC1D0: if ( TS_SPEC1D <= TSTEPS ) then

      call VecSet ( u1_1r1w, ZERO, ierr )
      call VecSet ( u2_1r1w, ZERO, ierr )
      call VecSet ( u3_1r1w, ZERO, ierr )

    end if IFSPEC1D0
  end if IFISOTROPIC0

  IFPART: if ( P_TRACK_PART == YES ) then

    !> u1 to u3: components of the velocity in real space
    call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for velocity &
                                         in real space. \n', ierr )
    call DMGetGlobalVector ( da3r, u1_3r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da3r, u2_3r, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da3r, u3_3r, ierr )
    CHKERRQ ( ierr )

    call VecSet ( u1_3r, ZERO, ierr )
    CHKERRQ ( ierr )
    call VecSet ( u2_3r, ZERO, ierr )
    CHKERRQ ( ierr )
    call VecSet ( u3_3r, ZERO, ierr )
    CHKERRQ ( ierr )

    !> Get local vector with halos
!    call  DMGetLocalVector( da3r, u1_3r_local, ierr )
!    CHKERRQ(ierr)
!    call  DMGetLocalVector( da3r, u2_3r_local, ierr )
!    CHKERRQ(ierr)
!    call  DMGetLocalVector( da3r, u3_3r_local, ierr )
!    CHKERRQ(ierr)

!    call  DMCreateLocalVector( da3r, u1_3r_local, ierr )
!    CHKERRQ(ierr)
!    call  DMCreateLocalVector( da3r, u2_3r_local, ierr )
!    CHKERRQ(ierr)
!    call  DMCreateLocalVector( da3r, u3_3r_local, ierr )
!    CHKERRQ(ierr)

!    call VecSet ( u1_3r_local, ZERO, ierr )
!    CHKERRQ ( ierr )
!    call VecSet ( u2_3r_local, ZERO, ierr )
!    CHKERRQ ( ierr )
!    call VecSet ( u3_3r_local, ZERO, ierr )
!    CHKERRQ ( ierr )

  end if IFPART

end subroutine f_ArraysVectorCreate
!----------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_ArraysVectorCreateRestart
!> Allocate vectors from DMDAs for restart file only.
!> @param ierr should return 0
subroutine f_ArraysVectorCreateRestart ( ierr )

  use g_parameters, only : V_ALLOCS, YES, MYRANK, MASTER, &
                           NODES_X, NODES_Y, NODES_Z, NHALO, &
                           P_TRACK_PART, YES, &
                           F_RESTART_SCALE
  use g_constants,  only : ZERO , PETSCONE, PETSCTWO, PETSCTHREE 
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, KJ_MAX

  implicit none

  PetscErrorCode,intent(inout)             :: ierr

  integer                           :: colourj, colouri
  integer                           :: wvncount, proccount, proccountrs

  PetscInt                          :: nkk, nkj, nki, nkpos, nkneg, &
                                       njpos, njneg 
  PetscInt                          :: npi3w, npj3w, npk3w
  PetscInt                          :: npi3wrs, npj3wrs, npk3wrs
  PetscInt                          :: nh
  PetscInt,dimension(:),allocatable :: lk3w, lj3w, li3w
  PetscInt,dimension(:),allocatable :: lk3wrs, lj3wrs, li3wrs

  PetscErrorCode        :: perr

  !> Total number of wave modes in the restart file in each direction.
  nkk = 2 * ( ( NODES_Z / F_RESTART_SCALE ) / 3 ) + 1
  nkj = 2 * ( ( NODES_Y / F_RESTART_SCALE ) / 3 ) + 1
  nki =     ( ( NODES_X / F_RESTART_SCALE ) / 3 ) + 1

  !> Number of positive and negative wave modes in k and j directions.
  !! Do not need this for i because of the conjugate symmetry.
  nkpos = ( ( NODES_Z / F_RESTART_SCALE ) / 3 ) + 1
  nkneg = ( ( NODES_Z / F_RESTART_SCALE ) / 3 ) 
  njpos = ( ( NODES_Y / F_RESTART_SCALE ) / 3 ) + 1
  njneg = ( ( NODES_Y / F_RESTART_SCALE ) / 3 )


  !> Get all the necessary info from the main wavespace DMDA 
  call DMDAGetInfo ( da3w, PETSC_NULL_INTEGER, &                                        ! DMDA, dimension of DMDA
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! global dimensions M, N, P
                     npk3w, npj3w, npi3w, &                                             ! number of procs m, n, p
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! number of dof per node, stencil width, ghost nodes bx
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr ) ! ghost nodes by, bz, stencil type, error code

  !> Get vectors with the number of nodes per process in each direction 
  allocate ( lk3w(npk3w), lj3w(npj3w), li3w(npi3w), stat=ierr )
  call DMDAGetOwnershipRanges ( da3w, lk3w, lj3w, li3w, ierr )

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,*) ' Array of k dimensions: ', lk3w
    write(*,*) ' Array of j dimensions: ', lj3w
    write(*,*) ' Array of i dimensions: ', li3w
    write(*,*) ' Dimensions of restart array: ', nkk, nkj, nki
    write(*,*) ' In k direction there are ', nkpos, ' positive wavemodes and ', nkneg, ' negative wavemodes. '
    write(*,*) ' In j direction there are ', njpos, ' positive wavemodes and ', njneg, ' negative wavemodes. '
  end if IFMASTER

  !> Find number of processes in each direction. Obvious only for k.
  npk3wrs = 1

  !> Initialise restart colours
  restartcolour = 0

  !> Determine j colour of each process for creating a communicator.
  !! If there is only one process in j direction, it will certainly contain
  !! the entire wavenumber range of the restart file.
  IF1J: if ( npj3w == 1 ) then

    colourj = 1
    npj3wrs = npj3w

  else

    colourj = 0

    IFJPOS: if ( j_min_3w <= KJ_MAX ) then
      !> Is the lowest j on this process smaller than the highest positive j
      !! of the restart file?
      IFJPOSWITHIN: if ( j_min_3w <= ( njpos - 1) ) then
!        write(*,*) ' j range on rank ', MYRANK, ': ', j_min_3w, j_max_3w
        colourj = 1
      end if IFJPOSWITHIN
    else 
      !> Is the lowest negative j mode on this process smaller than the highest
      !! negative j of the restart file?
      IFJNEGWITHIN: if ( ( j_max_3w - NODES_KJ ) >= - ( njneg - 1) ) then
!        write(*,*) ' j range on rank ', MYRANK, ': ', j_min_3w - NODES_KJ, j_max_3w - NODES_KJ
        colourj = 1
      end if IFJNEGWITHIN
    end if IFJPOS

    !> Sum up number of nodes per process and count how many processes in j
    !! direction contain the positive j range of the restart file.
    npj3wrs = 0
    wvncount  = 0
    proccount = 0
    DOJPOS: do while ( wvncount < njpos )
      proccount = proccount + 1
      wvncount = wvncount + lj3w(proccount)
      npj3wrs = npj3wrs + 1
!      write(*,*) 'j', proccount, wvncount, npj3wrs
    end do DOJPOS

    !> Same for the negative wave numbers.
    wvncount  = 0
    proccount = npj3w + 1
    DOJNEG: do while ( wvncount < njneg )
      proccount = proccount - 1
      wvncount = wvncount + lj3w(proccount)
      npj3wrs = npj3wrs + 1
!      write(*,*) 'j', proccount, wvncount, npj3wrs
    end do DOJNEG

  end if IF1J

  !> If there is only one process in i direction, it will certainly contain
  !! the entire wavenumber range of the restart file.
  IF1I: if ( npi3w == 1 ) then
    
    colouri = 1
    npi3wrs = npi3w

  else

    !> Determine i colour of each process for creating a communicator.
    colouri = 0
    IFI: if ( i_min_3w <= KI_MAX ) then
      !> Is the lowest i on this process smaller than the highest i 
      !! of the restart file?
      IFIWITHIN: if ( i_min_3w <= ( nki - 1) ) then
!        write(*,*) ' i range on rank ', MYRANK, ': ', i_min_3w, i_max_3w
        colouri = 1
      end if IFIWITHIN
    end if IFI

    !> Sum up number of nodes per process and count how many processes in i
    !! direction contain the i range of the restart file.
    npi3wrs = 0
    wvncount  = 0
    proccount = 0
    DOI: do while ( wvncount < nki )
      proccount = proccount + 1
      wvncount = wvncount + li3w(proccount)
      npi3wrs = npi3wrs + 1
!      write(*,*) 'i', proccount, wvncount, npi3wrs
    end do DOI

  end if IF1I

  !> Only include processes that contain the correct wavenumber range. Always true in k direction.
  !! Processes must own both the relevant i and j ranges.
  restartcolour = colourj * colouri
!  write(*,*) ' colour on rank ', MYRANK, ': ', restartcolour

  !> Create communicator for reading the restartfile.
  call MPI_Comm_Split ( PETSC_COMM_WORLD, restartcolour, 0, restartcomm, ierr )


!  write(*,*) 'Rank: ', MYRANK, ' Processes: ', npk3wrs, npj3wrs, npi3wrs

  !> Get vectors with the number of nodes per process in each direction 
  allocate ( lk3wrs(npk3wrs), lj3wrs(npj3wrs), li3wrs(npi3wrs), stat=ierr )

  !> Only one process in k direction, so it contains all k nodes.
  lk3wrs = [ nkk ]

  wvncount  = 0
  proccount = 0
  proccountrs = 0

  IF1JCOUNT: if ( npj3wrs == 1 ) then

    lj3wrs(1) = nkj

  else

    !> Go through processes from 0 to maximum wavenumber in j direction and
    !! sum up wave modes. If total lower than maximum take all, otherwise
    !! only take the remaining bit.
    DOJPOS2: do while ( wvncount < njpos )
      proccount = proccount + 1
      proccountrs = proccountrs + 1
      wvncount = wvncount + lj3w(proccount)
      IFJPOSLIMIT: if ( wvncount <= njpos ) then
        lj3wrs(proccountrs) = lj3w(proccount)
      else
        lj3wrs(proccountrs) = lj3w(proccount) - (wvncount - njpos)
      end if IFJPOSLIMIT
    end do DOJPOS2

    wvncount  = 0
    proccount = npj3w + 1
    proccountrs = npj3wrs + 1
    DOJNEG2: do while ( wvncount < njneg )
      proccount = proccount - 1
      proccountrs = proccountrs - 1
      wvncount = wvncount + lj3w(proccount)
      IFJNEGLIMIT: if ( wvncount <= njneg ) then
        lj3wrs(proccountrs) = lj3w(proccount)
      else
        lj3wrs(proccountrs) = lj3w(proccount) - (wvncount - njneg)
      end if IFJNEGLIMIT
    end do DOJNEG2

  end if IF1JCOUNT

!  li3wrs = [ 11 ]
  wvncount  = 0
  proccount = 0
  proccountrs = 0
  DOI2: do while ( wvncount < nki )
    proccount = proccount + 1
    proccountrs = proccountrs + 1
    wvncount = wvncount + li3w(proccount)
    IFILIMIT: if ( wvncount <= nki ) then
      li3wrs(proccountrs) = li3w(proccount)
    else
      li3wrs(proccountrs) = li3w(proccount) - (wvncount - nki)
    end if IFILIMIT
  end do DOI2

!  write(*,*) ' Array of k dimensions: ', lk3wrs
!  write(*,*) ' Array of j dimensions: ', lj3wrs
!  write(*,*) ' Array of i dimensions: ', li3wrs

!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  !> Create distributed array with three wave-space dimensions
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Creating distributed array for reading in &
                                            the scaled restart file. \n', ierr )

  IFRESTARTCOMM: if ( restartcolour == 1 ) then
    call DMDACreate3d ( restartcomm, &                     ! MPI communicator
                        DM_BOUNDARY_PERIODIC, &            ! periodic boundary conditions 
                        DM_BOUNDARY_PERIODIC, &
                        DM_BOUNDARY_PERIODIC, &
                        DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                        nkk, nkj, nki,  &   ! global dimension in each direction
                        npk3wrs, npj3wrs, npi3wrs, &   ! number of processors in each direction
                        PETSCTWO,                             &   ! degrees of freedom = array components, we use two for real and imaginary components
                        PETSC_NULL_INTEGER,                             &   ! stencil width = 0
                        lk3wrs, lj3wrs, li3wrs, & ! arrays containing number of nodes per cell - we don't use this
                        da3w_rs,                          &   ! DMDA to be created     
                        ierr )                             ! PETSc error code

    call DMSetUp ( da3w_rs, ierr )
  
    !> u1 to u3: components of the velocity in wave space
    call PetscPrintf ( restartcomm, '    | Allocating arrays for reading in &
                                         the scaled restart file. \n', ierr )
    call DMGetGlobalVector ( da3w_rs, u1_3w_rs, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da3w_rs, u2_3w_rs, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da3w_rs, u3_3w_rs, ierr )
    CHKERRQ ( ierr )


  end if IFRESTARTCOMM

  deallocate ( lk3w, lj3w, li3w, stat=ierr )
  deallocate ( lk3wrs, lj3wrs, li3wrs, stat=ierr )

end subroutine f_ArraysVectorCreateRestart
!----------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_ArraysVectorCreateRestartDeformed
!> Allocate vectors from DMDAs for restart file only.
!> @param ierr should return 0
subroutine f_ArraysVectorCreateRestartDeformed ( ierr )

  use g_parameters, only : V_ALLOCS, YES, MYRANK, MASTER, &
                           NODES_X, NODES_Y, NODES_Z, NHALO, &
                           P_TRACK_PART, YES, &
                           F_RESTART_SCALE, &
                           B110, B220, B330
  use g_constants,  only : ZERO , PETSCONE, PETSCTWO, PETSCTHREE 
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, KJ_MAX

  implicit none

  PetscErrorCode,intent(inout)             :: ierr

  integer                           :: colourj, colouri
  integer                           :: wvncount, wvncountscaled, proccount, proccountrs
  integer                           :: countthisproc

  integer                           :: i, j

  PetscInt                          :: nx, ny, nz, nodesrestart
  PetscInt                          :: nkk, nkj, nki, nkpos, nkneg, &
                                       njpos, njneg 
  PetscInt                          :: npi3w, npj3w, npk3w
  PetscInt                          :: npi3wrs, npj3wrs, npk3wrs
  PetscInt                          :: nh
  PetscInt,dimension(:),allocatable :: lk3w, lj3w, li3w
  PetscInt,dimension(:),allocatable :: lk3wrs, lj3wrs, li3wrs

  PetscErrorCode        :: perr

  !> Determine size of the restart file. We assume it is cubic with
  !! side length 2pi in all directions.
  nx = NODES_X / B110 / F_RESTART_SCALE
  ny = NODES_Y / B220 / F_RESTART_SCALE
  nz = NODES_Z / B330 / F_RESTART_SCALE
  nodesrestart = min ( nx, ny, nz ) 

  !> Total number of wave modes in the restart file in each direction.
  nkk = 2 * ( nodesrestart / 3 ) + 1
  nkj = 2 * ( nodesrestart / 3 ) + 1
  nki =     ( nodesrestart / 3 ) + 1

  !> Number of positive and negative wave modes in k and j directions.
  !! Do not need this for i because of the conjugate symmetry.
  nkpos = ( nodesrestart / 3 ) + 1
  nkneg = ( nodesrestart / 3 ) 
  njpos = ( nodesrestart / 3 ) + 1
  njneg = ( nodesrestart / 3 )


  !> Get all the necessary info from the main wavespace DMDA 
  call DMDAGetInfo ( da3w, PETSC_NULL_INTEGER, &                                        ! DMDA, dimension of DMDA
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! global dimensions M, N, P
                     npk3w, npj3w, npi3w, &                                             ! number of procs m, n, p
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! number of dof per node, stencil width, ghost nodes bx
                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr ) ! ghost nodes by, bz, stencil type, error code

  allocate ( lk3w(npk3w), lj3w(npj3w), li3w(npi3w), stat=ierr )

  call DMDAGetOwnershipRanges ( da3w, lk3w, lj3w, li3w, ierr )

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,*) ' Array of k dimensions: ', lk3w
    write(*,*) ' Array of j dimensions: ', lj3w
    write(*,*) ' Array of i dimensions: ', li3w
    write(*,*) ' Dimensions of restart array: ', nkk, nkj, nki
    write(*,*) ' In k direction there are ', nkpos, ' positive wavemodes and ', nkneg, ' negative wavemodes. '
    write(*,*) ' In j direction there are ', njpos, ' positive wavemodes and ', njneg, ' negative wavemodes. '
  end if IFMASTER

  !> Find number of processes in each direction. Obvious only for k.
  npk3wrs = 1

  !> Initialise restart colours
  restartcolour = 0

  !> Determine j colour of each process for creating a communicator.
  !! If there is only one process in j direction, it will certainly contain
  !! the entire wavenumber range of the restart file.
  IF1J: if ( npj3w == 1 ) then

    colourj = 1
    npj3wrs = npj3w

  else

    colourj = 0

    IFJPOS: if ( j_min_3w <= KJ_MAX ) then
      !> Is the lowest j on this process smaller than the highest positive j
      !! of the restart file?
      !! We need to scale j by B220 to get the index in the deformed system.
      IFJPOSWITHIN: if ( j_min_3w <= ( B220 * ( njpos - 1) ) ) then
!        write(*,*) ' j range on rank ', MYRANK, ': ', j_min_3w, j_max_3w
        !> Only every B220'th mode is counted, so make sure there actually is
        !! one on this process.
        IFJPOSWIDTH: if ( j_width_3w >= B220 ) then
          colourj = 1
        else
          DOJPOSFIND: do j = j_min_3w, j_max_3w, 1
            IFJPOSMOD: if ( modulo ( j, B220 ) == 0 ) then
              colourj = 1
            end if IFJPOSMOD
          end do DOJPOSFIND
        end if IFJPOSWIDTH
      end if IFJPOSWITHIN
    else 
      IFJNEGWITHIN: if ( ( j_max_3w - NODES_KJ ) >= - ( B220 * ( njneg - 1) ) ) then
!        write(*,*) ' j range on rank ', MYRANK, ': ', j_min_3w - NODES_KJ, j_max_3w - NODES_KJ
        !> Only every B220'th mode is counted, so make sure there actually is
        !! one on this process.
        IFJNEGWIDTH: if ( j_width_3w >= B220 ) then
          colourj = 1
        else
          DOJNEGFIND: do j = j_min_3w, j_max_3w, 1
            IFJNEGMOD: if ( modulo ( (j - NODES_KJ), B220 ) == 0 ) then
              colourj = 1
            end if IFJNEGMOD
          end do DOJNEGFIND
        end if IFJNEGWIDTH
      end if IFJNEGWITHIN
    end if IFJPOS

    !> This is the global part. Find out how wave modes in j direction are 
    !! distributed over the entire domain.
    npj3wrs = 0
    !> Start with (positive) mode 0.
    wvncount  = 0
    wvncountscaled  = 0
    proccount = 0
    !> Only stop when all wave modes from the restart file are counted.
    DOJPOS: do while ( wvncountscaled < njpos )
      countthisproc = 0
      DOJPOSCOUNT: do j = wvncount, ( wvncount + lj3w(proccount) ), 1 
!      write(*,*) 'j', proccount, wvncount, npj3wrs
        IFJPOSCOUNTED: if ( modulo ( j, B220 ) == 0 ) then
          wvncountscaled = wvncountscaled + 1
          countthisproc = 1
        end if IFJPOSCOUNTED
      end do DOJPOSCOUNT
      wvncount = wvncount + lj3w(proccount) + 1
      proccount = proccount + 1
      npj3wrs = npj3wrs + countthisproc
    end do DOJPOS

    !> Start with negative mode 1.
    wvncount  = 1
    proccount = npj3w + 1
    wvncountscaled  = 0
    !> Only stop when all wave modes from the restart file are counted.
    DOJNEG: do while ( wvncountscaled < njneg )
      countthisproc = 0
      DOJNEGCOUNT: do j = wvncount, ( wvncount + lj3w(proccount) ), 1 
!      write(*,*) 'j', proccount, wvncount, npj3wrs
        IFJNEGCOUNTED: if ( modulo ( j, B220 ) == 0 ) then
          wvncountscaled = wvncountscaled + 1
          countthisproc = 1
        end if IFJNEGCOUNTED
      end do DOJNEGCOUNT
      wvncount = wvncount + lj3w(proccount) + 1
      proccount = proccount - 1
      npj3wrs = npj3wrs + countthisproc
!      write(*,*) 'j', proccount, wvncount, npj3wrs
    end do DOJNEG

  end if IF1J


  IF1I: if ( npi3w == 1 ) then
    
    colouri = 1
    npi3wrs = npi3w

  else

    !> Determine i colour of each process for creating a communicator.
    colouri = 0
    IFI: if ( i_min_3w <= KI_MAX ) then
      IFIWITHIN: if ( i_min_3w <= ( B110 * ( nki - 1) ) ) then
!        write(*,*) ' i range on rank ', MYRANK, ': ', i_min_3w, i_max_3w
        IFIWIDTH: if ( i_width_3w >= B110 ) then
          colouri = 1
        else
          DOIFIND: do i = i_min_3w, i_max_3w, 1
            IFIMOD: if ( modulo ( j, B110 ) == 0 ) then
              colouri = 1
            end if IFIMOD
          end do DOIFIND
        end if IFIWIDTH
      end if IFIWITHIN
    end if IFI

    !> This is the global part. Find out how wave modes in i direction are 
    !! distributed over the entire domain.
    npi3wrs = 0
    !> Start with mode 0.
    wvncount  = 0
    wvncountscaled  = 0
    proccount = 0
    !> Only stop when all wave modes from the restart file are counted.
    DOI: do while ( wvncountscaled < nki )
      countthisproc = 0
      DOICOUNT: do i = wvncount, ( wvncount + li3w(proccount) ), 1 
!      write(*,*) 'j', proccount, wvncount, npj3wrs
        IFICOUNTED: if ( modulo ( i, B110 ) == 0 ) then
          wvncountscaled = wvncountscaled + 1
          countthisproc = 1
        end if IFICOUNTED
      end do DOICOUNT
      wvncount = wvncount + li3w(proccount) + 1
      proccount = proccount + 1
      npi3wrs = npi3wrs + countthisproc
!      write(*,*) 'i', proccount, wvncount, npi3wrs
    end do DOI

  end if IF1I

  !> Only include processes that contain the correct wavenumber range. (Always true in k direction.)
  restartcolour = colourj * colouri
!  write(*,*) ' colour on rank ', MYRANK, ': ', restartcolour

  !> Create communicator for reading the restartfile.
  call MPI_Comm_Split ( PETSC_COMM_WORLD, restartcolour, 0, restartcomm, ierr )


!  write(*,*) 'Rank: ', MYRANK, ' Processes: ', npk3wrs, npj3wrs, npi3wrs

  allocate ( lk3wrs(npk3wrs), lj3wrs(npj3wrs), li3wrs(npi3wrs), stat=ierr )

  !> Now determine correct distribution arrays.
  lk3wrs = [ nkk ]
!  lk3wrs = [ 21 ]

!  lj3wrs = [ 11, 10 ]
  IF1JCOUNT: if ( npj3wrs == 1 ) then

    lj3wrs(1) = nkj

  else

    !> Start with (positive) mode 0.
    wvncount  = 0
    wvncountscaled  = 0
    proccount = 0
    proccountrs = 0
    !> Only stop when all wave modes from the restart file are counted.
    DOJPOS2: do while ( wvncountscaled < njpos )
      countthisproc = 0
      lj3wrs(proccountrs) = 0 
      DOJPOS2COUNT: do j = wvncount, ( wvncount + lj3w(proccount) ), 1 
        IFJPOS2COUNTED: if ( modulo ( j, B220 ) == 0 ) then
          IFJPOS2LIMIT: if ( ( j / B220 ) <= njpos ) then
            wvncountscaled = wvncountscaled + 1
            lj3wrs(proccountrs) = lj3w(proccountrs) + 1
            countthisproc = 1
          end if IFJPOS2LIMIT
        end if IFJPOS2COUNTED
      end do DOJPOS2COUNT
      wvncount = wvncount + lj3w(proccount) + 1
      proccount = proccount + 1
      proccountrs = proccountrs + countthisproc
    end do DOJPOS2

    !> Start with negative mode 1.
    wvncount  = 1
    proccount = npj3w + 1
    proccountrs = npj3wrs + 1
    wvncountscaled  = 0
    !> Only stop when all wave modes from the restart file are counted.
    DOJNEG2: do while ( wvncountscaled < njneg )
      countthisproc = 0
      lj3wrs(proccountrs) = 0 
      DOJNEG2COUNT: do j = wvncount, ( wvncount + lj3w(proccount) ), 1 
        IFJNEG2COUNTED: if ( modulo ( j, B220 ) == 0 ) then
          IFJNEG2LIMIT: if ( ( j / B220 ) <= njneg ) then
            wvncountscaled = wvncountscaled + 1
            lj3wrs(proccountrs) = lj3w(proccountrs) + 1
            countthisproc = 1
          end if IFJNEG2LIMIT
        end if IFJNEG2COUNTED
      end do DOJNEG2COUNT
      wvncount = wvncount + lj3w(proccount) + 1
      proccount = proccount - 1
      proccountrs = proccountrs - countthisproc
    end do DOJNEG2

  end if IF1JCOUNT

!  li3wrs = [ 11 ]
  wvncount  = 0
  proccount = 0
  proccountrs = 0
  !> Only stop when all wave modes from the restart file are counted.
  DOI2: do while ( wvncountscaled < nki )
    proccount = proccount + 1
    proccountrs = proccountrs + 1
    wvncount = wvncount + li3w(proccount)
    IFILIMIT: if ( wvncount <= nki ) then
      li3wrs(proccountrs) = li3w(proccount)
    else
      li3wrs(proccountrs) = li3w(proccount) - (wvncount - nki)
    end if IFILIMIT
  end do DOI2

!  write(*,*) ' Array of k dimensions: ', lk3wrs
!  write(*,*) ' Array of j dimensions: ', lj3wrs
!  write(*,*) ' Array of i dimensions: ', li3wrs

!  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  !> Create distributed array with three wave-space dimensions
  call PetscPrintf ( PETSC_COMM_WORLD, '    | Creating distributed array for reading in &
                                            the scaled restart file. \n', ierr )

  IFRESTARTCOMM: if ( restartcolour == 1 ) then
    call DMDACreate3d ( restartcomm, &                     ! MPI communicator
                        DM_BOUNDARY_PERIODIC, &            ! periodic boundary conditions 
                        DM_BOUNDARY_PERIODIC, &
                        DM_BOUNDARY_PERIODIC, &
                        DMDA_STENCIL_BOX, &                ! need to choose a stencil even if not used
                        nkk, nkj, nki,  &   ! global dimension in each direction
                        npk3wrs, npj3wrs, npi3wrs, &   ! number of processors in each direction
                        PETSCTWO,                             &   ! degrees of freedom = array components, we use two for real and imaginary components
                        PETSC_NULL_INTEGER,                             &   ! stencil width = 0
                        lk3wrs, lj3wrs, li3wrs, & ! arrays containing number of nodes per cell - we don't use this
                        da3w_rs,                          &   ! DMDA to be created     
                        ierr )                             ! PETSc error code

    call DMSetUp ( da3w_rs, ierr )
  
    !> u1 to u3: components of the velocity in wave space
    call PetscPrintf ( restartcomm, '    | Allocating arrays for reading in &
                                         the scaled restart file. \n', ierr )
    call DMGetGlobalVector ( da3w_rs, u1_3w_rs, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da3w_rs, u2_3w_rs, ierr )
    CHKERRQ ( ierr )
    call DMGetGlobalVector ( da3w_rs, u3_3w_rs, ierr )
    CHKERRQ ( ierr )


  end if IFRESTARTCOMM

  deallocate ( lk3w, lj3w, li3w, stat=ierr )
  deallocate ( lk3wrs, lj3wrs, li3wrs, stat=ierr )

end subroutine f_ArraysVectorCreateRestartDeformed
!----------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_ArraysDMDADealloc
!> Deallocate fluid wave-space arrays.
!> @param ierr should return 0
subroutine f_ArraysDMDADealloc ( ierr )

  use g_parameters, only : V_ALLOCS, YES, MYRANK, MASTER
!> ZERO = 0.0 
  use g_constants,  only : ZERO

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  PetscErrorCode        :: perr


  deallocate ( j_allwidths_3w, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( i_allwidths_1r2w, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( z_allwidths_1r2w, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( y_allwidths_2r, STAT=ierr )
  CHKERRQ ( ierr )

end subroutine f_ArraysDMDADealloc
!----------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysRemeshBufferAllocate
!> Allocate buffer for remeshing.
!> @param ierr should return 0
subroutine f_ArraysRemeshBufferAllocate ( ierr )

  use g_parameters, only : V_ALLOCS, V_TEST, P_TWO_WAY, YES, NODES_Z
  use g_constants, only  : ZERO

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

!  allocate(arr_remesh(1:2,1:NODES_Z,i_min_3w:i_max_3w),STAT=ierr)
!  CHKERRQ(ierr)
  allocate(arr_u1_remesh(0:1,k_min_3w:k_max_3w),STAT=ierr)
  CHKERRQ(ierr)
  allocate(arr_u2_remesh(0:1,k_min_3w:k_max_3w),STAT=ierr)
  CHKERRQ(ierr)
  allocate(arr_u3_remesh(0:1,k_min_3w:k_max_3w),STAT=ierr)
  CHKERRQ(ierr)

end subroutine f_ArraysRemeshBufferAllocate
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysRemeshBufferDeallocate
!> Deallocate buffer for remeshing.
!> @param ierr should return 0
subroutine f_ArraysRemeshBufferDeallocate ( ierr )

  use g_parameters, only : V_ALLOCS, V_TEST, P_TWO_WAY, YES
  use g_constants, only  : ZERO

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

!  deallocate(arr_remesh,STAT=ierr)
!  CHKERRQ(ierr)
  deallocate(arr_u1_remesh,STAT=ierr)
  CHKERRQ(ierr)
  deallocate(arr_u2_remesh,STAT=ierr)
  CHKERRQ(ierr)
  deallocate(arr_u3_remesh,STAT=ierr)
  CHKERRQ(ierr)

end subroutine f_ArraysRemeshBufferDeallocate
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! function f_ArraysWavenumberPositive 
!> Compute a component of the wave vector for positive wavenumbers
!> @param k array index for the wave vector component
!> @param K0 minimum wavenumber for the wave vector component
pure function f_ArraysWavenumberPositive ( k, K0 ) result ( wvn )
  integer, intent(in)             :: k 
  real(kind=C_DOUBLE), intent(in) :: K0
  real(kind=C_DOUBLE)             :: wvn

  wvn = K0 * real ( k, kind=C_DOUBLE)

end function f_ArraysWavenumberPositive
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! function f_ArraysWavenumberNegative 
!> Compute a component of the wave vector for negative wavenumbers
!> @param k array index for the wave vector component
!> @param K0 minimum wavenumber for the wave vector component
pure function f_ArraysWavenumberNegative ( k, KMAX, K0 ) result ( wvn )
  integer, intent(in)             :: k, KMAX 
  real(kind=C_DOUBLE), intent(in) :: K0
  real(kind=C_DOUBLE)             :: wvn

  integer                         :: shiftedindex, nodes

  nodes = 2 * KMAX + 1
  shiftedindex = k - nodes
  wvn = K0 * real ( shiftedindex, kind=C_DOUBLE)

end function f_ArraysWavenumberNegative
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! function f_ArraysWavenumber 
!> Compute a component of the wave vector. Use f_ArraysWavenumberPositive if the 
!! component is positive and f_ArraysWavenumberNegative if the component is 
!! negative
!> @param k array index for the wave vector component
!> @param K0 minimum wavenumber for the wave vector component
!> @param KMAX maximum wavenumber for the wave vector component, needed to determine
!! if the wavenumber is positive ( <= KMAX ) or negative ( > KMAX ).
pure function f_ArraysWavenumber ( k, KMAX, K0 ) result ( wvn )
  integer, intent(in)              :: k, KMAX 
  real(kind=C_DOUBLE), intent(in)  :: K0
  real(kind=C_DOUBLE)              :: wvn

  IFPOSNEG: if ( k <= KMAX ) then 
    wvn = f_ArraysWavenumberPositive ( k, K0 )
  else
    wvn = f_ArraysWavenumberNegative ( k, KMAX, K0 )
  end if IFPOSNEG

end function f_ArraysWavenumber
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! function f_ArraysWavenumberPositiveLong
!> Same as f_ArraysWavenumberPositive, but in extended precision 
!> @param k array index for the wave vector component
!> @param K0 minimum wavenumber for the wave vector component
pure function f_ArraysWavenumberPositiveLong ( k, K0 ) result ( wvn )
  integer, intent(in)             :: k 
  real(kind=C_DOUBLE), intent(in) :: K0
  real(kind=C_LONG_DOUBLE)        :: wvn

  wvn = real( K0, kind=C_LONG_DOUBLE ) * real ( k, kind=C_LONG_DOUBLE )

end function f_ArraysWavenumberPositiveLong
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! function f_ArraysWavenumberNegativeLong
!> Same as f_ArraysWavenumberNegative, but in extended precision 
!> @param k array index for the wave vector component
!> @param K0 minimum wavenumber for the wave vector component
pure function f_ArraysWavenumberNegativeLong ( k, KMAX, K0 ) result ( wvn )
  integer, intent(in)             :: k, KMAX 
  real(kind=C_DOUBLE), intent(in) :: K0
  real(kind=C_LONG_DOUBLE)        :: wvn

  integer                         :: shiftedindex, nodes

  nodes = 2 * KMAX + 1
  shiftedindex = k - nodes
  wvn = real ( K0 , kind=C_LONG_DOUBLE ) * real ( shiftedindex, kind=C_LONG_DOUBLE )

end function f_ArraysWavenumberNegativeLong
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! function f_ArraysWavenumberLong
!> Same as f_ArraysWavenumber, but in extended precision 
!> @param k array index for the wave vector component
!> @param K0 minimum wavenumber for the wave vector component
!> @param KMAX maximum wavenumber for the wave vector component, needed to determine
!! if the wavenumber is positive ( <= KMAX ) or negative ( > KMAX ).
pure function f_ArraysWavenumberLong ( k, KMAX, K0 ) result ( wvn )
  integer, intent(in)              :: k, KMAX 
  real(kind=C_DOUBLE), intent(in)  :: K0
  real(kind=C_LONG_DOUBLE)         :: wvn

  IFPOSNEG: if ( k <= KMAX ) then 
    wvn = f_ArraysWavenumberPositiveLong ( k, K0 )
  else
    wvn = f_ArraysWavenumberNegativeLong ( k, KMAX, K0 )
  end if IFPOSNEG

end function f_ArraysWavenumberLong
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysDealiasSpherical
!> @param ierr should return 0
subroutine f_ArraysDealiasSpherical ( ierr )

  use g_constants,  only : ZERO, HALF
  use g_parameters, only : NODES_X, B110, B220, B330
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, bvalexist, KMAX

  implicit none

  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_DOUBLE),dimension(1:3)     :: wvnmov
  real(kind=C_DOUBLE),dimension(1:3)     :: wvnlab
  real(kind=C_DOUBLE)                    :: k_mag
!  integer                                :: kmax
  real(kind=C_DOUBLE)                    :: cutoff!, cutoffi, cutoffj, cutoffk

  integer                                :: i, j, k, m, n

!  cutoffi = K0I * real(KI_MAX,kind=C_DOUBLE) / real(B110,kind=C_DOUBLE)
!  cutoffj = K0J * real(KJ_MAX,kind=C_DOUBLE) / real(B220,kind=C_DOUBLE)
!  cutoffk = K0K * real(KK_MAX,kind=C_DOUBLE) / real(B330,kind=C_DOUBLE)
!  cutoff = min ( cutoffi, cutoffj, cutoffk ) 
  cutoff = KMAX

  DOI: do i = i_min_3w, i_max_3w, 1
  wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
    DOJ: do j = j_min_3w, j_max_3w, 1
    wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      DOK: do k = k_min_3w, k_max_3w, 1
      wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )


        wvnlab = ZERO

        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXIST
          end do DON
        end do DOM

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        IFALIASMOV: if ( k_mag > cutoff ) then
          arr_u1_3w(:,k,j,i) = ZERO 
          arr_u2_3w(:,k,j,i) = ZERO 
          arr_u3_3w(:,k,j,i) = ZERO 
        end if IFALIASMOV

      end do DOK
    end do DOJ
  end do DOI

end subroutine f_ArraysDealiasSpherical
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_ArraysLESCutoff
!> @param ierr should return 0
subroutine f_ArraysLESCutoff ( ierr )

  use g_constants,  only : ZERO, HALF
  use g_parameters, only : NODES_X, B110, B220, B330, LES_CUTOFF
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, bvalexist, KMAX

  implicit none

  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_DOUBLE),dimension(1:3)     :: wvnmov
  real(kind=C_DOUBLE),dimension(1:3)     :: wvnlab
  real(kind=C_DOUBLE)                    :: k_mag
!  integer                                :: kmax
  real(kind=C_DOUBLE)                    :: cutoff!, cutoffi, cutoffj, cutoffk

  integer                                :: i, j, k, m, n

!  cutoffi = K0I * real(KI_MAX,kind=C_DOUBLE) / real(B110,kind=C_DOUBLE)
!  cutoffj = K0J * real(KJ_MAX,kind=C_DOUBLE) / real(B220,kind=C_DOUBLE)
!  cutoffk = K0K * real(KK_MAX,kind=C_DOUBLE) / real(B330,kind=C_DOUBLE)
!  cutoff = min ( cutoffi, cutoffj, cutoffk ) 
  cutoff = LES_CUTOFF * KMAX

  DOI: do i = i_min_3w, i_max_3w, 1
  wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
    DOJ: do j = j_min_3w, j_max_3w, 1
    wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      DOK: do k = k_min_3w, k_max_3w, 1
      wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )


        wvnlab = ZERO

        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXIST
          end do DON
        end do DOM

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        IFALIASMOV: if ( k_mag > cutoff ) then
          arr_u1_3w(:,k,j,i) = ZERO 
          arr_u2_3w(:,k,j,i) = ZERO 
          arr_u3_3w(:,k,j,i) = ZERO 
        end if IFALIASMOV

      end do DOK
    end do DOJ
  end do DOI

end subroutine f_ArraysLESCutoff
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
end module f_arrays
