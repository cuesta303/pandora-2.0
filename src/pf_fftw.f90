!----------------------------------------------------------------------
! module f_fftw
!> Module containing all calls to FFTW for DFTs and transforms

module f_fftw

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


include 'fftw3-mpi.f03'

private

!> FFT_FORWARD is for the real-space to Fourier-space direction
!> FFT_BACKWARD is for the Fourier-space to real-space direction
integer,parameter                                       :: FFT_FORWARD = -1, FFT_BACKWARD = 1
!> FFTW plans
!> Fourier-space to real-space transform in k, j and i directions
!> Real-space to Fourier-space transform in k, j and i directions
type(C_PTR)                                             :: plan_k_fs2rs, plan_j_fs2rs, plan_i_fs2rs, &
                                                           plan_k_rs2fs, plan_j_rs2fs, plan_i_rs2fs, &
                                                           plan_transpose_kjjk, plan_transpose_jiij, &
                                                           plan_transpose_ijji, plan_transpose_jkkj
type(C_PTR)                                             :: p_k_in, p_k_out, &
                                                           p_j_in, p_j_out, &
                                                           p_i_complex, p_i_real
!> C pointers to transpose arrays if using same array for in and out
type(C_PTR)                                             :: p_k_transposebuffer, p_j_transposebuffer
!> C pointers to transpose arrays if using separate in and out arrays
type(C_PTR)                                             :: p_kj_transposebuffer, p_ji_transposebuffer
type(C_PTR)                                             :: p_jk_transposebuffer, p_ij_transposebuffer
integer(C_INTPTR_T)                                     :: ptrdifft_kjjk_n0, ptrdifft_kjjk_n1, & 
                                                           ptrdifft_jiij_n0, ptrdifft_jiij_n1, &
                                                           ptrdifft_kjjk_howmany, ptrdifft_kjjk_block0, &
                                                           ptrdifft_kjjk_block1, ptrdifft_jiij_howmany, &
                                                           ptrdifft_jiij_block0, ptrdifft_jiij_block1
!> The following are for the reverse direction and not strictly necessary, but make the naming convention more intuitive
integer(C_INTPTR_T)                                     :: ptrdifft_ijji_n0, ptrdifft_ijji_n1, & 
                                                           ptrdifft_jkkj_n0, ptrdifft_jkkj_n1, &
                                                           ptrdifft_ijji_howmany, ptrdifft_ijji_block0, &
                                                           ptrdifft_ijji_block1, ptrdifft_jkkj_howmany, &
                                                           ptrdifft_jkkj_block0, ptrdifft_jkkj_block1
complex(C_DOUBLE_COMPLEX),pointer,dimension(:,:)        :: arr_k_in, arr_k_out, &
                                                           arr_j_in, arr_j_out
!> Only give access to input arrays
real(C_DOUBLE),pointer,public,dimension(:,:,:)          :: arr_k_in_da, arr_j_in_da, arr_j_out_da, arr_k_out_da
!> Fortran pointer arrays if using same array for in and out
real(C_DOUBLE),pointer,dimension(:,:)                   :: arr_k_transposebuffer, arr_j_transposebuffer
!> Fortran pointer arrays if using separate in and out arrays
real(C_DOUBLE),pointer,dimension(:,:)                   :: arr_kj_transposebuffer, arr_jk_transposebuffer, &
                                                           arr_ji_transposebuffer, arr_ij_transposebuffer
real(C_DOUBLE),pointer,dimension(:,:,:)                 :: arr_kjjk_transposein, arr_jiij_transposein, arr_jiij_transposeout, &
                                                           arr_ijji_transposein, arr_ijji_transposeout, arr_jkkj_transposeout
!> Only give access to output arrays
real(C_DOUBLE),pointer,public,dimension(:,:,:)          :: arr_kjjk_transposeout, arr_jkkj_transposein
! These are complex arrays, but we want to access them directly, so the best way to do this is to treat them as real arrays
complex(C_DOUBLE_COMPLEX),pointer,dimension(:,:)        :: arr_i_complex
real(C_DOUBLE),pointer,public,dimension(:,:,:)          :: arr_i_complex_da
!> Only give access to output arrays
real(C_DOUBLE),pointer,public,dimension(:,:)            :: arr_i_real

integer                                                 :: size_k, size_j, size_i_complex, size_i_real, &
                                                           size_k_transpose, size_j_transpose
real(C_DOUBLE),save,private                             :: normalisationfactor                                     

public :: f_FftwInit
public :: f_FftwFinalise
public :: f_Fftw_k_Fs2Rs
public :: f_Fftw_ji_Fs2Rs
public :: f_Fftw_ji_Fs2RsPrepare
public :: f_Fftw_ji_Fs2RsFinalise
public :: f_Fftw_ji_Rs2Fs
public :: f_Fftw_k_Rs2Fs
public :: f_Fftw_FftOnly_k_Fs2Rs
public :: f_Fftw_FftOnly_k_Rs2Fs

contains
!---------------------------------------------------------------------------------
! subroutine f_FftwInit ( ierr )
!> Allocate buffer arrays and set up plans for Fourier transforms and global 
!! transpose
!> @param ierr should return 0
subroutine f_FftwInit ( ierr )

  use g_constants,   only : ZERO, ONE
  use g_domain,      only : NODES_KI, NODES_KJ, NODES_KK
  use g_parameters,  only : NODES_X, NODES_Y, NODES_Z, MYRANK
  use f_arrays,      only : j_width_3w, z_width_1r2w, i_width_1r2w, y_width_2r, &
                            j_maxwidth_3w, z_maxwidth_1r2w, i_maxwidth_1r2w, y_maxwidth_2r, & 
                            nprocsj_3w, nprocsz_1r2w, nprocsi_1r2w, nprocsy_2r, &
                            nodes_kj_transpose, nodes_z_transpose, &
                            nodes_ki_transpose, nodes_y_transpose, &
                            kjcomm, jicomm

  implicit none

  PetscErrorCode,intent(inout)                                  :: ierr
  integer,dimension(1)                                   :: fftsizek, fftsizej, fftsizei

! 1st transform in k direction
! need NODES_Z * NODES_KJ buffer array
! will then be transformed into a NODES_KJ * NODES_Z buffer array

! plan_many_dfts rather than the MPI version (no need for it)
! then fftw_mpi_transpose

! 2nd transform in j direction
! need NODES_Y * NODES_KI buffer aray
! will then be transformed into a NODES_KI * NODES_Y array

! However need 6 of the latter so can compute the cross product in real space
! Will also need a DMDA in 2D for the particles

!> Parameters for size of FFTs

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising f_fftw module \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Determining dimensions of Fourier transforms. \n', ierr )
  size_k = NODES_Z * j_width_3w
  fftsizek(1) = NODES_Z
  size_j = NODES_Y * i_width_1r2w
  fftsizej(1) = NODES_Y
  size_i_complex = (NODES_X/2 + 1) * y_width_2r
  size_i_real = NODES_X * y_width_2r
  fftsizei(1) = NODES_X

!> Parameters for transpose from wave-space to real-space

  call PetscPrintf( PETSC_COMM_WORLD, '    | Determining dimensions of global transpose. \n', ierr )
  size_k_transpose = z_maxwidth_1r2w * nodes_kj_transpose ! (nprocsj_3w * j_maxwidth_3w)
  ptrdifft_kjjk_n0 = nodes_kj_transpose ! NODES_KJ
  ptrdifft_kjjk_n1 = nodes_z_transpose ! NODES_Z
  ptrdifft_kjjk_howmany = 2
  ptrdifft_kjjk_block0 = j_maxwidth_3w
  ptrdifft_kjjk_block1 = z_maxwidth_1r2w

  size_j_transpose = y_maxwidth_2r * nodes_ki_transpose ! NODES_KI
  ptrdifft_jiij_n0 = nodes_ki_transpose ! NODES_KI 
  ptrdifft_jiij_n1 = nodes_y_transpose ! NODES_Y
  ptrdifft_jiij_howmany = 2
  ptrdifft_jiij_block0 = i_maxwidth_1r2w
  ptrdifft_jiij_block1 = y_maxwidth_2r

!> Transpose from real-space to wave-space has same parameters as the other way around, but in other direction (0 => 1, 1 => 0)

  ptrdifft_ijji_n0 = ptrdifft_jiij_n1 
  ptrdifft_ijji_n1 = ptrdifft_jiij_n0
  ptrdifft_ijji_howmany = ptrdifft_jiij_howmany 
  ptrdifft_ijji_block0 = ptrdifft_jiij_block1
  ptrdifft_ijji_block1 = ptrdifft_jiij_block0

  ptrdifft_jkkj_n0 = ptrdifft_kjjk_n1
  ptrdifft_jkkj_n1 = ptrdifft_kjjk_n0
  ptrdifft_jkkj_howmany = ptrdifft_kjjk_howmany
  ptrdifft_jkkj_block0 = ptrdifft_kjjk_block1
  ptrdifft_jkkj_block1 = ptrdifft_kjjk_block0

!> Allocate arrays for FFTs
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating buffer arrays for Fourier transforms. \n', ierr )
  p_k_in = fftw_alloc_complex ( int(size_k, C_SIZE_T) )
  p_k_out = fftw_alloc_complex ( int(size_k, C_SIZE_T) )
  p_j_in = fftw_alloc_complex ( int(size_j, C_SIZE_T) )
  p_j_out = fftw_alloc_complex ( int(size_j, C_SIZE_T) )
  p_i_complex = fftw_alloc_complex ( int(size_i_complex, C_SIZE_T) )
  p_i_real = fftw_alloc_real ( int(size_i_real, C_SIZE_T) )

  call c_f_pointer ( p_k_in, arr_k_in, [NODES_Z, j_width_3w])
  call c_f_pointer ( p_k_out, arr_k_out, [NODES_Z, j_width_3w])
  call c_f_pointer ( p_j_in, arr_j_in, [NODES_Y, i_width_1r2w])
  call c_f_pointer ( p_j_out, arr_j_out, [NODES_Y, i_width_1r2w])
  call c_f_pointer ( p_i_complex, arr_i_complex, [(NODES_X/2 + 1), y_width_2r])
  call c_f_pointer ( p_i_real, arr_i_real, [NODES_X, y_width_2r])

!> The following are pointers to the exact same arrays, only with the advantage of having direct access
  call c_f_pointer ( p_k_in, arr_k_in_da, [2, NODES_Z, j_width_3w])
  call c_f_pointer ( p_k_out, arr_k_out_da, [2, NODES_Z, j_width_3w])
  call c_f_pointer ( p_j_in, arr_j_in_da, [2, NODES_Y, i_width_1r2w])
  call c_f_pointer ( p_j_out, arr_j_out_da, [2, NODES_Y, i_width_1r2w])
  call c_f_pointer ( p_i_complex, arr_i_complex_da, [2, (NODES_X/2 + 1), y_width_2r])


  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning Fourier transforms in k direction. \n', ierr )

! FFTs in k direction
  plan_k_fs2rs = fftw_plan_many_dft ( &
  1,             &       !int rank, 1 for 1D transform 
  fftsizek,      &       !const int *n, length of each 1D transform (expressed as an array!)
  j_width_3w, &       !int howmany, number of 1D transforms
  arr_k_in,      &       !fftw_complex *in, input array
  fftsizek,      &       !const int *inembed, length of each 1D transform
  1,             &       !int istride, distance between two elements; must be 1 for contiguous memory
  NODES_Z,       &       !int idist, distance between first element of first array and first element of second array
  arr_k_out,     &       !fftw_complex *out, output array
  fftsizek,      &       !const int *onembed, length of each 1D transform
  1,             &       !int ostride, distance between two elements; must be 1 for contiguous memory
  NODES_Z,       &       !int odist, distance between first element of first array and first element of second array
  FFT_BACKWARD,  &       !int sign, FFT_FORWARD or FFT_BACKWARD
  FFTW_EXHAUSTIVE )      !unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  plan_k_rs2fs = fftw_plan_many_dft ( &
  1,             &       !int rank, 1 for 1D transform 
  fftsizek,      &       !const int *n, length of each 1D transform (expressed as an array!)
  j_width_3w, &       !int howmany, number of 1D transforms
  arr_k_in,      &       !fftw_complex *in, input array
  fftsizek,      &       !const int *inembed, length of each 1D transform
  1,             &       !int istride, distance between two elements; must be 1 for contiguous memory
  NODES_Z,       &       !int idist, distance between first element of first array and first element of second array
  arr_k_out,     &       !fftw_complex *out, output array
  fftsizek,      &       !const int *onembed, length of each 1D transform
  1,             &       !int ostride, distance between two elements; must be 1 for contiguous memory
  NODES_Z,       &       !int odist, distance between first element of first array and first element of second array
  FFT_FORWARD,   &       !int sign, FFT_FORWARD or FFT_BACKWARD
  FFTW_EXHAUSTIVE )      !unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning Fourier transforms in j direction. \n', ierr )

! FFTs in j direction
  plan_j_fs2rs = fftw_plan_many_dft ( &
  1,               &       !int rank, 1 for 1D transform 
  fftsizej,        &       !const int *n, length of each 1D transform (expressed as an array!)
  i_width_1r2w, &       !int howmany, number of 1D transforms
  arr_j_in,        &       !fftw_complex *in, input array
  fftsizej,        &       !const int *inembed, length of each 1D transform
  1,               &       !int istride, distance between two elements; must be 1 for contiguous memory
  NODES_Y,         &       !int idist, distance between first element of first array and first element of second array
  arr_j_out,       &       !fftw_complex *out, output array
  fftsizej,        &       !const int *onembed, length of each 1D transform
  1,               &       !int ostride, distance between two elements; must be 1 for contiguous memory
  NODES_Y,         &       !int odist, distance between first element of first array and first element of second array
  FFT_BACKWARD,    &       !int sign, FFT_FORWARD or FFT_BACKWARD
  FFTW_EXHAUSTIVE )      !unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  plan_j_rs2fs = fftw_plan_many_dft ( &
  1,               &       !int rank, 1 for 1D transform 
  fftsizej,        &       !const int *n, length of each 1D transform (expressed as an array!)
  i_width_1r2w, &       !int howmany, number of 1D transforms
  arr_j_in,        &       !fftw_complex *in, input array
  fftsizej,        &       !const int *inembed, length of each 1D transform
  1,               &       !int istride, distance between two elements; must be 1 for contiguous memory
  NODES_Y,         &       !int idist, distance between first element of first array and first element of second array
  arr_j_out,       &       !fftw_complex *out, output array
  fftsizej,        &       !const int *onembed, length of each 1D transform
  1,               &       !int ostride, distance between two elements; must be 1 for contiguous memory
  NODES_Y,         &       !int odist, distance between first element of first array and first element of second array
  FFT_FORWARD,     &       !int sign, FFT_FORWARD or FFT_BACKWARD
  FFTW_EXHAUSTIVE )      !unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning Fourier transforms in i direction. \n', ierr )

! FFTs in i direction
  plan_i_fs2rs = fftw_plan_many_dft_c2r ( &
  1,               &       !int rank, 1 for 1D transform 
  fftsizei,        &       !const int *n, length of each 1D transform (expressed as an array!)
  y_width_2r,   &       !int howmany, number of 1D transforms
  arr_i_complex,   &       !fftw_complex *in, input array
  fftsizei,        &       !const int *inembed, length of each 1D transform
  1,               &       !int istride, distance between two elements; must be 1 for contiguous memory
  NODES_X/2 + 1,   &       !int idist, distance between first element of first array and first element of second array
  arr_i_real,      &       !double *out, output array
  fftsizei,        &       !const int *onembed, length of each 1D transform
  1,               &       !int ostride, distance between two elements; must be 1 for contiguous memory
  NODES_X,         &       !int odist, distance between first element of first array and first element of second array
  FFTW_EXHAUSTIVE )      !unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  plan_i_rs2fs = fftw_plan_many_dft_r2c ( &
  1,               &       !int rank, 1 for 1D transform 
  fftsizei,        &       !const int *n, length of each 1D transform (expressed as an array!)
  y_width_2r,   &       !int howmany, number of 1D transforms
  arr_i_real,      &       !fftw_complex *in, input array
  fftsizei,        &       !const int *inembed, length of each 1D transform
  1,               &       !int istride, distance between two elements; must be 1 for contiguous memory
  NODES_X,         &       !int idist, distance between first element of first array and first element of second array
  arr_i_complex,   &       !fftw_complex *out, output array
  fftsizej,        &       !const int *onembed, length of each 1D transform
  1,               &       !int ostride, distance between two elements; must be 1 for contiguous memory
  NODES_X/2 + 1,   &       !int odist, distance between first element of first array and first element of second array
  FFTW_EXHAUSTIVE )      !unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  !> FFTW computes unnormalised transformed. Compute the normalisation factor here.
!  normalisationfactor = ONE / ( real ( NODES_X * NODES_Y * NODES_Z, kind=C_DOUBLE))
  normalisationfactor = ONE / ( real ( NODES_X, kind=C_DOUBLE) * &
                                real ( NODES_Y, kind=C_DOUBLE) * &
                                real ( NODES_Z, kind=C_DOUBLE) )


  !> @todo insert if condition here
  !> Transpose using the same buffer array.
  call f_FftwInitTransposeSameBuffer ( ierr )

!  !> Transpose using separate buffer arrays.
!  call f_FftwInitTransposeSeparateBuffers ( ierr )

end subroutine f_FftwInit
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_FftwInitTransposeSameBuffer 
!> Allocate buffer arrays and set up plans for Fourier transforms and global 
!! transpose
!> @param ierr should return 0
subroutine f_FftwInitTransposeSameBuffer ( ierr )

  use g_constants,   only : ZERO, ONE
  use g_domain,      only : NODES_KI, NODES_KJ, NODES_KK
  use g_parameters,  only : NODES_X, NODES_Y, NODES_Z, MYRANK
  use f_arrays,      only : j_width_3w, z_width_1r2w, i_width_1r2w, y_width_2r, &
                            j_maxwidth_3w, z_maxwidth_1r2w, i_maxwidth_1r2w, y_maxwidth_2r, & 
                            nprocsj_3w, nprocsz_1r2w, nprocsi_1r2w, nprocsy_2r, &
                            nodes_kj_transpose, nodes_z_transpose, &
                            nodes_ki_transpose, nodes_y_transpose, &
                            kjcomm, jicomm

  implicit none

  PetscErrorCode,intent(inout)                                  :: ierr

!> Allocate arrays for transpose.
!> We are using the same array for input and output.
!> For convenience we define an input pointer and an output pointer, so access becomes more natural.
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating buffer arrays for global transpose. \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Using same array for input and output. \n', ierr )
  size_k_transpose = z_maxwidth_1r2w * nodes_kj_transpose ! (nprocsj_3w * j_maxwidth_3w)
  ptrdifft_kjjk_n0 = nodes_kj_transpose ! NODES_KJ
  ptrdifft_kjjk_n1 = nodes_z_transpose ! NODES_Z
  ptrdifft_kjjk_howmany = 2
  ptrdifft_kjjk_block0 = j_maxwidth_3w
  ptrdifft_kjjk_block1 = z_maxwidth_1r2w

  size_j_transpose = y_maxwidth_2r * nodes_ki_transpose ! NODES_KI
  ptrdifft_jiij_n0 = nodes_ki_transpose ! NODES_KI 
  ptrdifft_jiij_n1 = nodes_y_transpose ! NODES_Y
  ptrdifft_jiij_howmany = 2
  ptrdifft_jiij_block0 = i_maxwidth_1r2w
  ptrdifft_jiij_block1 = y_maxwidth_2r

  p_k_transposebuffer = fftw_alloc_complex ( int(size_k_transpose, C_SIZE_T) )
  call c_f_pointer ( p_k_transposebuffer, arr_k_transposebuffer, [2, size_k_transpose])
  call c_f_pointer ( p_k_transposebuffer, arr_kjjk_transposein, [2, nodes_z_transpose, j_maxwidth_3w])
  call c_f_pointer ( p_k_transposebuffer, arr_kjjk_transposeout, [2, nodes_kj_transpose, z_maxwidth_1r2w])
  !> Swap in and out in real-to-Fourier direction. Just different names for the same arrays.
  call c_f_pointer ( p_k_transposebuffer, arr_jkkj_transposein, [2, nodes_kj_transpose, z_maxwidth_1r2w])
  call c_f_pointer ( p_k_transposebuffer, arr_jkkj_transposeout, [2, nodes_z_transpose, j_maxwidth_3w])
  p_j_transposebuffer = fftw_alloc_complex ( int(size_j_transpose, C_SIZE_T) )
  call c_f_pointer ( p_j_transposebuffer, arr_j_transposebuffer, [2, size_j_transpose])
  call c_f_pointer ( p_j_transposebuffer, arr_jiij_transposein, [2, nodes_y_transpose, i_maxwidth_1r2w])
  call c_f_pointer ( p_j_transposebuffer, arr_jiij_transposeout, [2, nodes_ki_transpose, y_maxwidth_2r])
  !> Swap in and out in real-to-Fourier direction. Just different names for the same arrays.
  call c_f_pointer ( p_j_transposebuffer, arr_ijji_transposein, [2, nodes_ki_transpose, y_maxwidth_2r])
  call c_f_pointer ( p_j_transposebuffer, arr_ijji_transposeout, [2, nodes_y_transpose, i_maxwidth_1r2w])


  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in kj direction. \n', ierr )

  plan_transpose_kjjk = fftw_mpi_plan_many_transpose ( &
  ptrdifft_kjjk_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_kjjk_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_kjjk_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_kjjk_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_kjjk_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_k_transposebuffer, &    !  double *in, input array
  arr_k_transposebuffer, &    !  double *out, output array, can be the same as input array
  kjcomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in ji direction. \n', ierr )

  plan_transpose_jiij = fftw_mpi_plan_many_transpose ( &
  ptrdifft_jiij_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_jiij_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_jiij_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_jiij_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_jiij_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_j_transposebuffer, &    !  double *in, input array
  arr_j_transposebuffer, &    !  double *out, output array, can be the same as input array
  jicomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in ij direction. \n', ierr )

  plan_transpose_ijji = fftw_mpi_plan_many_transpose ( &
  ptrdifft_ijji_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_ijji_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_ijji_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_ijji_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_ijji_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_j_transposebuffer, &    !  double *in, input array
  arr_j_transposebuffer, &    !  double *out, output array, can be the same as input array
  jicomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in jk direction. \n', ierr )

  plan_transpose_jkkj = fftw_mpi_plan_many_transpose ( &
  ptrdifft_jkkj_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_jkkj_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_jkkj_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_jkkj_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_jkkj_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_k_transposebuffer, &    !  double *in, input array
  arr_k_transposebuffer, &    !  double *out, output array, can be the same as input array
  kjcomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

end subroutine f_FftwInitTransposeSameBuffer 
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FftwInitTransposeSeparateBuffers 
!> Allocate buffer arrays and set up plans for Fourier transforms and global 
!! transpose
!> @param ierr should return 0
subroutine f_FftwInitTransposeSeparateBuffers ( ierr )

  use g_constants,   only : ZERO, ONE
  use g_domain,      only : NODES_KI, NODES_KJ, NODES_KK
  use g_parameters,  only : NODES_X, NODES_Y, NODES_Z, MYRANK
  use f_arrays,      only : j_width_3w, z_width_1r2w, i_width_1r2w, y_width_2r, &
                            j_maxwidth_3w, z_maxwidth_1r2w, i_maxwidth_1r2w, y_maxwidth_2r, & 
                            nprocsj_3w, nprocsz_1r2w, nprocsi_1r2w, nprocsy_2r, &
                            nodes_kj_transpose, nodes_z_transpose, &
                            nodes_ki_transpose, nodes_y_transpose, &
                            kjcomm, jicomm

  implicit none

  PetscErrorCode,intent(inout)                                  :: ierr

!> Allocate arrays for transpose.
!> We are using the same array for input and output.
!> For convenience we define an input pointer and an output pointer, so access becomes more natural.
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating buffer arrays for global transpose. \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Using separate input and output arrays. \n', ierr )
  size_k_transpose = z_maxwidth_1r2w * nodes_kj_transpose ! (nprocsj_3w * j_maxwidth_3w)
  ptrdifft_kjjk_n0 = nodes_kj_transpose ! NODES_KJ
  ptrdifft_kjjk_n1 = nodes_z_transpose ! NODES_Z
  ptrdifft_kjjk_howmany = 2
  ptrdifft_kjjk_block0 = j_maxwidth_3w
  ptrdifft_kjjk_block1 = z_maxwidth_1r2w

  size_j_transpose = y_maxwidth_2r * nodes_ki_transpose ! NODES_KI
  ptrdifft_jiij_n0 = nodes_ki_transpose ! NODES_KI 
  ptrdifft_jiij_n1 = nodes_y_transpose ! NODES_Y
  ptrdifft_jiij_howmany = 2
  ptrdifft_jiij_block0 = i_maxwidth_1r2w
  ptrdifft_jiij_block1 = y_maxwidth_2r

  p_kj_transposebuffer = fftw_alloc_complex ( int(size_k_transpose, C_SIZE_T) )
  p_jk_transposebuffer = fftw_alloc_complex ( int(size_k_transpose, C_SIZE_T) )
  !> Vectors as required by FFTW
  call c_f_pointer ( p_kj_transposebuffer, arr_kj_transposebuffer, [2, size_k_transpose])
  call c_f_pointer ( p_jk_transposebuffer, arr_jk_transposebuffer, [2, size_k_transpose])
  !> Aliases for accessing in a more natural way
  call c_f_pointer ( p_kj_transposebuffer, arr_kjjk_transposein, [2, nodes_z_transpose, j_maxwidth_3w])
  call c_f_pointer ( p_jk_transposebuffer, arr_kjjk_transposeout, [2, nodes_kj_transpose, z_maxwidth_1r2w])
  call c_f_pointer ( p_jk_transposebuffer, arr_jkkj_transposein, [2, nodes_kj_transpose, z_maxwidth_1r2w])
  call c_f_pointer ( p_kj_transposebuffer, arr_jkkj_transposeout, [2, nodes_z_transpose, j_maxwidth_3w])
  p_ji_transposebuffer = fftw_alloc_complex ( int(size_j_transpose, C_SIZE_T) )
  p_ij_transposebuffer = fftw_alloc_complex ( int(size_j_transpose, C_SIZE_T) )
  !> Vectors as required by FFTW
  call c_f_pointer ( p_ji_transposebuffer, arr_ji_transposebuffer, [2, size_j_transpose])
  call c_f_pointer ( p_ij_transposebuffer, arr_ij_transposebuffer, [2, size_j_transpose])
  !> Aliases for accessing in a more natural way
  call c_f_pointer ( p_ji_transposebuffer, arr_jiij_transposein, [2, nodes_y_transpose, i_maxwidth_1r2w])
  call c_f_pointer ( p_ij_transposebuffer, arr_jiij_transposeout, [2, nodes_ki_transpose, y_maxwidth_2r])
  call c_f_pointer ( p_ij_transposebuffer, arr_ijji_transposein, [2, nodes_ki_transpose, y_maxwidth_2r])
  call c_f_pointer ( p_ji_transposebuffer, arr_ijji_transposeout, [2, nodes_y_transpose, i_maxwidth_1r2w])


  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in kj direction. \n', ierr )

  plan_transpose_kjjk = fftw_mpi_plan_many_transpose ( &
  ptrdifft_kjjk_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_kjjk_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_kjjk_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_kjjk_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_kjjk_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_kj_transposebuffer, &    !  double *in, input array
  arr_jk_transposebuffer, &    !  double *out, output array, can be the same as input array
  kjcomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in ji direction. \n', ierr )

  plan_transpose_jiij = fftw_mpi_plan_many_transpose ( &
  ptrdifft_jiij_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_jiij_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_jiij_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_jiij_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_jiij_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_ji_transposebuffer, &   !  double *in, input array
  arr_ij_transposebuffer, &   !  double *out, output array, can be the same as input array
  jicomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in ij direction. \n', ierr )

  plan_transpose_ijji = fftw_mpi_plan_many_transpose ( &
  ptrdifft_ijji_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_ijji_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_ijji_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_ijji_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_ijji_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_ij_transposebuffer, &   !  double *in, input array
  arr_ji_transposebuffer, &   !  double *out, output array, can be the same as input array
  jicomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

  call PetscPrintf( PETSC_COMM_WORLD, '    | Planning global transpose in jk direction. \n', ierr )

  plan_transpose_jkkj = fftw_mpi_plan_many_transpose ( &
  ptrdifft_jkkj_n0,      &    !  ptrdiff_t n0, first array dimension
  ptrdifft_jkkj_n1,      &    !  ptrdiff_t n1, second array dimension
  ptrdifft_jkkj_howmany, &    !  ptrdiff_t howmany, 2 for complex numbers
  ptrdifft_jkkj_block0,  &    !  ptrdiff_t block0, n0 by n1 input distributed along n0 with block size block0
  ptrdifft_jkkj_block1,  &    !  ptrdiff_t block1, n1 by n0 output distributed along n1 with block size block1
  arr_jk_transposebuffer, &   !  double *in, input array
  arr_kj_transposebuffer, &   !  double *out, output array, can be the same as input array
  kjcomm,                &    !  MPI_Comm comm, should be the correct subset corresponding to the respective direction
  FFTW_EXHAUSTIVE )           !  unsigned flags Planner flag, reasonable choice FFTW_EXHAUSTIVE

end subroutine f_FftwInitTransposeSeparateBuffers 
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
! subroutine f_FftwFinalise ( ierr )
!> Deallocate buffer arrays and destroy all plans.
!> @param ierr should return 0
subroutine f_FftwFinalise ( ierr )

  implicit none

  PetscErrorCode,intent(inout)                                  :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Finalising f_fftw module \n', ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Destroying FFT plans \n', ierr )
  call fftw_destroy_plan ( plan_k_fs2rs )
  call fftw_destroy_plan ( plan_k_rs2fs )
  call fftw_destroy_plan ( plan_j_fs2rs )
  call fftw_destroy_plan ( plan_j_rs2fs )
  call fftw_destroy_plan ( plan_i_fs2rs )
  call fftw_destroy_plan ( plan_i_rs2fs )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Destroying transpose plans \n', ierr )
  call fftw_destroy_plan ( plan_transpose_kjjk )
  call fftw_destroy_plan ( plan_transpose_jiij )
  call fftw_destroy_plan ( plan_transpose_ijji )
  call fftw_destroy_plan ( plan_transpose_jkkj )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating buffer arrays \n', ierr )
  call fftw_free ( p_k_in )
  call fftw_free ( p_k_out )
  call fftw_free ( p_j_in )
  call fftw_free ( p_j_out )
  call fftw_free ( p_i_complex )
  call fftw_free ( p_i_real )

!> @todo this depends on which routine is chosen
!  call fftw_free ( p_k_transposebuffer )
!  call fftw_free ( p_j_transposebuffer )

end subroutine f_FftwFinalise

!---------------------------------------------------------------------------------
! subroutine f_Fftw_k_Fs2Rs
!> Perform an FFT in k direction on a kj slice and then transpose the slice from
!! kj to jk.
!> @param ierr should return 0
subroutine f_Fftw_k_Fs2Rs ( ierr )

  use f_arrays,     only : j_width_3w, nprocsz_1r2w, z_allwidths_1r2w, z_maxwidth_1r2w

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr
  integer                                                 :: proc_k, j, k, component, &
                                                             kfft, ktranspose

!> Perform a Fourier-space to real-space FFT on a kj slice in k direction

  call fftw_execute_dft ( plan_k_fs2rs, arr_k_in, arr_k_out )

!> Copy FFT output array to transpose input array

  DOJ: do j = 1, j_width_3w, 1
    kfft = 0
    DOPROCK: do proc_k = 0, nprocsz_1r2w - 1, 1
      DOK: do k = 1, z_allwidths_1r2w(proc_k), 1

        kfft = kfft + 1
        ktranspose = (proc_k * z_maxwidth_1r2w) + k

        DOCOMPONENT: do component = 1, 2, 1

          arr_kjjk_transposein(component, ktranspose, j) = arr_k_out_da(component, kfft, j)

        end do DOCOMPONENT

      end do DOK
    end do DOPROCK
  end do DOJ

!> Transpose the slice from kj to jk orientation
!  call fftw_mpi_execute_r2r ( plan_transpose_kjjk, arr_kj_transposebuffer, arr_jk_transposebuffer )
  call fftw_mpi_execute_r2r ( plan_transpose_kjjk, arr_k_transposebuffer, arr_k_transposebuffer )

end subroutine f_Fftw_k_Fs2Rs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_ji_Fs2Rs
!> Perform an FFT in j direction on a ji slice, then transpose the slice to ij orientation.
!! Zeropad the ij slice and perform an FFT in i direction. The 2D slice is now in real space.
!> @param ierr should return 0
subroutine f_Fftw_ji_Fs2Rs ( ierr )

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr

!> Perform an FFT in j direction on a ji slice, then transpose the slice to ij
!! orientation. Zeropad the ij slice.
  call f_Fftw_j_Fs2Rs ( ierr )

!  write(*,*) 'arr_i_complex_da(1,1,:): ', arr_i_complex_da(1,1,:)
!  write(*,*) 'arr_i_complex_da(2,1,:): ', arr_i_complex_da(2,1,:)

!> Perform a Fourier-space to real-space FFT on an ij slice in i direction
  call fftw_execute_dft_c2r ( plan_i_fs2rs, arr_i_complex, arr_i_real )

end subroutine f_Fftw_ji_Fs2Rs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_ji_Fs2RsPrepare
!> Perform an FFT in j direction on a ji slice, then transpose the slice to ij orientation.
!! Zeropad the ij slice.
!> @param ierr should return 0
subroutine f_Fftw_ji_Fs2RsPrepare ( ierr )

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr

!> Perform an FFT in j direction on a ji slice, then transpose the slice to ij
!! orientation. Zeropad the ij slice.
  call f_Fftw_j_Fs2Rs ( ierr )

end subroutine f_Fftw_ji_Fs2RsPrepare
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_ji_Fs2RsFinalise
!> Perform an FFT in i direction. The 2D slice is now in real space.
!> @param ierr should return 0
subroutine f_Fftw_ji_Fs2RsFinalise ( ierr )

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr

!  write(*,*) 'arr_i_complex_da(1,1,:): ', arr_i_complex_da(1,1,:)
!  write(*,*) 'arr_i_complex_da(2,1,:): ', arr_i_complex_da(2,1,:)

!> Perform a Fourier-space to real-space FFT on an ij slice in i direction
  call fftw_execute_dft_c2r ( plan_i_fs2rs, arr_i_complex, arr_i_real )

end subroutine f_Fftw_ji_Fs2RsFinalise
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_j_Fs2Rs
!> Perform an FFT in j direction on a ji slice, then transpose the slice to ij orientation.
!! Zeropad the ij slice. 
!> @param ierr should return 0
subroutine f_Fftw_j_Fs2Rs ( ierr )

  use g_parameters, only : NODES_X
  use g_constants,  only : ZERO
  use g_domain,     only : KI_MAX
  use f_arrays,     only : i_width_1r2w, nprocsy_2r, y_allwidths_2r, y_maxwidth_2r, &
                           y_width_2r, nprocsi_1r2w, i_allwidths_1r2w, i_maxwidth_1r2w
  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr
  integer                                                 :: proc_j, proc_i, i, j, component, &
                                                             jfft, jtranspose, ifft, itranspose

!> Perform a Fourier-space to real-space FFT on a ji slice in j direction
  call fftw_execute_dft ( plan_j_fs2rs, arr_j_in, arr_j_out )

!> Copy FFT output array to transpose input array

  DOI: do i = 1, i_width_1r2w, 1
    jfft = 0
    DOPROCJ: do proc_j = 0, nprocsy_2r - 1, 1
      DOJ: do j = 1, y_allwidths_2r(proc_j), 1

        jfft = jfft + 1
        jtranspose = (proc_j * y_maxwidth_2r) + j

        DOCOMPONENT: do component = 1, 2, 1

          arr_jiij_transposein(component, jtranspose, i) = arr_j_out_da(component, jfft, i)

        end do DOCOMPONENT
      end do DOJ
    end do DOPROCJ
  end do DOI

!> Transpose the slice from ji to ij orientation
!  call fftw_mpi_execute_r2r ( plan_transpose_jiij, arr_ji_transposebuffer, arr_ij_transposebuffer )
  call fftw_mpi_execute_r2r ( plan_transpose_jiij, arr_j_transposebuffer, arr_j_transposebuffer )

  !> j direction is now parallelised and in real space
  DOY: do j = 1, y_width_2r, 1

    ifft = 0
    !> Copy positive wavenumbers - there are no negative wavenumbers in the array!!
    DOPROCI: do proc_i = 0, nprocsi_1r2w - 1, 1
      DOIPOSITIVE: do i = 1, i_allwidths_1r2w(proc_i), 1
        !> Copy to buffer

        ifft = ifft + 1
        itranspose = (proc_i * i_maxwidth_1r2w) + i

        arr_i_complex_da(:,ifft,j) = arr_jiij_transposeout(:,itranspose,j)

      end do DOIPOSITIVE
    end do DOPROCI

    !> Zero padding for the highest 1/3 wavenumbers
    !! Start from (KI_MAX + 1) + 1 = KI_MAX + 2, because array indexing starts from 1 rather than 0

      DOIZEROS: do i = KI_MAX + 2, (NODES_X/2) + 1
        ifft = ifft + 1
        arr_i_complex_da(:,ifft,j) = ZERO 
      end do DOIZEROS

  end do DOY

end subroutine f_Fftw_j_Fs2Rs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_ji_Rs2Fs
!> Revert the transforms performed in f_Fftw_ji_Fs2Rs.
!> @param ierr should return 0
subroutine f_Fftw_ji_Rs2Fs ( ierr )

  use g_parameters, only : NODES_X
  use g_constants,  only : ZERO
  use g_domain,     only : KI_MAX
  use f_arrays,     only : i_width_1r2w, nprocsy_2r, y_allwidths_2r, y_maxwidth_2r, &
                           y_width_2r, nprocsi_1r2w, i_allwidths_1r2w, i_maxwidth_1r2w
  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr
  integer                                                 :: proc_j, proc_i, i, j, component, &
                                                             jfft, jtranspose, ifft, itranspose

!> Perform a real-space to Fourier-space FFT on an ij slice in i direction
  call fftw_execute_dft_r2c ( plan_i_rs2fs, arr_i_real, arr_i_complex )

!  write(*,*) 'arr_i_complex_da(1,1,:): ', arr_i_complex_da(1,1,:)
!  write(*,*) 'arr_i_complex_da(2,1,:): ', arr_i_complex_da(2,1,:)

  !> j direction is still parallelised and in real space
  DOY: do j = 1, y_width_2r, 1

    ifft = 0
    !> Copy positive wavenumbers - there are no negative wavenumbers in the array!!
    DOPROCI: do proc_i = 0, nprocsi_1r2w - 1, 1
      DOIPOSITIVE: do i = 1, i_allwidths_1r2w(proc_i), 1
        !> Copy to buffer

        ifft = ifft + 1
        itranspose = (proc_i * i_maxwidth_1r2w) + i

        !> We are going backwards. Swap output and input arrays.
        arr_ijji_transposein(:,itranspose,j) = arr_i_complex_da(:,ifft,j)

      end do DOIPOSITIVE
    end do DOPROCI

    !> Dealiasing for the highest 1/3 wavenumbers.
    !! There is nothing more to do here, just don't copy them.

  end do DOY

!> Transpose the slice from ij to ji orientation
!  call fftw_mpi_execute_r2r ( plan_transpose_ijji, arr_ij_transposebuffer, arr_ji_transposebuffer )
  call fftw_mpi_execute_r2r ( plan_transpose_ijji, arr_j_transposebuffer, arr_j_transposebuffer )

!> Copy transpose output array to FFT input array

  DOI: do i = 1, i_width_1r2w, 1
    jfft = 0
    DOPROCJ: do proc_j = 0, nprocsy_2r - 1, 1
      DOJ: do j = 1, y_allwidths_2r(proc_j), 1

        jfft = jfft + 1
        jtranspose = (proc_j * y_maxwidth_2r) + j

        DOCOMPONENT: do component = 1, 2, 1

          arr_j_in_da(component, jfft, i) = arr_ijji_transposeout (component, jtranspose, i)

        end do DOCOMPONENT
      end do DOJ
    end do DOPROCJ
  end do DOI

!> Perform a real-space to Fourier-space FFT on a ji slice in j direction
  call fftw_execute_dft ( plan_j_rs2fs, arr_j_in, arr_j_out )

end subroutine f_Fftw_ji_Rs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_k_Rs2Fs
!> Revert the transforms performed in f_Fftw_k_Fs2Rs.
!> @param ierr should return 0
subroutine f_Fftw_k_Rs2Fs ( ierr )

  use f_arrays,     only : j_width_3w, nprocsz_1r2w, z_allwidths_1r2w, z_maxwidth_1r2w

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr
  integer                                                 :: proc_k, j, k, component, &
                                                             kfft, ktranspose

!> Transpose the slice from jk to kj orientation
!  call fftw_mpi_execute_r2r ( plan_transpose_jkkj, arr_jk_transposebuffer, arr_kj_transposebuffer )
  call fftw_mpi_execute_r2r ( plan_transpose_jkkj, arr_k_transposebuffer, arr_k_transposebuffer )

!> Copy transpose output array to FFT input array

  DOJ: do j = 1, j_width_3w, 1
    kfft = 0
    DOPROCK: do proc_k = 0, nprocsz_1r2w - 1, 1
      DOK: do k = 1, z_allwidths_1r2w(proc_k), 1

        kfft = kfft + 1
        ktranspose = (proc_k * z_maxwidth_1r2w) + k

        DOCOMPONENT: do component = 1, 2, 1

          arr_k_in_da(component, kfft, j) = arr_kjjk_transposein(component, ktranspose, j)

        end do DOCOMPONENT
      end do DOK
    end do DOPROCK
  end do DOJ

!> Perform a Fourier-space to real-space FFT on a jk slice in k direction
  call fftw_execute_dft ( plan_k_rs2fs, arr_k_in, arr_k_out )

  arr_k_out = arr_k_out * normalisationfactor

end subroutine f_Fftw_k_Rs2Fs
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_Fftw_FftOnly_k_Fs2Rs
!> Perform an FFT in k direction. 
!> @param ierr should return 0
subroutine f_Fftw_FftOnly_k_Fs2Rs ( ierr )

  use f_arrays,     only : j_width_3w, nprocsz_1r2w, z_allwidths_1r2w, z_maxwidth_1r2w

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr
  integer                                                 :: proc_k, j, k, component, &
                                                             kfft, ktranspose

!> Perform a Fourier-space to real-space FFT on a kj slice in k direction

  call fftw_execute_dft ( plan_k_fs2rs, arr_k_in, arr_k_out )

end subroutine f_Fftw_FftOnly_k_Fs2Rs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_Fftw_FftOnly_k_Rs2Fs
!> Revert the transforms performed in f_Fftw_FftOnly_k_Fs2Rs.
!> @param ierr should return 0
subroutine f_Fftw_FftOnly_k_Rs2Fs ( ierr )

  use f_arrays,     only : j_width_3w, nprocsz_1r2w, z_allwidths_1r2w, z_maxwidth_1r2w
  use g_parameters, only : NODES_Z

  implicit none

  PetscErrorCode,intent(inout)                                   :: ierr
  integer                                                 :: proc_k, j, k, component, &
                                                             kfft, ktranspose

  !> Perform a Fourier-space to real-space FFT on a jk slice in k direction
  call fftw_execute_dft ( plan_k_rs2fs, arr_k_in, arr_k_out )

  !> Normalise the result by nodes in z direction ONLY
  arr_k_out = arr_k_out / ( real ( NODES_Z, kind=C_DOUBLE) )


end subroutine f_Fftw_FftOnly_k_Rs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
end module f_fftw

