! module g_random
!> Provides routines to handle random numbers.
module g_random

#include <g_petsc38.h>

use iso_c_binding
use g_petsc

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

public :: g_RandomNormal
public :: g_RandomInitFluid
public :: g_RandomInitParticle

contains

!-------------------------------------------------------------------------------
subroutine g_RandomInitFluid ( ierr )

  use g_parameters, only : MYRANK, MASTER, RANDOM_SEED_SIZE, RANDOM_SEED_FLUID, &
                           V_GENERAL, YES 

  implicit none

  PetscErrorCode, intent(inout) :: ierr

  !> Variables for setting random seed
  integer                                             :: seedsize 
  integer,dimension(:),allocatable                    :: seed

  !> Define random seed on each rank. 
  call random_seed ( size=seedsize )
  IFSEEDSIZE: if ( RANDOM_SEED_SIZE <= seedsize ) then
    allocate ( seed(seedsize) )
    seed(:) = 0
    seed(1:RANDOM_SEED_SIZE) = RANDOM_SEED_FLUID(1:RANDOM_SEED_SIZE)
  else
    allocate ( seed(RANDOM_SEED_SIZE) )
    seed(:) = 0
    seed(1:RANDOM_SEED_SIZE) = RANDOM_SEED_FLUID(1:RANDOM_SEED_SIZE)
  end if IFSEEDSIZE

  !> Set the random seed
  call random_seed(put=seed)

  !> Verbose information about random seed for debugging
  IFVERBOSE: if ( V_GENERAL == YES ) then
    write(*,*) ' Fluid random seed on rank ', MYRANK, ': ', seed
  end if IFVERBOSE

end subroutine g_RandomInitFluid
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine g_RandomInitParticle ( ierr )

  use g_parameters, only : MYRANK, MASTER, RANDOM_SEED_SIZE, RANDOM_SEED_PARTICLES, &
                           V_GENERAL, YES 

  implicit none

  PetscErrorCode, intent(inout) :: ierr

  !> Variables for setting random seed
  integer                                             :: seedsize 
  integer,dimension(:),allocatable                    :: seed

  !> Define random seed on each rank. 
  call random_seed ( size=seedsize )
  IFSEEDSIZE: if ( RANDOM_SEED_SIZE <= seedsize ) then
    allocate ( seed(seedsize) )
    seed(:) = 0
    seed(1:RANDOM_SEED_SIZE) = RANDOM_SEED_PARTICLES(1:RANDOM_SEED_SIZE)
  else
    allocate ( seed(RANDOM_SEED_SIZE) )
    seed(:) = 0
    seed(1:RANDOM_SEED_SIZE) = RANDOM_SEED_PARTICLES(1:RANDOM_SEED_SIZE)
  end if IFSEEDSIZE

  !> Set the random seed
  call random_seed(put=seed)

  !> Verbose information about random seed for debugging
  IFVERBOSE: if ( V_GENERAL == YES ) then
    write(*,*) ' Particle random seed on rank ', MYRANK, ': ', seed
  end if IFVERBOSE

end subroutine g_RandomInitParticle
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_RandomNormal
!> Generates pairs of independent random numbers that follow a normal distribution.
!! The implemented algorithm is the Box-Muller method as described e.g. in 
!! Press, William H. FORTRAN Numerical Recipes: Numerical recipes in FORTRAN 90. Cambridge University Press, 1996.
!! See also 
!! G. E. P. Box and Mervin E. Muller, A Note on the Generation of Random Normal Deviates, The Annals of Mathematical Statistics (1958), Vol. 29, No. 2 pp. 610–611,
!! J. Bell: 'Algorithm 334: Normal random deviates', Communications of the ACM, vol. 11, No. 7. 1968,
!! R. Knopp: 'Remark on algorithm 334 [G5]: normal random deviates', Communications of the ACM, vol. 12, No. 5. 1969.
!> @param mean mean of the normal distribution
!> @param variance variance of the normal distribution
!> @param r1 first random number
!> @param r2 second random number
!> @param ierr should return 0
subroutine g_RandomNormal ( mean, variance, r1, r2, ierr )

  use g_constants, only : ONE, TWO

  implicit none

  PetscErrorCode, intent(inout) :: ierr
  real(kind=C_DOUBLE),intent(in)      :: mean
  real(kind=C_DOUBLE),intent(in)      :: variance
  real(kind=C_DOUBLE),intent(out)     :: r1
  real(kind=C_DOUBLE),intent(out)     :: r2
  real(kind=C_DOUBLE)                 :: x1
  real(kind=C_DOUBLE)                 :: x2
  real(kind=C_DOUBLE)                 :: y1
  real(kind=C_DOUBLE)                 :: y2
  real(kind=C_DOUBLE)                 :: w

  r1 = 0.0
  r2 = 0.0
  x1 = 0.0
  x2 = 0.0
  y1 = 0.0
  y2 = 0.0

  w = TWO
    do while ( w >= ONE )
      call random_number(x1)
      call random_number(x2)
      x1 = two * x1 - one
      x2 = two * x2 - one
      w  = x1 * x1 + x2 * x2
    end do
      w  = sqrt ( ( - TWO * log ( w ) ) / w )
      y1 = x1 * w
      y2 = x2 * w

    r1 = mean + sqrt ( variance ) * y1
    r2 = mean + sqrt ( variance ) * y2

end subroutine g_RandomNormal
!---------------------------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!---------------------------------------------------------------------------------
end module
