!-----------------------------------------------------------------------------
!> The random routines in the original PANDORA code were sequential. Parallel
!! random number routines are faster, but are more likely to contain errors
!! that lead to non-random results. This module contains tests to ensure the 
!! random numbers from the parallel routines are identical to comparable 
!! sequential routines.
module tg_checkrandom

#include <g_petsc38.h>
use g_petsc
!use petscsys

use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>


private

public :: tg_CheckRandomSeed
public :: tg_CheckRandomNormal

contains

!-----------------------------------------------------------------------------
!  subroutine tg_CheckRandomSeed
!> Test if the random routines give identical results on all processes.
!> @param ierr should return 0
  subroutine tg_CheckRandomSeed ( ierr )

  use g_parameters, only : MYRANK, MASTER, NPROCS, RANDOM_SEED_SIZE, &
                           RANDOM_SEED_FLUID, RANDOM_SEED_PARTICLES
  use g_random,     only : g_RandomInitFluid

  implicit none

  PetscErrorCode,intent(inout)                            :: ierr
  real(kind=C_DOUBLE)                              :: x1
  integer                                          :: i, iproc
  integer                                          :: seedsize 
  integer,dimension(:),allocatable                 :: seed

  IFMASTER: if(myrank.eq.master)then
    write(*,100)
!> Find out the minimum size of the random seed and write it to standard out. 
!! The size varies for different compilers.
    call random_seed ( size=seedsize )
    write(*,110) seedsize
  endif IFMASTER

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOALLPROC: do iproc=0, NPROCS-1, 1
    IFRANK: if ( MYRANK == iproc ) then
      write(*,120) MYRANK 
!> Determine the random seed size on each rank.
      call random_seed ( size=seedsize )
      allocate ( seed(seedsize) )
!> Obtain the initial random seed and write it to standard out.
      call random_seed ( get=seed )
      write(*,130)
      write(*,*) seed
!> Obtain a random number and write it to standard out.
      call random_number(x1)
      write(*,140) x1
!> Obtain the random seed after getting the random number.
      call random_seed ( get=seed )
      write(*,150) 
      write(*,*) seed
!> Get 999 other random numbers and write the last one to standard out.
      DOGETRANDOM: do i=2,1000,1
        call random_number(x1)
      end do DOGETRANDOM
      write(*,160) x1
!> Obtain the random seed after 1000 random numbers and write it to standard out.
      call random_seed ( get=seed )
      write(*,170) 
      write(*,*) seed
!> Set random seed on each rank.
      deallocate ( seed )
      write(*,220) MYRANK 
      call g_RandomInitFluid ( ierr )
!      call random_seed ( size=RANDOM_SEED_SIZE )
!      allocate ( seed(RANDOM_SEED_SIZE) )
!      seed(:) = RANDOM_SEED_FLUID(:)
!      call random_seed(put=seed)
!      write(*,230)
!      write(*,*) seed
!> Obtain a random number and write it to standard out.
      call random_number(x1)
      write(*,240) x1
!!> Obtain the random seed after getting the random number.
!      call random_seed ( get=seed )
!      write(*,250) 
!      write(*,*) seed
!> Get 999 other random numbers and write the last one to standard out.
      DOGETRANDOMFLUID: do i=2,1000,1
        call random_number(x1)
      end do DOGETRANDOMFLUID
      write(*,260) x1
!!> Obtain the random seed after 1000 random numbers and write it to standard out.
!      call random_seed ( get=seed )
!      write(*,270) 
!      write(*,*) seed
    end if IFRANK
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  end do DOALLPROC

100 format('    | Checking required random seed size.')
110 format('    | The minimum seed size on this system is ', I4)
120 format('    | Getting initial seed on rank ', I5, '.')
130 format('    | The seed array is: ')
140 format('    | Random number: ', F16.10)
150 format('    | Seed after getting random number: ')
160 format('    | 1000th random number: ', F16.10)
170 format('    | Seed after getting 1000 random numbers: ')
220 format('    | Putting initial fluid seed on rank ', I5, '.')
230 format('    | The seed array is: ')
240 format('    | Random number: ', F16.10)
250 format('    | Seed after getting random number: ')
260 format('    | 1000th random number: ', F16.10)
270 format('    | Seed after getting 1000 random numbers: ')

end subroutine tg_CheckRandomSeed
!----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!  subroutine tg_CheckRandomNormal
!> Test if the random routines give identical results on all processes.
!> Neave, Henry R.: 
!! "On using the Box-Muller transformation with multiplicative congruential pseudo-random number generators." 
!! Applied Statistics (1973): 92-97.
!> @param ierr should return 0
  subroutine tg_CheckRandomNormal ( ierr )

  use g_parameters, only : MYRANK, MASTER, NPROCS, RANDOM_SEED_FLUID, RANDOM_SEED_PARTICLES
  use g_constants, only: ZERO, ONE
  use g_random, only : g_RandomNormal 

  implicit none

  PetscErrorCode,intent(inout)                            :: ierr
  real(kind=C_DOUBLE)                              :: x1, x2
  real(kind=C_DOUBLE),parameter                    :: HIGH = 3.6
  real(kind=C_DOUBLE),parameter                    :: L33 = -3.3, L32 = -3.2, L31 = -3.1, L30 = -3.0, &  
                                                      L29 = -2.9, L28 = -2.8
  integer                                          :: i, iproc
  integer                                          :: seedsize 
  integer,dimension(:),allocatable                 :: seed
  integer                                          :: x1gthigh, x2gthigh
  integer                                          :: x1l33, x1l32, x1l31, x1l30, x1l29, x1l28, &
                                                      x2l33, x2l32, x2l31, x2l30, x2l29, x2l28
  x1l33 = 0
  x1l32 = 0
  x1l31 = 0
  x1l30 = 0
  x1l29 = 0
  x1l28 = 0
  x2l33 = 0
  x2l32 = 0
  x2l31 = 0
  x2l30 = 0
  x2l29 = 0
  x2l28 = 0

  x1gthigh = 0
  x2gthigh = 0

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking accuracy of normally distributed random numbers. \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | for details see Neave, Henry R.:  \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | On using the Box-Muller transformation with multiplicative congruential &
                                             pseudo-random number generators. \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Applied Statistics (1973): 92-97. \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOALLPROC: do iproc=0, NPROCS-1, 1
    IFRANK: if ( MYRANK == iproc ) then
      write(*,110) MYRANK 
!> Determine the random seed size on each rank.
      call random_seed ( size=seedsize )
      allocate ( seed(seedsize) )
!> Obtain the initial random seed and write it to standard out.
      call random_seed ( get=seed )
!> Get 1000000 normally distributed random numbers and write the last one to standard out.
      DOGETRANDOM: do i=1,1000000,1
        call g_RandomNormal ( ZERO, ONE, x1, x2, ierr )

        IFX1L: if ( x1 < L28 ) then
          if ( x1 < L33 ) then
            x1l33 = x1l33 + 1
          else if ( x1 < L32 ) then
            x1l32 = x1l32 + 1
          else if ( x1 < L31 ) then
            x1l31 = x1l31 + 1
          else if ( x1 < L30 ) then
            x1l30 = x1l30 + 1
          else if ( x1 < L29 ) then
            x1l29 = x1l29 + 1
          else if ( x1 < L28 ) then
            x1l28 = x1l28 + 1
          end if
        end if IFX1L
        IFX1G: if ( x1 > HIGH ) then
          x1gthigh = x1gthigh + 1
        end if IFX1G
        IFX2L: if ( x2 < L28 ) then
          if ( x2 < L33 ) then
            x2l33 = x2l33 + 1
          else if ( x2 < L32 ) then
            x2l32 = x2l32 + 1
          else if ( x2 < L31 ) then
            x2l31 = x2l31 + 1
          else if ( x2 < L30 ) then
            x2l30 = x2l30 + 1
          else if ( x2 < L29 ) then
            x2l29 = x2l29 + 1
          else if ( x2 < L28 ) then
            x2l28 = x2l28 + 1
          end if
        end if IFX2L
        IFX2G: if ( x2 > HIGH ) then
          x2gthigh = x2gthigh + 1
        end if IFX2G

      end do DOGETRANDOM
      write(*,120) x1, x2
      write(*,130) x1l33, x2l33
      write(*,140) x1l32, x2l32
      write(*,150) x1l31, x2l31
      write(*,160) x1l30, x2l30
      write(*,170) x1l29, x2l29
      write(*,180) x1l28, x2l28
      write(*,190) x1gthigh, x2gthigh

    end if IFRANK
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOALLPROC

110 format('    | Getting initial seed on rank ', I5, '.')
120 format('    | 1000000th random number pair: ', F16.10, 2X, F16.10)
130 format('    | Numbers in each set below -3.3 (expected 483): ', I7, 2X, I7)
140 format('    | Numbers in each set between -3.3 and -3.2 (expected 204): ', I7, 2X, I7)
150 format('    | Numbers in each set between -3.2 and -3.1 (expected 281): ', I7, 2X, I7)
160 format('    | Numbers in each set between -3.1 and -3.0 (expected 382): ', I7, 2X, I7)
170 format('    | Numbers in each set between -3.0 and -2.9 (expected 516): ', I7, 2X, I7)
180 format('    | Numbers in each set between -2.9 and -2.8 (expected 689): ', I7, 2X, I7)
190 format('    | Numbers in each set above 3.6 (expected 159): ', I7, 2X, I7)

end subroutine tg_CheckRandomNormal
!----------------------------------------------------------------------------

end module tg_checkrandom
