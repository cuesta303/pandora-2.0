!--------------------------------------------------
! module tg_checkdomain
!> Test all domain parameters determined by g_domain. 
module tg_checkdomain
 
!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

!------- data section begins ----------------------
implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

public :: tg_CheckDomainParameters

!---------- data section ends ----------------------
contains

!---------------------------------------------------------------------------------
! subroutine tg_CheckDomainParameters 
!> Write the domain parameters on all processes to standard out to make sure they
!! are identical.
!> @param ierr should return 0
subroutine tg_CheckDomainParameters ( ierr )

  use g_parameters, only : MASTER, MYRANK, NPROCS
  use g_domain, only: LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME, DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME, &
                      lx_laboratoryframe, ly_laboratoryframe, lz_laboratoryframe, &
                      dx_laboratoryframe, dy_laboratoryframe, dz_laboratoryframe, &
                      bmat, amat, &
                      KI_MIN, KI_MAX, KJ_MIN, KJ_MAX, KK_MIN, KK_MAX

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  integer               :: iproc
 
  namelist /nl_movingframe/     MYRANK, LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME, &
                                DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME
  namelist /nl_laboratoryframe/ MYRANK, lx_laboratoryframe, ly_laboratoryframe, lz_laboratoryframe, &
                                dx_laboratoryframe, dy_laboratoryframe, dz_laboratoryframe
  namelist /nl_transform/       MYRANK, bmat, amat 
  namelist /nl_wavespace/       MYRANK, KI_MIN, KI_MAX, KJ_MIN, KJ_MAX, KK_MIN, KK_MAX

!> Using MPI_Barrier to make the output more likely to be organised by increasing process rank.
!! There is no guarantee for that though. 
  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking domain parameters on all processes... \n', ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking domain parameters in the moving frame... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOMOVING: do iproc=0, NPROCS-1, 1
    IFRANKMOVING: if ( MYRANK == iproc ) then
      write(*,nml=nl_movingframe) 
    end if IFRANKMOVING
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOMOVING

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking domain parameters in the laboratory frame... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOLAB: do iproc=0, NPROCS-1, 1
    IFRANKLAB: if ( MYRANK == iproc ) then
      write(*,nml=nl_laboratoryframe) 
    end if IFRANKLAB
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOLAB

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking transform parameters... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOTRANSFORM: do iproc=0, NPROCS-1, 1
    IFRANKTRANSFORM: if ( MYRANK == iproc ) then
      write(*,nml=nl_transform) 
    end if IFRANKTRANSFORM
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOTRANSFORM

  call PetscPrintf( PETSC_COMM_WORLD, '    | Checking domain parameters in wave space... \n', ierr )

  call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

  DOWAVESPACE: do iproc=0, NPROCS-1, 1
    IFRANKWAVESPACE: if ( MYRANK == iproc ) then
      write(*,nml=nl_wavespace) 
    end if IFRANKWAVESPACE
    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
  end do DOWAVESPACE


end subroutine tg_CheckDomainParameters
!---------------------------------------------------------------------------------

end module tg_checkdomain
