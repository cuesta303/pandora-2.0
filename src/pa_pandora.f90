!> @mainpage
!> PANDORA 2.0 is a pseudo-spectral DNS code that includes particle routines.
!! It is capable of simulating homogeneous isotropic turbulence, homogeneous
!! shear flows and other deformed flows that use the Rogallo transform.
!> The main program is pa_pandora.f90, which calls the initialisation, time-stepping and 
!! finalising subroutines provided by the module g_control.
!> Modules beginning with g are general modules, modules beginning with f are for the
!! fluid part and modules beginning with p are for the particle part.
!> Testing modules begin with a t.
!> @author Stephen Scott (author of the original code, 2005-2010)
!> @author Aditya Karnik (co-author of the original code, 2005-2010)
!> @author John Shrimpton 
!> @author Thorsten Wittemeier (since 2016)
!> @author David Scott (2017-2018)
!> @date 2005 - 2018
#include <g_petsc38.h>
!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

program pandora

  use g_petsc
  use g_control

  implicit none

  PetscErrorCode            :: ierr

!> Initialise PETSc and MPI. 
  call PetscInitialize ( PETSC_NULL_CHARACTER, ierr )

  call ProgramInit ( ierr )      !> initialise all modules
  call Calculate ( ierr )        !> run calculation
  call ProgramFinalise ( ierr )  !> finalise all modules

!> Finalise PETSc and MPI.
  call PetscFinalize ( ierr )

end program pandora

