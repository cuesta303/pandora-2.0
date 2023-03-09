!---------------------------------------------------------------------------------
! module g_headers
!> Write information to standard out.
module g_headers

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

public

contains
!---------------------------------------------------------------------------------
! subroutine g_HeadersInit
!> Write program header (first few lines in standard output file)
!> @param ierr should return 0
subroutine g_HeadersInit ( ierr )

  use g_parameters, only : NO, YES

  implicit none

  PetscErrorCode,intent(inout) :: ierr
  integer               :: rank 

  call PetscPrintf( PETSC_COMM_WORLD, ' \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '-------------------------------------------------------- \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' PANDORA 2.0: A pseudo-spectral DNS code with particles  \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' The original PANDORA code was developed by Stephen      \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' Scott and co-authored by Aditya Karnik, with the        \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' support of John Shrimpton and Andreas Kronenburg.       \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' The new PANDORA 2.0 has been further developed within   \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' an eCSE project by Thorsten Wittemeier, David Scott     \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' and John Shrimpton.                                     \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' < include detailed funding and license information >    \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '-------------------------------------------------------- \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, ' \n', ierr )

end subroutine g_HeadersInit

!---------------------------------------------------------------------------------
! subroutine g_HeadersSimulation
!> Writes information to standard out if a simulation is performed.
!> @param ierr should return 0
subroutine g_HeadersSimulation(ierr)

  use g_parameters, only : MASTER, MYRANK

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,140)
    write(*,150)
    write(*,160)
    write(*,170)
  end if IFMASTER

140   FORMAT('--------------------------------------------------------')
150   FORMAT('========>         Initialisation Phase         <========')
160   FORMAT('--------------------------------------------------------')
170   FORMAT('')

end subroutine g_HeadersSimulation

!---------------------------------------------------------------------------------
! subroutine g_HeadersTest
!> Writes information to standard out if the test routines are performed.
!> @param ierr should return 0
subroutine g_HeadersTest ( ierr )

  use g_parameters, only : MASTER, MYRANK

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,180)
    write(*,190)
    write(*,200)
    write(*,210)
  end if IFMASTER

180   FORMAT('--------------------------------------------------------')
190   FORMAT('========>              Test Phase              <========')
200   FORMAT('--------------------------------------------------------')
210   FORMAT('')
 
end subroutine g_HeadersTest

!---------------------------------------------------------------------------------
! subroutine g_HeadersCalc 
!> Writes header for calculation phase
!> @param ierr should return 0
subroutine g_HeadersCalc ( ierr )

  use g_parameters, only : MASTER, MYRANK

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,10)
    write(*,20)
    write(*,30)
    write(*,40)
    write(*,50)
    write(*,60)
  end if IFMASTER

10    format('')
20    format('')
30    format('--------------------------------------------------------')
40    format('========>           Calculation Phase          <========')
50    format('--------------------------------------------------------')
60    format('')

!   end if

end subroutine g_HeadersCalc

!---------------------------------------------------------------------------------
! subroutine g_HeadersTstep 
!> Write header for the current time-step
!> @param timestep number of the current time step 
!> @param simulationtime time in the simulation time frame
!> @param ierr should return 0
subroutine g_HeadersTstep ( timestep, simulationtime, ierr )

  use iso_c_binding
  use g_parameters, only: MYRANK, MASTER

  implicit none

  integer,intent(in)               :: timestep
  real(kind=C_DOUBLE),intent(in)   :: simulationtime
  PetscErrorCode,intent(inout)            :: ierr

  IFMASTER: if ( MYRANK == MASTER )then
    write(*,10)
    write(*,20)
    write(*,30) timestep, simulationtime
    write(*,40)
  end if IFMASTER

10    format('')
20    format('    ----------------------------------------------------')
30    format('      time step = 'i8',   time = 'f15.5)
40    format('    ----------------------------------------------------')

!   end if

end subroutine g_HeadersTstep

!---------------------------------------------------------------------------------
! subroutine g_HeadersFinal
!> Writes headers for the finalise phase.
!> @param ierr should return 0
subroutine g_HeadersFinal ( ierr )

  use g_parameters, only: MYRANK, MASTER

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  IFMASTER: if ( MYRANK == MASTER )then
    write(*,10)
    write(*,20)
    write(*,30)
    write(*,40)
    write(*,50)
    write(*,60)
  end if IFMASTER

10    format('')
20    format('')
30    format('--------------------------------------------------------')
40    format('========>             Finalise Phase           <========')
50    format('--------------------------------------------------------')
60    format('')

end subroutine g_HeadersFinal

!---------------------------------------------------------------------------------

end module g_headers
