module pff_arrays

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>

use g_petsc
use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

  interface DMDA3D
    subroutine pff_ArraysPrepare_da3w ( nkk, nkj, nki, ierr ) bind(C, name='pfc_ArraysPrepare_da3w')
      PetscInt,intent(in)        :: nkk, nkj, nki
!      DM,intent(out)             :: da3w
      PetscErrorCode,intent(out) :: ierr
    end subroutine pff_ArraysPrepare_da3w
  end interface DMDA3D

end module pff_arrays
