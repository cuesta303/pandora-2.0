#include <petscdm.h>
#include <petscdmda.h>
#include "xfc_arrays.h"

PetscErrorCode pfc_ArraysPrepare_da3w(PetscInt nkk, PetscInt nkj, PetscInt nki){

  PetscErrorCode ierr;
  DM da3w;

  /* Create parallel array */
  ierr = DMDACreate3d(PETSC_COMM_WORLD,               /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,           /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,               /* need to choose a stencil even if not used */
                      nkk, nkj, nki,                  /* global dimension in each direction */
                      1, PETSC_DECIDE, PETSC_DECIDE,  /* number of processors in each direction - let PETSc do this */
                      2,                              /* degrees of freedom = array components, we use two for real and imaginary components */
                      0,                              /* stencil width = 0 */
                      NULL, NULL, NULL,               /* number of nodes per cell - we don't use this */
                      &da3w);
  CHKERRQ(ierr);
  
  /* Set up parallel array */
  ierr = DMSetUp(da3w);
  CHKERRQ(ierr);

  return ierr;
}
