#include <petscdm.h>
#include <petscdmda.h>

/*! Prototypes for functions used to set up the DMDAs */
PetscErrorCode xfc_ArraysPrepare_da3w(PetscInt nkk, PetscInt nkj, PetscInt nki, DM *da3wout);
PetscErrorCode xfc_ArraysPrepare_da1r2w(PetscInt nkj, PetscInt nz, PetscInt nki, PetscInt npi3w, DM *da1r2wout);
PetscErrorCode xfc_ArraysPrepare_da3r(PetscInt nx, PetscInt nz, PetscInt ny, PetscInt npi3w, PetscInt nh, DM *da3rout);
PetscErrorCode xfc_ArraysPrepare_da3w_rs(DM da3w, MPI_Comm restartcomm, PetscInt nkk, PetscInt nkj, PetscInt nki, 
                                         PetscInt npk3wrs, PetscInt npj3wrs, PetscInt npi3wrs, DM *da3w_rsout);
PetscErrorCode xfc_ArraysPrepare_da2w(MPI_Comm kjcomm, PetscInt nkk, PetscInt nkj, DM *da2wout);
PetscErrorCode xfc_ArraysPrepare_da2r(MPI_Comm jicomm, PetscInt nx, PetscInt ny, DM *da2rout);
PetscErrorCode xfc_ArraysPrepare_da1r1w(MPI_Comm jicomm, PetscInt nki, PetscInt ny, DM *da1r1wout);

/*! Prototypes for functions used to obtain necessary info from DMDAs */
PetscErrorCode xfc_ArraysGetInfo_da3w(DM da3w, PetscInt *npk3w, PetscInt *npj3w, PetscInt *npi3w,
                                      PetscInt *kmin3w, PetscInt *jmin3w, PetscInt *imin3w, 
                                      PetscInt *kwidth3w, PetscInt *jwidth3w, PetscInt *iwidth3w, 
                                      PetscInt *kmax3w, PetscInt *jmax3w, PetscInt *imax3w);
PetscErrorCode xfc_ArraysGetInfo_da1r2w(DM da1r2w, PetscInt *npj1r2w, PetscInt *npz1r2w, PetscInt *npi1r2w,
                                        PetscInt *jmin1r2w, PetscInt *zmin1r2w, PetscInt *imin1r2w, 
                                        PetscInt *jwidth1r2w, PetscInt *zwidth1r2w, PetscInt *iwidth1r2w, 
                                        PetscInt *jmax1r2w, PetscInt *zmax1r2w, PetscInt *imax1r2w);
PetscErrorCode xfc_ArraysGetInfo_da3r(DM da3r, PetscInt *xmin3r, PetscInt *jmin3r, PetscInt *imin3r, 
                                      PetscInt *xwidth3r, PetscInt *jwidth3r, PetscInt *iwidth3r, 
                                      PetscInt *xmax3r, PetscInt *jmax3r, PetscInt *imax3r);
PetscErrorCode xfc_ArraysGetInfo_da2r(DM da2r, PetscInt *npx2r, PetscInt *npy2r,
                                        PetscInt *xmin2r, PetscInt *ymin2r, 
                                        PetscInt *xwidth2r, PetscInt *ywidth2r, 
                                        PetscInt *xmax2r, PetscInt *ymax2r);
PetscErrorCode xfc_ArraysGetInfo_da1r1w(DM da1r1w, PetscInt *imin1r1w, PetscInt *ymin1r1w, 
                                        PetscInt *iwidth1r1w, PetscInt *ywidth1r1w, 
                                        PetscInt *imax1r1w, PetscInt *ymax1r1w);
PetscErrorCode xfc_ArraysGetDomainCoordinates(PetscInt zmin, PetscInt zmax,
                                              PetscInt ymin, PetscInt ymax,
                                              PetscInt xmin, PetscInt xmax,
                                              PetscScalar dz, PetscScalar dy, PetscScalar dx,
                                              PetscScalar *zmincoord, PetscScalar *ymincoord, PetscScalar *xmincoord,
                                              PetscScalar *zmaxcoord, PetscScalar *ymaxcoord, PetscScalar *xmaxcoord);

/*! Prototypes for functions used to obtain global vectors from DMDAs */
PetscErrorCode xfc_ArraysGetGlobalVelocity_3w(DM da3w, Vec *u1, Vec *u2, Vec *u3);
PetscErrorCode xfc_ArraysGetGlobalVelocity_1r2w(DM da1r2w, Vec *u1, Vec *u2, Vec *u3);
PetscErrorCode xfc_ArraysGetGlobalVelocity_3r(DM da3r, Vec *u1, Vec *u2, Vec *u3);
PetscErrorCode xfc_ArraysGetGlobalVorticity_1r2w(DM da1r2w, Vec *u4, Vec *u5, Vec *u6);
