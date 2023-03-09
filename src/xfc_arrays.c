#include <petscdm.h>
#include <petscdmda.h>
#include "xfc_arrays.h"

/* Create distributed array with three wave-space dimensions */
PetscErrorCode xfc_ArraysPrepare_da3w(PetscInt nkk, PetscInt nkj, PetscInt nki, DM *da3wout){

  PetscErrorCode ierr;
  DM da3w;

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

  /* Return parallel array */
  *da3wout = da3w;
  return 0;
}


/* Create distributed array with one real-space and two wave-space dimensions */
PetscErrorCode xfc_ArraysPrepare_da1r2w(PetscInt nkj, PetscInt nz, PetscInt nki, PetscInt npi3w, DM *da1r2wout){

  PetscErrorCode ierr;
  DM da1r2w;

  ierr = DMDACreate3d(PETSC_COMM_WORLD,        /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,    /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,        /* need to choose a stencil even if not used */
                      nkj, nz, nki,            /* global dimension in each direction */
                      1, PETSC_DECIDE, npi3w,  /* number of processors in each direction - let PETSc do this */
                      2,                       /* degrees of freedom = array components, we use two for real and imaginary components */
                      0,                       /* stencil width = 0 */
                      NULL, NULL, NULL,        /* arrays containing number of nodes per cell - we don't use this */
                      &da1r2w);
  CHKERRQ(ierr);

  /* Set up parallel array */
  ierr = DMSetUp(da1r2w);
  CHKERRQ(ierr);

  /* Return parallel array */
  *da1r2wout = da1r2w;
  return 0;
}


/* Create distributed array with one real-space and two wave-space dimensions */
PetscErrorCode xfc_ArraysPrepare_da3r(PetscInt nx, PetscInt nz, PetscInt ny, PetscInt npi3w, PetscInt nh, DM *da3rout){

  PetscErrorCode ierr;
  DM da3r;

  ierr = DMDACreate3d(PETSC_COMM_WORLD,        /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,    /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,        /* need to choose a stencil even if not used */
                      nx, nz, ny,              /* global dimension in each direction */
                      1, PETSC_DECIDE, npi3w,  /* number of processors in each direction - let PETSc do this */
                      1,                       /* degrees of freedom = array components, we use two for real and imaginary components */
                      nh,                       /* stencil width = 0 */
                      NULL, NULL, NULL,        /* arrays containing number of nodes per cell - we don't use this */
                      &da3r);
  CHKERRQ(ierr);

  /* Set up parallel array */
  ierr = DMSetUp(da3r);
  CHKERRQ(ierr);

  /* Return parallel array */
  *da3rout = da3r;
  return 0;
}



/* Create distributed array with three wave-space dimensions for reading in from scaled restart file */
PetscErrorCode xfc_ArraysPrepare_da3w_rs(DM da3w, MPI_Comm restartcomm, PetscInt nkk, PetscInt nkj, PetscInt nki, 
                                         PetscInt npk3wrs, PetscInt npj3wrs, PetscInt npi3wrs, DM *da3w_rsout){

  const PetscInt *lk3wrs, *lj3wrs, *li3wrs; 
  DM da3w_rs;
  PetscErrorCode ierr;

  /* Need structure of the final wave-space array to create scaled array */
  ierr = DMDAGetOwnershipRanges(da3w, &lk3wrs, &lj3wrs, &li3wrs);
  CHKERRQ(ierr);

  ierr = DMDACreate3d(restartcomm,                    /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,           /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,               /* need to choose a stencil even if not used */
                      nkk, nkj, nki,                  /* global dimension in each direction */
                      npk3wrs, npj3wrs, npi3wrs,      /* number of processors in each direction - let PETSc do this */
                      2,                              /* degrees of freedom = array components, we use two for real and imaginary components */
                      0,                              /* stencil width = 0 */
                      lk3wrs, lj3wrs, li3wrs,         /* number of nodes per cell */
                      &da3w_rs);
  CHKERRQ(ierr);
  
  /* Set up parallel array */
  ierr = DMSetUp(da3w_rs);
  CHKERRQ(ierr);

  /* Return parallel array */
  *da3w_rsout = da3w_rs;
  return 0;
}


 
/*! Create distributed array with two wave-space dimensions.
    Need this for computing and storing the vorticity components before the first FFT */
PetscErrorCode xfc_ArraysPrepare_da2w(MPI_Comm kjcomm, PetscInt nkk, PetscInt nkj, DM *da2wout){

  PetscErrorCode ierr;
  DM da2w;

  ierr = DMDACreate2d(kjcomm,                  /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,    /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,        /* need to choose a stencil even if not used */
                      nkk, nkj,                /* global dimension in each direction */
                      1, PETSC_DECIDE,         /* number of processors in each direction - let PETSc do this */
                      2,                       /* degrees of freedom = array components, we use two for real and imaginary components */
                      0,                       /* stencil width = 0 */
                      NULL, NULL,              /* arrays containing number of nodes per cell - we don't use this */
                      &da2w);
  CHKERRQ(ierr);

  /* Set up parallel array */
  ierr = DMSetUp(da2w);
  CHKERRQ(ierr);

  /* Return parallel array */
  *da2wout = da2w;
  return 0;
}



/*! Create distributed array with two real-space dimensions
    Need this for computing the non-linear component and for the particle part */
PetscErrorCode xfc_ArraysPrepare_da2r(MPI_Comm jicomm, PetscInt nx, PetscInt ny, DM *da2rout){

  PetscErrorCode ierr;
  DM da2r;

  ierr = DMDACreate2d(jicomm,                  /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,    /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,        /* need to choose a stencil even if not used */
                      nx, ny,                  /* global dimension in each direction */
                      1, PETSC_DECIDE,         /* number of processors in each direction - let PETSc do this */
                      1,                       /* degrees of freedom = array components, we use two for real and imaginary components */
                      0,                       /* stencil width = 0 */
                      NULL, NULL,              /* arrays containing number of nodes per cell - we don't use this */
                      &da2r);
  CHKERRQ(ierr);

  /* Set up parallel array */
  ierr = DMSetUp(da2r);
  CHKERRQ(ierr);

  /* Return parallel array */
  *da2rout = da2r;
  return 0;
}



/*! Create distributed array with one wave-space and one real-space dimension */
PetscErrorCode xfc_ArraysPrepare_da1r1w(MPI_Comm jicomm, PetscInt nki, PetscInt ny, DM *da1r1wout){

  PetscErrorCode ierr;
  DM da1r1w;

  ierr = DMDACreate2d(jicomm,                  /* MPI communicator */
                      DM_BOUNDARY_PERIODIC,    /* periodic boundary conditions */
                      DM_BOUNDARY_PERIODIC, 
                      DMDA_STENCIL_BOX,        /* need to choose a stencil even if not used */
                      nki, ny,                 /* global dimension in each direction */
                      1, PETSC_DECIDE,         /* number of processors in each direction - let PETSc do this */
                      2,                       /* degrees of freedom = array components, we use two for real and imaginary components */
                      0,                       /* stencil width = 0 */
                      NULL, NULL,              /* arrays containing number of nodes per cell - we don't use this */
                      &da1r1w);

  /* Set up parallel array */
  ierr = DMSetUp(da1r1w);
  CHKERRQ(ierr);

  /* Return parallel array */
  *da1r1wout = da1r1w;
  return 0;
}



/* Obtain information needed from da3w */
PetscErrorCode xfc_ArraysGetInfo_da3w(DM da3w, PetscInt *npk3w, PetscInt *npj3w, PetscInt *npi3w,
                                      PetscInt *kmin3w, PetscInt *jmin3w, PetscInt *imin3w, 
                                      PetscInt *kwidth3w, PetscInt *jwidth3w, PetscInt *iwidth3w, 
                                      PetscInt *kmax3w, PetscInt *jmax3w, PetscInt *imax3w){

  PetscInt nproc_k_3w; 
  PetscInt nproc_j_3w; 
  PetscInt nproc_i_3w;
  PetscInt k_min_3w; 
  PetscInt j_min_3w; 
  PetscInt i_min_3w;
  PetscInt k_width_3w; 
  PetscInt j_width_3w; 
  PetscInt i_width_3w;
  PetscInt k_max_3w; 
  PetscInt j_max_3w; 
  PetscInt i_max_3w;
  PetscErrorCode ierr;

  /* Get number of processors from the main wavespace DMDA */
  ierr = DMDAGetInfo(da3w, NULL,                           /* DMDA, dimension of DMDA */
                    NULL, NULL, NULL,                      /* global dimensions M, N, P */
                    &nproc_k_3w, &nproc_j_3w, &nproc_i_3w, /* number of procs m, n, p */
                    NULL, NULL, NULL,                      /* number of dof per node, stencil width, ghost nodes bx */
                    NULL, NULL, NULL);                     /* ghost nodes by, bz, stencil type */

  /* Get local corners */
  ierr = DMDAGetCorners(da3w, &k_min_3w, &j_min_3w, &i_min_3w, &k_width_3w, &j_width_3w, &i_width_3w);

  /* Compute maxima */
  k_max_3w = k_min_3w + k_width_3w - 1;
  j_max_3w = j_min_3w + j_width_3w - 1;
  i_max_3w = i_min_3w + i_width_3w - 1;

  /* Return all values */
  *npk3w = nproc_k_3w;
  *npj3w = nproc_j_3w;
  *npi3w = nproc_i_3w;
  *kmin3w = k_min_3w;
  *jmin3w = j_min_3w;
  *imin3w = i_min_3w;
  *kwidth3w = k_width_3w;
  *jwidth3w = j_width_3w;
  *iwidth3w = i_width_3w;
  *kmax3w = k_max_3w;
  *jmax3w = j_max_3w;
  *imax3w = i_max_3w;
  return 0;
}


/* Obtain information needed from da1r2w */
PetscErrorCode xfc_ArraysGetInfo_da1r2w(DM da1r2w, PetscInt *npj1r2w, PetscInt *npz1r2w, PetscInt *npi1r2w,
                                        PetscInt *jmin1r2w, PetscInt *zmin1r2w, PetscInt *imin1r2w, 
                                        PetscInt *jwidth1r2w, PetscInt *zwidth1r2w, PetscInt *iwidth1r2w, 
                                        PetscInt *jmax1r2w, PetscInt *zmax1r2w, PetscInt *imax1r2w){

  PetscInt nproc_j_1r2w; 
  PetscInt nproc_z_1r2w; 
  PetscInt nproc_i_1r2w;
  PetscInt j_min_1r2w; 
  PetscInt z_min_1r2w; 
  PetscInt i_min_1r2w;
  PetscInt j_width_1r2w; 
  PetscInt z_width_1r2w; 
  PetscInt i_width_1r2w;
  PetscInt j_max_1r2w; 
  PetscInt z_max_1r2w; 
  PetscInt i_max_1r2w;
  PetscErrorCode ierr;

  /* Get number of processors from the main wavespace DMDA */
  ierr = DMDAGetInfo(da1r2w, NULL,                                /* DMDA, dimension of DMDA */
                     NULL, NULL, NULL,                            /* global dimensions M, N, P */
                     &nproc_j_1r2w, &nproc_z_1r2w, &nproc_i_1r2w, /* number of procs m, n, p */
                     NULL, NULL, NULL,                            /* number of dof per node, stencil width, ghost nodes bx */
                     NULL, NULL, NULL);                           /* ghost nodes by, bz, stencil type */

  /* Get local corners */
  ierr = DMDAGetCorners(da1r2w, &j_min_1r2w, &z_min_1r2w, &i_min_1r2w, &j_width_1r2w, &z_width_1r2w, &i_width_1r2w);

  /* Compute maxima */
  j_max_1r2w = j_min_1r2w + j_width_1r2w - 1;
  z_max_1r2w = z_min_1r2w + z_width_1r2w - 1;
  i_max_1r2w = i_min_1r2w + i_width_1r2w - 1;

  /* Return all values */
  *npj1r2w = nproc_j_1r2w;
  *npz1r2w = nproc_z_1r2w;
  *npi1r2w = nproc_i_1r2w;
  *jmin1r2w = j_min_1r2w;
  *zmin1r2w = z_min_1r2w;
  *imin1r2w = i_min_1r2w;
  *jwidth1r2w = j_width_1r2w;
  *zwidth1r2w = z_width_1r2w;
  *iwidth1r2w = i_width_1r2w;
  *jmax1r2w = j_max_1r2w;
  *zmax1r2w = z_max_1r2w;
  *imax1r2w = i_max_1r2w;
  return 0;
}


/* Obtain information needed from da3r */
PetscErrorCode xfc_ArraysGetInfo_da3r(DM da3r, PetscInt *xmin3r, PetscInt *jmin3r, PetscInt *imin3r, 
                                      PetscInt *xwidth3r, PetscInt *jwidth3r, PetscInt *iwidth3r, 
                                      PetscInt *xmax3r, PetscInt *jmax3r, PetscInt *imax3r){

  PetscInt x_min_3r; 
  PetscInt j_min_3r; 
  PetscInt i_min_3r;
  PetscInt x_width_3r; 
  PetscInt j_width_3r; 
  PetscInt i_width_3r;
  PetscInt x_max_3r; 
  PetscInt j_max_3r; 
  PetscInt i_max_3r;
  PetscErrorCode ierr;

  /* Get local corners */
  ierr = DMDAGetCorners(da3r, &x_min_3r, &j_min_3r, &i_min_3r, &x_width_3r, &j_width_3r, &i_width_3r);

  /* Compute maxima */
  x_max_3r = x_min_3r + x_width_3r - 1;
  j_max_3r = j_min_3r + j_width_3r - 1;
  i_max_3r = i_min_3r + i_width_3r - 1;

  /* Return all values */
  *xmin3r = x_min_3r;
  *jmin3r = j_min_3r;
  *imin3r = i_min_3r;
  *xwidth3r = x_width_3r;
  *jwidth3r = j_width_3r;
  *iwidth3r = i_width_3r;
  *xmax3r = x_max_3r;
  *jmax3r = j_max_3r;
  *imax3r = i_max_3r;
  return 0;
}


/*! Obtain information needed from da2r */
PetscErrorCode xfc_ArraysGetInfo_da2r(DM da2r, PetscInt *npx2r, PetscInt *npy2r,
                                        PetscInt *xmin2r, PetscInt *ymin2r, 
                                        PetscInt *xwidth2r, PetscInt *ywidth2r, 
                                        PetscInt *xmax2r, PetscInt *ymax2r){

  PetscInt nproc_x_2r; 
  PetscInt nproc_y_2r; 
  PetscInt x_min_2r; 
  PetscInt y_min_2r; 
  PetscInt x_width_2r; 
  PetscInt y_width_2r; 
  PetscInt x_max_2r; 
  PetscInt y_max_2r; 
  PetscErrorCode ierr;

  /*! Get number of processors from the main wavespace DMDA */
  ierr = DMDAGetInfo(da2r, NULL,                      /* DMDA, dimension of DMDA */
                     NULL, NULL, NULL,                /* global dimensions M, N, P */
                     &nproc_x_2r, &nproc_y_2r, NULL,  /* number of procs m, n, p */
                     NULL, NULL, NULL,                /* number of dof per node, stencil width, ghost nodes bx */
                     NULL, NULL, NULL);               /* ghost nodes by, bz, stencil type */

  /*! Get local corners */
  ierr = DMDAGetCorners(da2r, &x_min_2r, &y_min_2r, NULL, &x_width_2r, &y_width_2r, NULL);

  /*! Compute maxima */
  x_max_2r = x_min_2r + x_width_2r - 1;
  y_max_2r = y_min_2r + y_width_2r - 1;

  /*! Return all values */
  *npx2r = nproc_x_2r;
  *npy2r = nproc_y_2r;
  *xmin2r = x_min_2r;
  *ymin2r = y_min_2r;
  *xwidth2r = x_width_2r;
  *ywidth2r = y_width_2r;
  *xmax2r = x_max_2r;
  *ymax2r = y_max_2r;
  return 0;
}


/*! Obtain information needed from da1r1w */
PetscErrorCode xfc_ArraysGetInfo_da1r1w(DM da1r1w, PetscInt *imin1r1w, PetscInt *ymin1r1w, 
                                        PetscInt *iwidth1r1w, PetscInt *ywidth1r1w, 
                                        PetscInt *imax1r1w, PetscInt *ymax1r1w){

  PetscInt i_min_1r1w; 
  PetscInt y_min_1r1w; 
  PetscInt i_width_1r1w; 
  PetscInt y_width_1r1w; 
  PetscInt i_max_1r1w; 
  PetscInt y_max_1r1w; 
  PetscErrorCode ierr;

  /*! Get local corners */
  ierr = DMDAGetCorners(da1r1w, &i_min_1r1w, &y_min_1r1w, NULL, &i_width_1r1w, &y_width_1r1w, NULL);

  /*! Compute maxima */
  i_max_1r1w = i_min_1r1w + i_width_1r1w - 1;
  y_max_1r1w = y_min_1r1w + y_width_1r1w - 1;

  /*! Return all values */
  *imin1r1w = i_min_1r1w;
  *ymin1r1w = y_min_1r1w;
  *iwidth1r1w = i_width_1r1w;
  *ywidth1r1w = y_width_1r1w;
  *imax1r1w = i_max_1r1w;
  *ymax1r1w = y_max_1r1w;
  return 0;
}


PetscErrorCode xfc_ArraysGetDomainCoordinates(PetscInt zmin, PetscInt zmax,
                                              PetscInt ymin, PetscInt ymax,
                                              PetscInt xmin, PetscInt xmax,
                                              PetscScalar dz, PetscScalar dy, PetscScalar dx,
                                              PetscScalar *zmincoord, PetscScalar *ymincoord, PetscScalar *xmincoord,
                                              PetscScalar *zmaxcoord, PetscScalar *ymaxcoord, PetscScalar *xmaxcoord){

  PetscScalar zmn, ymn, xmn, zmx, ymx, xmx;
  PetscErrorCode ierr;

  /*! Compute domain coordinates */
  zmn = ((( PetscScalar ) zmin) * dz) - (0.5 * dz);
  ymn = ((( PetscScalar ) ymin) * dy) - (0.5 * dy);
  xmn = ((( PetscScalar ) xmin) * dx) - (0.5 * dx);
  zmx = ((( PetscScalar ) zmax) * dz) - (0.5 * dz);
  ymx = ((( PetscScalar ) ymax) * dy) - (0.5 * dy);
  xmx = ((( PetscScalar ) xmax) * dx) - (0.5 * dx);

  /*! Return domain coordinates */
  *zmincoord = zmn;
  *ymincoord = ymn;
  *xmincoord = xmn;
  *zmaxcoord = zmx;
  *ymaxcoord = ymx;
  *xmaxcoord = xmx;
  return 0;
}


/*! Get 3w global velocity vector */
PetscErrorCode xfc_ArraysGetGlobalVelocity_3w(DM da3w, Vec *u1, Vec *u2, Vec *u3){

  Vec u1_3w;
  Vec u2_3w;
  Vec u3_3w;
  PetscErrorCode ierr;
  
  /*! Get global vectors from DMDA */
  ierr = DMGetGlobalVector(da3w, &u1_3w);
  CHKERRQ(ierr);
  ierr = DMGetGlobalVector(da3w, &u2_3w);
  CHKERRQ(ierr);
  ierr = DMGetGlobalVector(da3w, &u3_3w);
  CHKERRQ(ierr);
  /* Initialise all vectors to zero */
  ierr = VecSet(u1_3w, 0.0);
  CHKERRQ(ierr);
  ierr = VecSet(u2_3w, 0.0);
  CHKERRQ(ierr);
  ierr = VecSet(u3_3w, 0.0);
  CHKERRQ(ierr);

  /* Return vectors */
  *u1 = u1_3w;
  *u2 = u2_3w;
  *u3 = u3_3w;
  return 0;
}

/*
PetscErrorCode xfc_ArraysGetGlobalVelocity_1r2w(DM da1r2w, Vec *u1, Vec *u2, Vec *u3);
PetscErrorCode xfc_ArraysGetGlobalVelocity_3r(DM da3r, Vec *u1, Vec *u2, Vec *u3);
PetscErrorCode xfc_ArraysGetGlobalVorticity_1r2w(DM da1r2w, Vec *u4, Vec *u5, Vec *u6); */


/*

  !> u1 to u3: components of the velocity in wave space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for velocity &
                                       in wave space. \n', ierr )
  !> u1 to u3: components of the velocity in 2d wave space 1d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for velocity &
                                       in real space (1d) / wave space (2d). \n', ierr )
  call DMGetGlobalVector ( da1r2w, u1_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u2_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u3_1r2w, ierr )
  CHKERRQ ( ierr )
  call VecSet ( u1_1r2w, ZERO, ierr )
  call VecSet ( u2_1r2w, ZERO, ierr )
  call VecSet ( u3_1r2w, ZERO, ierr )

  !> u4 to u6: components of the vorticity or non-linear term in 2d wave space 1d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for vorticity or non-linear &
                                       component in real space (1d) / wave space (2d). \n', ierr )
  call DMGetGlobalVector ( da1r2w, u4_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u5_1r2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da1r2w, u6_1r2w, ierr )
  CHKERRQ ( ierr )
  call VecSet ( u4_1r2w, ZERO, ierr )
  call VecSet ( u5_1r2w, ZERO, ierr )
  call VecSet ( u6_1r2w, ZERO, ierr )

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
    call VecSet ( u7_1r2w, ZERO, ierr )
    call VecSet ( u8_1r2w, ZERO, ierr )
    call VecSet ( u9_1r2w, ZERO, ierr )
    call VecSet ( u10_1r2w, ZERO, ierr )
    call VecSet ( u11_1r2w, ZERO, ierr )
    call VecSet ( u12_1r2w, ZERO, ierr )
  end if IFCONV1R2W

  !> Runge-Kutta wave-space arrays 
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating Runge-Kutta arrays in wave space. \n', ierr )
  call DMGetGlobalVector ( da3w, rk3_3w_1, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da3w, rk3_3w_2, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da3w, rk3_3w_3, ierr )
  CHKERRQ ( ierr )
    call VecSet ( rk3_3w_1, ZERO, ierr )
    call VecSet ( rk3_3w_2, ZERO, ierr )
    call VecSet ( rk3_3w_3, ZERO, ierr )

 end if IFCONV1R2WZERO

    call DMCreateGlobalVector ( da3r, u1_3r, ierr )
    CHKERRQ ( ierr )
    call DMCreateGlobalVector ( da3r, u2_3r, ierr )
    CHKERRQ ( ierr )
    call DMCreateGlobalVector ( da3r, u3_3r, ierr )
    CHKERRQ ( ierr )

    call VecSet ( u1_3r, ZERO, ierr )
    CHKERRQ ( ierr )
    call VecSet ( u2_3r, ZERO, ierr )
    CHKERRQ ( ierr )
    call VecSet ( u3_3r, ZERO, ierr )
    CHKERRQ ( ierr )



 !> Continuing with the 2D arrays

      call DMGetGlobalVector ( da2w, u1_2w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da2w, u2_2w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da2w, u3_2w, ierr )
      CHKERRQ ( ierr )
  IFTWOWAYINIT: if ( P_TWO_WAY == YES ) then
    call VecSet ( u1_2w, ZERO, ierr )
    call VecSet ( u2_2w, ZERO, ierr )
    call VecSet ( u3_2w, ZERO, ierr )
  end if IFTWOWAYINIT


  !> u4 to u6: components of the vorticity or non-linear term in 2d wave space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for vorticity or non-linear &
                                       component in wave space (2d). \n', ierr )
  call DMGetGlobalVector ( da2w, u4_2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2w, u5_2w, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2w, u6_2w, ierr )
  CHKERRQ ( ierr )
  call VecSet ( u4_2w, ZERO, ierr )
  call VecSet ( u5_2w, ZERO, ierr )
  call VecSet ( u6_2w, ZERO, ierr )

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
    call VecSet ( u7_2w, ZERO, ierr )
    call VecSet ( u8_2w, ZERO, ierr )
    call VecSet ( u9_2w, ZERO, ierr )
    call VecSet ( u10_2w, ZERO, ierr )
    call VecSet ( u11_2w, ZERO, ierr )
    call VecSet ( u12_2w, ZERO, ierr )
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
  call VecSet ( u1_2r, ZERO, ierr )
  call VecSet ( u2_2r, ZERO, ierr )
  call VecSet ( u3_2r, ZERO, ierr )


  !> u4 to u6: components of the vorticity or non-linear term in 2d real space
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating arrays for vorticity or non-linear &
                                       component in real space (2d). \n', ierr )
  call DMGetGlobalVector ( da2r, u4_2r, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2r, u5_2r, ierr )
  CHKERRQ ( ierr )
  call DMGetGlobalVector ( da2r, u6_2r, ierr )
  CHKERRQ ( ierr )
  call VecSet ( u4_2r, ZERO, ierr )
  call VecSet ( u5_2r, ZERO, ierr )
  call VecSet ( u6_2r, ZERO, ierr )
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
    call VecSet ( u7_2r, ZERO, ierr )
    call VecSet ( u8_2r, ZERO, ierr )
    call VecSet ( u9_2r, ZERO, ierr )
    call VecSet ( u10_2r, ZERO, ierr )
    call VecSet ( u11_2r, ZERO, ierr )
    call VecSet ( u12_2r, ZERO, ierr )
  end if IFCONV2R

  !> Arrays for 1D spectra if needed
      call DMGetGlobalVector ( da1r1w, u1_1r1w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da1r1w, u2_1r1w, ierr )
      CHKERRQ ( ierr )
      call DMGetGlobalVector ( da1r1w, u3_1r1w, ierr )
      CHKERRQ ( ierr )
    call VecSet ( u1_1r1w, ZERO, ierr )
    call VecSet ( u2_1r1w, ZERO, ierr )
    call VecSet ( u3_1r1w, ZERO, ierr )

    end if IFSPEC1D
  end if IFISOTROPIC




*/

