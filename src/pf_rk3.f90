!---------------------------------------------------------------------------------
! module f_rk3
!> Contains all subroutines necessary to solve the fluid equations of motion.
module f_rk3

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

private

public :: f_RK3Stage1
public :: f_RK3Stage2
public :: f_RK3Stage3

public :: f_RK3PrepareRealSpace
public :: f_RK3RealSpace
public :: f_RK3TransformToRealSpace
public :: f_RK3TransformAllArraysjiRs2Fs
public :: f_RK3TransformAllArrayskRs2Fs 

contains
!---------------------------------------------------------------------------------
! subroutine f_RK3Stage1 
!> First Runge-Kutta stage for the fluid part
!> @param ierr should return 0
subroutine f_RK3Stage1 ( tstep, time, dt, ierr )

  use f_fluidstats, only : f_FluidstatsRsZero, f_FluidstatsFs
  use g_parameters, only : P_TRACK_PART, YES, P_TS_RELEASE
  use p_stats,      only : p_StatsPrimary, p_StatsRelativeVelocity, p_StatsWrite


  implicit none
 
  integer,intent(in)                    :: tstep
  real(kind=C_DOUBLE),intent(in)        :: time
  real(kind=C_DOUBLE),intent(inout)     :: dt
  PetscErrorCode,intent(inout)                 :: ierr

!> insert timing code here

!> Remesh if shear flow
  call f_RK3Remesh ( tstep, ierr )

!> Do wave space fluid analysis
  call f_FluidstatsFs ( (tstep-1), time, ierr )

  IFPARTRK3: if ( P_TRACK_PART == YES ) then
    IFRELEASEDRK3: if ( tstep > P_TS_RELEASE ) then

      !> Do particle statistics 
!      call p_ControlStats ( tstep, ierr )
!      call PetscPrintf( PETSC_COMM_WORLD, 'Preparing &
!           particle statistics. \n', ierr )
      call p_StatsPrimary ( tstep, dt, ierr )

    end if IFRELEASEDRK3
  end if IFPARTRK3

!> Initialise all quantities for real-space fluid analysis
  call f_FluidstatsRsZero ( ierr )

!> Compute spatial derivatives and vorticity. Do first 1D FFT 
!! and transpose of velocity and vorticity components.
  call f_RK3PrepareRealSpace ( ierr )                          ! form wave space terms      

!> insert timing code here

!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
  call f_RK3RealSpace ( 1, tstep, dt, ierr )

  IFPARTWRITE: if ( P_TRACK_PART == YES ) then
    IFRELEASEDWRITE: if ( tstep > P_TS_RELEASE ) then

      !> Sum up particle relative velocity (computed during RK3)
!      call PetscPrintf( PETSC_COMM_WORLD, 'Finalising &
!           particle statistics. \n', ierr )
      call p_StatsRelativeVelocity ( ierr )

      !> Write particle statistics 
      call p_StatsWrite ( tstep, ierr )

    end if IFRELEASEDWRITE
  end if IFPARTWRITE

!!>  Determine the time step --> Moved to position between non-linear term and particle RK3
!  call f_RK3Tstep ( dt, ierr )

!> insert timing code here

!> if we actually use the spectrum arrays, set them to zero here
!  fluid_spec = zero

!> Solve fluid equation of motion
  call f_RK3FormG1U2 ( dt, tstep, ierr )                         ! form g1 and u2 
    
!> insert timing code here
 
end subroutine f_RK3Stage1

!---------------------------------------------------------------------------------
! subroutine f_RK3Stage2 
!> Second Runge-Kutta stage for the fluid part
!> @param ierr should return 0
subroutine f_RK3Stage2 ( tstep, time, dt, ierr )

!  use f_fluidstats, only: f_FluidstatsRsZero, f_FluidstatsRsCollect, f_FluidstatsFs

  implicit none
 
  integer,intent(in)                    :: tstep
  real(kind=C_DOUBLE),intent(in)        :: time
  real(kind=C_DOUBLE),intent(inout)     :: dt
  PetscErrorCode,intent(inout)                 :: ierr

!> insert timing code here

!!> Do wave space fluid analysis
!  call f_FluidstatsFs ( tstep, time, ierr )

!> Compute spatial derivatives and vorticity. Do first 1D FFT 
!! and transpose of velocity and vorticity components.
  call f_RK3PrepareRealSpace(ierr)                          ! form wave space terms      

!> insert timing code here

!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
  call f_RK3RealSpace ( 2, tstep, dt, ierr )

!> insert timing code here

!> Solve fluid equation of motion
  call f_RK3FormG2U3( dt, tstep, ierr)                         ! form g2 and u3 
    
!> insert timing code here
 
end subroutine f_RK3Stage2

!---------------------------------------------------------------------------------
! subroutine f_RK3Stage3 
!> Third Runge-Kutta stage for the fluid part
!> @param ierr should return 0
subroutine f_RK3Stage3 ( tstep, time, dt, ierr )

  use g_parameters, only: TSTEPS
  use f_fluidstats, only: f_FluidstatsFs 

  implicit none
 
  integer,intent(in)                    :: tstep
  real(kind=C_DOUBLE),intent(in)        :: time
  real(kind=C_DOUBLE),intent(inout)     :: dt
  PetscErrorCode,intent(inout)                 :: ierr

!> insert timing code here

!> Compute spatial derivatives and vorticity. Do first 1D FFT 
!! and transpose of velocity and vorticity components.
  call f_RK3PrepareRealSpace(ierr)                          ! form wave space terms      

!> insert timing code here

  
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
  call f_RK3RealSpace ( 3, tstep, dt, ierr )

!> insert timing code here

!> Solve fluid equation of motion
  call f_RK3FormG3Un ( dt, tstep, ierr )                         ! form g3 and un

!> If last time step do wave space fluid analysis of the final result

  IFLASTTS: if ( tstep == TSTEPS ) then
    call f_FluidstatsFs ( tstep, time, ierr )
  end if IFLASTTS

!> insert timing code here

end subroutine f_RK3Stage3

!---------------------------------------------------------------------------------
! subroutine f_RK3PrepareRealSpace
!> @param ierr should return 0
subroutine f_RK3PrepareRealSpace ( ierr )

  use g_parameters, only : F_TYPE, F_ISOTROPIC, DEALIASING, &
                           DEALIASING_SPHERICAL, DEALIASING_LES, &
                           NS_DEFORMED, NS_R_CONV
  use f_arrays,     only : da3w, da1r2w, da2w, i_min_3w, i_max_3w, &
                           u1_3w, u2_3w, u3_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           u1_1r2w, u2_1r2w, u3_1r2w, &
                           u4_1r2w, u5_1r2w, u6_1r2w, &
                           u7_1r2w, u8_1r2w, u9_1r2w, &
                           u10_1r2w, u11_1r2w, u12_1r2w, &
                           arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                           arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                           arr_u7_1r2w, arr_u8_1r2w, arr_u9_1r2w, &
                           arr_u10_1r2w, arr_u11_1r2w, arr_u12_1r2w, &
                           u4_2w, u5_2w, u6_2w, &
                           u7_2w, u8_2w, u9_2w, &
                           u10_2w, u11_2w, u12_2w, &
                           arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                           arr_u7_2w, arr_u8_2w, arr_u9_2w, &
                           arr_u10_2w, arr_u11_2w, arr_u12_2w, &
                           f_ArraysDealiasSpherical, f_ArraysLESCutoff

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: islice                 

  PetscErrorCode        :: perr

!> Get read access to all velocity components
  call DMDAVecGetArrayReadF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayReadF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayReadF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )

!> Get write access to all vorticity components
  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  IFCONVGET: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMDAVecGetArrayF90 ( da2w, u7_2w, arr_u7_2w, perr )
    call DMDAVecGetArrayF90 ( da2w, u8_2w, arr_u8_2w, perr )
    call DMDAVecGetArrayF90 ( da2w, u9_2w, arr_u9_2w, perr )
    call DMDAVecGetArrayF90 ( da2w, u10_2w, arr_u10_2w, perr )
    call DMDAVecGetArrayF90 ( da2w, u11_2w, arr_u11_2w, perr )
    call DMDAVecGetArrayF90 ( da2w, u12_2w, arr_u12_2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u7_1r2w, arr_u7_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u8_1r2w, arr_u8_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u9_1r2w, arr_u9_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u10_1r2w, arr_u10_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u11_1r2w, arr_u11_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u12_1r2w, arr_u12_1r2w, perr )
  end if IFCONVGET

  IFSPHER: if ( DEALIASING == DEALIASING_SPHERICAL ) then
  !> Truncate all wavemodes outside a 2/3 sphere
    call f_ArraysDealiasSpherical ( ierr )
  end if IFSPHER

  IFLES: if ( DEALIASING == DEALIASING_LES ) then
  !> Truncate all wavemodes outside cutoff
    call f_ArraysLESCutoff ( ierr )
  end if IFLES  

!Loop going through all planes
  DOSLICES: do islice = i_min_3w, i_max_3w, 1

    !> Transform the velocity components in k direction and compute vorticity while at it
    call f_RK3PrepareVelocity ( islice, ierr )
    !> Transform the vorticity components in k direction
    call f_RK3PrepareVorticity ( islice, ierr )

  end do DOSLICES

!  write(*,*) ' arr_u1_1r2w(:,:,:,0): ', arr_u1_1r2w(:,:,:,0)
!  write(*,*) ' arr_u2_1r2w(:,:,:,0): ', arr_u2_1r2w(:,:,:,0)
!  write(*,*) ' arr_u3_1r2w(:,:,:,0): ', arr_u3_1r2w(:,:,:,0)
!  write(*,*) ' arr_u4_1r2w(:,:,:,0): ', arr_u4_1r2w(:,:,:,0)
!  write(*,*) ' arr_u5_1r2w(:,:,:,0): ', arr_u5_1r2w(:,:,:,0)
!  write(*,*) ' arr_u6_1r2w(:,:,:,0): ', arr_u6_1r2w(:,:,:,0)


!> Return velocity components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u3_3w, arr_u3_3w, perr )

  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )

!> Return vorticity components to PETSc
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  IFCONVRESTORE: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMDAVecRestoreArrayF90 ( da2w, u7_2w, arr_u7_2w, perr )
    call DMDAVecRestoreArrayF90 ( da2w, u8_2w, arr_u8_2w, perr )
    call DMDAVecRestoreArrayF90 ( da2w, u9_2w, arr_u9_2w, perr )
    call DMDAVecRestoreArrayF90 ( da2w, u10_2w, arr_u10_2w, perr )
    call DMDAVecRestoreArrayF90 ( da2w, u11_2w, arr_u11_2w, perr )
    call DMDAVecRestoreArrayF90 ( da2w, u12_2w, arr_u12_2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u7_1r2w, arr_u7_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u8_1r2w, arr_u8_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u9_1r2w, arr_u9_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u10_1r2w, arr_u10_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u11_1r2w, arr_u11_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u12_1r2w, arr_u12_1r2w, perr )
  end if IFCONVRESTORE

end subroutine f_RK3PrepareRealSpace
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3PrepareVelocity
!> @param ierr should return 0
subroutine f_RK3PrepareVelocity ( islice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, bmat
  use g_parameters, only: NODES_Z, F_TYPE, F_ISOTROPIC, &
                          NS_DEFORMED, NS_BF_ROT, NS_R_CONV, NS_R_ROT
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w, k_max_3w, z_max, j_max_1r2w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          arr_u7_2w, arr_u8_2w, arr_u9_2w, &
                          arr_u10_2w, arr_u11_2w, arr_u12_2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w, &
                          f_ArraysWavenumber, f_ArraysWavenumberPositive, f_ArraysWavenumberNegative
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_k_in_da, arr_kjjk_transposeout

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: islice

  integer                                :: slicelocal
  integer                                :: j, jj, k, kk, proc_j, jtranspose, ktranspose, j1r2w, component
  real(kind=C_DOUBLE)                    :: wavenumber_k, wavenumber_j, wavenumber_i
  real(kind=C_DOUBLE),dimension(0:1)     :: velocity, dudx, dudy, dudz, &
                                            dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: vorticity

  type(C_PTR)                            :: copytoarray

!> wave number to be computed from array indices

! Initialise vorticity buffer arrays to zero.
! Go through all velocity components and gradually compute the vorticity.

! Copy each velocity component to the zero-padded velocity buffer array.
  !> Numbering of buffer array is local, numbering of velocity array is global
  !> PETSc uses C indexing, FFTW Fortran pointers use Fortran indexing!! 

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = islice - i_min_3w + 1

  !> Compute i component of wave vector
  wavenumber_i = f_ArraysWavenumber ( islice, KI_MAX, K0I )

  !> u1 component
  jj = 1
  DOU1J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Compute j component of wave vector
    wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )

    !> Copy positive wavenumbers
    DOU1KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u1_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      !> Compute contribution to vorticity
      !> Compute k component of wave vector
      wavenumber_k = f_ArraysWavenumberPositive ( k, K0K )

      IFU1POSITIVE: if ( F_TYPE == F_ISOTROPIC ) then
        call f_RK3Derivative ( velocity, wavenumber_j, dudy )
        call f_RK3Derivative ( velocity, wavenumber_k, dudz )
        arr_u4_2w(:,k,j) = ZERO
        arr_u5_2w(:,k,j) = dudz
        arr_u6_2w(:,k,j) = - dudy
      else
        select case ( NS_DEFORMED )
        !-------------------------------
          case ( NS_BF_ROT ) 
            call f_RK3VorticityBFRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 1, vorticity )
            arr_u4_2w(:,k,j) = vorticity(:,1)
            arr_u5_2w(:,k,j) = vorticity(:,2)
            arr_u6_2w(:,k,j) = vorticity(:,3)
          case ( NS_R_ROT ) 
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 2, dudy )
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 3, dudz )
            arr_u4_2w(:,k,j) = ZERO
            arr_u5_2w(:,k,j) = dudz
            arr_u6_2w(:,k,j) = - dudy
          case ( NS_R_CONV ) 
            call f_RK3Derivative ( velocity, wavenumber_i, dudx )
            call f_RK3Derivative ( velocity, wavenumber_j, dudy )
            call f_RK3Derivative ( velocity, wavenumber_k, dudz )
            arr_u4_2w(:,k,j) = dudx
            arr_u7_2w(:,k,j) = dudy
            arr_u10_2w(:,k,j) = dudz
        end select
      end if IFU1POSITIVE

!      write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)

      kk = kk + 1
    end do DOU1KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU1KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU1KZEROS

    !> Copy negative wavenumbers
    DOU1KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u1_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      !> Compute contribution to vorticity
      !> Compute k component of wave vector
      wavenumber_k = f_ArraysWavenumberNegative ( k, KK_MAX, K0K )

      IFU1NEGATIVE: if ( F_TYPE == F_ISOTROPIC ) then
        call f_RK3Derivative ( velocity, wavenumber_j, dudy )
        call f_RK3Derivative ( velocity, wavenumber_k, dudz )
        arr_u4_2w(:,k,j) = ZERO
        arr_u5_2w(:,k,j) = dudz
        arr_u6_2w(:,k,j) = - dudy
      else
        select case ( NS_DEFORMED )
        !-------------------------------
          case ( NS_BF_ROT ) 
            call f_RK3VorticityBFRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 1, vorticity )
            arr_u4_2w(:,k,j) = vorticity(:,1)
            arr_u5_2w(:,k,j) = vorticity(:,2)
            arr_u6_2w(:,k,j) = vorticity(:,3)
          case ( NS_R_ROT ) 
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 2, dudy)
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 3, dudz)
            arr_u4_2w(:,k,j) = ZERO
            arr_u5_2w(:,k,j) = dudz
            arr_u6_2w(:,k,j) = - dudy
          case ( NS_R_CONV ) 
            call f_RK3Derivative ( velocity, wavenumber_i, dudx )
            call f_RK3Derivative ( velocity, wavenumber_j, dudy )
            call f_RK3Derivative ( velocity, wavenumber_k, dudz )
            arr_u4_2w(:,k,j) = dudx
            arr_u7_2w(:,k,j) = dudy
            arr_u10_2w(:,k,j) = dudz
        end select
      end if IFU1NEGATIVE

!      write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)

      kk = kk + 1
    end do DOU1KNEGATIVE

    jj = jj + 1
  end do DOU1J3W

  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  copytoarray = c_loc(arr_u1_1r2w)
  call f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )

  !> u2 component
  jj = 1
  DOU2J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Compute j component of wave vector only for deformed flow
    IFFTYPE: if ( F_TYPE /= F_ISOTROPIC ) then
      wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
    end if IFFTYPE

    !> Copy positive wavenumbers
    DOU2KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u2_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      !> Compute contribution to vorticity
      !> Compute k component of wave vector
      wavenumber_k = f_ArraysWavenumberPositive ( k, K0K )

      IFU2POSITIVE: if ( F_TYPE == F_ISOTROPIC ) then
        call f_RK3Derivative ( velocity, wavenumber_i, dvdx )
        call f_RK3Derivative ( velocity, wavenumber_k, dvdz )
        arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) - dvdz
        arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + dvdx
      else
        select case ( NS_DEFORMED )
        !-------------------------------
          case ( NS_BF_ROT ) 
            call f_RK3VorticityBFRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 2, vorticity )
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + vorticity(:,1)
            arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) + vorticity(:,2)
            arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + vorticity(:,3)
          case ( NS_R_ROT ) 
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 1, dvdx)
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 3, dvdz)
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) - dvdz
            arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + dvdx
          case ( NS_R_CONV ) 
            call f_RK3Derivative ( velocity, wavenumber_i, dvdx )
            call f_RK3Derivative ( velocity, wavenumber_j, dvdy )
            call f_RK3Derivative ( velocity, wavenumber_k, dvdz )
            arr_u5_2w(:,k,j) = dvdx
            arr_u8_2w(:,k,j) = dvdy
            arr_u11_2w(:,k,j) = dvdz
        end select
      end if IFU2POSITIVE

!      write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)

      kk = kk + 1
    end do DOU2KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU2KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU2KZEROS

    !> Copy negative wavenumbers
    DOU2KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u2_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      !> Compute contribution to vorticity
      !> Compute k component of wave vector
      wavenumber_k = f_ArraysWavenumberNegative ( k, KK_MAX, K0K )

      IFU2NEGATIVE: if ( F_TYPE == F_ISOTROPIC ) then
        call f_RK3Derivative ( velocity, wavenumber_i, dvdx )
        call f_RK3Derivative ( velocity, wavenumber_k, dvdz )
        arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) - dvdz
        arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + dvdx
      else
        select case ( NS_DEFORMED )
        !-------------------------------
          case ( NS_BF_ROT ) 
            call f_RK3VorticityBFRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 2, vorticity )
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + vorticity(:,1)
            arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) + vorticity(:,2)
            arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + vorticity(:,3)
          case ( NS_R_ROT ) 
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 1, dvdx)
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 3, dvdz)
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) - dvdz
            arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + dvdx
          case ( NS_R_CONV ) 
            call f_RK3Derivative ( velocity, wavenumber_i, dvdx )
            call f_RK3Derivative ( velocity, wavenumber_j, dvdy )
            call f_RK3Derivative ( velocity, wavenumber_k, dvdz )
            arr_u5_2w(:,k,j) = dvdx
            arr_u8_2w(:,k,j) = dvdy
            arr_u11_2w(:,k,j) = dvdz
        end select
      end if IFU2NEGATIVE

!      write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)

      kk = kk + 1
    end do DOU2KNEGATIVE

  jj = jj + 1
  end do DOU2J3W

  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  copytoarray = c_loc(arr_u2_1r2w)
  call f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )


  !> u3 component
  jj = 1
  DOU3J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Compute j component of wave vector
    wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )

    !> Copy positive wavenumbers
    DOU3KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u3_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      !> Compute contribution to vorticity
      IFU3POSITIVE: if ( F_TYPE == F_ISOTROPIC ) then
        call f_RK3Derivative ( velocity, wavenumber_i, dwdx )
        call f_RK3Derivative ( velocity, wavenumber_j, dwdy )
        arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + dwdy
        arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) - dwdx
      else
        !> Compute k component of wave vector
        wavenumber_k = f_ArraysWavenumberPositive ( k, K0K )
        select case ( NS_DEFORMED )
        !-------------------------------
          case ( NS_BF_ROT ) 
            call f_RK3VorticityBFRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 3, vorticity )
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + vorticity(:,1)
            arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) + vorticity(:,2)
            arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + vorticity(:,3)
          case ( NS_R_ROT ) 
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 1, dwdx)
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 2, dwdy)
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + dwdy
            arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) - dwdx
          case ( NS_R_CONV ) 
            call f_RK3Derivative ( velocity, wavenumber_i, dwdx )
            call f_RK3Derivative ( velocity, wavenumber_j, dwdy )
            call f_RK3Derivative ( velocity, wavenumber_k, dwdz )
            arr_u6_2w(:,k,j) = dwdx
            arr_u9_2w(:,k,j) = dwdy
            arr_u12_2w(:,k,j) = dwdz
        end select
      end if IFU3POSITIVE

!      write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)

      kk = kk + 1
    end do DOU3KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU3KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU3KZEROS

    !> Copy negative wavenumbers
    DOU3KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u3_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      IFU3NEGATIVE: if ( F_TYPE == F_ISOTROPIC ) then
        call f_RK3Derivative ( velocity, wavenumber_i, dwdx )
        call f_RK3Derivative ( velocity, wavenumber_j, dwdy )
        arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + dwdy
        arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) - dwdx
      else
        !> Compute k component of wave vector
        wavenumber_k = f_ArraysWavenumberNegative ( k, KK_MAX, K0K )
        select case ( NS_DEFORMED )
        !-------------------------------
          case ( NS_BF_ROT ) 
            call f_RK3VorticityBFRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 3, vorticity )
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + vorticity(:,1)
            arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) + vorticity(:,2)
            arr_u6_2w(:,k,j) = arr_u6_2w(:,k,j) + vorticity(:,3)
          case ( NS_R_ROT ) 
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 1, dwdx)
            call f_RK3DerivativeRRot ( velocity, wavenumber_i, wavenumber_j, wavenumber_k, 2, dwdy)
            arr_u4_2w(:,k,j) = arr_u4_2w(:,k,j) + dwdy
            arr_u5_2w(:,k,j) = arr_u5_2w(:,k,j) - dwdx
          case ( NS_R_CONV ) 
            call f_RK3Derivative ( velocity, wavenumber_i, dwdx )
            call f_RK3Derivative ( velocity, wavenumber_j, dwdy )
            call f_RK3Derivative ( velocity, wavenumber_k, dwdz )
            arr_u6_2w(:,k,j) = dwdx
            arr_u9_2w(:,k,j) = dwdy
            arr_u12_2w(:,k,j) = dwdz
        end select
      end if IFU3NEGATIVE

!      write(*,*) islice, j, k, arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), arr_u6_2w(:,k,j)

      kk = kk + 1
    end do DOU3KNEGATIVE

    jj = jj + 1
  end do DOU3J3W

  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  copytoarray = c_loc(arr_u3_1r2w)
  call f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )

end subroutine f_RK3PrepareVelocity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3TransformVelocitykFs2Rs
!> @param ierr should return 0
  subroutine f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )

  use f_arrays,     only: i_min_3w, z_width_1r2w, j_width_1r2w, i_width_1r2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_kjjk_transposeout

  implicit none
 
  type(C_PTR), intent(in)                :: copytoarray
  integer,intent(in)                     :: slicelocal
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                          :: j, k, proc_j, jtranspose, &
                                                      ktranspose, j1r2w
  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:)   :: targetarray

  !> Perform FFT in k direction
  call f_Fftw_k_Fs2Rs ( ierr )

  !> Link targetarray to the C pointer of the velocity component
  call c_f_pointer(copytoarray, targetarray, [2,j_width_1r2w,z_width_1r2w,i_width_1r2w])

! Copy the result of the FFT to the buffer array
  DOZ: do k = 1, z_width_1r2w, 1                     ! parallelised in this direction
    j1r2w = 0
    DOPROCJ: do proc_j = 0, nprocsj_3w - 1, 1 
      DOJLOCAL: do j = 1, j_allwidths_3w(proc_j), 1

        j1r2w = j1r2w + 1
        jtranspose = (proc_j * j_maxwidth_3w) + j
        targetarray(:, j1r2w, k, slicelocal) = arr_kjjk_transposeout(:, jtranspose, k)

      end do DOJLOCAL
    end do DOPROCJ
  end do DOZ

end subroutine f_RK3TransformVelocitykFs2Rs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3PrepareVorticity
!> @param ierr should return 0
subroutine f_RK3PrepareVorticity ( islice, ierr )

  use g_constants,  only: ZERO, ONE
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, bmat
  use g_parameters, only: NODES_Z, F_TYPE, F_ISOTROPIC, &
                          NS_DEFORMED, NS_R_CONV
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w, k_max_3w, z_max, j_max_1r2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          arr_u7_2w, arr_u8_2w, arr_u9_2w, &
                          arr_u10_2w, arr_u11_2w, arr_u12_2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          arr_u7_1r2w, arr_u8_1r2w, arr_u9_1r2w, &
                          arr_u10_1r2w, arr_u11_1r2w, arr_u12_1r2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: islice

  integer                                :: slicelocal
  integer                                :: j, jj, k, kk, proc_j, jtranspose, ktranspose, j1r2w, component

  type(C_PTR)                            :: inputarray, copytoarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = islice - i_min_3w + 1

  !> u4 component
  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  inputarray = c_loc(arr_u4_2w)
  copytoarray = c_loc(arr_u4_1r2w)
  call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  !> u5 component
  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  inputarray = c_loc(arr_u5_2w)
  copytoarray = c_loc(arr_u5_1r2w)
  call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )


  !> u6 component
    !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  inputarray = c_loc(arr_u6_2w)
  copytoarray = c_loc(arr_u6_1r2w)
  call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )
  
  !> Additional transforms for convective formulation
  IFCONV: if ( NS_DEFORMED == NS_R_CONV ) then

    !> u7 component
    !> Do FFT in k direction and transpose from kj to jk
    !! Copy the result of the FFT to the buffer array
    inputarray = c_loc(arr_u7_2w)
    copytoarray = c_loc(arr_u7_1r2w)
    call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u8 component
    !> Do FFT in k direction and transpose from kj to jk
    !! Copy the result of the FFT to the buffer array
    inputarray = c_loc(arr_u8_2w)
    copytoarray = c_loc(arr_u8_1r2w)
    call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u9 component
    !> Do FFT in k direction and transpose from kj to jk
    !! Copy the result of the FFT to the buffer array
    inputarray = c_loc(arr_u9_2w)
    copytoarray = c_loc(arr_u9_1r2w)
    call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u10 component
    !> Do FFT in k direction and transpose from kj to jk
    !! Copy the result of the FFT to the buffer array
    inputarray = c_loc(arr_u10_2w)
    copytoarray = c_loc(arr_u10_1r2w)
    call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )
  
    !> u11 component
    !> Do FFT in k direction and transpose from kj to jk
    !! Copy the result of the FFT to the buffer array
    inputarray = c_loc(arr_u11_2w)
    copytoarray = c_loc(arr_u11_1r2w)
    call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u12 component
    !> Do FFT in k direction and transpose from kj to jk
    !! Copy the result of the FFT to the buffer array
    inputarray = c_loc(arr_u12_2w)
    copytoarray = c_loc(arr_u12_1r2w)
    call f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  end if IFCONV

end subroutine f_RK3PrepareVorticity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3TransformVorticitykFs2Rs
!> @param ierr should return 0
  subroutine f_RK3TransformVorticitykFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  use g_domain,     only: KK_MAX
  use g_parameters, only: NODES_Z
  use g_constants,  only: ZERO
  use f_arrays,     only: k_width_3w, j_width_3w, &
                          z_width_1r2w, j_width_1r2w, i_width_1r2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_k_in_da, arr_kjjk_transposeout

  implicit none
 
  type(C_PTR), intent(in)                :: inputarray, copytoarray
  integer,intent(in)                     :: slicelocal
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                          :: j, k, kk, proc_j, jtranspose, &
                                                      ktranspose, j1r2w
  real(kind=C_DOUBLE),pointer,dimension(:,:,:)     :: fromarray
  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:)   :: targetarray

  !> Link inputarray to the C pointer of the vorticity component
  !> This array uses Fortran indexing, which means all indices are one higher
  call c_f_pointer(inputarray, fromarray, [2, k_width_3w, j_width_3w])

  DOJ3W: do j = 1, j_width_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOKPOSITIVE: do k = 1, KK_MAX + 1, 1
      !> Copy to buffer
      arr_k_in_da(:,kk,j) = fromarray(:,k,j)
      kk = kk + 1
    end do DOKPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOKZEROS: do k = KK_MAX + 2, NODES_Z - KK_MAX 
      arr_k_in_da(:,kk,j) = ZERO 
      kk = kk + 1
    end do DOKZEROS

    !> Copy negative wavenumbers
    DOKNEGATIVE: do k = KK_MAX + 2, k_width_3w
      !> Copy to buffer
      arr_k_in_da(:,kk,j) = fromarray(:,k,j)
      kk = kk + 1
    end do DOKNEGATIVE

  end do DOJ3W

  !> Perform FFT in k direction and transpose from kj to jk
  call f_Fftw_k_Fs2Rs ( ierr )

  !> Link targetarray to the C pointer of the vorticity component
  call c_f_pointer(copytoarray, targetarray, [2,j_width_1r2w,z_width_1r2w,i_width_1r2w])

! Copy the result of the FFT to the buffer array
  DOZ: do k = 1, z_width_1r2w, 1                     ! parallelised in this direction
    j1r2w = 0
    DOPROCJ: do proc_j = 0, nprocsj_3w - 1, 1 
      DOJLOCAL: do j = 1, j_allwidths_3w(proc_j), 1

        j1r2w = j1r2w + 1
        jtranspose = (proc_j * j_maxwidth_3w) + j
        targetarray(:, j1r2w, k, slicelocal) = arr_kjjk_transposeout(:, jtranspose, k)

      end do DOJLOCAL
    end do DOPROCJ
  end do DOZ

end subroutine f_RK3TransformVorticitykFs2Rs
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RealSpace 
!> Transform velocity and vorticity into real space and compute the non-linear term.
!! Transform the non-linear term back into wave space.
!! If there are particles, compute the particle motion while in real space.
!> @param stage stage of the Runge-Kutta scheme
!> @param tstep current timestep
!> @param ierr should return 0
subroutine f_RK3RealSpace ( stage, tstep, deltat, ierr )

  use g_constants,  only: ONE, TWO
  use g_parameters, only: P_TRACK_PART, P_TWO_WAY, YES, P_TS_RELEASE, &
                          F_TYPE, F_ISOTROPIC, MYRANK, &
                          NS_DEFORMED, NS_BF_ROT, NS_R_CONV, NS_R_ROT, &
                          TS_SPEC1D, TSTEPS
  use g_domain,     only: bmat
  use f_arrays, only:     da1r2w, da2r, da3r, z_min, z_max, &
                          u1_1r2w, u2_1r2w, u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, &
                          u7_1r2w, u8_1r2w, u9_1r2w, &
                          u10_1r2w, u11_1r2w, u12_1r2w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          arr_u7_1r2w, arr_u8_1r2w, arr_u9_1r2w, &
                          arr_u10_1r2w, arr_u11_1r2w, arr_u12_1r2w, &
                          u1_2r, u2_2r, u3_2r, &
                          u4_2r, u5_2r, u6_2r, &
                          u7_2r, u8_2r, u9_2r, &
                          u10_2r, u11_2r, u12_2r, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          arr_u7_2r, arr_u8_2r, arr_u9_2r, &
                          arr_u10_2r, arr_u11_2r, arr_u12_2r, &
                          u1_3r, u2_3r, u3_3r, &
!                          u1_3r_local, u2_3r_local, u3_3r_local, &
                          arr_u1_3r, arr_u2_3r, arr_u3_3r, &
                          x_min_2r, x_max_2r, y_min_2r, y_max_2r, &
                          da1r1w, u1_1r1w, u2_1r1w, u3_1r1w, &
                          arr_u1_1r1w, arr_u2_1r1w, arr_u3_1r1w
  use f_fluidstats, only: f_FluidstatsFs2D, f_FluidstatsRs2D, &
                          f_FluidstatsRsCollect, f_FluidstatsFsCollect
  use p_control,    only: p_ControlRK3, p_ControlFinishTimeStep

!p_ControlFreeGhostParticles, &
                          
  implicit none
 
  integer,intent(in)                :: stage
  integer,intent(in)                :: tstep
  real(kind=C_DOUBLE),intent(inout) :: deltat
  PetscErrorCode,intent(inout)             :: ierr

  integer                           :: zslice                 
  integer                           :: yval

  PetscErrorCode                    :: perr

  !> Get write access to all velocity and vorticity components.
  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

  IFGETCONV: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMDAVecGetArrayF90 ( da1r2w, u7_1r2w, arr_u7_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u8_1r2w, arr_u8_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u9_1r2w, arr_u9_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u10_1r2w, arr_u10_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u11_1r2w, arr_u11_1r2w, perr )
    call DMDAVecGetArrayF90 ( da1r2w, u12_1r2w, arr_u12_1r2w, perr )
    call DMDAVecGetArrayF90 ( da2r, u7_2r, arr_u7_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u8_2r, arr_u8_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u9_2r, arr_u9_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u10_2r, arr_u10_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u11_2r, arr_u11_2r, perr )
    call DMDAVecGetArrayF90 ( da2r, u12_2r, arr_u12_2r, perr )
  end if IFGETCONV

  IFGETARRAYPART: if ( P_TRACK_PART == YES ) then
    IFGETARRAYRELEASED: if ( tstep >= P_TS_RELEASE ) then

!      call PetscPrintf( PETSC_COMM_WORLD, 'Getting arrays for &
!           velocity interpolation \n', ierr )
      call DMDAVecGetArrayF90 ( da3r, u1_3r, arr_u1_3r, perr )
      call DMDAVecGetArrayF90 ( da3r, u2_3r, arr_u2_3r, perr )
      call DMDAVecGetArrayF90 ( da3r, u3_3r, arr_u3_3r, perr )

    end if IFGETARRAYRELEASED
  end if IFGETARRAYPART

  !> Arrays for 1D spectra only
  IFISOTROPICGETARRAY: if ( F_TYPE /= F_ISOTROPIC ) then

    call DMDAVecGetArrayF90 ( da1r1w, u1_1r1w, arr_u1_1r1w, perr )
    call DMDAVecGetArrayF90 ( da1r1w, u2_1r1w, arr_u2_1r1w, perr )
    call DMDAVecGetArrayF90 ( da1r1w, u3_1r1w, arr_u3_1r1w, perr )

  else
    IFSPEC1DGETARRAY: if ( TS_SPEC1D <= TSTEPS ) then

      call DMDAVecGetArrayF90 ( da1r1w, u1_1r1w, arr_u1_1r1w, perr )
      call DMDAVecGetArrayF90 ( da1r1w, u2_1r1w, arr_u2_1r1w, perr )
      call DMDAVecGetArrayF90 ( da1r1w, u3_1r1w, arr_u3_1r1w, perr )

    end if IFSPEC1DGETARRAY
  end if IFISOTROPICGETARRAY

  !> Loop through slices in z direction.
  DOSLICES: do zslice = z_min, z_max, 1
               
    !> Transform velocity and vorticity components to real space.
    call f_RK3TransformToRealSpace ( zslice, stage, (tstep-1), ierr )

    !> compute real-space fluid statistics
    IFRSTATS: if ( stage == 1 ) then
      call f_FluidstatsFs2D ( (tstep-1), ierr )
      call f_FluidstatsRs2D ( (tstep-1), ierr )
    end if IFRSTATS

    !> Compute the non-linear term.
    IFNLTYPE: if ( F_TYPE == F_ISOTROPIC ) then
      call f_RK3ComputeNonlinearIsotropic ( ierr )
    else
      select case ( NS_DEFORMED )
      !-------------------------------
        case ( NS_BF_ROT ) 
          call f_RK3ComputeNonlinearBFRot ( ierr )
        case ( NS_R_ROT ) 
          call f_RK3ComputeNonlinearIsotropic ( ierr )
        case ( NS_R_CONV ) 
          call f_RK3ComputeNonlinearRConv ( ierr )
      end select
    end if IFNLTYPE

    !> Do a 2D FFT of the non-linear term, transpose it and store it 
    !! in the vorticity buffer array.
    call f_RK3TransformNonLinearjiRs2Fs ( zslice, ierr )
               
    !> If there are particles, call the particle routines.
    IFPART: if ( P_TRACK_PART == YES ) then
      IFRELEASED: if ( tstep >= P_TS_RELEASE ) then

!        call PetscPrintf( PETSC_COMM_WORLD, 'Copying from 2D arrays &
!             to 3D arrays. \n', ierr )
        !> @todo copy 2D velocity to 3D velocity        
!        write(*,*) 'rank: ', MYRANK, ' slice: ', zslice, &
!                   ' min/max: ', z_min, z_max
        DOYVAL: do yval = y_min_2r, y_max_2r, 1
!          arr_u1_3r(x_min_2r:x_max_2r,zslice,yval) = arr_u1_2r(:,yval)
!          arr_u2_3r(x_min_2r:x_max_2r,zslice,yval) = arr_u2_2r(:,yval)
!          arr_u3_3r(x_min_2r:x_max_2r,zslice,yval) = arr_u3_2r(:,yval)
          arr_u1_3r(:,zslice,yval) = arr_u1_2r(:,yval)
          arr_u2_3r(:,zslice,yval) = arr_u2_2r(:,yval)
          arr_u3_3r(:,zslice,yval) = arr_u3_2r(:,yval)
        end do DOYVAL

!        write(*,*) '2D u1 array on slice ', zslice, ':', arr_u1_2r(:,:), &
!                   '2D u2 array on slice ', zslice, ':', arr_u2_2r(:,:), &
!                   '2D u3 array on slice ', zslice, ':', arr_u3_2r(:,:)

      end if IFRELEASED
    end if IFPART

  end do DOSLICES

  !> If first stage get real-space fluid statistics and then compute dt.
  !! Need dt for particle RK3. If applicable, also collect 1D spectra in
  !! wave space
  IFSTAGE1STATS: if ( stage == 1 ) then
    !> Real-space statistics
    call f_FluidstatsRsCollect ( (tstep-1), ierr )
    call f_RK3Tstep ( deltat, ierr )

    !> Wave-space statistics
    IFNOTISOTROPICSTATS: if ( F_TYPE /= F_ISOTROPIC ) then
!      call PetscPrintf( PETSC_COMM_WORLD, 'Collecting FS stats \n', ierr )
      call f_FluidstatsFsCollect ( (tstep-1), ierr )
    else if ( TS_SPEC1D <= TSTEPS ) then
      IFWRITENOW: if ( modulo((tstep-1),TS_SPEC1D) == 0 ) then
        call f_FluidstatsFsCollect ( (tstep-1), ierr )
      end if IFWRITENOW
    end if IFNOTISOTROPICSTATS
  end if IFSTAGE1STATS

  IFPARTRK3: if ( P_TRACK_PART == YES ) then
    IFRELEASEDRK3: if ( tstep >= P_TS_RELEASE ) then

!      !> Need to normalise energy by volume of the domain.
!      No, it is already correct!!
!      IFENERGYNORM: if ( F_TYPE /= F_ISOTROPIC ) then
!        arr_u1_3r = sqrt( bmat(1,1) * bmat(2,2) * bmat(3,3) ) * arr_u1_3r
!        arr_u2_3r = sqrt( bmat(1,1) * bmat(2,2) * bmat(3,3) ) * arr_u2_3r
!        arr_u3_3r = sqrt( bmat(1,1) * bmat(2,2) * bmat(3,3) ) * arr_u3_3r
!        write(*,*) sqrt( bmat(1,1) * bmat(2,2) * bmat(3,3) )
!      end if IFENERGYNORM

      !> Need to return the arrays here, because will use local
      !! vectors for the particle routines.
!      call PetscPrintf( PETSC_COMM_WORLD, 'Restoring global fluid arrays &
!           to obtain local arrays with halos. \n', ierr )
      call DMDAVecRestoreArrayF90 ( da3r, u1_3r, arr_u1_3r, ierr )
      call DMDAVecRestoreArrayF90 ( da3r, u2_3r, arr_u2_3r, ierr )
      call DMDAVecRestoreArrayF90 ( da3r, u3_3r, arr_u3_3r, ierr )

      !> Solve the equation of motion for all particles
 !     call PetscPrintf( PETSC_COMM_WORLD, 'Solving the &
 !          particle equation. \n', ierr )
      call p_ControlRK3 ( stage, tstep, deltat, ierr )

      !> Update particle positions
 !     call PetscPrintf( PETSC_COMM_WORLD, 'Finishing the  &
 !          particle routines for this time step. \n', ierr )
      call p_ControlFinishTimeStep ( tstep, ierr )

    end if IFRELEASEDRK3
  end if IFPARTRK3

  !> Return velocity and vorticity components to PETSc.
  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

  IFRESTORECONV: if ( NS_DEFORMED == NS_R_CONV ) then
    call DMDAVecRestoreArrayF90 ( da1r2w, u7_1r2w, arr_u7_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u8_1r2w, arr_u8_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u9_1r2w, arr_u9_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u10_1r2w, arr_u10_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u11_1r2w, arr_u11_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da1r2w, u12_1r2w, arr_u12_1r2w, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u7_2r, arr_u7_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u8_2r, arr_u8_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u9_2r, arr_u9_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u10_2r, arr_u10_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u11_2r, arr_u11_2r, perr )
    call DMDAVecRestoreArrayF90 ( da2r, u12_2r, arr_u12_2r, perr )
  end if IFRESTORECONV

  !> Arrays for 1D spectra only
  IFISOTROPICRESTOREARRAY: if ( F_TYPE /= F_ISOTROPIC ) then

    call DMDAVecRestoreArrayF90 ( da1r1w, u1_1r1w, arr_u1_1r1w, perr )
    call DMDAVecRestoreArrayF90 ( da1r1w, u2_1r1w, arr_u2_1r1w, perr )
    call DMDAVecRestoreArrayF90 ( da1r1w, u3_1r1w, arr_u3_1r1w, perr )

  else
    IFSPEC1DRESTOREARRAY: if ( TS_SPEC1D <= TSTEPS ) then

      call DMDAVecRestoreArrayF90 ( da1r1w, u1_1r1w, arr_u1_1r1w, perr )
      call DMDAVecRestoreArrayF90 ( da1r1w, u2_1r1w, arr_u2_1r1w, perr )
      call DMDAVecRestoreArrayF90 ( da1r1w, u3_1r1w, arr_u3_1r1w, perr )

    end if IFSPEC1DRESTOREARRAY
  end if IFISOTROPICRESTOREARRAY


!> @todo implement transform of coupling term - go through slices again.
!  DOTWOWAYSLICES: do zslice = z_min, z_max, 1
!    !> If the flow is two-way coupled, store the coupling term in the
!    !! velocity arrays and proceed as with the vorticity.
!    IFPARTICLESTRANSFORM: if ( P_TRACK_PART == YES) then
!      IFTWOWAYTRANSFORM: if ( P_TWO_WAY == YES ) then
!        call f_RK3TransformCouplingTermjiRs2Fs ( zslice, ierr )
!      end if IFTWOWAYTRANSFORM
!    end if IFPARTICLESTRANSFORM

!  end do DOTWOWAYSLICES
end subroutine f_RK3RealSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformToRealSpace
!> @stage stage of RK3 for statistics
!> @tstep previous (!) timestep for statistics
!> @param ierr should return 0
subroutine f_RK3TransformToRealSpace ( zslice, stage, tstep, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X, &
                          NS_DEFORMED, NS_R_CONV, &
                          TS_SPEC1D, TSTEPS, F_TYPE, F_ISOTROPIC
  use f_arrays,     only: z_min, i_min_1r2w, j_min_1r2w, y_min_2r, &
                          i_max_1r2w, j_max_1r2w, y_max_2r, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          arr_u7_2r, arr_u8_2r, arr_u9_2r, &
                          arr_u10_2r, arr_u11_2r, arr_u12_2r, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          arr_u7_1r2w, arr_u8_1r2w, arr_u9_1r2w, &
                          arr_u10_1r2w, arr_u11_1r2w, arr_u12_1r2w, &
                          nprocsi_1r2w, i_allwidths_1r2w, i_maxwidth_1r2w, &
                          arr_u1_1r1w, arr_u2_1r1w, arr_u3_1r1w
  use f_fftw,       only: f_Fftw_ji_Fs2Rs, arr_j_in_da, arr_i_real

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: zslice
  integer,intent(in)                     :: stage 
  integer,intent(in)                     :: tstep 

  integer                                :: slicelocal
  integer                                :: i, ii, j, jj, proc_i, itranspose, jtranspose, i2r, component

  type(C_PTR)                            :: inputarray, copytoarray, spectrumarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = zslice - z_min + 1

  !> If 1D spectrum is evaluated, get the data here.
  IFSTAGE: if ( stage /= 1 ) then
    !> u1 component
    inputarray = c_loc(arr_u1_1r2w)
    copytoarray = c_loc(arr_u1_2r)
!    write(*,*) 'u1'
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u2 component
    inputarray = c_loc(arr_u2_1r2w)
    copytoarray = c_loc(arr_u2_2r)
!    write(*,*) 'u2'
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u3 component
    inputarray = c_loc(arr_u3_1r2w)
    copytoarray = c_loc(arr_u3_2r)
!    write(*,*) 'u3'
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )
  else
    IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
      IF1DSPEC: if ( TS_SPEC1D > TSTEPS ) then
        !> u1 component
        inputarray = c_loc(arr_u1_1r2w)
        copytoarray = c_loc(arr_u1_2r)
!    write(*,*) 'u1'
        call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )
  
        !> u2 component
        inputarray = c_loc(arr_u2_1r2w)
        copytoarray = c_loc(arr_u2_2r)
!    write(*,*) 'u2'
        call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

        !> u3 component
        inputarray = c_loc(arr_u3_1r2w)
        copytoarray = c_loc(arr_u3_2r)
!    write(*,*) 'u3'
        call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )
      else
        IFSPECNOW: if ( modulo(tstep, TS_SPEC1D) == 0 ) then
          !> u1 component
          inputarray = c_loc(arr_u1_1r2w)
          copytoarray = c_loc(arr_u1_2r)
          spectrumarray = c_loc(arr_u1_1r1w)
!    write(*,*) 'u1'
          call f_RK3TransformjiFs2RsGetSpectrum ( inputarray, copytoarray, &
                                                  spectrumarray, slicelocal, ierr )
  
          !> u2 component
          inputarray = c_loc(arr_u2_1r2w)
          copytoarray = c_loc(arr_u2_2r)
          spectrumarray = c_loc(arr_u2_1r1w)
!    write(*,*) 'u2'
          call f_RK3TransformjiFs2RsGetSpectrum ( inputarray, copytoarray, &
                                                  spectrumarray, slicelocal, ierr )

          !> u3 component
          inputarray = c_loc(arr_u3_1r2w)
          copytoarray = c_loc(arr_u3_2r)
          spectrumarray = c_loc(arr_u3_1r1w)
!    write(*,*) 'u3'
          call f_RK3TransformjiFs2RsGetSpectrum ( inputarray, copytoarray, &
                                                  spectrumarray, slicelocal, ierr )
        end if IFSPECNOW
      end if IF1DSPEC
    else
      !> u1 component
      inputarray = c_loc(arr_u1_1r2w)
      copytoarray = c_loc(arr_u1_2r)
      spectrumarray = c_loc(arr_u1_1r1w)
      call f_RK3TransformjiFs2RsGetSpectrum ( inputarray, copytoarray, &
                                              spectrumarray, slicelocal, ierr )

      !> u2 component
      inputarray = c_loc(arr_u2_1r2w)
      copytoarray = c_loc(arr_u2_2r)
      spectrumarray = c_loc(arr_u2_1r1w)
      call f_RK3TransformjiFs2RsGetSpectrum ( inputarray, copytoarray, &
                                              spectrumarray, slicelocal, ierr )

      !> u3 component
      inputarray = c_loc(arr_u3_1r2w)
      copytoarray = c_loc(arr_u3_2r)
      spectrumarray = c_loc(arr_u3_1r1w)
      call f_RK3TransformjiFs2RsGetSpectrum ( inputarray, copytoarray, &
                                              spectrumarray, slicelocal, ierr )
    end if IFISOTROPIC
  end if IFSTAGE

  !> u4 component
  inputarray = c_loc(arr_u4_1r2w)
  copytoarray = c_loc(arr_u4_2r)
!    write(*,*) 'u4'
  call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  !> u5 component
  inputarray = c_loc(arr_u5_1r2w)
  copytoarray = c_loc(arr_u5_2r)
!    write(*,*) 'u5'
  call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  !> u6 component
  inputarray = c_loc(arr_u6_1r2w)
  copytoarray = c_loc(arr_u6_2r)
!    write(*,*) 'u6'
  call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  IFCONV: if ( NS_DEFORMED == NS_R_CONV ) then

    !> u7 component
    inputarray = c_loc(arr_u7_1r2w)
    copytoarray = c_loc(arr_u7_2r)
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u8 component
    inputarray = c_loc(arr_u8_1r2w)
    copytoarray = c_loc(arr_u8_2r)
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u9 component
    inputarray = c_loc(arr_u9_1r2w)
    copytoarray = c_loc(arr_u9_2r)
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u10 component
    inputarray = c_loc(arr_u10_1r2w)
    copytoarray = c_loc(arr_u10_2r)
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u11 component
    inputarray = c_loc(arr_u11_1r2w)
    copytoarray = c_loc(arr_u11_2r)
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

    !> u12 component
    inputarray = c_loc(arr_u12_1r2w)
    copytoarray = c_loc(arr_u12_2r)
    call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  end if IFCONV

!  write(*,*) 'arr_u1_2r: ', arr_u1_2r
!  write(*,*) 'arr_u2_2r: ', arr_u2_2r
!  write(*,*) 'arr_u3_2r: ', arr_u3_2r
!  write(*,*) 'arr_u4_2r: ', arr_u4_2r
!  write(*,*) 'arr_u5_2r: ', arr_u5_2r
!  write(*,*) 'arr_u6_2r: ', arr_u6_2r

end subroutine f_RK3TransformToRealSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformjiFs2Rs
!> @param ierr should return 0
  subroutine f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use g_constants,  only: ZERO
  use f_arrays,     only: j_width_1r2w, z_width_1r2w, i_width_1r2w, x_width_2r, y_width_2r
  use f_fftw,       only: f_Fftw_ji_Fs2Rs, arr_j_in_da, arr_i_real

  implicit none
 
  type(C_PTR), intent(in)                :: inputarray, copytoarray
  integer,intent(in)                     :: slicelocal
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                          :: i, j, jj 
  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:)   :: fromarray
  real(kind=C_DOUBLE),pointer,dimension(:,:)       :: targetarray

  !> Link inputarray to the C pointer of the component
  call c_f_pointer(inputarray, fromarray, [2, j_width_1r2w,z_width_1r2w,i_width_1r2w])

  DOI1R2W: do i = 1, i_width_1r2w, 1
    jj = 1

    !> Copy positive wavenumbers
    DOJPOSITIVE: do j = 1, KJ_MAX + 1, 1
     !> Copy to buffer
      arr_j_in_da(:,jj,i) = fromarray(:,j,slicelocal,i)
      jj = jj + 1
    end do DOJPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOJZEROS: do j = KJ_MAX + 2, NODES_Y - KJ_MAX 
      arr_j_in_da(:,jj,i) = ZERO 
      jj = jj + 1
    end do DOJZEROS

   !> Copy negative wavenumbers
    DOJNEGATIVE: do j = KJ_MAX + 2, j_width_1r2w
    !> Copy to buffer
      arr_j_in_da(:,jj,i) = fromarray(:,j,slicelocal,i)
      jj = jj + 1
    end do DOJNEGATIVE
  end do DOI1R2W

! Once copying is finished for a whole plane, perform an FFT for the component
  !> Do FFT in j direction and transpose from ji to ij, then do FFT in i direction
  call f_Fftw_ji_Fs2Rs ( ierr )

  !> Link targetarray to the C pointer of the component
  call c_f_pointer(copytoarray, targetarray, [x_width_2r, y_width_2r])

! Copy the result of the FFT to the buffer array
  DOY: do j = 1, y_width_2r, 1                     ! parallelised in this direction
    DOI: do i = 1, NODES_X, 1

        targetarray(i, j) = arr_i_real(i, j)

    end do DOI
  end do DOY

end subroutine f_RK3TransformjiFs2Rs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformjiFs2RsGetSpectrum
!> @param ierr should return 0
  subroutine f_RK3TransformjiFs2RsGetSpectrum ( inputpointer, copytopointer, spectrumpointer, slicelocal, ierr )

  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use g_constants,  only: ZERO
  use f_arrays,     only: j_width_1r2w, z_width_1r2w, i_width_1r2w, x_width_2r, y_width_2r, &
                          i_width_1r1w, y_width_1r1w
  use f_fftw,       only: f_Fftw_ji_Fs2RsPrepare, f_Fftw_ji_Fs2RsFinalise, arr_j_in_da, &
                          arr_i_complex_da, arr_i_real

  implicit none
 
  type(C_PTR), intent(in)                :: inputpointer, copytopointer, spectrumpointer
  integer,intent(in)                     :: slicelocal
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                          :: i, j, jj 
  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:)   :: inputarray
  real(kind=C_DOUBLE),pointer,dimension(:,:)       :: copytoarray
  real(kind=C_DOUBLE),pointer,dimension(:,:,:)       :: spectrumarray

  !> Link inputarray to the C pointer of the component
  call c_f_pointer(inputpointer, inputarray, [2, j_width_1r2w,z_width_1r2w,i_width_1r2w])

  DOI1R2W: do i = 1, i_width_1r2w, 1
    jj = 1

    !> Copy positive wavenumbers
    DOJPOSITIVE: do j = 1, KJ_MAX + 1, 1
     !> Copy to buffer
      arr_j_in_da(:,jj,i) = inputarray(:,j,slicelocal,i)
      jj = jj + 1
    end do DOJPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOJZEROS: do j = KJ_MAX + 2, NODES_Y - KJ_MAX 
      arr_j_in_da(:,jj,i) = ZERO 
      jj = jj + 1
    end do DOJZEROS

   !> Copy negative wavenumbers
    DOJNEGATIVE: do j = KJ_MAX + 2, j_width_1r2w
    !> Copy to buffer
      arr_j_in_da(:,jj,i) = inputarray(:,j,slicelocal,i)
      jj = jj + 1
    end do DOJNEGATIVE
  end do DOI1R2W


! Once copying is finished for a whole plane, perform an FFT for the component
  !> Do FFT in j direction and transpose from ji to ij
  call f_Fftw_ji_Fs2RsPrepare ( ierr )

    !> Link spectrum array to the C pointer of the component
    call c_f_pointer(spectrumpointer, spectrumarray, [2, i_width_1r1w, y_width_1r1w])

!    write(*,*) 'Shape of spectrum array: ', shape(spectrumarray), &
!               'shape of FFTW array: ', shape(arr_i_complex_da)
    !> Get array for 1D spectrum
    DOYSPEC: do j = 1, y_width_1r1w, 1                     ! parallelised in this direction
      DOISPEC: do i = 1, i_width_1r1w, 1

!        write(*,*) i, i_width_1r1w, j, y_width_1r1w, arr_i_complex_da(:,i, j), spectrumarray(:,i, j)
        spectrumarray(:,i, j) = arr_i_complex_da(:,i, j) / &
                                real(NODES_X,kind=C_DOUBLE)

      end do DOISPEC
    end do DOYSPEC


  !> Do FFT in i direction
  call f_Fftw_ji_Fs2RsFinalise ( ierr )


  !> Link targetarray to the C pointer of the component
  call c_f_pointer(copytopointer, copytoarray, [x_width_2r, y_width_2r])

! Copy the result of the FFT to the buffer array
  DOY: do j = 1, y_width_2r, 1                     ! parallelised in this direction
    DOX: do i = 1, NODES_X, 1

        copytoarray(i, j) = arr_i_real(i, j)

    end do DOX
  end do DOY

end subroutine f_RK3TransformjiFs2RsGetSpectrum
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3ComputeNonlinearIsotropic
!> @param ierr should return 0
subroutine f_RK3ComputeNonlinearIsotropic ( ierr )

  use g_constants,  only : ZERO
  use f_arrays,     only: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          x_min_2r, x_max_2r, y_min_2r, y_max_2r 
 
  implicit none
 
  PetscErrorCode,intent(inout)   :: ierr
  real(kind=C_DOUBLE)     :: nonlinear1, nonlinear2, nonlinear3
  integer                 :: i,j

!  write(*,*) 'arr_u1_2r: ', arr_u1_2r
!  write(*,*) 'arr_u2_2r: ', arr_u2_2r
!  write(*,*) 'arr_u3_2r: ', arr_u3_2r
  !> calculate cross product of velocity and vorticity
  DOY: do j = y_min_2r, y_max_2r, 1
      DOX: do i = x_min_2r, x_max_2r, 1
!        write(*,*) i, j, arr_u1_2r(i,j), arr_u2_2r(i,j), arr_u3_2r(i,j)
        nonlinear1 = arr_u3_2r(i,j) * arr_u5_2r(i,j) - arr_u2_2r(i,j) * arr_u6_2r(i,j)        
        nonlinear2 = arr_u1_2r(i,j) * arr_u6_2r(i,j) - arr_u3_2r(i,j) * arr_u4_2r(i,j)        
        nonlinear3 = arr_u2_2r(i,j) * arr_u4_2r(i,j) - arr_u1_2r(i,j) * arr_u5_2r(i,j)
!        write(*,*) i, j, arr_u1_2r(i,j), arr_u4_2r(i,j), nonlinear1, &
!                   arr_u2_2r(i,j), arr_u5_2r(i,j), nonlinear2, arr_u3_2r(i,j), arr_u6_2r(i,j), nonlinear3 
        arr_u4_2r(i,j) = nonlinear1        
        arr_u5_2r(i,j) = nonlinear2        
        arr_u6_2r(i,j) = nonlinear3       
    end do DOX
  end do DOY

end subroutine f_RK3ComputeNonlinearIsotropic
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3ComputeNonlinearBFRot
!> @param ierr should return 0
subroutine f_RK3ComputeNonlinearBFRot ( ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : binv, bmat, bvalexist
  use f_arrays,     only : arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                           arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                           x_min_2r, x_max_2r, y_min_2r, y_max_2r 
 
  implicit none
 
  PetscErrorCode,intent(inout)              :: ierr

  real(kind=C_DOUBLE)                :: nonlinear1, nonlinear2, nonlinear3
  real(kind=C_DOUBLE),dimension(1:3) :: u
  real(kind=C_DOUBLE),dimension(1:3) :: omega
  real(kind=C_DOUBLE),dimension(1:3) :: nonlinear

  !> Loop variables for array
  integer                 :: i,j
  !> Loop variables for tensor
  integer                 :: m

  !> calculate cross product of velocity and vorticity
  DOY: do j = y_min_2r, y_max_2r, 1
    DOX: do i = x_min_2r, x_max_2r, 1

      nonlinear1 = ZERO
      nonlinear2 = ZERO
      nonlinear3 = ZERO
      u(1) = arr_u1_2r(i,j)
      u(2) = arr_u2_2r(i,j)
      u(3) = arr_u3_2r(i,j)
      omega(1) = arr_u4_2r(i,j)
      omega(2) = arr_u5_2r(i,j)
      omega(3) = arr_u6_2r(i,j)

      !> compute tensor product for array value (i,j)
      DOM: do m = 1, 3, 1
        IF1M: if ( bvalexist(1,m) .eqv. .TRUE. ) then
          nonlinear3 = nonlinear3 - ( omega(2) * binv(1,m) * u(m) )
          nonlinear2 = nonlinear2 + ( omega(3) * binv(1,m) * u(m) ) 
        end if IF1M
        IF2M: if ( bvalexist(2,m) .eqv. .TRUE. ) then
          nonlinear1 = nonlinear1 - ( omega(3) * binv(2,m) * u(m) )
          nonlinear3 = nonlinear3 + ( omega(1) * binv(2,m) * u(m) ) 
        end if IF2M
        IF3M: if ( bvalexist(3,m) .eqv. .TRUE. ) then
          nonlinear1 = nonlinear1 + ( omega(2) * binv(3,m) * u(m) ) 
          nonlinear2 = nonlinear2 - ( omega(1) * binv(3,m) * u(m) )
        end if IF3M
      end do DOM

      arr_u4_2r(i,j) = nonlinear1        
      arr_u5_2r(i,j) = nonlinear2        
      arr_u6_2r(i,j) = nonlinear3       

!      nonlinear(:) = ZERO

!      DOM2: do m = 1, 3, 1
!        IFBM1: if ( bvalexist(m,1) .eqv. .TRUE. ) then
!          nonlinear(m) = nonlinear(m) + bmat(m,1) * nonlinear1
!        end if IFBM1
!        IFBM2: if ( bvalexist(m,2) .eqv. .TRUE. ) then
!          nonlinear(m) = nonlinear(m) + bmat(m,2) * nonlinear2
!        end if IFBM2
!        IFBM3: if ( bvalexist(m,3) .eqv. .TRUE. ) then
!          nonlinear(m) = nonlinear(m) + bmat(m,3) * nonlinear3
!        end if IFBM3
!      end do DOM2

!      arr_u4_2r(i,j) = nonlinear(1)        
!      arr_u5_2r(i,j) = nonlinear(2)      
!      arr_u6_2r(i,j) = nonlinear(3)      

    end do DOX
  end do DOY

end subroutine f_RK3ComputeNonlinearBFRot
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3ComputeNonlinearRConv
!> @param ierr should return 0
subroutine f_RK3ComputeNonlinearRConv ( ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, bvalexist
  use f_arrays,     only: arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          arr_u7_2r, arr_u8_2r, arr_u9_2r, &
                          arr_u10_2r, arr_u11_2r, arr_u12_2r, &
                          x_min_2r, x_max_2r, y_min_2r, y_max_2r 
 
  implicit none
 
  PetscErrorCode,intent(inout)                   :: ierr
  real(kind=C_DOUBLE),dimension(1:3)      :: nonlinear, uj
  real(kind=C_DOUBLE),dimension(1:3,1:3)  :: duidk
  integer                                 :: i, j, k, x, y


!  Bkj duidk uj

  DOY: do y = y_min_2r, y_max_2r, 1
    DOX: do x = x_min_2r, x_max_2r, 1

      uj(1) = arr_u1_2r(x,y)
      uj(2) = arr_u2_2r(x,y)
      uj(3) = arr_u3_2r(x,y)

      duidk(1,1) = arr_u4_2r(x,y)
      duidk(1,2) = arr_u5_2r(x,y)
      duidk(1,3) = arr_u6_2r(x,y)
      duidk(2,1) = arr_u7_2r(x,y)
      duidk(2,2) = arr_u8_2r(x,y)
      duidk(2,3) = arr_u9_2r(x,y)
      duidk(3,1) = arr_u10_2r(x,y)
      duidk(3,2) = arr_u11_2r(x,y)
      duidk(3,3) = arr_u12_2r(x,y)

      nonlinear = ZERO

      DOI: do i = 1, 3, 1
        DOJ: do j = 1, 3, 1
          DOK: do k = 1, 3, 1
            IFBKJ: if ( bvalexist(k,j) .eqv. .TRUE. ) then
              nonlinear(i) = nonlinear(i) + &
                               ( bmat(k,j) * duidk(k,i) * uj(j) )
            end if IFBKJ
          end do DOK
        end do DOJ
      end do DOI

      arr_u4_2r(x,y) = nonlinear(1)      
      arr_u5_2r(x,y) = nonlinear(2)        
      arr_u6_2r(x,y) = nonlinear(3)     

    end do DOX
  end do DOY

end subroutine f_RK3ComputeNonlinearRConv
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformNonLinearjiRs2Fs
!> @param ierr should return 0
subroutine f_RK3TransformNonLinearjiRs2Fs ( zslice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use f_arrays,     only: z_min, i_min_1r2w, j_min_1r2w, y_min_2r, &
                          i_max_1r2w, j_max_1r2w, y_max_2r, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          nprocsi_1r2w, i_allwidths_1r2w, i_maxwidth_1r2w
  use f_fftw,       only: f_Fftw_ji_Fs2Rs, arr_j_out_da, arr_i_real

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: zslice

  integer                                :: slicelocal
  integer                                :: i, ii, j, jj, proc_i, itranspose, jtranspose, i2r, component

  type(C_PTR)                            :: inputarray, copytoarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = zslice - z_min + 1

  !> u4 component
  inputarray = c_loc(arr_u4_2r)
  copytoarray = c_loc(arr_u4_1r2w)
  call f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u5 component
  inputarray = c_loc(arr_u5_2r)
  copytoarray = c_loc(arr_u5_1r2w)
  call f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u6 component
  inputarray = c_loc(arr_u6_2r)
  copytoarray = c_loc(arr_u6_1r2w)
  call f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

end subroutine f_RK3TransformNonLinearjiRs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformCouplingTermjiRs2Fs
!> @param ierr should return 0
subroutine f_RK3TransformCouplingTermjiRs2Fs ( zslice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use f_arrays,     only: z_min, i_min_1r2w, j_min_1r2w, y_min_2r, &
                          i_max_1r2w, j_max_1r2w, y_max_2r, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          nprocsi_1r2w, i_allwidths_1r2w, i_maxwidth_1r2w
  use f_fftw,       only: f_Fftw_ji_Fs2Rs, arr_j_out_da, arr_i_real

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: zslice

  integer                                :: slicelocal
  integer                                :: i, ii, j, jj, proc_i, itranspose, jtranspose, i2r, component

  type(C_PTR)                            :: inputarray, copytoarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = zslice - z_min + 1

  !> u1 component
  inputarray = c_loc(arr_u1_2r)
  copytoarray = c_loc(arr_u1_1r2w)
  call f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u2 component
  inputarray = c_loc(arr_u2_2r)
  copytoarray = c_loc(arr_u2_1r2w)
  call f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u3 component
  inputarray = c_loc(arr_u3_2r)
  copytoarray = c_loc(arr_u3_1r2w)
  call f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

end subroutine f_RK3TransformCouplingTermjiRs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformjiRs2Fs
!> @param ierr should return 0
  subroutine f_RK3TransformjiRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use g_constants,  only: ZERO
  use f_arrays,     only: i_min_1r2w, j_min_1r2w, y_min_2r, &
                          i_max_1r2w, j_max_1r2w, y_max_2r, &
                          j_width_1r2w, z_width_1r2w, i_width_1r2w, x_width_2r, y_width_2r
  use f_fftw,       only: arr_j_out_da, arr_i_real, f_Fftw_ji_Rs2Fs

  implicit none
 
  type(C_PTR), intent(in)                :: inputarray, copytoarray
  integer,intent(in)                     :: slicelocal
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                          :: i, j, jj 
  real(kind=C_DOUBLE),pointer,dimension(:,:)       :: fromarray
  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:)   :: targetarray

  !> Link inputarray to the C pointer of the component
  call c_f_pointer(inputarray, fromarray, [x_width_2r, y_width_2r])

! Copy from PETSc array to FFTW buffer for FFT in i direction
  DOY: do j = 1, y_width_2r, 1                     ! parallelised in this direction
    DOI: do i = 1, NODES_X, 1

        arr_i_real(i, j) = fromarray(i, j)

    end do DOI
  end do DOY

! Once copying is finished for a whole plane, perform an FFT for the component
  !> Do FFT in j direction and transpose from ji to ij, then do FFT in i direction
  call f_Fftw_ji_Rs2Fs ( ierr )

  !> Link targetarray to the C pointer of the component
  call c_f_pointer(copytoarray, targetarray, [2, j_width_1r2w,z_width_1r2w,i_width_1r2w])

  DOI1R2W: do i = 1, i_width_1r2w, 1
    jj = 1

    !> Copy positive wavenumbers
    DOJPOSITIVE: do j = 1, KJ_MAX + 1, 1
     !> Copy to buffer
      targetarray(:,j,slicelocal,i) = arr_j_out_da(:,jj,i)

      jj = jj + 1
    end do DOJPOSITIVE

    !> Get rid of the highest 1/3 wavenumbers, both positive and negative
!    jj = NODES_Y - KJ_MAX + 2
    jj = NODES_Y - KJ_MAX + 1

   !> Copy negative wavenumbers
    DOJNEGATIVE: do j = KJ_MAX + 2, j_width_1r2w 
     !> Copy to buffer
      targetarray(:,j,slicelocal,i) = arr_j_out_da(:,jj,i)
      jj = jj + 1
    end do DOJNEGATIVE
  end do DOI1R2W

end subroutine f_RK3TransformjiRs2Fs
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
! subroutine f_RK3FormG1U2
!> @param ierr should return 0
subroutine f_RK3FormG1U2 ( dt, tstep, ierr )

  implicit none
 
  real(kind=C_DOUBLE),intent(in)   :: dt
  integer,intent(in)               :: tstep
  PetscErrorCode,intent(inout)            :: ierr

  call f_RK3SolveNavierStokes ( dt, tstep, 1, ierr )
 
end subroutine f_RK3FormG1U2
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3FormG2U3
!> @param ierr should return 0
subroutine f_RK3FormG2U3 ( dt, tstep, ierr )

  implicit none
 
  real(kind=C_DOUBLE),intent(in)   :: dt
  integer,intent(in)               :: tstep
  PetscErrorCode,intent(inout)            :: ierr

  call f_RK3SolveNavierStokes ( dt, tstep, 2, ierr )
 
end subroutine f_RK3FormG2U3
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3FormG3Un
!> @param ierr should return 0
subroutine f_RK3FormG3Un ( dt, tstep, ierr )

  implicit none
 
  real(kind=C_DOUBLE),intent(in)   :: dt
  integer,intent(in)               :: tstep
  PetscErrorCode,intent(inout)            :: ierr

  call f_RK3SolveNavierStokes ( dt, tstep, 3, ierr )

end subroutine f_RK3FormG3Un
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3SolveNavierStokes
!> @param ierr should return 0
subroutine f_RK3SolveNavierStokes ( dt, tstep, stage, ierr )

  use g_parameters, only: P_TRACK_PART, P_TWO_WAY, YES, MYRANK, MASTER, &
                          F_TYPE, F_ISOTROPIC, DEALIASING, &
                          DEALIASING_SPHERICAL, DEALIASING_LES, &
                          NS_DEFORMED, NS_BF_ROT, NS_R_CONV, NS_R_ROT
  use f_arrays,     only: da3w, da1r2w, da2w, i_min_3w, i_max_3w, &
                          u1_3w, u2_3w, u3_3w, &
                          rk3_3w_1, rk3_3w_2, rk3_3w_3, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          u1_1r2w, u2_1r2w, u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          u1_2w, u2_2w, u3_2w, &
                          u4_2w, u5_2w, u6_2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          f_ArraysDealiasSpherical, f_ArraysLESCutoff
  use f_force,      only: f_ForceCalc

  implicit none
 
  real(kind=C_DOUBLE),intent(in) :: dt
  integer,intent(in)             :: tstep, stage
  PetscErrorCode,intent(inout)          :: ierr

  integer                        :: islice                 
 
  PetscErrorCode        :: perr

!> Get read access to all non-linear components
  call DMDAVecGetArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
!> Get write access to the corresponding 2D arrays in wave space
  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> If two-way coupled, get read access to all coupling components
  IFPARTICLESGET: if ( P_TRACK_PART == YES) then
    IFTWOWAYGET: if ( P_TWO_WAY == YES ) then

      call DMDAVecGetArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
      call DMDAVecGetArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
      call DMDAVecGetArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
      !> Get write access to the corresponding 2D arrays in wave space
      call DMDAVecGetArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
      call DMDAVecGetArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
      call DMDAVecGetArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

    end if IFTWOWAYGET
  end if IFPARTICLESGET

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Get write access to all RK3 components
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecGetArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )

!> Compute forcing if applicable
  call f_ForceCalc ( tstep, dt, ierr )

  !> Loop through slices in i direction
  DOSLICES: do islice = i_min_3w, i_max_3w, 1

    !> Transform non-linear components.
    call f_RK3TransformNonLinearkRs2Fs ( islice, ierr )

    !> If the flow is two-way coupled, transform coupling components.
    IFPARTICLESTRANSFORM: if ( P_TRACK_PART == YES) then
      IFTWOWAYTRANSFORM: if ( P_TWO_WAY == YES ) then
        call f_RK3TransformCouplingTermkRs2Fs ( islice, ierr )
      end if IFTWOWAYTRANSFORM
    end if IFPARTICLESTRANSFORM
               
    !> Solve the Navier-Stokes equation.
    ! if isotropic: 
    IFFLOWTYPE: if ( F_TYPE == F_ISOTROPIC ) then
      call f_RK3FormFIsotropic ( islice, stage, dt, ierr )
    else
      select case ( NS_DEFORMED )
      !-------------------------------
        case ( NS_BF_ROT ) 
          call f_RK3FormFBFRot ( islice, stage, dt, ierr )
        case ( NS_R_ROT ) 
          call f_RK3FormFRRot ( islice, stage, dt, ierr )
        case ( NS_R_CONV ) 
          !> Once the non-linear term is computed, the two Rogallo
          !! formulations are identical, so no need for separate
          !! routines
          call f_RK3FormFRRot ( islice, stage, dt, ierr )
      end select
    end if IFFLOWTYPE 

  end do DOSLICES
  !end loop

  !> Anything here will be from the non-linear term
  IFSPHER: if ( DEALIASING == DEALIASING_SPHERICAL ) then
  !> Truncate all wavemodes outside a 2/3 sphere
    call f_ArraysDealiasSpherical ( ierr )
  end if IFSPHER

  IFLES: if ( DEALIASING == DEALIASING_LES ) then
  !> Truncate all wavemodes outside cutoff
    call f_ArraysLESCutoff ( ierr )
  end if IFLES  

!> Return all non-linear components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecRestoreArrayReadF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

!> If two-way coupled, return all coupling components to PETSc
  IFPARTICLESRESTORE: if ( P_TRACK_PART == YES) then
    IFTWOWAYRESTORE: if ( P_TWO_WAY == YES ) then
      call DMDAVecRestoreArrayReadF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
      call DMDAVecRestoreArrayReadF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
      call DMDAVecRestoreArrayReadF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )
      call DMDAVecRestoreArrayF90 ( da2w, u1_2w, arr_u1_2w, perr )
      call DMDAVecRestoreArrayF90 ( da2w, u2_2w, arr_u2_2w, perr )
      call DMDAVecRestoreArrayF90 ( da2w, u3_2w, arr_u3_2w, perr )

    end if IFTWOWAYRESTORE
  end if IFPARTICLESRESTORE

!> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

!> Return all RK3 components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_1, arr_rk3_3w_1, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_2, arr_rk3_3w_2, perr )
  call DMDAVecRestoreArrayF90 ( da3w, rk3_3w_3, arr_rk3_3w_3, perr )

end subroutine f_RK3SolveNavierStokes
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformNonLinearkRs2Fs
!> @param ierr should return 0
subroutine f_RK3TransformNonLinearkRs2Fs ( islice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use f_arrays,     only: i_min_3w, arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: islice

  integer                                :: slicelocal
  integer                                :: i, ii, j, jj, proc_i, itranspose, jtranspose, i2r, component

  type(C_PTR)                            :: inputarray, copytoarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
!  slicelocal = zslice - z_min + 1
  slicelocal = islice - i_min_3w + 1

  !> Transpose from jk to kj and do FFT in k direction.
  !! Copy the result of the FFT to the buffer array

  !> u4 component
  inputarray = c_loc(arr_u4_1r2w)
  copytoarray = c_loc(arr_u4_2w)
  call f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u5 component
  inputarray = c_loc(arr_u5_1r2w)
  copytoarray = c_loc(arr_u5_2w)
  call f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u6 component
  inputarray = c_loc(arr_u6_1r2w)
  copytoarray = c_loc(arr_u6_2w)
  call f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

end subroutine f_RK3TransformNonLinearkRs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformCouplingTermkRs2Fs
!> @param ierr should return 0
subroutine f_RK3TransformCouplingTermkRs2Fs ( islice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use f_arrays,     only: i_min_3w, arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: islice

  integer                                :: slicelocal
  integer                                :: i, ii, j, jj, proc_i, itranspose, jtranspose, i2r, component

  type(C_PTR)                            :: inputarray, copytoarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
!  slicelocal = zslice - z_min + 1
  slicelocal = islice - i_min_3w + 1

  !> Transpose from jk to kj and do FFT in k direction.
  !! Copy the result of the FFT to the buffer array

  !> u1 component
  inputarray = c_loc(arr_u1_1r2w)
  copytoarray = c_loc(arr_u1_2w)
  call f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u2 component
  inputarray = c_loc(arr_u2_1r2w)
  copytoarray = c_loc(arr_u2_2w)
  call f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  !> u3 component
  inputarray = c_loc(arr_u3_1r2w)
  copytoarray = c_loc(arr_u3_2w)
  call f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

end subroutine f_RK3TransformCouplingTermkRs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformkRs2Fs
!> @param ierr should return 0
  subroutine f_RK3TransformkRs2Fs ( inputarray, copytoarray, slicelocal, ierr )

  use g_domain,     only: KK_MAX
  use g_parameters, only: NODES_Z
  use g_constants,  only: ZERO
  use f_arrays,     only: k_width_3w, j_width_3w, &
                          z_width_1r2w, j_width_1r2w, i_width_1r2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w
  use f_fftw,       only: f_Fftw_k_Rs2Fs, arr_jkkj_transposein, arr_k_out_da

  implicit none
 
  type(C_PTR), intent(in)                :: inputarray, copytoarray
  integer,intent(in)                     :: slicelocal
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                          :: j, k, kk, proc_j, jtranspose, &
                                                      ktranspose, j1r2w
  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:)   :: fromarray
  real(kind=C_DOUBLE),pointer,dimension(:,:,:)     :: targetarray

  !> Link inputarray to the C pointer of the vorticity component
  call c_f_pointer(inputarray, fromarray, [2,j_width_1r2w,z_width_1r2w,i_width_1r2w])

!  write(*,*) 'arr_jkkj_transposein: ', shape(arr_jkkj_transposein)
!  write(*,*) 'fromarray:            ', shape(fromarray)

! Copy component from the buffer array to the transpose input array
  DOZ: do k = 1, z_width_1r2w, 1                     ! parallelised in this direction
    j1r2w = 0
    DOPROCJ: do proc_j = 0, nprocsj_3w - 1, 1 
      DOJLOCAL: do j = 1, j_allwidths_3w(proc_j), 1

        j1r2w = j1r2w + 1
        jtranspose = (proc_j * j_maxwidth_3w) + j

        arr_jkkj_transposein(:, jtranspose, k) = fromarray(:, j1r2w, k, slicelocal)

      end do DOJLOCAL
    end do DOPROCJ
  end do DOZ

  !> Perform FFT in k direction
  call f_Fftw_k_Rs2Fs ( ierr )

  !> Link targetarray to the C pointer of the vorticity component
  call c_f_pointer(copytoarray, targetarray, [2, k_width_3w, j_width_3w])

!  write(*,*) 'arr_k_out_da: ', shape(arr_k_out_da)
!  write(*,*) 'targetarray:  ', shape(targetarray)


  DOJ3W: do j = 1, j_width_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOKPOSITIVE: do k = 1, KK_MAX + 1, 1
      !> Copy to buffer
      targetarray(:,k,j) = arr_k_out_da(:,kk,j)
      kk = kk + 1
    end do DOKPOSITIVE

    !> Get rid of the highest 1/3 wavenumbers, both positive and negative
     kk = NODES_Z - KK_MAX + 1

    !> Copy negative wavenumbers
    DOKNEGATIVE: do k = KK_MAX + 2, k_width_3w
      !> Copy to buffer
      targetarray(:,k,j) = arr_k_out_da(:,kk,j)
      kk = kk + 1
    end do DOKNEGATIVE

  end do DOJ3W

end subroutine f_RK3TransformkRs2Fs
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3FormFIsotropic
!> @param ierr should return 0
subroutine f_RK3FormFIsotropic ( i, stage, dt, ierr )

  use g_constants, only : ZERO, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                          ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, NODES_KJ, NODES_KK
  use g_parameters, only: NODES_Z, F_TYPE, F_ISOTROPIC, F_NU, FORCE_KF, &
                          P_TRACK_PART, P_TWO_WAY, YES
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w, k_max_3w, z_max, j_max_1r2w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w, &
                          f_ArraysWavenumber
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_k_in_da, arr_kjjk_transposeout
  use f_force,      only: fluid_f1, fluid_f2, fluid_f3, FORCE_NK

  implicit none
 
  integer,intent(in)                     :: i, stage
  real(kind=C_DOUBLE),intent(in)         :: dt
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                :: j, k 
  real(kind=C_DOUBLE)                    :: wavenumber_k, wavenumber_j, wavenumber_i
  real(kind=C_DOUBLE)                    :: kmag2   
  real(kind=C_DOUBLE),dimension(0:1)     :: k_dot_b, k_dot_nl

  real(kind=C_DOUBLE),dimension(0:1,1:3) :: viscous
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: pressure
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: nonlinear
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: coupling 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: force
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: fluid_f
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: totalforce

  real(kind=C_DOUBLE)                    :: rkbeta, rkgamma

  coupling = ZERO

  !> Compute i component of wave vector
  wavenumber_i = f_ArraysWavenumber ( i, KI_MAX, K0I )

  DOU1J: do j = j_min_3w, j_max_3w, 1

    !> Compute j component of wave vector
    wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )

    !> Positive wavenumbers
    DOU1K: do k = k_min_3w, k_max_3w, 1

      !> Compute k component of wave vector
      wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )

      !> compute scalar wave number
      kmag2 = ( wavenumber_i * wavenumber_i ) + &
              ( wavenumber_j * wavenumber_j ) + &
              ( wavenumber_k * wavenumber_k )

      !> procede only if wave number greater than zero
      IFZERO: if ( kmag2 > ZERO ) then

!------------------ viscous term ---------------------------------------------------------------------  
        viscous(:,1)   = - F_NU * kmag2 * arr_u1_3w(:,k,j,i)    
        viscous(:,2)   = - F_NU * kmag2 * arr_u2_3w(:,k,j,i)    
        viscous(:,3)   = - F_NU * kmag2 * arr_u3_3w(:,k,j,i)     
!-----------------------------------------------------------------------------------------------------                   

!------------------ nonlinear term -------------------------------------------------------------------             
        nonlinear(:,1) =  arr_u4_2w(:,k,j)
        nonlinear(:,2) =  arr_u5_2w(:,k,j)
        nonlinear(:,3) =  arr_u6_2w(:,k,j)
!------------------------------------------------------------------------------------------------------
          
!------------------------- two-way coupling term -------------------------------------------------------
        IFPARTICLES: if ( P_TRACK_PART == YES ) then
          IFTWOWAY: if ( P_TWO_WAY == YES ) then
            coupling(:,1) = arr_u1_2w(:,k,j)
            coupling(:,2) = arr_u2_2w(:,k,j)
            coupling(:,3) = arr_u3_2w(:,k,j)
          end if IFTWOWAY
        end if IFPARTICLES
!-------------------------------------------------------------------------------------------------------

!> Force terms for isotropic turbulence
!------------------ pressure term ---------------------------------------------------------------------
        k_dot_nl   = (wavenumber_i * nonlinear(:,1)) + &
                     (wavenumber_j * nonlinear(:,2)) + &
                     (wavenumber_k * nonlinear(:,3))

        pressure(:,1)  = ( wavenumber_i / kmag2 ) * k_dot_nl
        pressure(:,2)  = ( wavenumber_j / kmag2 ) * k_dot_nl
        pressure(:,3)  = ( wavenumber_k / kmag2 ) * k_dot_nl

!----------------------- forcing terms ----------------------------------------------------------------
        !> Zero if not forced
        force = ZERO

        IFFORCED: if ( sqrt(kmag2) <= FORCE_KF ) then
          IFFNOTZERO: if ( sqrt(kmag2) > ZERO ) then

            !> Read the correct forcing from forcing array
            IFJ: if ( j <= FORCE_NK ) then
              IFJPOSK: if ( k <= FORCE_NK ) then
                fluid_f(:,1) = fluid_f1(:,k,j,i)
                fluid_f(:,2) = fluid_f2(:,k,j,i)
                fluid_f(:,3) = fluid_f3(:,k,j,i)
!                write(*,*) fluid_f1(:,k,j,i), fluid_f3(:,k,j,i), fluid_f3(:,k,j,i)
              else
                fluid_f(:,1) = fluid_f1(:,k-NODES_KK,j,i)
                fluid_f(:,2) = fluid_f2(:,k-NODES_KK,j,i)
                fluid_f(:,3) = fluid_f3(:,k-NODES_KK,j,i)
!                write(*,*) fluid_f1(:,k-NODES_KK,j,i), fluid_f2(:,k-NODES_KK,j,i), fluid_f3(:,k-NODES_KK,j,i)
              end if IFJPOSK      
            else
              IFJNEGK: if ( k <= FORCE_NK ) then
                fluid_f(:,1) = fluid_f1(:,k,j-NODES_KJ,i)
                fluid_f(:,2) = fluid_f2(:,k,j-NODES_KJ,i)
                fluid_f(:,3) = fluid_f3(:,k,j-NODES_KJ,i)
              else
                fluid_f(:,1) = fluid_f1(:,k-NODES_KK,j-NODES_KJ,i)
                fluid_f(:,2) = fluid_f2(:,k-NODES_KK,j-NODES_KJ,i)
                fluid_f(:,3) = fluid_f3(:,k-NODES_KK,j-NODES_KJ,i)
              end if IFJNEGK      
            end if IFJ

            k_dot_b  = ( wavenumber_i * fluid_f(:,1) + &
                         wavenumber_j * fluid_f(:,2) + &
                         wavenumber_k * fluid_f(:,3) )

!            write(*,*) 'Forcing before continuity fix on ', i, j, k, ': ', fluid_f(:,1), &
!                       fluid_f(:,2), fluid_f(:,3)

!            write(*,*) 'Wavenumbers: ', wavenumber_i, wavenumber_j, wavenumber_k, &
!                       'k_dot_b: ', k_dot_b

            force(:,1) = fluid_f(:,1) - ( (wavenumber_i / kmag2) * k_dot_b )
            force(:,2) = fluid_f(:,2) - ( (wavenumber_j / kmag2) * k_dot_b )
            force(:,3) = fluid_f(:,3) - ( (wavenumber_k / kmag2) * k_dot_b )
                 
!            write(*,*) 'Forcing on ', i, j, k, ': ', force

          end if IFFNOTZERO
        end if IFFORCED

        !> combine terms ( navier-stokes in spectral form )
        totalforce = viscous - nonlinear + pressure + force + coupling

        call f_RK3RungeKutta ( arr_rk3_3w_1(:,k,j,i), arr_rk3_3w_2(:,k,j,i), &
                               arr_rk3_3w_3(:,k,j,i), arr_u1_3w(:,k,j,i), &
                               arr_u2_3w(:,k,j,i), arr_u3_3w(:,k,j,i), &
                               totalforce, stage, dt )

        call f_RK3ContinuityFixIsotropic ( wavenumber_i, wavenumber_j, wavenumber_k, & 
                                           kmag2, arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                           arr_u3_3w(:,k,j,i), ierr )

      end if IFZERO

    end do DOU1K
  end do DOU1J

end subroutine f_RK3FormFIsotropic
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3FormFBFRot
!> @param ierr should return 0
subroutine f_RK3FormFBFRot (  i, stage, dt, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, NODES_KJ, NODES_KK
  use g_parameters, only: NODES_Z, F_TYPE, F_SHEAR_DNS, F_SHEAR_RDT_VISCOUS, &
                          F_SHEAR_RDT_INVISCID, F_ISOTROPIC_DEFORMED, F_NU, FORCE_KF, &
                          P_TRACK_PART, P_TWO_WAY, YES
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w, k_max_3w, z_max, j_max_1r2w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w, &
                          f_ArraysWavenumber
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_k_in_da, arr_kjjk_transposeout
  use f_force,      only: fluid_f1, fluid_f2, fluid_f3, FORCE_NK

  implicit none

  integer,intent(in)                     :: i, stage
  real(kind=C_DOUBLE),intent(in)         :: dt
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                :: j, k 
  integer                                :: ider, ider1, ider2, ii, jj, kk, ivel
  real(kind=C_DOUBLE),dimension(1:3)     :: wavenumber
  real(kind=C_DOUBLE)                    :: kmag2   
  real(kind=C_DOUBLE),dimension(0:1)     :: k_dot_b, k_dot_nl, rhsg, rhsn, pres

  real(kind=C_DOUBLE),dimension(0:1,1:3) :: localvelocity 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: viscous
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: pressure
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: nonlinear
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: coupling 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: homogeneous 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: forcing
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: totalforce

  real(kind=C_DOUBLE)                    :: rkbeta, rkgamma

  coupling = ZERO

  !> Compute i component of wave vector
  wavenumber(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )

  DOU1J: do j = j_min_3w, j_max_3w, 1

    !> Compute j component of wave vector
    wavenumber(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )

    !> Positive wavenumbers
    DOU1K: do k = k_min_3w, k_max_3w, 1

      !> Compute k component of wave vector
      wavenumber(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

      !> compute scalar wave number 
      kmag2 = ( wavenumber(1) * wavenumber(1) ) + &
              ( wavenumber(2) * wavenumber(2) ) + &
              ( wavenumber(3) * wavenumber(3) )

      !> procede only if wave number greater than zero
      IFZERO: if ( kmag2 > ZERO ) then

        IFNONLINEAR: if ( F_TYPE == F_SHEAR_DNS ) then
        !------------------ nonlinear term ---------------------------------------------------------
          call f_RK3BFRotNonlinear ( arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), &
                                       arr_u6_2w(:,k,j), nonlinear )
!          nonlinear(:,1) =  arr_u4_2w(:,k,j)
!          nonlinear(:,2) =  arr_u5_2w(:,k,j)
!          nonlinear(:,3) =  arr_u6_2w(:,k,j)
        else if ( F_TYPE == F_ISOTROPIC_DEFORMED ) then
          call f_RK3BFRotNonlinear ( arr_u4_2w(:,k,j), arr_u5_2w(:,k,j), &
                                     arr_u6_2w(:,k,j), nonlinear )
!          nonlinear(:,1) =  arr_u4_2w(:,k,j)
!          nonlinear(:,2) =  arr_u5_2w(:,k,j)
!          nonlinear(:,3) =  arr_u6_2w(:,k,j)
        else
          nonlinear(:,1) = ZERO  
          nonlinear(:,2) = ZERO
          nonlinear(:,3) = ZERO
        end if IFNONLINEAR

!        if ( i == 0 ) then

!          write(*,*) ' wavenumber(:), nonlinear(:,:): ', wavenumber(:), nonlinear(:,:)

!        end if 

!------------------ pressure term --------------------------------------------------------------------
        call f_RK3BFRotPressure ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                    arr_u3_3w(:,k,j,i), wavenumber, nonlinear, pressure )
!-----------------------------------------------------------------------------------------------------

!------------------ two-way coupling term ------------------------------------------------------------
        IFPARTICLES: if ( P_TRACK_PART == YES ) then
          IFTWOWAY: if ( P_TWO_WAY == YES ) then
            coupling(:,1) = arr_u1_2w(:,k,j)
            coupling(:,2) = arr_u2_2w(:,k,j)
            coupling(:,3) = arr_u3_2w(:,k,j)
          end if IFTWOWAY
        end if IFPARTICLES
!-------------------------------------------------------------------------------------------------------

!------------------ homogeneous source term ------------------------------------------------------------

        call f_RK3BFRotSource ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                  arr_u3_3w(:,k,j,i), homogeneous )
!-------------------------------------------------------------------------------------------------------

        !> the remaining terms are different depending on the method used, currently DNS, inviscid RDT
        !! or viscous RDT
        select case ( F_TYPE )
        !-------------------------------
          case ( F_SHEAR_DNS ) 
            !------------------ viscous term -----------------------------------------------------------
            call f_RK3BFRotViscous ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                       arr_u3_3w(:,k,j,i), wavenumber, viscous )
            !-------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------
            !> combine terms ( navier-stokes in spectral form )
            totalforce = viscous - nonlinear + pressure + coupling - homogeneous

        !-------------------------------
          case ( F_SHEAR_RDT_VISCOUS )
            !------------------ viscous term -----------------------------------------------------------
            call f_RK3BFRotViscous ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                       arr_u3_3w(:,k,j,i), wavenumber, viscous )
            !-------------------------------------------------------------------------------------------

            !> combine terms ( navier-stokes in spectral form )
            totalforce = viscous + pressure + coupling - homogeneous

        !-------------------------------
          case ( F_SHEAR_RDT_INVISCID )

            !> combine terms ( navier-stokes in spectral form )
!            totalforce = pressure + coupling - homogeneous
            totalforce = pressure - homogeneous

        !-------------------------------
          case ( F_ISOTROPIC_DEFORMED ) 
            !------------------ viscous term -----------------------------------------------------------
            call f_RK3BFRotViscous ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                       arr_u3_3w(:,k,j,i), wavenumber, viscous )
            !-------------------------------------------------------------------------------------------

            call f_RK3BFRotForcing ( wavenumber, forcing, ierr )
            !-------------------------------------------------------------------------------------------
            !> combine terms ( navier-stokes in spectral form )
            totalforce = viscous - nonlinear + pressure + coupling + forcing

        !-------------------------------
        end select
        !-------------------------------

        !> Now do the Runge-Kutta timestep
        call f_RK3RungeKutta ( arr_rk3_3w_1(:,k,j,i), arr_rk3_3w_2(:,k,j,i), &
                               arr_rk3_3w_3(:,k,j,i), arr_u1_3w(:,k,j,i), &
                               arr_u2_3w(:,k,j,i), arr_u3_3w(:,k,j,i), &
                               totalforce, stage, dt )

        !> the continuity equation in the moving coordinate system is identical to isotropic turbulence
!        call f_RK3ContinuityFixBFRot ( wavenumber(1), wavenumber(2), wavenumber(3), & 
!                                         arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
!                                         arr_u3_3w(:,k,j,i), ierr )
        call f_RK3ContinuityFixIsotropic ( wavenumber(1), wavenumber(2), wavenumber(3), & 
                                           kmag2, arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                           arr_u3_3w(:,k,j,i), ierr )


      end if IFZERO

    end do DOU1K
  end do DOU1J
                 
end subroutine f_RK3FormFBFRot
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3FormFRRot
!> @param ierr should return 0
subroutine f_RK3FormFRRot (  i, stage, dt, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, NODES_KJ, NODES_KK
  use g_parameters, only: NODES_Z, F_TYPE, F_SHEAR_DNS, F_SHEAR_RDT_VISCOUS, &
                          F_SHEAR_RDT_INVISCID, F_ISOTROPIC_DEFORMED, F_NU, FORCE_KF, &
                          P_TRACK_PART, P_TWO_WAY, YES
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w, k_max_3w, z_max, j_max_1r2w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_rk3_3w_1, arr_rk3_3w_2, arr_rk3_3w_3, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u1_2w, arr_u2_2w, arr_u3_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w, &
                          f_ArraysWavenumber
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_k_in_da, arr_kjjk_transposeout
  use f_force,      only: fluid_f1, fluid_f2, fluid_f3, FORCE_NK

  implicit none

  integer,intent(in)                     :: i, stage
  real(kind=C_DOUBLE),intent(in)         :: dt
  PetscErrorCode,intent(inout)                  :: ierr

  integer                                :: j, k 
  integer                                :: ider, ider1, ider2, ii, jj, kk, ivel
  real(kind=C_DOUBLE),dimension(1:3)     :: wavenumber
  real(kind=C_DOUBLE)                    :: kmag2   
  real(kind=C_DOUBLE),dimension(0:1)     :: k_dot_b, k_dot_nl, rhsg, rhsn, pres

  real(kind=C_DOUBLE),dimension(0:1,1:3) :: localvelocity 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: viscous
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: pressure
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: nonlinear
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: coupling 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: homogeneous 
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: forcing
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: totalforce

  real(kind=C_DOUBLE)                    :: rkbeta, rkgamma

  coupling = ZERO

  !> Compute i component of wave vector
  wavenumber(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )

  DOU1J: do j = j_min_3w, j_max_3w, 1

    !> Compute j component of wave vector
    wavenumber(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )

    !> Positive wavenumbers
    DOU1K: do k = k_min_3w, k_max_3w, 1

      !> Compute k component of wave vector
      wavenumber(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

      !> compute scalar wave number 
      kmag2 = ( wavenumber(1) * wavenumber(1) ) + &
              ( wavenumber(2) * wavenumber(2) ) + &
              ( wavenumber(3) * wavenumber(3) )

      !> procede only if wave number greater than zero
      IFZERO: if ( kmag2 > ZERO ) then

        IFNONLINEAR: if ( F_TYPE == F_SHEAR_DNS ) then
        !------------------ nonlinear term ---------------------------------------------------------
          nonlinear(:,1) =  arr_u4_2w(:,k,j)
          nonlinear(:,2) =  arr_u5_2w(:,k,j)
          nonlinear(:,3) =  arr_u6_2w(:,k,j)
        else if ( F_TYPE == F_ISOTROPIC_DEFORMED ) then
          nonlinear(:,1) =  arr_u4_2w(:,k,j)
          nonlinear(:,2) =  arr_u5_2w(:,k,j)
          nonlinear(:,3) =  arr_u6_2w(:,k,j)
        else
          nonlinear(:,1) = ZERO  
          nonlinear(:,2) = ZERO
          nonlinear(:,3) = ZERO
        end if IFNONLINEAR

!        if ( i == 0 ) then

!          write(*,*) ' wavenumber(:), nonlinear(:,:): ', wavenumber(:), nonlinear(:,:)

!        end if 

!------------------ pressure term --------------------------------------------------------------------
        call f_RK3RRotPressure ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                    arr_u3_3w(:,k,j,i), wavenumber, nonlinear, pressure )
!-----------------------------------------------------------------------------------------------------

!------------------ two-way coupling term ------------------------------------------------------------
        IFPARTICLES: if ( P_TRACK_PART == YES ) then
          IFTWOWAY: if ( P_TWO_WAY == YES ) then
            coupling(:,1) = arr_u1_2w(:,k,j)
            coupling(:,2) = arr_u2_2w(:,k,j)
            coupling(:,3) = arr_u3_2w(:,k,j)
          end if IFTWOWAY
        end if IFPARTICLES
!-------------------------------------------------------------------------------------------------------

!------------------ homogeneous source term ------------------------------------------------------------

        call f_RK3RRotSource ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                  arr_u3_3w(:,k,j,i), homogeneous )
!-------------------------------------------------------------------------------------------------------

        !> the remaining terms are different depending on the method used, currently DNS, inviscid RDT
        !! or viscous RDT
        select case ( F_TYPE )
        !-------------------------------
          case ( F_SHEAR_DNS ) 
            !------------------ viscous term -----------------------------------------------------------
            call f_RK3RRotViscous ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                       arr_u3_3w(:,k,j,i), wavenumber, viscous )
            !-------------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------------
            !> combine terms ( navier-stokes in spectral form )
            totalforce = viscous - nonlinear + pressure + coupling - homogeneous

        !-------------------------------
          case ( F_SHEAR_RDT_VISCOUS )
            !------------------ viscous term -----------------------------------------------------------
            call f_RK3RRotViscous ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                       arr_u3_3w(:,k,j,i), wavenumber, viscous )
            !-------------------------------------------------------------------------------------------

            !> combine terms ( navier-stokes in spectral form )
            totalforce = viscous + pressure + coupling - homogeneous

        !-------------------------------
          case ( F_SHEAR_RDT_INVISCID )

            !> combine terms ( navier-stokes in spectral form )
!            totalforce = pressure + coupling - homogeneous
            totalforce = pressure - homogeneous

        !-------------------------------
          case ( F_ISOTROPIC_DEFORMED ) 
            !------------------ viscous term -----------------------------------------------------------
            call f_RK3RRotViscous ( arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                       arr_u3_3w(:,k,j,i), wavenumber, viscous )
            !-------------------------------------------------------------------------------------------

            call f_RK3RRotForcing ( wavenumber, forcing, ierr )
            !-------------------------------------------------------------------------------------------
            !> combine terms ( navier-stokes in spectral form )

            totalforce = viscous - nonlinear + pressure + coupling + forcing

        !-------------------------------
        end select
        !-------------------------------

        !> Now do the Runge-Kutta timestep
        call f_RK3RungeKutta ( arr_rk3_3w_1(:,k,j,i), arr_rk3_3w_2(:,k,j,i), &
                               arr_rk3_3w_3(:,k,j,i), arr_u1_3w(:,k,j,i), &
                               arr_u2_3w(:,k,j,i), arr_u3_3w(:,k,j,i), &
                               totalforce, stage, dt )

        !> the continuity equation in the moving coordinate system is identical to isotropic turbulence
        call f_RK3ContinuityFixRRot ( wavenumber(1), wavenumber(2), wavenumber(3), & 
                                         arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
                                         arr_u3_3w(:,k,j,i), ierr )
!        call f_RK3ContinuityFixIsotropic ( wavenumber(1), wavenumber(2), wavenumber(3), & 
!                                           kmag2, arr_u1_3w(:,k,j,i), arr_u2_3w(:,k,j,i), &
!                                           arr_u3_3w(:,k,j,i), ierr )


      end if IFZERO

    end do DOU1K
  end do DOU1J
                 
end subroutine f_RK3FormFRRot
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3BFRotNonlinear
!> Compute the viscous term
!> @param ierr should return 0
pure subroutine f_RK3BFRotNonlinear ( nl1, nl2, nl3, nonlinearforce )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, bvalexist
  use g_parameters, only : F_NU

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: nl1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: nl2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: nl3
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: nonlinearforce

  integer :: i, r

  nonlinearforce(:,1) = ZERO 
  nonlinearforce(:,2) = ZERO
  nonlinearforce(:,3) = ZERO

  DOI: do i = 1, 3, 1
    IFBI1: if ( bvalexist(i,1) .eqv. .TRUE. ) then
      nonlinearforce(:,i) = nonlinearforce(:,i) + bmat(i,1) * nl1
    end if IFBI1
    IFBI2: if ( bvalexist(i,2) .eqv. .TRUE. ) then
      nonlinearforce(:,i) = nonlinearforce(:,i) + bmat(i,2) * nl2
    end if IFBI2
    IFBI3: if ( bvalexist(i,3) .eqv. .TRUE. ) then
      nonlinearforce(:,i) = nonlinearforce(:,i) + bmat(i,3) * nl3
    end if IFBI3
  end do DOI

end subroutine f_RK3BFRotNonlinear
!------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3BFRotViscous
!> Compute the viscous term
!> @param ierr should return 0
pure subroutine f_RK3BFRotViscous ( u1, u2, u3, wvn, viscousforce )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, bvalexist
  use g_parameters, only : F_NU

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u3
  real(kind=C_DOUBLE),dimension(1:3),intent(in)       :: wvn
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: viscousforce

  integer :: l,j,n

  viscousforce(:,1) = ZERO 
  viscousforce(:,2) = ZERO
  viscousforce(:,3) = ZERO

  DOL: do l = 1, 3, 1
    DOJ: do j = 1, 3, 1
      IFBLJ: if ( bvalexist(l,j) .eqv. .TRUE. ) then
        DON: do n = 1, 3, 1
          IFBNJ: if ( bvalexist(n,j) .eqv. .TRUE. ) then
            viscousforce(:,1) = viscousforce(:,1) - F_NU * bmat(l,j) * bmat(n,j) * wvn(l) * wvn(n) * u1
            viscousforce(:,2) = viscousforce(:,2) - F_NU * bmat(l,j) * bmat(n,j) * wvn(l) * wvn(n) * u2
            viscousforce(:,3) = viscousforce(:,3) - F_NU * bmat(l,j) * bmat(n,j) * wvn(l) * wvn(n) * u3
          end if IFBNJ
        end do DON
      end if IFBLJ
    end do DOJ
  end do DOL

end subroutine f_RK3BFRotViscous
!------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RRotViscous
!> Compute the viscous term
!> @param ierr should return 0
pure subroutine f_RK3RRotViscous ( u1, u2, u3, wvn, viscousforce )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, bvalexist
  use g_parameters, only : F_NU

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u3
  real(kind=C_DOUBLE),dimension(1:3),intent(in)       :: wvn
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: viscousforce

  integer :: l,j,n

  viscousforce(:,1) = ZERO 
  viscousforce(:,2) = ZERO
  viscousforce(:,3) = ZERO

  DOL: do l = 1, 3, 1
    DOJ: do j = 1, 3, 1
      IFBLJ: if ( bvalexist(l,j) .eqv. .TRUE. ) then
        DON: do n = 1, 3, 1
          IFBNJ: if ( bvalexist(n,j) .eqv. .TRUE. ) then
            viscousforce(:,1) = viscousforce(:,1) - F_NU * bmat(l,j) * bmat(n,j) * wvn(l) * wvn(n) * u1
            viscousforce(:,2) = viscousforce(:,2) - F_NU * bmat(l,j) * bmat(n,j) * wvn(l) * wvn(n) * u2
            viscousforce(:,3) = viscousforce(:,3) - F_NU * bmat(l,j) * bmat(n,j) * wvn(l) * wvn(n) * u3
          end if IFBNJ
        end do DON
      end if IFBLJ
    end do DOJ
  end do DOL

end subroutine f_RK3RRotViscous
!------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3BFRotPressure
!> Compute the pressure term
!> @param ierr should return 0
pure subroutine f_RK3BFRotPressure ( u1, u2, u3, wvn, nl, pressureforce )

  use g_constants, only : ZERO, TWO
  use g_domain,    only : amat, bmat, binv, avalexist, bvalexist

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u3
  real(kind=C_DOUBLE),dimension(1:3),intent(in)       :: wvn
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(in)   :: nl
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: pressureforce

  integer                            :: m, r, s, k, i, p, j

  !> pressureforce = ( p1 / p4 ) * ( p2 + p3 )
!  real(kind=C_DOUBLE),dimension(0:1,1:3) :: p2, p3
  real(kind=C_DOUBLE),dimension(0:1)     :: p2, p3
  real(kind=C_DOUBLE)                    :: p4
  real(kind=C_DOUBLE),dimension(1:3)     :: p1
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: u

  p1 = ZERO
  p2 = ZERO
  p3 = ZERO
  p4 = ZERO

  u(:,1) = u1
  u(:,2) = u2
  u(:,3) = u3

  DOM: do m = 1, 3, 1

!    p2(:,m) =  wvn(m) * nl(:,m) 
    p2(:) = p2(:) + ( wvn(m) * nl(:,m) )

    DOR: do r = 1, 3, 1
      IFBMR: if ( bvalexist(m,r) .eqv. .TRUE. ) then
        DOS: do s = 1, 3, 1
          IFARS: if ( avalexist(r,s) .eqv. .TRUE. ) then
            DOKP3: do k = 1, 3, 1
              IFBSK: if ( bvalexist(s,k) .eqv. .TRUE. ) then
!                p3(:,m) = TWO * ( wvn(m) * bmat(m,r) * amat(r,s) * binv(s,k) * u(:,k) )
                p3(:) = p3(:) + ( TWO * ( wvn(m) * bmat(m,r) * amat(r,s) * binv(s,k) * u(:,k) ) )
              end if IFBSK
            end do DOKP3
          end if IFARS
        end do DOS
        DOKP4: do k = 1, 3, 1
          IFBKR: if ( bvalexist(k,r) .eqv. .TRUE. ) then
            p4 = p4 + ( wvn(m) * bmat(m,r) * bmat(k,r) * wvn(k) )
          end if IFBKR
        end do DOKP4      
      end if IFBMR
    end do DOR

  end do DOM

  DOI: do i = 1, 3, 1
    DOP: do p = 1, 3, 1
      IFBIP: if ( bvalexist(i,p) .eqv. .TRUE. ) then
        DOJ: do j = 1, 3, 1
          IFBJP: if ( bvalexist(j,p) .eqv. .TRUE. ) then
            p1(i) = p1(i) + bmat(i,p) * bmat(j,p) * wvn(j)
          end if IFBJP
        end do DOJ
      end if IFBIP
    end do DOP
  end do DOI

  IFP4: if ( p4 > ZERO ) then
!    pressureforce(:,1) = ( p1(1) / p4 ) * ( p2(:,1) + p3(:,1) )
!    pressureforce(:,2) = ( p1(2) / p4 ) * ( p2(:,2) + p3(:,2) )
!    pressureforce(:,3) = ( p1(3) / p4 ) * ( p2(:,3) + p3(:,3) )
    pressureforce(:,1) = ( p1(1) / p4 ) * ( p2(:) + p3(:) )
    pressureforce(:,2) = ( p1(2) / p4 ) * ( p2(:) + p3(:) )
    pressureforce(:,3) = ( p1(3) / p4 ) * ( p2(:) + p3(:) )
  else
    pressureforce(:,1) = ZERO 
    pressureforce(:,2) = ZERO
    pressureforce(:,3) = ZERO
  end if IFP4

end subroutine f_RK3BFRotPressure
!------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RRotPressure
!> Compute the pressure term
!> @param ierr should return 0
pure subroutine f_RK3RRotPressure ( u1, u2, u3, wvn, nl, pressureforce )

  use g_constants, only : ZERO, TWO
  use g_domain,    only : amat, bmat, binv, avalexist, bvalexist

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u3
  real(kind=C_DOUBLE),dimension(1:3),intent(in)       :: wvn
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(in)   :: nl
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: pressureforce

  integer                            :: m, r, s, k, i, p, j

  !> pressureforce = ( p1 / p4 ) * ( p2 + p3 )
!  real(kind=C_DOUBLE),dimension(0:1,1:3) :: p2, p3
  real(kind=C_DOUBLE),dimension(0:1)     :: p2, p3
  real(kind=C_DOUBLE)                    :: p4
  real(kind=C_DOUBLE),dimension(1:3)     :: p1
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: u

  p1 = ZERO
  p2 = ZERO
  p3 = ZERO
  p4 = ZERO

  u(:,1) = u1
  u(:,2) = u2
  u(:,3) = u3

  DOM: do m = 1, 3, 1

!    p2(:,m) =  wvn(m) * nl(:,m) 

    DOR: do r = 1, 3, 1
      IFBMR: if ( bvalexist(m,r) .eqv. .TRUE. ) then
        p2(:) = p2(:) + ( wvn(m) * bmat(m,r) * nl(:,r) )
        DOKP3: do k = 1, 3, 1
          IFARK: if ( avalexist(r,k) .eqv. .TRUE. ) then
!                p3(:,m) = TWO * ( wvn(m) * bmat(m,r) * amat(r,s) * binv(s,k) * u(:,k) )
            p3(:) = p3(:) + ( &
!                       TWO * &
                    ( wvn(m) * bmat(m,r) * amat(r,k) * u(:,k) ) )
          end if IFARK
        end do DOKP3
        DOKP4: do k = 1, 3, 1
          IFBKR: if ( bvalexist(k,r) .eqv. .TRUE. ) then
            p4 = p4 + ( wvn(m) * bmat(m,r) * wvn(k) * bmat(k,r) )
          end if IFBKR
        end do DOKP4      
      end if IFBMR
    end do DOR

  end do DOM

  DOI: do i = 1, 3, 1
    DOP: do p = 1, 3, 1
      IFBPI: if ( bvalexist(p,i) .eqv. .TRUE. ) then
        p1(i) = p1(i) +  wvn(p) * bmat(p,i)
      end if IFBPI
    end do DOP
  end do DOI

  IFP4: if ( p4 > ZERO ) then
!    pressureforce(:,1) = ( p1(1) / p4 ) * ( p2(:,1) + p3(:,1) )
!    pressureforce(:,2) = ( p1(2) / p4 ) * ( p2(:,2) + p3(:,2) )
!    pressureforce(:,3) = ( p1(3) / p4 ) * ( p2(:,3) + p3(:,3) )
    pressureforce(:,1) = ( p1(1) / p4 ) * ( p2(:) + p3(:) )
    pressureforce(:,2) = ( p1(2) / p4 ) * ( p2(:) + p3(:) )
    pressureforce(:,3) = ( p1(3) / p4 ) * ( p2(:) + p3(:) )
  else
    pressureforce(:,1) = ZERO 
    pressureforce(:,2) = ZERO
    pressureforce(:,3) = ZERO
  end if IFP4

end subroutine f_RK3RRotPressure
!------------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------------

!---------------------- homogeneous source term --------------------------------------------------------
!---------------------------------------------------------------------------------
! subroutine f_RK3BFRotSource
!> Compute the homogeneous source term
!> @param ierr should return 0
pure subroutine f_RK3BFRotSource ( u1, u2, u3, homogeneousforce )

  use g_constants,  only : ZERO, TWO
  use g_domain,     only : amat, bmat, binv, avalexist, bvalexist

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u3
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: homogeneousforce

  integer :: i, r, j

  homogeneousforce = ZERO

  DOI: do i = 1, 3, 1
    DOR: do r = 1, 3, 1
      IFBIR: if ( bvalexist(i,r) .eqv. .TRUE. ) then
        DOJ: do j = 1, 3, 1
          IFARJ: if ( avalexist(r,j) .eqv. .TRUE. ) then
            IFBJ1: if ( bvalexist(j,1) .eqv. .TRUE. ) then
              homogeneousforce(:,i) = homogeneousforce(:,i) + &
                                      ( TWO * &
                                        ( bmat(i,r) * amat(r,j) * binv(j,1) * u1(:) ) )
            end if IFBJ1
            IFBJ2: if ( bvalexist(j,2) .eqv. .TRUE. ) then
              homogeneousforce(:,i) = homogeneousforce(:,i) + &
                                      ( TWO * &
                                        ( bmat(i,r) * amat(r,j) * binv(j,2) * u2(:) ) )
            end if IFBJ2
            IFBJ3: if ( bvalexist(j,3) .eqv. .TRUE. ) then
              homogeneousforce(:,i) = homogeneousforce(:,i) + &
                                      ( TWO * &
                                        ( bmat(i,r) * amat(r,j) * binv(j,3) * u3(:) ) )
            end if IFBJ3
          end if IFARJ
        end do DOJ
      end if IFBIR
    end do DOR
  end do DOI

end subroutine f_RK3BFRotSource
!-------------------------------------------------------------------------------------------------------


!---------------------- homogeneous source term --------------------------------------------------------
!---------------------------------------------------------------------------------
! subroutine f_RK3RRotSource
!> Compute the homogeneous source term
!> @param ierr should return 0
pure subroutine f_RK3RRotSource ( u1, u2, u3, homogeneousforce )

  use g_constants,  only : ZERO, TWO
  use g_domain,     only : amat, bmat, binv, avalexist, bvalexist

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(in)       :: u3
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out)  :: homogeneousforce

  integer :: i, r, j

  homogeneousforce = ZERO

  DOI: do i = 1, 3, 1
    IFAI1: if ( avalexist(i,1) .eqv. .TRUE. ) then
      homogeneousforce(:,i) = homogeneousforce(:,i) + ( &
                              ( amat(i,1) * u1(:) ) )
    end if IFAI1
    IFAI2: if ( avalexist(i,2) .eqv. .TRUE. ) then
      homogeneousforce(:,i) = homogeneousforce(:,i) + ( &
                              ( amat(i,2) * u2(:) ) )
    end if IFAI2
    IFAI3: if ( avalexist(i,3) .eqv. .TRUE. ) then
      homogeneousforce(:,i) = homogeneousforce(:,i) + ( &
                              ( amat(i,3) * u3(:) ) )
    end if IFAI3
  end do DOI

end subroutine f_RK3RRotSource
!-------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3BFRotForcing
!> @param ierr should return 0
subroutine f_RK3BFRotForcing ( wvn, force, ierr )

  use f_force,      only: fluid_f1, fluid_f2, fluid_f3, FORCE_NK
  use g_domain,     only: bmat, NODES_KK, NODES_KJ
  use g_parameters, only: B110, B220, B330, FORCE_KF
  use g_constants,  only: ZERO

  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in)      :: wvn
  PetscErrorCode,intent(inout)                              :: ierr
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out) :: force 

  integer                                            :: i, j, k 
  integer                                            :: ilab, jlab, klab
  integer                                            :: m
  real(kind=C_DOUBLE)                                :: kmag2   
  real(kind=C_DOUBLE),dimension(1:3)                 :: wvnlab
  real(kind=C_DOUBLE),dimension(0:1,1:3)             :: fluid_f
  real(kind=C_DOUBLE),dimension(0:1)                 :: k_dot_b

  !> Get indices from wavenumber. Do not use FFTW indexing.
  i = nint(wvn(1))
  j = nint(wvn(2))
  k = nint(wvn(3))

  !> Zero if not forced
  force = ZERO

  IFB110: if ( modulo(i,B110) == 0 ) then
    IFB220: if ( modulo(j,B220) == 0 ) then
      IFB330: if ( modulo(k,B330) == 0 ) then
        !> Wavenumber in the laboratory system.
        !! Currently only diagonal matrix allowed.
        wvnlab = ZERO

        DOM: do m = 1, 3, 1
          wvnlab(m) = wvnlab(m) + bmat(m,m) * wvn(m)
        end do DOM

        !> compute wave number magnitude in laboratory system
        kmag2 = ( wvnlab(1) * wvnlab(1) ) + &
                ( wvnlab(2) * wvnlab(2) ) + &
                ( wvnlab(3) * wvnlab(3) )

        IFFORCED: if ( sqrt(kmag2) <= FORCE_KF ) then
          IFFNOTZERO: if ( sqrt(kmag2) > ZERO ) then

            ilab = i / B110
            jlab = j / B220
            klab = k / B330

            !> Read the correct forcing from forcing array
            IFJ: if ( jlab <= FORCE_NK ) then
              IFJPOSK: if ( klab <= FORCE_NK ) then
                fluid_f(:,1) = fluid_f1(:,klab,jlab,ilab)
                fluid_f(:,2) = fluid_f2(:,klab,jlab,ilab)
                fluid_f(:,3) = fluid_f3(:,klab,jlab,ilab)
              else
                fluid_f(:,1) = fluid_f1(:,klab-NODES_KK,jlab,ilab)
                fluid_f(:,2) = fluid_f2(:,klab-NODES_KK,jlab,ilab)
                fluid_f(:,3) = fluid_f3(:,klab-NODES_KK,jlab,ilab)
              end if IFJPOSK      
            else
              IFJNEGK: if ( klab <= FORCE_NK ) then
                fluid_f(:,1) = fluid_f1(:,klab,jlab-NODES_KJ,ilab)
                fluid_f(:,2) = fluid_f2(:,klab,jlab-NODES_KJ,ilab)
                fluid_f(:,3) = fluid_f3(:,klab,jlab-NODES_KJ,ilab)
              else
                fluid_f(:,1) = fluid_f1(:,klab-NODES_KK,jlab-NODES_KJ,ilab)
                fluid_f(:,2) = fluid_f2(:,klab-NODES_KK,jlab-NODES_KJ,ilab)
                fluid_f(:,3) = fluid_f3(:,klab-NODES_KK,jlab-NODES_KJ,ilab)
              end if IFJNEGK      
            end if IFJ

            !> Enforce continuity in the laboratory system
            k_dot_b  = ( wvnlab(1) * fluid_f(:,1) + &
                         wvnlab(2) * fluid_f(:,2) + &
                         wvnlab(3) * fluid_f(:,3) )

            force(:,1) = fluid_f(:,1) - ( (wvnlab(1) / kmag2) * k_dot_b )
            force(:,2) = fluid_f(:,2) - ( (wvnlab(2) / kmag2) * k_dot_b )
            force(:,3) = fluid_f(:,3) - ( (wvnlab(3) / kmag2) * k_dot_b )
  
!            write(*,*) wvn, wvnlab, force
               
          end if IFFNOTZERO
        end if IFFORCED

      end if IFB330
    end if IFB220
  end if IFB110

end subroutine f_RK3BFRotForcing 
!------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RRotForcing
!> @param ierr should return 0
subroutine f_RK3RRotForcing ( wvn, force, ierr )

  use f_force,      only: fluid_f1, fluid_f2, fluid_f3, FORCE_NK
  use g_domain,     only: bmat, NODES_KK, NODES_KJ
  use g_parameters, only: B110, B220, B330, FORCE_KF
  use g_constants,  only: ZERO

  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in)      :: wvn
  PetscErrorCode,intent(inout)                              :: ierr
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out) :: force 

  integer                                            :: i, j, k 
  integer                                            :: ilab, jlab, klab
  integer                                            :: m
  real(kind=C_DOUBLE)                                :: kmag2   
  real(kind=C_DOUBLE),dimension(1:3)                 :: wvnlab
  real(kind=C_DOUBLE),dimension(0:1,1:3)             :: fluid_f
  real(kind=C_DOUBLE),dimension(0:1)                 :: k_dot_b

  !> Get indices from wavenumber. Do not use FFTW indexing.
  i = nint(wvn(1))
  j = nint(wvn(2))
  k = nint(wvn(3))

  !> Zero if not forced
  force = ZERO

  IFB110: if ( modulo(i,B110) == 0 ) then
    IFB220: if ( modulo(j,B220) == 0 ) then
      IFB330: if ( modulo(k,B330) == 0 ) then
        !> Wavenumber in the laboratory system.
        !! Currently only diagonal matrix allowed.
        wvnlab = ZERO

        DOM: do m = 1, 3, 1
          wvnlab(m) = wvnlab(m) + bmat(m,m) * wvn(m)
        end do DOM

        !> compute wave number magnitude in laboratory system
        kmag2 = ( wvnlab(1) * wvnlab(1) ) + &
                ( wvnlab(2) * wvnlab(2) ) + &
                ( wvnlab(3) * wvnlab(3) )

        IFFORCED: if ( sqrt(kmag2) <= FORCE_KF ) then
          IFFNOTZERO: if ( sqrt(kmag2) > ZERO ) then

            ilab = i / B110
            jlab = j / B220
            klab = k / B330

            !> Read the correct forcing from forcing array
            IFJ: if ( jlab <= FORCE_NK ) then
              IFJPOSK: if ( klab <= FORCE_NK ) then
                fluid_f(:,1) = fluid_f1(:,klab,jlab,ilab)
                fluid_f(:,2) = fluid_f2(:,klab,jlab,ilab)
                fluid_f(:,3) = fluid_f3(:,klab,jlab,ilab)
              else
                fluid_f(:,1) = fluid_f1(:,klab-NODES_KK,jlab,ilab)
                fluid_f(:,2) = fluid_f2(:,klab-NODES_KK,jlab,ilab)
                fluid_f(:,3) = fluid_f3(:,klab-NODES_KK,jlab,ilab)
              end if IFJPOSK      
            else
              IFJNEGK: if ( klab <= FORCE_NK ) then
                fluid_f(:,1) = fluid_f1(:,klab,jlab-NODES_KJ,ilab)
                fluid_f(:,2) = fluid_f2(:,klab,jlab-NODES_KJ,ilab)
                fluid_f(:,3) = fluid_f3(:,klab,jlab-NODES_KJ,ilab)
              else
                fluid_f(:,1) = fluid_f1(:,klab-NODES_KK,jlab-NODES_KJ,ilab)
                fluid_f(:,2) = fluid_f2(:,klab-NODES_KK,jlab-NODES_KJ,ilab)
                fluid_f(:,3) = fluid_f3(:,klab-NODES_KK,jlab-NODES_KJ,ilab)
              end if IFJNEGK      
            end if IFJ

            !> Enforce continuity in the laboratory system
            k_dot_b  = ( wvnlab(1) * fluid_f(:,1) + &
                         wvnlab(2) * fluid_f(:,2) + &
                         wvnlab(3) * fluid_f(:,3) )

            force(:,1) = fluid_f(:,1) - ( (wvnlab(1) / kmag2) * k_dot_b )
            force(:,2) = fluid_f(:,2) - ( (wvnlab(2) / kmag2) * k_dot_b )
            force(:,3) = fluid_f(:,3) - ( (wvnlab(3) / kmag2) * k_dot_b )
  
!            write(*,*) wvn, wvnlab, force
               
          end if IFFNOTZERO
        end if IFFORCED

      end if IFB330
    end if IFB220
  end if IFB110

end subroutine f_RK3RRotForcing 
!------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RungeKutta
!> Solve the Runge-Kutta step at a grid point.
!> @param ierr should return 0
pure subroutine f_RK3RungeKutta ( rk31, rk32, rk33, u1, u2, u3, force, rkstage, deltat )

  use g_constants, only : ZERO, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                          ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: rk31
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: rk32
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: rk33
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: u1
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: u2
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: u3
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(in) :: force
  integer, intent(in)                               :: rkstage
  real(kind=C_DOUBLE),intent(in)                    :: deltat

  real(kind=C_DOUBLE)                :: rkbeta
  real(kind=C_DOUBLE)                :: rkgamma

  IFSTAGE: if ( rkstage == 1 ) then
    rkbeta = ZERO
    rkgamma = THIRD

  else if ( rkstage == 2 ) then
    rkbeta = - FIVE_OVER_NINE 
    rkgamma = FIFTEEN_OVER_SIXTEEN

  else if ( rkstage == 3 ) then
    rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
    rkgamma = EIGHT_OVER_FIFTEEN
  end if IFSTAGE

  rk31(:) = ( rkbeta * rk31(:) ) + force(:,1)
  rk32(:) = ( rkbeta * rk32(:) ) + force(:,2)
  rk33(:) = ( rkbeta * rk33(:) ) + force(:,3)

  u1(:) = u1(:) + ( rkgamma * deltat * rk31(:) )
  u2(:) = u2(:) + ( rkgamma * deltat * rk32(:) )
  u3(:) = u3(:) + ( rkgamma * deltat * rk33(:) )

end subroutine f_RK3RungeKutta
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3ContinuityFixIsotropic 
!> @param ierr should return 0
pure subroutine f_RK3ContinuityFixIsotropic ( k1, k2, k3, kmag2, u1, u2, u3, ierr )

  use g_constants,  only : ZERO

  implicit none
 
  real(kind=C_DOUBLE),intent(in)                           :: k1, k2, k3
  real(kind=C_DOUBLE),intent(in)                           :: kmag2
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)         :: u1, u2, u3
  PetscErrorCode,intent(inout)                                    :: ierr

  real(kind=C_DOUBLE),dimension(0:1)                       :: kdotu

  kdotu = ZERO

  IFNOTZERO: if ( kmag2 > ZERO ) then 
               
    !> calculate dot product of wave vector and velocity
    kdotu = k1*u1 + k2*u2 + k3*u3
              
    !> enforce continuity in wave space
    u1 = u1 - ( k1 * kdotu / kmag2 )
    u2 = u2 - ( k2 * kdotu / kmag2 )
    u3 = u3 - ( k3 * kdotu / kmag2 )
                   
  end if IFNOTZERO
               
end subroutine f_RK3ContinuityFixIsotropic
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3ContinuityFixBFRot 
!> @param ierr should return 0
pure subroutine f_RK3ContinuityFixBFRot ( k1, k2, k3, u1, u2, u3, ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, binv, bvalexist

  implicit none
 
  real(kind=C_DOUBLE),intent(in)                    :: k1, k2, k3
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: u1, u2, u3
  PetscErrorCode,intent(inout)                             :: ierr

  real(kind=C_DOUBLE),dimension(1:3)                :: klab
  real(kind=C_DOUBLE),dimension(1:3)                :: kmov
  real(kind=C_DOUBLE),dimension(0:1,1:3)            :: ulab
  real(kind=C_DOUBLE),dimension(0:1,1:3)            :: umov
  real(kind=C_DOUBLE),dimension(0:1)                :: kdotu
  real(kind=C_DOUBLE)                               :: kmag2

  integer                                           :: m, n

  klab = ZERO 
  ulab = ZERO  

  umov(:,1) = u1
  umov(:,2) = u2
  umov(:,3) = u3

  kmov(1) = k1
  kmov(2) = k2
  kmov(3) = k3

  DOM1: do m = 1, 3, 1
    DON1: do n = 1, 3, 1
      IFN1: if ( bvalexist(m,n) .eqv. .TRUE. ) then
        klab(n) = klab(n) + ( bmat(m,n) * kmov(m) )
        ulab(:,m) = ulab(:,m) + ( binv(m,n) * umov(:,n) )
      end if IFN1
    end do DON1
  end do DOM1 

!  write(*,*) 'kmov, klab:', kmov, klab, 'umov, ulab:', umov, ulab

  !> compute scalar wave number 
  kmag2 = ( klab(1) * klab(1) ) + &
          ( klab(2) * klab(2) ) + &
          ( klab(3) * klab(3) )

  IFNOTZERO: if ( kmag2 > ZERO ) then 
               
    kdotu = ZERO

    !> calculate dot product of wave vector and velocity
    kdotu = klab(1)*ulab(:,1) + klab(2)*ulab(:,2) + klab(3)*ulab(:,3) 
              
    !> enforce continuity in wave space
    ulab(:,1) = ulab(:,1) - ( klab(1) * kdotu / kmag2 )
    ulab(:,2) = ulab(:,2) - ( klab(2) * kdotu / kmag2 )
    ulab(:,3) = ulab(:,3) - ( klab(3) * kdotu / kmag2 )
    
    umov = ZERO
               
    DOM2: do m = 1, 3, 1
      DON2: do n = 1, 3, 1
        IFN2: if ( bvalexist(m,n) .eqv. .TRUE. ) then
          umov(:,m) = umov(:,m) + ( bmat(m,n) * ulab(:,n) )
        end if IFN2
      end do DON2
    end do DOM2 

    u1 = umov(:,1)
    u2 = umov(:,2)
    u3 = umov(:,3)

  end if IFNOTZERO
               
end subroutine f_RK3ContinuityFixBFRot
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3ContinuityFixRRot 
!> @param ierr should return 0
pure subroutine f_RK3ContinuityFixRRot ( k1, k2, k3, u1, u2, u3, ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, binv, bvalexist

  implicit none
 
  real(kind=C_DOUBLE),intent(in)                    :: k1, k2, k3
  real(kind=C_DOUBLE),dimension(0:1),intent(inout)  :: u1, u2, u3
  PetscErrorCode,intent(inout)                             :: ierr

  real(kind=C_DOUBLE),dimension(1:3)                :: klab
  real(kind=C_DOUBLE),dimension(1:3)                :: kmov
  real(kind=C_DOUBLE),dimension(0:1,1:3)            :: ulab
  real(kind=C_DOUBLE),dimension(0:1,1:3)            :: umov
  real(kind=C_DOUBLE),dimension(0:1)                :: kdotu
  real(kind=C_DOUBLE)                               :: kmag2

  integer                                           :: m, n

  ulab(:,1) = u1
  ulab(:,2) = u2
  ulab(:,3) = u3

  umov = ZERO
               
  DOM2: do m = 1, 3, 1
    DON2: do n = 1, 3, 1
      IFN2: if ( bvalexist(m,n) .eqv. .TRUE. ) then
        umov(:,m) = umov(:,m) + ( bmat(m,n) * ulab(:,n) )
      end if IFN2
    end do DON2
  end do DOM2 

  kmov(1) = k1
  kmov(2) = k2
  kmov(3) = k3

!  write(*,*) 'kmov, klab:', kmov, klab, 'umov, ulab:', umov, ulab

  !> compute scalar wave number 
  kmag2 = ( kmov(1) * kmov(1) ) + &
          ( kmov(2) * kmov(2) ) + &
          ( kmov(3) * kmov(3) )

  IFNOTZERO: if ( kmag2 > ZERO ) then 
               
    kdotu = ZERO

    !> calculate dot product of wave vector and velocity
    kdotu = kmov(1)*umov(:,1) + kmov(2)*umov(:,2) + kmov(3)*umov(:,3) 
              
    !> enforce continuity in wave space
    umov(:,1) = umov(:,1) - ( kmov(1) * kdotu / kmag2 )
    umov(:,2) = umov(:,2) - ( kmov(2) * kdotu / kmag2 )
    umov(:,3) = umov(:,3) - ( kmov(3) * kdotu / kmag2 )
    
    ulab = ZERO  

    DOM1: do m = 1, 3, 1
      DON1: do n = 1, 3, 1
        IFN1: if ( bvalexist(m,n) .eqv. .TRUE. ) then
          ulab(:,m) = ulab(:,m) + ( binv(m,n) * umov(:,n) )
        end if IFN1
      end do DON1
    end do DOM1 

    u1 = ulab(:,1)
    u2 = ulab(:,2)
    u3 = ulab(:,3)

  end if IFNOTZERO
               
end subroutine f_RK3ContinuityFixRRot
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3TransformAllArraysjiRs2Fs
!> @param ierr should return 0
subroutine f_RK3TransformAllArraysjiRs2Fs ( zslice, ierr )

  implicit none
 
  integer,intent(in)                     :: zslice
  PetscErrorCode,intent(inout)                  :: ierr

  call f_RK3TransformCouplingTermjiRs2Fs ( zslice, ierr )
  call f_RK3TransformNonLinearjiRs2Fs ( zslice, ierr )

end subroutine f_RK3TransformAllArraysjiRs2Fs
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3TransformAllArrayskRs2Fs
!> @param ierr should return 0
subroutine f_RK3TransformAllArrayskRs2Fs ( islice, ierr )

  implicit none
 
  integer,intent(in)                     :: islice
  PetscErrorCode,intent(inout)                  :: ierr

  call f_RK3TransformCouplingTermkRs2Fs ( islice, ierr )
  call f_RK3TransformNonLinearkRs2Fs ( islice, ierr )

end subroutine f_RK3TransformAllArrayskRs2Fs
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3Derivative
!> @param u
!> @param k
!> @param dudx
pure subroutine f_RK3Derivative ( u, k, dudx )
  real(kind=C_DOUBLE), dimension(0:1), intent(in)  :: u 
  real(kind=C_DOUBLE), intent(in)                  :: k
  real(kind=C_DOUBLE), dimension(0:1), intent(out) :: dudx

  dudx(1) = k * u(0)      ! i * real part = imaginary 
  dudx(0) = - k * u(1)    ! i * imaginary part = i**2 * real = - real 

end subroutine f_RK3Derivative
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3DerivativeRRot 
!> @param k1
pure subroutine f_RK3DerivativeRRot ( uk, k1, k2, k3, j, dudx )

  use g_constants, only : ZERO
  use g_domain,    only : bmat, binv, bvalexist

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)      :: uk 
  real(kind=C_DOUBLE),intent(in)                     :: k1, k2, k3
  integer, intent(in)                                :: j 
  real(kind=C_DOUBLE),dimension(0:1),intent(out)     :: dudx 
  
  real(kind=C_DOUBLE),dimension(1:3)                 :: wvn 

  integer :: m

  wvn(1) = k1
  wvn(2) = k2
  wvn(3) = k3

  dudx = ZERO

  !> Complex derivative is explicitly written out for 0: real and 
  !! 1: imaginary components. See f_RK3_Derivative for details.
  DOM: do m = 1, 3, 1
    IFBMK: if ( bvalexist(m,j) .eqv. .TRUE. ) then
    end if IFBMK
      dudx(1) = dudx(1) + ( bmat(m,j) * wvn(m) * uk(0) )
      dudx(0) = dudx(0) - ( bmat(m,j) * wvn(m) * uk(1) )
  end do DOM

end subroutine f_RK3DerivativeRRot
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3VorticityBFRot 
!> @param k1
pure subroutine f_RK3VorticityBFRot ( un, k1, k2, k3, n, omega )

  use g_constants, only : ZERO
  use g_domain,    only : bmat, binv, bvalexist

  implicit none

  real(kind=C_DOUBLE),dimension(0:1),intent(in)      :: un 
  real(kind=C_DOUBLE),intent(in)                     :: k1, k2, k3
  integer, intent(in)                                :: n 
  real(kind=C_DOUBLE),dimension(0:1,1:3),intent(out) :: omega 
  
  real(kind=C_DOUBLE),dimension(1:3)                 :: k 

  integer :: i, m, j, p

  k(1) = k1
  k(2) = k2
  k(3) = k3

  omega = ZERO

  !> Complex derivative is explicitly written out for 0: real and 
  !! 1: imaginary components. See f_RK3_Derivative for details.
  DOM: do m = 1, 3, 1
    IFBM1: if ( bvalexist(m,1) .eqv. .TRUE. ) then
      IFBN13: if ( bvalexist(3,n) .eqv. .TRUE. ) then
        omega(1,2) = omega(1,2) - ( bmat(m,1) * binv(3,n) * k(m) * un(0) )
        omega(0,2) = omega(0,2) + ( bmat(m,1) * binv(3,n) * k(m) * un(1) )
      end if IFBN13
      IFBN12: if ( bvalexist(2,n) .eqv. .TRUE. ) then
        omega(1,3) = omega(1,3) + ( bmat(m,1) * binv(2,n) * k(m) * un(0) ) 
        omega(0,3) = omega(0,3) - ( bmat(m,1) * binv(2,n) * k(m) * un(1) ) 
      end if IFBN12
    end if IFBM1
    IFBM2: if ( bvalexist(m,2) .eqv. .TRUE. ) then
      IFBN23: if ( bvalexist(3,n) .eqv. .TRUE. ) then
        omega(1,1) = omega(1,1) + ( bmat(m,2) * binv(3,n) * k(m) * un(0) ) 
        omega(0,1) = omega(0,1) - ( bmat(m,2) * binv(3,n) * k(m) * un(1) ) 
      end if IFBN23
      IFBN21: if ( bvalexist(1,n) .eqv. .TRUE. ) then
        omega(1,3) = omega(1,3) - ( bmat(m,2) * binv(1,n) * k(m) * un(0) )
        omega(0,3) = omega(0,3) + ( bmat(m,2) * binv(1,n) * k(m) * un(1) )
      end if IFBN21
    end if IFBM2
    IFBM3: if ( bvalexist(m,3) .eqv. .TRUE. ) then
      IFBN32: if ( bvalexist(2,n) .eqv. .TRUE. ) then
        omega(1,1) = omega(1,1) - ( bmat(m,3) * binv(2,n) * k(m) * un(0) )
        omega(0,1) = omega(0,1) + ( bmat(m,3) * binv(2,n) * k(m) * un(1) )
      end if IFBN32
      IFBN31: if ( bvalexist(1,n) .eqv. .TRUE. ) then
        omega(1,2) = omega(1,2) + ( bmat(m,3) * binv(1,n) * k(m) * un(0) ) 
        omega(0,2) = omega(0,2) - ( bmat(m,3) * binv(1,n) * k(m) * un(1) ) 
      end if IFBN31
    end if IFBM3
  end do DOM

end subroutine f_RK3VorticityBFRot
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3Tstep
!> Determine the new time step resolution
!> @param rk3_dt \f$ \Delta t \f$ based on the CFL criterion or fixed when the 
!! Rogallo transform is used
!> @param ierr should return 0
subroutine f_RK3Tstep ( dt, ierr )

  use g_parameters, only : MYRANK, MASTER, REPORT, YES, RK3_TIMESTEP, RK3_TIMESTEP_DECIDE, &
                           RK3_TIMESTEP_FIXED, RK3_TIMESTEP_CFL, RK3_CFL, RK3_DT, &
                           F_TYPE, F_ISOTROPIC, F_ISOTROPIC_DEFORMED, F_REMESH_T, F_DEFORM_S, FORCE_TF
  use g_domain,     only : DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME
  use f_fluidstats, only : ucfl
  use g_constants,  only : ZERO, HALF, ONE

  implicit none

  PetscErrorCode,intent(inout)                   :: ierr
  real(kind=C_DOUBLE),intent(inout)       :: dt
  real(kind=C_DOUBLE)                     :: global_u_max

  select case ( RK3_TIMESTEP )
  !-------------------------------
    case ( RK3_TIMESTEP_FIXED ) 

      dt = RK3_DT 

  !-------------------------------
    case ( RK3_TIMESTEP_CFL ) 

      IFCFLMAXZERO: if ( ucfl == ZERO ) then

        dt = 0.02

      else

        !> CFL criterion for isotropic turbulence
        dt = RK3_CFL * DX_MOVINGFRAME / ucfl

      end if IFCFLMAXZERO

  !-------------------------------
    case ( RK3_TIMESTEP_DECIDE ) 
    !> Different procedure to determine step width for different flow types 
    IFFLOWTYPE: if( F_TYPE == F_ISOTROPIC ) then
      IFMAXZERO: if ( ucfl == ZERO ) then
        dt = 0.02
      else
        !> CFL criterion for isotropic turbulence
        dt = RK3_CFL * DX_MOVINGFRAME / ucfl

        IFLARGE: if ( dt > FORCE_TF ) then

          dt = FORCE_TF

        end if IFLARGE

      endif IFMAXZERO

    else if ( F_TYPE == F_ISOTROPIC_DEFORMED ) then

      IFMAXZERODEF: if ( ucfl == ZERO ) then
        dt = 0.02
      else
        !> CFL criterion for isotropic turbulence
        dt = RK3_CFL * min ( DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME ) / ucfl

        IFLARGEDEF: if ( dt > FORCE_TF ) then
  
          dt = FORCE_TF

        end if IFLARGEDEF

      endif IFMAXZERODEF

    else

      !> Fixed time step for Rogallo transform, based on remesh time
      dt = ( ( HALF / F_DEFORM_S ) / real ( F_REMESH_T, kind=C_DOUBLE ) )
!      dt = ONE / ( F_DEFORM_S * real(F_REMESH_T,kind=C_DOUBLE) )
!            should still check for CFL
!            could divide by 2 while CFL > 1
!            in that case would need to adapt number of time steps until remesh

    endif IFFLOWTYPE
 
  end select

  !> report \f$ \Delta t \f$
  IFMASTER: if ( MYRANK == MASTER ) then
    IFREPORT: if ( REPORT == YES ) then
      write(*,20) dt, RK3_CFL * min ( DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME ) / ucfl
      write(*,30)
    end if IFREPORT
  end if IFMASTER

!10  format('    | - setting tstep based on fluid velocity:')
!15  format('    | - setting tstep based on particle velocity:')
20  format('    | -  tstep = ', f22.19, ' ( cfl = ', f8.5, ' )')
30  format('    |')
!40  format('    | - fluid velocity is zero. initialising rk3_dt - ',f12.6)

end subroutine f_RK3Tstep
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3Remesh
!> If flow type is shear flow, see if remesh is required 
!! and if so call remesh subroutine.
!> @param tstep number of the current time step 
!> @param ierr should return 0
subroutine f_RK3Remesh ( tstep, ierr )

  use g_parameters,    only : F_REMESH_T, F_TYPE, F_ISOTROPIC, &
                              RK3_TIMESTEP, RK3_TIMESTEP_DECIDE, F_DEFORM_S, &
                              MYRANK, MASTER, DEALIASING, DEALIASING_SPHERICAL, &
                              DEALIASING_LES, P_TRACK_PART, YES, B110
  use g_constants,     only : ONE
  use f_fluidstats,    only : f_FluidstatsFs
  use f_arrays,        only : f_ArraysDealiasSpherical, f_ArraysLESCutoff
  use p_control,       only : p_ControlRemesh

  implicit none
  
  integer,intent(in)                 :: tstep
  PetscErrorCode,intent(inout)              :: ierr

  real(kind=C_DOUBLE)                :: dt

  integer                            :: a, b

  !> if shear flow, check if remeshing at current time step required
  IFSHEAR: if ( F_TYPE /= F_ISOTROPIC ) then

    !> Only remesh if timestepping method is chosen by PANDORA
    IFTSTEPDECIDE: if ( RK3_TIMESTEP == RK3_TIMESTEP_DECIDE ) then
  
      !> Correction for stretched domain @todo extend to general case
      a = tstep - ( B110 * F_REMESH_T ) - 1
!      a = tstep - F_REMESH_T - 1
      b = 2 * B110 * F_REMESH_T
!      b = 2 * F_REMESH_T
    
      !> call g_RK3Remesh if remesh required at this time step
      IFREMESH: if ( mod ( a, b ) == 0 ) then
!        !> Do wave space fluid analysis before remesh.
!        call f_FluidstatsFs ( tstep, ONE, ierr )

        IFMASTER: if ( MYRANK == MASTER ) then
          write(*,*) 'tstep: ', tstep, ' Remeshing now.'
        end if IFMASTER

        IFSPH1: if ( DEALIASING == DEALIASING_SPHERICAL ) then
          call f_ArraysDealiasSpherical ( ierr )
        end if IFSPH1

        IFLES1: if ( DEALIASING == DEALIASING_LES ) then
        !> Truncate all wavemodes outside cutoff
          call f_ArraysLESCutoff ( ierr )
        end if IFLES1

        !> If particles remesh particles
        IFPART: if ( P_TRACK_PART == YES ) then
          call PetscPrintf( PETSC_COMM_WORLD, 'Remeshing particles. \n', ierr )
          call p_ControlRemesh ( ierr )
        end if IFPART

        !> Remesh fluid grid
        call f_RK3RemeshWaveSpace ( ierr )

        IFSPH2: if ( DEALIASING == DEALIASING_SPHERICAL ) then
          call f_ArraysDealiasSpherical ( ierr )
        end if IFSPH2

        IFLES2: if ( DEALIASING == DEALIASING_LES ) then
          !> Truncate all wavemodes outside cutoff
          call f_ArraysLESCutoff ( ierr )
        end if IFLES2

      end if IFREMESH



!      !> At the beginning we do not need a remesh,
!      !! but we do need to calculate the phase shift.
!      IFINIT: if ( tstep > 1 ) then

!        !> Call subroutine to shift velocity planes to the laboratory system.
!        call f_RK3RemeshPhaseShift ( ierr )

!      else

!        !> dt is not yet calculated, so we need to do it here
!        dt = ONE / ( F_DEFORM_S * real(F_REMESH_T,kind=C_DOUBLE) )
!        call f_RK3RemeshPhaseShiftCalculate ( dt, ierr ) 

!      end if IFINIT 

    end if IFTSTEPDECIDE

  end if IFSHEAR
  
end subroutine f_RK3Remesh
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshWaveSpace
!> Remesh the flow field.
!> @param ierr should return 0
subroutine f_RK3RemeshWaveSpace ( ierr )
  
  use g_constants,  only : ZERO, HALF
  use g_domain,     only : KK_MIN, KK_MAX, bmat, binv
  use g_parameters, only : B110, NS_DEFORMED, NS_BF_ROT, NS_R_CONV, NS_R_ROT
  use f_arrays,     only : da3w, u1_3w, u2_3w, u3_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           arr_u1_remesh, arr_u2_remesh, arr_u3_remesh, &
                           i_min_3w, j_min_3w, k_min_3w, &
                           i_max_3w, j_max_3w, k_max_3w

  implicit none
   
  PetscErrorCode,intent(inout) :: ierr
  
  integer               :: cutoff, actualk, knew, nodes
  integer               :: i, j, k

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

  nodes = 2 * KK_MAX + 1

  cutoff = KK_MAX

!> @todo we are assuming equal number of wave modes in each direction here. might need modifying for longer domains

  !> B**(-1) for conversion to laboratory system before remesh (St=0.5).
  binv(1,3) = ( real(B110,kind=C_DOUBLE) * HALF ) / bmat(3,3)
!  binv(1,3) =  HALF  / bmat(3,3)
  !> B for conversion to moving system after remesh (St=-0.5).
  bmat(1,3) = - bmat(1,1) * ( - ( real(B110,kind=C_DOUBLE) * HALF ) )
!  bmat(1,3) = - bmat(1,1) * ( - HALF ) 

  DOI: do i = i_min_3w, i_max_3w, 1
    DOJ: do j = j_min_3w, j_max_3w, 1

      !> initialise remesh buffer 
      arr_u1_remesh = ZERO
      arr_u2_remesh = ZERO
      arr_u3_remesh = ZERO

      DOK: do k = k_min_3w, k_max_3w, 1

        IFOLDNEGATIVE: if ( k > KK_MAX ) then
          actualk = k - nodes
        else
          actualk = k
        end if IFOLDNEGATIVE

        ! i has no negative wavenumbers, so this is always true:
        knew = actualk - i
!        knew = actualk + i

        !> Only remesh the wave numbers that are not aliased
        IFNOTALIASED: if ( knew >= - KK_MAX ) then
!        IFNOTALIASED: if ( knew <= KK_MAX ) then
!        IFNOTALIASED: if ( knew <= (HALF*KK_MAX) ) then

          !> For negative wave numbers we need to go back to FFTW indexing
          IFNEWNEGATIVE: if ( knew < 0 ) then
            knew = knew + nodes
          end if IFNEWNEGATIVE

!!          write(*,*) 'i, j: ', i, j, ', k old: ', actualk, ', k new: ', knew

          arr_u1_remesh(:,knew) = arr_u1_3w(:,k,j,i)
          arr_u2_remesh(:,knew) = arr_u2_3w(:,k,j,i)
          arr_u3_remesh(:,knew) = arr_u3_3w(:,k,j,i)

        end if IFNOTALIASED

      enddo DOK

      !> Correct velocity from moving system before remesh to moving
      !! system after remesh. 

      !> Only u1 component changes.

!      write(*,*) ' i, j: ', i, j, &
!                 ', u1 before: ', arr_u1_3w(:,:,j,i), ', u1 after: ', &
!                 ( arr_u1_remesh(:,:) )!+ &
                ! ( binv(1,3) * arr_u3_remesh(:,:) ) ) + &
                ! ( bmat(1,3) * arr_u3_remesh(:,:) )
!      arr_u1_3w(:,:,j,i) = ( arr_u1_remesh(:,:) + &
!                           ( binv(1,3) * arr_u3_remesh(:,:) ) ) + &
!                           ( bmat(1,3) * arr_u3_remesh(:,:) ) 

!      arr_u1_3w(:,:,j,i) = ( arr_u1_remesh(:,:) - &
!                           ( binv(1,3) * arr_u3_remesh(:,:) ) ) - &
!                           ( bmat(1,3) * arr_u3_remesh(:,:) ) 
      select case ( NS_DEFORMED )
      !-------------------------------
        case ( NS_BF_ROT ) 
          arr_u1_3w(:,:,j,i) = arr_u1_remesh(:,:) + arr_u3_remesh(:,:)
          arr_u2_3w(:,:,j,i) = arr_u2_remesh(:,:)
          arr_u3_3w(:,:,j,i) = arr_u3_remesh(:,:)
        case ( NS_R_ROT ) 
          arr_u1_3w(:,:,j,i) = arr_u1_remesh(:,:) 
          arr_u2_3w(:,:,j,i) = arr_u2_remesh(:,:)
          arr_u3_3w(:,:,j,i) = arr_u3_remesh(:,:)
        case ( NS_R_CONV ) 
          arr_u1_3w(:,:,j,i) = arr_u1_remesh(:,:) 
          arr_u2_3w(:,:,j,i) = arr_u2_remesh(:,:)
          arr_u3_3w(:,:,j,i) = arr_u3_remesh(:,:)
      end select

    enddo DOJ
  enddo DOI

  !> Need to update binv for new mesh
  binv(1,3) = ( - ( real(B110,kind=C_DOUBLE) * HALF ) ) / bmat(3,3)

!  write(*,*) 'bmat(1,3): ', bmat(1,3), ', binv(1,3): ', binv(1,3)

!> Return velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

end subroutine f_RK3RemeshWaveSpace
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_ControlCalcRemesh
!> If flow type is shear flow, see if remesh is required 
!! and if so call remesh subroutine.
!> @param tstep number of the current time step 
!> @param ierr should return 0
!subroutine g_ControlCalcRemesh ( tstep, ierr )

!  use g_parameters,    only : F_REMESH_T, F_TYPE, F_ISOTROPIC, &
!                              RK3_TIMESTEP, RK3_TIMESTEP_DECIDE
!  use g_rk3,           only : g_RK3Remesh
!  use f_rk3,           only : f_RK3RemeshFluidRealSpace

!  implicit none
  
!  integer,intent(in)                 :: tstep
!  integer,intent(inout)              :: ierr

!  integer                            :: a, b

  !> if shear flow, check if remeshing at current time step required
!  IFSHEAR: if ( F_TYPE /= F_ISOTROPIC ) then

    !> Only remesh if timestepping method is chosen by PANDORA
!    IFTSTEPDECIDE: if ( RK3_TIMESTEP == RK3_TIMESTEP_DECIDE ) then
  
!      a = tstep - F_REMESH_T - 1
!      b = 2 * F_REMESH_T
    
      !> call g_RK3Remesh if remesh required at this time step
!      IFREMESH: if ( mod ( a, b ) == 0 ) then
!        call g_RK3Remesh ( ierr )
!        call f_RK3RemeshFluidRealSpace ( ierr )
!      end if IFREMESH

!    end if IFTSTEPDECIDE

!  end if IFSHEAR
  
!end subroutine g_ControlCalcRemesh
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshFluidRealSpace
!> @param ierr should return 0
subroutine f_RK3RemeshFluidRealSpace ( ierr )

  use f_arrays,     only: da3w, da1r2w, da2w, da2r, i_min_3w, i_max_3w, &
                          z_min, z_max, x_min_2r, x_max_2r, u1_3w, u2_3w, u3_3w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          u1_1r2w, u2_1r2w, u3_1r2w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          u4_1r2w, u5_1r2w, u6_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          u4_2w, u5_2w, u6_2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          u1_2r, u2_2r, u3_2r, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          u4_2r, u5_2r, u6_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r

  use g_parameters, only: NODES_X, NODES_Z
  use g_constants,  only: HALF
  use g_domain,     only: bmat, binv

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: islice                 
  integer               :: zslice                 

  integer               :: xold, xnew

  PetscErrorCode        :: perr

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

  call DMDAVecGetArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )

  call DMDAVecGetArrayF90 ( da1r2w, u4_1r2w, arr_u4_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u5_1r2w, arr_u5_1r2w, perr )
  call DMDAVecGetArrayF90 ( da1r2w, u6_1r2w, arr_u6_1r2w, perr )

  call DMDAVecGetArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )

  call DMDAVecGetArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecGetArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

  call DMDAVecGetArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecGetArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

  !Loop through slices in i direction.
  DOISLICES: do islice = i_min_3w, i_max_3w, 1

  !> Transform the velocity components in k direction and compute vorticity while at it
    call f_RK3RemeshPrepareVelocity ( islice, ierr )

  end do DOISLICES


  !> Loop through slices in z direction.
  DOZSLICES: do zslice = z_min, z_max, 1
               
    !> Transform velocity to real space.
    call f_RK3RemeshTransformToRealSpace ( zslice, ierr )

    !> @todo do the remesh here
    !! using 4, 5, 6 components as buffer for the remesh

!    write(*,*) 'z, zmin, zmax: ', zslice, z_min, z_max, ' NODES_Z: ', NODES_Z, &
!               'xmin, xmax: ', x_min_2r, x_max_2r

    DOX: do xold = x_min_2r, x_max_2r, 1
      
      xnew = xold - zslice !+ ( NODES_Z / 2 )

      IFSHIFT: if ( xnew < 0 ) then

        xnew = xnew + NODES_X

      else if ( xnew >= NODES_X ) then
       
        xnew = xnew - NODES_X

      end if IFSHIFT

!      write(*,*) 'x old, x new: ', xold, xnew, ' zslice: ', zslice

      arr_u4_2r(xnew,:) = arr_u1_2r(xold,:)
      arr_u5_2r(xnew,:) = arr_u2_2r(xold,:)
      arr_u6_2r(xnew,:) = arr_u3_2r(xold,:)

    end do DOX

    !> Transform velocity to wave space.
    call f_RK3TransformNonLinearjiRs2Fs ( zslice, ierr )
               
  end do DOZSLICES

  !> Transformation matrices
  binv(1,3) =               HALF / bmat(3,3)
  bmat(1,3) = - bmat(1,1) * ( - HALF )

  !> Loop through slices in i direction
  DOISLICESBACK: do islice = i_min_3w, i_max_3w, 1

    call f_RK3TransformNonLinearkRs2Fs ( islice, ierr )

!      write(*,*) ' i, j: ', i, j, &
!                 ', u1 before: ', arr_u1_3w(:,:,j,i), ', u1 after: ', &
!                 ( arr_u1_remesh(:,:) )!+ &
                ! ( binv(1,3) * arr_u3_remesh(:,:) ) ) + &
                ! ( bmat(1,3) * arr_u3_remesh(:,:) )

    arr_u1_3w(:,:,:,islice) = ( arr_u4_2w(:,:,:) + &
                              ( binv(1,3) * arr_u6_2w(:,:,:) ) ) + &
                              ( bmat(1,3) * arr_u6_2w(:,:,:) ) 

    arr_u2_3w(:,:,:,islice) = arr_u5_2w(:,:,:)
    arr_u3_3w(:,:,:,islice) = arr_u6_2w(:,:,:)
               
  end do DOISLICESBACK

!> Return velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da2w, u4_2w, arr_u4_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u5_2w, arr_u5_2w, perr )
  call DMDAVecRestoreArrayF90 ( da2w, u6_2w, arr_u6_2w, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u1_2r, arr_u1_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u2_2r, arr_u2_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u3_2r, arr_u3_2r, perr )

  call DMDAVecRestoreArrayF90 ( da2r, u4_2r, arr_u4_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u5_2r, arr_u5_2r, perr )
  call DMDAVecRestoreArrayF90 ( da2r, u6_2r, arr_u6_2r, perr )

  call DMDAVecRestoreArrayF90 ( da1r2w, u1_1r2w, arr_u1_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u2_1r2w, arr_u2_1r2w, perr )
  call DMDAVecRestoreArrayF90 ( da1r2w, u3_1r2w, arr_u3_1r2w, perr )

  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

end subroutine f_RK3RemeshFluidRealSpace
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshPrepareVelocity
!> @param ierr should return 0
subroutine f_RK3RemeshPrepareVelocity ( islice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, bmat
  use g_parameters, only: NODES_Z, F_TYPE, F_ISOTROPIC
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w, k_max_3w, z_max, j_max_1r2w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_2w, arr_u5_2w, arr_u6_2w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w, &
                          f_ArraysWavenumber, f_ArraysWavenumberPositive, f_ArraysWavenumberNegative
  use f_fftw,       only: f_Fftw_k_Fs2Rs, arr_k_in_da, arr_kjjk_transposeout

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: islice

  integer                                :: slicelocal
  integer                                :: j, jj, k, kk, proc_j, jtranspose, ktranspose, j1r2w, component
  real(kind=C_DOUBLE)                    :: wavenumber_k, wavenumber_j, wavenumber_i
  real(kind=C_DOUBLE),dimension(0:1)     :: velocity, dudy, dudz, dvdx, dvdz, dwdx, dwdy
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: vorticity

  type(C_PTR)                            :: copytoarray

!> wave number to be computed from array indices

! Copy each velocity component to the zero-padded velocity buffer array.
  !> Numbering of buffer array is local, numbering of velocity array is global
  !> PETSc uses C indexing, FFTW Fortran pointers use Fortran indexing!! 

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = islice - i_min_3w + 1

  !> u1 component
  jj = 1
  DOU1J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOU1KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u1_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU1KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU1KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU1KZEROS

    !> Copy negative wavenumbers
    DOU1KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u1_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU1KNEGATIVE

    jj = jj + 1
  end do DOU1J3W

  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  copytoarray = c_loc(arr_u1_1r2w)
  call f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )

  !> u2 component
  jj = 1
  DOU2J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOU2KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u2_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU2KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU2KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU2KZEROS

    !> Copy negative wavenumbers
    DOU2KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u2_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU2KNEGATIVE

  jj = jj + 1
  end do DOU2J3W

  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  copytoarray = c_loc(arr_u2_1r2w)
  call f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )


  !> u3 component
  jj = 1
  DOU3J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOU3KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u3_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU3KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU3KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU3KZEROS

    !> Copy negative wavenumbers
    DOU3KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u3_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU3KNEGATIVE

    jj = jj + 1
  end do DOU3J3W

  !> Do FFT in k direction and transpose from kj to jk
  !! Copy the result of the FFT to the buffer array
  copytoarray = c_loc(arr_u3_1r2w)
  call f_RK3TransformVelocitykFs2Rs ( copytoarray, slicelocal, ierr )

end subroutine f_RK3RemeshPrepareVelocity
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshTransformToRealSpace
!> @param ierr should return 0
subroutine f_RK3RemeshTransformToRealSpace ( zslice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KJ_MAX
  use g_parameters, only: NODES_Y, NODES_X
  use f_arrays,     only: z_min, i_min_1r2w, j_min_1r2w, y_min_2r, &
                          i_max_1r2w, j_max_1r2w, y_max_2r, &
                          arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                          arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                          arr_u1_1r2w, arr_u2_1r2w, arr_u3_1r2w, &
                          arr_u4_1r2w, arr_u5_1r2w, arr_u6_1r2w, &
                          nprocsi_1r2w, i_allwidths_1r2w, i_maxwidth_1r2w
  use f_fftw,       only: f_Fftw_ji_Fs2Rs, arr_j_in_da, arr_i_real

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: zslice

  integer                                :: slicelocal
  integer                                :: i, ii, j, jj, proc_i, itranspose, jtranspose, i2r, component

  type(C_PTR)                            :: inputarray, copytoarray

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = zslice - z_min + 1

  !> u1 component
  inputarray = c_loc(arr_u1_1r2w)
  copytoarray = c_loc(arr_u1_2r)
  call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  !> u2 component
  inputarray = c_loc(arr_u2_1r2w)
  copytoarray = c_loc(arr_u2_2r)
  call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

  !> u3 component
  inputarray = c_loc(arr_u3_1r2w)
  copytoarray = c_loc(arr_u3_2r)
  call f_RK3TransformjiFs2Rs ( inputarray, copytoarray, slicelocal, ierr )

end subroutine f_RK3RemeshTransformToRealSpace
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshPhaseShift
!> @param ierr should return 0
subroutine f_RK3RemeshPhaseShift ( ierr )

  use f_arrays,     only: da3w, i_min_3w, i_max_3w, &
                          u1_3w, u2_3w, u3_3w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w

  use g_parameters, only: NODES_X, NODES_Z, MYRANK
  use g_constants,  only: ZERO, HALF
  use g_domain,     only: bmat, binv

  implicit none
 
  PetscErrorCode,intent(inout) :: ierr
  integer               :: islice                 

  PetscErrorCode        :: perr

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

  !> Transform u1 into the laboratory system
!  arr_u1_3w(:,:,:,:) = ( binv(1,1) * arr_u1_3w(:,:,:,:) ) &
!                     + ( binv(1,3) * arr_u3_3w(:,:,:,:) )

!  write(*,*) MYRANK, ' bmat before: ', bmat, ' binv before: ', binv, &
!             'u1 before: ', arr_u1_3w(:,:,:,:), &
!             'u2 before: ', arr_u2_3w(:,:,:,:), &
!             'u3 before: ', arr_u3_3w(:,:,:,:)
             

  !Loop through slices in i direction.
  DOISLICES: do islice = i_min_3w, i_max_3w, 1

  !> Transform the velocity components in k direction and compute vorticity while at it
    call f_RK3RemeshPhaseShiftSlice ( islice, ierr )

  end do DOISLICES

!  !> Transform u1 into the moving system
!  arr_u1_3w(:,:,:,:) = bmat(1,1) * arr_u1_3w(:,:,:,:)

  !> Set the shear components of bmat to zero.
  bmat(1,3) = ZERO
  binv(1,3) = ZERO

!  write(*,*) MYRANK, ' bmat after: ', bmat, ' binv after: ', binv, &
!             'u1 after: ', arr_u1_3w(:,:,:,:), &
!             'u2 after: ', arr_u2_3w(:,:,:,:), &
!             'u3 after: ', arr_u3_3w(:,:,:,:)
             
  !> Return velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

end subroutine f_RK3RemeshPhaseShift
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshPhaseShiftSlice
!> @param ierr should return 0
subroutine f_RK3RemeshPhaseShiftSlice ( islice, ierr )

  use g_constants,  only: ZERO
  use g_domain,     only: KI_MAX, KJ_MAX, KK_MAX, K0I, K0J, K0K, bmat
  use g_parameters, only: NODES_Z, F_TYPE, F_ISOTROPIC
  use f_arrays,     only: i_min_3w, j_min_3w, k_min_3w, z_min, j_min_1r2w, &
                          j_max_3w,  k_max_3w, z_max, j_max_1r2w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          nprocsj_3w, j_allwidths_3w, j_maxwidth_3w, &
                          j_width_3w 
  use f_fftw,       only: f_Fftw_FftOnly_k_Fs2Rs, f_Fftw_FftOnly_k_Rs2Fs, &
                          arr_k_in_da, arr_k_out_da 

  implicit none
 
  PetscErrorCode,intent(inout)                  :: ierr
  integer,intent(in)                     :: islice

  integer                                :: slicelocal
  integer                                :: j, jj, k, kk, proc_j, jtranspose, ktranspose, j1r2w, component
  real(kind=C_DOUBLE),dimension(0:1)     :: velocity


!> wave number to be computed from array indices

! Copy each velocity component to the zero-padded velocity buffer array.
  !> Numbering of buffer array is local, numbering of velocity array is global
  !> PETSc uses C indexing, FFTW Fortran pointers use Fortran indexing!! 

  !> We will work in local coordinates when communicating with FFTW.
  !! Therefore we need the local coordinate of the slice we are working on.
  slicelocal = islice - i_min_3w + 1

  !> u1 component
  jj = 1
  DOU1J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOU1KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u1_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU1KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU1KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU1KZEROS

    !> Copy negative wavenumbers
    DOU1KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u1_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU1KNEGATIVE

    jj = jj + 1
  end do DOU1J3W

  !> Do FFT in k direction
  call f_Fftw_FftOnly_k_Fs2Rs ( ierr )

  call f_RK3RemeshPhaseShiftComponent ( islice, ierr )

  !> Do backwards FFT in k direction
  call f_Fftw_FftOnly_k_Rs2Fs ( ierr )

  !> Copy back to original array
  jj = 1
  DOJ3W1: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOKPOSITIVE1: do k = 0, KK_MAX, 1
      !> Copy to buffer
      arr_u1_3w(:,k,j,islice) = arr_k_out_da(:,kk,jj)
      kk = kk + 1
    end do DOKPOSITIVE1

    !> Get rid of the highest 1/3 wavenumbers, both positive and negative
     kk = NODES_Z - KK_MAX + 1

    !> Copy negative wavenumbers
    DOKNEGATIVE1: do k = KK_MAX + 1, k_max_3w
      !> Copy to buffer
      arr_u1_3w(:,k,j,islice) = arr_k_out_da(:,kk,jj)
      kk = kk + 1
    end do DOKNEGATIVE1

    jj = jj + 1
  end do DOJ3W1




  !> u2 component
  jj = 1
  DOU2J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOU2KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u2_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU2KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU2KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU2KZEROS

    !> Copy negative wavenumbers
    DOU2KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u2_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU2KNEGATIVE

    jj = jj + 1
  end do DOU2J3W

  !> Do FFT in k direction
  call f_Fftw_FftOnly_k_Fs2Rs ( ierr )

  call f_RK3RemeshPhaseShiftComponent ( islice, ierr )

  !> Do backwards FFT in k direction
  call f_Fftw_FftOnly_k_Rs2Fs ( ierr )

  !> Copy back to original array
  jj = 1
  DOJ3W2: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOKPOSITIVE2: do k = 0, KK_MAX, 1
      !> Copy to buffer
      arr_u2_3w(:,k,j,islice) = arr_k_out_da(:,kk,jj)
      kk = kk + 1
    end do DOKPOSITIVE2

    !> Get rid of the highest 1/3 wavenumbers, both positive and negative
     kk = NODES_Z - KK_MAX + 1

    !> Copy negative wavenumbers
    DOKNEGATIVE2: do k = KK_MAX + 1, k_max_3w
      !> Copy to buffer
      arr_u2_3w(:,k,j,islice) = arr_k_out_da(:,kk,jj)
      kk = kk + 1
    end do DOKNEGATIVE2

    jj = jj + 1
  end do DOJ3W2


  !> u3 component
  jj = 1
  DOU3J3W: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOU3KPOSITIVE: do k = k_min_3w, KK_MAX, 1

      !> Copy to buffer
      velocity = arr_u3_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU3KPOSITIVE

    !> Zero padding for the highest 1/3 wavenumbers, both positive and negative
    DOU3KZEROS: do k = KK_MAX + 1, NODES_Z - KK_MAX - 1

      arr_k_in_da(:,kk,jj) = ZERO 
      kk = kk + 1

    end do DOU3KZEROS

    !> Copy negative wavenumbers
    DOU3KNEGATIVE: do k = KK_MAX + 1, k_max_3w

      !> Copy to buffer
      velocity = arr_u3_3w(:,k,j,islice)
      arr_k_in_da(:,kk,jj) = velocity

      kk = kk + 1
    end do DOU3KNEGATIVE

    jj = jj + 1
  end do DOU3J3W

  !> Do FFT in k direction
  call f_Fftw_FftOnly_k_Fs2Rs ( ierr )

  call f_RK3RemeshPhaseShiftComponent ( islice, ierr )

  !> Do backwards FFT in k direction
  call f_Fftw_FftOnly_k_Rs2Fs ( ierr )

  !> Copy back to original array
  jj = 1
  DOJ3W3: do j = j_min_3w, j_max_3w, 1
    kk = 1

    !> Copy positive wavenumbers
    DOKPOSITIVE3: do k = 0, KK_MAX, 1
      !> Copy to buffer
      arr_u3_3w(:,k,j,islice) = arr_k_out_da(:,kk,jj)
      kk = kk + 1
    end do DOKPOSITIVE3

    !> Get rid of the highest 1/3 wavenumbers, both positive and negative
     kk = NODES_Z - KK_MAX + 1

    !> Copy negative wavenumbers
    DOKNEGATIVE3: do k = KK_MAX + 1, k_max_3w
      !> Copy to buffer
      arr_u3_3w(:,k,j,islice) = arr_k_out_da(:,kk,jj)
      kk = kk + 1
    end do DOKNEGATIVE3

    jj = jj + 1
  end do DOJ3W3

end subroutine f_RK3RemeshPhaseShiftSlice
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshPhaseShiftComponent
!> @param ierr should return 0
subroutine f_RK3RemeshPhaseShiftComponent ( i, ierr )

  use f_fftw,       only : arr_k_in_da, arr_k_out_da 
  use f_arrays,     only : arr_remesh, j_width_3w
  use g_parameters, only : NODES_Z

  implicit none
 
  integer,intent(in)                      :: i
  PetscErrorCode,intent(inout)                   :: ierr

  integer                                 :: j, z

  real(kind=C_LONG_DOUBLE),dimension(1:2) :: ushift
  real(kind=C_DOUBLE),dimension(1:2)      :: uout
  real(kind=C_DOUBLE),dimension(1:2)      :: uin

  arr_k_in_da = arr_k_out_da

  DOJ: do j = 1, j_width_3w, 1
    DOZ: do z = 1, NODES_Z, 1

      ushift = arr_remesh(:,z,i)
      uout = arr_k_out_da(:,z,j)

      call f_RK3RemeshPhaseShiftMultiply ( ushift, uout, uin, ierr )

      arr_k_in_da(:,z,j) = uin

!      write(*,*) 'before: ', uout, ' after: ', uin, ' phase shift: ', ushift, &
!                 ' i, z: ', i, z

    end do DOZ
  end do DOJ

end subroutine f_RK3RemeshPhaseShiftComponent
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshPhaseShiftMultiply
!> @param ierr should return 0
pure subroutine f_RK3RemeshPhaseShiftMultiply ( phaseshift, u_before, u_after, ierr )

  implicit none
 
  real(kind=C_LONG_DOUBLE),dimension(1:2),intent(in) :: phaseshift
  real(kind=C_DOUBLE),dimension(1:2),intent(in)      :: u_before
  real(kind=C_DOUBLE),dimension(1:2),intent(out)     :: u_after

  real(kind=C_LONG_DOUBLE),dimension(1:2)            :: u_before_long
  real(kind=C_LONG_DOUBLE),dimension(1:2)            :: u_after_long

  PetscErrorCode,intent(inout)                              :: ierr

!  phaseshift(1) = 1.0
!  phaseshift(2) = 0.0

  u_before_long(1) = real(u_before(1),kind=C_LONG_DOUBLE)
  u_before_long(2) = real(u_before(2),kind=C_LONG_DOUBLE)

  u_after_long(1) = ( phaseshift(1) * u_before_long(1) ) &
                  - ( phaseshift(2) * u_before_long(2) )

  u_after_long(2) = ( phaseshift(1) * u_before_long(2) ) &
                  + ( phaseshift(2) * u_before_long(1) )

  u_after(1) = real(u_after_long(1),kind=C_DOUBLE)
  u_after(2) = real(u_after_long(2),kind=C_DOUBLE)

end subroutine f_RK3RemeshPhaseShiftMultiply
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_RK3RemeshPhaseShiftCalculate
!> @param ierr should return 0
subroutine f_RK3RemeshPhaseShiftCalculate ( dt, ierr )

  use f_arrays,     only : arr_remesh, i_min_3w, i_max_3w, &
                           f_ArraysWavenumberPositiveLong
  use g_domain,     only : K0I
  use g_parameters, only : NODES_Z, F_DEFORM_S
  use g_constants,  only : LONGTWOPI

  implicit none
 
  real(kind=C_DOUBLE),intent(in)         :: dt
  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_LONG_DOUBLE)               :: wvn
  real(kind=C_LONG_DOUBLE)               :: zshift
  real(kind=C_LONG_DOUBLE)               :: shift
  real(kind=C_LONG_DOUBLE)               :: dtlong
  real(kind=C_LONG_DOUBLE)               :: slong

  integer                                :: i, z

  DOI: do i = i_min_3w, i_max_3w, 1
    DOZ: do z = 1, NODES_Z, 1

      wvn = f_ArraysWavenumberPositiveLong ( i, K0I )
      dtlong = real(dt,kind=C_LONG_DOUBLE)
      slong = real(F_DEFORM_S,kind=C_LONG_DOUBLE)
      zshift = real( (- z + ( NODES_Z / 2 )), kind=C_LONG_DOUBLE ) &
!      zshift = real( (+ z - ( NODES_Z / 2 )), kind=C_LONG_DOUBLE ) &
             / real( NODES_Z, kind=C_LONG_DOUBLE )
      zshift = zshift * LONGTWOPI
      shift = slong * dtlong * wvn * zshift
      arr_remesh(1,z,i) = cos ( shift )
      arr_remesh(2,z,i) = sin ( shift )

!      write(*,*) 'argument: ', shift, ' phase shift = ', &
!                 cos ( shift ), sin (shift ), &
!                 ' dt = ', dt, ' s = ', F_DEFORM_S, &
!                 ' z = ', z, ' dz = ', zshift, ' wvn = ', wvn
    end do DOZ
  end do DOI

!  write(*,*) ' Phase shift array: ', arr_remesh(:,:,:)   
!  allocate(arr_remesh(1:2,1:NODES_Z,i_min_3w:i_max_3w),STAT=ierr)

end subroutine f_RK3RemeshPhaseShiftCalculate
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
end module f_rk3
