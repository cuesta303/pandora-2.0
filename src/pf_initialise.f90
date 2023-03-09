!---------------------------------------------------------------------------------
! module f_initialise
!> Initialise the velocity field either from a restart file or according to a 
!! predefined spectrum.
module f_initialise

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

public :: f_InitialiseInit
public :: f_InitialiseInitTaylorGreen

contains
!---------------------------------------------------------------------------------
! subroutine f_InitialiseInit
!> Call subroutine to initialise the velocity either from a restart file or as
!! a solenoidal velocity field.
!> @param ierr should return 0
subroutine f_InitialiseInit ( ierr )

    use g_constants,  only : ZERO
    use g_parameters, only : F_INIT_TYPE, F_INIT_RESTART, F_INIT_TAYLOR_GREEN, READ_FILE, &
                             F_RESTART_SCALE, B110, B220, B330, DEALIASING, &
                             DEALIASING_SPHERICAL, DEALIASING_LES
    use g_restart,    only : g_RestartPetscFluid 
    use f_fluidstats, only : f_FluidstatsFs, f_FluidstatsInitialSpectrum
    use f_arrays,     only : f_ArraysDealiasSpherical, f_ArraysLESCutoff 
    implicit none

    PetscErrorCode,intent(inout)   :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising f_initialise module \n', ierr )

  IFRESTART: if ( F_INIT_TYPE == F_INIT_RESTART ) then
    !> Initialise fluid arrays from a restart file
    call PetscPrintf( PETSC_COMM_WORLD, '    | Initialising fluid arrays from restart files \n', ierr )
    call g_RestartPetscFluid ( READ_FILE, ierr )

    IFSCALED: if ( F_RESTART_SCALE > 1 ) then
      call PetscPrintf( PETSC_COMM_WORLD, '    | Copying fluid arrays from restart arrays \n', ierr )
      call f_InitialiseCopyRestart ( ierr )
    end if IFSCALED

  else if ( F_INIT_TYPE == F_INIT_TAYLOR_GREEN ) then
    !> Initialise fluid arrays with a Taylor-Green vortex
    call f_InitialiseInitTaylorGreen ( ierr )
  else
    !> Initialise fluid arrays as solenoidal
    call PetscPrintf( PETSC_COMM_WORLD, '    | Initialising fluid arrays as solenoidal \n', ierr )

    IFDEFORMED: if ( ( B110 + B220 + B330 ) > 3 ) then
      call f_InitialiseFluidDeformed(ierr)
    else
      call f_InitialiseFluid(ierr)
    end if IFDEFORMED

  end if IFRESTART
    
  IFSPHER: if ( DEALIASING == DEALIASING_SPHERICAL ) then
  !> Truncate all wavemodes outside a 2/3 sphere
    call f_ArraysDealiasSpherical ( ierr )
  end if IFSPHER

  IFLES: if ( DEALIASING == DEALIASING_LES ) then
  !> Truncate all wavemodes outside cutoff
    call f_ArraysLESCutoff ( ierr )
  end if IFLES

    !> Do wave space fluid analysis
!    call f_FluidstatsFs ( 0, ZERO, ierr )
!    call f_FluidstatsInitialSpectrum ( ierr )

end subroutine f_InitialiseInit
!---------------------------------------------------------------------------------
   

!---------------------------------------------------------------------------------
! subroutine f_InitialiseFluid
!> Initialise fluid velocity as a solenoidal field with a given energy spectrum.
!> @param ierr should return 0
subroutine f_InitialiseFluid ( ierr )

  use g_parameters, only : MYRANK, MASTER, NODES_X, NODES_Y, NODES_Z, &
                           RANDOM_SEED_SIZE, RANDOM_SEED_FLUID, &
                           B110, B220, B330 
  use g_random,     only : g_RandomInitFluid
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K
  use g_constants,  only : TWOPI, FOURPI, ZERO
  use f_arrays,     only : f_ArraysWavenumber, i_min_3w, j_min_3w, k_min_3w, &
                           i_max_3w, j_max_3w, k_max_3w, da3w, u1_3w, u2_3w, &
                           u3_3w, arr_u1_3w, arr_u2_3w, arr_u3_3w

  implicit none

  PetscErrorCode,intent(inout)                :: ierr
  real(kind=C_DOUBLE)                  :: rand1, rand2, rand3
  real(kind=C_DOUBLE)                  :: kmag, kmagsquared2d
  real(kind=C_DOUBLE)                  :: wavenumber_i,wavenumber_j,wavenumber_k
  real(kind=C_DOUBLE)                  :: factor, correction
  real(kind=C_DOUBLE)                  :: ek
  real(kind=C_DOUBLE),dimension(0:1)   :: alpha, beta
  real(kind=C_DOUBLE),dimension(0:1)   :: temp1, temp2, temp3
  real(kind=C_DOUBLE),dimension(0:1)   :: continuity
  integer                              :: i,j,k
  integer                              :: ii,jj,kk

  !> Variables for setting random seed
  integer                              :: seedsize 
  integer,dimension(:),allocatable     :: seed

  !> Set random seed on each rank. This is crucial for
  !! conjugate symmetry on the i=0 plane.
  call g_RandomInitFluid ( ierr )

!  call random_seed ( size=seedsize )
!  allocate ( seed(seedsize) )
!  IFSEEDSIZE: if ( RANDOM_SEED_SIZE <= seedsize ) then
!    seed(:) = 0
!    seed(1:RANDOM_SEED_SIZE) = RANDOM_SEED_FLUID(1:RANDOM_SEED_SIZE)
!  else
!    seed(:) = 0
!    seed(1:seedsize) = RANDOM_SEED_FLUID(1:seedsize)
!  end if IFSEEDSIZE

!  call random_seed(put=seed)
!  write(*,*) ' Random seed on rank ', MYRANK, ': ', seed

  !> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

  !> generate random number on all processes for every wave mode
  !! regardless of where that mode is stored.  this is to ensure that 
  !! conjugate symmetry across processes is maintained

  !> On i=0 plane go through the positive indices and derive
  !! conjugate symmetric values from there
  wavenumber_i = f_ArraysWavenumber ( 0, KI_MAX, K0I )
  DOI0J: do j = 0, KJ_MAX
    wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
    DOI0K: do k = 0, KK_MAX
      wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )
      ! generate rand1,rand2,rand3
      call random_number ( rand1 )
      call random_number ( rand2 )
      call random_number ( rand3 )

      ! multiply by 2pi
      rand1 = TWOPI * rand1
      rand2 = TWOPI * rand2
      rand3 = TWOPI * rand3

!      write(*,*) ' wavenumber_j, wavenumber_k, MYRANK, rand1, rand2, rand3: ', &
!                 wavenumber_j, wavenumber_k, MYRANK, rand1, rand2, rand3

      ! find magnitude of wave vector for current mode
      kmag = sqrt( (wavenumber_i * wavenumber_i) + &
                   (wavenumber_j * wavenumber_j) + &
                   (wavenumber_k * wavenumber_k) )

      kmagsquared2d = (wavenumber_i * wavenumber_i) + &
                      (wavenumber_j * wavenumber_j)

      ! find energy from prescribed energy spectrum
      call f_InitialiseSpec ( kmag, ek, ierr )
 
      ! calculate factor and ignore zero magnitude wave vector
      IFJKNOTZERO: if ( kmag > ZERO ) then
          factor = sqrt ( ek / (FOURPI * kmag * kmag) )

        ! form complex random variable alpha
        alpha(0) = factor * cos(rand1) * cos(rand3)
        alpha(1) = factor * sin(rand1) * cos(rand3)

        ! form complex random variable beta
        beta(0)  = factor * cos(rand2) * sin(rand3)
        beta(1)  = factor * sin(rand2) * sin(rand3)

        ! calculate random fourier mode according to rogallo method
        IFJKK1K2: if ( kmagsquared2d > ZERO ) then
          temp1 = ( alpha * kmag * wavenumber_j + beta * wavenumber_i * wavenumber_k ) / &
                  ( kmag * sqrt( kmagsquared2d ) ) 

          temp2 = ( beta * wavenumber_j * wavenumber_k - alpha * kmag * wavenumber_i ) / &
                  ( kmag * sqrt( kmagsquared2d ) )
   
          temp3 = - ( beta * sqrt( kmagsquared2d ) ) / kmag

        else

          temp1 = ZERO
          temp2 = ZERO
          temp3 = ZERO

        end if IFJKK1K2

      else

        temp1 = ZERO
        temp2 = ZERO
        temp3 = ZERO

      end if IFJKNOTZERO

      ! now check to see if the current process has the mode being broadcast
      ! note: because we only store half of the fourier modes there is no
      !       symmetry to enforce except on the i=0 plane

      ! complex conjugate symmetry does exist on the i=0 plane so we have to
      ! make sure that uk(0,j,k)=uk*(0,-j,-k)

      IFJJ0: if ( j /= 0 ) then
        jj = NODES_KJ - j !+ 1
      else
        jj = j
      end if IFJJ0
      IFKK0: if ( k /= 0 ) then
        kk = NODES_KK - k !+ 1
      else
        kk = k
      end if IFKK0

      !> If i = 0 is on this node, we know that i_min_3w =0
      IFJKKIMIN: if( i_min_3w == 0 ) then
        IFJKKJMIN: if ( j >= j_min_3w ) then
          IFJKKJMAX: if ( j <= j_max_3w ) then
            IFJKKKMIN: if ( k >= k_min_3w ) then
              IFJKKKMAX: if ( k <= k_max_3w ) then

!                arr_u1_3w(:,k,j,0) = ZERO
!                arr_u2_3w(:,k,j,0) = temp2
!                arr_u3_3w(:,k,j,0) = ZERO
                arr_u1_3w(:,k,j,0) = temp1
                arr_u2_3w(:,k,j,0) = temp2
                arr_u3_3w(:,k,j,0) = temp3

              end if IFJKKKMAX
            end if IFJKKKMIN
          end if IFJKKJMAX
        end if IFJKKJMIN
      end if IFJKKIMIN

      IFJJKKKIMIN: if( i_min_3w == 0 ) then
        IFJJKKKJMIN: if ( jj >= j_min_3w ) then
          IFJJKKKJMAX: if ( jj <= j_max_3w ) then
            IFJJKKKKMIN: if ( kk >= k_min_3w ) then
              IFJJKKKKMAX: if ( kk <= k_max_3w ) then

!                arr_u1_3w(0,kk,jj,0) = ZERO
!                arr_u2_3w(0,kk,jj,0) = temp2(0)
!                arr_u3_3w(0,kk,jj,0) = ZERO
!                arr_u1_3w(1,kk,jj,0) = - ZERO
!                arr_u2_3w(1,kk,jj,0) = - temp2(1)
!                arr_u3_3w(1,kk,jj,0) = - ZERO
                arr_u1_3w(0,kk,jj,0) = temp1(0)
                arr_u2_3w(0,kk,jj,0) = temp2(0)
                arr_u3_3w(0,kk,jj,0) = temp3(0)
                arr_u1_3w(1,kk,jj,0) = - temp1(1)
                arr_u2_3w(1,kk,jj,0) = - temp2(1)
                arr_u3_3w(1,kk,jj,0) = - temp3(1)

              end if IFJJKKKKMAX
            end if IFJJKKKKMIN
          end if IFJJKKKJMAX
        end if IFJJKKKJMIN
      end if IFJJKKKIMIN

      IFK0: if ( k /= 0 ) then
        !> Now if k /= 0 still need to find (-k, j) and (k, -j).
        !! Arbitrarily choose negative j.
        wavenumber_j = f_ArraysWavenumber ( jj, KJ_MAX, K0J )
        IFBCONJ: if ( B220 > 1 ) then
          wavenumber_j = wavenumber_j / real(B220,kind=C_DOUBLE)
        end if IFBCONJ

        ! generate rand1,rand2,rand3
        call random_number ( rand1 )
        call random_number ( rand2 )
        call random_number ( rand3 )

        ! multiply by 2pi
        rand1 = TWOPI * rand1
        rand2 = TWOPI * rand2
        rand3 = TWOPI * rand3

        ! calculate factor and ignore zero magnitude wave vector
        IFJJKNOTZERO: if ( kmag > ZERO ) then
          factor = sqrt ( ek / (FOURPI * kmag * kmag) )

          ! form complex random variable alpha
          alpha(0) = factor * cos(rand1) * cos(rand3)
          alpha(1) = factor * sin(rand1) * cos(rand3)

          ! form complex random variable beta
          beta(0)  = factor * cos(rand2) * sin(rand3)
          beta(1)  = factor * sin(rand2) * sin(rand3)

          ! calculate random fourier mode according to rogallo method
          IFJJKK1K2: if ( kmagsquared2d > ZERO ) then
            temp1 = ( alpha * kmag * wavenumber_j + beta * wavenumber_i * wavenumber_k ) / &
                    ( kmag * sqrt( kmagsquared2d ) ) 

            temp2 = ( beta * wavenumber_j * wavenumber_k - alpha * kmag * wavenumber_i ) / &
                    ( kmag * sqrt( kmagsquared2d ) )
   
            temp3 = - ( beta * sqrt( kmagsquared2d ) ) / kmag

          else

            temp1 = ZERO
            temp2 = ZERO
            temp3 = ZERO

          end if IFJJKK1K2

        else

          temp1 = ZERO
          temp2 = ZERO
          temp3 = ZERO

        end if IFJJKNOTZERO
  
        !> If i = 0 is on this node, we know that i_min_3w =0
        IFJJKKIMIN: if( i_min_3w == 0 ) then
          IFJJKKJMIN: if ( jj >= j_min_3w ) then
            IFJJKKJMAX: if ( jj <= j_max_3w ) then
              IFJJKKKMIN: if ( k >= k_min_3w ) then
                IFJJKKKMAX: if ( k <= k_max_3w ) then

!                  arr_u1_3w(:,k,jj,0) = ZERO
!                  arr_u2_3w(:,k,jj,0) = temp2
!                  arr_u3_3w(:,k,jj,0) = ZERO
                  arr_u1_3w(:,k,jj,0) = temp1
                  arr_u2_3w(:,k,jj,0) = temp2
                  arr_u3_3w(:,k,jj,0) = temp3

                end if IFJJKKKMAX
              end if IFJJKKKMIN
            end if IFJJKKJMAX
          end if IFJJKKJMIN
        end if IFJJKKIMIN

        IFJKKKIMIN: if( i_min_3w == 0 ) then
          IFJKKKJMIN: if ( j >= j_min_3w ) then
            IFJKKKJMAX: if ( j <= j_max_3w ) then
              IFJKKKKMIN: if ( kk >= k_min_3w ) then
                IFJKKKKMAX: if ( kk <= k_max_3w ) then

                  !> Conjugate symmetric equivalents
!                  arr_u1_3w(0,kk,j,0) = ZERO
!                  arr_u1_3w(1,kk,j,0) = - ZERO
!                  arr_u2_3w(0,kk,j,0) = temp2(0)
!                  arr_u2_3w(1,kk,j,0) = - temp2(1)
!                  arr_u3_3w(0,kk,j,0) = ZERO
!                  arr_u3_3w(1,kk,j,0) = - ZERO
                  arr_u1_3w(0,kk,j,0) = temp1(0)
                  arr_u1_3w(1,kk,j,0) = - temp1(1)
                  arr_u2_3w(0,kk,j,0) = temp2(0)
                  arr_u2_3w(1,kk,j,0) = - temp2(1)
                  arr_u3_3w(0,kk,j,0) = temp3(0)
                  arr_u3_3w(1,kk,j,0) = - temp3(1)

                end if IFJKKKKMAX
              end if IFJKKKKMIN
            end if IFJKKKJMAX
          end if IFJKKKJMIN
        end if IFJKKKIMIN
      end if IFK0

!      write(*,*) ' f_ArraysWavenumber ( j, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( jj, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( k, KJ_MAX, K0K ), &
!                   f_ArraysWavenumber ( kk, KK_MAX, K0K ): ', &
!                   f_ArraysWavenumber ( j, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( jj, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( k, KJ_MAX, K0K ), &
!                   f_ArraysWavenumber ( kk, KK_MAX, K0K ), &
!                 ' arr_u1_3w(:,k,j,0), arr_u1_3w(:,kk,jj,0), &
!                   arr_u1_3w(:,k,jj,0), arr_u1_3w(:,kk,j,0): ', &
!                   arr_u1_3w(:,k,j,0), arr_u1_3w(:,kk,jj,0), &
!                   arr_u1_3w(:,k,jj,0), arr_u1_3w(:,kk,j,0), &
!                 ' arr_u2_3w(:,k,j,0), arr_u2_3w(:,kk,jj,0), &
!                   arr_u2_3w(:,k,jj,0), arr_u2_3w(:,kk,j,0): ', &
!                   arr_u2_3w(:,k,j,0), arr_u2_3w(:,kk,jj,0), &
!                   arr_u2_3w(:,k,jj,0), arr_u2_3w(:,kk,j,0), &
!                 ' arr_u3_3w(:,k,j,0), arr_u3_3w(:,kk,jj,0), &
!                   arr_u3_3w(:,k,jj,0), arr_u3_3w(:,kk,j,0): ', &
!                   arr_u3_3w(:,k,j,0), arr_u3_3w(:,kk,jj,0), &
!                   arr_u3_3w(:,k,jj,0), arr_u3_3w(:,kk,j,0)


    end do DOI0K
  end do DOI0J



  !> Same for other nodes, where conjugate symmetry does not need to be enforced.
  DOI: do i = 1, NODES_KI
    wavenumber_i = f_ArraysWavenumber ( i, KI_MAX, K0I )
    DOJ: do j = 0, NODES_KJ
      wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      DOK: do k = 0, NODES_KK
        wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )
        ! generate rand1,rand2,rand3
        call random_number ( rand1 )
        call random_number ( rand2 )
        call random_number ( rand3 )

        ! multiply by 2pi
        rand1 = TWOPI * rand1
        rand2 = TWOPI * rand2
        rand3 = TWOPI * rand3

        ! find magnitude of wave vector for current mode
        kmag = sqrt( (wavenumber_i * wavenumber_i) + &
                     (wavenumber_j * wavenumber_j) + &
                     (wavenumber_k * wavenumber_k) )

        kmagsquared2d = (wavenumber_i * wavenumber_i) + &
                        (wavenumber_j * wavenumber_j)

        ! find energy from prescribed energy spectrum
        call f_InitialiseSpec ( kmag, ek, ierr )
 
        ! calculate factor and ignore zero magnitude wave vector
        IFNOTZERO: if ( kmag > ZERO ) then

          factor = sqrt ( ek / (FOURPI * kmag * kmag) )

          ! form complex random variable alpha
          alpha(0) = factor * cos(rand1) * cos(rand3)
          alpha(1) = factor * sin(rand1) * cos(rand3)

          ! form complex random variable beta
          beta(0)  = factor * cos(rand2) * sin(rand3)
          beta(1)  = factor * sin(rand2) * sin(rand3)

        ! calculate random fourier mode according to rogallo method
          IFK1K2: if ( kmagsquared2d > ZERO ) then
            temp1 = ( alpha * kmag * wavenumber_j + beta * wavenumber_i * wavenumber_k ) / &
                    ( kmag * sqrt( kmagsquared2d ) ) 

            temp2 = ( beta * wavenumber_j * wavenumber_k - alpha * kmag * wavenumber_i ) / &
                    ( kmag * sqrt( kmagsquared2d ) )
   
            temp3 = - ( beta * sqrt( kmagsquared2d ) ) / kmag

          else

            temp1 = ZERO
            temp2 = ZERO
            temp3 = ZERO

          end if IFK1K2

        else

          temp1 = ZERO
          temp2 = ZERO
          temp3 = ZERO

        end if IFNOTZERO

        ! now check to see if the current process owns this mode
        ! note: because we only store half of the fourier modes there is no
        !       symmetry to enforce except on the i=0 plane
        IFIKIMIN: if( i >= i_min_3w ) then
          IFIKIMAX: if ( i <= i_max_3w ) then
            IFIKJMIN: if ( j >= j_min_3w ) then
              IFIKJMAX: if ( j <= j_max_3w ) then
                IFIKKMIN: if ( k >= k_min_3w ) then
                  IFIKKMAX: if ( k <= k_max_3w ) then

                    arr_u1_3w(:,k,j,i) = temp1
                    arr_u2_3w(:,k,j,i) = temp2
                    arr_u3_3w(:,k,j,i) = temp3

                  end if IFIKKMAX
                end if IFIKKMIN
              end if IFIKJMAX
            end if IFIKJMIN
          end if IFIKIMAX
        end if IFIKIMIN

      end do DOK
    end do DOJ
  end do DOI

!> Return velocity components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, ierr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u2_3w, arr_u2_3w, ierr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u3_3w, arr_u3_3w, ierr )

end subroutine f_InitialiseFluid
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_InitialiseFluidDeformed
!> Initialise fluid velocity as a solenoidal field with a given energy spectrum.
!> @param ierr should return 0
subroutine f_InitialiseFluidDeformed ( ierr )

  use g_parameters, only : MYRANK, MASTER, NODES_X, NODES_Y, NODES_Z, &
                           RANDOM_SEED_SIZE, RANDOM_SEED_FLUID, &
                           B110, B220, B330 
  use g_random,     only : g_RandomInitFluid
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bvalexist, bmat
  use g_constants,  only : ZERO, ONE, HALF, TWOPI, FOURPI
  use f_arrays,     only : f_ArraysWavenumber, i_min_3w, j_min_3w, k_min_3w, &
                           i_max_3w, j_max_3w, k_max_3w, da3w, u1_3w, u2_3w, &
                           u3_3w, arr_u1_3w, arr_u2_3w, arr_u3_3w

  implicit none

  PetscErrorCode,intent(inout)                               :: ierr
  real(kind=C_DOUBLE)                                 :: rand1, rand2, rand3
  real(kind=C_DOUBLE)                                 :: kmag, kmaglab, kmagsquared2d
  real(kind=C_DOUBLE)                                 :: wavenumber_i,wavenumber_j,wavenumber_k
  real(kind=C_DOUBLE)                                 :: wvncentre_i,wvncentre_j,wvncentre_k
  real(kind=C_DOUBLE)                                 :: factor, correction
  real(kind=C_DOUBLE)                                 :: ek
  real(kind=C_DOUBLE),dimension(0:1)                  :: alpha, beta
  real(kind=C_DOUBLE),dimension(0:1,1:3)              :: templab
  real(kind=C_DOUBLE),dimension(0:1,1:3)              :: tempdef
  real(kind=C_DOUBLE),dimension(0:1)                  :: continuity
  integer                                             :: i,j,k,m,n
  integer                                             :: ii,jj,kk

  !> Arrays for energy on all wave modes in the volume
  real(kind=C_DOUBLE),dimension(0:B110)               :: wavenumbers_i
  real(kind=C_DOUBLE),dimension(0:B220)               :: wavenumbers_j
  real(kind=C_DOUBLE),dimension(0:B330)               :: wavenumbers_k
  real(kind=C_DOUBLE),dimension(0:B330,0:B220,0:B110) :: kmags,kmags2dsquared
  real(kind=C_DOUBLE),dimension(0:B330,0:B220,0:B110) :: energy

  !> Variables for setting random seed
  integer                                             :: seedsize 
  integer,dimension(:),allocatable                    :: seed

  !> Set random seed on each rank. This is crucial for
  !! conjugate symmetry on the i=0 plane.
  call g_RandomInitFluid ( ierr )

!  call random_seed ( size=seedsize )
!  allocate ( seed(seedsize) )
!  IFSEEDSIZE: if ( RANDOM_SEED_SIZE <= seedsize ) then
!    seed(:) = 0
!    seed(1:RANDOM_SEED_SIZE) = RANDOM_SEED_FLUID(1:RANDOM_SEED_SIZE)
!  else
!    seed(:) = 0
!    seed(1:seedsize) = RANDOM_SEED_FLUID(1:seedsize)
!  end if IFSEEDSIZE

!  call random_seed(put=seed)
!  write(*,*) ' Random seed on rank ', MYRANK, ': ', seed


  !> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

  !> generate random number on all processes for every wave mode
  !! regardless of where that mode is stored.  this is to ensure that 
  !! conjugate symmetry across processes is maintained


!  arr_u1_3w = ZERO
!  arr_u2_3w = ZERO
!  arr_u3_3w = ZERO

  !> On i=0 plane go through the positive indices and derive
  !! conjugate symmetric values from there
  wavenumber_i = f_ArraysWavenumber ( 0, KI_MAX, K0I )
  wavenumber_i = wavenumber_i / real(B110,kind=C_DOUBLE)
  DOI0J: do j = 0, KJ_MAX
    wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
    wavenumber_j = wavenumber_j / real(B220,kind=C_DOUBLE)
    DOI0K: do k = 0, KK_MAX
      wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )
      wavenumber_k = wavenumber_k / real(B330,kind=C_DOUBLE)
      !> generate rand1,rand2,rand3
      call random_number ( rand1 )
      call random_number ( rand2 )
      call random_number ( rand3 )

      !> multiply by 2pi
      rand1 = TWOPI * rand1
      rand2 = TWOPI * rand2
      rand3 = TWOPI * rand3

      !> find magnitude of wave vector for current mode
      kmag = sqrt( (wavenumber_i * wavenumber_i) + &
                   (wavenumber_j * wavenumber_j) + &
                   (wavenumber_k * wavenumber_k) )

      kmagsquared2d = (wavenumber_i * wavenumber_i) + &
                      (wavenumber_j * wavenumber_j)

      !> find energy from prescribed energy spectrum
      call f_InitialiseSpec ( kmag, ek, ierr )
 
      !> calculate factor and ignore zero magnitude wave vector
      IFJKNOTZERO: if ( kmag > ZERO ) then
        factor = sqrt ( ek / (FOURPI * kmag * kmag) )

        !> Correct for size of real domain vs computational domain
        factor = factor / ( real(B110,kind=C_DOUBLE) &
                 * real(B220,kind=C_DOUBLE) &
                 * real(B330,kind=C_DOUBLE) )**0.5

        ! form complex random variable alpha
        alpha(0) = factor * cos(rand1) * cos(rand3)
        alpha(1) = factor * sin(rand1) * cos(rand3)

        ! form complex random variable beta
        beta(0)  = factor * cos(rand2) * sin(rand3)
        beta(1)  = factor * sin(rand2) * sin(rand3)

        ! calculate random fourier mode according to rogallo method
        IFJKK1K2: if ( kmagsquared2d > ZERO ) then
          templab(:,1) = ( alpha * kmag * wavenumber_j + beta * wavenumber_i * wavenumber_k ) / &
                         ( kmag * sqrt( kmagsquared2d ) ) 

          templab(:,2) = ( beta * wavenumber_j * wavenumber_k - alpha * kmag * wavenumber_i ) / &
                         ( kmag * sqrt( kmagsquared2d ) )
   
          templab(:,3) = - ( beta * sqrt( kmagsquared2d ) ) / kmag

        else

          templab = ZERO

        end if IFJKK1K2

      else

        templab = ZERO

      end if IFJKNOTZERO

      ! now check to see if the current process has the mode being broadcast
      ! note: because we only store half of the fourier modes there is no
      !       symmetry to enforce except on the i=0 plane

      ! complex conjugate symmetry does exist on the i=0 plane so we have to
      ! make sure that uk(0,j,k)=uk*(0,-j,-k)

      IFJJ0: if ( j /= 0 ) then
        jj = NODES_KJ - j !+ 1
      else
        jj = j
      end if IFJJ0
      IFKK0: if ( k /= 0 ) then
        kk = NODES_KK - k !+ 1
      else
        kk = k
      end if IFKK0

      !> If i = 0 is on this node, we know that i_min_3w =0
      IFJKKIMIN: if( i_min_3w == 0 ) then
        IFJKKJMIN: if ( j >= j_min_3w ) then
          IFJKKJMAX: if ( j <= j_max_3w ) then
            IFJKKKMIN: if ( k >= k_min_3w ) then
              IFJKKKMAX: if ( k <= k_max_3w ) then

                tempdef = ZERO

                DOM1: do m = 1, 3, 1
                  DON1: do n = 1, 3, 1
                    IFBEXIST1: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                      tempdef(:,m) = tempdef(:,m) + bmat(m,n) * templab(:,n)
                    end if IFBEXIST1
                  end do DON1
                end do DOM1


                arr_u1_3w(:,k,j,0) = tempdef(:,1)
                arr_u2_3w(:,k,j,0) = tempdef(:,2)
                arr_u3_3w(:,k,j,0) = tempdef(:,3)

              end if IFJKKKMAX
            end if IFJKKKMIN
          end if IFJKKJMAX
        end if IFJKKJMIN
      end if IFJKKIMIN

      IFJJKKKIMIN: if( i_min_3w == 0 ) then
        IFJJKKKJMIN: if ( jj >= j_min_3w ) then
          IFJJKKKJMAX: if ( jj <= j_max_3w ) then
            IFJJKKKKMIN: if ( kk >= k_min_3w ) then
              IFJJKKKKMAX: if ( kk <= k_max_3w ) then

                tempdef = ZERO

                DOM2: do m = 1, 3, 1
                  DON2: do n = 1, 3, 1
                    IFBEXIST2: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                      tempdef(:,m) = tempdef(:,m) + bmat(m,n) * templab(:,n)
                    end if IFBEXIST2
                  end do DON2
                end do DOM2


                arr_u1_3w(0,kk,jj,0) = tempdef(0,1)
                arr_u2_3w(0,kk,jj,0) = tempdef(0,2)
                arr_u3_3w(0,kk,jj,0) = tempdef(0,3)
                arr_u1_3w(1,kk,jj,0) = - tempdef(1,1)
                arr_u2_3w(1,kk,jj,0) = - tempdef(1,2)
                arr_u3_3w(1,kk,jj,0) = - tempdef(1,3)

              end if IFJJKKKKMAX
            end if IFJJKKKKMIN
          end if IFJJKKKJMAX
        end if IFJJKKKJMIN
      end if IFJJKKKIMIN

      IFK0: if ( k /= 0 ) then
        !> Now still need to find (-k, j) and (k, -j).
        !! Arbitrarily choose negative j.
        wavenumber_j = f_ArraysWavenumber ( jj, KJ_MAX, K0J )
        wavenumber_j = wavenumber_j / real(B220,kind=C_DOUBLE)

        ! generate rand1,rand2,rand3
        call random_number ( rand1 )
        call random_number ( rand2 )
        call random_number ( rand3 )

        ! multiply by 2pi
        rand1 = TWOPI * rand1
        rand2 = TWOPI * rand2
        rand3 = TWOPI * rand3

        ! calculate factor and ignore zero magnitude wave vector
        IFJJKNOTZERO: if ( kmag > ZERO ) then
          factor = sqrt ( ek / (FOURPI * kmag * kmag) )

          !> Correct for size of real domain vs computational domain
          factor = factor / ( real(B110,kind=C_DOUBLE) &
                   * real(B220,kind=C_DOUBLE) &
                   * real(B330,kind=C_DOUBLE) )**0.5

          ! form complex random variable alpha
          alpha(0) = factor * cos(rand1) * cos(rand3)
          alpha(1) = factor * sin(rand1) * cos(rand3)

          ! form complex random variable beta
          beta(0)  = factor * cos(rand2) * sin(rand3)
          beta(1)  = factor * sin(rand2) * sin(rand3)

          ! calculate random fourier mode according to rogallo method
          IFJJKK1K2: if ( kmagsquared2d > ZERO ) then
            templab(:,1) = ( alpha * kmag * wavenumber_j + beta * wavenumber_i * wavenumber_k ) / &
                           ( kmag * sqrt( kmagsquared2d ) ) 

            templab(:,2) = ( beta * wavenumber_j * wavenumber_k - alpha * kmag * wavenumber_i ) / &
                           ( kmag * sqrt( kmagsquared2d ) )
   
            templab(:,3) = - ( beta * sqrt( kmagsquared2d ) ) / kmag

          else
  
            templab = ZERO

          end if IFJJKK1K2

        else

          templab = ZERO

        end if IFJJKNOTZERO

        !> If i = 0 is on this node, we know that i_min_3w =0
        IFJJKKIMIN: if( i_min_3w == 0 ) then
          IFJJKKJMIN: if ( jj >= j_min_3w ) then
            IFJJKKJMAX: if ( jj <= j_max_3w ) then
              IFJJKKKMIN: if ( k >= k_min_3w ) then
                IFJJKKKMAX: if ( k <= k_max_3w ) then

                  tempdef = ZERO
  
                  DOM3: do m = 1, 3, 1
                    DON3: do n = 1, 3, 1
                      IFBEXIST3: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                        tempdef(:,m) = tempdef(:,m) + bmat(m,n) * templab(:,n)
                      end if IFBEXIST3
                    end do DON3
                  end do DOM3


                  arr_u1_3w(:,k,jj,0) = tempdef(:,1)
                  arr_u2_3w(:,k,jj,0) = tempdef(:,2)
                  arr_u3_3w(:,k,jj,0) = tempdef(:,3)

                end if IFJJKKKMAX
              end if IFJJKKKMIN
            end if IFJJKKJMAX
          end if IFJJKKJMIN
        end if IFJJKKIMIN

        IFJKKKIMIN: if( i_min_3w == 0 ) then
          IFJKKKJMIN: if ( j >= j_min_3w ) then
            IFJKKKJMAX: if ( j <= j_max_3w ) then
              IFJKKKKMIN: if ( kk >= k_min_3w ) then
                IFJKKKKMAX: if ( kk <= k_max_3w ) then

                  tempdef = ZERO

                  DOM4: do m = 1, 3, 1
                    DON4: do n = 1, 3, 1
                      IFBEXIST4: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                        tempdef(:,m) = tempdef(:,m) + bmat(m,n) * templab(:,n)
                      end if IFBEXIST4
                    end do DON4
                  end do DOM4


                  !> Conjugate symmetric equivalents
                  arr_u1_3w(0,kk,j,0) = tempdef(0,1)
                  arr_u1_3w(1,kk,j,0) = - tempdef(1,1)
                  arr_u2_3w(0,kk,j,0) = tempdef(0,2)
                  arr_u2_3w(1,kk,j,0) = - tempdef(1,2)
                  arr_u3_3w(0,kk,j,0) = tempdef(0,3)
                  arr_u3_3w(1,kk,j,0) = - tempdef(1,3)

                end if IFJKKKKMAX
              end if IFJKKKKMIN
            end if IFJKKKJMAX
          end if IFJKKKJMIN
        end if IFJKKKIMIN
      end if IFK0

!      write(*,*) ' f_ArraysWavenumber ( j, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( jj, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( k, KJ_MAX, K0K ), &
!                   f_ArraysWavenumber ( kk, KK_MAX, K0K ): ', &
!                   f_ArraysWavenumber ( j, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( jj, KJ_MAX, K0J ), &
!                   f_ArraysWavenumber ( k, KJ_MAX, K0K ), &
!                   f_ArraysWavenumber ( kk, KK_MAX, K0K ), &
!                 ' arr_u1_3w(:,k,j,0), arr_u1_3w(:,kk,jj,0), &
!                   arr_u1_3w(:,k,jj,0), arr_u1_3w(:,kk,j,0): ', &
!                   arr_u1_3w(:,k,j,0), arr_u1_3w(:,kk,jj,0), &
!                   arr_u1_3w(:,k,jj,0), arr_u1_3w(:,kk,j,0), &
!                 ' arr_u2_3w(:,k,j,0), arr_u2_3w(:,kk,jj,0), &
!                   arr_u2_3w(:,k,jj,0), arr_u2_3w(:,kk,j,0): ', &
!                   arr_u2_3w(:,k,j,0), arr_u2_3w(:,kk,jj,0), &
!                   arr_u2_3w(:,k,jj,0), arr_u2_3w(:,kk,j,0), &
!                 ' arr_u3_3w(:,k,j,0), arr_u3_3w(:,kk,jj,0), &
!                   arr_u3_3w(:,k,jj,0), arr_u3_3w(:,kk,j,0): ', &
!                   arr_u3_3w(:,k,j,0), arr_u3_3w(:,kk,jj,0), &
!                   arr_u3_3w(:,k,jj,0), arr_u3_3w(:,kk,j,0)
    end do DOI0K
  end do DOI0J



  !> Same for other nodes, where conjugate symmetry does not need to be enforced.
  DOI: do i = 1, NODES_KI
    wavenumber_i = f_ArraysWavenumber ( i, KI_MAX, K0I )
    wavenumber_i = wavenumber_i / real(B110,kind=C_DOUBLE)
    DOJ: do j = 0, NODES_KJ
      wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      wavenumber_j = wavenumber_j / real(B220,kind=C_DOUBLE)
      DOK: do k = 0, NODES_KK
        wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )
        wavenumber_k = wavenumber_k / real(B330,kind=C_DOUBLE)
        ! generate rand1,rand2,rand3
        call random_number ( rand1 )
        call random_number ( rand2 )
        call random_number ( rand3 )

        ! multiply by 2pi
        rand1 = TWOPI * rand1
        rand2 = TWOPI * rand2
        rand3 = TWOPI * rand3

        ! find magnitude of wave vector for current mode
        kmag = sqrt( (wavenumber_i * wavenumber_i) + &
                     (wavenumber_j * wavenumber_j) + &
                     (wavenumber_k * wavenumber_k) )

        kmagsquared2d = (wavenumber_i * wavenumber_i) + &
                        (wavenumber_j * wavenumber_j)

        ! find energy from prescribed energy spectrum
        call f_InitialiseSpec ( kmag, ek, ierr )
 
        ! calculate factor and ignore zero magnitude wave vector
        IFNOTZERO: if ( kmag > ZERO ) then

          factor = sqrt ( ek / (FOURPI * kmag * kmag) )

          !> Correct for size of real domain vs computational domain
          factor = factor / ( real(B110,kind=C_DOUBLE) &
                   * real(B220,kind=C_DOUBLE) &
                   * real(B330,kind=C_DOUBLE) )**0.5

          ! form complex random variable alpha
          alpha(0) = factor * cos(rand1) * cos(rand3)
          alpha(1) = factor * sin(rand1) * cos(rand3)

          ! form complex random variable beta
          beta(0)  = factor * cos(rand2) * sin(rand3)
          beta(1)  = factor * sin(rand2) * sin(rand3)

        ! calculate random fourier mode according to rogallo method
          IFK1K2: if ( kmagsquared2d > ZERO ) then
            templab(:,1) = ( alpha * kmag * wavenumber_j + beta * wavenumber_i * wavenumber_k ) / &
                           ( kmag * sqrt( kmagsquared2d ) ) 

            templab(:,2) = ( beta * wavenumber_j * wavenumber_k - alpha * kmag * wavenumber_i ) / &
                           ( kmag * sqrt( kmagsquared2d ) )
   
            templab(:,3) = - ( beta * sqrt( kmagsquared2d ) ) / kmag

          else

            templab = ZERO

          end if IFK1K2

        else

          templab = ZERO

        end if IFNOTZERO

        ! now check to see if the current process owns this mode
        ! note: because we only store half of the fourier modes there is no
        !       symmetry to enforce except on the i=0 plane
        IFIKIMIN: if( i >= i_min_3w ) then
          IFIKIMAX: if ( i <= i_max_3w ) then
            IFIKJMIN: if ( j >= j_min_3w ) then
              IFIKJMAX: if ( j <= j_max_3w ) then
                IFIKKMIN: if ( k >= k_min_3w ) then
                  IFIKKMAX: if ( k <= k_max_3w ) then

                    tempdef = ZERO

                    DOM5: do m = 1, 3, 1
                      DON5: do n = 1, 3, 1
                        IFBEXIST5: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                          tempdef(:,m) = tempdef(:,m) + bmat(m,n) * templab(:,n)
                        end if IFBEXIST5
                      end do DON5
                    end do DOM5


                    arr_u1_3w(:,k,j,i) = tempdef(:,1)
                    arr_u2_3w(:,k,j,i) = tempdef(:,2)
                    arr_u3_3w(:,k,j,i) = tempdef(:,3)

                  end if IFIKKMAX
                end if IFIKKMIN
              end if IFIKJMAX
            end if IFIKJMIN
          end if IFIKIMAX
        end if IFIKIMIN

      end do DOK
    end do DOJ
  end do DOI

  !> Return velocity components to PETSc
  call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, ierr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u2_3w, arr_u2_3w, ierr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u3_3w, arr_u3_3w, ierr )

end subroutine f_InitialiseFluidDeformed
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
!> @param kmag wavenumber of the wave mode
!> @param ek energy on the wave mode
!> @param ierr should return 0
subroutine f_InitialiseSpec ( kmag, ek, ierr )

  use g_parameters, only : F_INIT_TYPE, F_INIT_NULL, F_INIT_GAUSS, &
                           F_INIT_KERR, F_INIT_PULSE, F_INIT_SQUARE, &
                           F_INIT_MINUSFIVETHIRD, F_INIT_MATSUMOTO, &
                           F_INIT_AHMED_ELGHOBASHI, &
                           F_INIT_KSTART, F_INIT_KEND, F_INIT_KV, &
                           F_INIT_KP, F_INIT_U0, F_INIT_EPSILON
  use g_constants,  only : ZERO, ONE, TWO, FOUR, THREE_OVER_TWO, TWO_OVER_THREE, &
                           FIVE_OVER_THREE, ELEVEN, ELEVEN_OVER_THREE, ONEPI, &
                           TWELVE, SIXTEEN

  implicit none

  PetscErrorCode,intent(inout)             :: ierr

  real(kind=C_DOUBLE),intent(in)    :: kmag
  real(kind=C_DOUBLE),intent(out)   :: ek

  select case ( F_INIT_TYPE )
    !-------------------------------
    ! init_dist not defined
    case ( F_INIT_NULL )
      ek = ZERO
      !-------------------------------
      ! standard gaussian
    case ( F_INIT_GAUSS )
      ek = ( F_INIT_U0 / ( F_INIT_KV * sqrt ( FOUR * asin ( ONE ) ) ) ) *          &
             exp ( -( ( kmag - F_INIT_KP )**2 ) / ( TWO * ( F_INIT_KV**2 ) ) )
    !-------------------------------
    ! kerr
    case ( F_INIT_KERR )
      ek = ( F_INIT_U0 / ( FOUR * asin ( ONE ) ) ) * ( kmag**4 ) *           &
             exp ( -( TWO * kmag**2 ) / ( F_INIT_KP**2 ) )
    !-------------------------------
    ! pulse
    case( F_INIT_PULSE )
      if ( sqrt ( ( kmag - F_INIT_KP )**2 ) < F_INIT_KV ) then
        ek = F_INIT_U0
      else
        ek = ZERO
      end if
    !-------------------------------
    case( F_INIT_SQUARE )
      if ( kmag > F_INIT_KSTART ) then
        if ( kmag < F_INIT_KEND ) then
          ek = ONE
        else
          ek = ZERO
        end if
      else
        ek = ZERO
      end if
    !-------------------------------
    case( F_INIT_MINUSFIVETHIRD )
      if ( kmag < F_INIT_KSTART ) then
        ek = THREE_OVER_TWO * ( F_INIT_EPSILON**TWO_OVER_THREE ) * &
             F_INIT_KSTART**(-FIVE_OVER_THREE) * &
             ( ( kmag / F_INIT_KSTART )**2 )
      else
        if ( kmag < F_INIT_KEND ) then
          ek = THREE_OVER_TWO * ( F_INIT_EPSILON**TWO_OVER_THREE ) * &
               F_INIT_KSTART**(-FIVE_OVER_THREE) * &
               ( ( kmag / F_INIT_KSTART )**(-FIVE_OVER_THREE) )
        else
          ek = ZERO
        end if
      end if
    !-------------------------------
    case( F_INIT_MATSUMOTO )

      if ( kmag < F_INIT_KSTART ) then
        ek = F_INIT_U0 * ( kmag**2 )
      else
        if ( kmag < F_INIT_KEND ) then
          ek = F_INIT_U0 * ( F_INIT_KSTART**ELEVEN_OVER_THREE ) * ( kmag **(-FIVE_OVER_THREE) )
        else
          ek = ZERO
        end if
      end if
    !-------------------------------
    case( F_INIT_AHMED_ELGHOBASHI )

      ek = SIXTEEN * ( ( TWO / ONEPI )**0.5 ) * THREE_OVER_TWO * &
           F_INIT_U0 * ( ( kmag**4 ) / ( F_INIT_KP**5 ) ) * &
           exp ( -( TWO * kmag**2 ) / ( F_INIT_KP**2 ) )

    !-------------------------------
    end select

    ierr = 0

end subroutine f_InitialiseSpec
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_InitialiseInitTaylorGreen 
!> Initialise velocity with a Taylor Green function in wave space
!> @param ierr should return 0
subroutine f_InitialiseInitTaylorGreen ( ierr )

  use g_constants,  only: ZERO, EIGHT
  use g_domain,     only: NODES_KJ, NODES_KK 
  use g_parameters, only: F_INIT_A, F_INIT_B, F_INIT_C, &
                          F_INIT_I, F_INIT_J, F_INIT_K
  use f_arrays,     only: da3w, u1_3w, u2_3w, u3_3w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          i_min_3w, j_min_3w, k_min_3w, &
                          i_max_3w, j_max_3w, k_max_3w

  implicit none
 
  PetscErrorCode,intent(inout)     :: ierr

  integer                   :: i, j, k

  PetscErrorCode            :: perr

  !> Initialise all arrays to zero
  call VecSet ( u1_3w, ZERO, perr )
  call VecSet ( u2_3w, ZERO, perr )
  call VecSet ( u3_3w, ZERO, perr )

!> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

  DOIINIT: do i = i_min_3w, i_max_3w, 1  
    DOJINIT: do j = j_min_3w, j_max_3w, 1  
      DOKINIT: do k = k_min_3w, k_max_3w, 1  

        IFI: if (i == F_INIT_I) then

          IFJPOS: if (j == F_INIT_J) then
            IFJPKPOS: if (k == F_INIT_K) then
              arr_u1_3w(0,k,j,i) = F_INIT_A / EIGHT           
              arr_u2_3w(0,k,j,i) = F_INIT_B / EIGHT              
              arr_u3_3w(0,k,j,i) = F_INIT_C / EIGHT             
!              write(*,*) 'i,j,k: ', i,j,k, 'u1, u1, u3: ', &
!                         arr_u1_3w(0,k,j,i), arr_u2_3w(0,k,j,i), arr_u3_3w(0,k,j,i)
            end if IFJPKPOS
            IFJPKNEG: if (k == NODES_KK - F_INIT_K) then
              arr_u1_3w(0,k,j,i) = - F_INIT_A / EIGHT              
              arr_u2_3w(0,k,j,i) = - F_INIT_B / EIGHT              
              arr_u3_3w(0,k,j,i) = F_INIT_C / EIGHT             
!              write(*,*) 'i,j,k: ', i,j,k, 'u1, u1, u3: ', &
!                         arr_u1_3w(0,k,j,i), arr_u2_3w(0,k,j,i), arr_u3_3w(0,k,j,i)
            end if IFJPKNEG
          end if IFJPOS

          IFJNEG: if (j == NODES_KJ - F_INIT_J) then
            IFJNKPOS: if (k == F_INIT_K) then
              arr_u1_3w(0,k,j,i) = - F_INIT_A / EIGHT              
              arr_u2_3w(0,k,j,i) = F_INIT_B / EIGHT              
              arr_u3_3w(0,k,j,i) = - F_INIT_C / EIGHT            
!              write(*,*) 'i,j,k: ', i,j,k, 'u1, u1, u3: ', &
!                         arr_u1_3w(0,k,j,i), arr_u2_3w(0,k,j,i), arr_u3_3w(0,k,j,i)
            end if IFJNKPOS
            IFJNKNEG: if (k == NODES_KK - F_INIT_K) then
              arr_u1_3w(0,k,j,i) = F_INIT_A / EIGHT              
              arr_u2_3w(0,k,j,i) = - F_INIT_B / EIGHT              
              arr_u3_3w(0,k,j,i) = - F_INIT_C / EIGHT            
!              write(*,*) 'i,j,k: ', i,j,k, 'u1, u1, u3: ', &
!                         arr_u1_3w(0,k,j,i), arr_u2_3w(0,k,j,i), arr_u3_3w(0,k,j,i)
            end if IFJNKNEG
          end if IFJNEG

        end if IFI

      end do DOKINIT
    end do DOJINIT
  end do DOIINIT

!> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

end subroutine f_InitialiseInitTaylorGreen
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_InitialiseCopyRestart 
!> Initialise velocity with a Taylor Green function in wave space
!> @param ierr should return 0
subroutine f_InitialiseCopyRestart ( ierr )

  use g_constants,  only: ZERO, EIGHT
  use g_parameters, only: F_RESTART_SCALE, NODES_Z, NODES_Y, NODES_X
  use g_domain,     only: NODES_KK, NODES_KJ
  use f_arrays,     only: da3w, u1_3w, u2_3w, u3_3w, &
                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                          i_min_3w, j_min_3w, k_min_3w, &
                          i_max_3w, j_max_3w, k_max_3w, &
                          restartcolour, &
                          da3w_rs, u1_3w_rs, u2_3w_rs, u3_3w_rs, &
                          arr_u1_3w_rs, arr_u2_3w_rs, arr_u3_3w_rs

  implicit none
 
  PetscErrorCode,intent(inout)     :: ierr

  integer                   :: i, j, k
  integer                   :: imin, jmin, kmin
  integer                   :: imax, jmax, kmaxneg
  integer                   :: imaxpos, nkj, nkk
  integer                   :: kmaxpos, nkneg, jmaxpos, njneg
  integer                   :: diffnj, diffnk, kshift, jshift

  PetscErrorCode            :: perr

  !> Initialise all arrays to zero
  call VecSet ( u1_3w, ZERO, perr )
  call VecSet ( u2_3w, ZERO, perr )
  call VecSet ( u3_3w, ZERO, perr )

  !> Get write access to all velocity components
  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

  IFRESTART: if ( restartcolour == 1 ) then

    !> Information about wave modes in restart file
    nkk = 2 * ( ( NODES_Z / F_RESTART_SCALE ) / 3 ) + 1
    nkj = 2 * ( ( NODES_Y / F_RESTART_SCALE ) / 3 ) + 1

    diffnk = NODES_KK - nkk
    !> shift k restart file indices for negative wavenumbers
    kshift = diffnk
    diffnj= NODES_KJ - nkj

!    write(*,*) ' k nodes in domain: ', NODES_KK, ' k nodes in restart file: ', nkk, ' difference: ', diffnk
!    write(*,*) ' j nodes in domain: ', NODES_KJ, ' j nodes in restart file: ', nkj, ' difference: ', diffnj

    imaxpos = ( ( NODES_X / F_RESTART_SCALE ) / 3 ) 
    jmaxpos = ( ( NODES_Y / F_RESTART_SCALE ) / 3 )
    kmaxpos = ( ( NODES_Z / F_RESTART_SCALE ) / 3 )
    nkneg = ( ( NODES_Z / F_RESTART_SCALE ) / 3 )
    njneg = ( ( NODES_Y / F_RESTART_SCALE ) / 3 ) 

    !> In i direction the maximum on this process is either the last element
    !! or the highest wavenumber in the restart file
    imin = i_min_3w
    imax = min ( i_max_3w, imaxpos )

    !> In j direction first find out if this process contains positive or
    !! negative wavenumbers. Then decide what to copy.
    IFJISPOS: if ( j_min_3w <= jmaxpos ) then
      jmin = j_min_3w
      jmax = min ( j_max_3w, jmaxpos )
      !> same indices for restart file
      jshift = 0
    else
!      write(*,*) ' Negative. ', j_min_3w, j_max_3w 
!      write(*,*) ' Wavenumbers on here: ', NODES_KJ - j_min_3w, NODES_KJ - j_max_3w
      jmin = max ( j_min_3w, ( NODES_KJ - njneg ) )
      jmax = j_max_3w
      !> different indices for restart file
      jshift = diffnj
    endif IFJISPOS

!    write(*,*) ' i range on this node: ', imin, imax
!    write(*,*) ' j range on this node: ', jmin, jmax

    !> Get read access to all velocity components
    call DMDAVecGetArrayReadF90 ( da3w_rs, u1_3w_rs, arr_u1_3w_rs, perr )
    call DMDAVecGetArrayReadF90 ( da3w_rs, u2_3w_rs, arr_u2_3w_rs, perr )
    call DMDAVecGetArrayReadF90 ( da3w_rs, u3_3w_rs, arr_u3_3w_rs, perr )

!    write(*,*) ' NODES_KK - nkneg: ', NODES_KK - nkneg, ' k_max_3w: ', k_max_3w
!    write(*,*) ' k_min_3w: ', k_min_3w, ' kmaxpos: ', kmaxpos

    DOI: do i = imin, imax, 1  
      DOJ: do j = jmin, jmax, 1  
        DOKPOS: do k = k_min_3w, kmaxpos, 1  
          arr_u1_3w(:,k,j,i) = arr_u1_3w_rs(:,k,(j-jshift),i)
          arr_u2_3w(:,k,j,i) = arr_u2_3w_rs(:,k,(j-jshift),i)
          arr_u3_3w(:,k,j,i) = arr_u3_3w_rs(:,k,(j-jshift),i)
        end do DOKPOS
        DOKNEG: do k = ( NODES_KK - nkneg ), k_max_3w, 1  
          arr_u1_3w(:,k,j,i) = arr_u1_3w_rs(:,(k-kshift),(j-jshift),i)
          arr_u2_3w(:,k,j,i) = arr_u2_3w_rs(:,(k-kshift),(j-jshift),i)
          arr_u3_3w(:,k,j,i) = arr_u3_3w_rs(:,(k-kshift),(j-jshift),i)
        end do DOKNEG
      end do DOJ
    end do DOI

    !> Return all velocity components to PETSc
    call DMDAVecRestoreArrayReadF90 ( da3w_rs, u1_3w_rs, arr_u1_3w_rs, perr )
    call DMDAVecRestoreArrayReadF90 ( da3w_rs, u2_3w_rs, arr_u2_3w_rs, perr )
    call DMDAVecRestoreArrayReadF90 ( da3w_rs, u3_3w_rs, arr_u3_3w_rs, perr )

  end if IFRESTART

  !> Return all velocity components to PETSc
  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, perr )

end subroutine f_InitialiseCopyRestart
!---------------------------------------------------------------------------------

end module
