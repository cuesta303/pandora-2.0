!---------------------------------------------------------------------------------
! module f_fluidstats
!> Contains subroutines for computing fluid statistics in both Fourier space and
!! real space.
module f_fluidstats

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

!------- data section begins ----------------------
!----------------------------------------------------------
! spectrum arrays 
!----------------------------------------------------------
!> kinetic energy spectrum local
real(kind=C_DOUBLE),allocatable :: temp_ken(:) 
!>kinetic energy spectrum global
real(kind=C_DOUBLE),allocatable :: spec_ken(:) 
!>kinetic energy spectrum global
real(kind=C_DOUBLE),allocatable :: spec_ken_initial(:) 
!> kinetic energy spectrum local
real(kind=C_DOUBLE),allocatable :: temp_ken1(:) 
!>kinetic energy spectrum global
real(kind=C_DOUBLE),allocatable :: spec_ken1(:) 
!> kinetic energy spectrum local
real(kind=C_DOUBLE),allocatable :: temp_ken2(:) 
!>kinetic energy spectrum global
real(kind=C_DOUBLE),allocatable :: spec_ken2(:) 
!> kinetic energy spectrum local
real(kind=C_DOUBLE),allocatable :: temp_ken3(:) 
!>kinetic energy spectrum global
real(kind=C_DOUBLE),allocatable :: spec_ken3(:) 
!> local array for kinetic energy in streamwise direction
real(kind=C_DOUBLE),allocatable :: temp_ken11(:) 
!> local array for streamwise energy spectrum
real(kind=C_DOUBLE),allocatable :: temp_spec_1i(:) 
!> global array for streamwise energy spectrum
real(kind=C_DOUBLE),allocatable :: spec_1i(:) 
!> local array for streamwise energy spectrum
real(kind=C_DOUBLE),allocatable :: temp_spec_2i(:) 
!> global array for streamwise energy spectrum
real(kind=C_DOUBLE),allocatable :: spec_2i(:) 
!!> global array for kinetic energy on i=0 mode
!real(kind=C_DOUBLE),allocatable :: spec_ken11(:) 
!> kinetic energy spectrum local
real(kind=C_DOUBLE),allocatable :: temp_ken21(:) 
!>kinetic energy spectrum global
real(kind=C_DOUBLE),allocatable :: spec_ken21(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken31(:) 
!!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken31(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken12(:) 
!!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken12(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken22(:) 
!!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken22(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken32(:) 
!!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken32(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken13(:) 
!!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken13(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken23(:) 
!!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken23(:) 
!!> kinetic energy spectrum local
!real(kind=C_DOUBLE),allocatable :: temp_ken33(:) 
!>kinetic energy spectrum global
!real(kind=C_DOUBLE),allocatable :: spec_ken33(:) 
!> Streamwise autocorrelation
real(kind=C_DOUBLE),allocatable :: temp_r1i(:) 
real(kind=C_DOUBLE),allocatable :: r1i(:) 
!>dissipation spectrum local
real(kind=C_DOUBLE),allocatable :: temp_dis(:) 
!>dissipation spectrum global
real(kind=C_DOUBLE),allocatable :: spec_dis(:) 

!>dissipation spectrum local
real(kind=C_DOUBLE),allocatable :: temp_disskew(:) 
!>dissipation spectrum global
real(kind=C_DOUBLE),allocatable :: spec_disskew(:) 


!real(kind=pandora_kind),allocatable :: local_rate(:) 
!real(kind=pandora_kind),allocatable :: global_rate(:) 

!>mode density
real(kind=C_DOUBLE),allocatable :: dens_mod(:) 
!>mode density for 1-D spectra
real(kind=C_DOUBLE),allocatable :: dens_mod_i(:) 
real(kind=C_DOUBLE),allocatable :: dens_mod_j(:) 
real(kind=C_DOUBLE),allocatable :: dens_mod_k(:) 

!real(kind=pandora_kind),allocatable :: rxlocal(:,:)
!real(kind=pandora_kind),allocatable :: rxglobal(:,:)

!real(kind=pandora_kind),allocatable :: exlocal(:,:)
!real(kind=pandora_kind),allocatable :: exglobal(:,:)

!real(kind=pandora_kind)   ,allocatable,save :: fluid_der (:,:,:,:)
!complex(kind=pandora_kind),allocatable,save :: fluid_derk(:,:,:,:)

!real(kind=pandora_kind)   ,allocatable,save :: tempr(:,:,:,:)
!complex(kind=pandora_kind),allocatable,save :: tempk(:,:,:,:)

!real(kind=pandora_kind),allocatable,public,save :: fluid_ensxz(:,:)

!----------------------------------------------------------
!flow quantities 
!----------------------------------------------------------
real(kind=C_DOUBLE),save                :: u1u1, u2u2, u3u3, u1u2, u1u3, u2u3, umax
real(kind=C_DOUBLE),save                :: rmsvort1, rmsvort2, rmsvort3 
!> Skewness and kurtosis
real(kind=C_DOUBLE),save                :: skew1, skew2, skew3, flat1, flat2, flat3
real(kind=C_DOUBLE),save                :: skewvort1, skewvort2, skewvort3, flatvort1, flatvort2, flatvort3
!> Maximum velocity for cfl criterion
real(kind=C_DOUBLE),save,public         :: ucfl
!> Local maximum vorticity
real(kind=C_DOUBLE),save,dimension(1:3) :: lomegamax
!> Maximum vorticity
real(kind=C_DOUBLE),save,dimension(1:3) :: omegamax
!> highest non-zero wave mode
integer,save,public                     :: kkmax                 
!> flow turbulent kinetic energy
real(kind=C_DOUBLE),save                :: f_ken,  fn_ken       
!> flow turbulent kinetic energy
real(kind=C_DOUBLE),save                :: fn_ken1, fn_ken2, fn_ken3
!> flow turbulent kinetic energy
real(kind=C_DOUBLE),save                :: f_dis,  fn_dis       
!> flow turbulent dissipation skew
real(kind=C_DOUBLE),save                :: f_disskew,  fn_disskew       

!> minimum wave vector
real(kind=C_DOUBLE),save                :: k_min                
!> maximum wave vector
real(kind=C_DOUBLE),save                :: k_max                

!real(kind=C_DOUBLE),save            :: l_i11,  ln_i11       ! integral scale l11
!real(kind=C_DOUBLE),save            :: l_i22,  ln_i22       ! integral scale l22
!real(kind=C_DOUBLE),save            :: l_i33,  ln_i33       ! integral scale l33

!> flow turbulent intensity
real(kind=C_DOUBLE),save            :: fn_int       
!> flow turbulent intensity component 1
real(kind=C_DOUBLE),save            :: fn_int1 
!> flow turbulent intensity component 2
real(kind=C_DOUBLE),save            :: fn_int2
!> flow turbulent intensity component 3
real(kind=C_DOUBLE),save            :: fn_int3
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int       
!> integral scale component l1 
real(kind=C_DOUBLE),save            :: ln_int1
!> integral scale component l2
real(kind=C_DOUBLE),save            :: ln_int2
!> integral scale component l3
real(kind=C_DOUBLE),save            :: ln_int3
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_11       
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_12       
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_13    
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_21
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_22    
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_23    
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_31    
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_32    
!> integral scale l
real(kind=C_DOUBLE),save            :: ln_int_33    
!> dimensionless shear rate 
real(kind=C_DOUBLE),save            :: Sstar
!> energy partition parameter (Lee/Kim/Moin 1990) 
real(kind=C_DOUBLE),save            :: Kstar
!> eddy elongation parameter (Lee/Kim/Moin 1990)
real(kind=C_DOUBLE),save            :: Lstar
!> turbulent length scale
real(kind=C_DOUBLE),save            :: ln_tur       
!> turbulent time scale
real(kind=C_DOUBLE),save            :: tn_tur
!> Kolmogorov length scale
real(kind=C_DOUBLE),public,save     :: ln_kol       
!> Taylor length scale
real(kind=C_DOUBLE),save            :: ln_tay       

!> eddy time scale
real(kind=C_DOUBLE),save,public     :: tn_edy       
!> Kolmogorov time scale
real(kind=C_DOUBLE),save,public     :: tn_kol       

!real(kind=C_DOUBLE),save            :: r_int,  rn_int       ! integral Reynolds number
! Taylor Reynolds number
real(kind=C_DOUBLE),save            :: rn_tay       
!> Dissipation skew
real(kind=C_DOUBLE),save            :: disskew       
!> Dissipation constant
real(kind=C_DOUBLE),save            :: C_epsilon     

!> u1u1 Reynolds stress in Fourier space
real(kind=C_DOUBLE),save            :: fs_u1u1              
!> u2u2 Reynolds stress in Fourier space
real(kind=C_DOUBLE),save            :: fs_u2u2              
!> u3u3 Reynolds stress in Fourier space
real(kind=C_DOUBLE),save            :: fs_u3u3              
!> u1u2 Reynolds stress in Fourier space
real(kind=C_DOUBLE),save            :: fs_u1u2              
!> u1u3 Reynolds stress in Fourier space
real(kind=C_DOUBLE),save            :: fs_u1u3              
!> u2u3 Reynolds stress in Fourier space
real(kind=C_DOUBLE),save            :: fs_u2u3              

!> u1u1 Reynolds stress in real space
real(kind=C_DOUBLE),save            :: rs_u1u1              
!> u2u2 Reynolds stress in real space
real(kind=C_DOUBLE),save            :: rs_u2u2              
!> u3u3 Reynolds stress in real space
real(kind=C_DOUBLE),save            :: rs_u3u3              
!> u1u2 Reynolds stress in real space
real(kind=C_DOUBLE),save            :: rs_u1u2              
!> u1u3 Reynolds stress in real space
real(kind=C_DOUBLE),save            :: rs_u1u3              
!> u2u3 Reynolds stress in real space
real(kind=C_DOUBLE),save            :: rs_u2u3              
!> Longitudinal integral length scale in real space
real(kind=C_DOUBLE),save            :: rs_lint11              
real(kind=C_DOUBLE),save            :: lint11              
!> Skewness and kurtosis
real(kind=C_DOUBLE),save            :: rs_skew1, rs_skew2, rs_skew3, &
                                       rs_flat1, rs_flat2, rs_flat3
real(kind=C_DOUBLE),save            :: rs_rmsvort1, rs_rmsvort2, rs_rmsvort3, &
                                       rs_skewvort1, rs_skewvort2, rs_skewvort3, &
                                       rs_flatvort1, rs_flatvort2, rs_flatvort3

!> turbulent kinetic energy
real(kind=C_DOUBLE),save            :: tke                  
!real(kind=C_DOUBLE),public,save :: rs_ens

!real(kind=C_DOUBLE),public,save :: rs_e1e1
!real(kind=C_DOUBLE),public,save :: rs_e2e2
!real(kind=C_DOUBLE),public,save :: rs_e3e3

!real(kind=C_DOUBLE),save        :: e110
!real(kind=C_DOUBLE),save        :: e220
!real(kind=C_DOUBLE),save        :: e330           

!real(kind=C_DOUBLE),save        :: fn_trn

!real(kind=C_DOUBLE),save        :: su, ku

!real(kind=C_DOUBLE),save        :: u1g, u2g, u3g

!real(kind=C_DOUBLE),save        :: global_ke

public :: f_FluidstatsRsZero
public :: f_FluidstatsRs2D
public :: f_FluidstatsRsCollect
public :: f_FluidstatsFs2D
public :: f_FluidstatsFsCollect
public :: f_FluidstatsFs
public :: f_FluidstatsInit
public :: f_FluidstatsFinalise
public :: f_FluidstatsInitialSpectrum
!public :: fluid_analysis
!------- data section ends ----------------------

contains
!---------------------------------------------------------------------------------
! subroutine f_FluidstatsInit
!> Initialise fluid statistics module.
!> @param ierr should return 0
subroutine f_FluidstatsInit ( ierr )

  implicit none

  PetscErrorCode,intent(inout)            :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising f_fluidstats module \n', ierr )

! allocate arrays required for calculating fluid quantities
  call PetscPrintf( PETSC_COMM_WORLD, '    | Allocating fluid spectrum arrays \n', ierr )
  call f_FluidstatsAlloc ( ierr )

  call PetscPrintf( PETSC_COMM_WORLD, '    | Setting initial values \n', ierr )
! calculate mode density for each node
  call f_FluidstatsModeDensity ( ierr )
!    call kwave_set(ierr)

end subroutine f_FluidstatsInit
!----------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_FluidstatsInitialSpectrum 
!> Initialise fluid statistics module.
!> @param ierr should return 0
subroutine f_FluidstatsInitialSpectrum ( ierr )

  implicit none

  PetscErrorCode,intent(inout)            :: ierr

  spec_ken_initial(:) = spec_ken(:)

end subroutine f_FluidstatsInitialSpectrum
!----------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_FluidstatsAlloc 
!> Allocate all arrays needed to compute the fluid statistics
!> @param ierr should return 0
subroutine f_FluidstatsAlloc ( ierr )

  use g_parameters, only : F_TYPE, F_ISOTROPIC, B110, B220, B330, F_TYPE, F_ISOTROPIC, NODES_X
  use g_domain,     only : KI_MAX, KJ_MAX, KK_MAX !myrank, master, v_allocs
  use g_constants,  only : ZERO, ONE, THREE_OVER_TWO !, c_zero

  implicit none

  PetscErrorCode,intent(inout)            :: ierr
  integer                          :: alloc_stat


  IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
    !> Determine lowest mode
    k_min = ONE

    !> Determine highest non-zero mode
    k_max = real ( ( KI_MAX * KI_MAX ) + &
                   ( KJ_MAX * KJ_MAX ) + & 
                   ( KK_MAX * KK_MAX ), kind=C_DOUBLE ) 

    k_max = sqrt ( k_max )
  else
    !> Determine lowest mode
!    k_min = min ( ( ONE / real(B110,kind=C_DOUBLE) ), &
!                  ( ONE / real(B220,kind=C_DOUBLE) ), &
!                  ( ONE / real(B330,kind=C_DOUBLE) ) )
    k_min = ONE

    !> Determine highest non-zero mode
    k_max = real ( ( KI_MAX * KI_MAX / real(B110*B110,kind=C_DOUBLE) ) + &
                   ( KJ_MAX * KJ_MAX / real(B220*B220,kind=C_DOUBLE) ) + & 
                   ( KK_MAX * KK_MAX / real(B330*B330,kind=C_DOUBLE) ), kind=C_DOUBLE ) 

    k_max = sqrt ( k_max )
  end if IFISOTROPIC

  !> Index corresponding to highest wave mode
  kkmax = nint ( k_max / k_min )

  !> allocate global energy spectrum arrays
  allocate ( temp_ken(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_ken = ZERO
  allocate ( spec_ken(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_ken = ZERO
  allocate ( spec_ken_initial(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_ken_initial = ZERO
  allocate ( temp_ken1(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_ken1 = ZERO
  allocate ( spec_ken1(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_ken1 = ZERO
  allocate ( temp_ken2(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_ken2 = ZERO
  allocate ( spec_ken2(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_ken2 = ZERO
  allocate ( temp_ken3(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_ken3 = ZERO
  allocate ( spec_ken3(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_ken3 = ZERO

  !> allocate energy spectrum arrays for longitudinal statistics
  allocate ( temp_ken11(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_ken11 = ZERO
  !> Need energy on all modes, in particular the i=0 mode
  allocate ( temp_spec_1i(0:KI_MAX), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_spec_1i = ZERO
  allocate ( spec_1i(0:KI_MAX), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_1i = ZERO
  allocate ( temp_spec_2i(0:KI_MAX), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_spec_2i = ZERO
  allocate ( spec_2i(0:KI_MAX), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_2i = ZERO
!  allocate ( spec_ken11(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken11 = ZERO
!  allocate ( temp_ken12(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken12 = ZERO
!  allocate ( spec_ken12(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken12 = ZERO
!  allocate ( temp_ken13(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken13 = ZERO
!  allocate ( spec_ken13(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken13 = ZERO
  allocate ( temp_ken21(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_ken21 = ZERO
  allocate ( spec_ken21(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken21 = ZERO
!  allocate ( temp_ken22(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken22 = ZERO
!  allocate ( spec_ken22(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken22 = ZERO
!  allocate ( temp_ken23(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken23 = ZERO
!  allocate ( spec_ken23(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken23 = ZERO
!  allocate ( temp_ken31(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken31 = ZERO
!  allocate ( spec_ken31(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken31 = ZERO
!  allocate ( temp_ken32(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken32 = ZERO
!  allocate ( spec_ken32(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken32 = ZERO
!  allocate ( temp_ken33(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  temp_ken33 = ZERO
!  allocate ( spec_ken33(1:kkmax), stat=alloc_stat )
!  CHKERRQ ( ierr )
!  spec_ken33 = ZERO

  !> allocate autocorrelation arrays
  allocate ( temp_r1i(0:(NODES_X/2)), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_r1i = ZERO
  allocate ( r1i(0:(NODES_X/2)), stat=alloc_stat )
  CHKERRQ ( ierr )
  r1i = ZERO
!--------------------------------------------------- 

  !>  allocate local and global dissipation spectrum arrays
  allocate ( temp_dis(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_dis = ZERO
  allocate ( spec_dis(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_dis = ZERO
  allocate ( temp_disskew(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  temp_disskew = ZERO
  allocate( spec_disskew(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  spec_disskew = ZERO

!------------------------------------------------------

  !> allocate mode density array
  allocate ( dens_mod(1:kkmax), stat=alloc_stat )
  CHKERRQ ( ierr )
  dens_mod = ZERO

  IF1DSPEC: if ( F_TYPE /= F_ISOTROPIC ) then
    allocate ( dens_mod_i(1:kkmax), stat=alloc_stat )
    CHKERRQ ( ierr )
    dens_mod = ZERO
    allocate ( dens_mod_j(1:kkmax), stat=alloc_stat )
    CHKERRQ ( ierr )
    dens_mod = ZERO
    allocate ( dens_mod_k(1:kkmax), stat=alloc_stat )
    CHKERRQ ( ierr )
    dens_mod = ZERO
  end if IF1DSPEC

!50  format('allocating temp_disipation array on rank:    'i2'  ('i3')')
!70  format('allocating spec_dis array on rank:           'i2'  ('i3')')

!    if(v_allocs.gt.0)write(*,90)myrank,nodes_x
!    if(alloc_stat/=0)write(*,95)

!------------------------------------------------

!------------------------------------------------------

!    allocate(local_rate(1:nodes_x),stat=alloc_stat)
!    local_rate = zero

!    allocate(global_rate(1:nodes_x),stat=alloc_stat)
!    global_rate = zero
 
!    allocate(rxlocal(1:6,1:((nodes_x/2)+1)),stat=alloc_stat)
!    rxlocal = zero

!    allocate(rxglobal(1:6,1:((nodes_x/2)+1)),stat=alloc_stat)
!    rxglobal = zero

!    allocate(exlocal(1:6,1:((nodes_x/2)+1)),stat=alloc_stat)
!    exlocal = zero

!    allocate(exglobal(1:6,1:((nodes_x/2)+1)),stat=alloc_stat)
!    exglobal = zero

!    allocate(tempr(ri_min:ri_max,rj_min:rj_max,rk_min:rk_max,3),stat=alloc_stat)
!    tempr = zero

!    allocate(tempk(ki_min:ki_max,kj_min:kj_max,kk_min:kk_max,3),stat=alloc_stat)
!    tempk = c_zero

!    if(myrank.eq.master)then
!       allocate( fluid_ensxz(1:nodes_x,1:nodes_z) )
!       fluid_ensxz = zero
!    endif

!   allocate(fluid_der(ri_min:ri_max,rj_min:rj_max,rk_min:rk_max,1:3))
!   allocate(fluid_derk(ki_min:ki_max,kj_min:kj_max,kk_min:kk_max,1:3))

end subroutine f_FluidstatsAlloc
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFinalise
!> Finalise the fluid statistics module.
!> @param ierr should return 0
subroutine f_FluidstatsFinalise ( ierr )

  use g_parameters, only : MYRANK, MASTER
 
  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Finalising f_fluidstats module \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Deallocating fluid spectrum arrays \n', ierr )

!> Call subroutine to deallocate PETSc arrays for the fluid components
!! in wave space.
  call f_FluidstatsDealloc ( ierr )

end subroutine f_FluidstatsFinalise
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsDealloc
!> Deallocate fluid statistics arrays.
!> @param ierr should return 0
subroutine f_FluidstatsDealloc ( ierr )

  use g_parameters, only : F_TYPE, F_ISOTROPIC

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  !> deallocate global energy spectrum arrays
  deallocate ( temp_ken, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_ken, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_ken_initial, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_ken1, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_ken1, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_ken2, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_ken2, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_ken3, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_ken3, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_dis, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_dis, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_disskew, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_disskew, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( dens_mod, STAT=ierr )
  CHKERRQ ( ierr )

  !> deallocate arrays for longitudinal statistics
  deallocate ( temp_ken11, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_spec_1i, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_1i, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( temp_spec_2i, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_2i, STAT=ierr )
  CHKERRQ ( ierr )
!  deallocate ( spec_ken11, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( temp_ken12, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken12, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( temp_ken13, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken13, STAT=ierr )
!  CHKERRQ ( ierr )
  deallocate ( temp_ken21, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( spec_ken21, STAT=ierr )
  CHKERRQ ( ierr )
!  deallocate ( temp_ken22, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken22, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( temp_ken23, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken23, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( temp_ken31, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken31, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( temp_ken32, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken32, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( temp_ken33, STAT=ierr )
!  CHKERRQ ( ierr )
!  deallocate ( spec_ken33, STAT=ierr )
!  CHKERRQ ( ierr )

  !> deallocate autocorrelation arrays
  deallocate ( temp_r1i, STAT=ierr )
  CHKERRQ ( ierr )
  deallocate ( r1i, STAT=ierr )
  CHKERRQ ( ierr )

  IF1DSPEC: if ( F_TYPE /= F_ISOTROPIC ) then
    deallocate ( dens_mod_i, stat=ierr )
    CHKERRQ ( ierr )
    deallocate ( dens_mod_j, stat=ierr )
    CHKERRQ ( ierr )
    deallocate ( dens_mod_k, stat=ierr )
    CHKERRQ ( ierr )
  end if IF1DSPEC

end subroutine f_FluidstatsDealloc
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsModeDensity 
!> Find out how many wave modes are in each bin.
!> @param ierr should return 0
subroutine f_FluidstatsModeDensity ( ierr )

!    use d_arrays,    only : nodes_x,nodes_y,nodes_z 
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K
  use g_constants,  only : HALF, THREE, FOUR, ONEPI, TWO !one, two
  use f_arrays,     only : i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           f_ArraysWavenumber
!    use f_arrays,    only : fluid_k, ki_min, kj_min, kk_min, &
!                                     ki_max, kj_max, kk_max

  implicit none

  PetscErrorCode,intent(inout)      :: ierr
  integer                    :: alloc_stat
  integer                    :: i,j,k,kk
  integer                    :: num_modes
  integer                    :: num_temp
  real(kind=C_DOUBLE)        :: wavenumber_i, wavenumber_j, wavenumber_k, k_mag
  real(kind=C_DOUBLE)        :: bin_upper
  real(kind=C_DOUBLE)        :: bin_lower

  ! find appropriate bin for current mode
  DOBIN: do kk = 1, kkmax, 1

    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    ! find bin limits
    bin_upper = real(kk,kind=C_DOUBLE) + ( k_min * HALF )
    bin_lower = real(kk,kind=C_DOUBLE) - ( k_min * HALF )

    ! find number of modes in each bin
    DOI: do i = i_min_3w, i_max_3w, 1
    wavenumber_i = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOJ: do j = j_min_3w, j_max_3w, 1
      wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
        DOK: do k = k_min_3w, k_max_3w, 1
        wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )

        !> Compute magnitude of wave vector
        k_mag = sqrt ( (wavenumber_i*wavenumber_i) + &
                       (wavenumber_j*wavenumber_j) + &
                       (wavenumber_k*wavenumber_k) )

          !> Check to see if k is in current bin
          IFLOWER: if ( k_mag >= bin_lower ) then
            IFUPPER: if ( k_mag < bin_upper ) then

              !> Need to count twice for conjugate symmetry except i=0
              IFI0: if ( i /= 0 ) then
                num_temp = num_temp + 2
              else
                num_temp = num_temp + 1
              end if IFI0

            end if IFUPPER
          end if IFLOWER

        end do DOK
      end do DOJ
    end do DOI

    ! sum counters from all processes
    call mpi_allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                 ( FOUR / THREE * ONEPI *( bin_upper**3 - bin_lower**3 ) )

!    !> Correct for size of real domain vs computational domain
!    dens_mod(kk) = dens_mod(kk) *  real(B330,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE)  

  end do DOBIN

end subroutine f_FluidstatsModeDensity
!----------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsRsZero 
!> Set real-space statistics to zero.
!> @param ierr should return 0
subroutine f_FluidstatsRsZero ( ierr )

  use g_constants, only: ZERO
 
  implicit none
 
  PetscErrorCode,intent(inout)            :: ierr

  u1u1  = ZERO
  u2u2  = ZERO
  u3u3  = ZERO
  u1u2  = ZERO
  u1u3  = ZERO
  u2u3  = ZERO
  skew1 = ZERO
  skew2 = ZERO
  skew3 = ZERO
  flat1 = ZERO
  flat2 = ZERO
  flat3 = ZERO
  rmsvort1 = ZERO
  rmsvort2 = ZERO
  rmsvort3 = ZERO
  skewvort1 = ZERO
  skewvort2 = ZERO
  skewvort3 = ZERO
  flatvort1 = ZERO
  flatvort2 = ZERO
  flatvort3 = ZERO
  lomegamax = ZERO
  omegamax = ZERO
  umax  = ZERO
  lint11 = ZERO
  rs_lint11 = ZERO
  temp_r1i = ZERO
  r1i = ZERO

end subroutine f_FluidstatsRsZero
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsRs2D 
!> Compute real-space statistics locally on 2-D slices.
!> @param ierr should return 0
subroutine f_FluidstatsRs2D ( tstep, ierr )

  use g_constants,  only : ZERO
  use g_parameters, only : TSTEPS, TS_CORR, NODES_X, F_TYPE, F_ISOTROPIC
  use g_domain,     only : dx_laboratoryframe, bmat, binv, bvalexist
  use f_arrays,     only : arr_u1_2r, arr_u2_2r, arr_u3_2r, &
                           arr_u4_2r, arr_u5_2r, arr_u6_2r, &
                           x_min_2r, x_max_2r, x_width_2r, y_min_2r, y_max_2r 
 
  implicit none
 
  integer,intent(in)                 :: tstep 
  PetscErrorCode,intent(inout)              :: ierr
  integer                            :: i, j, ii, m, n 
  real(kind=C_DOUBLE)                :: u1, u2, u3, u, &
                                        omega1, omega2, omega3
  real(kind=C_DOUBLE),dimension(1:3) :: ulab, umov 
  real(kind=C_DOUBLE),dimension(1:3) :: omegaabs

  DOY: do j = y_min_2r, y_max_2r, 1
    DOX: do i = x_min_2r, x_max_2r, 1

      !> Read local velocity
      IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
        u1 = arr_u1_2r(i,j)
        u2 = arr_u2_2r(i,j)
        u3 = arr_u3_2r(i,j)
      else

        !> @todo this is only for the Bardino/Ferziger type equations
        umov(1) = arr_u1_2r(i,j)
        umov(2) = arr_u2_2r(i,j)
        umov(3) = arr_u3_2r(i,j)

        ulab = ZERO

        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              ulab(m) = ulab(m) + binv(m,n) * umov(n)
            end if IFBEXIST
          end do DON
        end do DOM

        u1 = ulab(1)
        u2 = ulab(2)
        u3 = ulab(3)
      end if IFISOTROPIC

      !> Compute Reynolds stresses
      u1u1 = u1u1 + ( u1 * u1 )
      u2u2 = u2u2 + ( u2 * u2 )
      u3u3 = u3u3 + ( u3 * u3 )
      u1u2 = u1u2 + ( u1 * u2 )
      u1u3 = u1u3 + ( u1 * u3 )
      u2u3 = u2u3 + ( u2 * u3 )

      !> Compute skewness and kurtosis;
      !! need to normalise by Reynolds stresses
      skew1 = skew1 + ( u1**3 )
      skew2 = skew2 + ( u2**3 )
      skew3 = skew3 + ( u3**3 )
      flat1 = flat1 + ( u1**4 )
      flat2 = flat2 + ( u2**4 )
      flat3 = flat3 + ( u3**4 )

      !> Compute quantities for vorticity
      !> Read local vorticity
      omega1 = arr_u4_2r(i,j)
      omega2 = arr_u5_2r(i,j)
      omega3 = arr_u6_2r(i,j)

      !> Compute square vorticity for RMS 
      rmsvort1 = rmsvort1 + ( omega1 * omega1 )
      rmsvort2 = rmsvort2 + ( omega2 * omega2 )
      rmsvort3 = rmsvort3 + ( omega3 * omega3 )

      !> Compute skewness and kurtosis;
      !! need to normalise after MPI_Sum 
      skewvort1 = skewvort1 + ( omega1**3 )
      skewvort2 = skewvort2 + ( omega2**3 )
      skewvort3 = skewvort3 + ( omega3**3 )
      flatvort1 = flatvort1 + ( omega1**4 )
      flatvort2 = flatvort2 + ( omega2**4 )
      flatvort3 = flatvort3 + ( omega3**4 )

      !> Determine maximum velocity for CFL
      u = abs(u1) + abs(u2) + abs(u3)
      IFMAX: if ( u > umax ) then
        umax = u
      end if IFMAX

      !> Determine maximum vorticity
      omegaabs(1) = abs(omega1)
      omegaabs(2) = abs(omega2)
      omegaabs(3) = abs(omega3)
      IFMAX1: if ( omegaabs(1) > lomegamax(1) ) then
        lomegamax(1) = omegaabs(1)
      end if IFMAX1
      IFMAX2: if ( omegaabs(2) > lomegamax(2) ) then
        lomegamax(2) = omegaabs(2)
      end if IFMAX2
      IFMAX3: if ( omegaabs(3) > lomegamax(3) ) then
        lomegamax(3) = omegaabs(3)
      end if IFMAX3

      IFCORR: if ( TS_CORR <= TSTEPS ) then
        IFCORRTS: if ( modulo(tstep, TS_CORR) == 0 ) then
          !> Compute longitudinal autocorrelation and integral lengthscale
          DOII: do ii = i, ( i + ( x_max_2r/2 ) + 1 ), 1
            IFII: if ( ii <= x_max_2r ) then
              temp_r1i(ii-i) = temp_r1i(ii-i) + u1 * arr_u1_2r(ii,j)
              lint11 = lint11 + ( u1 * arr_u1_2r(ii,j) * dx_laboratoryframe )
            else
              temp_r1i(ii-i) = temp_r1i(ii-i) + u1 * arr_u1_2r(ii-x_width_2r,j)
              lint11 = lint11 + ( u1 * arr_u1_2r(ii-x_width_2r,j) * dx_laboratoryframe )
            end if IFII
          end do DOII
      end if IFCORRTS
      end if IFCORR

    end do DOX
  end do DOY

end subroutine f_FluidstatsRs2D
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsRsCollect 
!> Sum up Reynolds stresses from all processes and compute turbulent kinetic energy.
!! Write Reynolds stresses and turbulent kinetic energy to the real-space statistics
!! file. 
!> @param tstep current timestep
!> @param ierr should return 0
subroutine f_FluidstatsRsCollect ( tstep, ierr )

  use g_parameters, only: MYRANK, MASTER, NODES_X, NODES_Y, NODES_Z, &
                          F_TYPE, F_ISOTROPIC, B110, B220, B330, &
                          TSTEPS, TS_CORR
  use g_domain,     only: dx_laboratoryframe
  use g_files,      only: FUNIT_FSTATR, FUNIT_FCORRL
  use g_constants,  only: ZERO, HALF

  implicit none

  integer,intent(in)               :: tstep 
  PetscErrorCode,intent(inout)            :: ierr

  integer                          :: ii

  !> Determine total maximum velocity for CFL criterion
  call MPI_Allreduce ( umax, ucfl , 1, MPI_Double, MPI_MAX, PETSC_COMM_WORLD, ierr )

  !> Determine total maximum vorticity 
  call MPI_Allreduce ( lomegamax(1), omegamax(1) , 1, MPI_Double, MPI_MAX, PETSC_COMM_WORLD, ierr )
  !> Determine total maximum vorticity 
  call MPI_Allreduce ( lomegamax(2), omegamax(2) , 1, MPI_Double, MPI_MAX, PETSC_COMM_WORLD, ierr )
  !> Determine total maximum vorticity 
  call MPI_Allreduce ( lomegamax(3), omegamax(3) , 1, MPI_Double, MPI_MAX, PETSC_COMM_WORLD, ierr )

  !> Finish velocity statistics
  call MPI_Allreduce ( u1u1, rs_u1u1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2, rs_u2u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3, rs_u3u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2, rs_u1u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u3, rs_u1u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u3, rs_u2u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( skew1, rs_skew1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( skew2, rs_skew2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( skew3, rs_skew3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( flat1, rs_flat1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( flat2, rs_flat2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( flat3, rs_flat3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

  rs_u1u1 = rs_u1u1 / ( real ( NODES_X, kind=C_DOUBLE) * &
                        real ( NODES_Y, kind=C_DOUBLE) * &
                        real ( NODES_Z, kind=C_DOUBLE) ) 
  rs_u2u2 = rs_u2u2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                        real ( NODES_Y, kind=C_DOUBLE) * &
                        real ( NODES_Z, kind=C_DOUBLE) )
  rs_u3u3 = rs_u3u3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                        real ( NODES_Y, kind=C_DOUBLE) * &
                        real ( NODES_Z, kind=C_DOUBLE) )
  rs_u1u2 = rs_u1u2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                        real ( NODES_Y, kind=C_DOUBLE) * &
                        real ( NODES_Z, kind=C_DOUBLE) )
  rs_u1u3 = rs_u1u3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                        real ( NODES_Y, kind=C_DOUBLE) * &
                        real ( NODES_Z, kind=C_DOUBLE) )
  rs_u2u3 = rs_u2u3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                        real ( NODES_Y, kind=C_DOUBLE) * &
                        real ( NODES_Z, kind=C_DOUBLE) )
  rs_skew1 = rs_skew1 / ( real ( NODES_X, kind=C_DOUBLE) * &
                          real ( NODES_Y, kind=C_DOUBLE) * &
                          real ( NODES_Z, kind=C_DOUBLE) )
  rs_skew1 = rs_skew1 / ( rs_u1u1**1.5 )
  rs_skew2 = rs_skew2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                          real ( NODES_Y, kind=C_DOUBLE) * &
                          real ( NODES_Z, kind=C_DOUBLE) )
  rs_skew2 = rs_skew2 / ( rs_u2u2**1.5 )
  rs_skew3 = rs_skew3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                          real ( NODES_Y, kind=C_DOUBLE) * &
                          real ( NODES_Z, kind=C_DOUBLE) )
  rs_skew3 = rs_skew3 / ( rs_u3u3**1.5 )
  rs_flat1 = rs_flat1 / ( real ( NODES_X, kind=C_DOUBLE) * &
                          real ( NODES_Y, kind=C_DOUBLE) * &
                          real ( NODES_Z, kind=C_DOUBLE) )
  rs_flat1 = rs_flat1 / ( rs_u1u1**2 )
  rs_flat2 = rs_flat2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                          real ( NODES_Y, kind=C_DOUBLE) * &
                          real ( NODES_Z, kind=C_DOUBLE) )
  rs_flat2 = rs_flat2 / ( rs_u2u2**2 )
  rs_flat3 = rs_flat3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                          real ( NODES_Y, kind=C_DOUBLE) * &
                          real ( NODES_Z, kind=C_DOUBLE) )
  rs_flat3 = rs_flat3 / ( rs_u3u3**2 )

  !> Finish vorticity statistics
  call MPI_Allreduce ( rmsvort1, rs_rmsvort1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( rmsvort2, rs_rmsvort2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( rmsvort3, rs_rmsvort3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( skewvort1, rs_skewvort1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( skewvort2, rs_skewvort2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( skewvort3, rs_skewvort3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( flatvort1, rs_flatvort1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( flatvort2, rs_flatvort2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( flatvort3, rs_flatvort3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

  rs_rmsvort1 = rs_rmsvort1  
  rs_rmsvort2 = rs_rmsvort2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                real ( NODES_Y, kind=C_DOUBLE) * &
                                real ( NODES_Z, kind=C_DOUBLE) )
  rs_rmsvort3 = rs_rmsvort3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                real ( NODES_Y, kind=C_DOUBLE) * &
                                real ( NODES_Z, kind=C_DOUBLE) )
  rs_skewvort1 = rs_skewvort1 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                  real ( NODES_Y, kind=C_DOUBLE) * &
                                  real ( NODES_Z, kind=C_DOUBLE) )
  rs_skewvort2 = rs_skewvort2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                  real ( NODES_Y, kind=C_DOUBLE) * &
                                  real ( NODES_Z, kind=C_DOUBLE) )
  rs_skewvort3 = rs_skewvort3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                  real ( NODES_Y, kind=C_DOUBLE) * &
                                  real ( NODES_Z, kind=C_DOUBLE) )
  rs_flatvort1 = rs_flatvort1 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                  real ( NODES_Y, kind=C_DOUBLE) * &
                                  real ( NODES_Z, kind=C_DOUBLE) )
  rs_flatvort2 = rs_flatvort2 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                  real ( NODES_Y, kind=C_DOUBLE) * &
                                  real ( NODES_Z, kind=C_DOUBLE) )
  rs_flatvort3 = rs_flatvort3 / ( real ( NODES_X, kind=C_DOUBLE) * &
                                  real ( NODES_Y, kind=C_DOUBLE) * &
                                  real ( NODES_Z, kind=C_DOUBLE) )

  rs_rmsvort1 = rs_rmsvort1**0.5 
  rs_rmsvort2 = rs_rmsvort2**0.5 
  rs_rmsvort3 = rs_rmsvort3**0.5 
  rs_skewvort1 = rs_skewvort1 / ( rs_rmsvort1**3 )
  rs_skewvort2 = rs_skewvort2 / ( rs_rmsvort2**3 )
  rs_skewvort3 = rs_skewvort3 / ( rs_rmsvort3**3 )
  rs_flatvort1 = rs_flatvort1 / ( rs_rmsvort1**4 )
  rs_flatvort2 = rs_flatvort2 / ( rs_rmsvort2**4 )
  rs_flatvort3 = rs_flatvort3 / ( rs_rmsvort3**4 )

  IFCORR: if ( TS_CORR <= TSTEPS ) then
    IFCORRTS: if ( modulo(tstep, TS_CORR) == 0 ) then
      !> Finish longitudinal integral length scale
      call MPI_Allreduce ( lint11, rs_lint11, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
      call MPI_Allreduce ( temp_r1i, r1i, (NODES_X/2 + 1), MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
      r1i = r1i / ( rs_u1u1 * real ( NODES_X, kind=C_DOUBLE) * &
                                     real ( NODES_Y, kind=C_DOUBLE) * &
                                     real ( NODES_Z, kind=C_DOUBLE) )
      rs_lint11 = rs_lint11 / ( rs_u1u1 * real ( NODES_X, kind=C_DOUBLE) * &
                                          real ( NODES_Y, kind=C_DOUBLE) * &
                                          real ( NODES_Z, kind=C_DOUBLE) )
    end if IFCORRTS
  else
    rs_lint11 = ZERO
  end if IFCORR  

  IFMASTER: if ( MYRANK == MASTER ) then
    write(FUNIT_FSTATR,50) tstep, rs_u1u1, rs_u2u2, rs_u3u3, rs_u1u2, rs_u1u3, rs_u2u3, &
                                  rs_skew1, rs_skew2, rs_skew3, rs_flat1, rs_flat2, rs_flat3, &
                                  rs_rmsvort1, rs_rmsvort2, rs_rmsvort3, & 
                                  rs_skewvort1, rs_skewvort2, rs_skewvort3, & 
                                  rs_flatvort1, rs_flatvort2, rs_flatvort3, &
                                  omegamax(1), omegamax(2), omegamax(3), &
                                  rs_lint11
    !> if the correct timestep, write energy and dissipation spectrum to spectra file as a new record
    IFWRITECORR: if ( TS_CORR <= TSTEPS ) then 
      IFWRITECORRTS: if ( modulo(tstep, TS_CORR) == 0 ) then
        write ( FUNIT_FCORRL, 10 )
        write ( FUNIT_FCORRL, 11 ) tstep

        DOWRITENODES: do ii = 0, (NODES_X/2), 1
          write ( FUNIT_FCORRL, 20 ) ( ii * dx_laboratoryframe ), &
                                    r1i(ii)
        end do DOWRITENODES

      end if IFWRITECORRTS
    end if IFWRITECORR

  end if IFMASTER

10  format( '' )
11  format( ' ts: ',i6 )
20  format( 2e14.6 )

50  FORMAT('ts: ',i6,' u1u1: ',e14.6,' u2u2: ',e14.6,' u3u3: ',e14.6, &
           ' u1u2: ',e14.6,' u1u3: ',e14.6,' u2u3: ',e14.6, &
           ' skew1: ', e14.6, ' skew2: ', e14.6, ' skew3: ', e14.6, &
           ' flat1: ', e14.6, ' flat2: ', e14.6, ' flat3: ', e14.6, &
           ' rmsvort1: ',e14.6,' rmsvort2: ',e14.6,' rmsvort3: ',e14.6, &
           ' skewvort1: ', e14.6, ' skewvort2: ', e14.6, ' skewvort3: ', e14.6, &
           ' flatvort1: ', e14.6, ' flatvort2: ', e14.6, ' flatvort3: ', e14.6, & 
           ' vort1max: ', e14.6, ' vort2max: ', e14.6, ' vort3max: ', e14.6 , &
           ' lint11: ', e14.6 )

end subroutine f_FluidstatsRsCollect
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFs2D 
!> Compute 1D spectra locally on 2-D slices.
!> @param ierr should return 0
subroutine f_FluidstatsFs2D ( tstep, ierr )

  use g_constants,  only : ZERO, TWO
  use g_parameters, only : TSTEPS, TS_SPEC1D, NODES_X, F_TYPE, F_ISOTROPIC
  use g_domain,     only : dx_laboratoryframe, bmat, binv, bvalexist
  use f_arrays,     only : arr_u1_1r1w, arr_u2_1r1w, arr_u3_1r1w, &
                           i_min_1r1w, i_max_1r1w, y_min_1r1w, y_max_1r1w
 
  implicit none
 
  integer,intent(in)                                   :: tstep 
  PetscErrorCode,intent(inout)                                :: ierr

  integer                                              :: i, j, m, n
  real(kind=C_DOUBLE),dimension(0:1)                   :: u1, u2, u3
  real(kind=C_DOUBLE),dimension(0:1,1:3)               :: umov, ulab
  real(kind=C_DOUBLE)                                  :: phi11, phi22

  DOJ: do j = y_min_1r1w, y_max_1r1w, 1
    DOI: do i = i_min_1r1w, i_max_1r1w, 1

      !> Get velocity components
      IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
        u1(:) = arr_u1_1r1w(:,i,j)
        u2(:) = arr_u2_1r1w(:,i,j)
      else
        umov(:,1) = arr_u1_1r1w(:,i,j)
        umov(:,2) = arr_u2_1r1w(:,i,j)
        umov(:,3) = arr_u3_1r1w(:,i,j)

        !> @todo this is only for the Bardino/Ferziger type equations
        ulab = ZERO

        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              ulab(:,m) = ulab(:,m) + binv(m,n) * umov(:,n)
            end if IFBEXIST
          end do DON
        end do DOM

        u1(:) = ulab(:,1)
        u2(:) = ulab(:,2)

      end if IFISOTROPIC
! can test conjugate symmetry here - should give identical results with or
! without if condition, as the velocity components must be real on the i=0
! plane
!      IFZERO: if ( i /= 0 ) then              
        phi11 = u1(0)*u1(0) + u1(1)*u1(1)
        phi22 = u2(0)*u2(0) + u2(1)*u2(1)
!      else
!        phi11 = u1(0)*u1(0) 
!        phi22 = u2(0)*u2(0)
!      end if IFZERO

      !> Normalise energy on y-z plane by length in x direction
      IFNOTISOTROPIC: if ( F_TYPE /= F_ISOTROPIC ) then
        phi11 = phi11 * ( binv(1,1)**3 )
        phi22 = phi22 * ( binv(1,1)**3 )
      end if IFNOTISOTROPIC

      temp_spec_1i(i) = temp_spec_1i(i) + ( TWO * phi11 )
      temp_spec_2i(i) = temp_spec_2i(i) + ( TWO * phi22 )
    end do DOI
  end do DOJ

end subroutine f_FluidstatsFs2D
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsCollect 
!> Sum up and normalise 1D spectra and write them to file.
!> @param tstep current timestep
!> @param ierr should return 0
subroutine f_FluidstatsFsCollect ( tstep, ierr )

  use g_parameters, only : MYRANK, MASTER, NODES_Y, NODES_Z, TS_SPEC1D, TSTEPS, &
                           F_TYPE, F_ISOTROPIC, B110
  use g_domain,     only : NODES_KI, KI_MAX
  use g_files,      only : FUNIT_1DSPEC
  use g_constants,  only : ZERO, HALF

  implicit none

  integer,intent(in)               :: tstep 
  PetscErrorCode,intent(inout)            :: ierr

  integer                          :: i

  !> sum contributions from all processes 
  call MPI_Allreduce ( temp_spec_1i, spec_1i, NODES_KI, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_spec_2i, spec_2i, NODES_KI, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

!  !> Divide by number of 1D spectra to obtain average
!  spec_1i = spec_1i / ( real(NODES_Y,kind=C_DOUBLE) * &
!                        real(NODES_Z,kind=C_DOUBLE) )
!  spec_2i = spec_2i / ( real(NODES_Y,kind=C_DOUBLE) * &
!                        real(NODES_Z,kind=C_DOUBLE) )

  !> Write spectrum if required
  IFSPEC1D: if ( TS_SPEC1D <= TSTEPS ) then
    IFSPEC1DTS: if ( modulo(tstep, TS_SPEC1D) == 0 ) then
      IFMASTER: if ( MYRANK == MASTER ) then

        write ( FUNIT_1DSPEC, 10 )
        write ( FUNIT_1DSPEC, 11 ) tstep

        DOWRITENODES1D: do i = 0, KI_MAX, 1
          write ( FUNIT_1DSPEC, 21 ) ( real(i,kind=C_DOUBLE) / real(B110,kind=C_DOUBLE) ), &
                                     spec_1i(i) , &
                                     spec_2i(i)
        end do DOWRITENODES1D

      end if IFMASTER
    end if IFSPEC1DTS
  end if IFSPEC1D

  !> Derive further statistics
  IFNOTISOTROPIC: if ( F_TYPE /= F_ISOTROPIC ) then
    call f_FluidstatsHomogeneous ( tstep, ierr )
  end if IFNOTISOTROPIC

  10  format( '' )
  11  format( ' ts: ',i6 )
  20  format( 7e14.6 )
  21  format( 3e14.6 )

end subroutine f_FluidstatsFsCollect
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFs 
!> Get access to velocity arrays and call subroutines to compute wave-space 
!! statistics.
!> @param tstep current time step
!> @param time current simulation time
!> @param ierr should return 0
subroutine f_FluidstatsFs ( tstep, time, ierr )

  use f_arrays,     only : da3w, i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           u1_3w, u2_3w, u3_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w
  use g_parameters, only : F_TYPE, F_ISOTROPIC, &
                           NS_DEFORMED, NS_BF_ROT, NS_R_CONV, NS_R_ROT
 
  implicit none
 
  integer,intent(in)             :: tstep
  real(kind=C_DOUBLE),intent(in) :: time
  PetscErrorCode,intent(inout)          :: ierr

  PetscErrorCode        :: perr

!> Get read access to all velocity components
!  call DMDAVecGetArrayReadF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!  call DMDAVecGetArrayReadF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!  call DMDAVecGetArrayReadF90 ( da3w, u3_3w, arr_u3_3w, ierr )
  call DMDAVecGetArrayReadF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecGetArrayReadF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecGetArrayReadF90 ( da3w, u3_3w, arr_u3_3w, perr )

  call f_FluidstatsFsZeroQuants ( ierr )
  !> If the domain is isotropic, statistics can be computed directly.
  !! Otherwise a transformation into the laboratory system is needed.
  IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
    call f_FluidstatsFsPrimaryQuants ( tstep, ierr )
  else
    select case ( NS_DEFORMED )
    !-------------------------------
      case ( NS_BF_ROT ) 
        !> Count number of modes per bin, as this changes at each
        !! time step
        call f_FluidstatsModeDensityBFRot ( ierr )
        !> Compute primary quantities in the laboratory system
        call f_FluidstatsFsPrimaryQuantsBFRot ( tstep, ierr )
      case ( NS_R_ROT ) 
        !> Count number of modes per bin, as this changes at each
        !! time step
        call f_FluidstatsModeDensityRRot ( ierr )
        !> Compute primary quantities in the laboratory system
        call f_FluidstatsFsPrimaryQuantsRRot ( tstep, ierr )
      case ( NS_R_CONV ) 
        !> Count number of modes per bin, as this changes at each
        !! time step
        call f_FluidstatsModeDensityRRot ( ierr )
        !> Compute primary quantities in the laboratory system
        call f_FluidstatsFsPrimaryQuantsRRot ( tstep, ierr )
    end select
  end if IFISOTROPIC
  call f_FluidstatsFsSecondaryQuants ( ierr )
!  IFOTHER: if ( F_TYPE /= F_ISOTROPIC ) then
!    call f_FluidstatsHomogeneous ( tstep, ierr )
!  end if IFOTHER
  call f_FluidstatsFsWriteQuants ( tstep, time, ierr )
!  if ( tstep > 0 ) call f_FluidstatsAnalyticalViscous ( tstep, time, ierr )

!> Return velocity components to PETSc
!  call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!  call DMDAVecRestoreArrayReadF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!  call DMDAVecRestoreArrayReadF90 ( da3w, u3_3w, arr_u3_3w, ierr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, perr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u2_3w, arr_u2_3w, perr )
  call DMDAVecRestoreArrayReadF90 ( da3w, u3_3w, arr_u3_3w, perr )

end subroutine f_FluidstatsFs
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsZeroQuants 
!> Set wave-space quantities to zero.
!> @param ierr should return 0
subroutine f_FluidstatsFsZeroQuants ( ierr )

  use g_constants, only: ZERO
 
  implicit none
 
  PetscErrorCode,intent(inout)            :: ierr

  !> Initialise Reynolds stresses
  u1u1 = ZERO
  u2u2 = ZERO
  u3u3 = ZERO
  u1u2 = ZERO
  u1u3 = ZERO
  u2u3 = ZERO

  !> Initialise primary quantities
  f_ken = ZERO
  fn_ken = ZERO
  fn_ken1 = ZERO
  fn_ken2 = ZERO
  fn_ken3 = ZERO
  f_dis = ZERO
  fn_dis = ZERO
  f_disskew = ZERO
  fn_disskew = ZERO
  temp_ken = ZERO
  temp_ken1 = ZERO
  temp_ken2 = ZERO
  temp_ken3 = ZERO
  temp_spec_1i = ZERO
  temp_spec_2i = ZERO
  spec_1i = ZERO
  temp_spec_2i = ZERO
  spec_2i = ZERO
!  temp_ken12 = ZERO
!  temp_ken22 = ZERO
!  temp_ken32 = ZERO
!  temp_ken13 = ZERO
!  temp_ken23 = ZERO
!  temp_ken33 = ZERO
  temp_dis = ZERO
  temp_disskew = ZERO

end subroutine f_FluidstatsFsZeroQuants
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsPrimaryQuants 
!> Compute energy and dissipation spectra and get turbulent kinetic energy,
!! dissipation and dissipation skew from these.
!> @param tstep current time step
!> @param ierr should return 0
subroutine f_FluidstatsFsPrimaryQuants ( tstep, ierr )

  use g_constants,  only : ZERO, HALF, ONE, TWO, FOUR
  use g_parameters, only : MYRANK, MASTER, F_NU
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K
  use f_arrays,     only : da3w, i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           f_ArraysWavenumber
 
  implicit none
 
  integer,intent(in)                  :: tstep
  PetscErrorCode,intent(inout)        :: ierr

  real(kind=C_DOUBLE),dimension(0:1)  :: u1, u2, u3
  real(kind=C_DOUBLE)                 :: phi11, phi22, phi33, phi12, phi13, phi23
  integer                             :: i, j, k, kk
  real(kind=C_DOUBLE)                 :: bin_lower, bin_upper

  real(kind=C_DOUBLE)                 :: ek, dk, dsk
  real(kind=C_DOUBLE)                 :: wavenumber_i, wavenumber_j, wavenumber_k, k_mag

  !> Calculate energy spectrum on each node
  DOI: do i = i_min_3w, i_max_3w, 1
  wavenumber_i = f_ArraysWavenumber ( i, KI_MAX, K0I )
    DOJ: do j = j_min_3w, j_max_3w, 1
    wavenumber_j = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      DOK: do k = k_min_3w, k_max_3w, 1
      wavenumber_k = f_ArraysWavenumber ( k, KK_MAX, K0K )

        !> Compute magnitude of wave vector
        k_mag = sqrt ( (wavenumber_i*wavenumber_i) + &
                       (wavenumber_j*wavenumber_j) + &
                       (wavenumber_k*wavenumber_k) )

        u1(:) = arr_u1_3w(:,k,j,i)
        u2(:) = arr_u2_3w(:,k,j,i)
        u3(:) = arr_u3_3w(:,k,j,i)

        !> multiplication of two complex numbers expressed as real numbers
        !> u1*conj(u2) = real(u1)*real(u2) + imag(u1)*imag(u2); rest disappears due to conjugate symmetry
        phi11 = u1(0) * u1(0) + u1(1) * u1(1) 
        phi22 = u2(0) * u2(0) + u2(1) * u2(1) 
        phi33 = u3(0) * u3(0) + u3(1) * u3(1) 
        phi12 = u1(0) * u2(0) + u1(1) * u2(1) 
        phi13 = u1(0) * u3(0) + u1(1) * u3(1) 
        phi23 = u2(0) * u3(0) + u2(1) * u3(1) 

        !> double values due to conjugate symmetry
        IFNOTIZERO: if (i /= 0) then
          phi11 = phi11 * TWO
          phi22 = phi22 * TWO
          phi33 = phi33 * TWO
          phi12 = phi12 * TWO
          phi13 = phi13 * TWO
          phi23 = phi23 * TWO

        end if IFNOTIZERO

        !> calculate components of Reynolds stress tensor
          u1u1 = u1u1 + phi11
          u2u2 = u2u2 + phi22
          u3u3 = u3u3 + phi33
          u1u2 = u1u2 + phi12
          u1u3 = u1u3 + phi13
          u2u3 = u2u3 + phi23

          ! calculate energy spectrum (as defined in Tennekes & Lumley, 1972)
          ! see also Pope (6.188)
          ek = HALF * (phi11 + phi22 + phi33)

          ! calculate dissipation spectrum (Turbulent Flows, Pope, 2001, 6.191)
          dk = TWO * k_mag * k_mag * F_NU * ek
                
          dsk = FOUR * k_mag * k_mag * k_mag * k_mag * F_NU * F_NU * ek

                ! find appropriate bin for current mode
          DOBIN: do kk = 1, kkmax, 1

            ! find bin limits
            bin_upper = real(kk,kind=C_DOUBLE) + ( k_min * HALF )
            bin_lower = real(kk,kind=C_DOUBLE) - ( k_min * HALF )

            ! check to see if k is in current bin
            IFLOWER: if ( k_mag >= bin_lower ) then
              IFUPPER: if ( k_mag < bin_upper ) then

                temp_ken(kk) = temp_ken(kk) + ek ! add energy to appropriate bin
                temp_ken1(kk) = temp_ken1(kk) + ( HALF * phi11 )
                temp_ken2(kk) = temp_ken2(kk) + ( HALF * phi22 )
                temp_ken3(kk) = temp_ken3(kk) + ( HALF * phi33 )
                temp_dis(kk) = temp_dis(kk) + dk ! add dissipation to appropriate bin
                temp_disskew(kk) = temp_disskew(kk) + dsk

              end if IFUPPER
            end if IFLOWER

          end do DOBIN

      end do DOK
    end do DOJ
  end do DOI


    !> sum contributions from all processes
  call MPI_Allreduce ( temp_ken, spec_ken, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken1, spec_ken1, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken2, spec_ken2, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken3, spec_ken3, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_dis, spec_dis, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_disskew, spec_disskew, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

    !> adjust for wave vector density and calculate normalised 
    !! turbulent kinetic energy and dissipation
  DOPRIMARYSTATS: do kk = 1, kkmax, 1
    IFNOTEMPTY: if ( dens_mod(kk) /= ZERO ) then

      f_ken           = f_ken       + spec_ken(kk) ! un-normalised total kinetic energy
      f_dis           = f_dis       + spec_dis(kk) ! un-normalised total dissipation
      f_disskew       = f_disskew   + spec_disskew(kk) ! un-normalised total dissipation

      spec_ken(kk)     = spec_ken(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken1(kk)    = spec_ken1(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken2(kk)    = spec_ken2(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken3(kk)    = spec_ken3(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_dis(kk)     = spec_dis(kk) / dens_mod(kk) ! normalised dissipation spectrum
      spec_disskew(kk) = spec_disskew(kk) / dens_mod(kk)

      fn_ken          = fn_ken      + spec_ken(kk) ! normalised total kinetic energy
      fn_ken1         = fn_ken1     + spec_ken1(kk) ! normalised total kinetic energy
      fn_ken2         = fn_ken2     + spec_ken2(kk) ! normalised total kinetic energy
      fn_ken3         = fn_ken3     + spec_ken3(kk) ! normalised total kinetic energy
      fn_dis          = fn_dis      + spec_dis(kk) ! normalised total dissipation rate
      fn_disskew      = fn_disskew  + spec_disskew(kk)

    end if IFNOTEMPTY
  end do DOPRIMARYSTATS

!    if(myrank.eq.master)then
!    end if

  !> Determine Reynolds stresses
  call MPI_Allreduce ( u1u1, fs_u1u1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2, fs_u2u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3, fs_u3u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2, fs_u1u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u3, fs_u1u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u3, fs_u2u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

  !> Compute turbulent kinetic energy
  tke = HALF * ( fs_u1u1 + fs_u2u2 + fs_u3u3 )

end subroutine f_FluidstatsFsPrimaryQuants
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsPrimaryQuantsBFRot 
!> Compute energy and dissipation spectra and get turbulent kinetic energy,
!! dissipation and dissipation skew from these.
!> @param tstep current time step
!> @param ierr should return 0
subroutine f_FluidstatsFsPrimaryQuantsBFRot ( tstep, ierr )

  use g_constants,  only : ZERO, HALF, ONE, TWO, FOUR
  use g_parameters, only : MYRANK, MASTER, F_NU, B110, B220, B330
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, binv, bvalexist
  use f_arrays,     only : da3w, i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           f_ArraysWavenumber
 
  implicit none
 
  integer,intent(in)                     :: tstep
  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_DOUBLE)                    :: phi11, phi22, phi33, phi12, phi13, phi23
  integer                                :: i, j, k, kk
  integer                                :: m,n
  real(kind=C_DOUBLE)                    :: bin_lower, bin_upper

  real(kind=C_DOUBLE)                    :: ek, dk, dsk
  real(kind=C_DOUBLE),dimension(1:3)     :: wvnmov, wvnlab
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: umov, ulab
  real(kind=C_DOUBLE)                    :: k_mag

  !> Calculate energy spectrum on each node
  DOI: do i = i_min_3w, i_max_3w, 1
  wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
    DOJ: do j = j_min_3w, j_max_3w, 1
    wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      DOK: do k = k_min_3w, k_max_3w, 1
      wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

        !> Compute wavenumber and velocity in the laboratory system.
        umov(:,1) = arr_u1_3w(:,k,j,i)
        umov(:,2) = arr_u2_3w(:,k,j,i)
        umov(:,3) = arr_u3_3w(:,k,j,i)

        wvnlab = ZERO
        ulab = ZERO

        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
              ulab(:,m) = ulab(:,m) + binv(m,n) * umov(:,n)
            end if IFBEXIST
          end do DON
        end do DOM

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        !> multiplication of two complex numbers expressed as real numbers
        !> u1*conj(u2) = real(u1)*real(u2) + imag(u1)*imag(u2); rest disappears due to conjugate symmetry
        phi11 = ulab(0,1) * ulab(0,1) + ulab(1,1) * ulab(1,1) 
        phi22 = ulab(0,2) * ulab(0,2) + ulab(1,2) * ulab(1,2) 
        phi33 = ulab(0,3) * ulab(0,3) + ulab(1,3) * ulab(1,3) 
        phi12 = ulab(0,1) * ulab(0,2) + ulab(1,1) * ulab(1,2) 
        phi13 = ulab(0,1) * ulab(0,3) + ulab(1,1) * ulab(1,3) 
        phi23 = ulab(0,2) * ulab(0,3) + ulab(1,2) * ulab(1,3) 

        !> double values due to conjugate symmetry
        IFNOTIZERO: if (i /= 0) then
          phi11 = phi11 * TWO
          phi22 = phi22 * TWO
          phi33 = phi33 * TWO
          phi12 = phi12 * TWO
          phi13 = phi13 * TWO
          phi23 = phi23 * TWO

        end if IFNOTIZERO

        !> calculate components of Reynolds stress tensor
          u1u1 = u1u1 + phi11
          u2u2 = u2u2 + phi22
          u3u3 = u3u3 + phi33
          u1u2 = u1u2 + phi12
          u1u3 = u1u3 + phi13
          u2u3 = u2u3 + phi23

          ! calculate energy spectrum (as defined in Tennekes & Lumley, 1972)
          ek = HALF * (phi11 + phi22 + phi33)

          ! calculate dissipation spectrum (Turbulent Flows, Pope, 2001, 6.191)
          dk = TWO * k_mag * k_mag * F_NU * ek
                
          dsk = FOUR * k_mag * k_mag * k_mag * k_mag * F_NU * F_NU * ek

                ! find appropriate bin for current mode
          DOBIN: do kk = 1, kkmax, 1

            ! find bin limits
            bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
            bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )

            ! check to see if k is in current bin
            IFLOWER: if ( k_mag >= bin_lower ) then
              IFUPPER: if ( k_mag < bin_upper ) then

                temp_ken(kk) = temp_ken(kk) + ek ! add energy to appropriate bin
                temp_ken1(kk) = temp_ken1(kk) + ( HALF * phi11 )
                temp_ken2(kk) = temp_ken2(kk) + ( HALF * phi22 )
                temp_ken3(kk) = temp_ken3(kk) + ( HALF * phi33 )
                temp_dis(kk) = temp_dis(kk) + dk ! add dissipation to appropriate bin
                temp_disskew(kk) = temp_disskew(kk) + dsk

              end if IFUPPER
            end if IFLOWER

          end do DOBIN

      end do DOK
    end do DOJ
  end do DOI

!  write(*,*) 'wvn', wvnmov, wvnlab, 'bmat', bmat, 'binv', binv, 'umov', umov, 'ulab', ulab!,'ken', temp_ken, 'dis', temp_dis, 'disskew', temp_disskew

  !> Correct for domain length
!  temp_ken     = temp_ken * bmat(1,1) * bmat(2,2) * bmat(3,3)
!  temp_ken1    = temp_ken1 * bmat(1,1) * bmat(2,2) * bmat(3,3)
!  temp_ken2    = temp_ken2 * bmat(1,1) * bmat(2,2) * bmat(3,3)
!  temp_ken3    = temp_ken3 * bmat(1,1) * bmat(2,2) * bmat(3,3)
!  temp_dis     = temp_dis * bmat(1,1) * bmat(2,2) * bmat(3,3)
!  temp_disskew = temp_disskew * bmat(1,1) * bmat(2,2) * bmat(3,3)

    !> sum contributions from all processes
  call MPI_Allreduce ( temp_ken, spec_ken, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken1, spec_ken1, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken2, spec_ken2, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken3, spec_ken3, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_dis, spec_dis, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_disskew, spec_disskew, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

    !> adjust for wave vector density and calculate normalised 
    !! turbulent kinetic energy and dissipation
  DOPRIMARYSTATS: do kk = 1, kkmax, 1
    IFNOTEMPTY: if ( dens_mod(kk) /= ZERO ) then

      f_ken           = f_ken       + spec_ken(kk) ! un-normalised total kinetic energy
      f_dis           = f_dis       + spec_dis(kk) ! un-normalised total dissipation
      f_disskew       = f_disskew   + spec_disskew(kk) ! un-normalised total dissipation

      spec_ken(kk)     = spec_ken(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken1(kk)    = spec_ken1(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken2(kk)    = spec_ken2(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken3(kk)    = spec_ken3(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_dis(kk)     = spec_dis(kk) / dens_mod(kk) ! normalised dissipation spectrum
      spec_disskew(kk) = spec_disskew(kk) / dens_mod(kk)

      fn_ken          = fn_ken      + spec_ken(kk) ! normalised total kinetic energy
      fn_ken1         = fn_ken1     + spec_ken1(kk) ! normalised total kinetic energy
      fn_ken2         = fn_ken2     + spec_ken2(kk) ! normalised total kinetic energy
      fn_ken3         = fn_ken3     + spec_ken3(kk) ! normalised total kinetic energy
      fn_dis          = fn_dis      + spec_dis(kk) ! normalised total dissipation rate
      fn_disskew      = fn_disskew  + spec_disskew(kk)

    end if IFNOTEMPTY
  end do DOPRIMARYSTATS

!  write(*,*) 'fn_ken, 1, 2, 3: ', fn_ken, fn_ken1, fn_ken2, fn_ken3
!    if(myrank.eq.master)then
!    end if

  !> Determine Reynolds stresses
  call MPI_Allreduce ( u1u1, fs_u1u1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2, fs_u2u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3, fs_u3u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2, fs_u1u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u3, fs_u1u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u3, fs_u2u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

!  !> Correct for domain length
!  fs_u1u1 = fs_u1u1 * bmat(1,1) * bmat(2,2) * bmat(3,3) 
!  fs_u2u2 = fs_u2u2 * bmat(1,1) * bmat(2,2) * bmat(3,3) 
!  fs_u3u3 = fs_u3u3 * bmat(1,1) * bmat(2,2) * bmat(3,3) 
!  fs_u1u2 = fs_u1u2 * bmat(1,1) * bmat(2,2) * bmat(3,3) 
!  fs_u1u3 = fs_u1u3 * bmat(1,1) * bmat(2,2) * bmat(3,3) 
!  fs_u2u3 = fs_u2u3 * bmat(1,1) * bmat(2,2) * bmat(3,3) 

!  !> Correct Reynolds stresses for size of real domain vs computational domain
!  fs_u1u1 = fs_u1u1 * ( real(B110,kind=C_DOUBLE) &
!             * real(B220,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) )
!  fs_u2u2 = fs_u2u2 * ( real(B110,kind=C_DOUBLE) &
!             * real(B220,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) )
!  fs_u3u3 = fs_u3u3 * ( real(B110,kind=C_DOUBLE) &
!             * real(B220,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) )
!  fs_u1u2 = fs_u1u2 * ( real(B110,kind=C_DOUBLE) &
!             * real(B220,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) )
!  fs_u1u3 = fs_u1u3 * ( real(B110,kind=C_DOUBLE) &
!             * real(B220,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) )
!  fs_u2u3 = fs_u2u3 * ( real(B110,kind=C_DOUBLE) &
!             * real(B220,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) )

  !> Compute turbulent kinetic energy
  tke = HALF * ( fs_u1u1 + fs_u2u2 + fs_u3u3 )

end subroutine f_FluidstatsFsPrimaryQuantsBFRot
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsPrimaryQuantsRRot 
!> Compute energy and dissipation spectra and get turbulent kinetic energy,
!! dissipation and dissipation skew from these.
!> @param tstep current time step
!> @param ierr should return 0
subroutine f_FluidstatsFsPrimaryQuantsRRot ( tstep, ierr )

  use g_constants,  only : ZERO, HALF, ONE, TWO, FOUR
  use g_parameters, only : MYRANK, MASTER, F_NU
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, binv, bvalexist
  use f_arrays,     only : da3w, i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           f_ArraysWavenumber
 
  implicit none
 
  integer,intent(in)                     :: tstep
  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_DOUBLE)                    :: phi11, phi22, phi33, phi12, phi13, phi23
  integer                                :: i, j, k, kk
  integer                                :: m,n
  real(kind=C_DOUBLE)                    :: bin_lower, bin_upper

  real(kind=C_DOUBLE)                    :: ek, dk, dsk
  real(kind=C_DOUBLE),dimension(1:3)     :: wvnmov, wvnlab
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: umov, ulab
  real(kind=C_DOUBLE)                    :: k_mag

  !> Calculate energy spectrum on each node
  DOI: do i = i_min_3w, i_max_3w, 1
  wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
    DOJ: do j = j_min_3w, j_max_3w, 1
    wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
      DOK: do k = k_min_3w, k_max_3w, 1
      wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

        !> Compute wavenumber and velocity in the laboratory system.
        umov(:,1) = arr_u1_3w(:,k,j,i)
        umov(:,2) = arr_u2_3w(:,k,j,i)
        umov(:,3) = arr_u3_3w(:,k,j,i)

        wvnlab = ZERO
        ulab = umov

        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXIST
          end do DON
        end do DOM

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        !> multiplication of two complex numbers expressed as real numbers
        !> u1*conj(u2) = real(u1)*real(u2) + imag(u1)*imag(u2); rest disappears due to conjugate symmetry
        phi11 = ulab(0,1) * ulab(0,1) + ulab(1,1) * ulab(1,1) 
        phi22 = ulab(0,2) * ulab(0,2) + ulab(1,2) * ulab(1,2) 
        phi33 = ulab(0,3) * ulab(0,3) + ulab(1,3) * ulab(1,3) 
        phi12 = ulab(0,1) * ulab(0,2) + ulab(1,1) * ulab(1,2) 
        phi13 = ulab(0,1) * ulab(0,3) + ulab(1,1) * ulab(1,3) 
        phi23 = ulab(0,2) * ulab(0,3) + ulab(1,2) * ulab(1,3) 

        !> double values due to conjugate symmetry
        IFNOTIZERO: if (i /= 0) then
          phi11 = phi11 * TWO
          phi22 = phi22 * TWO
          phi33 = phi33 * TWO
          phi12 = phi12 * TWO
          phi13 = phi13 * TWO
          phi23 = phi23 * TWO

        end if IFNOTIZERO

        !> calculate components of Reynolds stress tensor
          u1u1 = u1u1 + phi11
          u2u2 = u2u2 + phi22
          u3u3 = u3u3 + phi33
          u1u2 = u1u2 + phi12
          u1u3 = u1u3 + phi13
          u2u3 = u2u3 + phi23

          ! calculate energy spectrum (as defined in Tennekes & Lumley, 1972)
          ek = HALF * (phi11 + phi22 + phi33)

          ! calculate dissipation spectrum (Turbulent Flows, Pope, 2001, 6.191)
          dk = TWO * k_mag * k_mag * F_NU * ek
                
          dsk = FOUR * k_mag * k_mag * k_mag * k_mag * F_NU * F_NU * ek

                ! find appropriate bin for current mode
          DOBIN: do kk = 1, kkmax, 1

            ! find bin limits
            bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
            bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )

            ! check to see if k is in current bin
            IFLOWER: if ( k_mag >= bin_lower ) then
              IFUPPER: if ( k_mag < bin_upper ) then

                temp_ken(kk) = temp_ken(kk) + ek ! add energy to appropriate bin
                temp_ken1(kk) = temp_ken1(kk) + ( HALF * phi11 )
                temp_ken2(kk) = temp_ken2(kk) + ( HALF * phi22 )
                temp_ken3(kk) = temp_ken3(kk) + ( HALF * phi33 )
                temp_dis(kk) = temp_dis(kk) + dk ! add dissipation to appropriate bin
                temp_disskew(kk) = temp_disskew(kk) + dsk

              end if IFUPPER
            end if IFLOWER

          end do DOBIN

      end do DOK
    end do DOJ
  end do DOI

!  write(*,*) 'wvn', wvnmov, wvnlab, 'bmat', bmat, 'binv', binv, 'umov', umov, 'ulab', ulab!,'ken', temp_ken, 'dis', temp_dis, 'disskew', temp_disskew

    !> sum contributions from all processes
  call MPI_Allreduce ( temp_ken, spec_ken, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken1, spec_ken1, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken2, spec_ken2, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_ken3, spec_ken3, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_dis, spec_dis, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( temp_disskew, spec_disskew, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

    !> adjust for wave vector density and calculate normalised 
    !! turbulent kinetic energy and dissipation
  DOPRIMARYSTATS: do kk = 1, kkmax, 1
    IFNOTEMPTY: if ( dens_mod(kk) /= ZERO ) then

      f_ken           = f_ken       + spec_ken(kk) ! un-normalised total kinetic energy
      f_dis           = f_dis       + spec_dis(kk) ! un-normalised total dissipation
      f_disskew       = f_disskew   + spec_disskew(kk) ! un-normalised total dissipation

      spec_ken(kk)     = spec_ken(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken1(kk)    = spec_ken1(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken2(kk)    = spec_ken2(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_ken3(kk)    = spec_ken3(kk) / dens_mod(kk) ! normalised energy spectrum
      spec_dis(kk)     = spec_dis(kk) / dens_mod(kk) ! normalised dissipation spectrum
      spec_disskew(kk) = spec_disskew(kk) / dens_mod(kk)

      fn_ken          = fn_ken      + spec_ken(kk) ! normalised total kinetic energy
      fn_ken1         = fn_ken1     + spec_ken1(kk) ! normalised total kinetic energy
      fn_ken2         = fn_ken2     + spec_ken2(kk) ! normalised total kinetic energy
      fn_ken3         = fn_ken3     + spec_ken3(kk) ! normalised total kinetic energy
      fn_dis          = fn_dis      + spec_dis(kk) ! normalised total dissipation rate
      fn_disskew      = fn_disskew  + spec_disskew(kk)

    end if IFNOTEMPTY
  end do DOPRIMARYSTATS

!  write(*,*) 'fn_ken, 1, 2, 3: ', fn_ken, fn_ken1, fn_ken2, fn_ken3
!    if(myrank.eq.master)then
!    end if

  !> Determine Reynolds stresses
  call MPI_Allreduce ( u1u1, fs_u1u1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2, fs_u2u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3, fs_u3u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2, fs_u1u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u3, fs_u1u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u3, fs_u2u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

  !> Compute turbulent kinetic energy
  tke = HALF * ( fs_u1u1 + fs_u2u2 + fs_u3u3 )

end subroutine f_FluidstatsFsPrimaryQuantsRRot
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine f_FluidstatsModeDensityBFRot 
!> Find out how many wave modes are in each bin.
!> @param ierr should return 0
subroutine f_FluidstatsModeDensityBFRot ( ierr )

!    use d_arrays,    only : nodes_x,nodes_y,nodes_z 
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, binv, bvalexist
  use g_parameters, only : B110, B220, B330
  use g_constants,  only : ZERO, HALF, TWO, THREE, FOUR, ONEPI
  use f_arrays,     only : i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           f_ArraysWavenumber
!    use f_arrays,    only : fluid_k, ki_min, kj_min, kk_min, &
!                                     ki_max, kj_max, kk_max

  implicit none

  PetscErrorCode,intent(inout)               :: ierr
  integer                             :: alloc_stat
  integer                             :: i,j,k,kk
  integer                             :: m,n
  integer                             :: num_modes
  integer                             :: num_temp
  real(kind=C_DOUBLE),dimension(1:3)  :: wvnmov, wvnlab
  real(kind=C_DOUBLE)                 :: k_mag
  real(kind=C_DOUBLE)                 :: bin_upper
  real(kind=C_DOUBLE)                 :: bin_lower

  ! find appropriate bin for current mode
  DOBIN: do kk = 1, kkmax, 1

    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    ! find bin limits
    bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
    bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )

  !> 3D mode density

    ! find number of modes in each bin
    DOI: do i = i_min_3w, i_max_3w, 1
    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOJ: do j = j_min_3w, j_max_3w, 1
      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
        DOK: do k = k_min_3w, k_max_3w, 1
        wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

          !> Compute wavenumber in the laboratory system.
          wvnlab = ZERO

          DOM: do m = 1, 3, 1
            DON: do n = 1, 3, 1
              IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
              end if IFBEXIST
            end do DON
          end do DOM

          !> Compute magnitude of wave vector
          k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                         ( wvnlab(2) * wvnlab(2) ) + &
                         ( wvnlab(3) * wvnlab(3) ) )

            ! check to see if k is in current bin
            IFLOWER: if ( k_mag >= bin_lower ) then
              IFUPPER: if ( k_mag < bin_upper ) then

                !> Need to count twice for conjugate symmetry except i=0
                IFI0: if ( i /= 0 ) then
                  num_temp = num_temp + 2
                else
                  num_temp = num_temp + 1
                end if IFI0

              end if IFUPPER
            end if IFLOWER

          end do DOK
        end do DOJ
      end do DOI

    ! sum counters from all processes
    call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                   ( FOUR / THREE * ONEPI *( bin_upper**3 - bin_lower**3 ) )

  
    !> Correct for size of real domain vs computational domain
    dens_mod(kk) = dens_mod(kk) / ( real(B110,kind=C_DOUBLE) &
             * real(B220,kind=C_DOUBLE) &
             * real(B330,kind=C_DOUBLE) )

  !> 1D mode density i direction
    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    IFI0MODE: if ( i_min_3w <= 0 ) then
      ! find number of modes in each bin
      wvnmov(1) = f_ArraysWavenumber ( 0, KI_MAX, K0I )

      DOJI: do j = j_min_3w, j_max_3w, 1
      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
        DOKI: do k = k_min_3w, k_max_3w, 1
        wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

          !> Compute wavenumber in the laboratory system.
          wvnlab = ZERO

          DOMI: do m = 1, 3, 1
            DONI: do n = 1, 3, 1
              IFBEXISTI: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
              end if IFBEXISTI
            end do DONI
          end do DOMI

          !> Compute magnitude of wave vector
          k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                         ( wvnlab(2) * wvnlab(2) ) + &
                         ( wvnlab(3) * wvnlab(3) ) )

          ! check to see if k is in current bin
          IFLOWERI: if ( k_mag >= bin_lower ) then
            IFUPPERI: if ( k_mag < bin_upper ) then

              num_temp = num_temp + 1

            end if IFUPPERI
          end if IFLOWERI

        end do DOKI
      end do DOJI
    end if IFI0MODE

    ! sum counters from all processes
    call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod_i(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                     ( ONEPI *( bin_upper**2 - bin_lower**2 ) )

    !> Correct for size of real domain vs computational domain
    dens_mod_i(kk) = dens_mod(kk) / ( real(B110,kind=C_DOUBLE) &
                                    * real(B220,kind=C_DOUBLE) &
                                    * real(B330,kind=C_DOUBLE) )

!    write(*,*) ' kk, dens_mod_i(kk): ', kk, dens_mod_i(kk) 



  !> 1D mode density j direction
    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

   ! find number of modes in each bin
     wvnmov(2) = f_ArraysWavenumber ( 0, KJ_MAX, K0J )
    DOIJ: do i = i_min_3w, i_max_3w, 1
    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOKJ: do k = k_min_3w, k_max_3w, 1
      wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

        !> Compute wavenumber in the laboratory system.
        wvnlab = ZERO

        DOMJ: do m = 1, 3, 1
          DONJ: do n = 1, 3, 1
            IFBEXISTJ: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXISTJ
          end do DONJ
        end do DOMJ

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        ! check to see if k is in current bin
        IFLOWERJ: if ( k_mag >= bin_lower ) then
          IFUPPERJ: if ( k_mag < bin_upper ) then

            num_temp = num_temp + 1

          end if IFUPPERJ
        end if IFLOWERJ

      end do DOKJ
    end do DOIJ

   ! sum counters from all processes
   call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

   ! calculate mode density
   dens_mod_j(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                    ( ONEPI *( bin_upper**2 - bin_lower**2 ) )




!  !> 1D mode density k direction
    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    ! find number of modes in each bin
    wvnmov(3) = f_ArraysWavenumber ( 0, KK_MAX, K0K )
    DOIK: do i = i_min_3w, i_max_3w, 1
    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOJK: do j = j_min_3w, j_max_3w, 1
      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )

        !> Compute wavenumber in the laboratory system.
        wvnlab = ZERO

        DOMK: do m = 1, 3, 1
          DONK: do n = 1, 3, 1
            IFBEXISTK: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXISTK
          end do DONK
        end do DOMK

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        ! check to see if k is in current bin
        IFLOWERK: if ( k_mag >= bin_lower ) then
          IFUPPERK: if ( k_mag < bin_upper ) then

            num_temp = num_temp + 1

          end if IFUPPERK
        end if IFLOWERK

      end do DOJK
    end do DOIK

    ! sum counters from all processes
    call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod_k(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                    ( ONEPI *( bin_upper**2 - bin_lower**2 ) )

  end do DOBIN

end subroutine f_FluidstatsModeDensityBFRot
!----------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsModeDensityRRot 
!> Find out how many wave modes are in each bin.
!> @param ierr should return 0
subroutine f_FluidstatsModeDensityRRot ( ierr )

!    use d_arrays,    only : nodes_x,nodes_y,nodes_z 
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, binv, bvalexist
  use g_constants,  only : ZERO, HALF, TWO, THREE, FOUR, ONEPI
  use f_arrays,     only : i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           f_ArraysWavenumber
!    use f_arrays,    only : fluid_k, ki_min, kj_min, kk_min, &
!                                     ki_max, kj_max, kk_max

  implicit none

  PetscErrorCode,intent(inout)               :: ierr
  integer                             :: alloc_stat
  integer                             :: i,j,k,kk
  integer                             :: m,n
  integer                             :: num_modes
  integer                             :: num_temp
  real(kind=C_DOUBLE),dimension(1:3)  :: wvnmov, wvnlab
  real(kind=C_DOUBLE)                 :: k_mag
  real(kind=C_DOUBLE)                 :: bin_upper
  real(kind=C_DOUBLE)                 :: bin_lower

  ! find appropriate bin for current mode
  DOBIN: do kk = 1, kkmax, 1

    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    ! find bin limits
    bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
    bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )

  !> 3D mode density

    ! find number of modes in each bin
    DOI: do i = i_min_3w, i_max_3w, 1
    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOJ: do j = j_min_3w, j_max_3w, 1
      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
        DOK: do k = k_min_3w, k_max_3w, 1
        wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

          !> Compute wavenumber in the laboratory system.
          wvnlab = ZERO

          DOM: do m = 1, 3, 1
            DON: do n = 1, 3, 1
              IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
              end if IFBEXIST
            end do DON
          end do DOM

          !> Compute magnitude of wave vector
          k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                         ( wvnlab(2) * wvnlab(2) ) + &
                         ( wvnlab(3) * wvnlab(3) ) )

            ! check to see if k is in current bin
            IFLOWER: if ( k_mag >= bin_lower ) then
              IFUPPER: if ( k_mag < bin_upper ) then

                !> Need to count twice for conjugate symmetry except i=0
                IFI0: if ( i /= 0 ) then
                  num_temp = num_temp + 2
                else
                  num_temp = num_temp + 1
                end if IFI0

              end if IFUPPER
            end if IFLOWER

          end do DOK
        end do DOJ
      end do DOI

    ! sum counters from all processes
    call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                   ( FOUR / THREE * ONEPI *( bin_upper**3 - bin_lower**3 ) )

!    !> Correct for size of real domain vs computational domain
!    dens_mod(kk) = dens_mod(kk) *  real(B330,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE) &
!             * real(B330,kind=C_DOUBLE)  


  !> 1D mode density i direction
    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    ! find number of modes in each bin
    IFI0MODE: if ( i_min_3w <= 0 ) then
      wvnmov(1) = f_ArraysWavenumber ( 0, KI_MAX, K0I )

      DOJI: do j = j_min_3w, j_max_3w, 1
      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
        DOKI: do k = k_min_3w, k_max_3w, 1
        wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

          !> Compute wavenumber in the laboratory system.
          wvnlab = ZERO

          DOMI: do m = 1, 3, 1
            DONI: do n = 1, 3, 1
              IFBEXISTI: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
                wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
              end if IFBEXISTI
            end do DONI
          end do DOMI

          !> Compute magnitude of wave vector
          k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                         ( wvnlab(2) * wvnlab(2) ) + &
                         ( wvnlab(3) * wvnlab(3) ) )

          ! check to see if k is in current bin
          IFLOWERI: if ( k_mag >= bin_lower ) then
            IFUPPERI: if ( k_mag < bin_upper ) then

              num_temp = num_temp + 1

            end if IFUPPERI
          end if IFLOWERI

        end do DOKI
      end do DOJI
    end if IFI0MODE

    ! sum counters from all processes
    call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod_i(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                     ( ONEPI *( bin_upper**2 - bin_lower**2 ) )

!    write(*,*) ' kk, dens_mod_i(kk): ', kk, dens_mod_i(kk) 


  !> 1D mode density j direction
    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

   ! find number of modes in each bin
     wvnmov(2) = f_ArraysWavenumber ( 0, KJ_MAX, K0J )
    DOIJ: do i = i_min_3w, i_max_3w, 1
    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOKJ: do k = k_min_3w, k_max_3w, 1
      wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

        !> Compute wavenumber in the laboratory system.
        wvnlab = ZERO

        DOMJ: do m = 1, 3, 1
          DONJ: do n = 1, 3, 1
            IFBEXISTJ: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXISTJ
          end do DONJ
        end do DOMJ

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        ! check to see if k is in current bin
        IFLOWERJ: if ( k_mag >= bin_lower ) then
          IFUPPERJ: if ( k_mag < bin_upper ) then

            num_temp = num_temp + 1

          end if IFUPPERJ
        end if IFLOWERJ

      end do DOKJ
    end do DOIJ

   ! sum counters from all processes
   call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

   ! calculate mode density
   dens_mod_j(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                    ( ONEPI *( bin_upper**2 - bin_lower**2 ) )




!  !> 1D mode density k direction
    ! reset counters to zero
    num_modes = 0
    num_temp  = 0

    ! find number of modes in each bin
    wvnmov(3) = f_ArraysWavenumber ( 0, KK_MAX, K0K )
    DOIK: do i = i_min_3w, i_max_3w, 1
    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
      DOJK: do j = j_min_3w, j_max_3w, 1
      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )

        !> Compute wavenumber in the laboratory system.
        wvnlab = ZERO

        DOMK: do m = 1, 3, 1
          DONK: do n = 1, 3, 1
            IFBEXISTK: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
            end if IFBEXISTK
          end do DONK
        end do DOMK

        !> Compute magnitude of wave vector
        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
                       ( wvnlab(2) * wvnlab(2) ) + &
                       ( wvnlab(3) * wvnlab(3) ) )

        ! check to see if k is in current bin
        IFLOWERK: if ( k_mag >= bin_lower ) then
          IFUPPERK: if ( k_mag < bin_upper ) then

            num_temp = num_temp + 1

          end if IFUPPERK
        end if IFLOWERK

      end do DOJK
    end do DOIK

    ! sum counters from all processes
    call MPI_Allreduce( num_temp, num_modes, 1, MPI_Integer, MPI_Sum, PETSC_COMM_WORLD, ierr)

    ! calculate mode density
    dens_mod_k(kk) = real ( num_modes, kind=C_DOUBLE ) / &
                    ( ONEPI *( bin_upper**2 - bin_lower**2 ) )

  end do DOBIN

end subroutine f_FluidstatsModeDensityRRot
!----------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsSecondaryQuants 
!> Compute length and time scales, the Taylor-scale Reynolds number and the 
!! dissipation constant.
!> @param ierr should return 0
subroutine f_FluidstatsFsSecondaryQuants ( ierr )

  use g_constants,  only : ZERO, ONE, TWO, THREE, FOUR, FIFTEEN, ONEPI, THIRTYFIVE
  use g_parameters, only : F_NU
!    use d_arrays,    only : nodes_x

  implicit none

!    integer,intent(in)                 :: tstep
  PetscErrorCode,intent(inout)             :: ierr
  integer                           :: kk
!    real(kind=pandora_kind)           :: kmin
!    real(kind=pandora_kind)           :: kmax
!    real(kind=pandora_kind)           :: temp

! turbulent kinetic energy - fn_ken
! mean dissipation rate    - fn_dis
! dissipation skewness    - fn_disskew

! r.m.s. velocity          - u'
  fn_int = sqrt ( TWO * fn_ken / THREE )
  fn_int1 = sqrt ( TWO * fn_ken1 )
  fn_int2 = sqrt ( TWO * fn_ken2 )
  fn_int3 = sqrt ( TWO * fn_ken3 )

! turbulent time scale     - t_{epsilon}
  tn_tur = fn_ken / fn_dis
  
! turbulent length scale   - l_{epsilon}
  ln_tur = ( fn_ken**( THREE / TWO ) ) / fn_dis
  
! Taylor microscale        - lambda
  ln_tay = sqrt( ( fifteen * F_NU * ( fn_int**2 ) ) / fn_dis )

! Taylor Reynolds number   - r_{lambda}
  rn_tay = ln_tay * fn_int / F_NU

! longitudinal integral length scale

  ln_int = ZERO
  ln_int1 = ZERO
  ln_int2 = ZERO
  ln_int3 = ZERO

!> Pope 6.225
  DOMODES: do kk = 1, kkmax
    ln_int = ln_int + ( ( ONEPI / TWO / fn_int**2 ) * spec_ken(kk) / real(kk, kind=C_DOUBLE) )
    ln_int1 = ln_int1 + ( ( ONEPI / TWO / fn_int**2 ) * spec_ken1(kk) / real(kk, kind=C_DOUBLE) )
    ln_int2 = ln_int2 + ( ( ONEPI / TWO / fn_int**2 ) * spec_ken2(kk) / real(kk, kind=C_DOUBLE) )
    ln_int3 = ln_int3 + ( ( ONEPI / TWO / fn_int**2 ) * spec_ken3(kk) / real(kk, kind=C_DOUBLE) )
  end do DOMODES

! eddy turnover time
  tn_edy = ln_int / fn_int

! Kolmogorov length scale
  ln_kol = ( F_NU**3 / fn_dis )**( ONE / FOUR ) 

! Kolmogorov time scale
  tn_kol = sqrt ( F_NU / fn_dis )

  disskew = fn_disskew * ( ( fifteen / fn_dis )**( THREE / TWO ) ) * sqrt ( F_NU ) / THIRTYFIVE

! Dissipation constant

  C_epsilon = fn_dis * ln_int / ( fn_int**3 )

end subroutine f_FluidstatsFsSecondaryQuants
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsHomogeneous 
!> Compute energy and dissipation spectra and get turbulent kinetic energy,
!! dissipation and dissipation skew from these.
!> @param tstep current time step
!> @param ierr should return 0
subroutine f_FluidstatsHomogeneous ( tstep, ierr )

  use g_constants,  only : ZERO, TWO, ONEPI, HALF
  use g_parameters, only : MYRANK, MASTER, F_DEFORM_S
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K, &
                           bmat, binv, bvalexist
  use g_files,      only : FUNIT_FSTATH
  use f_arrays,     only : da3w, i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           f_ArraysWavenumber
 
  implicit none
 
  integer,intent(in)                     :: tstep
  PetscErrorCode,intent(inout)                  :: ierr

  real(kind=C_DOUBLE)                    :: phi11_i, phi22_i, phi33_i
  real(kind=C_DOUBLE)                    :: phi11_j, phi22_j, phi33_j
  real(kind=C_DOUBLE)                    :: phi11_k, phi22_k, phi33_k
  real(kind=C_DOUBLE)                    :: u1u1_i, u2u2_i, u3u3_i
  real(kind=C_DOUBLE)                    :: u1u1_j, u2u2_j, u3u3_j
  real(kind=C_DOUBLE)                    :: u1u1_k, u2u2_k, u3u3_k
!  real(kind=C_DOUBLE)                    :: spec1_i, spec2_i, spec3_i
  real(kind=C_DOUBLE)                    :: spec1_j, spec2_j, spec3_j
  real(kind=C_DOUBLE)                    :: spec1_k, spec2_k, spec3_k
  real(kind=C_DOUBLE)                    :: ek1, ek2, ek3
  integer                                :: i, j, k, m, n, kk

  real(kind=C_DOUBLE)                    :: k_mag, bin_lower, bin_upper
  real(kind=C_DOUBLE),dimension(1:3)     :: wvnmov, wvnlab
  real(kind=C_DOUBLE),dimension(0:1,1:3) :: umov, ulab

  !> Pope 6.213 
!  ln_int_11 = ONEPI * spec_1i(0) / ( TWO * sum(spec_1i) )
  ln_int_11 = ONEPI * spec_1i(0) / ( TWO * fs_u1u1 )
  !> @todo need to use actual energy
!  ln_int_12 = ZERO!ONEPI * spec2_i / fs_u2u2
!  ln_int_13 = ZERO!ONEPI * spec3_i / fs_u3u3

!  !> Initialise to zero
!  u1u1_i = ZERO
!  u2u2_i = ZERO
!  u3u3_i = ZERO
!  u1u1_j = ZERO
!  u2u2_j = ZERO
!  u3u3_j = ZERO
!  u1u1_k = ZERO
!  u2u2_k = ZERO
!  u3u3_k = ZERO

!  !> Calculate streamwise 1-D energy spectrum
!  DOII: do i = i_min_3w, i_max_3w, 1
!    wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )

!    !> Initialise energy on this mode
!    temp_spec_1i(i) = ZERO
!    temp_ken11 = ZERO
!    temp_ken21 = ZERO
!    temp_ken31 = ZERO

!    DOJI: do j = j_min_3w, j_max_3w, 1
!      wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )
!      DOKI: do k = k_min_3w, k_max_3w, 1
!        wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

!        !> Compute wavenumber and velocity in the laboratory system.
!        umov(:,1) = arr_u1_3w(:,k,j,i)
!        umov(:,2) = arr_u2_3w(:,k,j,i)
!        umov(:,3) = arr_u3_3w(:,k,j,i)

!        wvnlab = ZERO
!        ulab = ZERO

!        DOMI: do m = 1, 3, 1
!          DONI: do n = 1, 3, 1
!            IFBEXISTI: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
!              ulab(:,m) = ulab(:,m) + binv(m,n) * umov(:,n)
!              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
!            end if IFBEXISTI
!          end do DONI
!        end do DOMI

!        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
!                       ( wvnlab(2) * wvnlab(2) ) + &
!                       ( wvnlab(3) * wvnlab(3) ) )

!        phi11_i = ulab(0,1) * ulab(0,1) + ulab(1,1) * ulab(1,1) 
!        phi22_i = ulab(0,2) * ulab(0,2) + ulab(1,2) * ulab(1,2) 
!        phi33_i = ulab(0,3) * ulab(0,3) + ulab(1,3) * ulab(1,3) 


!        ! calculate energy spectrum (see Pope 6.210) 
!        ek1 = TWO * phi11_i
!        ek2 = TWO * phi22_i

!        ! find appropriate bin for current mode
!        DOBINI: do kk = 1, kkmax, 1

!          ! find bin limits
!          bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
!          bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )

!          ! check to see if k is in current bin
!          IFLOWERI: if ( k_mag >= bin_lower ) then
!            IFUPPERI: if ( k_mag < bin_upper ) then

!              temp_ken11(kk) = temp_ken11(kk) + ek1
!              temp_ken21(kk) = temp_ken21(kk) + ek2
!              temp_ken31(kk) = temp_ken31(kk) + ek3

!            end if IFUPPERI
!          end if IFLOWERI

!        end do DOBINI

!!        temp_spec_1i(i) = temp_spec_1i(i) + ek1

!      end do DOKI
!    end do DOJI

!    !> Sum up over all j and k on this process and normalise by mode density
!    temp_spec_1i(i) = ZERO
!    temp_spec_2i(i) = ZERO
!    DOSUMKK: do kk = 1, kkmax, 1
!      IFNOTEMPTY: if ( dens_mod_i(kk) > 1.0e-8 ) then
    
!        temp_spec_1i(i) = temp_spec_1i(i) + ( temp_ken11(kk) / dens_mod_i(kk) ) 
!        temp_spec_2i(i) = temp_spec_2i(i) + ( temp_ken21(kk) / dens_mod_i(kk) ) 
!
!      end if IFNOTEMPTY
!    end do DOSUMKK

!  end do DOII

!  !> sum contributions from all processes - they are already normalised by mode density
!  call MPI_Allreduce ( temp_spec_1i, spec_1i, NODES_KI, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_spec_2i, spec_2i, NODES_KI, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

!!  write(*,*) ' temp_spec_1i: ', temp_spec_1i
!  write(*,*) ' spec_1i: ', spec_1i



!  !> Calculate 1-D energy spectrum on j = 0 plane
!  IFJ: if ( j_min_3w <= 0 ) then

!    wvnmov(2) = f_ArraysWavenumber ( 0, KJ_MAX, K0J )
!    DOIJ: do i = i_min_3w, i_max_3w, 1

!      wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
!      DOKJ: do k = k_min_3w, k_max_3w, 1
!        wvnmov(3) = f_ArraysWavenumber ( k, KK_MAX, K0K )

!        !> Compute wavenumber and velocity in the laboratory system.
!        umov(:,1) = arr_u1_3w(:,k,0,i)
!        umov(:,2) = arr_u2_3w(:,k,0,i)
!        umov(:,3) = arr_u3_3w(:,k,0,i)

!        wvnlab = ZERO
!        ulab = ZERO

!        DOMJ: do m = 1, 3, 1
!          DONJ: do n = 1, 3, 1
!            IFBEXISTJ: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
!              ulab(:,m) = ulab(:,m) + binv(m,n) * umov(:,n)
!              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
!            end if IFBEXISTJ
!          end do DONJ
!        end do DOMJ

!        !> Compute magnitude of wave vector
!        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
!                       ( wvnlab(2) * wvnlab(2) ) + &
!                       ( wvnlab(3) * wvnlab(3) ) )

!        phi11_j = ulab(0,1) * ulab(0,1) + ulab(1,1) * ulab(1,1) 
!        phi22_j = ulab(0,2) * ulab(0,2) + ulab(1,2) * ulab(1,2) 
!        phi33_j = ulab(0,3) * ulab(0,3) + ulab(1,3) * ulab(1,3) 

!        IFJNOTIZERO: if (i /= 0) then
!          phi11_j = phi11_j * TWO
!          phi22_j = phi22_j * TWO
!          phi33_j = phi33_j * TWO

!        end if IFJNOTIZERO

!        !> calculate components of Reynolds stress tensor
!        u1u1_j = u1u1_j + phi11_j
!        u2u2_j = u2u2_j + phi22_j
!        u3u3_j = u3u3_j + phi33_j

!        ! calculate energy spectrum 
!        ek1 = HALF * phi11_j
!        ek2 = HALF * phi22_j
!        ek3 = HALF * phi33_j
!
!
!        ! find appropriate bin for current mode
!        DOBINJ: do kk = 1, kkmax, 1
!
!          ! find bin limits
!          bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
!          bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )
!
!          ! check to see if k is in current bin
!          IFLOWERJ: if ( k_mag >= bin_lower ) then
!            IFUPPERJ: if ( k_mag < bin_upper ) then
!
!              temp_ken12(kk) = temp_ken12(kk) + ek1
!              temp_ken22(kk) = temp_ken22(kk) + ek2
!              temp_ken32(kk) = temp_ken32(kk) + ek3
!
!            end if IFUPPERJ
!          end if IFLOWERJ
!
!        end do DOBINJ

!      end do DOKJ
!    end do DOIJ

!  end if IFJ

!  !> Calculate 1-D energy spectrum on k = 0 plane
!  IFK: if ( k_min_3w <= 0 ) then

!    wvnmov(3) = f_ArraysWavenumber ( 0, KK_MAX, K0K )
!    DOIK: do i = i_min_3w, i_max_3w, 1

!      wvnmov(1) = f_ArraysWavenumber ( i, KI_MAX, K0I )
!      DOJK: do j = j_min_3w, j_max_3w, 1

!        wvnmov(2) = f_ArraysWavenumber ( j, KJ_MAX, K0J )

!        !> Compute wavenumber and velocity in the laboratory system.
!        umov(:,1) = arr_u1_3w(:,0,j,i)
!        umov(:,2) = arr_u2_3w(:,0,j,i)
!        umov(:,3) = arr_u3_3w(:,0,j,i)

!        wvnlab = ZERO
!        ulab = ZERO

!        DOMK: do m = 1, 3, 1
!          DONK: do n = 1, 3, 1
!            IFBEXISTK: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
!              ulab(:,m) = ulab(:,m) + binv(m,n) * umov(:,n)
!              wvnlab(n) = wvnlab(n) + bmat(m,n) * wvnmov(m)
!            end if IFBEXISTK
!          end do DONK
!        end do DOMK

!        !> Compute magnitude of wave vector
!        k_mag = sqrt ( ( wvnlab(1) * wvnlab(1) ) + &
!                       ( wvnlab(2) * wvnlab(2) ) + &
!                       ( wvnlab(3) * wvnlab(3) ) )

!        phi11_k = ulab(0,1) * ulab(0,1) + ulab(1,1) * ulab(1,1) 
!        phi22_k = ulab(0,2) * ulab(0,2) + ulab(1,2) * ulab(1,2) 
!        phi33_k = ulab(0,3) * ulab(0,3) + ulab(1,3) * ulab(1,3) 

!        IFKNOTIZERO: if (i /= 0) then
!          phi11_k = phi11_k * TWO
!          phi22_k = phi22_k * TWO
!          phi33_k = phi33_k * TWO

!        end if IFKNOTIZERO

!        !> calculate components of Reynolds stress tensor
!        u1u1_k = u1u1_k + phi11_k
!        u2u2_k = u2u2_k + phi22_k
!        u3u3_k = u3u3_k + phi33_k

!        ! calculate energy spectrum 
!        ek1 = HALF * phi11_k
!        ek2 = HALF * phi22_k
!        ek3 = HALF * phi33_k
!
!
!        ! find appropriate bin for current mode
!        DOBINK: do kk = 1, kkmax, 1
!
!          ! find bin limits
!          bin_upper = k_min * ( real(kk,kind=C_DOUBLE) + HALF )
!          bin_lower = k_min * ( real(kk,kind=C_DOUBLE) - HALF )
!
!          ! check to see if k is in current bin
!          IFLOWERK: if ( k_mag >= bin_lower ) then
!            IFUPPERK: if ( k_mag < bin_upper ) then
!
!              temp_ken13(kk) = temp_ken13(kk) + ek1
!              temp_ken23(kk) = temp_ken23(kk) + ek2
!              temp_ken33(kk) = temp_ken33(kk) + ek3
!
!            end if IFUPPERK
!          end if IFLOWERK
!
!        end do DOBINK

!      end do DOJK
!    end do DOIK

!  end if IFK

!  write(*,*) ' imin, jmin, kmin: ', i_min_3w, j_min_3w, k_min_3w, &
!             ' u1u1, ... :       ', u1u1_i, u2u2_i, u3u3_i, u1u1_j, &
!             u2u2_j, u3u3_j, u1u1_k, u2u2_k, u3u3_k

!  call MPI_Allreduce ( temp_ken21, spec_ken21, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken31, spec_ken31, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken12, spec_ken12, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken22, spec_ken22, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken32, spec_ken32, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken13, spec_ken13, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken23, spec_ken23, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( temp_ken33, spec_ken33, kkmax, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u1u1_i, spec1_i, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u2u2_i, spec2_i, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u3u3_i, spec3_i, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u1u1_j, spec1_j, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u2u2_j, spec2_j, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u3u3_j, spec3_j, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u1u1_k, spec1_k, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u2u2_k, spec2_k, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
!  call MPI_Allreduce ( u3u3_k, spec3_k, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

!  write(*,*) 'spec1_i: ', spec1_i, ', spec_ken11: ', spec_ken11, &
!             ' Mode density: ', dens_mod_i

!  spec_ken21 = spec_ken21 / dens_mod_i ! normalised energy spectrum
!  spec_ken31 = spec_ken31 / dens_mod_i ! normalised energy spectrum
!  spec_ken12 = spec_ken12 / dens_mod_j ! normalised energy spectrum
!  spec_ken22 = spec_ken22 / dens_mod_j ! normalised energy spectrum
!  spec_ken32 = spec_ken32 / dens_mod_j ! normalised energy spectrum
!  spec_ken13 = spec_ken13 / dens_mod_k ! normalised energy spectrum
!  spec_ken23 = spec_ken23 / dens_mod_k ! normalised energy spectrum
!  spec_ken33 = spec_ken33 / dens_mod_k ! normalised energy spectrum


!  Sstar = F_DEFORM_S * TWO * fn_ken / fn_dis
  Sstar = tn_edy * F_DEFORM_S
  Kstar = TWO * fs_u1u1 / ( fs_u2u2 + fs_u3u3 )
  Lstar = ln_int_11 / ( TWO * ln_int_21 )

  IFMASTER: if ( MYRANK == MASTER ) then
    write(FUNIT_FSTATH,60) tstep, &
                           Sstar, &
                           Kstar, &
                           Lstar, &
                           ln_int1, &
                           ln_int2, &
                           ln_int3, &
                           ln_int_11, &
                           ln_int_12, &
                           ln_int_13, &
                           fn_ken1, &
                           fn_ken2, &
                           fn_ken3, &
                           fn_int1, &
                           fn_int2, &
                           fn_int3
  end if IFMASTER

60  format ( ' ts: ',i6, ' S*: ', e14.6, ' K*: ', e14.6, &
             ' L*: ', e14.6' ln_int1: ', e14.6, ' ln_int2: ', e14.6, ' ln_int3: ', e14.6, &
             ' l11i: ', e14.6, ' l22i: ', e14.6, ' l33i: ', e14.6, &
             ' tke1: ', e14.6, ' tke2: ', e14.6, ' tke3: ', e14.6, &
             ' u1^\prime: ', e14.6, ' u2^\prime: ', e14.6, ' u3^\prime: ', e14.6)

end subroutine f_FluidstatsHomogeneous
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsAnalyticalViscous
!> Compute energy and dissipation spectra and get turbulent kinetic energy,
!! dissipation and dissipation skew from these.
!> @param tstep current time step
!> @param ierr should return 0
subroutine f_FluidstatsAnalyticalViscous ( tstep, t, ierr )

  use g_constants,  only : ZERO, HALF, ONE, TWO, THREE, ONEPI, TWOPI, FOURPI
  use g_parameters, only : MYRANK, MASTER, F_DEFORM_S, F_NU
  use g_domain,     only : NODES_KI, NODES_KJ, NODES_KK, KI_MAX, K0I, KJ_MAX, K0J, KK_MAX, K0K
  use f_arrays,     only : da3w, i_min_3w, i_max_3w, &
                           j_min_3w, j_max_3w, k_min_3w, k_max_3w, &
                           arr_u1_3w, arr_u2_3w, arr_u3_3w, &
                           f_ArraysWavenumber
  use g_files,      only : FUNIT_FSTATF, FUNIT_ESPEC, FUNIT_FSTATH
 
  implicit none
 
  integer,intent(in)             :: tstep
  real(kind=C_DOUBLE),intent(in) :: t
  PetscErrorCode,intent(inout)          :: ierr

  real(kind=C_DOUBLE)   :: alpha, phi
  real(kind=C_DOUBLE)   :: dalpha, dphi, dk0, dk
  real(kind=C_DOUBLE)   :: k0, k1, k2, k3
  real(kind=C_DOUBLE)   :: beta, psi, zeta, D, Q1, Q2
  real(kind=C_DOUBLE)   :: phi11_an, phi22_an, phi33_an, phi12_an
  real(kind=C_DOUBLE)   :: u1u1_an, u2u2_an, u3u3_an, u1u2_an
  real(kind=C_DOUBLE)   :: fs_u1u1_an, fs_u2u2_an, fs_u3u3_an, fs_u1u2_an
  integer               :: k, a, p

  !> Initialise to zero
  u1u1_an = ZERO
  u2u2_an = ZERO
  u3u3_an = ZERO
  u1u2_an = ZERO

  dalpha = ONEPI / real ( NODES_KI, kind=C_DOUBLE )
  dphi   = TWOPI / real ( NODES_KJ, kind=C_DOUBLE ) 
  dk0    = ONE

  !> Calculate spectrum 
  DOK: do k = 1, kkmax, 1
    DOA: do a = i_min_3w, i_max_3w, 1
      DOP: do p = j_min_3w, j_max_3w, 1

        alpha = real ( a, kind=C_DOUBLE ) / real ( NODES_KI, kind=C_DOUBLE ) * ONEPI
        phi   = real ( p, kind=C_DOUBLE ) / real ( NODES_KJ, kind=C_DOUBLE ) * TWOPI

        beta = F_DEFORM_S * t

        k0 = real ( k, kind=C_DOUBLE ) * ONE 
        k1 = k0 * cos ( alpha )
        k2 = k0 * sin ( alpha ) * sin ( phi )
        k3 = k0 * sin ( alpha ) * cos ( phi )

        dk = k0*k0 * sin(alpha) * dalpha * dphi * dk0

        psi = atan2 ( beta * k1 * ( k1*k1 + k3*k3 )**0.5, &
                      k0*k0 - beta * k1 * k2 )

        zeta = ( k0*k0 - ( TWO * k2*k2 + beta * k1 * k2 ) ) / &
               ( k0*k0 - ( TWO * beta * k1 * k2 + beta*beta * k1*k1 ) )

        D = exp ( - F_NU * t * ( k0*k0 - ( beta * k1 * k2 ) + &
                               ( beta*beta * k1*k1 / THREE ) ) )

        Q1 = ( - k3*k3 / ( k1 * ( k1*k1 + k3*k3 )**0.5 ) ) * psi &
             + beta * psi * ( k1*k1 / (k0*k0) )
        
        Q2 = ( k3 / ( k1*k1 + k3*k3 )**0.5 ) * psi + &
             ( beta * zeta * ( k1*k1 / ( k0*k0 ) ) )

        phi11_an = spec_ken(k) / ( FOURPI * k0*k0 ) * D*D * &
                   ( k0*k0 / ( k1*k1 + k3*k3 ) * Q1 * ( Q1 - &
                     ( ( TWO * k1 * k2 ) / ( k0*k0 ) ) ) &
                     + ( k2*k2 + k3*k3 ) / ( k0*k0) )

        phi22_an = spec_ken(k) / ( FOURPI * k0*k0 ) * D*D * &
                   ( k0*k0 * ( k1*k1 + k3*k3 ) ) / &
                   ( k0*k0 - TWO * beta * k1 * k2 + beta*beta * k1*k1 )**2

        phi33_an = spec_ken(k) / ( FOURPI * k0*k0 ) * D*D * &
                   ( k0*k0 / ( k1*k1 + k3*k3 ) * Q2 * ( Q2 - &
                     ( ( TWO * k2 * k3 ) / ( k0*k0 ) ) ) &
                    + ( k1*k1 + k2*k2 ) / ( k0*k0) )

        phi12_an = spec_ken(k) / ( FOURPI * k0*k0 ) * D*D * &
                   ( ( k0*k0 * Q1 ) - ( k1 * k2 ) ) / &
                   ( k0*k0 - TWO * beta * k1 * k2 + beta*beta * k1*k1 )**2

        !> calculate components of Reynolds stress tensor
        u1u1_an = u1u1_an + ( phi11_an * dk )
        u2u2_an = u2u2_an + ( phi22_an * dk )
        u3u3_an = u3u3_an + ( phi33_an * dk )
        u1u2_an = u1u2_an + ( phi12_an * dk )

      end do DOP
    end do DOA
  end do DOK


  !> Determine Reynolds stresses
  call MPI_Allreduce ( u1u1_an, fs_u1u1_an, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2_an, fs_u2u2_an, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3_an, fs_u3u3_an, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2_an, fs_u1u2_an, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

  IFMASTER: if ( MYRANK == MASTER ) then

    write(FUNIT_FSTATF,*) tstep, &
                          t, &
                          fs_u1u1_an, &
                          fs_u2u2_an, &
                          fs_u3u3_an, &
                          fs_u1u2_an

  end if IFMASTER

end subroutine f_FluidstatsAnalyticalViscous
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine f_FluidstatsFsWriteQuants 
!> Write fluid wave-space statistics and fluid spectra to the corresponding files.
!> @param tstep current time step
!> @param time current simulation time
!> @param ierr should return 0
subroutine f_FluidstatsFsWriteQuants ( tstep, time, ierr )

  use g_parameters, only : MYRANK, MASTER, F_TYPE, F_ISOTROPIC, TSTEPS, TS_SPEC3D, TS_SPEC1D
  use g_files,      only : FUNIT_FSTATF, FUNIT_ESPEC, FUNIT_FSTATH, FUNIT_1DSPEC
  use g_domain,     only : KI_MAX
!    use g_files,      only : funit_fstats, funit_test
!    use f_arrays, only : s

!    implicit none

  PetscErrorCode,intent(inout)          :: ierr
  integer,intent(in)             :: tstep
  real(kind=C_DOUBLE),intent(in) :: time

  integer                        :: kk, i
!    real(kind=pandora_kind) :: q2, nd_time

!    nd_time = s * time


!    endif

  IFMASTER: if ( MYRANK == MASTER ) then

   IFSPEC3D: if ( TS_SPEC3D <= TSTEPS ) then 
      !> if the correct timestep, write energy and dissipation spectrum to spectra file as a new record
      IFSPEC3DTS: if ( modulo(tstep, TS_SPEC3D) == 0 ) then
        write ( FUNIT_ESPEC, 10 )
        write ( FUNIT_ESPEC, 11 ) tstep

        DOWRITENODES: do kk = 1, kkmax, 1
          write ( FUNIT_ESPEC, 20 ) ( k_min * real(kk,kind=C_DOUBLE) ), &
                                    spec_ken(kk), &
                                    spec_dis(kk), &
                                    spec_disskew(kk), &
                                    spec_ken1(kk), &
                                    spec_ken2(kk), &
                                    spec_ken3(kk)
        end do DOWRITENODES
      end if IFSPEC3DTS
    end if IFSPEC3D

    ! write output quantities to fluid analysis file
    write(FUNIT_FSTATF,50) tstep,             & 
                           time,              & 
                           rn_tay,            &
                           fn_int,            &
                           fn_ken,            &
                           fn_dis,            &
                           ln_tur,            &
                           tn_tur,            &
                           ln_tay,            &
                           ln_int,            &
                           tn_edy,            &
                           ln_kol,            &
                           tn_kol,            &
!                           rs_ens,            &
                           disskew,           &
                           C_epsilon,         &
!                                global_ke,         &
!                                f_ken,             &
!                                su
                           fs_u1u1,           & 
                           fs_u2u2,           &
                           fs_u3u3,           &
                           fs_u1u2,           &
                           fs_u1u3,           &
                           fs_u2u3

  end if IFMASTER

10  format( '' )
11  format( ' ts: ',i6 )
20  format( 7e14.6 )
21  format( 3e14.6 )
50  format ( ' ts: ',i6, ' time: ', e14.6, ' Re_\lambda: ', e14.6, ' u^\prime: ', e14.6, &
             ' tke: ', e14.6, ' \epsilon: ', e14.6, ' ln_tur: ', e14.6, ' tn_tur: ', e14.6, &
             ' ln_tay: ', e14.6, ' ln_int: ', e14.6, ' tn_edy: ', e14.6, ' ln_kol: ', e14.6, &
             ' tn_kol: ', e14.6, ' disskew: ', e14.6, ' C_epsilon: ', e14.6, &
             ' u1u1: ', e14.6,' u2u2: ', e14.6,' u3u3: ', e14.6, &
             ' u1u2: ', e14.6,' u1u3: ', e14.6,' u2u3: ', e14.6)
!50  format(i6,18e14.6)

end subroutine f_FluidstatsFsWriteQuants
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
end module
