! module g_rk3
!> Contains all the routines that control the stages of each Runge-Kutta time step.
module g_rk3

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
use g_constants, only : HALF

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>


private

!---------------------------------------------------
!g_rk3 parameters 
!---------------------------------------------------

real(kind=C_DOUBLE),save        :: nd_remesh_time = HALF
real(kind=C_DOUBLE),save        :: remesh_time   = 0.0

!integer,save :: bfirstremesh, tlocal

public :: g_RK3Init
public :: g_RK3Calc
!public :: g_RK3Remesh
!public :: g_RK3Tstep

contains

!---------------------------------------------------------------------------------
! subroutine g_RK3Init
!> Initialise Runge-Kutta routines
!> @param ierr should return 0
subroutine g_RK3Init ( ierr )

!  use g_constants, only : TWO, ZERO
  use g_parameters, only : MYRANK, MASTER, F_DEFORM_S
!  use f_arrays, only : s

  implicit none
 
  PetscErrorCode,intent(inout)           :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,10)
    write(*,20)
  end if IFMASTER

  call g_RK3Check ( ierr )

!    tlocal = 0
!    time_local = zero
!    bfirstremesh = yes

!    remesh_time = length_x / ( two * s * length_y )
    remesh_time = nd_remesh_time / F_DEFORM_S

10  format('==> initialising g_rk3 module')
20  format('    | checking Runge-Kutta parameters')

end subroutine g_RK3Init

!---------------------------------------------------------------------------------
! subroutine g_RK3Calc 
!> Coordinate all four stages of the Runge-Kutta fluid and particle routines
!! and corresponding auxiliary routines.
!> @param tstep number of the current time step 
!> @param simtime time in the simulation time frame
!> @param ierr should return 0
subroutine g_RK3Calc ( tstep, simtime, ierr )

  use g_constants,    only : ONE, THREE, FOUR

  implicit none

  integer,intent(in)                      :: tstep
  real(kind=C_DOUBLE),intent(inout)       :: simtime
  PetscErrorCode,intent(inout)                   :: ierr
  real(kind=C_DOUBLE)                     :: tstart !,tl
  real(kind=C_DOUBLE)                     :: sim_dt = 99999.0
!    real(kind=8),dimension(1:12)            :: timingtags
!    integer                                 :: i
 
!    tlocal = tlocal + 1
!
  tstart = simtime

!> Determine the new time step resolution
!  call g_RK3Tstep ( rk3_dt, ierr )
  
!    call f_force_calc(tstep,rk3_dt,ierr)
  
!> Runge-Kutta stage 1

!> Compute all quantities, then update time, then update mesh.
  call g_RK3Stage1 ( tstep, simtime, sim_dt, ierr )
  simtime = tstart + (ONE / THREE) * sim_dt
!  call g_RK3UpdateB ( 1, sim_dt, simtime, ierr)
  call g_RK3UpdateB ( tstep, sim_dt, simtime, ierr)

!> Runge-Kutta stage 2

!> Compute all quantities, then update time, then update mesh.
  call g_RK3Stage2 ( tstep, simtime, sim_dt, ierr )
  simtime = tstart + (THREE/FOUR) * sim_dt
!  call g_RK3UpdateB ( 2, sim_dt, simtime, ierr)
  call g_RK3UpdateB ( tstep, sim_dt, simtime, ierr)

!> Runge-Kutta stage 3

!> Compute all quantities, then update time, then update mesh.
  call g_RK3Stage3 ( tstep, simtime, sim_dt, ierr )
  simtime = tstart + sim_dt
!  call g_RK3UpdateB ( 3, sim_dt, simtime, ierr)
  call g_RK3UpdateB ( tstep, sim_dt, simtime, ierr)

end subroutine g_RK3Calc

!---------------------------------------------------------------------------------
! subroutine g_RK3Remesh
!> Remesh the flow field.
!> @param ierr should return 0
!subroutine g_RK3Remesh ( ierr )
  
!  use g_constants, only : ZERO, HALF
!  use g_domain,    only : KK_MIN, KK_MAX, bmat, binv
!  use f_arrays,    only : da3w, u1_3w, u2_3w, u3_3w, &
!                          arr_u1_3w, arr_u2_3w, arr_u3_3w, &
!                          arr_u1_remesh, arr_u2_remesh, arr_u3_remesh, &
!                          i_min_3w, j_min_3w, k_min_3w, &
!                          i_max_3w, j_max_3w, k_max_3w

!  implicit none
!   
!  PetscErrorCode,intent(inout) :: ierr
!  
!  integer               :: cutoff, actualk, knew, nodes
!  integer               :: i, j, k

!> Get write access to all velocity components
!  call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!  call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!  call DMDAVecGetArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

!  nodes = 2 * KK_MAX + 1

!  cutoff = KK_MAX

!!> @todo we are assuming equal number of wave modes in each direction here. might need modifying for longer domains
!  DOI: do i = i_min_3w, i_max_3w, 1
!    DOJ: do j = j_min_3w, j_max_3w, 1

!      !> initialise remesh buffer 
!      arr_u1_remesh = ZERO
!      arr_u2_remesh = ZERO
!      arr_u3_remesh = ZERO

!      DOK: do k = k_min_3w, k_max_3w, 1
!
!        IFOLDNEGATIVE: if ( k > KK_MAX ) then
!          actualk = k - nodes
!        else
!          actualk = k
!        end if IFOLDNEGATIVE
!
!        ! i has no negative wavenumbers, so this is always true:
!        knew = actualk - i
!        knew = actualk + i
!
        !> Only remesh the wave numbers that are not aliased
!        IFNOTALIASED: if ( knew >= - KK_MAX ) then
!        IFNOTALIASED: if ( knew <= KK_MAX ) then
!!        IFNOTALIASED: if ( knew <= (HALF*KK_MAX) ) then
!
          !> For negative wave numbers we need to go back to FFTW indexing
!          IFNEWNEGATIVE: if ( knew < 0 ) then
!            knew = knew + nodes
!          end if IFNEWNEGATIVE

!!          write(*,*) 'i, j: ', i, j, ', k old: ', actualk, ', k new: ', knew

!          arr_u1_remesh(:,knew) = arr_u1_3w(:,k,j,i)
!          arr_u2_remesh(:,knew) = arr_u2_3w(:,k,j,i)
!          arr_u3_remesh(:,knew) = arr_u3_3w(:,k,j,i)

!        end if IFNOTALIASED

!      enddo DOK

      !> Correct velocity from moving system before remesh to moving
      !! system after remesh. 

      !> Only u1 component changes.
      !> B**(-1) for conversion to laboratory system before remesh (St=0.5).
!      binv(1,3) =               HALF / bmat(3,3)
      !> B for conversion to moving system after remesh (St=-0.5).
!      bmat(1,3) = - bmat(1,1) * ( - HALF )

!      write(*,*) ' i, j: ', i, j, &
!                 ', u1 before: ', arr_u1_3w(:,:,j,i), ', u1 after: ', &
!                 ( arr_u1_remesh(:,:) )!+ &
                ! ( binv(1,3) * arr_u3_remesh(:,:) ) ) + &
                ! ( bmat(1,3) * arr_u3_remesh(:,:) )
!      arr_u1_3w(:,:,j,i) = ( arr_u1_remesh(:,:) + &
!                             ( binv(1,3) * arr_u3_remesh(:,:) ) ) + &
!                           ( bmat(1,3) * arr_u3_remesh(:,:) ) 

!      arr_u1_3w(:,:,j,i) = arr_u1_remesh(:,:) 
      !> u2 and u3 stay the same.
!      arr_u2_3w(:,:,j,i) = arr_u2_remesh(:,:)
!      arr_u3_3w(:,:,j,i) = arr_u3_remesh(:,:)

!    enddo DOJ
!  enddo DOI

!> Return velocity components to PETSc
!  call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!  call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!  call DMDAVecRestoreArrayF90 ( da3w, u3_3w, arr_u3_3w, ierr )

!end subroutine g_RK3Remesh

!---------------------------------------------------------------------------------
! subroutine g_RK3Check(ierr)
!> Checks for the Courant-Friedrichs-Lewy criterion
!> @param ierr should return 0

subroutine g_RK3Check ( ierr )

  use g_parameters, only : MYRANK, MASTER, RK3_CFL
  use g_constants, only: ONE

  implicit none

  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: status

  status = 0

  IFMASTER: if ( MYRANK == MASTER ) then
    IFCFL: if ( RK3_CFL >= ONE ) then
      write(*,10)
        status = 1
      end if IFCFL
  end if IFMASTER

  call MPI_Bcast ( status, 1, MPI_Integer, 0, PETSC_COMM_WORLD, ierr )

  IFABORT: if ( status /= 0 ) then
    call MPI_Abort ( PETSC_COMM_WORLD, status, ierr )
  end if IFABORT

10  format('STOP: Courant number is greater than 1.0 !!!')

end subroutine g_RK3Check
!--------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
! subroutine g_RK3UpdateB
!> Updates the Rogallo transform tensor. Use local time for shear or simulation
!! time for other flow.
!> @param tstep number of the current time step 
!> @param rk3dt \f$ \Delta t \f$ based on the CFL criterion or fixed when the 
!! Rogallo transform is used
!> @param simtime time in the simulation time frame
!> @param ierr should return 0
subroutine g_RK3UpdateB ( tstep, rk3dt, simtime, ierr )
!subroutine g_RK3UpdateB ( rk3stage, rk3dt, simtime, ierr )

  use g_parameters, only : MYRANK, MASTER, REPORT, YES, RK3_TIMESTEP, &
                           RK3_TIMESTEP_DECIDE, &
                           F_TYPE, F_ISOTROPIC, F_REMESH_T, &
                           B110, B220, B330, &
                           F_DEFORM_A, F_DEFORM_B, F_DEFORM_C, F_DEFORM_S
  use g_constants,  only : ZERO, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN, &
                           ONE, THREE, FOUR
  use g_domain,     only : amat, rk3mat, bmat, binv

  implicit none

  integer,intent(in)                      :: tstep 
!  integer,intent(in)                      :: rk3stage
  real(kind=C_DOUBLE),intent(in)          :: rk3dt
  real(kind=C_DOUBLE),intent(in)          :: simtime
  PetscErrorCode,intent(inout)                   :: ierr

  integer                                 :: localtstep
  real(kind=C_DOUBLE)                     :: localtime
  real(kind=C_DOUBLE)                     :: a, b
!  real(kind=C_DOUBLE)                     :: rkbeta
!  real(kind=C_DOUBLE)                     :: rkgamma


  IFFTYPE: if ( F_TYPE /= F_ISOTROPIC ) then

    IFREMESH: if ( RK3_TIMESTEP == RK3_TIMESTEP_DECIDE ) then

!      !> Local time according to Runge-Kutta stage.

      IFFIRSTREMESH: if ( tstep <= ( B110 * F_REMESH_T ) ) then
!      IFFIRSTREMESH: if ( tstep <= F_REMESH_T ) then

        !> Nothing to do before first remesh, just use simulation time as is
        localtime = simtime

      else 

        !> Correct for stretched domain
        a = tstep - ( B110 * F_REMESH_T ) - 1
!        a = tstep - F_REMESH_T - 1
        b = 2 * ( B110 * F_REMESH_T )
!        b = 2 * F_REMESH_T

        !> Determine number of time steps after last remesh          
        !> Determine which remesh this is
        localtstep =  a / b
        !> Compute time passed since last remesh
        localtime = simtime - ( real ( localtstep + 1, kind=C_DOUBLE ) / &
                    F_DEFORM_S * real(B110,kind=C_DOUBLE) ) ! correction for stretched domain
!        localtime = simtime - ( real ( localtstep + 1, kind=C_DOUBLE ) / F_DEFORM_S )
!       !> In order to obtain the correct shear rate, time needs to
!       !! be negative after the remesh
!       localtime = localtime - HALF

      end if IFFIRSTREMESH
!      select case ( rk3stage )
!      !-------------------------------
!        !> Runge-Kutta stage 1
!        case ( 1 ) 
!          rkbeta = ZERO
!          rkgamma = THIRD

!          localtime = ( ONE / THREE ) * rk3dt
!        !-------------------------------
!        !> Runge-Kutta stage 2
!        case ( 2 )
!          rkbeta = - FIVE_OVER_NINE 
!          rkgamma = FIFTEEN_OVER_SIXTEEN
!          localtime = ( THREE / FOUR ) * rk3dt 
!        !-------------------------------
!        !> Runge-Kutta stage 3
!        case ( 3 ) 
!          rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
!          rkgamma = EIGHT_OVER_FIFTEEN
!          localtime = rk3dt 
        !-------------------------------
!      end select
!
    else
!
      localtime = simtime
!
    end if IFREMESH


!> Update shear component of transformation tensor
!    rk3mat(1,3) = ( rkbeta * rk3mat(1,3) ) - ( bmat(1,1) * amat(1,3) )
!    bmat(1,3)   = bmat(1,3) + ( rkgamma * rk3dt * rk3mat(1,3) )

    bmat(1,3) = - bmat(1,1) * F_DEFORM_S * localtime
    binv(1,3) =               F_DEFORM_S * localtime / bmat(3,3)

!    write(*,*) 'bmat: ', bmat(1,3), 'binv: ', binv(1,3)

!    write(*,*) 'tstep: ', tstep, ', a/b = localtstep: ', localtstep, ', local time: ', &
!               F_DEFORM_S * localtime, ' bmat, binv: ', bmat(1,3), binv(1,3), 'a: ', a, 'b: ', b

  end if IFFTYPE


end subroutine g_RK3UpdateB


!---------------------------------------------------------------------------------
! subroutine g_RK3Stage1
!> Coordinate stage 1 of the Runge-Kutta routines
!> @param ierr should return 0

subroutine g_RK3Stage1 ( tstep, time, sim_dt, ierr )

  use f_rk3, only: f_RK3Stage1

  implicit none

  integer,intent(in)                   :: tstep
  real(kind=C_DOUBLE),intent(in)       :: time
  real(kind=C_DOUBLE),intent(inout)    :: sim_dt
  PetscErrorCode,intent(inout)                :: ierr

!> insert timing code here

!> \f$ U_1 = U_n \f$

!> \f$ G_1 = F(U_n, t_n) \f$

!> \f$ U_2 = U_1 + \frac{1}{3} \Delta t G_1 \f$

!> Call the fluid Runge-Kutta routines. If there are particles, the particle
!! routines are called between the Fourier transforms in the fluid routines
!! to save memory.
  call f_RK3Stage1 ( tstep, time, sim_dt, ierr )

!> insert timing code here
     

end subroutine g_RK3Stage1

!---------------------------------------------------------------------------------
! subroutine g_RK3Stage2
!> Coordinate stage 2 of the Runge-Kutta routines
!> @param ierr should return 0

subroutine g_RK3Stage2 ( tstep, time, sim_dt, ierr )

  use f_rk3, only: f_RK3Stage2

  implicit none

  integer,intent(in)                   :: tstep
  real(kind=C_DOUBLE),intent(in)       :: time
  real(kind=C_DOUBLE),intent(inout)    :: sim_dt
  PetscErrorCode,intent(inout)                :: ierr

!> insert timing code here

!> \f$ G_2 = - \frac{5}{9} G_1 + F(U_1, t_n + \frac{1}{3} \Delta t) \f$

!> \f$ U_3 = U_2 + \frac{15}{16} \Delta t G_2 \f$

!> Call the fluid Runge-Kutta routines. If there are particles, the particle
!! routines are called between the Fourier transforms in the fluid routines
!! to save memory.
  call f_RK3Stage2 ( tstep, time, sim_dt, ierr )

!> insert timing code here
     

end subroutine g_RK3Stage2

!---------------------------------------------------------------------------------
! subroutine g_RK3Stage3
!> Coordinate stage 3 of the Runge-Kutta routines
!> @param ierr should return 0

subroutine g_RK3Stage3 ( tstep, time, sim_dt, ierr )

  use g_parameters, only : P_TRACK_PART, P_TS_RELEASE, YES
  use f_rk3,        only : f_RK3Stage3
!  use p_control,    only : p_ControlStats
  use p_stats,      only : p_StatsPrimary

  implicit none

  integer,intent(in)                   :: tstep
  real(kind=C_DOUBLE),intent(in)       :: time
  real(kind=C_DOUBLE),intent(inout)    :: sim_dt
  PetscErrorCode,intent(inout)                :: ierr

!> insert timing code here

!> \f$ G_3 = - \frac{153}{128} G_1 + F(U_3, t_n + \frac{3}{4} \Delta t) \f$

!> \f$ U_{n+1} = U_2 + \frac{8}{15} \Delta t G_3 \f$

!> Call the fluid Runge-Kutta routines. If there are particles, the particle
!! routines are called between the Fourier transforms in the fluid routines
!! to save memory.
  call f_RK3Stage3 ( tstep, time, sim_dt, ierr )

!> insert timing code here
     
!  IFPARTRK3: if ( P_TRACK_PART == YES ) then
!    IFRELEASEDRK3: if ( tstep >= P_TS_RELEASE ) then

!      !> Do particle statistics 
!      call p_ControlStats ( tstep, ierr )
!      call p_StatsPrimary ( tstep, ierr )

!    end if IFRELEASEDRK3
!  end if IFPARTRK3

end subroutine g_RK3Stage3
!------------------------------------------------------------------

end module g_rk3
 
