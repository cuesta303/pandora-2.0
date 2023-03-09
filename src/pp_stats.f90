!---------------------------------------------------------------------------------
! module p_stats
!> Contains routines for particle statistics 
module p_stats

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>


private

!----------------------------------------------------------------------
real(kind=C_DOUBLE),save                :: u1u1
real(kind=C_DOUBLE),save                :: u1u2
real(kind=C_DOUBLE),save                :: u1u3
real(kind=C_DOUBLE),save                :: u2u2
real(kind=C_DOUBLE),save                :: u2u3
real(kind=C_DOUBLE),save                :: u3u3
real(kind=C_DOUBLE),save                :: rs_u1u1
real(kind=C_DOUBLE),save                :: rs_u1u2
real(kind=C_DOUBLE),save                :: rs_u1u3
real(kind=C_DOUBLE),save                :: rs_u2u2
real(kind=C_DOUBLE),save                :: rs_u2u3
real(kind=C_DOUBLE),save                :: rs_u3u3
real(kind=C_DOUBLE),save                :: tke
real(kind=C_DOUBLE),save                :: pken
real(kind=C_DOUBLE),save                :: lpken
real(kind=C_DOUBLE),save                :: gpken
real(kind=C_DOUBLE),save                :: vrel_abs
real(kind=C_DOUBLE),save                :: gvrel_abs
real(kind=C_DOUBLE),save                :: vrel1
real(kind=C_DOUBLE),save                :: gvrel1
real(kind=C_DOUBLE),save                :: vrel2
real(kind=C_DOUBLE),save                :: gvrel2
real(kind=C_DOUBLE),save                :: vrel3
real(kind=C_DOUBLE),save                :: gvrel3
real(kind=C_DOUBLE),save                :: rep
real(kind=C_DOUBLE),save                :: grep
real(kind=C_DOUBLE),save                :: stk 
real(kind=C_DOUBLE),save                :: ste 
real(kind=C_DOUBLE),save,dimension(1:3) :: lxrel
real(kind=C_DOUBLE),save,dimension(1:3) :: lx2rel
real(kind=C_DOUBLE),save,dimension(1:3) :: gxrel
real(kind=C_DOUBLE),save,dimension(1:3) :: gx2rel
real(kind=C_DOUBLE),save,dimension(1:3) :: xrel
real(kind=C_DOUBLE),save                :: lx1x2
real(kind=C_DOUBLE),save                :: lx1x3
real(kind=C_DOUBLE),save                :: lx2x3
real(kind=C_DOUBLE),save                :: gx1x2
real(kind=C_DOUBLE),save                :: gx1x3
real(kind=C_DOUBLE),save                :: gx2x3
real(kind=C_DOUBLE),save                :: x1x2
real(kind=C_DOUBLE),save                :: x1x3
real(kind=C_DOUBLE),save                :: x2x3
real(kind=C_DOUBLE),save,dimension(1:3) :: xrelold
real(kind=C_DOUBLE),save,dimension(1:3) :: dfp 
real(kind=C_DOUBLE),save,dimension(1:3) :: vt 
real(kind=C_DOUBLE),save,dimension(1:3) :: vs 

public :: p_StatsPrimary
public :: p_StatsRelativeVelocity
public :: p_StatsWrite

contains
!---------------------------------------------------------------------------------
! subroutine p_StatsZero
!> Initialise particle  statistics.
!> @param ierr should return 0
subroutine p_StatsZero ( ierr )

  use p_arrays,     only : particle_dp, particle_up
  use p_control,    only : np_centre 
  use p_rk3,        only : lvrel_abs, lvrel1, lvrel2, lvrel3, lrep
  use g_parameters, only : MYRANK, MASTER, P_NP_GLOBAL, P_RHO
  use g_constants,  only : ZERO, HALF, ONEPI, SIX

  implicit none

  PetscErrorCode,intent(inout)            :: ierr
  
  u1u1      = ZERO 
  u2u2      = ZERO 
  u3u3      = ZERO 
  u1u2      = ZERO
  u1u3      = ZERO
  u2u3      = ZERO
  lpken     = ZERO
  lvrel_abs = ZERO
  lvrel1    = ZERO 
  lvrel2    = ZERO 
  lvrel3    = ZERO 
  lxrel     = ZERO
  lx2rel    = ZERO
  gxrel     = ZERO
  gx2rel    = ZERO
  dfp       = ZERO
  lrep      = ZERO
  lx1x2     = ZERO
  lx1x3     = ZERO
  lx2x3     = ZERO
 
end subroutine p_StatsZero
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_StatsPrimary
!> Compute particle  statistics.
!> @param ierr should return 0
subroutine p_StatsPrimary ( tstep, dt, ierr )

  use p_arrays,     only : particle_dp, particle_up
  use p_control,    only : np_centre 
  use p_rk3,        only : lvrel_abs, lvrel1, lvrel2, lvrel3, upgrav
  use g_parameters, only : MYRANK, MASTER, P_NP_GLOBAL, P_RHO, F_RHO, &
                           P_KEEP_INITIAL_POSITION, YES, P_DIA_VALS, F_NU, &
                           P_FOLLOW
  use g_constants,  only : ZERO, HALF, ONEPI, SIX, EIGHTEEN
  use f_fluidstats, only : tn_kol, tn_edy

  implicit none

  integer,intent(in)                 :: tstep
  real(kind=C_DOUBLE),intent(in)     :: dt
  PetscErrorCode,intent(inout)              :: ierr
  integer                            :: np
  real(kind=C_DOUBLE)                :: u1, u2, u3
  real(kind=C_DOUBLE)                :: dp
  real(kind=C_DOUBLE),dimension(1:3) :: xp
  real(kind=C_DOUBLE),dimension(1:3) :: xp0
  real(kind=C_DOUBLE)                :: taup

  call p_StatsZero ( ierr ) 

  DOPART: do np = 1, np_centre, 1

    dp = particle_dp(1,np)

    !> Read local velocity
    u1 = particle_up(1,np)
    u2 = particle_up(2,np)
    u3 = particle_up(3,np)

!    write(*,*) 'Velocity for stats: ', u1, u2, u3

    !> Compute Reynolds stresses
    u1u1 = u1u1 + u1 * u1
    u2u2 = u2u2 + u2 * u2
    u3u3 = u3u3 + u3 * u3
    u1u2 = u1u2 + u1 * u2
    u1u3 = u1u3 + u1 * u3
    u2u3 = u2u3 + u2 * u3

    !> Particle kinetic energy
    lpken = lpken + ( HALF * ((ONEPI / SIX) * (dp**3) * P_RHO) * \
            ( u1**2 + u2**2 + u3**2 ) )

    IFX0: if ( P_KEEP_INITIAL_POSITION == YES ) then

      call p_StatsDispersion ( np, ierr )

    end if IFX0

  end do DOPART

  !> Determine Reynolds stresses
  call MPI_Allreduce ( u1u1, rs_u1u1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u2, rs_u2u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u3u3, rs_u3u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u2, rs_u1u2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u1u3, rs_u1u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  call MPI_Allreduce ( u2u3, rs_u2u3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)

!  write(*,*), ' u1u1, rs_u1u1, P_NP_GLOBAL: ', u1u1, rs_u1u1, P_NP_GLOBAL

  rs_u1u1 = rs_u1u1 / real ( P_NP_GLOBAL, kind=C_DOUBLE )
  rs_u2u2 = rs_u2u2 / real ( P_NP_GLOBAL, kind=C_DOUBLE )
  rs_u3u3 = rs_u3u3 / real ( P_NP_GLOBAL, kind=C_DOUBLE )
  rs_u1u2 = rs_u1u2 / real ( P_NP_GLOBAL, kind=C_DOUBLE )
  rs_u1u3 = rs_u1u3 / real ( P_NP_GLOBAL, kind=C_DOUBLE )
  rs_u2u3 = rs_u2u3 / real ( P_NP_GLOBAL, kind=C_DOUBLE )

  !> Compute turbulent kinetic energy
  tke = HALF * ( rs_u1u1 + rs_u2u2 + rs_u3u3 )

  !> Particle kinetic energy
  call MPI_Allreduce ( lpken, gpken, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr)
  pken = gpken / REAL ( P_NP_GLOBAL, kind=C_DOUBLE )

  !> Compute particle response time
! do i = 1,DIA_LVLS
  taup = (P_DIA_VALS(1)**2) * (P_RHO/F_RHO) / (EIGHTEEN * F_NU)
  !> Compute Stokes number
  stk = taup / tn_kol 
  ste = taup / tn_edy 
! end do

  IFDISPERSION: if ( P_KEEP_INITIAL_POSITION == YES ) then
  !> Compute relative displacement

    !> Keep old displacement for diffusivity
    xrelold = xrel

    call MPI_Allreduce ( lxrel, gxrel, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
    call MPI_Allreduce ( lx2rel, gx2rel, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
    call MPI_Allreduce ( lx1x2, gx1x2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
    call MPI_Allreduce ( lx1x3, gx1x3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
    call MPI_Allreduce ( lx2x3, gx2x3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
!    gxrel = ( gxrel / real(P_NP_GLOBAL,kind=C_DOUBLE) )**2
!    gx2rel = gx2rel / real(P_NP_GLOBAL,kind=C_DOUBLE)
!    xrel = gx2rel - ( gxrel**2 )
    xrel = gx2rel / real(P_NP_GLOBAL,kind=C_DOUBLE)
    x1x2 = gx1x2 / real(P_NP_GLOBAL,kind=C_DOUBLE)
    x1x3 = gx1x3 / real(P_NP_GLOBAL,kind=C_DOUBLE)
    x2x3 = gx2x3 / real(P_NP_GLOBAL,kind=C_DOUBLE)

    !> Compute particle diffusivity
    dfp(1) = HALF * ( xrel(1) - xrelold(1) ) / dt
    dfp(2) = HALF * ( xrel(2) - xrelold(2) ) / dt
    dfp(3) = HALF * ( xrel(3) - xrelold(3) ) / dt

  end if IFDISPERSION

  !> Velocity in quiescent fluid
  vs = upgrav

end subroutine p_StatsPrimary
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_StatsRelativeVelocity
!> Compute particle relative velocity.
!> @param ierr should return 0
subroutine p_StatsRelativeVelocity ( ierr )

  use p_rk3,        only : lvrel_abs, lvrel1, lvrel2, lvrel3, lrep, vabs
  use g_constants,  only : ZERO, EIGHTEEN
  use g_parameters, only : MYRANK, MASTER, P_NP_GLOBAL, P_RHO, &
                           P_KEEP_INITIAL_POSITION, YES, &
                           F_RHO, P_DIA_VALS, F_NU, GRAVITY

  implicit none

  PetscErrorCode,intent(inout)              :: ierr

  real(kind=C_DOUBLE)                :: taup
  real(kind=C_DOUBLE)                :: g

  !> Compute relative velocity
  call MPI_Allreduce (lvrel_abs, gvrel_abs, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  call MPI_Allreduce (lvrel1, gvrel1, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  call MPI_Allreduce (lvrel2, gvrel2, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  call MPI_Allreduce (lvrel3, gvrel3, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )

!real(kind=C_DOUBLE),save,public,dimension(1:3)  :: lvabs 

  vrel_abs = sqrt(gvrel_abs / real(P_NP_GLOBAL,kind=C_DOUBLE))
  vrel1    = gvrel1 / real(P_NP_GLOBAL,kind=C_DOUBLE)
  vrel2    = gvrel2 / real(P_NP_GLOBAL,kind=C_DOUBLE)
  vrel3    = gvrel3 / real(P_NP_GLOBAL,kind=C_DOUBLE)

  !> Compute particle Re
  call  MPI_Allreduce ( lrep, grep, 1, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  rep = grep / REAL(P_NP_GLOBAL,kind=C_DOUBLE)

  !> Get settling velocity (need to copy value, because vabs will be updated
  !! during RK3
  vt = vabs

end subroutine p_StatsRelativeVelocity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_StatsDispersion
!> Compute dispersion.
!> @param ierr should return 0
subroutine p_StatsDispersion ( np, ierr )

  use p_arrays,     only : particle_xp, particle_xp0
  use p_control,    only : np_centre 
  use p_rk3,        only : lvrel_abs, lvrel1, lvrel2, lvrel3, xpgrav
  use g_parameters, only : MYRANK, MASTER, P_NP_GLOBAL, P_RHO, YES, &
                           P_KEEP_INITIAL_POSITION, F_TYPE, F_ISOTROPIC, &
                           P_GRAVITY, GRAVITY
  use g_constants,  only : ZERO, HALF, ONEPI, SIX
  use g_domain,     only : bvalexist, binv

  implicit none

  integer,intent(in)                 :: np
  PetscErrorCode,intent(inout)              :: ierr

  real(kind=C_DOUBLE),dimension(1:3) :: xp
  real(kind=C_DOUBLE),dimension(1:3) :: xp0
  real(kind=C_DOUBLE),dimension(1:3) :: xrel_moving
  real(kind=C_DOUBLE),dimension(1:3) :: xrel_local

  integer                            :: m, n
  integer                            :: i

  xp(:)  = particle_xp(:,np)
  xp0(:) = particle_xp0(:,np)

  xrel_moving(:) = xp(:) - xp0(:) 

  IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
    !> Distance is the same
    xrel_local(:)  = xrel_moving(:)
  else
    !> Need to transform into laboratory system
    xrel_local = ZERO
    DOM: do m = 1, 3, 1
      DON: do n = 1, 3, 1
        IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
          xrel_local(m) = xrel_local(m) + binv(m,n) * xrel_moving(n)
        end if IFBEXIST
      end do DON
    end do DOM
  end if IFISOTROPIC

!  write(*,*) ' Distance moving: ', xrel_moving, &
!             ' , distance lab: ', xrel_local

!  write(*,*) ' xp(:), xp0(:), ', xp(:), xp0(:), &
!             ' distance moving: ', xrel_moving, &
!             ' , distance lab: ', xrel_local, &
!             ' , distance gravity: ', xpgrav, xpgrav_local


!  IFGRAVITY: if ( P_GRAVITY == YES ) then
!    xrel_local = xrel_local - xpgrav
!  end if IFGRAVITY

  lxrel(:) = lxrel(:) + xrel_local(:)
  lx2rel(:) = lx2rel(:) + ( xrel_local(:)**2 )
  lx1x2 = lx1x2 + ( xrel_local(1) * xrel_local(2) )
  lx1x3 = lx1x3 + ( xrel_local(1) * xrel_local(3) )
  lx2x3 = lx2x3 + ( xrel_local(2) * xrel_local(3) )

!real(kind=C_DOUBLE),save,public,dimension(1:3)  :: xpgrav
end subroutine p_StatsDispersion
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
subroutine p_StatsWrite ( tstep, ierr )

  use g_parameters, only : MYRANK, MASTER, P_FOLLOW
  use g_files,      only : FUNIT_PSTATS, FUNIT_PSINGL

  use p_rk3,        only : xpsingle, upsingle, xp0single

  implicit none

  integer,intent(in)                  :: tstep
  PetscErrorCode,intent(inout)               :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(FUNIT_PSTATS,50) tstep, rs_u1u1, rs_u2u2, rs_u3u3, rs_u1u2, rs_u1u3, rs_u2u3, tke, &
                           pken, vrel1, vrel2, vrel3, vrel_abs, xrel, dfp, x1x2, x1x3, x2x3, &
                           vt, vs, rep, stk, ste

  end if IFMASTER


50  FORMAT('ts: ',i6,' u1u1: ',e14.6,' u2u2: ',e14.6,' u3u3: ',e14.6, &
           ' u1u2: ',e14.6,' u1u3: ',e14.6,' u2u3: ',e14.6, ' tke: ',e14.6, &
           ' ptke: ',e14.6,' vrel1: ',e14.6,' vrel2: ',e14.6,' vrel3: ',e14.6, &
           ' vrelabs: ',e14.6,' xrel1: ',e14.6,' xrel2: ',e14.6,' xrel3: ',e14.6, &
           ' dfp1: ',e14.6,' dfp2: ',e14.6,' dfp3: ',e14.6, &
           ' x1x2: ',e14.6,' x1x3: ',e14.6,' x2x3: ',e14.6, &
           ' vt: ', 3e14.6, ' vs: ', 3e14.6, ' Rep: ',e14.6,' Stk: ',e14.6,' Ste: ',e14.6)

end subroutine p_StatsWrite
!---------------------------------------------------------------------------------


!----------------------------------------------------------------
!subroutine p_quants(ierr)

!use d_arrays, only : length_x,length_y,length_z

!implicit none

!integer, intent(inout) :: ierr

!    eul_nd = real(np_global,kind=C_DOUBLE) /  (length_x*length_y*length_z)

!end subroutine
!---------------------------------------------------------------------------------

end module p_stats
