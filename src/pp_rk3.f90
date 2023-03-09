module p_rk3

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

!real(kind=C_DOUBLE),save,public :: vmin

!real(kind=C_DOUBLE),save,public :: sfdrag(3) = 0.0
!real(kind=C_DOUBLE),save,public :: sfelec(3) = 0.0
!real(kind=C_DOUBLE),save,public :: sfgrav(3) = 0.0
!real(kind=C_DOUBLE),save,public :: fdrag = 0.0
!real(kind=C_DOUBLE),save,public :: felec = 0.0
!real(kind=C_DOUBLE),save,public :: fgrav = 0.0
real(kind=C_DOUBLE),save,public                 :: lvrel_abs
real(kind=C_DOUBLE),save,public                 :: lvrel1
real(kind=C_DOUBLE),save,public                 :: lvrel2
real(kind=C_DOUBLE),save,public                 :: lvrel3
real(kind=C_DOUBLE),save,public                 :: lrep
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: xpgrav
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: upgrav
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: rk3xpgrav
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: rk3upgrav
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: lvabs 
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: gvabs 
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: vabs 
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: xpsingle 
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: upsingle
real(kind=C_DOUBLE),save,public,dimension(1:3)  :: xp0single

!real(kind=C_DOUBLE),save :: drag_accln(3), elec_accln(3), grav_accln(3)
!----------------------------------------------------------------------

public :: p_RK3SolveParticle
public :: p_RK3SolveParticleBF
public :: p_RK3SolveParticleR
public :: p_RK3AbsoluteVelocity
public :: p_RK3AbsoluteVelocityBF
public :: p_RK3AbsoluteVelocityR
public :: p_RK3Mov2Lab
public :: p_RK3Lab2Mov
public :: p_RK3UAbs

contains
!---------------------------------------------------------------------------------
subroutine p_RK3SolveParticle ( npin, stage, deltat, ierr )    

  use g_parameters, only : YES, P_TWO_WAY, P_GRAVITY, P_FOLLOW
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use p_arrays,     only : particle_xp, particle_up, particle_rk3, particle_up_rk3, particle_dp
  use p_interp,     only : p_InterpInterpolate3D
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  real(kind=C_DOUBLE)                        :: rkbeta
  real(kind=C_DOUBLE)                        :: rkgamma
  real(kind=C_DOUBLE)                        :: dp
  integer                                    :: np
!  real(kind=C_DOUBLE),dimension(1:3)         :: upold
  real(kind=C_DOUBLE),dimension(1:3)         :: particleuf
  real(kind=C_DOUBLE),dimension(1:3)         :: accln 
  real(kind=C_DOUBLE),dimension(1:3)         :: xrk3 

  !> Nothing to do if there are no particles
  IFPARTEXIST: if ( npin > 0 ) then
    IFSTAGE: if ( stage == 1 ) then
      rkbeta = ZERO
      rkgamma = THIRD

    else if (stage == 2 ) then
      rkbeta = - FIVE_OVER_NINE 
      rkgamma = FIFTEEN_OVER_SIXTEEN

    else if ( stage == 3 ) then
      rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
      rkgamma = EIGHT_OVER_FIFTEEN
    end if IFSTAGE

    DOPART: do np = 1, npin, 1

      dp = particle_dp(1,np)

      call p_InterpInterpolate3D ( particle_xp(:,np), particleuf, ierr )

    !> If inertial particle, solve Maxey-Riley equation for velocity
      IFINERTIAL: if ( dp > ZERO ) then

        !> If first stage compute relative velocity before updating
        IFSTAGE1RELVEL: if ( stage == 1 ) then
          
          call p_RK3RelativeVelocity ( particleuf, particle_up(:,np), particle_dp(1,np), ierr )

        end if IFSTAGE1RELVEL

        accln = p_RK3ComputeForce ( particleuf, particle_up(:,np), particle_dp(1,np) )

        particle_rk3(:,np) = ( rkbeta * particle_rk3(:,np) ) + accln(:)

        particle_up(:,np)  = particle_up(:,np) + ( rkgamma * deltat * particle_rk3(:,np))

      !> If fluid particle, the particle velocity is just the fluid velocity
      else

        !> Particle velocity is just the fluid velocity at this position
        particle_up(:,np)  = particleuf(:)

        !> Log velocity for quality control
        particleuf = ZERO
        call p_RK3RelativeVelocity ( particleuf, particle_up(:,np), particle_dp(1,np), ierr )

      end if IFINERTIAL

      !> Intermediate stage of position RK3
      particle_up_rk3(:,np) = ( rkbeta * particle_up_rk3(:,np) ) + particle_up(:,np)

      !> Compute new particle position
      particle_xp(:,np)  = particle_xp(:,np) + ( rkgamma * deltat * particle_up_rk3(:,np) )
            
  !    mp = ( onepi / six ) * (dp**3) * particle_rho

  !    if ( P_TWO_WAY == YES ) then
  !      lfdrag = lfdrag + sqrt( (mp*drag_accln(1))**2 + (mp*drag_accln(2))**2 + (mp*drag_accln(3))**2 )
  !      lfgrav = lfgrav + sqrt( (mp*grav_accln(1))**2 + (mp*grav_accln(2))**2 + (mp*grav_accln(3))**2 )
  !      lfelec = lfelec + sqrt( (mp*elec_accln(1))**2 + (mp*elec_accln(2))**2 + (mp*elec_accln(3))**2 )
           
  !      xp = part_data(np,1:3)
  !      call interp_find(xp,ii,jj,kk)

  !      force_partr(ii,jj,kk,:) = force_partr(ii,jj,kk,:) + (((-1.0 * accln(:))*particle_rho)*vol_particle(np))
  !    endif
    end do DOPART


  end if IFPARTEXIST

  !> If gravity, compute the absolute displacement of particles
  !! without turbulence for statistics
  IFGRAVITY: if ( P_GRAVITY == YES ) then
    call p_RK3AbsoluteDisplacement ( stage, deltat, ierr )    
  end if IFGRAVITY

  !> If single particle tracking on, evaluate quantities in laboratory system
  !! now.
!  IFFOLLOW: if ( P_FOLLOW > 0 ) then
!    call p_RK3Follow ( npin, stage, deltat, ierr )
!  end if IFFOLLOW

end subroutine p_RK3SolveParticle
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_RK3SolveParticleBF ( npin, stage, deltat, ierr )    

  use g_parameters, only : YES, P_TWO_WAY, P_GRAVITY, P_FOLLOW, P_KEEP_INITIAL_POSITION
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use p_arrays,     only : particle_xp, particle_up, particle_rk3, particle_up_rk3, particle_dp
  use p_interp,     only : p_InterpInterpolate3D
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  real(kind=C_DOUBLE)                        :: rkbeta
  real(kind=C_DOUBLE)                        :: rkgamma
  real(kind=C_DOUBLE)                        :: dp
  integer                                    :: np
!  real(kind=C_DOUBLE),dimension(1:3)         :: upold
  real(kind=C_DOUBLE),dimension(1:3)         :: ufmov
  real(kind=C_DOUBLE),dimension(1:3)         :: uflab
  real(kind=C_DOUBLE),dimension(1:3)         :: particleuf
  real(kind=C_DOUBLE),dimension(1:3)         :: xpmov
  real(kind=C_DOUBLE),dimension(1:3)         :: xplab
  real(kind=C_DOUBLE),dimension(1:3)         :: ufabs
  real(kind=C_DOUBLE),dimension(1:3)         :: accln 
  real(kind=C_DOUBLE),dimension(1:3)         :: xrk3 

  !> Nothing to do if there are no particles
  IFPARTEXIST: if ( npin > 0 ) then

    !> If gravity, compute the absolute displacement of particles
    !! without turbulence for statistics
    IFGRAVITY: if ( P_GRAVITY == YES ) then
      call p_RK3AbsoluteDisplacement(stage, deltat, ierr )    

          IFINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then 
            call p_RK3CorrectInitialPosition(npin, stage, deltat, ierr )
          end if IFINITIAL

    end if IFGRAVITY

    IFSTAGE: if ( stage == 1 ) then
      rkbeta = ZERO
      rkgamma = THIRD

    else if (stage == 2 ) then
      rkbeta = - FIVE_OVER_NINE 
      rkgamma = FIFTEEN_OVER_SIXTEEN

    else if ( stage == 3 ) then
      rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
      rkgamma = EIGHT_OVER_FIFTEEN
    end if IFSTAGE

    DOPART: do np = 1, npin, 1

      dp = particle_dp(1,np)
      xpmov = particle_xp(:,np)

      !> Find fluid velocity at particle position in moving system.
      call p_InterpInterpolate3D ( xpmov, ufmov, ierr )

      !> Transform coordinate to laboratory system
      call p_RK3Mov2Lab ( xpmov, xplab, ierr )

!      write(*,*) ' Old position: ', xpmov, xplab

      !> Transform fluid velocity to laboratory system
      call p_RK3Mov2Lab ( ufmov, uflab, ierr )

!      !> Find absolute velocity 
!      call p_RK3UAbs ( xplab, ufabs, ierr )

      !> Find total velocity
      particleuf = uflab !+ ufabs 

!      write(*,*) ' Velocity: ', ufmov, uflab, ufabs, particleuf

      !> If inertial particle, solve Maxey-Riley equation for velocity
      IFINERTIAL: if ( dp > ZERO ) then

        !> If first stage compute relative velocity before updating
        IFSTAGE1RELVEL: if ( stage == 1 ) then
          
          call p_RK3RelativeVelocity ( particleuf, particle_up(:,np), particle_dp(1,np), ierr )

        end if IFSTAGE1RELVEL


        accln = p_RK3ComputeForce ( particleuf, particle_up(:,np), particle_dp(1,np) ) - &
                !> Find homogeneous source term of acceleration
                p_RK3ComputeHomogeneousSourceTerm ( particle_up(:,np) )

        particle_rk3(:,np) = ( rkbeta * particle_rk3(:,np) ) + accln(:)

        particle_up(:,np)  = particle_up(:,np) + ( rkgamma * deltat * particle_rk3(:,np))

      !> If fluid particle, the particle velocity is just the fluid velocity
      else

        !> Particle velocity is just the fluid velocity at this position
        particle_up(:,np)  = particleuf(:)

        !> Log velocity for quality control
        particleuf = ZERO
        call p_RK3RelativeVelocity ( particleuf, particle_up(:,np), particle_dp(1,np), ierr )

      end if IFINERTIAL

      !> Intermediate stage of position RK3
      particle_up_rk3(:,np) = ( rkbeta * particle_up_rk3(:,np) ) + particle_up(:,np)

      !> Compute new particle position in laboratory system
      xplab(:)  = xplab(:) + ( rkgamma * deltat * particle_up_rk3(:,np) )
    
      !> Transform particle position to moving system
      call p_RK3Lab2Mov( xplab, xpmov, ierr ) 
      particle_xp(:,np) = xpmov
        
!      write(*,*) ' New position: ', xpmov, xplab

    end do DOPART


  end if IFPARTEXIST

  !> If single particle tracking on, evaluate quantities in laboratory system
  !! now.
!  IFFOLLOW: if ( P_FOLLOW > 0 ) then
!    call p_RK3Follow ( npin, stage, deltat, ierr )
!  end if IFFOLLOW

end subroutine p_RK3SolveParticleBF
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_RK3SolveParticleR ( npin, stage, deltat, ierr )    

  use g_parameters, only : YES, P_TWO_WAY, P_GRAVITY, P_FOLLOW
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use p_arrays,     only : particle_xp, particle_up, particle_rk3, particle_up_rk3, particle_dp
  use p_interp,     only : p_InterpInterpolate3D
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  real(kind=C_DOUBLE)                        :: rkbeta
  real(kind=C_DOUBLE)                        :: rkgamma
  real(kind=C_DOUBLE)                        :: dp
  integer                                    :: np
!  real(kind=C_DOUBLE),dimension(1:3)         :: upold
  real(kind=C_DOUBLE),dimension(1:3)         :: particleuf
  real(kind=C_DOUBLE),dimension(1:3)         :: accln 
  real(kind=C_DOUBLE),dimension(1:3)         :: xrk3 

  !> Nothing to do if there are no particles
  IFPARTEXIST: if ( npin > 0 ) then
    IFSTAGE: if ( stage == 1 ) then
      rkbeta = ZERO
      rkgamma = THIRD

    else if (stage == 2 ) then
      rkbeta = - FIVE_OVER_NINE 
      rkgamma = FIFTEEN_OVER_SIXTEEN

    else if ( stage == 3 ) then
      rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
      rkgamma = EIGHT_OVER_FIFTEEN
    end if IFSTAGE

    DOPART: do np = 1, npin, 1

      dp = particle_dp(1,np)

      call p_InterpInterpolate3D ( particle_xp(:,np), particleuf, ierr )

    !> If inertial particle, solve Maxey-Riley equation for velocity
      IFINERTIAL: if ( dp > ZERO ) then

        !> If first stage compute relative velocity before updating
        IFSTAGE1RELVEL: if ( stage == 1 ) then
          
          call p_RK3RelativeVelocity ( particleuf, particle_up(:,np), particle_dp(1,np), ierr )

        end if IFSTAGE1RELVEL

        accln = p_RK3ComputeForce ( particleuf, particle_up(:,np), particle_dp(1,np) )

        particle_rk3(:,np) = ( rkbeta * particle_rk3(:,np) ) + accln(:)

        particle_up(:,np)  = particle_up(:,np) + ( rkgamma * deltat * particle_rk3(:,np))

      !> If fluid particle, the particle velocity is just the fluid velocity
      else

        !> Particle velocity is just the fluid velocity at this position
        particle_up(:,np)  = particleuf(:)

        !> Log velocity for quality control
        particleuf = ZERO
        call p_RK3RelativeVelocity ( particleuf, particle_up(:,np), particle_dp(1,np), ierr )

      end if IFINERTIAL

      !> Intermediate stage of position RK3
      particle_up_rk3(:,np) = ( rkbeta * particle_up_rk3(:,np) ) + particle_up(:,np)

      !> Compute new particle position
      particle_xp(:,np)  = particle_xp(:,np) + ( rkgamma * deltat * particle_up_rk3(:,np) )
            
    end do DOPART


  end if IFPARTEXIST

  !> If gravity, compute the absolute displacement of particles
  !! without turbulence for statistics
  IFGRAVITY: if ( P_GRAVITY == YES ) then
    call p_RK3AbsoluteDisplacement ( stage, deltat, ierr )    
  end if IFGRAVITY

  !> If single particle tracking on, evaluate quantities in laboratory system
  !! now.
!  IFFOLLOW: if ( P_FOLLOW > 0 ) then
!    call p_RK3Follow ( npin, stage, deltat, ierr )
!  end if IFFOLLOW

end subroutine p_RK3SolveParticleR
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_RK3RelativeVelocity ( uf, up, dp, ierr )

  use g_constants,  only : ZERO
  use g_parameters, only : F_NU, F_TYPE, F_ISOTROPIC, P_GRAVITY, YES
  use g_domain,     only : bvalexist, binv

  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(in)      :: uf
  real(kind=C_DOUBLE),dimension(1:3),intent(in)      :: up
  real(kind=C_DOUBLE),intent(in)                     :: dp

  PetscErrorCode,intent(inout)                              :: ierr

  real(kind=C_DOUBLE),dimension(1:3)                 :: vrel_moving
  real(kind=C_DOUBLE),dimension(1:3)                 :: vrel_local
  real(kind=C_DOUBLE)                                :: vrel_abs

  integer                                            :: m, n

!  vrel_moving(:) = uf(:) - up(:)

!  IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
    !> Velocity is the same
!    vrel_local(:)  = vrel_moving

   vrel_local(:) = up(:) - uf(:)
!  else
!    !> Need to transform into laboratory system
!    vrel_local = ZERO
!    DOM: do m = 1, 3, 1
!      DON: do n = 1, 3, 1
!        IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
!          vrel_local(m) = vrel_local(m) + binv(m,n) * vrel_moving(n)
!        end if IFBEXIST
!      end do DON
!    end do DOM
!  end if IFISOTROPIC

!  write(*,*) ' Velocity moving: ', vrel_moving, &
!             ' , velocity lab: ', vrel_local

  IFGRAVITY: if ( P_GRAVITY == YES ) then
    vrel_local = vrel_local - upgrav
  end if IFGRAVITY

!  write(*,*) ' vrel_local, upgrav ', vrel_local, upgrav

  vrel_abs = (vrel_local(1)**2) + (vrel_local(2)**2) + (vrel_local(3)**2)
  lvrel_abs = lvrel_abs + vrel_abs
  lvrel1 = lvrel1 + vrel_local(1)**2 
  lvrel2 = lvrel2 + vrel_local(2)**2 
  lvrel3 = lvrel3 + vrel_local(3)**2 

  !> Compute particle Re
  lrep   = lrep + ( (sqrt(vrel_abs) * dp) / F_NU )

end subroutine p_RK3RelativeVelocity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_RK3AbsoluteDisplacement
!> Apply gravity to initial particle position.  
subroutine p_RK3AbsoluteDisplacement ( stage, deltat, ierr )    

  use g_parameters, only : P_DIA_VALS, GRAVITY, F_TYPE, F_ISOTROPIC
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use g_domain,     only : bvalexist, binv
  use p_arrays,     only : particle_up
 
  implicit none
 
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  real(kind=C_DOUBLE)                        :: rkbeta
  real(kind=C_DOUBLE)                        :: rkgamma

  real(kind=C_DOUBLE),dimension(1:3)         :: drag_accln
  real(kind=C_DOUBLE),dimension(1:3)         :: upold
  real(kind=C_DOUBLE),dimension(1:3)         :: uf

  integer                                    :: np
  real(kind=C_DOUBLE),dimension(1:3)         :: up
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_moving
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_local
  integer                                    :: m, n

  !> Compute movement of a particle with gravity, but no turbulence
  IFSTAGE: if ( stage == 1 ) then
    rkbeta = ZERO
    rkgamma = THIRD

  else if (stage == 2 ) then
    rkbeta = - FIVE_OVER_NINE 
    rkgamma = FIFTEEN_OVER_SIXTEEN

  else if ( stage == 3 ) then
    rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
    rkgamma = EIGHT_OVER_FIFTEEN
  end if IFSTAGE

  !> Need drag force to end up with terminal velocity
  uf = ZERO
  drag_accln = p_RK3SolidForce ( uf, upgrav, P_DIA_VALS(1) )

  !> Velocity
  rk3upgrav = ( rkbeta * rk3upgrav ) + GRAVITY + drag_accln  - &
              !> Find homogeneous source term of acceleration
              p_RK3ComputeHomogeneousSourceTerm ( upgrav )
  upgrav    = upgrav + ( rkgamma * deltat * rk3upgrav )
  rk3xpgrav = ( rkbeta * rk3xpgrav ) + upgrav
  xpgrav    = xpgrav + ( rkgamma * deltat * rk3xpgrav )

!  write(*,*) ' Velocity: ', upgrav, ', position: ', xpgrav

end subroutine p_RK3AbsoluteDisplacement
!!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_RK3CorrectInitialPosition ( npin, stage, deltat, ierr )    

  use g_parameters, only : YES, P_TWO_WAY, P_GRAVITY, P_FOLLOW
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use p_arrays,     only : particle_xp0, particle_dp
  use p_interp,     only : p_InterpInterpolate3D
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  real(kind=C_DOUBLE)                        :: rkbeta
  real(kind=C_DOUBLE)                        :: rkgamma
  real(kind=C_DOUBLE)                        :: dp
  integer                                    :: np
  real(kind=C_DOUBLE),dimension(1:3)         :: xp0mov
  real(kind=C_DOUBLE),dimension(1:3)         :: xp0lab

  !> Nothing to do if there are no particles
  IFPARTEXIST: if ( npin > 0 ) then

    IFSTAGE: if ( stage == 1 ) then
      rkbeta = ZERO
      rkgamma = THIRD

    else if (stage == 2 ) then
      rkbeta = - FIVE_OVER_NINE 
      rkgamma = FIFTEEN_OVER_SIXTEEN

    else if ( stage == 3 ) then
      rkbeta = - ONEFIVETHREE_OVER_ONETWOEIGHT
      rkgamma = EIGHT_OVER_FIFTEEN
    end if IFSTAGE

    DOPART: do np = 1, npin, 1

      dp = particle_dp(1,np)
      xp0mov = particle_xp0(:,np)

      !> Transform coordinate to laboratory system
      call p_RK3Mov2Lab ( xp0mov, xp0lab, ierr )

      !> Only correct initial position of inertial particles
      IFINERTIAL: if ( dp > ZERO ) then

        !> Compute new particle position in laboratory system.
        !! Acceleration is the same for all particles.
        !! @todo This needs to be extended for polydisperse
        !! simulations.
        xp0lab(:)  = xp0lab(:) + ( rkgamma * deltat * rk3xpgrav )
    
        !> Transform particle position to moving system
        call p_RK3Lab2Mov( xp0lab, xp0mov, ierr ) 
        particle_xp0(:,np) = xp0mov
        
      end if IFINERTIAL

    end do DOPART


  end if IFPARTEXIST

end subroutine p_RK3CorrectInitialPosition
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_RK3AbsoluteVelocity
!>  Compute absolute velocity of all particles. 
subroutine p_RK3AbsoluteVelocity ( npin, stage, deltat, ierr )    

  use g_parameters, only : P_DIA_VALS, GRAVITY, F_TYPE, F_ISOTROPIC, &
                           P_NP_GLOBAL, P_FOLLOW
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use g_domain,     only : bvalexist, binv
  use p_arrays,     only : particle_up
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: np
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_moving
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_local
  integer                                    :: m, n

  !> Compute average particle velocity for computing the correct dispersion
  !! and settling velocity
  lvabs = ZERO

  IFPARTEXIST: if ( npin > 0 ) then
    DOPART: do np = 1, npin, 1

      vabs_moving(:) =  particle_up(:,np)

      IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
        vabs_local(:)  = vabs_moving
      else
        !> Need to transform into laboratory system
        vabs_local = ZERO
        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              vabs_local(m) = vabs_local(m) + binv(m,n) * vabs_moving(n)
            end if IFBEXIST
          end do DON
        end do DOM
      end if IFISOTROPIC

      lvabs = lvabs + vabs_local 

    end do DOPART
  end if IFPARTEXIST

  !> Compute absolute velocity
  call MPI_Allreduce (lvabs, gvabs, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  vabs = gvabs / real(P_NP_GLOBAL,kind=C_DOUBLE)

end subroutine p_RK3AbsoluteVelocity
!!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_RK3AbsoluteVelocityBF
!>  Compute absolute velocity of all particles. 
subroutine p_RK3AbsoluteVelocityBF ( npin, stage, deltat, ierr )    

  use g_parameters, only : P_DIA_VALS, GRAVITY, F_TYPE, F_ISOTROPIC, &
                           P_NP_GLOBAL, P_FOLLOW
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use g_domain,     only : bvalexist, binv
  use p_arrays,     only : particle_up
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: np
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_moving
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_local
  integer                                    :: m, n

  !> Compute average particle velocity for computing the correct dispersion
  !! and settling velocity
  lvabs = ZERO

  IFPARTEXIST: if ( npin > 0 ) then
    DOPART: do np = 1, npin, 1

      vabs_moving(:) =  particle_up(:,np)

!      IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
        vabs_local(:)  = vabs_moving
!      else
!        !> Need to transform into laboratory system
!        vabs_local = ZERO
!        DOM: do m = 1, 3, 1
!          DON: do n = 1, 3, 1
!            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
!              vabs_local(m) = vabs_local(m) + binv(m,n) * vabs_moving(n)
!            end if IFBEXIST
!          end do DON
!        end do DOM
!      end if IFISOTROPIC

      lvabs = lvabs + vabs_local 

    end do DOPART
  end if IFPARTEXIST

  !> Compute absolute velocity
  call MPI_Allreduce (lvabs, gvabs, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  vabs = gvabs / real(P_NP_GLOBAL,kind=C_DOUBLE)

end subroutine p_RK3AbsoluteVelocityBF
!!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_RK3AbsoluteVelocityR
!>  Compute absolute velocity of all particles. 
subroutine p_RK3AbsoluteVelocityR ( npin, stage, deltat, ierr )    

  use g_parameters, only : P_DIA_VALS, GRAVITY, F_TYPE, F_ISOTROPIC, &
                           P_NP_GLOBAL, P_FOLLOW
  use g_constants,  only : ZERO, ONEPI, SIX, THIRD, FIVE_OVER_NINE, FIFTEEN_OVER_SIXTEEN, &
                           ONEFIVETHREE_OVER_ONETWOEIGHT, EIGHT_OVER_FIFTEEN
  use g_domain,     only : bvalexist, binv
  use p_arrays,     only : particle_up
 
  implicit none
 
  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: np
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_moving
  real(kind=C_DOUBLE),dimension(1:3)         :: vabs_local
  integer                                    :: m, n

  !> Compute average particle velocity for computing the correct dispersion
  !! and settling velocity
  lvabs = ZERO

  IFPARTEXIST: if ( npin > 0 ) then
    DOPART: do np = 1, npin, 1

      vabs_moving(:) =  particle_up(:,np)

      IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
        vabs_local(:)  = vabs_moving
      else
        !> Need to transform into laboratory system
        vabs_local = ZERO
        DOM: do m = 1, 3, 1
          DON: do n = 1, 3, 1
            IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
              vabs_local(m) = vabs_local(m) + binv(m,n) * vabs_moving(n)
            end if IFBEXIST
          end do DON
        end do DOM
      end if IFISOTROPIC

      lvabs = lvabs + vabs_local 

    end do DOPART
  end if IFPARTEXIST

  !> Compute absolute velocity
  call MPI_Allreduce (lvabs, gvabs, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
  vabs = gvabs / real(P_NP_GLOBAL,kind=C_DOUBLE)

end subroutine p_RK3AbsoluteVelocityR
!!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_RK3Follow
!> Compute particle relative velocity.
!> @param ierr should return 0
subroutine p_RK3Follow ( npin, stage, deltat, ierr )

  use g_constants,  only : ZERO, EIGHTEEN
  use g_parameters, only : MYRANK, MASTER, P_NP_GLOBAL, P_FOLLOW, &
                           P_KEEP_INITIAL_POSITION, YES, F_TYPE, &
                           F_ISOTROPIC
  use g_domain,     only : bvalexist, binv
  use g_files,      only : FUNIT_PSINGL
  use p_arrays,     only : particle_status, particle_xp, &
                           particle_up, particle_xp0

  implicit none

  integer,intent(in)                         :: npin
  integer,intent(in)                         :: stage
  real(kind=C_DOUBLE),intent(in)             :: deltat
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: np
  integer                                    :: nfollow

  real(kind=C_DOUBLE),dimension(1:3)         :: xp
  real(kind=C_DOUBLE),dimension(1:3)         :: up
  real(kind=C_DOUBLE),dimension(1:3)         :: xp0
  real(kind=C_DOUBLE),dimension(1:3)         :: xplocal
  real(kind=C_DOUBLE),dimension(1:3)         :: uplocal
  real(kind=C_DOUBLE),dimension(1:3)         :: xp0local

  integer                                    :: m, n

  IFSTAGE: if ( stage == 1 ) then

    xp = ZERO
    up = ZERO
    xp0 = ZERO

    IFLARGE: if ( P_FOLLOW > P_NP_GLOBAL ) then
      nfollow = P_NP_GLOBAL
    else
      nfollow = P_FOLLOW
    end if IFLARGE

    !> Nothing to do if there are no particles
    IFPARTEXIST: if ( npin > 0 ) then
      DOPART: do np = 1, npin, 1

        IFTHIS: if ( particle_status(1,np) == nfollow ) then
          xp = particle_xp(:,np)
          up = particle_up(:,np)

          IFINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then 
            xp0 = particle_xp0(:,np)
          end if IFINITIAL

        end if IFTHIS

      end do DOPART
    end if IFPARTEXIST

    IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
      !> Distance is the same
      xplocal  = xp
      uplocal  = up
      xp0local = xp0

    else
      !> Need to transform into laboratory system
      xplocal = ZERO
!      uplocal = ZERO
      uplocal  = up
      xp0local = ZERO
      DOM: do m = 1, 3, 1
        DON: do n = 1, 3, 1
          IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
            xplocal(m) = xplocal(m) + binv(m,n) * xp(n)
!            uplocal(m) = uplocal(m) + binv(m,n) * up(n)
            xp0local(m) = xp0local(m) + binv(m,n) * xp0(n)
          end if IFBEXIST
        end do DON
      end do DOM
    end if IFISOTROPIC

    !> Compute relative velocity
    call MPI_Allreduce (xplocal, xpsingle, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
    call MPI_Allreduce (uplocal, upsingle, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )
    call MPI_Allreduce (xp0local, xp0single, 3, MPI_Double, MPI_SUM, PETSC_COMM_WORLD, ierr )

  end if IFSTAGE

  IFPFOLLOW: if ( P_FOLLOW > 0 ) then
    write(FUNIT_PSINGL,60) deltat, xpsingle, upsingle, xp0single
  end if IFPFOLLOW

60  FORMAT('dt: ',e14.6,' xp: ',3e14.6,' up: ', 3e14.6, ' xp0: ', 3e14.6)

end subroutine p_RK3Follow
!---------------------------------------------------------------------------------

!--------------------------------------------------------------
pure function p_RK3ComputeForce ( uf, up, dp ) result ( acceleration )

    use g_constants,  only : ZERO
    use g_parameters, only : YES, P_CHARGE_PART, P_GRAVITY, GRAVITY, &
                             F_TYPE, F_ISOTROPIC
    use g_domain,     only : bvalexist, bmat, binv

  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(in) :: uf
  real(kind=C_DOUBLE),dimension(1:3),intent(in) :: up
  real(kind=C_DOUBLE),intent(in)                :: dp
  real(kind=C_DOUBLE),dimension(1:3)            :: acceleration
  real(kind=C_DOUBLE),dimension(1:3)            :: drag_accln 
  real(kind=C_DOUBLE),dimension(1:3)            :: elec_accln
  real(kind=C_DOUBLE),dimension(1:3)            :: grav_accln
    
  
  drag_accln = ZERO

  IFGRAVITY: if ( P_GRAVITY == YES ) then

! No need for subroutine, as only velocity in laboratory system is used
    grav_accln = GRAVITY
!    call p_RK3Gravity ( grav_accln )

  else

    grav_accln = ZERO

  end if IFGRAVITY

  drag_accln = p_RK3SolidForce ( uf, up, dp )

  !> Implement Lorentz force here
!  IFCHARGED: if(P_CHARGE_PART == YES) then
!    elec_accln = lorentz_force( npart )
!  else
  elec_accln = ZERO
!  end if IFCHARGED

  acceleration = drag_accln + elec_accln + grav_accln

end function p_RK3ComputeForce
!--------------------------------------------------------------

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
pure function p_RK3SolidForce ( uf, up, dp ) result ( fp ) 

  use g_constants,  only : ONE, EIGHTEEN, POINT_ONE_FIVE, POINT_SIX_EIGHT_SEVEN, ZERO
  use g_parameters, only : F_RHO, F_NU, P_RHO

  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in) :: uf
  real(kind=C_DOUBLE),dimension(1:3),intent(in) :: up
  real(kind=C_DOUBLE),intent(in)                :: dp
  real(kind=C_DOUBLE),dimension(1:3)            :: fp
  real(kind=C_DOUBLE)                           :: tau

  !> Particle relaxation time for small Re_p
  tau      = ( ( dp**2 ) / ( EIGHTEEN * F_NU ) ) * ( P_RHO / F_RHO ) 

  !> Drag acceleration 
  fp       = ( (uf - up) / tau ) 

end function p_RK3SolidForce
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
pure subroutine p_RK3Gravity ( gravityacc )

  use g_constants,  only : ZERO
  use g_parameters, only : GRAVITY, F_TYPE, F_ISOTROPIC
  use g_domain,     only : bvalexist, bmat, binv

  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(inout) :: gravityacc

  integer                                          :: m, n

  !> For isotropic turbulence apply gravity in the laboratory system
  IFISOTROPIC: if ( F_TYPE == F_ISOTROPIC ) then
    gravityacc = GRAVITY

  !> In the Rogallo domain 
  else

    gravityacc = ZERO
    DOM: do m = 1, 3, 1
      DON: do n = 1, 3, 1
        IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
          gravityacc(m) = gravityacc(m) + ( bmat(m,n) * GRAVITY(n) )
        end if IFBEXIST
      end do DON
    end do DOM

  end if IFISOTROPIC

end subroutine p_RK3Gravity
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_RK3Mov2Lab
!> Convert position or velocity from moving to laboratory system
!> @param ierr should return 0
pure subroutine p_RK3Mov2Lab ( xumov, xulab, ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : binv, bvalexist
 
  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in)  :: xumov
  real(kind=C_DOUBLE),dimension(1:3),intent(out) :: xulab
  PetscErrorCode,intent(inout)                          :: ierr

  integer                                        :: m,n

  xulab = ZERO

  DOM: do m = 1, 3, 1
    DON: do n = 1, 3, 1
      IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
        xulab(m) = xulab(m) + binv(m,n) * xumov(n)
      end if IFBEXIST
    end do DON
  end do DOM

end subroutine p_RK3Mov2Lab
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_RK3Lab2Mov
!> Convert position or velocity from moving to laboratory system
!> @param ierr should return 0
pure subroutine p_RK3Lab2Mov ( xulab, xumov, ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : bmat, bvalexist
 
  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in)  :: xulab
  real(kind=C_DOUBLE),dimension(1:3),intent(out) :: xumov
  PetscErrorCode,intent(inout)                          :: ierr

  integer                                        :: m,n

  xumov = ZERO

  DOM: do m = 1, 3, 1
    DON: do n = 1, 3, 1
      IFBEXIST: if ( bvalexist(m,n)  .eqv. .TRUE. ) then
        xumov(m) = xumov(m) + bmat(m,n) * xulab(n)
      end if IFBEXIST
    end do DON
  end do DOM

end subroutine p_RK3Lab2Mov
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_RK3UAbs
!> Obtain absolute particle velocity from particle position
!> @param ierr should return 0
pure subroutine p_RK3UAbs ( xp, upabs, ierr )

  use g_constants,  only : ZERO
  use g_domain,     only : amat, avalexist
 
  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in)  :: xp
  real(kind=C_DOUBLE),dimension(1:3),intent(out) :: upabs
  PetscErrorCode,intent(inout)                          :: ierr

  integer                                        :: m,n

  upabs = ZERO

  DOM: do m = 1, 3, 1
    DON: do n = 1, 3, 1
      IFBEXIST: if ( avalexist(m,n)  .eqv. .TRUE. ) then
        upabs(m) = upabs(m) + amat(m,n) * xp(n)
      end if IFBEXIST
    end do DON
  end do DOM

end subroutine p_RK3UAbs
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_RKComputeHomogeneousSourceTerm
!> Obtain pseudoforce from homogeneous shear (see Shotorban & Balachandar 2006,
!! eq. 25)
pure function p_RK3ComputeHomogeneousSourceTerm(up) result (accln)

  use g_constants,  only : ZERO
  use g_domain,     only : amat, avalexist
 
  implicit none
 
  real(kind=C_DOUBLE),dimension(1:3),intent(in)  :: up
  real(kind=C_DOUBLE),dimension(1:3)             :: accln

  integer                                        :: m,n

  accln = ZERO

  DOM: do m = 1, 3, 1
    DON: do n = 1, 3, 1
      IFBEXIST: if ( avalexist(m,n)  .eqv. .TRUE. ) then
        accln(m) = accln(m) + amat(m,n) * up(n)
      end if IFBEXIST
    end do DON
  end do DOM

end function p_RK3ComputeHomogeneousSourceTerm
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
!function lorentz_force(npart)
    
!    use g_constants, only : zero, onepi, six
!    use p_control,   only : part_data,part_charge
!    use p_interp,    only : interp_order, interp_type, field_el, interp_uf
!    use p_values,    only : particle_rho

!    implicit none

!    integer,intent(in) :: npart 
!    real(kind=C_DOUBLE) :: lorentz_force(3)
!    real(kind=C_DOUBLE) :: ef(3)
!    real(kind=C_DOUBLE) :: xp(3)
!    real(kind=C_DOUBLE) :: diam, mp
 
!    ef = zero

!    xp = part_data(npart,1:3)
!    diam = part_data(npart,10)

!    call interp_uf(interp_order,xp,ef,interp_type,field_el)
!
!    mp = particle_rho * ( (onepi/six) * (diam**3) )

!    lorentz_force = (ef * part_charge(npart)) / mp  

!end function
!--------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------------------------------------
!subroutine calc_efield(ierr)

!    use g_constants, only : epsilon, zero, c_imag, c_zero
!    use p_control, only : np_total, part_data, part_stat, part_chr, part_efieldr, part_efieldk, part_chk, part_charge, &
!                          part_potk, part_potr, part_ek_verify, part_er_verify
!    use p_values,  only : np_packet_ch
!    use f_dfts,    only : f_dfts_rsks, f_dfts_ksrs
!    use f_arrays,  only : ki_min, ki_max, kj_min, kj_max, kk_min, kk_max, fluid_k
!    use d_arrays,  only : domain_dx
!    use p_interp,  only : interp_find

!    implicit none
!    PetscErrorCode, intent(inout) :: ierr
!    integer :: np
!    real(kind=C_DOUBLE) :: xp(3), kmag2
!    integer :: i,j,k,ii,jj,kk
!    complex(kind=C_DOUBLE) :: kdotu

!    part_chr = zero

!    part_potk = c_zero
!    part_potr = zero

!    part_ek_verify = c_zero
!    part_er_verify = zero

!    do np = 1,np_total
!        if( part_stat(np).ne.0 )then
!            xp = part_data(np,1:3)
!            call interp_find(xp,ii,jj,kk)
!
!            part_chr(ii,jj,kk,1) = part_chr(ii,jj,kk,1) + part_charge(np)

!        endif
!    enddo

!    part_chr = part_chr * (np_packet_ch / ( (domain_dx**3) * epsilon ))
!
!    call f_dfts_rsks(part_chr,part_chk,1,1,ierr)      ! transform rs -> fs

!a0: do i=ki_min,ki_max
!a1:     do j=kj_min,kj_max
!a2:         do k=kk_min,kk_max

                ! calculate dot product
!                kdotu = fluid_k(i,j,k,1) * part_efieldk(i,j,k,1) + &
!                        fluid_k(i,j,k,2) * part_efieldk(i,j,k,2) + &
!                        fluid_k(i,j,k,3) * part_efieldk(i,j,k,3)

                ! calculate square of wave vector magnitude
!                kmag2 =        ( fluid_k(i,j,k,1) * fluid_k(i,j,k,1) ) + &
!                               ( fluid_k(i,j,k,2) * fluid_k(i,j,k,2) ) + &
!                               ( fluid_k(i,j,k,3) * fluid_k(i,j,k,3) )

!                if(kmag2.gt.zero)then 

                    ! enforce continuity in fs
!                    part_efieldk(i,j,k,1) = part_efieldk(i,j,k,1) - ( fluid_k(i,j,k,1) * (kdotu - (-1.0 * c_imag * part_chk(i,j,k,1)) ) / kmag2 )
!                    part_efieldk(i,j,k,2) = part_efieldk(i,j,k,2) - ( fluid_k(i,j,k,2) * (kdotu - (-1.0 * c_imag * part_chk(i,j,k,1)) ) / kmag2 )
!                    part_efieldk(i,j,k,3) = part_efieldk(i,j,k,3) - ( fluid_k(i,j,k,3) * (kdotu - (-1.0 * c_imag * part_chk(i,j,k,1)) ) / kmag2 )


!                    part_potk(i,j,k,1) = (part_chk(i,j,k,1) / kmag2)
!                end if

                
!            end do a2
!        end do a1
!    end do a0


! do i=ki_min,ki_max
!     do j=kj_min,kj_max
!         do k=kk_min,kk_max
!                part_ek_verify(i,j,k,1) = -1.0 * c_imag * fluid_k(i,j,k,1) * part_potk(i,j,k,1)
!                part_ek_verify(i,j,k,2) = -1.0 * c_imag * fluid_k(i,j,k,2) * part_potk(i,j,k,1)
!                part_ek_verify(i,j,k,3) = -1.0 * c_imag * fluid_k(i,j,k,3) * part_potk(i,j,k,1)
!         enddo
!     enddo
! enddo

!    call f_dfts_ksrs(part_efieldr,part_efieldk,1,3,ierr)

!    call f_dfts_ksrs(part_potr,part_potk,1,1,ierr)

!    call f_dfts_ksrs(part_er_verify,part_ek_verify,1,3,ierr)

!end subroutine
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
end module p_rk3

