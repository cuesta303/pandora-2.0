!----------------------------------------------------------------------
! module g_domain
!> Define the length and resolution of the domain both in real space
!! and in wave space
module g_domain
 
!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

use iso_c_binding
!------- data section begins ----------------------
 
implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

!----------------------------------------------------------------------
!domain parameters 
!----------------------------------------------------------------------
!> Pandora uses a moving frame for the computation. The domain 
!! length LX_MOVINGFRAME in x direction is fixed.
real(kind=C_DOUBLE),save,public                :: LX_MOVINGFRAME = 0.0
!> Pandora uses a moving frame for the computation. The domain 
!! length LY_MOVINGFRAME in y direction is fixed.
real(kind=C_DOUBLE),save,public                :: LY_MOVINGFRAME = 0.0
!> Pandora uses a moving frame for the computation. The domain 
!! length LZ_MOVINGFRAME in z direction is fixed.
real(kind=C_DOUBLE),save,public                :: LZ_MOVINGFRAME = 0.0
!> Pandora uses a moving frame for the computation. The domain 
!! resolution DX_MOVINGFRAME in x direction is fixed.
real(kind=C_DOUBLE),save,public                :: DX_MOVINGFRAME = 0.0
!> Pandora uses a moving frame for the computation. The domain 
!! resolution DY_MOVINGFRAME in y direction is fixed.
real(kind=C_DOUBLE),save,public                :: DY_MOVINGFRAME = 0.0
!> Pandora uses a moving frame for the computation. The domain 
!! resolution DZ_MOVINGFRAME in z direction is fixed.
real(kind=C_DOUBLE),save,public                :: DZ_MOVINGFRAME = 0.0
!> The domain length lx_laboratoryframe in the laboratory frame 
!! is variable (depending on the deformation tensor)
real(kind=C_DOUBLE),save,public                :: lx_laboratoryframe = 0.0
!> The domain length ly_laboratoryframe in the laboratory frame 
!! is variable (depending on the deformation tensor)
real(kind=C_DOUBLE),save,public                :: ly_laboratoryframe = 0.0
!> The domain length lz_laboratoryframe in the laboratory frame 
!! is variable (depending on the deformation tensor)
real(kind=C_DOUBLE),save,public                :: lz_laboratoryframe = 0.0
!> The domain resolution dx_laboratoryframe in the laboratory frame 
!! is variable (depending on the deformation tensor)
real(kind=C_DOUBLE),save,public                :: dx_laboratoryframe = 0.0
!> The domain resolution dy_laboratoryframe in the laboratory frame 
!! is variable (depending on the deformation tensor)
real(kind=C_DOUBLE),save,public                :: dy_laboratoryframe = 0.0
!> The domain resolution dz_laboratoryframe in the laboratory frame 
!! is variable (depending on the deformation tensor)
real(kind=C_DOUBLE),save,public                :: dz_laboratoryframe = 0.0

real(kind=C_DOUBLE),save,public                :: XMIN_MOVINGFRAME = 0.0
real(kind=C_DOUBLE),save,public                :: XMAX_MOVINGFRAME = 0.0
real(kind=C_DOUBLE),save,public                :: YMIN_MOVINGFRAME = 0.0
real(kind=C_DOUBLE),save,public                :: YMAX_MOVINGFRAME = 0.0
real(kind=C_DOUBLE),save,public                :: ZMIN_MOVINGFRAME = 0.0
real(kind=C_DOUBLE),save,public                :: ZMAX_MOVINGFRAME = 0.0

!> length of the halo in x direction
real(kind=C_DOUBLE),save,public                :: LX_HALO = 0.0
!> length of the halo in y direction
real(kind=C_DOUBLE),save,public                :: LY_HALO = 0.0
!> length of the halo in z direction
real(kind=C_DOUBLE),save,public                :: LZ_HALO = 0.0

!> transformation tensor
real(kind=C_DOUBLE),save,public,dimension(3,3) :: bmat = 0.0   
!> inverse of the transformation tensor     
real(kind=C_DOUBLE),save,public,dimension(3,3) :: binv = 0.0        
!> fluid velocity deformation tensor
real(kind=C_DOUBLE),save,public,dimension(3,3) :: amat = 0.0 
real(kind=C_DOUBLE),save,public,dimension(3,3) :: rk3mat = 0.0 
!real(kind=C_DOUBLE),save,public,dimension(3,3) :: FLUID_UGRAD = 0.0
!> boolean array for computations of homogeneous turbulence 
logical,save,public,dimension(3,3)             :: avalexist = .FALSE.
logical,save,public,dimension(3,3)             :: bvalexist = .FALSE.

integer,save,public                            :: KI_MIN = 0
integer,save,public                            :: KI_MAX = 0
integer,save,public                            :: NODES_KI = 0
integer,save,public                            :: KJ_MIN = 0
integer,save,public                            :: KJ_MAX = 0
integer,save,public                            :: NODES_KJ = 0
integer,save,public                            :: KK_MIN = 0
integer,save,public                            :: KK_MAX = 0
integer,save,public                            :: NODES_KK = 0
real(kind=C_DOUBLE),save,public                :: KMAX = 0

real(kind=C_DOUBLE),save,public                :: K0I = 0.0
real(kind=C_DOUBLE),save,public                :: K0J = 0.0
real(kind=C_DOUBLE),save,public                :: K0K = 0.0

public :: g_DomainInit

!------- data section ends ----------------------

contains
!---------------------------------------------------------------------------------
! subroutine g_DomainInit
!> Call g_DomainCalc to calculate the domain parameters and then call g_DomainHeader
!! to write information about the domain to the program header.
!> @param ierr should return 0
subroutine g_DomainInit ( ierr )

  use g_parameters, only : MYRANK, MASTER, P_TRACK_PART, YES !, read_file, write_file, v_test, no

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
    write(*,10)
  end if IFMASTER

  call g_DomainTypeCheck ( ierr )
  call g_DomainCalc ( ierr )

  IFPART: if ( P_TRACK_PART == YES ) then
    call PetscPrintf( PETSC_COMM_WORLD, '    | Setting number of halo nodes \n', ierr )
    call g_DomainInitHalo ( ierr )
  end if IFPART

  call g_DomainHeader ( ierr )

10  format('==> Initialising g_domain module')

end subroutine g_DomainInit

!---------------------------------------------------------------------------------
! subroutine g_DomainCalc 
!> Call g_DomainCheck to check that the number of nodes in real space in each
!! direction is reasonable.
!! Calculate the domain length and resolution in the moving and laboratory frames.
!! Calculate the minimum and maximum wave numbers for the wave space arrays.
!> @param ierr should return 0
subroutine g_DomainCalc ( ierr )

  use g_parameters, only : MYRANK, MASTER, &
                           NODES_X, NODES_Y, NODES_Z, B110, B220, B330, &
                           F_DEFORM_A, F_DEFORM_B, F_DEFORM_C, F_DEFORM_S, &
                           F_TYPE, F_ISOTROPIC_DEFORMED 
  use g_constants, only  : TWOPI, ZERO, ONE

  implicit none

  PetscErrorCode,intent(inout)   :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then
!    write(*,10)
!    call g_DomainCheck ( NODES_X, ierr )
!    call g_DomainCheck ( NODES_Y, ierr )
!    call g_DomainCheck ( NODES_Z, ierr )
    write(*,20)
  end if IFMASTER

  !> domain lengths in moving frame
  LX_MOVINGFRAME = TWOPI
  LY_MOVINGFRAME = TWOPI
  LZ_MOVINGFRAME = TWOPI

  !> grid spacing in moving frame
  DX_MOVINGFRAME = LX_MOVINGFRAME / real(NODES_X,kind=C_DOUBLE)
  DY_MOVINGFRAME = LY_MOVINGFRAME / real(NODES_Y,kind=C_DOUBLE)
  DZ_MOVINGFRAME = LZ_MOVINGFRAME / real(NODES_Z,kind=C_DOUBLE)

  !> domain lengths in laboratory frame
  lx_laboratoryframe = TWOPI * real(B110,kind=C_DOUBLE)
  ly_laboratoryframe = TWOPI * real(B220,kind=C_DOUBLE)
  lz_laboratoryframe = TWOPI * real(B330,kind=C_DOUBLE)

  !> grid spacing in laboratory frame
  dx_laboratoryframe = lx_laboratoryframe / real(NODES_X,kind=C_DOUBLE)
  dy_laboratoryframe = ly_laboratoryframe / real(NODES_Y,kind=C_DOUBLE)
  dz_laboratoryframe = lz_laboratoryframe / real(NODES_Z,kind=C_DOUBLE)

  !> initialise the deformation tensor
  amat = ZERO
  avalexist = .FALSE.

  !> initialise the inverse of the transformation tensor
  binv = ZERO
  bvalexist = .FALSE.

  binv(1,1) = real(B110,kind=C_DOUBLE)
  bvalexist(1,1) = .TRUE.
  binv(2,2) = real(B220,kind=C_DOUBLE)
  bvalexist(2,2) = .TRUE.
  binv(3,3) = real(B330,kind=C_DOUBLE)
  bvalexist(3,3) = .TRUE.

  IFNOTISOTROPICDEFORMED: if ( F_TYPE /= F_ISOTROPIC_DEFORMED ) then
  !> Shear also exists (in all currently allowed cases)
  bvalexist(1,3) = .TRUE.
  amat(1,3) = F_DEFORM_S
  avalexist(1,3) = .TRUE.
  end if IFNOTISOTROPICDEFORMED  

  !> initialise the transformation tensor
  bmat = ZERO

  bmat(1,1) = ONE/binv(1,1)
  bmat(2,2) = ONE/binv(2,2)
  bmat(3,3) = ONE/binv(3,3)

  !> wave space array boundaries,
  !! for N real-space nodes -N/2+1 to N/2 needed
  !! using conjugate symmetry in i direction: 0 to N/2
  !! dealiasing with 2/3 rule: only wave modes from -N/3 (not -N/3 +1!) to N/3
  KI_MIN  =  0
  KI_MAX  =  NODES_X / 3
  NODES_KI = KI_MAX + 1
  KJ_MIN  =  0
  KJ_MAX  =  NODES_Y / 3
  NODES_KJ = 2 * KJ_MAX + 1
  KK_MIN  =  0
  KK_MAX  =  NODES_Z / 3
  NODES_KK = 2 * KK_MAX + 1
  
  !> kmax for resolution
  KMAX = real(min( KI_MAX, KJ_MAX, KK_MAX ),kind=C_DOUBLE)

  !> calculate minimum wave number in each direction
  K0I      = ONE !TWOPI/lx_laboratoryframe 
  K0J      = ONE !TWOPI/ly_laboratoryframe 
  K0K      = ONE !TWOPI/lz_laboratoryframe 
!  K0I      = TWOPI/lx_laboratoryframe 
!  K0J      = TWOPI/ly_laboratoryframe 
!  K0K      = TWOPI/lz_laboratoryframe 

10  format('    | Checking domain size')
20  format('    | Computing domain parameters')
end subroutine g_DomainCalc


!---------------------------------------------------------------------------------
! subroutine g_DomainCheck 
!> Check if number of nodes in real space is positive. Check if number of nodes in
!! real space is a power of 2 in order to efficiently use FFTW.
!> @param nodes is the number of nodes supplied by the calling routine
!> @param ierr should return 0
subroutine g_DomainCheck ( nodes, ierr )

  use g_parameters, only : MYRANK, MASTER

  implicit none

  integer,intent(in)              :: nodes
  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: status
  integer                         :: poweroftwo
  integer                         :: factorthree 
  integer                         :: factorfive
  integer                         :: factorseven

  status = 0
  poweroftwo = 1

  IFMASTER: if ( MYRANK == MASTER ) then
    IFSMALL: if ( nodes <= poweroftwo ) then
      status = 1
      write(*,10)
    else
      WHILELESS: do while ( poweroftwo < nodes ) 
        poweroftwo = poweroftwo * 2
        IFNOTPOWEROFTWO: if ( poweroftwo > nodes ) then
          !> Multiplication by 3, 5 or 7 is now allowed.
          poweroftwo = poweroftwo / 4
          factorthree = poweroftwo * 3
          poweroftwo = poweroftwo / 2
          factorfive = poweroftwo * 5
          factorseven = poweroftwo * 7
          IFNOTFACTORTWO: if ( factorthree == nodes ) then
            !> factor 3 allowed
          else if ( factorfive == nodes ) then
            !> factor 5 allowed
          else if ( factorseven == nodes ) then
            !> factor 7 allowed
          else
            status = 1
            write(*,20)
          end if IFNOTFACTORTWO
        end if IFNOTPOWEROFTWO
      end do WHILELESS
    end if IFSMALL

    IFABORT: if ( status /= 0 ) then
      call MPI_Abort ( PETSC_COMM_WORLD, status, ierr )
    end if IFABORT
  end if IFMASTER

10  format('stop: domain size negative or zero!')
20  format('stop: number of nodes not power of two!')

end subroutine g_DomainCheck


!---------------------------------------------------------------------------------
! subroutine g_DomainTypeCheck 
!> Check if flow type exists. Write flow type to standard out or abort.
!> @param ierr should return 0
subroutine g_DomainTypeCheck ( ierr )

  use g_parameters, only : MYRANK, MASTER, F_TYPE, F_ISOTROPIC, F_SHEAR_DNS, &
                           F_SHEAR_RDT_INVISCID, F_SHEAR_RDT_VISCOUS, F_ISOTROPIC_DEFORMED

  implicit none

  PetscErrorCode,intent(inout)           :: ierr
  integer                         :: status

  status = 0

  IFMASTER: if ( MYRANK == MASTER ) then

!> Write to the header which flow type has been chosen.
!! If no flow type is chosen the simulation is aborted.
    IFFLOWTYPE: if ( F_TYPE == F_ISOTROPIC ) then
      write(*,10) 
    else if ( F_TYPE == F_SHEAR_DNS ) then
      write(*,20) 
    else if ( F_TYPE == F_SHEAR_RDT_INVISCID ) then
      write(*,30) 
    else if ( F_TYPE == F_SHEAR_RDT_VISCOUS ) then
      write(*,40) 
    else if ( F_TYPE == F_ISOTROPIC_DEFORMED ) then
      write(*,50) 
    else
      write (*,60)
      status = 1
      call MPI_Abort ( PETSC_COMM_WORLD, status, ierr )
    end if IFFLOWTYPE

  end if IFMASTER

10  format('    | Simulating homogeneous isotropic turbulence')
20  format('    | Simulating homogeneous shear flow (DNS)')
30  format('    | Simulating homogeneous shear flow (inviscid RDT)')
40  format('    | Simulating homogeneous shear flow (viscous RDT)')
50  format('    | Simulating homogeneous isotropic turbulence on a deformed mesh')
60  format('stop: flow type not defined!')

end subroutine g_DomainTypeCheck


!---------------------------------------------------------------------------------
! subroutine g_DomainInitHalo
!> Set number of fluid halo nodes needed for particle interpolation
!> @param ierr should return 0
subroutine g_DomainInitHalo ( ierr )

  use g_constants,  only : HALF
  use g_parameters, only : NHALO, P_INTERP_ORDER 

  implicit none

  PetscErrorCode, intent(inout) :: ierr

  !> Halo chosen such that there is enough information at the centre of an interpolation molecule.
  NHALO = P_INTERP_ORDER / 2

  !> Determine length of halo in each direction @todo might need to be updated
  LX_HALO = real(NHALO,kind=C_DOUBLE) * DX_MOVINGFRAME
  LY_HALO = real(NHALO,kind=C_DOUBLE) * DY_MOVINGFRAME
  LZ_HALO = real(NHALO,kind=C_DOUBLE) * DZ_MOVINGFRAME

  !> Needs to be corrected if even order
  IFEVEN: if ( modulo( P_INTERP_ORDER, 2 ) == 0 ) then
    LX_HALO = LX_HALO - ( HALF * DX_MOVINGFRAME )
    LY_HALO = LY_HALO - ( HALF * DY_MOVINGFRAME )
    LZ_HALO = LZ_HALO - ( HALF * DZ_MOVINGFRAME )
  end if IFEVEN

end subroutine g_DomainInitHalo
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_DomainHeader
!> Write the relevant information about the domain to the program header
!! depending on the flow type.
!> @param ierr should return 0
subroutine g_DomainHeader ( ierr )

  use g_parameters, only : MYRANK, MASTER, &
                           F_TYPE, F_ISOTROPIC, &
                           NODES_X, NODES_Y, NODES_Z
  use g_constants,  only : ONEPI

  implicit none

  PetscErrorCode,intent(inout)           :: ierr

  IFMASTER: if ( MYRANK == MASTER ) then

    IFDEFORMED: if ( F_TYPE /= F_ISOTROPIC ) then
      write(*,50)  bmat(1,1)
      write(*,60)  bmat(2,2)
      write(*,70)  bmat(3,3)
      write(*,80)
      write(*,90)  amat(1,1)
      write(*,100) amat(2,2)
      write(*,110) amat(3,3)
      write(*,120) amat(1,3)
      write(*,130) 
      write(*,140) lx_laboratoryframe / ONEPI, ly_laboratoryframe / ONEPI, lz_laboratoryframe / ONEPI
    end if IFDEFORMED

    write(*,150) NODES_KI
    write(*,160) NODES_KJ
    write(*,170) NODES_KK

  end if IFMASTER

50  format('    | B11(0) = ', F6.1)  
60  format('    | B22(0) = ', F6.1) 
70  format('    | B33(0) = ', F6.1) 
80  format('    | Fluid deformation tensor:')  
90  format('    | a = ', F6.1)  
100 format('    | b = ', F6.1) 
110 format('    | c = ', F6.1) 
120 format('    | s = ', F6.1) 
130 format('    | Initial domain dimensions in the laboratory system: ') 
140 format('    | lx = ', F6.1, 'pi,  ly = ', F6.1, 'pi,  lz = ', F6.1, 'pi') 
150 format('    | ', I6, ' wave modes in x direction') 
160 format('    | ', I6, ' wave modes in y direction') 
170 format('    | ', I6, ' wave modes in z direction') 

end subroutine g_DomainHeader
!---------------------------------------------------------------------------------
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!---------------------------------------------------------------------------------
end module g_domain
