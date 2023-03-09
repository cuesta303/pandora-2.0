!-----------------------------------------------------------------
! module g_constants
!> In this module the most common constants used throughout the
!! program are defined.
module g_constants

!------- data section begins ----------------------
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

!-----------------------------------------------------------------
! Petsc integers
!-----------------------------------------------------------------
PetscInt,parameter,public                         :: PETSCONE   = 1
PetscInt,parameter,public                         :: PETSCTWO   = 2
PetscInt,parameter,public                         :: PETSCTHREE = 3
!-----------------------------------------------------------------
!'real integers' 
!-----------------------------------------------------------------
real(kind=C_DOUBLE),parameter,public              :: ZERO                          =    0.0
real(kind=C_DOUBLE),parameter,public              :: ONE                           =    1.0
real(kind=C_DOUBLE),parameter,public              :: TWO                           =    2.0
real(kind=C_DOUBLE),parameter,public              :: THREE                         =    3.0
real(kind=C_DOUBLE),parameter,public              :: FOUR                          =    4.0
real(kind=C_DOUBLE),parameter,public              :: FIVE                          =    5.0
real(kind=C_DOUBLE),parameter,public              :: SIX                           =    6.0
real(kind=C_DOUBLE),parameter,public              :: SEVEN                         =    7.0
real(kind=C_DOUBLE),parameter,public              :: EIGHT                         =    8.0
real(kind=C_DOUBLE),parameter,public              :: NINE                          =    9.0
real(kind=C_DOUBLE),parameter,public              :: TEN                           =   10.0
real(kind=C_DOUBLE),parameter,public              :: ELEVEN                        =   11.0
real(kind=C_DOUBLE),parameter,public              :: TWELVE                        =   12.0
real(kind=C_DOUBLE),parameter,public              :: FIFTEEN                       =   15.0
real(kind=C_DOUBLE),parameter,public              :: SIXTEEN                       =   16.0
real(kind=C_DOUBLE),parameter,public              :: EIGHTEEN                      =   18.0
real(kind=C_DOUBLE),parameter,public              :: THIRTYFIVE                    =   35.0
real(kind=C_DOUBLE),parameter,public              :: ONEFIVETHREE                  =  153.0
real(kind=C_DOUBLE),parameter,public              :: ONETWOEIGHT                   =  128.0

real(kind=C_DOUBLE),parameter,public              :: HALF                          = ONE/TWO
real(kind=C_DOUBLE),parameter,public              :: THIRD                         = ONE/THREE
real(kind=C_DOUBLE),parameter,public              :: SIXTH                         = ONE/SIX
real(kind=C_DOUBLE),parameter,public              :: THREE_OVER_TWO                = THREE/TWO
real(kind=C_DOUBLE),parameter,public              :: TWO_OVER_THREE                = TWO/THREE
real(kind=C_DOUBLE),parameter,public              :: FIVE_OVER_TWO                 = FIVE/TWO
real(kind=C_DOUBLE),parameter,public              :: FIVE_OVER_THREE               = FIVE/THREE
real(kind=C_DOUBLE),parameter,public              :: FIVE_OVER_NINE                = FIVE/NINE
real(kind=C_DOUBLE),parameter,public              :: ELEVEN_OVER_SIX               = ELEVEN/SIX
real(kind=C_DOUBLE),parameter,public              :: ELEVEN_OVER_THREE             = ELEVEN/THREE
real(kind=C_DOUBLE),parameter,public              :: EIGHT_OVER_FIFTEEN            = EIGHT/FIFTEEN
real(kind=C_DOUBLE),parameter,public              :: FIFTEEN_OVER_SIXTEEN          = FIFTEEN/SIXTEEN
real(kind=C_DOUBLE),parameter,public              :: ONEFIVETHREE_OVER_ONETWOEIGHT = ONEFIVETHREE/ONETWOEIGHT
real(kind=C_DOUBLE),parameter,public              :: FOUR_OVER_THREE               = FOUR/THREE

real(kind=C_DOUBLE),parameter,public              :: POINT_ONE_FIVE                = 0.15
real(kind=C_DOUBLE),parameter,public              :: POINT_SIX_EIGHT_SEVEN         = 0.687

real(kind=C_DOUBLE),parameter,public              :: EPSILON         =  8.854e-12 

!-----------------------------------------------------------------
!'complex integers' 
!-----------------------------------------------------------------
!> @var C_IMAG
!> \f$ 0 + 1i \f$
!> @var C_ZERO
!> \f$ 0 + 0i \f$
!> @var C_ONE
!> \f$ 1 + 0i \f$
!> @var C_ONEONE
!> \f$ 1 + 1i \f$
complex(kind=C_DOUBLE),parameter,public           :: C_IMAG    = (0.0,1.0)
complex(kind=C_DOUBLE),parameter,public           :: C_ZERO    = (0.0,0.0)
complex(kind=C_DOUBLE),parameter,public           :: C_ONE     = (1.0,0.0)
complex(kind=C_DOUBLE),parameter,public           :: C_ONEONE  = (1.0,1.0)

!-----------------------------------------------------------------
!pi 
!-----------------------------------------------------------------
real(kind=C_DOUBLE),save,public                   :: ONEPI
real(kind=C_DOUBLE),save,public                   :: TWOPI
real(kind=C_DOUBLE),save,public                   :: FOURPI
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!> higher precision values of a few constants
real(kind=C_LONG_DOUBLE),save,public              :: LONGONEPI
real(kind=C_LONG_DOUBLE),save,public              :: LONGTWOPI
real(kind=C_LONG_DOUBLE),save,public              :: LONGFOURPI
real(kind=C_LONG_DOUBLE),parameter,public         :: LONGTWO  = 2.0
real(kind=C_LONG_DOUBLE),parameter,public         :: LONGFOUR = 4.0
!-----------------------------------------------------------------



public :: g_ConstantsInit
!------- data section ends ----------------------

contains
!---------------------------------------------------------------------------------
! subroutine g_ConstantsInit
!> Compute \f$ \pi \f$, \f$ 2\,\pi \f$ and \f$ 4\,\pi \f$.
!> @param ierr should return 0
subroutine g_ConstantsInit ( ierr )

  implicit none

  include 'mpif.h'

  PetscErrorCode,intent(inout) :: ierr

 
  !> A few digits of pi taken from
  !! http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
  LONGONEPI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534

  LONGTWOPI = LONGTWO * LONGONEPI
  LONGFOURPI = LONGFOUR * LONGONEPI

  ONEPI  = real(LONGONEPI,kind=C_DOUBLE)
  TWOPI  = real(LONGTWOPI,kind=C_DOUBLE)
  FOURPI = real(LONGFOURPI,kind=C_DOUBLE)

!  ONEPI =  TWO   * asin(ONE)
!  TWOPI =  FOUR  * asin(ONE)
!  FOURPI = EIGHT * asin(ONE)

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising g_constants module \n', ierr )

end subroutine g_ConstantsInit
!---------------------------------------------------------------------------------
end module g_constants
