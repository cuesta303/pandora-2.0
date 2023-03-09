!-----------------------------------------------------------------------
! module g_files
!> Handle all input and output files centrally.
module g_files

!------- data section begins ----------------------
!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
use g_petsc
!use petscsys

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>


private

!-----------------------------------------------------------------------
!file units 
!-----------------------------------------------------------------------
!> file unit for the configuration file
integer,parameter,public                   :: FUNIT_CONFIG   = 10
!> file unit for the input file 
integer,parameter,public                   :: FUNIT_INPUT    = 11
!> file unit for the log file 
integer,parameter,public                   :: FUNIT_OUTPUT   = 12

!> file unit for the fluid statistics file 
integer,parameter,public                   :: FUNIT_FSTATF   = 31
!> file unit for the fluid spectrum file 
integer,parameter,public                   :: FUNIT_ESPEC    = 32
!> file unit for the fluid enstrophy file 
integer,parameter,public                   :: FUNIT_FENSXZ   = 33
!> file unit for writing the forcing file 
integer,parameter,public                   :: FUNIT_FORCEW   = 34
!> file unit for reading the forcing file 
integer,parameter,public                   :: FUNIT_FORCER   = 35
!> file unit for the fluid correlations file 
integer,parameter,public                   :: FUNIT_FCORRL   = 36
!> file unit for the electric field correlations file 
integer,parameter,public                   :: FUNIT_ECORRL   = 37
!> file unit for the fluid VTK file 
integer,parameter,public                   :: FUNIT_FVTKXZ   = 38

!> file unit for the homogeneous fluid statistics file 
integer,parameter,public                   :: FUNIT_FSTATH   = 41
!> file unit for the fluid 1D spectrum file 
integer,parameter,public                   :: FUNIT_1DSPEC   = 42

!> file unit for the particle statistics file 
integer,parameter,public                   :: FUNIT_PSTATS   = 50
!> file unit for the particle preferential concentration file 
integer,parameter,public                   :: FUNIT_PCSTAT   = 51
integer,parameter,public                   :: FUNIT_CSTATS   = 52
!> file unit for the particle restart file 
integer,parameter,public                   :: FUNIT_PRESRT   = 53
integer,parameter,public                   :: FUNIT_PDSTAT   = 54
integer,parameter,public                   :: FUNIT_PSINGL   = 55
integer,parameter,public                   :: FUNIT_PFSTAT   = 56
integer,parameter,public                   :: FUNIT_PLOCXZ   = 57
integer,parameter,public                   :: FUNIT_EXZVTK   = 58
integer,parameter,public                   :: FUNIT_VXZVTK   = 59
integer,parameter,public                   :: FUNIT_PRMETA   = 60

integer,parameter,public                   :: FUNIT_FSTATR   = 71
integer,parameter,public                   :: FUNIT_GRIDXY   = 72
integer,parameter,public                   :: FUNIT_VECTXY   = 73
integer,parameter,public                   :: FUNIT_TIMING   = 74

!-----------------------------------------------------------------------
!FILE NAMES 
!-----------------------------------------------------------------------
character(len=100)                         :: INPUT_FILE
character(len=100),save,public             :: RESULTS_DIR
character(len=100),save,public             :: HOME_DIR
character(len=100),save,public             :: RANDOM_SEED_FILE

character(len=20),parameter                :: CONFIG_FILE     = '/pandora.conf       '

character(len=20),parameter                :: FNAME_OUTPUT    = '/pandora.out        '
character(len=20),parameter                :: FNAME_FSTATF    = '/fstatsfs.dat       '
character(len=20),parameter                :: FNAME_ESPEC     = '/fspect.dat         '
character(len=20),parameter                :: FNAME_1DSPEC    = '/fspec1d.dat        '
character(len=20),parameter                :: FNAME_FCORRL    = '/fcorrl.dat         '
character(len=20),parameter                :: FNAME_FENSXZ    = '/fensxz.vtk         '
character(len=20),parameter                :: FNAME_FVTKXZ    = '/fvelxz.vtk         '
character(len=20),parameter                :: FNAME_ECORRL    = '/ecorrl.dat         '

character(len=20),parameter                :: FNAME_FSTATH    = '/fstatsh.dat        '

character(len=20),parameter                :: FNAME_PSTATS    = '/pstats.dat         '
character(len=20),parameter                :: FNAME_PCSTAT    = '/pcstat.dat         '
character(len=20),parameter                :: FNAME_PDSTAT    = '/pd2val.dat         '
character(len=20),parameter                :: FNAME_CSTATS    = '/estats.dat         '
character(len=20),parameter                :: FNAME_PRESRT    = '/presrt.dat         '
character(len=20),parameter                :: FNAME_PRMETA    = '/prsrt.meta         '
character(len=20),parameter                :: FNAME_PSINGL    = '/psingl.dat         '
character(len=20),parameter                :: FNAME_PFSTAT    = '/pfstat.dat         '
character(len=20),parameter                :: FNAME_PLOCXZ    = '/plocxz.vtk         '
character(len=20),parameter                :: FNAME_EXZVTK    = '/exzvtk.vtk         '
character(len=20),parameter                :: FNAME_VXZVTK    = '/vxzvtk.vtk         '

character(len=20),parameter                :: FNAME_FSTATR    = '/fstatsrs.dat       '
character(len=20),parameter                :: FNAME_GRIDXY    = '/gridxy.vtk         '
character(len=20),parameter                :: FNAME_VECTXY    = '/vectxy.vtk         '
character(len=20),parameter                :: FNAME_TIMING    = '/cputime.dat        '

public :: g_FilesInput
public :: g_FilesOutput
public :: g_FilesForcing 
public :: g_FilesFinalise

!------- data section ends ----------------------


contains
!---------------------------------------------------------------------------------
! subroutine g_FilesInput
!> Read the name of the input file and results directory from the command line
!! and open the configuration, input and log files.
!> @param ierr should return 0
subroutine g_FilesInput ( ierr )

  use g_parameters, only : MYRANK, MASTER, NO, OPEN_FILE, READ_FILE, WRITE_FILE, &
                           F_FORCING, YES, F_INIT_TYPE, F_INIT_RESTART, &
                           F_INIT_RESTART_DEFORMED, F_FORCE_FILE_READ

  implicit none

  PetscErrorCode,intent(inout)    :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising g_files module \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Opening files: \n', ierr )

  IFMASTER: if ( MYRANK == MASTER ) then

!> read program arguments
    call getarg ( 1 , INPUT_FILE )
    call getarg ( 2 , RESULTS_DIR )

!> read HOME directory to find configuration file
    call get_environment_variable ( 'PANDORA_HOME', HOME_DIR ) 

!> create output directory if not existing
    call system('mkdir -p '//trim(RESULTS_DIR))

!> open configuration file $HOME/.pandora/pandora.conf
    open ( unit=FUNIT_CONFIG, status='OLD', file=trim(HOME_DIR)//'/.pandora'//trim(CONFIG_FILE), iostat=ierr )

    IFCONFIGERROR: if ( ierr /= 0 ) then
      write(*,*) '   | Configuration file ', trim(HOME_DIR)//'/.pandora'//trim(CONFIG_FILE), ' does not exist.'
    else
      write(*,*) '   | Opening configuration file ', trim(HOME_DIR)//'/.pandora'//trim(CONFIG_FILE), '.' 
    end if IFCONFIGERROR

!> open input file
    call g_FilesIO ( file_oc=OPEN_FILE, file_rw=READ_FILE, file_unit=FUNIT_INPUT, file_name=INPUT_FILE, ierr=ierr )

!> open log file
    call g_FilesIO ( file_oc=OPEN_FILE, file_rw=WRITE_FILE, file_unit=FUNIT_OUTPUT, file_name=FNAME_OUTPUT, ierr=ierr )

!> open forcing file if forced turbulence starting from restart file
    IFFORCED: if ( F_FORCING == YES ) then
      IFRESTART: if ( F_INIT_TYPE == F_INIT_RESTART ) then

        call g_FilesIO ( file_oc=OPEN_FILE, file_rw=READ_FILE, file_unit=FUNIT_FORCER, &
                         file_name=F_FORCE_FILE_READ, ierr=ierr )
 
      else if ( F_INIT_TYPE == F_INIT_RESTART_DEFORMED ) then

        call g_FilesIO ( file_oc=OPEN_FILE, file_rw=READ_FILE, file_unit=FUNIT_FORCER, &
                         file_name=F_FORCE_FILE_READ, ierr=ierr )
 
      end if IFRESTART
    end if IFFORCED

  end if IFMASTER

end subroutine g_FilesInput


!---------------------------------------------------------------------------------
! subroutine g_FilesOutput
!> Call g_FilesHandling to open all output files. This is delegated to a subroutine
!! to make sure all files that have been opened will also be closed later.
!> @param ierr should return 0
subroutine g_FilesOutput ( ierr )

  use g_parameters, only : MYRANK, MASTER, OPEN_FILE, F_RESTART_WRITE, F_FORCING, YES, &
                           F_FORCE_FILE_WRITE, WRITE_FILE

  implicit none

  PetscErrorCode,intent(inout)             :: ierr
    

  call PetscPrintf( PETSC_COMM_WORLD, '==> opening output files \n', ierr )

  IFMASTER: if ( MYRANK == MASTER ) then
 
!> create output directory if not existing
    call system('mkdir -p '//trim(RESULTS_DIR))

!> open all files
    call g_FilesHandling ( OPEN_FILE, ierr )

  end if IFMASTER

end subroutine g_FilesOutput
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine g_FilesForcing
!> Open forcing file for writing.
!> @param ierr should return 0
subroutine g_FilesForcing ( ierr )

  use g_parameters, only : MYRANK, MASTER, OPEN_FILE, F_FORCE_FILE_WRITE, F_FORCE_FILE_READ, &
                           WRITE_FILE, F_INIT_TYPE, F_INIT_RESTART, F_INIT_RESTART_DEFORMED, CLOSE_FILE 

  implicit none

  PetscErrorCode,intent(inout)             :: ierr
    
  IFMASTER: if ( MYRANK == MASTER ) then

    !> If forcing was read from a file, close it in case it is the same as the
    !! one for writing. In that case it needs to be overwritten. 
    IFRESTARTREAD: if ( F_INIT_TYPE == F_INIT_RESTART ) then

      call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_FORCER, &
                       file_name=F_FORCE_FILE_READ, ierr=ierr )
 
    else if ( F_INIT_TYPE == F_INIT_RESTART_DEFORMED ) then

      call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_FORCER, &
                       file_name=F_FORCE_FILE_READ, ierr=ierr )
 
    end if IFRESTARTREAD

    !> Now open the forcing file for writing.
    call g_FilesIO ( file_oc=OPEN_FILE, file_rw=WRITE_FILE, file_unit=FUNIT_FORCEW, &
                     file_name=F_FORCE_FILE_WRITE, ierr=ierr )
 
  end if IFMASTER

end subroutine g_FilesForcing
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine g_FilesFinalise 
!> Close configuration, input and log files and all output files that were opened during initialisation.
!> @param ierr should return 0
subroutine g_FilesFinalise ( ierr )

  use g_parameters, only : MYRANK, MASTER, CLOSE_FILE, &
                           F_FORCING, YES, F_INIT_TYPE, F_INIT_RESTART, &
                           F_INIT_RESTART_DEFORMED, F_FORCE_FILE_READ, &
                           F_RESTART_WRITE, F_FORCE_FILE_WRITE

  implicit none

  PetscErrorCode,intent(inout)    :: ierr

  call PetscPrintf( PETSC_COMM_WORLD, '==> Finalising g_files module \n', ierr )
  call PetscPrintf( PETSC_COMM_WORLD, '    | Closing files: \n', ierr )

  IFMASTER: if ( MYRANK == MASTER ) then

!> close configuration file
    close(unit=FUNIT_CONFIG, iostat=ierr)

    IFCONFIGERROR: if ( ierr /= 0 ) then
      write(*,30) trim(HOME_DIR)//'/.pandora'//trim(CONFIG_FILE)
    else
      write(*,40) trim(HOME_DIR)//'/.pandora'//trim(CONFIG_FILE)
    end if IFCONFIGERROR

!> close input file
    call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_INPUT, file_name=INPUT_FILE, ierr=ierr )

!> close log file
    call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_OUTPUT, file_name=FNAME_OUTPUT, ierr=ierr )

    IFFORCED: if ( F_FORCING == YES ) then
      !> close forcing file if forced turbulence starting from restart file
      IFRESTARTREAD: if ( F_INIT_TYPE == F_INIT_RESTART ) then

        call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_FORCER, &
                         file_name=F_FORCE_FILE_READ, ierr=ierr )
 
      else if ( F_INIT_TYPE == F_INIT_RESTART_DEFORMED ) then

        call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_FORCER, &
                         file_name=F_FORCE_FILE_READ, ierr=ierr )
 
      end if IFRESTARTREAD
      !> close forcing file if forced turbulence and writing to restart file
      IFRESTARTWRITE: if ( F_RESTART_WRITE == YES ) then

        call g_FilesIO ( file_oc=CLOSE_FILE, file_rw=0, file_unit=FUNIT_FORCEW, &
                         file_name=F_FORCE_FILE_WRITE, ierr=ierr )
 
      end if IFRESTARTWRITE
    end if IFFORCED
!> Call g_FilesHandling to close all output files. This is delegated to a subroutine
!! to make sure all files that have been opened will also be closed later.
    call g_FilesHandling ( CLOSE_FILE, ierr )

  end if IFMASTER

30  format('    | Error: cannot close ', 100A)
40  format('    | Closing ', 100A)

end subroutine g_FilesFinalise


!---------------------------------------------------------------------------------
! subroutine g_FilesHandling
!> Call g_FilesIO to either open or close all necessary output files
!> @param file_mode is either OPEN_FILE or CLOSE_FILE passed on by the calling subroutine 
!> @param ierr should return 0
subroutine g_FilesHandling( file_mode, ierr )

  use g_parameters, only : MYRANK, MASTER, WRITE_FILE, YES, &
                           TSTEPS, TS_PLANE, TS_CORR, TS_SPEC3D, TS_SPEC1D, RUNMODE, TIMINGMODE, &
                           P_TRACK_PART, P_CHARGE_PART, P_NP_GLOBAL, P_TS_D2VALS, P_FOLLOW, &
                           F_TYPE, F_ISOTROPIC

  implicit none

  integer,intent(in)            :: file_mode
  PetscErrorCode,intent(inout)         :: ierr
    

!> open or close (depending on the parameter file_mode) all files
  IFMASTER: if ( MYRANK == MASTER ) then
    call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_FSTATF, file_name=FNAME_FSTATF, ierr=ierr )
    IFSPEC3D: if ( TS_SPEC3D <= TSTEPS ) then
      call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_ESPEC, file_name=FNAME_ESPEC, ierr=ierr )
    end if IFSPEC3D

    IFSPEC1D: if ( TS_SPEC1D <= TSTEPS ) then
        call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_1DSPEC, file_name=FNAME_1DSPEC, ierr=ierr )
    end if IFSPEC1D

    call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_FSTATR  , file_name=FNAME_FSTATR  , ierr=ierr )
    call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_FSTATH, file_name=FNAME_FSTATH, ierr=ierr )

!> only open or close files that write information on the x-z plane if TS_PLANE, the time step at which
!! information is written, is less or equal the number of time steps in the simulation
    IFPLANE: if ( TS_PLANE <= TSTEPS ) then
      call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_FVTKXZ, file_name=FNAME_FVTKXZ, ierr=ierr )
      call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_FENSXZ, file_name=FNAME_FENSXZ, ierr=ierr )

!> particle locations on the x-z plane only if particles are introduced
      IFPLANEPART: if ( P_TRACK_PART == YES ) then
        call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_PLOCXZ, file_name=FNAME_PLOCXZ, ierr=ierr )

!> electric field and electric potential only if charged particles
        IFPLANECHARGE: if ( P_CHARGE_PART == YES ) then
          call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_EXZVTK, file_name=FNAME_EXZVTK, ierr=ierr )
          call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_VXZVTK, file_name=FNAME_VXZVTK, ierr=ierr )

        end if IFPLANECHARGE 
      end if IFPLANEPART
    end if IFPLANE

!> only open or close files that write correlations if TS_CORR, the time step at which
!! correlations are written, is less or equal the number of time steps in the simulation
    IFCORR: if ( TS_CORR <= TSTEPS ) then
      call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_FCORRL, file_name=FNAME_FCORRL, ierr=ierr )

!> electric field correlations only if charged particles
      IFCORRPART: if ( P_TRACK_PART == YES ) then
        IFCORRCHARGE: if ( P_CHARGE_PART == YES ) then
          call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_ECORRL, file_name=FNAME_ECORRL, ierr=ierr )
        end if IFCORRCHARGE
      end if IFCORRPART

    end if IFCORR

!> particle statistics only if particles are introduced
      IFPART: if ( P_TRACK_PART == YES ) then
        call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_PSTATS, file_name=FNAME_PSTATS, ierr=ierr )
        call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_PCSTAT, file_name=FNAME_PCSTAT, ierr=ierr )
        call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_PFSTAT, file_name=FNAME_PFSTAT, ierr=ierr )

!> only open D2 statistics file if time step for writing smaller than total number of time steps
        IFD2: if ( P_TS_D2VALS <= TSTEPS ) then
          call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_PDSTAT, file_name=FNAME_PDSTAT, ierr=ierr )
       
        end if IFD2

!> only open single particle file if particle number of followed particles is between 1 and total number of particles
        IFFOLLOWMIN: if ( P_FOLLOW > 0 ) then
            call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_PSINGL, file_name=FNAME_PSINGL, ierr=ierr )

        end if IFFOLLOWMIN

!> particle charge statistics only if charged particles
        IFCHARGE: if ( P_CHARGE_PART == YES ) then
          call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_CSTATS, file_name=FNAME_CSTATS, ierr=ierr )
        
        end if IFCHARGE
      end if IFPART

!> only open or close files for timing information if timing mode is activated
    IFTIMING: if ( RUNMODE == TIMINGMODE ) then
      call g_FilesIO ( file_oc=file_mode, file_rw=WRITE_FILE, file_unit=FUNIT_TIMING, file_name=FNAME_TIMING, ierr=ierr )

    end if IFTIMING
  end if IFMASTER

!> Comments:
!! The files FUNIT_GRIDXY and FUNIT_VECTXY are used only by g_output_gridxy, which is not called.
!! They could be useful for homogeneous shear.
!! The variable TS_FLUIDVTK seems to be unused.

end subroutine g_FilesHandling


!---------------------------------------------------------------------------------
! subroutine g_FilesIO
!> Open or close files for reading or writing.
!> @param file_oc is either OPEN_FILE or CLOSE_FILE
!> @param file_rw is either READ_FILE or WRITE_FILE, where reading is only handled
!! for the configuration and input files
!> @param file_unit is the unit of the file to be opened or closed
!> @param file_name is the name of the file to be opened or closed
!> @param ierr should return 0
subroutine g_FilesIO ( file_oc, file_rw, file_unit, file_name, ierr )

  use g_parameters, only : MYRANK, MASTER, OPEN_FILE, CLOSE_FILE, READ_FILE, WRITE_FILE

  implicit none

  integer,intent(in)            :: file_oc, file_rw, file_unit
  character(len=20),intent(in)  :: file_name
  PetscErrorCode,intent(inout)         :: ierr
    
  IFMASTER: if ( MYRANK == MASTER ) then

    IFOPENCLOSE: if ( file_oc == OPEN_FILE ) then
    ! open file
      IFREADWRITE: if ( file_rw == WRITE_FILE ) then
        open ( unit=file_unit, status='REPLACE', file=trim(RESULTS_DIR)//trim(file_name), iostat=ierr )
      else if ( file_rw == READ_FILE ) then
        open ( unit=file_unit, status='OLD', file=trim(file_name), iostat=ierr )
      end if IFREADWRITE

      IFOPENERROR: if ( ierr /= 0 ) then
        write(*,10) trim(file_name)
      else
        write(*,20) trim(file_name)
      end if IFOPENERROR

    else if ( file_oc == CLOSE_FILE ) then
    ! close file
      close(unit=file_unit, iostat=ierr)

      IFCLOSEERROR: if(ierr.ne.0) then
        write(*,30) trim(file_name)
      else
        write(*,40) trim(file_name)
      end if IFCLOSEERROR

    end if IFOPENCLOSE
  end if IFMASTER

10  format('    | Error: cannot open ', 40A)
20  format('    | Opening ', 40A)
30  format('    | Error: cannot close ', 40A)
40  format('    | Closing ', 40A)

end subroutine g_FilesIO 
!---------------------------------------------------------------------------------

end module g_files
