!---------------------------------------------------------------------------------
! module g_restart
!> Read and write restart files.
module g_restart

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

!> Displacement between the two fluid velocity components
integer(kind=MPI_OFFSET_KIND)    :: f_disp
!> Displacement between particle position and particle velocity 
integer(kind=MPI_OFFSET_KIND)    :: p_disp
!> Fluid restart file type defined by MPI Darray
integer                          :: FLUID_RESTART_TYPE
!> Particle restart file type defined by MPI Darray
integer                          :: PARTICLE_RESTART_TYPE

public :: g_RestartPetscFluid
public :: g_RestartPetscParticle

contains
!---------------------------------------------------------------------------------
subroutine g_RestartPetscFluid ( readorwrite, ierr )

  use f_arrays,     only : u1_3w, u2_3w, u3_3w, f_ArraysVectorCreateRestart, &
                           da3w_rs, u1_3w_rs, u2_3w_rs, u3_3w_rs, &
                           restartcolour, restartcomm
  use g_parameters, only : MYRANK, MASTER, WRITE_FILE, READ_FILE, &
                           F_RESTART_FILE_READ, F_RESTART_FILE_WRITE, &
                           F_RESTART_SCALE

  implicit none

  integer,intent(in)              :: readorwrite
  PetscErrorCode,intent(inout)           :: ierr
  integer(kind=MPI_OFFSET_KIND)   :: disp=0
  integer                         :: local_size
  integer                         :: fluid_restart_file

  PetscViewer                     :: fluidviewer

  PetscScalar                     :: sum_u1, sum_u2, sum_u3

  IFREADWRITE: if ( readorwrite == READ_FILE ) then

    IFSCALED: if ( F_RESTART_SCALE == 1 ) then

      call PetscPrintf( PETSC_COMM_WORLD, '    | creating viewer for parallel fluid I/O \n', ierr )

      !> Create viewer
      call PetscViewerCreate ( PETSC_COMM_WORLD, fluidviewer, ierr )

      !> Set type to binary
      call PetscViewerSetType ( fluidviewer, PETSCVIEWERBINARY, ierr )

      !> Don't use info file
      call PetscViewerBinarySetSkipInfo ( fluidviewer, PETSC_TRUE, ierr )

      !> Set file to read mode
      call PetscViewerFileSetMode ( fluidviewer, FILE_MODE_READ, ierr )

      !> Use MPI parallel I/O
      call PetscViewerBinarySetUseMPIIO ( fluidviewer, PETSC_TRUE, ierr )

      !> Set file name to the one specified in the input file
      call PetscViewerFileSetName ( fluidviewer, trim(F_RESTART_FILE_READ), ierr ) 

      !> Load u1 velocity component from file
      call VecLoad ( u1_3w, fluidviewer, ierr )

      !> Load u2 velocity component from file
      call VecLoad ( u2_3w, fluidviewer, ierr )

      !> Load u3 velocity component from file
      call VecLoad ( u3_3w, fluidviewer, ierr )

      !> Compute the sum of u1 and u2 as a quick check if we still have the same array in the restart file
      call VecSum ( u1_3w, sum_u1, ierr )  
      call VecSum ( u2_3w, sum_u2, ierr )  
      call VecSum ( u3_3w, sum_u3, ierr )  

    else

      !> Create distributed array for loading the restart file only on relevant processes.
      !! Define communicator to identify these.
      call f_ArraysVectorCreateRestart ( ierr )

      IFRESTARTCOMM: if ( restartcolour == 1 ) then

        call PetscPrintf( restartcomm, '    | creating viewer for scaled parallel fluid I/O \n', ierr )

        !> Create viewer
        call PetscViewerCreate ( restartcomm, fluidviewer, ierr )

        !> Set type to binary
        call PetscViewerSetType ( fluidviewer, PETSCVIEWERBINARY, ierr )

        !> Don't use info file
        call PetscViewerBinarySetSkipInfo ( fluidviewer, PETSC_TRUE, ierr )

        !> Set file to read mode
        call PetscViewerFileSetMode ( fluidviewer, FILE_MODE_READ, ierr )

        !> Use MPI parallel I/O
        call PetscViewerBinarySetUseMPIIO ( fluidviewer, PETSC_TRUE, ierr )

        !> Set file name to the one specified in the input file
        call PetscViewerFileSetName ( fluidviewer, trim(F_RESTART_FILE_READ), ierr ) 

        !> Load u1 velocity component from file
        call VecLoad ( u1_3w_rs, fluidviewer, ierr )

        !> Load u2 velocity component from file
        call VecLoad ( u2_3w_rs, fluidviewer, ierr )

        !> Load u3 velocity component from file
        call VecLoad ( u3_3w_rs, fluidviewer, ierr )

        !> Compute the sum of u1 and u2 as a quick check if we still have the same array in the restart file
        call VecSum ( u1_3w_rs, sum_u1, ierr )  
        call VecSum ( u2_3w_rs, sum_u2, ierr )  
        call VecSum ( u3_3w_rs, sum_u3, ierr )  

      end if IFRESTARTCOMM

    end if IFSCALED

    IFMASTERREAD: if ( MYRANK == MASTER ) then
      write(*,10) sum_u1, sum_u2, sum_u3
    end if IFMASTERREAD

  else if ( readorwrite == WRITE_FILE ) then

    call PetscPrintf( PETSC_COMM_WORLD, '    | creating viewer for parallel fluid I/O \n', ierr )

    !> Create viewer
    call PetscViewerCreate ( PETSC_COMM_WORLD, fluidviewer, ierr )

    !> Set type to binary
    call PetscViewerSetType ( fluidviewer, PETSCVIEWERBINARY, ierr )

    !> Don't use info file
    call PetscViewerBinarySetSkipInfo ( fluidviewer, PETSC_TRUE, ierr )

    !> Set file to write mode
    call PetscViewerFileSetMode ( fluidviewer, FILE_MODE_WRITE, ierr )

    !> Use MPI parallel I/O
    call PetscViewerBinarySetUseMPIIO ( fluidviewer, PETSC_TRUE, ierr )

    !> Set file name to the one specified in the input file
    call PetscViewerFileSetName ( fluidviewer, trim(F_RESTART_FILE_WRITE), ierr ) 

    !> Write u1 velocity component to file
    call VecView ( u1_3w, fluidviewer, ierr )

    !> Write u2 velocity component to file
    call VecView ( u2_3w, fluidviewer, ierr )

    !> Write u3 velocity component to file
    call VecView ( u3_3w, fluidviewer, ierr )

    !> Compute the sum of u1 and u2 as a quick check if we still have the same array in the restart file
    call VecSum ( u1_3w, sum_u1, ierr )  
    call VecSum ( u2_3w, sum_u2, ierr )  
    call VecSum ( u3_3w, sum_u3, ierr )  

    IFMASTERWRITE: if ( MYRANK == MASTER ) then
      write(*,10) sum_u1, sum_u2, sum_u3
    end if IFMASTERWRITE

  end if IFREADWRITE

  call PetscViewerDestroy ( fluidviewer, ierr )

  10  format('    | Sum of u1: ', f10.5, ' Sum of u2: ', f10.5, ' Sum of u3: ', f10.5)

end subroutine g_RestartPetscFluid
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine g_RestartPetscParticle ( readorwrite, ierr )

  use p_arrays,     only : p_restart_xp1, p_restart_xp2, p_restart_xp3, &
                           p_restart_xp01, p_restart_xp02, p_restart_xp03, &
                           p_restart_up1, p_restart_up2, p_restart_up3, &
                           p_restart_dp, p_restart_status
  use g_parameters, only : MYRANK, MASTER, WRITE_FILE, READ_FILE, &
                           P_RESTART_FILE_READ, P_RESTART_FILE_WRITE, &
                           P_KEEP_INITIAL_POSITION, YES

  implicit none

  integer,intent(in)              :: readorwrite
  PetscErrorCode,intent(inout)           :: ierr
  integer(kind=MPI_OFFSET_KIND)   :: disp=0
  integer                         :: local_size
  integer                         :: particle_restart_file

  PetscViewer                     :: particleviewer

  PetscScalar                     :: sum_u1, sum_u2, sum_u3

  IFREADWRITE: if ( readorwrite == READ_FILE ) then

    call PetscPrintf( PETSC_COMM_WORLD, '    | creating viewer for parallel particle I/O \n', ierr )

    !> Create viewer
    call PetscViewerCreate ( PETSC_COMM_WORLD, particleviewer, ierr )

    !> Set type to binary
    call PetscViewerSetType ( particleviewer, PETSCVIEWERBINARY, ierr )

    !> Don't use info file
    call PetscViewerBinarySetSkipInfo ( particleviewer, PETSC_TRUE, ierr )

    !> Set file to read mode
    call PetscViewerFileSetMode ( particleviewer, FILE_MODE_READ, ierr )

    !> Use MPI parallel I/O
    call PetscViewerBinarySetUseMPIIO ( particleviewer, PETSC_TRUE, ierr )

    !> Set file name to the one specified in the input file
    call PetscViewerFileSetName ( particleviewer, trim(P_RESTART_FILE_READ), ierr ) 

!    !> Load u1 velocity component from file
!    call VecLoad ( u1_3w, particleviewer, ierr )

!    !> Load u2 velocity component from file
!    call VecLoad ( u2_3w, particleviewer, ierr )

!    !> Load u3 velocity component from file
!    call VecLoad ( u3_3w, particleviewer, ierr )

!    !> Compute the sum of u1 and u2 as a quick check if we still have the same array in the restart file
!    call VecSum ( u1_3w, sum_u1, ierr )  
!    call VecSum ( u2_3w, sum_u2, ierr )  
!    call VecSum ( u3_3w, sum_u3, ierr )  


!    IFMASTERREAD: if ( MYRANK == MASTER ) then
!      write(*,10) sum_u1, sum_u2, sum_u3
!    end if IFMASTERREAD

  else if ( readorwrite == WRITE_FILE ) then

    call PetscPrintf( PETSC_COMM_WORLD, '    | creating viewer for parallel particle I/O \n', ierr )

    !> Create viewer
    call PetscViewerCreate ( PETSC_COMM_WORLD, particleviewer, ierr )

    !> Set type to binary
    call PetscViewerSetType ( particleviewer, PETSCVIEWERBINARY, ierr )

    !> Don't use info file
    call PetscViewerBinarySetSkipInfo ( particleviewer, PETSC_TRUE, ierr )

    !> Set file to write mode
    call PetscViewerFileSetMode ( particleviewer, FILE_MODE_WRITE, ierr )

    !> Use MPI parallel I/O
    call PetscViewerBinarySetUseMPIIO ( particleviewer, PETSC_TRUE, ierr )

    !> Set file name to the one specified in the input file
    call PetscViewerFileSetName ( particleviewer, trim(P_RESTART_FILE_WRITE), ierr ) 

    !> Write x1 position component to file
    call VecView ( p_restart_xp1, particleviewer, ierr )

    !> Write x2 position component to file
    call VecView ( p_restart_xp2, particleviewer, ierr )

    !> Write x3 position component to file
    call VecView ( p_restart_xp3, particleviewer, ierr )

    !> Write u1 velocity component to file
    call VecView ( p_restart_up1, particleviewer, ierr )

    !> Write u2 velocity component to file
    call VecView ( p_restart_up2, particleviewer, ierr )

    !> Write u3 velocity component to file
    call VecView ( p_restart_up3, particleviewer, ierr )

    !> Write particle diameter to file
    call VecView ( p_restart_dp, particleviewer, ierr )

    !> Write particle number to file
    call VecView ( p_restart_status, particleviewer, ierr )

    IFKEEPINITIAL: if ( P_KEEP_INITIAL_POSITION == YES ) then

      !> Write x01 position component to file
      call VecView ( p_restart_xp01, particleviewer, ierr )

      !> Write x02 position component to file
      call VecView ( p_restart_xp02, particleviewer, ierr )

      !> Write x03 position component to file
      call VecView ( p_restart_xp03, particleviewer, ierr )

    end if IFKEEPINITIAL

!    !> Compute the sum of u1 and u2 as a quick check if we still have the same array in the restart file
!    call VecSum ( u1_3w, sum_u1, ierr )  
!    call VecSum ( u2_3w, sum_u2, ierr )  
!    call VecSum ( u3_3w, sum_u3, ierr )  

!    IFMASTERWRITE: if ( MYRANK == MASTER ) then
!      write(*,10) sum_u1, sum_u2, sum_u3
!    end if IFMASTERWRITE

  end if IFREADWRITE

  call PetscViewerDestroy ( particleviewer, ierr )

  10  format('    | Sum of u1: ', f10.5, ' Sum of u2: ', f10.5, ' Sum of u3: ', f10.5)

end subroutine g_RestartPetscParticle
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!subroutine g_RestartInit ( ierr )
!
!!    use g_parameters, only : MYRANK, MASTER, READ_FILE, WRITE_FILE
!
!    implicit none
!
!    integer,intent(inout)    :: ierr
!
!  call PetscPrintf( PETSC_COMM_WORLD, '==> Initialising f_restart module \n', ierr )
!
!  call PetscPrintf( PETSC_COMM_WORLD, '    | setting displacements for parallel I/O \n', ierr )
!  call g_RestartSetDisplacements ( ierr )
!
!end subroutine g_RestartInit
!!---------------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------------
!subroutine g_RestartSetDisplacements ( ierr )
!
!  use g_domain, only     : NODES_KI, NODES_KJ, NODES_KK
!  use g_parameters, only : P_NP_GLOBAL, MYRANK, MASTER, V_MPI, YES
!
!  implicit none
!
!  integer,intent(inout)  :: ierr
!
!  real(kind=C_DOUBLE)           :: testvariable
!  integer(kind=C_SIZE_T)        :: sizeoftestvariable
!  !> Variable for testing if MPI_OFFSET_KIND is sufficient for MPI I/O
!  integer(kind=MPI_OFFSET_KIND) :: offsettest
!
!  !> Get the size of the test variable in bytes
!  sizeoftestvariable = c_sizeof ( testvariable )
!
!  !> Calculate the global size of one fluid velocity component in bytes
!  f_disp = int(sizeoftestvariable,kind=MPI_OFFSET_KIND) * int(2,kind=MPI_OFFSET_KIND) * &
!           int(NODES_KI,kind=MPI_OFFSET_KIND) * int(NODES_KJ,kind=MPI_OFFSET_KIND) * &
!           int(NODES_KK,kind=MPI_OFFSET_KIND) 
!
!  !> Calculate the global size of the particle position arrays in bytes
!  p_disp = int(sizeoftestvariable,kind=MPI_OFFSET_KIND) * int(3,kind=MPI_OFFSET_KIND) * &
!           int(P_NP_GLOBAL,kind=MPI_OFFSET_KIND)
!
!  IFMASTER: if ( MYRANK == MASTER ) then
!    IFVERBOSE: if ( V_MPI == YES ) then
!      write(*,*) 'size of test variable: ', sizeoftestvariable, 'bytes.'
!      write(*,*) 'size of fluid displacement: ', f_disp, 'bytes.'
!      write(*,*) 'size of particle displacement: ', p_disp, 'bytes.'
!
!      write(*,*) 'Testing if MPI_OFFSET_KIND is sufficient on the system'
!      offsettest = int(2048,kind=MPI_OFFSET_KIND)**3 * int(2,kind=MPI_OFFSET_KIND) * &
!               int(sizeoftestvariable,kind=MPI_OFFSET_KIND)
!      write(*,*) 'size for 2048**3 grid: ', offsettest, 'bytes.'
!      offsettest = int(4096,kind=MPI_OFFSET_KIND)**3 * int(2,kind=MPI_OFFSET_KIND) * &
!               int(sizeoftestvariable,kind=MPI_OFFSET_KIND)
!      write(*,*) 'size for 4096**3 grid: ', offsettest, 'bytes.'
!      offsettest = int(8192,kind=MPI_OFFSET_KIND)**3 * int(2,kind=MPI_OFFSET_KIND) * &
!               int(sizeoftestvariable,kind=MPI_OFFSET_KIND)
!      write(*,*) 'size for 8192**3 grid: ', offsettest, 'bytes.'
!      offsettest = int(16384,kind=MPI_OFFSET_KIND)**3 * int(2,kind=MPI_OFFSET_KIND) * &
!               int(sizeoftestvariable,kind=MPI_OFFSET_KIND)
!      write(*,*) 'size for 16384**3 grid: ', offsettest, 'bytes.'
!    end if IFVERBOSE
!  end if IFMASTER
!
!end subroutine g_RestartSetDisplacements
!!---------------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------------
!subroutine g_RestartCreateDarrayFluid ( ierr )
!
!  use f_arrays,     only : da3w, i_width_3w, j_width_3w, k_width_3w
!  use g_parameters, only : MYRANK, NPROCS, P_TRACK_PART, YES, V_ALLOCS
!  use g_domain,     only : NODES_KK, NODES_KJ, NODES_KI
!
!  implicit none
!
!  integer,intent(inout)  :: ierr
!
!  integer                :: cdof
!  integer                :: arraydimension 
!  integer,dimension(1:4) :: gsizes
!  integer,dimension(1:4) :: distribs
!  integer,dimension(1:4) :: dargs
!  integer,dimension(1:4) :: psizes
!
!  !> real and imaginary part of one velocity component
!  cdof     = 2 
!  !> Define global array for fluid restart file
!  gsizes   = [ cdof, NODES_KK, NODES_KJ, NODES_KI ]
!
!  !> Describe how the global array is distributed
!  distribs(1) = MPI_Distribute_None
!  distribs(2) = MPI_Distribute_None
!  !> Distributed in j direction
!  distribs(3) = MPI_Distribute_Block
!  !> Distributed in i direction
!  distribs(4) = MPI_Distribute_Block
!  !> Distribution argument - just use the default
!  dargs       = MPI_Distribute_Dflt_Darg
!
!  !> Fluid array: 4-dimensional  
!  arraydimension = 4
!
!  !> Find out the size of the process grid in each dimension
!  psizes(1) = 1
!  call DMDAGetInfo ( da3w, PETSC_NULL_INTEGER, &                                        ! DMDA, dimension of DMDA
!                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! global dimensions M, N, P
!                     psizes(2), psizes(3), psizes(4), &                                 ! number of procs m, n, p
!                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &      ! number of dof per node, stencil width, ghost nodes bx
!                     PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr ) ! ghost nodes by, bz, stencil type, error code
!
!!  IFVERBOSE: if ( V_ALLOCS == YES ) then
!    write(*,*) 'local sizes: ', k_width_3w, j_width_3w, i_width_3w
!    write(*,*) 'gsizes: ', gsizes, 'psizes: ', psizes, 'arraydimension: ', arraydimension
!!  end if IFVERBOSE
!
!  call MPI_Type_Create_Darray ( NPROCS, &
!                                MYRANK, &            
!                                arraydimension, &     ! 4 for fluid, 2 for particle
!                                gsizes, &             ! array of sizes in each direction
!                                distribs, &           ! distribution of the array in each dimension
!                                dargs, &              ! distribution argument, can use default
!                                psizes, &             ! size of process grid in each dimension
!                                MPI_Order_Fortran, &  ! Fortran-type array ordering - we are using Fortran arrays
!                                MPI_Double, &         ! variable type
!                                FLUID_RESTART_TYPE, & ! Name of the type we are about to commit 
!                                ierr )
!
!  call MPI_Type_Commit ( FLUID_RESTART_TYPE, ierr )
!
!end subroutine g_RestartCreateDarrayFluid
!!---------------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------------
!subroutine g_RestartCreateDarrayParticles ( ierr )
!
!  use f_arrays,     only : da3w, i_width_3w, j_width_3w, k_width_3w
!  use g_parameters, only : MYRANK, NPROCS, P_TRACK_PART, YES
!  use p_control,    only : np_centre
!
!  implicit none
!
!  integer,intent(inout)  :: ierr
!
!  integer                :: cdof
!  integer                :: arraydimension 
!  integer,dimension(1:2) :: gsizes
!  integer,dimension(1:2) :: distribs
!  integer,dimension(1:2) :: dargs
!  integer,dimension(1:2) :: psizes
!
!  !> Define global array for fluid restart file
!  !> position and velocity components: total of 6
!  cdof     = 2 * 3
!  !> @todo All particles need to be in the centre arrays for this to be true
!  gsizes   = [ cdof, np_centre ]
!
!  !> Describe how the global array is distributed
!  distribs(1) = MPI_Distribute_None
!  distribs(2) = MPI_Distribute_Block
!  dargs       = MPI_Distribute_Dflt_Darg
!
!  !> We don't need detailed information about the exact 2-D processor grid.
!  !! For purposes of I/O the particles are simply distributed over
!  !! NPROCS processes on a 1-D array.
!  psizes = [ 1, NPROCS ]
!
!  !> Particle array: 2-dimensional  
!  arraydimension = 4
!
!  call MPI_Type_Create_Darray ( NPROCS, &
!                                MYRANK, &            
!                                arraydimension, &        ! 4 for fluid, 2 for particle
!                                gsizes, &                ! array of sizes in each direction
!                                distribs, &              ! distribution of the array in each dimension
!                                dargs, &                 ! distribution argument, can use default
!                                psizes, &                ! size of process grid in each dimension
!                                MPI_Order_Fortran, &     ! Fortran-type array ordering - we are using Fortran arrays
!                                MPI_Double, &            ! variable type 
!                                PARTICLE_RESTART_TYPE, & ! Name of the type we are about to commit 
!                                ierr )
!
!  call MPI_Type_Commit ( PARTICLE_RESTART_TYPE, ierr )
!
!end subroutine g_RestartCreateDarrayParticles
!!---------------------------------------------------------------------------------
!
!!---------------------------------------------------------------------------------
!subroutine g_RestartFluid ( readorwrite, ierr )
!
!  use f_arrays,     only : da3w, i_width_3w, j_width_3w, k_width_3w, &
!                           u1_3w, u2_3w, arr_u1_3w, arr_u2_3w, &
!                           i_min_3w, j_min_3w, k_min_3w, &
!                           i_max_3w, j_max_3w, k_max_3w
!  use g_parameters, only : MYRANK, NPROCS, MASTER, WRITE_FILE, READ_FILE, &
!                           F_RESTART_FILE_READ, F_RESTART_FILE_WRITE
!!    use f_force,      only : f_force_restart
!
!  implicit none
!
!  integer,intent(in)              :: readorwrite
!  integer,intent(inout)           :: ierr
!  integer(kind=MPI_OFFSET_KIND)   :: disp=0
!  integer                         :: local_size
!  integer                         :: fluid_restart_file
!!    integer                         :: pandora_comm
!  integer                         :: status(mpi_status_size)
!
!  integer                         :: i
!
!!  real(kind=C_DOUBLE),pointer,dimension(:,:,:,:) :: niftyarray
!
!  call PetscPrintf( PETSC_COMM_WORLD, '    | creating distributed array for parallel fluid I/O \n', ierr )
!  call g_RestartCreateDarrayFluid ( ierr )
!
!  IFREADWRITE: if ( readorwrite == READ_FILE ) then
!
!    !> open file and read arrays
!    call PetscPrintf( PETSC_COMM_WORLD, '    | reading fluid restart file \n', ierr )
!    call MPI_File_Open ( PETSC_COMM_WORLD, trim(F_RESTART_FILE_READ), &
!                         MPI_MODE_RDONLY, MPI_INFO_NULL, fluid_restart_file, ierr)
!
!    !> size of the local array identical for both velocity components
!    local_size = 2 * i_width_3w * j_width_3w * k_width_3w
!
!    !> Set displacement for first velocity component
!    disp = int(0,kind=MPI_OFFSET_KIND)
!
!    !> Set view for first velocity component
!    call MPI_File_Set_View ( fluid_restart_file, &
!                             disp, &
!                             MPI_Double, &
!                             FLUID_RESTART_TYPE, &
!                             "native", &
!                             MPI_INFO_NULL, &
!                             ierr )
!
!    !> Get write access to first velocity component
!    call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!
!    call mpi_file_read_all ( fluid_restart_file, &
!                             arr_u1_3w, &
!                             local_size, &
!                             MPI_Double, &
!                             status, &
!                             ierr )
!
!    !> Return velocity component to PETSc
!    call DMDAVecRestoreArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!
!    !> Set displacement for second velocity component
!    disp = f_disp 
!
!    !> Set view for second velocity component
!    call MPI_File_Set_View ( fluid_restart_file, &
!                             disp, &
!                             MPI_Double, &
!                             fluid_restart_type, &
!                             "native", &
!                             MPI_INFO_NULL, &
!                             ierr )
!
!    !> Get write access to second velocity component
!    call DMDAVecGetArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!
!    call mpi_file_read_all ( fluid_restart_file, &
!                             arr_u2_3w, &
!                             local_size, &
!                             MPI_Double, &
!                             status, &
!                             ierr )
!
!    !> Return velocity component to PETSc
!    call DMDAVecRestoreArrayF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!
!    call MPI_File_Close ( fluid_restart_file, ierr )
!    ! read forcing coefficients from file
!!    call f_force_restart(force_file_read,read_file,ierr)
!        !
!        !-----------------------------------------------
!
!  else if ( readorwrite == WRITE_FILE ) then
!
!    !> open file and write arrays
!    call PetscPrintf( PETSC_COMM_WORLD, '    | writing fluid restart files \n', ierr )
!    call MPI_File_Open ( PETSC_COMM_WORLD, trim(F_RESTART_FILE_WRITE)//'.u1', &
!                         MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fluid_restart_file, ierr)
!
!    !> size of the local array identical for both velocity components
!    local_size = 2 * i_width_3w * j_width_3w * k_width_3w
!
!    !> Set displacement for first velocity component
!    disp = 0!int(0,kind=MPI_OFFSET_KIND)
!
!    !> Set view for first velocity component
!    call MPI_File_Set_View ( fluid_restart_file, &
!                             disp, &
!                             MPI_Double, &
!                             FLUID_RESTART_TYPE, &
!                             "native", &
!                             MPI_INFO_NULL, &
!                             ierr )
!
!    !> Get read access to first velocity component
!    call DMDAVecGetArrayF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!
!!    write(*,*) 'arr_u1_3w: ', size(arr_u1_3w), shape(arr_u1_3w), local_size, arr_u1_3w(1,1,1,1), arr_u1_3w(2,1,1,1)
!!    write(*,*) 'arr_u1_3w: ', arr_u1_3w(0,k_min_3w,j_min_3w,i_min_3w), arr_u1_3w(1,k_max_3w,j_max_3w,i_max_3w), arr_u1_3w
!
!!    write(*,*) arr_u1_3w(0:1,k_min_3w:k_max_3w,j_min_3w:j_max_3w,i_min_3w:i_max_3w)
!    call mpi_file_write_all ( fluid_restart_file, &
!!    do i = 1, NPROCS, 1
! !     if ( MYRANK == i - 1 ) then
!  !      call mpi_file_write ( fluid_restart_file, &
!                             arr_u1_3w(0:1,k_min_3w:k_max_3w,j_min_3w:j_max_3w,i_min_3w:i_max_3w), &
!                             local_size, &
!                             MPI_Double, &
!                             status, &
!                             ierr )
!!        write(*,*) 'Written on rank ', MYRANK
!!        call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
!!      end if
!!    end do
!
!!    write(*,*) status
!
!!    call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
!
!    !> Return velocity component to PETSc
!!    call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!
!    call MPI_File_Close ( fluid_restart_file, ierr )
!    call MPI_File_Open ( PETSC_COMM_WORLD, trim(F_RESTART_FILE_WRITE)//'.u2', &
!                         MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, fluid_restart_file, ierr)
!
!
!    !> Set displacement for second velocity component
!!    disp = f_disp 
!    write(*,*) 'Displacement: ', disp
!
!    !> Set view for second velocity component
!    call MPI_File_Set_View ( fluid_restart_file, &
!                             disp, &
!                             MPI_Double, &
!                             fluid_restart_type, &
!                             "native", &
!                             MPI_INFO_NULL, &
!                             ierr )
!
!    !> Get read access to first velocity component
!    call DMDAVecGetArrayReadF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!
!!    niftyarray = [ arr_u1_3w, arr_u2_3w ]
!!    write(*,*) 'Joined velocity: ', niftyarray
!
!    call mpi_file_write_all ( fluid_restart_file, &
!                              arr_u2_3w, &
!                              local_size, &
!                              MPI_Double, &
!                              status, &
!                              ierr )
!
!    !> Return velocity component to PETSc
!    call DMDAVecRestoreArrayReadF90 ( da3w, u1_3w, arr_u1_3w, ierr )
!    call DMDAVecRestoreArrayReadF90 ( da3w, u2_3w, arr_u2_3w, ierr )
!
!    call MPI_File_Close ( fluid_restart_file, ierr )
!
!        ! write forcing coefficients to file
!!        call f_force_restart(force_file_write,write_file,ierr)
!        !
!        !-----------------------------------------------
!  end if IFREADWRITE
!
!
!end subroutine g_RestartFluid
!!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
end module g_restart
