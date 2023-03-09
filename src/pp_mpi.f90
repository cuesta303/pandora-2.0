module p_mpi

!> Empty if PETSc version <= 3.7
#include <g_petsc38.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
use g_petsc
!use petscsys
!use petscvec

use iso_c_binding

implicit none

!> Empty if PETSc version >= 3.8
#include <g_petsc37.h>

private

!-----------------------------------------------------------------
! variables 
!-----------------------------------------------------------------

integer,dimension(1:8)                                       :: p_sendneighbours
integer,dimension(1:8)                                       :: p_recvneighbours

integer                                                      :: p_neighbourminhalo
integer                                                      :: p_neighbourmaxhalo

integer,dimension(1:8)                                       :: p_yindexminneighbours
integer,dimension(1:8)                                       :: p_yindexmaxneighbours
integer,dimension(1:8)                                       :: p_zindexminneighbours
integer,dimension(1:8)                                       :: p_zindexmaxneighbours

integer                                                      :: p_yindexmin
integer                                                      :: p_yindexmax
integer                                                      :: p_zindexmin
integer                                                      :: p_zindexmax

real(kind=C_DOUBLE),dimension(1:8)                           :: p_yminneighbours
real(kind=C_DOUBLE),dimension(1:8)                           :: p_ymaxneighbours
real(kind=C_DOUBLE),dimension(1:8)                           :: p_zminneighbours
real(kind=C_DOUBLE),dimension(1:8)                           :: p_zmaxneighbours

real(kind=C_DOUBLE),save,allocatable,dimension(:,:),public   :: localbuffer
real(kind=C_DOUBLE),save,allocatable,dimension(:,:,:),public :: sendbuffer
real(kind=C_DOUBLE),save,allocatable,dimension(:,:,:),public :: recvbuffer

PetscInt,save,allocatable,dimension(:,:),public              :: intlocalbuffer
PetscInt,save,allocatable,dimension(:,:,:),public            :: intsendbuffer
PetscInt,save,allocatable,dimension(:,:,:),public            :: intrecvbuffer

integer,parameter,public                                     :: PART_IN_CENTRE     = 0  
integer,parameter,public                                     :: PART_INSIDE        = 4
integer,parameter,public                                     :: PART_OUTSIDE       = 3  
integer,parameter,public                                     :: PARTICLE_NOT_FOUND = -1

integer,parameter,public                                     :: HALO_UPDATE_REAL = 1
integer,parameter,public                                     :: HALO_UPDATE_STATUS = 2
integer,parameter,public                                     :: PART_MIGRATE_REAL = 3
integer,parameter,public                                     :: PART_MIGRATE_STATUS = 4

integer,parameter,public                                     :: PART_RESTART_XP     = 1
integer,parameter,public                                     :: PART_RESTART_UP     = 2
integer,parameter,public                                     :: PART_RESTART_DP     = 3
integer,parameter,public                                     :: PART_RESTART_STATUS = 4
integer,parameter,public                                     :: PART_RESTART_XP0    = 5

real(kind=C_DOUBLE),save                                     :: x1min
real(kind=C_DOUBLE),save                                     :: x1max
real(kind=C_DOUBLE),save                                     :: x2min
real(kind=C_DOUBLE),save                                     :: x2max
real(kind=C_DOUBLE),save                                     :: x3min
real(kind=C_DOUBLE),save                                     :: x3max

public :: p_MPIFindNeighbours
public :: p_MPIInitParticleDomain
public :: p_MPICountIncomingParticles
public :: p_MPIAllocateBufferMigrate
public :: p_MPIDeallocateBufferMigrate
public :: p_MPIFindLocalParticles
public :: p_MPIFindLocalParticlesLong
!public :: p_MPIFindParticlesInside
public :: p_MPIFindParticlesOutside
!public :: p_MPIFindParticlesInsideNeighbour
public :: p_MPICopyToSendBuffer
public :: p_MPIMigrateParticles
public :: p_MPICopyFromRecvBuffer
public :: p_MPICorrectPosition
public :: p_MPIRemeshPosition
public :: p_MPIRemeshInitial
public :: p_MPICopyToRestart 
 
contains
!---------------------------------------------------------------------------------
subroutine p_MPIFindNeighbours ( ierr )

  use g_domain,     only : YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                           YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME, &
                           LZ_MOVINGFRAME, LZ_HALO
  use g_constants,  only : ZERO
  use g_parameters, only : MYRANK, V_ALLOCS, YES, NODES_Y, NODES_Z
  use f_arrays,     only : da1r2w, z_min, z_max, y_min_2r, y_max_2r

  implicit none

  PetscErrorCode,intent(inout)      :: ierr

  integer                    :: i
  integer,dimension(1:27)    :: fulllistofneighbours

  integer                    :: mpi_status ( MPI_STATUS_SIZE )
  integer                    :: sendtag
  integer                    :: recvtag
  integer                    :: sendcount
  integer                    :: recvcount

  call DMDAGetNeighbors ( da1r2w, fulllistofneighbours, ierr )

  !> The array of neighbours is ordered in x, z, y directions.
  !! Our array is not parallelised in x direction, so we there
  !! are triplets of identical processes in this array. We need
  !! only the 2D description. 

  ! Ordering:
  ! (19 20 21) (22 23 24) (25 26 27)
  ! (10 11 12) (13 14 15) (16 17 18)
  ! ( 1  2  3) ( 4  5  6) ( 7  8  9)

  ! We use counterclockwise ordering:
  !  7  6  5
  !  8     4
  !  1  2  3

  !> Define counterclockwise send array
  p_sendneighbours(1) = fulllistofneighbours(1)
  p_sendneighbours(2) = fulllistofneighbours(4)
  p_sendneighbours(3) = fulllistofneighbours(7)
  p_sendneighbours(4) = fulllistofneighbours(16)
  p_sendneighbours(5) = fulllistofneighbours(25)
  p_sendneighbours(6) = fulllistofneighbours(22)
  p_sendneighbours(7) = fulllistofneighbours(19)
  p_sendneighbours(8) = fulllistofneighbours(10)

  !> Define counterclockwise receive array (simply shifted by 4 elements)
  p_recvneighbours(1) = p_sendneighbours(5)
  p_recvneighbours(2) = p_sendneighbours(6)
  p_recvneighbours(3) = p_sendneighbours(7)
  p_recvneighbours(4) = p_sendneighbours(8)
  p_recvneighbours(5) = p_sendneighbours(1)
  p_recvneighbours(6) = p_sendneighbours(2)
  p_recvneighbours(7) = p_sendneighbours(3)
  p_recvneighbours(8) = p_sendneighbours(4)

!  DOFINDDIMENSIONS: do i = 1, 8, 1

!    sendtag = MYRANK
!    recvtag = p_recvneighbours(i)
!    sendcount = 1
!    recvcount = 1

!    call MPI_Sendrecv ( YMIN_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_sendneighbours(i), sendtag, &
!                        p_yminneighbours(i), recvcount, MPI_DOUBLE, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( ZMIN_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_sendneighbours(i), sendtag, &
!                        p_zminneighbours(i), recvcount, MPI_DOUBLE, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( YMAX_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_sendneighbours(i), sendtag, &
!                        p_ymaxneighbours(i), recvcount, MPI_DOUBLE, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( ZMAX_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_sendneighbours(i), sendtag, &
!                        p_zmaxneighbours(i), recvcount, MPI_DOUBLE, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( y_min_2r, sendcount, MPI_INTEGER, p_sendneighbours(i), sendtag, &
!                        p_yindexminneighbours(i), recvcount, MPI_INTEGER, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( z_min,    sendcount, MPI_INTEGER, p_sendneighbours(i), sendtag, &
!                        p_zindexminneighbours(i), recvcount, MPI_INTEGER, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( y_max_2r, sendcount, MPI_INTEGER, p_sendneighbours(i), sendtag, &
!                        p_yindexmaxneighbours(i), recvcount, MPI_INTEGER, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!    call MPI_Sendrecv ( z_max,    sendcount, MPI_INTEGER, p_sendneighbours(i), sendtag, &
!                        p_zindexmaxneighbours(i), recvcount, MPI_INTEGER, p_recvneighbours(i), recvtag, &
!                        PETSC_Comm_World, mpi_status, ierr )

!  end do DOFINDDIMENSIONS

  !> Need to swap send and receive numbers in order to get the correct 
  !! coordinates of the neighbour nodes
  DOFINDDIMENSIONS: do i = 1, 8, 1

    sendtag = MYRANK
    recvtag = p_sendneighbours(i)
    sendcount = 1
    recvcount = 1

    call MPI_Sendrecv ( YMIN_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_recvneighbours(i), sendtag, &
                        p_yminneighbours(i), recvcount, MPI_DOUBLE, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( ZMIN_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_recvneighbours(i), sendtag, &
                        p_zminneighbours(i), recvcount, MPI_DOUBLE, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( YMAX_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_recvneighbours(i), sendtag, &
                        p_ymaxneighbours(i), recvcount, MPI_DOUBLE, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( ZMAX_MOVINGFRAME,    sendcount, MPI_DOUBLE, p_recvneighbours(i), sendtag, &
                        p_zmaxneighbours(i), recvcount, MPI_DOUBLE, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( y_min_2r, sendcount, MPI_INTEGER, p_recvneighbours(i), sendtag, &
                        p_yindexminneighbours(i), recvcount, MPI_INTEGER, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( z_min,    sendcount, MPI_INTEGER, p_recvneighbours(i), sendtag, &
                        p_zindexminneighbours(i), recvcount, MPI_INTEGER, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( y_max_2r, sendcount, MPI_INTEGER, p_recvneighbours(i), sendtag, &
                        p_yindexmaxneighbours(i), recvcount, MPI_INTEGER, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

    call MPI_Sendrecv ( z_max,    sendcount, MPI_INTEGER, p_recvneighbours(i), sendtag, &
                        p_zindexmaxneighbours(i), recvcount, MPI_INTEGER, p_sendneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

  end do DOFINDDIMENSIONS

  !> Find index corresponding to the upper physical limit (next node to the 'right').
  p_yindexmin = y_min_2r
  p_yindexmax = y_max_2r + 1
  p_zindexmin = z_min
  p_zindexmax = z_max + 1

  !> If same as length, set to zero.
  IFYEDGE: if ( p_yindexmax == NODES_Y ) then
    p_yindexmax = 0
  end if IFYEDGE

  IFZEDGE: if ( p_zindexmax == NODES_Z ) then
    p_zindexmax = 0
  end if IFZEDGE

  !> Find index corresponding to the upper physical limit (next node to the 'right').
  DOSHIFTNEIGHBOURS: do i = 1, 8, 1
    
    p_yindexmaxneighbours(i) = p_yindexmaxneighbours(i) + 1
    p_zindexmaxneighbours(i) = p_zindexmaxneighbours(i) + 1

    !> If same as length, set to zero.
    IFYEDGEN: if ( p_yindexmaxneighbours(i) == NODES_Y ) then
      p_yindexmaxneighbours(i) = 0
    end if IFYEDGEN

    IFZEDGEN: if ( p_zindexmaxneighbours(i) == NODES_Z ) then
      p_zindexmaxneighbours(i) = 0
    end if IFZEDGEN
  end do DOSHIFTNEIGHBOURS

  !> Find the neighbours for exchanging particles in the halos
!  DOFINDHALOSEND: do i = 1, 8, 1
!    IFSAMEY: if ( p_yindexminneighbours(i) == p_yindexmin ) then
!      IFZMINHALO: if ( p_zindexmaxneighbours(i) == p_zindexmin ) then
!        p_neighbourminhalo = p_sendneighbours(i)
!        p_neighbourmaxhalo = p_sendneighbours(i)
!      end if IFZMINHALO
!      IFZMAXHALO: if ( p_zindexminneighbours(i) == p_zindexmax ) then
!        p_neighbourmaxhalo = p_sendneighbours(i)
!        p_neighbourminhalo = p_sendneighbours(i)
!      end if IFZMAXHALO
!    end if IFSAMEY
!  end do DOFINDHALOSEND

  IFVERBOSE: if ( V_ALLOCS == YES ) then

    write(*,*) 'Send:    ', MYRANK, p_sendneighbours
    write(*,*) 'Receive: ', MYRANK, p_recvneighbours

    write(*,*) ' ymin: ', MYRANK, p_yminneighbours
    write(*,*) ' zmin: ', MYRANK, p_zminneighbours
    write(*,*) ' ymax: ', MYRANK, p_ymaxneighbours
    write(*,*) ' zmax: ', MYRANK, p_zmaxneighbours

    write(*,*) ' ymin: ', MYRANK, p_yindexminneighbours
    write(*,*) ' zmin: ', MYRANK, p_zindexminneighbours
    write(*,*) ' ymax: ', MYRANK, p_yindexmaxneighbours
    write(*,*) ' zmax: ', MYRANK, p_zindexmaxneighbours

  end if IFVERBOSE

end subroutine p_MPIFindNeighbours
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_MPIInitParticleDomain ( ierr )

  use g_constants, only : HALF
  use g_domain,    only : LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME, &
                          DX_MOVINGFRAME, DY_MOVINGFRAME, DZ_MOVINGFRAME

  implicit none

  PetscErrorCode,intent(inout) :: ierr

  x1min = - HALF * DX_MOVINGFRAME
  x1max = x1min + LX_MOVINGFRAME
  x2min = - HALF * DY_MOVINGFRAME
  x2max = x2min + LY_MOVINGFRAME
  x3min = - HALF * DZ_MOVINGFRAME
  x3max = x3min + LZ_MOVINGFRAME

end subroutine p_MPIInitParticleDomain
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_MPICountIncomingParticles ( npsendcentre, nprecvcentre, ierr )

  use g_parameters, only : MYRANK, NPROCS

  implicit none

  integer,dimension(1:8),intent(in)    :: npsendcentre
  integer,dimension(1:8),intent(out)   :: nprecvcentre
  PetscErrorCode,intent(inout)                :: ierr

  integer                              :: mpi_status ( MPI_STATUS_SIZE )
  integer                              :: i

  integer                              :: sendtag
  integer                              :: recvtag
  integer                              :: sendcount
  integer                              :: recvcount

!  write(*,*) ' My rank: ', MYRANK, ' Sending to: ', p_sendneighbours, &
!             ' Receiving from: ', p_recvneighbours

  DOEXCHANGE: do i = 1, 8, 1

    sendcount = 1
    recvcount = 1

    sendtag = MYRANK
    recvtag = p_recvneighbours(i)

    call MPI_Sendrecv ( npsendcentre(i), sendcount, MPI_INTEGER, p_sendneighbours(i), sendtag, &
                        nprecvcentre(i), recvcount, MPI_INTEGER, p_recvneighbours(i), recvtag, &
                        PETSC_Comm_World, mpi_status, ierr )

  end do DOEXCHANGE

end subroutine p_MPICountIncomingParticles
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_MPIAllocateBufferMigrate ( buffertype, dof, size_local, size_send, size_recv, ierr )

  implicit none

  integer,intent(in)                   :: buffertype
  integer,intent(in)                   :: dof
  integer,intent(in)                   :: size_local
  integer,dimension(1:8),intent(in)    :: size_send
  integer,dimension(1:8),intent(in)    :: size_recv
  integer                              :: size_send_max
  integer                              :: size_recv_max
  PetscErrorCode,intent(inout)                :: ierr

  select case ( buffertype )
    !-------------------------------
    case ( PART_MIGRATE_REAL )
      IFRL: if ( size_local > 0 ) then
        allocate ( localbuffer(1:dof,1:size_local), stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFRL

      size_send_max = maxval ( size_send)
      IFRS: if ( size_send_max > 0 ) then
        allocate ( sendbuffer(1:dof,1:size_send_max,1:8), stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFRS
        
      size_recv_max = maxval ( size_recv)
      IFRR: if ( size_recv_max > 0 ) then
        allocate ( recvbuffer(1:dof,1:size_recv_max,1:8), stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFRR

    !-------------------------------
    case ( PART_MIGRATE_STATUS )
      IFIL: if ( size_local > 0 ) then
        allocate ( intlocalbuffer(dof,size_local), stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFIL

      size_send_max = maxval ( size_send)
      IFIS: if ( size_send_max > 0 ) then
        allocate ( intsendbuffer(1:dof,1:size_send_max,1:8), stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFIS
        
      size_recv_max = maxval ( size_recv)
      IFIR: if ( size_recv_max > 0 ) then
        allocate ( intrecvbuffer(1:dof,1:size_recv_max,1:8), stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFIR
        
  end select

end subroutine p_MPIAllocateBufferMigrate
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_MPIDeallocateBufferMigrate ( buffertype, ierr )

  implicit none

  integer,intent(in)                   :: buffertype
  PetscErrorCode,intent(inout)                :: ierr

  select case ( buffertype )
    !-------------------------------
    case ( PART_MIGRATE_REAL )
      IFRL: if ( allocated ( localbuffer ) ) then
        deallocate ( localbuffer, stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFRL

      IFRS: if ( allocated ( sendbuffer ) ) then
        deallocate ( sendbuffer, stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFRS
        
      IFRR: if ( allocated ( recvbuffer ) ) then
        deallocate ( recvbuffer, stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFRR

    !-------------------------------
    case ( PART_MIGRATE_STATUS )
      IFIL: if ( allocated ( intlocalbuffer ) ) then
        deallocate ( intlocalbuffer, stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFIL

      IFIS: if ( allocated ( intsendbuffer ) ) then
        deallocate ( intsendbuffer, stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFIS
        
      IFIR: if ( allocated ( intrecvbuffer ) ) then
        deallocate ( intrecvbuffer, stat=ierr ) 
!        CHKERRQ ( ierr )
      end if IFIR
        
  end select

end subroutine p_MPIDeallocateBufferMigrate
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_MPICopyToSendBuffer
!> Distribute particles from an existing array to buffer arrays. Sort particles
!! by their owner ( MPI process ) and their location in either the centre or one
!! of the halos. This is achieved using the status array.
subroutine p_MPICopyToSendBuffer ( buffertype, dof_array, dof_status, &
                                   size_centre, size_send_centre, &
                                   ptr_centre, ptr_status_centre, size_centre_old, ierr )

  implicit none

  integer,intent(in)                         :: buffertype
  integer,intent(in)                         :: dof_array
  integer,intent(in)                         :: dof_status
  integer,intent(in)                         :: size_centre
  integer,dimension(1:8),intent(in)          :: size_send_centre
  integer                                    :: size_send_max
  type(C_PTR),intent(in)                     :: ptr_centre
  type(C_PTR),intent(in)                     :: ptr_status_centre
  integer,intent(in)                         :: size_centre_old
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: np
  integer                                    :: owner
  integer                                    :: npcentre
  integer,dimension(1:8)                     :: npsendcentre

  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particlearray_centre
  integer,pointer,dimension(:,:)             :: particlestatus_centre

  !> Set particle index to start position, sorted by halo_min, centre and halo_max, for both centre and send buffers
  npcentre      = 0
  npsendcentre  = 0

!  write(*,*) 'Sizes for sections: ', size_halo_min, size_centre, size_halo_max
!  write(*,*) 'Initial particle index for local particles: ', nphalomin, npcentre, nphalomax !, &
!             ' Initial particle indices for send particles: ', npsendhalomin, npsendcentre, npsendhalomax

!> @todo The same loop is repeated three times. Calling a subroutine for this would be more transparent.

  select case ( buffertype )
    !-------------------------------
    case ( PART_MIGRATE_REAL )

!      write(*,*) 'sendbuffer before copying: ', sendbuffer

      IFRCENTREEXIST: if ( size_centre_old > 0 ) then

        !> Associate pointer arrays with C pointers
        call c_f_pointer ( ptr_centre, particlearray_centre, [dof_array, size_centre_old])
        call c_f_pointer ( ptr_status_centre, particlestatus_centre, [dof_status, size_centre_old])
!        if ( dof_array == 1) write(*,*) 'centre before: ', particlearray_centre

        DORC: do np = 1, size_centre_old
          owner = particlestatus_centre(2,np)
!          write(*,*) 'Centre - owner: ', owner
          IFRCL: if ( owner == 0 ) then
            npcentre = npcentre + 1
            localbuffer(:,npcentre) = particlearray_centre(:,np)
          else
            npsendcentre(owner) = npsendcentre(owner) + 1
            sendbuffer(:,npsendcentre(owner),owner) = particlearray_centre(:,np)
!            write(*,*) 'Copying ', particlearray_centre(:,np), npsendcentre(owner), owner
          end if IFRCL
        end do DORC
      end if IFRCENTREEXIST


!      if ( dof_array == 1) write(*,*) 'sendbuffer after copying: ', sendbuffer
!      write(*,*) 'sendbuffer after copying: ', sendbuffer

    !-------------------------------
    case ( PART_MIGRATE_STATUS )
      IFICENTREEXIST: if ( size_centre_old > 0 ) then
        !> Associate pointer arrays with C pointers
        call c_f_pointer ( ptr_status_centre, particlestatus_centre, [dof_status, size_centre_old])
        DOIC: do np = 1, size_centre_old
          owner = particlestatus_centre(2,np)
          IFICL: if ( owner == 0 ) then
            npcentre = npcentre + 1
            intlocalbuffer(:,npcentre) = particlestatus_centre(:,np)
          else
            npsendcentre(owner) = npsendcentre(owner) + 1
            intsendbuffer(:,npsendcentre(owner),owner) = particlestatus_centre(:,np)
          end if IFICL
        end do DOIC
      end if IFICENTREEXIST

  end select

!  write(*,*) 'Particle index for local particles after copying: ', nphalomin, npcentre, nphalomax !, &

end subroutine p_MPICopyToSendBuffer
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine p_MPIMigrateParticles ( buffertype, dof, size_send, size_recv, ierr )

  use g_parameters, only : MYRANK, V_64BIT, YES

  implicit none

  integer,intent(in)                   :: buffertype
  integer,intent(in)                   :: dof
  integer,dimension(1:8),intent(in)    :: size_send
  integer,dimension(1:8),intent(in)    :: size_recv
  PetscErrorCode,intent(inout)                :: ierr

  integer                              :: i
  integer                              :: sendsize
  integer                              :: recvsize
  integer                              :: request
  integer                              :: status(MPI_STATUS_SIZE)

  request = MPI_REQUEST_NULL

  select case ( buffertype )
    !-------------------------------
    case ( PART_MIGRATE_REAL )
    ! loop to do the MPI
      DORMIGRATE: do i = 1, 8, 1
        !> send if there are particles
        IFRSEND: if ( size_send(i) > 0 ) then
          sendsize = dof * size_send(i)
          !> non-blocking send of the relevant buffer section
!          if ( dof == 1 ) write(*,*) 'Sending ', sendsize, ' from ', MYRANK, ' to ', p_sendneighbours(i)
!          if ( dof == 1 ) write(*,*) 'Rank ', MYRANK, 'Send buffer: ', sendbuffer(:,1:size_send(i),i) 
          call MPI_Isend(sendbuffer(:,1:size_send(i),i), sendsize, MPI_Double, &
                         p_sendneighbours(i), MYRANK, PETSC_COMM_WORLD, request, ierr)
        end if IFRSEND
    
        !> receive if there are particles expected
        IFRRECV: if ( size_recv(i) > 0 ) then
          recvsize = dof * size_recv(i)
          !> blocking receive of the relevant buffer section
!          if ( dof == 1 ) write(*,*) 'Receiving ', recvsize, ' from ', p_recvneighbours(i), ' to ', MYRANK
          call MPI_Recv(recvbuffer(:,1:size_recv(i),i), recvsize, MPI_Double, &
                        p_recvneighbours(i), p_recvneighbours(i), PETSC_COMM_WORLD, status, ierr)
!          if ( dof == 1 ) write(*,*) 'Rank ', MYRANK, 'Receive buffer: ', recvbuffer(:,1:size_recv(i),i)
          !> MPI_Wait just in case. Hopefully we don't need this.
          call MPI_Wait(request,status,ierr)
        end if IFRRECV

      end do DORMIGRATE

    !-------------------------------
    case ( PART_MIGRATE_STATUS )
    ! loop to do the MPI

      IF64BIT: if ( V_64BIT == YES ) then
        DOIMIGRATE64: do i = 1, 8, 1
          !> send if there are particles
          IFISEND64: if ( size_send(i) > 0 ) then
            sendsize = dof * size_send(i)
            !> non-blocking send of the relevant buffer section
            call MPI_Isend(intsendbuffer(:,1:size_send(i),i), sendsize, MPI_Integer8, &
                           p_sendneighbours(i), MYRANK, PETSC_COMM_WORLD, request, ierr)
          end if IFISEND64
    
          !> receive if there are particles expected
          IFIRECV64: if ( size_recv(i) > 0 ) then
            recvsize = dof * size_recv(i)
            !> blocking receive of the relevant buffer section
            call MPI_Recv(intrecvbuffer(:,1:size_recv(i),i), recvsize, MPI_Integer8, &
                          p_recvneighbours(i), p_recvneighbours(i), PETSC_COMM_WORLD, status, ierr)
            !> MPI_Wait just in case. Hopefully we don't need this.
            call MPI_Wait(request,status,ierr)
          end if IFIRECV64

        end do DOIMIGRATE64
      else
        DOIMIGRATE32: do i = 1, 8, 1
          !> send if there are particles
          IFISEND32: if ( size_send(i) > 0 ) then
            sendsize = dof * size_send(i)
            !> non-blocking send of the relevant buffer section
            call MPI_Isend(intsendbuffer(:,1:size_send(i),i), sendsize, MPI_Integer, &
                           p_sendneighbours(i), MYRANK, PETSC_COMM_WORLD, request, ierr)
          end if IFISEND32
    
          !> receive if there are particles expected
          IFIRECV32: if ( size_recv(i) > 0 ) then
            recvsize = dof * size_recv(i)
            !> blocking receive of the relevant buffer section
            call MPI_Recv(intrecvbuffer(:,1:size_recv(i),i), recvsize, MPI_Integer, &
                          p_recvneighbours(i), p_recvneighbours(i), PETSC_COMM_WORLD, status, ierr)
            !> MPI_Wait just in case. Hopefully we don't need this.
            call MPI_Wait(request,status,ierr)
          end if IFIRECV32

        end do DOIMIGRATE32
      end if IF64BIT

  end select

end subroutine p_MPIMigrateParticles
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
! subroutine p_MPICopyFromRecvBuffer
!> Distribute particles from an existing array to buffer arrays. Sort particles
!! by their owner ( MPI process ) and their location in either the centre or one
!! of the halos. This is achieved using the status array.
subroutine p_MPICopyFromRecvBuffer ( buffertype, dof, size_centre, size_recv_centre, &
                                     ptr_centre, size_centre_new, ierr )

  implicit none

  integer,intent(in)                         :: buffertype
  integer,intent(in)                         :: dof
  integer,intent(in)                         :: size_centre
  integer,dimension(1:8),intent(in)          :: size_recv_centre
  integer                                    :: size_recv_max
  type(C_PTR),intent(in)                     :: ptr_centre
  integer,intent(in)                         :: size_centre_new
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: i
  integer                                    :: np
  integer                                    :: nparray
  integer                                    :: npmin
  integer                                    :: npmax

  real(kind=C_DOUBLE),pointer,dimension(:,:) :: particlearray_centre
  integer,pointer,dimension(:,:)             :: particlestatus_centre

  select case ( buffertype )
    !-------------------------------
    case ( PART_MIGRATE_REAL )
    ! Copy from local buffer to halo_min
    ! Copy from recv buffer to halo_min
      nparray = 0

!      if ( dof == 1) write(*,*) 'recvbuffer: ', recvbuffer

    ! Copy from local buffer to centre
    ! Copy from recv buffer to centre
      nparray = 0

      IFRCENTREEXIST: if ( size_centre_new > 0 ) then
        call c_f_pointer ( ptr_centre, particlearray_centre, [dof, size_centre_new])

        IFRCENTRELOCALEXIST: if ( size_centre > 0 ) then
          !> Copy from local buffer
          npmin = 1 
          npmax = size_centre
          DORCOPYCENTRELOCAL: do np = npmin, npmax, 1
            nparray = nparray + 1
            particlearray_centre(:,nparray) = localbuffer(:,np)
          end do DORCOPYCENTRELOCAL
        end if IFRCENTRELOCALEXIST
        DORCENTRERECVBUFFER: do i = 1, 8, 1
          IFRCENTRERECVEXIST: if ( size_recv_centre(i) > 0 ) then
            !> Copy from recv buffer
            npmin = 1
            npmax = size_recv_centre(i)
            DORCOPYCENTRERECV: do np = npmin, npmax, 1
              nparray = nparray + 1
              particlearray_centre(:,nparray) = recvbuffer(:,np, i)
            end do DORCOPYCENTRERECV
          end if IFRCENTRERECVEXIST
        end do DORCENTRERECVBUFFER
!        if ( dof == 1) write(*,*) 'centre after: ', particlearray_centre
      end if IFRCENTREEXIST


    !-------------------------------
    case ( PART_MIGRATE_STATUS )
    ! Copy from local buffer to halo_min
    ! Copy from recv buffer to halo_min
      nparray = 0

      IFICENTREEXIST: if ( size_centre_new > 0 ) then
        call c_f_pointer ( ptr_centre, particlestatus_centre, [dof, size_centre_new])
        IFICENTRELOCALEXIST: if ( size_centre > 0 ) then
         !> Copy from local buffer
          npmin = 1
          npmax = size_centre
          DOICOPYCENTRELOCAL: do np = npmin, npmax, 1
            nparray = nparray + 1
            particlestatus_centre(:,nparray) = intlocalbuffer(:,np)
!             write(*,*) 'particlestatus_centre: ', particlestatus_centre(:,nparray)
          end do DOICOPYCENTRELOCAL
        end if IFICENTRELOCALEXIST
        DOICENTRERECVBUFFER: do i = 1, 8, 1
          IFICENTRERECVEXIST: if ( size_recv_centre(i) > 0 ) then
            !> Copy from recv buffer
            npmin = 1
            npmax = size_recv_centre(i)
            DOICOPYCENTRERECV: do np = npmin, npmax, 1
              nparray = nparray + 1
              particlestatus_centre(:,nparray) = intrecvbuffer(:,np, i)
            end do DOICOPYCENTRERECV
          end if IFICENTRERECVEXIST
        end do DOICENTRERECVBUFFER
      end if IFICENTREEXIST

  end select

end subroutine p_MPICopyFromRecvBuffer
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
! subroutine p_MPICopyToRestart
!> Distribute particles from an existing array to buffer arrays. Sort particles
!! by their owner ( MPI process ) and their location in either the centre or one
!! of the halos. This is achieved using the status array.
subroutine p_MPICopyToRestart ( restartarray, size_centre, size_recv_centre, &
                                size_total_new, ierr )

  use g_parameters, only : V_MPI, YES, P_KEEP_INITIAL_POSITION

  use p_arrays, only :     p_restart_xp1, p_restart_xp2, p_restart_xp3, &
                           p_restart_xp01, p_restart_xp02, p_restart_xp03, &
                           p_restart_up1, p_restart_up2, p_restart_up3, &
                           p_restart_dp, p_restart_status, &
                           arr_p_restart_xp1, arr_p_restart_xp2, &
                           arr_p_restart_xp3, arr_p_restart_xp01, &
                           arr_p_restart_xp02, arr_p_restart_xp03, &
                           arr_p_restart_up1, arr_p_restart_up2, &
                           arr_p_restart_up3, arr_p_restart_dp, &
                           arr_p_restart_status

  implicit none

  integer,intent(in)                         :: restartarray
  integer,intent(in)                         :: size_centre
  integer,dimension(1:8),intent(in)          :: size_recv_centre
  integer                                    :: size_recv_max
  integer,intent(in)                         :: size_total_new
  PetscErrorCode,intent(inout)                      :: ierr

  integer                                    :: i
  integer                                    :: np
  integer                                    :: nparray
  integer                                    :: npmin
  integer                                    :: npmax

  IFLOCALEXIST: if ( size_total_new > 0 ) then

    !> Get arrays from PETSc vectors
    select case ( restartarray )
      !-------------------------------
      case ( PART_RESTART_XP )
        !> Get arrays for copying from buffer
        call VecGetArrayF90 ( p_restart_xp1, arr_p_restart_xp1, ierr )
        call VecGetArrayF90 ( p_restart_xp2, arr_p_restart_xp2, ierr )
        call VecGetArrayF90 ( p_restart_xp3, arr_p_restart_xp3, ierr )
      !-------------------------------
      case ( PART_RESTART_XP0 )
        !> Get arrays for copying from buffer
        call VecGetArrayF90 ( p_restart_xp01, arr_p_restart_xp01, ierr )
        call VecGetArrayF90 ( p_restart_xp02, arr_p_restart_xp02, ierr )
        call VecGetArrayF90 ( p_restart_xp03, arr_p_restart_xp03, ierr )
      !-------------------------------
      case ( PART_RESTART_UP )
        !> Get arrays for copying from buffer
        call VecGetArrayF90 ( p_restart_up1, arr_p_restart_up1, ierr )
        call VecGetArrayF90 ( p_restart_up2, arr_p_restart_up2, ierr )
        call VecGetArrayF90 ( p_restart_up3, arr_p_restart_up3, ierr )
      !-------------------------------
      case ( PART_RESTART_DP )
        !> Get arrays for copying from buffer
        call VecGetArrayF90 ( p_restart_dp, arr_p_restart_dp, ierr )
      !-------------------------------
      case ( PART_RESTART_STATUS )
        !> Get arrays for copying from buffer
        call VecGetArrayF90 ( p_restart_status, arr_p_restart_status, ierr )
    end select

    nparray = 0

    IFRCENTRELOCALEXIST: if ( size_centre > 0 ) then
      !> Copy from local buffer
      npmin = 1 
      npmax = size_centre
      select case ( restartarray )
        !-------------------------------
        case ( PART_RESTART_XP )
          DORCOPYCENTRELOCALX: do np = npmin, npmax, 1
            nparray = nparray + 1
            arr_p_restart_xp1(nparray) = localbuffer(1,np)
            arr_p_restart_xp2(nparray) = localbuffer(2,np)
            arr_p_restart_xp3(nparray) = localbuffer(3,np)
          end do DORCOPYCENTRELOCALX
        !-------------------------------
        case ( PART_RESTART_XP0 )
          DORCOPYCENTRELOCALX0: do np = npmin, npmax, 1
            nparray = nparray + 1
            arr_p_restart_xp01(nparray) = localbuffer(1,np)
            arr_p_restart_xp02(nparray) = localbuffer(2,np)
            arr_p_restart_xp03(nparray) = localbuffer(3,np)
          end do DORCOPYCENTRELOCALX0
        !-------------------------------
        case ( PART_RESTART_UP )
          DORCOPYCENTRELOCALU: do np = npmin, npmax, 1
            nparray = nparray + 1
            arr_p_restart_up1(nparray) = localbuffer(1,np)
            arr_p_restart_up2(nparray) = localbuffer(2,np)
            arr_p_restart_up3(nparray) = localbuffer(3,np)
          end do DORCOPYCENTRELOCALU
        !-------------------------------
        case ( PART_RESTART_DP )
          DORCOPYCENTRELOCALD: do np = npmin, npmax, 1
            nparray = nparray + 1
            arr_p_restart_dp(nparray) = localbuffer(1,np)
          end do DORCOPYCENTRELOCALD
        !-------------------------------
        case ( PART_RESTART_STATUS )
          DORCOPYCENTRELOCALS: do np = npmin, npmax, 1
            nparray = nparray + 1
            arr_p_restart_status(nparray) = real ( intlocalbuffer(1,np), kind=C_DOUBLE )
          end do DORCOPYCENTRELOCALS
      end select
    end if IFRCENTRELOCALEXIST
    DORCENTRERECVBUFFER: do i = 1, 8, 1
      IFRCENTRERECVEXIST: if ( size_recv_centre(i) > 0 ) then
        !> Copy from recv buffer
        npmin = 1
        npmax = size_recv_centre(i)
        select case ( restartarray )
          !-------------------------------
          case ( PART_RESTART_XP )
            DORCOPYCENTRERECVX: do np = npmin, npmax, 1
              nparray = nparray + 1
              arr_p_restart_xp1(nparray) = recvbuffer(1,np, i)
              arr_p_restart_xp2(nparray) = recvbuffer(2,np, i)
              arr_p_restart_xp3(nparray) = recvbuffer(3,np, i)
            end do DORCOPYCENTRERECVX
          !-------------------------------
          case ( PART_RESTART_XP0 )
            DORCOPYCENTRERECVX0: do np = npmin, npmax, 1
              nparray = nparray + 1
              arr_p_restart_xp01(nparray) = recvbuffer(1,np, i)
              arr_p_restart_xp02(nparray) = recvbuffer(2,np, i)
              arr_p_restart_xp03(nparray) = recvbuffer(3,np, i)
            end do DORCOPYCENTRERECVX0
          !-------------------------------
          case ( PART_RESTART_UP )
            DORCOPYCENTRERECVU: do np = npmin, npmax, 1
              nparray = nparray + 1
              arr_p_restart_up1(nparray) = recvbuffer(1,np, i)
              arr_p_restart_up2(nparray) = recvbuffer(2,np, i)
              arr_p_restart_up3(nparray) = recvbuffer(3,np, i)
            end do DORCOPYCENTRERECVU
          !-------------------------------
          case ( PART_RESTART_DP )
            DORCOPYCENTRERECVD: do np = npmin, npmax, 1
              nparray = nparray + 1
              arr_p_restart_dp(nparray) = recvbuffer(1,np, i)
            end do DORCOPYCENTRERECVD
          !-------------------------------
          case ( PART_RESTART_STATUS )
            DORCOPYCENTRERECVS: do np = npmin, npmax, 1
              nparray = nparray + 1
              arr_p_restart_status(nparray) = real ( intrecvbuffer(1,np, i), kind=C_DOUBLE )
            end do DORCOPYCENTRERECVS
        end select
      end if IFRCENTRERECVEXIST
    end do DORCENTRERECVBUFFER
!     if ( dof == 1) write(*,*) 'centre after: ', particlearray_centre

    !> Return arrays to PETSc vectors
    select case ( restartarray )
      !-------------------------------
      case ( PART_RESTART_XP )
        !> Return arrays
        call VecRestoreArrayF90 ( p_restart_xp1, arr_p_restart_xp1, ierr )
        call VecRestoreArrayF90 ( p_restart_xp2, arr_p_restart_xp2, ierr )
        call VecRestoreArrayF90 ( p_restart_xp3, arr_p_restart_xp3, ierr )
      !-------------------------------
      case ( PART_RESTART_XP0 )
        !> Return arrays
        call VecRestoreArrayF90 ( p_restart_xp01, arr_p_restart_xp01, ierr )
        call VecRestoreArrayF90 ( p_restart_xp02, arr_p_restart_xp02, ierr )
        call VecRestoreArrayF90 ( p_restart_xp03, arr_p_restart_xp03, ierr )
      !-------------------------------
      case ( PART_RESTART_UP )
        !> Return arrays
        call VecRestoreArrayF90 ( p_restart_up1, arr_p_restart_up1, ierr )
        call VecRestoreArrayF90 ( p_restart_up2, arr_p_restart_up2, ierr )
        call VecRestoreArrayF90 ( p_restart_up3, arr_p_restart_up3, ierr )
      !-------------------------------
      case ( PART_RESTART_DP )
        !> Return arrays
        call VecRestoreArrayF90 ( p_restart_dp, arr_p_restart_dp, ierr )
      !-------------------------------
      case ( PART_RESTART_STATUS )
        !> Return arrays
        call VecRestoreArrayF90 ( p_restart_status, arr_p_restart_status, ierr )
    end select

  end if IFLOCALEXIST

  IFVERBOSE: if ( V_MPI == YES ) then

    !> View PETSc vectors
    select case ( restartarray )
    !-------------------------------
      case ( PART_RESTART_XP )
        !> Return arrays
        call PetscPrintf( PETSC_COMM_WORLD, '    | Viewing xp vectors \n', ierr )
        call VecView ( p_restart_xp1, PETSC_VIEWER_STDOUT_WORLD, ierr )
        call VecView ( p_restart_xp2, PETSC_VIEWER_STDOUT_WORLD, ierr )
        call VecView ( p_restart_xp3, PETSC_VIEWER_STDOUT_WORLD, ierr )
      !-------------------------------
      case ( PART_RESTART_XP0 )
        !> Return arrays
        call PetscPrintf( PETSC_COMM_WORLD, '    | Viewing xp0 vectors \n', ierr )
        call VecView ( p_restart_xp01, PETSC_VIEWER_STDOUT_WORLD, ierr )
        call VecView ( p_restart_xp02, PETSC_VIEWER_STDOUT_WORLD, ierr )
        call VecView ( p_restart_xp03, PETSC_VIEWER_STDOUT_WORLD, ierr )
      !-------------------------------
      case ( PART_RESTART_UP )
        !> Return arrays
        call PetscPrintf( PETSC_COMM_WORLD, '    | Viewing up vectors \n', ierr )
        call VecView ( p_restart_up1, PETSC_VIEWER_STDOUT_WORLD, ierr )
        call VecView ( p_restart_up2, PETSC_VIEWER_STDOUT_WORLD, ierr )
        call VecView ( p_restart_up3, PETSC_VIEWER_STDOUT_WORLD, ierr )
      !-------------------------------
      case ( PART_RESTART_DP )
        !> Return arrays
        call PetscPrintf( PETSC_COMM_WORLD, '    | Viewing dp vector \n', ierr )
        call VecView ( p_restart_dp, PETSC_VIEWER_STDOUT_WORLD, ierr )
      !-------------------------------
      case ( PART_RESTART_STATUS )
        !> Return arrays
        call PetscPrintf( PETSC_COMM_WORLD, '    | Viewing status vector \n', ierr )
        call VecView ( p_restart_status, PETSC_VIEWER_STDOUT_WORLD, ierr )
    end select

  end if IFVERBOSE

end subroutine p_MPICopyToRestart
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
pure subroutine p_MPICorrectPosition ( xp )

  use g_domain,    only : LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME, &
                          XMIN_MOVINGFRAME, YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                          XMAX_MOVINGFRAME, YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME
  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(inout) :: xp

  !> If particle position is < min or >= max, shift.
  !! The >= max condition is necessary to be consistent
  !! with the routines that find particles. Otherwise
  !! with increasing number of particles it becomes more
  !! likely that some particles are exactly at the max
  !! position and can't be found.
  IFMIN1: if ( xp(1) < x1min ) then
!    DOMIN1: do while ( xp(1) < x1min )
      xp(1) = xp(1) + LX_MOVINGFRAME
!    end do DOMIN1
  end if IFMIN1
  IFMAX1: if ( xp(1) >= x1max ) then
!    DOMAX1: do while ( xp(1) >= x1max )
      xp(1) = xp(1) - LX_MOVINGFRAME
!    end do DOMAX1
  end if IFMAX1
  IFMIN2: if ( xp(2) < x2min ) then
!    DOMIN2: do while ( xp(2) < x2min )
      xp(2) = xp(2) + LY_MOVINGFRAME
!    end do DOMIN2
  end if IFMIN2
  IFMAX2: if ( xp(2) >= x2max ) then
!    DOMAX2: do while ( xp(2) >= x2max )
      xp(2) = xp(2) - LY_MOVINGFRAME
!    end do DOMAX2
  end if IFMAX2
  IFMIN3: if ( xp(3) < x3min ) then
!    DOMIN3: do while ( xp(3) < x3min )
      xp(3) = xp(3) + LZ_MOVINGFRAME
!    end do DOMIN3
  end if IFMIN3
  IFMAX3: if ( xp(3) > x3max ) then
!    DOMAX3: do while ( xp(3) >= x3max )
      xp(3) = xp(3) - LZ_MOVINGFRAME
!    end do DOMAX3
  end if IFMAX3

end subroutine p_MPICorrectPosition
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
pure subroutine p_MPIRemeshPosition ( xp )

  use g_parameters, only : B110
  use g_domain,     only : LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME, &
                           XMIN_MOVINGFRAME, YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                           XMAX_MOVINGFRAME, YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME
  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(inout) :: xp

  !> If particle position is < min or >= max, shift.
  !! The >= max condition is necessary to be consistent
  !! with the routines that find particles. Otherwise
  !! with increasing number of particles it becomes more
  !! likely that some particles are exactly at the max
  !! position and can't be found.

  xp(1) = xp(1) + xp(3) 

  IFMIN1: if ( xp(1) < x1min ) then
!    DOMIN1: do while ( xp(1) < x1min )
      xp(1) = xp(1) + LX_MOVINGFRAME
!    end do DOMIN1
  end if IFMIN1
  IFMAX1: if ( xp(1) >= x1max ) then
!    DOMAX1: do while ( xp(1) >= x1max )
      xp(1) = xp(1) - LX_MOVINGFRAME
!    end do DOMAX1
  end if IFMAX1

end subroutine p_MPIRemeshPosition
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
pure subroutine p_MPIRemeshInitial ( xp0, xp )

  use g_parameters, only : B110
  use g_domain,     only : LX_MOVINGFRAME, LY_MOVINGFRAME, LZ_MOVINGFRAME, &
                           XMIN_MOVINGFRAME, YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                           XMAX_MOVINGFRAME, YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME
  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(inout) :: xp0
  real(kind=C_DOUBLE),dimension(1:3),intent(inout) :: xp

  !> If particle position is < min or >= max, shift.
  !! The >= max condition is necessary to be consistent
  !! with the routines that find particles. Otherwise
  !! with increasing number of particles it becomes more
  !! likely that some particles are exactly at the max
  !! position and can't be found.

  xp0(1) = xp0(1) + xp0(3) 
  xp(1)  = xp(1)  + xp(3) 

  IFMIN1: if ( xp(1) < x1min ) then
!    DOMIN1: do while ( xp(1) < x1min )
      xp0(1) = xp0(1) + LX_MOVINGFRAME
!    end do DOMIN1
  end if IFMIN1
  IFMAX1: if ( xp(1) >= x1max ) then
!    DOMAX1: do while ( xp(1) >= x1max )
      xp0(1) = xp0(1) - LX_MOVINGFRAME
!    end do DOMAX1
  end if IFMAX1

end subroutine p_MPIRemeshInitial
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
pure function p_MPIFindLocalParticlesLong ( xp ) result ( particleposition )

  use g_domain,     only : YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                           YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME
  implicit none

  real(kind=C_LONG_DOUBLE),dimension(3),intent(in) :: xp

  integer                                          :: particleposition

  IFYMIN: if ( xp(2) >= real(YMIN_MOVINGFRAME,kind=C_LONG_DOUBLE) ) then
!  IFYMIN: if ( xp(2) >= x2min ) then
    IFYMAX: if ( xp(2) < real(YMAX_MOVINGFRAME,kind=C_LONG_DOUBLE) ) then
!    IFYMAX: if ( xp(2) < x2max ) then
      !> Check if on this node. Ignore halos.
      IFZMIN: if ( xp(3) >= real(ZMIN_MOVINGFRAME,kind=C_LONG_DOUBLE) ) then
!      IFZMIN: if ( xp(3) >= x3min ) then
        IFZMAX: if ( xp(3) < real(ZMAX_MOVINGFRAME,kind=C_LONG_DOUBLE) ) then
!        IFZMAX: if ( xp(3) < x3max ) then
          particleposition = PART_INSIDE
        else
          particleposition = PART_OUTSIDE
        end if IFZMAX
      end if IFZMIN
    end if IFYMAX
  end if IFYMIN

end function p_MPIFindLocalParticlesLong
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
pure function p_MPIFindLocalParticles ( xp ) result ( particleposition )

  use g_domain,     only : YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
                           YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME
  implicit none

  real(kind=C_DOUBLE),dimension(3),intent(in) :: xp

  integer                                     :: particleposition

  IFYMIN: if ( xp(2) >= YMIN_MOVINGFRAME ) then
!  IFYMIN: if ( xp(2) >= x2min ) then
    IFYMAX: if ( xp(2) < YMAX_MOVINGFRAME ) then
!    IFYMAX: if ( xp(2) < x2max ) then
      !> Check if on this node. Ignore halos.
      IFZMIN: if ( xp(3) >= ZMIN_MOVINGFRAME ) then
!      IFZMIN: if ( xp(3) >= x3min ) then
        IFZMAX: if ( xp(3) < ZMAX_MOVINGFRAME ) then
!        IFZMAX: if ( xp(3) < x3max ) then
          particleposition = PART_INSIDE
        else
          particleposition = PART_OUTSIDE
        end if IFZMAX
      end if IFZMIN
    end if IFYMAX
  end if IFYMIN

end function p_MPIFindLocalParticles
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
function p_MPIFindParticlesOutside ( xp ) result ( particleowner )

  implicit none

  real(kind=C_DOUBLE),dimension(1:3),intent(in) :: xp

  integer                                       :: particleowner
  integer                                       :: i

  !> Particleposition not yet known
  particleowner = PARTICLE_NOT_FOUND

  DONEIGHBOURS: do i = 1, 8, 1
    IFYMIN: if ( xp(2) >= p_yminneighbours(i) ) then
      IFYMAX: if ( xp(2) < p_ymaxneighbours(i) ) then
        IFZMIN: if ( xp(3) >= p_zminneighbours(i) ) then
          IFZMAX: if ( xp(3) < p_zmaxneighbours(i) ) then
            particleowner = i
          end if IFZMAX
        end if IFZMIN
      end if IFYMAX
    end if IFYMIN
  end do DONEIGHBOURS

  IFNOTFOUND: if ( particleowner == PARTICLE_NOT_FOUND ) then
    write(*,*) 'Particle lost.'
    DOLOST:  do i = 1, 8, 1
      write(*,*) 'Particle position: ', xp, ' y range: ', &
      p_yminneighbours(i), p_ymaxneighbours(i), &
      ' z range: ', p_zminneighbours(i), p_zmaxneighbours(i), &
      xp(2), '>', p_yminneighbours(i), (xp(2) > p_yminneighbours(i)), &
      xp(2), '<', p_ymaxneighbours(i), (xp(2) < p_ymaxneighbours(i)), &
      xp(3), '>', p_zminneighbours(i), (xp(3) > p_zminneighbours(i)), &
      xp(3), '<', p_zmaxneighbours(i), (xp(3) < p_zmaxneighbours(i))
     
    end do DOLOST
  end if IFNOTFOUND

end function p_MPIFindParticlesOutside
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!pure function p_MPIFindParticlesInside ( xp ) result ( particleposition )

!  use g_domain,     only : YMIN_MOVINGFRAME, ZMIN_MOVINGFRAME, &
!                           YMAX_MOVINGFRAME, ZMAX_MOVINGFRAME, LZ_HALO

!  implicit none

!  real(kind=C_DOUBLE),dimension(3),intent(in) :: xp

!  integer                                     :: particleposition

  !> Particle outside unless proven otherwise
!  particleposition = PART_OUTSIDE

!  IFYMIN: if ( xp(2) >= YMIN_MOVINGFRAME ) then
!    IFYMAX: if ( xp(2) < YMAX_MOVINGFRAME ) then
      !> check if in either of the three z arrays
!      IFZMIN: if ( xp(3) >= ZMIN_MOVINGFRAME ) then
!        IFZMAX: if ( xp(3) < ZMIN_MOVINGFRAME + LZ_HALO ) then
!          particleposition = PART_IN_HALO_MIN
!        else if ( xp(3) < ZMAX_MOVINGFRAME - LZ_HALO ) then
!          particleposition = PART_IN_CENTRE
!        else if ( xp(3) < ZMAX_MOVINGFRAME ) then
!          particleposition = PART_IN_HALO_MAX
!        end if IFZMAX
!      end if IFZMIN
!    end if IFYMAX
!  end if IFYMIN

!end function p_MPIFindParticlesInside
!---------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
end module p_mpi
