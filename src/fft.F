!
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
!******************************************************************************
! MODULE m_fft
! 3-D fast complex Fourier transform
! Written by J.D.Gale (July 1999)
! Modified by Rogeli Grima (June 2022)
!******************************************************************************
!
!   PUBLIC procedures available from this module:
! fft_init   : Initiallization of the FFT
! getFFTLims : Returns the limits of the arrays in the Transformed domain
! fft        : 3-D complex FFT
!
!   PUBLIC parameters, types, and variables available from this module:
! none
!
!******************************************************************************
!
!   USED module procedures:
! use sys,          only: die           ! Termination routine
! use alloc,        only: de_alloc      ! De-allocation routine
! use alloc,        only: re_alloc      ! Re-allocation routine
! use fft1d,        only: gpfa          ! 1-D FFT routine
! use fft1d,        only: setgpfa       ! Sets gpfa routine
! use m_timer,      only: timer_start   ! Start counting CPU time
! use m_timer,      only: timer_stop    ! Stop counting CPU time
!
!   USED MPI procedures
! use mpi_siesta
!
!   USED module parameters:
! use precision,    only: dp            ! Real double precision type
! use precision,    only: gp=>grid_p    ! Real precision type of mesh arrays
! use parallel,     only: Node, Nodes, ProcessorY
! use mesh,         only: nsm
!
!   EXTERNAL procedures used:
! gpfa    : 1D complex FFT
! setgpfa : Initializes gpfa
! timer   : CPU time counter
!
!******************************************************************************
#include "mpi_macros.f"

      MODULE m_fft

      use precision,    only : dp, grid_p
      use, intrinsic :: iso_c_binding, only : C_PTR, C_DOUBLE_COMPLEX
      use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
      use parallel,     only : Node, Nodes, processorY
      use moreMeshSubs, only : UNIFORM, getMeshBox
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
      use mesh,         only : nsm
      use gpfa_fft,     only : gpfa        ! 1-D FFT routine
      use gpfa_fft,     only : setgpfa=>setgpfa_check     ! Sets gpfa routine
      use m_timer,      only : timer_start ! Start counting CPU time
      use m_timer,      only : timer_stop ! Stop counting CPU time
#ifdef MPI
      use mpi_siesta
#endif
#ifdef SIESTA__FFTW
      use fftw3_siesta_wrapper
#endif

      implicit none

      PUBLIC :: fft
      PUBLIC :: fft_init
      PUBLIC :: getFFTLims
      PUBLIC :: fft_set_fftw
      PUBLIC :: fft_trigs_reset
      PUBLIC :: fft_serial_global

      PRIVATE ! Nothing is declared public beyond this point
      logical :: frsttime = .TRUE.

#ifdef SIESTA__FFTW
      type(C_PTR), target :: plan_fw(3)
      type(C_PTR), target :: plan_bw(3)
      logical             :: use_fftw = .true.
#endif
      integer :: maxtrigs = 256
      real(dp), pointer :: trigs(:,:) => null()
      integer :: OldMesh(3) = (/ 0, 0, 0 /)

      integer :: gmesh(3)
      integer :: ProcessorZ
      integer :: PY
      integer :: PZ
#ifdef MPI
      MPI_COMM_TYPE :: COMM_Y
      MPI_COMM_TYPE :: COMM_Z
#endif
      logical :: k_pin_x_distribution = .false.
      integer, pointer :: xdispl(:) => null()
      integer, pointer :: ydispl(:) => null()
      integer, pointer :: zdispl(:) => null()
      integer, pointer :: xzdisp(:) => null()
      integer, pointer :: yzdisp(:) => null()

      CONTAINS

!******************************************************************************
      subroutine fft_trigs_reset()
        implicit none

        call de_alloc( trigs, 'trigs', 'fft' )
      end subroutine fft_trigs_reset

      subroutine fft_set_fftw( use_fftw_in )
        !! Sets FFTW usage at runtime.
        implicit none
        logical, intent(in) :: use_fftw_in
#ifdef SIESTA__FFTW
        use_fftw = use_fftw_in
        if ( use_fftw ) then
          if ( node == 0 ) write(6,'(A)') "FFT: Using FFTW."
        else
          if ( node == 0 ) write(6,'(A)') "FFT: Using regular FFT."
        endif
#else
        if ( use_fftw_in .and. (node == 0)) write(6,'(A)')
     &     "FFT: FFTW support is not enabled in this compilation."
#endif
      end subroutine fft_set_fftw

      subroutine fft_init( nMesh, size1, size2,  pin_x_distribution )
C **********************************************************************
C Initiallize data structures to compute the FFT:
C   INTEGER PROCESSORZ : NODES/PROCESSORY
C   INTEGER PY, PZ     : INDICES of the current process in a 2-D grid of
C                        processors of size PROCESSORY by PROCESSORZ.
C   INTEGER COMM_Y     : global communicator for Y dimension
C   INTEGER COMM_Z     : global communicator for Z dimension
C   INTEGER XDISPL()   : (ProcessorY) cuts in the X direction
C   INTEGER YDISPL()   : (ProcessorY) cuts in the Y direction
C   INTEGER ZDISPL()   : (ProcessorZ) cuts in the Z direction
C   INTEGER XZDISP()   : (ProcessorZ) cuts in the X direction (if needed)
C   INTEGER YZDISP()   : (ProcessorZ) cuts in the Y direction (if needed)
C Initiallize the FFT plans
C Rogeli Grima (June 2022)
C ************ INPUT ***************************************************
C LOGICAL  pin_x_distribution :
C   true: To keep the original distribution (Pencil in X direction)
C   false: If we swap from X to Z (forward) or from Z to X (backward)
C INTEGER nMesh(3)  : Global mesh dimension
C INTEGER SIZE1     : returns the size needed for the first array
C INTEGER SIZE2     :  returns the size needed for the second array
C **********************************************************************
      implicit none
      logical :: pin_x_distribution                      ! Keep original distribution
      integer :: nMesh(3), size1, size2
      integer, pointer :: box(:,:,:)
      integer :: i, di, mo, s1, s2, s3
      integer :: ierr, maxmaxtrigs, ntrigs
      logical :: lredimension
#ifdef SIESTA__FFTW
      integer :: rank, n(1), howmany, inembed(1), istride, idist
      complex(C_DOUBLE_COMPLEX) :: in(1)
#endif
      gmesh(1:3) = nMesh(1:3)

      if (frsttime) then
#ifdef SIESTA__FFTW
        if ( .not. use_fftw ) then
          nullify(trigs)
          call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
        endif
#else
        nullify(trigs)
        call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
#endif
        ProcessorZ = NODES/ProcessorY
        PY = mod(Node,ProcessorY)+1
        PZ = Node/ProcessorY+1
#ifdef MPI
        call MPI_Comm_Split( MPI_Comm_World, PZ, PY, COMM_Y, ierr )
        call MPI_Comm_Split( MPI_Comm_World, PY, PZ, COMM_Z, ierr )
#endif
        call re_alloc( xdispl, 1, ProcessorY+1, 'xdispl', 'fft' )
        call re_alloc( ydispl, 1, ProcessorY+1, 'ydispl', 'fft' )
        call re_alloc( zdispl, 1, ProcessorZ+1, 'zdispl', 'fft' )
        call re_alloc( xzdisp, 1, ProcessorZ+1, 'xzdisp', 'fft' )
        yzdisp => xzdisp

        frsttime = .false.
      endif
#ifdef MPI
#else
        size1 = gMesh(1)*gMesh(2)*gMesh(3)
        size2 = 0
#endif
C
C     Recompute plans and displacements if there is a change in the grid size or
C     in the value of k_pin_x_distribution
C
      if ( any(OldMesh/=gmesh) .or.
     $     (pin_x_distribution .neqv. k_pin_x_distribution) ) then
C
C       Compute the limits of the differents pencils in each possible direction
C
        call getMeshBox( UNIFORM, box )
        di = gmesh(1)/ProcessorY
        mo = mod(gmesh(1),ProcessorY)
        xdispl(1) = 1
        do i= 1, ProcessorY
          xdispl(i+1) = xdispl(i) + di + merge(1,0,i<=mo)
        enddo
        ydispl(1) = 1
        do i= 1, ProcessorY
          ydispl(i+1) = box(2,2,i)*NSM+1
        enddo
        zdispl(1) = 1
        do i= 1, ProcessorZ
          zdispl(i+1) = box(2,3,(i-1)*ProcessorY+1)*NSM+1
        enddo
        if (pin_x_distribution) then
          di = gmesh(1)/ProcessorZ
          mo = mod(gmesh(1),ProcessorZ)
          xzdisp(1) = 1
          do i= 1, ProcessorZ
            xzdisp(i+1) = xzdisp(i) + di + merge(1,0,i<=mo)
          enddo
        else
          di = gmesh(2)/ProcessorZ
          mo = mod(gmesh(2),ProcessorZ)
          yzdisp(1) = 1
          do i= 1, ProcessorZ
            yzdisp(i+1) = yzdisp(i) + di + merge(1,0,i<=mo)
          enddo
        endif
#ifdef SIESTA__FFTW
        if ( use_fftw ) then

        if (any(OldMesh/=0)) then
          call fftw_destroy_plan( plan_fw(1) )
          call fftw_destroy_plan( plan_fw(2) )
          call fftw_destroy_plan( plan_fw(3) )
          call fftw_destroy_plan( plan_bw(1) )
          call fftw_destroy_plan( plan_bw(2) )
          call fftw_destroy_plan( plan_bw(3) )
        endif
        rank = 1
        n = (/ gmesh(1) /)
        howmany = (ydispl(py+1)-ydispl(py))*(zdispl(pz+1)-zdispl(pz))
        inembed = (/ gmesh(1) /)
        istride = 1
        idist   = gmesh(1)
        plan_fw(1) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_FORWARD, FFTW_ESTIMATE )
        plan_bw(1) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_BACKWARD, FFTW_ESTIMATE )
        n = (/ gmesh(2) /)
        howmany = (xdispl(py+1)-xdispl(py))
        inembed = (/ gmesh(2) /)
        istride = (xdispl(py+1)-xdispl(py))
        idist   = 1
        plan_fw(2) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_FORWARD, FFTW_ESTIMATE )
        plan_bw(2) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_BACKWARD, FFTW_ESTIMATE )

        n = (/ gmesh(3) /)
        if (pin_x_distribution) then
          howmany = (xzdisp(pz+1)-xzdisp(pz))*(ydispl(py+1)-ydispl(py))
        else
          howmany = (xdispl(py+1)-xdispl(py))*(yzdisp(pz+1)-yzdisp(pz))
        endif
        istride = howmany
        idist   = 1
        inembed = (/ gmesh(3) /)
        plan_fw(3) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_FORWARD, FFTW_ESTIMATE )
        plan_bw(3) = fftw_plan_many_dft( rank, n, howmany,
     &                                    in, inembed, istride, idist,
     &                                    in, inembed, istride, idist,
     &                                    FFTW_BACKWARD, FFTW_ESTIMATE )
        OldMesh = gmesh

      else ! use_fftw = false
#endif
C       Initialise the tables for the FFT if the mesh has changed
        do
          lredimension = .false.
          maxmaxtrigs = maxtrigs
          do i=1,3
            if (OldMesh(i).ne.gmesh(i)) then
              call setgpfa( trigs(:,i), maxtrigs, ntrigs, gmesh(i) )
              if (ntrigs.gt.maxmaxtrigs) then
                lredimension = .true.
                maxmaxtrigs = ntrigs
              else
                OldMesh(i) = gmesh(i)
              endif
            endif
          enddo
          if (.not. lredimension) exit
C
C         Resize FFT array for trig factors and set OldMesh to 0 to force recalculation
C
          maxtrigs = maxmaxtrigs
          call re_alloc( trigs, 1, maxtrigs, 1, 3, 'trigs', 'fft' )
          OldMesh(1:3) = 0
        enddo

#ifdef SIESTA__FFTW
      endif ! use_fftw
#endif
      endif
#ifdef MPI
C
C       Compute the size of the local pencils before
C       the allocation of input arrays
C
        s1 = gMesh(1)*(ydispl(py+1)-ydispl(py))*
     &                (zdispl(pz+1)-zdispl(pz))
        s2 = gMesh(2)*(xdispl(py+1)-xdispl(py))*
     &                (zdispl(pz+1)-zdispl(pz))
      if (pin_x_distribution) then
        s3 = gMesh(3)*(ydispl(py+1)-ydispl(py))*
     &                (xzdisp(pz+1)-xzdisp(pz))
        size1 = s1
        size2 = max(s2,s3)
      else
        s3 = gMesh(3)*(xdispl(py+1)-xdispl(py))*
     &       (yzdisp(pz+1)-yzdisp(pz))
        size1 = max(s1,s3)
        size2 = s2
      endif
#else
        size1 = gMesh(1)*gMesh(2)*gMesh(3)
        size2 = 0
#endif
      k_pin_x_distribution = pin_x_distribution
      end subroutine fft_init

      subroutine getFFTLims( lbox )
C **********************************************************************
C Get the limits of the local pencil in the transformed domain
C Rogeli Grima (June 2022)
C ************ INPUT ***************************************************
C   INTEGER LBOX(2,3) : Local limits
C **********************************************************************
      implicit none
      integer :: lbox(2,3)
#ifdef MPI
      if (k_pin_x_distribution) then
        lbox(1,1) = 1
        lbox(2,1) = gmesh(1)
        lbox(1,2) = ydispl(py)
        lbox(2,2) = ydispl(py+1)-1
        lbox(1,3) = zdispl(pz)
        lbox(2,3) = zdispl(pz+1)-1
      else
        lbox(1,1) = xdispl(py)
        lbox(2,1) = xdispl(py+1)-1
        lbox(1,2) = yzdisp(pz)
        lbox(2,2) = yzdisp(pz+1)-1
        lbox(1,3) = 1
        lbox(2,3) = gmesh(3)
      endif
#else
      lbox(1,1:3) = 1
      lbox(2,1:3) = gmesh(1:3)
#endif
      end subroutine getFFTLims

      subroutine fft( f1, f2, isn )
C **********************************************************************
C Computes the 3-D FFT
C In order to compute the 3D FFT we have to reorder our data in 3 different
C directions: pencil-X, pencil-Y and pencil-Z.
C The user can set two modes for k_pin_x_distribution
C   true: To keep the original distribution (Pencil in X direction)
C         FFT-X MOV(X=>Y) FFT-Y MOV(Y=>X) MOV(X=>Z) FFT-Z MOV(Z=>X)
C   false: If we swap from X to Z (forward) or from Z to X (backward)
C         FFT-X MOV(X=>Y) FFT-Y MOV(Y=>Z) FFT-Z (forward mode)
C         FFT-Z MOV(Z=>Y) FFT-Y MOV(Y=>X) FFT-X (backward mode)
C When we set k_pin_x_distribution to false we can achieve the 3-D FFT with less
C communications, but we have to adapt the calling code in order to work
C with the pencil-Z data distribution
C Rogeli Grima (June 2022)
C ************ INPUT ***************************************************
C   real f1(*) : Input array
C   real f2(*) : Auxiliar array
C   integer isn : FFT direction
C **********************************************************************
      implicit none
      real(grid_p), target :: f1(*), f2(*)
      integer :: isn
      integer :: ierr, n, n1, n2, n3, n2a, n3a, n1b, n3b, n1c, n2c
      integer :: i, j, k, ng, ioffset
      real(dp) :: scale
#ifdef SIESTA__FFTW
      complex(grid_p), pointer :: in(:)
      type(C_PTR), pointer :: plan(:)
#endif
#ifdef DEBUG
      call write_debug( '    PRE fft' )
#endif
      call timer( 'fft', 1 )

      n1 = gmesh(1)
      n2 = gmesh(2)
      n3 = gmesh(3)
#ifdef SIESTA__FFTW
      if (isn<0) then
        plan => plan_fw
      else
        plan => plan_bw
      endif
#endif
#ifdef MPI
      if (k_pin_x_distribution) then
        n2a = ydispl(py+1)-ydispl(py)
        n3a = zdispl(pz+1)-zdispl(pz)
        n1b = xdispl(py+1)-xdispl(py)
        n3b = n3a
        n1c = xzdisp(pz+1)-xzdisp(pz)
        n2c = n2a
C
C       FFT in X direction
C
        n = n1*n2a*n3a
        if (n>0) then
#ifdef SIESTA__FFTW
          if ( use_fftw ) then
            call c_f_pointer( c_loc(f1), in, shape=[n] )
            call fftw_execute_dft( plan(1), in, in )
          else
#endif
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &                 n2a*n3a, -isn )
#ifdef SIESTA__FFTW
          endif ! use_fftw
#endif
        endif
C
C       Transpose data from X to Y direction
C
        call redistributePencil( 1, n1, n2a, n3a, f1,
     &                           2, n1b, n2, n3b, f2,
     &                           py, ProcessorY, xdispl, ydispl,
     &                           COMM_Y )
C
C       FFT in Y direction
C
        n = n1b*n2*n3b
        if (n>0) then
          IOffSet = 1
          do i= 1, n3b
#ifdef SIESTA__FFTW
            if ( use_fftw ) then
              call c_f_pointer( c_loc(f2(IOffSet:IOffSet+n1b*n2-1)),
     &                          in, shape=[n1b*n2] )
              call fftw_execute_dft( plan(2), in, in )
            else
#endif
              call gpfa(f2(IOffSet:2*n),f2(IOffSet+1:2*n), trigs(:,2),
     .              2*n1b, 2, n2, n1b, -isn )
#ifdef SIESTA__FFTW
            endif ! use_fftw
#endif
            IOffSet = Ioffset + 2*n1b*n2
          enddo
        endif
C
C       Transpose data from Y to X and then to Z direction
C
        call redistributePencil( 2, n1b, n2, n3b, f2,
     &                           1, n1, n2a, n3a, f1,
     &                           py, ProcessorY, ydispl, xdispl,
     &                           COMM_Y )
        call redistributePencil( 1, n1, n2a, n3a, f1,
     &                           3, n1c, n2c, n3, f2,
     &                           pz, ProcessorZ, xzdisp, zdispl,
     &                           COMM_Z )
C
C       FFT in Z direction
C
        n = n3*n1c*n2c
        if (n>0) then
#ifdef SIESTA__FFTW
          if ( use_fftw ) then
            call c_f_pointer( c_loc(f2), in, shape=[n] )
            call fftw_execute_dft( plan(3), in, in )
          else
#endif
            call gpfa( f2(1:2*n), f2(2:2*n), trigs(:,3), 2*n1c*n2c, 2,
     &                 n3, n1c*n2c,-isn )
#ifdef SIESTA__FFTW
          endif ! use_fftw
#endif
        endif
C
C       Transpose data back from Z to X
C
        call redistributePencil( 3, n1c, n2c, n3, f2,
     &                           1, n1, n2a, n3a, f1,
     &                           pz, ProcessorZ, zdispl, xzdisp,
     &                           COMM_Z )
C
C       Scale values
C
        if (isn.gt.0) then
          ng = n1*n2*n3
          scale=1.0_dp/dble(ng)
          do i=1,2*n
            f1(i)=f1(i)*scale
          enddo
        endif
      else ! .not. k_pin_x_distribution
        n2a = ydispl(py+1)-ydispl(py)
        n3a = zdispl(pz+1)-zdispl(pz)
        n1b = xdispl(py+1)-xdispl(py)
        n3b = n3a
        n1c = n1b
        n2c = yzdisp(pz+1)-yzdisp(pz)
        if (isn==-1) then
C
C         FFT in X direction
C
          n = n1*n2a*n3a
          if (n>0) then
#ifdef SIESTA__FFTW
            if ( use_fftw ) then
              call c_f_pointer( c_loc(f1), in, shape=[n] )
              call fftw_execute_dft( plan(1), in, in )
            else
#endif
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &                  n2a*n3a, -isn )
#ifdef SIESTA__FFTW
            endif ! use_fftw
#endif
           endif
C
C         Transpose data from X to Y direction
C
          call redistributePencil( 1, n1, n2a, n3a, f1,
     &                             2, n1b, n2, n3b, f2,
     &                             py, ProcessorY, xdispl, ydispl,
     &                             COMM_Y )
C
C         FFT in Y direction
C
          n = n1b*n2*n3b
          if (n>0) then
            IOffSet = 1
            do i= 1, n3b
#ifdef SIESTA__FFTW
              if ( use_fftw ) then
                call c_f_pointer( c_loc(f2(IOffSet:IOffSet+n1b*n2-1)),
     &                            in, shape=[n1b*n2] )
                call fftw_execute_dft( plan(2), in, in )
              else
#endif
              call gpfa(f2(IOffSet:2*n),f2(IOffSet+1:2*n), trigs(:,2),
     .              2*n1b, 2, n2, n1b, -isn )
#ifdef SIESTA__FFTW
              endif ! use_fftw
#endif
              IOffSet = Ioffset + 2*n1b*n2
            enddo
          endif
C
C         Transpose data from Y to Z direction
C
          call redistributePencil( 2, n1b, n2, n3b, f2,
     &                             3, n1c, n2c, n3, f1,
     &                             pz, ProcessorZ, yzdisp, zdispl,
     &                             COMM_Z )
C
C         FFT in Z direction
C
          n = n3*n1c*n2c
          if (n>0) then
#ifdef SIESTA__FFTW
            if ( use_fftw ) then
              call c_f_pointer( c_loc(f1), in, shape=[n] )
              call fftw_execute_dft( plan(3), in, in )
            else
#endif
              call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,3), 2*n1c*n2c, 2,
     &                   n3, n1c*n2c,-isn )
#ifdef SIESTA__FFTW
            endif ! use_fftw
#endif
          endif
        else ! isn==1
C
C         FFT in Z direction
C
          n = n3*n1c*n2c
          if (n>0) then
#ifdef  SIESTA__FFTW
            if ( use_fftw ) then
              call c_f_pointer( c_loc(f1), in, shape=[n] )
              call fftw_execute_dft( plan(3), in, in )
            else
#endif
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,3), 2*n1c*n2c, 2,
     &      n3, n1c*n2c,-isn )
#ifdef SIESTA__FFTW
            endif ! use_fftw
#endif
          endif
C
C         Transpose data from Z to Y direction
C
          call redistributePencil( 3, n1c, n2c, n3, f1,
     &                             2, n1b, n2, n3b, f2,
     &                             pz, ProcessorZ, zdispl, yzdisp,
     &                             COMM_Z )
C
C         FFT in Y direction
C
          n = n1b*n2*n3b
          if (n>0) then
            IOffSet = 1
            do i= 1, n3b
#ifdef  SIESTA__FFTW
            if ( use_fftw ) then
              call c_f_pointer( c_loc(f2(IOffSet:IOffSet+n1b*n2-1)),
     &                  in, shape=[n1b*n2] )
              call fftw_execute_dft( plan(2), in, in )
            else
#endif
              call gpfa(f2(IOffSet:2*n),f2(IOffSet+1:2*n), trigs(:,2),
     .              2*n1b, 2, n2, n1b, -isn )
#ifdef SIESTA__FFTW
            endif ! use_fftw
#endif
              IOffSet = Ioffset + 2*n1b*n2
            enddo
          endif
C
C         Transpose data from Y to X direction
C
          call redistributePencil( 2, n1b, n2, n3b, f2,
     &                             1, n1, n2a, n3a, f1,
     &                             py, ProcessorY, ydispl, xdispl,
     &                             COMM_Y )
C
C         FFT in X direction
C
          n = n1*n2a*n3a
          if (n>0) then
#ifdef SIESTA__FFTW
            if ( use_fftw ) then
              call c_f_pointer( c_loc(f1), in, shape=[n] )
              call fftw_execute_dft( plan(1), in, in )
            else
#endif
            call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &                  n2a*n3a, -isn )
#ifdef SIESTA__FFTW
            endif ! use_fftw
#endif
C
C           SCALE the result
C
            ng = n1*n2*n3
            scale=1.0_dp/dble(ng)
            do i=1,2*n
              f1(i)=f1(i)*scale
            enddo
          endif
        endif
      endif
#else /* NOT MPI */

        n = n1*n2*n3
C
C       FFT in X direction
C
#ifdef SIESTA__FFTW
        if ( use_fftw ) then
          call c_f_pointer( c_loc(f1), in, shape=[n] )
          call fftw_execute_dft( plan(1), in, in )
        else
#endif
        call gpfa( f1(1:2*n), f1(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &             n2*n3, -isn )
#ifdef SIESTA__FFTW
        endif ! use_fftw
#endif
C
C  FFT in Y direction
C
      IOffSet = 1
      do i=0,n3-1
#ifdef  SIESTA__FFTW
        if ( use_fftw ) then
          call c_f_pointer( c_loc(f1(IOffSet:IOffSet+n1*n2-1)),
     &              in, shape=[n1*n2] )
          call fftw_execute_dft( plan(2), in, in )
        else
#endif
        call gpfa(f1(IOffSet:2*n),f1(IOffSet+1:2*n),trigs(:,2),
     .             2*n1,2,n2,n1,-isn)
#ifdef SIESTA__FFTW
        endif ! use_fftw
#endif
        IOffSet = Ioffset + 2*n1*n2
      enddo
C
C  FFT in Z direction
C
#ifdef  SIESTA__FFTW
      if ( use_fftw ) then
        call c_f_pointer( c_loc(f1), in, shape=[n] )
        call fftw_execute_dft( plan(3), in, in )
      else
#endif
      call gpfa(f1(1:2*n),f1(2:2*n),trigs(:,3),2*n1*n2,2,n3,n1*n2,-isn)
#ifdef SIESTA__FFTW
        endif ! use_fftw
#endif
C
C  Scale values
C
      if (isn.gt.0) then
        scale=1.0_dp/dble(n)
        do i=1,2*n
          f1(i)=f1(i)*scale
        enddo
      endif
#endif /* MPI */
      call timer( 'fft', 2 )
#ifdef DEBUG
      call write_debug( '    POS fft' )
#endif
      end subroutine fft

#ifdef MPI
      subroutine redistributePencil( src, s1, s2, s3, fsrc,
     &                               dst, d1, d2, d3, fdst,
     &                               p, np, srcd, dstd, comm )
      implicit none
      integer :: src, dst, s1, s2, s3, d1, d2, d3
      integer :: p, np, srcd(np+1), dstd(np+1)
      real*8 :: fsrc(2,s1,s2,s3), fdst(2,d1,d2,d3)
      integer :: i, j, k, l, ii, jj, kk, pi, pl, pr, nsend, nrecv
      integer :: ierr
      MPI_REQUEST_TYPE :: r_recv
      MPI_REQUEST_TYPE :: r_send
      MPI_STATUS_TYPE :: Status
      MPI_COMM_TYPE   :: COMM
      integer :: bsize_s, bsize_r
      integer :: sbox(2,3), dbox(2,3)
      real(grid_p), pointer  :: sbuf(:,:) => null(), rbuf(:,:) => null()
      sbox(1,:) = 1
      sbox(2,:) = (/ s1, s2, s3 /)
      dbox(1,:) = 1
      dbox(2,:) = (/ d1, d2, d3 /)
C     Compute the buffer sizes.  Use the MAXIMUM pencil extent over all
C     nodes (not node 1's), so the buffers are large enough when the mesh
C     divides unevenly across processors (else later, larger nodes overrun).
      bsize_s = (s1*s2*s3/sbox(2,src))*maxval(srcd(2:np+1)-srcd(1:np))
      bsize_r = (d1*d2*d3/dbox(2,dst))*maxval(dstd(2:np+1)-dstd(1:np))
      call re_alloc( sbuf, 1, 2, 1, bsize_s, 'sbuf', 'fft' )
      call re_alloc( rbuf, 1, 2, 1, bsize_r, 'rbuf', 'fft' )
C
C     Handle transfer of terms which are purely local
C
      sbox(1,src) = srcd(p)
      sbox(2,src) = srcd(p+1)-1
      dbox(1,dst) = dstd(p)
      dbox(2,dst) = dstd(p+1)-1
      ii = dbox(1,3)
      do i= sbox(1,3), sbox(2,3)
        jj = dbox(1,2)
        do j= sbox(1,2), sbox(2,2)
          kk = dbox(1,1)
          do k= sbox(1,1), sbox(2,1)
            fdst(:,kk,jj,ii) = fsrc(:,k,j,i)
            kk = kk+1
          enddo
          jj = jj+1
        enddo
        ii =ii+1
      enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
      do pi=1, NP-1
        ! Let's receive from left node and send to the right node
        pl = mod(p-1+NP-pi,NP)+1
        pr = mod(p-1+pi,NP)+1
        ! Destination boxes
        dbox(1,dst) = dstd(pl)
        dbox(2,dst) = dstd(pl+1)-1
        sbox(1,src) = srcd(pr)
        sbox(2,src) = srcd(pr+1)-1
        ! Number of elements to send/receive
        nsend = 2*PRODUCT( sbox(2,:)-sbox(1,:)+1 )
        nrecv = 2*PRODUCT( dbox(2,:)-dbox(1,:)+1 )
C
C       Collect data to send
C
        l = 1
        do i= sbox(1,3), sbox(2,3)
          do j= sbox(1,2), sbox(2,2)
            do k= sbox(1,1), sbox(2,1)
              sbuf(:,l) = fsrc(:,k,j,i)
              l = l+1
            enddo
          enddo
        enddo
C
C       Exchange data - send to right and receive from left
C
        call MPI_IRecv( rbuf(1,1), nrecv,
     .    MPI_grid_real, pl-1, 1, COMM, r_recv, ierr )
        call MPI_ISend( sbuf(1,1), nsend,
     .    MPI_grid_real, pr-1, 1, COMM, r_send, ierr )
C
C       Wait for receive to complete
C
        call MPI_Wait( r_recv, Status, ierr )
C
C       Place received data into correct array
C
        l = 1
        do i= dbox(1,3), dbox(2,3)
          do j= dbox(1,2), dbox(2,2)
            do k= dbox(1,1), dbox(2,1)
              fdst(:,k,j,i) = rbuf(:,l)
              l = l+1
            enddo
          enddo
        enddo
C
C  Wait for send to complete
C
        call MPI_Wait( r_send, Status, ierr )
      enddo

      call de_alloc( sbuf, 'sbuf', 'fft' )
      call de_alloc( rbuf, 'rbuf', 'fft' )
      end subroutine redistributePencil
#endif /* MPI */

!******************************************************************************
!
!     Exact serial global 3-D FFT (no MPI).  Operates on a packed
!     complex array f(2,*) in GLOBAL natural order
!        index = 1 + i1 + n1*i2 + n1*n2*i3  (0-based global coords).
!     dir = -1 (FORWARD) / +1 (BACKWARD), matching the convention of
!     the distributed fft() (isn).  Reuses the module trigs(:,1:3) which
!     were built by fft_init for the full dims gmesh(1:3).
!
      subroutine fft_serial_global( f, n1, n2, n3, dir )
      implicit none
      integer :: n1, n2, n3, dir
      real(dp), target :: f(2,*)      ! packed complex, global natural order
      integer :: n, i, ioff, isn
      real(dp), pointer :: a(:)
      isn = dir
      n = n1*n2*n3
      call c_f_pointer( c_loc(f(1,1)), a, shape=[2*n] )
C
C     FFT in X direction
C
      call gpfa( a(1:2*n), a(2:2*n), trigs(:,1), 2, 2*n1, n1,
     &           n2*n3, -isn )
C
C     FFT in Y direction
C
      ioff = 1
      do i = 0, n3-1
        call gpfa( a(ioff:2*n), a(ioff+1:2*n), trigs(:,2), 2*n1, 2,
     &             n2, n1, -isn )
        ioff = ioff + 2*n1*n2
      enddo
C
C     FFT in Z direction
C
      call gpfa( a(1:2*n), a(2:2*n), trigs(:,3), 2*n1*n2, 2, n3,
     &           n1*n2, -isn )
C
C     Scale on the backward transform
C
      if (isn > 0) then
        do i = 1, 2*n
          a(i) = a(i) / dble(n)
        enddo
      endif
      end subroutine fft_serial_global

!******************************************************************************

      END MODULE m_fft
