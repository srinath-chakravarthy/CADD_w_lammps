!*==feap.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**    PROGRAM FEAP:
!**      executes the subprograms "feap", "macros" or "stop"
!**
!**   Sub-Programs:
!**      feap: reads global settings, allocates storage and calls pmesh
!**            where material, grain, model and constitutive info is rea
!**      macros: calls the solution macros through the control
!**              subroutine pmacr1
!--
!
      PROGRAM FEAP
      USE MOD_MAIN
      USE MOD_GLOBAL
      USE MOD_BOUNDARY
      USE MOD_DYNAMO
      USE MOD_PARALLEL
      USE MOD_FILE
      use lammps
      IMPLICIT NONE
!*--FEAP20
      CHARACTER*4 wd(7)
      CHARACTER(len=80) :: arg
      CHARACTER(len=:), allocatable :: input_fname
      CHARACTER*80 input
      DATA wd/'feap' , 'macr' , 'stop' , 'xxxx' , 'xxxx' , 'xxxx' , &
     &     'xxxx'/
!
 
      INTEGER neqad , l
      INTEGER i
      DOUBLE PRECISION rcutsq , CUTOFFR2
      integer :: arg_count
      type(c_ptr) :: lmp
!
!      call MPI_INIT(ierr)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

!>    New stuff to read command line arguments
      arg_count = command_argument_count()
      if (arg_count /= 2) then
         call lammps_close(lmp)
         call error_handler("Need exactly 2 command line arguments")
      else
         print *, 'The command line arguments are'
         do i = 1, arg_count
            CALL get_command_argument(i, arg)
            if (i == 1) then
               if (trim(arg) /= "-f" ) then
                  print *, "The format must be -f filename"
                  call lammps_close(lmp)
               end if                  
            end if
            if (i == 2) then
               allocate(character(len=LEN(TRIM(arg))) :: input_fname)
               input_fname = trim(arg)
            end if
            print *, trim(arg)
         end do
      end if
      !> TODO exit on error
      open(unit = input_file_unit, file=input_fname, status='old')
      
      NPRocs = 1
      RANk = 0
      IERr = 0
      call lammps_open_no_mpi('lmp -log log.CADD', lmp)
 
!
!---- variable definition
!
!---- numnp  = number of nodal points
!---- numel  = number of elements
!---- ndm    = dimension of ambient space
!---- nxdm   = dimension of x array (=3)
!---- ndf    = number of degrees of freedom per node
!---- nen    = number of nodes per element
!---- nen1    = number of nodes per element + 1
!---- nad    = added size to element matrices in excess of ndf*nen
!---- nsdm   = dimension of stress array
!---- nquad  = number of quadrature points per element
!---- neq    = number of equations of system
!---- nq     = dimension of element internal variable array
!---- nshp   = dimension of element shape function array
!---- nxsj   = dimension of element jacobian array
!
!
!---- write banner
      WRITE (6,&
     &'('' f i n i t e   e l e m e n t   a n a l y s i s'',       ''   p&
      &r o g r a m'')')
      DO
!
!---- read a card and compare first 4 columns with macro list
         READ (input_file_unit,'(a80)') input
         IF ( input(1:4)==wd(1) ) THEN
!
!---- macro 'feap'
!---- initialize variables
            NEWmesh = .TRUE.
! for neighbor finding:
            NMEth = 2
            NEWlst = 1
            NGTlst = 0
            PERub(1:2) = 1.E30
            PERlb(1:2) = -1.E30
            PERlen(1:2) = 2E30
!
            DD_set = .FALSE.
            NCE = 0
            ALLOCATE (ELIst(2,NCEmax))
!---- read and print control information
            HEAd = input(5:)
            CALL GLOBALSETTINGS
!
            WRITE (6,&
     &'(/'' '',a80//                                             5x,''nu&
     &mber of nodal points (max) ='',i6/                         5x,''nu&
     &mber of elements     (max) ='',i6/                         5x,''di&
     &mension of coordinate space='',i6/                         5x,''de&
     &grees of freedom/node      ='',i6/                         5x,''no&
     &des/element (maximum)      ='',i6/                         5x,''di&
     &mension of stress array    ='',i6/                         5x,''nu&
     &mber of quad pts/element   ='',i6/                         5x,''nu&
     &mber of local nodes/element='',i6)') HEAd , MAXnp , MAXel , NDM , &
     &NDF , NEN , NSDm , NQUad , NAD
 
            WRITE (*,*)
 
            NEQ = MAXnp*NDF
            neqad = (MAXnp+MAXel*NAD)*NDF
 
            ALLOCATE (DR(NDF*MAXnp))
            ALLOCATE (X(NXDm*(MAXnp+MAXel*NAD)))
            ALLOCATE (ISRelaxed(MAXnp))
            ALLOCATE (IDTemp(NDF,MAXnp))
            ALLOCATE (ENErgy(MAXnp))
            ALLOCATE (F(NDF*MAXnp))
            ALLOCATE (DB(neqad))
            ALLOCATE (ID(NDF*MAXnp))
            ALLOCATE (IX(NEN1*MAXel))
            ALLOCATE (ITX(3,MAXel))
            ALLOCATE (B(neqad))
 
            ALLOCATE (B0(NDF,MAXnp+MAXel*NAD))
!C--Jun Song: allocate neighflag for each atom
!     c--DampForce for Marder and langevin Thermostats
            ALLOCATE (UPDateneigh(MAXnp))
            ALLOCATE (DAMpforce(NDF,MAXnp))
 
!
! dynamo variables
!
            ALLOCATE (DIS(3,MAXnp),RDYn(MAXnp),KNBr(MAXnp))
            ALLOCATE (ROLd(3,MAXnp))
            ALLOCATE (NNIndx(0:MAXnp),NNLst(MAXnp,NEIMAX))
!c--JS      nnindx(0)=0
!c--JS      nnindx(1)=0
            NNIndx = 0
 
!     Initializing dynamical allocated variables
            ISRelaxed = 1
            IDTemp = .FALSE.
            ENErgy = 0.D0
            DR = 0.D0
            X = 0.D0
            F = 0.D0
            B = 0.D0
            B0 = 0.D0
            DB = 0.D0
            ID = 0
            IX = 0
            UPDateneigh = 0
 
!     Initialize crack tip position to 0,0
!     If getc macro is not called then crack tip is assumed to be 0,0
!     if getc is called it will initialize the crack tip position
            XTIp_init = 0.0D0
!---- call mesh input subroutine to read and print all mesh data
 
 
            DO i = 1 , NPRocs
               IF ( i==RANk+1 ) THEN
                  CALL PMESH(ID,X,IX,F,B,DR,ITX,lmp)
                  CALL CHKPER
                  NEQ = NUMnp*NDF
                  rcutsq = CUTOFFR2(1)
                  IF ( CUTfact<=1.D0 ) THEN
                     DRAdn = 0.1*DSQRT(rcutsq)
                  ELSE
                     DRAdn = (CUTfact-1.D0)*DSQRT(rcutsq)
                  ENDIF
                  RCTsqn = (SQRT(rcutsq)+DRAdn)**2
                  CALL ECHOSETTINGS
                  CALL LATTICECHECK(X)
                  CALL STATUSCALC(X,IX,.FALSE.)
                  IF ( NQC/=0 ) CALL GETDETECTIONBAND(IX,NEN1,NUMel,X,&
     &                 NUM2dnode,NXDm,ITX,ISRelaxed)
                  CALL PLOTESI(X,IX)
                  CALL PLOTESI_VTK(X,IX)
 
               ENDIF
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            ENDDO
         ELSEIF ( input(1:4)==wd(2) ) THEN
!
!---- macro 'macr'
!---- set up macro program for execution
!---- call macro solution module for establishing solution algorithm
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            CALL PMACR(ID,X,IX,F,B,DR,DB,ITX,lmp)
         ELSEIF ( input(1:4)==wd(3) ) THEN
!
!---- macro 'stop'
            DO l = 7 , 99
               CLOSE (l)
            ENDDO
            EXIT
         ELSEIF ( input(1:4)/=wd(4) ) THEN
            IF ( input(1:4)/=wd(5) ) THEN
               IF ( input(1:4)/=wd(6) ) THEN
                  IF ( input(1:4)==wd(7) ) THEN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
!      call MPI_FINALIZE(ierr)
      END PROGRAM FEAP
 
