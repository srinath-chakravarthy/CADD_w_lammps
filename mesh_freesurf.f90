!*==mp01.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
! Remove 3 layers of atoms around crack faces, for
!
! Define detection band rings around crack tip
! mesh generator for atomically sharp crack tip
!        1         2         3         4         5         6         7
!23456789012345678901234567890123456789012345678901234567890123456789012
! build a blunt *center, crack
      SUBROUTINE MP01(Id,X,Ix,F,B,Itx,lmp)
      USE MOD_GRAIN
      USE MOD_GLOBAL
      USE MOD_FILE
      USE MOD_BOUNDARY
      USE MOD_CRACK
      USE MOD_MATERIAL
      USE MOD_DD_SLIP
      USE MOD_DISL_PARAMETERS
      use lammps
      IMPLICIT NONE
!*--MP0119
      !     Input variables
      type(c_ptr) :: lmp
      INTEGER Id , Ix , Itx(3,*)
      DOUBLE PRECISION X , F , B
      DIMENSION Id(NDF,1) , X(NXDm,1) , Ix(NEN1,1) , F(NDF,1) , B(NDF,1)

!     common declarations
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
      DOUBLE PRECISION tol

!     local variables
      LOGICAL n1 , n2 , n3 , m1 , m2 , m3 , usedetectionband
      DOUBLE PRECISION , POINTER :: nodeangle(:)
      INTEGER nxregions , nyregions , icell , ndxdy , numnodes , numx , &
     &        numy , ispace , ixrem , iyrem , igrain , coincidentnodes ,&
     &        nodestart , inode , nsort , np1 , numnp0 , countlast , i ,&
     &        j , k , node1 , node2 , node3
      INTEGER nr1 , nr2 , nr3 , nr4
      DOUBLE PRECISION xmax(0:20) , ymax(0:20) , nodesite(3) , dx , dy ,&
     &                 dxdy , xxmin , xxmax , yymin , yymax , &
     &                 yyminorig , delx , dely , xx , yy , mindb , &
     &                 rcutmesh , DIST2 , large
      DATA tol/1.D-6/
      LOGICAL placenode , top , bot , left , right , mirror

      TYPE REGION
         DOUBLE PRECISION XMIN , xmax , YMIN , ymax
      END TYPE REGION
      TYPE (REGION) atomregion , detectionband , innerregion , mirroratomregion , simulationcell
      LOGICAL INSIDEREGION , inside

      INTEGER temp_slip , ii , jj , islp , iii
      DOUBLE PRECISION xslip_start , xslip_end , yslip_start , yslip_end
      DOUBLE PRECISION slip_angle(3) , xxx1 , yyy1 , dxslip , dyslip
      DOUBLE PRECISION xslp1 , xendslp1 , zbqlu01 , rr1 , lnuc , sn , smax
      INTEGER*4 timearray(3)
      DOUBLE PRECISION :: aangle

!! VBS added this to read in tolerances
      COMMON /ELEMCONV/ NUMelold
      INTEGER lcrack , NUMelold , numnpc
!c--JS: Specially for crack asymmetry
      INTEGER dwnumx , dwnumy , dwfactorx , dwfactory
      INTEGER logic
      CHARACTER*80 filename
!c--JS: dwfactor 2 for asym and 1 for sym
      dwfactorx = 1
      dwfactory = 1
!c
      WRITE (*,*)
      WRITE (*,*) 'Generating mesh containing an embedded crack'
      WRITE (*,*)


      X0Crack = 0.01
      Y0Crack = 0.01
!
      xmax(0) = -1.
      ymax(0) = -1.


! Read mesh data
      READ (input_file_unit,*) nxregions , nyregions
      READ (input_file_unit,*) (xmax(i),i=1,nxregions)
      READ (input_file_unit,*) (ymax(i),i=1,nyregions)
      READ (input_file_unit,*) mindb
      READ (input_file_unit,*) rcutmesh
      READ (input_file_unit,*) X0Crack , Y0Crack
      READ (input_file_unit,*) PAD_width

      mirror = .FALSE.

!       if(rcutmesh.lt.0.d0) then
!          rcutmesh=abs(rcutmesh)
! 	write (6,*) 'Mirror not allowed'
! 	stop
!       endif


      DO k = MIN(nxregions,nyregions) + 1 , MAX(nxregions,nyregions)
         IF ( nxregions<nyregions ) THEN
            xmax(k) = xmax(nxregions)
         ELSE
            ymax(k) = ymax(nyregions)
         ENDIF
      ENDDO


!     Extract lattice data
!     Normally, the second atom is dx away from the first due to the
!     way cell is sorted.
!     if there is only 1 atom on the lowest y-plane of cell, then dx is
!     cell width in the x dirn.
      dx = GRAins(1)%CELL(1,2)
      IF ( GRAins(1)%CELL(2,2)>tol ) dx = GRAins(1)%DCELL(1)
!
      DO icell = 1 , GRAins(1)%NCELL
         IF ( GRAins(1)%CELL(2,icell)>tol ) THEN
            dy = GRAins(1)%CELL(2,icell)
            dxdy = GRAins(1)%CELL(1,icell)
            EXIT
         ENDIF
      ENDDO


      IF ( dxdy/=0. ) THEN
         ndxdy = NINT(dx/dxdy)
      ELSE
         ndxdy = 1
      ENDIF

      PRINT * , 'Generating nodes'
! Generate the nodes in each box, the mesh is coarsened with increasing
! distance from the center of the box.
      numx = INT(xmax(nxregions)/dx) + 1
      numy = INT(ymax(nyregions)/dy) + 1
      numnodes = 0
      DO i = -numx , numx

         xx = i*dx
         IF ( ABS(xx)<=xmax(nxregions) ) THEN
!
            DO j = -numy , 0
!
               yy = j*dy
               IF ( ABS(yy)<=ymax(nyregions) ) THEN
!
!! Determine the region in which the node is being placed
                  DO k = MAX(nxregions,nyregions) - 1 , 1 , -1
!c--JS: Symmetric geometry-Commented if ther is crack
!               if ( (abs(xx).gt.XMax(k))
!     $              .or.(abs(yy).gt.YMax(k)) ) go to 21
!c--JS: Asymetry due to the presence of crack
                     IF ( (xx>xmax(k)) .OR. (yy>ymax(k)/dwfactory) )&
     &                    EXIT
                     IF ( (xx<-xmax(k)/dwfactorx) .OR. &
     &                    (yy>ymax(k)/dwfactory) ) EXIT
                     IF ( (xx>xmax(k)) .OR. (yy<-ymax(k)) ) EXIT
                     IF ( (xx<-xmax(k)/dwfactorx) .OR. (yy<-ymax(k)) )&
     &                    EXIT
                  ENDDO
                  ispace = 2**k
                  IF ( ispace==0 ) ispace = 1

!! Decide if a node should be placed in this region
                  ixrem = MOD(ABS(i),ispace)
                  iyrem = MOD(ABS(j),ispace)
                  placenode = (ixrem+iyrem)==0
                  IF ( placenode ) THEN

! Assign the node a position
                     numnodes = numnodes + 1
                     X(1,numnodes) = xx + MOD(j,ndxdy)*dxdy
                     X(2,numnodes) = yy
                  ENDIF
               ENDIF
! Qu modification to remove extra layers of atoms around crack faces
!$$$            if((j.eq.0.and.x(1,numNodes).lt.-0.01*dx).or.
!$$$     &           (j.eq.1.and.x(1,numNodes).lt.-0.01*dx) .or.
!$$$     &      (j.eq.-1.and.x(1,numNodes).lt.-dx) )then
!$$$c     **************************************************************
!$$$c     Add this line for a sharp crack
!$$$c     Comment this line and uncomment to the last 3 lines to
!$$$c     for blunt crack with 3 layers missing
!$$$c     **************************************************************
!$$$c            if(j.eq.0.and.x(1,numNodes).lt.-0.1*dx) then
!$$$              numNodes=numNodes-1
!$$$            endif

! FOR NOW ONLY, assign nodes to the crack faces. Leo mentions
! this is a hack. But this appears legitimate.
!     Add back nodes according to desired crack shape
!$$$  *************************************************
!$$$  Chakravarthy mods for no crack
!$$$  *************************************************
!$$$            if(j.eq.0.and.x(1,numNodes).lt.-0.01*dx) then
!$$$              numNodes = numNodes+2
!$$$              x(1,numNodes-1) = xx+mod(j,ndxdy)*dxdy
!$$$!              x(2,numNodes-1) = yy+2*dy
!$$$              x(2,numNodes-1) = yy+dy
!$$$              x(1,numNodes) = xx+mod(j,ndxdy)*dxdy
!$$$!              x(2,numNodes) = yy-2*dy
!$$$              x(2,numNodes) = yy
!$$$!              print *, i, x(1, numNodes)
!$$$              if(i.eq.-320) then
!$$$                countLast = numNodes-1
!$$$              endif
!$$$            endif
! Qu modification ends
!$$$  *************************************************
!$$$  End Chakravarthy mods for no crack
!$$$  *************************************************

            ENDDO
         ENDIF
      ENDDO


      PRINT * , 'Moving nodes to the nearest atomic sites'
! Move nodes to the nearest atomic sites
      large = 1.E30
      xxmax = -large
      xxmin = large
      yymax = -large
      yymin = large
      DO i = 1 , numnodes
!
         nodesite(1) = X(1,i)
         nodesite(2) = X(2,i)
         IF ( i/=countlast ) CALL NEARESTBSITE(nodesite,1,.FALSE.,X(1,i)&
     &        ,igrain)

!! find xxmax, xxmin, etc.
         xxmax = MAX(X(1,i),xxmax)
         xxmin = MIN(X(1,i),xxmin)
         yymax = MAX(X(2,i),yymax)
         yymin = MIN(X(2,i),yymin)
!                print *, i, x(1,i), x(2,i)
         IF ( X(1,i)==0.0D0 ) THEN
            IF ( X(2,i)==0.0D0 ) temp_slip = i
         ENDIF
      ENDDO

      xslip_start = X(1,temp_slip+1)
      yslip_start = (X(2,temp_slip)+X(2,temp_slip+1))/2.D0

      xxx1 = X(1,temp_slip+1) - X(1,temp_slip)
      yyy1 = X(2,temp_slip+1) - X(2,temp_slip)
      slip_angle(1) = ATAN2(yyy1,xxx1)
      dxslip = ABS(xxx1)*2.D0
      dyslip = ABS(yyy1)
      PRINT * , xxx1 , yyy1 , slip_angle(1)
      PRINT * , 'Slip plane start and end' , temp_slip
      PRINT * , 'x_start' , xslip_start , X(1,temp_slip+1)
      PRINT * , 'y_start' , yslip_start


      xxx1 = X(1,temp_slip-1) - X(1,temp_slip)
      yyy1 = ABS(X(2,temp_slip-1)-X(2,temp_slip))
      slip_angle(2) = ATAN2(yyy1,xxx1)
      PRINT * , xxx1 , yyy1 , slip_angle(2)
      slip_angle(3) = 0.0D0
      PRINT * , 'DX' , dxslip
      PRINT * , 'DY' , dyslip
      slip_angle(3) = 0.0D0


      simulationcell%XMIN = xxmin
      simulationcell%xmax = xxmax
      simulationcell%YMIN = yymin
      simulationcell%ymax = yymax
!     hard coded for 1 grain

      IF ( mirror ) THEN
         NUMnp = numnodes
         DO i = 1 , NUMnp
!
            IF ( X(2,i)<yymin+tol ) THEN
               nodesite(1) = X(1,i)
               nodesite(2) = -136.85
               CALL NEARESTBSITE(nodesite,1,.FALSE.,X(1,i),igrain)
            ENDIF
!
            numnodes = numnodes + 1
            nodesite(1) = X(1,i)
            nodesite(2) = 2*yymin - X(2,i)
            CALL NEARESTBSITE(nodesite,1,.FALSE.,X(1,numnodes),igrain)
!
         ENDDO
!
         yyminorig = yymin
         yymin = 2*yymin - yymax
!
         simulationcell%YMIN = yymin
      ENDIF


! Remove coincident nodes
      PRINT * , 'Removing coincident nodes'
      coincidentnodes = 0
      NUMnp = numnodes
      DO i = NUMnp , 2 , -1
         IF ( i<=numnodes ) THEN
            DO j = 1 , i - 1
               delx = ABS(X(1,i)-X(1,j))
               dely = ABS(X(2,i)-X(2,j))
               IF ( delx+dely<2.*tol ) THEN
                  X(1,j) = X(1,numnodes)
                  X(2,j) = X(2,numnodes)
                  numnodes = numnodes - 1
                  coincidentnodes = coincidentnodes + 1
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF ( coincidentnodes/=0 ) WRITE (6,*) coincidentnodes , &
     &                                 ' coincident nodes removed'

      NUMnp = numnodes
      IF ( NUMnp>MAXnp ) THEN
         WRITE (6,*) '***ERROR: Insufficient storage for nodes'
         WRITE (6,*) '         numnp = ' , NUMnp
         STOP
      ENDIF
      WRITE (6,*) 'Total nodes: numnp = ' , NUMnp


!     Apply boundary conditions and define the boundary for the
!     triangulator

      NCE = 0
      nodestart = 0

!     find all external boundary nodes
      simulationcell%XMIN = simulationcell%XMIN + 10.D0
      simulationcell%xmax = simulationcell%xmax - 10.D0
      simulationcell%YMIN = simulationcell%YMIN + 10.D0
      simulationcell%ymax = simulationcell%ymax - 0.1D0

      PRINT * , 'Detecting the outer cell boundary' , NUMnp
      PRINT * , simulationcell%XMIN , simulationcell%xmax
      PRINT * , simulationcell%YMIN , simulationcell%ymax


      DO i = 1 , NUMnp

         top = (X(2,i)>simulationcell%ymax)
         bot = (X(2,i)<simulationcell%YMIN)
         left = (X(1,i)<simulationcell%XMIN)
         right = (X(1,i)>simulationcell%xmax)

!     store all boundary points, but put crack faces at the beginning of
!     elist.  While you are at it, apply the b.c.s

         IF ( bot .OR. right .OR. left ) THEN
            NCE = NCE + 1
            IF ( NCE>NCEmax ) THEN
               IF ( NCE>NCEmax ) CALL INCREASEELIST(100)
            ENDIF
            ELIst(1,NCE) = i

! apply the b.c's
            IF ( bot .OR. right .OR. left ) THEN
!$$$            if (bot) then
               Id(1,i) = 1
               Id(2,i) = 1
               PRINT * , 'BCs on  node' , i
            ENDIF
!
         ENDIF
      ENDDO

! sort the boundary so that is goes CCW.
      ALLOCATE (nodeangle(NUMnp))
      nodeangle = 0.
      DO i = nodestart + 1 , NCE
         inode = ELIst(1,i)
!     YET ANOTHER HACK
         nodeangle(inode) = DATAN2(X(2,inode)-1.0D0,X(1,inode))
      ENDDO
      nsort = NCE - nodestart
      CALL QSORTR(nsort,ELIst(1,nodestart+1),nodeangle,1,1,1.D0)
      DEALLOCATE (nodeangle)

      NCE = NCE - 1
      nodestart = NCE
!     Now add top boundary
      DO i = 1 , NUMnp
         top = (X(2,i)>simulationcell%ymax)
         IF ( top ) THEN
            NCE = NCE + 1
            IF ( NCE>NCEmax ) THEN
               IF ( NCE>NCEmax ) CALL INCREASEELIST(100)
            ENDIF
            ELIst(1,NCE) = i
            IF ( X(1,i)>-3.0*dx .AND. X(1,i)<3.0*dx ) THEN
               Id(2,i) = 1
               PRINT * , "BC's on atom" , i
            ENDIF
         ENDIF
!$$$         if (top) then
!$$$            id(i,1) = 1
!$$$            id(i,2) = 1
!$$$         end if
      ENDDO

! sort the upper ledge so that is goes CW (from left to right).
      ALLOCATE (nodeangle(NUMnp))
      nodeangle = 0.
      DO i = nodestart + 1 , NCE
         inode = ELIst(1,i)
         IF ( X(2,inode)<tol ) THEN
            nodeangle(inode) = -3.14159265358979323846*X(1,inode)&
     &                         /100.0D0
         ELSE
            nodeangle(inode) = DATAN2(X(2,inode)-tol,X(1,inode))
         ENDIF


      ENDDO
      nsort = NCE - nodestart
      CALL QSORTR(nsort,ELIst(1,nodestart+1),nodeangle,1,1,1.D0)
      DEALLOCATE (nodeangle)
      NCE = NCE - 1

! finish defining the boundary
      DO i = 1 , NCE - 1
         inode = ELIst(1,i)
         ELIst(2,i) = ELIst(1,i+1)
         aangle = -DATAN2(X(2,inode),X(1,inode))
         WRITE (*,'(A7,3I7,1X,3E15.7)') 'Elist ' , i , ELIst(1,i) , &
     &                                  ELIst(2,i) , X(1,inode) , &
     &                                  X(2,inode) , aangle
      ENDDO
      ELIst(2,NCE) = ELIst(1,1)
      PRINT * , 'Elist ' , NCE , ELIst(1,NCE) , ELIst(2,NCE)

      PRINT * , 'No. of outer boundary edges = ' , NCE

      NCB = NCE

!     Triangulate, sets all elements to material 1 for this mesh
      PRINT * , 'Triangulating'
      numnpc = NUMnp
      CALL DELAUNAY(Id,X,Ix,F,B,Itx)
      PRINT * , 'Done'

      WRITE (*,*) 'BEFORE adding overlap'
      WRITE (6,*) 'Number of nodes: numnp = ' , NUMnp
      WRITE (6,*) 'Number of elements: numel = ' , NUMel
      IF ( NUMel>MAXel ) STOP 'too many elements'
      IF ( NUMnp>MAXnp ) STOP 'too many nodes'




!     Find max/min of atomistic region
      innerregion%XMIN = -xmax(1)/dwfactorx
      innerregion%xmax = xmax(1)
      innerregion%YMIN = -ymax(1)
      innerregion%ymax = ymax(1)/dwfactory
      CALL FINDATOMREGIONSIZE(X,dx,dy,tol,innerregion,atomregion)


      DO i = 1 , NUMel
         Ix(NEN1,i) = 1
      ENDDO



!     Find continuum region elements
      IF ( mirror ) THEN
         mirroratomregion%XMIN = atomregion%XMIN
         mirroratomregion%xmax = atomregion%xmax
         mirroratomregion%YMIN = 2*yyminorig - atomregion%YMIN
         mirroratomregion%ymax = 2*yyminorig - atomregion%ymax
      ENDIF
!

!	Check if each node of any continuum element is in the
!	atomistic region
      DO i = 1 , NUMel
!

         node1 = Ix(1,i)
         node2 = Ix(2,i)
         node3 = Ix(3,i)
!
!! Determine if any node is in the atomistic region
         n1 = .NOT.INSIDEREGION(X(1:2,node1),atomregion)
         n2 = .NOT.INSIDEREGION(X(1:2,node2),atomregion)
         n3 = .NOT.INSIDEREGION(X(1:2,node3),atomregion)

!! Determine if any node is in the mirrored atomistic region
         IF ( mirror ) THEN
            m1 = .NOT.INSIDEREGION(X(1:2,node1),mirroratomregion)
            m2 = .NOT.INSIDEREGION(X(1:2,node2),mirroratomregion)
            m3 = .NOT.INSIDEREGION(X(1:2,node3),mirroratomregion)
         ELSE
            m1 = .TRUE.
            m2 = .TRUE.
            m3 = .TRUE.
         ENDIF
!
         IF ( (n1 .AND. n2 .AND. n3) .AND. (m1 .AND. m2 .AND. m3) )&
     &        Ix(NEN1,i) = 0
      ENDDO




!     Add pad atoms in the interface region. Note that elements are
!     not needed in this region, so just add atoms.
!      go to 1234
      PRINT * , 'Adding pad atoms' , PAD_width
      numnodes = NUMnp
!      XMax(2)=XMax(1)+2.d0*rcutmesh
!      YMax(2)=YMax(1)+2.d0*rcutmesh
      xmax(2) = xmax(1) + PAD_width
      ymax(2) = ymax(1) + PAD_width
      numx = INT(xmax(2)/dx) + 1
      numy = INT(ymax(2)/dy) + 1
      dwnumx = INT(xmax(2)/dwfactorx/dx) + 1
      dwnumy = INT(ymax(2)/dwfactory/dy) + 1
      dwnumy = 0
      DO i = -dwnumx , numx
         xx = i*dx
         DO j = -numy , dwnumy
            yy = j*dy
!
            numnodes = numnodes + 1
!
            nodesite(1) = xx + MOD(j,ndxdy)*dxdy
            nodesite(2) = yy
!
            CALL NEARESTBSITE(nodesite,1,.FALSE.,X(1,numnodes),igrain)

!!          Skip this node if it is the atomistic region
            IF ( ABS(j)<0 .AND. i<0 ) THEN
               numnodes = numnodes - 1
            ELSEIF ( INSIDEREGION(X(1:2,numnodes),atomregion) ) THEN
               numnodes = numnodes - 1
            ENDIF
!
            IF ( mirror ) THEN
               nodesite(1) = X(1,numnodes)
               nodesite(2) = 2*yyminorig - X(2,numnodes)
               numnodes = numnodes + 1
               CALL NEARESTBSITE(nodesite,1,.FALSE.,X(1,numnodes),&
     &                           igrain)
            ENDIF
         ENDDO
      ENDDO
      PRINT * , 'Done adding pad atoms'

      np1 = NUMnp + 1
      numnp0 = NUMnp
      NUMnp = numnodes
      IF ( NUMel>MAXel ) STOP 'Too many elements'
      IF ( NUMnp>MAXnp ) STOP 'Too many nodes'



!     Determine nodal character -- continuum, interface, atomistic etc.
      NUM2dnode = NUMnp
      CALL STATUSCALC(X,Ix,.TRUE.)


!     remove the coincident nodes on the interface
      PRINT * , 'Removing coincident nodes near the interface'
      DO i = NUMnp , np1 , -1
         DO j = 1 , numnp0
            IF ( ISRelaxed(j)/=0 ) THEN
               IF ( DIST2(X(1:2,i),X(1:2,j),2)<1.E-6 ) THEN
                  F(1:NDF,i) = F(1:NDF,NUMnp)
                  X(1:NXDm,i) = X(1:NXDm,NUMnp)
                  NUMnp = NUMnp - 1
                  EXIT
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      NUM2dnode = NUMnp
      IF ( NUMperiodz>1 ) CALL INCATOMS(X)

      IF ( NUMnp>MAXnp ) STOP 'Too many nodes'

!      allocate(atomSpecie(numnp))
!      do i=1,numnp
!         atomSpecie(i)=1
!      enddo

!--- Qu added on 08/26/2005
!----- For multimaterils purpose

      numnp0 = NUMnp
!#######################################################################
!------- insert H atoms in the crack tip region
      PRINT * , 'nmaterials' , NMAterials
      IF ( NMAterials>1 ) THEN
         CALL INSERTHATOM(X)
         PRINT * , ' Adding interstitial atom'
      ENDIF

      IF ( NUMnp>MAXnp ) STOP 'Too many nodes'

      ALLOCATE (ATOmspecie(NUMnp))
      DO i = 1 , numnp0
         ATOmspecie(i) = 1
      ENDDO
      DO i = numnp0 + 1 , NUMnp
         PRINT * , '# of H atom' , i
         ATOmspecie(i) = 2
         WRITE (6,*)  , 'int' , X(1,i) , '' , X(2,i) , '' , X(3,i)
      ENDDO


      IHNumber = i
!M modif for removing one atom
      PRINT * , 'done'


      WRITE (*,*) 'Final mesh size'
      WRITE (6,*) 'Total nodes: numnp = ' , NUMnp
      WRITE (6,*) 'Total elements: numel = ' , NUMel
      WRITE (6,*) 'rcutmesh = ' , rcutmesh
      NUMnpp1 = -1


      call write_lammps_data(Id, X, Ix, F, B, Itx, -(xmax(1)+pad_width+10.0), xmax(1)+pad_width+10.0,-(ymax(1) + pad_width+10.0), 0.0)

!     Create a detection band
!     Identify a path, defined by a closed polygon, ccw around the
!     vertices, along which the detection band elements will be placed.
!     If the atomistic continuum interface is not a closed path, define
!     this polygon to extend outside the mesh and surround the atomistic
!     region.


      usedetectionband = .FALSE.
      IF ( .NOT.usedetectionband ) THEN
         NDBpoly = 0
         GOTO 100
      ENDIF

      PRINT * , 'Using a detection band'
!     There is one detection band with 4 points
! Qu modified detection band rings starts
      detectionband%XMIN = atomregion%XMIN + rcutmesh
      detectionband%xmax = atomregion%xmax - rcutmesh
      detectionband%YMIN = atomregion%YMIN + rcutmesh
      detectionband%ymax = atomregion%ymax - rcutmesh
      mindb = INT(mindb/dy)*dy
      nr1 = INT((ABS(detectionband%XMIN)-mindb)/dx) + 1
      nr2 = INT((ABS(detectionband%xmax)-mindb)/dx) + 1
      nr2 = nr1
      nr3 = INT((ABS(detectionband%YMIN)-mindb)/dy) + 1
      nr4 = INT((ABS(detectionband%ymax)-mindb)/dy) + 1

      NDBpoly = MAX(nr1,nr2,nr3,nr4) + 1
      NDBpoly = 2
      PRINT * , 'Detection band' , nr1 , nr2 , nr3 , nr4 , NDBpoly
      CALL ALLOCATE_DB

      DO i = 1 , NDBpoly
         detectionband%XMIN = atomregion%XMIN + rcutmesh*(i-1)*dx
!$$$         detectionBand%xmax = atomRegion%xmax - 2.0*rcutmesh
!$$$         detectionBand%ymin = atomRegion%ymin + 2.0*rcutmesh
!$$$         detectionBand%ymax = atomRegion%ymax - 2.0*rcutmesh

         detectionband%xmax = atomregion%xmax - rcutmesh - (i-1)*dx
         detectionband%YMIN = atomregion%YMIN + rcutmesh + (i-1)*dy
         detectionband%ymax = atomregion%ymax - rcutmesh

!$$$
!$$$         detectionBand%xmax = min( minDb+(i-1)*dx
!$$$     $                            ,atomRegion%xmax - rcutmesh)
!$$$         detectionBand%ymin = max(-minDb-(i-1)*dy
!$$$     $                            ,atomRegion%ymin + rcutmesh)
!$$$         detectionBand%ymax = min( minDb+(i-1)*dy
!$$$     $                            ,atomRegion%ymax - rcutmesh)


         DBPoly(1,1,i) = detectionband%XMIN
         DBPoly(2,1,i) = detectionband%YMIN

         DBPoly(1,2,i) = detectionband%xmax
         DBPoly(2,2,i) = detectionband%YMIN

         DBPoly(1,3,i) = detectionband%xmax
         DBPoly(2,3,i) = detectionband%ymax

         DBPoly(1,4,i) = detectionband%XMIN
         DBPoly(2,4,i) = detectionband%ymax

!$$$     For this particular detection band the first layer is closest
!$$$     to the boundary
!$$$     Any other logic can be used ...
         IF ( i==1 ) DBBoundnear(i) = .TRUE.


!$$$         print *, 'Detection band region', i,
!$$$     $        dbpoly(1,1,i), dbpoly(2,1,i)
!$$$         print *,
!$$$     $        'Detection band region', i,dbpoly(1,2,i), dbpoly(2,2,i
!$$$         print
!$$$     $        *, 'Detection band region', i,dbpoly(1,3,i), dbpoly(2,
!$$$         print
!$$$     $        *, 'Detection band region', i,dbpoly(1,4,i), dbpoly(2,
         PRINT * , 'Detection Band Region' , i
         WRITE (*,'(8f10.3)') DBPoly(1,1,i) , DBPoly(2,1,i) , &
     &                        DBPoly(1,2,i) , DBPoly(2,2,i) , &
     &                        DBPoly(1,3,i) , DBPoly(2,3,i) , &
     &                        DBPoly(1,4,i) , DBPoly(2,4,i)

      ENDDO
! Qu modified detection band rings ends

      WRITE (6,*) 'width of the detection band: rcutmesh = ' , rcutmesh
 100  CALL GEN_SLIP_PLANES(simulationcell%XMIN,simulationcell%xmax,&
     &                     simulationcell%YMIN,simulationcell%ymax,&
     &                     atomregion%XMIN,atomregion%xmax,&
     &                     atomregion%YMIN,atomregion%ymax,slip_angle,&
     &                     GRAins(1)%DCELL(1),xslip_start,yslip_start,&
     &                     dxslip,dyslip,PAD_width)

      X_Move_mesh = 5.D0*dxslip
      MOVemesh = .FALSE.
      MOVed = .FALSE.
      call initialize_lammps(Id,X,Ix,F,B,Itx,lmp)
      
      END SUBROUTINE MP01
!*==insideregion.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015



!
!***********************************************************************
!	Checks to see if a 2D point x is located inside thisRegion.
      LOGICAL FUNCTION INSIDEREGION(X,Thisregion)
      IMPLICIT NONE
!*--INSIDEREGION741

      TYPE REGION
         DOUBLE PRECISION XMIN , XMAX , YMIN , YMAX
      END TYPE REGION
      TYPE (REGION) Thisregion
      DOUBLE PRECISION X(2)

      INSIDEREGION = (X(1)>Thisregion%XMIN .AND. &
     &               X(1)<Thisregion%XMAX .AND. &
     &               X(2)>Thisregion%YMIN .AND. X(2)<Thisregion%YMAX)

      END FUNCTION INSIDEREGION
!*==qsortr.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015

!***********************************************************************
      SUBROUTINE QSORTR(N,List,Xkey,Nxdm,Ind,Sign)
      IMPLICIT NONE
!*--QSORTR759
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION guess , Sign , Xkey
      INTEGER Ind , Nxdm
!*** End of declarations inserted by SPAG
!
      DIMENSION Xkey(Nxdm,1)
      INTEGER List(2,*) , N , ll , lr , lm , nl , nr , ltemp , stktop , &
     &        MAXSTK
!
      PARAMETER (MAXSTK=32)
!
      INTEGER lstack(MAXSTK) , rstack(MAXSTK)
!
      ll = 1
      lr = N
      stktop = 0
      DO
         IF ( ll<lr ) THEN
            nl = ll
            nr = lr
            lm = (ll+lr)/2
            guess = Sign*Xkey(Ind,List(1,lm))
!
!     Find xkeys for exchange
!
 20         DO WHILE ( Sign*Xkey(Ind,List(1,nl))<guess )
               nl = nl + 1
            ENDDO
            DO
               IF ( guess<Sign*Xkey(Ind,List(1,nr)) ) THEN
                  nr = nr - 1
                  CYCLE
               ENDIF
               IF ( nl<(nr-1) ) THEN
                  ltemp = List(1,nl)
                  List(1,nl) = List(1,nr)
                  List(1,nr) = ltemp
                  nl = nl + 1
                  nr = nr - 1
                  GOTO 20
               ENDIF
!
!     Deal with crossing of pointers
!
               IF ( nl<=nr ) THEN
                  IF ( nl<nr ) THEN
                     ltemp = List(1,nl)
                     List(1,nl) = List(1,nr)
                     List(1,nr) = ltemp
                  ENDIF
                  nl = nl + 1
                  nr = nr - 1
               ENDIF
!
!     Select sub-list to be processed next
!
               stktop = stktop + 1
               IF ( nr<lm ) THEN
                  lstack(stktop) = nl
                  rstack(stktop) = lr
                  lr = nr
               ELSE
                  lstack(stktop) = ll
                  rstack(stktop) = nr
                  ll = nl
               ENDIF
               GOTO 100
            ENDDO
         ENDIF
!
!     Process any stacked sub-lists
!
         IF ( stktop/=0 ) THEN
            ll = lstack(stktop)
            lr = rstack(stktop)
            stktop = stktop - 1
            CYCLE
         ENDIF
         EXIT
 100  ENDDO
!
      END SUBROUTINE QSORTR
!*==findatomregionsize.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct



!     Find max/min of atomistic region
      SUBROUTINE FINDATOMREGIONSIZE(Atomcoord,Dx,Dy,Tol,Innerregion,Atomregion)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--FINDATOMREGIONSIZE851

      TYPE REGION
         DOUBLE PRECISION xmin , xmax , ymin , ymax
      END TYPE REGION
      TYPE (REGION) Atomregion , Innerregion
      DOUBLE PRECISION Atomcoord(NDF,*) , Dx , Dy , Tol

!	local variables
      DOUBLE PRECISION xmin , xmax , ymin , ymax , x , y , large
      INTEGER inode
!	functions
      LOGICAL INSIDEREGION

      large = 1.E30

      xmax = -large
      xmin = large

      ymax = -large
      ymin = large

      DO inode = 1 , NUMnp

         x = Atomcoord(1,inode)
         y = Atomcoord(2,inode)
!
         IF ( INSIDEREGION(Atomcoord(1:2,inode),Innerregion) ) THEN
            xmax = MAX(xmax,x)
            ymax = MAX(ymax,y)
            xmin = MIN(xmin,x)
            ymin = MIN(ymin,y)
         ENDIF
!
      ENDDO

      xmax = xmax - Dx + Tol
      ymax = ymax + Tol
      xmin = xmin + Dx - Tol
      ymin = ymin - Tol

      xmax = xmax - Dx
      xmin = xmin + Dx
      ymax = ymax - Dy
      ymin = ymin + Dy

      Atomregion%xmin = xmin
      Atomregion%xmax = xmax
      Atomregion%ymin = ymin
      Atomregion%ymax = ymax

      END SUBROUTINE FINDATOMREGIONSIZE
!*==incatoms.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015

! Qu added begin
!
      SUBROUTINE INCATOMS(X)
      USE MOD_GLOBAL
      USE MOD_GRAIN
      IMPLICIT NONE
!*--INCATOMS911

      DOUBLE PRECISION X
      DIMENSION X(NXDm,1)

      INTEGER i , j , numnp0

      numnp0 = NUMnp
      DO i = 2 , NUMperiodz
         DO j = 1 , numnp0
            IF ( ISRelaxed(j)/=0 ) THEN
               NUMnp = NUMnp + 1
               ISRelaxed(NUMnp) = ISRelaxed(j)
               X(1:2,NUMnp) = X(1:2,j)
               X(3,NUMnp) = X(3,j) + (i-1)*GRAins(1)%DCELL(3)
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE INCATOMS
!*==inserthatom.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
!
! Qu added ends

!	c*********************************************************************
      SUBROUTINE INSERTHATOM(X)
      USE MOD_GLOBAL
      USE MOD_GRAIN
      IMPLICIT NONE
!*--INSERTHATOM940

      DOUBLE PRECISION X
      DIMENSION X(NXDm,1)

!----Local variables
      INTEGER i , j , atomhnum , inhatoms
      CHARACTER*80 input_file

      input_file = 'interstitial.inp'
      OPEN (UNIT=10,FILE=input_file,STATUS='old')
      READ (10,*) inhatoms
!c--  JS
      NUMtoth = inhatoms
      DO i = 1 , inhatoms
         NUMnp = NUMnp + 1
         PRINT * , '# of H atom' , NUMnp
         READ (10,*) atomhnum , (X(j,NUMnp),j=1,3)
      ENDDO
      END SUBROUTINE INSERTHATOM
!*==perturbmesh.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015

!***********************************************************************

!***********************************************************************
      SUBROUTINE PERTURBMESH(X,Nxdm,Numnp,Perturb)
      IMPLICIT NONE
!*--PERTURBMESH967
!*** Start of declarations inserted by SPAG
      INTEGER Numnp , Nxdm
      REAL Perturb , X
!*** End of declarations inserted by SPAG
      END SUBROUTINE PERTURBMESH
!*==getrandomnumber1.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 20


! **********************************************************************
!       logical function inside(x,x1,x2,y1,y2)
!       implicit none
!       double precision x(2),x1,x2,y1,y2
!       inside=(x(1).gt.x1).and.(x(1).lt.x2).and.(x(2).gt.y1).and.(x(2).
!      $     .y2)
!       end
!
      DOUBLE PRECISION FUNCTION GETRANDOMNUMBER1()
      IMPLICIT NONE
!*--GETRANDOMNUMBER1986
      DOUBLE PRECISION random , ZBQLU01
!
      random = -1.D0 + 2.D0*ZBQLU01(0.0D0)
!
      IF ( random>0 ) THEN
         random = 1.D0
      ELSEIF ( random<0 ) THEN
         random = -1.D0
      ELSEIF ( DABS(random)<1.0D-6 ) THEN
         random = 0.D0
      ENDIF
!
      GETRANDOMNUMBER1 = random
      END FUNCTION GETRANDOMNUMBER1


!$$$      subroutine move_mesh(id, x, ix, f, b, dr, db)
!$$$!     Interpolates the displacements and forcs at the current node
!$$$!     x from x+xtip(1), xtip(1) is the currect crack tip position
!$$$!     1) Calculate element in which x(1) + xtip(1) exists
!$$$!     2) Calculate the tri-coord of the x_new
!$$$!     3) Interpolate Displacments and the forces to x_new
!$$$!     4) b(x(1)) = b(x_new(1))
!$$$!     5) dr(x(1)) = dr(x_new(1))
!$$$!     6) db(x(1)) = db(x_new(1))
!$$$!     6) id, ix, f, itx all stay the same
!$$$      use mod_grain
!$$$      use mod_global
!$$$      use mod_file
!$$$      use mod_boundary
!$$$      use mod_crack
!$$$      use mod_material
!$$$      use mod_dd_slip
!$$$      implicit none
!$$$      use "./Disl/mod_disl_parameters"
!$$$      use "./Disl/mod_fem_paramters"
!$$$
!$$$!     Input variables
!$$$      double precision b(ndf,*), x(nxdm,*),f(ndf,*),dr(*),db(*)
!$$$      integer id(ndf,*),ix(nen1,*)
!$$$
!$$$!     common declarations
!$$$      common/debugger/debug
!$$$      logical debug
!$$$      double precision tol
!$$$
!$$$!     Local variables
!$$$      integer i, j, iNode, iel, nel
!$$$      logical MoveMesh
!$$$      logical, allocatable :: examined(:)
!$$$      double precision, allocatable :: x_new(:,:)
!$$$      double precision :: coord(3), el_coord(3,3), disp(3), force(3)
!$$$      double precision fe_locate
!$$$!     x_new = x + xtip
!$$$
!$$$
!$$$
!$$$!!!   Operations are performed only if xtip(1) > x_move_mesh
!$$$!     Need to make this value an input parameter or
!$$$!     a fixed value of the atomistic box size
!$$$      if (xtip(1) > x_move_mesh) then
!$$$         MoveMesh = .true.
!$$$         allocate(x_new(nxdm,numnp))
!$$$         allocate(examined(numnp))
!$$$         examined = .true.
!$$$         x_new(1,:) = x(1,1:numnp) + xtip(1)
!$$$         x_new(2:nxdm,:) = x(2:nxdm,1:numnp)
!$$$      else
!$$$         MoveMesh = .false.
!$$$         print *, 'No moving mesh operations performed'
!$$$         return
!$$$      end if
!$$$
!$$$!     Loop through all the FE elements
!$$$      do iel = 1, numel
!$$$!        Check if the element is a continuum element
!$$$         if (ix(4,iel) .eq. 0) then
!$$$!           Loop through all the nodes of the element
!$$$            do i = 1,nen1-1
!$$$               inode = ix(i,iel)
!$$$               if (isRelaxed(inode) .eq. 1) then
!$$$                  if (examined(inode)) then
!$$$                     ! Locate the element in which x_new exists
!$$$                     nel = fe_locate(x_new(1:nxdm,inode),iel)
!$$$                     ! Store coordinates of the 3 nodes
!$$$                     el_coord = 0.0d0
!$$$                     do j = 1,nen1-1
!$$$                        el_coord(j,1:nxdm) = x(:,ix(j,nel))
!$$$                     end do
!$$$                     examined(inode) = .false.
!$$$                     call fe_tricoord(el_coord(1,:), el_coord(2,:),
!$$$     $                    el_coord(3,:), coord)
!$$$                  end if
!$$$               end if
!$$$            end do
!$$$         end if
!$$$      end do
!$$$
      !$$$      end subroutine move_mesh


