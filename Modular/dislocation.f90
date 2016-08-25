!*==mod_dislocation.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
!//   Qu modified on 01/11/2005 to check if the dected dislocation
!//   related to
!//   thermal fluctuation. One subroutine findNearburgers() is added
!//   into.
!
!//   Qu modified on 05/18/2006 to use detection band rings in the
!//   atomistic region.
!//   when dislocation is detected, the corresponding ring is removed,
!//   which makes
!//   the # of the detection band rings decreased.
!
!     this file contains routines that manage the passing of
!     dislocations to and from the atomistic region
!
!     global module variables:
!
!     nburger - number of possible burgers vectors for this crystal
!     structure.  Note that this only uses those whose slip plane
!     normals or in the x-y plane in 2D CADD.
!     nslip - number of detection band elements
!     utilde - current u~ field superimposed of all the existing DD's
!     epsloc(1:3,1:3,i) - strain in detection band (DB) element i
!     amat(1:3,1:2,i) - matrix of shape function derivatives in DB
!     element i
!     enorm - for a given DB element, contains the L2 norm of epsloc
!     -epslib
!     imap(i) - global element number of DB element i
!     iburg(i) - the burgers vector type found in DB element i
!     possible(1:nburger,1:nslip) - .true. if a particular burgers
!     vector is possible in a particular DB element.  This means that
!     the burgers vector is not parallel to the "entry side" of the
!     element.  The entry side is the side of the element facing the
!     atomistics - into which we expect dislocations to enter the DB.
!     newslip - flag that signals when new DD's have been added,
!     necessitating an update of utilde.
!     burg(1:3,i) - library entry of the burger's vector for dislocation
!     i
!     normal(1:3,i) - library entry of the slip plane normal for
!     dislocation i
!     epslib(1:3,1:3,i) - eigenstrain for dislocation i
!     flib(1:3,1:3,i) - deformation gradient for dislocation i
!
      MODULE MOD_DISLOCATION
      IMPLICIT NONE
!*--MOD_DISLOCATION46
      INTEGER nburger , nslip
      DOUBLE PRECISION , ALLOCATABLE :: utilde(:,:) , epsloc(:,:,:) , amat(:,:,:) , enorm(:)
      INTEGER , ALLOCATABLE :: imap(:) , iburg(:)
      LOGICAL , ALLOCATABLE :: possible(:,:)
      LOGICAL newslip
      DOUBLE PRECISION , POINTER :: burg(:,:) , normal(:,:)
      DOUBLE PRECISION , POINTER :: epslib(:,:,:) , flib(:,:,:)
      END MODULE MOD_DISLOCATION
!*==newdislocation.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**********************************************************************
!
!     NewDisocation: provides a means to insert continuum dislocations
!     "by
!     hand" using the macro "newd":
!
!     newd,direct,xpos,ypos,bx,by,bz,itheta_e,theta_s
!
!     where:
!     xpos,ypos - dislocation coordinates
!     bx,by,bz  - burgers vector
!     itheta_e  - angle of branch cut of edge field, in multiples
!     of PI, measured from the direction of b
!     theta_s   - angle of branch cut of the screw field, measured
!     from pos. x axis.
!
!     or:
!
!     newd,file,filename
!
!     to put multiple dislocations in a separate file.
!
      SUBROUTINE NEWDISLOCATION(Input,X,B,Isrelaxed,Numnp,Nxdm,Ndf)
      USE MOD_FILE
      IMPLICIT NONE
!*--NEWDISLOCATION82
!
!     transferred variables
!
      INTEGER Ndf , Numnp , Nxdm
      INTEGER Isrelaxed(Numnp) , X(Nxdm,Numnp) , B(Ndf,Numnp)
      CHARACTER Input*80
!
!     local variables
!
      DOUBLE PRECISION disx(2) , disb(3) , theta_e , theta_s , PI
!---  PI is deliberately less than full 3.14159 to be PI-epsilon
      PARAMETER (PI=3.1415)
 
      INTEGER upper , lower , NEXT , idum , logic , nnewdis , inew , itheta_e
      CHARACTER*80 filename
      CHARACTER*4 key
      key = Input(1:4)
 
      lower = 4
      upper = NEXT(lower,Input)
      IF ( key=='dire' ) THEN
         CALL FREEIN(Input,lower,upper,idum,disx(1),2)
         lower = upper
         upper = NEXT(lower,Input)
         CALL FREEIN(Input,lower,upper,idum,disx(2),2)
         lower = upper
         upper = NEXT(lower,Input)
         CALL FREEIN(Input,lower,upper,idum,disb(1),2)
         lower = upper
         upper = NEXT(lower,Input)
         CALL FREEIN(Input,lower,upper,idum,disb(2),2)
         lower = upper
         upper = NEXT(lower,Input)
         CALL FREEIN(Input,lower,upper,idum,disb(3),2)
         lower = upper
         upper = NEXT(lower,Input)
         CALL FREEIN(Input,lower,upper,itheta_e,idum,1)
         lower = upper
         upper = NEXT(lower,Input)
         CALL FREEIN(Input,lower,upper,idum,theta_s,2)
         nnewdis = 1
      ELSE
         filename = Input(1:(upper-1))
         IF ( .NOT.FILEEXISTS(filename,.FALSE.) ) THEN
            WRITE (*,*) '** WARNING: no dislocation file found'
            RETURN
         ENDIF
         CALL IOFILE(filename,'formatted  ',logic,.TRUE.)
         READ (logic,*) nnewdis
         READ (logic,*) disx , disb , itheta_e , theta_s
      ENDIF
      DO inew = 1 , nnewdis
         IF ( inew>1 ) READ (logic,*) disx , disb , itheta_e , theta_s
         theta_e = PI*MOD(itheta_e,2)
 
         CALL DISL_PASS(disx,disx,disb,theta_e,theta_s,X,B,Isrelaxed,&
     &                  Numnp,.FALSE.,.TRUE.)
         WRITE (*,*) '    New Dislocation Added at:' , disx
         WRITE (*,*) '         with Burgers vector:' , disb
         WRITE (*,*) '                     theta_e:' , theta_e
         WRITE (*,*) '                 and theta_s:' , theta_s
      ENDDO
      IF ( key/='dire' ) CLOSE (logic)
      END SUBROUTINE NEWDISLOCATION
!!$!*==dislcheck.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!!$!**********************************************************************
!!$!
!!$!     dislcheck:
!!$!     main point of contact for checking if dislocations want to pass
!!$!     either way.
!!$!
      LOGICAL FUNCTION DISLCHECK(Checkslip,Lostslip,Addedslip,Movedisl,&
     &                           Ix,X,B,Itx,Isrelaxed,Numnp,Ndf,Nxdm,&
     &                           Numel,Nen1,Newmesh,Plottime, dislpass, npass)
      IMPLICIT NONE
!!$!*--DISLCHECK160
      INTEGER Npass
      LOGICAL Checkslip , Lostslip , Addedslip , Movedisl , Newmesh, dislpass
      INTEGER Numnp , Ndf , Nxdm , Numel , Nen1
      INTEGER Ix(Nen1,Numel) , Itx(3,Numel) , Isrelaxed(Numnp)
      DOUBLE PRECISION X(Nxdm,Numnp) , B(Ndf,Numnp) , Plottime

      DISLCHECK = .FALSE.
!!$!
!!$!     if we are checking for dislocation passings, proceed, otherwise
!!$!     just return
!!$!
      IF ( Checkslip ) THEN
!
!!$!     only check for dislocations leaving the continuum if the
!!$!     dislocations are mobile.  Lostslipcheck looks for dislocations
!!$!     that want to go from continuum to atomistic
         !
         IF ( Movedisl ) THEN
            CALL LOSTSLIPCHECK(Lostslip,Ix,X,B)
            IF ( Lostslip ) THEN
               DISLCHECK = .TRUE.
               RETURN
            ENDIF
         ENDIF
!!$!     check detection band for dislocations that want to go from
!!$!     atomistic to continumm


         IF ( .NOT.Lostslip ) THEN
            CALL SLIPCHECK(X,B,Ix,Itx,Isrelaxed,Numnp,Ndf,Nxdm,Numel,Nen1,Newmesh,Addedslip,Plottime, dislpass)
            Lostslip = .FALSE.
 
            IF ( Addedslip ) THEN
               DISLCHECK = .TRUE.
               RETURN
            ENDIF
         ENDIF
!!$         write(*, '(A,3I7)')'IX after slip_check ==============================', ix(1,1), ix(2,1), ix(3,1)

      ENDIF
      END FUNCTION DISLCHECK
!!$!*==slipcheck.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!!$!*********************************************************************
!!$!
!!$!     slipcheck:  checks detection band for dislocations leaving the
!!$!     atomistic region.
!!$!
      SUBROUTINE SLIPCHECK(X,B,Ix,Itx,Isrelaxed,Numnp,Ndf,Nxdm,Numel,&
     &                     Nen1,Newmesh,Addedslip,Plottime, dislpass)
      USE MOD_DISLOCATION
      USE MOD_FILE
      USE MOD_BOUNDARY
      USE MOD_PARALLEL
      IMPLICIT NONE
!!$!*--SLIPCHECK214
      LOGICAL :: image_flag
      LOGICAL Newmesh , Addedslip
      INTEGER Numnp , Ndf , Numel , Nen1 , Nxdm , i , n1 , n2 , n3
      DOUBLE PRECISION X(Nxdm,Numnp) , B(Ndf,Numnp) , det , Plottime
      INTEGER Ix(Nen1,Numel) , Isrelaxed(Numnp) , Itx(3,*)
      LOGICAL Dislpass
 
      INTEGER j1 , j2 , j3 , i1 , i2 , i3 , node1 , node2 , j , k1 , k2 , iel , k , kmod
      INTEGER itotal , itheta , kp1 , kp2 , ndisl , numnew , ifactor , idb
      INTEGER , VOLATILE :: idbcopy , jcopy
      LOGICAL found(3) , incontinuum , inatoms
      CHARACTER*80 filename
      DOUBLE PRECISION dvec(2,3) , delp(2) , delm(2) , xav(2) , x3(2) , &
     &                 theta_e , theta_s , xd(3) , x0(3) , xi(3) , &
     &                 cross , cvec(2) , bvec(3)
      DOUBLE PRECISION LENGTHTOL , LENGTHTOL2 , PI
!     Qu modification begins
!
      INTEGER nelenear , nelenext , nsidenext , nsideright , nsideleft , ibnext , nsidenearr , nsidenearl , nside , kk
      DOUBLE PRECISION amatnext(3,2) , epslocnext(3,3) , coordside(2) , &
     &                 coordright(2) , coordleft(2) , vecright(2) , &
     &                 vecleft(2) , veclength , dotright , dotleft
      LOGICAL possiblenext(NBUrger)
      LOGICAL , SAVE , ALLOCATABLE :: examined(:,:)
      INTEGER , SAVE :: oidb
      DOUBLE PRECISION :: s_dis
 
      INTEGER , VOLATILE :: islp
 
!!$!     Qu modification ends
!!$!
!!$!---  PI is deliberately less than full 3.14159 to be PI-epsilon
      PARAMETER (LENGTHTOL=1.E-4,LENGTHTOL2=1.E-8,PI=3.1415)
!!$!
!!$!     allocate and initialize data structure
!!$!
      Addedslip = .FALSE.
      IF ( Newmesh ) THEN
         oidb = 0
!!$!     write(*,*) '** First Entry to slipcheck: initializing'
         Newmesh = .FALSE.
         NEWslip = .TRUE.
         IF ( ALLOCATED(EPSloc) ) DEALLOCATE (EPSloc,IMAp,UTIlde,AMAt,&
     &        IBUrg,POSsible,ENOrm,ELIdb,examined)
!!$!--   count the slip detection elements
!!$!     Qu modified detection band rings starts
         NSLip = 0
         DO idb = NDBpoly , 1 , -1
            DO i = 1 , NELidb(idb)
               NSLip = NSLip + 1
            ENDDO
         ENDDO
         PRINT * , 'Number of detection band slip elements =' , NSLip
!     do i=1,numel
!     if(ix(nen1,i).lt.0) then
!     nslip=nslip+1
!     endif
!     enddo
!--   allocate
         IF ( .NOT.ALLOCATED(EPSloc) ) THEN
            ALLOCATE (EPSloc(3,3,NSLip),IMAp(NSLip),UTIlde(3,Numnp),&
     &                AMAt(3,2,NSLip),IBUrg(NSLip),&
     &                POSsible(NBUrger,NSLip),ENOrm(NBUrger),&
     &                ELIdb(NSLip),examined(NBUrger,Numel))
         ELSE
            EPSloc = 0.0D0
            IMAp = 0
            UTIlde = 0.0D0
            AMAt = 0.0D0
            IBUrg = 0
            ENOrm = 0.0D0
            ELIdb = 0
            examined = 0
         ENDIF
 
!--   compute imap
         NSLip = 0
         DO idb = NDBpoly , 1 , -1
            DO i = 1 , NELidb(idb)
               NSLip = NSLip + 1
               IMAp(NSLip) = ELDb(i,idb)
               ELIdb(NSLip) = idb
            ENDDO
         ENDDO
!     do i=1,numel
!     if(ix(nen1,i).lt.0) then
!     nslip=nslip+1
!     imap(nslip)=i
!     endif
!     enddo
!     Qu modified detection band rings ends
!
!--   pre-process the slip elements (compute amat for each)
!
         DO i = 1 , NSLip
            n1 = Ix(1,IMAp(i))
            n2 = Ix(2,IMAp(i))
            n3 = Ix(3,IMAp(i))
            det = (X(1,n1)-X(1,n3))*(X(2,n2)-X(2,n3))&
     &            - (X(1,n2)-X(1,n3))*(X(2,n1)-X(2,n3))
            AMAt(1,1,i) = (X(2,n2)-X(2,n3))/det
            AMAt(2,1,i) = (X(2,n3)-X(2,n1))/det
            AMAt(3,1,i) = (X(2,n1)-X(2,n2))/det
            AMAt(1,2,i) = (X(1,n3)-X(1,n2))/det
            AMAt(2,2,i) = (X(1,n1)-X(1,n3))/det
            AMAt(3,2,i) = (X(1,n2)-X(1,n1))/det
!
!     only allow burgers vectors that are not parallel to entry side of
!     the element.
!
            POSsible(NBUrger,i) = .TRUE.
            k = ABS(Ix(Nen1,IMAp(i)))
            IF ( k/=0 ) THEN
               kp1 = MOD(k,3) + 1
               k = Ix(k,IMAp(i))
               kp1 = Ix(kp1,IMAp(i))
               DO j = 1 , NBUrger - 1
                  cross = BURg(1,j)*(X(2,k)-X(2,kp1)) - BURg(2,j)&
     &                    *(X(1,k)-X(1,kp1))
                  POSsible(j,i) = (ABS(cross)>LENGTHTOL)
!     dw hack
!$$$  if (j.ne.10) possible(j,i)=.false.
!     end hack
               ENDDO
            ELSE
               POSsible(1:NBUrger-1,i) = .FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!     subtract all previously passed dislocations for computation of DB
!     elemental strains
!
 
!     print *, 'In slipcheck and adding displ_contribution'
      IF ( NEWslip ) THEN
         UTIlde = 0
         DO i = 1 , Numnp
            CALL DISL_DISPL(X(1:3,i),UTIlde(1:3,i))
         ENDDO
         NEWslip = .FALSE.
      ENDIF
      B(1:3,1:Numnp) = B(1:3,1:Numnp) - UTIlde(1:3,1:Numnp)
!
!     get current strain in each DB element
!
      DO i = 1 , NSLip
         n1 = Ix(1,IMAp(i))
         n2 = Ix(2,IMAp(i))
         n3 = Ix(3,IMAp(i))
         CALL GETELEMENTSTRAIN(B(1,n1),B(1,n2),B(1,n3),AMAt(1:3,1:2,i),&
     &                         EPSloc(1:3,1:3,i))
      ENDDO
!
!     find nearest possible slip vector for each element
!
      ndisl = 0
      DO i = 1 , NSLip
!     if(elidb(nslip).ne.1) then
         CALL FINDBURGERS(EPSloc(1:3,1:3,i),POSsible(1:NBUrger,i),&
     &                    IBUrg(i),ndisl,ENOrm)
!     endif
      ENDDO
!     **** Finished finding Burger's vector in elements **
      DO i = 1 , NSLip
         DO j = 1 , NBUrger
            examined(j,IMAp(i)) = .FALSE.
         ENDDO
      ENDDO
 
!$$$c     Qu's modification begins
!$$$c     check if the detected slip is related to thermal fluctuation
!$$$      if(ndisl.gt.0) then
!$$$         do i=1,nslip
!$$$
!$$$c     Dw's mod begin
!$$$            if(examined(iburg(i),imap(i))) then
!$$$               ndisl=ndisl-1
!$$$               iburg(i)=nburger
!$$$               goto 21
!$$$            endif
!$$$c$$$            print *, 'XXXX Dislocations in Atomistic = ', ndisl
!$$$
!$$$c     Dw's mod end
!$$$
!$$$            if(iburg(i).ne.nburger) then
!$$$               call checkburgers(epsloc(1:3,1:3,i),possible(1:nburge
!$$$     &              ,iburg(i),x,b,ix,imap(i),enorm)
!$$$c     find the neighbor element of the detection band element i
!$$$               neleNear=itx(abs(ix(nen1,imap(i))),imap(i))
!$$$c--   pre-process the element (compute amat)
!$$$c
!$$$               if(neleNear.ne.0) then
!$$$                  possibleNext(1:nburger)=possible(1:nburger,i)
!$$$                  k=abs(ix(nen1,imap(i)))
!$$$c     vector from centroid to midside of edge k
!$$$                  kp1=mod(k,3)+1
!$$$                  kp2=mod(kp1,3)+1
!$$$                  k=ix(k,imap(i))
!$$$                  kp1=ix(kp1,imap(i))
!$$$                  kp2=ix(kp2,imap(i))
!$$$                  cvec=(x(1:2,k)+x(1:2,kp1)-2*x(1:2,kp2))/6.d0
!$$$                  cross=normal(1,iburg(i))*cvec(2)-normal(2,iburg(i)
!$$$     $                 *cvec(1)
!$$$                  if(cross.lt.0.d0) then
!$$$                     bvec=-burg(1:3,iburg(i))
!$$$                  else
!$$$                     bvec=burg(1:3,iburg(i))
!$$$                  endif
!$$$c
!$$$                  do j=1,3
!$$$                     if(itx(j,neleNear).eq.imap(i))then
!$$$                        nside=j
!$$$                     endif
!$$$                  enddo
!$$$                  nsideRight=mod(nside,3)+1
!$$$                  nsideLeft=mod(nsideRight,3)+1
!$$$
!$$$                  n1=ix(nside,neleNear)
!$$$                  n2=ix(nsideRight,neleNear)
!$$$                  n3=ix(nsideLeft,neleNear)
!$$$c     find the coord of the mid point of nside
!$$$                  coordSide(1:2)=(x(1:2,n1)+x(1:2,n2))/2.d0
!$$$c     find the coord of the mid point of nsideRight
!$$$                  coordRight(1:2)=(x(1:2,n2)+x(1:2,n3))/2.d0
!$$$c     find the coord of the mid point of nsideLeft
!$$$                  coordLeft(1:2)=(x(1:2,n3)+x(1:2,n1))/2.d0
!$$$c
!$$$c     find vector from mid isideRight to mid j
!$$$                  vecRight(1:2)=coordRight(1:2)-coordSide(1:2)
!$$$                  vecLength=dot_product(vecRight(1:2),vecRight(1:2))
!$$$                  vecRight(1:2)=vecRight(1:2)/vecLength
!$$$c     find vector from mid isideLeft to mid j
!$$$                  vecLeft(1:2)=coordLeft(1:2)-coordSide(1:2)
!$$$                  vecLength=dot_product(vecLeft(1:2),vecLeft(1:2))
!$$$                  vecLeft(1:2)=vecLeft(1:2)/vecLength
!$$$
!$$$                  dotRight=dot_product(bvec(1:2),vecRight(1:2))
!$$$                  dotLeft=dot_product(bvec(1:2),vecLeft(1:2))
!$$$                  if(dotRight.gt.dotLeft)then
!$$$                     nsideNext=nsideRight
!$$$                  else
!$$$                     nsideNext=nsideLeft
!$$$                  endif
!$$$                  neleNext=neleNear
!$$$ 10               n1=ix(1,neleNext)
!$$$                  n2=ix(2,neleNext)
!$$$                  n3=ix(3,neleNext)
!$$$
!$$$                  det =(x(1,n1)-x(1,n3))*(x(2,n2)-x(2,n3)) -
!$$$     &                 (x(1,n2)-x(1,n3))*(x(2,n1)-x(2,n3))
!$$$                  amatNext(1,1)=(x(2,n2)-x(2,n3))/det
!$$$                  amatNext(2,1)=(x(2,n3)-x(2,n1))/det
!$$$                  amatNext(3,1)=(x(2,n1)-x(2,n2))/det
!$$$                  amatNext(1,2)=(x(1,n3)-x(1,n2))/det
!$$$                  amatNext(2,2)=(x(1,n1)-x(1,n3))/det
!$$$                  amatNext(3,2)=(x(1,n2)-x(1,n1))/det
!$$$c
!$$$c     only allow burgers vectors that are not parallel to entry side
!$$$c     the element.
!$$$c
!$$$                  possibleNext(nburger)=.true.
!$$$                  kk=abs(ix(nen1,neleNext))
!$$$                  if(kk.ne.0)then
!$$$                     k=ix(nsideNext,neleNext)
!$$$                     kp1=mod(nsideNext,3)+1
!$$$                     kp1=ix(kp1,neleNext)
!$$$                     do j=1,nburger-1
!$$$                        cross=burg(1,j)*(x(2,k)-x(2,kp1))-burg(2,j)
!$$$     &                       *(x(1,k)-x(1,kp1))
!$$$                        possibleNext(j)=(abs(cross).gt.LENGTHTOL)
!$$$                     enddo
!$$$                  else
!$$$                     possibleNext(1:nburger-1)=.false.
!$$$                  endif
!$$$c
!$$$c     get current strain in each DB element
!$$$c
!$$$                  n1=ix(1,neleNext)
!$$$                  n2=ix(2,neleNext)
!$$$                  n3=ix(3,neleNext)
!$$$                  call GetElementStrain(b(1,n1),b(1,n2),b(1,n3),
!$$$     &                 amatNext(1:3,1:2),epslocNext(1:3,1:3))
!$$$     $
!$$$c     find nearest possible slip vector for each element
!$$$                  call findNextburgers(epslocNext,possibleNext,ibNex
!$$$                  if(ibNext.ne.nburger)then
!$$$                     if(itx(nsideNext,neleNext).eq.0)then
!$$$	                goto 20
!$$$                     endif
!$$$                     if(ibNext.eq.iburg(i)) then
!$$$                        examined(iburg(i),neleNext)=.true.
!$$$                     endif
!$$$                     neleNear=neleNext
!$$$                     nsideNearR=nsideRight
!$$$                     nsideNearL=nsideLeft
!$$$                     neleNext=itx(nsideNext,neleNear)
!$$$                     do j=1,3
!$$$                        if(itx(j,neleNext).eq.neleNear)then
!$$$                           nside=j
!$$$                        endif
!$$$                     enddo
!$$$                     nsideRight=mod(nside,3)+1
!$$$                     nsideLeft=mod(nsideRight,3)+1
!$$$                     if(nsideNext.eq.nsideNearR)then
!$$$                        nsideNext=nsideLeft
!$$$                     else
!$$$                        nsideNext=nsideRight
!$$$                     endif
!$$$                     goto 10
!$$$                  else
!$$$                     ndisl=ndisl-1
!$$$                     iburg(i)=nburger
!$$$                  endif
!$$$               else
!$$$                  write(*,*)'!Warning: No past Dislocation path'
!$$$               endif
!$$$ 20            continue
!$$$            endif
!$$$ 21         continue
!$$$c$$$            print *, 'NSLIP = ', nslip, nburger
!$$$	 enddo
!$$$      end if
!$$$c     Qu's modification ends
 
!     Qu modified detection band rings starts
!     if(ndisl.gt.1)then
!     write(*,*)'ERROR----, more than one disl'
!     stop
!     endif
      IF ( ndisl>0 ) THEN
!     find outermost detection ring that is triggered
         j = 0
         DO i = 1 , NSLip
            IF ( j==0 .AND. IBUrg(i)/=NBUrger ) THEN
               idb = ELIdb(i)
               j = 1
               idbcopy = idb
               jcopy = j
            ENDIF
         ENDDO
 
         IF ( idb<=NDBpoly ) THEN
            DO i = 1 , NSLip
               IF ( IBUrg(i)/=NBUrger .AND. idb/=oidb ) THEN
                  IF ( ndisl>1 ) THEN
                     DO j = 1 , NSLip
                        IF ( j/=i ) THEN
!     AddedSlip=.true.
!     ndisl=0
                           IF ( ELIdb(j)==3 .AND. IBUrg(j)/=NBUrger )&
     &                          WRITE (*,*)&
     &                                  '!!!!!more than one disl on ' , &
     &                                 RANk
                        ENDIF
                     ENDDO
                     WRITE (*,*) 'ERROR----, more than one disl'
!     stop
                  ENDIF
!
!     recompute enorms
!
                  CALL FINDBURGERS(EPSloc(1:3,1:3,i),&
     &                             POSsible(1:NBUrger,i),IBUrg(i),ndisl,&
     &                             ENOrm)
!
!     resolve degeneracies - sometimes 2 or more dislocations produce
!     the
!     same strain matrix.  checkburgers chooses the b that
!     produces the best fit with the rotation of the element
!
                  CALL CHECKBURGERS(EPSloc(1:3,1:3,i),&
     &                              POSsible(1:NBUrger,i),IBUrg(i),X,B,&
     &                              Ix,IMAp(i),ENOrm)
 
 
 
!$$$                  write(*,*)
!$$$                  write(*,*) 'slip found in element',imap(i),abs(ix(
!$$$     $                 ,imap(i))), elidb(i), ndbpoly, ' :'m
!$$$
 
!     write(*,*) x(1:2,ix(1,imap(i))),b(1:3,ix(1,imap(i)))
!     write(*,*) x(1:2,ix(2,imap(i))),b(1:3,ix(2,imap(i)))
!     write(*,*) x(1:2,ix(3,imap(i))),b(1:3,ix(3,imap(i)))
!     write(*,*)
!     write(*,*) 'strain matrix:'
!     write(*,'(3e15.6)') (epsloc(j,1:3,i),j=1,3)
 
!$$$                  write(*,*) 'burgers vector:',iburg(i)
!$$$                  write(*,*) burg(1:3,iburg(i))
!$$$                  x0=0.d0
!$$$                  do k=1,3
!$$$                     x0(1:2)=x0(1:2)+x(1:2,ix(k,imap(i)))
!$$$                  enddo
!$$$                  x0(1:2)=x0(1:2)/3.d0
!$$$
!$$$                  write(*,*) 'dislocation at ',x0(1:2)
!$$$                  write(*,*)'time = ', plottime,'ps'
!$$$                  x0(3)=dsqrt(x0(2)*x0(2)+x0(1)*x0(1))
!$$$                  write(*,107) rank,iburg(i),plottime,x0(3),x0(1),x0
!$$$ 107              format (I4,I4,' disdata ',4f10.3)
!$$$                  write(*,*)
                  oidb = idb
 
               ENDIF
            ENDDO
!     write(*,*) 'dislocation at ',x0
!     write(*,*)'time = ', plottime,'ps'
!     write(*,*)
!     nslip=nslip-nelidb(idb)
!            ndisl=0
         ENDIF
      ENDIF
      IF ( ndisl>0 ) PRINT * , 'Outermost Detection Band' , idb
 
!     Qu modified detection band rings ends
!
!     put in the new dislocations
!
      IF ( ndisl>0 ) THEN
         WRITE (*,*) '   ***** New Slip Detected *****   '
         numnew = 0
!$$$         filename='out/detection.plt'
!$$$         call plottrigger(x,b,ix,numel,numnp,nen1,ndf,nxdm,iburg,ima
!$$$     $        ,nslip,utilde,filename,nburger)
!$$$         call plotdisp(x,b,ix,numel,numnp,nen1,ndf,nxdm,IsRelaxed
!$$$     $        ,filename)
         NEWslip = .TRUE.
         Addedslip = .TRUE.
!
!     assume core is at center of each slipped element.
!
!     find x0, xd, theta
!
 
 
!     dw hack
!     goto 505
 
         DO i = 1 , NSLip
            IF ( IBUrg(i)/=NBUrger ) THEN
!
!     recompute enorms
!
               CALL FINDBURGERS(EPSloc(1:3,1:3,i),POSsible(1:NBUrger,i),&
     &                          IBUrg(i),ndisl,ENOrm)
!
!     resolve degeneracies - sometimes 2 or more dislocations produce
!     the
!     same strain matrix.  checkburgers chooses the b that
!     produces the best fit with the rotation of the element
!
               CALL CHECKBURGERS(EPSloc(1:3,1:3,i),POSsible(1:NBUrger,i)&
     &                           ,IBUrg(i),X,B,Ix,IMAp(i),ENOrm)
               WRITE (*,'(A,4I7,A)') 'slip found in element' , IMAp(i) ,  ABS(Ix(Nen1,IMAp(i))) , ELIdb(i) , NDBpoly , ' :'
!$$$               write(*,*) x(1:2,ix(1,imap(i))),b(1:3,ix(1,imap(i)))
!$$$               write(*,*) x(1:2,ix(2,imap(i))),b(1:3,ix(2,imap(i)))
!$$$               write(*,*) x(1:2,ix(3,imap(i))),b(1:3,ix(3,imap(i)))
!$$$               write(*,*)
!$$$               write(*,*) 'strain matrix:'
!$$$               write(*,'(3e15.6)') (epsloc(j,1:3,i),j=1,3)
!$$$               write(*,*)
               WRITE (*,'(A,I5,3F16.6)') 'burgers vector:' , IBUrg(i) , BURg(1:3,IBUrg(i))
 
!
!     process the dislocation:  determine its location, the direction of
!     its
!     branch cut
!     x0: initial location (centroid of DB element)
!     xd: place in the continuum to which disl will be moved.
!     xi: location of the "image" dislocation as far from any continuum
!     region as possible but still on the slip plane.
!
!--   ix(nen1,i) stores the side of the element into which dislocations
!--   may
!     pass
               k = ABS(Ix(Nen1,IMAp(i)))
               IF ( k/=0 ) THEN
!     vector from centroid to midside of edge k
!     Check if dislocation is in the last detection band polygon
                  idb = ELIdb(i)
                  idbcopy = idb
!     Pass dislocations only if the stored Outer polygon is detected
                  kp1 = MOD(k,3) + 1
                  kp2 = MOD(kp1,3) + 1
                  k = Ix(k,IMAp(i))
                  kp1 = Ix(kp1,IMAp(i))
                  kp2 = Ix(kp2,IMAp(i))
                  cvec = (X(1:2,k)+X(1:2,kp1)-2*X(1:2,kp2))/6.D0
                  cross = NORmal(1,IBUrg(i))*cvec(2)&
     &                    - NORmal(2,IBUrg(i))*cvec(1)
                  IF ( cross<0.D0 ) THEN
                     bvec = -BURg(1:3,IBUrg(i))
                  ELSE
                     bvec = BURg(1:3,IBUrg(i))
                  ENDIF
                  IF ( DOT_PRODUCT(bvec(1:2),cvec(1:2))<0.D0 ) THEN
                     itheta = 1
                     theta_s = DATAN2(-bvec(2),-bvec(1))
                  ELSE
                     itheta = 0
                     theta_s = DATAN2(bvec(2),bvec(1))
                  ENDIF
                  x0 = 0.D0
                  DO k = 1 , 3
                     x0(1:2) = x0(1:2) + X(1:2,Ix(k,IMAp(i)))
                  ENDDO
                  x0(1:2) = x0(1:2)/3.D0
!
!     compute location for the continuum core.  FindEntryPoint finds the
!     place where the slip plane intersects the atom/continuum
!     interface.
!     Then we move it 1 burgers vector further along to get it out of
!     the continuum detection region.
!
                  numnew = numnew + 1
!
!     Check of Detection band polygon is a boundary (dbboundnear))
!     pass it ony if it is in this polygon
                  IF ( DBBoundnear(idb) ) THEN
                     CALL FINDENTRYPOINT(bvec,x0,xd)
                     IF ( itheta==0 ) THEN
                        ifactor = -1
                     ELSE
                        ifactor = 1
                     ENDIF
                     xd(1:2)=xd(1:2)+15.0d0*ifactor*bvec(1:2)
                     !! Try to place the dilsocation 

!!$                     CALL FINDSLIPPLANE(bvec,x0,ifactor,islp,theta_s,&
!!$     &                                  s_dis,xd,xi)

                     !
!--   choose location for the "image" dislocation as far from any
!--   detection bands as possible.
!
                     call FindImageLocation(xi,ifactor,x0,bvec(1:3),ix,x,nxdm,numnp,numel,nen1)
                     theta_e = itheta*PI
                     WRITE (*,'(A,3F15.6)', advance = 'no') 'dislocation at ' , x0
                     WRITE (*,'(A,3F15.6)') ' passed to ' , xd
                     WRITE (*,'(A,3F15.6)') 'Image Location' , xi
                     WRITE (*,'(A,3F15.6)') 'with b=' , bvec(1:3)
                     WRITE (*,'(A,F15.6)', advance = 'no') ' and theta_e=' , theta_e
                     WRITE (*,'(A,F15.6)') ' and theta_s=' , theta_s
		     image_flag = .false.
                     CALL DISL_PASS(x0,xd,bvec(1:3),theta_e,theta_s,X,B,Isrelaxed,Numnp,.TRUE.,.TRUE.)
                     !! --- Number of dislocation passed into atomistics
                     itheta=mod(itheta+1,2)
                     theta_e=itheta*PI
                     image_flag = .true. 
                     call disl_pass(xi,xi,-bvec(1:3),theta_e,theta_s,X,B,IsRelaxed,numnp,.true.,.true.)
                     image_flag = .false. 
                     
                     CALL DISL_PRINT(0)
99001                FORMAT (5E15.6)
                     Dislpass = .TRUE.
                  ENDIF
               ELSE
                  WRITE (*,*) &
     &                  'warning: slip in an element on the interior of'
                  WRITE (*,*) '       the detection band'
               ENDIF
            ENDIF
         ENDDO
         IF ( numnew==0 ) STOP 'ERROR: couldn''t resolve dislocations'
         filename = 'out/detection.plt'
         CALL PLOTDISP(X,B,Ix,Numel,Numnp,Nen1,Ndf,Nxdm,Isrelaxed,&
     &                 filename)
      ENDIF
 
!
!     restore b-vector
!
      B(1:3,1:Numnp) = B(1:3,1:Numnp) + UTIlde(1:3,1:Numnp)
 
      END SUBROUTINE SLIPCHECK
!*==getelementstrain.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
!*******************************************************************
!
!     GetElementStrain:  compute strain in an element
!
      SUBROUTINE GETELEMENTSTRAIN(B1,B2,B3,Amat,Eps)
      IMPLICIT NONE
!*--GETELEMENTSTRAIN802
      DOUBLE PRECISION Amat(3,2) , B1(3) , B2(3) , B3(3) , u(3,3) , Eps(3,3) , ua(3,3)
      u(1:3,1) = B1
      u(1:3,2) = B2
      u(1:3,3) = B3
      ua(1:3,1:2) = MATMUL(u,Amat)
      ua(1:3,3) = 0.D0
      Eps = MATMUL(TRANSPOSE(ua),ua)
      Eps = Eps + ua + TRANSPOSE(ua)
      END SUBROUTINE GETELEMENTSTRAIN
!*==checkburgers.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     checkburgers:
!     resolve degeneracies in the dislocation strain library.
!
      SUBROUTINE CHECKBURGERS(Eps,P,Ib,X,B,Ix,Iel,Test)
      USE MOD_GLOBAL
      USE MOD_DISLOCATION
      IMPLICIT NONE
!*--CHECKBURGERS823
      DOUBLE PRECISION Eps(3,3) , X(NXDm,*) , B(NDF,*) , Test(*)
      INTEGER Iel , Ix(NEN1,*) , Ib , ibest , i
      LOGICAL P(*)
!
      DOUBLE PRECISION dx1(3) , dx2(3) , dy1(3) , dy2(3) , norm , del(3), normmin , testmin , dz1(3) , dz2(3)
      dx1 = X(1:3,Ix(2,Iel)) - X(1:3,Ix(1,Iel)) + B(1:3,Ix(2,Iel)) - B(1:3,Ix(1,Iel))
      dx2 = X(1:3,Ix(3,Iel)) - X(1:3,Ix(1,Iel)) + B(1:3,Ix(3,Iel)) - B(1:3,Ix(1,Iel))
      dy1 = X(1:3,Ix(2,Iel)) - X(1:3,Ix(1,Iel))
      dy2 = X(1:3,Ix(3,Iel)) - X(1:3,Ix(1,Iel))
      normmin = 1.E30
      testmin = 1.001*Test(Ib)
      DO i = 1 , NBUrger - 1
         IF ( P(i) .AND. Test(i)<testmin ) THEN
            dz1 = MATMUL(FLIb(1:3,1:3,i),dy1)
            dz2 = MATMUL(FLIb(1:3,1:3,i),dy2)
            del = dx1 - dz1
            norm = DOT_PRODUCT(del,del)
            del = dx2 - dz2
            norm = norm + DOT_PRODUCT(del,del)
            IF ( norm<normmin ) THEN
               normmin = norm
               ibest = i
            ENDIF
         ENDIF
      ENDDO
!//   Qu modification begins
!     write(*,*) '**NOTICE: burgers vector changed in checkburgers'
!//   Qu modification ends
      IF ( ibest/=Ib ) Ib = ibest
      END SUBROUTINE CHECKBURGERS
!*==findburgers.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     findburgers:
!     determine what burgers vector lies in a DB element.
!
      SUBROUTINE FINDBURGERS(Eps,P,Ib,Ndisl,Test)
      USE MOD_DISLOCATION
      IMPLICIT NONE
!*--FINDBURGERS866
      DOUBLE PRECISION Eps(3,3) , Test(*) , xmin , del(3,3)
      INTEGER Ib , Ndisl , i , k , j
      LOGICAL P(NBUrger)
      xmin = 1.E30
      DO i = 1 , NBUrger
         IF ( P(i) ) THEN
            del = Eps - EPSlib(1:3,1:3,i)
            Test(i) = 0
            DO j = 1 , 3
               DO k = j , 3
                  Test(i) = Test(i) + (del(j,k))**2
               ENDDO
            ENDDO
            IF ( Test(i)<xmin ) THEN
               Ib = i
               xmin = Test(i)
            ENDIF
         ELSE
            Test(i) = 1.E30
         ENDIF
      ENDDO
      IF ( Ib==NBUrger ) RETURN
      Ndisl = Ndisl + 1
      END SUBROUTINE FINDBURGERS
!*==findnextburgers.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
!*********************************************************************
!     Qu modification begins
!
!     findNearburgers:
!     determine what burgers vector lies in a DB element.
!
      SUBROUTINE FINDNEXTBURGERS(Eps,P,Ib)
      USE MOD_DISLOCATION
      IMPLICIT NONE
!*--FINDNEXTBURGERS901
      DOUBLE PRECISION Eps(3,3) , test(NBUrger) , xmin , del(3,3)
      INTEGER Ib , ndislnear , i , k , j
      LOGICAL P(NBUrger)
      xmin = 1.E30
      DO i = 1 , NBUrger
         IF ( P(i) ) THEN
            del = Eps - EPSlib(1:3,1:3,i)
            test(i) = 0
            DO j = 1 , 3
               DO k = j , 3
                  test(i) = test(i) + (del(j,k))**2
               ENDDO
            ENDDO
            IF ( test(i)<xmin ) THEN
               Ib = i
               xmin = test(i)
            ENDIF
         ELSE
            test(i) = 1.E30
         ENDIF
      ENDDO
      END SUBROUTINE FINDNEXTBURGERS
!*==rotateburgers.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**********************************************************************
!     Qu modification ends
!
!     rotateburgers:
!     for a given crystal structure and orientation, find all the
!     possible burgers vectors and rotate them to the global coord
!     system.  Also build the strain and defm gradient libraries.
!
      SUBROUTINE ROTATEBURGERS(Q,Rotmat,A0,Struct)
      USE MOD_DISLOCATION
      IMPLICIT NONE
!*--ROTATEBURGERS936
      DOUBLE PRECISION Q(3,3) , Rotmat(2) , q2(3,3) , a6 , a62 , x(3) , &
     &                 A0 , b(3) , b1(3) , m(3) , bx , bz , rt6 , rt3 , &
     &                 d , TOL , rt2 , b2(3)
      PARAMETER (TOL=1.E-6)
      INTEGER i1 , i2 , i3 , i , ib , ip , ipp , is , j , ip2 , ip3 , &
     &        is2
      CHARACTER*3 Struct
      rt6 = DSQRT(6.D0)
      rt3 = DSQRT(3.D0)
      rt2 = DSQRT(2.D0)
      IF ( Struct=='fcc' ) THEN
         NBUrger = 25
      ELSEIF ( Struct=='hex' ) THEN
         NBUrger = 7
      ELSEIF ( Struct=='bcc' ) THEN
         NBUrger = 25
      ELSE
         WRITE (*,*) 'structure type:' , Struct
         STOP 'not recognized in rotateburgers'
      ENDIF
      ALLOCATE (EPSlib(3,3,NBUrger),FLIb(3,3,NBUrger),BURg(4,NBUrger),NORmal(3,NBUrger))
      q2 = 0.D0
      q2(1,1) = Rotmat(1)
      q2(2,2) = Rotmat(1)
      q2(1,2) = -Rotmat(2)
      q2(2,1) = Rotmat(2)
      q2(3,3) = 1.D0
      q2 = MATMUL(q2,Q)
!
      NBUrger = 0
      IF ( Struct=='fcc' ) THEN
         a6 = A0/6.D0
         a62 = 2*a6
         d = A0/rt3
         DO ip = 0 , 3
            m(1:3) = 1.D0/rt3
            IF ( ip/=0 ) m(ip) = -m(ip)
            m = MATMUL(q2,m)
            IF ( ABS(m(3))<TOL ) THEN
               DO ib = 1 , 3
                  CALL GETB(Struct,ip,ib,b1,a6,a62)
                  DO is = -1 , 1 , 2
                     NBUrger = NBUrger + 1
                     NORmal(1:3,NBUrger) = m
                     BURg(1:3,NBUrger) = is*MATMUL(q2,b1)
                     BURg(4,NBUrger) = SQRT(DOT_PRODUCT(b1,b1))
                     CALL GETSTRAIN(BURg(1:3,NBUrger),d,m, EPSlib(1:3,1:3,NBUrger),FLIb(1:3,1:3,NBUrger))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ELSEIF ( Struct=='bcc' ) THEN
         a6 = rt3*A0/2.D0
         d = A0/rt2
         DO ip = 1 , 3
            DO is2 = -1 , 1 , 2
               m(1:3) = 1.D0/rt2
               m(ip) = 0
               ip2 = MOD(ip,3) + 1
               ip3 = MOD(ip2,3) + 1
               m(ip2) = m(ip2)*is2
               b1(1:3) = A0/2.D0
               b2(1:3) = A0/2.D0
               b1(ip2) = -b1(ip2)*(m(ip2)/ABS(m(ip2)))
                                                     !sgn fcn
               b1(ip3) = b1(ip3)*(m(ip3)/ABS(m(ip3)))
                                                     !sgn fcn
               b2(ip2) = -b1(ip2)
               b2(ip3) = -b1(ip3)
               m = MATMUL(q2,m)
               IF ( ABS(m(3))<TOL ) THEN
                  DO is = -1 , 1 , 2
                     NBUrger = NBUrger + 1
                     NORmal(1:3,NBUrger) = m
                     BURg(1:3,NBUrger) = is*MATMUL(q2,b1)
                     BURg(4,NBUrger) = SQRT(DOT_PRODUCT(b1,b1))
                     CALL GETSTRAIN(BURg(1:3,NBUrger),d,m,EPSlib(1:3,1:3,NBUrger),FLIb(1:3,1:3,NBUrger))
                     NBUrger = NBUrger + 1
                     NORmal(1:3,NBUrger) = m
                     BURg(1:3,NBUrger) = is*MATMUL(q2,b2)
                     BURg(4,NBUrger) = SQRT(DOT_PRODUCT(b2,b2))
                     CALL GETSTRAIN(BURg(1:3,NBUrger),d,m,EPSlib(1:3,1:3,NBUrger),FLIb(1:3,1:3,NBUrger))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
!!$--   Qu modification begins
         d = A0/rt6
         DO ip = 1 , 3
            ip2 = MOD(ip,3) + 1
            ip3 = MOD(ip2,3) + 1
            DO ipp = 0 , 3
               m(ip) = 2.D0/rt6
               m(ip2) = 1.D0/rt6
               m(ip3) = 1.D0/rt6
               IF ( ipp/=0 ) m(ipp) = -m(ipp)
               b1(1:3) = A0/2.D0
               DO is = 1 , 3
                  b1(is) = b1(is)*(m(is)/ABS(m(is)))
               ENDDO
               b1(ip) = -b1(ip)
               m = MATMUL(q2,m)
               IF ( ABS(m(3))<TOL ) THEN
                  DO is = -1 , 1 , 2
                     NBUrger = NBUrger + 1
                     NORmal(1:3,NBUrger) = m
                     BURg(1:3,NBUrger) = is*MATMUL(q2,b1)
                     BURg(4,NBUrger) = SQRT(DOT_PRODUCT(b1,b1))
                     CALL GETSTRAIN(BURg(1:3,NBUrger),d,m,EPSlib(1:3,1:3,NBUrger), FLIb(1:3,1:3,NBUrger))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
!!$--   Qu modification ends
      ELSEIF ( Struct=='hex' ) THEN
         a6 = A0/2.D0
         a62 = rt3*a6
         d = a62
         DO ip = 0 , 2
            m(1:3) = 0
            IF ( ip==0 ) THEN
               m(2) = 1
            ELSEIF ( ip==1 ) THEN
               m(1) = -rt3/2.D0
               m(2) = 0.5D0
            ELSE
               m(1) = -rt3/2.D0
               m(2) = -0.5D0
            ENDIF
            m = MATMUL(q2,m)
            IF ( ABS(m(3))<TOL ) THEN
               CALL GETB(Struct,ip,1,b1,a6,a62)
               DO is = -1 , 1 , 2
                  NBUrger = NBUrger + 1
                  NORmal(1:3,NBUrger) = m
                  BURg(1:3,NBUrger) = is*MATMUL(q2,b1)
                  BURg(4,NBUrger) = SQRT(DOT_PRODUCT(b1,b1))
                  WRITE (*,'(i3,3e15.6,2x,4e15.6)') NBUrger , b1 , BURg(1:4,NBUrger)
                  CALL GETSTRAIN(BURg(1:3,NBUrger),d,m,EPSlib(1:3,1:3,NBUrger),FLIb(1:3,1:3,NBUrger))
                  WRITE (*,'(10x,3e15.6)') (EPSlib(j,1:3,NBUrger),j=1,3)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      NBUrger = NBUrger + 1
      EPSlib(1:3,1:3,NBUrger) = 0.D0
      FLIb(1:3,1:3,NBUrger) = 0.D0
      DO i = 1 , 3
         FLIb(i,i,NBUrger) = 1.D0
      ENDDO
      b1(1:3) = 0.D0
      NORmal(1:3,NBUrger) = 0.D0
      BURg(1:4,NBUrger) = 0.D0
      WRITE (*,*) '------burgers vectors-------'
      DO i = 1 , NBUrger
         WRITE (*,'(i3,3e15.6,5x,4e15.6)') i , NORmal(1:3,i) ,  BURg(1:4,i)
         WRITE (*,'(2(10x,3e15.6))') (EPSlib(j,1:3,i),FLIb(j,1:3,i),j=1,3)
      ENDDO
      WRITE (*,*) '-----------------------'
      END SUBROUTINE ROTATEBURGERS
!*==getstrain.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     getstrain:
!     for a given dislocation find the eigenstrain and defm gradient.
!
      SUBROUTINE GETSTRAIN(B,D,M,Eps,F)
      IMPLICIT NONE
!*--GETSTRAIN1119
      DOUBLE PRECISION B(3) , D , M(3) , Eps(3,3) , F(3,3)
      INTEGER i , j
!
!     store dudx in f:
!
      DO i = 1 , 3
         DO j = 1 , 3
            F(i,j) = B(i)*M(j)/D
         ENDDO
      ENDDO
!
!     get strain:
!
      Eps = MATMUL(TRANSPOSE(F),F)
      Eps = (Eps+F+TRANSPOSE(F))
!
!     add delta to f:
!
      DO i = 1 , 3
         F(i,i) = F(i,i) + 1.D0
      ENDDO
      END SUBROUTINE GETSTRAIN
!*==getb.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     getb:
!     just a lookup table of burgers vectors for different crystal
!     structures.
!
      SUBROUTINE GETB(Struct,Ip,Ib,B,A6,A62)
      IMPLICIT NONE
!*--GETB1151
      INTEGER Ip , Ib
      DOUBLE PRECISION B(3) , A6 , A62
      CHARACTER*3 Struct
      IF ( Struct=='fcc' ) THEN
         IF ( Ip==0 ) THEN
            IF ( Ib==1 ) THEN
               B(1) = -A6
               B(2) = A62
               B(3) = -A6
            ELSEIF ( Ib==2 ) THEN
               B(1) = A62
               B(2) = -A6
               B(3) = -A6
            ELSE
               B(1) = -A6
               B(2) = -A6
               B(3) = A62
            ENDIF
         ELSEIF ( Ip==1 ) THEN
            IF ( Ib==1 ) THEN
               B(1) = -A6
               B(2) = A6
               B(3) = -A62
            ELSEIF ( Ib==2 ) THEN
               B(1) = A62
               B(2) = A6
               B(3) = A6
            ELSE
               B(1) = -A6
               B(2) = -A62
               B(3) = A6
            ENDIF
         ELSEIF ( Ip==2 ) THEN
            IF ( Ib==1 ) THEN
               B(1) = A6
               B(2) = -A6
               B(3) = -A62
            ELSEIF ( Ib==2 ) THEN
               B(1) = -A62
               B(2) = -A6
               B(3) = A6
            ELSE
               B(1) = A6
               B(2) = A62
               B(3) = A6
            ENDIF
         ELSEIF ( Ib==1 ) THEN
            B(1) = -A62
            B(2) = A6
            B(3) = -A6
         ELSEIF ( Ib==2 ) THEN
            B(1) = A6
            B(2) = A6
            B(3) = A62
         ELSE
            B(1) = A6
            B(2) = -A62
            B(3) = -A6
         ENDIF
      ELSEIF ( Struct=='bcc' ) THEN
         IF ( Ip==-3 ) THEN
 
         ELSEIF ( Ip==-2 ) THEN
         ELSEIF ( Ip==-1 ) THEN
         ELSEIF ( Ip==1 ) THEN
         ELSEIF ( Ip==2 ) THEN
         ELSEIF ( Ip==3 ) THEN
         ENDIF
      ELSEIF ( Struct=='hex' ) THEN
         IF ( Ip==0 ) THEN
            B(1) = 2.0*A6
            B(2) = 0.D0
            B(3) = 0.D0
         ELSEIF ( Ip==1 ) THEN
            B(1) = A6
            B(2) = A62
            B(3) = 0.D0
         ELSE
            B(1) = -A6
            B(2) = A62
            B(3) = 0.D0
         ENDIF
      ENDIF
      END SUBROUTINE GETB
!*==plottrigger.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     plottrigger:
!     mostly for debugging, plots the DB element that triggered a
!     dislocation pass.
!
      SUBROUTINE PLOTTRIGGER(X,B,Ix,Numel,Numnp,Nen1,Ndf,Nxdm,Iburg,&
     &                       Imap,Nslip,Utilde,Filename,Nburger)
      USE MOD_FILE
      IMPLICIT NONE
!*--PLOTTRIGGER1247
      INTEGER Numel , Numnp , Nen1 , Ndf , Nxdm , Nslip , idisfile
      INTEGER Imap(Nslip) , Ix(Nen1,Numel) , Iburg(Nslip) , Nburger
      DOUBLE PRECISION X(Nxdm,Numnp) , B(Ndf,Numnp) , Utilde(Ndf,Numnp)
      INTEGER i , j , k , im , nplot
      CHARACTER*80 Filename
      CALL IOFILE(Filename,'formatted  ',idisfile,.FALSE.)
      nplot = 0
      DO i = 1 , Nslip
         IF ( Iburg(i)/=Nburger ) nplot = nplot + 1
      ENDDO
      WRITE (idisfile,99001) 'zone, f=fepoint, et=triangle, n=' , &
     &                       3*nplot , ', e=' , nplot
99001 FORMAT (a,i10,a,i10)
      DO i = 1 , Nslip
         DO j = 1 , 3
            IF ( Iburg(i)/=Nburger ) THEN
               DO k = 1 , 3
                  im = Ix(k,Imap(i))
                  WRITE (idisfile,99002) X(1:3,im) + B(1:3,im) , &
     &                   Iburg(i) , Iburg(i) , Iburg(i)
99002             FORMAT (3E15.6,3I4)
               ENDDO
               EXIT
            ENDIF
         ENDDO
      ENDDO
      DO i = 1 , nplot
         WRITE (idisfile,99003) 3*i - 2 , 3*i - 1 , 3*i
99003    FORMAT (3I6)
      ENDDO
      CALL FLUSH(idisfile)
      END SUBROUTINE PLOTTRIGGER
!*==plotdisp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     plotdisp:
!     mostly for debugging, plots the current state of affairs in the
!     atomistic region.
!
      SUBROUTINE PLOTDISP(X,B,Ix,Numel,Numnp,Nen1,Ndf,Nxdm,Isrelaxed,&
     &                    Filename)
      USE MOD_FILE
      IMPLICIT NONE
!*--PLOTDISP1291
      INTEGER Numel , Numnp , Nen1 , Ndf , Nxdm , i
      INTEGER Ix(Nen1,Numel) , numel1 , Isrelaxed(Numnp) , idisfile
      DOUBLE PRECISION X(Nxdm,Numnp) , B(Ndf,Numnp)
      CHARACTER*80 Filename
      CALL IOFILE(Filename,'formatted  ',idisfile,.FALSE.)
      numel1 = 0
      DO i = 1 , Numel
         IF ( Ix(Nen1,i)/=0 ) numel1 = numel1 + 1
      ENDDO
      WRITE (idisfile,99001) 'zone, f=fepoint, et=triangle, n=' , &
     &                       Numnp , ', e=' , numel1
99001 FORMAT (a,i10,a,i10)
      DO i = 1 , Numnp
         IF ( Isrelaxed(i)==0 ) THEN
            WRITE (idisfile,99002) 1000 , 1000 , 0 , 0 , 0 , 0
99002       FORMAT (2I5,4I2)
         ELSE
            WRITE (idisfile,99003) X(1:3,i) + B(1:3,i) , B(1:3,i)
99003       FORMAT (6E15.6)
         ENDIF
      ENDDO
      DO i = 1 , Numel
         IF ( Ix(Nen1,i)/=0 ) WRITE (idisfile,99004) Ix(1:3,i)
99004    FORMAT (3I6)
      ENDDO
      CALL FLUSH(idisfile)
      END SUBROUTINE PLOTDISP
!*==incontinuum.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     InContinuum: check if the point xd is in one of the the continuum
!     region elements
!
      LOGICAL FUNCTION INCONTINUUM(Xd,Ix,X,Nxdm,Numnp,Numel,Nen1,Inatoms)
      IMPLICIT NONE
!*--INCONTINUUM1328
      INTEGER Nxdm , Numnp , Numel , Nen1
      INTEGER Ix(Nen1,Numel) , i
      DOUBLE PRECISION Xd(2) , X(Nxdm,Numnp) , s(3)
      LOGICAL INTRI , ontri , in , Inatoms
      INCONTINUUM = .FALSE.
      Inatoms = .FALSE.
      DO i = 1 , Numel
         in = INTRI(X(1:2,Ix(1,i)),X(1:2,Ix(2,i)),X(1:2,Ix(3,i)),Xd,s, ontri)
         IF ( in .OR. ontri ) THEN
            IF ( Ix(Nen1,i)==0 ) THEN
               INCONTINUUM = .TRUE.
            ELSE
               Inatoms = .TRUE.
            ENDIF
            RETURN
         ENDIF
      ENDDO
      END FUNCTION INCONTINUUM
!*==findimagelocation.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2
!*********************************************************************
!
!     FindImageLocation:  figure out where to put the image dislocation
!     so that it is far away from any continuum regions.
!
      SUBROUTINE FINDIMAGELOCATION(Xi,Ifactor,X0,Bmat,Ix,X,Nxdm,Numnp,Numel,Nen1)
      IMPLICIT NONE
!*--FINDIMAGELOCATION1357
      INTEGER Ifactor , Numel , Nen1 , Numnp , Nxdm
      INTEGER Ix(Nen1,Numel)
      DOUBLE PRECISION X(Nxdm,Numnp)
      DOUBLE PRECISION Bmat(3) , X0(2) , Xi(3) , xnew(2)
      LOGICAL ic , ia , INCONTINUUM , towardscontinuum
      INTEGER istep
!
!     make sure that x0 is in the atomistic region, as it should be.
!
      ic = INCONTINUUM(X0,Ix,X,Nxdm,Numnp,Numel,Nen1,ia)
      IF ( .NOT.ia ) STOP 'ERROR: bug in findimagelocation'
!
!     search for the "other end" of the slip plane, in initial steps of
!     16b, quit refining when step is down to 1b.  March away from the
!     place where the dislocation is passing across the interface until
!     you are in the continuum (or free space, which is even better),
!     then march back with smaller steps until you narrow down the
!     "other end" of the slip plane.
!
      istep = -16
      xnew = X0(1:2)
      DO
         towardscontinuum = istep<0
         IF ( ABS(istep)<1 ) THEN
!
!     at this point, we have either found the point where the slip plane
!     re-enters the continuum or we have found a free surface
!
!--   re-entry: put the image half-way between the two points where slip
!     plane meets continuum.
!--   free:surface: put image just outside the mesh near the free
!--   surface.
!
            IF ( ic ) xnew = 0.5D0*(xnew+X0(1:2))
            Xi(1:2) = xnew
            Xi(3) = 0.D0
            EXIT
         ELSE
            DO
               xnew = xnew + Ifactor*istep*Bmat(1:2)
               ic = INCONTINUUM(xnew,Ix,X,Nxdm,Numnp,Numel,Nen1,ia)
               IF ( towardscontinuum ) THEN
                  IF ( ic .OR. (.NOT.ic .AND. .NOT.ia) ) THEN
                     istep = istep/(-2)
                     EXIT
                  ENDIF
               ELSEIF ( ia ) THEN
                  istep = istep/(-2)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE FINDIMAGELOCATION
!*==getdetectionband.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
!*********************************************************************
!
!     GetDetectionBand:  Takes the user-defined DB polygons (in
!     mod_boundary) and finds the elements that make up the DB.  One
!     should check "esi.tec" to make sure the DB makes sense.  The DB
!     must only hit atomistic elements (ix(4,i)=1) or this routine will
!     complain.
!
      SUBROUTINE GETDETECTIONBAND(Ix,Nen1,Numel,X,Numnp,Nxdm,Itx,&
     &                            Isrelaxed)
      USE MOD_BOUNDARY
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--GETDETECTIONBAND1426
      INTEGER Nen1 , Numel , Numnp , Nxdm
      INTEGER Ix(Nen1,Numel) , Itx(3,Numel) , Isrelaxed(Numnp)
      DOUBLE PRECISION X(Nxdm,Numnp)
 
 
      INTEGER next , NMAX , ifound , iel , j , i , n1 , n2 , jp1 , iring
      PARAMETER (NMAX=1600)
      LOGICAL found
      DOUBLE PRECISION d1 , dmin , vec(2)
      INTEGER , ALLOCATABLE :: nn(:) , neigh(:,:)
      INTEGER ip(2) , idb
      ALLOCATE (nn(Numnp),neigh(DBNMAX,Numnp))
!$$$      allocate(eldb(NMAX,ndbpoly))
!$$$      allocate(dbbound(NMAX,ndbpoly))
!$$$      allocate(nelidb(ndbpoly))
!$$$      dbbound = .false.
!
!     using ix and itx, compute the possible paths between nodes
!
      WRITE (*,*) ' *** computing detection band'
      nn = 0
      DO iel = Numel , 1 , -1
         DO j = 1 , 3
            IF ( Itx(j,iel)<iel ) THEN
               jp1 = MOD(j,3) + 1
               n1 = Ix(j,iel)
               n2 = Ix(jp1,iel)
               nn(n1) = nn(n1) + 1
               IF ( nn(n1)>NMAX ) STOP 'increase NMAX'
               neigh(nn(n1),n1) = n2
               nn(n2) = nn(n2) + 1
               IF ( nn(n2)>NMAX ) STOP 'increase NMAX'
               neigh(nn(n2),n2) = n1
            ENDIF
         ENDDO
      ENDDO
      DO idb = 1 , NDBpoly
         iring = 0
!
!     find nearest node to the first vertex of the dbpoly
!
         dmin = 1.E20
         DO i = 1 , Numnp
            IF ( Isrelaxed(i)==1 ) THEN
               vec = X(1:2,i) - DBPoly(1:2,1,idb)
               d1 = DOT_PRODUCT(vec,vec)
               IF ( d1<dmin ) THEN
                  dmin = d1
                  ip(1) = i
               ENDIF
            ENDIF
         ENDDO
!
!     error:
!
         IF ( nn(ip(1))==0 ) THEN
            WRITE (*,*) ip(1) , X(1:2,ip(1))
            DO i = 1 , NDBvtx(idb)
               WRITE (*,*) DBPoly(1:2,i,idb)
            ENDDO
            DO i = 1 , Numel
               DO j = 1 , 3
                  IF ( Ix(j,i)==ip(1) ) WRITE (*,*) Ix(j,i)
               ENDDO
            ENDDO
            STOP 'ERROR: detection band hits atom pad'
         ENDIF
!
!     choose minimum path to next vertex by moving from one node to the
!     next, always choosing the path that minimizes distance to the
!     next DB polygon vertex.
!
         next = 2
         DO
            vec = X(1:2,ip(1)) - DBPoly(1:2,next,idb)
            dmin = DOT_PRODUCT(vec,vec)
            found = .FALSE.
            DO i = 1 , nn(ip(1))
               vec = X(1:2,neigh(i,ip(1))) - DBPoly(1:2,next,idb)
               d1 = DOT_PRODUCT(vec,vec)
               IF ( d1<dmin ) THEN
                  found = .TRUE.
                  ifound = neigh(i,ip(1))
                  dmin = d1
               ENDIF
            ENDDO
            IF ( found ) THEN
               ip(2) = ifound
!!$	       print *, 'vertex ', ip(1), ip(2), next
               DO iel = 1 , Numel
                  DO j = 1 , 3
                     jp1 = MOD(j,3) + 1
                     IF ( Ix(j,iel)==ip(2) .AND. Ix(jp1,iel)==ip(1) )  THEN
!!$			print *, 'Detection band element = ', iel, -j, next
                        Ix(Nen1,iel) = -j
!     Qu modified detection band ring starts
                        iring = iring + 1
                        IF ( iring>NMAX ) THEN
                           WRITE (*,*) 'error---NMX should be increased'
                           STOP
                        ENDIF
                        ELDb(iring,idb) = iel
                        DBBound(iring,idb) = DBBoundnear(idb)
!     Qu modified detection band ring ends
                        GOTO 10
                     ENDIF
                  ENDDO
               ENDDO
 10            ip(1) = ip(2)
               CYCLE
            ELSE
               next = MOD(next,NDBvtx(idb)) + 1
               IF ( next/=1 ) CYCLE
            ENDIF
            NELidb(idb) = iring
            PRINT * , 'DDDD' , NELidb(idb) , idb , iring
            EXIT
         ENDDO
      ENDDO
!     Find the boundary detection band elements
!     Found by distance to atom boundary such that the minimum
!       DB elements at minimum distance to atom boundary are
!       chosen to be the boundary elements
      DEALLOCATE (nn,neigh)
      END SUBROUTINE GETDETECTIONBAND
!*==lostslipinit.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
      SUBROUTINE LOSTSLIPINIT(Lostslip)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--LOSTSLIPINIT1557
      LOGICAL Lostslip
      Lostslip = .FALSE.
      R_Old(1:3,1:NDIsl) = R_Disl(1:3,1:NDIsl)
      END SUBROUTINE LOSTSLIPINIT
!*==lostslipcheck.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     lostslipcheck:
!     checks for dislocations that want to leave the continuum.  Each
!     dislocation has a range in which it can live (disl_range) that
!     spans the continuum region but for a little gap near the
!     atom/continuum interface and new free surfaces.
!
 
      SUBROUTINE LOSTSLIPCHECK(Lostslip,Ix,X,B)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--LOSTSLIPCHECK1575
      DOUBLE PRECISION X , B , r , r1 , r2
      INTEGER Ix , i , Npass
      LOGICAL Lostslip , pass
      Lostslip = .FALSE.
      Npass = 0
!print *, 'In lostslipcheck'
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i) > 0 ) THEN
            pass = .FALSE.
            r1 = R_Disl(1,i)
            r2 = R_Disl(2,i)
            r = SQRT(r1**2+r2**2)
!$$$            print *, 'In lostslipcheck', i, r1, disl_range(1,i),
!$$$     $           r2, disl_range(2,i)
	    if (abs(r1) <= abs(disl_residence(1,1,i)) .and. abs(r1) <= abs(disl_residence(2,1,i))) then 
		if (abs(r2) <= abs(disl_residence(1,2,i)) .and. abs(r2) <= abs(disl_residence(2,2,i))) then 
		pass = .TRUE. 
		Npass = Npass + 1
		end if
	    end if	

!!$            IF ( ABS(r1)<ABS(DISl_range(1,i)) ) THEN
!!$               IF ( ABS(r2)<ABS(DISl_range(2,i)) ) pass = .TRUE.
!!$               Npass = Npass + 1
!!$            ENDIF
            IF ( pass ) THEN
	       write(*,'(A,I4,4(1X,E15.6),I7,4(1X,E15.6))') 'Entering PasstoAtomistic' , i, &
	       disl_residence(1:2,1,i), disl_residence(1:2,2,i), & 
		ELEM_disl(i), R_old(1:2,i), R_Disl(1:2,i)
!           if(r.lt.disl_range(1,i).or.r.gt.disl_range(2,i)) then
               CALL PASSTOATOMISTIC(R_Disl(1,i),R_Old(1,i),BURgers(1,i),&
     &                              THEta_e(i),THEta_s(i),Ix,X,B,&
     &                              Lostslip,NDIsl_dd(i))
               ELEm_disl(i) = 0
               Lostslip = .TRUE.
!               exit
            ENDIF
         ENDIF
      ENDDO
      R_Old(1:3,1:NDIsl) = R_Disl(1:3,1:NDIsl)
      END SUBROUTINE LOSTSLIPCHECK
!*==passtoatomistic.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
!*********************************************************************
!
!     PassToAtomistic.  As the name suggests, come here when we have
!     found a DD that wants to leave the continuum.
!
      SUBROUTINE PASSTOATOMISTIC(R,Rold,Burgers,Te,Ts,Ix,X,B,Lostslip)
      USE MOD_DISLOCATION
      USE MOD_GLOBAL
      USE MOD_BOUNDARY
      IMPLICIT NONE
!*--PASSTOATOMISTIC1622
      DOUBLE PRECISION R(3) , Rold(3) , Burgers(3) , Te , Ts , vec(2) , &
     &                 xi(3) , xd(3) , X(NXDm,NUMnp) , s(3) , &
     &                 B(NDF,NUMnp)
      INTEGER ifactor , Ix(NEN1,NUMel) , idb , iel , Idis_slip
      INTEGER , STATIC :: nshift = 1
      CHARACTER*80 filename
      LOGICAL POINTINPOLY , Lostslip , notin , INTRI , ontri , ntin
      DO iel = 1 , NUMel
         IF ( INTRI(X(1:2,Ix(1,iel)),X(1:2,Ix(2,iel)),X(1:2,Ix(3,iel)),&
     &        R,s,ontri) ) GOTO 100
      ENDDO
!
!     dislocation is outside of any element, assume it passed to free
!     space.
!
      RETURN
!
!     figure out where to pass the dislocation
!
!
!     use last known location of the dislocation to figure out which way
!     it is moving relative to b
!
 100  vec(1:2) = R(1:2) - Rold(1:2)
      IF ( DOT_PRODUCT(vec(1:2),Burgers(1:2))<0.D0 ) THEN
         ifactor = -1
      ELSE
         ifactor = 1
      ENDIF
      xd(1:2) = R(1:2) + ifactor*Burgers(1:2)
!
!     if disl is still in the continuum, march along by 'b' until it is
!     either in free space or the atomistics
!
 200  IF ( Ix(NEN1,iel)/=0 ) THEN
         DO
!
!     disl is definitely moving to atomistic region.  march along until
!     you are inside a detection band polygon.
!
            PRINT * , 'dislocation in atomistic region' , xd
            notin = .TRUE.
            DO idb = 1 , NDBpoly
!     Make sure that the dislocation is just outside the outermost
!      Detection band polygon
               IF ( DBBoundnear(idb) ) THEN
                  notin = notin .AND. &
     &                    (.NOT.POINTINPOLY(xd,NDBvtx(idb),DBPoly(1,1,&
     &                    idb)))
                  PRINT * , 'Checking dislocation in ' , idb
               ENDIF
            ENDDO
            IF ( notin ) THEN
               xd(1:2) = xd(1:2) + ifactor*Burgers(1:2)
               CYCLE
            ENDIF
            xd(1:2) = xd(1:2) + nshift*ifactor*Burgers(1:2)
            PRINT * , 'Disl. Position in atomistics' , ifactor , xd
            ifactor = -ifactor
            CALL FINDIMAGELOCATION(xi,ifactor,xd,Burgers,Ix,X,NXDm,NUMnp,NUMel,NEN1)
!
!     puts the discrete dislocation as far from the detection bands as
!     possible
!
            R = xi
!
!     adjusts the atomistic displacements so that the new atomistic core
!     is just inside the detection band
!
 
            PRINT * , 'Dislocation initially at' , Rold(1:2)
            CALL DISL_PASS(Rold,xd,Burgers,Te,Ts,X,B,ISRelaxed,NUMnp,.TRUE.,.FALSE.)
            NEWslip = .TRUE.
            Lostslip = .TRUE.
            WRITE (*,*) 'Dislocation passed to Atomistics'
            WRITE (*,*)  , 'With burgers vector' , Burgers(1:2)
            WRITE (*,*) 'initially at:' , Rold(1:2)
            WRITE (*,*) 'moved to:    ' , R(1:2)
            WRITE (*,*) 'atomistic core at:    ' , xd(1:2)
            filename = 'out/final.plt'
            CALL PLOTDISP(X,B,Ix,NUMel,NUMnp,NEN1,NDF,NXDm,ISRelaxed,&
     &                    filename)
            GOTO 99999
         ENDDO
      ENDIF
!
!     disl is still in the continuum - nudge it along by b:
!
      R(1:2) = R(1:2) + ifactor*Burgers(1:2)
      DO iel = 1 , NUMel
         IF ( INTRI(X(1:2,Ix(1,iel)),X(1:2,Ix(2,iel)),X(1:2,Ix(3,iel)),&
     &        R,s,ontri) ) GOTO 200
      ENDDO
!
!     disl is in free space, return without passing to atomistics0
!
      RETURN
99999 END SUBROUTINE PASSTOATOMISTIC
!*==disl_restart.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!*********************************************************************
!
!     disl_restart:
!     write a restart file during minimization, just in case something
!     goes wrong.
!
      SUBROUTINE DISL_RESTART(Id,X,Ix,Itx,F,B)
      USE MOD_GLOBAL
      USE MOD_FILE
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--DISL_RESTART1735
      INTEGER Id(*) , Ix(*) , Itx(*) , logic
      CHARACTER key*4 , filename*80
      DOUBLE PRECISION X(*) , F(*) , B(*)
      key = 'writ'
      filename = 'disloc.res'
      CALL IOFILE(filename,'unformatted',logic,.FALSE.)
      CALL REST(Id,X,Ix,Itx,F,B,logic,key)
      CLOSE (logic)
      WRITE (*,*) '** updating disloc.res ' , NDIsl
      END SUBROUTINE DISL_RESTART
