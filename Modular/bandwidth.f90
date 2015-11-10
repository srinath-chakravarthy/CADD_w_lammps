!*==bandnl.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
      SUBROUTINE BANDNL(Id,X,Ix,F,B)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--BANDNL6
!
!---- optimize  bandwidth/profile
!
      INTEGER NGRa , IDPth , IDEg
      COMMON /GRA   / NGRa , IDPth , IDEg
      DOUBLE PRECISION X(NXDm,*) , F(NDF,*) , B(NDF,*)
      INTEGER Ix(NEN1,*) , Id(NDF,*)
      INTEGER , POINTER :: ndstk(:,:) , ndeg(:) , iold(:) , lvl(:) , &
     &                     renum(:) , ccstor(:) , lvls1(:) , lvls2(:)
      LOGICAL deglimit
      INTEGER maxdeg , i , ibw2 , ipf2
 
      PRINT * , '** Bandwidth Optimization **'
 
      maxdeg = 600
      ALLOCATE (ndstk(maxdeg,NUMnp),ndeg(NUMnp),iold(NUMnp),&
     &          renum(NUMnp+1),lvl(NUMnp),lvls1(NUMnp),lvls2(NUMnp),&
     &          ccstor(NUMnp))
      NGRa = NUMnp
!
      renum = 0
      lvl = 0
      lvls1 = 0
      lvls2 = 0
      ccstor = 0
      DO i = 1 , NUMnp
         iold(i) = i
      ENDDO
      DO
!
         CALL SETCONGR(Ix,ndstk,ndeg,maxdeg,deglimit,B,X)
         IF ( deglimit ) THEN
            maxdeg = maxdeg*2
            DEALLOCATE (ndstk)
            ALLOCATE (ndstk(maxdeg,NUMnp))
            WRITE (*,*) '**WARNING: maxdeg increased to ' , maxdeg
            WRITE (*,*) '           in bandwidth optimization'
            CYCLE
         ENDIF
         CALL REDUCE(ndstk,maxdeg,iold,renum,ndeg,lvl,lvls1,lvls2,&
     &               ccstor,ibw2,ipf2)
         CALL SWAPALLGR(renum,Id,X,Ix,F,B)
         DEALLOCATE (ndstk,ndeg,iold,renum,lvl,lvls1,lvls2,ccstor)
         EXIT
      ENDDO
      END SUBROUTINE BANDNL
!*==setcongr.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!---VBS------------------------------
!-------------------
      SUBROUTINE SETCONGR(Ix,Ndstk,Ndeg,Maxdeg,Deglimit,B,X)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--SETCONGR61
!
      INTEGER Maxdeg
      INTEGER Ix(NEN1,*) , Ndstk(Maxdeg,*) , Ndeg(*)
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*)
      LOGICAL Deglimit , needlist
!
      INTEGER NGRa , IDPth , IDEg
      COMMON /GRA   / NGRa , IDPth , IDEg
      INTEGER nloclist(40)
      INTEGER node , j , irep , jpoint , nnlist , ilist , jlist , ii , &
     &        l , iel , natms , igrain , i
      Deglimit = .FALSE.
 
!
      DO node = 1 , NUMnp
         Ndeg(node) = 0
         DO j = 1 , Maxdeg
            Ndstk(j,node) = 0
         ENDDO
      ENDDO
!
      !Now Handle Local Elements
 
      DO iel = 1 , NUMel
         nnlist = 3
         DO ilist = 1 , nnlist
            node = Ix(ilist,iel)
            DO jlist = 1 , nnlist
               ii = Ix(jlist,iel)
               IF ( ii/=node ) THEN
                  IF ( Ndeg(node)>0 ) THEN
                     DO l = 1 , Ndeg(node)
                        IF ( Ndstk(l,node)==ii ) GOTO 20
                     ENDDO
                  ENDIF
                  IF ( Ndeg(node)<Maxdeg ) THEN
                     Ndeg(node) = Ndeg(node) + 1
                     Ndstk(Ndeg(node),node) = ii
                  ELSE
                     Deglimit = .TRUE.
                     RETURN
                  ENDIF
               ENDIF
 20         ENDDO
         ENDDO
      ENDDO
 
      IDEg = 0
      DO i = 1 , NUMnp
         IF ( Ndeg(i)>IDEg ) IDEg = Ndeg(i)
         Ndeg(i) = 0
      ENDDO
      PRINT * , 'Maximum node degree = ' , IDEg
      END SUBROUTINE SETCONGR
!*==swapallgr.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!-------------------VBS---------------------
!-------------------------------------------
      SUBROUTINE SWAPALLGR(Num,Id,X,Ix,F,B)
      USE MOD_GRAIN
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      USE MOD_BOUNDARY
      IMPLICIT NONE
!*--SWAPALLGR127
!
      INTEGER Num(*) , Id(NDF,*) , Ix(NEN1,*)
      DOUBLE PRECISION X(NXDm,*) , F(NDF,*) , B(NDF,*)
!
      INTEGER i , j , irep , jpoint , nnlist , jlist , ii
      INTEGER , DIMENSION(:) , POINTER :: itemp
      DOUBLE PRECISION , DIMENSION(:) , POINTER :: temp
 
!
!  force a new neighbor list
!
      NEWlst = 1
!
!     allocate
!
      ALLOCATE (temp(NUMnp),itemp(NUMnp))
!
      DO j = 1 , NDF
         DO i = 1 , NUMnp
            itemp(Num(i)) = Id(j,i)
         ENDDO
         DO i = 1 , NUMnp
            Id(j,i) = itemp(i)
         ENDDO
      ENDDO
!
      DO j = 1 , NXDm
         DO i = 1 , NUMnp
            temp(Num(i)) = X(j,i)
         ENDDO
         DO i = 1 , NUMnp
            X(j,i) = temp(i)
         ENDDO
      ENDDO
!
      DO j = 1 , NDF
         DO i = 1 , NUMnp
            temp(Num(i)) = F(j,i)
         ENDDO
         DO i = 1 , NUMnp
            F(j,i) = temp(i)
         ENDDO
      ENDDO
!
      DO j = 1 , NDF
         DO i = 1 , NUMnp
            temp(Num(i)) = B(j,i)
         ENDDO
         DO i = 1 , NUMnp
            B(j,i) = temp(i)
         ENDDO
      ENDDO
!
      DO j = 1 , NEN
         DO i = 1 , NUMel
            Ix(j,i) = Num(Ix(j,i))
         ENDDO
      ENDDO
!
!
! for constrained delaunay
!
      IF ( NCE>0 ) THEN
         DO i = 1 , NCE
            DO j = 1 , 2
               ii = ELIst(j,i)
               ELIst(j,i) = Num(ii)
            ENDDO
         ENDDO
      ENDIF
 
      DO i = 1 , NUMnp
         itemp(Num(i)) = ISRelaxed(i)
      ENDDO
      DO i = 1 , NUMnp
         ISRelaxed(i) = itemp(i)
      ENDDO
      DO i = 1 , NUMnp
         itemp(Num(i)) = ATOmspecie(i)
      ENDDO
      DO i = 1 , NUMnp
         ATOmspecie(i) = itemp(i)
      ENDDO
!
      DO i = 1 , NUMnp
         temp(Num(i)) = ENErgy(i)
      ENDDO
      DO i = 1 , NUMnp
         ENErgy(i) = temp(i)
      ENDDO
!
      DEALLOCATE (temp,itemp)
      END SUBROUTINE SWAPALLGR
!*==dgree.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE DGREE(Ndstk,Nr,Ndeg,Iold,Ibw1,Ipf1)
      IMPLICIT NONE
!*--DGREE227
!*** Start of declarations inserted by SPAG
      INTEGER i , Ibw1 , IDEg , idif , IDPth , Iold , Ipf1 , irw , &
     &        itst , j , N , Ndeg , Nr
!*** End of declarations inserted by SPAG
!
!
!  dgree computes the degree of each node in ndstk and stores
!  it in the array ndeg.  The bandwidth and profile for the original
!  or input renumbering of the graph is computed also.
!  Use integer*2 ndstk  with an ibm 360 or 370.
!
!
      INTEGER Ndstk
!
      COMMON /GRA   / N , IDPth , IDEg
!
      DIMENSION Ndstk(Nr,*) , Ndeg(*) , Iold(*)
!
      Ibw1 = 0
      Ipf1 = 0
!
      DO i = 1 , N
         Ndeg(i) = 0
         irw = 0
         DO j = 1 , IDEg
            itst = Ndstk(j,i)
            IF ( itst<=0 ) EXIT
            Ndeg(i) = Ndeg(i) + 1
            idif = Iold(i) - Iold(itst)
            IF ( irw<idif ) irw = idif
         ENDDO
         Ipf1 = Ipf1 + irw
         IF ( irw>Ibw1 ) Ibw1 = irw
      ENDDO
!
!
      END SUBROUTINE DGREE
!*==fndiam.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FNDIAM(Snd1,Snd2,Ndstk,Nr,Ndeg,Lvl,Lvls1,Lvls2,Iwk,&
     &                  Idflt)
      IMPLICIT NONE
!*--FNDIAM272
!*** Start of declarations inserted by SPAG
      INTEGER i , IDEg , Idflt , IDPth , Iwk , Lvl , lvlbot , lvln , &
     &        Lvls1 , Lvls2 , lvlwth , maxlw , mtw1 , mtw2 , N , Ndeg , &
     &        NDLst , ndxl , ndxn , Nr
!*** End of declarations inserted by SPAG
!
!
!  fndiam is the control procedure for finding the pseudo-diameter of
!  ndstk as well as the level structure from each end
!  snd1-        on input this is the node number of the first
!               attempt at finding a diameter.  on output it
!               contains the actual number used.
!  snd2-        on output contains other end of diameter
!  lvls1-       array containing level structure with snd1 as root
!  lvls2-       array containing level structure with snd2 as root
!  idflt-       flag used in picking final level structure, set
!               =1 if width of lvls1 .le. width of lvls2, otherwise =2
!  lvl,iwk-     working storage
!  Use integer*2 ndstk  with an ibm 360 or 370.
!
      INTEGER Ndstk
      INTEGER flag , snd , Snd1 , Snd2
!
      COMMON /GRA   / N , IDPth , IDEg
!
!  It is assumed that the last level has at most 'maxlvl' nodes.
!
      INTEGER MAXLVL
      PARAMETER (MAXLVL=2000)
      COMMON /CC    / NDLst(MAXLVL)
      DIMENSION Ndstk(Nr,*) , Ndeg(*) , Lvl(*) , Lvls1(*) , Lvls2(*) , &
     &          Iwk(*)
!
      flag = 0
      mtw2 = N
      snd = Snd1
 100  DO
!
!  Zero lvl to indicate all nodes are available to tree.
!
         DO i = 1 , N
            Lvl(i) = 0
         ENDDO
         lvln = 1
!
!  Drop a tree from snd.
!
         CALL TREE(snd,Ndstk,Nr,Lvl,Iwk,Ndeg,lvlwth,lvlbot,lvln,maxlw,&
     &             mtw2)
         IF ( flag<1 ) THEN
!
            flag = 1
            EXIT
         ELSEIF ( IDPth>=lvln-1 ) THEN
            IF ( maxlw<mtw2 ) THEN
               mtw2 = maxlw
               Snd2 = snd
!
!  Store narrowest reverse level structure in lvls2.
!
               DO i = 1 , N
                  Lvls2(i) = Lvl(i)
               ENDDO
            ENDIF
 
            IF ( ndxn==ndxl ) THEN
               Idflt = 1
               IF ( mtw2<=mtw1 ) Idflt = 2
               GOTO 99999
            ELSE
!
!  Try next node in ndlst.
!
               ndxn = ndxn + 1
               IF ( ndxn>MAXLVL ) THEN
                  PRINT * , &
     &            '***ERROR: Insufficient storage in fndiam for ndlst()'
                  STOP
               ENDIF
               snd = NDLst(ndxn)
            ENDIF
         ELSE
!
!  Start again with new starting node.
!
            Snd1 = snd
            EXIT
         ENDIF
      ENDDO
      IDPth = lvln - 1
      mtw1 = maxlw
!
!  Copy level structure into lvls1.
!
      DO i = 1 , N
         Lvls1(i) = Lvl(i)
      ENDDO
      ndxn = 1
      ndxl = 0
      mtw2 = N
!
!  Sort last level by degree  and store in ndlst.
!
      CALL SORTDG(NDLst,Iwk(lvlbot),ndxl,lvlwth,Ndeg)
      snd = NDLst(1)
      GOTO 100
!
!
99999 END SUBROUTINE FNDIAM
!*==number.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE NUMBER(Snd,Num,Ndstk,Lvls2,Ndeg,Renum,Lvlst,Lstpt,Nr,&
     &                  Nflg,Ibw2,Ipf2,Ipfa,Isdir)
      IMPLICIT NONE
!*--NUMBER389
!*** Start of declarations inserted by SPAG
      INTEGER i , Ibw2 , IDEg , IDPth , inx , Ipf2 , Ipfa , ipro , &
     &        Isdir , j , lnd , lst , Lstpt , lvln , Lvls2 , Lvlst , &
     &        max , N , nbw , Ndeg
      INTEGER Nflg , Nr , nstpt , Num
!*** End of declarations inserted by SPAG
!
!  Number produces the numbering of the graph for min bandwidth
!  snd-         on input the node to begin numbering on
!  num-         on input and output, the next available number
!  lvls2-       the level structure to be used in numbering
!  renum-       the array used to store the new numbering
!  lvlst-       on output contains level structure
!  lstpt(i)-    on output, index into lvlst to first node in ith lvl
!               lstpt(i+1) - lstpt(i) = number of nodes in ith lvl
!  nflg-        =+1 if snd is forward end of pseudo-diam
!               =-1 if snd is reverse end of pseudo-diam
!  ibw2-        bandwidth of new numbering computed by number
!  ipf2-        profile of new numbering computed by number
!  ipfa-        working storage used to compute profile and bandwidth
!  isdir-       indicates step direction used in numbering(+1 or -1)
!  Use integer*2 ndstk  with an ibm 360 or 370.
!
      INTEGER Ndstk
      INTEGER Snd , STKa , STKb , STKc , STKd , xa , xb , xc , xd , cx ,&
     &        end , Renum , test
!
      COMMON /GRA   / N , IDPth , IDEg
!
!  The storage in common blocks cc and lvlw is now free and can
!  be used for stacks.
!
      INTEGER MAXLVL
      PARAMETER (MAXLVL=2000)
      COMMON /LVLW  / STKa(MAXLVL) , STKb(MAXLVL) , STKc(MAXLVL)
      COMMON /CC    / STKd(MAXLVL)
!
      DIMENSION Ipfa(*)
      DIMENSION Ndstk(Nr,*) , Lvls2(*) , Ndeg(*) , Renum(*) , Lvlst(*) ,&
     &          Lstpt(*)
!
!  Set up lvlst and lstpt from lvls2.
!
      DO i = 1 , N
         Ipfa(i) = 0
      ENDDO
      nstpt = 1
      DO i = 1 , IDPth
         Lstpt(i) = nstpt
         DO j = 1 , N
            IF ( Lvls2(j)==i ) THEN
               Lvlst(nstpt) = j
               nstpt = nstpt + 1
            ENDIF
         ENDDO
      ENDDO
      Lstpt(IDPth+1) = nstpt
!
!  stka, stkb, stkc and stkd are stacks with pointers
!  xa,xb,xc, and xd.  cx is a special pointer into stkc which
!  indicates the particular node being processed.
!  lvln keeps track of the level we are working at.
!  initially stkc contains only the initial node, snd.
!
      lvln = 0
      IF ( Nflg<0 ) lvln = IDPth + 1
      xc = 1
      STKc(xc) = Snd
 100  cx = 1
      xd = 0
      lvln = lvln + Nflg
      lst = Lstpt(lvln)
      lnd = Lstpt(lvln+1) - 1
!
!  Begin processing node stkc(cx).
!
 200  ipro = STKc(cx)
      Renum(ipro) = Num
      Num = Num + Isdir
      end = Ndeg(ipro)
      xa = 0
      xb = 0
!
!  Check all adjacent nodes.
!
      DO i = 1 , end
         test = Ndstk(i,ipro)
         inx = Renum(test)
!
!  Only nodes not numbered or already on a stack are added.
!
         IF ( inx==0 ) THEN
            Renum(test) = -1
!
!  Put nodes on same level on stka, all others on stkb.
!
            IF ( Lvls2(test)==Lvls2(ipro) ) THEN
               xa = xa + 1
               IF ( xa>MAXLVL ) THEN
                  PRINT * , &
     &             '***ERROR: Insufficient storage in subroutine number'
                  PRINT * , '           for stka()'
                  STOP
               ENDIF
               STKa(xa) = test
            ELSE
               xb = xb + 1
               IF ( xb>MAXLVL ) THEN
                  PRINT * , &
     &             '***ERROR: Insufficient storage in subroutine number'
                  PRINT * , '           for stkb()'
                  STOP
               ENDIF
               STKb(xb) = test
            ENDIF
         ELSEIF ( inx>=0 ) THEN
!
!  Do preliminary bandwidth and profile calculations.
!
            nbw = (Renum(ipro)-inx)*Isdir
            IF ( Isdir>0 ) inx = Renum(ipro)
            IF ( Ipfa(inx)<nbw ) Ipfa(inx) = nbw
         ENDIF
      ENDDO
!
!  Sort stka and stkb into increasing degree and add stka to stkc
!  and stkb to stkd.
!
      IF ( xa/=0 ) THEN
         IF ( xa==1 ) THEN
            xc = xc + 1
            IF ( xc>MAXLVL ) THEN
               PRINT * , &
     &             '***ERROR: Insufficient storage in subroutine number'
               PRINT * , '           for stkc()'
               STOP
            ENDIF
            STKc(xc) = STKa(xa)
         ELSE
            CALL SORTDG(STKc,STKa,xc,xa,Ndeg)
         ENDIF
      ENDIF
      IF ( xb/=0 ) THEN
         IF ( xb==1 ) THEN
            xd = xd + 1
            IF ( xd>MAXLVL ) THEN
               PRINT * , &
     &             '***ERROR: Insufficient storage in subroutine number'
               PRINT * , '           for stkd()'
               STOP
            ENDIF
            STKd(xd) = STKb(xb)
         ELSE
            CALL SORTDG(STKd,STKb,xd,xb,Ndeg)
         ENDIF
      ENDIF
!
!  Be sure to process all nodes in stkc.
!
      cx = cx + 1
      IF ( cx>MAXLVL ) THEN
         PRINT * , '***ERROR: Insufficient storage in subroutine number'
         PRINT * , '           for stkc() (index cx)'
         STOP
      ENDIF
      IF ( xc>=cx ) GOTO 200
!
!  When stkc is exhausted look for min degree node in same level
!  which has not been processed.
!
      max = IDEg + 1
      Snd = N + 1
      DO i = lst , lnd
         test = Lvlst(i)
         IF ( Renum(test)==0 ) THEN
            IF ( Ndeg(test)<max ) THEN
               Renum(Snd) = 0
               Renum(test) = -1
               max = Ndeg(test)
               Snd = test
            ENDIF
         ENDIF
      ENDDO
      IF ( Snd/=N+1 ) THEN
         xc = xc + 1
         IF ( xc>MAXLVL ) THEN
            PRINT * , &
     &            '***ERROR: Insufficient storage in subroutine number'
            PRINT * , '           for stkc()'
            STOP
         ENDIF
         STKc(xc) = Snd
         GOTO 200
!
!  If stkd is empty we are done, otherwise copy stkd onto stkc
!  and begin processing new stkc.
!
      ELSEIF ( xd==0 ) THEN
!
!  Do final bandwidth and profile calculations.
!
         DO i = 1 , N
 
            IF ( Ipfa(i)>Ibw2 ) Ibw2 = Ipfa(i)
            Ipf2 = Ipf2 + Ipfa(i)
         ENDDO
      ELSE
         IF ( xd>MAXLVL ) THEN
            PRINT * , &
     &            '***ERROR: Insufficient storage in subroutine number'
            PRINT * , '           for stkc()'
            STOP
         ENDIF
         DO i = 1 , xd
            STKc(i) = STKd(i)
         ENDDO
         xc = xd
         GOTO 100
      ENDIF
!
!
      END SUBROUTINE NUMBER
!*==piklvl.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE PIKLVL(Lvls1,Lvls2,Ccstor,Idflt,Isdir)
      IMPLICIT NONE
!*--PIKLVL618
!*** Start of declarations inserted by SPAG
      INTEGER i , IDEg , Idflt , IDPth , inode , Isdir , it , j , k , &
     &        lvlnh , lvlnl , Lvls1 , Lvls2 , max1 , max2 , MAXCMP , N ,&
     &        NACum , NHIgh , NLOw
!*** End of declarations inserted by SPAG
!
! piklvl chooses the level structure  used in numbering graph
! lvls1-       on input contains forward leveling info
! lvls2-       on input contains reverse leveling info
!              on output the final level structure chosen
! ccstor-      on input contains connected component info
! idflt-       on input =1 if wdth lvls1.le.wdth lvls2, =2 otherwise
! nhigh        keeps track of level widths for high numbering
! nlow-        keeps track of level widths for low numbering
! nacum-       keeps track of level widths for chosen level structure
! xc-          number of connected components
! size(i)-     size of ith connected component
! stpt(i)-     index into ccstore of 1st node in ith con compt
! isdir-       flag which indicates which way the largest connected
!              component fell.  =+1 if low and -1 if high
!
      INTEGER Ccstor , SIZe , STPt , XC , end
!
      COMMON /GRA   / N , IDPth , IDEg
!
!  It is assumed that the graph has at most 'maxcmp' (was 50) components
!  that there are at most 'maxlvl' (was 100) levels.
!
      INTEGER MAXLVL
      PARAMETER (MAXLVL=2000,MAXCMP=400) ! Note: maxcmp must be less tha
                                         !       (maxlvl-1)/2
      COMMON /LVLW  / NHIgh(MAXLVL) , NLOw(MAXLVL) , NACum(MAXLVL)
      COMMON /CC    / XC , SIZe(MAXCMP) , STPt(MAXCMP)
!     common / ccc  / xc, size(maxcmp), stpt(maxcmp)
!
      DIMENSION Lvls1(*) , Lvls2(*) , Ccstor(*)
!
!  For each connected component do.
!
      IF ( XC>MAXCMP ) THEN
         PRINT * , '***ERROR: Insufficient storage in piklvl for'
         PRINT * , '          stpt() and size()'
         STOP
      ENDIF
      DO i = 1 , XC
         j = STPt(i)
         end = SIZe(i) + j - 1
!
!  Set nhigh and nlow equal to nacum.
!
         IF ( IDPth>MAXLVL ) THEN
            PRINT * , '***ERROR: Insufficient storage in piklvl for'
            PRINT * , '          nacum(), nhigh() and nlow()'
            STOP
         ENDIF
         DO k = 1 , IDPth
            NHIgh(k) = NACum(k)
            NLOw(k) = NACum(k)
         ENDDO
!
!  Update nhigh and nlow for each node in connected component.
!
         DO k = j , end
            inode = Ccstor(k)
            lvlnh = Lvls1(inode)
            IF ( lvlnh>MAXLVL ) THEN
               PRINT * , '***ERROR: pointer in lvls1 out-of-bounds'
               STOP
            ENDIF
            NHIgh(lvlnh) = NHIgh(lvlnh) + 1
            lvlnl = Lvls2(inode)
            IF ( lvlnl>MAXLVL ) THEN
               PRINT * , '***ERROR: pointer in lvls2 out-of-bounds'
               STOP
            ENDIF
            NLOw(lvlnl) = NLOw(lvlnl) + 1
         ENDDO
         max1 = 0
         max2 = 0
!
!  Set max1=largest new number in nhigh.
!  Set max2=largest new number in nlow.
!
         DO k = 1 , IDPth
            IF ( 2*NACum(k)/=NLOw(k)+NHIgh(k) ) THEN
               IF ( NHIgh(k)>max1 ) max1 = NHIgh(k)
               IF ( NLOw(k)>max2 ) max2 = NLOw(k)
            ENDIF
         ENDDO
!
!  Set it= number of level structure to be used.
!
         it = 1
         IF ( max1>max2 ) it = 2
         IF ( max1==max2 ) it = Idflt
         IF ( it==2 ) THEN
!
!  Update nacum to be the same as nlow.
!
            DO k = 1 , IDPth
               NACum(k) = NLOw(k)
            ENDDO
         ELSE
            IF ( i==1 ) Isdir = -1
!
!  Copy lvls1 into lvls2 for each node in connected component.
!
            DO k = j , end
               inode = Ccstor(k)
               Lvls2(inode) = Lvls1(inode)
            ENDDO
!
!  Update nacum to be the same as nhigh.
!
            DO k = 1 , IDPth
               NACum(k) = NHIgh(k)
            ENDDO
         ENDIF
      ENDDO
!
!
      END SUBROUTINE PIKLVL
!*==reduce.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE REDUCE(Ndstk,Nr,Iold,Renum,Ndeg,Lvl,Lvls1,Lvls2,Ccstor,&
     &                  Ibw2,Ipf2)
      IMPLICIT NONE
!*--REDUCE748
!*** Start of declarations inserted by SPAG
      REAL dmy
      INTEGER i , ibw1 , Ibw2 , IDEg , idflt , IDPth , Iold , ipf1 , &
     &        Ipf2 , isdir , lowdg , lroot , Lvl , lvlbot , lvln , &
     &        Lvls1 , Lvls2 , lvlwth , MAXCMP , MAXLVL
      INTEGER maxlw , N , NACum , Ndeg , nflg , NHIgh , NLOw , Nr , num
!*** End of declarations inserted by SPAG
!
!  Subroutine reduce determines a row and column permutation which,
!  when applied to a given sparse matrix, produces a permuted
!  matrix with a smaller bandwidth and profile.
!  The input array is a connection table which represents the
!  indices of the nonzero elements of the matrix, a.  The algo-
!  rithm is described in terms of the adjacency graph which
!  has the characteristic that there is an edge (connection)
!  between nodes i and j if a(i,j) .ne. 0 and i .ne. j.
!
!  Dimensioning information:
!  The following integer arrays must be dimensioned in the calling routi
!    ndstk(nr,d1)        nr is .ge. maximum degree of all nodes.
!    iold(d2)            d2 and d1 are .ge. the total number of
!    renum(d2+1)         nodes in the graph.
!    ndeg(d2)            storage requirements can be significantly
!    lvl(d2)             decreased for ibm 360 and 370 computers
!    lvls1(d2)           by replacing integer ndstk by
!    lvls2(d2)           integer*2 ndstk in subroutines reduce,
!    ccstor(d2)          dgree, fndiam, tree and number.
!
!  Common information:
!  The following common block must be in the calling routine.
!
!    common / gra / n, idpth, ideg
!
!  Explanation of input variables:
!    ndstk-     connection table representing graph.
!               ndstk(i,j)=node number of ith connection to node
!               number-j.  a connection of a node to itself is not
!               listed.  extra positions must have zero fill.
!    nr-        row dimension assigned ndstk in calling program.
!               this is the maximum number of node point connections
!               allowed in the graph.
!    iold(i)-   numbering of ith node upon input.
!               if no numbering exists then iold(i)=i.
!    n-         number of nodes in graph (equal to order of matrix).
!    ideg-      maximum degree of any node in the graph.
!
!  Explanation of output variables:
!    renum(i)-  the new number for the ith node.
!    ndeg(i)-   the degree of the ith node.
!    ibw2-      the bandwidth after renumbering.
!    ipf2-      the profile after renumbering.
!    idpth-     number of levels in reduce level structure.
!  The following only have meaning if the graph was connected--
!    lvl(i)-    index into lvls1 to the first node in level i.
!               lvl(i+1)-lvl(i)= number of nodes in ith level
!    lvls1(i)-  node numbers listed by level.
!    lvls2(i)-  the level assigned to node i by reduce.
!
!  Working storage variable:
!    ccstor
!
!  Local storage:
!    common/cc/-subroutines reduce, sort2 and piklvl assume that
!               the graph has at most 100 (was 50) connected components.
!               subroutine fndiam assumes that there are at most
!               500 (was 100) nodes in the last level.
!    common/lvlw/-subroutines setup and piklvl assume that there
!               are at most 500 (was 100) levels.
!
!  Use integer*2 ndstk  with an ibm 360 or 370.
!
      INTEGER Ndstk
      INTEGER stnode , rvnode , Renum , XC , SORT2 , stnum , Ccstor , &
     &        SIZe , STPt , sbnum
!
      COMMON /GRA   / N , IDPth , IDEg
!
!  It is assumed that the graph has at most 'maxcmp' (was 50) connected
!  components.
!
      PARAMETER (MAXLVL=2000,MAXCMP=400) ! Note: maxcmp must be less tha
                                         !       (maxlvl-1)/2
      COMMON /CC    / XC , SIZe(MAXCMP) , STPt(MAXCMP)
!     common / ccc  / xc, size(maxcmp), stpt(maxcmp)
      COMMON /LVLW  / NHIgh(MAXLVL) , NLOw(MAXLVL) , NACum(MAXLVL)
!
      DIMENSION Ccstor(*) , Iold(*)
      DIMENSION Ndstk(Nr,*) , Lvl(*) , Lvls1(*) , Lvls2(*) , Renum(*) , &
     &          Ndeg(*)
!
      Ibw2 = 0
      Ipf2 = 0
!
!  Set renum(i)=0 for all i to indicate node i is unnumbered.
!
      DO i = 1 , N
         Renum(i) = 0
      ENDDO
!
!  Compute degree of each node and original bandwidth and profile.
!
      CALL DGREE(Ndstk,Nr,Ndeg,Iold,ibw1,ipf1)
!
!  Display the original bandwidth and profile.
!
!     write(6,'(x,''..original bandwidth  = '',i10,
!    +            /,x,''..original profile    = '',i10)')
!    +            ibw1, ipf1
!
!  sbnum = low end of available numbers for renumbering
!  stnum = high end of available numbers for renumbering
!
      sbnum = 1
      stnum = N
!
!  Number the nodes of degree zero.
!
      DO i = 1 , N
         IF ( Ndeg(i)<=0 ) THEN
            Renum(i) = stnum
            stnum = stnum - 1
         ENDIF
      ENDDO
!
!  Find an unnumbered node of min degree to start on.
!
 100  lowdg = IDEg + 1
      nflg = 1
      isdir = 1
      DO i = 1 , N
         IF ( Ndeg(i)<lowdg ) THEN
            IF ( Renum(i)<=0 ) THEN
               lowdg = Ndeg(i)
               stnode = i
            ENDIF
         ENDIF
      ENDDO
!
!  Find pseudo-diameter and associated level structures.
!  stnode and rvnode are the ends of the diam and lvls1 and lvls2
!  are the respective level structures.
!
      CALL FNDIAM(stnode,rvnode,Ndstk,Nr,Ndeg,Lvl,Lvls1,Lvls2,Ccstor,&
     &            idflt)
      IF ( Ndeg(stnode)>Ndeg(rvnode) ) THEN
!
!  nflg indicates the end to begin numbering on
!
         nflg = -1
         stnode = rvnode
      ENDIF
      CALL SETUP(Lvl,Lvls1,Lvls2)
!
!  Find all the connected components  (xc counts them).
!
      XC = 0
      lroot = 1
      lvln = 1
      DO i = 1 , N
         IF ( Lvl(i)==0 ) THEN
            XC = XC + 1
            IF ( XC>MAXCMP ) THEN
               PRINT * , &
     &             '***ERROR: Insufficient storage in reduce for stpt()'
               STOP
            ENDIF
            STPt(XC) = lroot
            CALL TREE(i,Ndstk,Nr,Lvl,Ccstor,Ndeg,lvlwth,lvlbot,lvln,&
     &                maxlw,N)
            SIZe(XC) = lvlbot + lvlwth - lroot
            lroot = lvlbot + lvlwth
            lvln = lroot
         ENDIF
      ENDDO
      IF ( SORT2(dmy)/=0 ) CALL PIKLVL(Lvls1,Lvls2,Ccstor,idflt,isdir)
!
!  On return from piklvl, isdir indicates the direction the largest
!  component fell.  isdir is modified now to indicate the numbering
!  direction.  num is set to the proper value for this direction.
!
      isdir = isdir*nflg
      num = sbnum
      IF ( isdir<0 ) num = stnum
      CALL NUMBER(stnode,num,Ndstk,Lvls2,Ndeg,Renum,Lvls1,Lvl,Nr,nflg,&
     &            Ibw2,Ipf2,Ccstor,isdir)
!
!  Update stnum or sbnum after numbering.
!
      IF ( isdir<0 ) stnum = num
      IF ( isdir>0 ) sbnum = num
      IF ( sbnum<=stnum ) GOTO 100
!
!     write(6,'(x,''..minimized bandwidth = '',i10,
!    +            /,x,''..minimized profile   = '',i10)')
!    +            ibw2, ipf2
!
!  If the original bandwidth is the same as the minimum bandwidth
!  keep the old labeling and return.
!
!      write(6,'(x,''..original bandwidth is minimized'')')
      IF ( Ibw2<=ibw1 ) RETURN
!
!  If original numbering is better than new one, set up to return it.
!
      DO i = 1 , N
         Renum(i) = Iold(i)
      ENDDO
      Ibw2 = ibw1
      Ipf2 = ipf1
!
!
      END SUBROUTINE REDUCE
!*==setup.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SETUP(Lvl,Lvls1,Lvls2)
      IMPLICIT NONE
!*--SETUP967
!*** Start of declarations inserted by SPAG
      INTEGER i , IDEg , IDPth , itemp , Lvl , Lvls1 , Lvls2 , MAXLVL , &
     &        N , NACum , NHIgh , NLOw
!*** End of declarations inserted by SPAG
!
! Setup computes the reverse leveling info from lvls2 and stores
! it into lvls2.  nacum(i) is initialized to nodes/ith level for nodes
! on the pseudo-diameter of the graph.  lvl is initialized to non-
! zero for nodes on the pseudo-diam and nodes in a different
! component of the graph.
!
      COMMON /GRA   / N , IDPth , IDEg
!
! It is assumed that there are at most 'maxlvl' (was 100) levels.
!
      PARAMETER (MAXLVL=2000)
      COMMON /LVLW  / NHIgh(MAXLVL) , NLOw(MAXLVL) , NACum(MAXLVL)
!
      DIMENSION Lvl(*) , Lvls1(*) , Lvls2(*)
!
      IF ( IDPth>MAXLVL ) THEN
         PRINT * , '***ERROR: Insufficient storage in setup for nacum()'
         STOP
      ENDIF
      DO i = 1 , IDPth
         NACum(i) = 0
      ENDDO
!
      DO i = 1 , N
         Lvl(i) = 1
         Lvls2(i) = IDPth + 1 - Lvls2(i)
         itemp = Lvls2(i)
         IF ( itemp<=IDPth ) THEN
            IF ( itemp/=Lvls1(i) ) THEN
               Lvl(i) = 0
            ELSE
               NACum(itemp) = NACum(itemp) + 1
            ENDIF
         ENDIF
      ENDDO
!
!
      END SUBROUTINE SETUP
!*==sort2.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      INTEGER FUNCTION SORT2(Dmy)
      IMPLICIT NONE
!*--SORT21016
!*** Start of declarations inserted by SPAG
      REAL Dmy
      INTEGER i , ind , itest , j , MAXCMP , MAXLVL
!*** End of declarations inserted by SPAG
!
! Sort2 sorts size and stpt into descending order according to
! values of size. xc=number of entries in each array.
!
      INTEGER temp , ccstor , SIZe , STPt , XC
!
! It is assumed that the graph has at most 'maxcmp' (was 50) connected
! components.
!
      PARAMETER (MAXLVL=2000,MAXCMP=400) ! Note: maxcmp must be less tha
                                         !       (maxlvl-1)/2
      COMMON /CC    / XC , SIZe(MAXCMP) , STPt(MAXCMP)
!     common / ccc  / xc, size(maxcmp), stpt(maxcmp)
!
      SORT2 = 0
      IF ( XC==0 ) RETURN
      IF ( XC>MAXCMP ) THEN
         PRINT * , '***ERROR: Insufficient storage in sort2 for size()'
         PRINT * , '          and stpt()'
         STOP
      ENDIF
!
      SORT2 = 1
      ind = XC
      DO
         itest = 0
         ind = ind - 1
         IF ( ind<1 ) RETURN
!
         DO i = 1 , ind
            j = i + 1
            IF ( SIZe(i)<SIZe(j) ) THEN
               itest = 1
               temp = SIZe(i)
               SIZe(i) = SIZe(j)
               SIZe(j) = temp
               temp = STPt(i)
               STPt(i) = STPt(j)
               STPt(j) = temp
            ENDIF
         ENDDO
!
         IF ( itest/=1 ) EXIT
      ENDDO
!
!
      END FUNCTION SORT2
!*==sortdg.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SORTDG(Stk1,Stk2,X1,X2,Ndeg)
      IMPLICIT NONE
!*--SORTDG1074
!*** Start of declarations inserted by SPAG
      INTEGER i , IDEg , IDPth , ind , istk2 , itest , j , jstk2 , N , &
     &        Ndeg
!*** End of declarations inserted by SPAG
! sortdg sorts stk2 by degree of the node and adds it to the end
! of stk1 in order of lowest to highest degree.  x1 and x2 are the
! number of nodes in stk1 and stk2 respectively.
      INTEGER X1 , X2 , Stk1 , Stk2 , temp
      COMMON /GRA   / N , IDPth , IDEg
      DIMENSION Ndeg(1) , Stk1(1) , Stk2(1)
      ind = X2
      DO
         itest = 0
         ind = ind - 1
         IF ( ind<1 ) EXIT
         DO i = 1 , ind
            j = i + 1
            istk2 = Stk2(i)
            jstk2 = Stk2(j)
            IF ( Ndeg(istk2)>Ndeg(jstk2) ) THEN
               itest = 1
               temp = Stk2(i)
               Stk2(i) = Stk2(j)
               Stk2(j) = temp
            ENDIF
         ENDDO
         IF ( itest/=1 ) EXIT
      ENDDO
      DO i = 1 , X2
         X1 = X1 + 1
         Stk1(X1) = Stk2(i)
      ENDDO
      END SUBROUTINE SORTDG
!*==tree.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TREE(Iroot,Ndstk,Nr,Lvl,Iwk,Ndeg,Lvlwth,Lvlbot,Lvln,&
     &                Maxlw,Ibort)
      IMPLICIT NONE
!*--TREE1115
!*** Start of declarations inserted by SPAG
      INTEGER Ibort , inow , Iroot , itest , itop , Iwk , iwknow , j , &
     &        Lvl , Lvlbot , Lvln , lvltop , Lvlwth , Maxlw , Ndeg , &
     &        ndrow , Nr
!*** End of declarations inserted by SPAG
!
!  tree drops a tree in ndstk from iroot.
!  lvl-         array indicating available nodes in ndstk with zero
!               entries. tree enters level numbers assigned
!               during execution of this procedure
!  iwk-         on output contains node numbers used in tree
!               arranged by levels (iwk(lvln) contains iroot
!               and iwk(lvlbot+lvlwth-1) contains last node entered)
!  lvlwth-      on output contains width of last level
!  lvlbot-      on output contains index into iwk of first
!               node in last level
!  maxlw-       on output contains the maximum level width
!  lvln-        on input the first available location in iwk
!               usually one but if iwk is used to store previous
!               connected components, lvln is next available location.
!               on output the total number of levels + 1
!  ibort-       input param which triggers early return if
!               maxlw becomes .ge. ibort
! use integer*2 ndstk  with an ibm 360 or 370.
!
      INTEGER Ndstk
!
      DIMENSION Ndstk(Nr,*) , Lvl(*) , Iwk(*) , Ndeg(*)
!
      Maxlw = 0
      itop = Lvln
      inow = Lvln
      Lvlbot = Lvln
      lvltop = Lvln + 1
      Lvln = 1
      Lvl(Iroot) = 1
      Iwk(itop) = Iroot
 100  Lvln = Lvln + 1
      DO
         iwknow = Iwk(inow)
         ndrow = Ndeg(iwknow)
!
         DO j = 1 , ndrow
            itest = Ndstk(j,iwknow)
            IF ( Lvl(itest)==0 ) THEN
               Lvl(itest) = Lvln
               itop = itop + 1
               Iwk(itop) = itest
            ENDIF
         ENDDO
!
         inow = inow + 1
         IF ( inow>=lvltop ) THEN
            Lvlwth = lvltop - Lvlbot
            IF ( Maxlw<Lvlwth ) Maxlw = Lvlwth
            IF ( Maxlw>=Ibort ) RETURN
            IF ( itop<lvltop ) RETURN
            Lvlbot = inow
            lvltop = itop + 1
            GOTO 100
         ENDIF
      ENDDO
!
!
      END SUBROUTINE TREE
 
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
