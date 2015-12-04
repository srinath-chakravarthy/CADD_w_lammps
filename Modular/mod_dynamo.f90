!*==mod_dynamo.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
      MODULE MOD_DYNAMO
      IMPLICIT NONE
!*--MOD_DYNAMO4
!
!     This use is for use with the "dynamo" solver only.
!     It is the part of dyn84_2.inc that is needed by the dynamo solver,
!     but is not already used in poten_eam.inc
!
!.....$Modified: Wed Jan 1810:27:091995 by foiles $
      INTEGER NEIMAX
      PARAMETER (NEIMAX=300)
      
      !> NEIMAX -> Max number of neighbors
      DOUBLE PRECISION rneigh(NEIMAX) , dneigh(3,NEIMAX)
      !> rneigh --> distance squared from atom i to neighbor
      !> dneigh --> Vector from atom i to neighbor

      INTEGER jneigh(NEIMAX)
      !> jneigh ---> atom numbers in neighbor list
      
      DOUBLE PRECISION , ALLOCATABLE :: rold(:,:) , dis(:,:) , rdyn(:)

      INTEGER , ALLOCATABLE :: nnindx(:) , nnlst(:,:) , knbr(:) ,  neighborlist(:,:) , numneighbors(:)
      DOUBLE PRECISION perub(3) , perlb(3) , perlen(3) , alat , xbound ,  ybound(2) , zbound(2)
      LOGICAL onlyoi , twoi , threei , fouri
      INTEGER ngtlst , nneimx , mxlstu , mxnnei , nforce
      DOUBLE PRECISION rctsqn , dradn
      INTEGER nneigh , nneips , nmeth , newlst , nlstmx
      CONTAINS
 
!**********************************************************************
      DOUBLE PRECISION FUNCTION RV2(I,J,B,X)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--RV230
      INTEGER I , J
      DOUBLE PRECISION B(NDF,NUMnp) , X(NXDm,NUMnp)
      RV2 = B(I,J) + X(I,J)
      END FUNCTION RV2
!**********************************************************************
 
!.....
!.....this routine determines the neighbors of a given atom
!.....
      SUBROUTINE GNEIGH(I,B,X,Needlist)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--GNEIGH43
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION B , CUTOFFR2 , disz , rcutsq , rcutsq2 , rim , X

      INTEGER I , iatom , j , jend , jtmp , kcoord
!*** End of declarations inserted by SPAG
      DIMENSION B(NDF,*) , X(NXDm,*)
      INTEGER totneighbors
      LOGICAL Needlist , firsttime
      DATA firsttime/.TRUE./
 
!--   Vijay
      IF ( firsttime==.TRUE. ) THEN
         firsttime = .FALSE.
         ALLOCATE (NEIghborlist(MAXneighbors,NUMnp))
         ALLOCATE (NUMneighbors(NUMnp))
 
         DO iatom = 1 , NUMnp
            !print*, 'iAtom: ', iAtom
            DO j = 1 , MAXneighbors
               NEIghborlist(j,iatom) = 0
            ENDDO
            NUMneighbors(iatom) = 0
         ENDDO
      ENDIF
!--   Vijay
 
 
!      print*, 'In gNeigh: iAtom: ', i
      IF ( I==NUMnpp1 ) THEN
         rcutsq = CUTOFFR2(-1)
      ELSE
         rcutsq = CUTOFFR2(1)
      ENDIF
      RCTsqn = (SQRT(rcutsq)+DRAdn)**2
 
 
      totneighbors = 0
!.....
!.....beginning of nmeth = 2
!.....storage of neighbor indices
!.....
!.....newlst specifies whether the neighbor list is being updated
!.....
      IF ( SIMstep==1 ) THEN
         IF ( NEWlst/=0 ) GOTO 100
      ENDIF
 
!     C--Jun Song New: only update those atoms sepcified in chkdis
      IF ( SIMstep>1 ) THEN
         IF ( (NEWlst/=0) .AND. (UPDateneigh(I)==1) ) GOTO 100
      ENDIF
!     c-end
 
!.....
!.....if here then use the neighbors found earlier
!.....
      IF ( .NOT.Needlist ) RETURN
!     c-JS         jend = nnindx(i) - nnindx(i-1)
      jend = NNIndx(I)
      DO jtmp = 1 , jend
!     c-JS            knbr(jtmp) = nnlst(jtmp+nnindx(i-1))
         KNBr(jtmp) = NNLst(I,jtmp)
 
         DIS(1,jtmp) = RV2(1,I,B,X) - RV2(1,KNBr(jtmp),B,X)
         DIS(1,jtmp) = DIS(1,jtmp) - PERlen(1)&
     &                 *NINT(DIS(1,jtmp)/PERlen(1))
         RDYn(jtmp) = DIS(1,jtmp)**2
         DIS(2,jtmp) = RV2(2,I,B,X) - RV2(2,KNBr(jtmp),B,X)
         DIS(2,jtmp) = DIS(2,jtmp) - PERlen(2)&
     &                 *NINT(DIS(2,jtmp)/PERlen(2))
         RDYn(jtmp) = RDYn(jtmp) + DIS(2,jtmp)**2
         DIS(3,jtmp) = RV2(3,I,B,X) - RV2(3,KNBr(jtmp),B,X)
         DIS(3,jtmp) = DIS(3,jtmp) - PERlen(3)&
     &                 *NINT(DIS(3,jtmp)/PERlen(3))
         RDYn(jtmp) = RDYn(jtmp) + DIS(3,jtmp)**2
      ENDDO
      NNEigh = 0
      DO jtmp = 1 , jend
!.....
!     determine which pairs are separated by less than rcut
!     and store the needed information about these pairs
!.....
!     if nearest periodic image is out of range, then all images
!     will be
!.....
         IF ( KNBr(jtmp)==NUMnpp1 ) THEN
            rcutsq2 = INDradsq
         ELSE
            rcutsq2 = rcutsq
         ENDIF
         IF ( RDYn(jtmp)<=rcutsq2 ) THEN
            j = KNBr(jtmp)
            NNEigh = NNEigh + 1
            RNEigh(NNEigh) = RDYn(jtmp)
            JNEigh(NNEigh) = j
!--   Vijay
            totneighbors = totneighbors + 1
            NEIghborlist(totneighbors,I) = j
!--   Vijay
 
            DO kcoord = 1 , 3
               DNEigh(kcoord,NNEigh) = DIS(kcoord,jtmp)
            ENDDO
!     check periodic images in z direction
            IF ( .NOT.(ONLyoi) ) THEN
               disz = -SIGN(PERlen(3),DIS(3,jtmp))
               rim = RDYn(jtmp) + 2.*DIS(3,jtmp)*disz + disz**2
 
!     if next nearest image is out of range, subsequent ones will also b
               IF ( rim<=rcutsq2 ) THEN
                  NNEigh = NNEigh + 1
                  RNEigh(NNEigh) = rim
                  JNEigh(NNEigh) = j
!--   Vijay
                  totneighbors = totneighbors + 1
                  NEIghborlist(totneighbors,I) = j
!--   Vijay
 
                  DNEigh(1,NNEigh) = DIS(1,jtmp)
                  DNEigh(2,NNEigh) = DIS(2,jtmp)
                  DNEigh(3,NNEigh) = DIS(3,jtmp) + disz
                  IF ( THReei ) THEN
                     disz = -disz
                     rim = RDYn(jtmp) + 2.*DIS(3,jtmp)*disz + disz**2
                     IF ( rim<=rcutsq2 ) THEN
                        NNEigh = NNEigh + 1
                        RNEigh(NNEigh) = rim
                        JNEigh(NNEigh) = j
!--   Vijay
                        totneighbors = totneighbors + 1
                        NEIghborlist(totneighbors,I) = j
!--   Vijay
                        DNEigh(1,NNEigh) = DIS(1,jtmp)
                        DNEigh(2,NNEigh) = DIS(2,jtmp)
                        DNEigh(3,NNEigh) = DIS(3,jtmp) + disz
                        IF ( FOUri ) THEN
                           disz = -2.*SIGN(PERlen(3),DIS(3,jtmp))
                           rim = RDYn(jtmp) + 2.*DIS(3,jtmp)&
     &                           *disz + disz**2
                           IF ( rim<=rcutsq2 ) THEN
                              NNEigh = NNEigh + 1
                              RNEigh(NNEigh) = rim
                              JNEigh(NNEigh) = j
!--   Vijay
                              totneighbors = totneighbors + 1
                              NEIghborlist(totneighbors,I) = j
!--   Vijay
                              DNEigh(1,NNEigh) = DIS(1,jtmp)
                              DNEigh(2,NNEigh) = DIS(2,jtmp)
                              DNEigh(3,NNEigh) = DIS(3,jtmp) + disz
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!.....
!.....now do diagonal (i=j) term for three or four images of self
!.....both cases produce two images
!.....
      NNEips = NNEigh
      IF ( THReei .AND. I/=NUMnpp1 ) THEN
!.....first image of self
         NNEips = NNEips + 1
         disz = PERlen(3)
         rim = disz**2
!     don't need to check against range here
         RNEigh(NNEips) = rim
         JNEigh(NNEips) = I
!--   Vijay
         totneighbors = totneighbors + 1
         NEIghborlist(totneighbors,I) = I
!--   Vijay
 
         DNEigh(1,NNEips) = 0.0
         DNEigh(2,NNEips) = 0.0
         DNEigh(3,NNEips) = disz
!.....second image of self
         NNEips = NNEips + 1
         disz = -disz
!.....rim = disz**2    done above
         RNEigh(NNEips) = rim
         JNEigh(NNEips) = I
!--   Vijay
         totneighbors = totneighbors + 1
         NEIghborlist(totneighbors,I) = I
!--   Vijay
         DNEigh(1,NNEips) = 0.0
         DNEigh(2,NNEips) = 0.0
         DNEigh(3,NNEips) = disz
      ENDIF
 
      IF ( NNEips>NEIMAX ) THEN
         WRITE (*,99001) NNEips , NEIMAX
99001    FORMAT (' number of neighbors',i5,' exceeds array bound ',i5)
         STOP
      ENDIF
!.....
!.....compute the maximum number of neighbors seen so far
!.....
      NNEimx = MAX0(NNEimx,NNEips)
!.....
!.....compute the timing
!.....
!.....call timeused(icpu,io,isys,imem)
!.....gnetim = (1.e-6)*float(icpu) - timin
!.....gnetmx = amax1(gnetmx,gnetim)
!.....gnetmn = amin1(gnetmn,gnetim)
 
!--   Vijay
!     print*, 'Atom: ', i, ' Neighbors: ', NumNeighbors(i)
      NUMneighbors(I) = totneighbors
!--   Vijay
      RETURN
!.....
!.....end of neighbor finding when using old neighbor list
!.....
!-------------------------------------------------------------
!.....
!.....determine new neighbor list while getting the neighbors
!.....
!.....
!     c-JS            nnindx(i) = nnindx(i-1)
 100  NNIndx(I) = 0
 
      NNEigh = 0
      NNEips = 0
      IF ( Needlist ) THEN
         DO j = 1 , NUMnp
            IF ( I/=j .AND. ISRelaxed(j)/=0 ) THEN
!.....
!     compute the square of the distance to the closest periodic image
!.....
               DIS(1,j) = RV2(1,I,B,X) - RV2(1,j,B,X)
               DIS(1,j) = DIS(1,j) - PERlen(1)*NINT(DIS(1,j)/PERlen(1))
               RDYn(j) = DIS(1,j)**2
               IF ( RDYn(j)<=RCTsqn ) THEN
                  DIS(2,j) = RV2(2,I,B,X) - RV2(2,j,B,X)
                  DIS(2,j) = DIS(2,j) - PERlen(2)&
     &                       *NINT(DIS(2,j)/PERlen(2))
                  RDYn(j) = RDYn(j) + DIS(2,j)**2
                  IF ( RDYn(j)<=RCTsqn ) THEN
                     DIS(3,j) = RV2(3,I,B,X) - RV2(3,j,B,X)
                     DIS(3,j) = DIS(3,j) - PERlen(3)&
     &                          *NINT(DIS(3,j)/PERlen(3))
                     RDYn(j) = RDYn(j) + DIS(3,j)**2
!.....
!.....determine if these particles are within the storage distance
!.....
                     IF ( RDYn(j)<=RCTsqn ) THEN
!.....
!.....store the index of the particle
!.....
!     c--JS
                        NNIndx(I) = NNIndx(I) + 1
!$$$         print *, 'Atom = ', i, nnindx(i), UpdateNeigh(i), newlst, j
!$$$     $        IsRelaxed(j), rdyn(j)
                        NNLst(I,NNIndx(I)) = j
!.....
!     determine which pairs are separated by less than rcut
!     and store the needed information about these pairs
!.....
!     if nearest periodic image is out of range, then all images
!     will be
!.....
                        IF ( RDYn(j)<=rcutsq ) THEN
                           NNEigh = NNEigh + 1
                           RNEigh(NNEigh) = RDYn(j)
                           JNEigh(NNEigh) = j
!--   Vijay
                           totneighbors = totneighbors + 1
                           NEIghborlist(totneighbors,I) = jtmp
!--   Vijay
                           DO kcoord = 1 , 3
                              DNEigh(kcoord,NNEigh) = DIS(kcoord,j)
                           ENDDO
!     check periodic images in z direction
                           IF ( .NOT.(ONLyoi) ) THEN
                              disz = -SIGN(PERlen(3),DIS(3,j))
                              rim = RDYn(j) + 2.*DIS(3,j)*disz + disz**2
!     if next nearest image is out of range, subsequent ones will also b
                              IF ( rim<=rcutsq ) THEN
                                 NNEigh = NNEigh + 1
                                 RNEigh(NNEigh) = rim
                                 JNEigh(NNEigh) = j
!--   Vijay
                                 totneighbors = totneighbors + 1
                                 NEIghborlist(totneighbors,I) = jtmp
!--   Vijay
                                 DNEigh(1,NNEigh) = DIS(1,j)
                                 DNEigh(2,NNEigh) = DIS(2,j)
                                 DNEigh(3,NNEigh) = DIS(3,j) + disz
                                 IF ( THReei ) THEN
                                    disz = -disz
                                    rim = RDYn(j) + 2.*DIS(3,j)&
     &                                 *disz + disz**2
                                    IF ( rim<=rcutsq ) THEN
                                       NNEigh = NNEigh + 1
                                       RNEigh(NNEigh) = rim
                                       JNEigh(NNEigh) = j
!--   Vijay
                                       totneighbors = totneighbors + 1
                                       NEIghborlist(totneighbors,I)&
     &                                    = jtmp
!--   Vijay
                                       DNEigh(1,NNEigh) = DIS(1,j)
                                       DNEigh(2,NNEigh) = DIS(2,j)
                                       DNEigh(3,NNEigh) = DIS(3,j)&
     &                                    + disz
                                       IF ( FOUri ) THEN
                                         disz = -&
     &                                      2.*SIGN(PERlen(3),DIS(3,j))
                                         rim = RDYn(j) + 2.*DIS(3,j)&
     &                                      *disz + disz**2
                                         IF ( rim<=rcutsq ) THEN
                                         NNEigh = NNEigh + 1
                                         RNEigh(NNEigh) = rim
                                         JNEigh(NNEigh) = j
!--   Vijay
                                         totneighbors = totneighbors + 1
                                         NEIghborlist(totneighbors,I)&
     &                                      = jtmp
!--   Vijay
                                         DNEigh(1,NNEigh) = DIS(1,j)
                                         DNEigh(2,NNEigh) = DIS(2,j)
                                         DNEigh(3,NNEigh) = DIS(3,j)&
     &                                      + disz
                                         ENDIF
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
!
!     add indenter
!
         IF ( NUMnpp1>NUMnp .AND. I/=NUMnpp1 ) THEN
            NNIndx(I) = NNIndx(I) + 1
            NNLst(I,NNIndx(I)) = NUMnpp1
            NNEigh = NNEigh + 1
            JNEigh(NNEigh) = NUMnpp1
            RNEigh(NNEigh) = 0.D0
            DO kcoord = 1 , 3
               DNEigh(kcoord,NNEigh) = RV2(kcoord,I,B,X)&
     &                                 - RV2(kcoord,NUMnpp1,B,X)
               RNEigh(NNEigh) = RNEigh(NNEigh) + DNEigh(kcoord,NNEigh)&
     &                          **2
            ENDDO
            IF ( RNEigh(NNEigh)>INDradsq ) NNEigh = NNEigh - 1
         ENDIF
!.....
!.....now do diagonal (i=j) term for three or four images of self
!.....both cases produce two images
!.....
         NNEips = NNEigh
         IF ( THReei .AND. I/=NUMnpp1 ) THEN
!.....first image of self
            NNEips = NNEips + 1
            disz = PERlen(3)
            rim = disz**2
!     don't need to check against range here
            RNEigh(NNEips) = rim
            JNEigh(NNEips) = I
!--   Vijay
            totneighbors = totneighbors + 1
            NEIghborlist(totneighbors,I) = jtmp
!--   Vijay
            DNEigh(1,NNEips) = 0.0
            DNEigh(2,NNEips) = 0.0
            DNEigh(3,NNEips) = disz
!.....second image of self
            NNEips = NNEips + 1
            disz = -disz
!.....rim = disz**2    done above
            RNEigh(NNEips) = rim
            JNEigh(NNEips) = I
!--   Vijay
            totneighbors = totneighbors + 1
            NEIghborlist(totneighbors,I) = jtmp
!--   Vijay
            DNEigh(1,NNEips) = 0.0
            DNEigh(2,NNEips) = 0.0
            DNEigh(3,NNEips) = disz
         ENDIF
         IF ( NNEips>NEIMAX ) THEN
            WRITE (*,99002) NNEips , NEIMAX
99002       FORMAT (' number of neighbors',i5,' exceeds array bound ',&
     &              i5)
            STOP
         ENDIF
!.....
!.....compute the maximum number of neighbors seen so far
!.....
         NNEimx = MAX0(NNEimx,NNEips)
         MXNnei = MAX0(MXNnei,NNIndx(I))
      ENDIF
!     c-JS               mxlstu = max0(mxlstu,nnindx(i))
 
!.....
!.....if this is the last call to gneigh, then set newlst=0 to indicate
!.....that the current list can be used and also increment ngtlst
!.....
      IF ( (I==NUMnpp1) .OR. (NUMnpp1<0 .AND. I==NUMnp) ) THEN
         NEWlst = 0
         NGTlst = NGTlst + 1
         WRITE (*,*) 'UPDATING NEIGHBORS'
      ENDIF
!.....
!.....compute the timing
!.....
!.....call timeused(icpu,io,isys,imem)
!.....gnetim = (1.e-6)*float(icpu) - timin
!.....gnetmx = amax1(gnetmx,gnetim)
!.....gnetmn = amin1(gnetmn,gnetim)
 
!--   Vijay
      NUMneighbors(I) = totneighbors
!     print*, 'Atom: ', i, ' Neighbors: ', NumNeighbors(i)
!!$      if (ISRelaxed(i) >= 1) then 
!!$	write(*,fmt='(A6, I7, 1X, 2I5)', advance = "no") 'Atom =', i, ISRelaxed(i), numneighbors(i) 
!!$	write(*, fmt='(26(1X,I7))') (neighborlist(j,i), j=1,numneighbors(i))
!!$      end if

!--   Vijay
!.....
!.....end of nmeth = 2
!.....
      END SUBROUTINE GNEIGH
!*********************************************************************
!
!     this subroutine computes the displacement of each particle since
!     the last update of the neighbor list.
!     if the maximum displacement is more than 1/2 of dradn, then newlst
!     is set to flag the creation of a new neighbor list by gneigh
!
      SUBROUTINE CHKDIS(B,X)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--CHKDIS485
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION B , drmax , drmax1 , drmax2, tmp , X
      INTEGER i , j , lstold
!*** End of declarations inserted by SPAG
      DIMENSION lstold(NEIMAX)
      DIMENSION B(NDF,*) , X(NXDm,*)
      LOGICAL firsttime
      DATA firsttime/.TRUE./
!     C--Jun Song: User defined distance-test whether atom moves too muc
      DOUBLE PRECISION defdist(NUMnp)
      INTEGER nloop
 
      defdist = 0.0D0
!     C--JS: Reset UpdateNeigh(:)
      UPDateneigh = 0
!
!     branch to appropriate neighbor method
!
!     Ron's Changes
      IF ( firsttime ) THEN
         NEWlst = 1
         firsttime = .FALSE.
      ELSE
         NEWlst = 0
 
!     c--		return
      ENDIF
 
      SELECT CASE (NMEth)
      CASE (2,4)
!
!     nmeth = 2 by default
!
!
!     treat the second call the same as the first in case defects have b
!     added
!
!     also if newlst has been set to 1, do not override even if the part
!     have not moved
!
         IF ( NEWlst/=1 ) THEN
!
!     compare the new positions with the old ones
!
            DO i = 1 , NUMnp
               IF ( ISRelaxed(i)/=0 ) THEN
                  DIS(1,i) = ROLd(1,i) - RV2(1,i,B,X)
                  DIS(1,i) = DIS(1,i) - PERlen(1)&
     &                       *NINT(DIS(1,i)/PERlen(1))
                  DIS(2,i) = ROLd(2,i) - RV2(2,i,B,X)
                  DIS(2,i) = DIS(2,i) - PERlen(2)&
     &                       *NINT(DIS(2,i)/PERlen(2))
                  DIS(3,i) = ROLd(3,i) - RV2(3,i,B,X)
                  DIS(3,i) = DIS(3,i) - PERlen(3)&
     &                       *NINT(DIS(3,i)/PERlen(3))
                  RDYn(i) = DIS(1,i)**2 + DIS(2,i)**2 + DIS(3,i)**2
               ENDIF
            ENDDO
!
!     determine the maximum displacement and compare with  dradn
!
            drmax1 = 0.0
            drmax2 = 0.0
 
            DO i = 1 , NUMnp
               IF ( ISRelaxed(i)/=0 ) THEN
                  tmp = MIN(drmax1,RDYn(i))
                  drmax1 = MAX(drmax1,RDYn(i))
                  drmax2 = MAX(drmax2,tmp)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     C--JS: set UpdateNeigh for i atom if exceed max allowable distance
!     c--this is to find the exact atom that leads to neighbr updating
!     c--Note that here large neighbor list is used based on rctsq
                  defdist(i) = 2.0D0*DSQRT(RDYn(i))
                  IF ( defdist(i)>DRAdn ) THEN
                     UPDateneigh(i) = 1
!     c--also set Update flag for all i's neighbors
                     nloop = NNIndx(i)
                     DO j = 1 , NNIndx(i)
                        UPDateneigh(NNLst(i,j)) = 1
                     ENDDO
                  ENDIF
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               ENDIF
            ENDDO
            drmax = DSQRT(drmax1) + DSQRT(drmax2)
            IF ( drmax<=DRAdn ) THEN
!
!     if here the old neighbor list can be used
!
               NEWlst = 0
               RETURN
            ENDIF
         ENDIF
!-----------------
!
!     if here, a new neighbor list is needed so store the current coordi
!
 
         NEWlst = 1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     C--Jun Song, only update rold for specific atoms
         IF ( SIMstep==1 ) THEN
            DO j = 1 , 3
               DO i = 1 , NUMnp
                  ROLd(j,i) = RV2(j,i,B,X)
               ENDDO
            ENDDO
         ENDIF
 
         IF ( SIMstep>1 ) THEN
            DO i = 1 , NUMnp
               IF ( UPDateneigh(i)==1 ) THEN
                  DO j = 1 , 3
                     ROLd(j,i) = RV2(j,i,B,X)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
 
         RETURN
      CASE (3)
      CASE DEFAULT
         RETURN
      END SELECT
!
!     nmeth = 3
!
      STOP 'method three not here.'
      END SUBROUTINE CHKDIS
 
 
!     Get the neighborlist
 
 
 
 
!.....&Modified: Wed Jan 1810:10:471995 by foiles &
!***************************************************************
      SUBROUTINE CHKPER
      IMPLICIT NONE
!*--CHKPER627
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION CUTOFFR2 , permin , rcut , rcutsq
      INTEGER i , istop , ncalls , nx , nx0 , ny , ny0 , nz , nz0
!*** End of declarations inserted by SPAG
      DATA ncalls/0/
      SAVE nx0 , ny0 , nz0
!
!     check periodicity:
!     allow only one periodic image in the x and y directions
!     allow up to four in the z direction
!
!     x and y:
!     if perlen is lt 2*dsqrt(rcutsq), then have troubles with periodici
!
      rcutsq = CUTOFFR2(1)
      rcut = DSQRT(rcutsq)
      permin = 2.*rcut
      istop = 0
      DO i = 1 , 2
         IF ( PERlen(i)<permin ) THEN
            WRITE (*,99002) i
            WRITE (*,99003) permin , PERlen(i)
            istop = 1
         ENDIF
      ENDDO
      IF ( istop==1 ) STOP
!
!     z direction:
!     four levels
!
!.....initial guess
      ONLyoi = .TRUE.
      TWOi = .FALSE.
      THReei = .FALSE.
      FOUri = .FALSE.
!
      IF ( PERlen(3)<2.*rcut ) THEN
         ONLyoi = .FALSE.
         TWOi = .TRUE.
      ENDIF
      IF ( .NOT.(ONLyoi) ) THEN
         IF ( PERlen(3)<rcut ) THReei = .TRUE.
         IF ( PERlen(3)<2.*rcut/3. ) FOUri = .TRUE.
         IF ( PERlen(3)<0.5*rcut ) THEN
            i = 3
            permin = 0.5*rcut
            WRITE (*,99002) i
            WRITE (*,99003) permin , PERlen(i)
            STOP
         ENDIF
      ENDIF
      nx = 1
      ny = 1
      nz = 1
      IF ( FOUri ) nz = 4
      IF ( THReei .AND. .NOT.FOUri ) nz = 3
      IF ( TWOi .AND. .NOT.THReei ) nz = 2
      IF ( ONLyoi ) nz = 1
!!$      nz = 1
!!$      Threei = .false. 
!!$      TWoi = .false. 
!!$      Fouri = .false. 
!
!     print out the number of images used on the first call only
!
      IF ( ncalls==0 ) THEN
         WRITE (*,99004) nx , ny , nz
         nx0 = nx
         ny0 = ny
         nz0 = nz
      ENDIF
      ncalls = 1
!
!     print out notice if # of periodic images changes
      IF ( nx/=nx0 .OR. ny/=ny0 .OR. nz/=nz0 ) THEN
         nx0 = nx
         ny0 = ny
         nz0 = nz
         WRITE (*,99001)
99001    FORMAT (//' ****** changing # of periodic images ')
         WRITE (*,99004) nx , ny , nz
      ENDIF
99002 FORMAT ('   periodicity is too short in the ',i2,'  direction')
99003 FORMAT ('   permin = ',e15.5,'  periodicity = ',e15.5)
99004 FORMAT (' # of periodic images   ',i10,5x,i10,5x,i10)
      END SUBROUTINE CHKPER
 
      END MODULE MOD_DYNAMO
