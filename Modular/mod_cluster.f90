!*==mod_cluster.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!******************************************************************
!**
!**   MODULE mod_cluster: routines pertaining to building, storing, and
!**   manipulating clusters.
!**
!**   Variable TYPES Defined:
!**   ---------------------
!**
!**   type(cluster) variables:
!**   %natoms          - number of atoms in the cluster
!**   %x(3,natoms)     - coordinates of the atoms
!**   %spec(natoms)    - species of each atom
!**   %lattice(natoms) - lattice from which the atom is derived (for
!**   complex lattices)
!**
!**   Contains Routines:
!**   ------------------
!**
!**   BuildCluster         - Builds a cluster of given radius from a
!**   given material
!**   deallocate_cluster   - Deallocates a type(cluster) variable
!**   SortBravais          - pre-processing for buildcluster
!**
!********************************************************************
 
      MODULE MOD_CLUSTER
      IMPLICIT NONE
!*--MOD_CLUSTER29
 
!     * Type Defintions
      TYPE CLUSTER
         INTEGER NATOMS
         DOUBLE PRECISION , DIMENSION(:,:) , POINTER :: X , RLAT
         INTEGER , DIMENSION(:) , POINTER :: SPEC , LATTICE
      END TYPE CLUSTER
 
      CONTAINS
 
!------------------------------------------------------------------
! BuildCluster -- Build a cluster of atoms
!
!      Passed Parameters :
!            mat  (in) : type(bravaismat), contains the material
!                        from which to build the cluster
!            rcut (in) : radius of cluste
!            atoms(out): type(cluster), contains the built cluster
!
!      Algorithm :
!         Span the directions of the three Bravais vectors,
!         but try to be clever about minimizing the range of
!         each loop (important for non-orthogonal bravais lattice
!         vectors).
!
!         Add to the cluster all atoms associated with bravais
!         sites less than rcut from the origin, making the origin
!         site atoms first in the list.
!
!      Notes :
!         Needs to be rigourously tested for speed and optimized
!         if possible - this will be called for every call to
!     "cauchyborn"
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      SUBROUTINE BUILDCLUSTER(Mat,Rcut,Atoms)
      USE MOD_MATERIAL
      IMPLICIT NONE
!*--BUILDCLUSTER74
 
!     variables that are passed
      TYPE (BRAVAISMAT) Mat
      DOUBLE PRECISION Rcut
      TYPE (CLUSTER) Atoms
 
!--   Local Variables
!
      TYPE (CLUSTER) at
      INTEGER ii(3) , i1 , i2 , i3 , is1 , is2 , is3 , if1 , if2 , if3 ,&
     &        ic1 , ic2 , ic3 , nguess , i , j , ia , ib , ib2
      DOUBLE PRECISION , PARAMETER :: PI43 = 4.18879D0
      DOUBLE PRECISION latlen(3) , org(2) , per3(2) , bperp , rlen2 , &
     &                 r(3) , rb(3) , per2 , cperp , rcut2 , vec(3) , &
     &                 rtol
      LOGICAL notzero1 , notzero2 , notzero3
!
 
      rcut2 = Rcut*Rcut
!
!     estimate cluster size using (volume of sphere)/(volume per atom),
!     with
!     a safty factor of 2
!
      nguess = 2*NINT(Mat%NBASIS*(PI43*Rcut*rcut2)/(Mat%VOLUME))
      ALLOCATE (at%X(3,nguess))
      ALLOCATE (at%RLAT(Mat%NBASIS,nguess))
      ALLOCATE (at%SPEC(nguess))
      ALLOCATE (at%LATTICE(nguess))
!     C@
      at%X = 0.
      at%RLAT = 0.
      at%SPEC = 0
      at%LATTICE = 0
      DO
!
!     compute lengths of bravais lattice vectors
!
         DO j = 1 , 3
            latlen(j) = DOT_PRODUCT(Mat%BVEC(1:3,j),Mat%BVEC(1:3,j))
            latlen(j) = SQRT(latlen(j))
         ENDDO
!
!     to optimize speed, sort (using sorting map ii) the bvec.
!     Call the bvec vectors "a" "b" and "c".  "a" is chosen to be
!     the shortest of the three vectors.  Of the remaining 2, "b" is
!     the one that is "least perpendicular" to "a".  Finally, "c" is the
!     remaining vector.  map these using ii(1)=a, ii(2)=b, ii(3)=c.
!
!     bperp is the length of b perp. to a.  Similarly cperp.
!     "per3" is the projection of c in the plane defined by a and b
!     "per2" is the projection of b along a.
!
         CALL SORTBRAVAIS(Mat%BVEC,latlen,ii,bperp,cperp,per3,per2)
!
!     Put in the center site atoms (bravais lattice (0,0,0))
!     now, so that they appear first in list.
!
         at%NATOMS = 0
         DO ib = 1 , Mat%NBASIS
            DO j = 1 , 3
               rb(j) = Mat%BASIS(ii(1),ib)*Mat%BVEC(j,ii(1))&
     &                 + Mat%BASIS(ii(2),ib)*Mat%BVEC(j,ii(2))&
     &                 + Mat%BASIS(ii(3),ib)*Mat%BVEC(j,ii(3))
            ENDDO
            at%NATOMS = at%NATOMS + 1
            at%X(1:3,at%NATOMS) = rb
            at%SPEC(at%NATOMS) = Mat%ISPEC(ib)
            at%LATTICE(at%NATOMS) = ib
         ENDDO
!     update the rlat
         DO ia = 1 , at%NATOMS
            DO ib = 1 , at%NATOMS
               IF ( ia==ib ) THEN
                  at%RLAT(ib,ia) = 0
               ELSE
                  vec = at%X(1:3,ia) - at%X(1:3,ib)
                  at%RLAT(ib,ia) = SQRT(DOT_PRODUCT(vec,vec))
               ENDIF
            ENDDO
         ENDDO
 
!
!     decide how big to make the loop over the "c" direction
!
         if3 = INT(Rcut/cperp)
         is3 = INT(-Rcut/cperp)
         IF ( if3>is3 ) THEN
            ic3 = 1
         ELSE
            ic3 = -1
         ENDIF
!
!     loop in the c direction
!
         DO i3 = is3 , if3 , ic3
            notzero3 = i3/=0
            org(2) = i3*per3(2)
!
!     decide how big to make the loop over the "b" direction
!
            if2 = INT((Rcut-org(2))/bperp)
            is2 = INT((-Rcut-org(2))/bperp)
            IF ( if2>is2 ) THEN
               ic2 = 1
            ELSE
               ic2 = -1
            ENDIF
!
!     loop in the b direction
!
            DO i2 = is2 , if2 , ic2
               notzero2 = i2/=0
               org(1) = i3*per3(1) + i2*per2
!
!     decide how big to make the loop over the "a" direction
!
               if1 = INT((Rcut-org(1))/latlen(ii(1)))
               is1 = INT((-Rcut-org(1))/latlen(ii(1)))
               IF ( if1>is1 ) THEN
                  ic1 = 1
               ELSE
                  ic1 = -1
               ENDIF
!
!     loop in the a direction
!
               DO i1 = is1 , if1 , ic1
                  notzero1 = i1/=0
                  IF ( notzero1 .OR. notzero2 .OR. notzero3 ) THEN
                     rlen2 = 0.
                     DO j = 1 , 3
                        r(j) = i1*Mat%BVEC(j,ii(1))&
     &                         + i2*Mat%BVEC(j,ii(2))&
     &                         + i3*Mat%BVEC(j,ii(3))
                        rlen2 = rlen2 + r(j)*r(j)
                     ENDDO
                     IF ( rlen2<=rcut2 ) THEN
!
!     add all atoms associated with this BL site
!
                        DO ib = 1 , Mat%NBASIS
                           DO j = 1 , 3
                              rb(j) = r(j) + Mat%BASIS(ii(1),ib)&
     &                                *Mat%BVEC(j,ii(1))&
     &                                + Mat%BASIS(ii(2),ib)&
     &                                *Mat%BVEC(j,ii(2))&
     &                                + Mat%BASIS(ii(3),ib)&
     &                                *Mat%BVEC(j,ii(3))
                           ENDDO
!     add an atom
                           at%NATOMS = at%NATOMS + 1
!
!     make sure the assumed array size is still okay.  If not, start
!     again.
!
                           IF ( at%NATOMS>nguess ) THEN
                              WRITE (*,*) &
     &                     '***WARNING: BuildCluster had to reallocate!'
                              nguess = nguess*2
                              ALLOCATE (at%X(3,nguess))
                              ALLOCATE (at%RLAT(Mat%NBASIS,nguess))
                              ALLOCATE (at%SPEC(nguess))
                              ALLOCATE (at%LATTICE(nguess))
!C@
                              at%X = 0.
                              at%RLAT = 0.
                              at%SPEC = 0
                              at%LATTICE = 0
 
                              GOTO 100
                           ENDIF
                           at%X(1:3,at%NATOMS) = rb
                           at%SPEC(at%NATOMS) = Mat%ISPEC(ib)
                           at%LATTICE(at%NATOMS) = ib
! compute rlat
                           DO ib2 = 1 , Mat%NBASIS
                              vec = at%X(1:3,at%NATOMS) - at%X(1:3,ib2)
                              at%RLAT(ib2,at%NATOMS)&
     &                           = SQRT(DOT_PRODUCT(vec,vec))
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
!
! Now allocate the array to be returned, and copy the data from
!  the temporary storage
!
         Atoms%NATOMS = at%NATOMS
         ALLOCATE (Atoms%X(3,at%NATOMS))
!C@
         Atoms%X = 0.
         Atoms%X(1:3,1:at%NATOMS) = at%X(1:3,1:at%NATOMS)
         ALLOCATE (Atoms%RLAT(Mat%NBASIS,at%NATOMS))
         Atoms%RLAT = 0.
         Atoms%RLAT(1:Mat%NBASIS,1:at%NATOMS)&
     &      = at%RLAT(1:Mat%NBASIS,1:at%NATOMS)
         ALLOCATE (Atoms%SPEC(at%NATOMS))
         Atoms%SPEC = 0
         Atoms%SPEC(1:at%NATOMS) = at%SPEC(1:at%NATOMS)
         ALLOCATE (Atoms%LATTICE(at%NATOMS))
         Atoms%LATTICE = 0
         Atoms%LATTICE(1:at%NATOMS) = at%LATTICE(1:at%NATOMS)
!
! deallocate the temporary storage
!
         CALL DEALLOCATE_CLUSTER(at)
         EXIT
 100  ENDDO
      END SUBROUTINE BUILDCLUSTER
 
   !------------------------------------------------------------------
   ! deallocate_cluster -- deallocate a type(cluster) variable
   !
   !      Passed Parameters :
   !            c  (in) : type(cluster) to be deallocated
   !
   !      Notes :
   !         It is important to deallocate all temporary clusters after
   !         use!
   !
   !      Author :
   !            R. Miller (01/17/98)
   !
   !      Revisions :
   !              none
   !
   !--
 
      SUBROUTINE DEALLOCATE_CLUSTER(C)
      IMPLICIT NONE
!*--DEALLOCATE_CLUSTER309
      TYPE (CLUSTER) C
      DEALLOCATE (C%X)
      DEALLOCATE (C%RLAT)
      DEALLOCATE (C%SPEC)
      DEALLOCATE (C%LATTICE)
      END SUBROUTINE DEALLOCATE_CLUSTER
 
   !------------------------------------------------------------------
   ! SortBravais -- Pre-processing for buildcluster
   !
   !      Passed Parameters :
   !            lat   (in): bravais lattice vectors
   !            latlen(in): lengths of lat
   !            ii   (out): sort map for lat a=1,b=2,c=3
   !            bperp(out): component of b perp. to a
   !            bperp(out): component of c perp. to a
   !            per3 (out): projection of c on a-b plane
   !            per2 (out): projection of b on a
   !
   !      Algorithm :
   !         make the shortest vector "a", make "b" the one
   !         that is least perpendicular to "a" and "c" is the leftover
   !
   !      Notes :
   !         Needs to be rigourously tested for speed and optimized
   !         if possible - this will be called for every call to "cauchy
   !
   !      Author :
   !            R. Miller (01/17/98)
   !
   !      Revisions :
   !              none
   !
   !--
 
      SUBROUTINE SORTBRAVAIS(Lat,Latlen,Ii,Bperp,Cperp,Per3,Per2)
      IMPLICIT NONE
!*--SORTBRAVAIS347
      INTEGER Ii(3) , j , i , itemp
      DOUBLE PRECISION Lat(3,3) , Bperp , Per3(2) , Latlen(3) , t3(3) , &
     &                 Per2 , Cperp , t2(3) , bperp2 , bperp3 , temp2 , &
     &                 temp3 , per22 , per23
!
! make the shortest vector "a"
!
      Ii(1) = 1
      DO i = 2 , 3
         IF ( Latlen(i)<Latlen(Ii(1)) ) Ii(1) = i
      ENDDO
      Ii(2) = MOD(Ii(1),3) + 1
      Ii(3) = MOD(Ii(2),3) + 1
!
! bperp is the length of the component of b perp. to a.
! bperp2 is using vector 2 as b,
! bperp3 is using vector 3 as b.
!
! choose the shorter of bperp2,bperp3 to be bperp
!
      per22 = DOT_PRODUCT(Lat(1:3,Ii(1)),Lat(1:3,Ii(2)))/Latlen(Ii(1))
      temp2 = per22/Latlen(Ii(1))
      per23 = DOT_PRODUCT(Lat(1:3,Ii(1)),Lat(1:3,Ii(3)))/Latlen(Ii(1))
      temp3 = per23/Latlen(Ii(1))
      bperp2 = 0.
      bperp3 = 0.
      DO j = 1 , 3
         t2(j) = Lat(j,Ii(2)) - temp2*Lat(j,Ii(1))
         bperp2 = bperp2 + t2(j)**2
         t3(j) = Lat(j,Ii(3)) - temp3*Lat(j,Ii(1))
         bperp3 = bperp3 + t3(j)**2
      ENDDO
      IF ( bperp2<bperp3 ) THEN
         Bperp = DSQRT(bperp2)
         Per2 = per22
         Per3(1) = DOT_PRODUCT(Lat(1:3,Ii(3)),Lat(1:3,Ii(1)))&
     &             /Latlen(Ii(1))
         Per3(2) = DOT_PRODUCT(Lat(1:3,Ii(3)),t2)/Bperp
      ELSE
         Bperp = DSQRT(bperp3)
         Per2 = per23
         itemp = Ii(2)
         Ii(2) = Ii(3)
         Ii(3) = itemp
         Per3(1) = DOT_PRODUCT(Lat(1:3,Ii(3)),Lat(1:3,Ii(1)))&
     &             /Latlen(Ii(1))
         Per3(2) = DOT_PRODUCT(Lat(1:3,Ii(3)),t3)/Bperp
      ENDIF
!
      CALL CROSS_PRODUCT(Lat(1:3,Ii(1)),Lat(1:3,Ii(2)),t3)
      Cperp = DOT_PRODUCT(t3,Lat(1:3,Ii(3)))/SQRT(DOT_PRODUCT(t3,t3))
      END SUBROUTINE SORTBRAVAIS
      END MODULE MOD_CLUSTER
