!*==mod_grain.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!****************************************************************
!**
!**   MODULE mod_grain : contains definition of model grains and related
!**   routines.
!**
!**
!**   Variable Definitions:
!**   ---------------------
!**
!**   integer variables:
!**   ngrains          - number of grains
!**   type(graintype) variables:
!**   grains(ngrains) -  All data relevant to the grain structure
!**   ncell -  number of atoms in the orthogonal repeating cell
!**   cell(3,ncell) -  coordinates of repeating cell atoms
!**   dcell(3) -  length of periodic cell in xyz directions
!**
!**   matgrain -  material in the grain
!**   xlatvect(3,3) -  lattice vectors defining natural coordinates of
!**   the grain
!**   rotation -  rotation (degrees) of the natural coordinates defining
!**   the grain
!**   rotmat -  rotation matrix elements
!**   numvrts -  number of vertices defining the grain
!**   grainarr -  polygon defining grain
!**   rotgrain -  polygon defining grain rotated to grain's natural
!**   coordinate system
!**   refatom -  reference atom coordinates
!**   fudgevec -  fudge vector
!**   rfudgvec -  fudge vector rotated to grain's natural coordinate
!**   system
!**
!**   Contains Routines:
!**   ------------------
!**
!**   ReadGrainData        - Reads in the number of grains and their
!**   properties
!**   OutputGrainData      - Echoes all grain data, either read in or
!**   computed
!**   ProcessGrains        - Initializes grain data based on that which
!**   is read in
!**   ProcessLatVect       - given 2 lattice vectors, computes third,
!**   and computes rotmat
!**   ProcessGrainGeometry - Computes rotgrain and rfudgvec
!**
!***************************************************************
 
      MODULE MOD_GRAIN
      IMPLICIT NONE
!*--MOD_GRAIN51
 
!     * Type Defintions
      TYPE GRAINTYPE
!     stuff for repeating cell
         INTEGER NCELL
         DOUBLE PRECISION , DIMENSION(:,:) , POINTER :: CELL
         DOUBLE PRECISION , DIMENSION(3) :: DCELL
!     stuff for grain structure and energy
         DOUBLE PRECISION ROTATION , XLATVECT(3,3) , ROTMAT(2) , &
     &                    REFATOM(3)
         INTEGER MATGRAIN , NUMVRTS
         DOUBLE PRECISION , DIMENSION(:) , POINTER :: FUDGEVEC , &
     &                                RFUDGVEC
         DOUBLE PRECISION , DIMENSION(:,:) , POINTER :: GRAINARR , &
     &                                ROTGRAIN
      END TYPE GRAINTYPE
 
!     * Variable Definintions
      INTEGER ngrains
      TYPE (GRAINTYPE)  , DIMENSION(:) , POINTER::grains
 
      CONTAINS
!------------------------------------------------------------------
!     ReadGrainData -- Read in th grain data from .geo file
!
!     Passed Parameters :
!     geomfile  (in):  name of geometry file
!
!     Module Parameters :
!     ngrains         (out) :defined above
!     grains(ngrains) (out)
!
!     Algorithm :
!     obvious
!
!     Notes :
!     may not be general for 3D
!
!     Author :
!     R. Miller (01/17/98)
!
!     Revisions :
!     none
!
!--
      SUBROUTINE READGRAINDATA(Geomfile,Key)
      USE MOD_GLOBAL
      USE MOD_FILE
      USE MOD_MATERIAL
      IMPLICIT NONE
!*--READGRAINDATA102
 
!--   Variables transferred
      CHARACTER(LEN=80) :: Geomfile
      CHARACTER(LEN=4) :: Key
!--   Local variables
      INTEGER i , j , k , iunit
!
!--   Check that the materials are ready
      IF ( NMAterials<1 ) THEN
         WRITE (*,*) &
     &'***ERROR: Material definitions must                come before gr&
     &ain definitions'
         STOP
      ENDIF
!--   Open the specified file
      IF ( Key=='dire' ) THEN
         iunit = 5
      ELSEIF ( FILEEXISTS(Geomfile,.FALSE.) ) THEN
         CALL IOFILE(Geomfile,'formatted  ',iunit,.TRUE.)
      ELSE
         WRITE (*,*) '**ERROR: Grain geometry file not found.**'
         WRITE (*,*) '         Filename:' , Geomfile
         STOP
      ENDIF
 
!--   Readin data
      READ (iunit,*) NGRains
      ALLOCATE (GRAins(NGRains))
      DO i = 1 , NGRains
         ALLOCATE (GRAins(i)%FUDGEVEC(NDM))
         ALLOCATE (GRAins(i)%RFUDGVEC(NDM))
      ENDDO
 
      DO i = 1 , NGRains
         READ (iunit,*) GRAins(i)%MATGRAIN
         IF ( GRAins(i)%MATGRAIN<1 .OR. GRAins(i)%MATGRAIN>NMAterials )&
     &        THEN
            WRITE (*,*) '***ERROR: An undefined material type ' , &
     &                  GRAins(i)%MATGRAIN
            WRITE (*,*) '          has been assigned to grain number ' ,&
     &                  i
            STOP
         ENDIF
         READ (iunit,*) (GRAins(i)%REFATOM(j),j=1,NXDm)
         READ (iunit,*) (GRAins(i)%FUDGEVEC(j),j=1,NDM)
         READ (iunit,*) ((GRAins(i)%XLATVECT(j,k),k=1,3),j=1,2)
         READ (iunit,*) GRAins(i)%ROTATION
         READ (iunit,*) GRAins(i)%NUMVRTS
 
         ALLOCATE (GRAins(i)%GRAINARR(NDM,GRAins(i)%NUMVRTS))
         ALLOCATE (GRAins(i)%ROTGRAIN(NDM,GRAins(i)%NUMVRTS))
 
         READ (iunit,*) ((GRAins(i)%GRAINARR(k,j),k=1,NDM),j=1,GRAins(i)&
     &                  %NUMVRTS)
      ENDDO
!--   Close file
      IF ( iunit/=5 ) CLOSE (iunit)
 
      END SUBROUTINE READGRAINDATA
!------------------------------------------------------------------
! OutputGrainData : Print out all the grain data
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            ngrains         (out) :defined above
!            grains(ngrains) (out)
!
!      Algorithm :
!            obvious
!
!      Notes :
!
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      SUBROUTINE OUTPUTGRAINDATA()
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--OUTPUTGRAINDATA188
!--   Local variables
      INTEGER i , j , k
      WRITE (*,'(a)') 'GRAIN INFORMATION'
      WRITE (6,99001) NGRains
 
!--   Format statements
99001 FORMAT ('    Number of Grains   = ',i3)
      DO i = 1 , NGRains
         WRITE (6,99002) i
99002    FORMAT (&
     &'    ================================================       == Gra&
     &in ',i4)
         WRITE (6,99003) GRAins(i)%MATGRAIN
99003    FORMAT ('        Material ',i4)
         WRITE (6,99004) (GRAins(i)%REFATOM(j),j=1,NXDm)
99004    FORMAT ('        Reference Atom at '/8x,3F10.5)
         WRITE (6,99005) ((GRAins(i)%XLATVECT(j,k),k=1,3),j=1,3)
99005    FORMAT ('        Lattice Vectors',3(/8x,3F10.5))
         WRITE (6,99006) GRAins(i)%ROTATION
99006    FORMAT ('        Rotation (degrees) ',f10.5)
         WRITE (6,99007) GRAins(i)%NUMVRTS
99007    FORMAT ('        Number of vertices = ',i4)
         WRITE (6,99008) ((GRAins(i)%GRAINARR(k,j),k=1,NDM),j=1,&
     &                   GRAins(i)%NUMVRTS)
99008    FORMAT (8x,2E15.5)
         WRITE (6,99009)
99009    FORMAT ('        Cell Structure ')
         WRITE (6,99010) (GRAins(i)%DCELL(j),j=1,3)
99010    FORMAT ('         Cell dimensions',/8x,3F10.5)
         WRITE (6,99011) GRAins(i)%NCELL
99011    FORMAT ('         Cell Bravais Lattice sites :',i4)
         WRITE (6,99012) ((GRAins(i)%CELL(k,j),k=1,3),j=1,GRAins(i)&
     &                   %NCELL)
99012    FORMAT (8x,3E15.5)
      ENDDO
 
      END SUBROUTINE OUTPUTGRAINDATA
!------------------------------------------------------------------
! ProcessGrains -- Given the grain data that was read in, compute other
!                  grain related data.
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            none
!
!      Algorithm :
!            obvious
!
!      Notes :
!          GrainCry call should soon be obsolete
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      SUBROUTINE PROCESSGRAINS
      IMPLICIT NONE
!*--PROCESSGRAINS251
      INTEGER i
!--   Process Lattice Vectors
      CALL PROCESSLATVECT()
 
!--   Check grain geometry
!--   Can we relax the requirement of convex grains?  I think so.
      CALL PROCESSGRAINGEOMETRY()
 
!--   Make Cell Structure and representative crystallite.
      DO i = 1 , NGRains
         CALL GETCELLDATA(GRAins(i))
      ENDDO
      END SUBROUTINE PROCESSGRAINS
 
!------------------------------------------------------------------
! ProcessLatVect : Form three unit lattice vectors
!                  for each grain.
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            grains%xlatvect  (in/out)
!            grains%rotmat    (out)
!
!      Algorithm :
!        For each grain do
!           1)  Find third lattice vector by cross product
!           2)  Normalise the three vectors
!           3)  Check if determinent of resulting matrix is unity
!           4)  Make cosine and sine of rotation angle
!        end do
!
!      Notes :
!            the only variables in the grain data that are changed
!            are xlatvect and rotmat
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      SUBROUTINE PROCESSLATVECT()
      IMPLICIT NONE
!*--PROCESSLATVECT298
 
!--   Local Variables
      INTEGER i , j , k
      DOUBLE PRECISION detq , qmag , q(3,3)
      DOUBLE PRECISION , PARAMETER :: EPS = 1.D-09 , &
     &                                PI = 3.1415926535898
 
 
!--   Loop over all grains
      DO i = 1 , NGRains
 
!Load in lattice vectors to local arrays
         DO j = 1 , 2
            DO k = 1 , 3
               q(j,k) = GRAins(i)%XLATVECT(j,k)
            ENDDO
         ENDDO
 
!Find third vector by cross product
         q(3,1) = q(1,2)*q(2,3) - q(1,3)*q(2,2)
         q(3,2) = q(1,3)*q(2,1) - q(1,1)*q(2,3)
         q(3,3) = q(1,1)*q(2,2) - q(1,2)*q(2,1)
 
!Normalise the  vectors
         DO j = 1 , 3
            qmag = DSQRT(1.D0/(q(j,1)*q(j,1)+q(j,2)*q(j,2)+q(j,3)*q(j,3)&
     &             ))
            DO k = 1 , 3
               q(j,k) = q(j,k)*qmag
            ENDDO
         ENDDO
 
 
!Check if determinant is unity
         detq = q(1,1)*(q(2,2)*q(3,3)-q(2,3)*q(3,2)) - q(1,2)&
     &          *(q(2,1)*q(3,3)-q(2,3)*q(3,1)) + q(1,3)&
     &          *(q(2,1)*q(3,2)-q(2,2)*q(3,1))
         IF ( DABS(detq-1.D0)>EPS ) THEN
            PRINT *
            PRINT * , '*** ERROR: Lattice Vector Error in grain ' , i
            PRINT * , 'Determinant = ' , detq
            PRINT *
            STOP
         ENDIF
 
!Put back the lattice vector into grainstrc
         DO j = 1 , 3
            DO k = 1 , 3
               GRAins(i)%XLATVECT(j,k) = q(j,k)
            ENDDO
         ENDDO
 
!Make cosine and sine of the rotation angle
         GRAins(i)%ROTMAT(1) = COS(PI*GRAins(i)%ROTATION/180.0)
         GRAins(i)%ROTMAT(2) = SIN(PI*GRAins(i)%ROTATION/180.0)
      ENDDO
      END SUBROUTINE PROCESSLATVECT
!------------------------------------------------------------------
! ProcessGrainGeometry -- perform rotations on grain data
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            ngrains        (in)
!            grains%refatom (in)
!            grains%numvrts (in)
!            grains%rotmat  (in)
!            grains%grainarr(in)
!            grains%rotgrain(out)
!            grains%fudgevec(in)
!            grains%rfudgvec(out)
!
!      Algorithm :
!       Algorithm :
!          1) check if reference atoms are in the interior of grain
!          1) compute rotgrain and rfudgvec
!
!      Notes :
!            notes
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      SUBROUTINE PROCESSGRAINGEOMETRY
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--PROCESSGRAINGEOMETRY390
!--   Local variables
      INTEGER i , j , k , nvt
      LOGICAL POINTINGRAIN
      DOUBLE PRECISION u(NDM)
 
      DO i = 1 , NGRains
         IF ( .NOT.POINTINGRAIN(GRAins(i)%REFATOM,i) ) THEN
            WRITE (*,*) '***ERROR: The reference atom of grain' , i
            WRITE (*,*) '          is outside of the grain polygon.'
            STOP
         ENDIF
         nvt = GRAins(i)%NUMVRTS
         DO j = 1 , nvt
            DO k = 1 , NDM
               u(k) = GRAins(i)%GRAINARR(k,j) - GRAins(i)%REFATOM(k)
            ENDDO
            GRAins(i)%ROTGRAIN(1,j) = u(1)*GRAins(i)%ROTMAT(1) + u(2)&
     &                                *GRAins(i)%ROTMAT(2)
            GRAins(i)%ROTGRAIN(2,j) = -u(1)*GRAins(i)%ROTMAT(2) + u(2)&
     &                                *GRAins(i)%ROTMAT(1)
         ENDDO
         GRAins(i)%RFUDGVEC(1) = GRAins(i)%FUDGEVEC(1)*GRAins(i)&
     &                           %ROTMAT(1) + GRAins(i)%FUDGEVEC(2)&
     &                           *GRAins(i)%ROTMAT(2)
         GRAins(i)%RFUDGVEC(2) = -GRAins(i)%FUDGEVEC(1)*GRAins(i)&
     &                           %ROTMAT(2) + GRAins(i)%FUDGEVEC(2)&
     &                           *GRAins(i)%ROTMAT(1)
      ENDDO
      END SUBROUTINE PROCESSGRAINGEOMETRY
!***********************************************************************
      SUBROUTINE GETCELLDATA(Grain)
      USE MOD_GLOBAL
      USE MOD_MATERIAL
      USE MOD_CLUSTER
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--GETCELLDATA427
 
!** Transferred Variables **!
      TYPE (GRAINTYPE) Grain
 
!** Parameter Declarations **!
      DOUBLE PRECISION , PARAMETER :: TOL = 1.D-6
 
!** Local Variables **!
      DOUBLE PRECISION tmp(3) , brot(3,3) , rclust , rtol , cx , cy , cz
      DOUBLE PRECISION , POINTER :: cell2(:,:)
      INTEGER icry , i , j , k , maxcell , nlast , nsort
      LOGICAL xzero , yzero , zzero , break
      TYPE (BRAVAISMAT) matl
      TYPE (CLUSTER) clust
!     c
!     C     Analyze rotated crystal structure
!     c
!     c Find the repeat cell of bravais sites in the natural coord
!     system of the grain
!     c by building a temporary cluster of bravais sites and analysing
!     it.  Note that
!     c if the assumed size of the cluster is too small, it increases it
!     and tries
!     c again.
!     c
!     ctry something simple: find atoms closest to origin along
!     c each axis and take this as the unit cell.
 
! Rotate BL vectors into slip c.s. (note sexy f90 function call)
      brot = MATMUL(Grain%XLATVECT,MATerial(Grain%MATGRAIN)%BVEC)
      rclust = SQRT(DOT_PRODUCT(brot(1:3,1),brot(1:3,1)))
      rtol = TOL*rclust
      rclust = rclust*5.D0
      maxcell = 100
      DO
!     c
         ALLOCATE (cell2(3,maxcell))
         matl%BVEC = brot
         matl%NBASIS = 1
         ALLOCATE (matl%BASIS(3,1))
         matl%BASIS(1:3,1) = 0.0D0
         ALLOCATE (matl%ISPEC(1))
         matl%ISPEC(1) = 1
         matl%VOLUME = MATerial(Grain%MATGRAIN)%VOLUME
 
         CALL BUILDCLUSTER(matl,rclust,clust)
 
         DEALLOCATE (matl%BASIS)
         DEALLOCATE (matl%ISPEC)
         cx = 1.D30
         cy = 1.D30
         cz = 1.D30
         DO i = 1 , clust%NATOMS
            IF ( (clust%X(1,i)>=-rtol) .AND. (clust%X(2,i)>=-rtol) .AND.&
     &           (clust%X(3,i)>=-rtol) ) THEN
               xzero = DABS(clust%X(1,i))<rtol
               yzero = DABS(clust%X(2,i))<rtol
               zzero = DABS(clust%X(3,i))<rtol
               IF ( xzero .AND. yzero .AND. (.NOT.zzero) ) THEN
                  IF ( clust%X(3,i)<cz ) cz = clust%X(3,i)
               ELSEIF ( xzero .AND. zzero .AND. (.NOT.yzero) ) THEN
                  IF ( clust%X(2,i)<cy ) cy = clust%X(2,i)
               ELSEIF ( yzero .AND. zzero .AND. (.NOT.xzero) ) THEN
                  IF ( clust%X(1,i)<cx ) cx = clust%X(1,i)
               ENDIF
            ENDIF
         ENDDO
         IF ( cx>1.D29 ) THEN
            rclust = 2*rclust
            PRINT * , &
     &           '***WARNING: Repeating pattern not located in x-dir''n'
            PRINT * , '          Increased computed cluster size:' , &
     &            rclust
            CYCLE
         ELSEIF ( cy>1.D29 ) THEN
            rclust = 2*rclust
            PRINT * , &
     &           '***WARNING: Repeating pattern not located in y-dir''n'
            PRINT * , '          Increased computed cluster size:' , &
     &            rclust
            CYCLE
         ELSEIF ( cz>1.D29 ) THEN
            rclust = 2*rclust
            PRINT * , &
     &           '***WARNING: Repeating pattern not located in z-dir''n'
            PRINT * , '          Increased computed cluster size:' , &
     &            rclust
            CYCLE
         ENDIF
         Grain%NCELL = 0
         DO i = 1 , clust%NATOMS
            xzero = (clust%X(1,i)<cx-rtol) .AND. (clust%X(1,i)>-rtol)
            yzero = (clust%X(2,i)<cy-rtol) .AND. (clust%X(2,i)>-rtol)
            zzero = (clust%X(3,i)<cz-rtol) .AND. (clust%X(3,i)>-rtol)
            IF ( xzero .AND. yzero .AND. zzero ) THEN
               Grain%NCELL = Grain%NCELL + 1
               IF ( Grain%NCELL>maxcell ) THEN
                  maxcell = maxcell*2
                  PRINT * , &
     &                 '***WARNING: Insufficient storage in cell matrix'
                  PRINT * , '            increased maxcell to:' , &
     &                  maxcell
                  GOTO 100
               ENDIF
               DO j = 1 , 3
                  cell2(j,Grain%NCELL) = clust%X(j,i)
               ENDDO
            ENDIF
         ENDDO
         CALL DEALLOCATE_CLUSTER(clust)
         ALLOCATE (Grain%CELL(3,Grain%NCELL))
         Grain%CELL(1:3,1:Grain%NCELL) = cell2(1:3,1:Grain%NCELL)
!     c
!     c sort into ascending y-layers.  In each layer sort into ascending
!     c x coord.  If more than one with same x coord, sort into
!     ascending z
!     c coord.
!     c
!     c  First sort, y coord as key:
!     c
         CALL SORT(Grain%CELL,3,Grain%NCELL,2,Grain%NCELL)
!     c
!     c  Second sort, within each y-layer sort by x coord
!     c
         nlast = 0
         DO i = 1 , Grain%NCELL
            break = (i==Grain%NCELL)
            IF ( .NOT.break ) break = (ABS(Grain%CELL(2,i+1)-Grain%CELL(&
     &                                2,i))>rtol)
            IF ( break ) THEN
               nsort = i - nlast
               DO j = 1 , nsort
                  DO k = 1 , 3
                     cell2(k,j) = Grain%CELL(k,nlast+j)
                  ENDDO
               ENDDO
               CALL SORT(cell2,3,maxcell,1,nsort)
               DO j = 1 , nsort
                  DO k = 1 , 3
                     Grain%CELL(k,nlast+j) = cell2(k,j)
                  ENDDO
               ENDDO
               nlast = i
            ENDIF
         ENDDO
!     c
!     c third sort, within each line of atoms with the same x and y,
!     c sort by z-coord.
!     c
         nlast = 0
         DO i = 1 , Grain%NCELL
            break = (i==Grain%NCELL)
            IF ( .NOT.break ) break = (ABS(Grain%CELL(2,i+1)-Grain%CELL(&
     &                                2,i))>rtol) .OR. &
     &                                (ABS(Grain%CELL(1,i+1)&
     &                                -Grain%CELL(1,i))>rtol)
            IF ( break ) THEN
               nsort = i - nlast
               DO j = 1 , nsort
                  DO k = 1 , 3
                     cell2(k,j) = Grain%CELL(k,nlast+j)
                  ENDDO
               ENDDO
               CALL SORT(cell2,3,maxcell,3,nsort)
               DO j = 1 , nsort
                  DO k = 1 , 3
                     Grain%CELL(k,nlast+j) = cell2(k,j)
                  ENDDO
               ENDDO
               nlast = i
            ENDIF
         ENDDO
         DEALLOCATE (cell2)
!     c
!     c at this point, cell is sorted.
!     c
!     c     Compute projected bravais site area in the global xy plane
         Grain%DCELL(1) = cx
         Grain%DCELL(2) = cy
         Grain%DCELL(3) = cz
         Z_Length = NUMperiodz*cz
         PERlen(3) = NUMperiodz*cz
         PERub(3) = NUMperiodz*cz
         PERlb(3) = 0.D0
         EXIT
 100  ENDDO
      END SUBROUTINE GETCELLDATA
 
      SUBROUTINE SORT(Ra,Mra,Nra,M,N)
!     C     Based on a Numerical Recipes routine.
      IMPLICIT NONE
!*--SORT619
      INTEGER , INTENT(IN) :: Nra , Mra , N , M
      DOUBLE PRECISION , INTENT(INOUT) :: Ra(Mra,Nra)
 
      INTEGER , PARAMETER :: MAXM = 10
      DOUBLE PRECISION rra(MAXM)
      INTEGER ir , i , j , k , l
      IF ( N<=1 ) RETURN
      IF ( Mra>MAXM ) THEN
         PRINT * , '***ERROR: MRA too large in SORT'
         STOP
      ENDIF
      l = N/2 + 1
      ir = N
 100  IF ( l>1 ) THEN
         l = l - 1
         DO k = 1 , Mra
            rra(k) = Ra(k,l)
         ENDDO
      ELSE
         DO k = 1 , Mra
            rra(k) = Ra(k,ir)
         ENDDO
         DO k = 1 , Mra
            Ra(k,ir) = Ra(k,1)
         ENDDO
         ir = ir - 1
         IF ( ir==1 ) THEN
            DO k = 1 , Mra
               Ra(k,1) = rra(k)
            ENDDO
            RETURN
         ENDIF
      ENDIF
      i = l
      j = l + l
      DO
         IF ( j<=ir ) THEN
            IF ( j<ir ) THEN
               IF ( Ra(M,j)<Ra(M,j+1) ) j = j + 1
            ENDIF
            IF ( rra(M)<Ra(M,j) ) THEN
               DO k = 1 , Mra
                  Ra(k,i) = Ra(k,j)
               ENDDO
               i = j
               j = j + j
            ELSE
               j = ir + 1
            ENDIF
            CYCLE
         ENDIF
         DO k = 1 , Mra
            Ra(k,i) = rra(k)
         ENDDO
         GOTO 100
      ENDDO
      END SUBROUTINE SORT
 
 
      END MODULE MOD_GRAIN
