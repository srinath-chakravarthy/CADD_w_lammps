!*==nearestbsite.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**------------------------------------------------------------------
!** NearestBSite : get nearest Bravais site to the point in the model
!**
!**    Parmeters :
!**            u   (in) : Coordinates of the point
!**        nnear   (in) : Number of nearest sites wanted
!**   CheckModel   (in) : Flag, true if point must be in model
!**           xa  (out) : coordinates of the B-site
!**       igrain  (out) : Grain in which these sites are
!**
!**    Algorithm :
!**          for each grain
!**              find that B-site that is closest
!**              is it the smallest
!**          end do
!**
!--
      SUBROUTINE NEARESTBSITE(U,Nnear,Checkmodel,Xa,Igrain)
 
      USE MOD_GRAIN
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--NEARESTBSITE24
 
!--Variables Transferred
      DOUBLE PRECISION U(NDM) , Xa(3,1)
!90 ==> This variable should be removed and dynamic allocation used
      INTEGER NNEAMAX
      PARAMETER (NNEAMAX=9)
      INTEGER Igrain(1) , Nnear
      LOGICAL Checkmodel
 
!--Local Variables
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
 
      INTEGER i , j , k , k1 , kt , ispt1
      DOUBLE PRECISION x(3,NNEAMAX) , rmin(NNEAMAX) , r(NNEAMAX) , &
     &                 rtmp , xtmp(3)
 
      LOGICAL change
!
 
      IF ( Nnear>NNEAMAX ) THEN
         PRINT * , 'NNEAMAX smaller than requested number neighs.'
         PRINT * , 'Increase NNEAMAX in subroutine NearestBSite'
         PRINT * , 'NNEAMAX =' , NNEAMAX
         STOP
      ENDIF
 
      !Start with the first grain
      IF ( Nnear==1 ) THEN
         CALL NEARBSITEINGRAIN1(U,1,Checkmodel,x(1,1),r(1))
      ELSE
         CALL NEARBSITEINGRAIN(U,Nnear,1,Checkmodel,x(1,1),r(1))
      ENDIF
      DO k = 1 , Nnear
         Igrain(k) = 1
         rmin(k) = r(k)
         DO j = 1 , 3
            Xa(j,k) = x(j,k)
         ENDDO
      ENDDO
 
      !now for rest of the grains
      DO i = 2 , NGRains
         change = .FALSE.
         IF ( Nnear==1 ) THEN
            CALL NEARBSITEINGRAIN1(U,i,Checkmodel,x(1,1),r(1))
         ELSE
            CALL NEARBSITEINGRAIN(U,Nnear,i,Checkmodel,x(1,1),r(1))
         ENDIF
         DO k = 1 , Nnear
            rtmp = rmin(k)
            DO k1 = 1 , Nnear
               IF ( r(k1)<rtmp ) THEN
                  change = .TRUE.
                  kt = k1
                  rtmp = r(k1)
               ENDIF
            ENDDO
 
            IF ( change ) THEN
               DO j = 1 , 3
                  xtmp(j) = Xa(j,k)   !flip atoms
                  Xa(j,k) = x(j,kt)
                  x(j,kt) = xtmp(j)
               ENDDO
               r(kt) = rmin(k)        !flip distances
               rmin(k) = rtmp
               Igrain(k) = i          !set grain
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE NEARESTBSITE
!*==nearbsiteingrain.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
 
!**------------------------------------------------------------------
!** NearBSiteInGrain : Get Bravais site nearest to given point in the gr
!**
!**      Parameters :
!**                u  (in) : coordinates of the point
!**            nnear  (in) : number of neartest B-sites requested
!**           igrain  (in) : the grain that is required
!**       CheckModel  (in) : Flag true if point must lie in the model
!**               xa (out) : coordinates of the B-sites
!**               rd (out) : distance square of near B-site
!**
!**      Algorithm :
!**           check if point is in the grain
!**           transform point to grain coordinates.
!**           ff (exterior point) then
!**              Find the nearest interior point
!**           end if
!**           find Cell and its eight neighbours
!**           for all B-sites in the cell,
!**              if (not interior) discard
!**              find the smallest
!**           end do
!**
!**      Notes :
!**           Does not change polygrain.
!**           This is a very messy sub.
!**
!--
      SUBROUTINE NEARBSITEINGRAIN(U,Nnear,Igrain,Checkmodel,Xa,Rd)
 
      USE MOD_GRAIN
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--NEARBSITEINGRAIN132
 
      INTEGER NNEAMAX
      PARAMETER (NNEAMAX=9)
 
!--Variables Transferred
      DOUBLE PRECISION U(NDM) , Xa(3,NNEAMAX) , Rd(NNEAMAX)
      INTEGER Igrain , Nnear
      LOGICAL Checkmodel
 
!--Local Varibles
 
      LOGICAL POINTINGRAIN , foundit , POINTINPOLY , pim
      DOUBLE PRECISION xx , yy , une(NDM) , TOL , cx , cy , ul(NDM) , &
     &                 ymin , zmin , uat(3) , umin(3,NNEAMAX) , xcorg , &
     &                 ycorg
      DOUBLE PRECISION , POINTER :: cell(:,:)
      DOUBLE PRECISION rmin(NNEAMAX) , r , xc , yc , rtmp , utmp(3) , &
     &                 xtmp(3) , xat(3)
      INTEGER i , j , n , m , k , k1 , ncell , nvt , icell , idummy
      PARAMETER (TOL=1.0D-6)
 
      IF ( NDM==3 ) STOP 'NearBSiteInGrain not ready for 3D'
      ALLOCATE (cell(3,GRAins(Igrain)%ncell))
      DO i = 1 , NDM
         une(i) = U(i)
      ENDDO
 
      IF ( .NOT.POINTINGRAIN(U,Igrain) )&
     &     CALL FINDNEARESTPOINT(U,GRAins(Igrain)%NUMVRTS,GRAins(Igrain)&
     &     %GRAINARR,une,idummy)
 
 
      !load required data to local arrays
      nvt = GRAins(Igrain)%NUMVRTS
      ncell = GRAins(Igrain)%ncell
      cx = GRAins(Igrain)%DCELL(1)
      cy = GRAins(Igrain)%DCELL(2)
      DO i = 1 , ncell
         DO j = 1 , 3
            cell(j,i) = GRAins(Igrain)%cell(j,i)
         ENDDO
      ENDDO
 
 
      !find corrdinates of the point in local system
      ul(1) = GRAins(Igrain)%ROTMAT(1)*(U(1)-GRAins(Igrain)%REFATOM(1))&
     &        + GRAins(Igrain)%ROTMAT(2)&
     &        *(U(2)-GRAins(Igrain)%REFATOM(2))
      ul(2) = -GRAins(Igrain)%ROTMAT(2)*(U(1)-GRAins(Igrain)%REFATOM(1))&
     &        + GRAins(Igrain)%ROTMAT(1)&
     &        *(U(2)-GRAins(Igrain)%REFATOM(2))
 
 
      !shift origin to reference atom
      DO i = 1 , NDM
         une(i) = une(i) - GRAins(Igrain)%REFATOM(i)
      ENDDO
 
      !get into grain coordinates
      xx = une(1)*GRAins(Igrain)%ROTMAT(1) + une(2)*GRAins(Igrain)&
     &     %ROTMAT(2)
      yy = -une(1)*GRAins(Igrain)%ROTMAT(2) + une(2)*GRAins(Igrain)&
     &     %ROTMAT(1)
 
      !find which repeating cell
      m = INT(xx*(1.0+TOL)/cx)
      n = INT(yy*(1.0+TOL)/cy)
      IF ( xx<TOL .AND. DABS(DBLE(m)*cx-xx)>TOL ) m = m - 1
      IF ( yy<-TOL .AND. DABS(DBLE(n)*cy-yy)>TOL ) n = n - 1
      xcorg = DBLE(m)*cx
      ycorg = DBLE(n)*cy
 
 
      DO i = 1 , Nnear
         rmin(i) = 1.D16
      ENDDO
      DO i = -2 , 2 , 1
         DO j = -2 , 2 , 1
            xc = xcorg + DBLE(i)*cx
            yc = ycorg + DBLE(j)*cy
            DO icell = 1 , ncell
               uat(1) = xc + cell(1,icell) + GRAins(Igrain)%RFUDGVEC(1)
               uat(2) = yc + cell(2,icell) + GRAins(Igrain)%RFUDGVEC(2)
               uat(3) = cell(3,icell)
               IF ( POINTINPOLY(uat,nvt,GRAins(Igrain)%ROTGRAIN(1,1)) )&
     &              THEN
                  uat(1) = uat(1) - GRAins(Igrain)%RFUDGVEC(1)
                  uat(2) = uat(2) - GRAins(Igrain)%RFUDGVEC(2)
                  xat(1) = GRAins(Igrain)%ROTMAT(1)*uat(1)&
     &                     - GRAins(Igrain)%ROTMAT(2)*uat(2)
                  xat(2) = GRAins(Igrain)%ROTMAT(2)*uat(1)&
     &                     + GRAins(Igrain)%ROTMAT(1)*uat(2)
                  xat(1) = xat(1) + GRAins(Igrain)%REFATOM(1)
                  xat(2) = xat(2) + GRAins(Igrain)%REFATOM(2)
                  xat(3) = uat(3)
                  r = (uat(1)-ul(1))**2 + (uat(2)-ul(2))**2
                  DO k = 1 , Nnear
                     pim = .TRUE.
                     IF ( (r<=rmin(k)) .AND. pim ) THEN
                        rtmp = rmin(k)
                        rmin(k) = r
                        r = rtmp
                        DO k1 = 1 , 3
                           utmp(k1) = umin(k1,k)
                           umin(k1,k) = uat(k1)
                           uat(k1) = utmp(k1)
                           xtmp(k1) = Xa(k1,k)
                           Xa(k1,k) = xat(k1)
                           xat(k1) = xtmp(k1)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
 
      foundit = .TRUE.
      DO i = 1 , Nnear
         IF ( rmin(k)>0.9E16 ) foundit = .FALSE.
      ENDDO
 
      IF ( foundit ) THEN
         DO i = 1 , Nnear
            Rd(i) = rmin(i)
         ENDDO
      ELSE
         PRINT * , '***Error : Could not locate all near atoms '
         STOP
      ENDIF
      DEALLOCATE (cell)
      END SUBROUTINE NEARBSITEINGRAIN
!*==nearbsiteingrain1.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2
 
!**------------------------------------------------------------------
!** NearBSiteInGrain : Get Bravais site nearest to given point in the gr
!**
!**      Parameters :
!**                u  (in) : coordinates of the point
!**            nnear  (in) : number of neartest B-sites requested
!**           igrain  (in) : the grain that is required
!**       CheckModel  (in) : Flag true if point must lie in the model
!**               xa (out) : coordinates of the B-sites
!**               rd (out) : distance square of near B-site
!**
!**      Algorithm :
!**           check if point is in the grain
!**           transform point to grain coordinates.
!**           ff (exterior point) then
!**              Find the nearest interior point
!**           end if
!**           find Cell and its eight neighbours
!**           for all B-sites in the cell,
!**              if (not interior) discard
!**              find the smallest
!**           end do
!**
!**      Notes :
!**           Does not change polygrain.
!**           This is a very messy sub.
!**
!--
      SUBROUTINE NEARBSITEINGRAIN1(U,Igrain,Checkmodel,Xa,Rd)
 
      USE MOD_GRAIN
      USE MOD_MATERIAL
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--NEARBSITEINGRAIN1301
 
!--Variables Transferred
      DOUBLE PRECISION U(NDM) , Xa(3) , Rd
      INTEGER Igrain , nnear
      LOGICAL Checkmodel
 
!--Local Varibles
 
      LOGICAL POINTINGRAIN
      DOUBLE PRECISION une(NDM) , TOL
      INTEGER idummy
      DOUBLE PRECISION a(3,3) , b(3,3) , bvec(3,3)
      PARAMETER (TOL=1.0D-6)
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
 
      IF ( NDM==3 ) STOP 'NearBSiteInGrain not ready for 3D'
 
      a = GRAins(Igrain)%XLATVECT
      b = MATerial(GRAins(Igrain)%MATGRAIN)%bvec
      bvec = MATMUL(a,b)
      a = 0.D0
      a(3,3) = 1.D0
      a(1,1) = GRAins(Igrain)%ROTMAT(1)
      a(2,2) = GRAins(Igrain)%ROTMAT(1)
      a(1,2) = -GRAins(Igrain)%ROTMAT(2)
      a(2,1) = GRAins(Igrain)%ROTMAT(2)
      bvec = MATMUL(a,bvec)
      IF ( POINTINGRAIN(U,Igrain) ) THEN
         CALL NEARESTBRAVAIS(bvec,GRAins(Igrain)%REFATOM,U,Xa,Rd,&
     &                       Checkmodel,Igrain)
      ELSE
         CALL FINDNEARESTPOINT(U,GRAins(Igrain)%NUMVRTS,GRAins(Igrain)&
     &                         %GRAINARR,une,idummy)
         CALL NEARESTBRAVAIS(bvec,GRAins(Igrain)%REFATOM,une,Xa,Rd,&
     &                       Checkmodel,Igrain)
         Rd = SQRT((U(1)-Xa(1))**2+(U(2)-Xa(2))**2)
      ENDIF
 
      END SUBROUTINE NEARBSITEINGRAIN1
!*==nearestbravais.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE NEARESTBRAVAIS(Bvec,Org,U,X,Rd,Checkmodel,Igrain)
      IMPLICIT NONE
!*--NEARESTBRAVAIS346
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
      DOUBLE PRECISION Bvec(3,3) , Org(3) , U(2) , X(3) , ZERO , ONE , &
     &                 BIG , Rd , r2 , r2best , u3(3) , binv(3,3) , &
     &                 vec(3) , p(3)
      INTEGER i , j , l(3) , i1 , i2 , i3 , Igrain
      PARAMETER (ZERO=0.D0,ONE=1.0D0,BIG=1.E6)
      LOGICAL Checkmodel , POINTINGRAIN , pim
      u3(1) = U(1) - Org(1)
      u3(2) = U(2) - Org(2)
      u3(3) = ZERO - Org(3)
      CALL INV33(Bvec,binv)
      DO i = 1 , 3
         vec(i) = 0
         DO j = 1 , 3
            vec(i) = vec(i) + binv(i,j)*u3(j)
         ENDDO
         IF ( vec(i)<ZERO ) THEN
            l(i) = INT(vec(i)-ONE)
         ELSE
            l(i) = INT(vec(i))
         ENDIF
      ENDDO
      r2best = BIG
      DO i1 = l(1) , l(1) + 1
         DO i2 = l(2) , l(2) + 1
            DO i3 = l(3) , l(3) + 1
               p = Org
               DO i = 1 , 3
                  p(i) = p(i) + Bvec(i,1)*i1 + Bvec(i,2)*i2 + Bvec(i,3)&
     &                   *i3
               ENDDO
               pim = .TRUE.
               IF ( pim ) pim = POINTINGRAIN(p,Igrain)
               IF ( pim ) THEN
                  r2 = (p(1)-U(1))**2 + (p(2)-U(2))**2
                  IF ( r2<r2best ) THEN
                     r2best = r2
                     X = p
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      Rd = SQRT(r2best)
      END SUBROUTINE NEARESTBRAVAIS
!*==findnearestpoint.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
 
 
!**------------------------------------------------------------------
!** FindNearestPoint : Find the nearest point in the polygon to the
!**                    given point.
!**
!**      Parameters :
!**                uin   (in) : the given point
!**                nvt   (in) : number of vertices in the polygon
!**           nvt_phys   (in) : physical array dimension
!**               poly   (in) : the polygon
!**                  v  (out) : the nearest point
!**               iedge (out) : edge on which the nearest point is locat
!**
!**      Algorithm :
!**             Scan over all edges in the polygon,
!**             find the nearest point on each edge,
!**             take the nearest of the nearest
!**
!**      Notes :
!**             Works only for 2D.
!**
!**      Changes :
!**            VBS, Aug 24, 1996 : Returns the edge on which the nearest
!**                                point sits
!**
!--
      SUBROUTINE FINDNEARESTPOINT(Uin,Nvt,Poly,V,Iedge)
 
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--FINDNEARESTPOINT425
 
!--Variables Transferred
      INTEGER Nvt , Iedge
      DOUBLE PRECISION Uin(NDM) , Poly(NDM,*) , V(NDM)
 
!--Local Variables
      INTEGER ivt , ivt1 , i
      DOUBLE PRECISION x , RINTERSECTRATIO , vtmp(3) , DIST2 , dbest , &
     &                 dtemp
!
      dbest = 1.E30
      DO i = 1 , 3
         vtmp(i) = 0.D0
      ENDDO
!
      DO ivt = 1 , Nvt
         ivt1 = ivt + 1
         IF ( ivt1>Nvt ) ivt1 = 1
!
!     intersection with edge(ivt)
!
         x = RINTERSECTRATIO(Poly(1,ivt),Poly(1,ivt1),Uin)
         IF ( x<0. ) THEN
            vtmp(1) = Poly(1,ivt)
            vtmp(2) = Poly(2,ivt)
         ELSEIF ( x>1. ) THEN
            vtmp(1) = Poly(1,ivt1)
            vtmp(2) = Poly(2,ivt1)
         ELSE
            vtmp(1) = Poly(1,ivt) + x*(Poly(1,ivt1)-Poly(1,ivt))
            vtmp(2) = Poly(2,ivt) + x*(Poly(2,ivt1)-Poly(2,ivt))
         ENDIF
         dtemp = DIST2(vtmp,Uin,2)
         IF ( dtemp<dbest ) THEN
            Iedge = ivt
            dbest = dtemp
            V(1) = vtmp(1)
            V(2) = vtmp(2)
         ENDIF
      ENDDO
      END SUBROUTINE FINDNEARESTPOINT
!*==rintersectratio.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
 
 
!**------------------------------------------------------------------
!** RIntersectRatio : Given two points and a third point, return the
!**                  intersection ratio of the perpendicular
!**                  line form the third point to the line fromed
!**                  by the first two points
!**
!**     Parameters :
!**                x1 (in) : coordinates of the first point
!**                x2 (in) : coordinates of the second point
!**                 a (in) : coordinates of the third point
!**
!**
!**     Algorithm :
!**                  o x2
!**                  |                u  = x2 - x1
!**                  |                u1 = a  - x1
!**                  |
!**                  |                                  u1 . u
!**                 -+-------o a   RInterSectRatio =  ---------
!**                  |      /                           |u|^2
!**                  |     /
!**                  |    /
!**                  |   /
!**                  |  /
!**                  | /
!**                  |/
!**                  o
!**                  x1
!**
!**
!--
      DOUBLE PRECISION FUNCTION RINTERSECTRATIO(X1,X2,A)
 
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--RINTERSECTRATIO505
 
!--Variables Transferred
      DOUBLE PRECISION X1(NDM) , X2(NDM) , A(NDM)
 
!--Local Variables
      DOUBLE PRECISION u(NDM) , u1(NDM) , d , dot
      INTEGER i
 
      DO i = 1 , NDM
         u(i) = X2(i) - X1(i)
         u1(i) = A(i) - X1(i)
      ENDDO
 
      d = 0.0
      dot = 0.0
 
      DO i = 1 , NDM
         d = d + u(i)*u(i)
         dot = dot + u(i)*u1(i)
      ENDDO
 
      RINTERSECTRATIO = dot/d
 
      END FUNCTION RINTERSECTRATIO
!*==pointingrain.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!**----------------------------------------------------------------
!** PointInGrain : Checks if given point is in given grain
!**
!**    Parameters :
!**               u (in)  : coordinates of the point
!**          igrain (in)  : grain number
!**
!**    Algorithm :
!**              Check using fudge vector
!**
!--
      LOGICAL FUNCTION POINTINGRAIN(U,Igrain)
      USE MOD_GRAIN
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--POINTINGRAIN547
 
!--Variables transferred
      DOUBLE PRECISION U(NDM)
      INTEGER Igrain
 
!--Local Variables
      LOGICAL POINTINPOLY
      DOUBLE PRECISION utmp(NDM)
      INTEGER i
 
      DO i = 1 , NDM
         utmp(i) = U(i) + GRAins(Igrain)%FUDGEVEC(i)
      ENDDO
      POINTINGRAIN = POINTINPOLY(utmp(1:NDM),GRAins(Igrain)%NUMVRTS,&
     &               GRAins(Igrain)%GRAINARR(1,1))
 
      END FUNCTION POINTINGRAIN
 
