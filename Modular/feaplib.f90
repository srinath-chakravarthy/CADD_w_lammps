!*==pcomp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!** FEAPLIB:  These are a lot of handy utilities and other routines that
!**   don't fit neatly elsewhere in the code.  Most are self
!**   explanatory.  Here is a brief table of contents, but note that the
!**   actual order of the routines is different.
!**
!** data dumping routines:
!**    putout          : dump a double precision array to std. output
!**    intout          : dump an integer array to std. output
!** vector and matrix manipulations:
!**    distance        : 2-D distance between 2 points
!**    Dist2           : n-D distance-squared between 2 points
!**    dot             : n-D dot product
!**    cross_product   : 3-D cross product
!**    norm            : n-D vector normalization
!**    promul          : computes c=c+a*b, a is stiffness, c and b are v
!**    eq_inte         : equate two integer vectors
!**    addb            : do x=x+b
!**    subb            : do x=x-b
!**    maxsearch       : find the maximum absolute value in a vector
!**    store           : store a vector in another vector
!** handy geometric routines:
!**    FindArea        : find area of a quadrilateral
!**    PolygonArea     : find area of an n-sided polygon
!**    PointInPoly     : check if a point is in a polygon
!**    getaspect       : get the aspect ratio of a triangle
!**    area2x          : convert area coordinates to cartesian
!**    between0        : check if the point zero is between two points
!**    areaix          : get the area of an element
!** character string manipulations:
!**    freein          : parse free-form input lines
!**    downcase        : make a string all lower case
!**    pcomp           : compare to 4-character strings
!**    next            : find the next delimiter (comma or space) in a s
!** Feap things:
!**    blank           : blanks unconstrained components of a vector
!**    simplex         : shape function calculations
!**    dot11           : dot product of unconstrained components
!** QC specific, but called in many places:
!**    MapToRef        : map a point from deformed to reference configur
!**    MapToCurrent    : map a point from reference to deformed configur
!**    TwoMaxSum       : find sum of two maximum displacements
!**    MustUpdateStatus: check if status must be updated due to large di
!**    RINFEIGEN       : find (rcut)*(max. eigenvalue of F)
!**    BALANC          : balance a matrix
!**    ELMHES          : get hessian form of matrix
!**    HQR             : get eigenvalues of a matrix
!***********************************************************************
      LOGICAL FUNCTION PCOMP(A,B)
      IMPLICIT NONE
!*--PCOMP51
!
      CHARACTER*4 A , B
!
      PCOMP = .FALSE.
      IF ( A==B ) PCOMP = .TRUE.
      END FUNCTION PCOMP
!*==next.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      FUNCTION NEXT(Last,Input)
      IMPLICIT NONE
!*--NEXT62
!*** Start of declarations inserted by SPAG
      INTEGER i , Last , NEXT
!*** End of declarations inserted by SPAG
!
      CHARACTER*80 Input
!
      DO i = Last + 1 , 80
         IF ( (Input(i:i)==',') .OR. (Input(i:i)==' ') ) THEN
            NEXT = i
            RETURN
         ENDIF
      ENDDO
      NEXT = 81
      END FUNCTION NEXT
!*==freein.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE FREEIN(Input,Ll,Lr,Int,Real,Isw)
      IMPLICIT NONE
!*--FREEIN81
!*** Start of declarations inserted by SPAG
      INTEGER Int , Isw , ldif , Ll , Lr
      DOUBLE PRECISION Real
!*** End of declarations inserted by SPAG
      CHARACTER*80 Input
      CHARACTER*1 range1
      CHARACTER*2 range2
      CHARACTER*80 string
!
      string = Input((Ll+1):(Lr-1))
      ldif = Lr - Ll - 1
      IF ( ldif==0 ) THEN
         IF ( Isw==1 ) Int = 0
         IF ( Isw==2 ) Real = 0.
      ELSEIF ( ldif<10 ) THEN
         ENCODE (1,'(i1)',RANGE1) ldif
         IF ( Isw==1 ) DECODE (ldif,'(i'//range1//')',STRING) Int
         IF ( Isw==2 ) DECODE (ldif,'(f'//range1//'.0)',STRING) Real
      ELSE
         ENCODE (2,'(i2)',RANGE2) ldif
         IF ( Isw==1 ) DECODE (ldif,'(i'//range2//')',STRING) Int
         IF ( Isw==2 ) DECODE (ldif,'(f'//range2//'.0)',STRING) Real
      ENDIF
      END SUBROUTINE FREEIN
!*==dist2.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------------
!** Dist2 : return distance square of between two points
!**
!**      Parameters :
!**              x, y (in) : two points
!**                 n (in) : dimensions
!--
      DOUBLE PRECISION FUNCTION DIST2(X,Y,N)
      IMPLICIT NONE
!*--DIST2116
 
      INTEGER N
      DOUBLE PRECISION X(N) , Y(N)
 
      DIST2 = DOT_PRODUCT(X-Y,X-Y)
      END FUNCTION DIST2
!*==pointinpoly.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!***********************************************************************
!** PointInPoly : Checks if given point is in given polygon
!**
!**    Parameters :
!**               u (in)  : coordinates of the point
!**            numv (in)  : number of vertices
!**         grainar (in)  : grain array (only one polygon)
!**
!**    Algorithm :
!**              The one by Joseph O'Rourke.
!**              Refer his book (look up in lib).
!**
!--
 
      LOGICAL FUNCTION POINTINPOLY(U,Numv,Grainar)
      IMPLICIT NONE
!*--POINTINPOLY141
 
 
!--Variables transferred
      DOUBLE PRECISION U(2) , Grainar(2,*)
      INTEGER Numv
 
!--Local variables
      COMMON /DEBUGGER/ DEBug
      LOGICAL DEBug
      INTEGER NVMAX
      PARAMETER (NVMAX=100)
      INTEGER icross , i , i1 , j
      DOUBLE PRECISION x , grain(2,NVMAX)
      LOGICAL BETWEEN0 , onvertex
      DOUBLE PRECISION TOLE
      PARAMETER (TOLE=1.0D-6)
 
!--Check data
      IF ( Numv>NVMAX ) THEN
         PRINT * , 'ERROR ** Too many vertices in PointInPoly'
         PRINT * , 'Increase NVMAX' , NVMAX , Numv
         STOP
      ENDIF
 
 
!--Set up local variables
      DO i = 1 , Numv
         DO j = 1 , 2
            grain(j,i) = Grainar(j,i) - U(j)
         ENDDO
      ENDDO
 
 
!--For each edge, check if it crosses the ray
      icross = 0
      DO i = 1 , Numv
         i1 = MOD((i+Numv-2),Numv) + 1
         ! if edge straddles the xaxis ...
         IF ( ((grain(2,i)>0.0) .AND. (grain(2,i1)<=0.0)) .OR. &
     &        ((grain(2,i1)>0.0) .AND. (grain(2,i)<=0.0)) ) THEN
            ! edge straddles ray, compute intersection
            x = (grain(1,i)*grain(2,i1)-grain(1,i1)*grain(2,i))&
     &          /(grain(2,i1)-grain(2,i))
            IF ( x>0.0 ) icross = icross + 1
!--A point on the boundary is considered in the polygon.
         ELSEIF ( BETWEEN0(grain(1,i),grain(1,i1)) ) THEN
            POINTINPOLY = .TRUE.
            RETURN
         ENDIF
      ENDDO
 
 
      !u is inside for odd number of crossings
      IF ( MOD(icross,2)==1 ) THEN
         POINTINPOLY = .TRUE.
         RETURN
      ELSE
!--If all else fails make sure that the point is not one of the vertices
         onvertex = .FALSE.
         DO i = 1 , Numv
            IF ( DABS(grain(1,i))<TOLE .AND. DABS(grain(2,i))<TOLE )&
     &           onvertex = .TRUE.
         ENDDO
         IF ( onvertex ) THEN
            POINTINPOLY = .TRUE.
         ELSE
            POINTINPOLY = .FALSE.
         ENDIF
         RETURN
      ENDIF
 
      END FUNCTION POINTINPOLY
!*==between0.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!***********************************************************************
      LOGICAL FUNCTION BETWEEN0(X,Y)
! checks if the point zero fall on the closed segment between x and y
      IMPLICIT NONE
!*--BETWEEN0220
      DOUBLE PRECISION X(2) , Y(2) , area2 , ZERO , TOL
      PARAMETER (ZERO=0.D0,TOL=1.E-9)
      BETWEEN0 = .FALSE.
!--make sure the three points are colinear.
      IF ( ABS(X(2)*Y(1)-X(1)*Y(2))>TOL ) RETURN
!--check betweenness.
      IF ( X(1)/=Y(1) ) THEN
         BETWEEN0 = (((X(1)>=ZERO) .AND. (Y(1)<=ZERO)) .OR. &
     &              ((X(1)<=ZERO) .AND. (Y(1)>=ZERO)))
      ELSE
         BETWEEN0 = (((X(2)>=ZERO) .AND. (Y(2)<=ZERO)) .OR. &
     &              ((X(2)<=ZERO) .AND. (Y(2)>=ZERO)))
      ENDIF
      END FUNCTION BETWEEN0
!*==between.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      LOGICAL FUNCTION BETWEEN(X1,X2,Y)
      IMPLICIT NONE
!*--BETWEEN239
      DOUBLE PRECISION X1(2) , X2(2) , Y(2) , RINTERSECTRATIO , r1
      r1 = RINTERSECTRATIO(X1,X2,Y)
      BETWEEN = (r1>=0.D0) .AND. (r1<=1.D0)
      END FUNCTION BETWEEN
!*==getaspect.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!***********************************************************************
!** get the aspect ratio of a triangle
      SUBROUTINE GETASPECT(Xtri,Aspect)
      IMPLICIT NONE
!*--GETASPECT250
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION area , Aspect , det , equi , perim , Xtri
!*** End of declarations inserted by SPAG
!
      DIMENSION Xtri(2,3)
!
      equi = 36.D0/DSQRT(3.D0)
      det = Xtri(1,2)*Xtri(2,3) - Xtri(1,3)*Xtri(2,2) + Xtri(1,3)&
     &      *Xtri(2,1) - Xtri(1,1)*Xtri(2,3) + Xtri(1,1)*Xtri(2,2)&
     &      - Xtri(1,2)*Xtri(2,1)
      area = det/2.D0
      perim = DSQRT((Xtri(1,2)-Xtri(1,1))**2+(Xtri(2,2)-Xtri(2,1))**2)&
     &        + DSQRT((Xtri(1,3)-Xtri(1,2))**2+(Xtri(2,3)-Xtri(2,2))**2)&
     &        + DSQRT((Xtri(1,1)-Xtri(1,3))**2+(Xtri(2,1)-Xtri(2,3))**2)
      Aspect = equi*area/(perim*perim)
      END SUBROUTINE GETASPECT
!*==areaix.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      DOUBLE PRECISION FUNCTION AREAIX(Kount,Ixnew,Xnew,Nxdm)
      IMPLICIT NONE
!*--AREAIX271
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION det , Xnew , xtri
      INTEGER i , ii , Ixnew , j , Kount , Nxdm
!*** End of declarations inserted by SPAG
      DIMENSION Ixnew(3,1) , Xnew(Nxdm,1) , xtri(2,3)
      DO i = 1 , 3
         DO j = 1 , 2
            ii = Ixnew(i,Kount)
            xtri(j,i) = Xnew(j,ii)
         ENDDO
      ENDDO
      det = xtri(1,2)*xtri(2,3) - xtri(1,3)*xtri(2,2) + xtri(1,3)&
     &      *xtri(2,1) - xtri(1,1)*xtri(2,3) + xtri(1,1)*xtri(2,2)&
     &      - xtri(1,2)*xtri(2,1)
      AREAIX = det/2.D0
      END FUNCTION AREAIX
!*==cross_product.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE CROSS_PRODUCT(X1,X2,X3)
!
!     finds cross prod of x1,x2 --> x3
!
      IMPLICIT NONE
!*--CROSS_PRODUCT295
      DOUBLE PRECISION X1(3) , X2(3) , X3(3)
      X3(1) = X1(2)*X2(3) - X1(3)*X2(2)
      X3(2) = X1(3)*X2(1) - X1(1)*X2(3)
      X3(3) = X1(1)*X2(2) - X1(2)*X2(1)
      END SUBROUTINE CROSS_PRODUCT
!*==intri.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
!     Logical function intri returns true when point p is inside
!     triangle with vertices x1,x2,x3 and transforms p coordinates
!     to area coordinates in s. ontri is set to true if point p
!     is on the triangle boundary.
      LOGICAL FUNCTION INTRI(X1,X2,X3,P,S,Ontri)
      IMPLICIT NONE
!*--INTRI309
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION a , a1 , a2 , a3 , P , S , TOL , tola , TRAREA , &
     &                 X1 , X2 , X3
!*** End of declarations inserted by SPAG
      LOGICAL Ontri
      PARAMETER (TOL=1.D-9)
      DIMENSION X1(2) , X2(2) , X3(2) , P(2) , S(3)
!     Calculate area of triangle
      a = DABS(TRAREA(X1,X2,X3))
      IF ( a==0.D0 ) THEN
         WRITE (6,*) '***ERROR: Model contains coincident nodes.'
         WRITE (6,*) X1(1) , X1(2)
         WRITE (6,*) X2(1) , X2(2)
         WRITE (6,*) X3(1) , X3(2)
         STOP
      ENDIF
!     Calculate area of 3 new triangles formed by connecting p with
!     with vertices
      a1 = DABS(TRAREA(P,X2,X3))
      a2 = DABS(TRAREA(X1,P,X3))
      a3 = DABS(TRAREA(X1,X2,P))
!     If p is inside triangle we should have a1+a2+a3 = a
      tola = TOL*a
      INTRI = (DABS(a1+a2+a3-a)<=tola)
      Ontri = .FALSE.
      IF ( INTRI ) THEN
         Ontri = (a1<=tola .OR. a2<=tola .OR. a3<=tola)
         S(1) = a1/a
         S(2) = a2/a
         S(3) = a3/a
      ENDIF
      END FUNCTION INTRI
!*==trarea.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
!     Area of a triangle with vertices at x1,x2,x3
      FUNCTION TRAREA(X1,X2,X3)
      IMPLICIT NONE
!*--TRAREA347
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION TRAREA , X1 , X2 , X3
!*** End of declarations inserted by SPAG
      DIMENSION X1(2) , X2(2) , X3(2)
      TRAREA = 0.5*(X1(1)*(X2(2)-X3(2))-X1(2)*(X2(1)-X3(1))+X2(1)*X3(2)&
     &         -X3(1)*X2(2))
      END FUNCTION TRAREA
!*==inv33.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE INV33(A,B)
      IMPLICIT NONE
!*--INV33359
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION A , B , det , DET33
!*** End of declarations inserted by SPAG
      DIMENSION A(3,3) , B(3,3)
      det = DET33(A)
      B(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
      B(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3))/det
      B(1,3) = (A(1,2)*A(2,3)-A(2,2)*A(1,3))/det
      B(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3))/det
      B(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
      B(2,3) = (A(2,1)*A(1,3)-A(1,1)*A(2,3))/det
      B(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
      B(3,2) = (A(3,1)*A(1,2)-A(1,1)*A(3,2))/det
      B(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
      END SUBROUTINE INV33
!*==det33.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      DOUBLE PRECISION FUNCTION DET33(A)
      IMPLICIT NONE
!*--DET33379
!*** Start of declarations inserted by SPAG
      DOUBLE PRECISION A
!*** End of declarations inserted by SPAG
      DIMENSION A(3,3)
      DET33 = A(1,1)*A(2,2)*A(3,3) + A(2,1)*A(3,2)*A(1,3) + A(1,2)&
     &        *A(2,3)*A(3,1) - A(3,1)*A(2,2)*A(1,3) - A(1,2)*A(2,1)&
     &        *A(3,3) - A(3,2)*A(2,3)*A(1,1)
      END FUNCTION DET33
!*==linecoef.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE LINECOEF(X1,X2,Abc)
      IMPLICIT NONE
!*--LINECOEF392
!
!     given two points on a line x1 and x2, find the equation of the
!     line:   abc(1)*x + abc(2)*y + abc(3) = 0
!
      DOUBLE PRECISION X1(2) , X2(2) , Abc(3) , bot
      IF ( X1(1)==0.D0 .AND. X2(1)==0.D0 ) THEN
         Abc(1) = 1.D0
         Abc(2) = 0.D0
         Abc(3) = 0.D0
      ELSEIF ( X1(2)==0.D0 .AND. X2(2)==0.D0 ) THEN
         Abc(1) = 0.D0
         Abc(2) = 1.D0
         Abc(3) = 0.D0
      ELSE
         bot = X1(1)*X2(2) - X2(1)*X1(2)
         Abc(1) = (X1(2)-X2(2))/bot
         Abc(2) = (-X1(1)+X2(1))/bot
         Abc(3) = 1.D0
      ENDIF
      END SUBROUTINE LINECOEF
 
