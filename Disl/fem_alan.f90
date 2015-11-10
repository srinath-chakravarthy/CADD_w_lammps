!*==fe_makemat.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
! $Id: fem_alan.f,v 1.22004/04/2114:29:39 shastry Exp $
!
      SUBROUTINE FE_MAKEMAT()
!
!12345678901234567890123456789012345678901234567890123456789012345678901
!         1         2         3         4         5         6         7
!
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_MAKEMAT12
      INTEGER MAXSTS , MAXSTN
      PARAMETER (MAXSTS=4)
      PARAMETER (MAXSTN=3)
      DOUBLE PRECISION d(MAXSTS,MAXSTN)
      DOUBLE PRECISION xl(KNODE*NDOF) , s(KNODE*NDOF,KNODE*NDOF)
!
      INTEGER nodedf , mdif , i1 , i2 , mtest , i , j
      CHARACTER*80 error_message
!
      CALL FE_DMAT(d)
!
! compute bandwidth
!
      NEQu = NDOF*NNOdes
!
      nodedf = 0
      DO i = 1 , NELm
         mdif = 0
         DO i1 = 1 , KNODE
            DO i2 = 1 , KNODE
               mtest = ABS(ICOnn(i1,i)-ICOnn(i2,i))
               IF ( mtest>=mdif ) mdif = mtest
            ENDDO
         ENDDO
         IF ( mdif>=nodedf ) nodedf = mdif
      ENDDO
      MBAndw = NDOF*(nodedf+1)
      WRITE (*,99001) NEQu , MBAndw
99001 FORMAT (' no. equations, bandwidth= ',2I8)
      IF ( MBAndw>MAXBND ) THEN
         error_message = ' *** not enough storage *** '
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , NEQu
         DO j = 1 , MBAndw
            A_Stiff(j,i) = 0.0D0
            AD_stiff(j,i) = 0.0D0
         ENDDO
      ENDDO
!
      DO i = 1 , NELm
         CALL FE_FIND(i,xl)
         CALL FE_STIFF(xl,s,d)
         CALL FE_PLACE(i,s)
      ENDDO
      END SUBROUTINE FE_MAKEMAT
!*==fe_find.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE FE_FIND(Lmn,Xl)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_FIND66
      INTEGER Lmn
      DOUBLE PRECISION Xl(KNODE*NDOF)
      INTEGER j , k
!
      DO j = 1 , KNODE
         DO k = 1 , NDOF
            Xl(NDOF*(j-1)+k) = X0(k,ICOnn(j,Lmn))
         ENDDO
      ENDDO
      END SUBROUTINE FE_FIND
!*==fe_dmat.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_DMAT(D)
      IMPLICIT NONE
!*--FE_DMAT83
      INTEGER MAXSTS , MAXSTN
      PARAMETER (MAXSTS=4)
      PARAMETER (MAXSTN=3)
      DOUBLE PRECISION D(MAXSTS,MAXSTN)
      CHARACTER*80 error_message
!
      DOUBLE PRECISION CC(6,6)
      DOUBLE PRECISION XE , XNU , XLAmbda , XMU
      INTEGER I_Elas
!
      COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
!
!  1<=> eps(1,1)
!  2<=> eps(2,2)
!  3<=> 2*eps(1,2)
!
!  1<=> s(1,1)
!  2<=> s(2,2)
!  3<=> s(1,2)
!  4<=> s(3,3)
!
      IF ( I_Elas/=1 ) THEN
         error_message = 'fe_dmat: call fe_elastic first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      D(1,1) = CC(1,1)
      D(1,2) = CC(1,2)
      D(1,3) = CC(1,6)
      D(2,1) = CC(2,1)
      D(2,2) = CC(2,2)
      D(2,3) = CC(2,6)
      D(3,1) = CC(6,1)
      D(3,2) = CC(6,2)
      D(3,3) = CC(6,6)
      D(4,1) = CC(3,1)
      D(4,2) = CC(3,2)
      D(4,3) = CC(3,6)
      END SUBROUTINE FE_DMAT
!*==fe_stiff.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_STIFF(Xl,S,D)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_STIFF130
      INTEGER MAXSTS , MAXSTN
      PARAMETER (MAXSTS=4)
      PARAMETER (MAXSTN=3)
      DOUBLE PRECISION D(MAXSTS,MAXSTN)
      DOUBLE PRECISION Xl(KNODE*NDOF) , S(KNODE*NDOF,KNODE*NDOF)
      DOUBLE PRECISION pn(KNODE) , qn(KNODE)
      DOUBLE PRECISION bpq(NDOF,NDOF,KNODE*NDOF) , b(MAXSTN,KNODE*NDOF)
      DOUBLE PRECISION xjac(NDOF,NDOF) , xji(NDOF,NDOF)
!
      INTEGER i , j , i1 , i2 , j1 , j2
      DOUBLE PRECISION det , area
!
! Everything is hard-coded, since I don't know where I may need
! the derivatives of shape functions
!
!  1<=> eps(1,1)
!  2<=> eps(2,2)
!  3<=> 2*eps(1,2)
!
! Shape functions are: N_1 = 1 - p - q; N_2 = p; N_3 = q;
!
      pn(1) = -1.0D0
      pn(2) = 1.0D0
      pn(3) = 0.0D0
      qn(1) = -1.0D0
      qn(2) = 0.0D0
      qn(3) = 1.0D0
!
      DO i = 1 , KNODE*NDOF
         DO i2 = 1 , NDOF
            DO i1 = 1 , NDOF
               bpq(i1,i2,i) = 0.0D0
            ENDDO
         ENDDO
      ENDDO
!
      DO i = 1 , KNODE
         bpq(1,1,NDOF*(i-1)+1) = pn(i)
         bpq(1,2,NDOF*(i-1)+1) = qn(i)
         bpq(2,1,NDOF*(i-1)+NDOF) = pn(i)
         bpq(2,2,NDOF*(i-1)+NDOF) = qn(i)
      ENDDO
!
      DO i = 1 , NDOF
         DO j = 1 , NDOF
            xjac(i,j) = 0.0D0
         ENDDO
      ENDDO
!
      DO i = 1 , KNODE
         xjac(1,1) = xjac(1,1) + pn(i)*Xl((i-1)*NDOF+1)
         xjac(1,2) = xjac(1,2) + qn(i)*Xl((i-1)*NDOF+1)
         xjac(2,1) = xjac(2,1) + pn(i)*Xl((i-1)*NDOF+2)
         xjac(2,2) = xjac(2,2) + qn(i)*Xl((i-1)*NDOF+2)
      ENDDO
!
      det = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
!
      xji(1,1) = xjac(2,2)/det
      xji(2,2) = xjac(1,1)/det
      xji(1,2) = -xjac(1,2)/det
      xji(2,1) = -xjac(2,1)/det
!
      area = ABS(det)/2.0D0
!
      DO i1 = 1 , KNODE*NDOF
         b(1,i1) = bpq(1,1,i1)*xji(1,1) + bpq(1,2,i1)*xji(2,1)
         b(2,i1) = bpq(2,1,i1)*xji(1,2) + bpq(2,2,i1)*xji(2,2)
         b(3,i1) = bpq(1,1,i1)*xji(1,2) + bpq(1,2,i1)*xji(2,2)&
     &             + bpq(2,1,i1)*xji(1,1) + bpq(2,2,i1)*xji(2,1)
      ENDDO
!
      DO i1 = 1 , NDOF*KNODE
         DO i2 = 1 , NDOF*KNODE
            S(i1,i2) = 0.0D0
            DO j1 = 1 , MAXSTN
               DO j2 = 1 , MAXSTN
                  S(i1,i2) = S(i1,i2) + b(j1,i1)*D(j1,j2)*b(j2,i2)
               ENDDO
            ENDDO
            S(i1,i2) = S(i1,i2)*area
         ENDDO
      ENDDO
!
      END SUBROUTINE FE_STIFF
!*==fe_place.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_PLACE(Lmn,S)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_PLACE223
      INTEGER Lmn
      DOUBLE PRECISION S(KNODE*NDOF,KNODE*NDOF)
      INTEGER i1 , i2 , i , k , j1 , j2 , j , l
!
      DO i1 = 1 , NDOF
         DO i2 = 1 , KNODE
            i = NDOF*(i2-1) + i1
            k = NDOF*(ICOnn(i2,Lmn)-1) + i1
            DO j1 = 1 , NDOF
               DO j2 = 1 , KNODE
                  j = NDOF*(j2-1) + j1
                  l = NDOF*(ICOnn(j2,Lmn)-1) + j1
                  IF ( k>=l ) THEN
                     A_Stiff(k-l+1,l) = A_Stiff(k-l+1,l) + S(i,j)
                     AD_stiff(k-l+1,l) = AD_stiff(k-l+1,l) + S(i,j)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE FE_PLACE
!*==fe_fixdsp.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE FE_FIXDSP()
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_FIXDSP251
      INTEGER i , ic , j
      DOUBLE PRECISION hold_stiff
!
      hold_stiff = 0.0D0
      DO i = 1 , NEQu
         IF ( ABS(AD_stiff(1,i))>hold_stiff ) THEN
            hold_stiff = ABS(AD_stiff(1,i))
         ENDIF
      ENDDO
      DO i = 1 , NFIxed
         ic = IFIxed(i)
         DO j = 1 , NEQu
            IF ( (j<ic) .AND. (ic-j+1<=MBAndw) ) THEN
               AD_stiff(ic-j+1,j) = 0.0D0
            ELSEIF ( (j>ic) .AND. (j-ic+1<=MBAndw) ) THEN
               AD_stiff(j-ic+1,ic) = 0.0D0
            ENDIF
         ENDDO
         AD_stiff(1,ic) = hold_stiff
      ENDDO
      END SUBROUTINE FE_FIXDSP
!*==fe_substitute.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE FE_SUBSTITUTE(Rhs,Presv)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_SUBSTITUTE278
      DOUBLE PRECISION Rhs(*) , Presv(*)
      INTEGER i , ic , j , k , kmin , kmax
!
      DO i = 1 , NFIxed
         ic = IFIxed(i)
         DO j = 1 , NEQu
            IF ( (j<ic) .AND. (ic-j+1<=MBAndw) ) THEN
               Rhs(j) = Rhs(j) - Presv(i)*A_Stiff(ic-j+1,j)
            ELSEIF ( (j>ic) .AND. (j-ic+1<=MBAndw) ) THEN
               Rhs(j) = Rhs(j) - Presv(i)*A_Stiff(j-ic+1,ic)
            ENDIF
         ENDDO
      ENDDO
      DO i = 1 , NFIxed
         ic = IFIxed(i)
         Rhs(ic) = Presv(i)*AD_stiff(1,ic)
      ENDDO
!
! substitution
!
      DO i = 2 , NEQu
         kmin = i + 1 - MBAndw
         IF ( kmin<1 ) kmin = 1
         DO k = kmin , (i-1)
            Rhs(i) = Rhs(i) - AD_stiff(i-k+1,k)*Rhs(k)
         ENDDO
      ENDDO
!
! divide by diagonals
!
      DO i = 1 , NEQu
         Rhs(i) = Rhs(i)/AD_stiff(1,i)
      ENDDO
!
! back substitution
!
      DO i = (NEQu-1) , 1 , -1
         kmax = i - 1 + MBAndw
         IF ( kmax>NEQu ) kmax = NEQu
         DO k = (i+1) , kmax
            Rhs(i) = Rhs(i) - AD_stiff(k-i+1,i)*Rhs(k)
         ENDDO
      ENDDO
      END SUBROUTINE FE_SUBSTITUTE
!*==fe_chol.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_CHOL(Ier)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_CHOL330
      INTEGER Ier , i , k , j , kmin , jmax
      DOUBLE PRECISION piv , test
!
      Ier = 0
      test = 0.0D0
      DO i = 1 , NEQu
         IF ( AD_stiff(1,i)>test ) test = AD_stiff(1,i)
      ENDDO
      test = 1.D-8*test
!
! factorization
!
      IF ( AD_stiff(1,1)<test ) THEN
         Ier = 1
         RETURN
      ENDIF
      piv = 1.0D0/AD_stiff(1,1)
      DO k = 2 , MBAndw
         AD_stiff(k,1) = piv*AD_stiff(k,1)
      ENDDO
      DO i = 2 , NEQu
         kmin = i + 1 - MBAndw
         IF ( kmin<1 ) kmin = 1
         DO k = kmin , (i-1)
            AD_stiff(1,i) = AD_stiff(1,i) - AD_stiff(i-k+1,k)&
     &                      *AD_stiff(1,k)*AD_stiff(i-k+1,k)
         ENDDO
         IF ( i<NEQu ) THEN
            IF ( AD_stiff(1,i)<test ) THEN
               Ier = i
               RETURN
            ENDIF
            piv = 1.0D0/AD_stiff(1,i)
            jmax = i - 1 + MBAndw
            IF ( jmax>NEQu ) jmax = NEQu
            DO j = (i+1) , jmax
               kmin = j + 1 - MBAndw
               IF ( kmin<1 ) kmin = 1
               DO k = kmin , (i-1)
                  AD_stiff(j-i+1,i) = AD_stiff(j-i+1,i)&
     &                                - AD_stiff(i-k+1,k)*AD_stiff(1,k)&
     &                                *AD_stiff(j-k+1,k)
               ENDDO
               AD_stiff(j-i+1,i) = piv*AD_stiff(j-i+1,i)
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE FE_CHOL
!*==fe_strain.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_STRAIN(Lmn,Rhs,Ee_out)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_STRAIN386
      INTEGER MAXSTS , MAXSTN
      PARAMETER (MAXSTS=4)
      PARAMETER (MAXSTN=3)
!
      INTEGER Lmn
      DOUBLE PRECISION Rhs(*) , Ee_out(MAXSTN)
!
      DOUBLE PRECISION xl(KNODE*NDOF) , s(KNODE*NDOF,KNODE*NDOF)
      DOUBLE PRECISION pn(KNODE) , qn(KNODE)
      DOUBLE PRECISION bpq(NDOF,NDOF,KNODE*NDOF) , b(MAXSTN,KNODE*NDOF)
      DOUBLE PRECISION xjac(NDOF,NDOF) , xji(NDOF,NDOF)
!
      INTEGER i , j , i1 , i2 , j1 , j2
      DOUBLE PRECISION det
!
      CHARACTER*80 error_message
!
      IF ( (Lmn<1) .OR. (Lmn>NELm) ) THEN
         PRINT * , Lmn
         error_message = 'fe_strain: element is out of range!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
 
      CALL FE_FIND(Lmn,xl)
!
! Everything is hard-coded, since I don't know where I may need
! the derivatives of shape functions
!
!  1<=> eps(1,1)
!  2<=> eps(2,2)
!  3<=> 2*eps(1,2)
!
! Shape functions are: N_1 = 1 - p - q; N_2 = p; N_3 = q;
!
      pn(1) = -1.0D0
      pn(2) = 1.0D0
      pn(3) = 0.0D0
      qn(1) = -1.0D0
      qn(2) = 0.0D0
      qn(3) = 1.0D0
!
      DO i = 1 , KNODE*NDOF
         DO i2 = 1 , NDOF
            DO i1 = 1 , NDOF
               bpq(i1,i2,i) = 0.0D0
            ENDDO
         ENDDO
      ENDDO
!
      DO i = 1 , KNODE
         bpq(1,1,NDOF*(i-1)+1) = pn(i)
         bpq(1,2,NDOF*(i-1)+1) = qn(i)
         bpq(2,1,NDOF*(i-1)+NDOF) = pn(i)
         bpq(2,2,NDOF*(i-1)+NDOF) = qn(i)
      ENDDO
!
      DO i = 1 , NDOF
         DO j = 1 , NDOF
            xjac(i,j) = 0.0D0
         ENDDO
      ENDDO
!
      DO i = 1 , KNODE
         xjac(1,1) = xjac(1,1) + pn(i)*xl((i-1)*NDOF+1)
         xjac(1,2) = xjac(1,2) + qn(i)*xl((i-1)*NDOF+1)
         xjac(2,1) = xjac(2,1) + pn(i)*xl((i-1)*NDOF+2)
         xjac(2,2) = xjac(2,2) + qn(i)*xl((i-1)*NDOF+2)
      ENDDO
!
      det = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
!
      xji(1,1) = xjac(2,2)/det
      xji(2,2) = xjac(1,1)/det
      xji(1,2) = -xjac(1,2)/det
      xji(2,1) = -xjac(2,1)/det
!
      DO i = 1 , KNODE*NDOF
         b(1,i) = bpq(1,1,i)*xji(1,1) + bpq(1,2,i)*xji(2,1)
         b(2,i) = bpq(2,1,i)*xji(1,2) + bpq(2,2,i)*xji(2,2)
         b(3,i) = bpq(1,1,i)*xji(1,2) + bpq(1,2,i)*xji(2,2) + bpq(2,1,i)&
     &            *xji(1,1) + bpq(2,2,i)*xji(2,1)
      ENDDO
!
      CALL FE_FIND_RHS(Lmn,xl,Rhs)
      DO i = 1 , MAXSTN
         Ee_out(i) = 0.0D0
         DO j = 1 , KNODE*NDOF
            Ee_out(i) = Ee_out(i) + b(i,j)*xl(j)
         ENDDO
      ENDDO
      END SUBROUTINE FE_STRAIN
!*==fe_find_rhs.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_FIND_RHS(Lmn,Xl,Rhs)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_FIND_RHS485
      INTEGER Lmn
      DOUBLE PRECISION Xl(KNODE*NDOF) , Rhs(*)
      INTEGER j , k
!
      DO j = 1 , KNODE
         DO k = 1 , NDOF
            Xl(NDOF*(j-1)+k) = Rhs((ICOnn(j,Lmn)-1)*NDOF+k)
         ENDDO
      ENDDO
      END SUBROUTINE FE_FIND_RHS
!*==fe_stress.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FE_STRESS(Lmn,Rhs,S_out)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_STRESS503
      INTEGER MAXSTS , MAXSTN
      PARAMETER (MAXSTS=4)
      PARAMETER (MAXSTN=3)
!
      INTEGER Lmn
      DOUBLE PRECISION Rhs(*) , S_out(3)
!
      DOUBLE PRECISION d(MAXSTS,MAXSTN) , ee_out(MAXSTN)
      INTEGER i , j
      CHARACTER*80 error_message
!
!  1<=> s(1,1)
!  2<=> s(2,2)
!  3<=> s(1,2)
!
      IF ( (Lmn<1) .OR. (Lmn>NELm) ) THEN
         PRINT * , Lmn
         error_message = 'fe_strain: element is out of range!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      CALL FE_DMAT(d)
      CALL FE_STRAIN(Lmn,Rhs,ee_out)
      DO i = 1 , 3
         S_out(i) = 0
         DO j = 1 , MAXSTN
            S_out(i) = S_out(i) + d(i,j)*ee_out(j)
         ENDDO
      ENDDO
!
      END SUBROUTINE FE_STRESS
!*==fe_locate.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      FUNCTION FE_LOCATE(R,Iel0)
!
!  Returns the element number corresponding to position r or zero.
!  Breadth first search in FORTRAN
!
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_LOCATE546
      INTEGER FE_LOCATE , Iel0
      DOUBLE PRECISION R(NDOF+1)
!
      INTEGER Q_MAX
      PARAMETER (Q_MAX=512)
      INTEGER queue(Q_MAX)
      INTEGER idummy , iel , ptr_put , ptr_get , i , lmn
      LOGICAL todo(MAXLMN) , FE_IN_TRI
!
      CHARACTER*80 error_message
!
      IF ( I_Flag/=1 ) THEN
         error_message = 'fe_locate: call fem_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      iel = Iel0
      IF ( ABS(iel)>NELm ) THEN
         PRINT * , 'iel: ' , iel , ' nelm: ' , NELm
         error_message = 'fe_locate: initial element out of boundary'
         CALL ERROR_HANDLER(error_message)
      ENDIF
      IF ( iel==0 ) iel = 1
      IF ( iel<0 ) iel = -iel
!
!  Checking iel in order not to waste time on array initialization
!
      IF ( FE_IN_TRI(X0(1,ICOnn(1,iel)),X0(1,ICOnn(2,iel)),&
     &     X0(1,ICOnn(3,iel)),R) ) THEN
         FE_LOCATE = iel
         RETURN
      ENDIF
!
! initialize the fringe
!
      DO i = 1 , NELm
         todo(i) = .TRUE.
      ENDDO
!
! Breadth-first search checking first iel again to simplify the code
! Circular queue of size Q_MAX
! ptr_put pointer to where to put an element
! ptr_get where the last element was taken (i.e., needs to be incremente
! before taking something out of the queue).
!
      ptr_get = Q_MAX
      ptr_put = 2
      queue(1) = iel
      todo(iel) = .FALSE.
!
      DO idummy = 1 , NELm
         ptr_get = ptr_get + 1
         IF ( ptr_get==Q_MAX+1 ) ptr_get = 1
!  Empty queue. The search over the connected component is exhasted
!  Probably we are done
         IF ( ptr_get==ptr_put ) EXIT
!
         iel = queue(ptr_get)
         IF ( FE_IN_TRI(X0(1,ICOnn(1,iel)),X0(1,ICOnn(2,iel)),&
     &        X0(1,ICOnn(3,iel)),R) ) THEN
            FE_LOCATE = iel
            RETURN
         ENDIF
         DO i = 1 , 3
            lmn = IADj(i,iel)
            IF ( (lmn/=0) .AND. todo(lmn) ) THEN
               IF ( ptr_put==ptr_get ) THEN
                  PRINT * , 'Q_MAX: ' , Q_MAX
                  error_message = 'fe_locate: queue overflow'
                  CALL ERROR_HANDLER(error_message)
               ENDIF
               queue(ptr_put) = lmn
               todo(lmn) = .FALSE.
               ptr_put = ptr_put + 1
               IF ( ptr_put==Q_MAX+1 ) ptr_put = 1
            ENDIF
         ENDDO
      ENDDO
!
      FE_LOCATE = 0
      END FUNCTION FE_LOCATE
!*==fe_in_tri.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      FUNCTION FE_IN_TRI(R1,R2,R3,P)
!
!  Wrapper for intri from feaplib
!
      IMPLICIT NONE
!*--FE_IN_TRI637
      DOUBLE PRECISION R1(3) , R2(3) , R3(3) , P(3)
      LOGICAL FE_IN_TRI
!
      DOUBLE PRECISION x1(2) , x2(2) , x3(2) , pt(2) , s(3)
      LOGICAL INTRI , ontri
      INTEGER i
!
      DO i = 1 , 2
         x1(i) = R1(i)
         x2(i) = R2(i)
         x3(i) = R3(i)
         pt(i) = P(i)
      ENDDO
      FE_IN_TRI = INTRI(x1,x2,x3,pt,s,ontri)
      END FUNCTION FE_IN_TRI
 
 
!
! $Log: fem_alan.f,v $
! Revision 1.2  2004/04/2114:29:39  shastry
! vijay-    mod_fem_paramters fem_alan.f: increased storage.
!
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.7  2001/12/1307:31:24  shilkrot
! Implemented breadth first search to find the element number for
! a dislocation. Changed the interface of fe_locate to use the starting
! element for the search. Old fe_locate is in fem_services.
! Changed the interface of fem_setup. Now two arrays used as temp space 
! passed from outside as the last two parameters.
!
! Revision 1.6  2001/11/1303:32:44  shilkrot
! Added functions computing strain and stress for a given element.
! Wrote a function fe_locate which must be replaced in the future.
!
! Revision 1.5  2001/08/2203:18:35  shilkrot
! Fixed the expression for the energy and polished fem_alan a little bit
! This wersion works with dislocation passing.
!
! Revision 1.4  2001/07/1206:36:49  shilkrot
! Elastic constants are now passed to fem_setup and further to fe_elasti
!
! Revision 1.3  2001/06/2518:55:53  shilkrot
! Changed stiffness matrix for the transposed one
! Otherwise, no major changes
!
! Revision 1.2  2001/06/1801:19:17  shilkrot
! a_stiff becomes ad_stiff and vice versa
!
! Revision 1.1  2001/06/1800:17:25  shilkrot
! The routines initially written by A. Needleman
!
!
