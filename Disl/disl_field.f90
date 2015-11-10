!*==disl_u.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!     $Id: disl_field.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE DISL_U(R0,Burgers,Th_e,Th_s,R,U_out)
!
!     calculate displacement field for dislocation using elastic
!     anisotropy
!
!     1234567890123456789012345678901234567890123456789012345678901234
!     1         2         3         4         5         6        7
!
      IMPLICIT NONE
!*--DISL_U14
      DOUBLE PRECISION PI , TOL
      PARAMETER (PI=3.14159265358979323844D0)
      PARAMETER (TOL=1.0D-3)
      INTEGER i , j
      LOGICAL :: aniso = .FALSE.
      LOGICAL , SAVE :: calledf = .FALSE.
      COMPLEX*16 , SAVE :: ad(3,3) , p(3)
      COMPLEX*16 eta
      DOUBLE PRECISION R0(3) , Burgers(3) , Th_e , Th_s , R(3) , &
     &                 U_out(3) , u_temp(3) , factor
!
      DOUBLE PRECISION CC(6,6)
      DOUBLE PRECISION XE , XNU , XLAmbda , XMU
      INTEGER I_Elas
!
      COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
!
      CHARACTER*80 error_message
      DOUBLE PRECISION b2 , dx , dy , r2 , b_ort_x , b_ort_y , x , y , &
     &                 theta , ux , uy , dxr , dyr
!
      IF ( I_Elas/=1 ) THEN
         error_message = 'dislocations: call fe_elastic first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      b2 = Burgers(1)*Burgers(1) + Burgers(2)*Burgers(2)
      dx = R(1) - R0(1)
      dy = R(2) - R0(2)
      r2 = dx*dx + dy*dy
      IF ( r2<TOL ) THEN
         U_out(1:3) = 0.D0
         RETURN
      ENDIF
!
!
      IF ( aniso==.FALSE. ) THEN
!
         b_ort_x = -Burgers(2)
         b_ort_y = Burgers(1)
         x = dx*Burgers(1) + dy*Burgers(2)
         y = dx*b_ort_x + dy*b_ort_y
!
         theta = DATAN2(y,x)
         IF ( theta<=Th_e ) THEN
            theta = theta + (PI-Th_e)
         ELSE
            theta = theta - (PI+Th_e)
         ENDIF
!
         ux = (theta+x*y/b2/2.0D0/(1.0D0-XNU)/r2)
         uy = -((1.0D0-2.0D0*XNU)/4.0D0/(1.0D0-XNU)*DLOG(r2)+(x*x-y*y)&
     &        /b2/4.0D0/(1.0D0-XNU)/r2)
!
         U_out(1) = (ux*Burgers(1)-uy*Burgers(2))/2.0D0/PI
         U_out(2) = (ux*Burgers(2)+uy*Burgers(1))/2.0D0/PI
!
         dxr = -COS(Th_s)*dx - SIN(Th_s)*dy
         dyr = SIN(Th_s)*dx - COS(Th_s)*dy
         U_out(3) = Burgers(3)/2.0D0/PI*DATAN2(dyr,dxr)
!
!
      ELSEIF ( aniso ) THEN
!
         IF ( calledf==.FALSE. ) THEN
            calledf = .TRUE.
            CALL READEIGFILE(p,ad)
         ENDIF
!
!
         IF ( Burgers(1)*COS(Th_s)+Burgers(2)*SIN(Th_s)<0.0 ) THEN
            factor = -1.0
         ELSE
            factor = 1.0
         ENDIF
!
         x = COS(Th_s)*dx + SIN(Th_s)*dy
         y = -SIN(Th_s)*dx + COS(Th_s)*dy
!
         U_out(1:3) = 0.0
         DO i = 1 , 3
            DO j = 1 , 3
               eta = x + p(j)*y
               U_out(i) = U_out(i) - (0.5D0/PI)&
     &                    *AIMAG(factor*ad(i,j)*(LOG(-eta)))
            ENDDO
         ENDDO
!
         u_temp(1) = COS(Th_s)*U_out(1) - SIN(Th_s)*U_out(2)
         u_temp(2) = SIN(Th_s)*U_out(1) + COS(Th_s)*U_out(2)
         U_out(1) = u_temp(1)
         U_out(2) = u_temp(2)
!
!
      ENDIF
!
      END SUBROUTINE DISL_U
!*==readeigfile.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!
!
!
!
!
      SUBROUTINE READEIGFILE(P,Ad)
      IMPLICIT NONE
!*--READEIGFILE121
      COMPLEX*16 P(3) , Ad(3,3)
      INTEGER i , j
      DOUBLE PRECISION aa1 , aa2
      OPEN (UNIT=12,FILE='anisoDis.inp',STATUS='old')
      DO i = 1 , 3
         READ (12,*) aa1 , aa2
         P(i) = DCMPLX(aa1,aa2)
!     write(*,*) p(i)
      ENDDO
      DO i = 1 , 3
         DO j = 1 , 3
            READ (12,*) aa1 , aa2
            Ad(i,j) = DCMPLX(aa1,aa2)
!     write(*,*) ad(i,j)
         ENDDO
      ENDDO
      CLOSE (12)
      END SUBROUTINE READEIGFILE
!*==readeigfile2.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!
!
!
      SUBROUTINE READEIGFILE2(P,Sfact)
      IMPLICIT NONE
!*--READEIGFILE2147
      COMPLEX*16 P(3) , ad(3,3) , Sfact(3,3,3)
      INTEGER i , j , n
      DOUBLE PRECISION aa1 , aa2
      OPEN (UNIT=12,FILE='anisoDis.inp',STATUS='old')
      DO i = 1 , 3
         READ (12,*) aa1 , aa2
         P(i) = DCMPLX(aa1,aa2)
!     write(*,*) p(i)
      ENDDO
      DO i = 1 , 3
         DO j = 1 , 3
            READ (12,*) aa1 , aa2
            ad(i,j) = DCMPLX(aa1,aa2)
!     write(*,*) ad(i,j)
         ENDDO
      ENDDO
      DO i = 1 , 3
         DO j = 1 , 3
            DO n = 1 , 3
               READ (12,*) aa1 , aa2
               Sfact(i,j,n) = DCMPLX(aa1,aa2)
            ENDDO
         ENDDO
      ENDDO
      CLOSE (12)
      END SUBROUTINE READEIGFILE2
!*==disl_s.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
 
!
!
!
      SUBROUTINE DISL_S(R0,Burgers,R,S_out,Th_s)
!
!     2D Stress tensor
!
      IMPLICIT NONE
!*--DISL_S185
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323844D0)
      DOUBLE PRECISION R0(3) , Burgers(3) , R(3) , S_out(3)
      DOUBLE PRECISION CC(6,6)
      DOUBLE PRECISION XE , XNU , XLAmbda , XMU
      INTEGER I_Elas
      COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
      DOUBLE PRECISION dx , dy , dxr , dyr , r2 , factor , vp , sp
      CHARACTER*80 error_message
      INTEGER i , j , n , v
      LOGICAL :: aniso = .FALSE.
      LOGICAL , SAVE :: calledf = .FALSE.
      COMPLEX*16 , SAVE :: sfact(3,3,3) , p(3)
      COMPLEX*16 eta , dw
      DOUBLE PRECISION Th_s , s_temp(3) , rcore , b2 , b1 , r1
!
!     1<=> s(1,1)
!     2<=> s(2,2)
!     3<=> s(1,2)
!
      IF ( I_Elas/=1 ) THEN
         error_message = 'dislocations: call fe_elastic first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
 
      b2 = Burgers(2)*Burgers(2) + Burgers(1)*Burgers(1)
      b1 = DSQRT(b2)
 
      rcore = 2.D0*b1
 
      dx = R(1) - R0(1)
      dy = R(2) - R0(2)
      r2 = dx*dx + dy*dy
!     Cutoff dislocation stresses at 2*burgers magnitude
!     *************************************************
!     This is not strictly correct, at the pad interface
!     but this is ok until the template is put in
      IF ( r2<rcore**2 ) THEN
         r1 = DSQRT(r2)
         dx = dx/r1*rcore
         dy = dy/r1*rcore
         r2 = rcore**2
      ENDIF
!     *************************************************
!
!
      IF ( aniso==.FALSE. ) THEN
!
!
         factor = XMU/PI/2.0D0/(1.0D0-XNU)
!
         vp = Burgers(2)*dx - Burgers(1)*dy
         sp = Burgers(1)*dx + Burgers(2)*dy
!
         S_out(1) = (vp-2.0D0*dx*dy*sp/r2)/r2*factor
         S_out(2) = (vp+2.0D0*dx*dy*sp/r2)/r2*factor
         S_out(3) = sp*(dx*dx-dy*dy)/r2/r2*factor
!
      ELSEIF ( aniso ) THEN
!
         IF ( calledf==.FALSE. ) THEN
            calledf = .TRUE.
            CALL READEIGFILE2(p,sfact)
         ENDIF
!
         IF ( Burgers(1)*COS(Th_s)+Burgers(2)*SIN(Th_s)<0.0 ) THEN
            factor = -1.0
         ELSE
            factor = 1.0
         ENDIF
!
         S_out(1:3) = 0.D0
         dxr = COS(Th_s)*dx + SIN(Th_s)*dy
         dyr = -SIN(Th_s)*dx + COS(Th_s)*dy
!
         DO v = 1 , 3
            IF ( v==1 ) THEN
               i = 1
               j = 1
            ELSEIF ( v==2 ) THEN
               i = 2
               j = 2
            ELSEIF ( v==3 ) THEN
               i = 1
               j = 2
            ENDIF
            DO n = 1 , 3
               eta = dxr + p(n)*dyr
               S_out(v) = S_out(v) - (0.5/PI)&
     &                    *AIMAG(factor*sfact(i,j,n)/eta)
            ENDDO
         ENDDO
!
         s_temp(1) = COS(Th_s)*COS(Th_s)*S_out(1) + SIN(Th_s)*SIN(Th_s)&
     &               *S_out(2) - 2*SIN(Th_s)*COS(Th_s)*S_out(3)
         s_temp(2) = SIN(Th_s)*SIN(Th_s)*S_out(1) + COS(Th_s)*COS(Th_s)&
     &               *S_out(2) + 2*SIN(Th_s)*COS(Th_s)*S_out(3)
         s_temp(3) = COS(Th_s)*SIN(Th_s)*(S_out(1)-S_out(2))&
     &               + (COS(Th_s)*COS(Th_s)-SIN(Th_s)*SIN(Th_s))&
     &               *S_out(3)
!
         S_out(1:3) = s_temp(1:3)
!
      ENDIF
      END SUBROUTINE DISL_S
!*==disl_displ.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!
!
!
!
!
      SUBROUTINE DISL_DISPL(R,U_out)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--DISL_DISPL301
      DOUBLE PRECISION R(3) , U_out(3) , u11(3)
!
      CHARACTER*80 error_message
      DOUBLE PRECISION u_temp(3)
      INTEGER i , j
!
      IF ( I_Disl/=1 ) THEN
         error_message = 'disl_displ: call disl_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , 3
         U_out(i) = 0.0D0
      ENDDO
!
!     print *, 'Disl displ' , ndisl
      DO i = 1 , NDIsl
!$$$         if (elem_disl(i) .eq. 0) then
!$$$            print *, "Calculating bucket dislocation displacment"
!$$$            print *, elem_disl(i), r_disl(1:2,i)
!$$$         end if
         CALL DISL_U(R_Disl(1,i),BURgers(1,i),THEta_e(i),THEta_s(i),R,&
     &               u_temp)
         DO j = 1 , 3
            U_out(j) = U_out(j) + u_temp(j)
         ENDDO
      ENDDO
      END SUBROUTINE DISL_DISPL
!*==disl_stress.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!
!
      SUBROUTINE DISL_STRESS(R,S_out)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--DISL_STRESS337
      DOUBLE PRECISION R(3) , S_out(3)
!
      DOUBLE PRECISION s_temp(3)
      CHARACTER*80 error_message
      INTEGER i , j
!
      IF ( I_Disl/=1 ) THEN
         error_message = 'disl_stress: call disl_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , 3
         S_out(i) = 0.0D0
      ENDDO
!
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) THEN
            CALL DISL_S(R_Disl(1,i),BURgers(1,i),R,s_temp,THEta_s(i))
            DO j = 1 , 3
               S_out(j) = S_out(j) + s_temp(j)
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE DISL_STRESS
 
 
