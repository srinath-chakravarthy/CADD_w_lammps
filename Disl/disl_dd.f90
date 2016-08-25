!*==disl_setup.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!     $Id: disl_dd.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE DISL_SETUP()
      USE MOD_DD_SLIP
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--DISL_SETUP9
      CHARACTER*80 error_message
      I_Disl = 0
      irmdisl = .false.
!
!      DATA I_Disl/0/
!
      IF ( I_Disl==0 ) THEN
         I_Disl = 1
!          CALL GEN_SOURCES
!          CALL GEN_OBSTACLES
!!$         CALL PLOT_SLIP
!$$$  stop
!     call gen_slip_planes
      ELSE
         error_message = 'Second call to disl_setup'
         CALL ERROR_HANDLER(error_message)
      ENDIF
      NDIsl = 0
      END SUBROUTINE DISL_SETUP
!*==disl_accept.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
      SUBROUTINE DISL_REMOVE(idisl)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
      integer :: idisl, j
      !!> Shift all arrays using idisl position 
      !!> Arays to update are burgers, burg_length, theta_e, theta_s
      !!>   r_disl, pk_stress, pk_force, disl_range, r_old, disl_residence
      do j = idisl + 1, ndisl
         burgers(1:2,j-1) = burgers(1:2, j)
         BURg_length(j-1) = burg_length(j)
         theta_e(j-1) = theta_e(j)
         theta_s(j-1) = theta_s(j)
         r_disl(:,j-1) = r_disl(:,j)
         pk_stress(:,j-1) = pk_stress(:,j)
         pk_force(:,j-1) = pk_force(:,j)
         disl_range(:,j-1) = disl_range(:,j)
         r_old(:,j-1) = r_old(:,j)
         disl_residence(:,:,j-1) = disl_residence(:,:,j)
         irmdisl(j-1) = irmdisl(j)
         disl_timer(j-1) = disl_timer(j)
      end do
      ndisl = ndisl -1
!
      END SUBROUTINE DISL_REMOVE
 
      SUBROUTINE DISL_ACCEPT(R0,Bv,Th_e,Th_s)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--DISL_ACCEPT33
      DOUBLE PRECISION R0(3) , Bv(3) , Th_e , Th_s , rr , RTOL
      LOGICAL :: image_flag
      PARAMETER (RTOL=1.D0)
      INTEGER ishift
      LOGICAL moved
!
      INTEGER N_Total
      COMMON N_Total
      DATA N_Total/0/
!
      CHARACTER*80 error_message
      INTEGER FE_LOCATE , i , j
!
      IF ( I_Disl/=1 ) THEN
         error_message = 'disl_accept: call disl_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      N_Total = N_Total + 1
!
!     forced shut-off of lumping:
!
      N_Total = 1
!
      IF ( (N_Total==1) .OR. (N_Total==2) ) THEN
         NDIsl = NDIsl + 1
         PRINT * , '***************  Holding ' , NDIsl , &
     &         ' dislocations *****************'
 
         IF ( NDIsl>MAX_DISL ) THEN
            error_message = 'max_disl needs to be increased'
            CALL ERROR_HANDLER(error_message)
         ENDIF
!     if(ndisl.ge.2) then
!     stop 'Amit nucleation stop'
!     endif
         ELEm_disl(NDIsl) = FE_LOCATE(R0,1)
!
 
	    
         DO i = 1 , 3
	    if (R0(2) > -2.0d0) then
		Elem_disl(Ndisl) = 0
	    end if	
            R_Disl(i,NDIsl) = R0(i)
         ENDDO
         IF ( ELEm_disl(NDIsl)>0 ) THEN
            moved = .FALSE.
            ishift = -NINT(COS(Th_e))
            DO
               DO i = 1 , NDIsl - 1
                  rr = 0.D0
                  DO j = 1 , 2
                     rr = rr + (R_Disl(j,i)-R0(j))**2
                  ENDDO
                  IF ( rr<RTOL ) THEN
                     R0(1) = R0(1) + ishift*Bv(1)
                     R0(2) = R0(2) + ishift*Bv(2)
                     moved = .TRUE.
                     GOTO 20
                  ENDIF
               ENDDO
               IF ( moved ) THEN
                  ELEm_disl(NDIsl) = FE_LOCATE(R0,ELEm_disl(NDIsl))
                  WRITE (*,*) NDIsl , ' overlap found.  moved from'
                  WRITE (*,*) R_Disl(1:3,NDIsl) , 'to'
                  WRITE (*,*) R0(1:3)
               ENDIF
               EXIT
 20         ENDDO
         ENDIF
!
         DO i = 1 , 3
            BURgers(i,NDIsl) = Bv(i)
            R_Disl(i,NDIsl) = R0(i)
         ENDDO
         BURg_length(NDIsl) = DSQRT(BURgers(1,NDIsl)*BURgers(1,NDIsl)&
     &                        +BURgers(2,NDIsl)*BURgers(2,NDIsl))
         THEta_e(NDIsl) = Th_e
         THEta_s(NDIsl) = Th_s
         irmdisl(NDisl) = .false. 
         disl_timer(NDISL) = 0
         PRINT * , 'In element: ' , ELEm_disl(NDIsl)
         IF ( ELEm_disl(NDIsl)/=0 ) then 
            CALL SLIPRANGE(Bv,R0,DISl_range(1,NDIsl),DISl_index(NDIsl), disl_residence(1:2,1:2,ndisl))
         END IF
!
      ELSE
         PRINT * , 'Lumping with disl # ' , NDIsl + N_Total - 4
         DO i = 1 , 3
            BURgers(i,NDIsl+N_Total-4) = BURgers(i,NDIsl+N_Total-4)&
     &         + Bv(i)
         ENDDO
         BURg_length(NDIsl+N_Total-4) = SQRT(BURgers(1,NDIsl+N_Total-4)&
     &                                  **2+BURgers(2,NDIsl+N_Total-4)&
     &                                  **2)
         IF ( N_Total==4 ) N_Total = N_Total - 4
      ENDIF
      END SUBROUTINE DISL_ACCEPT
!*==disl_pass.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE DISL_PASS(R0,Rd,Burg,Th_e,Th_s,X,B,Is_relaxed,Numnp, Subtract,Store)
 
!
!     given location x and xd and the burg/theta of a dislocation,
!     compute
!     and return b=b-utilde(x)+utilde(xd)
!     and add this new disl to the d.d. side.
!
!      USE MOD_DD_SLIP
      IMPLICIT NONE
!*--DISL_PASS141
      INTEGER :: disl_num !> @ number of dislocation accepted
      DOUBLE PRECISION R0(3) , Rd(3) , Burg(3) , B(3,*) , X(3,*)
      DOUBLE PRECISION Th_e , Th_s
      INTEGER Is_relaxed(numnp) , Numnp
      LOGICAL Subtract , Store , nucl, image_flag
!
      DOUBLE PRECISION u(3) , ud(3)
      INTEGER i , j
!!$      INTEGER , OPTIONAL :: Idis_slip
!!$      INTEGER , OPTIONAL :: Islp
!!$      DOUBLE PRECISION , OPTIONAL :: S_dis
!

!!$      PRINT * , 'Image Locations'
!!$      PRINT * , 'Subtract =' , R0(1:2)
!!$      PRINT * , 'Image = ' , Rd(1:2)
	if (store) then 
	    call disl_accept(rd, burg, th_e, th_s, disl_num, image_flag)
	end if
!!$      IF ( Store ) THEN
!!$         IF ( nucl ) THEN
!!$            CALL NUCLEATE_ATOMISTIC_DISL(Islp,S_dis,Th_e,Th_s,Burg)
!!$         ELSE
!!$            PRINT * , 'Slip plane num and distance not given to pass'
!!$         ENDIF
!!$      ENDIF
      DO i = 1 , Numnp
         IF ( Is_relaxed(i)/=0 ) THEN
            IF ( Subtract ) THEN
               CALL DISL_U(R0,Burg,Th_e,Th_s,X(1,i),u)
            ELSE
               DO j = 1 , 3
                  u(j) = 0.D0
               ENDDO
            ENDIF
            CALL DISL_U(Rd,Burg,Th_e,Th_s,X(1,i),ud)
            DO j = 1 , 3
               B(j,i) = B(j,i) - u(j) + ud(j)
            ENDDO
         ENDIF
      ENDDO
!     if (.not. store) then
!     endif
      END SUBROUTINE DISL_PASS
!*==fd_no_disl.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      FUNCTION FD_NO_DISL()
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FD_NO_DISL198
      LOGICAL FD_NO_DISL
!
      FD_NO_DISL = (NDIsl==0)
      END FUNCTION FD_NO_DISL
!*==fd_peach_koeller_nuc.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oc
 
      SUBROUTINE FD_PEACH_KOELLER_NUC(Rhs)
      USE MOD_DD_SLIP
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FD_PEACH_KOELLER_NUC209
      INTEGER islp , i , j , k , l , jslp , li , ii
      DOUBLE PRECISION Rhs(*) , sout(3) , s(2)
      PRINT * , 'Cos and Sin and Burger in fd_nuc' , BLEn
 
      ii = 0
      DO islp = 1 , NSLp
         li = LOCphi(islp)
         IF ( NNUc(islp)>0 ) THEN
            DO i = 1 , NNUc(islp)
               ii = ii + 1
               sout = 0.0D0
               s = 0.0D0
               IF ( ELEm_source(i,islp)>0 )&
     &              CALL FE_STRESS(ELEm_source(i,islp),Rhs,sout)
!$$$  do jslp = i, nslp
!$$$  if (nnuc(jslp) > 0) then
!$$$  if (ndis(jslp) > 0) then
!$$$  do j = 1, ndis(islp)
!$$$  call disl_s(r_disl(1,j), burgers(1,j)
!$$$  $                            r_nuc(3,i,islp), burgers
!$$$  end do
!$$$  endif
!$$$  endif
!$$$  end do
               s(2) = s(1) + sout(1)*BLEn*COSphi(li) + sout(3)&
     &                *BLEn*SINphi(li)
               s(1) = s(2) + sout(3)*BLEn*COSphi(li) + sout(2)&
     &                *BLEn*SINphi(li)
               TAUi(i,islp) = (s(1)*COSphi(li)+s(2)*SINphi(li))/BLEn
               WRITE (*,FMT='(3I7,1X,8(1X,E15.7))') ii , i , islp , &
     &                sout(1) , sout(2) , sout(3) , TAUi(i,islp) , &
     &                T_Fr(i,islp) , TNLaps(i,islp)
 
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE FD_PEACH_KOELLER_NUC
!*==fd_peach_koeller.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
 
      SUBROUTINE FD_PEACH_KOELLER(Rhs)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FD_PEACH_KOELLER252
      DOUBLE PRECISION Rhs(*)
!
      INTEGER i , j , k
      DOUBLE PRECISION s_out(3)
      CHARACTER*80 error_message
!
      IF ( I_Disl/=1 ) THEN
         error_message = 'fd_peach_koeller: call disl_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) THEN
!!$	    call print_element(ELEm_disl(i), rhs)
            CALL FE_STRESS(ELEm_disl(i),Rhs,PK_stress(1,i))
!!$            write(*,'(A,I7,6E15.6)') 'Stress on disl = ', i, pk_stress(1:3,i), pk_stress(1:3,i)/1.602176/1.d-5
            DO j = 1 , NDIsl
               IF ( j/=i ) THEN
                  if (elem_disl(j) > 0) then 
!!$ 		     Only include contribution of actual dislocations that can move 
!!$		     if (disl_timer(j) >= time_static) then 
			CALL DISL_S(R_Disl(1,j),BURgers(1,j),R_Disl(1,i),s_out,THEta_s(j))
			DO k = 1 , 3
			    PK_stress(k,i) = PK_stress(k,i) + s_out(k)
			ENDDO
!!$		    end if
		  end if
               ENDIF
            ENDDO
!!$            write(*,'(A,I7,6E15.6)') 'Stress on disl total = ', i, pk_stress(1:3,i), pk_stress(1:3,i)/1.602176/1.d-5
         ENDIF
      ENDDO
 
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) THEN
            PK_force(2,i) = -(PK_stress(1,i)*BURgers(1,i)+PK_stress(3,i)&
     &                      *BURgers(2,i))
            PK_force(1,i) = PK_stress(3,i)*BURgers(1,i) + PK_stress(2,i)&
     &                      *BURgers(2,i)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     if(.not.NumericalPK) then
            PK_f(i) = PK_force(1,i)*BURgers(1,i) + PK_force(2,i)&
     &                *BURgers(2,i)
            PK_f(i) = PK_f(i)/BURg_length(i)
            PRINT * , 'Pk-force =' , i , PK_f(i)
!     if(abs(pk_f(i)).lt.PEIERLS) pk_f(i)=0.d0
!     endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ELSE
            PK_force(1,i) = 0.0D0
            PK_force(2,i) = 0.0D0
            PK_f(i) = 0.0D0
         ENDIF
      ENDDO
!
      END SUBROUTINE FD_PEACH_KOELLER
 
 
 
