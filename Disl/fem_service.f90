!*==fe_plot.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
! $Id: fem_service.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE FE_PLOT(X,B,Ix,Rhs,Numel,Avevirst)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_PLOT8
!
      DOUBLE PRECISION X(3,*) , B(3,*) , Rhs(*) , Avevirst(3,3,*)
      INTEGER Ix(4,*) , Numel
      DOUBLE PRECISION xl(KNODE*NDOF) , xu(KNODE*NDOF) , s_out(3)
      INTEGER i , j , k
      CHARACTER*80 filename
      CHARACTER*3 cnt
      INTEGER , SAVE :: icount1
      DATA icount1/0/
      WRITE (cnt,'(i3)') icount1
      IF ( icount1<10 ) cnt(1:2) = '00'
      IF ( icount1<100 ) cnt(1:1) = '0'
      filename(1:7) = 'Stress_'
      filename(8:10) = cnt
      filename(11:14) = '.plt'
!
      OPEN (903,FILE=filename)
      WRITE (903,*) 'ZONE N=' , 3*Numel , ' E=' , Numel , &
     &              ' , F=FEPOINT, ET=TRIANGLE'
      DO i = 1 , NELm
         CALL FE_FIND(i,xl)
         CALL FE_FIND_RHS(i,xu,Rhs)
         CALL FE_STRESS(i,Rhs,s_out)
         DO j = 1 , KNODE
            WRITE (903,'(7e16.8)') (xl(NDOF*(j-1)+k),k=1,NDOF) , &
     &                             (xl(NDOF*(j-1)+k)+xu(NDOF*(j-1)+k),&
     &                             k=1,NDOF) , (s_out(k),k=1,3)
         ENDDO
      ENDDO
      DO i = 1 , Numel
         IF ( Ix(4,i)/=0 ) THEN
            DO j = 1 , KNODE
               WRITE (903,'(7e16.8)') (X(k,Ix(j,i)),k=1,NDOF) , &
     &                                (X(k,Ix(j,i))+B(k,Ix(j,i)),k=1,&
     &                                NDOF) , Avevirst(1,1,Ix(j,i)) , &
     &                                Avevirst(2,2,Ix(j,i)) , &
     &                                Avevirst(1,2,Ix(j,i))
            ENDDO
         ENDIF
      ENDDO
      WRITE (903,*)
      DO i = 1 , Numel
         WRITE (903,*) (KNODE*(i-1)+k,k=1,KNODE)
      ENDDO
      CLOSE (903)
      icount1 = icount1 + 1
!
      END SUBROUTINE FE_PLOT
!*==fe_pk_print.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE FE_PK_PRINT()
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FE_PK_PRINT63
      DOUBLE PRECISION pk
      INTEGER i , j
!
      DO i = 1 , NDIsl
         IF ( ELEm_disl(i)>0 ) THEN
            PRINT * , 'Disl: ' , i
            PRINT * , 'Element: ' , i , ELEm_disl(i)
            PRINT * , 'Stress: ' , i , (PK_stress(j,i),j=1,3)
            pk = PK_force(1,i)*BURgers(1,i) + PK_force(2,i)*BURgers(2,i)
            pk = pk/BURg_length(i)
            PRINT * , 'Force: ' , i , (PK_force(j,i),j=1,2) , pk
            PRINT * , 'Gliding force: ' , i , PK_f(i)
         ENDIF
      ENDDO
!
      END SUBROUTINE FE_PK_PRINT
!*==disl_print.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE DISL_PRINT(Iflag)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--DISL_PRINT87
      INTEGER i , j , Iflag
!
      OPEN (UNIT=401,FILE='dislocations.plt',STATUS='unknown')
      WRITE (401,*) 'zone'
      DO i = 1 , NDIsl
         WRITE (401,'(5e15.6,2i5)') (R_Disl(j,i),j=1,2) , &
     &                              (BURgers(j,i),j=1,3) , i , Iflag
      ENDDO
      CALL FLUSH(401)
!
      END SUBROUTINE DISL_PRINT
!*==fd_move.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FD_MOVE(Id,Dr)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FD_MOVE106
      INTEGER Id , FE_LOCATE
      DOUBLE PRECISION Dr(2)
!
      INTEGER elem_old
!
      R_Disl(1,Id) = R_Disl(1,Id) + Dr(1)
      R_Disl(2,Id) = R_Disl(2,Id) + Dr(2)
      elem_old = ELEm_disl(Id)
      ELEm_disl(Id) = FE_LOCATE(R_Disl(1,Id),ELEm_disl(Id))
      IF ( ELEm_disl(Id)==0 ) THEN
         IF ( elem_old>0 ) THEN
            ELEm_disl(Id) = -elem_old
         ELSE
            ELEm_disl(Id) = elem_old
         ENDIF
      ENDIF
      END SUBROUTINE FD_MOVE
!*==fe_stress_check.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
 
 
      SUBROUTINE FE_STRESS_CHECK(Rhs,E0)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_STRESS_CHECK130
      DOUBLE PRECISION Rhs(*)
      DOUBLE PRECISION e(3) , s(3) , E0 , area , FE_AREA
      INTEGER i
!
      E0 = 0.0D0
      DO i = 1 , NELm
         CALL FE_STRAIN(i,Rhs,e)
         CALL FE_STRESS(i,Rhs,s)
         area = FE_AREA(i)
         E0 = E0 + (s(1)*e(1)+s(2)*e(2)+s(3)*e(3))*area
      ENDDO
      E0 = E0/2.D0
!        print *, 'stress_check: ', e0
      END SUBROUTINE FE_STRESS_CHECK
!*==fe_area.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      FUNCTION FE_AREA(Lmn)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_AREA152
      INTEGER Lmn
      DOUBLE PRECISION FE_AREA
!
      DOUBLE PRECISION xl(KNODE*NDOF)
      DOUBLE PRECISION pn(KNODE) , qn(KNODE)
      DOUBLE PRECISION xjac(NDOF,NDOF) , det , area
      INTEGER i , j
!
! Everything is hard-coded, since I don't know where I may need
! the derivatives of shape functions
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
      CALL FE_FIND(Lmn,xl)
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
      FE_AREA = ABS(det)/2.0D0
!
      END FUNCTION FE_AREA
!*==fe_locate_old.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      FUNCTION FE_LOCATE_OLD(R,Idummy)
!
!  Returns the element number corresponding to position r or zero.
!  Needs to be replaced by something smarter than this.
!
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_LOCATE_OLD202
      INTEGER FE_LOCATE_OLD , Idummy
      DOUBLE PRECISION R(3)
!
      INTEGER i , i1 , i2 , i3
      LOGICAL FE_IN_TRI
      CHARACTER*80 error_message
!
      IF ( I_Flag/=1 ) THEN
         error_message = 'fe_locate_old: call fem_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      FE_LOCATE_OLD = 0
      DO i = 1 , NELm
         i1 = ICOnn(1,i)
         i2 = ICOnn(2,i)
         i3 = ICOnn(3,i)
         IF ( FE_IN_TRI(X0(1,i1),X0(1,i2),X0(1,i3),R) ) THEN
            FE_LOCATE_OLD = i
            RETURN
         ENDIF
      ENDDO
      END FUNCTION FE_LOCATE_OLD
 
 
!
! $Log: fem_service.f,v $
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.2  2001/12/1307:31:24  shilkrot
! Implemented breadth first search to find the element number for
! a dislocation. Changed the interface of fe_locate to use the starting
! element for the search. Old fe_locate is in fem_services.
! Changed the interface of fem_setup. Now two arrays used as temp space 
! passed from outside as the last two parameters.
!
! Revision 1.1  2001/11/1304:30:36  shilkrot
! Service routines like plotting, printing etc.
!
