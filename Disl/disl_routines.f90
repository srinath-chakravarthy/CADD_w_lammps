!*==fd_update_u_tilde_bc.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oc
!
! $Id: disl_routines.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE FD_UPDATE_U_TILDE_BC(U_bc)
!
!12345678901234567890123456789012345678901234567890123456789012345678901
!         1         2         3         4         5         6         7
!
      USE MOD_DD_SLIP
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_UPDATE_U_TILDE_BC13
      DOUBLE PRECISION U_bc(*)
!
      DOUBLE PRECISION u_out(3) , u11(3)
      CHARACTER*80 error_message
      INTEGER i
!
      IF ( I_Flag/=1 ) THEN
         error_message = 'fd_update_u_tilde_bc: call fem_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , NFIxed
         CALL DISL_DISPL(X0(1,IFIx_hold(2,i)),u_out)
!          call dislp(x0(1, ifix_hold(2,i)), u11)
!          write(*, fmt='(A6,I7,1X,4(1X,E15.8))') 'dislp ',i,
!     $         u11(1),u_out(1)
!     $         , u11(2), u_out(2)
         U_bc(i) = u_out(IFIx_hold(1,i))
!          print *, 'Bc = ', i, u_bc(i)
      ENDDO
      END SUBROUTINE FD_UPDATE_U_TILDE_BC
!*==fd_update_f_tilde_bc.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oc
 
 
 
      SUBROUTINE FD_UPDATE_F_TILDE_BC(F_segm)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_UPDATE_F_TILDE_BC42
      INTEGER NINP
      PARAMETER (NINP=2)
      DOUBLE PRECISION F_segm(NDOF,*)
!
      DOUBLE PRECISION s(3) , s_out(3) , xnorm(NDOF) , x(NDOF+1)
      DOUBLE PRECISION xinp(NINP) , w(NINP)
      INTEGER i , j , k
      CHARACTER*80 error_message
      DATA xinp/ - 0.577350269189625731D0 , 0.577350269189625731D0/
      DATA w/1.D0 , 1.D0/
!
      IF ( I_Flag/=1 ) THEN
         error_message = 'fd_update_f_tilde_bc: call fem_setup first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      DO i = 1 , NSEgm
         xnorm(1) = X0(2,ISEgm(2,i)) - X0(2,ISEgm(1,i))
         xnorm(2) = -(X0(1,ISEgm(2,i))-X0(1,ISEgm(1,i)))
!
!  Integral over the segment using 2 Gaussian integration points
!
         DO j = 1 , 3
            s(j) = 0.D0
         ENDDO
!
         DO j = 1 , NINP
            DO k = 1 , NDOF
               x(k) = (1.D0+xinp(j))/2.0D0*X0(k,ISEgm(1,i))&
     &                + (1.D0-xinp(j))/2.0D0*X0(k,ISEgm(2,i))
            ENDDO
            x(NDOF+1) = (1.D0+xinp(j))/2.0D0*X0(NDOF+1,ISEgm(1,i))&
     &                  + (1.D0-xinp(j))/2.0D0*X0(NDOF+1,ISEgm(2,i))
!
!            print *, 'Nsegm', i, x(1),x(2)
            CALL DISL_STRESS(x,s_out)
            DO k = 1 , 3
               s(k) = s(k) + s_out(k)*w(j)/2.0D0
            ENDDO
         ENDDO
!
         F_segm(1,i) = s(1)*xnorm(1) + s(3)*xnorm(2)
         F_segm(2,i) = s(2)*xnorm(2) + s(3)*xnorm(1)
      ENDDO
      END SUBROUTINE FD_UPDATE_F_TILDE_BC
!*==fd_corrective_force_bc.spg  processed by SPAG 6.70Rc at 12:39 on 29 
 
 
 
      SUBROUTINE FD_CORRECTIVE_FORCE_BC(Rhs,F_segm)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_CORRECTIVE_FORCE_BC95
      DOUBLE PRECISION Rhs(*) , F_segm(NDOF,*)
!
      INTEGER i , j , k
!
! Puts a force on constrained nodes. This is overwritten by
! fe_substitute()
!
      DO i = 1 , NSEgm
         DO j = 1 , 2
            DO k = 1 , NDOF
               Rhs((ISEgm(j,i)-1)*NDOF+k) = Rhs((ISEgm(j,i)-1)*NDOF+k)&
     &            - F_segm(k,i)/2.0D0
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE FD_CORRECTIVE_FORCE_BC
!*==fd_augment_force_energy.spg  processed by SPAG 6.70Rc at 12:39 on 29
 
 
 
      SUBROUTINE FD_AUGMENT_FORCE_ENERGY(Rhs,Forces,F_segm,E0)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_AUGMENT_FORCE_ENERGY119
      DOUBLE PRECISION Rhs(*) , Forces(*) , F_segm(NDOF,*) , E0
!
      INTEGER i , j , k
!
!  f~ * u^ term
!
      DO i = 1 , NSEgm
         DO j = 1 , 2
            DO k = 1 , NDOF
               Forces((ISEgm(j,i)-1)*NDOF+k)&
     &            = Forces((ISEgm(j,i)-1)*NDOF+k) + F_segm(k,i)/2.0D0
               E0 = E0 + F_segm(k,i)/2.0D0*Rhs((ISEgm(j,i)-1)*NDOF+k)
            ENDDO
         ENDDO
      ENDDO
!
      END SUBROUTINE FD_AUGMENT_FORCE_ENERGY
!*==fd_corrective_displ_bc.spg  processed by SPAG 6.70Rc at 12:39 on 29 
 
 
 
      SUBROUTINE FD_CORRECTIVE_DISPL_BC(Presv,U_bc)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_CORRECTIVE_DISPL_BC144
      DOUBLE PRECISION Presv(*) , U_bc(*)
!
      INTEGER i
!
      DO i = 1 , NFIxed
         Presv(i) = Presv(i) - U_bc(i)
      ENDDO
!
      END SUBROUTINE FD_CORRECTIVE_DISPL_BC
!*==fd_restore_displ_bc.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct
 
 
 
      SUBROUTINE FD_RESTORE_DISPL_BC(Rhs,U_bc)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_RESTORE_DISPL_BC161
      DOUBLE PRECISION Rhs(*) , U_bc(*)
!
      INTEGER i , ic
!
      DO i = 1 , NFIxed
         ic = IFIxed(i)
         Rhs(ic) = Rhs(ic) + U_bc(i)
      ENDDO
      END SUBROUTINE FD_RESTORE_DISPL_BC
!*==fd_full_field.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FD_FULL_FIELD(Rhs,B)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_FULL_FIELD178
      DOUBLE PRECISION Rhs(*) , B(3,*)
!
      DOUBLE PRECISION u(3)
      INTEGER i , j
!
      DO i = 1 , NNOdes
         CALL DISL_DISPL(X0(1,i),u)
         DO j = 1 , NDOF
            Rhs((i-1)*NDOF+j) = Rhs((i-1)*NDOF+j) + u(j)
            B(j,IMAp(i)) = Rhs((i-1)*NDOF+j)
         ENDDO
         B(NDOF+1,IMAp(i)) = u(3)
      ENDDO
      END SUBROUTINE FD_FULL_FIELD
!*==fd_rescale_peach_koeller.spg  processed by SPAG 6.70Rc at 12:39 on 2
 
 
 
      SUBROUTINE FD_RESCALE_PEACH_KOELLER(Cz)
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FD_RESCALE_PEACH_KOELLER200
      DOUBLE PRECISION Cz
      INTEGER i , j
!
      DO i = 1 , NDIsl
         DO j = 1 , 2
            PK_force(j,i) = PK_force(j,i)*Cz
         ENDDO
         PK_f(i) = PK_f(i)*Cz
      ENDDO
!
      END SUBROUTINE FD_RESCALE_PEACH_KOELLER
 
 
!
! $Log: disl_routines.f,v $
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.3  2001/11/1303:34:47  shilkrot
! Added the routines that compute the P.-K. force.
!
! Revision 1.2  2001/08/2203:18:35  shilkrot
! Fixed the expression for the energy and polished fem_alan a little bit
! This wersion works with dislocation passing.
!
! Revision 1.1  2001/07/1206:30:20  shilkrot
! The routines used to apply the tilde field and to handle the array of
! dislocations.
!
!
