!*==fem_setup.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!
!       $Id: fem_routines.f,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      SUBROUTINE FEM_SETUP(Numnp,Numel,X,Id,Is_relaxed,Ix,Itx,Cc,&
     &                     Inverse_map,Inverse_el_map)
!
!       3456789012345678901234567890123456789012345678901234567890123
!       1         2         3         4         5         6         7
!
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FEM_SETUP13
!
      INTEGER Numnp , Numel
      DOUBLE PRECISION X(3,*)
      INTEGER Id(3,*) , Is_relaxed(*) , Ix(4,*) , Itx(3,*)
      DOUBLE PRECISION Cc(6,6)
      INTEGER Inverse_map(*) , Inverse_el_map(*)
!
      INTEGER FE_LOCATE
!
      INTEGER i , j , ij , k , ier
      CHARACTER*80 error_message
!
      I_Flag = 0
!
      IF ( I_Flag==0 ) THEN
         I_Flag = 1
      ELSE
         error_message = 'Second call to fem_setup'
         CALL ERROR_HANDLER(error_message)
      ENDIF
!
      CALL FE_ELASTIC(Cc)
!
      NNOdes = 0
      NFIxed = 0
      NPAd = 0
      DO i = 1 , Numnp
         IF ( (Is_relaxed(i)==2) .OR. (Is_relaxed(i)==0) ) THEN
            NNOdes = NNOdes + 1
            IF ( NNOdes>MAXGEO ) THEN
               error_message = 'maxgeo needs to be increased (is_relaxed=(2/0))'
               CALL ERROR_HANDLER(error_message)
            ENDIF
            IMAp(NNOdes) = i
            Inverse_map(i) = NNOdes
            DO j = 1 , NDOF
               X0(j,NNOdes) = X(j,i)
            ENDDO
!       z-component (not needed)
            X0(NDOF+1,NNOdes) = X(NDOF+1,i)
!
            DO j = 1 , NDOF
               IF ( (Is_relaxed(i)==2) .OR. (Id(j,i)==1) ) THEN
                  NFIxed = NFIxed + 1
                  IF ( NFIxed>MAXFIXED ) THEN
                     error_message = 'maxfixed too small (id=1 + 3*n_int)'
                     CALL ERROR_HANDLER(error_message)
                  ENDIF
                  IFIxed(NFIxed) = NDOF*(NNOdes-1) + j
                  IFIx_hold(1,NFIxed) = j
                  IFIx_hold(2,NFIxed) = NNOdes
               ENDIF
            ENDDO
         ELSE
            Inverse_map(i) = 0
!
            IF ( Is_relaxed(i)==-1 ) THEN
               NPAd = NPAd + 1
               IF ( NPAd>MAXPAD ) THEN
                  error_message = 'maxpad too small'
                  CALL ERROR_HANDLER(error_message)
               ENDIF
               PADmap(NPAd) = i
            ENDIF
!
         ENDIF
      ENDDO
!
      NELm = 0
      NSEgm = 0
!
      DO i = 1 , Numel
         IF ( Ix(4,i)==0 ) THEN
            IF ( (Inverse_map(Ix(1,i))==0) .OR. (Inverse_map(Ix(2,i))==0) .OR. (Inverse_map(Ix(3,i))==0) ) THEN
               PRINT * , 'Element number: ' , i
               error_message = 'Inconsistency between ix and id in this element'
               CALL ERROR_HANDLER(error_message)
            ENDIF
            NELm = NELm + 1
            IF ( NELm>MAXLMN ) THEN
               error_message = 'maxlmn (els with ix=0) needs to be increased'
               CALL ERROR_HANDLER(error_message)
            ENDIF
            Inverse_el_map(i) = NELm
            DO j = 1 , KNODE
               ICOnn(j,NELm) = Inverse_map(Ix(j,i))
               IADj(j,NELm) = Itx(j,i)
            ENDDO
!
            DO j = 1 , KNODE
               IF ( (Itx(j,i)==0) .OR. (Ix(4,Itx(j,i))/=0) ) THEN
                  NSEgm = NSEgm + 1
                  IF ( NSEgm>MAXSEGM ) THEN
                     error_message = 'maxsegm needs to be increased'
                     CALL ERROR_HANDLER(error_message)
                  ENDIF
                  ISEgm(1,NSEgm) = ICOnn(j,NELm)
                  ij = j + 1
                  IF ( ij>KNODE ) ij = ij - KNODE
                  ISEgm(2,NSEgm) = ICOnn(ij,NELm)
                  CALL LINECOEF(X0(1,ISEgm(1,NSEgm)),X0(1,ISEgm(2,NSEgm)),ABSegm(1,NSEgm))
               ENDIF
            ENDDO
         ELSE
            Inverse_el_map(i) = 0
         ENDIF
      ENDDO
!
      DO i = 1 , NELm
         DO j = 1 , KNODE
            IF ( IADj(j,i)/=0 ) IADj(j,i) = Inverse_el_map(IADj(j,i))
         ENDDO
      ENDDO
!
      IF ( NPAd/=0 ) THEN
         j = FE_LOCATE(X(1,PADmap(1)),0)
         IF ( j==0 ) THEN
            PRINT * , 'Atom #' , PADmap(1)
            PRINT * , (X(j,PADmap(1)),j=1,NDOF)
            error_message = 'fem_setup: Pad atom outside the continuum region'
            CALL ERROR_HANDLER(error_message)
         ENDIF
         PADelmnt(1) = j
         CALL FE_TRICOORD(X0(1,ICOnn(1,j)),X0(1,ICOnn(2,j)),&
     &                    X0(1,ICOnn(3,j)),X(1,PADmap(1)),&
     &                    PADtricoord(1,1))
         DO i = 2 , NPAd
            j = FE_LOCATE(X(1,PADmap(i)),PADelmnt(i-1))
            IF ( j==0 ) THEN
               PRINT * , 'Atom #' , PADmap(i)
               PRINT * , (X(j,PADmap(i)),j=1,NDOF)
               error_message = &
     &                'fem_setup: Pad atom outside the continuum region'
               CALL ERROR_HANDLER(error_message)
            ENDIF
            PADelmnt(i) = j
            CALL FE_TRICOORD(X0(1,ICOnn(1,j)),X0(1,ICOnn(2,j)),&
     &                       X0(1,ICOnn(3,j)),X(1,PADmap(i)),&
     &                       PADtricoord(1,i))
         ENDDO
      ENDIF
!
      CALL FE_MAKEMAT()
      CALL FE_FIXDSP()
      CALL FE_CHOL(ier)
      PRINT * , 'ier: ' , ier
      PRINT *
      IF ( ier/=0 ) THEN
         error_message = 'Matrix inversion failed!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
      END SUBROUTINE FEM_SETUP
!*==fd_solve.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
      SUBROUTINE FD_SOLVE(B,Bc,Prop,Cz,Id,Is_relaxed,Rhs,Forces,E0)
!
!       1234567890123456789012345678901234567890123456789012345
!       1         2         3         4         5         6        7
!
!       if the node is on the interface u=b
!       if the node is on the boundary u=prop*bc
!
!       Returns rhs, forces, e0
!       This subroutine should never be called!!!
!
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FD_SOLVE191
      DOUBLE PRECISION Bc(3,*) , B(3,*)
      INTEGER Id(3,*) , Is_relaxed(*)
      DOUBLE PRECISION Prop , Cz
!
      DOUBLE PRECISION Rhs(*) , Forces(*) , E0
      DOUBLE PRECISION u_bc(MAXFIXED) , f_segm(NDOF,MAXSEGM)
      DOUBLE PRECISION presv(MAXFIXED)
      INTEGER i , j
!
      CALL FD_UPDATE_U_TILDE_BC(u_bc)
      CALL FD_UPDATE_F_TILDE_BC(f_segm)
!
      DO i = 1 , NEQu
         Rhs(i) = 0.0D0
      ENDDO
!
!       body forces and boundary tractions
!
      DO i = 1 , NNOdes
         DO j = 1 , NDOF
!
!       add body forces to all unconstrained d.o.f. including the
!       interface
!       fe_substitute() overwrites interfacial nodes
!
            IF ( Id(j,IMAp(i))==0 ) THEN
               Rhs((i-1)*NDOF+j) = Prop*Bc(j,IMAp(i))/Cz
            ENDIF
         ENDDO
      ENDDO
!
      DO i = 1 , NFIxed
         IF ( Is_relaxed(IMAp(IFIx_hold(2,i)))==2 ) THEN
            presv(i) = B(IFIx_hold(1,i),IMAp(IFIx_hold(2,i)))
         ELSE
            presv(i) = Prop*Bc(IFIx_hold(1,i),IMAp(IFIx_hold(2,i)))
         ENDIF
      ENDDO
!
!       Set up the b.c. for the ^ field
!
      CALL FD_CORRECTIVE_DISPL_BC(presv,u_bc)
      CALL FD_CORRECTIVE_FORCE_BC(Rhs,f_segm)
!
      CALL FE_SUBSTITUTE(Rhs,presv)
      CALL FE_FORCE_ENERGY(Rhs,Forces,E0)
!
!       Add f~*u^ term to the energy and the appropriate forces
!
      CALL FD_AUGMENT_FORCE_ENERGY(Rhs,Forces,f_segm,E0)
!
      END SUBROUTINE FD_SOLVE
!*==fe_force_energy.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 201
 
 
 
      SUBROUTINE FE_FORCE_ENERGY(Rhs,Forces,E0)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--FE_FORCE_ENERGY250
      DOUBLE PRECISION Rhs(*) , Forces(*) , E0
      INTEGER i , jmin , jmax , j
!
      DO i = 1 , NEQu
         Forces(i) = A_Stiff(1,i)*Rhs(i)
         jmax = MBAndw
         IF ( jmax+i-1>NEQu ) jmax = NEQu - i + 1
         DO j = 2 , jmax
            Forces(i) = Forces(i) + A_Stiff(j,i)*Rhs(i+j-1)
         ENDDO
!
         jmin = MBAndw - 1
         IF ( (i+1-MBAndw)<=0 ) jmin = i - 1
         DO j = 1 , jmin
            Forces(i) = Forces(i) + A_Stiff(j+1,i-j)*Rhs(i-j)
         ENDDO
      ENDDO
!
      E0 = 0.0D0
      DO i = 1 , NEQu
         E0 = E0 + Rhs(i)*Forces(i)/2.0D0
      ENDDO
      END SUBROUTINE FE_FORCE_ENERGY
!*==error_handler.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE ERROR_HANDLER(Message)
      IMPLICIT NONE
!*--ERROR_HANDLER279
      CHARACTER*80 Message
      PRINT '(A,A)' , 'Error: ' , Message
      STOP
      END SUBROUTINE ERROR_HANDLER
!*==findslipplane.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**********************************************************************
      SUBROUTINE FINDSLIPPLANE(B,X0,Ifactor,Islp,Th_s,S_dis,Xd,Xi)
      USE MOD_DD_SLIP
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FINDSLIPPLANE290
      DOUBLE PRECISION :: B(3) , r(2) , Xd(3) , xd1(2) , X0(3) , Xi(3)
      INTEGER , OPTIONAL :: Islp
      INTEGER i , Ifactor , li
      LOGICAL :: between , upper
      DOUBLE PRECISION :: bot , TOL , dtol , S_dis , dist , rmin , Th_s
      DOUBLE PRECISION , VOLATILE :: a1 , a2 , a3 , xst , angl , bmag
      DOUBLE PRECISION :: pi
      PARAMETER (TOL=1.E-6)
      dtol = 1.D-3
      pi = 4.D0*ATAN(1.0)
!
!       eqn of line on which disl moves: a1*x+a2*y+a3=0
!
      r(1:2) = X0(1:2)
      IF ( B(1)==0.D0 .AND. r(1)==0.D0 ) THEN
         a1 = 1.D0
         a2 = 0.D0
         a3 = 0.D0
      ELSEIF ( B(2)==0.D0 .AND. r(2)==0.D0 ) THEN
         a1 = 0.D0
         a2 = 1.D0
         a3 = 0.D0
      ELSE
         a1 = B(2)/(B(1)*r(2)-B(2)*r(1))
         a2 = -B(1)/(B(1)*r(2)-B(2)*r(1))
         a3 = 1.D0
      ENDIF
 
!       Now put it on the slip plane 'islp' at the nearest location
!       Move the dislocation 20 burgers vector away from current locatio
!       This is to ensure that the dislocation is not in the Detection
!       Band
!       Checking for dislocation in the detection band
      IF ( Xd(2)>0.0D0 ) THEN
         upper = .TRUE.
      ELSE
         upper = .FALSE.
      ENDIF
 
      angl = DATAN2(B(2),B(1))
      IF ( angl<0.0D0 ) angl = pi + angl
      xst = -a3/a1
      li = 0
      DO i = 1 , 3
         IF ( DABS(angl-PHIslp(i))<dtol ) THEN
            li = i
            EXIT
         ENDIF
      ENDDO
      IF ( li<=0 .OR. li>3 ) THEN
         PRINT * , 'Slip system not found'
         STOP
      ENDIF
      PRINT * , 'Slip system number = ' , li
      rmin = 1.D30
      DO i = 1 , NSLp
         IF ( LOCphi(i)==li ) THEN
            dist = ABS(xst-XSLp(i))
            IF ( dist<DX_slip ) THEN
               PRINT * , 'Slip planes considered' , i , XSLp(i) , xst , &
     &               dist
               IF ( upper .AND. YENdslp(i)>0 ) THEN
                  IF ( dist<rmin ) THEN
                     rmin = dist
                     Islp = i
                  ENDIF
               ENDIF
               IF ( .NOT.upper .AND. YENdslp(i)<0 ) THEN
                  IF ( dist<rmin ) THEN
                     rmin = dist
                     Islp = i
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      PRINT * , 'Nearest slip plane is =' , Islp
 
 
 
      bmag = DSQRT(B(1)**2+B(2)**2)
      S_dis = DSQRT((Xd(1)-XSLp(Islp))**2+(Xd(2)-YSLp(Islp))**2)
      S_dis = SDIs_out1(1,Islp) + DISl_range_tol + 2.D0*Ifactor*bmag
 
 
      IF ( YENdslp(Islp)>0.0D0 ) THEN
         xd1(1) = XSLp(Islp) + S_dis*COSphi(li)
         xd1(2) = YSLp(Islp) + S_dis*SINphi(li)
      ELSE
         IF ( li==1 ) THEN
            xd1(1) = XSLp(Islp) - S_dis*COSphi(li)
         ELSE
            xd1(1) = XSLp(Islp) + ABS(S_dis*COSphi(li))
         ENDIF
         xd1(2) = -(YSLp(Islp)+S_dis*SINphi(li))
      ENDIF
      PRINT * , 'Slip plane coord = ' , Xd(1:2) , xd1
      Xd(1:2) = xd1
      Xi(1) = XSLp(Islp) + SDIs_out(1,Islp)*COSphi(li)
      Xi(2) = YSLp(Islp) + SDIs_out(1,Islp)*SINphi(li)
      Xi(3) = 0.0D0
 
 
      END SUBROUTINE FINDSLIPPLANE
!*==findentrypoint.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**********************************************************************
      SUBROUTINE FINDENTRYPOINT(B,R,Xd,Ifactor)
      USE MOD_DD_SLIP
      USE MOD_FEM_PARAMETERS
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--FINDENTRYPOINT402
      DOUBLE PRECISION :: B(3) , R(2) , Xd(3) , xint(2) , dist2 , &
     &                    rdist , xd1(3)
      INTEGER i , Ifactor
      LOGICAL :: BETWEEN , upper
      DOUBLE PRECISION bot , TOL , dtol , s_dis
      DOUBLE PRECISION , VOLATILE :: a1 , a2 , a3 , xst , angl , bmag
      PARAMETER (TOL=1.E-6)
      INTEGER , VOLATILE :: ism , li
      dtol = 1.D-3
!
!       eqn of line on which disl moves: a1*x+a2*y+a3=0
!
      IF ( B(1)==0.D0 .AND. R(1)==0.D0 ) THEN
         a1 = 1.D0
         a2 = 0.D0
         a3 = 0.D0
      ELSEIF ( B(2)==0.D0 .AND. R(2)==0.D0 ) THEN
         a1 = 0.D0
         a2 = 1.D0
         a3 = 0.D0
      ELSE
         a1 = B(2)/(B(1)*R(2)-B(2)*R(1))
         a2 = -B(1)/(B(1)*R(2)-B(2)*R(1))
         a3 = 1.D0
      ENDIF
!
!       find all intersections between boundary segments and this line,
!       keep the closest one as the entry point
!
      rdist = 1.E30
      DO ism = 1 , NSEgm
         bot = (a1*ABSegm(2,ism)-ABSegm(1,ism)*a2)
         IF ( ABS(bot)>TOL ) THEN
            xint(1) = (ABSegm(3,ism)*a2-ABSegm(2,ism)*a3)/bot
            xint(2) = (-ABSegm(3,ism)*a1+ABSegm(1,ism)*a3)/bot
!
!       found an intersection:
!
            IF ( BETWEEN(X0(1,ISEgm(1,ism)),X0(1,ISEgm(2,ism)),xint) )&
     &           THEN
               dist2 = (R(1)-xint(1))**2 + (R(2)-xint(2))**2
               IF ( dist2<rdist ) THEN
                  rdist = dist2
                  Xd(1:2) = xint(1:2)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!	Entry point in continuum/pad region is found
      Xd(1:2) = Xd(1:2) + 2.0D0*Ifactor*B(1:2)
!
      END SUBROUTINE FINDENTRYPOINT
!*==sliprange.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE SLIPRANGE(B,R,Range,Index)
      USE MOD_FEM_PARAMETERS
      IMPLICIT NONE
!*--SLIPRANGE460
      DOUBLE PRECISION B(3) , R(2) , Range(2) , xint(2) , RANGETOL
      PARAMETER (RANGETOL=0.5)
      INTEGER Index , i
      LOGICAL BETWEEN
      DOUBLE PRECISION a1 , a2 , a3 , bot , TOL
      PARAMETER (TOL=1.E-6)
      INTEGER ism
!
!       eqn of line on which disl moves: a1*x+a2*y+a3=0
!
      IF ( B(1)==0.D0 .AND. R(1)==0.D0 ) THEN
         a1 = 1.D0
         a2 = 0.D0
         a3 = 0.D0
      ELSEIF ( B(2)==0.D0 .AND. R(2)==0.D0 ) THEN
         a1 = 0.D0
         a2 = 1.D0
         a3 = 0.D0
      ELSE
         a1 = B(2)/(B(1)*R(2)-B(2)*R(1))
         a2 = -B(1)/(B(1)*R(2)-B(2)*R(1))
         a3 = 1.D0
      ENDIF
!
!       depending on orientation of b, use x or y for comparisons
!
      IF ( ABS(B(1))>ABS(B(2)) ) THEN
         Index = 1
      ELSE
         Index = 2
      ENDIF
!
!       find all intersections between boundary segments and this line,
!       limit dislocation motion accordingly
!
      Range(1) = -1.E30
      Range(2) = 1.E30
      DO ism = 1 , NSEgm
         bot = (a1*ABSegm(2,ism)-ABSegm(1,ism)*a2)
         IF ( ABS(bot)>TOL ) THEN
            xint(1) = (ABSegm(3,ism)*a2-ABSegm(2,ism)*a3)/bot
            xint(2) = (-ABSegm(3,ism)*a1+ABSegm(1,ism)*a3)/bot
!
!       found an intersection:
!
            IF ( BETWEEN(X0(1,ISEgm(1,ism)),X0(1,ISEgm(2,ism)),xint) )&
     &           THEN
               IF ( R(Index)<xint(Index) ) THEN
                  Range(2) = MIN(xint(Index),Range(2))
               ELSE
                  Range(1) = MAX(xint(Index),Range(1))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!
!       set the closeness to the interface allowed.
!
      Range(2) = Range(2) - RANGETOL
      Range(1) = Range(1) + RANGETOL
      IF ( Index==1 ) THEN
         WRITE (*,*) 'X-range for the dislocation' , Range(1) , Range(2)
      ELSE
         WRITE (*,*) 'Y-range for the dislocation' , Range(1) , Range(2)
      ENDIF
      END SUBROUTINE SLIPRANGE
 
!
!       $Log: fem_routines.f,v $
!       Revision 1.1.1.1  2003/03/1220:09:00  shastry
!       vijay-   Initial import.
!
!       Revision 1.10  2002/06/0420:31:44  shilkrot
!       1. Rewrote fem solve (changed commons in fem_parameters and
!       energy
!       and numerical force computation.
!       2. Introduced negative element # and penalty.
!       3. Add flag MoveDisl to fem_solve.
!
!       Revision 1.9  2001/12/1307:31:24  shilkrot
!       Implemented breadth first search to find the element number for
!       a dislocation. Changed the interface of fe_locate to use the
!       starting
!       element for the search. Old fe_locate is in fem_services.
!       Changed the interface of fem_setup. Now two arrays used as temp
!       space are
!       passed from outside as the last two parameters.
!
!       Revision 1.8  2001/11/1304:09:40  shilkrot
!       Added computation of the P.-K. force.
!
!       Revision 1.7  2001/08/2523:51:12  shilkrot
!       Fixed the bug by dividing input forces by cz, moved fd_update_
!       *_bc
!       form disl_pass here, eliminated the call to
!       fd_restore_displ_bc()
!
!       Revision 1.6  2001/08/2203:18:35  shilkrot
!       Fixed the expression for the energy and polished fem_alan a
!       little bit.
!       This wersion works with dislocation passing.
!
!       Revision 1.5  2001/07/1215:45:29  shilkrot
!       Changed error messages to make it more convenient for Ron.
!
!       Revision 1.4  2001/07/1206:36:49  shilkrot
!       Elastic constants are now passed to fem_setup and further to
!       fe_elastic
!
!       Revision 1.3  2001/07/1206:33:41  shilkrot
!       Updating vector b is moved into the subroutine fd_full_field
!       where
!       the tilde field is added.
!
!       Revision 1.2  2001/06/2520:56:46  shilkrot
!       Added the part correcting for the tilde field
!
!       Revision 1.1  2001/06/1800:22:07  shilkrot
!       Interface routines between QC and FEM to minimize the total
!       energy
!
!
