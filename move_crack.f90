!*==get_crack_tip.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 2015
!**********************************************************************
      SUBROUTINE GET_CRACK_TIP(X,B)
!     Calculate the crack tip position from atomistics
!     Atomistic crack tip position is computed using the Neighbor list
!     of the atoms
!     This is used rather than a coordination number to avoid any
!     other calculations
!     In this instance Surface atoms are identified using the number of
!     the near neighbors < 23
!     Perfect atoms = 26
!     Dislocations > 26
!     The crack tip position is max(x) of the surface atoms
 
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--GET_CRACK_TIP18
!*** Start of declarations inserted by SPAG
      INTEGER i , imin , j , npad , npatoms , ntot
      DOUBLE PRECISION pe , startpad , x1 , y1
!*** End of declarations inserted by SPAG
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*)
      DOUBLE PRECISION box_max(2) , box_min(2)
      DOUBLE PRECISION xtip_min(2)
      xtip_min = -1.E30
      DO j = 1 , 2
         box_max(j) = -1E30
         box_min(j) = 1E30
      ENDDO
      npatoms = 0
      npad = 0
      ntot = 0
!     ********************************************************
!     ---- Calculate the Box size ignoring the interface atoms
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. ISRelaxed(i)&
     &        ==-1 ) THEN
            IF ( ISRelaxed(i)==-1 ) THEN
               IF ( npad==0 ) startpad = i
               npad = npad + 1
            ELSE
               pe = pe + ENErgy(i)
               npatoms = npatoms + 1
            ENDIF
            ntot = ntot + 1
            DO j = 1 , 2
               IF ( ISRelaxed(i)/=-1 ) THEN
                  IF ( X(j,i)>box_max(j) ) box_max(j) = X(j,i)
                  IF ( X(j,i)<box_min(j) ) box_min(j) = X(j,i)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      startpad = npatoms - npad + 1
!     ********************************************************
!     Compute surface atoms and therefore crack tip atom pos
!     Surface atoms are classified as those with less than the
!     26 nearest neighbors
!     ********************************************************
      imin = -1000000
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 ) THEN
            IF ( X(1,i)<box_max(1) .AND. X(1,i)>box_min(1) ) THEN
               IF ( X(2,i)<box_max(2) .AND. X(2,i)>box_min(2) ) THEN
                  IF ( NUMneighbors(i)<=23 ) THEN
                     IF ( X(1,i)>xtip_min(1) ) THEN
                        xtip_min(1) = X(1,i)
                        xtip_min(2) = X(2,i)
                        imin = i
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF ( imin>0 ) THEN
         XTIp(1) = xtip_min(1)
         XTIp(2) = xtip_min(2)
         x1 = (X(1,imin)+box_min(1))/(box_max(1)-box_min(1))
         y1 = (X(2,imin)+box_min(2))/(box_max(2)-box_min(2))
         WRITE (*,'(A20,2I8,5F16.11)') 'Crack Tip Position = ' , imin , &
     &                                 NUMneighbors(imin) , xtip_min(1)&
     &                                 , X(1,imin) , X(2,imin) , x1 , y1
         WRITE (*,'(A20,2F16.11)') 'New Crack Tip = ' , X(1,imin)&
     &                             + B(1,imin) , X(2,imin) + B(2,imin)
 
      ELSE
         PRINT * , 'Crack not in atomistic region **********'
         PRINT * , 'Assuming positio for crack tip'
         xtip_min(1) = -1.E-30
         xtip_min(2) = 0.0D0
      ENDIF
 
!$$$      xtip = xtip_min
      END SUBROUTINE GET_CRACK_TIP
!*==move_atomistic_crack.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oc
!**********************************************************************
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Move crack tip back to (x0,y0)
!     This means moving the AtomDispl from current atom position
!     to the new one
!     1) xtip(1) contains current x position of crack
!     2)
!
      SUBROUTINE MOVE_ATOMISTIC_CRACK(X,Ix,B)
      USE MOD_GLOBAL
      USE MOD_FILE
      USE MOD_DYNAMO
      USE mod_disl_parameters
      IMPLICIT NONE
!*--MOVE_ATOMISTIC_CRACK113
 
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*)
      INTEGER Ix(NEN1,*) , iq
 
      INTEGER iatom , j , i , ifem
 
      INTEGER logic
      INTEGER , ALLOCATABLE :: mmap(:)
      DOUBLE PRECISION :: xnew(NXDm) , xold(NXDm) , rr , box_min(2) , &
     &                    box_max(2) , umag , bb(NXDm)
      DOUBLE PRECISION :: xtip_move
      INTEGER ntot , npatoms
      DOUBLE PRECISION xdef , ydef , zdef
      CHARACTER*80 input , filename , temp
      CHARACTER*4 key
      CHARACTER*6 cnt
      INTEGER numatoms
      DOUBLE PRECISION , ALLOCATABLE :: new_b(:,:)
      INTEGER , ALLOCATABLE :: new_numneigh(:) , new_neighlist(:,:)
      DOUBLE PRECISION , ALLOCATABLE :: utilde(:,:) , new_utilde(:,:)
 
      ALLOCATE (new_b(NDF,NUMnp))
      ALLOCATE (new_numneigh(NUMnp))
      ALLOCATE (new_neighlist(MAXneighbors,NUMnp))
      umag = 1.0
      iq = 0
      filename = "out/before_crack_atoms.cfg"
      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
      CALL DUMP_ATOM(X,B,logic)
      CLOSE (logic)
 
 
      filename = "out/before_mm.vtk"
      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
      CALL DUMP_MESH(X,B,Ix,logic,iq)
      CLOSE (logic)
 
      ALLOCATE (mmap(NUMnp))
      mmap = 0
      PRINT * , 'Crack tip motion = ' , XTIp_actual*X_Tip_dir
      xtip_move = INT(XTIp_actual(1)/X_Move_mesh)*X_Move_mesh
      PRINT * , 'Actual crack tip motion = ' , &
     &      INT(XTIp_actual/X_Move_mesh)*X_Move_mesh
 
      DO i = 1 , 2
         CRAck_motion(i) = CRAck_motion(i) + X_Tip_dir(i)*xtip_move
      ENDDO
!     Calculate the utilde vector everywhere
      ALLOCATE (utilde(3,NUMnp),new_utilde(3,NUMnp))
      utilde = 0.0D0
      new_utilde = 0.0D0
      DO i = 1 , NUMnp
         CALL DISL_DISPL(X(1:3,i),utilde(1:3,i))
      ENDDO
 
!     Create mapping from new coordinates to old coordinates
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)==1 .OR. ISRelaxed(iatom)==2 .OR. &
     &        ISRelaxed(iatom)==-1 ) THEN
            DO j = 1 , 1
!!               xold(j) = x(j,iAtom) + x_tip_dir(j)*xtip_actual(j)
               xold(j) = X(j,iatom) + xtip_move
            ENDDO
            xold(2:NXDm) = X(2:NXDm,iatom)
            DO i = 1 , NUMnp
               IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. &
     &              ISRelaxed(i)==-1 ) THEN
                  xnew(1:NXDm) = X(1:NXDm,i)
!     Calculate distance between 2 atoms
                  rr = 0.D0
                  DO j = 1 , NXDm
                     rr = rr + (xnew(j)-xold(j))**2
                  ENDDO
                  rr = DSQRT(rr)
                  IF ( rr<1.D-3 ) THEN
                     mmap(iatom) = i
                     new_b(:,i) = B(:,i)
                     new_numneigh(i) = NUMneighbors(i)
                     new_neighlist(:,i) = NEIghborlist(:,i)
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
!$$$  write(*,*) 'MMap', iAtom, mmap(iAtom)
      ENDDO
      DO iatom = 1 , NUMnp
         IF ( ISRelaxed(iatom)==1 .OR. ISRelaxed(iatom)==2 .OR. &
     &        ISRelaxed(iatom)==-1 ) THEN
!     Copy displacements from atom i to iAtom to move the crack tip
!     Also copy velocities and neighbour lists ???
            i = mmap(iatom)
!$$$            write(*,'(A10,2I7,3(1X,E15.8),2I7)'),
!$$$     $           'Mapping = ',iAtom, i
!$$$     $           , x(1,iAtom),x(1:nxdm,i), b(1,iAtom), b(1,i)
!$$$               write(*,'(A13,2E15.8,A4,2E15.8, A8, 4(1X,E15.8)),2I7'
!$$$     $              'Mapping from ', x(1:2,i), ' to ', x(1:2,iAtom),
!$$$     $              ' Disp = ', b(1:2,iAtom), b(1:2, i),
!$$$     $              NumNeighbors(iAtom), NumNeighbors(i)
!     *****************************************************************
!     Updates the neighbor list as well
!$$$               NumNeighbors(iAtom) = new_numneigh(i)
!$$$               NeighborList(:,iAtom) = new_neighlist(:,i)
!     *****************************************************************
            IF ( i/=0 ) B(1:NDF,iatom) = new_b(1:NDF,i)
!     End copy displacements
         ENDIF
      ENDDO
 
 
      CALL MOVE_SLIP_PLANES(xtip_move)
 
      DO i = 1 , NUMnp
         CALL DISL_DISPL(X(1:3,i),new_utilde(1:3,i))
      ENDDO
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==0 ) B(1:3,i) = B(1:3,i) - utilde(1:3,i)&
     &        + new_utilde(1:3,i)
      ENDDO
 
      DEALLOCATE (new_b)
      DEALLOCATE (mmap)
      DEALLOCATE (new_numneigh)
      DEALLOCATE (new_neighlist)
 
 
!$$$c     Just use the finite deformation quantities of the atomistics
!$$$c     b(:,i) = b(:,i) - xtip_actual
!$$$c     Since this will be rigid body translation
!$$$c     The deformation gradient within the atomistics will not change
!$$$c     Right now pad atoms are also moved
!$$$c     They should be moved
!$$$      do iAtom = 1, numnp
!$$$         if (isRelaxed(iAtom) .eq. 1 .or. isRelaxed(iAtom) .eq. 2
!$$$     $        .or. isRelaxed(iAtom) .eq. -1) then
!$$$           do j = 1, 1
!$$$              b(j,iAtom) = b(j,iAtom) - x_tip_dir(j)*xtip_actual(j)
!$$$           end do
!$$$         end if
!$$$      end do
 
 
      filename = "out/after_crack_atoms.cfg"
      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
      CALL DUMP_ATOM(X,B,logic)
      CLOSE (logic)
 
      filename = "out/after_mm.vtk"
      CALL IOFILE(filename,'formatted  ',logic,.FALSE.)
      CALL DUMP_MESH(X,B,Ix,logic)
      PRINT * , 'Total Crack tip motion (deltaA) = ' , CRAck_motion
!     Set newmesh parameter to ensure that the DB strains are recalculat
      NEWmesh = .TRUE.
 
      CLOSE (logic)
 
      END SUBROUTINE MOVE_ATOMISTIC_CRACK
!*==move_slip_planes.spg  processed by SPAG 6.70Rc at 12:37 on 29 Oct 20
!**********************************************************************
 
      SUBROUTINE MOVE_SLIP_PLANES(Xtip_move)
      USE MOD_GLOBAL
      USE MOD_DD_SLIP
      IMPLICIT NONE
!*--MOVE_SLIP_PLANES278
      INTEGER islp
      DOUBLE PRECISION :: Xtip_move
      DO islp = 1 , NSLp
         XSLp(islp) = XSLp(islp) - Xtip_move
      ENDDO
      CALL GEN_SLIP_ENDS
 
      CALL ASSIGN_NEW_GLOBAL(Xtip_move)
      END SUBROUTINE MOVE_SLIP_PLANES
