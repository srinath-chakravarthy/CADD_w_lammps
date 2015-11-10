!*==pdel.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**   pdel : write a load-delta point to the pdel file
!**
      SUBROUTINE PDEL(F,Dr,Id,Logic,Ndf,X,Time)
      IMPLICIT NONE
!*--PDEL7
      DOUBLE PRECISION F , Dr , X , Time
      INTEGER Id , Logic , Ndf
!
      DOUBLE PRECISION force
      DOUBLE PRECISION TOTener , STRener
      COMMON /ENERDAT/ TOTener , STRener
!
      CALL PDELCALC(F,Dr,Id,force,X)
      WRITE (Logic,'(4e15.5)') Time , force , TOTener , STRener
      CALL FLUSH(Logic)
      END SUBROUTINE PDEL
!*==initialiseenergy.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 20
!**---------------------------------------------------------------------
!** InitialiseEnergy
!**
!--
      SUBROUTINE INITIALISEENERGY(E_flag,F_flag,Id,Dr,F)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--INITIALISEENERGY27
 
!--passed variables
      LOGICAL E_flag , F_flag
      INTEGER Id(NDF,*)
      DOUBLE PRECISION Dr(NDF,*) , F(NDF,*)
!--Local Variables
      INTEGER i , j
 
!     Initialize out-of-balance force vector
      IF ( F_flag ) THEN
         DO i = 1 , NUMnp
            DO j = 1 , NDF
               IF ( Id(j,i)==0 ) THEN
                  Dr(j,i) = TIMe*F(j,i)
               ELSE
                  Dr(j,i) = 0.D0
               ENDIF
            ENDDO
         ENDDO
         IF ( NUMnpp1>NUMnp ) Dr(1:NDF,NUMnpp1) = 0.D0
      ENDIF
 
      IF ( E_flag ) THEN
         DO i = 1 , NUMnp
            ENErgy(i) = 0.0
         ENDDO
      ENDIF
      END SUBROUTINE INITIALISEENERGY
!*==statuscalc.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!**
      SUBROUTINE STATUSCALC(X,Ix,Silent)
      USE MOD_FILE
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--STATUSCALC63
!
      INTEGER Ix(NEN1,NUMel)
      DOUBLE PRECISION X(NXDm,NUMnp) , s(3)
      CHARACTER*80 filename
      INTEGER i , j , k , icheck , icheck2 , logic
      LOGICAL first , Silent , intri , ontri
 
      IF ( ALLOCATED(ISRelaxed) ) DEALLOCATE (ISRelaxed)
      IF ( .NOT.ALLOCATED(ISRelaxed) ) ALLOCATE (ISRelaxed(NUMnp))
!--- only atoms in the atomistic region (abs(ix(nen1,i)).eq.1) are "rela
!--- in CG.
!-- IsRelaxed=0: continuum
!-- IsRelaxed=1: atomistic
!-- IsRelaxed=2: interface
!
      ISRelaxed(1:NUM2dnode) = -1
!
      DO i = 1 , NUM2dnode
         first = .TRUE.
         DO j = 1 , NUMel
            DO k = 1 , 3
               IF ( Ix(k,j)==i ) THEN
                  IF ( first ) THEN
                     icheck = Ix(NEN1,j)
                     IF ( icheck<0 ) icheck = 1
                     first = .FALSE.
                     IF ( Ix(NEN1,j)==0 ) THEN
                        ISRelaxed(i) = 0
                     ELSE
                        ISRelaxed(i) = 1
                     ENDIF
                     EXIT
                  ELSE
                     icheck2 = Ix(NEN1,j)
                     IF ( icheck2<0 ) icheck2 = 1
                     IF ( icheck==icheck2 ) THEN
                        IF ( Ix(NEN1,j)==0 ) THEN
                           ISRelaxed(i) = 0
                        ELSE
                           ISRelaxed(i) = 1
                        ENDIF
                        EXIT
                     ELSE
                        ISRelaxed(i) = 2
                        GOTO 100
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
 100  ENDDO
      IF ( NUMnpp1>NUMnp ) ISRelaxed(NUMnpp1) = 1
 
      IF ( Silent ) RETURN
      NQC = 0
      NSPring = 0
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==1 ) THEN
            NQC = NQC + 1
         ELSEIF ( ISRelaxed(i)==2 ) THEN
            NSPring = NSPring + 1
         ENDIF
      ENDDO
      WRITE (*,*) 'Number of nodes in the atomistic region:' , NQC
      WRITE (*,*) 'Number of nodes in the elastic region:' , &
     &            NUMnp - NQC - NSPring
      WRITE (*,*) 'Number of nodes on the interface:' , NSPring
      IF ( NEN1/=4 ) STOP 'ERROR: nen must be 3!'
      filename = 'out/check.plt'
      OPEN (UNIT=123,FILE=filename,STATUS='unknown')
      WRITE (123,*) 'VARIABLES = X Y Z ISRELAXED'
      WRITE (123,*) 'zone, t=continuum'
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==0 ) WRITE (123,99001) X(1:3,i) , 0
      ENDDO
      WRITE (123,*) 'zone, t=atomistic'
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==1 ) WRITE (123,99001) X(1:3,i) , 1
      ENDDO
      WRITE (123,*) 'zone, t=interface'
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==2 ) WRITE (123,99001) X(1:3,i) , 2
      ENDDO
      WRITE (123,*) 'zone, t=Pad'
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==-1 ) WRITE (123,99001) X(1:3,i) , -1
      ENDDO
      CLOSE (123)
      filename = 'out/check.vtk'
      CALL IOFILE(filename,'formatted  ',logic,.TRUE.)
      WRITE (logic,FMT='(A)') '# vtk DataFile Version 2.0'
      WRITE (logic,FMT='(A)') 'Atoms/Nodes colored by region'
      WRITE (logic,FMT='(A)') 'ASCII'
      WRITE (logic,FMT='(A)') 'DATASET POLYDATA'
      WRITE (logic,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , NUMnp , 'float'
      DO i = 1 , NUMnp
         WRITE (logic,'(3e13.5)') X(1,i) , X(2,i) , 0.0
      ENDDO
 
      WRITE (logic,*) 'POINT_DATA' , NUMnp
 
      WRITE (logic,*) 'SCALARS atoms int  1'
      WRITE (logic,*) 'LOOKUP_TABLE default'
      DO i = 1 , NUMnp
         WRITE (logic,*) ISRelaxed(i)
      ENDDO
99001 FORMAT (3E15.6,i10)
 
 
      END SUBROUTINE STATUSCALC
!*==latticecheck.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE LATTICECHECK(X)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--LATTICECHECK179
      DOUBLE PRECISION X(NXDm,NUMnp)
      INTEGER i , igrain(3)
      DOUBLE PRECISION pm(3)
 
! make sure nodes are on atomic sites
      DO i = 1 , NUM2dnode
         CALL NEARESTBSITE(X(1,i),1,.TRUE.,pm,igrain(1))
         IF ( (ABS(pm(1)-X(1,i))>0.0001) .OR. (ABS(pm(2)-X(2,i))>0.0001)&
     &        ) THEN
            WRITE (*,*) '***WARNING: node is not on atomic site!'
            WRITE (*,*) i , X(1,i) , X(2,i) , X(3,i)
            WRITE (*,*) igrain(1) , pm(1) , pm(2) , pm(3)
            WRITE (*,*) pm(1) - X(1,i) , pm(2) - X(2,i) , pm(3) - X(3,i)
         ENDIF
      ENDDO
      END SUBROUTINE LATTICECHECK
!*==plotesi.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
!**---------------------------------------------------------------
!**     PlotEsi : makes the esi.tec file
!**
!--
      SUBROUTINE PLOTESI(X,Ix)
      USE MOD_GLOBAL
      USE MOD_FILE
      IMPLICIT NONE
!*--PLOTESI206
 
      DOUBLE PRECISION X(NXDm,*)
      INTEGER Ix(NEN1,*)
!
      CHARACTER*80 vars , filename
      INTEGER logic , n , i1 , i2 , i3 , i , iregion
      DOUBLE PRECISION x1 , x2 , x3 , y1 , y2 , y3
 
!     Write file header
      filename = 'out/esi.tec'
      CALL IOFILE(filename,'formatted  ',logic,.TRUE.)
      WRITE (logic,'('' VARIABLES = X Y IREGION'')')
      WRITE (logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'', J = '',i5,     '', F = FEP&
     &OINT'')') NUMel*3 , NUMel
 
      DO n = 1 , NUMel
         iregion = Ix(NEN1,n)
         i1 = Ix(1,n)
         i2 = Ix(2,n)
         i3 = Ix(3,n)
         x1 = X(1,i1)
         x2 = X(1,i2)
         x3 = X(1,i3)
         y1 = X(2,i1)
         y2 = X(2,i2)
         y3 = X(2,i3)
         WRITE (logic,'(2e13.5,i4)') x1 , y1 , iregion
         WRITE (logic,'(2e13.5,i4)') x2 , y2 , iregion
         WRITE (logic,'(2e13.5,i4)') x3 , y3 , iregion
      ENDDO
      DO n = 1 , NUMel
         i = (n-1)*3
         WRITE (logic,'(4i6)') i + 1 , i + 2 , i + 3 , i + 3
      ENDDO
      CLOSE (logic)
      END SUBROUTINE PLOTESI
!*==plotesi_vtk.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
 
 
 
 
!***********************************************************************
!**---------------------------------------------------------------
!**     PlotEsi : makes the esi.tec file
!**
!--
      SUBROUTINE PLOTESI_VTK(X,Ix)
      USE MOD_GLOBAL
      USE MOD_FILE
      IMPLICIT NONE
!*--PLOTESI_VTK259
 
      DOUBLE PRECISION X(NXDm,*)
      INTEGER Ix(NEN1,*)
!
      CHARACTER*80 vars , filename
      INTEGER logic , n , i1 , i2 , i3 , i , iregion
      DOUBLE PRECISION x1 , x2 , x3 , y1 , y2 , y3
      INTEGER scal(NUMel)
 
!     Write file header
      filename = 'out/esi.vtk'
      CALL IOFILE(filename,'formatted  ',logic,.TRUE.)
!$$$      write(logic,'('' VARIABLES = X Y IREGION'')')
!$$$      write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'', J = '',i5,
!$$$     1 '', F = FEPOINT'')') numel*3,numel
 
      WRITE (logic,FMT='(A)') '# vtk DataFile Version 2.0'
      WRITE (logic,FMT='(A)') 'Detection Band Check from CADD'
      WRITE (logic,FMT='(A)') 'ASCII'
      WRITE (logic,FMT='(A)') 'DATASET UNSTRUCTURED_GRID'
      WRITE (logic,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , NUMel*3 , 'float'
      DO n = 1 , NUMel
         iregion = Ix(NEN1,n)
         i1 = Ix(1,n)
         i2 = Ix(2,n)
         i3 = Ix(3,n)
         x1 = X(1,i1)
         x2 = X(1,i2)
         x3 = X(1,i3)
         y1 = X(2,i1)
         y2 = X(2,i2)
         y3 = X(2,i3)
         WRITE (logic,'(3e13.5)') x1 , y1 , 0.0
         WRITE (logic,'(3e13.5)') x2 , y2 , 0.0
         WRITE (logic,'(3e13.5)') x3 , y3 , 0.0
      ENDDO
      WRITE (logic,*)
      WRITE (logic,'(A5,1X,I7,1X,I7)') 'CELLS' , NUMel , 4*NUMel
 
      DO n = 1 , NUMel
         i = (n-1)*3 - 1
         WRITE (logic,'(4i6)') 3 , i + 1 , i + 2 , i + 3
      ENDDO
      WRITE (logic,'(A10,1X,I7)') 'CELL_TYPES' , NUMel
 
      DO n = 1 , NUMel
         WRITE (logic,FMT='(5(1x,I7))') 5
      ENDDO
 
      WRITE (logic,*) 'CELL_DATA' , NUMel
      WRITE (logic,*) 'SCALARS DB integer'
      WRITE (logic,*) 'LOOKUP_TABLE default'
 
      DO n = 1 , NUMel
         WRITE (logic,'(i5)') Ix(NEN1,n)
      ENDDO
 
      CLOSE (logic)
      END SUBROUTINE PLOTESI_VTK
!*==checktime.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE CHECKTIME
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--CHECKTIME324
      WRITE (*,*) '**** Current time is:' , TIMe
      WRITE (*,*) '    and time step is:' , DT
      END SUBROUTINE CHECKTIME
