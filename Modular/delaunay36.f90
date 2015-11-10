!*==delaunay.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**    delaunay  :  generate the mesh
!**
!**   Algorithm :-
!**        allocate temp storage
!**        call contri
!**        write a mesh file for debugging purposes
!**        deallocate temp storage
!--
      SUBROUTINE DELAUNAY(Id,X,Ix,F,B,Itx)
      USE MOD_GLOBAL
      USE MOD_BOUNDARY
      IMPLICIT NONE
!*--DELAUNAY15
!
      INTEGER Id(NDF,*) , Ix(NEN1,*)
      INTEGER Itx(3,*)
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*) , F(NDF,*)
!
      INTEGER , POINTER :: list(:) , w(:) , ixnew(:,:)
      INTEGER i , j
      CHARACTER*80 filename
 
      ! Allocate temporary storage
      ALLOCATE (list(NUMnp))
      ALLOCATE (w((NUMnp+3)*2))
      ALLOCATE (ixnew(3,(2*NUMnp+1)))
 
      ! create plot file of mesh
      filename = 'out/before.tec'
      CALL MESHOUT(filename,NEN1,NXDm,NDF,X,B,Ix,NUMel,NUMnp)
 
      DO i = 1 , NUMnp
         list(i) = i
      ENDDO
      CALL CONTRI(NUMnp,NUMnp,NCE,NCB,ELIst,X,NXDm,list,w,ixnew,Itx,&
     &            NUMel)
 
      IF ( NUMel>MAXel ) THEN
         WRITE (6,&
     &'(///''  **error** detected by subroutine delaunay''           /''&
     &  **number of elements exceeds maximum''                       /''&
     &  **elements requested = '',i5                                 /''&
     &  **elements maximum   = '',i5)') NUMel , MAXel
         STOP
      ENDIF
!
!     this is necessary because contri expects ix to be (3,numel) but we
!     have (4,numel).
 
      DO i = 1 , NUMel
         DO j = 1 , 3
            Ix(j,i) = ixnew(j,i)
         ENDDO
      ENDDO
      NEQ = NDF*NUMnp
 
      ! create plot file of mesh
      filename = 'out/after.tec'
      CALL MESHOUT(filename,NEN1,NXDm,NDF,X,B,Ix,NUMel,NUMnp)
 
      ! Deallocate storage
      DEALLOCATE (w,ixnew,list)
 
      END SUBROUTINE DELAUNAY
!*==meshout.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**     meshout  :  dump a mesh to a tecplot file
!**
!**   Non-Obvious Parameters :-
!**     filename   (in) : output file name
!**           ne   (in) : number of elements
!**           np   (in) : number of nodes
!--
      SUBROUTINE MESHOUT(Filename,Nen1,Nxdm,Ndf,X,B,Ix,Ne,Np)
      USE MOD_BOUNDARY
      USE MOD_FILE
      IMPLICIT NONE
!*--MESHOUT80
!
      INTEGER Nen1 , Nxdm , Ndf , Ix(Nen1,*) , Ne , Np
      DOUBLE PRECISION X(Nxdm,*) , B(Ndf,*)
      CHARACTER*80 Filename
!
      INTEGER logic , i , j
      CALL IOFILE(Filename,'formatted  ',logic,.TRUE.)
      IF ( Ne>0 .AND. Np>0 ) THEN
         WRITE (logic,*) 'zone n=' , Np , ', e=' , Ne , &
     &                   ',f=fepoint,et=triangle'
      ELSEIF ( Np>0 ) THEN
         WRITE (logic,*) 'zone'
      ENDIF
      DO i = 1 , Np
         WRITE (logic,'(4e14.5)') X(1,i) , X(2,i) , X(1,i) + B(1,i) , &
     &                            X(2,i) + B(2,i)
      ENDDO
      DO i = 1 , Ne
         WRITE (logic,'(3i6)') (Ix(j,i),j=1,3)
      ENDDO
      IF ( NCE==0 ) RETURN
      WRITE (logic,*) 'zone n=' , NCE*2 , ', e=' , NCE , &
     &                ',f=fepoint,et=triangle'
      DO i = 1 , NCE
         j = ELIst(1,i)
         WRITE (logic,'(4e14.5)') X(1,j) , X(2,j) , X(1,j) + B(1,j) , &
     &                            X(2,j) + B(2,j)
         j = ELIst(2,i)
         WRITE (logic,'(4e14.5)') X(1,j) , X(2,j) , X(1,j) + B(1,j) , &
     &                            X(2,j) + B(2,j)
      ENDDO
      DO i = 1 , NCE
         WRITE (logic,'(3i6)') 2*i , 2*i - 1 , 2*i
      ENDDO
      CLOSE (logic)
      END SUBROUTINE MESHOUT
