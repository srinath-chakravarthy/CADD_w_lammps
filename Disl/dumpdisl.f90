!*==mod_disl_files.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!$$$  subroutine dumpdisl(scale,logic,umag,b,ndf)
!$$$  use mod_dd_slip
!$$$  implicit none
!$$$  use mod_disl_parameters
!$$$  use mod_fem_parameters
!$$$  double precision xdef(2),u(3),umag,scale,bhat(3),s(3)
!$$$  integer logic,i,j,k1,k2,iel,node,ndf
!$$$  double precision xl(2,3),bl(3,3),b(ndf,*)
!$$$  logical intri,on,flag
!$$$  if(ndisl.lt.1) return
!$$$  write(logic,'('' VARIABLES = X Y UX UY UZ BX BY BZ'')')
!$$$  write(logic,*) 'ZONE'
!$$$  do i=1,ndisl
!$$$  iel=elem_disl(i)
!$$$  if(iel.gt.0) then
!$$$  do j=1,3
!$$$  node=iconn(j,iel)
!$$$  xl(1:2,j)=x0(1:2,node)
!$$$  call disl_displ(xl(1,j),u)
!$$$  bl(1:3,j)=b(1:3,imap(node))-u(1:3)
!$$$  enddo
!$$$  flag=intri(xl(1,1),xl(1,2),xl(1,3),r_disl(1,i),s,on)
!$$$  do k1=1,3
!$$$  bhat(k1)=0.
!$$$  do k2=1,3
!$$$  bhat(k1)=bhat(k1)+s(k2)*bl(k1,k2)
!$$$  enddo
!$$$  enddo
!$$$  else
!$$$  bhat=0
!$$$  endif
!$$$  call disl_displ(r_disl(1,i),u)
!$$$  u=u+bhat
!$$$  xdef(1:2)=r_disl(1:2,i)+umag*u(1:2)
!$$$  write(logic,1000) xdef,scale*u(1:3),burgers(1:3,i)
!$$$  end do
!$$$  1000   format(8e15.6)
!$$$  end

MODULE MOD_DISL_FILES
  IMPLICIT NONE
  !*--MOD_DISL_FILES43
  TYPE DISL_FILES
     CHARACTER(LEN=80) :: FNAME
  END TYPE DISL_FILES

END MODULE MOD_DISL_FILES
!*==dumpdisl_vtk.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015


SUBROUTINE DUMPDISL_VTK1(Scale,Logic,Umag,B,Ndf)
  USE MOD_DD_SLIP
  USE MOD_DISL_PARAMETERS
  USE MOD_FEM_PARAMETERS
  IMPLICIT NONE

  DOUBLE PRECISION xdef(2) , u(3) , Umag , Scale , bhat(3) , s(3)
  INTEGER Logic , i , j , k1 , k2 , iel , node , Ndf , k
  INTEGER idisl , li , n , l
  DOUBLE PRECISION xl(2,3) , bl(3,3) , B(Ndf,*)
  LOGICAL intri , on , flag
  DOUBLE PRECISION :: sd , xd(4,3)
  INTEGER , DIMENSION(:,:,:) , ALLOCATABLE :: lines
  DOUBLE PRECISION :: cosi , sini , ds , ds2 , bsign
  ds = 5.0D0
  ds2 = ds/2.D0


  !! VTK files for dislocations
  WRITE (Logic,FMT='(A)') '# vtk DataFile Version 2.0'
  WRITE (Logic,FMT='(A)') 'Dislocations from CADD'
  WRITE (Logic,FMT='(A)') 'ASCII'
  WRITE (Logic,FMT='(A)') 'DATASET POLYDATA'
  WRITE (Logic,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , ndisl*4 , 'float'
  ALLOCATE (lines(2,3,ndisl))  
  lines = 0
  n = 0
  do idisl = 1, ndisl
     xd = 0.0
     xd(1,1:2) = r_disl(1:2, idisl)
     xd(3,1:2) = xd(1,1:2) - ds2*burgers(1:2,idisl)
     xd(4,1:2) = xd(1,1:2) + ds2*burgers(1:2,idisl)
     xd(2,1) = xd(1,1) + ds2*burgers(2,idisl)
     xd(2,2) = xd(1,2) - ds2*burgers(1,idisl)
     n = n + 1
     lines(1,1,n) = 2
     lines(1,2,n) = 4*(n-1) + 1 - 1
     lines(1,3,n) = 4*(n-1) + 2 - 1
     lines(2,1,n) = 2
     lines(2,2,n) = 4*(n-1) + 3 - 1
     lines(2,3,n) = 4*(n-1) + 4 - 1
     DO k = 1 , 4
        WRITE (Logic,FMT='(3(1X,E15.8))') (xd(k,l),l=1,3)
     ENDDO
     WRITE (Logic,*)     
  end do

  WRITE (Logic,*)
  WRITE (Logic,FMT='(A5,1X,I7,1X,I7)') 'LINES' , 2*ndisl ,6*ndisl
  DO i = 1 , ndisl
     DO j = 1 , 2
        DO k = 1 , 3
           WRITE (Logic,FMT='(1X,I7)',ADVANCE='no') lines(j,k,i)
        ENDDO
        WRITE (Logic,*)
     ENDDO
  ENDDO
  DEALLOCATE (lines)

END SUBROUTINE DUMPDISL_VTK1



SUBROUTINE DUMPDISL_VTK(Scale,Logic,Umag,B,Ndf)
  USE MOD_DD_SLIP
  USE MOD_DISL_PARAMETERS
  USE MOD_FEM_PARAMETERS
  IMPLICIT NONE
  !*--DUMPDISL_VTK56

  DOUBLE PRECISION xdef(2) , u(3) , Umag , Scale , bhat(3) , s(3)
  INTEGER Logic , i , j , k1 , k2 , iel , node , Ndf , k
  INTEGER islp , li , n , l
  DOUBLE PRECISION cphi , sphi , x , y , xstart , xend , ystart , &
       &                 yend
  DOUBLE PRECISION xl(2,3) , bl(3,3) , B(Ndf,*)
  LOGICAL intri , on , flag
  DOUBLE PRECISION :: sd , xd(4,3)
  INTEGER , DIMENSION(:,:,:) , ALLOCATABLE :: lines
  DOUBLE PRECISION :: cosi , sini , ds , ds2 , bsign
  ds = 50.0D0
  ds2 = ds/2.D0

  !! VTK files for dislocations
  WRITE (Logic,FMT='(A)') '# vtk DataFile Version 2.0'
  WRITE (Logic,FMT='(A)') 'Dislocations in DD bending'
  WRITE (Logic,FMT='(A)') 'ASCII'
  WRITE (Logic,FMT='(A)') 'DATASET POLYDATA'
  WRITE (Logic,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , TOT_disl*4 , &
       &       'float'
  ALLOCATE (lines(2,3,TOT_disl))

  lines = 0
  n = 0
  DO islp = i , NSLp
     if (ndis(islp) > 0) then 
        DO i = 1 , NDIs(islp)
           li = LOCphi(islp)
           cosi = COSphi(li)
           sini = SINphi(li)
           sd = SDIs(i,islp)
           IF ( YENdslp(islp)<0.0D0 ) sd = -sd
           IF ( B_Dd(i,islp)<0.0D0 ) THEN
              bsign = -1.0D0
           ELSE
              bsign = 1.0D0
           ENDIF
           xd = 0.0
           xd(1,1) = XSLp(islp) + cosi*sd
           xd(1,2) = YSLp(islp) + sini*sd
           xd(2,1) = xd(1,1) + bsign*ds*sini
           xd(2,2) = xd(1,2) - bsign*ds*cosi
           xd(3,1) = xd(1,1) - cosi*ds
           xd(3,2) = xd(1,2) - sini*ds
           xd(4,1) = xd(1,1) + cosi*ds
           xd(4,2) = xd(1,2) + sini*ds
           n = n + 1
           lines(1,1,n) = 2
           lines(1,2,n) = 4*(n-1) + 1 - 1
           lines(1,3,n) = 4*(n-1) + 2 - 1
           lines(2,1,n) = 2
           lines(2,2,n) = 4*(n-1) + 3 - 1
           lines(2,3,n) = 4*(n-1) + 4 - 1
           DO k = 1 , 4
              WRITE (Logic,FMT='(3(1X,E15.8))') (xd(k,l),l=1,3)
           ENDDO
           WRITE (Logic,*)
        ENDDO
     end if
  ENDDO
  WRITE (Logic,*)
  WRITE (Logic,FMT='(A5,1X,I7,1X,I7)') 'LINES' , 2*TOT_disl , &
       &       6*TOT_disl
  DO i = 1 , TOT_disl
     DO j = 1 , 2
        DO k = 1 , 3
           WRITE (Logic,FMT='(1X,I7)',ADVANCE='no') lines(j,k,i)
        ENDDO
        WRITE (Logic,*)
     ENDDO
  ENDDO
  DEALLOCATE (lines)


END SUBROUTINE DUMPDISL_VTK
!*==dumpdisl.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015


SUBROUTINE DUMPDISL(Scale,Logic,Umag,B,Ndf)
  USE MOD_DD_SLIP
  USE MOD_DISL_PARAMETERS
  USE MOD_FEM_PARAMETERS
  IMPLICIT NONE
  !*--DUMPDISL140
  DOUBLE PRECISION xdef(2) , u(3) , Umag , Scale , bhat(3) , s(3)
  INTEGER Logic , i , j , k1 , k2 , iel , node , Ndf
  INTEGER islp , li
  DOUBLE PRECISION cphi , sphi , x , y , xstart , xend , ystart , &
       &                 yend
  DOUBLE PRECISION xl(2,3) , bl(3,3) , B(Ndf,*)
  LOGICAL intri , on , flag
  IF ( NDIsl<1 ) RETURN
  WRITE (Logic,'('' VARIABLES = X Y'')')

  WRITE (Logic,*) 'ZONE T= pos1'
  DO islp = 1 , NSLp
     li = LOCphi(islp)
     cphi = COSphi(li)
     sphi = SINphi(li)
     xstart = XSLp(islp)
     xend = XENdslp(islp)
     ystart = YSLp(islp)
     yend = YENdslp(islp)
     IF ( li==1 ) THEN
        DO i = 1 , NDIs(islp)
           IF ( B_Dd(i,islp)>0 ) THEN
              IF ( yend>0 ) THEN
                 x = xstart + SDIs(i,islp)*cphi
                 y = ystart + SDIs(i,islp)*sphi
              ELSE
                 x = xstart - SDIs(i,islp)*cphi
                 y = ystart - SDIs(i,islp)*sphi
              ENDIF
              WRITE (Logic,FMT='(2(1x,E16.9))') x , y
           ENDIF
        ENDDO
     ENDIF
  ENDDO

  WRITE (Logic,*) 'ZONE T= neg1'
  DO islp = 1 , NSLp
     li = LOCphi(islp)
     cphi = COSphi(li)
     sphi = SINphi(li)
     xstart = XSLp(islp)
     xend = XENdslp(islp)
     ystart = YSLp(islp)
     yend = YENdslp(islp)
     IF ( li==1 ) THEN
        DO i = 1 , NDIs(islp)
           IF ( B_Dd(i,islp)<0 ) THEN
              IF ( yend>0 ) THEN
                 x = xstart + SDIs(i,islp)*cphi
                 y = ystart + SDIs(i,islp)*sphi
              ELSE
                 x = xstart - SDIs(i,islp)*cphi
                 y = ystart - SDIs(i,islp)*sphi
              ENDIF
              WRITE (Logic,FMT='(2(1x,E16.9))') x , y
           ENDIF
        ENDDO
     ENDIF
  ENDDO


  WRITE (Logic,*) 'ZONE T= pos2'
  DO islp = 1 , NSLp
     li = LOCphi(islp)
     cphi = COSphi(li)
     sphi = SINphi(li)
     xstart = XSLp(islp)
     xend = XENdslp(islp)
     ystart = YSLp(islp)
     yend = YENdslp(islp)
     IF ( li==3 ) THEN
        DO i = 1 , NDIs(islp)
           IF ( B_Dd(i,islp)>0 ) THEN
              IF ( yend>0 ) THEN
                 x = xstart + SDIs(i,islp)*cphi
                 y = ystart + SDIs(i,islp)*sphi
              ELSE
                 x = xstart - SDIs(i,islp)*cphi
                 y = ystart - SDIs(i,islp)*sphi
              ENDIF
              WRITE (Logic,FMT='(2(1x,E16.9))') x , y
           ENDIF
        ENDDO
     ENDIF
  ENDDO


  WRITE (Logic,*) 'ZONE T= neg2'
  DO islp = 1 , NSLp
     li = LOCphi(islp)
     cphi = COSphi(li)
     sphi = SINphi(li)
     xstart = XSLp(islp)
     xend = XENdslp(islp)
     ystart = YSLp(islp)
     yend = YENdslp(islp)
     IF ( li==2 ) THEN
        DO i = 1 , NDIs(islp)
           IF ( B_Dd(i,islp)<0 ) THEN
              IF ( yend>0 ) THEN
                 x = xstart + SDIs(i,islp)*cphi
                 y = ystart + SDIs(i,islp)*sphi
              ELSE
                 x = xstart - SDIs(i,islp)*cphi
                 y = ystart - SDIs(i,islp)*sphi
              ENDIF
              WRITE (Logic,FMT='(2(1x,E16.9))') x , y
           ENDIF
        ENDDO
     ENDIF
  ENDDO


  WRITE (Logic,*) 'ZONE T= pos3'
  DO islp = 1 , NSLp
     li = LOCphi(islp)
     cphi = COSphi(li)
     sphi = SINphi(li)
     xstart = XSLp(islp)
     xend = XENdslp(islp)
     ystart = YSLp(islp)
     yend = YENdslp(islp)
     IF ( li==3 ) THEN
        DO i = 1 , NDIs(islp)
           IF ( B_Dd(i,islp)>0 ) THEN
              x = xstart + SDIs(i,islp)*cphi
              y = ystart + SDIs(i,islp)*sphi
           ENDIF
           WRITE (Logic,FMT='(2(1x,E16.9))') x , y
        ENDDO
     ENDIF
  ENDDO

  WRITE (Logic,*) 'ZONE T= neg3'
  DO islp = 1 , NSLp
     li = LOCphi(islp)
     cphi = COSphi(li)
     sphi = SINphi(li)
     xstart = XSLp(islp)
     xend = XENdslp(islp)
     ystart = YSLp(islp)
     yend = YENdslp(islp)
     IF ( li==3 ) THEN
        DO i = 1 , NDIs(islp)
           IF ( B_Dd(i,islp)<0 ) THEN
              x = xstart + SDIs(i,islp)*cphi
              y = ystart + SDIs(i,islp)*sphi
           ENDIF
           WRITE (Logic,FMT='(2(1x,E16.9))') x , y
        ENDDO
     ENDIF
  ENDDO
END SUBROUTINE DUMPDISL
