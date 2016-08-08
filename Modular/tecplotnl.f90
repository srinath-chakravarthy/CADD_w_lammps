!s!*==ma02.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------
!**     ma02  :  produces tecplot read-able files of the mesh and
!**              various contour variables.
!**
!**   Calling format:
!**   ma02,key,filename,index,scale,umag, nprint
!**
!**   key - type of plot - I only trust 'bcon', 'disp', 'disl' (see belo
!**   filename - name of output file (routine adds a numerical counter)
!**   index - not used
!**   scale - magnifies the output variable (scale=1.0 for no magnificat
!**   umag  - exagerates the displacements on the deformed mesh
!**           (umag=1.0 for true deformed mesh). Read from input
!**	      by the subroutine freein()
!**
!**   all above comments are pre 9/1/06
!**
!**   things have been updated 9/28/06 now key='stra' plots strain value
!**   and k='atom' plots an atom eye file
!**
!**   thing have been updated 10/4/06 so that key='stre' plots stress
!**   values
!**
!**
!--
      SUBROUTINE MA02(Id,X,Ix,F,B,Dr,Db,Input,Flag)
      USE MOD_FILE
      USE MOD_DYNAMO
      USE MOD_GLOBAL
      USE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--MA0233
!
      DOUBLE PRECISION X(NXDm,*) , F(NDF,*) , B(NDF,*) , Dr(*) , Db(*)
      INTEGER Ix(NEN1,*) , Id(NDF,*)
!
!---- tecplot file
!
      CHARACTER*80 Input , filename , temp , fname_d1 , fname_d2
      CHARACTER*4 key
      CHARACTER*6 cnt
      LOGICAL Flag
!c--JS: update icount entries if adding new ikey
      INTEGER , SAVE :: icount(13) = 0 , nprint
      INTEGER i , lower , upper , NEXT , ikey , j , logic , index , idum
      DOUBLE PRECISION dum , scale , umag
!
      key = Input(1:4)
      lower = 4
      upper = NEXT(lower,Input)
      filename((lower-3):(upper-5)) = Input((lower+1):(upper-1))
      DO i = upper - 4 , 80
         filename(i:i) = ' '
      ENDDO
!
      ikey = 0
      IF ( key=='noda' ) ikey = 1
      IF ( key=='disp' ) ikey = 2
                                ! plot of displacements
      IF ( key=='bcon' ) ikey = 3
                                ! plot of boundary conditions
      IF ( key=='ener' ) ikey = 4
      IF ( key=='disl' ) ikey = 5
                                ! plot of discrete dislocations
      IF ( key=='stre' ) ikey = 6
                                ! plot of stresses is now working
      IF ( key=='stra' ) ikey = 7
                                ! plot of strains this is now working
      IF ( key=='viri' ) ikey = 8
                                ! plot of averaged virial stresses
      IF ( key=='atom' ) ikey = 9
                                ! plot of atoms to be read with atomeye
      IF ( key=='virH' ) ikey = 10
                                 ! plot the viri stress for H atoms
      IF ( key=='ovit' ) ikey = 11
      ! Write Cfg file for ovito

      
      IF ( key=='lamp' ) ikey = 12
      ! Write lammps dump file

      IF ( key=='olmp' ) ikey = 13

      IF ( ikey==0 ) THEN
         WRITE (*,*) 'ERROR: unknown key'
         RETURN
      ENDIF
!
      icount(ikey) = icount(ikey) + 1
      WRITE (cnt,'(i6)') icount(ikey)
      IF ( icount(ikey)<10 ) cnt(1:5) = '00000'
      IF ( icount(ikey)<100 ) cnt(1:4) = '0000'
      IF ( icount(ikey)<1000 ) cnt(1:3) = '000'
      IF ( icount(ikey)<1000 ) cnt(1:2) = '00'
      IF ( icount(ikey)<10000 ) cnt(1:1) = '0'
 
      i = upper - 4
      DO j = 1 , upper - 5
         IF ( filename(j:j)=='.' ) i = j
      ENDDO
      temp = filename(i:)
      filename(i:i+5) = cnt
      filename(i+6:) = temp
 
!      if (ndisl > 0 ) then
!$$$       if (ikey .ne. 5) then
      CALL IOFILE(filename,'formatted  ',logic,.TRUE.)
!$$$      else
!$$$         if (ndisl > 0) then
!$$$            call iofile(filename,'formatted  ',logic,.true.)
!$$$         endif
!$$$      endif
 
      lower = upper
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,index,dum,1)
      lower = upper
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,idum,scale,2)
      lower = upper
      upper = NEXT(lower,Input)
      CALL FREEIN(Input,lower,upper,idum,umag,2)
      lower = upper
!$$$      upper = next(lower,input)
!$$$      call freein(input,lower,upper,idum,nprint,2)
 
 
 
      PRINT * , 'NPRINT = ' , nprint , umag , scale
!
!
      IF ( ikey==5 ) THEN
         PRINT * , 'Total no. of dislocations = ' , NDIsl
!        if (ndisl > 0 .or. nprint .eq. 0) then
!         if (ndisl > 0) then
         CALL DUMPDISL_VTK(scale,logic,umag,B,NDF)
!         else
!           print *, 'no dislocations to print in output file'
!        endif
      ELSE
         IF ( ikey<6 .OR. ikey>8 ) CALL DUMPIT(X,Ix,B,Db,Id,F,Dr,scale,&
     &        logic,key,index,umag)
         IF ( ikey>=6 .OR. ikey<=8 )&
     &        CALL DUMPIT_VTK(X,Ix,B,Db,Id,F,Dr,scale,logic,key,index,&
     &        umag)
      ENDIF
      CLOSE (logic)
      END SUBROUTINE MA02
!*==getstresses.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!******************************************************************
      SUBROUTINE GETSTRESSES(Strain,Stress)
      IMPLICIT NONE
!*--GETSTRESSES152
      DOUBLE PRECISION Strain(3,3) , Stress(3,3)
      DOUBLE PRECISION vstrain(6) , vstress(6)
      DOUBLE PRECISION t1 , t2 , theta , pi
      INTEGER i , j
 
      CHARACTER*80 error_message
      DOUBLE PRECISION CC(6,6)
      DOUBLE PRECISION XE , XNU , XLAmbda , XMU
      INTEGER I_Elas
      COMMON /ELASTIC/ XE , XNU , XLAmbda , XMU , CC , I_Elas
      DOUBLE PRECISION ev_convert
      ev_convert = 1.602176462
 
      IF ( I_Elas/=1 ) THEN
         error_message = 'fe_dmat: call fe_elastic first!'
         CALL ERROR_HANDLER(error_message)
      ENDIF
 
 
 
      vstrain(1) = Strain(1,1)
      vstrain(2) = Strain(2,2)
      vstrain(3) = Strain(3,3)
      vstrain(4) = Strain(2,3) + Strain(3,2)
      vstrain(5) = Strain(1,3) + Strain(3,1)
      vstrain(6) = Strain(1,2) + Strain(2,1)
 
      DO i = 1 , 6
         vstress(i) = 0.0
         DO j = 1 , 6
!          vstress(i)=vstress(i)+cc(i,j)*vstrain(j)*100e9/1e6
            vstress(i) = (vstress(i)+CC(i,j)*vstrain(j))
         ENDDO
      ENDDO
 
      Stress(1,1) = vstress(1)
      Stress(2,2) = vstress(2)
      Stress(3,3) = vstress(3)
      Stress(2,3) = vstress(4)
      Stress(1,3) = vstress(5)
      Stress(1,2) = vstress(6)
      Stress(3,2) = vstress(4)
      Stress(3,1) = vstress(5)
      Stress(2,1) = vstress(6)
 
!      section if we want to plot stress on slip system
!      pi=3.141592653
!      theta=-70.0*pi/180.0
!      t1=cos(theta+pi/2)*stress(1,1)+sin(theta+pi/2)*stress(1,2);
!      t2=cos(theta+pi/2)*stress(1,2)+sin(theta+pi/2)*stress(2,2);
!     // get shear stress //
!      stress(2,3)=t1*cos(theta)+t2*sin(theta);
 
 
 
      END SUBROUTINE GETSTRESSES
!*==dump_atom.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE DUMP_ATOM(X,B,Dr,Logic)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--DUMP_ATOM216
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*), DR(NDF,*)
      INTEGER Logic , ntot , numpnts
      INTEGER n1 , n2 , n3 , npatoms
 
      DOUBLE PRECISION dev(6) , xdef , ydef , zdef , hyd , sigeff , y ,  epseff
      INTEGER numtri , i , j , n , ndfmax , ii , jj , nout , npad ,  startpad
      DOUBLE PRECISION dwx1 , dwx2
      DOUBLE PRECISION box_max(2) , box_min(2)
      DOUBLE PRECISION pe , ke
      INTEGER isurf, atom_type
      
      character(len=1024) :: outstr
 
      DOUBLE PRECISION :: umag
      CHARACTER tatom*2
      DOUBLE PRECISION mass

 
      umag = 1.0D0
 
      pe = 0.0D0
      DO j = 1 , 2
         box_max(j) = -1.0D30
         box_min(j) = 1.0D30
      ENDDO
      npatoms = 0
      npad = 0
      ntot = 0
      DO i = 1 , NUMnp
         IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. ISRelaxed(i) ==-1 ) THEN
            IF ( ISRelaxed(i)==-1 ) THEN
               IF ( npad==0 ) startpad = i
               npad = npad + 1
            ELSE
               pe = pe + ENErgy(i)
               npatoms = npatoms + 1
            ENDIF
            ntot = ntot + 1
            DO j = 1 , 2
               IF ( X(j,i)+umag*B(j,i)>box_max(j) ) box_max(j) = X(j,i) + umag*B(j,i)
               IF ( X(j,i)+umag*B(j,i)<box_min(j) ) box_min(j) = X(j,i) + umag*B(j,i)
            ENDDO
         ENDIF
      ENDDO
      startpad = npatoms - npad + 1

      do i = 1, numnp
         if (isrelaxed(i) == -1) then
            atom_type = 2
         else if (isrelaxed(i) == 2) then
            atom_type = 3
         else
            atom_type = 1
         end if

         if (i == 1) then
            write(outstr, *) 'ITEM: TIMESTEP'
            write(logic, *) adjustl(trim(outstr))

            write(logic, '(I1)') 0

            write(outstr, *) 'ITEM: NUMBER OF ATOMS'
            write(logic, *) adjustl(trim(outstr))

            write(outstr, '(I7)') ntot
            write(logic, *) adjustl(trim(outstr))

            write(outstr, *) 'ITEM: BOX BOUNDS ss ss pp'
            write(logic, *) adjustl(trim(outstr))

            write(outstr, '(2F16.9)') box_min(1), box_max(1)
            write(logic, *) adjustl(trim(outstr))
            write(outstr, '(2F16.9)') box_min(2), box_max(2)
            write(logic, *) adjustl(trim(outstr))
            write(outstr, '(2F16.9)') -0.5, 0.5
            write(logic, *) adjustl(trim(outstr))

            write(outstr, '(A300)') 'ITEM: ATOMS id type x y z fx fy fz'
            write(logic, *) adjustl(trim(outstr))
         end if

         xdef = x(1,i) + b(1,i)
         ydef = x(2,i) + b(2,i)
         zdef = x(3,i) + b(3,i)
         if (isrelaxed(i) /= 0) then 
            write(outstr,'(I7,1X,I2,1X,6(F16.9))') &
                 i, atom_type, xdef, ydef, zdef, dr(1,i), dr(2,i), dr(3,i)
            write(logic, '(A107)') adjustl(trim(outstr))
         end if
      end do
      
!!$      PRINT * , 'Start pad atoms = ' , startpad
!!$      PRINT * , 'Total Potential Energy = ' , pe
!!$      PRINT * , 'Average Potential Energy = ' , DBLE(pe/npatoms)
 
!!$      DO i = 1 , NUMnp
!!$         IF ( i==1 ) THEN
!!$             WRITE (Logic,'(''Number of particles = '',i5)') ntot
!!$            WRITE (Logic,'(''A = 1.0 Angstrom '')')
!!$            WRITE (Logic,'(''H0(1,1) = '',f10.4,'' A'')') box_max(1) - box_min(1)
!!$            WRITE (Logic,'(''H0(1,2) = 0 A'')')
!!$            WRITE (Logic,'(''H0(1,3) = 0 A'')')
!!$            WRITE (Logic,'(''H0(2,1) = 0 A'')')
!!$            WRITE (Logic,'(''H0(2,2) = '',f10.4,'' A'')') box_max(2) - box_min(2)
!!$            WRITE (Logic,'(''H0(2,3) = 0 A'')')
!!$            WRITE (Logic,'(''H0(3,1) = 0 A'')')
!!$            WRITE (Logic,'(''H0(3,2) = 0 A'')')
!!$            WRITE (Logic,'(''H0(3,3) = '',f10.3,'' A'')') Z_Length
!!$            
!!$         ENDIF
!!$         
!!$ 
!!$         xdef = (X(1,i)+umag*B(1,i)+box_min(1))/(box_max(1)-box_min(1))
!!$         ydef = (X(2,i)+umag*B(2,i)+box_max(1))/(box_max(2)-box_min(2))
!!$         zdef = (X(3,i)+umag*B(3,i))/Z_Length
!!$
!!$
!!$         IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 ) THEN
!!$            tatom = "Al"
!!$            mass = 13.0
!!$         ENDIF
!!$         IF ( ISRelaxed(i)==-1 ) THEN
!!$            tatom = "Ni"
!!$            mass = 14.0
!!$         ENDIF
!!$
!!$         !C           zdef = b(3,i)
!!$         IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. ISRelaxed(i)&
!!$              &           ==-1 ) WRITE (Logic,'(f4.0,1X,A2, 1X, 6f16.11)') mass ,&
!!$              &                         tatom , xdef , ydef , zdef , 0.0 , 0.0 , &
!!$              &                         0.0
!!$         end DO

     END SUBROUTINE DUMP_ATOM
!*==dumpit.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
!***********************************************************************
      SUBROUTINE DUMPIT(X,Ix,B,Db,Id,F,Dr,Scale,Logic,Key,Index,Umag)
      USE MOD_GLOBAL
      USE MOD_DYNAMO
      IMPLICIT NONE
!*--DUMPIT331
!
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*) , Db(NDF,*) , F(NDF,*) , &
     &                 Dr(NDF,*) , Scale , Umag
      INTEGER Id(NDF,*) , Ix(NEN1,*) , Index , Logic , numpnts , ntot
      CHARACTER*4 Key
!
      INTEGER numtri , i , j , n , NDFMAX , ii , jj , nout , npad , &
     &        startpad
      CHARACTER out*100
      CHARACTER tatom*2
      DOUBLE PRECISION mass
      PARAMETER (NDFMAX=18)
      DOUBLE PRECISION dev(6) , xdef , ydef , zdef , hyd , sigeff , y , &
     &                 epseff
 
      INTEGER n1 , n2 , n3 , npatoms
      DOUBLE PRECISION det , amat(3,2) , epsloc(3,3) , tarray(3,3)
      DOUBLE PRECISION , ALLOCATABLE :: nstrain(:,:,:)
      DOUBLE PRECISION , ALLOCATABLE :: avgnum(:)
      DOUBLE PRECISION dwx1 , dwx2
      DOUBLE PRECISION box_max(2) , box_min(2)
      DOUBLE PRECISION pe , ke
      INTEGER isurf, atom_type
      character(len = 1024) outstr
      pe = 0.0
!
      numtri = NUMel
      IF ( Key/='atom' .AND. Key/='virH' .AND. Key/='lamp') THEN
         IF ( Key/='viri' .AND. Key/='stra' .AND. Key/='stre' ) THEN
!!$            write(logic,'('' TITLE = " '',a4,'' "'')') key
         ENDIF
      ENDIF
 
 
!     if we are doing strains or stresses then calculate them at each no
      IF ( Key=='stra' .OR. Key=='stre' .OR. Key=='viri' ) THEN
 
         ALLOCATE (nstrain(3,3,NUMnp))
         ALLOCATE (avgnum(NUMnp))
         DO i = 1 , NUMnp
            avgnum(i) = 0.0
            DO j = 1 , 3
               DO ii = 1 , 3
                  nstrain(ii,j,i) = 0.0
               ENDDO
            ENDDO
         ENDDO
 
         DO i = 1 , NUMel
            n1 = Ix(1,i)
            n2 = Ix(2,i)
            n3 = Ix(3,i)
            det = (X(1,n1)-X(1,n3))*(X(2,n2)-X(2,n3))&
     &            - (X(1,n2)-X(1,n3))*(X(2,n1)-X(2,n3))
            amat(1,1) = (X(2,n2)-X(2,n3))/det
            amat(2,1) = (X(2,n3)-X(2,n1))/det
            amat(3,1) = (X(2,n1)-X(2,n2))/det
            amat(1,2) = (X(1,n3)-X(1,n2))/det
            amat(2,2) = (X(1,n1)-X(1,n3))/det
            amat(3,2) = (X(1,n2)-X(1,n1))/det
!           ****compute the green strain time 2*****
            CALL GETELEMENTSTRAIN(B(1,n1),B(1,n2),B(1,n3),amat(1:3,1:2),epsloc(1:3,1:3))
            epsloc(1:3,1:3) = epsloc(1:3,1:3)*0.5
 
            IF ( Key=='stre' .OR. Key=='viri' ) THEN
!           ***** compute stresses in units of eV/A^3
               CALL GETSTRESSES(epsloc,tarray)
            ELSE
               tarray(1:3,1:3) = epsloc(1:3,1:3)
!!$         convert to engineering strains
               tarray(2,3) = tarray(2,3) + tarray(3,2)
               tarray(1,3) = tarray(1,3) + tarray(3,1)
               tarray(1,2) = tarray(1,2) + tarray(2,1)
            ENDIF
 
            DO j = 1 , 3
               DO jj = 1 , 3
                  nstrain(j,jj,n1) = nstrain(j,jj,n1) + tarray(j,jj)
                  nstrain(j,jj,n2) = nstrain(j,jj,n2) + tarray(j,jj)
                  nstrain(j,jj,n3) = nstrain(j,jj,n3) + tarray(j,jj)
               ENDDO
            ENDDO
 
            avgnum(n1) = avgnum(n1) + 1.0
            avgnum(n2) = avgnum(n2) + 1.0
            avgnum(n3) = avgnum(n3) + 1.0
 
         ENDDO
 
         DO i = 1 , NUMnp
            IF ( avgnum(i)>0 ) THEN
               DO j = 1 , 3
                  DO jj = 1 , 3
                     nstrain(j,jj,i) = nstrain(j,jj,i)/avgnum(i)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
 
      ENDIF
 
 
!     find max and min of coordintes so that a box can be made.
      IF ( Key=='atom' .OR. Key=='ovit' .or. Key == 'lamp' .or. Key == 'olmp') THEN
         pe = 0.0D0
         DO j = 1 , 2
            box_max(j) = -1.0D30
            box_min(j) = 1.0D30
         ENDDO
         npatoms = 0
         npad = 0
         ntot = 0
         DO i = 1 , NUMnp
            IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. ISRelaxed(i)==-1 ) THEN
               IF ( ISRelaxed(i)==-1 ) THEN
                  IF ( npad==0 ) startpad = i
                  npad = npad + 1
               ELSE
                  pe = pe + ENErgy(i)
                  npatoms = npatoms + 1
               ENDIF
               ntot = ntot + 1
               DO j = 1 , 2
                  IF ( X(j,i)+Umag*B(j,i)>box_max(j) ) box_max(j) = X(j,i) + Umag*B(j,i)
                  IF ( X(j,i)+Umag*B(j,i)<box_min(j) ) box_min(j) = X(j,i) + Umag*B(j,i)
               ENDDO
            ENDIF
         ENDDO
         startpad = npatoms - npad + 1
         PRINT * , 'Start pad atoms = ' , startpad
         PRINT * , 'Total Potential Energy = ' , pe
         PRINT * , 'Average Potential Energy = ' , DBLE(pe/npatoms)
      ENDIF
 
!     begin loop over nodes for all cases
!      numpnts=max(numnp,numnpp1)
      numpnts = MAX(NUMnp,NUMnpp1)
      DO i = 1 , numpnts
         IF ( Key=='noda' ) THEN
            dev(1) = Scale*B(1,i)
            dev(2) = Scale*B(2,i)
            IF ( i==1 ) THEN
               WRITE (Logic,'('' VARIABLES = X Y UX UY VX VY AX AY'')')
               WRITE (Logic,'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''              ,i5,'', F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
            WRITE (Logic,'(8e14.6)') xdef , ydef , (dev(j),j=1,6)
         ENDIF
         IF ( Key=='disp' ) THEN
            dev(1) = Scale*B(1,i)
            dev(2) = Scale*B(2,i)
            IF ( NDF>=3 ) dev(3) = Scale*B(3,i)
            IF ( NDF==4 ) dev(4) = Scale*B(4,i)
            IF ( i==1 ) THEN
               IF ( NDF>NDFMAX ) THEN
                  WRITE (*,*) &
     &                    '***WARNING: tecplot file can only show first'&
     &                    , NDFMAX
                  WRITE (*,*) '            degrees of freedom'
               ENDIF
               IF ( NDF==2 ) THEN
                  WRITE (Logic,'('' VARIABLES = X Y UX UY '')')
               ELSEIF ( NDF==3 ) THEN
                  WRITE (Logic,'('' VARIABLES = X Y Z UX UY UZ'')')
               ELSEIF ( NDF==4 ) THEN
                  WRITE (Logic,'('' VARIABLES = X Y Z UX UY UZ EZ'')')
               ELSE
                  out = 'VARIABLES = X Y Z UX UY UZ'
                  nout = 26
                  jj = 1
                  DO ii = 4 , MIN(NDF,NDFMAX)
                     IF ( jj==1 ) THEN
                        out(nout+1:nout+3) = ' SX'
                        jj = 2
                     ELSEIF ( jj==2 ) THEN
                        out(nout+1:nout+3) = ' SY'
                        jj = 3
                     ELSEIF ( jj==3 ) THEN
                        out(nout+1:nout+3) = ' SZ'
                        jj = 1
                     ENDIF
                     nout = nout + 3
                  ENDDO
                  WRITE (Logic,*) out
               ENDIF
               WRITE (Logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''              ,i5,'',&
     & F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
!!$C            print*, 'spatial dimension', ndf
            IF ( NDF>=3 ) zdef = Umag*B(3,i) + X(3,i)
            IF ( NDF==2 ) THEN
               WRITE (Logic,'(4e14.6)') xdef , ydef , dev(1) , dev(2)
            ELSEIF ( NDF==3 ) THEN
               WRITE (Logic,'(6e14.6)') xdef , ydef , zdef , dev(1) , dev(2) , dev(3)
            ELSEIF ( NDF==4 ) THEN
               WRITE (Logic,'(7e14.6)') xdef , ydef , zdef , dev(1) , dev(2) , dev(3) , dev(4)
            ELSE
               WRITE (Logic,'(21e14.6)') xdef , ydef , zdef , dev(1) , dev(2) , dev(3) , (B(ii,i),ii=4,MIN(NDF,NDFMAX))
            ENDIF
         ENDIF
 
         IF ( Key=='bcon' ) THEN
!!$     dev(1) = scale*f(1,i)
!!$     dev(2) = scale*f(2,i)
!!$     if (ndf.ge.3) dev(3) = scale*f(3,i)
!!$            dev(1) = scale*dr(1,i)
!!$            dev(2) = scale*dr(2,i)
!!$            if (ndf.ge.3) dev(3) = scale*dr(3,i)
            dev(1) = DBLE(Id(1,i))
            dev(2) = DBLE(Id(2,i))
            IF ( NDF>=3 ) dev(3) = DBLE(Id(3,i))
            IF ( i==1 ) THEN
               IF ( NDF==2 ) THEN
                  WRITE (Logic,'('' VARIABLES = X Y FX FY'')')
               ELSE
                  WRITE (Logic,'('' VARIABLES = X Y FX FY FZ'')')
               ENDIF
               WRITE (Logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'', J =''              ,i5,'',&
     & F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
            IF ( NDF==2 ) THEN
               WRITE (Logic,'(5e14.6)') xdef , ydef , dev(1) , dev(2)
            ELSE
               WRITE (Logic,'(6e14.6)') xdef , ydef , dev(1) , dev(2) , dev(3)
            ENDIF
         ENDIF
 
 
         IF ( Key=='ener' ) THEN
            dev(1) = 1.D0*ENErgy(i)
            IF ( i==1 ) THEN
               WRITE (Logic,'('' VARIABLES = X Y Z E_atom'')')
               WRITE (Logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''              ,i5,'',&
     & F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
!            zdef =  b(3,i)
            zdef = X(3,i) + Umag*B(3,i)
            WRITE (Logic,'(6e14.6)') xdef , ydef , zdef , dev(1)
         ENDIF
 
 
         IF ( Key=='atom' ) THEN
 
            IF ( i==1 ) THEN
 
               WRITE (Logic,'(''Number of particles = '',i5)') ntot
               WRITE (Logic,'(''A = 1.0 Angstrom '')')
               WRITE (Logic,'(''H0(1,1) = '',f10.4,'' A'')') box_max(1) - box_min(1)
!               write(logic,'(''H0(1,2) = '',f10.4,'' A'')') box_min(2)
               WRITE (Logic,'(''H0(1,2) = 0 A'')')
               WRITE (Logic,'(''H0(1,3) = 0 A'')')
               WRITE (Logic,'(''H0(2,1) = 0 A'')')
!              write(logic,'(''H0(2,1) = '',f10.4,'' A'')') box_min(1)
               WRITE (Logic,'(''H0(2,2) = '',f10.4,'' A'')') box_max(2) - box_min(2)
               WRITE (Logic,'(''H0(2,3) = 0 A'')')
               WRITE (Logic,'(''H0(3,1) = 0 A'')')
               WRITE (Logic,'(''H0(3,2) = 0 A'')')
!              write(logic,'(''H0(3,1) = '',f10.4,'' A'')') box_min(1)
!              write(logic,'(''H0(3,2) = '',f10.4,'' A'')') box_min(2)
               WRITE (Logic,'(''H0(3,3) = '',f10.3,'' A'')') Z_Length
               WRITE (Logic,'(''.NO_VELOCITY. '')')
               WRITE (Logic,'(''entry_count = 6 '')')
               WRITE (Logic,'(''auxiliary[0] = pote [eV] '')')
               WRITE (Logic,'(''auxiliary[1] = coord [num] '')')
               WRITE (Logic,'(''auxiliary[2] = Surface [num] '')')
               WRITE (Logic,'(''1.000000 '')')
               WRITE (Logic,'(''Al '')')
            ENDIF
 
            xdef = (X(1,i)+Umag*B(1,i)+box_min(1)) /(box_max(1)-box_min(1))
            ydef = (X(2,i)+Umag*B(2,i)+box_max(1))/(box_max(2)-box_min(2))
            zdef = (X(3,i)+Umag*B(3,i))/Z_Length
!!$           zdef = b(3,i)
            IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 ) THEN
               IF ( NUMneighbors(i)==26 ) isurf = 0
               IF ( (NUMneighbors(i)>22 .AND. NUMneighbors(i)<=26) .OR. NUMneighbors(i)>26 ) isurf = 1
 
               IF ( NUMneighbors(i)<=22 ) isurf = 2
                WRITE (Logic,'(4f16.11,1x,2I4)') xdef , ydef , zdef ,  1.D0*ENErgy(i) , NUMneighbors(i) , isurf
            ENDIF
 
         ENDIF
 
         IF ( Key=='ovit' ) THEN
 
            IF ( i==1 ) THEN
 
               WRITE (Logic,'(''Number of particles = '',i5)') ntot
               WRITE (Logic,'(''A = 1.0 Angstrom '')')
               WRITE (Logic,'(''H0(1,1) = '',f10.4,'' A'')') box_max(1) - box_min(1)
               WRITE (Logic,'(''H0(1,2) = 0 A'')')
               WRITE (Logic,'(''H0(1,3) = 0 A'')')
               WRITE (Logic,'(''H0(2,1) = 0 A'')')
               WRITE (Logic,'(''H0(2,2) = '',f10.4,'' A'')') box_max(2) - box_min(2)
               WRITE (Logic,'(''H0(2,3) = 0 A'')')
               WRITE (Logic,'(''H0(3,1) = 0 A'')')
               WRITE (Logic,'(''H0(3,2) = 0 A'')')
               WRITE (Logic,'(''H0(3,3) = '',f10.3,'' A'')') Z_Length
            ENDIF
 
            xdef = (X(1,i)+Umag*B(1,i)-box_min(1))/(box_max(1)-box_min(1))
            ydef = (X(2,i)+Umag*B(2,i)-box_min(2))/(box_max(2)-box_min(2))
            zdef = (X(3,i)+Umag*B(3,i))/Z_Length
            IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 ) THEN
               tatom = "Al"
               mass = 13.0
            ENDIF
            IF ( ISRelaxed(i)==-1 ) THEN
               tatom = "Ni"
               mass = 14.0
            ENDIF
 
!C           zdef = b(3,i)
            IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. ISRelaxed(i)==-1 ) THEN
               WRITE (Logic,'(f4.0,1X,A2, 1X, 6f16.11)') mass , tatom , xdef , ydef , zdef , 0.0 , 0.0 , 0.0
            END IF
         ENDIF

         if (key == 'lamp') then
            if (i == 1) then
               if (isrelaxed(i) == -1) then
                  atom_type = 1
               else if (isrelaxed(i) == 2) then
                  atom_type = 3
               else
                  atom_type = 1
               end if                 
            end if
            write(outstr, *) 'ITEM: TIMESTEP'
            write(logic, *) adjustl(trim(outstr))
            
            write(logic, '(I1)') 0

            write(outstr, *) 'ITEM: NUMBER OF ATOMS'
            write(logic, *) adjustl(trim(outstr))

            write(outstr, '(I7)') ntot
            write(logic, *) adjustl(trim(outstr))

            write(outstr, *) 'ITEM: BOX BOUNDS ss ss pp'
            write(logic, *) adjustl(trim(outstr))
            
            write(outstr, '(2F16.9)') box_min(1), box_max(1)
            write(logic, *) adjustl(trim(outstr))
            write(outstr, '(2F16.9)') box_min(2), box_max(2)
            write(logic, *) adjustl(trim(outstr))
            write(outstr, '(2F16.9)') -0.5, 0.5
            write(logic, *) adjustl(trim(outstr))
            
            write(outstr, '(A300)') 'ITEM: ATOMS id type x y z fx fy fz'
            write(logic, *) adjustl(trim(outstr))
            xdef = x(1,i) + b(1,i)
            ydef = x(2,i) + b(2,i)
            zdef = x(3,i) + b(3,i)
            if (isrelaxed(i) /= 0) then 
               write(outstr,'(I7,1X,I2,1X,6(F16.9))') &
                    i, atom_type, xdef, ydef, zdef, dr(1,i), dr(2,i), dr(3,i)
               write(logic, '(A107)') adjustl(trim(outstr))
            end if
            
         end if

        if (key.eq.'olmp') then

          if (i.eq.1) then

            write(logic,'(''ITEM: TIMESTEP'')')
            write(logic,'(i1)') 0
            write(logic,'(''ITEM: NUMBER OF ATOMS'')')
            write(logic,'(i5)') ntot
            write(logic,'(''ITEM: BOX BOUNDS pp pp pp'')')
            write(logic,'(f10.4,'' '',f10.4)') 0.0,& 
     &            box_max(1)-box_min(1)
            write(logic,'(f10.4,'' '',f10.4)') 0.0,& 
     &            box_max(2)-box_min(2)
            write(logic,'(f10.4,'' '',f10.4)') -z_length/2, z_length/2
            write(logic,'(''ITEM: ATOMS type x y z'')')
          end if 


          xdef =(x(1,i)+umag*b(1,i)-box_min(1))
          ydef =(x(2,i)+umag*b(2,i)-box_min(2))
          zdef = (x(3,i) + umag*b(3,i))
          if (IsRelaxed(i) == 1) then 
            TAtom = "1"
          end if
          if (IsRelaxed(i) == 2) then 
            TAtom = "2"
          end if
          if (IsRelaxed(i) == -1) then 
            TAtom = "3"
          end if

          if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 .or. &
     &       IsRelaxed(i) == -1) then
            write(logic,'(A1, 3f16.8)') TAtom,xdef,ydef,zdef
          endif         


        endif         
         
  
         IF ( Key=='stra' .OR. Key=='stre' ) THEN
            IF ( i==1 ) THEN
               WRITE (Logic,&
     &               '('' VARIABLES = X Y Z E11 E22 E33 E23 E13 E12 '')'&
     &               )
               WRITE (Logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''              ,i5,'',&
     & F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
            zdef = X(3,i) + Umag*B(3,i)
            WRITE (Logic,'(9e14.6)') xdef , ydef , zdef , nstrain(1,1,i)&
     &                               , nstrain(2,2,i) , nstrain(3,3,i) ,&
     &                               nstrain(2,3,i) , nstrain(1,3,i) , &
     &                               nstrain(1,2,i)
!            if(i.eq.numpnts) then
!               deallocate(nstrain,avgnum)
!            endif
         ENDIF
 
         IF ( Key=='inte' ) THEN
            IF ( i==1 ) THEN
               WRITE (Logic,'('' VARIABLES = X Y SED'')')
               WRITE (Logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'', J= ''              ,i5,'',&
     & F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
            WRITE (Logic,'(3e15.6)') xdef , ydef
         ENDIF
 
 
 
         IF ( Key=='viri' ) THEN
            IF ( i==1 ) THEN
               WRITE (Logic,&
     &               '('' VARIABLES = X Y Z V11 V22 V33 V23 V13 V12 '')'&
     &               )
               WRITE (Logic,&
     &'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''              ,i5,'',&
     & F = FEPOINT'')') numpnts , numtri
            ENDIF
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
            zdef = X(3,i) + Umag*B(3,i)
            IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 ) THEN
               WRITE (Logic,'(9e14.6)') xdef , ydef , zdef , &
     &                                  AVEvirst(1,1,i) , &
     &                                  AVEvirst(2,2,i) , &
     &                                  AVEvirst(3,3,i) , &
     &                                  AVEvirst(2,3,i) , &
     &                                  AVEvirst(1,3,i) , &
     &                                  AVEvirst(1,2,i)
            ELSE
               WRITE (Logic,'(6e14.6)') xdef , ydef , zdef , &
     &                                  nstrain(1,1,i) , nstrain(2,2,i)&
     &                                  , nstrain(3,3,i) , &
     &                                  nstrain(2,3,i) , nstrain(1,3,i)&
     &                                  , nstrain(1,2,i)
            ENDIF
         ENDIF
 
!c--Now output H position only!
         IF ( Key=='virH' ) THEN
            IF ( ISRelaxed(i)==1 .AND. ATOmspecie(i)==2 ) THEN
               xdef = X(1,i) + Umag*B(1,i)
               ydef = X(2,i) + Umag*B(2,i)
               zdef = X(3,i) + Umag*B(3,i)
 
               WRITE (Logic,'(6e14.6)') xdef , ydef , zdef
!     &            avevirst(1,1,i)
!     $          ,avevirst(2,2,i),avevirst(3,3,i),avevirst(2,3,i)
!     $          ,avevirst(1,3,i),avevirst(1,2,i)
            ENDIF
         ENDIF
 
 
      ENDDO
      DO i = 1 , NUMnp
         IF ( Key=='atom' ) THEN
            IF ( i==1 ) WRITE (Logic,'(''Ni '')')
 
            xdef = (X(1,i)+Umag*B(1,i)+box_min(1))&
     &             /(box_max(1)-box_min(1))
            ydef = (X(2,i)+Umag*B(2,i)+box_max(1))&
     &             /(box_max(2)-box_min(2))
            zdef = (X(3,i)+Umag*B(3,i))/Z_Length
!     C           zdef = b(3,i)
            IF ( ISRelaxed(i)==-1 ) WRITE (Logic,'(4f16.11,1X,2I4)')&
     &           xdef , ydef , zdef , ENErgy(i) , NUMneighbors(i) , 0
         ENDIF
      ENDDO
 
      IF ( (Key/='atom') .AND. (Key/='pdbf') .AND. (Key/='virH') ) THEN
         DO n = 1 , NUMel
            WRITE (Logic,'(4i6)') Ix(1,n) , Ix(2,n) , Ix(3,n) , Ix(3,n)
         ENDDO
      ENDIF
 
      END SUBROUTINE DUMPIT
!*==dumpit_vtk.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!***********************************************************************
      SUBROUTINE DUMPIT_VTK(X,Ix,B,Db,Id,F,Dr,Scale,Logic,Key,Index,&
     &                      Umag)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--DUMPIT_VTK821
!
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*) , Db(NDF,*) , F(NDF,*) , &
     &                 Dr(NDF,*) , Scale , Umag
      INTEGER Id(NDF,*) , Ix(NEN1,*) , Index , Logic , numpnts , ntot
      CHARACTER*4 Key
!
      INTEGER numtri , i , j , n , NDFMAX , ii , jj , nout , npad , &
     &        startpad
      CHARACTER out*100
      PARAMETER (NDFMAX=18)
      DOUBLE PRECISION dev(6) , xdef , ydef , zdef , hyd , sigeff , y , &
     &                 epseff
 
      INTEGER n1 , n2 , n3 , npatoms
      DOUBLE PRECISION det , amat(3,2) , epsloc(3,3) , tarray(3,3)
      DOUBLE PRECISION , ALLOCATABLE :: nstrain(:,:,:)
      DOUBLE PRECISION , ALLOCATABLE :: avgnum(:)
      DOUBLE PRECISION , ALLOCATABLE :: virist(:,:) , nstr1(:,:)
      DOUBLE PRECISION dwx1 , dwx2
      DOUBLE PRECISION box_max(2) , box_min(2)
      DOUBLE PRECISION :: rx1(3) , ud1(3) , b1(3) , b2(3) , b3(3) , &
     &                    ud2(3) , ud3(3) , rx2(3) , rx3(3)
      DOUBLE PRECISION :: sout1(3) , sout2(3) , sout3(3)
!     New Variable for vtk
      INTEGER ndf22
      DOUBLE PRECISION :: ev_convert , fact
      ev_convert = 1.602176462
!     Converts from ev/A^3 to MPa
      fact = 1.D0/(ev_convert)/1.D-5
 
!
      numtri = NUMel
!$$$      if(key.ne.'atom'.and.key.ne.'virH') then
!$$$c$$$        write(logic,'('' TITLE = " '',a4,'' "'')') key
!$$$      endif
 
 
!     if we are doing strains or stresses then calculate them at each no
      IF ( Key=='stra' .OR. Key=='stre' .OR. Key=='viri' ) THEN
 
         ALLOCATE (nstrain(3,3,NUMnp))
         ALLOCATE (nstr1(3,3))
         ALLOCATE (virist(3,3))
         ALLOCATE (avgnum(NUMnp))
         DO i = 1 , NUMnp
            avgnum(i) = 0.0
            DO j = 1 , 3
               DO ii = 1 , 3
                  nstrain(ii,j,i) = 0.0
               ENDDO
            ENDDO
         ENDDO
 
         DO i = 1 , NUMel
            n1 = Ix(1,i)
            n2 = Ix(2,i)
            n3 = Ix(3,i)
            b1(1:3) = B(1:3,n1)
            b2(1:3) = B(1:3,n2)
            b3(1:3) = B(1:3,n3)
!     ****************************************************
!     Adding back displacement contribution to satisfy superposition
            rx1(1:3) = X(1:3,n1)
            rx2(1:3) = X(1:3,n2)
            rx3(1:3) = X(1:3,n3)
 
            CALL DISL_DISPL(rx1,ud1)
            CALL DISL_DISPL(rx2,ud2)
            CALL DISL_DISPL(rx3,ud3)
 
!     Can also call disl_stress to calculate stress directly at the
!     nodal points to make sure that everything is ok
!     Currently the plots show some anamalies due to interpolation and
!     slip planes passing through elements etc.
!     Interpolation from FE displacments is not quite right
            b1 = b1 - ud1
            b2 = b2 - ud2
            b3 = b3 - ud3
!     ****************************************************
            det = (X(1,n1)-X(1,n3))*(X(2,n2)-X(2,n3))&
     &            - (X(1,n2)-X(1,n3))*(X(2,n1)-X(2,n3))
            amat(1,1) = (X(2,n2)-X(2,n3))/det
            amat(2,1) = (X(2,n3)-X(2,n1))/det
            amat(3,1) = (X(2,n1)-X(2,n2))/det
            amat(1,2) = (X(1,n3)-X(1,n2))/det
            amat(2,2) = (X(1,n1)-X(1,n3))/det
            amat(3,2) = (X(1,n2)-X(1,n1))/det
!           ****compute the green strain time 2*****
            CALL GETELEMENTSTRAIN(b1,b2,b3,amat(1:3,1:2),epsloc(1:3,1:3)&
     &                            )
 
!$$$            call GetElementStrain(b1,b2,b3,amat(1:3,1:2),epsloc(1:3
!$$$     $           ,1:3))
 
            epsloc(1:3,1:3) = epsloc(1:3,1:3)*0.5
 
            IF ( Key=='stre' .OR. Key=='viri' ) THEN
!           ***** compute stresses in units of eV/A^3
               CALL GETSTRESSES(epsloc,tarray)
            ELSE
               tarray(1:3,1:3) = epsloc(1:3,1:3)
!***  convert to engineering strains
               tarray(2,3) = tarray(2,3) + tarray(3,2)
               tarray(1,3) = tarray(1,3) + tarray(3,1)
               tarray(1,2) = tarray(1,2) + tarray(2,1)
            ENDIF
 
 
            DO j = 1 , 3
               DO jj = 1 , 3
                  nstrain(j,jj,n1) = nstrain(j,jj,n1) + tarray(j,jj)
                  nstrain(j,jj,n2) = nstrain(j,jj,n2) + tarray(j,jj)
                  nstrain(j,jj,n3) = nstrain(j,jj,n3) + tarray(j,jj)
               ENDDO
            ENDDO
            IF ( Key=='stre' .OR. Key=='viri' ) THEN
!           ***** Compute superposition stress due to dislocations
               sout1 = 0.0D0
               sout2 = 0.0D0
               sout3 = 0.0D0
               CALL DISL_STRESS(rx1,sout1)
               CALL DISL_STRESS(rx2,sout2)
               CALL DISL_STRESS(rx3,sout3)
 
               nstrain(1,1,n1) = nstrain(1,1,n1) + sout1(1)
               nstrain(1,2,n1) = nstrain(1,2,n1) + sout1(3)
               nstrain(2,1,n1) = nstrain(1,2,n1) + sout1(3)
               nstrain(2,2,n1) = nstrain(2,2,n1) + sout1(2)
 
               nstrain(1,1,n2) = nstrain(1,1,n2) + sout2(1)
               nstrain(1,2,n2) = nstrain(1,2,n2) + sout2(3)
               nstrain(2,1,n2) = nstrain(1,2,n2) + sout2(3)
               nstrain(2,2,n2) = nstrain(2,2,n2) + sout2(2)
 
               nstrain(1,1,n3) = nstrain(1,1,n3) + sout3(1)
               nstrain(1,2,n3) = nstrain(1,2,n3) + sout3(3)
               nstrain(2,1,n3) = nstrain(1,2,n3) + sout3(3)
               nstrain(2,2,n3) = nstrain(2,2,n3) + sout3(2)
            ENDIF
 
!$$$            call disl_stress(rx1, ud1)
!$$$            call disl_stress(rx2, ud2)
!$$$            call disl_stress(rx3, ud3)
!$$$
!$$$            do j = 1,3
!$$$               do jj = 1, 3
!$$$                  if (j==1 .and. jj==1) then
!$$$                     nstrain(j,j,n1) = nstrain(j,j,n1) + ud1(1)
!$$$                     nstrain(j,j,n2) = nstrain(j,j,n2) + ud2(1)
!$$$                     nstrain(j,j,n3) = nstrain(j,j,n3) + ud3(1)
!$$$                  end if
!$$$                  if (j==1 .and. jj==2) then
!$$$                     nstrain(j,j,n1) = nstrain(j,j,n1) + ud1(1)
!$$$                     nstrain(j,j,n2) = nstrain(j,j,n2) + ud2(1)
!$$$                     nstrain(j,j,n3) = nstrain(j,j,n3) + ud3(1)
!$$$                  end if
!$$$
!$$$               end do
!$$$            end do
            avgnum(n1) = avgnum(n1) + 1.0
            avgnum(n2) = avgnum(n2) + 1.0
            avgnum(n3) = avgnum(n3) + 1.0
         ENDDO
 
         DO i = 1 , NUMnp
            IF ( avgnum(i)>0 ) THEN
               DO j = 1 , 3
                  DO jj = 1 , 3
                     nstrain(j,jj,i) = nstrain(j,jj,i)/avgnum(i)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
 
      ENDIF
 
 
!     find max and min of coordintes so that a box can be made.
      IF ( Key=='atom' ) THEN
         DO j = 1 , 2
            box_max(j) = -1.0D30
            box_min(j) = 1.0D30
         ENDDO
         npatoms = 0
         npad = 0
         ntot = 0
         DO i = 1 , NUMnp
            IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 .OR. ISRelaxed(i)&
     &           ==-1 ) THEN
               IF ( ISRelaxed(i)==-1 ) THEN
                  IF ( npad==0 ) startpad = i
                  npad = npad + 1
               ELSE
                  npatoms = npatoms + 1
               ENDIF
               ntot = ntot + 1
               DO j = 1 , 2
                  IF ( X(j,i)+Umag*B(j,i)>box_max(j) ) box_max(j)&
     &                 = X(j,i) + Umag*B(j,i)
                  IF ( X(j,i)+Umag*B(j,i)<box_min(j) ) box_min(j)&
     &                 = X(j,i) + Umag*B(j,i)
               ENDDO
            ENDIF
         ENDDO
         startpad = npatoms - npad + 1
         PRINT * , 'Start pad atoms = ' , startpad
      ENDIF
 
!     begin loop over nodes for all cases
!      numpnts=max(numnp,numnpp1)
      numpnts = MAX(NUMnp,NUMnpp1)
      IF ( Key=='stra' .OR. Key=='stre' .OR. Key=='viri' ) THEN
         WRITE (Logic,FMT='(A)') '# vtk DataFile Version 2.0'
         IF ( Key=='stra' ) WRITE (Logic,FMT='(A)') 'Strains from CADD'
         IF ( Key=='stre' ) WRITE (Logic,FMT='(A)') 'Stresses from CADD'
         IF ( Key=='viri' ) WRITE (Logic,FMT='(A)')&
     &                              'Virial Stresses from CADD'
         WRITE (Logic,FMT='(A)') 'ASCII'
         WRITE (Logic,FMT='(A)') 'DATASET UNSTRUCTURED_GRID'
         WRITE (Logic,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , numpnts , &
     &          'float'
         DO i = 1 , numpnts
 
            xdef = X(1,i) + Umag*B(1,i)
            ydef = X(2,i) + Umag*B(2,i)
            zdef = X(3,i) + Umag*B(3,i)
            WRITE (Logic,'(3(1x,e14.6))') xdef , ydef , zdef
         ENDDO
         WRITE (Logic,*)
         WRITE (Logic,'(A5,1X,I7,1X,I7)') 'CELLS' , NUMel , 4*numtri
 
         DO n = 1 , NUMel
            WRITE (Logic,'(4i6)') 3 , Ix(1,n) - 1 , Ix(2,n) - 1 , &
     &                            Ix(3,n) - 1
         ENDDO
         WRITE (Logic,*)
 
         WRITE (Logic,'(A10,1X,I7)') 'CELL_TYPES' , NUMel
 
         DO n = 1 , NUMel
            WRITE (Logic,FMT='(5(1x,I7))') 5
         ENDDO
 
         WRITE (Logic,*) 'POINT_DATA' , numpnts
 
         IF ( Key=='stra' ) THEN
            WRITE (Logic,*) 'SCALARS strain float 3'
         ELSE
            WRITE (Logic,*) 'TENSORS stress float'
         ENDIF
 
!$$$         write(logic, *) 'LOOKUP_TABLE default'
 
         DO i = 1 , numpnts
!$$$               , nstrain(2,3,i), nstrain(1
!$$$     $              ,3,i), nstrain(1,2,i)
            IF ( Key=='stra' ) WRITE (Logic,'(6(1x,e14.6))')&
     &                                nstrain(1,1,i) , nstrain(2,2,i) , &
     &                                nstrain(1,2,i)
            virist(:,:) = AVEvirst(:,:,i)*fact
            nstr1(:,:) = nstrain(:,:,i)*fact
            IF ( Key=='viri' ) THEN
               IF ( ISRelaxed(i)==1 .OR. ISRelaxed(i)==2 ) THEN
                  WRITE (Logic,'(3e14.6)') virist(1,1) , virist(1,2) , &
     &                   virist(1,3)
                  WRITE (Logic,'(3e14.6)') virist(1,2) , virist(2,2) , &
     &                   virist(2,3)
                  WRITE (Logic,'(3e14.6)') virist(1,2) , virist(2,2) , &
     &                   virist(2,3)
                  WRITE (Logic,*)
               ELSE
                  WRITE (Logic,'(3e14.6)') nstr1(1,1) , nstr1(1,2) , &
     &                   nstr1(1,3)
                  WRITE (Logic,'(3e14.6)') nstr1(1,2) , nstr1(2,2) , &
     &                   nstr1(2,3)
                  WRITE (Logic,'(3e14.6)') nstr1(1,3) , nstr1(2,3) , &
     &                   nstr1(3,3)
                  WRITE (Logic,*)
               ENDIF
            ENDIF
         ENDDO
         WRITE (Logic,*) 'VECTORS F float'
         DO i = 1 , numpnts
            WRITE (Logic,'(3e14.6)') DBLE(Id(1,i)) , DBLE(Id(2,i)) , 0.0
         ENDDO
         WRITE (Logic,*) 'CELL_DATA' , NUMel
         WRITE (Logic,*) 'SCALARS DB integer'
         WRITE (Logic,*) 'LOOKUP_TABLE default'
 
         DO n = 1 , NUMel
            WRITE (Logic,'(i5)') Ix(NEN1,n)
         ENDDO
 
 
      ENDIF
      END SUBROUTINE DUMPIT_VTK
!*==dump_mesh.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***************************************************************
      SUBROUTINE DUMP_MESH(X,B,Ix,Logic,Iq)
      USE MOD_GLOBAL
      IMPLICIT NONE
!*--DUMP_MESH1123
!
      DOUBLE PRECISION X(NXDm,*) , B(NDF,*) , scale , umag
      INTEGER Ix(NEN1,*) , index , Logic , numpnts , ntot
      INTEGER , OPTIONAL :: Iq
      CHARACTER*4 key
!
      INTEGER numtri , i , j , n , NDFMAX , ii , jj , nout , npad , &
     &        startpad
      CHARACTER out*100
      PARAMETER (NDFMAX=18)
      DOUBLE PRECISION dev(6) , xdef , ydef , zdef , hyd , sigeff , y , &
     &                 epseff
 
      INTEGER n1 , n2 , n3 , npatoms
      DOUBLE PRECISION det , amat(3,2) , epsloc(3,3) , tarray(3,3)
      DOUBLE PRECISION , ALLOCATABLE :: nstrain(:,:,:)
      DOUBLE PRECISION , ALLOCATABLE :: avgnum(:)
      DOUBLE PRECISION dwx1 , dwx2
      DOUBLE PRECISION box_max(2) , box_min(2)
      DOUBLE PRECISION :: rx1(3) , ud1(3) , b1(3) , b2(3) , b3(3) , &
     &                    ud2(3) , ud3(3) , rx2(3) , rx3(3)
 
      umag = 1.0D0
      numpnts = NUMnp
      numtri = NUMel
 
      WRITE (Logic,FMT='(A)') '# vtk DataFile Version 2.0'
      WRITE (Logic,FMT='(A)') 'Deformed mesh from CADD'
      WRITE (Logic,FMT='(A)') 'ASCII'
      WRITE (Logic,FMT='(A)') 'DATASET UNSTRUCTURED_GRID'
      WRITE (Logic,FMT='(A6,1x,I7,1x,A5)') 'POINTS' , numpnts , 'float'
      DO i = 1 , numpnts
         IF ( PRESENT(Iq) ) THEN
            xdef = X(1,i) - XTIp_actual(1) + umag*B(1,i)
            ydef = X(2,i) - XTIp_actual(2) + umag*B(2,i)
            zdef = X(3,i) + umag*B(3,i)
         ELSE
            xdef = X(1,i) + umag*B(1,i)
            ydef = X(2,i) + umag*B(2,i)
            zdef = X(3,i) + umag*B(3,i)
         ENDIF
         WRITE (Logic,'(3(1x,e14.6))') xdef , ydef , zdef
      ENDDO
      WRITE (Logic,*)
      WRITE (Logic,'(A5,1X,I7,1X,I7)') 'CELLS' , NUMel , 4*numtri
 
      DO n = 1 , NUMel
         WRITE (Logic,'(4i6)') 3 , Ix(1,n) - 1 , Ix(2,n) - 1 , Ix(3,n)&
     &                         - 1
      ENDDO
      WRITE (Logic,*)
 
      WRITE (Logic,'(A10,1X,I7)') 'CELL_TYPES' , NUMel
 
      DO n = 1 , NUMel
         WRITE (Logic,FMT='(5(1x,I7))') 5
      ENDDO
      END SUBROUTINE DUMP_MESH
 
!***************************************************************
