!*==pmacr.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
      SUBROUTINE PMACR(Id,X,Ix,F,B,Dr,Db,Itx,lmp)
      USE MOD_FILE
      USE MOD_GLOBAL
      USE MOD_BOUNDARY
      USE MOD_DYNAMO
      USE MOD_GRAIN
      USE MOD_MATERIAL
      use lammps
      use mod_lammps
      IMPLICIT NONE

!*--PMACR10
      !

      type(c_ptr) :: lmp
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*) , F(*) , Dr(NDF,*) , Db(*)
      INTEGER Id(*) , Ix(NEN1,*)
      INTEGER Itx(3,*)
 
!
!---- macro instruction interpreter
!---- controls solution algorithms
!
      INTEGER NWD
      PARAMETER (NWD=56)
      LOGICAL PCOMP , check , BOUndary , fladapt , statrue
      CHARACTER*80 input , filename , ctemp
      CHARACTER*4 key
      REAL*4 wd(NWD)
      INTEGER ICRit , ISIgn , NODecrit , IDFcrit
      DOUBLE PRECISION ARCrate
      COMMON /CSTEP / BOUndary , ICRit , ARCrate , ISIgn , NODecrit , &
     &                IDFcrit
 
      INTEGER lvs(9) , lve(9) , neqmax , i , j , k , l , lv , lx , ll , &
     &        lmax , ii , lower , upper , NEXT , ldif , i1 , i2 , &
     &        loops , idum , n1 , n2 , n , logic , l0 , numelsve , &
     &        loopshow , nloop , ngt0
      DOUBLE PRECISION dum , bnew , varstp , dtstep , tmax , dttol , &
     &                 dt0 , qrot(3,3) , cc2(6,6) , strainenergy , &
     &                 dtstp_orig
      LOGICAL flag02
      REAL*4 , POINTER :: ct(:,:) , ctpass(:,:)
      INTEGER , ALLOCATABLE :: itemp1(:) , itemp2(:)
      DATA wd/'tole' , 'dtim' , 'dump' , 'stre' , 'xxxx' , 'xxxx' , &
     &     'loop' , 'next' , 'xxxx' , 'chec' , 'time' , 'xxxx' , &
     &     'chti' , 'xxxx' , 'xxxx' , 'mesh' , 'xxxx' , 'end ' , &
     &     'pdel' , 'xxxx' , 'xxxx' , 'xxxx' , 'prec' , 'ddse' , &
     &     'xxxx' , 'newd' , 'ma05' , 'ma06' , 'getc' , 'xxxx' , &
     &     'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , 'rest' , &
     &     'xxxx' , 'ma01' , 'ma02' , 'ma03' , 'xxxx' , 'xxxx' , &
     &     'clea' , 'xxxx' , 'stat' , 'xxxx' , 'xxxx' , 'xxxx' , &
     &     'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , 'xxxx' , &
     &     'xxxx' , 'xxxx'/
!
!---- set initial values of parameters
      neqmax = MAXnp*NDF
      BOUndary = .TRUE.
      ICRit = 0
      NITer = 0
      RNMax = 0.D0
      TIMeol = 0.D0
      TIMe = 0.D0
      flag02 = .TRUE.
      CONvergetol = 1.D-9
!---- read macro program
      WRITE (6,*)
      WRITE (6,&
     &'('' '',a80//5x,''macro instructions''//1x,                ''nesti&
     &ng  statement'',10x,''variables'')') HEAd
      WRITE (6,'(1x,''-------  ---------'',10x,''---------'')')
      nloop = 0
      loopshow = 0
      ll = 1
      lmax = 16
      ALLOCATE (ct(20,lmax),ctpass(20,lmax))
      ct(1,1) = wd(7)
      ct(3,1) = 1.0
 100  ll = ll + 1
      IF ( ll>=lmax ) THEN
         lmax = lmax + 16
         ctpass = ct
         DEALLOCATE (ct)
         ALLOCATE (ct(20,lmax))
         ct(1:20,1:lmax) = ctpass
         DEALLOCATE (ctpass)
         ALLOCATE (ctpass(20,lmax))
      ENDIF
      DO
         READ (input_file_unit,'(a80)') input
!---- added to allow indentation with spaces in input file
         i = 1
         DO WHILE ( (input(i:i)==' ') .AND. (i<=80) )
            i = i + 1
         ENDDO
         IF ( i<=80 ) THEN
            IF ( input(i:i)=='%' ) CYCLE
         ENDIF
         EXIT
      ENDDO
      upper = i - 1
      DO ii = 1 , 2
         lower = upper
         upper = NEXT(lower,input)
         ldif = upper - lower - 1
         IF ( ldif==0 ) THEN
            ct(ii,ll) = '    '
         ELSE
            CALL EQUAL(ct(ii,ll),input,lower+1,upper-1)
         ENDIF
      ENDDO
      IF ( PCOMP(ct(1,ll),wd(7)) ) THEN
!-- added to show loop nesting
         nloop = nloop + 1
         loopshow = loopshow + nloop*10**(nloop-1)
         i1 = upper
         i2 = NEXT(i1,input)
         CALL FREEIN(input,i1,i2,loops,dum,1)
         ct(3,ll) = loops
      ELSE
         CALL EQUAL(ct(3,ll),input,upper+1,80)
      ENDIF
      IF ( (loopshow>0) ) THEN
         WRITE (6,'(i8,2x,a4,1x,a4,9x,a40)') loopshow , (ct(j,ll),j=1,2)&
     &          , input(upper+1:)
      ELSE
         WRITE (6,'(10x,a4,1x,a4,9x,a40)') (ct(j,ll),j=1,2) , &
     &          input(upper+1:)
      ENDIF
!--- added to show loop nesting
      IF ( PCOMP(ct(1,ll),wd(8)) ) THEN
         loopshow = loopshow - nloop*10**(nloop-1)
         nloop = nloop - 1
      ENDIF
      IF ( .NOT.PCOMP(ct(1,ll),wd(18)) ) GOTO 100
      ct(1,ll) = wd(8)
!---- set loop markers
      lx = ll - 1
      DO l = 1 , lx
         IF ( .NOT.PCOMP(ct(1,l),wd(7)) ) CYCLE
         j = 1
         k = l + 1
         DO i = k , ll
            IF ( PCOMP(ct(1,i),wd(7)) ) j = j + 1
            IF ( j>9 ) THEN
               WRITE (6,'(/'' **error** loops nested deeper than 8'')')
               STOP
            ENDIF
            IF ( PCOMP(ct(1,i),wd(8)) ) j = j - 1
            IF ( j==0 ) GOTO 150
         ENDDO
         WRITE (6,'(/'' **error** unbalanced loop/next macros'')')
         STOP
 150     ct(4,i) = l
         ct(4,l) = i
      ENDDO
      j = 0
      DO l = 1 , ll
         IF ( PCOMP(ct(1,l),wd(7)) ) j = j + 1
         IF ( PCOMP(ct(1,l),wd(8)) ) j = j - 1
      ENDDO
      IF ( j/=0 ) THEN
         WRITE (6,'(/'' **error** unbalanced loop/next macros'')')
         STOP
      ENDIF
!---- execute macro instruction program
      lv = 0
      l = 1
 200  DO j = 1 , NWD
         IF ( PCOMP(ct(1,l),wd(j)) ) GOTO 300
      ENDDO
      GOTO 500
 300  i = l - 1
      IF ( l/=1 .AND. l/=ll ) THEN
         WRITE (6,&
     &'(2x,''**macro instruction'',i4,'' ** '',                  2(a4,2x&
     &))') i , (ct(k,l),k=1,2)
         CALL FLUSH(6)
      ENDIF
      IF ( j<=24 ) THEN
         SELECT CASE (j)
         CASE (2)
!
!---- macro 'dtim'
!---- set time increment
            lower = 0
            upper = NEXT(lower,ct(3,l))
            CALL FREEIN(ct(3,l),lower,upper,idum,DT,2)
            WRITE (6,'(2x,''**time step set to '',e15.5)') DT
            GOTO 500
         CASE (3)
!
!---- macro 'dump'
            CALL LEODUMP(NUMnp,NUMel,NDF,NXDm,NEN1,X,Ix,Id,ISRelaxed,F,&
     &                   Itx)
            GOTO 500
         CASE (4)
!
!---- macro 'stre'
            CALL HSTRESS(NUMnp,ISRelaxed,AVEvirst)
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            GOTO 500
         CASE (5)
            GOTO 500
         CASE (6)
            GOTO 500
         CASE (7)
!
!---- macro 'loop'
!---- set loop start indicators
            lv = lv + 1
            lx = ct(4,l)
            lvs(lv) = l
            lve(lv) = lx
            ct(3,lx) = 1.
            GOTO 500
         CASE (8)
!
!---- macro 'next'
!---- loop terminator control
            n = ct(4,l)
            ct(3,l) = ct(3,l) + 1.0
            IF ( ct(3,l)>ct(3,n) ) lv = lv - 1
            IF ( ct(3,l)<=ct(3,n) ) l = n
!
!---- macro 'xxxx'
            GOTO 500
         CASE (9)
            GOTO 500
         CASE (10)
!
!---- macro 'chec'
!
            WRITE (*,*) " chec no longer in this version of CADD "
!      call derivcheck(id,b,x,ix,f,dr,itx)
            STOP
         CASE (11)
!
!---- macro 'time'
!---- increment time
            TIMeol = TIMe
            upper = 0
            dtstep = DT
            MOVED = .FALSE.
            IF ( dtstep/=0.D0 ) B0(1:NDF,1:NUMnp) = B(1:NDF,1:NUMnp)
!!$            IF ( MOVed ) THEN
!!$               dtstep = 0.0D0
!!$               IF ( XTIp_init(1)<0.0D0 ) THEN
!!$                  CALL GET_CRACK_TIP(X,B)
!!$                  XTIp_init(1:2) = XTIp(1:2)
!!$               ENDIF
!!$               PRINT * , &
!!$     &             'Moving mesh, repeat calculation without load change'
!!$               MOVed = .FALSE.
!!$!     Reset moved flag to .false.
!!$            ENDIF
            TIMe = TIMe + dtstep
            Db(1:MAXnp*NDF) = 0.D0
            WRITE (6,'(2x,''**time set to '',e15.5)') TIMe
            lower = upper
            upper = NEXT(lower,ct(3,l))
            CALL FREEIN(ct(3,l),lower,upper,idum,tmax,2)
            IF ( tmax>0.0 ) THEN
               IF ( TIMe>tmax ) THEN
                  WRITE (6,'(2x,''**maximum time attained**'')')
                  WRITE (6,'(2x,''**end of macro execution**'')')
                  RETURN
               ENDIF
            ENDIF
            RNMax = 0.0D0
            NITer = 0
!
!---- macro 'xxxx'
            GOTO 500
         CASE (12)
            GOTO 500
         CASE (13)
!
!---- macro 'chti'
            CALL CHECKTIME
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            GOTO 500
         CASE (14)
            GOTO 500
         CASE (15)
            GOTO 500
         CASE (16)
 
!
!---- macro 'mesh'
!---- modify mesh data
            CALL PMESH(Id,X,Ix,F,B,Dr,Itx,lmp)
 
!
!---- macro 'xxxx'
!----
!
!---- macro 'xxxx'
            GOTO 500
         CASE (17)
            GOTO 500
         CASE (18)
            GOTO 500
         CASE (19)
!
!---- macro 'pdel'
            lower = 0
            upper = NEXT(lower,ct(3,l))
            CALL EQUAL(filename,ct(3,l),lower+1,upper-1)
            CALL IOFILE(filename,'formatted  ',logic,.TRUE.)
            CALL PDEL(F,Dr,Id,logic,NDF,X,TIMe)
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            GOTO 500
         CASE (20)
            GOTO 500
         CASE (21)
            GOTO 500
         CASE (22)
            GOTO 500
         CASE (23)
!
!---- macro 'prec'
!---- initializes the displacement to the K-field
            CALL PRECRACK(Id,X,B,F,ct(2,l))
            call update_lammps_coords(X, B, .true., .true., lmp)
            call lammps_command(lmp, "run 0 pre yes post no")
            GOTO 500
         CASE (24)
!
!---- macro 'ddse'
!---- set up the discrete dislocations solving routine
            IF ( DD_set ) RETURN
            WRITE (*,*)
            WRITE (*,*) 'setting up the discrete dislocation solver'
            WRITE (*,*)
            CALL BANDNL(Id,X,Ix,F,B)
            call map_from_lammps(x, id,lmp)
            qrot = 0.D0
            IF ( NGRains/=1 ) STOP 'hardwired here for 1 grain'
            qrot(1,1) = GRAins(1)%ROTMAT(1)
            qrot(2,2) = GRAins(1)%ROTMAT(1)
            qrot(1,2) = -GRAins(1)%ROTMAT(2)
            qrot(2,1) = GRAins(1)%ROTMAT(2)
            qrot(3,3) = 1.D0
            qrot = MATMUL(qrot,GRAins(1)%XLATVECT(1:3,1:3))
            CALL ROTATEVOIGT(qrot,MATerial(1)%CC,cc2)
            WRITE (*,*)
            WRITE (*,*) 'Rotated elastic matrix:'
            DO i = 1 , 6
               WRITE (*,'(6e13.4)') cc2(i,1:6)
            ENDDO
            ALLOCATE (itemp1(NUMnp),itemp2(NUMel))
            CALL FEM_SETUP(NUMnp,NUMel,X,Id,ISRelaxed,Ix,Itx,cc2,itemp1,&
     &                     itemp2)
 
            DEALLOCATE (itemp1,itemp2)
            CALL DISL_SETUP
            ! ---- Updates any initial config of atoms
            ! ---- What i think is happening is that the K_field is not in the correc
            !      units for the atoms
            call update_lammps_coords(X, B, .true., .true., lmp)
            call lammps_command(lmp, "run 0 pre yes post no")
            DD_set = .TRUE.
!
!---- macro 'xxxx'
            GOTO 500
         CASE DEFAULT
         END SELECT
      ELSEIF ( j<=48 ) THEN
         SELECT CASE ((j-24))
         CASE (1)
            GOTO 500
         CASE (2)
!
!---- macro 'newd'
!---  add a new discrete dislocation
            PRINT * , 'Entering new dislocation'
            CALL NEWDISLOCATION(ct(2,l),X,B,ISRelaxed,NUMnp,NXDm,NDF)
            GOTO 500
         CASE (3)
!
!---- macro 'ma05'
!---- user supplied macro
            IF ( .NOT.DD_set ) STOP 'ERROR: must initialize D.D.'
            ngt0 = NGTlst
 
            PRINT * , 'did you mean ma06?'
            STOP
 
            WRITE (*,*) 'NUMBER NEIGHBOR UPDATES:' , NGTlst - ngt0
            GOTO 500
         CASE (4)
!
!---- macro 'xxxx'
!---- user supplied macro to call MD/DD routine
 
            IF ( .NOT.DD_set ) STOP 'ERROR: must initialize D.D.'
            ngt0 = NGTlst
 
            CALL MA06(Id,X,Ix,F,B,Dr,Db,ct(2,l),Itx, lmp)
            WRITE (*,*) 'NUMBER NEIGHBOR UPDATES:' , NGTlst - ngt0
!!$            CALL GET_CRACK_TIP(X,B)
!$$$      if (abs(xtip(1)-xtip_init(1)) > x_move_mesh) then
!$$$         do j = 1, 2
!$$$            xtip_actual(j)=abs(xtip(j)-xtip_init(j))
!$$$            if (xtip(j) > xtip_init(j)) then
!$$$               x_tip_dir(j) = 1.0d0
!$$$            else
!$$$               x_tip_dir(j) = -1.0d0
!$$$            end if
!$$$         end do
!$$$         call move_atomistic_crack(x,ix,b)
!$$$      end if
 
            GOTO 500
         CASE (5)
!
!---- macro 'getc'
!---- macro to set the initial crack tip position
            CALL GET_CRACK_TIP(X,B)
            XTIp_init(1:2) = XTIp(1:2)
            CRAck_motion = 0.0
            MOVemesh = .TRUE.
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            GOTO 500
         CASE (6)
            GOTO 500
         CASE (7)
            GOTO 500
         CASE (8)
            GOTO 500
         CASE (9)
            GOTO 500
         CASE (10)
            GOTO 500
         CASE (11)
            GOTO 500
         CASE (12)
!
!---- macro 'rest'
!---- read/write restart files
            WRITE (*,*) " (rest macro) not in this version! "
!
!---- macro 'xxxx'
            STOP
         CASE (13)
            GOTO 500
         CASE (14)
!
!---- macro 'ma01'
!---- user supplied macro - boundary conditions/loading
            CALL MA01(Id,X,Ix,F,B,Dr,Db,ct(2,l))
            GOTO 500
         CASE (15)
!
!---- macro 'ma02'
!---- user supplied macro - tecplot files
            CALL MA02(Id,X,Ix,F,B,Dr,Db,ct(2,l),flag02)
            GOTO 500
         CASE (16)
!
!---- macro 'ma03'
!---- user supplied macro - energy calculation
            WRITE (*,*) "ma03 not in this version of CADD"
 
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            STOP
         CASE (17)
            GOTO 500
         CASE (18)
            GOTO 500
         CASE (19)
!
!---- macro 'clea'
!---- Reinitialize variables within same time step for a different
!---- solution procedure.
            B(1:NDF,1:NUMnp) = B0(1:NDF,1:NUMnp)
            Db(1:MAXnp*NDF) = 0.D0
            WRITE (6,*) ' **Displacement increments initialized'
            RNMax = 0.0D0
            NITer = 0
!
!---- macro 'xxxx'
            GOTO 500
         CASE (20)
            GOTO 500
         CASE (21)
!
!---- macro 'stat'
!---- Recompute nonlocal element status
            CALL LATTICECHECK(X)
            CALL STATUSCALC(X,Ix,.FALSE.)
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
!
!---- macro 'xxxx'
            GOTO 500
         CASE (22)
            GOTO 500
         CASE (23)
            GOTO 500
         CASE (24)
            GOTO 500
         CASE DEFAULT
         END SELECT
      ELSE
         SELECT CASE ((j-48))
         CASE (1)
         CASE (2)
         CASE (3)
         CASE (4)
         CASE (5)
         CASE (6)
         CASE (7)
         CASE (8)
         CASE DEFAULT
            GOTO 400
         END SELECT
         GOTO 500
      ENDIF
!
!---- macro 'tole'
!---- set solution tolerance
 400  lower = 0
      upper = NEXT(lower,ct(3,l))
      CALL FREEIN(ct(3,l),lower,upper,idum,CONvergetol,2)
      WRITE (6,'(2x,''**tolerance set to '',e15.5)') CONvergetol
!
 500  l = l + 1
      IF ( l<=ll ) GOTO 200
      WRITE (6,'(2x,''**end of macro execution**'')')
      END SUBROUTINE PMACR
!*==equal.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE EQUAL(Char1,Char2,L1,L2)
      IMPLICIT NONE
!*--EQUAL568
!*** Start of declarations inserted by SPAG
      INTEGER L1 , L2
!*** End of declarations inserted by SPAG
!
      CHARACTER*80 Char1
      CHARACTER*80 Char2
!
      Char1 = Char2(L1:L2)
      END SUBROUTINE EQUAL
!*==keyequal.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
      SUBROUTINE KEYEQUAL(Char1,Char2,L1,L2)
      IMPLICIT NONE
!*--KEYEQUAL583
!*** Start of declarations inserted by SPAG
      INTEGER L1 , L2
!*** End of declarations inserted by SPAG
!
      CHARACTER*4 Char1
      CHARACTER*80 Char2
!
      Char1 = Char2(L1:L2)
      END SUBROUTINE KEYEQUAL
!*==rotatevoigt.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
!***********************************************************************
      SUBROUTINE ROTATEVOIGT(Q,C,C2)
      IMPLICIT NONE
!*--ROTATEVOIGT600
      DOUBLE PRECISION Q(3,3) , C(6,6) , C2(6,6) , cc(3,3,3,3) , &
     &                 cc2(3,3,3,3)
      INTEGER i1 , i2 , i3 , i4 , j1 , j2 , j3 , j4
      CALL UNVOIGT(C,cc)
      DO i1 = 1 , 3
         DO i2 = 1 , 3
            DO i3 = 1 , 3
               DO i4 = 1 , 3
                  cc2(i1,i2,i3,i4) = 0.
                  DO j1 = 1 , 3
                     DO j2 = 1 , 3
                        DO j3 = 1 , 3
                           DO j4 = 1 , 3
                              cc2(i1,i2,i3,i4) = cc2(i1,i2,i3,i4)&
     &                           + Q(i1,j1)*Q(i2,j2)*Q(i3,j3)*Q(i4,j4)&
     &                           *cc(j1,j2,j3,j4)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL VOIGT(cc2,C2)
      END SUBROUTINE ROTATEVOIGT
!*==voigt.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE VOIGT(C,Cv)
      IMPLICIT NONE
!*--VOIGT631
      DOUBLE PRECISION C(3,3,3,3) , Cv(6,6)
      INTEGER i , j
      Cv(1,1) = C(1,1,1,1)
      Cv(1,2) = C(1,1,2,2)
      Cv(1,3) = C(1,1,3,3)
      Cv(1,4) = C(1,1,2,3)
      Cv(1,5) = C(1,1,1,3)
      Cv(1,6) = C(1,1,1,2)
      Cv(2,2) = C(2,2,2,2)
      Cv(2,3) = C(2,2,3,3)
      Cv(2,4) = C(2,2,2,3)
      Cv(2,5) = C(2,2,1,3)
      Cv(2,6) = C(2,2,1,2)
      Cv(3,3) = C(3,3,3,3)
      Cv(3,4) = C(3,3,2,3)
      Cv(3,5) = C(3,3,1,3)
      Cv(3,6) = C(3,3,1,2)
      Cv(4,4) = C(2,3,2,3)
      Cv(4,5) = C(2,3,1,3)
      Cv(4,6) = C(2,3,1,2)
      Cv(5,5) = C(1,3,1,3)
      Cv(5,6) = C(1,3,1,2)
      Cv(6,6) = C(1,2,1,2)
      DO i = 1 , 6
         DO j = i , 6
            Cv(j,i) = Cv(i,j)
         ENDDO
      ENDDO
      END SUBROUTINE VOIGT
!*==unvoigt.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE UNVOIGT(Cv,C)
      IMPLICIT NONE
!*--UNVOIGT665
      DOUBLE PRECISION C(3,3,3,3) , Cv(6,6)
      INTEGER i1 , i2 , j1 , j2 , j3 , j4
      DO i1 = 1 , 6
         IF ( i1==1 ) THEN
            j1 = 1
            j2 = 1
         ELSEIF ( i1==2 ) THEN
            j1 = 2
            j2 = 2
         ELSEIF ( i1==3 ) THEN
            j1 = 3
            j2 = 3
         ELSEIF ( i1==4 ) THEN
            j1 = 2
            j2 = 3
         ELSEIF ( i1==5 ) THEN
            j1 = 1
            j2 = 3
         ELSEIF ( i1==6 ) THEN
            j1 = 1
            j2 = 2
         ENDIF
         DO i2 = i1 , 6
            IF ( i2==1 ) THEN
               j3 = 1
               j4 = 1
            ELSEIF ( i2==2 ) THEN
               j3 = 2
               j4 = 2
            ELSEIF ( i2==3 ) THEN
               j3 = 3
               j4 = 3
            ELSEIF ( i2==4 ) THEN
               j3 = 2
               j4 = 3
            ELSEIF ( i2==5 ) THEN
               j3 = 1
               j4 = 3
            ELSEIF ( i2==6 ) THEN
               j3 = 1
               j4 = 2
            ENDIF
            C(j1,j2,j3,j4) = Cv(i1,i2)
            C(j2,j1,j3,j4) = Cv(i1,i2)
            C(j1,j2,j4,j3) = Cv(i1,i2)
            C(j2,j1,j4,j3) = Cv(i1,i2)
            C(j3,j4,j1,j2) = Cv(i1,i2)
            C(j3,j4,j2,j1) = Cv(i1,i2)
            C(j4,j3,j1,j2) = Cv(i1,i2)
            C(j4,j3,j2,j1) = Cv(i1,i2)
         ENDDO
      ENDDO
      END SUBROUTINE UNVOIGT
!*==leodump.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
      SUBROUTINE LEODUMP(Numnp,Numel,Ndf,Nxdm,Nen1,X,Ix,Id,Isrelaxed,F,&
     &                   Itx)
      IMPLICIT NONE
!*--LEODUMP724
      INTEGER Numnp , Numel , Ndf , Nxdm , Nen1 , i
      INTEGER Ix(Nen1,Numel) , Id(Ndf,Numnp) , Isrelaxed(Numnp) , &
     &        Itx(3,Numel)
      DOUBLE PRECISION X(Nxdm,Numnp) , b(Ndf,Numnp) , F(Ndf,Numnp)
      WRITE (*,*) Numnp , Numel
      WRITE (*,*) 'x:'
      DO i = 1 , Numnp
         WRITE (*,*) X(1:3,i)
      ENDDO
      WRITE (*,*) 'f:'
      DO i = 1 , Numnp
         WRITE (*,*) F(1:3,i)
      ENDDO
      WRITE (*,*) 'id:'
      DO i = 1 , Numnp
         WRITE (*,*) Id(1:3,i)
      ENDDO
      WRITE (*,*) 'IsRelaxed:'
      DO i = 1 , Numnp
         WRITE (*,*) Isrelaxed(i)
      ENDDO
      WRITE (*,*) 'ix:'
      DO i = 1 , Numel
         WRITE (*,*) Ix(1:4,i)
      ENDDO
      WRITE (*,*) 'itx:'
      DO i = 1 , Numel
         WRITE (*,*) Itx(1:3,i)
      ENDDO
 
 
      END SUBROUTINE LEODUMP
!*==hstress.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE HSTRESS(Numnp,Isrelaxed,Virst)
      IMPLICIT NONE
!*--HSTRESS763
      INTEGER i , j , k , natoms
      INTEGER Numnp
      INTEGER Isrelaxed(Numnp)
      DOUBLE PRECISION Virst(3,3,Numnp)
      DOUBLE PRECISION hydro_stress
      DOUBLE PRECISION avg_stress(3,3)
      DOUBLE PRECISION tstress
 
      hydro_stress = 0.0
      avg_stress(1:3,1:3) = 0.0
      natoms = 0
      DO i = 1 , Numnp
         IF ( Isrelaxed(i)==1 ) THEN
            natoms = natoms + 1
            DO j = 1 , 3
               hydro_stress = hydro_stress + Virst(j,j,i)
               DO k = 1 , 3
                  avg_stress(j,k) = avg_stress(j,k) + Virst(j,k,i)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
 
      avg_stress(1:3,1:3) = avg_stress(1:3,1:3)/natoms
      hydro_stress = hydro_stress*1.0/3.0/natoms
      PRINT * , ' '
      PRINT * , &
     &      ' **** Average Stresses in Atomistic Region (eV/A^3) ****'
!      print*, avg_stress(1:3,1:3)
      PRINT * , 'hydrostatic stress = ' , hydro_stress
      tstress = 0.5*(avg_stress(1,1)+avg_stress(3,3))
      PRINT * , 's(1,1) s(2,2) s(3,3) '
      PRINT * , 's(2,3) s(1,3) s(1,2) '
      PRINT * , 'strse=' , avg_stress(1,1) , avg_stress(2,2) , &
     &      avg_stress(3,3)
      PRINT * , 'strss=' , avg_stress(2,3) , avg_stress(1,3) , &
     &      avg_stress(1,2)
      PRINT * , ' '
 
      END SUBROUTINE HSTRESS
