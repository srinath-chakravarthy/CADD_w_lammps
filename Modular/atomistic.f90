!*==processclump.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!**---------------------------------------------------------------------
!** atomistic.f : atomistic routines
!--
 
!**---------------------------------------------------------------------
!** ProcessClump : process clump
!**
!**      Non-Obvious Parameters :
!**          dr (out) : out-of-balance forces
!**            fls (in)  : true = compute forces
!**            flw (in)  : true = compute energy
!**      Algorithm :
!**              Loop through each atom
!**                   And add its contribution
!**              end do
!--
 
      SUBROUTINE PROCESSCLUMP(Id,X,Ix,F,B,Dr,Fls,Flw)
      USE MOD_GLOBAL
      USE MOD_POTEN
      USE MOD_DYNAMO
      USE MOD_FILE
      IMPLICIT NONE
!*--PROCESSCLUMP25
 
!--Variables Transferred
 
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*) , F(NDF,*) , Dr(NDF,*)
      DOUBLE PRECISION total_energy , xdef , ydef , zdef , xini , yini ,zini
      DOUBLE PRECISION bb1 , bb2 , bb3 , aa1 , aa2 , aa3
      INTEGER Id(NDF,*) , Ix(NEN1,*) , no_mdatoms
 
      LOGICAL Fls , Flw , needlist
 
 
!--Local Variables
      INTEGER irep , iel , i , j , inode , idf
      INTEGER logic
      CHARACTER*80 filename
      character*80 str
      LOGICAL STRessflag
      COMMON /FLAG  / STRessflag
!
      IF ( .NOT.STRessflag ) THEN
         ALLOCATE (VIRst(3,3,NUMnp))
         STRessflag = .TRUE.
      ENDIF
      VIRst(1:3,1:3,1:NUMnp) = 0.D0
 
 
!     chkdis checks whether neighbor lists need updating
!
      CALL CHKDIS(B,X)
      total_energy = 0.D0
      no_mdatoms = 0
      TOTal_energymd = 0.D0
      ENErgyh = 0.D0
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
!C--Jun Song Screen Print if updating neighborlist
      IF ( NEWlst==1 ) THEN
         WRITE (*,*) "ccccccccccccccccccccccccccccccc"
         WRITE (*,*) "Updating NeighborList at " , SIMstep
      ENDIF
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!      print *, 'Total no. of nodes =', numnp
      DO irep = 1 , NUMnp
         NeedList=IsRelaxed(irep).ne.0
!         needlist = ISRelaxed(irep)>=1
         CALL GNEIGH(irep,B,X,needlist)
         IF ( needlist ) THEN
	    if (irep == 4082) then 
	      write(*,*) 'Number of neighbors of ', 4082, ' = ', nneips
	      write(*,'(I6,1X,30I7)') 4082, (jneigh(i), i = 1,nneips) 
	      write(*,'(I6,1X,30F15.9)') 4082, (Dneigh(3,i), i = 1,nneips) 
	      write(*,'(I6,1X,30F15.9)') 4082, (rneigh(i), i = 1,nneips) 
	    end if
!     NonLocal adds the energy and forces due to the E_i term (energy of
!     atom i) in the total energy functional.
!
 
            CALL NONLOCAL(Id,X,Ix,F,B,Dr,Fls,Flw,irep)
           if (irep == 4082) then 
	    write(*, '(A25,I5,1X,2(F15.8))') 'Force on atom in nei 1 ', 4082, dr(1:2, 4082)
	   end if
            IF ( ISRelaxed(irep)==1 ) THEN
               TOTal_energymd = TOTal_energymd + ENErgy(irep)
 
!C--MS: Get H displacements and E
               IF ( ATOmspecie(irep)==2 ) THEN
                  bb1 = B(1,irep)
                  bb2 = B(2,irep)
                  bb3 = B(3,irep)
                  aa1 = X(1,irep)
                  aa2 = X(2,irep)
                  aa3 = X(3,irep)
                  ENErgyh = ENErgyh + ENErgy(irep)
               ENDIF
               no_mdatoms = no_mdatoms + 1
            ENDIF
         ENDIF
      ENDDO
!       print*, h_no
!          xdef =(x(1,h_no)+ b(1,h_no))/95.5171
!           ydef =(x(2,h_no)+ b(2,h_no))/98.2759
!           zdef = (x(3,h_no) + b(3,h_no))/4.999
!          xini =(x(1,h_no))/95.5171
!           yini =(x(2,h_no))/98.2759
!           zini = (x(3,h_no))/4.999
!       print*,'poz_H',xdef,  ydef, zdef
!       print*,'poz_H1',xini, xini,xini
!       print*,'poz_b',bb1, bb2,bb3
!       print*,'poz_a',aa1, aa2,aa3
!      print*,'total_energy', total_energy
!      print*,'total_energy',total_energy
!
!	print*,'total_energyMD',total_energyMD
!	print*,'no of MD atoms', no_MDatoms
!	stop
 
!      print*, 'energy of H atom', energyH
!
!     anything involving "numnpp1" is hard-wired to the brinell indenter
!     and is written so that it will be "turned off" if numnpp1=-1,
!     which is the initial value.
!
      write(*, '(A25,I5,1X,2(F15.8))') 'Force on atom in nei 2', 4082, dr(1:2, 4082)

      IF ( NUMnpp1<NUMnp ) RETURN
      CALL GNEIGH(NUMnpp1,B,X,.TRUE.)
      CALL NONLOCAL(Id,X,Ix,F,B,Dr,Fls,Flw,NUMnpp1)
      

      
      
      END SUBROUTINE PROCESSCLUMP
!*==processclumpwrapper.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct
 
 
      SUBROUTINE PROCESSCLUMPWRAPPER(Id,X,Ix,F,B,Dr,Eatom)
      USE MOD_GLOBAL
      USE MOD_POTEN
      USE MOD_DYNAMO
      USE MOD_FILE
      IMPLICIT NONE
!*--PROCESSCLUMPWRAPPER136
 
      DOUBLE PRECISION B(NDF,*) , X(NXDm,*) , F(NDF,*) , Dr(NDF,*) , &
     &                 Eatom(*)
      INTEGER Id(NDF,*) , Ix(NEN1,*) , iatom
      LOGICAL fls , flw
 
      fls = .TRUE.
      flw = .TRUE.
 
      CALL PROCESSCLUMP(Id,X,Ix,F,B,Dr,fls,flw)
      DO iatom = 1 , NUMnp
         Eatom(iatom) = ENErgy(iatom)
      ENDDO
 
 
      END SUBROUTINE PROCESSCLUMPWRAPPER
 
 
