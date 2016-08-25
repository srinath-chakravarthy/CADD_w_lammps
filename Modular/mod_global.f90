!!$*==mod_global.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!!$*******************************************************************
!!$**
!!$     <**  MODULE mod_global: contains routines which set all relevent
!!$     global
!!$     !!$**                        flags and dimensions at the start of a
!!$     run.
!!$     !!$**
!!$     !!$**  Variable Definitions:
!!$     !!$**  ---------------------
!!$**
!!$     !!$**  NOTE: variables marked by (????) are ones that I believe
!!$     should be
!!$**   eliminated eventually (R.E. Miller (jul 02))
!!$**
!!$     <**  Integer Variables:
!!$     <**    nxdm                 -  Dimension of the nodal coord array
!!$     x()
!!$     <**    ndm                  -  Number of spatial dimensions
!!$     <**    ndf                  -  Number of degrees of freedom per
!!$     node
!!$     <**    nen                  -  Number of nodes per element
!!$     <**    nen1                 -  Number of nodes per element plus 1
!!$     <**    nsdm                 -  Number of stress components
!!$     <**    numnp                -  current number of nodes in the mesh
!!$     <**    numel                -  current number of elements in the
!!$     mesh
!!$     <**    neq                  -  number of equations (ndf*numnp)
!!$     <**    nquad         (????) -  number of quadrature points
!!$     <**    nad           (????) -  number of additional degrees of
!!$     freedom per element
!!$     <**    nshpdm        (????) -
!!$     <**    nstad         (????) -
!!$     <**    nst0          (????) -
!!$     <**    nentot        (????) -
!!$     <**    maxnp                -  SIZE of allocated storage for nodes
!!$     <**    maxel                -  SIZE of allocated storage for
!!$     elements
!!$     <**    nxsj                 -  dimension of shape function array,
!!$     xsj
!!$     <**    nshp                 -  number of shape functions,
!!$     <**    nstr                 -  dimension of stress and strain
!!$     arrays, str,eps
!!$     <**    nst           (????) -
!!$     <**   CUTFACT:
!!$     !!$**   EFFECTIVE CUTOFF FACTOR (in multiples of rcut gives the
!!$     radius in
!!$     !!$**   the deformed configuration contributing to an atom's
!!$     neighbor lis).
!!$     !!$**   The larger this factor the fewer updates to the neighbor
!!$     lists will
!!$     !!$**   be necessary during relaxation, but more atoms will be
!!$     sampled at
!!$     !!$**   each step.
!!$**
!!$     <**   New CADD variables:
!!$**
!!$     <**   nqc       -        number of atoms in the atomistic region
!!$     <**   nspring   -        number of interface atoms
!!$     <**   numnpp1   -        if its -1, all reference to the brinell
!!$**   indenter is ignored.  Otherwise it should be
!!$**   set to numnp+1, the location where the brinell indenter is stored
!!$**   in
!!$**   the atom list
!!$     <**   z_length  - periodic length in the z direction
!!$     <**   DD_set    - flag to show that leo's stuff is initialized
!!$     <**   newmesh   - flag to tell the detection routine to re
!!$     -initialize
!!$     <**   b0        - stores the displacements at the start of the
!!$     time
!!$     !!$**               step, before any relaxation starts.
!!$     <**   INDRAD    - brinell indenter radius
!!$     <**   INDRADSQ  - INDRAD**2
!!$     <**   energy    - site energy of each atom
!!$     <**   IsRelaxed - identifies each atom type:  -1 pad
!!$     !!$**                                            0 continuum
!!$     !!$**                                            1 atom
!!$     !!$**                                            2 interface
!!$     <**   idtemp    - used to fix (set force to zero) constrained
!!$     degrees of freedom.
!!$**
!!$**
!!$**
!!$**
!!$**
!!$**   Contains Routines:
!!$**   ------------------
!!$**
!!$**   GlobalSettings -- Reads in global feap settings and sets others to
!!$**   default values
!!$**   EchoSettings   -- Echoes all settings.
!!$**
!!$**************************************************************
 
      MODULE MOD_GLOBAL
      IMPLICIT NONE
!!$*--MOD_GLOBAL97
 
!!$     * Variable Definintions
      INTEGER nxdm , ndm , ndf , nen , nen1 , nsdm , numnp , numel , &
     &        neq , nquad , nad , nshpdm , nstad , nst0 , nentot , &
     &        maxnp , maxel , nxsj , nshp , nstr , nst , niter , nqc , &
     &        nspring
      INTEGER :: numnpp1 = -1
      DOUBLE PRECISION z_length , convergetol , rnmax
      LOGICAL dd_set , newmesh
      DOUBLE PRECISION , POINTER :: b0(:,:)
      DOUBLE PRECISION :: cutfact , indrad = -1.D0 , indradsq
      DOUBLE PRECISION , ALLOCATABLE :: energy(:)
!!$     double precision, allocatable:: amass(:)
      INTEGER , ALLOCATABLE :: isrelaxed(:)
      LOGICAL , ALLOCATABLE :: idtemp(:,:)
      INTEGER i_initial , i_final , indextimeh
      INTEGER :: numtoth = 0
!!$     Jun Song Modification
!!$     NHrescale means how many steps to do H stablizer
!!$     MaxTemp, do H stablizer if exceed this Temp
!!$     NVEFlag determines whether to use H in thermostat
!!$     Number of MD steps doing temp rescaling
      INTEGER nummdrescale , debugflag
      LOGICAL userescale
      INTEGER nhrescale , hnveflag
      INTEGER numperiodz , num2dnode
      DOUBLE PRECISION twindow , maxhtemp
!!$     *Qu modification begins related to virial stress
      INTEGER , ALLOCATABLE :: imaterial(:)
                                           !!$!!$, boundary(:)
      DOUBLE PRECISION , ALLOCATABLE :: virst(:,:,:) , rssgb(:)
      DOUBLE PRECISION , ALLOCATABLE :: avevirst(:,:,:)
      DOUBLE PRECISION , ALLOCATABLE :: rsatomstress(:,:,:)
      DOUBLE PRECISION , ALLOCATABLE :: rsatomstressf(:,:,:)
      INTEGER , ALLOCATABLE :: atomspecie(:)
!!$     Jun Song: neighbor flag-determine if updating an atom's neighbor
!!$     DampForce store the damp forces for Marder&langevin, allocate in
!!$     pmain.f
      INTEGER , ALLOCATABLE :: updateneigh(:)
      DOUBLE PRECISION , ALLOCATABLE :: dampforce(:,:)
!!$     Jun Song: total_energyMD is total energy for system
!!$     energyH is energy for H atom
!!$     timestep1 is to implement 2 timesteps
!!$     SimStep is current # of simlation steps
      DOUBLE PRECISION total_energymd , energyh , timestep1
      INTEGER dim , simstep
      DOUBLE PRECISION perthick
!!$     Global parameter for Tempset-set in dosteps. Default 0.01
!!$     Global parameter for NoseHoover, from input
      DOUBLE PRECISION :: systemp = 0.01D0 , nhdampcoeff
!!$     JS-Scale Rationfor Langevin (e.g., Use the coefficient in Marder
!!$     as base
!!$     If 0.5, use half the coefficient in Marder)
      DOUBLE PRECISION lvscaleratio
!!$     *Qu modification ends
      CHARACTER*80 head
      DOUBLE PRECISION timeol , time , dt , boltzmannconst
      INTEGER indexpad , indexinterface , indexcontinuum , indexatom
      INTEGER maxneighbors , ihnumber
      LOGICAL initplot
!!$     ------Hack -SC
!!$     This is set in mesh.f to lattice spacing in x.
!!$      double precision, parameter :: x_move_mesh = 20.0
      DOUBLE PRECISION :: x_move_mesh
      LOGICAL :: movemesh , moved
!!$     Current Crack tip position and other crack tip related parameters
      DOUBLE PRECISION :: xtip(2) , xtip_old(2) , xtip_init(2)
      DOUBLE PRECISION :: xtip_actual(2) , x_tip_dir(2) , crack_motion(2)
!!$     -----End Hack -SC
!!$     Since the xtip is located in the atomistic region and is always an
!!$     atom
!!$     position in the reference coordinate it will always be a multiple
!!$     of the
!!$     lattice spacing in x
      DOUBLE PRECISION :: pad_width
!!$     Adjustable parameter read from input file to control the width of
!!$     the pad atom region. Needs to be bigger to accomodate more
!!$     dislocations on the same slip plane entering the atomistic region.
!!$     Dislocations can be in the pad region provided the template method
!!$     is used.
 
!!$     As the pad region gets bigger, there is a jump in the dislocation
!!$     passing algorithm that makes the dislocation move past the pad
!!$     atoms going either way Atom-> continuum or Continuum-> atom.
 
!!$     Template method not used here .... need implementation
!!$ ******************************************************************************************
!!$      MD Parameters
      !!$ ******************************************************************************************
      integer :: exclude_top, exclude_bot, exclude_right, exclude_left !!$> @var parameters to control particle identification

!!$     Stadium Langevin parameters
      double precision :: stadium_xmin, stadium_xmax, stadium_ymin, stadium_ymax, stadium_width
!!$     Other MD parameters, currently local to md.f90 
      double precision :: damp_coeff           !!$> @var Langevin damping coeff inverse time units --> gets converted to time internally
      double precision :: lammps_temperature   !!$> @var Desired temperature 
      double precision :: damp_ratio           !!$> @var Damping ratio ---> legacy parameter in case we need to scale ratio
      double precision :: lammps_timestep      !!$> @var lammps_timestep --> if not set = 1.e-3 (metal units in lammps)
      integer :: fem_update_steps              !!$> @var number of md steps before fem update  (10 - 50) typically
      integer :: lammps_output_steps           !!$> @var number of md steps to output dump files (1 or fem_update_steps or num_md_steps)
      integer :: num_md_steps                  !!$> @var number of md steps to perform before increasing fem load
      integer :: num_restart_steps             !!$> @var number of md steps before writing restart files
      integer :: num_initial_equil             !!$> @var number of initial steps for temperature equilibriation
      double precision :: tstart, tstop        !!$> @var start and stop of temperature for thermostat (currently equal)
      
      !!$     Particle Parameters
      double precision :: particle_velocity !!$> @var particle impact velocity (Angstrom/picosecond)
      double precision :: particle_radius   !!$> @var particle radius (Angstrom), <= (xmax - xmin)/2.0 and (ymax-ymin)/2.0
      double precision :: particle_height   !!$> @var particle initial height above top surface (Angstrom)
      double precision :: particle_rotation !!$> @var particle rotation for different impact orientations
      double precision :: impact_angle      !!$> @var impact angle of particle with substrate



      CONTAINS
 
!!$------------------------------------------------------------------
!!$ GlobalSettings -- Reads in global feap settings and sets others to
!!$                   default values
!!$
!!$      Author :
!!$            R. Miller (01/13/98)
!!$
!!$      Revisions :
!!$            R. Miller (10/10/98)
!!$            S. Chakravarthy (2012-2014)
!!$---
      SUBROUTINE GLOBALSETTINGS
        USE MOD_FILE
        IMPLICIT NONE
!!$*--GLOBALSETTINGS203
!!$
!!$**   Store FEAP values
!!$
      READ (input_file_unit,*) MAXnp , MAXel , NDM , NDF , NEN , NSDm , NQUad , NAD
      READ (input_file_unit,*) NUMperiodz
      NXDm = 3
      NEN1 = NEN + 1
!
      NSDm = MAX0(NSDm,1)
      NQUad = MAX0(NQUad,1)
      NAD = MAX0(NAD,0)
!
      NEN1 = NEN + 1
      NST0 = NEN*NDF
      NSTad = NAD*NDF
      NST = NST0 + NSTad
      NENtot = NEN + NAD
      NSTr = NQUad*NSDm
      NSHpdm = (NDM+1)*(NEN+NAD)
      NSHp = NQUad*NSHpdm
      NXSj = NQUad
 
      INDexcontinuum = 0
      INDexpad = -1
      INDexinterface = 2
      INDexatom = 1
      INItplot = .FALSE.
 
!**   default
      CUTfact = 1.4D0
      MAXneighbors = 100
 
!**   physical constants
      BOLtzmannconst = 8.629064E-05 ! eV/K
      END SUBROUTINE GLOBALSETTINGS
!**   ---------------------------------------------------------------
!**   EchoSetting  :  Echo all global settings of note to std out.
!**
      SUBROUTINE ECHOSETTINGS
 
      IMPLICIT NONE
!*--ECHOSETTINGS245
!**   echo
      WRITE (*,*)
      WRITE (*,*) '*** Global Parameter Settings ***'
      WRITE (*,*)
      WRITE (*,99002) 'Number of spatial dimensions: ' , NDM
      WRITE (*,99002) 'Number of degrees of freedom per node: ' , NDF
      WRITE (*,99002) 'Number of nodes per element: ' , NEN
      WRITE (*,99002) 'Number of stress components: ' , NSDm
      WRITE (*,*)
      WRITE (*,99001) 'Effective cutoff factor: ' , CUTfact
99001 FORMAT (5x,a50,g14.5)
      WRITE (*,*)
99002 FORMAT (5x,a50,i7)
      END SUBROUTINE ECHOSETTINGS
 
      END MODULE MOD_GLOBAL
!*==mod_parallel.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
      MODULE MOD_PARALLEL
      IMPLICIT NONE
!*--MOD_PARALLEL265
 
!     use '/opt/hpmpi/use/mpif.h'
      INTEGER nprocs , rank , ierr
 
      END MODULE MOD_PARALLEL
!*==mod_timming.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
 
 
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!module with timming variables!!!!!!!!!!!!!!!!!
      MODULE MOD_TIMMING
      IMPLICIT NONE
!*--MOD_TIMMING280
!
!     ct1 - time at call to ma06
!     ct2 - local start time
!     ct3 - local end time
!     ct4 - MD time
!     ct5 - detection band time
!     ct6 - time in ma06
!     ct7 - time in fem
 
      REAL ct1 , ct2 , ct3 , ct4 , ct5 , ct6 , ct7 , ct8
 
      END MODULE MOD_TIMMING
