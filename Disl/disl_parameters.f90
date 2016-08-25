!*==mod_disl_parameters.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct
! $Id: mod_disl_parameters,v 1.1.1.12003/03/1220:09:00 shastry Exp $
!
      MODULE MOD_DISL_PARAMETERS
      IMPLICIT NONE
!*--MOD_DISL_PARAMETERS6
 
      LOGICAL,PARAMETER ::  NUMERICALPK = .TRUE.
      DOUBLE PRECISION, PARAMETER ::  PEIERLS = 0.0D0
 
      INTEGER, PARAMETER :: MAX_DISL = 5000

!!$
!!$>ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$>
!!$> Common blocks
!!$> /arrays_disl/:
!!$>        burgers(3,max_disl) Burgers vectors
!!$>        burg_length(max_disl) length of the inplane component of
!!$>                               Burgers vectors
!!$>        theta_e(max_disl) direction of the cut for the edge cmponent
!!$>                           (w.r.t the Burgers vector ( [-P, Pi) )
!!$>        theta_s(max_disl) direction of the cut fot the screw component
!!$>                           ( w.r.t the X-axis, any walue )
!!$>        r_disl(3, max_disl) current positions
!!$>        r_old(3, max_disl) positions at last check for lost disl
!!$>        pk_stress(3, max_disl) stress from acting on a dislocation
!!$>        pk_force(2, max_disl) Peach-Koeller force acting on a
!!$>                               dislocation
!!$>        pk_f(max_disl)  Peach-Koeller force in the direction of the
!!$>                            Burgers vector
!!$>        elem_disl(max_disl) the number of an element containing a
!!$>                            dislocation
!!$>
!!$> /contrl_disl/:
!!$>        ndisl  number of dislocations
!!$>        i_disl flag
!!$>
!!$>cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$23456789012345678901234567890123456789012345678901234567890123456789012
!!$         1         2         3         4         5         6         7
!!$
      INTEGER :: ndisl , i_disl , ndisl_dd(MAX_DISL)
      INTEGER :: elem_disl(MAX_DISL) , disl_index(MAX_DISL)
      DOUBLE PRECISION :: burgers(3,MAX_DISL) , burg_length(MAX_DISL)
      DOUBLE PRECISION :: theta_e(MAX_DISL) , theta_s(MAX_DISL)
      DOUBLE PRECISION :: r_disl(3,MAX_DISL)
      DOUBLE PRECISION :: pk_stress(3,MAX_DISL)
      DOUBLE PRECISION :: pk_force(2,MAX_DISL) , pk_f(MAX_DISL)
      DOUBLE PRECISION :: disl_range(2,MAX_DISL) , r_old(3,MAX_DISL)
      DOUBLE PRECISION :: disl_residence(2,2,MAX_disl) !> @var containing both x and y coordinates for range of dislocation
      DOUBLE PRECISION :: disl_line(3,MAX_DISL), disl_image(MAX_disl)
      INTEGER :: disl_timer(MAX_DISL)
      LOGICAL :: irmdisl(MAX_DISL)
      
      LOGICAL MOVedisl1
      COMMON /KMVDIS/ MOVedisl1
      INTEGER, PARAMETER :: time_static = 20
 
      END MODULE MOD_DISL_PARAMETERS
 
 
!
! $Log: mod_disl_parameters,v $
! Revision 1.1.1.1  2003/03/1220:09:00  shastry
! vijay-   Initial import.
!
! Revision 1.4  2002/03/0503:00:45  shilkrot
! Moved pk_b out of dislocations.par into the cg common block.
!
! Revision 1.3  2002/02/2520:52:26  shilkrot
! Added an array burg_length holding the length of the in-plane componen
! Burgers vectors.
!
! Revision 1.2  2001/11/0623:34:16  shilkrot
! Added variables to store the P. - K. force and to do c. g.
!
! Revision 1.1  2001/07/1206:30:20  shilkrot
! The routines used to apply the tilde field and to handle the array of
! dislocations.
!
