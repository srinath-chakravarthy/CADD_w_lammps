!*==mod_file.spg  processed by SPAG 6.70Rc at 12:39 on 29 Oct 2015
!***********************************************************************
!**
!**   MODULE mod_file : contains routines relevant to file handling
!**
!**
!**   Variable Definitions:
!**   ---------------------
!**
!**   none
!**
!**
!**   Contains Routines:
!**   ------------------
!**
!**   FileExists   - check for the existence of a file
!**   GetFreeUnit  - get a unit number not being currently used
!**   iofile       - opens files and returns a unit number (uses getfree
!**
!***********************************************************************
 
MODULE MOD_FILE
  IMPLICIT NONE
    INTEGER :: input_file_unit

!*--MOD_FILE24
 
      CONTAINS
 
!------------------------------------------------------------------
! FileExists -- returns .true. if file 'filename' exists in the current
!               directory
!
!      Passed Parameters :
!          filename (in) : filename
!             print (in) : logical variable, set to .true. if function
!                          should print a message if file is not found
!
!      Module Parameters :
!                      none
!
!      Algorithm :
!            self explanatory
!
!      Notes :
!          none
!
!      Author :
!            E.B.Tadmor (12/31/97)
!
!      Revisions :
!              none
!
!--
      LOGICAL FUNCTION FILEEXISTS(Filename,Print)
 
      IMPLICIT NONE
!*--FILEEXISTS56
 
!** Transferred Variables **|
      CHARACTER(LEN=*) , INTENT(IN) :: Filename
      LOGICAL , INTENT(IN) :: Print
 
      INQUIRE (FILE=Filename,EXIST=FILEEXISTS)
      IF ( .NOT.FILEEXISTS .AND. Print ) PRINT * , '***ERROR: File ' , &
     &     TRIM(Filename) , ' not found.'
      END FUNCTION FILEEXISTS
 
 
!------------------------------------------------------------------
! GetFreeUnit -- return the number of a unit number currently not
!                being used.
!
!      Passed Parameters :
!                      none
!
!      Module Parameters :
!                      none
!
!      Algorithm :
!            Loop over all possible unit numbers until an unallocated
!            one is found, otherwise report an error and stop.
!
!      Notes :
!          none
!
!      Author :
!            E.B.Tadmor (12/31/97)
!
!      Revisions :
!              none
!
!--
      INTEGER FUNCTION GETFREEUNIT()
 
      IMPLICIT NONE
!*--GETFREEUNIT95
 
!** Local Variables **|
      LOGICAL inuse
!
!     reserved file numbers: 5 = std in
!     6 = std out
!     99 = abort.dat
!
      DO GETFREEUNIT = 7 , 98
         INQUIRE (UNIT=GETFREEUNIT,OPENED=inuse)
         IF ( .NOT.inuse ) RETURN
      ENDDO
      PRINT * , '***ERROR: Could not obtain a free unit handle.'
      STOP
      END FUNCTION GETFREEUNIT
!------------------------------------------------------------------
! iofile -- if filename is open, returns the unit number, otherwise
!           gets a free unit number and opens filename.
!
!      Passed Parameters :
!          filename (in) : filename
!            format (in) : character string for file format
!            logic (out) : unit number
!
!      Module Parameters :
!                      uses GetFreeUnit
!
!      Algorithm :
!            self explanatory
!
!      Notes :
!          none
!
!      Author :
!            R.E.Miller and Feap (09/25/95)
!
!      Revisions :
!              none
!
!--
      SUBROUTINE IOFILE(Filename,Format,Logic,Verbose)
      IMPLICIT NONE
!*--IOFILE138
!
!---- io file management: open file and return unit number
!
      CHARACTER*80 Filename
      CHARACTER*11 Format
      INTEGER Logic
      LOGICAL ofl , Verbose
!     :
      INQUIRE (FILE=Filename,OPENED=ofl)
      IF ( .NOT.ofl ) THEN
         Logic = GETFREEUNIT()
         OPEN (Logic,FILE=Filename,STATUS='unknown',FORM=Format)
         IF ( Verbose ) THEN
            WRITE (*,*) '  ** Opening ' , Format , ' file:'
            WRITE (*,*) '        ' , Filename
         ENDIF
      ELSE
         INQUIRE (FILE=Filename,NUMBER=Logic)
         IF ( Verbose ) THEN
            WRITE (*,*) '  ** Accessing previously opened file:'
            WRITE (*,*) '        ' , Filename
         ENDIF
      ENDIF
      END SUBROUTINE IOFILE
 
      END MODULE MOD_FILE
