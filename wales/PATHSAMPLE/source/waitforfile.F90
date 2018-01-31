!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2015 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

SUBROUTINE WAIT_FOR_FILE(FILENAME,SLEEPTIME,NTRIES)

! nagfor needs a module to recognise SLEEP
#ifdef NAGFOR
  USE F90_UNIX_PROC
#endif

! subroutine to check if the FILENAME file exists NTRIES times, SLEEPing SLEEPTIME seconds between checks
! if the file does not exist after all checks are done, STOP the code with an error.
IMPLICIT NONE
INTEGER, INTENT(IN)          :: NTRIES, SLEEPTIME
INTEGER                      :: TRY
CHARACTER(LEN=*), INTENT(IN) :: FILENAME
LOGICAL                      :: THERE

! set the existence flag to .FALSE.
THERE=.FALSE.

! check to see if the file has been created NTRIES times 
DO TRY=1, NTRIES
   INQUIRE(FILE=TRIM(ADJUSTL(FILENAME)), EXIST=THERE)
   IF(THERE) THEN
      ! if the file exists, exit the loop
      ! DEBUG PRINTING
      !PRINT *,"wait_for_file> ",TRIM(ADJUSTL(FILENAME))," present, moving on"
      EXIT
   ELSE
      ! if the file doesn't exist yet, SLEEP and check again
      ! DEBUG PRINTING
      !PRINT *,"wait_for_file> ",TRIM(ADJUSTL(FILENAME))," not present, SLEEPing for ",SLEEPTIME
      CALL SLEEP(SLEEPTIME)
   ENDIF
ENDDO

! if the file is still not present, STOP the program with a helpful error
IF (.NOT.THERE) THEN 
   PRINT *,"wait_for_file> ERROR: file ",TRIM(ADJUSTL(FILENAME))," not created correctly"
   STOP
ENDIF

RETURN
END SUBROUTINE WAIT_FOR_FILE
