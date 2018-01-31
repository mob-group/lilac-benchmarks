!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
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

! SAT  Wed Dec  8 13:55:12 GMT 2004
module utils
     implicit none
     contains

     subroutine openiostat(i,fname,subrid)
          implicit none
          
          integer,intent(in) :: i
          character(len=*),intent(in) :: fname
          character(len=*),intent(in),optional :: subrid

          if (i/=0) then
               print *, 'ERROR while opening file "'//trim(adjustl(fname))//'"'
               if (present(subrid)) print *, 'STOP in '//trim(adjustl(subrid))
               stop
          endif
end subroutine openiostat

INTEGER FUNCTION GETUNIT()
IMPLICIT NONE
LOGICAL :: INUSE
!
! start checking for available units > 10, to avoid system default units
!
INTEGER :: UNITNUM

INUSE=.TRUE.
UNITNUM=11

DO WHILE (INUSE)
   INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
   IF (.NOT.INUSE) THEN
      GETUNIT=UNITNUM 
   ELSE     
      UNITNUM=UNITNUM+1 
   ENDIF
ENDDO
END FUNCTION GETUNIT

end module utils
