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

MODULE DOCKMODULE
USE PORFUNCS
implicit none
character(len=30), allocatable  ::  PATHSTRINGARRAY(:)
double precision, allocatable :: COORDSLIGAND(:,:), COORDSCOMPLEX(:,:), COORDSPROTEIN(:,:)
double precision, allocatable :: ENERGIESLIGAND(:), ENERGIESCOMPLEX(:), ENERGIESPROTEIN(:)
double precision, allocatable :: ENERGIESLIGANDSAVE(:), ENERGIESCOMPLEXSAVE(:), ENERGIESPROTEINSAVE(:)
double precision, allocatable :: FELIGAND(:), FECOMPLEX(:), FEPROTEIN(:)
double precision, allocatable :: RELFELIGAND(:), RELFECOMPLEX(:), RELFEPROTEIN(:)
double precision, allocatable :: WLIGAND(:), WCOMPLEX(:), WPROTEIN(:) ! arrays storing Boltzmann weights
double precision, allocatable :: FVIBMINLIGAND(:), FVIBMINCOMPLEX(:), FVIBMINPROTEIN(:), &
                              & LPLIGAND(:),LPCOMPLEX(:),LPPROTEIN(:)
character(len=4), allocatable :: ZSYMLIGAND(:),ZSYMCOMPLEX(:),ZSYMPROTEIN(:)
integer :: natomsligand, natomscomplex, natomsprotein,totalmin(3)
integer, allocatable :: indexarrayl(:),indexarrayc(:),indexarrayp(:)
integer, allocatable :: indexarraylfe(:),indexarraycfe(:),indexarraypfe(:)
double precision :: gminligand,gmincomplex,gminprotein,gminbinding,meane(3),totalw(3),rmat(3,3)
contains

SUBROUTINE DOCK
USE COMMONS
USE PORFUNCS
implicit none
integer :: istat

if(dstage(1)) then
   write(*,*) 'dock> STAGE 1 -- calling preparegminfiles'
   CALL PREPAREGMINFILES
end if
if(dstage(2)) then
   write(*,*) 'dock> STAGE 2 -- calling dockcycle'
   CALL DOCKCYCLE
   CALL FLUSH(6)
end if
if(dstage(3)) then
   write(*,*) 'dock> STAGE 3 -- calling readminima'
   CALL READMINIMA
end if
if(dstage(4)) then
   write(*,*) 'dock> STAGE 4 -- calling prepareoptimfiles'
   CALL PREPAREOPTIMFILES 
end if
if(dstage(5)) then
   write(*,*) 'dock> STAGE 5 -- calling dockcycle2'
   CALL DOCKCYCLE2
end if
if(dstage(6)) then
   write(*,*) 'dock> STAGE 6 -- calling analyseresults'
   CALL ANALYSERESULTS
end if





! analytical frequencies from NAB (with AMBER) will be rubbish if cutoff is used!!!!
! we'll need to see how much using no cutoffs would slow down the runs.
! Otherwise numerical Hessian calculation should also work, the question is
! how much time one evaluation takes.
WRITE(*,*) 'master process exiting here'
CALL FLUSH(6)

CALL EXIT(ISTAT)
END SUBROUTINE DOCK

SUBROUTINE DOCKCYCLE
! submit GMIN jobs when CPUs become available
use COMMONS
USE PORFUNCS
implicit none

integer :: totalgminjobs,J1,J2,J3,PID(NCPU),ISTAT,PIDDONE,J1SAVE
logical :: killed,parent,child
integer :: temp,ncycles,newjob,status
character(len=8) :: J1CHAR

totalgminjobs = 3*parallel ! total number of GMIN jobs
 
WRITE(*,*) 'totalgminjobs in dockcycle =', totalgminjobs, NCPU
CALL FLUSH(6)
PID(:)=0
STATUS=0
DO J1=1,TOTALGMINJOBS
        INQUIRE(FILE=TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest',exist=YESNO)
        J1SAVE=0
      IF(YESNO) THEN
        DO
         J1SAVE=J1SAVE+1
         WRITE(J1CHAR,'(I8)') J1SAVE
         INQUIRE(FILE=TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest.'//TRIM(ADJUSTL(J1CHAR)),exist=YESNO)
         IF(.NOT.YESNO) EXIT
        END DO
        CALL MYSYSTEM(ISTAT,DEBUG,'mv '//TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest '//&
        &     TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest.'//TRIM(ADJUSTL(J1CHAR)))
      END IF
END DO

J1=0
  DO J3=1,NCPU
    J1=J1+1

   WRITE(*,*) 'J1 at the beginning of dockcycle=',J1

    CALL SYSTEM('sleep 1') ! to prevent us running out of source ports
    CALL FLUSH(6) ! the child process may duplicate output without this line!
    CALL FORK_SUBR(PID(J3)) ! PID is zero in the child, non-zero in the parent
!    IF (DEBUG.AND.(PID(J3).NE.0)) WRITE(*,'(A,I16)') 'dockcycle> forked connect run process id=',PID(J3)
!    IF (PID(J3).NE.0) WRITE(*,'(A,I16)') 'dockcycle> forked connect run process id=',PID(J3)
!    IF (PID(J3).EQ.0) THEN
!       WRITE(*,'(A,2I16)') 'dockcycle> I am the child! PID=',PID(J3), J3
!    END IF  
   IF (PID(J3).EQ.0) CALL SUBMITGMINJOBS(J3,DEBUG,PATHSTRINGARRAY(J1))
        ! child processes exit at the end of SUBMITGMINJOBS
   CALL FLUSH(6)
  ENDDO
PIDDONE=0
STATUS=0
NCYCLES=0
cycles: DO !NCYCLES=1,TOTALGMINJOBS-NCPU        ! NCPU jobs are already submitted
   J1=J1+1
   NCYCLES=NCYCLES+1
   KILLED=.FALSE.
   WRITE(*,*) 'calling WAIT, J1=',J1
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   CALL WAIT_SUBR(PIDDONE,STATUS)
!  PRINT '(A,5I8)','cycle2> PIDDONE,STATUS,PID=',PIDDONE,STATUS,PID(1:NCPU)
11 CONTINUE
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   IF (PIDDONE.GT.0) THEN
      IF (DEBUG) PRINT '(A,I8,A,I6)','dockcycle> PID ',PIDDONE,' has finished with exit status ',STATUS
      DO J2=1,NCPU
         IF (PIDDONE.EQ.PID(J2)) THEN
            IF (STATUS.NE.0) KILLED=.TRUE. ! INCOMPLETE OPTIM JOBS WOULD IDEALLY RETURN A NON-ZERO EXIT CODE 
            NEWJOB=J2
            GOTO 10
         ENDIF
      ENDDO
      PRINT*,'ERROR - PID of completed child process not recognised: ',PIDDONE
      CALL EXIT(ISTAT) 
   ELSE
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      PRINT '(A,I8)','dockcycle> ERROR - WAIT returned system error code ',-PIDDONE
!
! Try calling wait again to see if this fixes things.
!
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I8))','dockcycle> on calling wait again pid=',PIDDONE,' status=',STATUS
      IF (PIDDONE.GT.0) GOTO 11
      CALL EXIT(ISTAT)
   ENDIF
10 WRITE(*,'(3(A,I8))') 'dockcycle> forked GMIN run ',NCYCLES,' on CPU ',NEWJOB,' process id ',PID(NEWJOB)

   CALL FLUSH(6) ! the child process may duplicate output without this line!
   WRITE(*,'(3(A,I8))') 'dockcycle> analysing result of GMIN job ',NCYCLES,' on CPU ',NEWJOB,' for process id ',PID(NEWJOB)
!   WRITE(CONNSTR,'(I10)') PID(NEWJOB)
   IF (KILLED) WRITE(*,'(3(A,I8))') 'dockcycle> GMIN run ',NCYCLES,' on CPU ',NEWJOB,' was unsuccessful'

   CALL FLUSH(6) ! the child process may duplicate output without this line!
   CALL FORK_SUBR(PID(NEWJOB))
   IF (PID(NEWJOB).EQ.0) CALL SUBMITGMINJOBS(NEWJOB,DEBUG,PATHSTRINGARRAY(J1))
   IF (MOD(NCYCLES,NCPU).EQ.0) THEN 
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
      WRITE(*,'(A,I8)') 'dockcycle> end of cycle ',NCYCLES/NCPU
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
   END IF
  
      IF(J1==totalgminjobs) then
         EXIT
      END IF
END DO cycles

CALL FLUSH(6)

RETURN

END SUBROUTINE DOCKCYCLE


SUBROUTINE READMINIMA
! we have a bunch of subdirectories for each of the three components of the system. We need
! to read in the coordinates of the minima before we submit the OPTIM jobs, also compare them, 
! so that minima that are the same are not being added to the database.
! some 2-dimensional allocatable arrays might be needed:
! COORDSLIGAND(%NMINIMA%,3*NATOMS), COORDSPROTEIN and COORDSCOMPLEX, and
! ENERGIESLIGAND(%NMINIMA%), ENERGIESPROTEIN and ENERGIESCOMPLEX.
! If we will be dumping min.data.info at the end of the OPTIM run, we'll need to be able to 
! read the coordinates and frequencies into the min.data file from there directly, so the
! above arrays are just for submitting OPTIM jobs and perhaps some on-the-fly analysis.
use COMMONS, only : PARALLEL,DEBUG,EDIFFTOL
implicit none
integer :: J1,J1SAVE,J2,J3,J4,J5,ISTAT,iostatus,nlines(3*PARALLEL),nminima(3*PARALLEL),dummyint,dummyint1
logical :: YESNO,sameminimum,badminimum
character(len=15) :: check,J1CHAR
character(len=80) :: check1, check2
! count the number of minima found. The actual number of minima that we are going to read in 
! will be smaller than this. Allocate arrays afterwards.

! ligand

nminima(:) = 0
totalmin(1) = 0

DO J1=1,3*PARALLEL
        INQUIRE(FILE=TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest',exist=YESNO)
        J1SAVE=0
      IF(YESNO) THEN
        DO
         J1SAVE=J1SAVE+1
         WRITE(J1CHAR,'(I8)') J1SAVE
         INQUIRE(FILE=TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest.'//TRIM(ADJUSTL(J1CHAR)),exist=YESNO)
         IF(.NOT.YESNO) EXIT
        END DO
        CALL MYSYSTEM(ISTAT,DEBUG,'mv '//TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest '//&
        &     TRIM(ADJUSTL(PATHSTRINGARRAY(J1)))//'/lowest.'//TRIM(ADJUSTL(J1CHAR)))
      END IF
END DO


DO J1=1,PARALLEL
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
   IF(YESNO) CALL MYSYSTEM(ISTAT,.TRUE.,'rm '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save' )
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd '//trim(adjustl(PATHSTRINGARRAY(J1)))//'; cat lowest.[1-9]* >lowest.save' )
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
     IF(YESNO) THEN
        OPEN(unit=311,file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',status='OLD')
         nlines(J1) = 0
         do
           read (unit=311,iostat=iostatus,fmt='(A15)') check
             if (iostatus<0) then
                if (debug) write(*,'(A)') 'readminima> finished reading file '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save'
                rewind(311)
                read (unit=311,iostat=iostatus,fmt=*) natomsligand
                if (natomsligand==0) then
                    write(*,'(A)') 'readminima> FATAL ERROR, natoms is zero in '//&
                                    &trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save'
                    call flush(6)
                    stop
                end if
                exit
             end if
           nlines(J1) = nlines(J1) + 1
         end do
         nminima(J1) =  nlines(J1)/(natomsligand+2)     
         totalmin(1) = totalmin(1) + nminima(J1)
         write(*,*) 'nminima:',nminima(J1),J1
         REWIND(311)
         CLOSE(311)
     ELSE
        IF(DEBUG) WRITE(*,'(A)') 'readminima> ' // trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest' // &
                  & ' could not be found, GMIN job was probably killed before completion'
     END IF
END DO

IF(DEBUG) WRITE(*,'(A,I8,A,I8)') 'readminima> allocating arrays for ligand, NATOMS = ',NATOMSLIGAND,' TOTALMIN = ',totalmin(1)
IF(ALLOCATED(COORDSLIGAND)) DEALLOCATE(COORDSLIGAND)
IF(ALLOCATED(ENERGIESLIGAND)) DEALLOCATE(ENERGIESLIGAND)
IF(ALLOCATED(ZSYMLIGAND)) DEALLOCATE(ZSYMLIGAND)
ALLOCATE(COORDSLIGAND(totalmin(1),3*natomsligand),ENERGIESLIGAND(totalmin(1)),ZSYMLIGAND(natomsligand))
J4=1
DO J1=1,PARALLEL
        write(*,*) 'in this loop'
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
   IF(YESNO) CALL MYSYSTEM(ISTAT,.TRUE.,'rm '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save' )
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd '//trim(adjustl(PATHSTRINGARRAY(J1)))//'; cat lowest.[1-9]* >lowest.save' )
    IF(YESNO) THEN
        OPEN(unit=311,file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',status='OLD')
         do J2=1,nminima(J1)
           badminimum=.false.
           read (unit=311,iostat=iostatus,fmt='(A15)') check ! first line
           read (unit=311,iostat=iostatus,fmt='(A18,I6,A1,F20.10)') check1,dummyint,check2,ENERGIESLIGAND(J4)!,&
!                & dummyint,check2!,ENERGIESLIGAND(J4),check,dummyint1 ! second line
           do J3=1,natomsligand
                read (unit=311,fmt='(A4,3F20.10)') ZSYMLIGAND(J3),COORDSLIGAND(J4,3*J3-2), &
                     & COORDSLIGAND(J4,3*J3-1),COORDSLIGAND(J4,3*J3)
!                write (*,'(A4,3F20.10)') ZSYMLIGAND(J3),COORDSLIGAND(J4,3*J3-2), &
!                     & COORDSLIGAND(J4,3*J3-1),COORDSLIGAND(J4,3*J3)
!                        WRITE(*,*) COORDSLIGAND(J4,1:3),J4,J3,3*J3-2
! for backwards compatibility, when the 'lowest' file is filled with zeros.
                if(COORDSLIGAND(J4,3*J3-2)==0.0D0.AND.COORDSLIGAND(J4,3*J3)==0.0D0) badminimum=.true.
           end do
!        WRITE(*,*) 'COORDSLIGAND:', COORDSLIGAND(1,1:3),J4
                if(badminimum) exit
!          check energy difference between read minimum and all previous ones
           sameminimum=.false.
           if(J4>1) then
             do J5=1,J4-1
               if(ABS(ENERGIESLIGAND(J4)-ENERGIESLIGAND(J5))<EDIFFTOL) then
!                write(*,*) 'same minimum', ENERGIESLIGAND(J4),ENERGIESLIGAND(J5)
                  sameminimum=.true.
               end if
             end do
           end if
           if(.not.sameminimum) J4 = J4 + 1
         end do
        CLOSE(311)
     END IF
!    write(*,*) 'in loop', J1
END DO
totalmin(1) = J4 - 1
!WRITE(*,*) 'energies'
!DO J1=1,totalmin
! IF(ENERGIESLIGAND(J1)==0.0D0) THEN
!    totalmin = J1-1    
!    exit
! END IF
! WRITE(*,*) ENERGIESLIGAND(J1)
!END DO
WRITE(*,'(A,I8)') 'readminima> total number of ligand minima found: ',totalmin(1)

! determine lowest energy minimum  
gminligand=ENERGIESLIGAND(1)
DO J1=2,totalmin(1)
       IF(ENERGIESLIGAND(J1)<gminligand.AND.ENERGIESLIGAND(J1)>-100000.0D0) gminligand=ENERGIESLIGAND(J1) 
END DO
WRITE(*,'(A,F20.10)') 'readminima> lowest energy of the ligand minima is: ', gminligand

! complex
nminima(:) = 0
totalmin(2) = 0
DO J1=PARALLEL+1,2*PARALLEL
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
   IF(YESNO) CALL MYSYSTEM(ISTAT,.TRUE.,'rm '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save' )
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd '//trim(adjustl(PATHSTRINGARRAY(J1)))//'; cat lowest.[1-9]* >lowest.save' )
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
     IF(YESNO) THEN
        OPEN(unit=311,file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',status='OLD')
         nlines(J1) = 0
         do
           read (unit=311,iostat=iostatus,fmt='(A15)') check
             if (iostatus<0) then
                if (debug) write(*,'(A)') 'readminima> finished reading file '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save'
                rewind(311)
                read (unit=311,iostat=iostatus,fmt=*) natomscomplex
                if (natomscomplex==0) then
                    write(*,'(A)') 'readminima> FATAL ERROR, natoms is zero in '//&
                                    &trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save'
                    call flush(6)
                    stop
                end if
                exit
             end if
           nlines(J1) = nlines(J1) + 1
         end do
         nminima(J1) =  nlines(J1)/(natomscomplex+2)     
         totalmin(2) = totalmin(2) + nminima(J1)
        CLOSE(311)
     ELSE
        IF(DEBUG) WRITE(*,'(A)') 'readminima> ' // trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save' // &
                  & ' could not be found, GMIN job was probably killed before completion'
     END IF
END DO

IF(DEBUG) WRITE(*,'(A,I8,A,I8)') 'readminima> allocating arrays for complex, NATOMS = ',NATOMSCOMPLEX,' TOTALMIN = ',totalmin(2)
IF(ALLOCATED(COORDSCOMPLEX)) DEALLOCATE(COORDSCOMPLEX)
IF(ALLOCATED(ENERGIESCOMPLEX)) DEALLOCATE(ENERGIESCOMPLEX)
IF(ALLOCATED(ZSYMCOMPLEX)) DEALLOCATE(ZSYMCOMPLEX)
ALLOCATE(COORDSCOMPLEX(totalmin(2),3*natomscomplex),ENERGIESCOMPLEX(totalmin(2)),ZSYMCOMPLEX(natomscomplex))

J4=1
DO J1=PARALLEL+1,2*PARALLEL
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
     IF(YESNO) THEN
        OPEN(unit=311,file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',status='OLD')
         do J2=1,nminima(J1)
           badminimum=.false.
           read (unit=311,iostat=iostatus,fmt='(A15)') check ! first line
           read (unit=311,iostat=iostatus,fmt='(A18,I6,A1,F20.10)') check1,dummyint,check2,ENERGIESCOMPLEX(J4)
!           read (unit=311,iostat=iostatus,fmt='(A18,I7,A1,F20.10,A21,I8)') check1,&
!                & dummyint,check2,ENERGIESCOMPLEX(J4),check1,dummyint1 ! second line
           do J3=1,natomscomplex
                read (unit=311,fmt='(A4,3F20.10)') ZSYMCOMPLEX(J3),COORDSCOMPLEX(J4,3*J3-2), &
                     & COORDSCOMPLEX(J4,3*J3-1),COORDSCOMPLEX(J4,3*J3)
                if(COORDSCOMPLEX(J4,3*J3-1)==0.0D0.AND.COORDSCOMPLEX(J4,3*J3)==0.0D0) badminimum=.true.
           end do
! for backwards compatibility, when the 'lowest' file is filled with zeros.
                if(badminimum) exit
!          check energy difference between read minimum and all previous ones
           sameminimum=.false.
           if(J4>1) then
             do J5=1,J4-1
               if(ABS(ENERGIESCOMPLEX(J4)-ENERGIESCOMPLEX(J5))<EDIFFTOL) then
!                write(*,*) 'same minimum', ENERGIESCOMPLEX(J4),ENERGIESCOMPLEX(J5)
                  sameminimum=.true.
               end if
             end do
           end if
           if(.not.sameminimum) J4 = J4 + 1
         end do
        CLOSE(311)
     END IF
!    write(*,*) 'in loop', J1
END DO
totalmin(2) = J4 - 1
WRITE(*,'(A,I8)') 'readminima> total number of complex minima found: ',totalmin(2)

! determine lowest energy minimum  
gmincomplex=ENERGIESCOMPLEX(1)
DO J1=2,totalmin(2)
       IF(ENERGIESCOMPLEX(J1)<gmincomplex.AND.ENERGIESCOMPLEX(J1)>-100000.0D0) gmincomplex=ENERGIESCOMPLEX(J1) 
END DO
WRITE(*,'(A,F20.10)') 'readminima> lowest energy of the complex minima is: ', gmincomplex


! protein
nminima(:) = 0
totalmin(3) = 0
DO J1=2*PARALLEL+1,3*PARALLEL
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
   IF(YESNO) CALL MYSYSTEM(ISTAT,.TRUE.,'rm '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save' )
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd '//trim(adjustl(PATHSTRINGARRAY(J1)))//'; cat lowest.[1-9]* >lowest.save' )
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
         nlines(J1) = 0
     IF(YESNO) THEN
        OPEN(unit=311,file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',status='OLD')
         do
           read (unit=311,iostat=iostatus,fmt='(A15)') check
             if (iostatus<0) then
                if (debug) write(*,'(A)') 'readminima> finished reading file '//trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save'
                rewind(311)
                read (unit=311,iostat=iostatus,fmt=*) natomsprotein
                if (natomsprotein==0) then
                    write(*,'(A)') 'readminima> FATAL ERROR, natoms is zero in '//&
                                    &trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save'
                    call flush(6)
                    stop
                end if
                exit
             end if
           nlines(J1) = nlines(J1) + 1
         end do
         nminima(J1) =  nlines(J1)/(natomsprotein+2)
!        write(*,*) 'natomsprotein',natomsprotein, nlines(J1), nminima(J1),totalmin 
         totalmin(3) = totalmin(3) + nminima(J1)
        CLOSE(311)
     ELSE
        IF(DEBUG) WRITE(*,'(A)') 'readminima> ' // trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save' // &
                  & ' could not be found, GMIN job was probably killed before completion'
     END IF
END DO

IF(DEBUG) WRITE(*,'(A,I8,A,I8)') 'readminima> allocating arrays for protein, NATOMS = ',NATOMSPROTEIN,' TOTALMIN = ',totalmin(3)
IF(ALLOCATED(COORDSPROTEIN)) DEALLOCATE(COORDSPROTEIN)
IF(ALLOCATED(ENERGIESPROTEIN)) DEALLOCATE(ENERGIESPROTEIN)
IF(ALLOCATED(ZSYMPROTEIN)) DEALLOCATE(ZSYMPROTEIN)
ALLOCATE(COORDSPROTEIN(totalmin(3),3*natomsprotein),ENERGIESPROTEIN(totalmin(3)),ZSYMPROTEIN(natomsprotein))

J4=1
DO J1=2*PARALLEL+1,3*PARALLEL
   INQUIRE(file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',exist=YESNO)
     IF(YESNO) THEN
        OPEN(unit=311,file=trim(adjustl(PATHSTRINGARRAY(J1)))//'/lowest.save',status='OLD')
         do J2=1,nminima(J1)
           badminimum=.false.
           read (unit=311,iostat=iostatus,fmt='(A15)') check ! first line
           read (unit=311,iostat=iostatus,fmt='(A18,I6,A1,F20.10)') check1,dummyint,check2,ENERGIESPROTEIN(J4)
!           read (unit=311,iostat=iostatus,fmt='(A18,I7,A1,F20.10,A21,I8)') check1,&
!                & dummyint,check2,ENERGIESPROTEIN(J4),check1,dummyint1 ! second line
           do J3=1,natomsprotein
                read (unit=311,fmt='(A4,3F20.10)') ZSYMPROTEIN(J3),COORDSPROTEIN(J4,3*J3-2), &
                     & COORDSPROTEIN(J4,3*J3-1),COORDSPROTEIN(J4,3*J3)
                if(COORDSPROTEIN(J4,3*J3-1)==0.0D0.AND.COORDSPROTEIN(J4,3*J3)==0.0D0) badminimum=.true.
           end do
! for backwards compatibility, when the 'lowest' file is filled with zeros.
                if(badminimum) exit
!          check energy difference between read minimum and all previous ones
           sameminimum=.false.
           if(J4>1) then
             do J5=1,J4-1
               if(ABS(ENERGIESPROTEIN(J4)-ENERGIESPROTEIN(J5))<EDIFFTOL) then
!                write(*,*) 'same minimum', ENERGIESPROTEIN(J4),ENERGIESPROTEIN(J5)
                  sameminimum=.true.
               end if
             end do
           end if
           if(.not.sameminimum) J4 = J4 + 1
         end do
        CLOSE(311)
     END IF
!    write(*,*) 'in loop', J1
END DO
totalmin(3) = J4 - 1
WRITE(*,'(A,I8)') 'readminima> total number of protein minima found: ',totalmin(3)

! determine lowest energy minimum  
gminprotein=ENERGIESPROTEIN(1)
DO J1=2,totalmin(3)
       IF(ENERGIESPROTEIN(J1)<gminprotein.AND.ENERGIESPROTEIN(J1)>-100000.0D0) gminprotein=ENERGIESPROTEIN(J1)
END DO
WRITE(*,'(A,F20.10)') 'readminima> lowest energy of the protein minima is: ', gminprotein
gminbinding = gmincomplex-gminprotein-gminligand

WRITE(*,'(A,F20.10)') 'readminima> binding energy based on lowest energies found is: ', &
        &       gminbinding 

! now we have all coordinates in their arrays, can exit this routine

RETURN

END SUBROUTINE READMINIMA

SUBROUTINE PREPAREGMINFILES
use COMMONS, only : COPYFILES, COORDSLIGANDSTR, COORDSPROTEINSTR, COORDSCOMPLEXSTR, AMBERT, CHARMMT, PARALLEL, NCPU, DEBUG
USE PORFUNCS
implicit none

logical YESNOA, YESNOB, YESNOC, YESNOD
integer ISTAT, STATUS, J1, J3
CHARACTER(LEN=5) ISTR

YESNOA=.FALSE.
YESNOB=.FALSE.
YESNOC=.FALSE.
YESNOD=.FALSE.

INQUIRE(FILE='data.ligand',EXIST=YESNOA)
INQUIRE(FILE='data.protein',EXIST=YESNOB)
INQUIRE(FILE='data.complex',EXIST=YESNOC)
INQUIRE(FILE='odata.freq',EXIST=YESNOD)

IF(.NOT.(YESNOA.AND.YESNOB.AND.YESNOC)) THEN
   WRITE(*,'(A)') ' dock> ERROR - one of the files data.complex, data.ligand or data.protein is missing'
   STOP
END IF

IF(.NOT.YESNOD) THEN
   WRITE(*,'(A)') ' dock> ERROR - odata.freq is missing! How can we run frequency calculations without it?'
   WRITE(*,'(A)') " dock> Remember, it has to contain DUMPDATA and NAB keywords for analytical frequencies with AMBER"
   STOP
END IF


IF(AMBERT) THEN
   INQUIRE(FILE='ligand/coords.ligand',EXIST=YESNOA)
   INQUIRE(FILE='protein/coords.protein',EXIST=YESNOB)
   INQUIRE(FILE='complex/coords.complex',EXIST=YESNOC)

   IF(.NOT.(YESNOA.AND.YESNOB.AND.YESNOC)) THEN
      WRITE(*,'(A)') ' dock> ERROR - one of the files complex/coords.complex, '// &
&                'ligand/coords.ligand or protein/coords.protein is missing'
      STOP
   END IF
   WRITE(COORDSLIGANDSTR,'(A)') 'coords.ligand'
   WRITE(COORDSPROTEINSTR,'(A)') 'coords.protein'
   WRITE(COORDSCOMPLEXSTR,'(A)') 'coords.complex'
ELSE IF (CHARMMT) THEN
   INQUIRE(FILE='ligand/ligand.crd',EXIST=YESNOA)
   INQUIRE(FILE='protein/protein.crd',EXIST=YESNOB)
   INQUIRE(FILE='complex/complex.crd',EXIST=YESNOC)

   IF(.NOT.(YESNOA.AND.YESNOB.AND.YESNOC)) THEN
      WRITE(*,'(A)') ' dock> ERROR - one of the files complex/complex.crd, '// &
&        'ligand/ligand.crd or protein/protein.crd is missing'
      STOP
   END IF
   WRITE(COORDSLIGANDSTR,'(A)') 'ligand.crd'
   WRITE(COORDSPROTEINSTR,'(A)') 'protein.crd'
   WRITE(COORDSCOMPLEXSTR,'(A)') 'complex.crd'
END IF

! subdirectores complex, ligand and protein should already exist before the start of the run

!CALL SYSTEM('mkdir -p complex ligand protein', ISTAT)

! PARALLEL=10

IF(.NOT.ALLOCATED(PATHSTRINGARRAY)) ALLOCATE(PATHSTRINGARRAY(3*PARALLEL))

 CALL PREPAREGMINDIRS(PARALLEL,'ligand')
 CALL PREPAREGMINDIRS(PARALLEL,'protein')
 CALL PREPAREGMINDIRS(PARALLEL,'complex')

 CALL MYSYSTEM(STATUS,.TRUE.,'cp data.ligand ligand/data' )   
 CALL MYSYSTEM(STATUS,.TRUE.,'cp data.complex complex/data' )  
 CALL MYSYSTEM(STATUS,.TRUE.,'cp data.protein protein/data' )  

DO J3=1,PARALLEL
   CALL PREPAREDATAFILE(J3,COORDSLIGANDSTR, 'ligand')
   CALL PREPAREDATAFILE(J3,COORDSPROTEINSTR, 'protein')
   CALL PREPAREDATAFILE(J3,COORDSCOMPLEXSTR, 'complex')
   CALL FLUSH(6) ! the child process may duplicate output without this line!
ENDDO

! directories are made now, we'll need an array of strings to point to these paths
  J1=1
  DO J3=1,PARALLEL
     WRITE(ISTR,'(I5)') J3
     PATHSTRINGARRAY(J1)='ligand/'//TRIM(ADJUSTL(ISTR))
     J1=J1+1
  END DO
  DO J3=1,PARALLEL
     WRITE(ISTR,'(I5)') J3
     PATHSTRINGARRAY(J1)='complex/'//TRIM(ADJUSTL(ISTR))
     J1=J1+1
  END DO
  DO J3=1,PARALLEL
     WRITE(ISTR,'(I5)') J3
     PATHSTRINGARRAY(J1)='protein/'//TRIM(ADJUSTL(ISTR))
     J1=J1+1
  END DO

!  STOP

END SUBROUTINE PREPAREGMINFILES

SUBROUTINE SUBMITGMINJOBS(ICPU,DEBUG,DIRSTRING)
! submit jobs from each subdirectory
USE PORFUNCS
USE NODES, ONLY: SSHSUBMIT, SSHPARALLEL
USE COMMONS, ONLY: DUMMYRUNT, AMBERT, NCPU, CHARMMT, UNRST,EXECGMIN
IMPLICIT NONE
INTEGER, INTENT(IN) :: ICPU
INTEGER J3, ISTAT, PID(NCPU)
CHARACTER(LEN=10) CONNSTR1, CONNSTR2
CHARACTER(LEN=80) OUTPUTSTRING
CHARACTER(LEN=*) DIRSTRING
CHARACTER(LEN=256) MYJOBSTRING
INTEGER :: CHILDPID, CONNID
LOGICAL :: DEBUG

call getpid_subr(CHILDPID)
! PRINT '(A,2I8)','in SUBMITOPTIMJOB, CHILDPID,CONNID=',CHILDPID,CONNID
WRITE(CONNSTR1,'(I10)') CHILDPID
!WRITE(*,*) 'CONNSTR1 in submitgminjobs = ', CONNSTR1, SSHPARALLEL,DEBUG
!WRITE(CONNSTR2,'(I10)') CONNID
!STOP

MYJOBSTRING=TRIM(ADJUSTL(EXECGMIN))
IF (SSHPARALLEL) then
     call SYSTEM('sleep 10')
!     WRITE(*,*) 'submitting job', icpu,istat,trim(adjustl(myjobstring)),trim(adjustl(CONNSTR1)),trim(adjustl(DIRSTRING))
     call SSHsubmitgmin(icpu,istat,trim(adjustl(myjobstring)),trim(adjustl(CONNSTR1)),trim(adjustl(DIRSTRING)))
ELSE
     CALL MYSYSTEM(ISTAT,DEBUG,'cd ' // trim(adjustl(DIRSTRING)) // ';' //trim(adjustl(myjobstring)))
ENDIF
IF (ISTAT.NE.0) PRINT '(A,I8)','submitgminjobs> WARNING - '//trim(adjustl(MYJOBSTRING))//' exit status=',ISTAT

!call getpid_subr(CHILDPID)
!! PRINT '(A,2I8)','in SUBMITOPTIMJOB, CHILDPID,CONNID=',CHILDPID,CONNID

!WRITE(CONNSTR1,'(I10)') CHILDPID
!WRITE(CONNSTR2,'(I10)') CONNID

!CALL MYSYSTEM(ISTAT,DEBUG,'mv complex/data.' // TRIM(ADJUSTL(CONNSTR2)) // ' complex/data.' // TRIM(ADJUSTL(CONNSTR1)) )
!CALL MYSYSTEM(ISTAT,DEBUG,'mv ligand/data.' // TRIM(ADJUSTL(CONNSTR2)) // ' ligand/data.' // TRIM(ADJUSTL(CONNSTR1)) )
!CALL MYSYSTEM(ISTAT,DEBUG,'mv protein/data.' // TRIM(ADJUSTL(CONNSTR2)) // ' protein/data.' // TRIM(ADJUSTL(CONNSTR1)) )
!CALL EXIT_SUBR(ISTAT)
STOP
END SUBROUTINE SUBMITGMINJOBS

SUBROUTINE PREPAREGMINDIRS(IPARALLEL,TYPESTRING)
USE COMMONS
USE PORFUNCS
IMPLICIT NONE

CHARACTER(LEN=80) COORDSSTRING,DIRSTRING
CHARACTER(LEN=*) TYPESTRING
CHARACTER(LEN=5) ISTR
INTEGER IPARALLEL, J1, ISTAT

DO J1=1,IPARALLEL
   WRITE(ISTR,'(I5)') J1
   CALL MYSYSTEM(ISTAT,DEBUG,'mkdir -p ' // TRIM(ADJUSTL(TYPESTRING)) // '/' // TRIM(ADJUSTL(ISTR)))
END DO 

END SUBROUTINE PREPAREGMINDIRS

SUBROUTINE PREPAREDATAFILE(ICPU,COORDSSTRING,TYPESTRING)
USE COMMONS
USE PORFUNCS
IMPLICIT NONE

INTEGER ICPU, ISTAT
CHARACTER(LEN=*) TYPESTRING
CHARACTER(LEN=80) COORDSSTRING,DIRSTRING
CHARACTER(LEN=5) ICPUSTR


WRITE(ICPUSTR,'(I5)') ICPU

! copy files to subdirectories
CALL MYSYSTEM(ISTAT,DEBUG,'cp ' // TRIM(ADJUSTL(TYPESTRING)) // '/' // TRIM(ADJUSTL(COORDSSTRING)) // ' ' // &
&                TRIM(ADJUSTL(TYPESTRING)) // '/' // TRIM(ADJUSTL(ICPUSTR)) // '/coords.inpcrd')
CALL MYSYSTEM(ISTAT,DEBUG,'cp ' // TRIM(ADJUSTL(TYPESTRING)) // '/data ' // &
&                TRIM(ADJUSTL(TYPESTRING)) // '/' // TRIM(ADJUSTL(ICPUSTR)) // '/')
CALL MYSYSTEM(ISTAT,DEBUG,'cd ' // TRIM(ADJUSTL(TYPESTRING)) // '; cp ' // TRIM(ADJUSTL(COPYFILES)) // ' ' &
&                 // TRIM(ADJUSTL(ICPUSTR)) // '/; cd ..')

END SUBROUTINE PREPAREDATAFILE

SUBROUTINE SSHSUBMITGMIN(ICPU,STAT,JOBSTRING,CONNSTR1,DIRSTRING)
          USE NODES
          USE PORFUNCS
          USE COMMONS, ONLY: CHARMMT, ZSYM, DEBUG, COPYFILES, COPYOPTIMT, BHINTERPT
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: ICPU
          INTEGER,INTENT(OUT) :: STAT
          CHARACTER(LEN=*) :: JOBSTRING,CONNSTR1,DIRSTRING

          INTEGER :: STATUS(5),ISTAT
          CHARACTER(LEN=80) :: NODE

          NODE=TRIM(ADJUSTL(NODENAME(ICPU)))

          CALL MYSYSTEM(STATUS(1),DEBUG,'rsh ' // TRIM(node) // ' mkdir /scratch/' //&
                &       TRIM(ADJUSTL(USERNAME)) // '/' // connstr1)
          CALL FLUSH(6)

          CALL MYSYSTEM(STATUS(2),DEBUG,'rsync -ae rsh ' // DIRSTRING // '/* ' // TRIM(NODE) //&
                        &  ':/scratch/'// TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/')
          CALL FLUSH(6)
          CALL MYSYSTEM(STATUS(3),DEBUG,'rsh ' // TRIM(node) // ' "cd /scratch/' // TRIM(ADJUSTL(USERNAME)) // &
   &                           '/' // connstr1 // ' ; ' // jobstring // '"')
          CALL FLUSH(6)
          
!          IF (DEBUG.OR.(VERIFY('tssearch',JOBSTRING).EQ.0)) THEN ! copy everything
             CALL MYSYSTEM(STATUS(4),DEBUG,'rsync -ae rsh ' // TRIM(node) // &
                &   ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/* '// DIRSTRING //'/.')


!          ELSEIF (COPYOPTIMT) THEN ! copy path.info, OPTIM, odata and finish
!             CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                        '/' // connstr1 // '/OPTIM* .',status(4))
!             CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                        '/' // connstr1 // '/odata* .',status(4))
!             CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                        '/' // connstr1 // '/finish* .',status(4))
!             IF (BHINTERPT) THEN
!                CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                           '/' // connstr1 // '/min.data.info* .',status(4))
!             ELSE
!                CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                           '/' // connstr1 // '/path.info* .',status(4))
!             ENDIF
!          ELSE ! we only really need path.info or min.data.info 
!             IF (BHINTERPT) THEN
!                CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                           '/' // connstr1 // '/min.data.info* .',status(4))
!             ELSE
!                CALL SYSTEM('rsync -ae rsh ' // TRIM(node) // ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // &
!   &                           '/' // connstr1 // '/path.info* .',status(4))
!             ENDIF
!          ENDIF

          CALL MYSYSTEM(status(5),DEBUG,'rsh ' // TRIM(node) // ' rm -r /scratch/' // TRIM(ADJUSTL(USERNAME)) // &
   &                        '/' // connstr1)
     
          STAT=STATUS(4)
END SUBROUTINE SSHSUBMITGMIN

SUBROUTINE PREPAREOPTIMFILES
! this should dump all coordinates of the minima into ligand/start.x, coords/start.x etc.,
! also copy odata.freq file to ligand/odata.x etc.
implicit none

integer J1,J2,J3,ISTAT
character(len=10) ISTR

! copy odata.freq
DO J1=1,totalmin(1)
 WRITE(ISTR,'(I8)') J1
 CALL SYSTEM('cp odata.freq ligand/odata.'//TRIM(ADJUSTL(ISTR)),ISTAT)
 OPEN(UNIT=332,FILE='ligand/start.'//TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
 DO J2=1,natomsligand
  WRITE(332,'(3F20.10)') COORDSLIGAND(J1,3*J2-2), COORDSLIGAND(J1,3*J2-1), COORDSLIGAND(J1,3*J2)
 END DO
 CLOSE(332)
END DO
DO J1=1,totalmin(2)
 WRITE(ISTR,'(I8)') J1
 CALL SYSTEM('cp odata.freq complex/odata.'//TRIM(ADJUSTL(ISTR)),ISTAT)
 OPEN(UNIT=332,FILE='complex/start.'//TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
 DO J2=1,natomscomplex
  WRITE(332,'(3F20.10)') COORDSCOMPLEX(J1,3*J2-2), COORDSCOMPLEX(J1,3*J2-1), COORDSCOMPLEX(J1,3*J2)
 END DO
 CLOSE(332)
END DO
DO J1=1,totalmin(3)
 WRITE(ISTR,'(I8)') J1
 CALL SYSTEM('cp odata.freq protein/odata.'//TRIM(ADJUSTL(ISTR)),ISTAT)
 OPEN(UNIT=332,FILE='protein/start.'//TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
 DO J2=1,natomsprotein
  WRITE(332,'(3F20.10)') COORDSPROTEIN(J1,3*J2-2), COORDSPROTEIN(J1,3*J2-1), COORDSPROTEIN(J1,3*J2)
 END DO
 CLOSE(332)
END DO

END SUBROUTINE PREPAREOPTIMFILES

SUBROUTINE DOCKCYCLE2
! submit OPTIM frequency calculations when CPUs become available
USE COMMONS
USE PORFUNCS
IMPLICIT NONE

integer :: totaljobs,J1,J2,J3,tempj,PID(NCPU),ISTAT,piddone
logical :: killed,parent,child
integer :: temp,ncycles,newjob,status
character(len=20) :: directory

totaljobs = totalmin(1)+totalmin(2)+totalmin(3) ! total number of jobs
 
WRITE(*,*) 'total jobs in dockcycle2 =', totaljobs
CALL FLUSH(6)
!PID(:)=0
J1=0
  DO J3=1,NCPU
!  WRITE(*,*) 'j1,totaljobs',J1,totaljobs,totalmin(1:3),NCPU
    J1=J1+1
    IF (J1<totalmin(1)) THEN
        directory='ligand'
        tempj=J1
    ELSE IF(J1>totalmin(1).AND.J1<=totalmin(1)+totalmin(2)) THEN
        directory='complex'
        tempj=J1-totalmin(1)
    ELSE IF(J1>totalmin(1)+totalmin(2)) THEN
        directory='protein'
        tempj=J1-totalmin(1)-totalmin(2)
    END IF

    CALL SYSTEM('sleep 5') ! to prevent us running out of source ports
    CALL FLUSH(6) ! the child process may duplicate output without this line!
    CALL FORK_SUBR(PID(J3)) ! PID is zero in the child, non-zero in the parent
!    IF (DEBUG.AND.(PID(J3).NE.0)) WRITE(*,'(A,I16)') 'dockcycle> forked connect run process id=',PID(J3)
!    IF (PID(J3).NE.0) WRITE(*,'(A,I16)') 'dockcycle> forked connect run process id=',PID(J3)
!    IF (PID(J3).EQ.0) THEN
!       WRITE(*,'(A,2I16)') 'dockcycle> I am the child! PID=',PID(J3), J3
!    END IF  
   IF (PID(J3).EQ.0) CALL SUBMITFRQJOBS(J3,DEBUG,trim(adjustl(directory)),tempj)
        ! child processes exit at the end of SUBMITFRQJOBS
   CALL FLUSH(6)
  ENDDO

NCYCLES=0
STATUS=0
cycles: DO
  NCYCLES=NCYCLES+1
  J1=J1+1
    IF (J1<totalmin(1)) THEN
        directory='ligand'
        tempj=J1
    ELSE IF(J1>totalmin(1).AND.J1<=totalmin(1)+totalmin(2)) THEN
        directory='complex'
        tempj=J1-totalmin(1)
    ELSE IF(J1>totalmin(1)+totalmin(2)) THEN
        directory='protein'
        tempj=J1-totalmin(1)-totalmin(2)
    END IF
   KILLED=.FALSE.
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   CALL WAIT_SUBR(PIDDONE,STATUS)
!  PRINT '(A,5I8)','cycle2> PIDDONE,STATUS,PID=',PIDDONE,STATUS,PID(1:NCPU)
11 CONTINUE
   CALL FLUSH(6) ! the child process may duplicate output without this line!
   IF (PIDDONE.GT.0) THEN
      IF (DEBUG) PRINT '(A,I8,A,I6)','dockcycle2> PID ',PIDDONE,' has finished with exit status ',STATUS
      DO J2=1,NCPU
         IF (PIDDONE.EQ.PID(J2)) THEN
            IF (STATUS.NE.0) KILLED=.TRUE. ! INCOMPLETE OPTIM JOBS WOULD IDEALLY RETURN A NON-ZERO EXIT CODE 
            NEWJOB=J2
            GOTO 10
         ENDIF
      ENDDO
      PRINT*,'ERROR - PID of completed child process not recognised: ',PIDDONE
      CALL EXIT(ISTAT) 
   ELSE
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      PRINT '(A,I8)','dockcycle2> ERROR - WAIT returned system error code ',-PIDDONE
!
! Try calling wait again to see if this fixes things.
!
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I8))','dockcycle2> on calling wait again pid=',PIDDONE,' status=',STATUS
      IF (PIDDONE.GT.0) GOTO 11
      CALL EXIT(ISTAT)
   ENDIF
10 WRITE(*,'(3(A,I8))') 'dockcycle2> forked OPTIM run ',NCYCLES,' on CPU ',NEWJOB,' process id ',PID(NEWJOB)

   CALL FLUSH(6) ! the child process may duplicate output without this line!
   WRITE(*,'(3(A,I8))') 'dockcycle2> analysing result of OPTIM job ',NCYCLES,' on CPU ',NEWJOB,' for process id ',PID(NEWJOB)
!   WRITE(CONNSTR,'(I10)') PID(NEWJOB)
   IF (KILLED) WRITE(*,'(3(A,I8))') 'dockcycle2> run ',NCYCLES,' on CPU ',NEWJOB,' was unsuccessful'

   CALL FLUSH(6) ! the child process may duplicate output without this line!
   CALL FORK_SUBR(PID(NEWJOB))
   IF (PID(NEWJOB).EQ.0) CALL SUBMITFRQJOBS(NEWJOB,DEBUG,trim(adjustl(directory)),tempj)

   IF (MOD(NCYCLES,NCPU).EQ.0) THEN 
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
      WRITE(*,'(A,I8)') 'dockcycle2> end of cycle ',NCYCLES/NCPU
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
   END IF
  
   IF (J1==totaljobs) STOP
END DO cycles
 
CALL FLUSH(6)

RETURN

END SUBROUTINE DOCKCYCLE2

SUBROUTINE SUBMITFRQJOBS(ICPU,DEBUG,DIRSTRING,J1)
USE PORFUNCS
USE NODES, ONLY: SSHSUBMIT, SSHPARALLEL
USE COMMONS, ONLY: DUMMYRUNT, AMBERT, NCPU, CHARMMT, UNRST,EXEC
implicit none

integer :: J1,ICPU,ISTAT
CHARACTER(LEN=*) DIRSTRING
CHARACTER(LEN=256) MYJOBSTRING
INTEGER :: CHILDPID, CONNID
LOGICAL :: DEBUG
CHARACTER(LEN=10) CONNSTR1, CONNSTR2, J1CHAR

call getpid_subr(CHILDPID)
! PRINT '(A,2I8)','in SUBMITOPTIMJOB, CHILDPID,CONNID=',CHILDPID,CONNID
WRITE(CONNSTR1,'(I10)') CHILDPID
!WRITE(*,*) 'CONNSTR1 in submitgminjobs = ', CONNSTR1, SSHPARALLEL,DEBUG
!WRITE(CONNSTR2,'(I10)') CONNID
!STOP
WRITE(J1CHAR,'(I6)') J1
MYJOBSTRING=TRIM(ADJUSTL(EXEC))//' '//TRIM(ADJUSTL(J1CHAR))// ' > OPTIM.dock.' //TRIM(ADJUSTL(J1CHAR))

IF (SSHPARALLEL) then
!    CALL SLEEP(10)  ! this does not compile with NAG - DJW
     CALL SYSTEM('sleep 10') ! to prevent us running out of source ports
!     WRITE(*,*) 'submitting job', icpu,istat,trim(adjustl(myjobstring)),trim(adjustl(CONNSTR1)),trim(adjustl(DIRSTRING))
     call SSHsubmitfrq(ICPU,istat,trim(adjustl(myjobstring)),trim(adjustl(CONNSTR1)),trim(adjustl(DIRSTRING)),trim(adjustl(J1CHAR)))
ELSE
     CALL MYSYSTEM(ISTAT,DEBUG,'cd ' // trim(adjustl(DIRSTRING)) // ';' //trim(adjustl(myjobstring)))
ENDIF
IF (ISTAT.NE.0) PRINT '(A,I8)','submitfrqjobs> WARNING - '//trim(adjustl(MYJOBSTRING))//' exit status=',ISTAT

!call getpid_subr(CHILDPID)
!! PRINT '(A,2I8)','in SUBMITOPTIMJOB, CHILDPID,CONNID=',CHILDPID,CONNID

!WRITE(CONNSTR1,'(I10)') CHILDPID
!WRITE(CONNSTR2,'(I10)') CONNID

!CALL MYSYSTEM(ISTAT,DEBUG,'mv complex/data.' // TRIM(ADJUSTL(CONNSTR2)) // ' complex/data.' // TRIM(ADJUSTL(CONNSTR1)) )
!CALL MYSYSTEM(ISTAT,DEBUG,'mv ligand/data.' // TRIM(ADJUSTL(CONNSTR2)) // ' ligand/data.' // TRIM(ADJUSTL(CONNSTR1)) )
!CALL MYSYSTEM(ISTAT,DEBUG,'mv protein/data.' // TRIM(ADJUSTL(CONNSTR2)) // ' protein/data.' // TRIM(ADJUSTL(CONNSTR1)) )
CALL EXIT(ISTAT)






END SUBROUTINE SUBMITFRQJOBS

SUBROUTINE SSHSUBMITFRQ(ICPU,STAT,JOBSTRING,CONNSTR1,DIRSTRING,CONNSTR2)
          USE NODES
          USE PORFUNCS
          USE COMMONS, ONLY: CHARMMT, ZSYM, DEBUG, COPYFILES, COPYOPTIMT, BHINTERPT, AMBERT
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: ICPU
          INTEGER,INTENT(OUT) :: STAT
          CHARACTER(LEN=*) :: JOBSTRING,CONNSTR1,CONNSTR2,DIRSTRING

          INTEGER :: STATUS(5),ISTAT
          CHARACTER(LEN=80) :: NODE

          NODE=TRIM(ADJUSTL(NODENAME(ICPU)))

          CALL MYSYSTEM(STATUS(1),DEBUG,'rsh ' // TRIM(node) // ' mkdir /scratch/' //&
                &       TRIM(ADJUSTL(USERNAME)) // '/' // connstr1)
          CALL FLUSH(6)
          
          CALL MYSYSTEM(STATUS(2),DEBUG,'rsync -ae rsh ' // DIRSTRING // '/odata.' //TRIM(ADJUSTL(CONNSTR2)) // &
                        & ' '// TRIM(NODE) // ':/scratch/'// TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/')
          CALL MYSYSTEM(STATUS(2),DEBUG,'rsync -ae rsh ' // DIRSTRING // '/start.' //TRIM(ADJUSTL(CONNSTR2)) // &
                        & ' '// TRIM(NODE) // ':/scratch/'// TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/')
       IF(AMBERT) THEN
          CALL MYSYSTEM(STATUS(2),DEBUG,'rsync -ae rsh ' // DIRSTRING // '/coords.* ' // TRIM(NODE) //&
                        &  ':/scratch/'// TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/')
          CALL MYSYSTEM(STATUS(2),DEBUG,'rsync -ae rsh ' // DIRSTRING // '/min.in '  // TRIM(NODE) //&
                        &  ':/scratch/'// TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/')
          CALL MYSYSTEM(STATUS(2),DEBUG,'rsync -ae rsh ' // DIRSTRING // '/refc '  // TRIM(NODE) //&
                        &  ':/scratch/'// TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/')
       ELSE IF(CHARMMT) THEN
! not working yet for CHARMM runs!!!! Need to figure out what files to copy over.
       END IF
          CALL FLUSH(6)
          CALL MYSYSTEM(STATUS(3),DEBUG,'rsh ' // TRIM(node) // ' "cd /scratch/' // TRIM(ADJUSTL(USERNAME)) // &
   &                           '/' // connstr1 // ' && ' // jobstring // '"')
          CALL FLUSH(6)
          
!          IF (DEBUG.OR.(VERIFY('tssearch',JOBSTRING).EQ.0)) THEN ! copy everything

! copy back everything for now
          CALL MYSYSTEM(STATUS(4),DEBUG,'rsync -ae rsh ' // TRIM(node) // &
                &   ':/scratch/' // TRIM(ADJUSTL(USERNAME)) // '/' // connstr1 // '/* '// DIRSTRING //'/.')

          CALL SYSTEM('rsh ' // TRIM(node) // ' rm -r /scratch/' // TRIM(ADJUSTL(USERNAME)) // &
   &                        '/' // connstr1,status(5))
     
          STAT=STATUS(4)
END SUBROUTINE SSHSUBMITFRQ

SUBROUTINE ANALYSERESULTS
use porfuncs
use COMMONS, only : DEBUG,TEMPERATURE,BOXLX,BOXLY,BOXLZ,PERMDIST,PERMISOMER,EDIFFTOL,GEOMDIFFTOL,ZSYM
implicit none

LOGICAL :: YESNO
INTEGER :: J1, J2, J3, nminima(3),nlines(3),ISTAT,natomsdummy,iostatus,HORDERMIN,&
           & minindex(3),prmlines,nlines2
DOUBLE PRECISION :: IXMIN,IYMIN,IZMIN,ORDERMIN,MINE(3),DISTANCE,DIST2
DOUBLE PRECISION, ALLOCATABLE :: FVIBMIN(:),EMIN(:)
CHARACTER(LEN=10) :: DIRSTR,PRMFORMAT
CHARACTER(LEN=80) :: check,dummyline
CHARACTER(len=4), allocatable :: atomlabels(:)
! concatenate min.data.info files in the three directories

   INQUIRE(file='ligand/min.data.info.save',exist=YESNO)
   IF(YESNO) THEN 
        CALL MYSYSTEM(ISTAT,.TRUE.,'rm ligand/min.data.info.save')
   END IF
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd ligand; cat min.data.info.[1-9]* >min.data.info.save')
   INQUIRE(file='complex/min.data.info.save',exist=YESNO)
   IF(YESNO) CALL MYSYSTEM(ISTAT,.TRUE.,'rm complex/min.data.info.save')
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd complex; cat min.data.info.[1-9]* >min.data.info.save')
   INQUIRE(file='protein/min.data.info.save',exist=YESNO)
   IF(YESNO) CALL MYSYSTEM(ISTAT,.TRUE.,'rm protein/min.data.info.save')
   CALL MYSYSTEM(ISTAT,.TRUE.,'cd protein; cat min.data.info.[1-9]* >min.data.info.save')

! deallocate ENERGIES... arrays
  if(allocated(energiesligand)) deallocate(energiesligand)
  if(allocated(energiescomplex)) deallocate(energiescomplex)
  if(allocated(energiesprotein)) deallocate(energiesprotein)
  if(allocated(energiesligandsave)) deallocate(energiesligandsave)
  if(allocated(energiescomplexsave)) deallocate(energiescomplexsave)
  if(allocated(energiesproteinsave)) deallocate(energiesproteinsave)
  if(allocated(feligand)) deallocate(feligand)
  if(allocated(fecomplex)) deallocate(fecomplex)
  if(allocated(feprotein)) deallocate(feprotein)



! determine number of minima in the min.data.info.save file
DO j1=1,3
 if(j1==1) then
        natomsdummy=natomsligand
        dirstr='ligand'
 else if(j1==2) then
        natomsdummy=natomscomplex
        dirstr='complex'
 else
        natomsdummy=natomsprotein
        dirstr='protein'
 end if
         OPEN(unit=311,file=trim(adjustl(dirstr))//'/min.data.info.save',status='OLD')
         nlines(j1) = 0
         do
           read (unit=311,iostat=iostatus,fmt='(A15)') check
             if (iostatus<0) then
                if (debug) write(*,'(A)') 'analyseresults> finished reading file ' // &
                        &        trim(adjustl(dirstr)) //'/min.data.info.save'
                rewind(311)
                exit
             end if
           nlines(j1) = nlines(j1) + 1
         end do
         nminima(j1) =  nlines(j1)/(natomsdummy+1)     
         CLOSE(311)
END DO
! allocate energies and vibrational contribution arrays
  allocate(energiesligand( nminima(1)))
  allocate(energiescomplex(nminima(2)))
  allocate(energiesprotein(nminima(3)))
  allocate(energiesligandsave( nminima(1)))
  allocate(energiescomplexsave(nminima(2)))
  allocate(energiesproteinsave(nminima(3)))
  allocate( fvibminligand( nminima(1)))
  allocate( fvibmincomplex(nminima(2)))
  allocate( fvibminprotein(nminima(3)))
  allocate(lpligand(3*natomsligand),lpcomplex(3*natomscomplex),lpprotein(3*natomsprotein))

! now open min.data and points.min files

! LIGAND
  OPEN(unit=311,file='ligand/min.data.info.save',status='OLD')
  OPEN(UNIT=312,FILE='min.data.ligand',STATUS='UNKNOWN')
  OPEN(UNIT=313,FILE='points.min.ligand',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*natomsligand)

CALL determinelines(312,nlines(1))
ALLOCATE(FVIBMIN(nlines(1)),EMIN(nlines(1)))
!write(*,*) 'nlines determined is', nlines(:)
!write(*,*) 'nminima determined is', nminima(:)
do j2=1,nlines(1)
  READ(312,*) EMIN(j2),fvibmin(j2),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
end do

nlines2=0
do j1=1,nminima(1)
  READ(311,*) energiesligand(j1),fvibminligand(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
  READ(311,*) (COORDSLIGAND(1,J2),J2=1,3*natomsligand)
!  WRITE(*,*) 'checking whether minima are in the database', energiesligand(1:10),edifftol,EMIN(1)
! check whether minima are already in the database
  DO j2=1,nlines(1)
     DISTANCE=1.0D100
     IF (ABS(energiesligand(j1)-EMIN(J2)).LT.EDIFFTOL) THEN
        READ(313,REC=J2) (LPLIGAND(J3),J3=1,3*NATOMSLIGAND)
        CALL MINPERMDIST(COORDSLIGAND(1,:),LPLIGAND,NATOMSLIGAND,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,.FALSE.,DISTANCE,DIST2,.FALSE., &
  &                      RMAT,.FALSE.)
     ENDIF
     IF ((ABS(energiesligand(j1)-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
        PRINT '(A,I6)','analyseresults> ligand minimum is database minimum ',J2
        IF (ABS(fvibminligand(J1)-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
           WRITE(*,'(A,F15.5,A,F15.5)') 'analyseresults> WARNING, NEWFVIBMIN=',fvibminligand(J1),' should be ',FVIBMIN(J2)
        ENDIF
        GOTO 140
     ENDIF
  END DO
  nlines2=nlines2+1
!  WRITE(312,'(2F20.10,I6,4F20.10)') energiesligand(j1),fvibminligand(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
!  WRITE(313,REC=nlines(1)+nlines2) (COORDSLIGAND(1,J2),J2=1,3*natomsligand)
140 CONTINUE
end do
  close(311)
  close(312)
  close(313)
DEALLOCATE(FVIBMIN,EMIN)
! COMPLEX
  OPEN(unit=311,file='complex/min.data.info.save',status='OLD')
  OPEN(UNIT=312,FILE='min.data.complex',STATUS='UNKNOWN')
  OPEN(UNIT=313,FILE='points.min.complex',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*natomscomplex)

CALL determinelines(312,nlines(2))
ALLOCATE(FVIBMIN(nlines(2)),EMIN(nlines(2)))
!write(*,*) 'nlines determined is', nlines(:)
!write(*,*) 'nminima determined is', nminima(:)
do j2=1,nlines(2)
  READ(312,*) EMIN(j2),fvibmin(j2),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
end do

nlines2=0
do j1=1,nminima(2)
  READ(311,*) energiescomplex(j1),fvibmincomplex(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
  READ(311,*) (COORDSCOMPLEX(1,J2),J2=1,3*natomscomplex)
!  WRITE(*,*) 'checking whether minima are in the database', energiesligand(1:10),edifftol,EMIN(1)
! check whether minima are already in the database
  DO j2=1,nlines(2)
     DISTANCE=1.0D100
     IF (ABS(energiescomplex(j1)-EMIN(J2)).LT.EDIFFTOL) THEN
        READ(313,REC=J2) (LPCOMPLEX(J3),J3=1,3*NATOMSCOMPLEX)
        CALL MINPERMDIST(COORDSCOMPLEX(1,:),LPCOMPLEX,NATOMSCOMPLEX,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,.FALSE.,DISTANCE, &
  &                      DIST2,.FALSE., &
  &                      RMAT,.FALSE.)
     ENDIF
     IF ((ABS(energiescomplex(j1)-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
        PRINT '(A,I6)','analyseresults> complex minimum is database minimum ',J2
        IF (ABS(fvibmincomplex(J1)-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
           WRITE(*,'(A,F15.5,A,F15.5)') 'analyseresults> WARNING, NEWFVIBMIN=',fvibmincomplex(J1),' should be ',FVIBMIN(J2)
        ENDIF
        GOTO 150
     ENDIF
  END DO
  nlines2=nlines2+1
!  WRITE(312,'(2F20.10,I6,4F20.10)') energiescomplex(j1),fvibmincomplex(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
!  WRITE(313,REC=nlines(2)+nlines2) (COORDSCOMPLEX(1,J2),J2=1,3*natomscomplex)
150 CONTINUE
end do
  close(311)
  close(312)
  close(313)
DEALLOCATE(FVIBMIN,EMIN)
! PROTEIN 
  OPEN(unit=311,file='protein/min.data.info.save',status='OLD')
  OPEN(UNIT=312,FILE='min.data.protein',STATUS='UNKNOWN')
  OPEN(UNIT=313,FILE='points.min.protein',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*natomsprotein)


CALL determinelines(312,nlines(3))
ALLOCATE(FVIBMIN(nlines(3)),EMIN(nlines(3)))
!write(*,*) 'nlines determined is', nlines(:)
!write(*,*) 'nminima determined is', nminima(:)
do j2=1,nlines(3)
  READ(312,*) EMIN(j2),fvibmin(j2),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
end do

nlines2=0
do j1=1,nminima(3)
  READ(311,*) energiesprotein(j1),fvibminprotein(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
  READ(311,*) (COORDSPROTEIN(1,J2),J2=1,3*natomsprotein)
! check whether minima are already in the database
  DO j2=1,nlines(3)
     DISTANCE=1.0D100
     IF (ABS(energiesprotein(j1)-EMIN(J2)).LT.EDIFFTOL) THEN
        READ(313,REC=J2) (LPPROTEIN(J3),J3=1,3*NATOMSPROTEIN)
        CALL MINPERMDIST(COORDSPROTEIN(1,:),LPCOMPLEX,NATOMSCOMPLEX,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,.FALSE., &
  &                      DISTANCE,DIST2,.FALSE.,RMAT,.FALSE.)
     ENDIF
     IF ((ABS(energiesprotein(j1)-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
        PRINT '(A,I6)','analyseresults> protein minimum is database minimum ',J2
        IF (ABS(fvibminprotein(J1)-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
           WRITE(*,'(A,F15.5,A,F15.5)') 'analyseresults> WARNING, NEWFVIBMIN=',fvibminprotein(J1),' should be ',FVIBMIN(J2)
        ENDIF
        GOTO 160
     ENDIF
  END DO
  nlines2=nlines2+1
!  WRITE(312,'(2F20.10,I6,4F20.10)') energiesprotein(j1),fvibminprotein(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
!  WRITE(313,REC=nlines(3)+nlines2) (COORDSPROTEIN(1,J2),J2=1,3*natomsprotein)
160 CONTINUE
end do
  close(311)
  close(312)
  close(313)
deallocate(energiesligand,fvibminligand,energiescomplex,fvibmincomplex,energiesprotein,fvibminprotein)
DEALLOCATE(FVIBMIN,EMIN)
100 continue
OPEN(unit=211,file='min.data.ligand',status='old')
OPEN(unit=212,file='min.data.complex',status='old')
OPEN(unit=213,file='min.data.protein',status='old')

nlines(:)=0

! can't determine lines right after opening a file!!! seems to need some time until the file is being opened,
! behaviour probably due to NFS files
CALL SYSTEM('sleep 2')
CALL determinelines(211,nlines(1))
CALL SYSTEM('sleep 2')
CALL determinelines(212,nlines(2))
CALL SYSTEM('sleep 2')
CALL determinelines(213,nlines(3))

write(*,*) 'lines determined:',nlines(:)
allocate(energiesligand(nlines(1)),fvibminligand(nlines(1)))
allocate(energiescomplex(nlines(2)),fvibmincomplex(nlines(2)))
allocate(energiesprotein(nlines(3)),fvibminprotein(nlines(3)))
allocate(feligand( nlines(1)),relfeligand(nlines(1)),wligand(nlines(1)))
allocate(fecomplex(nlines(2)),relfecomplex(nlines(2)),wcomplex(nlines(2)))
allocate(feprotein(nlines(3)),relfeprotein(nlines(3)),wprotein(nlines(3)))
allocate(indexarrayl(nlines(1)),indexarrayc(nlines(2)),indexarrayp(nlines(3))  )
allocate(indexarraylfe(nlines(1)),indexarraycfe(nlines(2)),indexarraypfe(nlines(3))  )

! read in energies and log products from the min.data.info files

do j1=1,nlines(1)
   READ(211,*) energiesligand(j1),fvibminligand(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
end do

do j1=1,nlines(2)
   READ(212,*) energiescomplex(j1),fvibmincomplex(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
end do

do j1=1,nlines(3)
   READ(213,*) energiesprotein(j1),fvibminprotein(j1),HORDERMIN,IXMIN,IYMIN,IZMIN!,ORDERMIN
end do

close(211)
close(212)
close(213)

MINE(:)=1.0D10
minindex(:)=1
! determine lowest minima once again, save indices

do j1=1,nlines(1)
   feligand(j1)=energiesligand(j1)+TEMPERATURE*fvibminligand(j1) / 2.0D0
   if(energiesligand(j1)<MINE(1)) then
        MINE(1)=energiesligand(j1)
        minindex(1)=j1
   end if
end do
do j1=1,nlines(2)
   fecomplex(j1)=energiescomplex(j1)+TEMPERATURE*fvibmincomplex(j1) / 2.0D0  
   if(energiescomplex(j1)<MINE(2)) then
        MINE(2)=energiescomplex(j1)
        minindex(2)=j1
   end if
end do
do j1=1,nlines(3)
   feprotein(j1)=energiesprotein(j1)+TEMPERATURE*fvibminprotein(j1) / 2.0D0   
   if(energiesprotein(j1)<MINE(3)) then 
        MINE(3)=energiesprotein(j1)
        minindex(3)=j1
   end if
end do


! now calculate relative free energies for each minima, using the harmonic approximation
! all free energies are relative to the global energy minimum for each system 
! A = E + TEMP*FVIBMIN/2 - FE0, where TEMP is in units of kT

! relative free energies and Boltzmann weights
meane(:)=0.0D0
totalw(:)=0.0D0
do j1=1,nlines(1)
   relfeligand(j1)=feligand(j1)-feligand(minindex(1))
   wligand(j1)=EXP(-relfeligand(j1)/TEMPERATURE)
   totalw(1)=totalw(1)+wligand(j1)
end do
do j1=1,nlines(2)
   relfecomplex(j1)=fecomplex(j1)-fecomplex(minindex(2))
   wcomplex(j1)=EXP(-relfecomplex(j1)/TEMPERATURE)
   totalw(2)=totalw(2)+wcomplex(j1)
end do
do j1=1,nlines(3)
   relfeprotein(j1)=feprotein(j1)-feprotein(minindex(3))
   wprotein(j1)=EXP(-relfeprotein(j1)/TEMPERATURE)
   totalw(3)=totalw(3)+wprotein(j1)
end do
do j1=1,nlines(1)
   meane(1)=meane(1)+wligand(j1)*energiesligand(j1)/totalw(1)
end do
do j1=1,nlines(2)
   meane(2)=meane(2)+wcomplex(j1)*energiescomplex(j1)/totalw(2)
end do
do j1=1,nlines(3)
   meane(3)=meane(3)+wprotein(j1)*energiesprotein(j1)/totalw(3)
end do

! dump relative free energies and energies (to create histograms)

OPEN(unit=331,file='ligand.csv',status='unknown')
OPEN(unit=332,file='complex.csv',status='unknown')
OPEN(unit=333,file='protein.csv',status='unknown')
OPEN(unit=334,file='ligand_fe.csv',status='unknown')
OPEN(unit=335,file='complex_fe.csv',status='unknown')
OPEN(unit=336,file='protein_fe.csv',status='unknown')
energiesligandsave(:)=energiesligand(:)
energiescomplexsave(:)=energiescomplex(:)
energiesproteinsave(:)=energiesprotein(:)
CALL Insrt(energiesligand,indexarrayl)
CALL Insrt(energiescomplex,indexarrayc)
CALL Insrt(energiesprotein,indexarrayp)

WRITE(331,'(2F20.10,I6)') (relfeligand(indexarrayl(J1)),energiesligand(j1),indexarrayl(j1),J1=1,nminima(1))
WRITE(332,'(2F20.10,I6)') (relfecomplex(indexarrayc(J1)),energiescomplex(j1),indexarrayc(j1),J1=1,nminima(2))
WRITE(333,'(2F20.10,I6)') (relfeprotein(indexarrayp(J1)),energiesprotein(j1),indexarrayp(j1),J1=1,nminima(3))

CALL Insrt(relfeligand,indexarraylfe)
CALL Insrt(relfecomplex,indexarraycfe)
CALL Insrt(relfeprotein,indexarraypfe)

WRITE(334,'(2F20.10,I6)') (relfeligand(j1),energiesligandsave(indexarraylfe(J1)),indexarraylfe(j1),J1=1,nminima(1))
WRITE(335,'(2F20.10,I6)') (relfecomplex(j1),energiescomplexsave(indexarraycfe(J1)),indexarraycfe(j1),J1=1,nminima(2))
WRITE(336,'(2F20.10,I6)') (relfeprotein(j1),energiesproteinsave(indexarraypfe(J1)),indexarraypfe(j1),J1=1,nminima(3))

CLOSE(331)
CLOSE(332)
CLOSE(333)
CLOSE(334)
CLOSE(335)
CLOSE(336)

WRITE(*,*) 'Free energies:',mine(3),nminima(:),minindex(3),feprotein(minindex(3))
!write(*,*) feprotein(:)
!WRITE(*,*) 'Weights:'
!write(*,*) wligand(:),wcomplex(:),wprotein(:)
!fecomplex(:),feprotein(:)

WRITE(*,'(A,3F20.10)') 'Binding energy based on Boltzmann-weighted free energies: ', MEANE(2)-MEANE(1)-MEANE(3)
!WRITE(*,*) 'Lowest free energies for ligand, complex and protein: ',MINFE(:)
Write(*,'(A,3F20.10)') 'Weighted mean energies: ', meane(:)

! dump lowest energy structures for the complex into xyz format
ALLOCATE(ATOMLABELS(NATOMSCOMPLEX))
! get the atom labels
       PRMLINES=(NATOMSCOMPLEX-MOD(NATOMSCOMPLEX,20))/20
       WRITE(PRMFORMAT,'(A,I2,A)') "(",MAX(20,MOD(NATOMSCOMPLEX,20)),"A4)"
       OPEN(UNIT=4,FILE='complex/coords.prmtop',STATUS='OLD')
       DO
        READ(4,'(A)',IOSTAT=ISTAT) DUMMYLINE
        IF(TRIM(ADJUSTL(DUMMYLINE))=='%FLAG ATOM_NAME') THEN
           READ(4,'(A)') DUMMYLINE
           EXIT
        END IF
        IF(ISTAT<0) EXIT
       END DO
       DO J1=1,PRMLINES
        READ(4,'(20A4)') ATOMLABELS(20*(J1-1)+1:20*J1)
       END DO
       IF(MOD(NATOMSCOMPLEX,20)/=0) THEN
       WRITE(PRMFORMAT,'(A,I2,A)') "(",MOD(NATOMSCOMPLEX,20),"A4)"
        WRITE(*,*) 'prmformat=',prmformat
       READ(4,TRIM(ADJUSTL(PRMFORMAT))) ATOMLABELS(PRMLINES*20+1:PRMLINES*20+MOD(NATOMSCOMPLEX,20))
        write(*,*) ATOMLABELS(PRMLINES+1:PRMLINES+MOD(NATOMSCOMPLEX,20))
!       READ(4,TRIM(ADJUSTL(PRMFORMAT))) ATOMLABELS(PRMLINES+1:PRMLINES+MOD(NATOMS,20))
       END IF
       WRITE(*,*) 'ATOMLABELS='
       WRITE(*,*) ATOMLABELS(1:NATOMSCOMPLEX)

  OPEN(UNIT=312,FILE='min.data.complex',STATUS='UNKNOWN')
  OPEN(UNIT=313,FILE='points.min.complex',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*natomscomplex)
  OPEN(UNIT=314,FILE='complex.xyz',STATUS='UNKNOWN')
  OPEN(UNIT=315,FILE='complex_fe.xyz',STATUS='UNKNOWN')

!write(*,*) 'indexarrayc: ', indexarrayc(:)
J2=1000
IF(nlines(2)<1000) J2=nlines(2)
DO J1=1,J2
   READ(313,rec=indexarrayc(J1)) (LPCOMPLEX(J3),J3=1,3*NATOMSCOMPLEX) 
   WRITE(314,*) natomscomplex
   WRITE(314,'(A,I6)') 'minimum ', indexarrayc(J1)
   WRITE(314,'(A4,3F20.10)') (ATOMLABELS(J3),LPCOMPLEX(3*J3-2),LPCOMPLEX(3*J3-1),LPCOMPLEX(3*J3),J3=1,NATOMSCOMPLEX)
END DO
J2=1000
IF(nlines(2)<1000) J2=nlines(2)
DO J1=1,J2
   READ(313,rec=indexarraycfe(J1)) (LPCOMPLEX(J3),J3=1,3*NATOMSCOMPLEX) 
   WRITE(315,*) natomscomplex
   WRITE(315,'(A,I6)') 'minimum ', indexarraycfe(J1)
   WRITE(315,'(A4,3F20.10)') (ATOMLABELS(J3),LPCOMPLEX(3*J3-2),LPCOMPLEX(3*J3-1),LPCOMPLEX(3*J3),J3=1,NATOMSCOMPLEX)
END DO


close(312)
close(313)
close(314)
close(315)
DEALLOCATE(ATOMLABELS)


END SUBROUTINE ANALYSERESULTS

SUBROUTINE DETERMINELINES(nunit,nlines)
implicit none
integer nunit, nlines, iostatus
character(len=10) check

REWIND(nunit)

nlines=0
do
  IF(iostatus<0) EXIT
  nlines = nlines + 1
  READ(nunit,*,iostat=iostatus) check
!  write(*,*) check,nunit
end do
  nlines = nlines - 1
  REWIND(nunit)
RETURN


END SUBROUTINE DETERMINELINES

! sorting subroutine follows, taken from 
! http://en.wikibooks.org/wiki/Algorithm_implementation/Sorting/Insertion_sort#Fortran_90.2F95

! ***********************************
! *
  Subroutine Insrt(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascending order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. 
! ***********************************
 
    Double precision, Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Double precision :: Rtmp
    Integer :: I, J
 
 
    If (Present(Ipt)) Then
       Forall (I=1:Size(X)) Ipt(I) = I
 
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
                CALL InsrtSwap(Ipt, J, J+1)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    Else
       Do I = 2, Size(X)
          Rtmp = X(I)
          Do J = I-1, 1, -1
             If (Rtmp < X(J)) Then
                X(J+1) = X(J)
             Else
                Exit
             End If
          End Do
          X(J+1) = Rtmp
       End Do
    End If
 
    Return
  End Subroutine Insrt
 
! ***********************************
! *
  Subroutine InsrtSwap(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Integer :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine InsrtSwap

END MODULE DOCKMODULE
