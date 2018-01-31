!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE NEWCONNECTMODULE
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE NEWCONNECT(NA,EII,QQ,EFF,FINFIN,ENDPOINTSEP,PTEST,REDOPATH,REDOPATHXYZ)
          USE CONNECTDATA
          USE KEYCONNECT
          USE CONNECTUTILS
          USE DECIDEWHATTOCONNECT
          USE TRYCONNECTMODULE
          USE IDOMODULE
          USE KEY,ONLY : BHDISTTHRESH, BHINTERPT, BHDEBUG
          IMPLICIT NONE

          INTEGER,INTENT(IN)              :: NA
          DOUBLE PRECISION,INTENT(IN)           :: ENDPOINTSEP,EII,EFF,QQ(3*NA),FINFIN(3*NA)
          LOGICAL,INTENT(IN)              :: PTEST
          DOUBLE PRECISION,POINTER              :: EI,EF
          DOUBLE PRECISION,POINTER,DIMENSION(:) :: Q,FIN

          INTEGER :: JS,JF,JS2,J1,J2,NSTART,POSITION
          CHARACTER(LEN=132) :: STR
          LOGICAL REDOPATH, REDOPATHXYZ, YESNO, SUCCESS, MINNEW, PERMUTE
          DOUBLE PRECISION TSREDO(3*NA), DSTART, DFINISH
          DOUBLE PRECISION, POINTER :: PINTERPCOORDS(:), PENERGY
          DOUBLE PRECISION INTERPCOORDS(3*NA), ENERGY
          DOUBLE PRECISION CSTART(3*NA), CFINISH(3*NA)
          INTEGER, ALLOCATABLE :: TEMPDIJPAIR(:,:)
          DOUBLE PRECISION, ALLOCATABLE :: TEMPDIJPAIRDIST(:)
          INTEGER INVERT, INDEX(NA), IMATCH

          FINISHED=.FALSE.
          ALLOCATE(EI,EF,Q(3*NA),FIN(3*NA))
          EI=EII;EF=EFF;Q=QQ;FIN=FINFIN;

          MOREPRINTING=PTEST
          IF (MOREPRINTING) THEN
               CALL ALLKEYCONNECTPRINT
               PRINT*
          ENDIF
          CALL INITIALISE(NA,EI,Q,EF,FIN,ENDPOINTSEP)
          INQUIRE(FILE='redopoints',EXIST=YESNO)
          IF (YESNO) THEN
             IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                PRINT '(A)',' newconnect> Transition state coordinates will be read from file redopoints'
!            ELSE
!               PRINT '(A)',' newconnect> WARNING - redopoints file present, but no REDOPATH keyword'
             ENDIF
          ELSE
             IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                PRINT '(A)',' newconnect> WARNING - REDOPATH keyword was specified, but no redopoints file'
                REDOPATH=.FALSE.
             ENDIF
          ENDIF

          IF (REDOPATHXYZ) PRINT '(A)',' newconnect> Redo run will use available path.<n>.xyz files'
          IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) OPEN(99,FILE='redopoints',STATUS='OLD')

          DO NCONDONE=1,NCONMAX
               WRITE(CHR,'(i5)') NConDone
               WRITE(STR,'(3a)') '>>>>>>>>>>>>>>>>>>>>> CONNECT CYCLE ',trim(adjustl(chr)),' >>>>>>>>>>>>>>>>>>>>>'
               WRITE(CHR,'(i5)') Nmin
               WRITE(STR,'(a)')  trim(str)//' '//trim(adjustl(chr))//' minima and'
               WRITE(CHR,'(i5)') Nts
               WRITE(STR,'(a)')  trim(str)//' '//trim(adjustl(chr))//' ts are known'
               WRITE(*,'(1x,a)') trim(str)//' '//repeat('>',107-len_trim(str))
!
!  Read in transition state coordinates from redopoints if available.
!
               IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                  READ(99,*,END=32) (TSREDO(J1),J1=1,3*NA)
                  WRITE(*,'(A)') ' newconnect> transition state coordinates read from file redopoints'
                  GOTO 33
32                WRITE(*,'(A)') ' newconnect> no more transition state coordinates in redopoints'
                  REDOPATH=.FALSE.
                  CLOSE(99)
33                CONTINUE
               ENDIF

               IF (REDOPATH) THEN
                  CALL TRYCONNECT(1,2,TSREDO,REDOPATH,REDOPATHXYZ)
               ELSE
                  CALL DECIDE

!    call tryconnect for each pair of minima specified in the shortest path returned by DIJKSTRA
!    Call BHINTERP to fill in addtional minima if required.

                  IF (BHINTERPT) THEN
                     NSTART=1
34                   CONTINUE
                     DO J1=NSTART,NDIJPAIRS
!                       PRINT '(A)',' newconnect> DIJPAIR distances and minima:'
!                       PRINT '(I8,G20.10,2I8)',(J2,DIJPAIRDIST(J2),DIJPAIR(J2,1),DIJPAIR(J2,2),J2=1,NDIJPAIRS)
                        PRINT '(A,I8,A,I8,A,I8,A,G20.10)',' newconnect> gap ',J1,' BH interpolation for minima   ',DIJPAIR(J1,1), &
  &                                    ' and ',DIJPAIR(J1,2),' dist=',DIJPAIRDIST(J1)
                        IF (DIJPAIRDIST(J1).GT.BHDISTTHRESH) THEN
                           CSTART(1:3*NA)=MI(DIJPAIR(J1,1))%DATA%X(1:3*NA)
                           CFINISH(1:3*NA)=MI(DIJPAIR(J1,2))%DATA%X(1:3*NA)
                           CALL BHINTERP(CSTART,CFINISH,3*NA,NA,INTERPCOORDS,SUCCESS,DSTART,DFINISH,ENERGY)
                           IF ((DSTART.GT.DIJPAIRDIST(J1)).OR.(DFINISH.GT.DIJPAIRDIST(J1))) THEN
                              IF (SUCCESS.AND.BHDEBUG) PRINT '(A)',' newconnect> BH interpolated distance is larger - reject'
                              SUCCESS=.FALSE.
                           ENDIF
                           IF (SUCCESS) THEN
!
! If the extra minimum is new add it to the stack along with all its distances.
! Must also add it to the DIJPAIR list.
!
                              NULLIFY(PINTERPCOORDS,PENERGY)
                              ALLOCATE(PINTERPCOORDS(3*NA),PENERGY)
!
! Here we are reading data into the pointers PENERGY and PINTERPCOORDS so that when
! they are subsequently nullified the data in MI%DATA is preserved. 
! This is a truly horrible hack caused by the inappropriate pointer attribute for static data.
!
                              OPEN(UNIT=781,FILE='BHscratch',STATUS='UNKNOWN')
                              WRITE(781,*) ENERGY,INTERPCOORDS
                              REWIND(781)
                              READ(781,*) PENERGY,PINTERPCOORDS
                              CLOSE(781)

                              CALL ISNEWMIN(PENERGY,PINTERPCOORDS,POSITION,MINNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
                              IF (MINNEW) THEN
                                 WRITE(*,'(A,I7)') ' newconnect> Interpolated minimum is new minimum ',POSITION
                                 CALL ADDNEWMIN(PENERGY,PINTERPCOORDS)
                              ELSE
                                 WRITE(*,'(A,I7)') ' newconnect> Interpolated minimum is old minimum ',POSITION
                              ENDIF
                              NULLIFY(PINTERPCOORDS,PENERGY) 

                              ALLOCATE(TEMPDIJPAIR(NDIJPAIRS,2),TEMPDIJPAIRDIST(NDIJPAIRS))
                              TEMPDIJPAIR(1:NDIJPAIRS,1:2)=DIJPAIR(1:NDIJPAIRS,1:2)
                              TEMPDIJPAIRDIST(1:NDIJPAIRS)=DIJPAIRDIST(1:NDIJPAIRS)
                              DEALLOCATE(DIJPAIR,DIJPAIRDIST)
                              NDIJPAIRS=NDIJPAIRS+1
                              ALLOCATE(DIJPAIR(NDIJPAIRS,2),DIJPAIRDIST(NDIJPAIRS))
                              DO J2=1,J1-1
                                 DIJPAIR(J2,1:2)=TEMPDIJPAIR(J2,1:2)
                                 DIJPAIRDIST(J2)=TEMPDIJPAIRDIST(J2)
                              ENDDO
                              DIJPAIR(J1,1)=TEMPDIJPAIR(J1,1)
                              DIJPAIR(J1,2)=POSITION
                              DIJPAIRDIST(J1)=DSTART
                              DIJPAIR(J1+1,1)=POSITION
                              DIJPAIR(J1+1,2)=TEMPDIJPAIR(J1,2)
                              DIJPAIRDIST(J1+1)=DFINISH
                              DO J2=J1+1,NDIJPAIRS-1
                                 DIJPAIR(J2+1,1:2)=TEMPDIJPAIR(J2,1:2)
                                 DIJPAIRDIST(J2+1)=TEMPDIJPAIRDIST(J2)
                              ENDDO
                              NSTART=J1
                              DEALLOCATE(TEMPDIJPAIR,TEMPDIJPAIRDIST)
                              GOTO 34 ! go back, since NDIJPAIRS has increased by one
                           ELSE
                              NSTART=J1+1 ! to avoid an infinite loop
                           ENDIF
                        ENDIF
                     ENDDO
                     WRITE(*,'(A)') ' newconnect> The unconnected minima in the chain and their distances are now:'
                     DO J1=1,NDIJPAIRS
                        WRITE(*,'(I6,F12.2,I6)',ADVANCE="NO") DIJPAIR(J1,1),DIJPAIRDIST(J1),DIJPAIR(J1,2)
                     ENDDO
                  ENDIF
                  WRITE(*,'(A)') ' ' ! to advance to the next line

                  DO J1=1,NDIJPAIRS
                     JS=MAX(DIJPAIR(J1,1),DIJPAIR(J1,2))
                     JF=MIN(DIJPAIR(J1,1),DIJPAIR(J1,2))
                     WRITE(*,'(A,I5,A,2I5)') ' newconnect> trying DIJKSTRA pair number ',J1,' minima ',js,jf
                     MI(JS)%DATA%NTRIES(JF) = MI(JS)%DATA%NTRIES(JF) + 1
!    Set edge weight to infinity if we have reached the maximum number of tries for this pair.
                     CALL TRYCONNECT(JF,JS,TSREDO,REDOPATH,REDOPATHXYZ)
                     IF ((MI(JS)%DATA%NTRIES(JF).GE.NTRIESMAX).AND.  &
                         (MI(JS)%DATA%D(JF).GT.1.0D-5)) MI(JS)%DATA%D(JF) = HUGE(MI(JS)%DATA%D(JF))
                  ENDDO

                  DEALLOCATE(DIJPAIR,DIJPAIRDIST)
               ENDIF

               IF (FINISHED) EXIT
          ENDDO

          CALL OUTPUT
          CALL DEINITIALISE
          DEALLOCATE(Q,FIN,EI,EF)
     END SUBROUTINE NEWCONNECT

END MODULE NEWCONNECTMODULE
