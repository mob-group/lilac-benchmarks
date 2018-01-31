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

!    CALL       NEWCONNECT(NATOMS,EINITIAL,Q,EFINAL,FIN,DIST,.TRUE.,REDOPATH,REDOPATHXYZ)
     SUBROUTINE NEWCONNECT(NA,EII,QQ,EFF,FINFIN,ENDPOINTSEP,PTEST,REDOPATH,REDOPATHXYZ,LNOPT)

          USE CONNECTDATA
          USE KEYCONNECT
          USE CONNECTUTILS
          USE MODCHARMM
          USE DECIDEWHATTOCONNECT
          USE TRYCONNECTMODULE
          USE IDOMODULE
          USE PORFUNCS
          USE AMHGLOBALS, ONLY : NMRES
          USE KEY,ONLY : BHDISTTHRESH, BHINTERPT, BHDEBUG, DIJKSTRALOCAL, DUMPDATAT, REOPTIMISEENDPOINTS, &
  &                      AMHT, SEQ, MIN1REDO, MIN2REDO, REDOE1, REDOE2, BULKT, TWOD, PERMDIST, RIGIDBODY, &
  &                      INTCONSTRAINTT, INTLJT, INTTST, INTNTRIESMAX, QCICYCLEST, QCICYCN, QCICYCDIST,   &
  &                      QCIDNEBT
          USE COMMONS,ONLY : PARAM1, PARAM2, PARAM3, ZSYM
          USE MODNEB,ONLY :  NEWCONNECTT
          USE GENRIGID
          IMPLICIT NONE

          INTEGER LNOPT
          INTEGER,INTENT(IN)              :: NA
          DOUBLE PRECISION           :: ENDPOINTSEP,EII,EFF
          DOUBLE PRECISION           QQ(LNOPT),FINFIN(LNOPT)
          LOGICAL,INTENT(IN)              :: PTEST
          DOUBLE PRECISION,POINTER              :: EI,EF
          DOUBLE PRECISION,POINTER,DIMENSION(:) :: Q,FIN

          INTEGER :: JS,JF,J1,J2,NSTART,POSITION,J3,NCOUNT
          CHARACTER(LEN=132) :: STR
          LOGICAL REDOPATH, REDOPATHXYZ, YESNO, SUCCESS, MINNEW, PERMUTE, CHANGED, BIGGERGAP, NOPRINT
          DOUBLE PRECISION TSREDO(LNOPT), DSTART, DFINISH
          DOUBLE PRECISION, POINTER :: PINTERPCOORDS(:), PENERGY
          DOUBLE PRECISION INTERPCOORDS(LNOPT), ENERGY, OLDDISTS, OLDDISTF
          DOUBLE PRECISION CSTART(3*NA), CFINISH(LNOPT), ESTART, EFINISH, RMSINITIAL, RMSFINAL, LGDUMMY(3*NA), ETS, RMSTS
          INTEGER, ALLOCATABLE :: TEMPDIJPAIR(:,:)
          DOUBLE PRECISION, ALLOCATABLE :: TEMPDIJPAIRDIST(:)
          INTEGER INVERT, INDEX(NA), IMATCH, NMINSAVE, NMINSAVE2, ISTAT, MYJS, MYJF
          CHARACTER(LEN=80) FNAMEF
          CHARACTER(LEN=20) EFNAME
          DOUBLE PRECISION BHENERGY
          COMMON /BHINTE/ BHENERGY
          DOUBLE PRECISION DIST2, RMAT(3,3)

          FINISHED=.FALSE.
          ALLOCATE(EI,EF,Q(LNOPT),FIN(LNOPT))
          EI=EII;EF=EFF;Q=QQ;FIN=FINFIN;
          TSREDO(1:LNOPT)=-1.0D0 ! can be uninitialised otherwise

          MOREPRINTING=PTEST
          IF (MOREPRINTING) THEN
               CALL ALLKEYCONNECTPRINT
               PRINT*
          ENDIF
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
          IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
!
!  It is important to use EII, EFF, QQ and FINFIN here because EI, EF, Q and FIN are
!  declared as pointers. Trying to asign static data to the pointers cannot work!!
!
             OPEN(99,FILE='redopoints',STATUS='OLD')
             ALLOCATE(MIN1REDO(3*NA),MIN2REDO(3*NA))
             IF (AMHT) THEN
                NCOUNT=0
                DO J2=1,NMRES
                   IF (SEQ(J2).EQ.8) THEN
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                   ELSE
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                   ENDIF
                ENDDO
             ELSE
                READ(99,*) (QQ(J1),J1=1,3*NA)
             ENDIF
             NREDO=1
41           CONTINUE
             IF (AMHT) THEN
                NCOUNT=0
                DO J2=1,NMRES
                   IF (SEQ(J2).EQ.8) THEN
                      NCOUNT=NCOUNT+1
                      READ(99,*,END=42) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                   ELSE
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                   ENDIF
                ENDDO
             ELSE
                READ(99,*,END=42) (FINFIN(J1),J1=1,3*NA)
             ENDIF
             NREDO=NREDO+1
             GOTO 41
42           CONTINUE
             IF (RIGIDINIT) THEN
                CALL GENRIGID_POTENTIAL(QQ,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
             ELSE
                CALL POTENTIAL(QQ,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
             ENDIF
             IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
             IF (RIGIDINIT) THEN
                CALL GENRIGID_POTENTIAL(FINFIN,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
             ELSE
                CALL POTENTIAL(FINFIN,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
             ENDIF
             EI=EII;EF=EFF;Q=QQ;FIN=FINFIN;
             PRINT '(3(A,I6))',' newconnect> redopoints file contains ',(NREDO-1)/2,' transition states and ',(NREDO+1)/2,' minima'
             IF (MOD(NREDO,2).EQ.0) THEN
                PRINT '(A,I6,A)',' newconnect> WARNING *** # of coordinate sets read from redopoints is ', &
  &                          NREDO,' which is not an odd number!'
             ENDIF
 
             PRINT '(A,F20.10)',' newconnect> start and finish minima overwritten from redopoints file'
             WRITE(*,'(a,2(g20.10,a))') ' newconnect> Initial energy=',EI,' RMS force=',RMSINITIAL
             WRITE(*,'(a,2(g20.10,a))') ' newconnect> Final   energy=',EF,' RMS force=',RMSFINAL
             CALL ALIGN_DECIDE(QQ,FINFIN,NA,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ENDPOINTSEP,DIST2,RIGIDBODY,RMAT)
             WRITE(*,'(A,G20.10)') ' newconnect> separation=',ENDPOINTSEP
          ENDIF

          CALL INITIALISE(NA,EI,Q,EF,FIN,ENDPOINTSEP)

          IF (REDOPATHXYZ) PRINT '(A)',' newconnect> Redo run will use available path.<n>.xyz files'
          IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
             REWIND(99)
             IF (AMHT) THEN
                NCOUNT=0
                DO J2=1,NMRES
                   IF (SEQ(J2).EQ.8) THEN
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                   ELSE
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                   ENDIF
                ENDDO
             ELSE
                READ(99,*) (MIN2REDO(J1),J1=1,3*NA)
             ENDIF
          ENDIF

          DO NCONDONE=1,NCONMAX
               IF (QCICYCLEST) THEN
                  IF (QCICYCN.LT.NCONDONE) THEN
                     INTCONSTRAINTT=.FALSE.
                     INTLJT=.FALSE.
                     QCIDNEBT=.TRUE.
                     WRITE(*,'(A)') ' newconnect> turned QCI interpolation off'
                  ENDIF
               ENDIF
               IF (REDOPATH) THEN
                  IF (NCONDONE.EQ.(NREDO-1)/2+1) THEN
                     WRITE(*,'(A)') ' newconnect> no more coordinates in redopoints'
                     REDOPATH=.FALSE.
                     CLOSE(99)
                  ENDIF
               ENDIF
               WRITE(CHR,'(i5)') NConDone
               WRITE(STR,'(3a)') '>>>>>>>>>>>>>>>>>>>>> CONNECT CYCLE ',trim(adjustl(chr)),' >>>>>>>>>>>>>>>>>>>>>'
               WRITE(CHR,'(i5)') Nmin
               WRITE(STR,'(a)')  trim(str)//' '//trim(adjustl(chr))//' minima and'
               WRITE(CHR,'(i5)') NTSGOOD
               WRITE(STR,'(a)')  trim(str)//' '//trim(adjustl(chr))//' ts are known'
               WRITE(*,'(1x,a)') trim(str)//' '//repeat('>',107-len_trim(str))
!
!  Read in transition state coordinates from redopoints if available.
!  Must NOT use the pointers EI and EF here!!!
!
               IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                  IF (AMHT) THEN
                     MIN1REDO(1:3*NA)=MIN2REDO(1:3*NA)
                     NCOUNT=0
                     DO J2=1,NMRES
                        IF (SEQ(J2).EQ.8) THEN
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                        ELSE
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                        ENDIF
                     ENDDO
                     NCOUNT=0
                     DO J2=1,NMRES
                        IF (SEQ(J2).EQ.8) THEN
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                        ELSE
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                        ENDIF
                     ENDDO
                     IF (DEBUG) PRINT '(A,I6,A)',' newconnect> ',NCOUNT,' ts coordinates read from redopoints'
                     IF (RIGIDINIT) THEN
                        CALL GENRIGID_POTENTIAL(MIN1REDO,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
                        CALL GENRIGID_POTENTIAL(TSREDO,ETS,LGDUMMY,.TRUE.,.FALSE.,RMSTS,.FALSE.,.FALSE.)
                        CALL GENRIGID_POTENTIAL(MIN2REDO,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
                     ELSE
                        CALL POTENTIAL(MIN1REDO,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
                        CALL POTENTIAL(TSREDO,ETS,LGDUMMY,.TRUE.,.FALSE.,RMSTS,.FALSE.,.FALSE.)
                        CALL POTENTIAL(MIN2REDO,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
                     ENDIF
                     WRITE(*,'(a,2(g20.10,a))') ' newconnect> initial energy=',EII,' RMS force=',RMSINITIAL
                     WRITE(*,'(a,2(g20.10,a))') ' newconnect> ts energy     =',ETS,' RMS force=',RMSTS
                     WRITE(*,'(a,2(g20.10,a))') ' newconnect> final energy  =',EFF,' RMS force=',RMSFINAL
                     REDOE1=EII
                     REDOE2=EFF
                  ELSE
                     MIN1REDO(1:3*NA)=MIN2REDO(1:3*NA)
                     READ(99,*) (TSREDO(J1),J1=1,3*NA)
                     READ(99,*) (MIN2REDO(J1),J1=1,3*NA)
                     IF (RIGIDINIT) THEN
                        CALL GENRIGID_POTENTIAL(MIN1REDO,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
                        CALL GENRIGID_POTENTIAL(TSREDO,ETS,LGDUMMY,.TRUE.,.FALSE.,RMSTS,.FALSE.,.FALSE.)
                        IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
                        CALL GENRIGID_POTENTIAL(MIN2REDO,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
                     ELSE
                        CALL POTENTIAL(MIN1REDO,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
                        CALL POTENTIAL(TSREDO,ETS,LGDUMMY,.TRUE.,.FALSE.,RMSTS,.FALSE.,.FALSE.)
                        IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
                        CALL POTENTIAL(MIN2REDO,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
                     ENDIF
                     WRITE(*,'(a,2(g20.10,a))') ' newconnect> initial energy=',EII,' RMS force=',RMSINITIAL
                     WRITE(*,'(a,2(g20.10,a))') ' newconnect> ts energy     =',ETS,' RMS force=',RMSTS
                     WRITE(*,'(a,2(g20.10,a))') ' newconnect> final energy  =',EFF,' RMS force=',RMSFINAL
                     REDOE1=EII
                     REDOE2=EFF
                  ENDIF
                  WRITE(*,'(A)') ' newconnect> minimum and transition state coordinates read from file redopoints'
               ENDIF

               IF (REDOPATH) THEN
                  CALL TRYCONNECT(1,2,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
               ELSE
                  CALL DECIDE

!    call tryconnect for each pair of minima specified in the shortest path returned by DIJKSTRA
!    Call BHINTERP to fill in addtional minima if required.

                  IF (BHINTERPT) THEN
                     NSTART=1
                     CHANGED=.FALSE.
34                   CONTINUE
                     NOPRINT=.FALSE.
                     DO J1=NSTART,NDIJPAIRS
                        IF(.NOT.NOPRINT) PRINT '(A)',' newconnect> DIJPAIR distances and minima:'
                        IF(.NOT.NOPRINT) PRINT '(I8,G20.10,I8,G20.10,I8,G20.10)',(J2,DIJPAIRDIST(J2),DIJPAIR(J2,1), &
                                         MI(DIJPAIR(J2,1))%DATA%E,DIJPAIR(J2,2),MI(DIJPAIR(J2,2))%DATA%E,J2=1,NDIJPAIRS)
                        DO J2=1,NDIJPAIRS
                           IF (DIJPAIRDIST(J2).LT.GEOMDIFFTOL) THEN
!                             PRINT '(A,I8,2G20.10)','J2,DIJPAIRDIST(J2),GEOMDIFFTOL=',J2,DIJPAIRDIST(J2),GEOMDIFFTOL
                              PRINT '(A,G20.10,A,G20.10)',' newconnect> WARNING - distance ',DIJPAIRDIST(J2), &
  &                                       ' is less than GEOMDIFFTOL=',GEOMDIFFTOL 
!                             STOP
                           ENDIF
                        ENDDO
!
!  Only try interpolation at the first connection attempt for each pair.
!
                        IF ((BHDEBUG).AND.(.NOT.NOPRINT)) PRINT '(A,I8,A,I8,A,I8)', &
  &                      ' bhinterp> tries for DIJPAIR min ',DIJPAIR(J1,1),' and min ',DIJPAIR(J1,2),' = ', &
  &                        MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2)))
                        IF ((MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2))).EQ.0).AND. &
                            (DIJPAIRDIST(J1).GT.BHDISTTHRESH)) THEN
                           PRINT '(A,I8,A,I8,A,I8,A,G20.10)',' newconnect> gap ',J1,' trying BH interpolation for minima   ', &
  &                                DIJPAIR(J1,1),' and ',DIJPAIR(J1,2),' dist=',DIJPAIRDIST(J1)
!
!  Count this as a connection attempt. Otherwise we may try to interpolate twice.
!
                           MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2)))= &
  &                          MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2)))+1
                           CSTART(1:3*NA)=MI(DIJPAIR(J1,1))%DATA%X(1:3*NA)
                           CFINISH(1:3*NA)=MI(DIJPAIR(J1,2))%DATA%X(1:3*NA)
                           ESTART=MI(DIJPAIR(J1,1))%DATA%E
                           EFINISH=MI(DIJPAIR(J1,2))%DATA%E
                           CALL BHINTERP(CSTART,CFINISH,3*NA,NA,INTERPCOORDS,SUCCESS,DSTART,DFINISH,ENERGY,ESTART,EFINISH, &
  &                                      DIJPAIRDIST(J1))
                           BIGGERGAP=.FALSE.
! bs360                    IF ((DSTART.LT.(1.5D0*DIJPAIRDIST(J1))).AND.(DFINISH.LT.(1.5D0*DIJPAIRDIST(J1)))) THEN
                           IF ((DSTART.LT.DIJPAIRDIST(J1)).OR.(DFINISH.LT.DIJPAIRDIST(J1))) THEN
                              IF (SUCCESS.AND.BHDEBUG) PRINT '(A,2G20.10,A)', &
  &                              ' newconnect> at least BH one interpolated distance is shorter: ',DSTART,DFINISH,' accept'
!
! If we are going to allow this, then we skip bisection of either new gap.
!
                               IF ((DSTART.GT.DIJPAIRDIST(J1)).OR.(DFINISH.GT.DIJPAIRDIST(J1))) BIGGERGAP=.TRUE.
                           ELSE IF ((DSTART.GT.DIJPAIRDIST(J1)).AND.(DFINISH.GT.DIJPAIRDIST(J1))) THEN
                              IF (SUCCESS.AND.BHDEBUG) PRINT '(A,2G20.10,A)', &
  &                              ' newconnect> BH both interpolated distances are larger: ',DSTART,DFINISH,' reject'
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
!                                DO J2=1,NMIN-1
!                                   PRINT '(A,2I8,G20.10)','min1 min2 d=',NMIN,J2,MI(NMIN)%DATA%D(J2)
!                                ENDDO
                              ELSE
                                 WRITE(*,'(A,I7)') ' newconnect> Interpolated minimum is old minimum ',POSITION
                                  
                                 IF (DIJPAIR(J1,1).GE.POSITION) THEN
                                    IF (MI(DIJPAIR(J1,1))%DATA%NTRIES(POSITION).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' bhinterp> maximum connection attempts for ', &
  &                                      DIJPAIR(J1,1),' and ',POSITION,' is ', &
  &                                      MI(DIJPAIR(J1,1))%DATA%NTRIES(POSITION), &
  &                                      ' skip; distance=',MI(DIJPAIR(J1,1))%DATA%D(POSITION)
                                       GOTO 2
                                    ENDIF
                                 ELSE
                                    IF (MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,1)).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' newconnect> maximum connection attempts for ', &
  &                                      DIJPAIR(J1,1),' and ',POSITION,' is ', &
  &                                      MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,1)), &
  &                                      ' skip; distance=',MI(POSITION)%DATA%D(DIJPAIR(J1,1))
                                       GOTO 2
                                    ENDIF
                                 ENDIF
                                 IF (DIJPAIR(J1,2).GE.POSITION) THEN
                                    IF (MI(DIJPAIR(J1,2))%DATA%NTRIES(POSITION).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' newconnect> maximum connection attempts for ', &
  &                                      DIJPAIR(J1,2),' and ',POSITION,' is ', &
  &                                      MI(DIJPAIR(J1,2))%DATA%NTRIES(POSITION), &
  &                                      ' skip; distance=',MI(DIJPAIR(J1,2))%DATA%D(POSITION)
                                       GOTO 2
                                    ENDIF
                                 ELSE
                                    IF (MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,2)).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' newconnect> maximum connection attempts for ', &
  &                                      DIJPAIR(J1,2),' and ',POSITION,' is ', &
  &                                      MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,2)), &
  &                                      ' skip; distance=',MI(POSITION)%DATA%D(DIJPAIR(J1,2))
                                       GOTO 2
                                    ENDIF
                                 ENDIF
!
!  Reject interpolated minimum if it is already in the list, to avoid an
!  infinite loop.
!
                                 DO J3=1,NDIJPAIRS
                                    IF ((POSITION.EQ.DIJPAIR(J3,1)).OR.(POSITION.EQ.DIJPAIR(J3,2))) THEN
                                       PRINT '(A)',' newconnect> interpolated minimum is already in the list - reject'
                                       GOTO 2
                                    ENDIF
                                 ENDDO
                              ENDIF

                              CHANGED=.TRUE.

                              NULLIFY(PINTERPCOORDS,PENERGY) 
!
!  If the interpolated minmum is old there are four posibilities. It could already
!  be connected to both endpoints, to one of them, or to neither. The endpoints should
!  not have been chosen in the first place if a connection exists via the interpolated old
!  minimum, so that leaves three possibilities.
!
                              IF (.NOT.MINNEW) THEN
                                 IF (POSITION.GT.DIJPAIR(J1,1)) THEN
                                    OLDDISTS=MI(POSITION)%DATA%D(DIJPAIR(J1,1))
                                 ELSE
                                    OLDDISTS=MI(DIJPAIR(J1,1))%DATA%D(POSITION)
                                 ENDIF
                                 IF (POSITION.GT.DIJPAIR(J1,2)) THEN
                                    OLDDISTF=MI(POSITION)%DATA%D(DIJPAIR(J1,2))
                                 ELSE
                                    OLDDISTF=MI(DIJPAIR(J1,2))%DATA%D(POSITION)
                                 ENDIF
                                 PRINT '(2(A,I8),2(A,G20.10))','newconnect> old distance between ',POSITION,' and ',DIJPAIR(J1,1), &
  &                                                       ' is ',OLDDISTS,' new dist=',DSTART
                                 PRINT '(2(A,I8),2(A,G20.10))','newconnect> old distance between ',POSITION,' and ',DIJPAIR(J1,2), &
  &                                                       ' is ',OLDDISTF,' new dist=',DFINISH
                                 DSTART=OLDDISTS
                                 DFINISH=OLDDISTF
                              ENDIF

                              IF (DSTART.LT.GEOMDIFFTOL) THEN ! replace entry by old minimum/second minimum pair
                                 DIJPAIR(J1,1)=POSITION
                                 DIJPAIRDIST(J1)=DFINISH
                                 NSTART=J1
                                 IF (BIGGERGAP) NSTART=J1+1
                                 GOTO 34 ! go back
                              ELSE IF (DFINISH.LT.GEOMDIFFTOL) THEN ! replace entry by first minimum/old minimum pair
                                 DIJPAIR(J1,2)=POSITION
                                 DIJPAIRDIST(J1)=DSTART
                                 NSTART=J1
                                 IF (BIGGERGAP) NSTART=J1+1
                                 GOTO 34 ! go back
                              ELSE ! add an entry to the list - neither minimum connected to the minimum in question
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
                                 DIJPAIR(J1+1,1)=POSITION
                                 DIJPAIR(J1+1,2)=TEMPDIJPAIR(J1,2)
                                 DIJPAIRDIST(J1)=DSTART
                                 DIJPAIRDIST(J1+1)=DFINISH
                                 DO J2=J1+1,NDIJPAIRS-1
                                    DIJPAIR(J2+1,1:2)=TEMPDIJPAIR(J2,1:2)
                                    DIJPAIRDIST(J2+1)=TEMPDIJPAIRDIST(J2)
                                 ENDDO
                                 NSTART=J1
                                 IF (BIGGERGAP) NSTART=J1+2
                                 DEALLOCATE(TEMPDIJPAIR,TEMPDIJPAIRDIST)
                                 GOTO 34 ! go back, since NDIJPAIRS has increased by one
                              ENDIF
                           ELSE
                              NSTART=J1+1 ! to avoid an infinite loop
                           ENDIF
                        ELSE
                           NOPRINT=.TRUE.
                        ENDIF
2                       CONTINUE
                     ENDDO
                     IF (CHANGED) THEN
                        WRITE(*,'(A)') ' newconnect> The unconnected minima in the chain and their distances are now:'
                        DO J1=1,NDIJPAIRS
                           WRITE(*,'(I6,F12.2,I6)',ADVANCE="NO") DIJPAIR(J1,1),DIJPAIRDIST(J1),DIJPAIR(J1,2)
                        ENDDO
                     ENDIF
                  ENDIF
                  WRITE(*,'(A)') ' ' ! to advance to the next line

                  IF (.NOT.NEWCONNECTT) THEN ! this must be a BHINTERP only run
                     PRINT '(A,I8)',' newconnect> nmin=',nmin
                     IF (DUMPDATAT) THEN ! call geopt to generate min.data.info entries
                        FNAMEF='points.final.bh'
                        EFNAME='  '
                        DO J1=3,NMIN ! don't include the end point minima
!                          write(*,*) 'in newconnect.f90, bhenergy=',bhenergy
                           CALL FLUSH(6)
                           REOPTIMISEENDPOINTS=.FALSE.
                           BHENERGY=MI(J1)%DATA%E
                           CALL GEOPT(FNAMEF,EFNAME,MI(J1)%DATA%X(1:3*NA)) ! file names are dummies here
                           CALL FLUSH(881)    ! flush min.data.info unit again
                        ENDDO
                     ENDIF
                     STOP
                  ENDIF
!
! If intconstraint is specified use the intlbfgs routine to generate a series
! of images, minimise them, and use the minima in order for separate DNEB runs with the
! true potential.
! Or, use the interpolation routine to generate transition states directly.
!
                  DO J1=1,NDIJPAIRS
                     NMINSAVE=NMIN
                     JS=MAX(DIJPAIR(J1,1),DIJPAIR(J1,2))
                     JF=MIN(DIJPAIR(J1,1),DIJPAIR(J1,2))
!
! Check that we haven't already connected this pair in a previous search in this cycle.
!
                     IF (DIJPAIRDIST(J1).EQ.0.0D0) THEN
                        PRINT '(A,I8,A,2I8,A,G20.10)',' newconnect> skipping DIJKSTRA pair number ',J1,' minima ', &
  &                            JS,JF,' distance in now',DIJPAIRDIST(J1)
                        CYCLE
                     ENDIF
                     IF (MI(JS)%DATA%D(JF).EQ.0.0D0) THEN
                        PRINT '(A,I8,A,2I8,A,G20.10)',' newconnect> skipping DIJKSTRA pair number ',J1,' minima ', &
  &                            JS,JF,' distance is now ',MI(JS)%DATA%D(JF)
                        CYCLE
                     ENDIF
!                    WRITE(*,'(A,I5,A,2I5)') ' newconnect> trying DIJKSTRA pair number ',J1,' minima ',js,jf
!
! Do count constrained potential interpolations as connection attempts.
!
                     IF (INTCONSTRAINTT.OR.INTLJT) THEN
                        IF (QCICYCLEST.AND.(DIJPAIRDIST(J1).LT.QCICYCDIST)) THEN
                           MI(JS)%DATA%NTRIES(JF)=MI(JS)%DATA%NTRIES(JF)+1
                           WRITE(*,'(A)') ' newconnect> distance below threshold, use DNEB'
                           QCIDNEBT=.TRUE.
                        ELSE  
                           MI(JS)%DATA%INTNTRIES(JF)=MI(JS)%DATA%INTNTRIES(JF)+1
                           QCIDNEBT=.FALSE.
                        ENDIF
                     ELSE
                        MI(JS)%DATA%NTRIES(JF)=MI(JS)%DATA%NTRIES(JF)+1
                     ENDIF
!
! Minimum JF becomes image number one and minimum JS becomes image NIMAGE+2
! Just to be confusing, tryconnect has the arguments as JS,JF, so the minima
! get swapped around in making the DNEB interpolation. This shouldn't matter,
! except that the reconnection procedure for INTCONSTRAINTT/INTLJT below needs to
! know this order!
!
                     IF ((INTCONSTRAINTT.OR.INTLJT).AND.(QCICYCLEST.AND.(DIJPAIRDIST(J1).LT.QCICYCDIST))) THEN
                        CALL TRYCONNECT(JF,JS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
                     ELSE
                        CALL TRYCONNECT(JF,JS,TSREDO,REDOPATH,REDOPATHXYZ,INTCONSTRAINTT,INTLJT)
                     ENDIF
!                    IF (INTCONSTRAINTT.OR.INTLJT) THEN
!                       IF ((MI(JS)%DATA%INTNTRIES(JF).GE.INTNTRIESMAX).AND.(MI(JS)%DATA%D(JF).GT.1.0D-5)) THEN
!                          MI(JS)%DATA%D(JF)=HUGE(1.0D0)
!                          MI(JS)%DATA%INTERP(JF)=HUGE(1.0D0)
!                       ENDIF
!                    ELSE
                        IF ((MI(JS)%DATA%NTRIES(JF).GE.NTRIESMAX).AND.(MI(JS)%DATA%D(JF).GT.1.0D-5)) THEN
                           MI(JS)%DATA%D(JF)=HUGE(1.0D0)
                           IF (INTERPCOSTFUNCTION) MI(JS)%DATA%INTERP(JF)=HUGE(1.0D0)
                        ENDIF
!                    ENDIF

!                    IF (DEBUG) PRINT '(A,I6,A,I6)',' newconnect> After TRYCONNECT for minima ',JS,' and ',JF
!                    IF (DEBUG) PRINT '(A,2I6)',' newconnect> INTNTRIES,NTRIES=',MI(JS)%DATA%INTNTRIES(JF),MI(JS)%DATA%NTRIES(JF)
!                    IF (DEBUG) PRINT '(A,G20.10)',' newconnect> D=',MI(JS)%DATA%D(JF)
!                    IF (DEBUG.AND.INTERPCOSTFUNCTION) PRINT '(A,G20.10)',' newconnect> INTERP=',MI(JS)%DATA%INTERP(JF)
!
! Reweight distance metric between minima JS, JF and all the new minima found in the
! most recent connection attempt. This is to encourage trying new connection attempts
! within this set of minima, even if the distances would otherwise cause some other pairs
! to be tried.
!
                     ! PRINT '(A,2I8)',' newconnect> NMINSAVE,NMIN=',NMINSAVE,NMIN
                     IF (DIJKSTRALOCAL.NE.1.0D0) THEN
                        DO J2=NMINSAVE+1,NMIN
                           MI(J2)%DATA%D(JS)=MI(J2)%DATA%D(JS)*DIJKSTRALOCAL
                           MI(J2)%DATA%D(JF)=MI(J2)%DATA%D(JF)*DIJKSTRALOCAL
                           PRINT '(A,3I8,2G20.10)','newconnect> J2,JS,JF=',J2,JS,JF
                           PRINT '(A,3I8,2G20.10)','newconnect> J2,JS,JF,d(JS),d(JF)=',J2,JS,JF,MI(J2)%DATA%D(JS),MI(J2)%DATA%D(JF)
                           DO J3=J2+1,NMIN
                              MI(J3)%DATA%D(J2)=MI(J3)%DATA%D(J2)*DIJKSTRALOCAL
                              PRINT '(A,2I8,G20.10)','newconnect> J3,J2=',J2,J3,MI(J3)%DATA%D(J2)
                           ENDDO
                        ENDDO
                     ENDIF
!                    PRINT '(A)','newconnect distances summary:'
!                    DO J2=1,NMIN
!                       DO J3=J2+1,NMIN
!                          PRINT '(A,2I8,G20.10)','newconnect> J3,J2=',J2,J3,MI(J3)%DATA%D(J2)
!                       ENDDO
!                    ENDDO
!
! Do proper DNEB searches for new minima found in the previous run in order.
!
!                      IF ((INTCONSTRAINTT.OR.INTLJT).AND.(.NOT.INTTST)) THEN
!                         FCD=.FALSE. ! turn off special first connection image specification if set
!                         NMINSAVE2=NMIN ! because NMIN could increase in the next tryconnect calls! 
!                         IF (NMINSAVE2-NMINSAVE.GT.0) THEN
!                            MYJS=MAX(JF,NMINSAVE+1)
!                            MYJF=MIN(JF,NMINSAVE+1)
! !
! ! Don't try the pair if they have already been connected by a previous run!
! !
!                            MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF)+1
!                            IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
!                               WRITE(*,'(A,2I6)') ' newconnect> trying minima from interpolation run: ',MYJS,MYJF
!                               CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
!                               IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
!                                   (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF)=HUGE(MI(MYJS)%DATA%D(MYJF))
!                            ENDIF
!                            DO J2=1,NMINSAVE2-NMINSAVE-1
!                               MYJF=NMINSAVE+J2
!                               MYJS=NMINSAVE+J2+1
!                               IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
!                                  WRITE(*,'(A,2I6)') ' newconnect> trying minima from interpolation run: ',MYJS,MYJF
!                                  MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF)+1
! !
! !    Set edge weight to infinity if we have reached the maximum number of tries for this pair.
! !
!                                 CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
!                                 IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
!                                     (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF) = HUGE(MI(MYJS)%DATA%D(MYJF))
!                               ENDIF
!                            ENDDO
!                            MYJS=NMINSAVE2
!                            MYJF=JS
!                            IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
!                               WRITE(*,'(A,2I6)') ' newconnect> trying minima from interpolation run: ',MYJS,MYJF
!                               MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF)+1
!                               CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
!                               IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
!                                   (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF) = HUGE(MI(MYJS)%DATA%D(MYJF))
!                            ENDIF
!                         ELSE
!                            MYJS=JS
!                            MYJF=JF
!                            IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
!                               WRITE(*,'(A,2I6)') ' newconnect> trying minima from constrained potential run: ',MYJS,MYJF
!                               MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF)+1
!                               CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
!                               IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
!                                   (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF) = HUGE(MI(MYJS)%DATA%D(MYJF))
!                            ENDIF
!                         ENDIF
!                         FCD=.TRUE.
!                      ENDIF
                  ENDDO
                  DEALLOCATE(DIJPAIR,DIJPAIRDIST)
               ENDIF
               IF (FINISHED) EXIT
          ENDDO

          CALL OUTPUT
          CALL DEINITIALISE
          DEALLOCATE(Q,FIN,EI,EF)
          IF (ALLOCATED(MIN1REDO)) DEALLOCATE(MIN1REDO)
          IF (ALLOCATED(MIN2REDO)) DEALLOCATE(MIN2REDO)
     END SUBROUTINE NEWCONNECT

END MODULE NEWCONNECTMODULE
