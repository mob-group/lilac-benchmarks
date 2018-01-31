C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE CONNECT(NCDONE,Q)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODNEB
      USE modcharmm
      USE MODUNRES
      IMPLICIT NONE
      INTEGER NMAXMIN, NMAXTS, MSTART, MFINISH, MSTEP
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      INTEGER J1, J2, NMIN, JDOING, NIMAGESAVE, WHICHTS(NMAXMIN), JUSE, J3,
     1        NTSUSED, NMINCON, IGUESS
      LOGICAL MFLAG, CON(NMAXMIN), TSUSED(NMAXTS), USEOLD, PTEST, REVERSET, YESNO
      DOUBLE PRECISION ENERGY, VNEW(3*NATOMS), VECS(3*NATOMS), RMS, EVALMIN, QSTART(3*NATOMS), CMX1, CMY1, CMZ1,
     1 QPLUS(3*NATOMS), DISTPF, EFIN, ESTART, DUMMY(3*NATOMS),OVEC(3), CMX2, CMY2, CMZ2, 
     2 RMSFIN, RMSSTART, ETS, TSEN(NMAXTS), PLUSEN(NMAXTS), H1VEC(3),H2VEC(3), Q(3*NATOMS),
     3 QSAVETS(3*NATOMS,NMAXTS), MINUSEN(NMAXTS), QSAVE(3*NATOMS,NMAXMIN), EMIN(NMAXMIN),
     4 PATHLENGTH(NMAXMIN),DISP(NMAXMIN),GAMMA(NMAXMIN),NTILDE(NMAXMIN), QMINUS(3*NATOMS),
     5 DIAG(3*NATOMS), FRQSTS(3*NATOMS), FRQSPLUS(3*NATOMS), FRQSMINUS(3*NATOMS),
     6 FSAVETS(3*NATOMS,NMAXTS), FSAVEMIN(3*NATOMS,NMAXMIN), MXSTPSAVE, Q1(3*NATOMS), SEPARATION, RMAT(3,3)
      INTEGER NCDONE, NMOVE, NTS, NTSREDO
      CHARACTER(LEN=80) ITSTRING
      CHARACTER(LEN=80) EOFSSTRING
      CHARACTER(LEN=80) DSTRING
      COMMON /HS/ NMOVE
      COMMON /NEBRMS/ RMS,ESTART,EFIN,SEPARATION
      LOGICAL CONNECTT, DUMPPATH, TSFRQDONE, READPATH, CALCRATES, STOPFIRST, GT10, MINFRQDONE
      DOUBLE PRECISION TEMPERATURE, HRED
      INTEGER NCONNECT
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
C     COMMON /CONREAL/ QSAVETS, QSAVE, FSAVETS, FSAVEMIN
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
      DOUBLE PRECISION STOPDISP
      LOGICAL STOPDISPT
      COMMON /STOPD/ STOPDISP, STOPDISPT
      DOUBLE PRECISION ALLANG ! JMC
      LOGICAL PATHFAILT, ORDERSTOP ! jmc

      ORDERSTOP=.FALSE.
      IF (UNRST) ALLOCATE(MYQMINSAVE(3*NATOMS,50))

      INQUIRE(FILE='redopoints',EXIST=YESNO)
      IF (YESNO) THEN
         IF (REDOPATH) THEN
            PRINT '(A)','connect> Transition state coordinates will be read from file redopoints'
!        ELSE
!           PRINT '(A)','connect> WARNING - redopoints file present, but no REDOPATH keyword'
         ENDIF
      ELSE
         IF (REDOPATH) THEN
            PRINT '(A)','connect> WARNING - REDOPATH keyword was specified, but no redopoints file'
            REDOPATH=.FALSE.
         ENDIF
      ENDIF

!
!  Changed so that only transition state coordinates are read in.
!  The points in finish must therefore be correct.
!
      NTSREDO=1
      IF (REDOPATH) THEN
         OPEN(UNIT=99,FILE='redopoints',STATUS='OLD')
33       READ(99,*) (DUMMY(J1),J1=1,NOPT)    ! minimum
         READ(99,*,END=32) (QSAVETS(J1,NTSREDO),J1=1,NOPT) ! ts
         NTSREDO=NTSREDO+1
         GOTO 33
32       WRITE(*,'(I5,A)') NTSREDO,' transition state coordinates read from file redopoints'
         CLOSE(99)
      ENDIF
C
C IGUESS is used if GUESSTST is set (for CHARMM and UNRES)
C
      IGUESS=1
      IF (GUESSTST) TRYNEB=.FALSE. ! declared in modcharmm

      MXSTPSAVE=MXSTP
      NIMAGESAVE=NIMAGE
      TSFRQDONE=.TRUE.
      IF (BFGSTST) TSFRQDONE=.FALSE.
      MINFRQDONE=.TRUE.
      IF (BFGSMINT.OR.RKMIN.OR.BSMIN) MINFRQDONE=.FALSE.
      DO J1=2,NATOMS
         IF (ATMASS(J1).NE.ATMASS(J1-1)) THEN
            TSFRQDONE=.FALSE.
            MINFRQDONE=.FALSE.
            GOTO 31
         ENDIF 
      ENDDO
31    PTEST=.FALSE.
C
C  Even if a second derivative method was used the frequencies won;t have been
C  done for the first and last minima, and the zeros may be shifted. Let;s just
C  set the FRQDONE variables to false!
C
      TSFRQDONE=.FALSE.
      MINFRQDONE=.FALSE.

      IF (DEBUG) PTEST=.TRUE.
      NMIN=2
      DO J1=1,NOPT
         QSAVE(J1,1)=Q(J1)
         QSAVE(J1,2)=FIN(J1)
      ENDDO

      IF (UNRST) THEN
         DO J1=1,nres
            c(1,J1)=FIN(6*(J1-1)+1)
            c(2,J1)=FIN(6*(J1-1)+2)
            c(3,J1)=FIN(6*(J1-1)+3)
            c(1,J1+nres)=FIN(6*(J1-1)+4)
            c(2,J1+nres)=FIN(6*(J1-1)+5)
            c(3,J1+nres)=FIN(6*(J1-1)+6)
         ENDDO
         CALL UPDATEDC
         CALL int_from_cart(.true.,.false.)
         CALL chainbuild
      ENDIF

C jmc note can't use both calcdihe and stopdisp keywords, as using the refppsangle array here for finish's internal coordinates
C COULD CAUSE PROBLEMS FOR PATHSAMPLE....
C NEED TO RETURN TO THIS ISSUE!!!!!!!!!!!
C Also note that refcoord and refppsangle don't correspond to the same structure!
      IF (STOPDISPT.AND.UNRST) THEN
         DO J1=1,nres-3
            UREFPPSANGLE(J1)=phi(J1+3)
            UREFPPSANGLE(J1+nres-3)=omeg(J1+1)
         END DO
         UREFPPSANGLE(nres-2+nres-3)=omeg(nres-2+1)
      END IF

      CALL POTENTIAL(FIN,EFIN,  VNEW,.FALSE.,.FALSE.,RMSFIN,.FALSE.,.FALSE.)

      IF (UNRST) THEN
         DO J1=1,nres
            c(1,J1)=Q(6*(J1-1)+1)
            c(2,J1)=Q(6*(J1-1)+2)
            c(3,J1)=Q(6*(J1-1)+3)
            c(1,J1+nres)=Q(6*(J1-1)+4)
            c(2,J1+nres)=Q(6*(J1-1)+5)
            c(3,J1+nres)=Q(6*(J1-1)+6)
         ENDDO
         CALL UPDATEDC
         CALL int_from_cart(.true.,.false.)
         CALL chainbuild
      ENDIF
      CALL POTENTIAL(Q,  ESTART,VNEW,.FALSE.,.FALSE.,RMSSTART,.FALSE.,.FALSE.)
      EMIN(1)=ESTART
      EMIN(2)=EFIN
      DO J1=1,NMAXMIN
         WHICHTS(J1)=-1
         CON(J1)=.FALSE.
      ENDDO
      NCDONE=0
      NTS=0
      REVERSET=.FALSE.
      IF (STOPDISP.LT.0.0D0) REVERSET=.TRUE.
C
C---------------------------------------start of main loop---------------------------------------------------
C
10    NCDONE=NCDONE+1
C
C  Print current state, check for completion and choose next pair of minima to connect, moving their
C  coordinates into Q and FIN and energies into ESTART and EFIN
C
      IF (REDOPATH.AND.(NTS.EQ.NTSREDO)) REDOPATH=.FALSE.
      CALL STARTUP(MFLAG,NMIN,NTS,CON,EMIN,WHICHTS,TSEN,NTSUSED,NMINCON,NCDONE,NCONNECT,MSTART,MFINISH,MSTEP,
     1                   REVERSET,STOPFIRST,NOPT,Q,QSAVE,QSTART,FIN,JDOING,ESTART,EFIN,TSUSED)
      IF (MFLAG) THEN
         CALL CLEANUP(WHICHTS,TSEN,PLUSEN,MINUSEN,PATHLENGTH,DISP,GAMMA,NTILDE,EMIN,NMIN,DUMPPATH,QSAVE,QSAVETS,
     1                FSAVETS,FSAVEMIN,MINFRQDONE,TSFRQDONE,BULKT,FILTH,FILTHSTR,TWOD,UNRST)
         IF (UNRST) DEALLOCATE(MYQMINSAVE)
         RETURN
      ENDIF
      IF (STOPDISPT) THEN
         IF (.NOT.TAGT) THEN
            PRINT '(A)','connect> ERROR - STOPDISPT is true but TAGT is false'
            STOP
         ENDIF
         IF (STOPDISP.GT.0.0D0) THEN
            DO J1=2,NMIN
               IF (CON(J1-1)) THEN
                  IF (UNRST) THEN
                     DO J2=1,3*NATOMS
                        Q1(J2)=QSAVE(J2,J1)
                     ENDDO
                     CALL UNRESCALCDIHESEC(DISTPF,ALLANG,Q1,ORDERSTOP)
                  ELSE
                     CMX1=0.0D0
                     CMY1=0.0D0
                     CMZ1=0.0D0
                     CMX2=0.0D0
                     CMY2=0.0D0
                     CMZ2=0.0D0
                     DO J2=1,NATOMS
                        CMX1=CMX1+QSAVE(3*(J2-1)+1,1)
                        CMY1=CMY1+QSAVE(3*(J2-1)+2,1)
                        CMZ1=CMZ1+QSAVE(3*(J2-1)+3,1)
                        CMX2=CMX2+QSAVE(3*(J2-1)+1,J1)
                        CMY2=CMY2+QSAVE(3*(J2-1)+2,J1)
                        CMZ2=CMZ2+QSAVE(3*(J2-1)+3,J1)
                     ENDDO
                     CMX1=CMX1/NATOMS
                     CMY1=CMY1/NATOMS
                     CMZ1=CMZ1/NATOMS
                     CMX2=CMX2/NATOMS
                     CMY2=CMY2/NATOMS
                     CMZ2=CMZ2/NATOMS
                     DISTPF=SQRT((QSAVE(3*(TAGNUM(1)-1)+1,1)-CMX1+CMX2-QSAVE(3*(TAGNUM(1)-1)+1,J1))**2+
     1                           (QSAVE(3*(TAGNUM(1)-1)+2,1)-CMY1+CMY2-QSAVE(3*(TAGNUM(1)-1)+2,J1))**2+
     2                           (QSAVE(3*(TAGNUM(1)-1)+3,1)-CMZ1+CMZ2-QSAVE(3*(TAGNUM(1)-1)+3,J1))**2)
C                    WRITE(*,'(A,I5,A,F12.5)') ' displacement between minimum 1 and minimum ',J1,' is ',DISTPF
                  ENDIF

                  IF (((DISTPF.GT.STOPDISP).AND.(.NOT.UNRST)).OR.ORDERSTOP) THEN
                     NMIN=J1
                     NTS=NMIN-1
                     MFLAG=.TRUE.
                     IF (ORDERSTOP) THEN
                         WRITE(*,'(A,I5,A,F12.5,A)') ' order parameter condition satisfied for finish and minimum ',J1,
     1                                           ' - stopping'
                     ELSE
                         WRITE(*,'(A,I5,A,F12.5,A)') ' order parameter condition satisfied for minimum 1 and minimum ',J1,
     1                                           ' - stopping'
                     ENDIF
                     CALL CLEANUP(WHICHTS,TSEN,PLUSEN,MINUSEN,PATHLENGTH,DISP,GAMMA,NTILDE,EMIN,NMIN,DUMPPATH,QSAVE,QSAVETS,
     1                            FSAVETS,FSAVEMIN,MINFRQDONE,TSFRQDONE,BULKT,FILTH,FILTHSTR,TWOD,UNRST)
                     IF (UNRST) DEALLOCATE(MYQMINSAVE)
                     RETURN
                  ELSE
                     GOTO 222  ! don't want to exit do loop until have tested all connected minima
                  ENDIF
               ELSE
                  GOTO 111
               ENDIF
222            CONTINUE
            ENDDO
         ELSE
            DO J1=NMIN-1,1,-1
               IF (CON(J1)) THEN
C                 DO J2=1,3*NATOMS
C                    Q1(J2)=QSAVE(J2,NMIN)
C                    Q2(J2)=QSAVE(J2,J1+1)
C                 ENDDO
                  CMX1=0.0D0
                  CMY1=0.0D0
                  CMZ1=0.0D0
                  CMX2=0.0D0
                  CMY2=0.0D0
                  CMZ2=0.0D0
                  DO J2=1,NATOMS
                     CMX1=CMX1+QSAVE(3*(J2-1)+1,NMIN)
                     CMY1=CMY1+QSAVE(3*(J2-1)+2,NMIN)
                     CMZ1=CMZ1+QSAVE(3*(J2-1)+3,NMIN)
                     CMX2=CMX2+QSAVE(3*(J2-1)+1,J1)
                     CMY2=CMY2+QSAVE(3*(J2-1)+2,J1)
                     CMZ2=CMZ2+QSAVE(3*(J2-1)+3,J1)
                  ENDDO
                  CMX1=CMX1/NATOMS
                  CMY1=CMY1/NATOMS
                  CMZ1=CMZ1/NATOMS
                  CMX2=CMX2/NATOMS
                  CMY2=CMY2/NATOMS
                  CMZ2=CMZ2/NATOMS
                  DISTPF=SQRT((QSAVE(3*(TAGNUM(1)-1)+1,NMIN)-CMX1+CMX2-QSAVE(3*(TAGNUM(1)-1)+1,J1))**2+
     1                        (QSAVE(3*(TAGNUM(1)-1)+2,NMIN)-CMY1+CMY2-QSAVE(3*(TAGNUM(1)-1)+2,J1))**2+
     2                        (QSAVE(3*(TAGNUM(1)-1)+3,NMIN)-CMZ1+CMZ2-QSAVE(3*(TAGNUM(1)-1)+3,J1))**2)
C                 WRITE(*,'(A,I5,A,I5,A,F12.5)') ' displacement between minimum ',NMIN,' and minimum ',J1,' is ',DISTPF
                  IF (DISTPF.GT.-STOPDISP) THEN
                     DO J2=1,NMIN-(J1)+1
                        DO J3=1,NOPT
                           QSAVE(J3,J2)=QSAVE(J3,(J1)+J2-1)
                        ENDDO
                        CON(J2)=CON((J1)+J2-1)
                        EMIN(J2)=EMIN((J1)+J2-1)
                        WHICHTS(J2)=WHICHTS((J1)+J2-1)
                     ENDDO
                     NMIN=NMIN-(J1)+1
                     NTS=NMIN-1
                     MFLAG=.TRUE.
                     WRITE(*,'(A,I5,A,I5,A,F12.5,A)') 
     1                       ' displacement between minimum ',NMIN-J1,' and minimum ',(J1),' is ',DISTPF,' stopping'
                     CALL CLEANUP(WHICHTS,TSEN,PLUSEN,MINUSEN,PATHLENGTH,DISP,GAMMA,NTILDE,EMIN,NMIN,DUMPPATH,QSAVE,QSAVETS,
     1                            FSAVETS,FSAVEMIN,MINFRQDONE,TSFRQDONE,BULKT,FILTH,FILTHSTR,TWOD,UNRST)
                     IF (UNRST) DEALLOCATE(MYQMINSAVE)
                     RETURN
                  ELSE
                     GOTO 333
                  ENDIF
               ELSE
                  GOTO 111
               ENDIF
333            CONTINUE  
            ENDDO
         ENDIF
      ENDIF
111   CONTINUE
C
C  Try to find a transition state between minima JDOING and JDOING+1
C
22    CALL TSSEARCH(JDOING,JDOING+1,Q,FIN,BULKT,DTHRESH,FIXD,NATOMS,NIMAGESAVE,PTEST,NSTEPS,FILTH,FILTHSTR,MUPDATE,NSTEPMIN,
     1      MAXNEBBFGS,VNEW,NOSHIFT,NOPT,QSAVE,QSAVETS,CON,EMIN,WHICHTS,GT10,DIAG,STPMAX,MXSTPSAVE,NTS,
     2      NMIN,TSUSED,REVERSET,STOPFIRST,BFGSTST,INR,VECS,ENERGY,EVALMIN,RMS,JDOING,TWOD,REDOPATH,IGUESS,UNRST,NINTS)
      IF (GT10) GOTO 10
C
C  Is it a new transition state?
C
      USEOLD=.FALSE.
      DO J1=1,NTS
         IF (ABS(ENERGY-TSEN(J1)).LT.EDIFFTOL) THEN
            DO J2=1,NOPT
               DUMMY(J2)=QSAVETS(J2,J1)
            ENDDO
            CALL NEWMINDIST(DUMMY,Q,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            IF (DISTPF.GT.GEOMDIFFTOL) GOTO 11
            IF (.NOT.TSUSED(J1)) THEN
               WRITE(*,'(A,I4,A,I4,A,I4,A)') ' Same as transition state ',J1,' currently unused'
C
C  Trying to reuse transition states seems to cause cycling. Could introduce a logical ALLOWREUSE
C  and try reusing every other time...?
C
C              IF (ZSYMSAVE(1:1).EQ.'W') THEN
C              USEOLD=.TRUE.
C              JUSE=J1
C              DO J2=1,NOPT
C                 Q(J2)=QSAVETS(J2,J1)
C              ENDDO
C              GOTO 30
C              ELSE
               IF (GUESSTST.AND.(.NOT.TRYNEB)) THEN
                  WRITE(*,'(A,I4,A,I4,A)') ' Same as transition state ',J1,' try again with a different twist'
                  IGUESS=IGUESS+1
                  IF (IGUESS.GT.ABS(NGUESS)) THEN
                     WRITE(*,'(A,I5)') ' Too many ts guesses - try neb'
C                    CALL REMMIN(JDOING,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,
C 1                              REVERSET,STOPFIRST)
                     IGUESS=1
                     WRITE(*,'(A,I4)') 'IGUESS=',IGUESS
                     TRYNEB=.TRUE.
                     GOTO 10
                  ENDIF
                  DO J2=1,NOPT
                     Q(J2) =QSAVE(J2,JDOING)
                     QSTART(J2) =QSAVE(J2,JDOING)
                     FIN(J2)=QSAVE(J2,JDOING+1)
                  ENDDO
                  GOTO 22
               ELSE
                  IF (NIMAGE+4.GT.NIMAGEMAX) THEN
                     WRITE(*,'(A,I5,A,I5)') ' Same as transition state ',J1,
     1                                         ' too many neb images requested - remove minimum ',JDOING
                     CALL REMMIN(JDOING,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                     NIMAGE=NIMAGESAVE
                  ELSE
                     WRITE(*,'(A,I4,A,I4,A)') ' Same as transition state ',J1,' try again with ',NIMAGE+4,' images'
                     NIMAGE=NIMAGE+4
                     DO J2=1,NOPT
                        Q(J2) =QSAVE(J2,JDOING)
                        QSTART(J2) =QSAVE(J2,JDOING)
                        FIN(J2)=QSAVE(J2,JDOING+1)
                     ENDDO
                     GOTO 22
                  ENDIF
                  GOTO 10
               ENDIF
            ELSE
               DO J2=1,NOPT
                  STPMAX(J2)=MXSTPSAVE
               ENDDO
C
C  If JDOING is 1 or NMINMAX then we shouldn;t really try to
C  access WHICHTS(JDOING-or+1). Seems OK with the additional
C  tests added before accessing WHICHTS, but the tests should
C  really be nested?
C
               IF ((JDOING.GT.1).AND.(J1.EQ.WHICHTS(JDOING-1))) THEN
                  WRITE(*,'(A,I4,A,I4,A,I4)') 
     1                    ' Same as transition state ',J1,' linking minima ',JDOING-1,' and ',JDOING
                  WRITE(*,'(A,I4,A,I4,A)') ' Try linking minima ',JDOING-1,' and ',JDOING+1,' instead'
20                ESTART=EMIN(JDOING-1)
                  EFIN=EMIN(JDOING+1)
                  DO J2=1,NOPT
                     Q(J2) =QSAVE(J2,JDOING-1)
                     QSTART(J2) =QSAVE(J2,JDOING-1)
                     FIN(J2)=QSAVE(J2,JDOING+1)
                  ENDDO
C
C  Try to find a transition state between minima JDOING-1 and JDOING+1 instead
C
                 CALL TSSEARCH(JDOING-1,JDOING+1,Q,FIN,BULKT,DTHRESH,FIXD,NATOMS,NIMAGESAVE,PTEST,NSTEPS,FILTH,FILTHSTR,
     1           MUPDATE,NSTEPMIN,
     2           MAXNEBBFGS,VNEW,NOSHIFT,NOPT,QSAVE,QSAVETS,CON,EMIN,WHICHTS,GT10,DIAG,STPMAX,MXSTPSAVE,NTS,
     3           NMIN,TSUSED,REVERSET,STOPFIRST,BFGSTST,INR,VECS,ENERGY,EVALMIN,RMS,JDOING,TWOD,REDOPATH,IGUESS,UNRST,NINTS)
                  IF (GT10) GOTO 10
C
C  Is it a new transition state?
C
                  CALL ISNEWTS(ENERGY,TSEN,NOPT,QSAVETS,Q,NATOMS,BULKT,STPMAX,MXSTPSAVE,NIMAGE,NIMAGEMAX,JDOING,
     1            NMIN,QSAVE,CON,EMIN,WHICHTS,TSUSED,REVERSET,STOPFIRST,GT10,NIMAGESAVE,NTS,TWOD,IGUESS)
                  IF (GT10) THEN
                     IF ((NIMAGE.GT.NIMAGESAVE).OR.(IGUESS.GT.ABS(NGUESS))) GOTO 20 ! jmc added abs for unres related purposes...
                     GOTO 10
                  ENDIF

                  TSEN(NTS+1)=ENERGY
                  TSUSED(NTS+1)=.TRUE.
                  DO J2=1,NOPT
                     QSAVETS(J2,NTS+1)=Q(J2)
                  ENDDO
        
                  NTS=NTS+1

                  CALL MKFNAMES(NTS,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)

                  CALL PATH(Q,ENERGY,VNEW,RMS,EVALMIN,VECS,.FALSE.,QPLUS,QMINUS,PTEST,ETS,PLUSEN(NTS),MINUSEN(NTS),
     1         PATHLENGTH(NTS),DISP(NTS),GAMMA(NTS),NTILDE(NTS),FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
C jmc added logical PATHFAILT, which  is set to true if one side of the path fails to converge - try another ts guess 
C rather than stopping the entire run

                  IF (TSFRQDONE) THEN
                     DO J2=1,NOPT
                        FSAVETS(J2,NTS)=FRQSTS(J2)
                     ENDDO
                  ENDIF
C
C  Insert the new minima into the existing set, or reject.
C
                  CALL ADDMIN(JDOING-1,JDOING+1,FIN,QPLUS,QSTART,Q,NATOMS,BULKT,NTS,NOPT,STPMAX,MXSTPSAVE,REVERSET,
     1                        STOPFIRST,CON,WHICHTS,NMIN,QSAVE,MINFRQDONE,FSAVEMIN,GT10,TSUSED,EMIN,NIMAGE,NIMAGESAVE,
     2                        FRQSPLUS,FRQSMINUS,PLUSEN,MINUSEN,NIMAGEMAX,NTS,JDOING,TWOD,REDOPATH,QSAVETS,IGUESS,
     3                        UNRST,PATHFAILT,BFGSTST,ITSTRING)
                  IF (GT10) THEN
                     IF ((NIMAGE.GT.NIMAGESAVE).OR.(IGUESS.GT.ABS(NGUESS))) GOTO 20 ! jmc abs
                     GOTO 10
                  ENDIF

               ELSE IF ((JDOING.LT.NMAXMIN).AND.(J1.EQ.WHICHTS(JDOING+1))) THEN
                  WRITE(*,'(A,I4,A,I4,A,I4)') ' Same as transition state ',J1,' linking minima ',JDOING+1,' and ',JDOING+2
                  WRITE(*,'(A,I4,A,I4,A)') ' Try linking minima ',JDOING,' and ',JDOING+2,' instead'
21                ESTART=EMIN(JDOING)
                  EFIN=EMIN(JDOING+2)
                  DO J2=1,NOPT
                     Q(J2) =QSAVE(J2,JDOING)
                     QSTART(J2) =QSAVE(J2,JDOING)
                     FIN(J2)=QSAVE(J2,JDOING+2)
                  ENDDO
C
C  Try to find a transition state between minima JDOING and JDOING+2 instead
C
             CALL TSSEARCH(JDOING,JDOING+2,Q,FIN,BULKT,DTHRESH,FIXD,NATOMS,NIMAGESAVE,PTEST,NSTEPS,FILTH,FILTHSTR,MUPDATE,NSTEPMIN,
     1                     MAXNEBBFGS,VNEW,NOSHIFT,NOPT,QSAVE,QSAVETS,CON,EMIN,WHICHTS,GT10,DIAG,STPMAX,MXSTPSAVE,NTS,
     2                  NMIN,TSUSED,REVERSET,STOPFIRST,BFGSTST,INR,VECS,ENERGY,EVALMIN,RMS,JDOING,TWOD,REDOPATH,IGUESS,UNRST,NINTS)
                  IF (GT10) GOTO 10
C
C  Is it a new transition state?
C
                  CALL ISNEWTS(ENERGY,TSEN,NOPT,QSAVETS,Q,NATOMS,BULKT,STPMAX,MXSTPSAVE,NIMAGE,NIMAGEMAX,JDOING,
     1                         NMIN,QSAVE,CON,EMIN,WHICHTS,TSUSED,REVERSET,STOPFIRST,GT10,NIMAGESAVE,NTS,TWOD,IGUESS)
                  IF (GT10) THEN
                     IF ((NIMAGE.GT.NIMAGESAVE).OR.(IGUESS.GT.ABS(NGUESS))) GOTO 21 ! jmc abs
                     GOTO 10
                  ENDIF

                  TSEN(NTS+1)=ENERGY
                  TSUSED(NTS+1)=.TRUE.
                  DO J2=1,NOPT
                     QSAVETS(J2,NTS+1)=Q(J2)
                  ENDDO
        
                  NTS=NTS+1
                  CALL MKFNAMES(NTS,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
                  CALL PATH(Q,ENERGY,VNEW,RMS,EVALMIN,VECS,.FALSE.,QPLUS,QMINUS,PTEST,ETS,
     1                      PLUSEN(NTS),MINUSEN(NTS),PATHLENGTH(NTS),DISP(NTS),GAMMA(NTS),NTILDE(NTS),
     2                      FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)

                  IF (TSFRQDONE) THEN
                     DO J2=1,NOPT
                        FSAVETS(J2,NTS)=FRQSTS(J2)
                     ENDDO
                  ENDIF
C
C  Insert the new minima into the existing set, or reject.
C
                  CALL ADDMIN(JDOING,JDOING+2,FIN,QPLUS,QSTART,Q,NATOMS,BULKT,NTS,NOPT,STPMAX,MXSTPSAVE,REVERSET,
     1           STOPFIRST,CON,WHICHTS,NMIN,QSAVE,MINFRQDONE,FSAVEMIN,GT10,TSUSED,EMIN,NIMAGE,NIMAGESAVE,
     2           FRQSPLUS,FRQSMINUS,PLUSEN,MINUSEN,NIMAGEMAX,NTS,JDOING,TWOD,REDOPATH,QSAVETS,IGUESS,
     3           UNRST,PATHFAILT,BFGSTST,ITSTRING)
                  IF (GT10) THEN
                     IF ((NIMAGE.GT.NIMAGESAVE).OR.(IGUESS.GT.ABS(NGUESS))) GOTO 21 ! jmc abs
                     GOTO 10
                  ENDIF

               ELSE
                  IF (GUESSTST.AND.(.NOT.TRYNEB)) THEN
                     WRITE(*,'(A,I4,A,I4,A)') ' Same as transition state ',J1,' try again with a different twist'
                     IGUESS=IGUESS+1
C jmc                     IF (IGUESS.GT.NGUESS) THEN
                     IF (IGUESS.GT.ABS(NGUESS)) THEN
                        WRITE(*,'(A,I5)') ' Too many ts guesses - try neb'
C                       CALL REMMIN(JDOING,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,NATOMS,TSUSED,NMAXTS,
C    1                              REVERSET,STOPFIRST)
                        IGUESS=1
                        WRITE(*,'(A,I4)') 'IGUESS=',IGUESS
                        TRYNEB=.TRUE.
                        GOTO 10
                     ENDIF
                     DO J2=1,NOPT
                        Q(J2) =QSAVE(J2,JDOING)
                        QSTART(J2) =QSAVE(J2,JDOING)
                        FIN(J2)=QSAVE(J2,JDOING+1)
                     ENDDO
                     GOTO 22
                  ELSE
                     IF (NIMAGE+4.GT.NIMAGEMAX) THEN
                        WRITE(*,'(A,I5,A,I5)') ' Same as transition state ',J1,
     1                                         ' too many neb images requested - remove minimum ',JDOING
                        CALL REMMIN(JDOING,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                        NIMAGE=NIMAGESAVE
                     ELSE
                        WRITE(*,'(A,I4,A,I4,A)') ' Same as transition state ',J1,' try again with ',NIMAGE+4,' images'
                        NIMAGE=NIMAGE+4
                        DO J2=1,NOPT
                           Q(J2) =QSAVE(J2,JDOING)
                           QSTART(J2) =QSAVE(J2,JDOING)
                           FIN(J2)=QSAVE(J2,JDOING+1)
                        ENDDO
                        GOTO 22
                     ENDIF
                     GOTO 10
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
11       CONTINUE
      ENDDO
   
      TSEN(NTS+1)=ENERGY 
      TSUSED(NTS+1)=.TRUE.
      DO J1=1,NOPT
         QSAVETS(J1,NTS+1)=Q(J1)
      ENDDO
C
C  The two ends of the path are returned in QPLUS and Q (in commons.h)
C     
      NTS=NTS+1
      CALL MKFNAMES(NTS,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
      CALL PATH(Q,ENERGY,VNEW,RMS,EVALMIN,VECS,.FALSE.,QPLUS,QMINUS,PTEST,ETS,
     1          PLUSEN(NTS),MINUSEN(NTS),PATHLENGTH(NTS),DISP(NTS),GAMMA(NTS),NTILDE(NTS),
     2          FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
 
      IF (TSFRQDONE) THEN
         DO J1=1,NOPT
            FSAVETS(J1,NTS)=FRQSTS(J1)
         ENDDO
      ENDIF
   
      IF (USEOLD) THEN
C
C  For water the paths files are in xyz format. We could convert instead of rerunning the path.
C
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
            CALL MKFNAMES(JUSE,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
C           CALL PATH(Q,ENERGY,VNEW,RMS,EVALMIN,VECS,.FALSE.,QPLUS,QMINUS,PTEST,ETS,
C    1             PLUSEN(JUSE),MINUSEN(JUSE),PATHLENGTH(JUSE),DISP(JUSE),GAMMA(JUSE),NTILDE(JUSE),
C    2             FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
            TSUSED(JUSE)=.TRUE.
            OPEN(UNIT=45,FILE=ITSTRING,STATUS='OLD')
            READ(45,*)
            READ(45,*)
C           DO J1=1,NATOMS  !  WCOMMENT
            DO J1=1,NATOMS/2
               READ(45,'(A2,4X,3F20.10)') DSTRING,OVEC(1),OVEC(2),OVEC(3)
               READ(45,'(A2,4X,3F20.10)') DSTRING,H1VEC(1),H1VEC(2),H1VEC(3)
               READ(45,'(A2,4X,3F20.10)') DSTRING,H2VEC(1),H2VEC(2),H2VEC(3)
               CALL CONVERT2(OVEC,H1VEC,H2VEC,QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3),
C  !  WCOMMENT
C    1                       QPLUS(3*(NATOMS+J1-1)+1),QPLUS(3*(NATOMS+J1-1)+2),QPLUS(3*(NATOMS+J1-1)+3))
     1                       QPLUS(3*(NATOMS/2+J1-1)+1),QPLUS(3*(NATOMS/2+J1-1)+2),QPLUS(3*(NATOMS/2+J1-1)+3))
            ENDDO
C
C  Read coordinates of final minimum from end of path.
C
44          READ(45,*,END=45)
            READ(45,*)
C           DO J1=1,NATOMS
            DO J1=1,NATOMS  !  WCOMMENT
               READ(45,'(A2,4X,3F20.10)') DSTRING,OVEC(1),OVEC(2),OVEC(3)
               READ(45,'(A2,4X,3F20.10)') DSTRING,H1VEC(1),H1VEC(2),H1VEC(3)
               READ(45,'(A2,4X,3F20.10)') DSTRING,H2VEC(1),H2VEC(2),H2VEC(3)
               CALL CONVERT2(OVEC,H1VEC,H2VEC,Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),
C    1                       Q(3*(NATOMS+J1-1)+1),Q(3*(NATOMS+J1-1)+2),Q(3*(NATOMS+J1-1)+3))
     1                       Q(3*(NATOMS/2+J1-1)+1),Q(3*(NATOMS/2+J1-1)+2),Q(3*(NATOMS/2+J1-1)+3))
            ENDDO
            GOTO 44
45          CONTINUE
            CLOSE(45)
         ELSE
C
C  Put the old end points back in the right places.
C
            CALL MKFNAMES(JUSE,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
            TSUSED(JUSE)=.TRUE.
   
            OPEN(UNIT=45,FILE=ITSTRING,STATUS='OLD')
            READ(45,*)
            READ(45,*)
            READ(45,'(A2,4X,3F20.10)') (DSTRING,QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3),J1=1,NATOMS)
C
C  Read coordinates of final minimum from end of path.
C
34          READ(45,*,END=35)
            READ(45,*)
            READ(45,'(A2,4X,3F20.10)') (DSTRING,Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),J1=1,NATOMS)
            GOTO 34
35          CONTINUE
            CLOSE(45)
         ENDIF
      ELSE
         JUSE=NTS
      ENDIF
C
C  Insert the new minima into the existing set, or reject.
C
      CALL ADDMIN(JDOING,JDOING+1,FIN,QPLUS,QSTART,Q,NATOMS,BULKT,NTS,NOPT,STPMAX,MXSTPSAVE,REVERSET,
     1           STOPFIRST,CON,WHICHTS,NMIN,QSAVE,MINFRQDONE,FSAVEMIN,GT10,TSUSED,EMIN,NIMAGE,NIMAGESAVE,
     2           FRQSPLUS,FRQSMINUS,PLUSEN,MINUSEN,NIMAGEMAX,JUSE,JDOING,TWOD,REDOPATH,QSAVETS,IGUESS,
     3           UNRST,PATHFAILT,BFGSTST,ITSTRING)
      IF (GT10) THEN
         IF ((NIMAGE.GT.NIMAGESAVE).OR.(IGUESS.GT.ABS(NGUESS))) THEN ! jmc abs
            ESTART=EMIN(JDOING)
            EFIN=EMIN(JDOING+2)
            DO J2=1,NOPT
               Q(J2) =QSAVE(J2,JDOING)
               QSTART(J2) =QSAVE(J2,JDOING)
               FIN(J2)=QSAVE(J2,JDOING+1)
            ENDDO
            GOTO 22
         ENDIF
         GOTO 10
      ENDIF

      IF (UNRST) DEALLOCATE(MYQMINSAVE)

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Reverse all the lines in a file.
C
      SUBROUTINE REVERSE(IUNIT)
      IMPLICIT NONE
      INTEGER MAX, J, N, IUNIT
      PARAMETER (MAX=50000)
      CHARACTER(LEN=100) LSTRING(MAX)

      N=0
      DO J=1,MAX
         READ(IUNIT,'(A100)',END=20) LSTRING(J)      
         N=N+1
      ENDDO
      PRINT*,'WARNING: file on unit ',IUNIT,' too long for reverse'

20    REWIND(IUNIT)
      DO J=1,N
         WRITE(IUNIT,'(A100)') LSTRING(N-J+1)
      ENDDO
 
      REWIND(IUNIT)

      RETURN
      END

      SUBROUTINE REVERSEBLOCK(IUNIT,NBLOCK)
      IMPLICIT NONE
      INTEGER MAX, J, N, NBLOCK, IUNIT, NBLOCKS
      PARAMETER (MAX=1000)
      CHARACTER(LEN=100) LSTRING(MAX)

      IF (NBLOCK.GT.MAX) THEN
         PRINT*,'WARNING: block size requested in REVERSEBLOCK, ',NBLOCK,' exceeds MAX=',MAX
         RETURN
      ENDIF
C
C  Move to end of file.
C
      N=0
10    READ(IUNIT,'(A100)',END=20) LSTRING(1)
C     PRINT*,'N,LSTRING=',N,LSTRING(1)
      N=N+1
      GOTO 10

20    OPEN(UNIT=35,STATUS='SCRATCH')

      DO J=1,NBLOCK+1
         BACKSPACE(IUNIT)
      ENDDO

      NBLOCKS=N/NBLOCK
C     PRINT*,'N,NBLOCK,NBLOCKS=',N,NBLOCK,NBLOCKS

      N=0
30    DO J=1,NBLOCK
         READ(IUNIT,'(A100)') LSTRING(J)
C        PRINT*,'J,LSTRING=',J,LSTRING(J)
      ENDDO
      N=N+1
      DO J=1,NBLOCK
         WRITE(35,'(A100)') LSTRING(J)
      ENDDO
      IF (N.LT.NBLOCKS) THEN
         DO J=1,2*NBLOCK
            BACKSPACE(IUNIT)
         ENDDO
         GOTO 30
      ENDIF

      REWIND(IUNIT)
      REWIND(35)
      
      DO J=1,NBLOCKS*NBLOCK
         READ(35,'(A100)') LSTRING(1)
C        PRINT*,'J,LSTRING=',J,LSTRING(1)
         WRITE(IUNIT,'(A100)') LSTRING(1)
      ENDDO

      REWIND(IUNIT)
      CLOSE(35)

      RETURN
      END

      SUBROUTINE SHORTCUT(NMIN,QSAVE,NTEST,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,
     1                    MOVE,BULKT,TWOD,SHORT)
      USE KEY,ONLY : RIGIDBODY, GEOMDIFFTOL
      USE COMMONS,ONLY : DEBUG
      IMPLICIT NONE
      DOUBLE PRECISION RMAT(3,3)
      INTEGER J2,J3,J4,NMIN,NOPT,NATOMS,NMAXMIN,WHICHTS(NMAXMIN),NTEST,NMAXTS,MOVE(NMAXMIN)
      LOGICAL CON(NMAXMIN), TSUSED(NMAXTS), BULKT, TWOD, SHORT
      DOUBLE PRECISION QSAVE(NOPT,NMAXMIN), Q1(NOPT), Q2(NOPT), DISTPF, EMIN(NMAXMIN)
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
C
C  Check for any shortcuts in the list of minima.
C
      SHORT=.FALSE.
      DO J2=1,NMIN
         MOVE(J2)=J2
      ENDDO
      DO J2=1,NTEST-1
         DO J3=1,NOPT
            Q1(J3)=QSAVE(J3,NTEST)
            Q2(J3)=QSAVE(J3,J2)
         ENDDO
         CALL NEWMINDIST(Q1,Q2,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         IF (DISTPF.LT.GEOMDIFFTOL) THEN
            WRITE(*,'(A,I5,A,I5,A,F15.5)')' Minima ',NTEST,' and ',J2,' are the same: minimum distance=',DISTPF
            SHORT=.TRUE.
            DO J3=J2,NTEST-1
               IF (CON(J3)) THEN
                  TSUSED(WHICHTS(J3))=.FALSE.
               ENDIF
            ENDDO
            DO J3=J2,NTEST-1
               MOVE(J3)=-1
            ENDDO
C           CON(J2)=CON(NTEST)
C           WHICHTS(J2)=WHICHTS(NTEST)
C           MOVE(NTEST)=J2
C           MOVE(J2)=-1
C           DO J3=J2+1,NMIN-(NTEST-J2)  
            DO J3=J2,NMIN-(NTEST-J2)  
C              MOVE(J3)=-1
               MOVE(J3+NTEST-J2)=J3
               PRINT*,'moving minimum ',J3+NTEST-J2,' to minimum ',J3
               DO J4=1,NOPT
                  QSAVE(J4,J3)=QSAVE(J4,J3+NTEST-J2)
               ENDDO
               CON(J3)=CON(J3+NTEST-J2)
               EMIN(J3)=EMIN(J3+NTEST-J2)
               WHICHTS(J3)=WHICHTS(J3+NTEST-J2)
            ENDDO
            NMIN=NMIN-(NTEST-J2)
            RETURN
         ENDIF
      ENDDO

      DO J2=NTEST+1,NMIN
         DO J3=1,NOPT
            Q1(J3)=QSAVE(J3,NTEST)
            Q2(J3)=QSAVE(J3,J2)
         ENDDO
         CALL NEWMINDIST(Q1,Q2,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         IF (DISTPF.LT.GEOMDIFFTOL) THEN
            WRITE(*,'(A,I5,A,I5,A,F15.5)')' Minima ',NTEST,' and ',J2,' are the same: minimum distance=',DISTPF
            DO J3=NTEST,J2-1
               IF (CON(J3)) TSUSED(WHICHTS(J3))=.FALSE.
            ENDDO
            DO J3=NTEST,J2-1
               MOVE(J3)=-1
            ENDDO
C           CON(NTEST)=CON(J2)
C           WHICHTS(NTEST)=WHICHTS(J2)
C           MOVE(J2)=NTEST
C           MOVE(NTEST)=-1
C           DO J3=NTEST+1,NMIN-(J2-NTEST)  
            DO J3=NTEST,NMIN-(J2-NTEST)  
               MOVE(J3+J2-NTEST)=J3
               PRINT*,'moving minimum ',J3+J2-NTEST,' to minimum ',J3
               DO J4=1,NOPT
                  QSAVE(J4,J3)=QSAVE(J4,J3+J2-NTEST)
               ENDDO
               CON(J3)=CON(J3+J2-NTEST)
               EMIN(J3)=EMIN(J3+J2-NTEST)
               WHICHTS(J3)=WHICHTS(J3+J2-NTEST)
            ENDDO
            NMIN=NMIN-(J2-NTEST)
            RETURN
         ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE REMMIN(MIN,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
      use porfuncs
      IMPLICIT NONE
      INTEGER NMAXMIN, NMAXTS
      INTEGER MIN,NMIN,NOPT,WHICHTS(NMAXMIN), J2, J3
      DOUBLE PRECISION QSAVE(NOPT,NMAXMIN), EMIN(NMAXMIN)
      LOGICAL CON(NMAXMIN), TSUSED(NMAXTS), REVERSET, STOPFIRST
      DOUBLE PRECISION STOPDISP
      LOGICAL STOPDISPT
      COMMON /STOPD/ STOPDISP, STOPDISPT

      IF (NMIN.EQ.2) THEN
         PRINT*,'Cannot reduce minima below two - quit'
         STOP
      ENDIF
      IF (MIN.EQ.1) THEN
         PRINT*,'Cannot remove first minimum - remove the second one instead'
         MIN=2
      ENDIF
      IF (MIN.EQ.NMIN) THEN
         PRINT*,'Cannot remove last minimum - remove the penultimate one instead'
         MIN=NMIN-1
      ENDIF
      IF (CON(MIN)) THEN
         CON(MIN)=.FALSE.
         TSUSED(WHICHTS(MIN))=.FALSE.
      ENDIF
      IF (CON(MIN-1)) THEN
         TSUSED(WHICHTS(MIN-1))=.FALSE.
         CON(MIN-1)=.FALSE.
      ENDIF
      DO J2=MIN,NMIN-1
         DO J3=1,NOPT
            QSAVE(J3,J2)=QSAVE(J3,J2+1)
         ENDDO
         CON(J2)=CON(J2+1)
         EMIN(J2)=EMIN(J2+1)
         WHICHTS(J2)=WHICHTS(J2+1)
      ENDDO
      NMIN=NMIN-1
      IF (.NOT.STOPDISPT) THEN
         IF (REVERSET) THEN
            REVERSET=.FALSE.
         ELSE
            IF (.NOT.STOPFIRST) REVERSET=.TRUE.
         ENDIF
      ENDIF

      CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)

      RETURN
      END

      SUBROUTINE REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
      IMPLICIT NONE
      INTEGER J2,J3,J4,NMIN,NMAXMIN,WHICHTS(NMAXMIN),NOPT
      LOGICAL CON(NMAXMIN)
      DOUBLE PRECISION QSAVE(NOPT,NMAXMIN), EMIN(NMAXMIN)
C
C  If we have a minimum that is not connected to the one before it or the one after it, then
C  remove it.
C
10    CONTINUE
      IF (NMIN.GT.3) THEN
         DO J2=2,NMIN-1
            IF ((.NOT.CON(J2-1)).AND.(.NOT.CON(J2))) THEN
               WRITE(*,'(A,I5)') ' Removing isolated minimum ',J2
               DO J3=J2,NMIN-1
                  DO J4=1,NOPT
                     QSAVE(J4,J3)=QSAVE(J4,J3+1)
                  ENDDO
                  CON(J3)=CON(J3+1)
                  EMIN(J3)=EMIN(J3+1)
                  WHICHTS(J3)=WHICHTS(J3+1)
               ENDDO
               NMIN=NMIN-1
               GOTO 10 
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END

      SUBROUTINE CLEANUP(WHICHTS,TSEN,PLUSEN,MINUSEN,PATHLENGTH,DISP,GAMMA,NTILDE,EMIN,NMIN,DUMPPATH,
     1   QSAVE,QSAVETS,FSAVETS,FSAVEMIN,MINFRQDONE,TSFRQDONE,BULKT,FILTH,FILTHSTR,TWOD,UNRST)
      use porfuncs
      USE KEY,ONLY : RIGIDBODY, EDIFFTOL, GEOMDIFFTOL, NOFRQS, SDT
      USE COMMONS
      USE SYMINF
      USE MODHESS
      USE MODCHARMM
      USE MODUNRES
      IMPLICIT NONE
      DOUBLE PRECISION RMAT(3,3)
      INTEGER J1, J2, J3, J4, NMAXMIN, NMAXTS, NMIN, HORDER, FILTH, NATOMSSAVE
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      INTEGER WHICHTS(NMAXMIN), INFO, IPOT
      LOGICAL DUMPPATH, MINFRQDONE, TSFRQDONE, BULKT, TWOD, DOREVERSE
      CHARACTER(LEN=80) ITSTRING,ITSTRING3
      CHARACTER(LEN=80) EOFSSTRING
      CHARACTER(LEN=100) DUMMYST,DUMMYST2
      CHARACTER(LEN=80) PINFOSTRING
      CHARACTER(LEN=*) FILTHSTR
      DOUBLE PRECISION TSEN(NMAXTS), PLUSEN(NMAXTS),SSHIFT,DUMMY1S,MINUSEN(NMAXTS),PATHLENGTH(NMAXMIN),DISP(NMAXMIN),RMS,
     1                 GAMMA(NMAXMIN),NTILDE(NMAXMIN),DUMMY1,DUMMY2,EMIN(NMAXMIN),GRAD(3*NATOMS),TEMPA(9*NATOMS),
     2                 DISTPF, QSAVETS(3*NATOMS,NMAXTS), QSAVE(3*NATOMS,NMAXMIN), INERTIA(3,3),
     3                 FSAVETS(3*NATOMS,NMAXTS), FSAVEMIN(3*NATOMS,NMAXMIN),DIAG(3*NATOMS),
     4                 QMIN(3*NATOMS), QTS(3*NATOMS), Q(3*NATOMS), H1VEC(3), H2VEC(3), OVEC(3), DIHE, ALLANG
      INTEGER NATOMSIMUL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QW, COORDS, PCOORDS, PEPCOORDS
      DOUBLE PRECISION CMXA, CMYA, CMZA
      CHARACTER(LEN=5) ZSYMSAVE
      CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:) :: ZSYML
      COMMON /SYS/ ZSYMSAVE
      INTEGER K1,K2,KD,NNZ,NINTB ! jmc for unres
      LOGICAL UNRST ! jmc for unres

C     ZSYMSAVE=ZSYM(1)
C     IMUL=1   !  WCOMMENT
C     IF (ZSYMSAVE(1:1).EQ.'W') IMUL=3   !  WCOMMENT
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
         NATOMSIMUL=(NATOMS/2)*3  !  number of atoms for water
         ALLOCATE(COORDS(9*(NATOMS/2)),PCOORDS(9*(NATOMS/2)),QW(9*(NATOMS/2)),ZSYML(3*(NATOMS/2)))
      ELSE
         NATOMSIMUL=NATOMS
         ALLOCATE(COORDS(3*NATOMS),PCOORDS(3*NATOMS),ZSYML(NATOMS))
         IF (UNRST) ALLOCATE(PEPCOORDS(3*NATOMS))
      ENDIF

      WRITE(*,'(A)') 'Connected path found'
      WRITE(*,'(A)')
     1'  ts        E+         Ets - E+          Ets       Ets - E-          E-          S       D' //
     2 '      gamma   ~N'
      IF (FILTH.EQ.0) THEN 
         WRITE(ITSTRING,  '(A8)') 'path.xyz'
         IF (UNRST) WRITE(ITSTRING3,  '(A12)') 'path.unr.xyz'
         WRITE(EOFSSTRING,'(A4)')    'EofS'
      ELSE
         WRITE(ITSTRING,  '(A)') 'path.xyz.' // TRIM(ADJUSTL(FILTHSTR))
         WRITE(EOFSSTRING,'(A)') 'EofS.' // TRIM(ADJUSTL(FILTHSTR))
         IF (UNRST) WRITE(ITSTRING3,  '(A)') 'path.unr.xyz.' // TRIM(ADJUSTL(FILTHSTR))
      ENDIF
      OPEN(UNIT=82,FILE=EOFSSTRING,STATUS='UNKNOWN')
      OPEN(UNIT=83,FILE=ITSTRING,STATUS='UNKNOWN')
      IF (UNRST) OPEN(UNIT=86,FILE=ITSTRING3,STATUS='UNKNOWN')
      SSHIFT=0.0D0
C
C  Set PCOORDS to the coordinates of the starting minimum. This enables us
C  to test if the first pathway needs to be reversed. Comparing the energy
C  alone is insufficient for degenerate rearrangements.
C
      IF (ZSYMSAVE(1:1).NE.'W') THEN
         PCOORDS(1:3*NATOMSIMUL)=QSAVE(1:3*NATOMSIMUL,1) 
      ELSE ! must convert to Cartesians for TIP
         DO J2=1,NATOMS/2
            CALL CONVERT(QSAVE(3*(J2-1)+1,1),QSAVE(3*(J2-1)+2,1),QSAVE(3*(J2-1)+3,1),
     1                   QSAVE(3*(NATOMS/2+J2-1)+1,1),QSAVE(3*(NATOMS/2+J2-1)+2,1),QSAVE(3*(NATOMS/2+J2-1)+3,1),
     2                   OVEC,H1VEC,H2VEC)
            PCOORDS(9*(J2-1)+1)=OVEC(1)
            PCOORDS(9*(J2-1)+2)=OVEC(2)
            PCOORDS(9*(J2-1)+3)=OVEC(3)
            PCOORDS(9*(J2-1)+4)=H1VEC(1)
            PCOORDS(9*(J2-1)+5)=H1VEC(2)
            PCOORDS(9*(J2-1)+6)=H1VEC(3)
            PCOORDS(9*(J2-1)+7)=H2VEC(1)
            PCOORDS(9*(J2-1)+8)=H2VEC(2)
            PCOORDS(9*(J2-1)+9)=H2VEC(3)
         ENDDO
      ENDIF
C
C  Main loop over minima.
C
      DO J1=1,NMIN-1 
         J2=WHICHTS(J1)
         CALL MKFNAMES(J2,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
         OPEN(UNIT=84,FILE=EOFSSTRING,STATUS='OLD')
         OPEN(UNIT=85,FILE=ITSTRING,STATUS='OLD')
         READ(85,'(A100)') DUMMYST
         READ(85,'(A100)') DUMMYST2
         DO J4=1,NATOMSIMUL
            READ(85,'(A2,4X,3F20.10)') ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),COORDS(3*(J4-1)+3)
         ENDDO
C
C  Cannot just set PCOORDS = COORDS because we may need to reverse the path.
C
C  Call MIND to align the final minimum of the previous path with the first minimum of the next path; these
C  two minima should be identical of course! To put the other frames in the path into the same orientation
C  we need to save the overall rotation matrix from mind and apply the same rotation - this is done in
C  ROTGEOM. 
C
C  For W4/W3/W2/W1 COORDS contains the Cartesian coordinates, so we don;t want
C  to convert them again in MIND. However, ZSYM(1) has been overwritten at this
C  point, so it won;t contain W any more. Need to turn RIGIDBODY off in the call for 
C  water.
C
         DOREVERSE=.TRUE.
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
            CALL NEWMINDIST(PCOORDS,COORDS,NATOMSIMUL,DISTPF,BULKT,TWOD,ZSYML(1),.FALSE.,.FALSE.,DEBUG,RMAT)
         ELSE
            CALL NEWMINDIST(PCOORDS,COORDS,NATOMSIMUL,DISTPF,BULKT,TWOD,ZSYML(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
         ENDIF
C        PRINT '(A,G20.10)','DISTPF=',DISTPF
C        PRINT '(A,9F15.10)','RMAT=',RMAT(1:3,1:3)
C        IF (DEBUG) PRINT '(A,G20.10)','Initial DISTPF=',DISTPF
         IF (DISTPF.LE.GEOMDIFFTOL) THEN ! there should be no need to reverse the path!
                                         ! unless both distances are less than DISTPF, sigh
            DOREVERSE=.FALSE.
            IF (ABS(PLUSEN(J2)-EMIN(J1)).GT.EDIFFTOL) THEN
               PRINT'(A)','ERROR - distance and energy tolerances are not consistent'
               PRINT'(A,5G20.10)','DISTPF,GEOMDIFFTOL,PLUSEN,EMIN,EDIFFTOL=',DISTPF,GEOMDIFFTOL,PLUSEN(J2),EMIN(J1),EDIFFTOL
               PRINT '(A)','try reversing the path anyway'
C              STOP
               DOREVERSE=.TRUE.
               GOTO 22
            ENDIF
            WRITE(*,60) J2,PLUSEN(J2),TSEN(J2)-PLUSEN(J2),TSEN(J2),TSEN(J2)-MINUSEN(J2),MINUSEN(J2),
     1                  PATHLENGTH(J2),DISP(J2),GAMMA(J2),NTILDE(J2)
60          FORMAT(I4,1X,F16.7,G12.5,F16.7,G12.5,F16.7,F8.3,F8.3,F8.3,F8.3)
            DO J3=1,50000
               READ(84,*,END=22) DUMMY1, DUMMY2
               IF (J3.EQ.1) THEN
                  DUMMY1S=DUMMY1
               ELSE
                  SSHIFT=SSHIFT+DUMMY1-DUMMY1S
                  DUMMY1S=DUMMY1
               ENDIF
               WRITE(82,'(2F20.10)') SSHIFT, DUMMY2
            ENDDO
22          CONTINUE
         ENDIF
         IF (DOREVERSE) THEN
            CALL REVERSE(84)
            REWIND(85)
            CALL REVERSEBLOCK(85,NATOMSIMUL+2)
            READ(85,'(A100)') DUMMYST
            READ(85,'(A100)') DUMMYST2
            DO J4=1,NATOMSIMUL
               READ(85,'(A2,4X,3F20.10)') ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),COORDS(3*(J4-1)+3)
            ENDDO

            IF (ZSYMSAVE(1:1).EQ.'W') THEN
               CALL NEWMINDIST(PCOORDS,COORDS,NATOMSIMUL,DISTPF,BULKT,TWOD,ZSYML(1),.FALSE.,.FALSE.,DEBUG,RMAT)
            ELSE
               CALL NEWMINDIST(PCOORDS,COORDS,NATOMSIMUL,DISTPF,BULKT,TWOD,ZSYML(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
            ENDIF
C           PRINT '(A,G20.10)','DISTPF=',DISTPF
C           PRINT '(A,9F15.10)','RMAT=',RMAT(1:3,1:3)

C           IF (DEBUG) PRINT '(A,G20.10,A)','DISTPF for reversed path=',DISTPF ! DJW
            IF ((DISTPF.GT.GEOMDIFFTOL).OR.(ABS(MINUSEN(J2)-EMIN(J1)).GT.EDIFFTOL)) THEN
               PRINT'(A)','ERROR - distance and energy tolerances are not consistent for reversed path'
               PRINT'(A,5G20.10)','DISTPF,GEOMDIFFTOL,MINUSEN,EMIN,EDIFFTOL=',DISTPF,GEOMDIFFTOL,MINUSEN(J2),EMIN(J1),EDIFFTOL
               STOP
            ENDIF

            WRITE(*,60) J2,MINUSEN(J2),TSEN(J2)-MINUSEN(J2),TSEN(J2),TSEN(J2)-PLUSEN(J2),PLUSEN(J2),
     1                  PATHLENGTH(J2),DISP(J2),GAMMA(J2),NTILDE(J2)
            DO J3=1,50000
               READ(84,*,END=24) DUMMY1, DUMMY2
               IF (J3.EQ.1) THEN
                  DUMMY1S=DUMMY1
               ELSE
                  SSHIFT=SSHIFT-DUMMY1+DUMMY1S
                  DUMMY1S=DUMMY1
               ENDIF
               WRITE(82,'(2F20.10)') SSHIFT, DUMMY2
            ENDDO
24          CONTINUE
         ENDIF
         IF (J1.EQ.1) THEN ! jmc write coords of first minimum
            WRITE(83,'(A100)') DUMMYST
            WRITE(83,'(A100)') DUMMYST2

            CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
            DO J2=1,NATOMSIMUL
               CMXA=CMXA+COORDS(3*(J2-1)+1)
               CMYA=CMYA+COORDS(3*(J2-1)+2)
               CMZA=CMZA+COORDS(3*(J2-1)+3)
            ENDDO
            CMXA=CMXA/NATOMSIMUL; CMYA=CMYA/NATOMSIMUL; CMZA=CMZA/NATOMSIMUL
            DO J2=1,NATOMSIMUL
               COORDS(3*(J2-1)+1)=COORDS(3*(J2-1)+1)-CMXA
               COORDS(3*(J2-1)+2)=COORDS(3*(J2-1)+2)-CMYA
               COORDS(3*(J2-1)+3)=COORDS(3*(J2-1)+3)-CMZA
            ENDDO
            IF (UNRST) THEN
               WRITE(86,'(I4)') 2*NATOMS-2 ! jmc this should be the number of atoms, 
                                           ! so for carbons plus dummy peptide atoms, want 2*N-2...
               WRITE(86,'(A100)') DUMMYST2
            ENDIF
            DO J4=1,NATOMSIMUL
               WRITE(83,'(A2,4X,3F25.15)') ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),COORDS(3*(J4-1)+3)
               IF (UNRST) WRITE(86,'(A2,4X,3F20.10)') ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),
     &                                                COORDS(3*(J4-1)+3)
            ENDDO
            IF (UNRST) THEN
               DO K1=1,(NATOMS/2)-1 ! jmc add peptide atoms...
                  DO K2=1,3
                     PEPCOORDS(6*(K1-1)+K2)=(2.0D0*COORDS(6*(K1-1)+K2)+COORDS(6*K1+K2))/3.0D0
                     PEPCOORDS(6*(K1-1)+K2+3)=(COORDS(6*(K1-1)+K2)+2.0D0*COORDS(6*K1+K2))/3.0D0
                  END DO
               END DO
               WRITE(86,'(A2,4X,3F20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3)
     &                                   ,K1=1,(NATOMS/2)-1)
               WRITE(86,'(A2,4X,3F20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6)
     &                                   ,K1=1,(NATOMS/2)-1)
            ENDIF
         ENDIF

         DO J3=1,50000
            READ(85,'(A100)',END=23) DUMMYST
            WRITE(83,'(A100)') DUMMYST
            READ(85,'(A100)') DUMMYST
            WRITE(83,'(A100)') DUMMYST
            IF (UNRST) THEN
               WRITE(86,'(I4)') 2*NATOMS-2
               WRITE(86,'(A100)') DUMMYST
            ENDIF
            DO J4=1,NATOMSIMUL
               READ(85,'(A2,4X,3F20.10)',END=23) ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),COORDS(3*(J4-1)+3)
            ENDDO
            CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
            DO J2=1,NATOMSIMUL
               CMXA=CMXA+COORDS(3*(J2-1)+1)
               CMYA=CMYA+COORDS(3*(J2-1)+2)
               CMZA=CMZA+COORDS(3*(J2-1)+3)
            ENDDO
            CMXA=CMXA/NATOMSIMUL; CMYA=CMYA/NATOMSIMUL; CMZA=CMZA/NATOMSIMUL
C
C  Because PCOORDS could contain a rotated version of the first minimum the first
C  RMAT might not be the identity. Must therefore call NEWROTGEOM.
C
C           IF ((J1.GT.1).AND.(.NOT.UNRST)) CALL NEWROTGEOM(NATOMSIMUL,COORDS,RMAT,CMXA,CMYA,CMZA)
            IF (.NOT.UNRST) CALL NEWROTGEOM(NATOMSIMUL,COORDS,RMAT,CMXA,CMYA,CMZA)
            DO J2=1,NATOMSIMUL
               COORDS(3*(J2-1)+1)=COORDS(3*(J2-1)+1)-CMXA
               COORDS(3*(J2-1)+2)=COORDS(3*(J2-1)+2)-CMYA
               COORDS(3*(J2-1)+3)=COORDS(3*(J2-1)+3)-CMZA
            ENDDO
            DO J4=1,NATOMSIMUL
               WRITE(83,'(A2,4X,3F25.15)') ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),COORDS(3*(J4-1)+3)
               IF (UNRST) WRITE(86,'(A2,4X,3F20.10)') ZSYML(J4),COORDS(3*(J4-1)+1),COORDS(3*(J4-1)+2),
     &                                                COORDS(3*(J4-1)+3)
            ENDDO
            IF (UNRST) THEN
               DO K1=1,(NATOMS/2)-1 ! jmc add peptide atoms...
                  DO K2=1,3
                     PEPCOORDS(6*(K1-1)+K2)=(2.0D0*COORDS(6*(K1-1)+K2)+COORDS(6*K1+K2))/3.0D0
                     PEPCOORDS(6*(K1-1)+K2+3)=(COORDS(6*(K1-1)+K2)+2.0D0*COORDS(6*K1+K2))/3.0D0
                  ENDDO
               ENDDO
               WRITE(86,'(A2,4X,3F20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3)
     &                                   ,K1=1,(NATOMS/2)-1)
               WRITE(86,'(A2,4X,3F20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6)
     &                                   ,K1=1,(NATOMS/2)-1)
            ENDIF
         ENDDO
23       CONTINUE
C
C  Save the last frame to align the next path.
C
         CLOSE(84)
         CLOSE(85)
         DO J2=1,3*NATOMSIMUL
            PCOORDS(J2)=COORDS(J2)
         ENDDO
      ENDDO
      CLOSE(82)
      CLOSE(83)
C
C  Finished writing path.xyz and EofS
C
      IF (UNRST) CLOSE(86)

      IF (DUMPPATH) THEN
         CLOSE(88) ! path.info was opened on unit 88 in keywords.f for DUMPALLPATHS
         IF (FILTH.EQ.0) THEN
            WRITE(PINFOSTRING,'(A9)') 'path.info'
         ELSE
            WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
         ENDIF
         OPEN(UNIT=83,FILE=PINFOSTRING,STATUS='UNKNOWN')
         DO J1=1,NMIN
            IF (UNRST.AND.CALCDIHE) THEN
               DO J2=1,3*NATOMS
                  QMIN(J2)=QSAVE(J2,J1)
               ENDDO
               CALL UNRESCALCDIHEREF(DIHE,ALLANG,QMIN)
C               WRITE(83,'(3F25.15)') EMIN(J1), DIHE, ALLANG
               WRITE(83,'(3F25.15)') EMIN(J1), DIHE
            ELSE
               WRITE(83,'(F25.15)') EMIN(J1)
            ENDIF
            IF (ZSYMSAVE(1:1).EQ.'W') THEN
               IF (ZSYMSAVE.EQ.'W4') IPOT=4
               IF (ZSYMSAVE.EQ.'W3') IPOT=3
               IF (ZSYMSAVE.EQ.'W2') IPOT=2
               IF (ZSYMSAVE.EQ.'W1') IPOT=1
C              DO J2=1,NATOMS  !  WCOMMENT
               DO J2=1,NATOMS/2
                  CALL CONVERT(QSAVE(3*(J2-1)+1,J1),QSAVE(3*(J2-1)+2,J1),QSAVE(3*(J2-1)+3,J1),
C ! WCOMMENT
C    1                         QSAVE(3*(NATOMS+J2-1)+1,J1),QSAVE(3*(NATOMS+J2-1)+2,J1),QSAVE(3*(NATOMS+J2-1)+3,J1),
     1                         QSAVE(3*(NATOMS/2+J2-1)+1,J1),QSAVE(3*(NATOMS/2+J2-1)+2,J1),QSAVE(3*(NATOMS/2+J2-1)+3,J1),
     2                         OVEC,H1VEC,H2VEC)
                  QW(9*(J2-1)+1)=OVEC(1)
                  QW(9*(J2-1)+2)=OVEC(2)
                  QW(9*(J2-1)+3)=OVEC(3)
                  QW(9*(J2-1)+4)=H1VEC(1)
                  QW(9*(J2-1)+5)=H1VEC(2)
                  QW(9*(J2-1)+6)=H1VEC(3)
                  QW(9*(J2-1)+7)=H2VEC(1)
                  QW(9*(J2-1)+8)=H2VEC(2)
                  QW(9*(J2-1)+9)=H2VEC(3)
               ENDDO
C              NATOMS=NATOMS*3  ! WCOMMENT
               NATOMSSAVE=NATOMS
               NATOMS=(NATOMS/2)*3
               CALL SYMMETRY(HORDER,.FALSE.,QW,INERTIA)
C              NATOMS=NATOMS/3  !  WCOMMENT
               NATOMS=2*(NATOMS/3)
               NATOMS=NATOMSSAVE
               WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
C              DO J2=1,6*NATOMS  !  WCOMMENT
               DO J2=1,3*NATOMS
                  Q(J2)=QSAVE(J2,J1)
                  QMIN(J2)=Q(J2)
               ENDDO
C              CALL H2OMODES(NATOMS,IPOT,Q,DIAG)  !  WCOMMENT
C              WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,6*NATOMS)  !  WCOMMENT
               IF (.NOT.NOFRQS) THEN
                  CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
                  WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)  !  WCOMMENT
               ENDIF
            ELSE IF ((.NOT.MINFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST) THEN
               DO J2=1,3*NATOMS
                  Q(J2)=QSAVE(J2,J1)
                  QMIN(J2)=Q(J2)
               ENDDO
               IF (UNRST) THEN ! jmc update coords
                  DO K1=1,nres
                      c(1,K1)=Q(6*(K1-1)+1)
                      c(2,K1)=Q(6*(K1-1)+2)
                      c(3,K1)=Q(6*(K1-1)+3)
                      c(1,K1+nres)=Q(6*(K1-1)+4)
                      c(2,K1+nres)=Q(6*(K1-1)+5)
                      c(3,K1+nres)=Q(6*(K1-1)+6)
                   ENDDO
                   CALL UPDATEDC
                   CALL int_from_cart(.true.,.false.)
                   CALL chainbuild
               ENDIF
!              CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               IF (CHRMMT) THEN
                  HORDER=1
                  FPGRP='C1'
                  WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                  IF (.NOT.NOFRQS) THEN
                     IF (ENDNUMHESS) THEN
                        CALL MAKENUMHESS(Q,NATOMS)
                     ELSE
                        CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                     ENDIF
                     CALL MASSWT2(NATOMS,ATMASS,COORDS,GRAD,.TRUE.)
                     CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                     if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
                     WRITE(83,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                  ENDIF
               ELSEIF (UNRST) THEN
                  HORDER=1
                  FPGRP='C1'
                  WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                  IF (.NOT.NOFRQS) THEN
                     IF (ENDNUMHESS) THEN
                        CALL MAKENUMINTHESS(NINTS,NATOMS)
                        CALL GETSTUFF(KD,NNZ,NINTB)
                        CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,DIAG)
                     ELSE
                        CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                     ENDIF
                     WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                  ENDIF
               ELSE
                  CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
                  WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                  IF (.NOT.NOFRQS) THEN
                     IF (ENDNUMHESS) THEN
                        CALL MAKENUMHESS(Q,NATOMS)
                     ELSE
                        CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                     ENDIF
                     CALL MASSWT(NATOMS,ATMASS,COORDS,GRAD,.TRUE.)
                     CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                     if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
                     IF (SDT) THEN
                        WRITE(83,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                     ELSE
                        WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               DO J2=1,3*NATOMS
                  Q(J2)=QSAVE(J2,J1)
                  QMIN(J2)=Q(J2)
               ENDDO
               CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
               WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
               IF (.NOT.NOFRQS) WRITE(83,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
            ENDIF
            IF (J1.GT.1) CALL NEWMINDIST(QTS,QMIN,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            WRITE(83,'(3F25.15)') (QMIN(J2),J2=1,NOPT)
            IF (J1.LT.NMIN) THEN
               J3=WHICHTS(J1)
               WRITE(83,'(F25.15)') TSEN(J3)
               IF (ZSYMSAVE(1:1).EQ.'W') THEN
                  IF (ZSYMSAVE.EQ.'W4') IPOT=4
                  IF (ZSYMSAVE.EQ.'W3') IPOT=3
                  IF (ZSYMSAVE.EQ.'W2') IPOT=2
                  IF (ZSYMSAVE.EQ.'W1') IPOT=1
C                 DO J2=1,NATOMS  !  WCOMMENT
                  DO J2=1,NATOMS/2
                     CALL CONVERT(QSAVETS(3*(J2-1)+1,J3),QSAVETS(3*(J2-1)+2,J3),QSAVETS(3*(J2-1)+3,J3),
C    1                            QSAVETS(3*(NATOMS+J2-1)+1,J3),QSAVETS(3*(NATOMS+J2-1)+2,J3),QSAVETS(3*(NATOMS+J2-1)+3,J3),
     1                            QSAVETS(3*(NATOMS/2+J2-1)+1,J3),QSAVETS(3*(NATOMS/2+J2-1)+2,J3),QSAVETS(3*(NATOMS/2+J2-1)+3,J3),
     2                            OVEC,H1VEC,H2VEC)
                     QW(9*(J2-1)+1)=OVEC(1)
                     QW(9*(J2-1)+2)=OVEC(2)
                     QW(9*(J2-1)+3)=OVEC(3)
                     QW(9*(J2-1)+4)=H1VEC(1)
                     QW(9*(J2-1)+5)=H1VEC(2)
                     QW(9*(J2-1)+6)=H1VEC(3)
                     QW(9*(J2-1)+7)=H2VEC(1)
                     QW(9*(J2-1)+8)=H2VEC(2)
                     QW(9*(J2-1)+9)=H2VEC(3)
                  ENDDO
C                 NATOMS=NATOMS*3 ! WCOMMENT
                  NATOMSSAVE=NATOMS
                  NATOMS=(NATOMS/2)*3
                  CALL SYMMETRY(HORDER,.FALSE.,QW,INERTIA)
C                 NATOMS=NATOMS/3 ! WCOMMENT
                  NATOMS=2*(NATOMS/3)
                  NATOMS=NATOMSSAVE
                  WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
C                 DO J2=1,6*NATOMS  !  WCOMMENT
                  DO J2=1,3*NATOMS
                     Q(J2)=QSAVETS(J2,J3)
                     QTS(J2)=Q(J2)
                  ENDDO
C                 CALL H2OMODES(NATOMS,IPOT,Q,DIAG) ! WCOMMENT
                  IF (.NOT.NOFRQS) THEN
                     CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
                     WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                  ENDIF
               ELSE IF ((.NOT.TSFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST) THEN
                  DO J2=1,3*NATOMS
                     Q(J2)=QSAVETS(J2,J3)
                     QTS(J2)=Q(J2)
                  ENDDO
                  IF (UNRST) THEN
                    DO K1=1,nres ! jmc update coords
                        c(1,K1)=Q(6*(K1-1)+1)
                        c(2,K1)=Q(6*(K1-1)+2)
                        c(3,K1)=Q(6*(K1-1)+3)
                        c(1,K1+nres)=Q(6*(K1-1)+4)
                        c(2,K1+nres)=Q(6*(K1-1)+5)
                        c(3,K1+nres)=Q(6*(K1-1)+6)
                     ENDDO
                     CALL UPDATEDC
                     CALL int_from_cart(.true.,.false.)
                     CALL chainbuild
                  ENDIF
!                 CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                  IF (CHRMMT) THEN
                     HORDER=1
                     FPGRP='C1'
                     WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                     IF (.NOT.NOFRQS) THEN
                        IF (ENDNUMHESS) THEN
                           CALL MAKENUMHESS(Q,NATOMS)
                        ELSE
                           CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                        ENDIF
                        CALL MASSWT2(NATOMS,ATMASS,COORDS,GRAD,.TRUE.)
                        CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                        if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
                        WRITE(83,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                     ENDIF
                  ELSEIF (UNRST) THEN
                     HORDER=1
                     FPGRP='C1'
                     WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                     IF (.NOT.NOFRQS) THEN
                        IF (ENDNUMHESS) THEN
                           CALL MAKENUMINTHESS(NINTS,NATOMS)
                           CALL GETSTUFF(KD,NNZ,NINTB)
                           CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,DIAG)
                        ELSE
                           CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                        ENDIF
                        DO J2=1,NINTS-1
                           IF (DIAG(J2).LT.0.0D0) PRINT *,'Higher order saddle found in pathway - ts ',J3,'eigenvalue ',DIAG(J2)
                        END DO
                        WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                     ENDIF
                  ELSE
                     CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
                     WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                     IF (.NOT.NOFRQS) THEN
                        IF (ENDNUMHESS) THEN
                           CALL MAKENUMHESS(Q,NATOMS)
                        ELSE
                           CALL POTENTIAL(Q,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                        ENDIF
                        CALL MASSWT(NATOMS,ATMASS,COORDS,GRAD,.TRUE.)
                        CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                        if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
                        IF (SDT) THEN
                           WRITE(83,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                        ELSE
                           WRITE(83,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
                  DO J2=1,3*NATOMS
                     Q(J2)=QSAVETS(J2,J3)
                     QTS(J2)=Q(J2)
                  ENDDO
                  CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
                  WRITE(83,'(I6,1X,A4)') HORDER,FPGRP
                  IF (.NOT.NOFRQS) WRITE(83,'(3G20.10)') (FSAVETS(J2,WHICHTS(J1)),J2=1,3*NATOMS)
               ENDIF
               CALL NEWMINDIST(QMIN,QTS,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
               WRITE(83,'(3F25.15)') (QTS(J2),J2=1,NOPT)
            ENDIF
         ENDDO
         CLOSE(83)
      ENDIF

      IF (ALLOCATED(COORDS)) DEALLOCATE(COORDS)
      IF (ALLOCATED(PCOORDS)) DEALLOCATE(PCOORDS)
      IF (ALLOCATED(QW)) DEALLOCATE(QW)
      IF (ALLOCATED(ZSYML)) DEALLOCATE(ZSYML)
      IF (ALLOCATED(PEPCOORDS)) DEALLOCATE(PEPCOORDS)

      RETURN
      END

      SUBROUTINE RATES(NATOMS,NINTS)
      use porfuncs
      USE KEY
      IMPLICIT NONE
      INTEGER J1, J2, NATOMS, NMAXMIN, NMAXTS, NMIN
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      DOUBLE PRECISION TSEN(NMAXTS),EMIN(NMAXMIN),QSAVETS(3*NATOMS,NMAXTS),QSAVE(3*NATOMS,NMAXMIN),
     1                 FSAVETS(3*NATOMS,NMAXTS),FSAVEMIN(3*NATOMS,NMAXMIN), DUMMYF, DUMMYB, FTS, FMINF, FMINB,
     2                 LNKF(NMAXTS), LNKB(NMAXTS), ETA(NMAXTS), THETA(NMAXTS), WAIT(NMAXTS),
     3                 CLASSF, CLASSB, LNKFQM(NMAXTS), LNKBQM(NMAXTS)
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST
      DOUBLE PRECISION TEMPERATURE, HRED
      CHARACTER(LEN=14) PINFOSTRING
      INTEGER NCONNECT, HMIN(NMAXMIN), HTS(NMAXTS)
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
C     COMMON /CONREAL/ QSAVETS, QSAVE, FSAVETS, FSAVEMIN
      INTEGER NINTS ! jmc

      IF (FILTH.EQ.0) THEN
         WRITE(PINFOSTRING,'(A9)') 'path.info'
      ELSE
         WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF
      OPEN(UNIT=83,FILE=PINFOSTRING,STATUS='OLD')
      NMIN=0
      DO J1=1,NMAXMIN
         READ(83,*) EMIN(J1)
         NMIN=NMIN+1
         READ(83,*) HMIN(J1)
         IF (.NOT.NOFRQS) READ(83,'(3F20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
         READ(83,'(3F25.15)') (QSAVE(J2,J1),J2=1,3*NATOMS)
         READ(83,*,END=10) TSEN(J1)
         READ(83,*) HTS(J1)
         IF (.NOT.NOFRQS) READ(83,'(3F20.10)') (FSAVETS(J2,J1),J2=1,3*NATOMS)
         READ(83,'(3F25.15)') (QSAVETS(J2,J1),J2=1,3*NATOMS)
      ENDDO
      WRITE(*,'(A)') ' ERROR - path contains too many minima '
      STOP

10    CLOSE(83)
      WRITE(*,'(A,I3,A,I3,A)') ' path.info contains details of ',NMIN,' minima linked by ',NMIN-1,' transition states'
      WRITE(*,'(A,F15.5,A,F15.5)') ' Calculated classical canonical rates in reduced units at T*=',TEMPERATURE,' h*=',HRED
      WRITE(*,'(A,F15.5)') ' Rotational degrees of freedom ignored'
C
C  How many zero eigenvalues?
C
      IF (RTEST) THEN
         NZERO=2
         IF (JZ.NE.0.0D0) THEN
            NZERO=4
         ENDIF
      ELSE IF (MIEFT) THEN
         NZERO=0
      ELSE IF (NOTRANSROTT) THEN
         NZERO=0
      ELSE IF (BULKT) THEN
         NZERO=3
         IF (TWOD) THEN
            NZERO=2
         ENDIF
C     ELSE IF ((FPGRP.EQ.'DXh'.OR.FPGRP.EQ.'CXv').AND.(ZSYMSAVE(1:1).NE.'W')) THEN
C        NZERO=5
      ELSE IF (VARIABLES) THEN
      ELSE IF (FIELDT) THEN
         NZERO=3
      ELSE IF (TWOD) THEN
         NZERO=3
      ELSE IF (UNRST) THEN
         NZERO=3*NATOMS-NINTS
      ELSE
         NZERO=6
      ENDIF

      WRITE(*,'(A,I3)') ' Number of zero eigenvalues=',NZERO
C     WRITE(*,'(A)') ' Vibrational and Boltzmann contributions to k+ and k-:'

      DO J1=1,NMIN-1
         DUMMYF=0.0D0
         DUMMYB=0.0D0
         CLASSF=0.0D0
         CLASSB=0.0D0
         DO J2=1,3*NATOMS-NZERO-1
            FMINF=0.159154943D0*SQRT(FSAVEMIN(J2,J1))*HRED
            FMINB=0.159154943D0*SQRT(FSAVEMIN(J2,J1+1))*HRED
            FTS=  0.159154943D0*SQRT(FSAVETS(J2,J1))*HRED
            DUMMYF=DUMMYF+FMINF/(2.0D0*TEMPERATURE)-FTS/(2.0D0*TEMPERATURE)
     1                 +LOG((1.0D0-EXP(-FMINF/TEMPERATURE))/(1.0D0-EXP(-FTS/TEMPERATURE)))
            DUMMYB=DUMMYB+FMINB/(2.0D0*TEMPERATURE)-FTS/(2.0D0*TEMPERATURE)
     1                 +LOG((1.0D0-EXP(-FMINB/TEMPERATURE))/(1.0D0-EXP(-FTS/TEMPERATURE)))
            CLASSF=CLASSF+LOG(FMINF/FTS)
            CLASSB=CLASSB+LOG(FMINB/FTS)
         ENDDO
         FMINF=0.159154943D0*SQRT(FSAVEMIN(3*NATOMS-NZERO,J1))*HRED
         FMINB=0.159154943D0*SQRT(FSAVEMIN(3*NATOMS-NZERO,J1+1))*HRED
         CLASSF=CLASSF+LOG(FMINF/HRED)-(TSEN(J1)-EMIN(J1))/TEMPERATURE
         CLASSB=CLASSB+LOG(FMINB/HRED)-(TSEN(J1)-EMIN(J1+1))/TEMPERATURE
         DUMMYF=DUMMYF+LOG(TEMPERATURE/HRED)
     1       +FMINF/(2.0D0*TEMPERATURE)+LOG(1.0D0-EXP(-FMINF/TEMPERATURE))-(TSEN(J1)-EMIN(J1))/TEMPERATURE
         DUMMYB=DUMMYB+LOG(TEMPERATURE/HRED)
     1       +FMINB/(2.0D0*TEMPERATURE)+LOG(1.0D0-EXP(-FMINB/TEMPERATURE))-(TSEN(J1)-EMIN(J1+1))/TEMPERATURE
C        WRITE(*,'(I4,4F20.10)') J1,LOG10(CLASSF),LOG10(DUMMYF),LOG10(CLASSB),LOG10(DUMMYB)
         LNKF(J1)=CLASSF
         LNKB(J1)=CLASSB
         LNKFQM(J1)=DUMMYF
         LNKBQM(J1)=DUMMYB
      ENDDO
C
C  Mean first passage times according to Weiss, Adv. Chem. Phys., 13, 1, 1967
C
      ETA(1)=0.0D0
      THETA(1)=1.0D0
      DO J1=2,NMIN
         ETA(J1)=EXP(LNKF(J1-1)-LNKB(J1-1))*ETA(J1-1)+EXP(-LNKB(J1-1))
         THETA(J1)=THETA(J1-1)*EXP(LNKF(J1-1)-LNKB(J1-1))
      ENDDO
      DUMMYF=0.0D0
      DUMMYB=0.0D0
      WRITE(*,'(A)') ' ts           log10(k+cl)         log10(k-cl)         log10(k+qm)         log10(k-qm)   mfp time(classical)'
      DO J1=2,NMIN
         DUMMYF=DUMMYF+THETA(J1-1)
         DUMMYB=DUMMYB+ETA(J1-1)
         WAIT(J1)=ETA(J1)*DUMMYF/THETA(J1)-DUMMYB
         WRITE(*,'(I3,2X,4F20.10,E20.10)') J1-1,LNKF(J1-1)*0.434294481D0,LNKB(J1-1)*0.434294481D0,
     1                                     LNKFQM(J1-1)*0.434294481D0,LNKBQM(J1-1)*0.434294481D0,WAIT(J1)
      ENDDO

      RETURN
      END

      SUBROUTINE MKFNAMES(NTS,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
      IMPLICIT NONE
      CHARACTER(LEN=*) ITSTRING
      CHARACTER(LEN=*) EOFSSTRING
      CHARACTER(LEN=80) FILTHSTR, NTSSTR
      INTEGER NTS, FILTH

      WRITE(NTSSTR,'(I8)') NTS
      IF (FILTH.EQ.0) THEN
         WRITE(ITSTRING,  '(A)') 'path.'//TRIM(ADJUSTL(NTSSTR))//'.xyz'
         WRITE(EOFSSTRING,'(A)')    'EofS.'//TRIM(ADJUSTL(NTSSTR))
      ELSE
         WRITE(ITSTRING,  '(A)') 'path.'//TRIM(ADJUSTL(NTSSTR))//'.xyz.'//TRIM(ADJUSTL(FILTHSTR))
         WRITE(EOFSSTRING,'(A)')    'EofS.'//TRIM(ADJUSTL(NTSSTR))//'.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Initial printing and termination check.
C
      SUBROUTINE STARTUP(MFLAG,NMIN,NTS,CON,EMIN,WHICHTS,TSEN,NTSUSED,NMINCON,NCDONE,NCONNECT,MSTART,MFINISH,MSTEP,
     1                   REVERSET,STOPFIRST,NOPT,Q,QSAVE,QSTART,FIN,JDOING,ESTART,EFIN,TSUSED)
      use porfuncs
      IMPLICIT NONE
      INTEGER NMAXTS,NMAXMIN,NCDONE,NCONNECT,NOPT,JDOING,J2
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      INTEGER J1,NMIN,NTS,WHICHTS(NMAXMIN),NTSUSED,NMINCON,MSTART,MFINISH,MSTEP
      LOGICAL MFLAG,CON(NMAXMIN),TSUSED(NMAXTS),REVERSET,STOPFIRST
      DOUBLE PRECISION EMIN(NMAXMIN),TSEN(NMAXTS),Q(NOPT),QSAVE(NOPT,NMAXMIN),QSTART(NOPT),FIN(NOPT),
     1                 ESTART,EFIN

      WRITE(*,'(A)') '**********************************************************************************************'
      WRITE(*,'(A)') ' minimum      energy      forward connection?    which ts      ts         energy         used?'
      DO J1=1,MAX(NMIN,NTS)
         IF ((J1.LE.NMIN-1).AND.(J1.LE.NTS)) THEN
            IF (CON(J1)) THEN
                WRITE(*,'(I5,F20.10,L14,I15,I11,F20.10,L7)') J1,EMIN(J1),CON(J1),WHICHTS(J1),J1,TSEN(J1),TSUSED(J1)
            ELSE
                WRITE(*,'(I5,F20.10,L14,15X,I11,F20.10,L7)') J1,EMIN(J1),CON(J1),J1,TSEN(J1),TSUSED(J1)
            ENDIF
         ELSE IF (J1.EQ.NMIN) THEN
            IF (J1.LE.NTS) THEN               
               WRITE(*,'(I5,F20.10,29X,I11,F20.10,L7)') J1,EMIN(J1),J1,TSEN(J1),TSUSED(J1)
            ELSE
               WRITE(*,'(I5,F20.10,L14)') J1,EMIN(J1)
            ENDIF
         ELSE IF (J1.LE.NMIN-1) THEN
            IF (CON(J1)) THEN
               WRITE(*,'(I5,F20.10,L14,I15)') J1,EMIN(J1),CON(J1),WHICHTS(J1)
            ELSE
               WRITE(*,'(I5,F20.10,L14)') J1,EMIN(J1),CON(J1)
            ENDIF
         ELSE
             WRITE(*,'(54X,I11,F20.10,L7)') J1,TSEN(J1),TSUSED(J1)
         ENDIF
      ENDDO
      WRITE(*,'(A)') '**********************************************************************************************'
      NTSUSED=0
      DO J1=1,NTS
         IF (TSUSED(J1)) NTSUSED=NTSUSED+1
      ENDDO
      NMINCON=0
      DO J1=1,NMIN
         IF (CON(J1)) THEN
            NMINCON=NMINCON+1
            IF (.NOT.TSUSED(WHICHTS(J1))) THEN
               PRINT*,' Failed consistency check, minimum ',J1,' ts ',WHICHTS(J1),' is not true'
               STOP
            ENDIF
         ENDIF
      ENDDO
      IF (NMINCON.NE.NTSUSED) THEN
         PRINT*,' Failed consistency check, NTSUSED, NMINCON=',NTSUSED, NMINCON
         STOP
      ENDIF
      IF (NTS.EQ.NMAXTS) THEN
         WRITE(*,'(A)') 'Maximum transition states found - exit connect'
         STOP
      ELSE IF (NCDONE.GT.NCONNECT) THEN
         PRINT*,'Number of allowed connection steps reached'

C        STOP
c jmc stop program stopping here         STOP
C jmc make a path.info file now for the sequence of connected minima (if it exists!)
C
         DO J1=MSTART,MFINISH,MSTEP
            IF (.NOT.CON(J1)) STOP
         ENDDO   
         MFLAG=.TRUE.
         RETURN
C jmc end new stuff
      ENDIF

      MSTART=1
      MFINISH=NMIN-1
      MSTEP=1
      IF (REVERSET) THEN
         MSTART=NMIN-1
         MFINISH=1
         MSTEP=-1
      ENDIF
C
C  IF STOPFIRST is true then we want to stop as soon as the first minimum is connected.
C
      MFLAG=.FALSE.
      IF (STOPFIRST.AND.CON(1)) THEN
         NMIN=2
      ELSE
         DO J1=MSTART,MFINISH,MSTEP
            IF (.NOT.CON(J1)) THEN
               DO J2=1,NOPT
                  Q(J2) =QSAVE(J2,J1)
                  QSTART(J2) =QSAVE(J2,J1)
                  FIN(J2)=QSAVE(J2,J1+1)
               ENDDO
               JDOING=J1
               ESTART=EMIN(JDOING)
               EFIN=EMIN(JDOING+1)
               RETURN
            ENDIF
         ENDDO
      ENDIF
      MFLAG=.TRUE.

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Transition state search between minima JSTART and JFINISH
C
      SUBROUTINE TSSEARCH(JSTART,JFINISH,Q,FIN,BULKT,DTHRESH,FIXD,NATOMS,NIMAGESAVE,PTEST,NSTEPS,FILTH,FILTHSTR,MUPDATE,NSTEPMIN,
     1     MAXNEBBFGS,VNEW,NOSHIFT,NOPT,QSAVE,QSAVETS,CON,EMIN,WHICHTS,GT10,DIAG,STPMAX,MXSTPSAVE,NTS,
     2                   NMIN,TSUSED,REVERSET,STOPFIRST,BFGSTST,INR,VECS,ENERGY,EVALMIN,RMS,JREMOVE,TWOD,REDOPATH,
     3                   IGUESS,UNRST,NINTS)
      use porfuncs
      USE KEY,ONLY : RIGIDBODY, PERMDIST, EDIFFTOL, GEOMDIFFTOL
      USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3,DEBUG
      USE MODNEB
      use KEYNEB, ONLY: NNNIMAGE=>NIMAGE, OLDCONNECT
      USE NEWNEBMODULE
      USE MODCHARMM
      USE MODUNRES
      USE MODEFOL
      IMPLICIT NONE
      DOUBLE PRECISION RMAT(3,3)

      INTEGER NMAXTS,NMAXMIN,NATOMS,NIMAGESAVE,ITDONE,INR,INRSAVE,NOPT,J1,NMIN,JREMOVE,FILTH,MUPDATE,NSTEPMIN
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      CHARACTER(LEN=*) FILTHSTR
      INTEGER JSTART,JFINISH,NSTEPS,WHICHTS(NMAXMIN),J3,NTS, IGUESS, ITMAX, NINTS
      LOGICAL BULKT,PERM,FIXD,PTEST,BFGSTST,NOSHIFT,CON(NMAXMIN),GT10,MFLAG,TSUSED(NMAXTS),REVERSET,STOPFIRST,FIXDSAVE,TWOD,
     1        REDOPATH, MorePrinting2
      LOGICAL GUESSFAIL, FAILCHECK, GT20
      DOUBLE PRECISION Q(3*NATOMS),FIN(3*NATOMS),DISTPF,EMIN(NMAXMIN),DTHRESH,FIXDIR(3*NATOMS),STEP,DPRAND,
     1                 VNEW(3*NATOMS),ENERGY,RMS,EVALMIN,EVALMAX,VECS(3*NATOMS),DIST2,EMAX,MAXNEBBFGS,
     2 TIME,TIME0,QSAVE(3*NATOMS,NMAXMIN),DIAG(3*NATOMS),STPMAX(3*NATOMS),MXSTPSAVE,QSAVETS(3*NATOMS,NMAXTS)
      DOUBLE PRECISION QSTART(3*NATOMS),RMSTEMP,EINITIAL,EFINAL,TWISTFRAC,SEPARATION,DCOORDS(3*NATOMS)
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
      COMMON /NEBRMS/ RMSTEMP,EINITIAL,EFINAL,SEPARATION
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      DOUBLE PRECISION DIHE,ALLANG
      LOGICAL UNRST
      COMMON /OLDC/ EMAX

      DO J1=1,NOPT ! in case VECS would otherwise be unset
         VECS(J1)=DPRAND()*2-1.0D0
      ENDDO
C 
C The declared dimension of VECS is 3*NATOMS, so the above DO loop is OK but we need to get the dimension right below.
C
      IF (UNRST) THEN
         CALL VECNORM(VECS,NINTS)
      ELSE
         CALL VECNORM(VECS,NOPT)
      ENDIF
      FIXDSAVE=FIXD
      CALL ALIGN_DECIDE(Q,FIN,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DISTPF,DIST2,RIGIDBODY,RMAT)
      EINITIAL=EMIN(JSTART)
      EFINAL=EMIN(JFINISH)
      DO J1=1,NOPT
         QSTART(J1)=Q(J1)
      ENDDO
      WRITE(*,'(A,I5,A,I5,A,F15.5)')' Minima ',JSTART,' and ',JFINISH,' unconnected: minimum distance=',DISTPF
      IF (DISTPF.LT.GEOMDIFFTOL) THEN
         WRITE(*,'(A,I6,A,I6,A)') ' Minima ',JSTART,' and ',JFINISH,' appear to be identical - quit'
         STOP
      ENDIF
      IF (UNRST) THEN
         CALL UNRESCALCDIHE(DIHE,ALLANG,Q,FIN)
         PRINT *,'DIFF ',DIHE,ALLANG
      ENDIF
      
      PERM=.FALSE.
      IF (ABS(EMIN(JSTART)-EMIN(JFINISH)).LT.EDIFFTOL) PERM=.TRUE.
20    IF (REDOPATH) THEN
         DO J1=1,NOPT
            Q(J1)=QSAVETS(J1,NTS+1)
         ENDDO
         KNOWE=.FALSE.
         KNOWG=.FALSE.
         KNOWH=.FALSE.
      ELSE IF (FIXD.AND.(DISTPF.GT.DTHRESH)) THEN
C        IF (.NOT.FIXD) THEN
C           WRITE(*,'(A)') ' These minima appear to be permutational isomers - switching to hard sphere move guess'
C           T12FAC=1.1D0
C        ENDIF
         DO J1=1,NOPT
            FIXDIR(J1)=FIN(J1)-Q(J1)
         ENDDO 
         CALL HSMOVE(Q,Q,FIXDIR,STEP,.TRUE.)
      ELSE
         IF (PERM) THEN 
C           WRITE(*,'(A)') ' These minima appear to be permutational isomers - using two neb images'
            WRITE(*,'(A,I4,A)') ' These minima appear to be permutational isomers - using ',NIMAGE,'  neb images'
C           NIMAGE=3
         ENDIF
         IF (CHRMMT) THEN
            IF (GUESSTST.AND.(.NOT.TRYNEB)) THEN
C last guess is now midway around the long arc from start to finish; all others are equally spaced along shorter arc
               IF (IGUESS.EQ.NGUESS) THEN
                  TWISTFRAC = -1.D0
               ELSE
                  TWISTFRAC=1.0D0*IGUESS/NGUESS
               ENDIF
C              print *,'IGUESS TWISTFRAC',IGUESS,TWISTFRAC
               GUESSFAIL=.FALSE.
               CALL CHGUESSTS(Q,.FALSE.,PTEST,TWISTFRAC,GUESSFAIL,DISTPF)
               FAILT=.FALSE.
               IF (GUESSFAIL) THEN
                  CALL OLDNEB(.FALSE.,PTEST,ENERGY,MAXNEBBFGS,PERM,Q)
C                 CALL NEB(MAXNEBBFGS,MUPDATE,NIMAGE,NSTEPMIN,NSTEPNEB,RMSNEB,NEBMAG,NATOMS,EINITIAL,Q,EFINAL,FIN,PTEST,ENERGY,
C    1                     FILTH)
                  TRYNEB=.TRUE.
                  IGUESS=1
                  WRITE(*,'(A,I4)') 'IGUESS=',IGUESS
               ENDIF
            ELSE
               IF (NEBT) THEN
                  CALL OLDNEB(.FALSE.,PTEST,ENERGY,MAXNEBBFGS,PERM,Q)
                  CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ELSEIF (NEWNEBT) THEN
                  OldConnect=.TRUE. ! changes behaviour of NEWNEB, ts guess should come back in Q
                  CALL NEWNEB(.FALSE.,DCOORDS,EINITIAL,Q,EFINAL,FIN,.TRUE.)
               ENDIF
            ENDIF
         ELSEIF (UNRST) THEN
            IF (GUESSTST.AND.(.NOT.TRYNEB)) THEN
C last guess is now midway around the long arc from start to finish; all others are equally spaced along shorter arc
               IF (IGUESS.EQ.NGUESS) THEN
                  TWISTFRAC = -1.D0
               ELSE
                  TWISTFRAC=1.0D0*IGUESS/NGUESS
               ENDIF
C Apply Hammond's postulate i.e. start nearest to highest energy minimum
C              IF (EMIN(JSTART).LT.EMIN(JFINISH)) TWISTFRAC = 1.0D0 - TWISTFRAC
               print *,'IGUESS TWISTFRAC',IGUESS,TWISTFRAC
               GUESSFAIL=.FALSE.
               IF (CONSECT) THEN
                  CALL UNRESGUESSTSSEC(Q,.FALSE.,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
               ELSE
                  CALL UNRESGUESSTS(Q,.FALSE.,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
               END IF
               FAILT=.FALSE.
               IF (GUESSFAIL) THEN
                  PRINT *,'guessfail true, calling newneb'
                  IF (NEBT) THEN
                     CALL OLDNEB(.FALSE.,PTEST,ENERGY,MAXNEBBFGS,PERM,Q)
                  ELSEIF (NEWNEBT) THEN
                     OldConnect=.TRUE. ! changes behaviour of NEWNEB, ts guess should come back in Q
                     CALL NEWNEB(.FALSE.,DCOORDS,EINITIAL,Q,EFINAL,FIN,.TRUE.)
                  ENDIF
                  TRYNEB=.TRUE.
                  IGUESS=1
                  WRITE(*,'(A,I4)') 'IGUESS=',IGUESS
               ELSE
                  KNOWE=.FALSE. ! jmc testing, to fix rather large bug with intbfgsts.....
                  KNOWG=.FALSE.
               ENDIF
            ELSE
               IF (NEBT) THEN
                  CALL OLDNEB(.FALSE.,PTEST,ENERGY,MAXNEBBFGS,PERM,Q)
               ELSEIF (NEWNEBT) THEN
                  OldConnect=.TRUE. ! changes behaviour of NEWNEB, ts guess should come back in Q
                  CALL NEWNEB(.FALSE.,DCOORDS,EINITIAL,Q,EFINAL,FIN,.TRUE.)
               ENDIF
            ENDIF
         ELSE
            IF (NEBT) THEN
               CALL OLDNEB(.FALSE.,PTEST,ENERGY,MAXNEBBFGS,PERM,Q)
            ELSEIF (NEWNEBT) THEN
               OldConnect=.TRUE. ! changes behaviour of NEWNEB, ts guess should come back in Q
               CALL NEWNEB(.FALSE.,DCOORDS,EINITIAL,Q,EFINAL,FIN,.TRUE.)
            ENDIF
         ENDIF
C        IF (PERM) NIMAGE=NIMAGESAVE
      ENDIF
!
! Must set the energy of the highest image for use in secdiag routine with bfgsts.
!
      IF (NEWNEBT) ENERGY=EMAX
      FIXD=.FALSE.
C
C  Transition state search from points in Q.
C
      CALL MYCPU_TIME(TIME0,.TRUE.)
      IF (BFGSTST) THEN
         IF (UNRST) THEN 
            CALL INTBFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PTEST)
         ELSE
C
C DAE Make number of ts state search steps allowed proportional to the distance between minima
C
            IF (CHRMMT) ITMAX=NSTEPS+NINT(DISTPF*EXTRASTEPS)
            CALL BFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PTEST)
         ENDIF
      ELSE
         INRSAVE=INR
         INR=2
C        NOSHIFT=.FALSE.
         CALL EFOL(Q,MFLAG,NSTEPS,ENERGY,ITDONE,EVALMIN,PTEST,DIAG,2)
         INR=INRSAVE
      ENDIF
      CALL MYCPU_TIME(TIME,.FALSE.)
C
C DAE If EVALMIN large in magnitude, this TS is likely to be bogus, and cause problems
C when then the connected minima have to be found
C
C     IF ((CHRMMT.OR.UNRST).AND.(EVALMIN.LT.-100.D0)) THEN
C        MFLAG=.FALSE.
C        WRITE(*,'(A,F20.10,A)') 'Eigenvalue ',EVALMIN,' too negative, TS search failed'
C     ENDIF
C
C DAE for CHARMM check this transition state to see if its geometry has become unfeasible
C
      IF (CHRMMT) THEN
         CALL CHECKPOINT(Q,FAILCHECK)
         IF (FAILCHECK) THEN
            WRITE(*,'(A)') ' WARNING Transition state may have unphysical geometry'
C           MFLAG=.FALSE.
         ENDIF
      ENDIF

      IF (MFLAG) THEN
         WRITE(*,'(A,I5,A,F19.10,A,F15.5,A,F11.2)')
     1 ' ts search converged in ',ITDONE,' steps. Energy=',ENERGY,' eigenvalue=',EVALMIN,' time=',TIME-TIME0
         GT10=.FALSE.
      ELSE
         WRITE(*,'(A,I5,A)') ' Transition state search failed to converge in ',ITDONE,' steps'
         DO J3=1,NOPT
            STPMAX(J3)=MXSTPSAVE
         ENDDO
         FIXD=FIXDSAVE
C
C  Next block was introduced by DAE to replace the commented part.
C
         CALL TRYAGAINRESET(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1        JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST,
     2        Q,QSTART,GT20)
         IF (GT20) GOTO 20
!         IF (NIMAGE+4.GT.NIMAGEMAX) THEN
!            WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
!            CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
!            NIMAGE=NIMAGESAVE
!         ELSE
!            WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
!            NIMAGE=NIMAGE+4
!C
!C  We have to reset the starting geometry in Q!
!C
!            DO J3=1,NOPT
!               Q(J3)=QSTART(J3)
!            ENDDO
!            GOTO 20
!         ENDIF
         GT10=.TRUE.
         RETURN
      ENDIF

      FIXD=FIXDSAVE
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Subroutine to incorporate new min-saddle-min triple into the path, or reject it.
C
      SUBROUTINE ADDMIN(JSTART,JFINISH,FIN,QPLUS,QSTART,Q,NATOMS,BULKT,NTS,NOPT,STPMAX,MXSTPSAVE,REVERSET,
     1                  STOPFIRST,CON,WHICHTS,NMIN,QSAVE,MINFRQDONE,FSAVEMIN,GT10,TSUSED,EMIN,NIMAGE,NIMAGESAVE,
     2                  FRQSPLUS,FRQSMINUS,PLUSEN,MINUSEN,NIMAGEMAX,JUSE,JREMOVE,TWOD,REDOPATH,QSAVETS,
     3                  IGUESS,UNRST,PATHFAILT,BFGSTST,ITSTRING)
      use porfuncs
      USE MODCHARMM
      use modamber9, only : NOCISTRANSRNA,NOCISTRANSDNA,GOODSTRUCTURE1,GOODSTRUCTURE2,CISARRAY1,CISARRAY2
      use commons, only : ZSYM, DEBUG
      USE KEY,ONLY : RIGIDBODY, EDIFFTOL, GEOMDIFFTOL, AMBERT, NABT, AMBER12T
      USE CHIRALITY
      IMPLICIT NONE
      DOUBLE PRECISION RMAT(3,3)
      INTEGER NMAXMIN, J3, NLOST, J1, JUSE, NMAXTS, IGUESS
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      INTEGER JSTART,JFINISH,NATOMS,NTS,NOPT,WHICHTS(NMAXMIN),NMIN,J2,NIMAGE,NIMAGESAVE,MOVE(NMAXMIN),NIMAGEMAX,NTEMP,
     1        JREMOVE, INVERT, INDEX(NATOMS)
      DOUBLE PRECISION FIN(3*NATOMS),QPLUS(3*NATOMS),Q(3*NATOMS),DISTPS,DISTPF,DISTMS,DISTMF,STPMAX(3*NATOMS),
     1                 MXSTPSAVE,QSTART(3*NATOMS),QSAVE(3*NATOMS,NMAXMIN),FSAVEMIN(3*NATOMS,NMAXMIN),EMIN(NMAXMIN),
     2                 FRQSPLUS(3*NATOMS),FRQSMINUS(3*NATOMS),PLUSEN(NMAXTS),MINUSEN(NMAXTS),DUMMY1,QLOCAL(3*NATOMS),
     3                 QSAVETS(3*NATOMS,NMAXTS), QTEMP(3*NATOMS)
      LOGICAL BULKT,PEQS,MEQF,PEQF,MEQS,REVERSET,STOPFIRST,GT10,CON(NMAXMIN),MINFRQDONE,TSUSED(NMAXTS),TWOD,SHORT,REDOPATH
      LOGICAL CISTRANST,UNRST,PATHFAILT,BFGSTST
      CHARACTER(LEN=5) ZSYMSAVE
      CHARACTER(LEN=5) ZTEMP(NATOMS)
      CHARACTER(LEN=80) DSTRING, TMPSTRING
      CHARACTER(LEN=*) ITSTRING
      COMMON /SYS/ ZSYMSAVE
      DOUBLE PRECISION STOPDISP
      LOGICAL STOPDISPT, SUCCESS
      COMMON /STOPD/ STOPDISP, STOPDISPT

      IF (CHRMMT.AND.NOCISTRANS) THEN
         CISTRANST=.FALSE.
         CALL CHECKCISTRANS(Q,QPLUS,CISTRANST)
         IF (CISTRANST) THEN
            WRITE(*,*) ' connect> Cis-trans isomerisation detected, rejecting this min-sad-min insertion'
            NTS=NTS-1
            CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                   JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
            GT10=.TRUE.
            RETURN
         ENDIF
      ELSE IF ((AMBERT.OR.NABT).AND.NOCISTRANS) THEN
         GOODSTRUCTURE1=.TRUE.
         IF (NOCISTRANSRNA) THEN
            CALL CHECK_CISTRANS_RNA(Q,NATOMS,ZSYM,GOODSTRUCTURE1)
            CALL CHECK_CISTRANS_RNA(QPLUS,NATOMS,ZSYM,GOODSTRUCTURE2)
            IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
               WRITE(*,*) 
     &          ' connect> Cis-trans isomerisation detected in the RNA ring wrt. the original RNA ribose structure, rejecting'
               NTS=NTS-1
               CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                       JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
               GT10=.TRUE.
               RETURN
            ENDIF
         ELSE IF (NOCISTRANSDNA) THEN
            CALL CHECK_CISTRANS_DNA(Q,NATOMS,ZSYM,GOODSTRUCTURE1)
            CALL CHECK_CISTRANS_DNA(QPLUS,NATOMS,ZSYM,GOODSTRUCTURE2)
            IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
               WRITE(*,*) 
     &          ' connect> Cis-trans isomerisation detected in the DNA ring wrt. the original DNA ribose structure, rejecting'
               NTS=NTS-1
               CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                       JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
               GT10=.TRUE.
               RETURN
            ENDIF
C
C      check cis peptide bonds for proteins
C
         ELSE
C        instead of looking for cis-trans isomerisation processes, we now check if the new minima are in a different 
C        cis/trans orientation than the STARTING point. In theory, if we start from a clean database, this should be
C        enough.
            CALL CHECK_CISTRANS_PROTEIN(Q,NATOMS,GOODSTRUCTURE1,MINOMEGA,CISARRAY1)
            CALL CHECK_CISTRANS_PROTEIN(QSTART,NATOMS,GOODSTRUCTURE1,MINOMEGA,CISARRAY2)
            CISARRAY1=CISARRAY1-CISARRAY2
            GOODSTRUCTURE1=.TRUE.
            DO J1=1,NATOMS
               IF (CISARRAY1(J1)/=0) THEN
                  GOODSTRUCTURE1=.FALSE.
                  WRITE(*,'(A,I6)') ' connect> cis-trans isomerisation of a peptide bond detected involving atom ', J1 
               ENDIF
            ENDDO

            IF (.NOT.GOODSTRUCTURE1) THEN
               WRITE(*,'(A)') ' connect> Cis-trans isomerisation of a peptide bond detected, rejecting'
               NTS=NTS-1
               CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                       JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
               GT10=.TRUE.
               RETURN
            ENDIF

            CALL CHECK_CISTRANS_PROTEIN(QPLUS,NATOMS,GOODSTRUCTURE2,MINOMEGA,CISARRAY1)
            CALL CHECK_CISTRANS_PROTEIN(QSTART,NATOMS,GOODSTRUCTURE2,MINOMEGA,CISARRAY2)
            CISARRAY1=CISARRAY1-CISARRAY2
            GOODSTRUCTURE1=.TRUE.
            DO J1=1,NATOMS
               IF (CISARRAY1(J1)/=0) THEN
                  GOODSTRUCTURE1=.FALSE.
                  WRITE(*,'(A,I6)') ' connect> cis-trans isomerisation of a peptide bond detected involving atom ', J1 
               ENDIF
            ENDDO

            IF (.NOT.GOODSTRUCTURE1) THEN
               WRITE(*,'(A)') ' connect> Cis-trans isomerisation of a peptide bond detected, rejecting'
               NTS=NTS-1
               CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                       JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
               GT10=.TRUE.
               RETURN
            ENDIF
          ENDIF
      ELSE IF (AMBER12T .AND. NOCISTRANS) THEN
          GOODSTRUCTURE1 = .TRUE.
          GOODSTRUCTURE2 = .TRUE.
! Check for chiral inversions in DNA/RNA, since cis/trans refers to the substituents on the ribose rings.
          IF (NOCISTRANSRNA .OR. NOCISTRANSDNA) THEN
            CALL CHIRALITY_CHECK(Q, GOODSTRUCTURE1)
            CALL CHIRALITY_CHECK(QPLUS, GOODSTRUCTURE2)
            IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
               WRITE(*,*) 
     &          ' connect> Possible cis-trans isomerisation detected wrt. the starting ribose structure, rejecting'
               NTS=NTS-1
               CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                       JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
               GT10=.TRUE.
               RETURN
            ENDIF
! Otherwise, check for cis/trans flips of peptide bonds
          ELSE
!            CALL INIT_CIS_TRANS(QSTART) ! Already called in keywords
            CALL CIS_TRANS_CHECK(Q, GOODSTRUCTURE1)
            CALL CIS_TRANS_CHECK(QPLUS, GOODSTRUCTURE2)
            IF (.NOT. (GOODSTRUCTURE1 .OR. GOODSTRUCTURE2)) THEN
               WRITE(*,'(A)') ' connect> Cis-trans isomerisation of a peptide bond detected, rejecting'
               NTS=NTS-1
               CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                       JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
               GT10=.TRUE.
               RETURN
            ENDIF
          END IF
      ENDIF

! Check chirality for AMBER aminoacid residues 
      IF (CHECKCHIRALT.AND.(AMBERT.OR.NABT)) THEN
         CALL CHECK_CHIRALITY(Q,NATOMS,GOODSTRUCTURE1)
         CALL CHECK_CHIRALITY(QPLUS,NATOMS,GOODSTRUCTURE2)
         IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
            WRITE(*,'(A)') ' connect> Chirality inversion detected in at least one of the carbon centres, rejecting'
            NTS=NTS-1
            CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                   JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
            GT10=.TRUE.
            RETURN
         ENDIF
! Equivalent for AMBER12
      ELSE IF (CHECKCHIRALT .AND. AMBER12T) THEN
         GOODSTRUCTURE1 = .TRUE.
         GOODSTRUCTURE2 = .TRUE.
         CALL CHIRALITY_CHECK(Q, GOODSTRUCTURE1)
         CALL CHIRALITY_CHECK(QPLUS, GOODSTRUCTURE2)
         IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
            WRITE(*,'(A)') ' connect> Chirality inversion detected in at least one of the carbon centres, rejecting'
            NTS=NTS-1
            CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                   JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
            GT10=.TRUE.
            RETURN
         ENDIF
      ENDIF

C jmc gets here in the case that one side of the path has failed to converge.  Want to try a different ts, rather than 
C stopping completely...
      IF (PATHFAILT) THEN
         PATHFAILT=.FALSE.
C        BFGSTST=.TRUE.
         WRITE(*,*) 'One side of the path failed to converge, rejecting this min-sad-min insertion'
         NTS=NTS-1
         CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
         GT10=.TRUE.
         RETURN
      ENDIF

333   CALL NEWMINDIST(FIN,QPLUS,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
      CALL NEWMINDIST(QSTART,QPLUS,NATOMS,DISTPS,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
      CALL NEWMINDIST(FIN,Q,NATOMS,DISTMF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
      CALL NEWMINDIST(QSTART,Q,NATOMS,DISTMS,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)

      PRINT*
      WRITE(*,'(A,I6)') ' Minimum distances of Starting and Finishing geometries from +/- minima for ts ',JUSE
      PRINT*
      WRITE(*,'(A)')  '            S              F'
      WRITE(*,'(A,2F15.6)')  ' + ',DISTPS,DISTPF
      WRITE(*,'(A,2F15.6)')  ' - ',DISTMS,DISTMF
            
      PEQS=.FALSE.
      MEQF=.FALSE.
      PEQF=.FALSE.
      MEQS=.FALSE.
      IF (DISTPS.LT.GEOMDIFFTOL) PEQS=.TRUE.
      IF (DISTMF.LT.GEOMDIFFTOL) MEQF=.TRUE.
      IF (DISTPF.LT.GEOMDIFFTOL) PEQF=.TRUE.
      IF (DISTMS.LT.GEOMDIFFTOL) MEQS=.TRUE.
C
C  Allow redopath runs to do a permutation-inversion operation to put each new path into
C  correspondence with the existing segment. This means that we don't have to have the right
C  permutation-inversion isomers of the transition states, which is useful for comparison
C  with old master equation results.
C
      IF (REDOPATH.AND.((.NOT.PEQS).AND.(.NOT.MEQS))) THEN
         DO J2=1,3*NATOMS
            QLOCAL(J2)=QSAVE(J2,NMIN-1)
         ENDDO
         IF (DABS(EMIN(NMIN-1)-PLUSEN(NTS)).LT.EDIFFTOL) THEN
            CALL GETPERM(QLOCAL,QPLUS,INVERT,INDEX,SUCCESS)
         ELSE IF (DABS(EMIN(NMIN-1)-MINUSEN(NTS)).LT.EDIFFTOL) THEN
            CALL GETPERM(QLOCAL,Q,INVERT,INDEX,SUCCESS)
         ELSE
            PRINT '(A)',' Problem in redopath run, the last path cannot be properly linked'
            PRINT*,'NMIN,NTS,EMIN(NMIN-1),PLUSEN(NTS),MINUSEN(NTS)=',NMIN,NTS,EMIN(NMIN-1),PLUSEN(NTS),MINUSEN(NTS)
            PRINT '(A)',' Continuing anyway'
!           STOP
            GOTO 999
         ENDIF
         PRINT '(A)',' Trying a permutational isomerisation to fit this segment in'
         DO J2=1,NATOMS
            QLOCAL(3*(J2-1)+1)=INVERT*QPLUS(3*(INDEX(J2)-1)+1)
            IF (TWOD) THEN
               QLOCAL(3*(J2-1)+2)=QPLUS(3*(INDEX(J2)-1)+2)
               QLOCAL(3*(J2-1)+3)=QPLUS(3*(INDEX(J2)-1)+3)
            ELSE
               QLOCAL(3*(J2-1)+2)=INVERT*QPLUS(3*(INDEX(J2)-1)+2)
               QLOCAL(3*(J2-1)+3)=INVERT*QPLUS(3*(INDEX(J2)-1)+3)
            ENDIF
         ENDDO
         DO J2=1,3*NATOMS
            QPLUS(J2)=QLOCAL(J2)
         ENDDO
         DO J2=1,NATOMS
            QLOCAL(3*(J2-1)+1)=INVERT*Q(3*(INDEX(J2)-1)+1)
            IF (TWOD) THEN
               QLOCAL(3*(J2-1)+2)=Q(3*(INDEX(J2)-1)+2)
               QLOCAL(3*(J2-1)+3)=Q(3*(INDEX(J2)-1)+3)
            ELSE
               QLOCAL(3*(J2-1)+2)=INVERT*Q(3*(INDEX(J2)-1)+2)
               QLOCAL(3*(J2-1)+3)=INVERT*Q(3*(INDEX(J2)-1)+3)
            ENDIF
         ENDDO
         DO J2=1,3*NATOMS
            Q(J2)=QLOCAL(J2)
         ENDDO
         DO J2=1,NATOMS
            QLOCAL(3*(J2-1)+1)=INVERT*QSAVETS(3*(INDEX(J2)-1)+1,NTS)
            IF (TWOD) THEN
               QLOCAL(3*(J2-1)+2)=QSAVETS(3*(INDEX(J2)-1)+2,NTS)
               QLOCAL(3*(J2-1)+3)=QSAVETS(3*(INDEX(J2)-1)+3,NTS)
            ELSE
               QLOCAL(3*(J2-1)+2)=INVERT*QSAVETS(3*(INDEX(J2)-1)+2,NTS)
               QLOCAL(3*(J2-1)+3)=INVERT*QSAVETS(3*(INDEX(J2)-1)+3,NTS)
            ENDIF
         ENDDO
         DO J2=1,3*NATOMS
            QSAVETS(J2,NTS)=QLOCAL(J2)
         ENDDO
C
C  We need to fix all the points in the path.n.xyz file as well, or we will have
C  problems when we try to reconstruct the full path.
C
         OPEN(UNIT=81,FILE=ITSTRING,STATUS='OLD')
         TMPSTRING=TRIM(ADJUSTL("temp." // ITSTRING))
         OPEN(UNIT=82,FILE=TMPSTRING,STATUS='UNKNOWN')
222      READ(81,'(A80)',END=444) DSTRING
         WRITE(82,'(A80)') DSTRING
         READ(81,'(A80)') DSTRING
         WRITE(82,'(A80)') DSTRING
         DO J2=1,NATOMS
            READ(81,*) ZTEMP(J2),QTEMP(3*(J2-1)+1),QTEMP(3*(J2-1)+2),QTEMP(3*(J2-1)+3) 
         ENDDO
         DO J2=1,NATOMS
            QLOCAL(1)=INVERT*QTEMP(3*(INDEX(J2)-1)+1)
            IF (TWOD) THEN
               QLOCAL(2)=QTEMP(3*(INDEX(J2)-1)+2)
               QLOCAL(3)=QTEMP(3*(INDEX(J2)-1)+3)
            ELSE
               QLOCAL(2)=INVERT*QTEMP(3*(INDEX(J2)-1)+2)
               QLOCAL(3)=INVERT*QTEMP(3*(INDEX(J2)-1)+3)
            ENDIF
            WRITE(82,'(A2,4X,3F20.10)') ZTEMP(J2),QLOCAL(1),QLOCAL(2),QLOCAL(3)
         ENDDO
         GOTO 222
444      CLOSE(81)
         CLOSE(82)
         CALL SYSTEM('mv ' // TMPSTRING // ' ' // ITSTRING)
C
C  Can we finish the connection to the last minimum?
C
         IF ((DABS(EMIN(NMIN)-PLUSEN(NTS)).LT.EDIFFTOL).OR.(DABS(EMIN(NMIN)-MINUSEN(NTS)).LT.EDIFFTOL)) THEN
            DO J2=1,3*NATOMS
               QLOCAL(J2)=QSAVE(J2,NMIN)
            ENDDO
            IF (DABS(EMIN(NMIN)-PLUSEN(NTS)).LT.EDIFFTOL) THEN
               CALL GETPERM(QPLUS,QLOCAL,INVERT,INDEX,SUCCESS)
            ELSE IF (DABS(EMIN(NMIN)-MINUSEN(NTS)).LT.EDIFFTOL) THEN
               CALL GETPERM(Q,QLOCAL,INVERT,INDEX,SUCCESS)
            ENDIF
            DO J2=1,NATOMS
               QLOCAL(3*(J2-1)+1)=INVERT*QSAVE(3*(INDEX(J2)-1)+1,NMIN)
               IF (TWOD) THEN
                  QLOCAL(3*(J2-1)+2)=QSAVE(3*(INDEX(J2)-1)+2,NMIN)
                  QLOCAL(3*(J2-1)+3)=QSAVE(3*(INDEX(J2)-1)+3,NMIN)
               ELSE
                  QLOCAL(3*(J2-1)+2)=INVERT*QSAVE(3*(INDEX(J2)-1)+2,NMIN)
                  QLOCAL(3*(J2-1)+3)=INVERT*QSAVE(3*(INDEX(J2)-1)+3,NMIN)
               ENDIF
            ENDDO
            DO J2=1,3*NATOMS
               QSAVE(J2,NMIN)=QLOCAL(J2)
            ENDDO
         ENDIF
         GOTO 333
      ENDIF
999   CONTINUE
C
C  End of redopath permutation stuff.
C
      DO J2=1,NOPT
         STPMAX(J2)=MXSTPSAVE
      ENDDO
      IF ((PEQS.AND.MEQF).OR.(PEQF.AND.MEQS)) THEN
         WRITE(*,'(A,I5,A,I5)') ' Pathway links the required minima ',JSTART,' and ',JFINISH
C        IF (REVERSET) THEN
C           REVERSET=.FALSE.
C        ELSE
C           IF (.NOT.STOPFIRST) REVERSET=.TRUE.
C        ENDIF
         DO J3=JSTART,JFINISH-1
            IF (CON(J3)) TSUSED(WHICHTS(J3))=.FALSE.
         ENDDO
         WHICHTS(JSTART)=JUSE
         CON(JSTART)=.TRUE.
         NLOST=JFINISH-JSTART-1
         NMIN=NMIN-NLOST
         IF (NLOST.GT.0) THEN
            DO J2=JSTART+1,NMIN
               DO J3=1,NOPT
                  QSAVE(J3,J2)=QSAVE(J3,J2+NLOST)
               ENDDO
               IF (MINFRQDONE) THEN
                  DO J3=1,NOPT
                     FSAVEMIN(J3,J2)=FSAVEMIN(J3,J2+NLOST)
                  ENDDO
               ENDIF
               CON(J2)=CON(J2+NLOST)
               EMIN(J2)=EMIN(J2+NLOST)
               WHICHTS(J2)=WHICHTS(J2+NLOST)
            ENDDO
         ENDIF
         NIMAGE=NIMAGESAVE
         IF (UNRST) THEN
            TRYNEB=.FALSE.
            IGUESS=1
         ENDIF
         GT10=.TRUE.
         RETURN
      ELSE IF (PEQS.OR.MEQS) THEN
         WRITE(*,'(A,I5)') ' Pathway connected to first minimum ',JSTART
C
C  If it is an alternative transition state linking JSTART to the other neighbouring minimum
C  that is already connected then we want to try the current connection again, otherwise we
C  can go round in circles. JSTART is always < JFINISH, so we only need to compare with the
C  minimum at JSTART-1.
C
         IF (JSTART.GT.1) THEN
            IF (MEQS) THEN ! minus minimum of path corresponds to JSTART
               IF (ABS(PLUSEN(JUSE)-EMIN(JSTART-1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JSTART-1)
                  ENDDO
                  CALL NEWMINDIST(QPLUS,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the other end of the path is the previous minimum
                     IF (CON(JSTART-1)) THEN   ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JSTART-1))=.FALSE.
                     ENDIF
                     WHICHTS(JSTART-1)=JUSE
                     CON(JSTART-1)=.TRUE.
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JSTART,' to minimum ',JSTART-1
                     CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                      JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                     GT10=.TRUE.
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
            IF (PEQS) THEN ! plus minimum of path corresponds to JSTART
               IF (ABS(MINUSEN(JUSE)-EMIN(JSTART-1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JSTART-1)
                  ENDDO
                  CALL NEWMINDIST(Q,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the other end of the path is the previous minimum
                     IF (CON(JSTART-1)) THEN   ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JSTART-1))=.FALSE.
                     ENDIF
                     WHICHTS(JSTART-1)=JUSE
                     CON(JSTART-1)=.TRUE.
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JSTART,' to minimum ',JSTART-1
                     CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                      JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                     GT10=.TRUE.
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (NMIN.EQ.NMAXMIN) THEN
            WRITE(*,'(A)') 'Too many minima - quit'
            STOP
         ENDIF
         DO J2=NMIN,JSTART+1,-1
            DO J3=1,NOPT
               QSAVE(J3,J2+1)=QSAVE(J3,J2)
            ENDDO
            IF (MINFRQDONE) THEN
               DO J3=1,NOPT
                  FSAVEMIN(J3,J2+1)=FSAVEMIN(J3,J2)
               ENDDO
            ENDIF
            CON(J2+1)=CON(J2)
            EMIN(J2+1)=EMIN(J2)
            WHICHTS(J2+1)=WHICHTS(J2)
         ENDDO
         IF (MEQS) THEN
            DO J2=1,NOPT
               QSAVE(J2,JSTART+1)=QPLUS(J2)
            ENDDO
            IF (MINFRQDONE) THEN
               DO J2=1,NOPT
                  FSAVEMIN(J2,JSTART+1)=FRQSPLUS(J2)
               ENDDO
            ENDIF
            EMIN(JSTART+1)=PLUSEN(JUSE)
         ELSE
            DO J2=1,NOPT
               QSAVE(J2,JSTART+1)=Q(J2)
            ENDDO
            IF (MINFRQDONE) THEN
               DO J2=1,NOPT
                  FSAVEMIN(J2,JSTART+1)=FRQSMINUS(J2)
               ENDDO
            ENDIF
            EMIN(JSTART+1)=MINUSEN(JUSE)
         ENDIF
         IF (CON(JSTART)) TSUSED(WHICHTS(JSTART))=.FALSE.
         WHICHTS(JSTART)=JUSE
         WHICHTS(JSTART+1)=-1
         CON(JSTART)=.TRUE.
         CON(JSTART+1)=.FALSE.
         NMIN=NMIN+1
C
C  Check for shortcuts to new minimum that now lies at at position JSTART+1
C
         CALL SHORTCUT(NMIN,QSAVE,JSTART+1,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,
     1                 NMAXTS,MOVE,BULKT,TWOD,SHORT)
         DO J2=1,NOPT
            STPMAX(J2)=MXSTPSAVE
         ENDDO
         NIMAGE=NIMAGESAVE
         IF (UNRST) THEN
            TRYNEB=.FALSE.
            IGUESS=1
         ENDIF
         CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
         GT10=.TRUE.
         RETURN
      ELSE IF (PEQF.OR.MEQF) THEN
         WRITE(*,'(A,I5)') ' Pathway connected to second minimum ',JFINISH
C
C  If it is an alternative transition state linking JFINISH to the other neighbouring minimum
C  that is already connected then we want to try the current connection again, otherwise we
C  can go round in circles. JSTART is always < JFINISH, so we only need to compare with the
C  minimum at JFINISH+1.
C
         IF (JFINISH.LT.NMIN) THEN
            IF (MEQF) THEN ! minus minimum of path corresponds to JFINISH
               IF (ABS(PLUSEN(JUSE)-EMIN(JFINISH+1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JFINISH+1)
                  ENDDO
                  CALL NEWMINDIST(QPLUS,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the other end of the path is the next minimum
                     IF (CON(JFINISH)) THEN  ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JFINISH))=.FALSE.
                     ENDIF
                     WHICHTS(JFINISH)=JUSE
                     CON(JFINISH)=.TRUE.
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JFINISH,' to minimum ',JFINISH+1
                     CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                      JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                     GT10=.TRUE.
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
            IF (PEQF) THEN ! plus minimum of path corresponds to JFINISH
               IF (ABS(MINUSEN(JUSE)-EMIN(JFINISH+1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JFINISH+1)
                  ENDDO
                  CALL NEWMINDIST(Q,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the other end of the path is the next minimum
                     IF (CON(JFINISH)) THEN  ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JFINISH))=.FALSE.
                     ENDIF
                     WHICHTS(JFINISH)=JUSE
                     CON(JFINISH)=.TRUE.
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JFINISH,' to minimum ',JFINISH+1
                     CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                      JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                     GT10=.TRUE.
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (NMIN.EQ.NMAXMIN) THEN
            WRITE(*,'(A)') 'Too many minima - quit'
            STOP
         ENDIF
         DO J2=NMIN,JFINISH,-1
            DO J3=1,NOPT
               QSAVE(J3,J2+1)=QSAVE(J3,J2)
            ENDDO
            IF (MINFRQDONE) THEN
               DO J3=1,NOPT
                  FSAVEMIN(J3,J2+1)=FSAVEMIN(J3,J2)
               ENDDO
            ENDIF
            CON(J2+1)=CON(J2)
            EMIN(J2+1)=EMIN(J2)
            WHICHTS(J2+1)=WHICHTS(J2)
         ENDDO
         IF (MEQF) THEN
            DO J2=1,NOPT
               QSAVE(J2,JFINISH)=QPLUS(J2)
            ENDDO
            IF (MINFRQDONE) THEN
               DO J2=1,NOPT
                  FSAVEMIN(J2,JFINISH)=FRQSPLUS(J2)
               ENDDO
            ENDIF
            EMIN(JFINISH)=PLUSEN(JUSE)
         ELSE
            DO J2=1,NOPT
               QSAVE(J2,JFINISH)=Q(J2)
            ENDDO
            IF (MINFRQDONE) THEN
               DO J2=1,NOPT
                  FSAVEMIN(J2,JFINISH)=FRQSMINUS(J2)
               ENDDO
            ENDIF
            EMIN(JFINISH)=MINUSEN(JUSE)
         ENDIF
         IF (CON(JFINISH-1)) TSUSED(WHICHTS(JFINISH-1))=.FALSE.
         CON(JSTART+1)=.FALSE.
         CON(JFINISH)=.TRUE.
         WHICHTS(JSTART+1)=-1
         WHICHTS(JFINISH)=JUSE
         NMIN=NMIN+1
C
C  Check for shortcuts to first minimum at position JFINISH
C
         CALL SHORTCUT(NMIN,QSAVE,JFINISH,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,
     1                 NMAXTS,MOVE,BULKT,TWOD,SHORT)
         DO J2=1,NOPT
            STPMAX(J2)=MXSTPSAVE
         ENDDO
         NIMAGE=NIMAGESAVE
         IF (UNRST) THEN
            TRYNEB=.FALSE.
            IGUESS=1
         ENDIF
         GT10=.TRUE.
         CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
         RETURN
      ELSE
C
C  This is the block where we proceed differently if we were trying to connect adjacent versus
C  non-adjacent minima. The new pair of minima are added (order depending upon distances from
C  the minima we were trying to connect) if the search was for an adjacent pair, but not otherwise.
C  For non-adjacent minima we try increasing the number of images and ultimately remove minimum
C  JDOING if the upper limit for images has been reached.
C
         WRITE(*,'(A,F20.10,A,I4,A)') ' Pathway is not connected to either minimum'
         IF ((JFINISH-JSTART.GT.1).AND.(.NOT.REDOPATH)) THEN
            NTS=NTS-1
            CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1            JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
            GT10=.TRUE.
            RETURN
         ELSE
C
C  Does it link the same minima as an existing ts?
C
            DO J1=1,NMIN-1
               IF (CON(J1)) THEN
                  IF (((ABS(PLUSEN(JUSE)-EMIN(J1)).LT.EDIFFTOL).AND.(ABS(MINUSEN(JUSE)-EMIN(J1+1)).LT.EDIFFTOL)).OR.
     1                ((ABS(MINUSEN(JUSE)-EMIN(J1)).LT.EDIFFTOL).AND.(ABS(PLUSEN(NTS)-EMIN(J1+1)).LT.EDIFFTOL))) THEN
                     DO J3=1,NOPT
                        STPMAX(J3)=MXSTPSAVE
                     ENDDO
                     WRITE(*,'(A,I5,A,I5,A,I5,A)') ' Pathway links minima ',J1,' and ',J1+1
                     TSUSED(JUSE)=.FALSE.
                     CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1                     JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
                     GT10=.TRUE.
                     RETURN
                  ENDIF
               ENDIF
            ENDDO
C
C  If one of the minima corresponds to minimum JSTART-1 or JFINISH+1 then try replacing 
C  JSTART or JFINISH, as appropriate, with the other one.
C
C
C  If it is an alternative transition state linking JSTART to the other neighbouring minimum
C  that is already connected then we want to try the current connection again, otherwise we
C  can go round in circles. JSTART is always < JFINISH, so we only need to compare with the
C  minimum at JSTART-1 etc.
C
            IF (JSTART.GT.1) THEN
               IF (ABS(PLUSEN(JUSE)-EMIN(JSTART-1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JSTART-1)
                  ENDDO
                  CALL NEWMINDIST(QPLUS,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the plus end of the path is the previous minimum, replace JSTART by minus minimum
                     IF (CON(JSTART-1)) THEN   ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JSTART-1))=.FALSE.
                     ENDIF
                     WHICHTS(JSTART-1)=JUSE
                     CON(JSTART-1)=.TRUE.
                     EMIN(JSTART)=MINUSEN(JUSE)
                     DO J3=1,NOPT
                        QSAVE(J3,JSTART)=Q(J3)
                     ENDDO
                     IF (MINFRQDONE) THEN
                        DO J2=1,NOPT
                           FSAVEMIN(J2,JSTART)=FRQSMINUS(J2)
                        ENDDO
                     ENDIF
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JSTART-1,' to new minimum ',JSTART
C
C  Check for shortcuts to the new minimum at position JSTART. Since JSTART is being replaced we should
C  probably go back to NIMAGESAVE images.
C
                     CALL SHORTCUT(NMIN,QSAVE,JSTART,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,
     1                             MOVE,BULKT,TWOD,SHORT)
                     NIMAGE=NIMAGESAVE
                     IF (UNRST) THEN
                        TRYNEB=.FALSE.
                        IGUESS=1
                     ENDIF
C                    IF (SHORT) THEN
C                       NIMAGE=NIMAGESAVE
C                    ELSE IF (NIMAGE+4.GT.NIMAGEMAX) THEN
C                       WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
C                       CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
C                       NIMAGE=NIMAGESAVE
C                    ELSE
C                       WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
C                       NIMAGE=NIMAGE+4
C                    ENDIF
                     GT10=.TRUE.
                     RETURN
                  ENDIF
               ELSE IF (ABS(MINUSEN(JUSE)-EMIN(JSTART-1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JSTART-1)
                  ENDDO
                  CALL NEWMINDIST(Q,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the minus end of the path is the previous minimum, replace JSTART by minus minimum
                     IF (CON(JSTART-1)) THEN   ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JSTART-1))=.FALSE.
                     ENDIF
                     WHICHTS(JSTART-1)=JUSE
                     CON(JSTART-1)=.TRUE.
                     EMIN(JSTART)=PLUSEN(JUSE)
                     DO J3=1,NOPT
                        QSAVE(J3,JSTART)=QPLUS(J3)
                     ENDDO
                     IF (MINFRQDONE) THEN
                        DO J2=1,NOPT
                           FSAVEMIN(J2,JSTART)=FRQSPLUS(J2)
                        ENDDO
                     ENDIF
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JSTART-1,' to new minimum ',JSTART
C
C  Check for shortcuts to the new minimum at position JSTART.
C
                     CALL SHORTCUT(NMIN,QSAVE,JSTART,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,
     1                             MOVE,BULKT,TWOD,SHORT)
                     NIMAGE=NIMAGESAVE
                     IF (UNRST) THEN
                        TRYNEB=.FALSE.
                        IGUESS=1
                     ENDIF
C                    IF (SHORT) THEN
C                       NIMAGE=NIMAGESAVE
C                    ELSE IF (NIMAGE+4.GT.NIMAGEMAX) THEN
C                       WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
C                       CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
C                       NIMAGE=NIMAGESAVE
C                    ELSE
C                       WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
C                       NIMAGE=NIMAGE+4
C                    ENDIF
                     GT10=.TRUE.
                     CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
                     RETURN
                  ENDIF
               ENDIF
            ENDIF

            IF (JFINISH.LT.NMIN) THEN
               IF (ABS(PLUSEN(JUSE)-EMIN(JFINISH+1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JFINISH+1)
                  ENDDO
                  CALL NEWMINDIST(QPLUS,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the plus end of the path is the next minimum
                     IF (CON(JFINISH)) THEN    ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JFINISH))=.FALSE.
                     ENDIF
                     WHICHTS(JFINISH)=JUSE
                     CON(JFINISH)=.TRUE.
                     EMIN(JFINISH)=MINUSEN(JUSE)
                     DO J3=1,NOPT
                        QSAVE(J3,JFINISH)=Q(J3)
                     ENDDO
                     IF (MINFRQDONE) THEN
                        DO J2=1,NOPT
                           FSAVEMIN(J2,JFINISH)=FRQSMINUS(J2)
                        ENDDO
                     ENDIF
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JFINISH+1,' to new minimum ',JFINISH
C
C  Check for shortcuts to the new minimum at position JFINISH.
C
                     CALL SHORTCUT(NMIN,QSAVE,JFINISH,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,
     1                             MOVE,BULKT,TWOD,SHORT)
                     NIMAGE=NIMAGESAVE
                     IF (UNRST) THEN
                        TRYNEB=.FALSE.
                        IGUESS=1
                     ENDIF
C                    IF (SHORT) THEN
C                       NIMAGE=NIMAGESAVE
C                    ELSE IF (NIMAGE+4.GT.NIMAGEMAX) THEN
C                       WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
C                       CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
C                       NIMAGE=NIMAGESAVE
C                    ELSE
C                       WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
C                       NIMAGE=NIMAGE+4
C                    ENDIF
                     GT10=.TRUE.
                     CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
                     RETURN
                  ENDIF
               ELSE IF (ABS(MINUSEN(JUSE)-EMIN(JFINISH+1)).LT.EDIFFTOL) THEN
                  DO J3=1,NOPT
                     QLOCAL(J3)=QSAVE(J3,JFINISH+1)
                  ENDDO
                  CALL NEWMINDIST(Q,QLOCAL,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
                  IF (DISTPF.LT.GEOMDIFFTOL) THEN  ! the minus end of the path is the next minimum
                     IF (CON(JFINISH)) THEN  ! if there was already a connection replace it with this ts
                        TSUSED(WHICHTS(JFINISH))=.FALSE.
                     ENDIF
                     WHICHTS(JFINISH)=JUSE
                     CON(JFINISH)=.TRUE.
                     EMIN(JFINISH)=PLUSEN(JUSE)
                     DO J3=1,NOPT
                        QSAVE(J3,JFINISH)=QPLUS(J3)
                     ENDDO
                     IF (MINFRQDONE) THEN
                        DO J2=1,NOPT
                           FSAVEMIN(J2,JFINISH)=FRQSPLUS(J2)
                        ENDDO
                     ENDIF
                     WRITE(*,'(A,I4,A,I4)') ' This transition state connects minimum ',JFINISH+1,' to new minimum ',JFINISH
C
C  Check for shortcuts to the new minimum at position JFINISH.
C
                     CALL SHORTCUT(NMIN,QSAVE,JFINISH,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,
     1                             MOVE,BULKT,TWOD,SHORT)
                     NIMAGE=NIMAGESAVE
                     IF (UNRST) THEN
                        TRYNEB=.FALSE.
                        IGUESS=1
                     ENDIF
C                    IF (SHORT) THEN
C                       NIMAGE=NIMAGESAVE
C                    ELSE IF (NIMAGE+4.GT.NIMAGEMAX) THEN
C                       WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
C                       CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
C                       NIMAGE=NIMAGESAVE
C                    ELSE
C                       WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
C                       NIMAGE=NIMAGE+4
C                    ENDIF
                     GT10=.TRUE.
                     CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
   
            IF (.NOT.STOPDISPT) THEN
               IF (REVERSET) THEN
                  REVERSET=.FALSE.
               ELSE
                  IF (.NOT.STOPFIRST) REVERSET=.TRUE.
               ENDIF
            ENDIF

            IF (NMIN.EQ.NMAXMIN-1) THEN
               WRITE(*,'(A)') 'Too many minima - quit'
               STOP
            ENDIF
            DO J2=NMIN,JFINISH,-1
               DO J3=1,NOPT
                  QSAVE(J3,J2+2)=QSAVE(J3,J2)
               ENDDO
               IF (MINFRQDONE) THEN
                  DO J3=1,NOPT
                     FSAVEMIN(J3,J2+2)=FSAVEMIN(J3,J2)
                  ENDDO
               ENDIF
               CON(J2+2)=CON(J2)
               EMIN(J2+2)=EMIN(J2)
               WHICHTS(J2+2)=WHICHTS(J2)
            ENDDO
            DUMMY1=MIN(MIN(DISTPS,DISTMS),MIN(DISTPF,DISTMF))
            IF (DUMMY1.EQ.MIN(DISTPS,DISTMS)) THEN
               IF (DISTPS.LT.DISTMS) THEN
                  DO J1=1,NOPT
                     QSAVE(J1,JFINISH)=QPLUS(J1)
                     QSAVE(J1,JFINISH+1)=Q(J1)
                     EMIN(JFINISH)=PLUSEN(JUSE)
                     EMIN(JFINISH+1)=MINUSEN(JUSE)
                  ENDDO
                  IF (MINFRQDONE) THEN
                     DO J1=1,NOPT
                        FSAVEMIN(J1,JFINISH)=FRQSPLUS(J1)
                        FSAVEMIN(J1,JFINISH+1)=FRQSMINUS(J1)
                     ENDDO
                  ENDIF
               ELSE
                  DO J1=1,NOPT
                     QSAVE(J1,JFINISH)=Q(J1)
                     QSAVE(J1,JFINISH+1)=QPLUS(J1)
                     EMIN(JFINISH+1)=PLUSEN(JUSE)
                     EMIN(JFINISH)=MINUSEN(JUSE)
                  ENDDO
                  IF (MINFRQDONE) THEN
                     DO J1=1,NOPT
                        FSAVEMIN(J1,JFINISH)=FRQSMINUS(J1)
                        FSAVEMIN(J1,JFINISH+1)=FRQSPLUS(J1)
                     ENDDO
                  ENDIF
               ENDIF
            ELSE
               IF (DISTPF.LT.DISTMF) THEN
                  DO J1=1,NOPT
                     QSAVE(J1,JFINISH)=Q(J1)
                     QSAVE(J1,JFINISH+1)=QPLUS(J1)
                     EMIN(JFINISH)=MINUSEN(JUSE)
                     EMIN(JFINISH+1)=PLUSEN(JUSE)
                  ENDDO
                  IF (MINFRQDONE) THEN
                     DO J1=1,NOPT
                        FSAVEMIN(J1,JFINISH)=FRQSMINUS(J1)
                        FSAVEMIN(J1,JFINISH+1)=FRQSPLUS(J1)
                     ENDDO
                  ENDIF
               ELSE
                  DO J1=1,NOPT
                     QSAVE(J1,JFINISH+1)=Q(J1)
                     QSAVE(J1,JFINISH)=QPLUS(J1)
                     EMIN(JFINISH)=PLUSEN(JUSE)
                     EMIN(JFINISH+1)=MINUSEN(JUSE)
                  ENDDO
                  IF (MINFRQDONE) THEN
                     DO J1=1,NOPT
                        FSAVEMIN(J1,JFINISH+1)=FRQSMINUS(J1)
                        FSAVEMIN(J1,JFINISH)=FRQSPLUS(J1)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
            WHICHTS(JFINISH)=JUSE
            CON(JSTART)=.FALSE.
            CON(JFINISH)=.TRUE.
            CON(JFINISH+1)=.FALSE.
            NMIN=NMIN+2
C
C  Check for shortcuts to new minima at positions JSTART+1 and JSTART+2. If the first call
C  finds a shortcut then be careful on the second call!
C
            CALL SHORTCUT(NMIN,QSAVE,JSTART+1,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,
     1                    MOVE,BULKT,TWOD,SHORT)
            NTEMP=MOVE(JSTART+2)
            IF (NTEMP.GT.0) THEN
               CALL SHORTCUT(NMIN,QSAVE,NTEMP,NOPT,NATOMS,NMAXMIN,CON,EMIN,WHICHTS,TSUSED,NMAXTS,MOVE,BULKT,TWOD,SHORT)
            ENDIF
    
            DO J1=1,NOPT
               STPMAX(J1)=MXSTPSAVE
            ENDDO
            NIMAGE=NIMAGESAVE
            IF (UNRST) THEN
               print *,'dae setting TRYNEB F' ! DAE If we have introduced two new minima, set TRYNEB to FALSE
               TRYNEB=.FALSE.
               IGUESS=1
            ENDIF
            GT10=.TRUE.
            CALL REMMIN2(NMIN,CON,NMAXMIN,QSAVE,EMIN,WHICHTS,NOPT)
            RETURN
         ENDIF
      ENDIF
      
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C
      SUBROUTINE ISNEWTS(ENERGY,TSEN,NOPT,QSAVETS,Q,NATOMS,BULKT,STPMAX,MXSTPSAVE,NIMAGE,NIMAGEMAX,JREMOVE,
     1                   NMIN,QSAVE,CON,EMIN,WHICHTS,TSUSED,REVERSET,STOPFIRST,GT10,NIMAGESAVE,NTS,TWOD,IGUESS)
      USE MODCHARMM
      USE KEY,ONLY : RIGIDBODY, EDIFFTOL, GEOMDIFFTOL
      USE COMMONS,ONLY : DEBUG
      IMPLICIT NONE
      DOUBLE PRECISION RMAT(3,3)
      INTEGER NMAXMIN, NMAXTS
      PARAMETER (NMAXTS=1000,NMAXMIN=1001)
      LOGICAL GT10,BULKT,CON(NMAXMIN),TSUSED(NMAXTS),REVERSET,STOPFIRST,TWOD
      INTEGER J2,J3,NOPT,NATOMS,NIMAGE,NIMAGEMAX,NIMAGESAVE,JREMOVE,NMIN,WHICHTS(NMAXMIN),NTS, IGUESS
      DOUBLE PRECISION ENERGY,TSEN(NMAXTS),QSAVETS(3*NATOMS,NMAXTS),Q(3*NATOMS),STPMAX(3*NATOMS),MXSTPSAVE,QSAVE(3*NATOMS,NMAXMIN),
     1                 DUMMY(3*NATOMS),EMIN(NMAXMIN),DISTPF
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

      GT10=.FALSE.
      DO J2=1,NTS
         IF (ABS(ENERGY-TSEN(J2)).LT.EDIFFTOL) THEN
            DO J3=1,NOPT
               DUMMY(J3)=QSAVETS(J3,J2)
            ENDDO
            CALL NEWMINDIST(DUMMY,Q,NATOMS,DISTPF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            IF (DISTPF.GT.GEOMDIFFTOL) GOTO 10
            DO J3=1,NOPT
               STPMAX(J3)=MXSTPSAVE
            ENDDO
            WRITE(*,'(A,I4)') ' Same as transition state ',J2
            CALL TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGEMAX,NIMAGESAVE,
     1           JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,REVERSET,STOPFIRST)
            GT10=.TRUE.
            RETURN
         ENDIF
10       CONTINUE
      ENDDO

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Subroutine to determine the PI operation that will map the coordinates in TPOINTS2
C  onto those in TARGET, and check that tagged atoms are not permuted.
C  Once we have INDEX and INVERT the initial coordinates in
C  TARGET will be in correspondence with those in OLDPOINTSMIN(x,J1) after applying
C  TARGET(3*(INDEX(b)-1)+a)*INVERT = TPOINTS2(3*(b-1)+a,J1)
C  Note that TARGET coordinates are multiplied by -1 for inversion, so INVERT is only
C  for use on other coordinates. The TARGET coordinates are uninverted for OPTIM!!!
C
C
      SUBROUTINE GETPERM(TPOINTS2,TARGET,INVERT,INDEX,SUCCESS)
      use porfuncs
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      INTEGER INDEX(NATOMS), INVERT, J5, J4, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
      DOUBLE PRECISION TPOINTS1(3*NATOMS), TARGET(3*NATOMS), 
     1                 TPOINTS2(3*NATOMS), DIST, DMIN 
      LOGICAL SUCCESS

      SUCCESS=.FALSE.
      INVERT=1
      IF (TWOD) THEN
         CALL STD_ORIENT2D(TPOINTS2,TPOINTS1,NORBIT1,1,NORBIT2,1)
      ELSE IF (BULKT) THEN
         IF (.NOT.FREEZE .AND..NOT.MIEFT .AND. (.NOT. NOTRANSROTT)) THEN
            CALL CENTRE(TARGET)
            CALL CENTRE(TPOINTS2)
         ENDIF
      ELSE
         CALL STD_ORIENT(TPOINTS2,TPOINTS1,NORBIT1,1,NORBIT2,1)
      ENDIF
      NCHOOSE1=1
      NCHOOSE2=1
20    CONTINUE
      IF (TWOD) THEN
         CALL STD_ORIENT2D(TARGET,TPOINTS2,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2)
      ELSE IF (BULKT) THEN
         DO J4=1,3*NATOMS
            TPOINTS1(J4)=TARGET(J4)
         ENDDO
      ELSE
         CALL STD_ORIENT(TARGET,TPOINTS2,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2)
      ENDIF
      DO J4=1,NATOMS
         DMIN=1.0D90
         DO J5=1,NATOMS
            DIST=(TPOINTS1(3*(J4-1)+1)-TPOINTS2(3*(J5-1)+1))**2+
     1           (TPOINTS1(3*(J4-1)+2)-TPOINTS2(3*(J5-1)+2))**2+
     2           (TPOINTS1(3*(J4-1)+3)-TPOINTS2(3*(J5-1)+3))**2
            IF (DIST.LT.DMIN) THEN
               DMIN=DIST
               INDEX(J4)=J5
               IF ((TAGFAC(J4).NE.1.0D0).AND.(J4.NE.J5)) DMIN=1.0D0 ! must not permute different tagged atoms
            ENDIF
         ENDDO
         IF (DEBUG.AND.(J4.NE.INDEX(J4))) WRITE(*,'(A,6I5,F20.10)') 'J4,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,INDEX,DMIN=',
     1                                         J4,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,INDEX(J4),DMIN
         IF (DMIN.GT.0.1D0) GOTO 21
      ENDDO
21    IF (DMIN.GT.0.1D0) THEN
         IF (.NOT.BULKT) THEN
            IF (NORBIT2.GT.1) THEN
               NCHOOSE2=NCHOOSE2+1
               IF (NCHOOSE2.LE.NORBIT2) GOTO 20
            ENDIF
            NCHOOSE2=1
            IF (NORBIT1.GT.1) THEN
               NCHOOSE1=NCHOOSE1+1
               IF (NCHOOSE1.LE.NORBIT1) GOTO 20
            ENDIF
            NCHOOSE1=1
            IF ((INVERT.EQ.1).AND.(.NOT.NOINVERSION)) THEN
               INVERT=-1
               IF (DEBUG.AND.TWOD) WRITE(*,'(A)') 'trying to align reflected permutational isomer'
               IF (DEBUG.AND.(.NOT.TWOD)) WRITE(*,'(A)') 'trying to align inverted permutational isomer'
               DO J5=1,NATOMS
                  TARGET(3*(J5-1)+1)=-TARGET(3*(J5-1)+1)
                  IF (.NOT.TWOD) THEN
                     TARGET(3*(J5-1)+2)=-TARGET(3*(J5-1)+2)
                     TARGET(3*(J5-1)+3)=-TARGET(3*(J5-1)+3)
                  ENDIF
               ENDDO
               GOTO 20
            ELSE
               WRITE(*,'(A)') 'getperm> unable to align permutational isomers, giving up'
               PRINT*,NATOMS
               PRINT*,'TARGET:'
               DO J4=1,NATOMS
                  WRITE(*,'(A2,2X,3F20.10)') ZSYM(J4),TARGET(3*(J4-1)+1),TARGET(3*(J4-1)+2),TARGET(3*(J4-1)+3)
               ENDDO
               PRINT*,NATOMS
               PRINT*,'TPOINTS1:'
               DO J4=1,NATOMS
                  WRITE(*,'(A2,2X,3F20.10)') ZSYM(J4),TPOINTS1(3*(J4-1)+1),TPOINTS1(3*(J4-1)+2),TPOINTS1(3*(J4-1)+3)
               ENDDO
               PRINT*,NATOMS
               PRINT*,'TPOINTS2:'
               DO J4=1,NATOMS
                  WRITE(*,'(A2,2X,3F20.10)') ZSYM(J4),TPOINTS2(3*(J4-1)+1),TPOINTS2(3*(J4-1)+2),TPOINTS2(3*(J4-1)+3)
               ENDDO
               GOTO 55
C              STOP
            ENDIF
         ELSE
            WRITE(*,'(A)') 'getperm> unable to align permutational isomers, giving up'
            PRINT*,NATOMS
            PRINT*,'TARGET:'
            DO J4=1,NATOMS
               WRITE(*,'(A2,2X,3F20.10)') ZSYM(J4),TARGET(3*(J4-1)+1),TARGET(3*(J4-1)+2),TARGET(3*(J4-1)+3)
            ENDDO
            PRINT*,NATOMS
            PRINT*,'TPOINTS1:'
            DO J4=1,NATOMS
               WRITE(*,'(A2,2X,3F20.10)') ZSYM(J4),TPOINTS1(3*(J4-1)+1),TPOINTS1(3*(J4-1)+2),TPOINTS1(3*(J4-1)+3)
            ENDDO
            PRINT*,NATOMS
            PRINT*,'TPOINTS2:'
            DO J4=1,NATOMS
               WRITE(*,'(A2,2X,3F20.10)') ZSYM(J4),TPOINTS2(3*(J4-1)+1),TPOINTS2(3*(J4-1)+2),TPOINTS2(3*(J4-1)+3)
            ENDDO
            GOTO 55
C           STOP
         ENDIF
      ENDIF

      SUCCESS=.TRUE.
      IF (DEBUG) THEN
         WRITE(*,'(A,2I4)') 'alignment successful, INVERT=',INVERT
!        PRINT*,'transformed points:'
!        DO J5=1,NATOMS
!           IF (TWOD) WRITE(*,'(3F20.10)') INVERT*TARGET(3*(INDEX(J5)-1)+1),TARGET(3*(INDEX(J5)-1)+2),TARGET(3*(INDEX(J5)-1)+3)
!           IF (.NOT.TWOD) WRITE(*,'(3F20.10)') INVERT*TARGET(3*(INDEX(J5)-1)+1),INVERT*TARGET(3*(INDEX(J5)-1)+2),
!    1                                          INVERT*TARGET(3*(INDEX(J5)-1)+3)
!        ENDDO
      ENDIF
C
C  Undo the inversion of the TARGET coordinates!
C
55    DO J5=1,NATOMS
         TARGET(3*(J5-1)+1)=INVERT*TARGET(3*(J5-1)+1)
         IF (.NOT.TWOD) THEN
            TARGET(3*(J5-1)+2)=INVERT*TARGET(3*(J5-1)+2)
            TARGET(3*(J5-1)+3)=INVERT*TARGET(3*(J5-1)+3)
         ENDIF
      ENDDO

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Put permutational isomers into a standard orientation.
C
      SUBROUTINE STD_ORIENT2D(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION Q1(3*NATOMS), DIST(NATOMS), DMAX, RVEC(3), T1(3*NATOMS), CMX, CMY, CMZ,
     1                 COST, SINT, RDOTN, TX, TY, TZ
      INTEGER J1, I, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
C
C  Move centre of mass to the origin.
C
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+Q1(3*(I-1)+1)
         CMY=CMY+Q1(3*(I-1)+2)
         CMZ=CMZ+Q1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
      DO I=1,NATOMS
         Q1(3*(I-1)+1)=Q1(3*(I-1)+1)-CMX
         Q1(3*(I-1)+2)=Q1(3*(I-1)+2)-CMY
         Q1(3*(I-1)+3)=Q1(3*(I-1)+3)-CMZ
      ENDDO
      
      NORBIT1=1

      DO J1=1,3*NATOMS
         T1(J1)=Q1(J1)
      ENDDO

C
C  Find the atom with the largest distance from the z axis.
C
      DMAX=-1.0D0
      NORBIT2=1
      DO J1=1,NATOMS
         DIST(J1)=SQRT(T1(3*(J1-1)+1)**2+T1(3*(J1-1)+2)**2)
         IF (ABS(DIST(J1)-DMAX).LT.1.0D-3) THEN
            NORBIT2=NORBIT2+1
            IF (NORBIT2.EQ.NCHOOSE2) THEN
               JMAX2=J1
            ENDIF
         ELSE IF (DIST(J1).GT.DMAX) THEN
            DMAX=DIST(J1)
            NORBIT2=1
            JMAX2=J1
         ENDIF
      ENDDO
C     IF (DEBUG) PRINT*,'NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX2=',NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX2
C
C  and rotate it into the xz plane.
C
C     IF (DEBUG) PRINT*,'atom ',JMAX2,' is now the furthest from the  z axis, distance=',DMAX,' rotate into xz plane'
      IF (T1(3*(JMAX2-1)+2).EQ.0.0D0) GOTO 20
      COST=T1(3*(JMAX2-1)+1)/DMAX
      SINT=T1(3*(JMAX2-1)+2)/DMAX
C     PRINT*,'COST,SINT=',COST,SINT
      RVEC(1)=0.0D0
      RVEC(2)=0.0D0
      RVEC(3)=1.0D0
      DO J1=1,NATOMS
         IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=T1(3*(J1-1)+1)*RVEC(1)+T1(3*(J1-1)+2)*RVEC(2)+T1(3*(J1-1)+3)*RVEC(3)
            TX=T1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)+(T1(3*(J1-1)+2)*RVEC(3)-T1(3*(J1-1)+3)*RVEC(2))*SINT
            TY=T1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)+(T1(3*(J1-1)+3)*RVEC(1)-T1(3*(J1-1)+1)*RVEC(3))*SINT
            TZ=T1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)+(T1(3*(J1-1)+1)*RVEC(2)-T1(3*(J1-1)+2)*RVEC(1))*SINT
            T1(3*(J1-1)+1)=TX
            T1(3*(J1-1)+2)=TY
            T1(3*(J1-1)+3)=TZ
         ENDIF
      ENDDO

20    CONTINUE

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Put permutational isomers into a standard orientation.
C
      SUBROUTINE STD_ORIENT(Q1,T1,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION Q1(3*NATOMS), DIST(NATOMS), DMAX, RVEC(3), T1(3*NATOMS), CMX, CMY, CMZ,
     1                 COST, SINT, RDOTN, TX, TY, TZ
      INTEGER J1, I, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
C
C  Move centre of mass to the origin.
C
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+Q1(3*(I-1)+1)
         CMY=CMY+Q1(3*(I-1)+2)
         CMZ=CMZ+Q1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS
      DO I=1,NATOMS
         Q1(3*(I-1)+1)=Q1(3*(I-1)+1)-CMX
         Q1(3*(I-1)+2)=Q1(3*(I-1)+2)-CMY
         Q1(3*(I-1)+3)=Q1(3*(I-1)+3)-CMZ
      ENDDO
      
      DMAX=-1.0D0
      NORBIT1=1
      DO J1=1,NATOMS
         DIST(J1)=SQRT(Q1(3*(J1-1)+1)**2+Q1(3*(J1-1)+2)**2+Q1(3*(J1-1)+3)**2)
         IF (ABS(DIST(J1)-DMAX).LT.1.0D-3) THEN
            NORBIT1=NORBIT1+1
            IF (NORBIT1.EQ.NCHOOSE1) THEN
               JMAX1=J1
            ENDIF
         ELSE IF (DIST(J1).GT.DMAX) THEN
            DMAX=DIST(J1)
            NORBIT1=1
            JMAX1=J1
         ENDIF
      ENDDO
C
C  For tagged atoms the choice of the first atom matters if it belongs to an orbit of size > 1.
C
C     PRINT*,'atom ',JMAX1,' will be moved onto the z axis, distance=',DMAX
      IF ((ABS(Q1(3*(JMAX1-1)+1)).LT.1.0D-8).AND.(ABS(Q1(3*(JMAX1-1)+2)).LT.1.0D-8)) THEN
         DO J1=1,3*NATOMS
            T1(J1)=Q1(J1)
         ENDDO
         GOTO 10
      ENDIF
      COST=Q1(3*(JMAX1-1)+3)/DMAX
      SINT=SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)/DMAX
C
C  Rotate atom JMAX1 onto the z axis.
C  Rotate all the atoms through ANGLE about RVEC. Use rotation formula
C  from Goldstein p. 165.
C
      RVEC(1)= Q1(3*(JMAX1-1)+2)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
      RVEC(2)=-Q1(3*(JMAX1-1)+1)/SQRT(Q1(3*(JMAX1-1)+1)**2+Q1(3*(JMAX1-1)+2)**2)
      RVEC(3)=0.0D0
      DO J1=1,NATOMS
C        IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=Q1(3*(J1-1)+1)*RVEC(1)+Q1(3*(J1-1)+2)*RVEC(2)+Q1(3*(J1-1)+3)*RVEC(3)
            T1(3*(J1-1)+1)=Q1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)-(Q1(3*(J1-1)+2)*RVEC(3)-Q1(3*(J1-1)+3)*RVEC(2))*SINT
            T1(3*(J1-1)+2)=Q1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)-(Q1(3*(J1-1)+3)*RVEC(1)-Q1(3*(J1-1)+1)*RVEC(3))*SINT
            T1(3*(J1-1)+3)=Q1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)-(Q1(3*(J1-1)+1)*RVEC(2)-Q1(3*(J1-1)+2)*RVEC(1))*SINT
C        ENDIF
      ENDDO

10    CONTINUE
C
C  Find the atom with the largest distance from the z axis.
C
      DMAX=-1.0D0
      NORBIT2=1
      DO J1=1,NATOMS
         DIST(J1)=SQRT(T1(3*(J1-1)+1)**2+T1(3*(J1-1)+2)**2)
C        WRITE(*,'(A,I5,5F20.10)') 'J1,DIST,DMAX,ABS(diff),X,Y=',J1,DIST(J1),DMAX,ABS(DIST(J1)-DMAX),T1(3*(J1-1)+1),T1(3*(J1-1)+2)
         IF (ABS(DIST(J1)-DMAX).LT.1.0D-3) THEN
            NORBIT2=NORBIT2+1
            IF (NORBIT2.EQ.NCHOOSE2) THEN
               JMAX2=J1
            ENDIF
         ELSE IF (DIST(J1).GT.DMAX) THEN
            DMAX=DIST(J1)
            NORBIT2=1
            JMAX2=J1
         ENDIF
      ENDDO
C     PRINT*,'NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX1,JMAX2=',NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,JMAX1,JMAX2
C
C  and rotate it into the xz plane.
C
C     PRINT*,'atom ',JMAX2,' is now the furthest from the  z axis, distance=',DMAX,' rotate into xz plane'
      IF (T1(3*(JMAX2-1)+2).EQ.0.0D0) GOTO 20
      COST=T1(3*(JMAX2-1)+1)/DMAX
      SINT=T1(3*(JMAX2-1)+2)/DMAX
C     PRINT*,'COST,SINT=',COST,SINT
      RVEC(1)=0.0D0
      RVEC(2)=0.0D0
      RVEC(3)=1.0D0
      DO J1=1,NATOMS
         IF (DIST(J1).NE.0.0D0) THEN
            RDOTN=T1(3*(J1-1)+1)*RVEC(1)+T1(3*(J1-1)+2)*RVEC(2)+T1(3*(J1-1)+3)*RVEC(3)
            TX=T1(3*(J1-1)+1)*COST + RVEC(1)*RDOTN*(1-COST)+(T1(3*(J1-1)+2)*RVEC(3)-T1(3*(J1-1)+3)*RVEC(2))*SINT
            TY=T1(3*(J1-1)+2)*COST + RVEC(2)*RDOTN*(1-COST)+(T1(3*(J1-1)+3)*RVEC(1)-T1(3*(J1-1)+1)*RVEC(3))*SINT
            TZ=T1(3*(J1-1)+3)*COST + RVEC(3)*RDOTN*(1-COST)+(T1(3*(J1-1)+1)*RVEC(2)-T1(3*(J1-1)+2)*RVEC(1))*SINT
            T1(3*(J1-1)+1)=TX
            T1(3*(J1-1)+2)=TY
            T1(3*(J1-1)+3)=TZ
         ENDIF
      ENDDO

20    CONTINUE

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Move the centre of mass to the origin.
C
      SUBROUTINE CENTRE(T1)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION CMX, CMY, CMZ, T1(3*NATOMS)
      INTEGER I
C
C  Move centre of mass to the origin.
C

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO I=1,NATOMS
         CMX=CMX+T1(3*(I-1)+1)
         CMY=CMY+T1(3*(I-1)+2)
         CMZ=CMZ+T1(3*(I-1)+3)
      ENDDO
      CMX=CMX/NATOMS
      CMY=CMY/NATOMS
      CMZ=CMZ/NATOMS

      DO I=1,NATOMS
         T1(3*(I-1)+1)=T1(3*(I-1)+1)-CMX
         T1(3*(I-1)+2)=T1(3*(I-1)+2)-CMY
         T1(3*(I-1)+3)=T1(3*(I-1)+3)-CMZ
      ENDDO

      RETURN
      END

      SUBROUTINE TRYAGAIN(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGELIMIT,NIMAGESAVE,
     1                   JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,
     2                   REVERSET,STOPFIRST)
      IMPLICIT NONE

      INTEGER NATOMS,IGUESS,NGUESS,NIMAGE,NIMAGELIMIT,NIMAGESAVE,JREMOVE,NMIN,NOPT,NMAXMIN,NMAXTS
      INTEGER WHICHTS(NMAXMIN)
      DOUBLE PRECISION QSAVE(3*NATOMS,NMAXMIN),EMIN(NMAXMIN)
      LOGICAL GUESSTST,TRYNEB,REVERSET,STOPFIRST,CON(NMAXMIN),TSUSED(NMAXTS)

      IF (GUESSTST.AND.(.NOT.TRYNEB)) THEN
         IGUESS=IGUESS+1
         WRITE(*,'(A)') ' Try again with next twist'
         IF (IGUESS.GT.ABS(NGUESS)) THEN ! jmc
            WRITE(*,'(A)') ' Too many ts guesses - try neb'
            IGUESS=1
            WRITE(*,'(A,I4)') 'IGUESS=',IGUESS
            TRYNEB=.TRUE.
         ENDIF
      ELSE
         WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
         IF (NIMAGE+4.GT.NIMAGELIMIT) THEN
            WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
            CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,
     1                 REVERSET,STOPFIRST)
            NIMAGE=NIMAGESAVE
            TRYNEB=.FALSE.
         ELSE
            NIMAGE=NIMAGE+4
         ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE TRYAGAINRESET(NATOMS,GUESSTST,TRYNEB,IGUESS,NGUESS,NIMAGE,NIMAGELIMIT,NIMAGESAVE,
     1                   JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,
     2                   REVERSET,STOPFIRST,Q,QSTART,GT20)
      IMPLICIT NONE

      INTEGER NATOMS,IGUESS,NGUESS,NIMAGE,NIMAGELIMIT,NIMAGESAVE,JREMOVE,NMIN,NOPT,NMAXMIN,NMAXTS,J3
      INTEGER WHICHTS(NMAXMIN)
      DOUBLE PRECISION QSAVE(3*NATOMS,NMAXMIN),EMIN(NMAXMIN),Q(3*NATOMS),QSTART(3*NATOMS)
      LOGICAL GUESSTST,TRYNEB,REVERSET,STOPFIRST,CON(NMAXMIN),TSUSED(NMAXTS),GT20
 
      GT20=.FALSE.
      IF (GUESSTST.AND.(.NOT.TRYNEB)) THEN
         IGUESS=IGUESS+1
         WRITE(*,'(A)') ' Try again with next twist'
         IF (IGUESS.GT.ABS(NGUESS)) THEN ! jmc
            WRITE(*,'(A)') ' Too many ts guesses - try neb'
            IGUESS=1
            WRITE(*,'(A,I4)') 'IGUESS=',IGUESS
            TRYNEB=.TRUE.
         ENDIF
C
C  We have to reset the starting geometry in Q!
C
         DO J3=1,NOPT
            Q(J3)=QSTART(J3)
         ENDDO
         GT20=.TRUE.
      ELSE
         WRITE(*,'(A,I4,A)') ' Try again with ',NIMAGE+4,' images'
         IF (NIMAGE+4.GT.NIMAGELIMIT) THEN
            WRITE(*,'(A,I5)') ' Too many neb images requested - remove minimum ',JREMOVE
            CALL REMMIN(JREMOVE,NMIN,NOPT,QSAVE,CON,EMIN,WHICHTS,NMAXMIN,TSUSED,NMAXTS,
     1                 REVERSET,STOPFIRST)
            NIMAGE=NIMAGESAVE
            TRYNEB=.FALSE.
         ELSE
            NIMAGE=NIMAGE+4
C
C  We have to reset the starting geometry in Q!
C
            DO J3=1,NOPT
               Q(J3)=QSTART(J3)
            ENDDO
            GT20=.TRUE.
         ENDIF
      ENDIF

      RETURN
      END
