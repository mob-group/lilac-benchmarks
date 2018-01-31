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

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  This subroutine analyses a path.info in the new min-sad-min min-sad-min format
C  as generated with OPTIM keyword DUMPALLPATHS.
C
      SUBROUTINE GETALLPATHS
      USE PORFUNCS
      USE COMMONS
      USE UTILS,ONLY : GETUNIT
      IMPLICIT NONE

      INTEGER J1, J2, NMINOLD, TSNUMBER, J3, NCOUNT, NMINSAVE, NTSSAVE, J4, J5, PLUSMIN, LUNIT, JM, JN, NPOSITION
      DOUBLE PRECISION LOCALPOINTS(NOPT), ENERGY, NEWEMIN, NEWETS, DISTANCE, RMAT(3,3),
     1                 LPOINTSTS(NOPT), LPLUS(NOPT), LMINUS(NOPT), LOCALPOINTS2(NOPT)
      DOUBLE PRECISION DUMMY, DIST2, ELAPSED, TNEW
      DOUBLE PRECISION NEWFVIBMIN, NEWFVIBTS, NEWNEGEIG, NEWPOINTSMIN(NOPT), NEWPOINTSMINPLUS(NOPT), EPLUS,
     1                 NEWPOINTSTS(NOPT), NEWIXMIN,  NEWIYMIN, NEWIZMIN, IXPLUS, IYPLUS, IZPLUS,
     2                 NEWIXTS,  NEWIYTS, NEWIZTS, IXMINUS, IYMINUS, IZMINUS, FRICTIONFAC, TEMPD(PAIRDISTMAX)
      INTEGER NEWHORDERMIN, NEWHORDERTS, NEWMIN, NEWTS, NTRIPLES, TEMPL(PAIRDISTMAX), NFRQS, NVARS
      LOGICAL TSISOLD, FAILED, MINPOLD, MINMOLD, BADTRIPLE, TESTOP, TESTNAME
      CHARACTER(LEN=1) DUMMYSTRING

      NMINOLD=NMIN
      IF (MACHINE) THEN
           OPEN(1,FILE='path.info',STATUS='OLD',form='unformatted')
      ELSE
           OPEN(1,FILE='path.info',STATUS='OLD')
      ENDIF
      NFRQS=3*(NATOMS-NGLY)
      NVARS=NOPT
      IF (PHI4MODT) NFRQS=NATOMS
      IF (PHI4MODT) NVARS=NATOMS
      IF (MLP3T.OR.MLPB3T) NFRQS=NATOMS
      IF (MLP3T.OR.MLPB3T) NVARS=NATOMS
C
C  Nasty things can happen if an OPTIM job crashes and leaves a path.info file incomplete.
C  A transition state could be written without the necessary minima. 
C  To avoid this problem, we now parse the path.info file first and only read in complete triples.
C
C  To avoid non-Morse points and higher index saddles getting into the database we now skip
C  a triple where any stationary point has an inappropriate normal mode frequency.
C  We must therefore avoid writing to the open min.data and ts.data files until we
C  have parsed the data!
C
      NCOUNT=0
      DO
         READ(1,*,END=123) DUMMYSTRING
         NCOUNT=NCOUNT+1
      ENDDO
123   REWIND(1)
      IF (NOFRQS) THEN
         NTRIPLES=NCOUNT/(3*(NATOMS+2))
         J1=NTRIPLES*3*(NATOMS+2)
         IF (DEBUG) PRINT '(2(A,I8))','getallpaths> number of triples=',NTRIPLES,' number of trailing lines=',NCOUNT-J1
      ELSEIF (PHI4MODT.OR.MLP3T.OR.MLPB3T) THEN
         IF (MOD(NATOMS,3).EQ.0) THEN
            NTRIPLES=NCOUNT/(3*(2*(NATOMS/3)+2))
            J1=NTRIPLES*(3*(2*(NATOMS/3)+2))
         ELSE
            NTRIPLES=NCOUNT/(3*(2*(NATOMS/3+1)+2))
            J1=NTRIPLES*(3*(2*(NATOMS/3+1)+2))
         ENDIF
         IF (DEBUG) PRINT '(2(A,I8))','getallpaths> number of triples=',NTRIPLES,' number of trailing lines=',NCOUNT-J1
! hk286
      ELSEIF ( (.NOT. NOFRQS) .AND. RIGIDINIT) THEN
         NTRIPLES=NCOUNT/(NOPT+DEGFREEDOMS+6)
         J1=NTRIPLES*(NOPT+DEGFREEDOMS+6)
         IF (DEBUG) PRINT '(2(A,I8))','getallpaths> number of triples=',NTRIPLES,' number of trailing lines=',NCOUNT-J1
      ELSE
         NTRIPLES=NCOUNT/(3*(2*NATOMS-NGLY+2))
         J1=NTRIPLES*3*(2*NATOMS-NGLY+2)
         IF (DEBUG) PRINT '(2(A,I8))','getallpaths> number of triples=',NTRIPLES,' number of trailing lines=',NCOUNT-J1
      ENDIF

      TSISOLD=.TRUE.
      DO J1=1,NTRIPLES 
         IF (DEBUG) PRINT '(A,I6,A,2I10)','getallpaths> doing triple number ',J1,' number of minima and ts=',NMIN,NTS
         IF (DEBUG) CALL FLUSH(6)
         BADTRIPLE=.FALSE.
C
C  NMIN and NTS can be incremented locally within the loop, but are reset to
C  NMINSAVE and NTSSAVE if we diagnose a bad triple due to bad frequencies or the
C  threshold test on the transition state fails.
C
         NMINSAVE=NMIN
         NTSSAVE=NTS
C
C  Read data for plus minimum.
C
         IF (MACHINE) THEN
            READ(1,END=110) ENERGY
            READ(1) HORDER
         ELSE
            READ(1,*,END=110) ENERGY, HORDER
         ENDIF
         NEWEMIN=ENERGY
         IF (NOFRQS) THEN
!            DUMMY=4.675754133D0 ! 2 ln(2pi) +1
            DUMMY = 1.D0
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,NFRQS)
            ELSE
               IF (RIGIDINIT) THEN ! hk286
                  READ(1,*) (FRQS(J2),J2=1,DEGFREEDOMS)
                  FRQS(DEGFREEDOMS+1:NOPT) = 1.0D0
               ELSE
                  READ(1,*) (FRQS(J2),J2=1,NFRQS)
               ENDIF
            ENDIF
            DUMMY=0.0D0
            DO J2=NFSTART,NFFINISH
               IF (FRQS(J2).LE.EVCUT) THEN
                  PRINT '(A,I8,A,G20.10)','getallpaths> SKIPPING - vibrational frequency ',J2,' of + minimum is ',FRQS(J2)
                  BADTRIPLE=.TRUE.
               ELSE
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ENDIF
            ENDDO
         ENDIF
         NEWHORDERMIN=HORDER
         NEWFVIBMIN=DUMMY

         IF (MACHINE) THEN
            READ(1) (NEWPOINTSMIN(J2),J2=1,NVARS)  
         ELSE
            READ(1,*) (NEWPOINTSMIN(J2),J2=1,NVARS)  
         ENDIF
         LOCALPOINTS(1:NVARS)=NEWPOINTSMIN(1:NVARS)
         CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,NEWIXMIN,NEWIYMIN,NEWIZMIN)
         MINPOLD=.TRUE.
         DO J2=1,NMIN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NVARS)  
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF

            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               NEWMIN=J2
               IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path minimum ',2*(J1-1)+1,' is database minimum ',J2
               IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-3) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
               ENDIF
               IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                  NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                  WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERMIN
               ENDIF
               GOTO 130
            ENDIF
         ENDDO
         MINPOLD=.FALSE.
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         NEWMIN=NMIN
         EMIN(NMIN)=NEWEMIN
         FVIBMIN(NMIN)=NEWFVIBMIN
         HORDERMIN(NMIN)=NEWHORDERMIN
         IF (ENSEMBLE.EQ.'T') THEN
            PFMIN(NMIN) = -EMIN(NMIN)/TEMPERATURE - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
         ELSEIF (ENSEMBLE.EQ.'E') THEN
            IF (TOTALE.GT.EMIN(NMIN)) THEN
               PFMIN(NMIN) = (KAPPA-1)*LOG(TOTALE-EMIN(NMIN)) - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
            ELSE
               PFMIN(NMIN) = -1.0D250
            ENDIF
         ENDIF
         PFMIN(NMIN)=PFMIN(NMIN)-PFMEAN
         IXMIN(NMIN)=NEWIXMIN
         IYMIN(NMIN)=NEWIYMIN
         IZMIN(NMIN)=NEWIZMIN
         GPFOLD(NMIN)=0.0D0
         IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallpaths> new minimum ',NMIN
         TOPPOINTER(NMIN)=0
!        KSUM(NMIN)=0.0D0
C
C  We must delay this write until we know that all the stationary points are OK.
C  Writing the points to points.min doesn't matter - this is a direct access file!
C
C        WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
C        CALL FLUSH(UMINDATA)
C

!     INQUIRE(UNIT=UMIN,OPENED=TESTOP,NAMED=TESTNAME)
!     PRINT *,'UNIT UMIN is ',UMIN
!     PRINT *,'opened is ',TESTOP
!     PRINT *,'named is ',TESTNAME

         WRITE(UMIN,REC=NMIN) (NEWPOINTSMIN(J2),J2=1,NVARS)
!        PRINT *,'getallpaths>  calling flush for UMIN'
         CALL FLUSH(UMIN)

130      CONTINUE
         PLUSMIN=NEWMIN ! save in case we need it for a duplicate transition state
         NEWPOINTSMINPLUS(1:NVARS)=NEWPOINTSMIN(1:NVARS)
         EPLUS=NEWEMIN
C
C  Read TS data.
C
         IF (MACHINE) THEN
              READ(1) ENERGY
              READ(1) HORDER
         ELSE
              READ(1,*) ENERGY, HORDER
         ENDIF
         NEWETS=ENERGY
         IF (NOFRQS) THEN
            DUMMY=1.0D0
            NEWNEGEIG=-1.0D0
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,NFRQS)
            ELSE
               IF (RIGIDINIT) THEN ! hk286
                  READ(1,*) (FRQS(J2),J2=1,DEGFREEDOMS)
                  FRQS(DEGFREEDOMS+1:NVARS) = 1.0D0
               ELSE
                  READ(1,*) (FRQS(J2),J2=1,NFRQS)
               ENDIF
            ENDIF
            DUMMY=0.0D0
!
! The current tests do not allow higher index saddles to get into ts.data.
! We also need to test for the case where the negative eigenvalue itself is
! actually below the EVCUT threshold in magnitude. For 2D short-range Morse
! this condition happens a lot!
!
            IF (ABS(FRQS(NFRQS)).LT.EVCUT .AND. (.NOT.RIGIDINIT)) THEN
               WRITE(*,'(A,G20.10)') 'getallpaths> SKIPPING - negative eigenvalue of this transition state is only ',
     &                                FRQS(NFRQS)
               BADTRIPLE=.TRUE.
            ELSEIF (ABS(FRQS(DEGFREEDOMS)).LT.EVCUT .AND. (RIGIDINIT)) THEN
               WRITE(*,'(A,G20.10)') 'getallpaths> SKIPPING - negative eigenvalue of this transition state is only ',
     &                                FRQS(DEGFREEDOMS)
               BADTRIPLE=.TRUE.
            ELSE 
               DO J2=NFSTART,NFFINISH-1
                  IF (FRQS(J2).LT.EVCUT) THEN
                     WRITE(*,'(A,I6,A,G20.10)') 'getallpaths> SKIPPING - eigenvalue ',J2,
     &                                          ' of this transition state is only ',FRQS(J2)
                     BADTRIPLE=.TRUE.
                     EXIT
                  ELSE
                     DUMMY=DUMMY+LOG(FRQS(J2))
                  ENDIF
               ENDDO
               IF (RIGIDINIT) THEN
                  NEWNEGEIG=FRQS(DEGFREEDOMS)
               ELSE
                  NEWNEGEIG=FRQS(NFRQS)
               ENDIF
            ENDIF
         ENDIF
C
C  Now we store the transition state coordinates.
C
         IF (MACHINE) THEN
            READ(1) (NEWPOINTSTS(J2),J2=1,NVARS)  
         ELSE
            READ(1,*) (NEWPOINTSTS(J2),J2=1,NVARS)
         ENDIF
         LOCALPOINTS(1:NVARS)=NEWPOINTSTS(1:NVARS)
         CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,ANGLEAXIS,NEWIXTS,NEWIYTS,NEWIZTS)
         NEWFVIBTS=DUMMY
         NEWHORDERTS=HORDER

         IF (NEWETS.LT.NEWEMIN) PRINT '(2(A,G20.10),A)', 'getallpaths> WARNING *** New TS (energy = ',NEWETS,') is lower in energy 
     &than connected minimum (energy = ',NEWEMIN,')'
C
C  Now check for new transition states.
C
         TSISOLD=.FALSE.
!        IF ((DIJINITT.OR.DIJINITFLYT).AND.(NEWETS.GT.TSTHRESH)) THEN ! reject this ts
!
! We should also reject if both barriers are > MAXBARRIER, but we don't know the minus energy
! yet. We can allow such transition states into the database, but exclude them from
! pathway analysis.
!
         IF (NEWETS.GT.TSTHRESH) THEN ! reject this ts
            PRINT '(A,I8,A,G20.10,A)','getallpaths> Transition state ',J1,' energy ',NEWETS,' lies above the threshold'
            IF (.NOT.MINPOLD) PRINT '(A,I8,A,G20.10,A)','getallpaths> New connected plus minimum will not be added'
            BADTRIPLE=.TRUE.
            GOTO 120
         ENDIF
         DO J2=1,NTS
            DISTANCE=1.0D100
            IF (ABS(NEWETS-ETS(J2)).LT.EDIFFTOL) THEN
               READ(UTS,REC=J2) (LOCALPOINTS2(J3),J3=1,NVARS)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF
!           PRINT '(A,I8,5G20.10)','J2,NEWETS,ETS(J2),diff,tol,DISTANCE=',J2,NEWETS,ETS(J2),ABS(NEWETS-ETS(J2)),EDIFFTOL,DISTANCE
            IF ((ABS(NEWETS-ETS(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path ts ',J1,' energy matches database ts ',J2
               IF (ABS(NEWFVIBTS-FVIBTS(J2))/FVIBTS(J2).GT.1.0D-3) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBTS=',NEWFVIBTS,' should be ',FVIBTS(J2)
               ENDIF
               IF (NEWHORDERTS.NE.HORDERTS(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERTS=',NEWHORDERTS,' should be ',HORDERTS(J2)
                  NEWHORDERTS=MAX(NEWHORDERTS,HORDERTS(J2))
                  WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERTS
               ENDIF
               TSISOLD=.TRUE.
               TSNUMBER=J2
               GOTO 120
            ENDIF
         ENDDO
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         NEWTS=NTS
         ETS(NTS)=NEWETS
         FVIBTS(NTS)=NEWFVIBTS
         HORDERTS(NTS)=NEWHORDERTS
         IXTS(NTS)=NEWIXTS
         IYTS(NTS)=NEWIYTS
         IZTS(NTS)=NEWIZTS
         NEGEIG(NTS)=NEWNEGEIG
         IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
         IF (DEBUG) THEN
            WRITE(*,'(A,I6,A)') 'getallpaths> new intermediate ts ',NTS
         ENDIF
         PLUS(NTS)=NEWMIN

!     INQUIRE(UNIT=UTS,OPENED=TESTOP,NAMED=TESTNAME)
!     PRINT *,'UNIT UTS is ',UTS
!     PRINT *,'opened is ',TESTOP
!     PRINT *,'named is ',TESTNAME

         WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NVARS)
         CALL FLUSH(UTS)
120      CONTINUE
C
C  Read data for minus minimum.
C
         IF (MACHINE) THEN
              READ(1) ENERGY
              READ(1) HORDER
         ELSE
              READ(1,*) ENERGY, HORDER
         ENDIF
         NEWEMIN=ENERGY
         IF (NOFRQS) THEN
!            DUMMY=4.675754133D0 ! 2 ln(2pi) +1
            DUMMY= 1.D0
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,NFRQS)
            ELSE
               IF (RIGIDINIT) THEN ! hk286
                  READ(1,*) (FRQS(J2),J2=1,DEGFREEDOMS)
                  FRQS(DEGFREEDOMS+1:NOPT) = 1.0D0
               ELSE
                  READ(1,*) (FRQS(J2),J2=1,NFRQS)
               ENDIF
            ENDIF
            DUMMY=0.0D0
            DO J2=NFSTART,NFFINISH
               IF (FRQS(J2).LE.EVCUT) THEN
                  PRINT '(A,I8,A,G20.10)','getallpaths> SKIPPING - vibrational frequency ',J2,' of - minimum is ',FRQS(J2)
                  BADTRIPLE=.TRUE.
               ELSE
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ENDIF
            ENDDO
         ENDIF
         NEWHORDERMIN=HORDER
         NEWFVIBMIN=DUMMY
         IF (NEWETS.LT.NEWEMIN) PRINT '(2(A,G20.10),A)', 'getallpaths> WARNING *** New TS (energy = ',NEWETS,') is lower in energy 
     &than connected minimum (energy = ',NEWEMIN,')'

         IF (MACHINE) THEN
              READ(1) (NEWPOINTSMIN(J2),J2=1,NVARS)  
         ELSE
              READ(1,*) (NEWPOINTSMIN(J2),J2=1,NVARS)  
         ENDIF
         LOCALPOINTS(1:NVARS)=NEWPOINTSMIN(1:NVARS)
         CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,angleAxis,NEWIXMIN,NEWIYMIN,NEWIZMIN)
         MINMOLD=.TRUE.
         DO J2=1,NMIN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NVARS)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF

            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               NEWMIN=J2
               IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path minimum ',2*J1,' is database minimum ',J2
               IF (FVIBMIN(J2).NE.0.0D0) THEN
                  IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-3) THEN
                     WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
                  ENDIF
               ENDIF
               IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                  NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                  WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERMIN
               ENDIF
               GOTO 140
            ENDIF
         ENDDO
         MINMOLD=.FALSE.
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         NEWMIN=NMIN
         EMIN(NMIN)=NEWEMIN
         FVIBMIN(NMIN)=NEWFVIBMIN
         HORDERMIN(NMIN)=NEWHORDERMIN
         IF (ENSEMBLE.EQ.'T') THEN
            PFMIN(NMIN) = -EMIN(NMIN)/TEMPERATURE - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
         ELSEIF (ENSEMBLE.EQ.'E') THEN
            IF (TOTALE.GT.EMIN(NMIN)) THEN
               PFMIN(NMIN) = (KAPPA-1)*LOG(TOTALE-EMIN(NMIN)) - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
            ELSE
               PFMIN(NMIN) = -1.0D250
            ENDIF
         ENDIF
         PFMIN(NMIN)=PFMIN(NMIN)-PFMEAN
         IXMIN(NMIN)=NEWIXMIN
         IYMIN(NMIN)=NEWIYMIN
         IZMIN(NMIN)=NEWIZMIN
         GPFOLD(NMIN)=0.0D0
         IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallpaths> new minimum ',NMIN
         TOPPOINTER(NMIN)=0
!        KSUM(NMIN)=0.0D0
C
C  Delay writing data for new minimum until transition state threshold check is passed.
C  It is OK to write points data to the direct access file.
C
C        WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
C        CALL FLUSH(UMINDATA)
C
         WRITE(UMIN,REC=NMIN) (NEWPOINTSMIN(J2),J2=1,NVARS)
!        PRINT *,'getallpaths>  calling flush for UMIN'
         CALL FLUSH(UMIN)

140      CONTINUE
C
C  Now we know the id of the MINUS minimum we can finish the bookkeeping for this ts.
C
         IF (TSISOLD) THEN
            IF (ALLTST) THEN ! need to check all the ts, as there could be more than one entry already
               LOCALPOINTS(1:NVARS)=NEWPOINTSTS(1:NVARS)
               DO J2=1,NTS
                  DISTANCE=1.0D100
                  IF (ABS(NEWETS-ETS(J2)).LT.EDIFFTOL) THEN
                     READ(UTS,REC=J2) (LOCALPOINTS2(J3),J3=1,NVARS)
                     CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
                  ENDIF
!           PRINT '(A,I8,5G20.10)','J2,NEWETS,ETS(J2),diff,tol,DISTANCE=',J2,NEWETS,ETS(J2),ABS(NEWETS-ETS(J2)),EDIFFTOL,DISTANCE
                  IF ((ABS(NEWETS-ETS(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
                     IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path ts ',J1,' energy matches database ts ',J2
                     IF (ABS(NEWFVIBTS-FVIBTS(J2))/FVIBTS(J2).GT.1.0D-3) THEN
                        WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBTS=',NEWFVIBTS,' should be ',FVIBTS(J2)
                     ENDIF
                     IF (NEWHORDERTS.NE.HORDERTS(J2)) THEN
                        WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERTS=',NEWHORDERTS,' should be ',HORDERTS(J2)
                        NEWHORDERTS=MAX(NEWHORDERTS,HORDERTS(J2))
                        WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERTS
                     ENDIF
                     TSNUMBER=J2
                     IF (ABS(EPLUS-EMIN(PLUS(TSNUMBER))).LT.EDIFFTOL) THEN
                        IF (ABS(NEWEMIN-EMIN(MINUS(TSNUMBER))).LT.EDIFFTOL) THEN
                           CALL INERTIAWRAPPER(NEWPOINTSMINPLUS,NOPT,angleAxis,IXPLUS,IYPLUS,IZPLUS)
                           CALL INERTIAWRAPPER(NEWPOINTSMIN,NOPT,angleAxis,IXMINUS,IYMINUS,IZMINUS)
                           IF (DEBUG) THEN
                              PRINT '(A)','getallpaths> assigning + to + and - to -, energies are:'
                              PRINT '(A,2G25.15)','getallpaths> new: ',EPLUS,NEWEMIN
                              PRINT '(A,2G25.15)','getallpaths> old: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                           ENDIF
                           FAILED=.TRUE.
                           IF (DEBUG) WRITE(*,'(A,I10)') 'getallpaths> path matches ts ',J2
                           WRITE(*,'(A,I10)') 'getallpaths> path matches ts ',J2
                           GOTO 753
                        ENDIF
                     ELSEIF (ABS(EPLUS-EMIN(MINUS(TSNUMBER))).LT.EDIFFTOL) THEN
                        IF (ABS(NEWEMIN-EMIN(PLUS(TSNUMBER))).LT.EDIFFTOL) THEN
                           CALL INERTIAWRAPPER(NEWPOINTSMIN,NOPT,angleAxis,IXPLUS,IYPLUS,IZPLUS)
                           CALL INERTIAWRAPPER(NEWPOINTSMINPLUS,NOPT,angleAxis,IXMINUS,IYMINUS,IZMINUS)
                           IF (DEBUG) THEN
                              PRINT '(A)','getallpaths> assigning + to - and - to +, energies are:'
                              PRINT '(A,2G25.15)','getallpaths> new: ',EPLUS,NEWEMIN
                              PRINT '(A,2G25.15)','getallpaths> old: ',EMIN(MINUS(TSNUMBER)),EMIN(PLUS(TSNUMBER))
                           ENDIF
                           FAILED=.TRUE.
                           IF (DEBUG) WRITE(*,'(A,I10)') 'getallpaths> path matches ts ',J2
                           WRITE(*,'(A,I10)') 'getallpaths> path matches ts ',J2
                           GOTO 753
                        ENDIF
                     ENDIF
                PRINT '(A,I6,A,G20.10)','getallpaths> failed consistency check for old ts ',TSNUMBER,' energy=',ETS(TSNUMBER)
                PRINT '(A,2G20.10,A,I8,A,2G20.10)','getallpaths> +/- energies=',EPLUS,NEWEMIN
                PRINT '(A,I8,A,2G25.15)','getallpaths>  ts ',TSNUMBER,' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                  ENDIF
               ENDDO
               PRINT '(A,I6,A,G20.10)','getallpaths> failed consistency check for all old ts '
               TSISOLD=.FALSE.
               NTS=NTS+1
               IF (NTS.GT.MAXTS) CALL TSDOUBLE
               NEWTS=NTS
               ETS(NTS)=NEWETS
               FVIBTS(NTS)=NEWFVIBTS
               HORDERTS(NTS)=NEWHORDERTS
               IXTS(NTS)=NEWIXTS
               IYTS(NTS)=NEWIYTS
               IZTS(NTS)=NEWIZTS
               NEGEIG(NTS)=NEWNEGEIG
               IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
!              IF (DEBUG) THEN
                  WRITE(*,'(A,I6,A)') 'getallpaths> new intermediate ts recounted for different connections ',NTS
!              ENDIF
               PLUS(NTS)=PLUSMIN
               WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NVARS)
               CALL FLUSH(UTS)
               GOTO 140
753            CONTINUE
            ELSE IF (TSNUMBER.GT.0) THEN
               FAILED=.FALSE.
C
C  Consistency check for minima.
C
               IF (ABS(EPLUS-EMIN(PLUS(TSNUMBER))).LT.EDIFFTOL) THEN
                  IF (ABS(NEWEMIN-EMIN(MINUS(TSNUMBER))).LT.EDIFFTOL) THEN
                     CALL INERTIAWRAPPER(NEWPOINTSMINPLUS,NOPT,angleAxis,IXPLUS,IYPLUS,IZPLUS)
                     CALL INERTIAWRAPPER(NEWPOINTSMIN,NOPT,angleAxis,IXMINUS,IYMINUS,IZMINUS)
                     IF (DEBUG) THEN
                        PRINT '(A)','getallpaths> assigning + to + and - to -, energies are:'
                        PRINT '(A,2G25.15)','getallpaths> new: ',EPLUS,NEWEMIN
                        PRINT '(A,2G25.15)','getallpaths> old: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                     ENDIF
                  ELSE
                     PRINT '(A,I6,A,G20.10)','getallpaths> WARNING failed consistency check for old ts ',
     &                                        TSNUMBER,' energy=',ETS(TSNUMBER)
                     PRINT '(A,2G25.15)','getallpaths> +/- energies=',EPLUS,NEWEMIN
                     PRINT '(A,I8,A,2G25.15)','getallpaths> ts ',TSNUMBER,
     &                                      ' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
!
! Allow multiple ts with different connections if ALLTST is true.
!
                     IF (ALLTST) THEN
                        TSISOLD=.FALSE.
                        NTS=NTS+1
                        IF (NTS.GT.MAXTS) CALL TSDOUBLE
                        NEWTS=NTS
                        ETS(NTS)=NEWETS
                        FVIBTS(NTS)=NEWFVIBTS
                        HORDERTS(NTS)=NEWHORDERTS
                        IXTS(NTS)=NEWIXTS
                        IYTS(NTS)=NEWIYTS
                        IZTS(NTS)=NEWIZTS
                        NEGEIG(NTS)=NEWNEGEIG
                        IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
!                       IF (DEBUG) THEN
                           WRITE(*,'(A,I6,A)') 'getallpaths> new intermediate ts recounted for different connections ',NTS
!                       ENDIF
                        PLUS(NTS)=PLUSMIN
                        WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NVARS)
                        CALL FLUSH(UTS)
                        GOTO 140
                     ELSE
                        FAILED=.TRUE.
!                       IF (DEBUG) WRITE(*,'(A)') 'getallpaths> not recounting this ts for different connections'
                        WRITE(*,'(A)') 'getallpaths> not recounting this ts for different connections'
                     ENDIF
                  ENDIF
               ELSEIF (ABS(EPLUS-EMIN(MINUS(TSNUMBER))).LT.EDIFFTOL) THEN
                  IF (ABS(NEWEMIN-EMIN(PLUS(TSNUMBER))).LT.EDIFFTOL) THEN
                     CALL INERTIAWRAPPER(NEWPOINTSMIN,NOPT,angleAxis,IXPLUS,IYPLUS,IZPLUS)
                     CALL INERTIAWRAPPER(NEWPOINTSMINPLUS,NOPT,angleAxis,IXMINUS,IYMINUS,IZMINUS)
                     IF (DEBUG) THEN
                        PRINT '(A)','getallpaths> assigning + to - and - to +, energies are:'
                        PRINT '(A,2G25.15)','getallpaths> new: ',EPLUS,NEWEMIN
                        PRINT '(A,2G25.15)','getallpaths> old: ',EMIN(MINUS(TSNUMBER)),EMIN(PLUS(TSNUMBER))
                     ENDIF
                  ELSE
                     PRINT '(A,I6,A,G20.10)','getallpaths> WARNING failed consistency check for old ts ',
     &                                        TSNUMBER,' energy=',ETS(TSNUMBER)
                     PRINT '(A,2G20.10,A,I8,A,2G20.10)','getallpaths> +/- energies=',EPLUS,NEWEMIN
                     PRINT '(A,I8,A,2G25.15)','getallpaths>  ts ',TSNUMBER,
     &                                         ' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                     IF (ALLTST) THEN
                        TSISOLD=.FALSE.
                        NTS=NTS+1
                        IF (NTS.GT.MAXTS) CALL TSDOUBLE
                        NEWTS=NTS
                        ETS(NTS)=NEWETS
                        FVIBTS(NTS)=NEWFVIBTS
                        HORDERTS(NTS)=NEWHORDERTS
                        IXTS(NTS)=NEWIXTS
                        IYTS(NTS)=NEWIYTS
                        IZTS(NTS)=NEWIZTS
                        NEGEIG(NTS)=NEWNEGEIG
                        IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
!                       IF (DEBUG) THEN
                           WRITE(*,'(A,I6,A)') 'getallpaths> new intermediate ts recounted for different connections ',NTS
!                       ENDIF
                        PLUS(NTS)=PLUSMIN
                        WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NVARS)
                        CALL FLUSH(UTS)
                        GOTO 140
                     ELSE
                        FAILED=.TRUE.
!                       IF (DEBUG) WRITE(*,'(A)') 'getallpaths> not recounting this ts for different connections'
                        WRITE(*,'(A)') 'getallpaths> not recounting this ts for different connections'
                     ENDIF
                  ENDIF
               ELSE
                  PRINT '(A,I6,A,G20.10)','getallpaths> WARNING failed consistency check for old ts ',
     &                                     TSNUMBER,' energy=',ETS(TSNUMBER)
                  PRINT '(A,2G25.15)','getallpaths> +/- energies=',EPLUS,NEWEMIN
                  PRINT '(A,I8,A,2G25.15)','getallpaths> ts ',TSNUMBER,
     &                                     ' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                  IF (ALLTST) THEN
                     TSISOLD=.FALSE.
                     NTS=NTS+1
                     IF (NTS.GT.MAXTS) CALL TSDOUBLE
                     NEWTS=NTS
                     ETS(NTS)=NEWETS
                     FVIBTS(NTS)=NEWFVIBTS
                     HORDERTS(NTS)=NEWHORDERTS
                     IXTS(NTS)=NEWIXTS
                     IYTS(NTS)=NEWIYTS
                     IZTS(NTS)=NEWIZTS
                     NEGEIG(NTS)=NEWNEGEIG
                     IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
!                    IF (DEBUG) THEN
                        WRITE(*,'(A,I6,A)') 'getallpaths> new intermediate ts recounted for different connections ',NTS
!                    ENDIF
                     PLUS(NTS)=PLUSMIN
                     WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NVARS)
                     CALL FLUSH(UTS)
                     GOTO 140
                  ELSE
                     FAILED=.TRUE.
!                    IF (DEBUG) WRITE(*,'(A)') 'getallpaths> not recounting this ts for different connections'
                     WRITE(*,'(A)') 'getallpaths> not recounting this ts for different connections'
                  ENDIF
               ENDIF
               IF ((.NOT.FAILED).AND.(.NOT.BULKT)) THEN
                  IF (ABS(IXPLUS-IXMIN(PLUS(TSNUMBER))).GT.IDIFFTOL) THEN
                     PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia +x: ',
     &                                    IXPLUS,IXMIN(PLUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IYPLUS-IYMIN(PLUS(TSNUMBER))).GT.IDIFFTOL) THEN
                     PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia +y: ',
     &                                    IYPLUS,IYMIN(PLUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IZPLUS-IZMIN(PLUS(TSNUMBER))).GT.IDIFFTOL) THEN
                     PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia +z: ',
     &                                    IZPLUS,IZMIN(PLUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IXMINUS-IXMIN(MINUS(TSNUMBER))).GT.IDIFFTOL) THEN
                    PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: -x',
     &                                    IXMINUS,IXMIN(MINUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IYMINUS-IYMIN(MINUS(TSNUMBER))).GT.IDIFFTOL) THEN
                    PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: -y',
     &                                    IYMINUS,IYMIN(MINUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IZMINUS-IZMIN(MINUS(TSNUMBER))).GT.IDIFFTOL) THEN
                    PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: -z',
     &                                    IZMINUS,IZMIN(MINUS(TSNUMBER))
                  ENDIF
               ENDIF
            ENDIF
C
C  Old ts might link new minima. This is inconsistent and is ignored unless ALLTST is true.
C  Any new minima are not written to min.data, so we should reset to the saved NMIN value,
C  just to be safe.
C
C  For small systems we see lines of control characters occasionally written to
C  min.data and ts.data. Probably due to a linux nfs problem. 
C  Addressed this issue by closing and opening files.
C
            NMIN=NMINSAVE
         ELSE IF (.NOT.BADTRIPLE) THEN
            MINUS(NTS)=NEWMIN
            IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
            IF (IMFRQT) THEN
               WRITE(UTSDATA,'(2F25.15,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     &                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
            ELSE
               WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     &                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
            END IF
            CALL FLUSH(UTSDATA)
            IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
            PRINT '(A,I6,A)','getallpaths> writing data for new ts to ts.data'
            IF (NEWCONNECTIONST) THEN
               MINCONN(PLUS(NTS))=MINCONN(PLUS(NTS))+1
               MINCONN(MINUS(NTS))=MINCONN(MINUS(NTS))+1
            ENDIF
C           PRINT '(A,2L5,2I6)','MINPOLD,MINMOLD,NMINSAVE,NMIN=',MINPOLD,MINMOLD,NMINSAVE,NMIN
            IF (NMIN-NMINSAVE.GT.0) THEN
               PRINT '(A,I6,A)','getallpaths> writing data for ',NMIN-NMINSAVE,' new min to min.data'
               IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
               DO J2=NMINSAVE+1,NMIN
                  WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(J2),FVIBMIN(J2),HORDERMIN(J2),IXMIN(J2),IYMIN(J2),IZMIN(J2)
                  CALL FLUSH(UMINDATA)
               ENDDO
               IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
            ENDIF
!
! (1) Construct pairlist and pairdist entries for any new minima.
! (2) Insert zero entry into PAIRDIST for connected pair. 
! Must remove any previous entry for this pair first.
! Put the zero entry at position one, and move the others
! along, unless we hit the same minimum, in which case we
! can overwrite it and leave the rest of the sorted list
! unchanged.
! For INITIALDIST we maintain a full list of distances and not pairlist and pairdist.
! Need to insert zero entry into ALLPAIRS
!
            IF (DIJINITT) THEN
               CALL CPU_TIME(ELAPSED)
               CALL GETMETRIC(NMINSAVE+1,NMIN) ! NMINSAVE+1 is the first new minimum, NMIN is the last.
               CALL CPU_TIME(TNEW)
               IF (INITIALDIST) THEN
                  JM=MIN(PLUS(NTS),MINUS(NTS))
                  JN=MAX(PLUS(NTS),MINUS(NTS))
                  IF (JM.NE.JN) THEN
                     NPOSITION=((JN-2)*(JN-1))/2+JM
                     ALLPAIRS(NPOSITION)=0.0D0
!
! Append distances to file allpairs. Need to rewrite the file if we have found
! a connection that changes connections in the set of minima previously known.
!
                     LUNIT=GETUNIT()
                     IF (JM.LE.NMINSAVE) THEN
                        OPEN(UNIT=LUNIT,FILE='allpairs',STATUS='OLD')
                        WRITE(LUNIT,'(G20.10)') ALLPAIRS(1:(NMIN*(NMIN-1))/2)
                     ELSE
                        OPEN(UNIT=LUNIT,FILE='allpairs',POSITION='APPEND',ACTION='WRITE',STATUS='OLD')
                        WRITE(LUNIT,'(G20.10)') ALLPAIRS((NMINSAVE*(NMINSAVE-1))/2+1:(NMIN*(NMIN-1))/2)  
                     ENDIF
                     CLOSE(LUNIT)
                  ENDIF
               ELSE
                  J2=MAX(PLUS(NTS),MINUS(NTS))
                  J3=MIN(PLUS(NTS),MINUS(NTS))
!                 PRINT '(A,3I6)','getallpaths> J2,J3,NTS=',J2,J3,NTS
                  IF (J2.NE.J3) THEN
!                    PRINT '(A)','pairlist old J2:'
!                    PRINT '(10I10)',PAIRLIST(J2,1:PAIRDISTMAX)
                     TEMPL(1:PAIRDISTMAX)=PAIRLIST(J2,1:PAIRDISTMAX)
                     TEMPD(1:PAIRDISTMAX)=PAIRDIST(J2,1:PAIRDISTMAX)
                     DO J5=2,PAIRDISTMAX
                        IF (PAIRLIST(J2,J5-1).EQ.J3) EXIT
                        TEMPL(J5)=PAIRLIST(J2,J5-1)
                        TEMPD(J5)=PAIRDIST(J2,J5-1)
                     ENDDO
                     PAIRLIST(J2,2:PAIRDISTMAX)=TEMPL(2:PAIRDISTMAX)
                     PAIRDIST(J2,2:PAIRDISTMAX)=TEMPD(2:PAIRDISTMAX)
!                    PRINT '(A)','pairlist new J2:'
!                    PRINT '(10I10)',PAIRLIST(J2,1:PAIRDISTMAX)
                     PAIRDIST(J2,1)=0.0D0
                     PAIRLIST(J2,1)=J3
!                    PRINT '(A)','pairlist J2,1:'
!                    PRINT '(10I10)',PAIRLIST(J2,1)
!                    PRINT '(A)','pairlist old J3:'
!                    PRINT '(10I10)',PAIRLIST(J3,1:PAIRDISTMAX)
                     TEMPL(1:PAIRDISTMAX)=PAIRLIST(J3,1:PAIRDISTMAX)
                     TEMPD(1:PAIRDISTMAX)=PAIRDIST(J3,1:PAIRDISTMAX)
                     DO J5=2,PAIRDISTMAX
                        IF (PAIRLIST(J3,J5-1).EQ.J2) EXIT
                        TEMPL(J5)=PAIRLIST(J3,J5-1)
                        TEMPD(J5)=PAIRDIST(J3,J5-1)
                     ENDDO
                     PAIRLIST(J3,2:PAIRDISTMAX)=TEMPL(2:PAIRDISTMAX)
                     PAIRDIST(J3,2:PAIRDISTMAX)=TEMPD(2:PAIRDISTMAX)
!                    PRINT '(A)','pairlist new J3:'
!                    PRINT '(10I10)',PAIRLIST(J3,1:PAIRDISTMAX)
                     PAIRDIST(J3,1)=0.0D0
                     PAIRLIST(J3,1)=J2
!                    PRINT '(A)','pairlist J3,1:'
!                    PRINT '(10I10)',PAIRLIST(J3,1)
                     IF (DEBUG) THEN
                        PRINT '(A,2I8)','getallpaths> Changed pair distance list for minima ',J2,J3
                        PRINT '(10G13.2)',PAIRDIST(J2,1:PAIRDISTMAX)
                        PRINT '(10I13)',PAIRLIST(J2,1:PAIRDISTMAX)
                        PRINT '(10G13.2)',PAIRDIST(J3,1:PAIRDISTMAX)
                        PRINT '(10I13)',PAIRLIST(J3,1:PAIRDISTMAX)
                     ENDIF
!
! Since entries have changed we'd better rewrite these files rather
! than just append to them. 
!
                     LUNIT=GETUNIT()
                     OPEN(UNIT=LUNIT,FILE='pairdist',STATUS='UNKNOWN')
                     DO J3=1,NMIN
                        WRITE(LUNIT,'(10G20.10)') (PAIRDIST(J3,J4),J4=1,PAIRDISTMAX)
                     ENDDO
                     CALL FLUSH(LUNIT)
                     CLOSE(LUNIT)
                     OPEN(UNIT=LUNIT,FILE='pairlist',STATUS='UNKNOWN')
                     DO J3=1,NMIN
                        WRITE(LUNIT,'(10I10)') (PAIRLIST(J3,J4),J4=1,PAIRDISTMAX)
                     ENDDO
                     CALL FLUSH(LUNIT)
                     CLOSE(LUNIT)
                  ENDIF
               ENDIF
            ENDIF
C
C  Update ts pointers.
C
            POINTERP(NTS)=-1
            POINTERM(NTS)=-1
            IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
            IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
            TOPPOINTER(PLUS(NTS))=NTS
            TOPPOINTER(MINUS(NTS))=NTS

C           PRINT '(A,7I8)','NTS,PLUS(NTS),MINUS(NTS),TOP+,TOP-,PP,PM=',
C    &               NTS,PLUS(NTS),MINUS(NTS),TOPPOINTER(PLUS(NTS)),TOPPOINTER(MINUS(NTS)),POINTERP(NTS),POINTERM(NTS)
   
C
C  Calculate rates.
C
            IF (ENSEMBLE.EQ.'T') THEN
               KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2.0D0 * PI*HORDERTS(NTS))) +
     1          (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(PLUS(NTS)))/TEMPERATURE
               IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))

               KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2.0D0 * PI*HORDERTS(NTS))) +
     1          (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(MINUS(NTS)))/TEMPERATURE
               IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
            ELSE
               IF (TEMPERATURE.GT.ETS(NTS)) THEN
                  KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2*PI*HORDERTS(NTS))) +
     1               (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(PLUS(NTS))))
                  KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1              (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(MINUS(NTS))))
               ELSE
                  KPLUS(NTS)=-1.0D250
                  KMINUS(NTS)=-1.0D250
               ENDIF
            ENDIF
            IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
            IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)
            IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)
C
C  Don't update sum of rates out of the connected minima.
C  Assume KSUM will be calculated when needed.
C
            DO J2=1,NVARS
               LPOINTSTS(J2)=NEWPOINTSTS(J2)
               LPLUS(J2)=NEWPOINTSMINPLUS(J2)
               LMINUS(J2)=NEWPOINTSMIN(J2)
            ENDDO
            IF (ADDPT) CALL ADDPERM(LPOINTSTS,LPLUS,LMINUS) 
            IF (ADDPT2) CALL ADDPERM2(NTS,LPOINTSTS,LPLUS,LMINUS) 
            IF (ADDPT3) CALL ADDPERM3(NTS,LPOINTSTS,LPLUS,LMINUS) 
         ELSE
C
C  Old ts or bad triple encountered. Either way, resetting to saved NTS and NMIN values should be safe.
C
            NTS=NTSSAVE
            NMIN=NMINSAVE
         ENDIF
      ENDDO
110   CLOSE(1)
C
C  Find pathways out of the new intermediate minimum. We can;t do this before
C  because we need all the data for new minima and transition states to be assigned in the
C  right places.
C
C  This has to be turned off to avoid a recursive call to getallpaths, which
C  is now called from tssearch.
C
C     IF ((CONNECTIONS.GT.2).AND.(.NOT.CHECKCONNECTIONST)) THEN
C        WRITE(*,'(A,I6,A)') 'getallpaths> checking for at least ',CONNECTIONS,' connections per minimum'
C        DO J1=NMINOLD+1,NMIN
CC           CALL TSSEARCH2(NEWMIN(J1),0)
C           CALL TSSEARCH(J1,0)
C        ENDDO
C     ENDIF

      RETURN
      END
