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
C  This subroutine analyses a min.data.info file in the format generated 
C  the OPTIM keyword BHINTERP or BISECTT.
C
      SUBROUTINE GETALLMIN(MINS,MINF)
      USE PORFUNCS
      USE COMMONS
      USE UTILS,ONLY : GETUNIT
      IMPLICIT NONE

      INTEGER J1, J2, J3, MIN1, MIN2, J4, LUNIT
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS), NEWEMIN, DISTANCE, RMAT(3,3), LOCALPOINTS2(3*NATOMS)
      DOUBLE PRECISION DIST2 
      DOUBLE PRECISION NEWFVIBMIN, NEWPOINTSMIN(3*NATOMS), NEWIXMIN,  NEWIYMIN, NEWIZMIN, TSEGUESS, NEWMINCURVE, NEWMINFRQ2,
     &                 TSFVIBGUESS, ETSDUMMY, DUMMY, FRICTIONFAC
      INTEGER NEWHORDERMIN, NEWMIN, MINS, MINF, NMINSAVE
      LOGICAL MINISOLD, MATCHED, REWRITE

      LUNIT=GETUNIT()
      IF (MACHINE) THEN
         OPEN(LUNIT,FILE='min.data.info',STATUS='OLD',form='unformatted')
      ELSE
         OPEN(LUNIT,FILE='min.data.info',STATUS='OLD')
      ENDIF
      NMINSAVE=NMIN

      J1=0
      DO 
         J1=J1+1

         IF (DUMMYTST.AND.LOWESTFRQT) THEN
            READ(LUNIT,*,END=110) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN,NEWMINCURVE,NEWMINFRQ2
         ELSE
            READ(LUNIT,*,END=110) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN
         ENDIF
         READ(LUNIT,*) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)  
         MINISOLD=.TRUE.
         DO J2=1,NMIN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)  
               CALL MINPERMDIST(NEWPOINTSMIN,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2, 
     &                          RIGIDBODY,RMAT,.FALSE.)
            ENDIF

!           PRINT '(A,I8,5G15.5)',' getallmin> J2,NEWEMIN,EMIN(J2),EDIFFTOL.DISTANCE,GEOMDIFFTOL=',
!    &                              J2,NEWEMIN,EMIN(J2),EDIFFTOL,DISTANCE,GEOMDIFFTOL
            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               NEWMIN=J2
               IF (DEBUG) PRINT '(2(A,I6))','getallmin> minimum ',2*(J1-1)+1,' is database minimum ',J2
               IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getallmin> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
               ENDIF
               IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallmin> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                  NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                  WRITE(*,'(A,I6)') 'getallmin> using maximum value: ',NEWHORDERMIN
               ENDIF
               GOTO 130
            ENDIF
         ENDDO
         MINISOLD=.FALSE.
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         NEWMIN=NMIN
         EMIN(NMIN)=NEWEMIN
         FVIBMIN(NMIN)=NEWFVIBMIN
         HORDERMIN(NMIN)=NEWHORDERMIN
         IXMIN(NMIN)=NEWIXMIN
         IYMIN(NMIN)=NEWIYMIN
         IZMIN(NMIN)=NEWIZMIN
         IF (DUMMYTST.AND.LOWESTFRQT) MINCURVE(NMIN)=NEWMINCURVE
         IF (DUMMYTST.AND.LOWESTFRQT) MINFRQ2(NMIN)=NEWMINFRQ2
         GPFOLD(NMIN)=0.0D0
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

         IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallmin> new minimum ',NMIN,
     &                ' writing parameters to file min.data and points to points.min'
         TOPPOINTER(NMIN)=0
         IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
         IF (DUMMYTST.AND.LOWESTFRQT) THEN
            WRITE(UMINDATA,'(2F25.15,I6,5F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN), MINCURVE(NMIN), MINFRQ2(NMIN)
         ELSE
            WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
         ENDIF
         CALL FLUSH(UMINDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
         WRITE(UMIN,REC=NMIN) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)
!        PRINT '(A,I8)','writing these coords to record ',NMIN
!        PRINT '(3G20.10)',(NEWPOINTSMIN(J2),J2=1,3*NATOMS)
C
C  Set partition functions and transition states for new minimum if DUMMYTST.AND.BHINTERPT.
C
         IF (DUMMYTST.AND.BHINTERPT) THEN
            MINDISTMIN(NMIN)=HUGE(1.0D0)
            DO J3=1,NMIN-1
               READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                          DIST2,RIGIDBODY,RMAT,.FALSE.)
               IF ((DISTANCE.LT.MINDISTMIN(NMIN)).OR.(DISTANCE.LT.MINDISTMIN(J3))) THEN
                  IF (DISTANCE.LT.MINDISTMIN(NMIN)) MINDISTMIN(NMIN)=DISTANCE
                  IF (DISTANCE.LT.MINDISTMIN(J3)) MINDISTMIN(J3)=DISTANCE
C
C  Must create an entry in ts.data in this case.
C  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
C  
                  IF (IMFRQT) THEN
                     PRINT '(A)',"getallmin> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                     STOP
                  ENDIF
                  IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
                  WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') TSEGUESS(EMIN(NMIN),EMIN(J3),MINCURVE(NMIN),MINCURVE(J3),DISTANCE), 
     &                           TSFVIBGUESS(EMIN(NMIN),EMIN(J3),FVIBMIN(NMIN),FVIBMIN(J3),MINFRQ2(NMIN),MINFRQ2(J3),NATOMS),
     &                           1,J3,NMIN,1.0D0,1.0D0,1.0D0
                  CALL FLUSH(UTSDATA)
                  IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
                  NTS=NTS+1
                  IF (NTS.GT.MAXTS) CALL TSDOUBLE
                  ETS(NTS)=TSEGUESS(EMIN(NMIN),EMIN(J3),MINCURVE(NMIN),MINCURVE(J3),DISTANCE)
                  FVIBTS(NTS)=TSFVIBGUESS(EMIN(NMIN),EMIN(J3),FVIBMIN(NMIN),FVIBMIN(J3),MINFRQ2(NMIN),MINFRQ2(J3),NATOMS)
                  HORDERTS(NTS)=1
                  IXTS(NTS)=1.0D0
                  IYTS(NTS)=1.0D0
                  IZTS(NTS)=1.0D0
                  PLUS(NTS)=J3
                  MINUS(NTS)=NMIN
                  IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
                  IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallmin> dummy ts ',NTS,' writing parameters to file ts.data'
C
C  Update ts pointers.
C  
                  POINTERP(NTS)=-1
                  POINTERM(NTS)=-1
                  IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
                  IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
                  TOPPOINTER(PLUS(NTS))=NTS
                  TOPPOINTER(MINUS(NTS))=NTS
C
C  Fill in KPLUS and KMINUS
C
                  IF (ENSEMBLE.EQ.'T') THEN
                     KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2.0D0 * PI*HORDERTS(NTS))) +
     &                  (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(PLUS(NTS)))/TEMPERATURE
                     IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))

                     KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2.0D0 * PI*HORDERTS(NTS))) +
     &                  (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(MINUS(NTS)))/TEMPERATURE
                     IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
                  ELSE
                     IF (TEMPERATURE.GT.ETS(NTS)) THEN
                        KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2*PI*HORDERTS(NTS))) +
     &                   (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(PLUS(NTS))))
                        KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     &                  (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(MINUS(NTS))))
                     ELSE
                        KPLUS(NTS)=-1.0D250
                        KMINUS(NTS)=-1.0D250
                     ENDIF
                  ENDIF
                  IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
                  IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
                  IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)
                  IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)
               ENDIF
            ENDDO
         ENDIF

         IF (DIJINITT) THEN
            CALL GETMETRIC(NMIN,NMIN)

            PRINT '(A,I8)','setup> Finished pair distance calculation for new minimum ',NMIN
          
            IF (INITIALDIST) THEN
               PRINT '(10G13.2)',ALLPAIRS(((NMIN-1)*(NMIN-2))/2+1:(NMIN*(NMIN-1))/2)  
            ELSE
               PRINT '(10G13.2)',PAIRDIST(NMIN,1:PAIRDISTMAX)
               PRINT '(10I13)',PAIRLIST(NMIN,1:PAIRDISTMAX)
            ENDIF

         ENDIF

         CALL FLUSH(UMIN)
130      CONTINUE
      ENDDO
110   CLOSE(LUNIT)
      IF (DIJINITT) THEN
         IF (INITIALDIST) THEN
!
! Append distances to file allpairs. Do the whole file for the moment.
!
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='allpairs',STATUS='UNKNOWN')
            WRITE(LUNIT,'(G20.10)') ALLPAIRS(1:(NMIN*(NMIN-1))/2)
            CLOSE(LUNIT)
         ELSE
            IF (NMIN.GT.NMINSAVE) THEN ! write new pairdist and pairlist files
               LUNIT=GETUNIT()
               OPEN(UNIT=LUNIT,FILE='pairdist',STATUS='UNKNOWN',POSITION='APPEND')
               DO J3=NMINSAVE+1,NMIN
                  WRITE(LUNIT,'(10G20.10)') (PAIRDIST(J3,J4),J4=1,PAIRDISTMAX)
               ENDDO
               CALL FLUSH(LUNIT)
               CLOSE(LUNIT)
               OPEN(UNIT=LUNIT,FILE='pairlist',STATUS='UNKNOWN',POSITION='APPEND')
               DO J3=NMINSAVE+1,NMIN
                  WRITE(LUNIT,'(10I10)') (PAIRLIST(J3,J4),J4=1,PAIRDISTMAX)
               ENDDO
               CALL FLUSH(LUNIT)
               CLOSE(LUNIT)
            ENDIF
         ENDIF
      ENDIF
!
!  If we found new minima between the original ones then adjust the barrier for MINS and MINF
!  according to the sum of distances.
!
      IF ((NMIN-NMINSAVE.GT.0).AND.(BISECTT)) THEN
         MATCHED=.FALSE.
         DO J3=1,NTS
            IF (((PLUS(J3).EQ.MINS).AND.(MINUS(J3).EQ.MINF)).OR.((PLUS(J3).EQ.MINF).AND.(MINUS(J3).EQ.MINS))) THEN
               MATCHED=.TRUE.
               EXIT
            ENDIF
         ENDDO
         IF (.NOT.MATCHED) THEN
            PRINT '(A)','getallmin> no previous dummy ts links the end minima'
         ENDIF
!
!  Add new ts.data items according to the order of the new minima. Assume connections in this order!
!
         READ(UMIN,REC=MINS) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)
         DUMMY=0.0D0
         DO J3=NMINSAVE+1,NMIN+1
            IF (J3.EQ.NMIN+1) THEN
               READ(UMIN,REC=MINF) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ELSE
               READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDIF
            CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                       RMAT,.FALSE.)
            DUMMY=DUMMY+DISTANCE
            IF (J3.EQ.NMINSAVE+1) THEN
               PRINT '(A,I6,A,I6,A,F15.5,A,F15.5)','getallmin> minima ',MINS,' and ',J3,' distance=',DISTANCE,' sum=',DUMMY
               MIN1=MINS; MIN2=J3
            ELSEIF (J3.LE.NMIN) THEN
               PRINT '(A,I6,A,I6,A,F15.5,A,F15.5)','getallmin> minima ',J3-1,' and ',J3,' distance=',DISTANCE,' sum=',DUMMY
               MIN1=J3-1; MIN2=J3
            ELSE
               PRINT '(A,I6,A,I6,A,F15.5,A,F15.5)','getallmin> minima ',J3-1,' and ',MINF,' distance=',DISTANCE,' sum=',DUMMY
               MIN1=J3-1; MIN2=MINF
            ENDIF
!
!  create an entry in ts.data
!  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
!

            IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
            WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') TSEGUESS(EMIN(MIN1),EMIN(MIN2),MINCURVE(MIN1),MINCURVE(MIN2),DISTANCE),
     &                     TSFVIBGUESS(EMIN(MIN1),EMIN(MIN2),FVIBMIN(MIN1),FVIBMIN(MIN2),MINFRQ2(MIN1),MINFRQ2(MIN2),NATOMS),
     &                     1,MIN1,MIN2,1.0D0,1.0D0,1.0D0
            CALL FLUSH(UTSDATA)
            IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
            NTS=NTS+1
            IF (NTS.GT.MAXTS) CALL TSDOUBLE
            ETS(NTS)=TSEGUESS(EMIN(MIN1),EMIN(MIN2),MINCURVE(MIN1),MINCURVE(MIN2),DISTANCE)
            FVIBTS(NTS)=TSFVIBGUESS(EMIN(MIN1),EMIN(MIN2),FVIBMIN(MIN1),FVIBMIN(MIN2),MINFRQ2(MIN1),MINFRQ2(MIN2),NATOMS)
            HORDERTS(NTS)=1
            IXTS(NTS)=1.0D0
            IYTS(NTS)=1.0D0
            IZTS(NTS)=1.0D0
            PLUS(NTS)=MIN1
            MINUS(NTS)=MIN2
            IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
            IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallmin> dummy ts ',NTS,' writing parameters to file ts.data'
!
!  Update ts pointers.
!
            POINTERP(NTS)=-1
            POINTERM(NTS)=-1
            IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
            IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
            TOPPOINTER(PLUS(NTS))=NTS
            TOPPOINTER(MINUS(NTS))=NTS
!
!  Fill in KPLUS and KMINUS
!
            IF (ENSEMBLE.EQ.'T') THEN
               KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2.0D0 * PI*HORDERTS(NTS))) +
     &            (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(PLUS(NTS)))/TEMPERATURE
               IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))

               KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2.0D0 * PI*HORDERTS(NTS))) +
     &            (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(MINUS(NTS)))/TEMPERATURE
               IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
            ELSE
               IF (TEMPERATURE.GT.ETS(NTS)) THEN
                  KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2*PI*HORDERTS(NTS))) +
     &             (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(PLUS(NTS))))
                  KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     &            (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(MINUS(NTS))))
               ELSE
                  KPLUS(NTS)=-1.0D250
                  KMINUS(NTS)=-1.0D250
               ENDIF
            ENDIF
            IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
            IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)
            IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)

            NEWPOINTSMIN(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
         ENDDO

         ETSDUMMY=TSEGUESS(EMIN(MINS),EMIN(MINF),MINCURVE(MINS),MINCURVE(MINF),DUMMY)
         DO J3=1,NTS
            IF (((PLUS(J3).EQ.MINS).AND.(MINUS(J3).EQ.MINF)).OR.((PLUS(J3).EQ.MINF).AND.(MINUS(J3).EQ.MINS))) THEN
               PRINT '(A,F15.5,A,F15.5)','getallmin> new guess for ts energy=',ETSDUMMY,' compared to previous ',ETS(J3)
               IF (ETSDUMMY.GT.ETS(J3)) THEN
                  ETS(J3)=ETSDUMMY
!
!  Need to change the rate constants as well
!
                  IF (ENSEMBLE.EQ.'T') THEN
                     KPLUS(J3)  = LOG(1.0D0 * HORDERMIN(PLUS(J3))  / (2.0D0 * PI*HORDERTS(J3))) +
     1             (FVIBMIN(PLUS(J3))  - FVIBTS(J3)) / 2.0D0 - (ETS(J3) - EMIN(PLUS(J3)) )/TEMPERATURE
                     IF (FRICTIONT) KPLUS(J3)=KPLUS(J3)+LOG(FRICTIONFAC(NEGEIG(J3)))
                     KMINUS(J3) = LOG(1.0D0 * HORDERMIN(MINUS(J3)) / (2.0D0 * PI*HORDERTS(J3))) +
     1             (FVIBMIN(MINUS(J3)) - FVIBTS(J3)) / 2.0D0 - (ETS(J3) - EMIN(MINUS(J3)))/TEMPERATURE
                     IF (FRICTIONT) KMINUS(J3)=KMINUS(J3)+LOG(FRICTIONFAC(NEGEIG(J3)))
                     IF (ZSYM(1:2).EQ.'CA') KPLUS(J3)=KPLUS(J3)+30.66356D0
                     IF (ZSYM(1:2).EQ.'CA') KMINUS(J3)=KMINUS(J3)+30.66356D0
                     IF (PLUS(J3).EQ.MINUS(J3)) KPLUS(J3)=KPLUS(J3)+LOG(2.0D0)
                     IF (PLUS(J3).EQ.MINUS(J3)) KMINUS(J3)=KMINUS(J3)+LOG(2.0D0)
                  ELSE
                     IF (TOTALE.GT.ETS(J3)) THEN
                        KPLUS(J3)  = LOG(1.0D0 * HORDERMIN(PLUS(J3))  / (2*PI*HORDERTS(J3))) +
     1                (FVIBMIN(PLUS(J3))  - FVIBTS(J3))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J3))/(TOTALE-EMIN(PLUS(J3))))
                        KMINUS(J3) = LOG(1.0D0 * HORDERMIN(MINUS(J3)) / (2*PI*HORDERTS(J3))) +
     1                (FVIBMIN(MINUS(J3)) - FVIBTS(J3))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J3))/(TOTALE-EMIN(MINUS(J3))))
                        IF (ZSYM(1:2).EQ.'CA') KPLUS(J3)=KPLUS(J3)+30.66356D0
                        IF (ZSYM(1:2).EQ.'CA') KMINUS(J3)=KMINUS(J3)+30.66356D0
                        IF (PLUS(J3).EQ.MINUS(J3)) KPLUS(J3)=KPLUS(J3)+LOG(2.0D0)
                        IF (PLUS(J3).EQ.MINUS(J3)) KMINUS(J3)=KMINUS(J3)+LOG(2.0D0)
                     ELSE
                        KPLUS(J3)=-1.0D250
                        KMINUS(J3)=-1.0D250
                     ENDIF
                  ENDIF
!
!  KSUM is assumed to be recalculated as needed.
!
                  REWRITE=.TRUE.
               ENDIF
            ENDIF
         ENDDO
         IF (REWRITE) THEN
            PRINT '(A)','getallmin> rewriting ts.data with new barrier data'
            IF (.NOT.CLOSEFILEST) CLOSE(UTSDATA)
            OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='REPLACE')
            DO J3=1,NTS
               WRITE(UTSDATA,'(2F25.15,3I10,4F20.10)') ETS(J3),FVIBTS(J3),1,PLUS(J3),MINUS(J3),1.0,1.0,1.0,1.0
            ENDDO
            CALL FLUSH(UTSDATA)
            IF (CLOSEFILEST) CLOSE(UTSDATA)
         ENDIF
      ENDIF
140   CONTINUE

      RETURN
      END
