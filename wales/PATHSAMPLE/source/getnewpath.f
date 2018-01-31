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
C  This subroutine analyses a new path segment from an OPTIM connect run. 
C
      SUBROUTINE GETNEWPATH(MINS,MINF)
      USE PORFUNCS
      USE COMMONS
      IMPLICIT NONE

      INTEGER J1, J2, ISTAT, MINS, MINF, NMINOLD, J3
      DOUBLE PRECISION LOCALPOINTS(NOPT), ENERGY, NEWEMIN, NEWETS, DISTANCE, RMAT(3,3),
     1                 LPOINTSTS(NOPT), LPLUS(NOPT), LMINUS(NOPT), LOCALPOINTS2(NOPT), DIST2
      DOUBLE PRECISION DUMMY
      DOUBLE PRECISION NEWFVIBMIN, NEWFVIBTS, NEWPOINTSMIN(NOPT), NEWPOINTSMINPREV(NOPT),
     1                 NEWPOINTSTS(NOPT), NEWIXMIN,  NEWIYMIN, NEWIZMIN,
     2                 NEWIXTS,  NEWIYTS, NEWIZTS, NEWNEGEIG, FRICTIONFAC
      INTEGER NEWHORDERMIN, NEWHORDERTS, NEWMIN, NEWTS
      LOGICAL TSISOLD

      NMINOLD=NMIN

      if (machine) then
           OPEN(1,FILE='path.info',STATUS='OLD',form='unformatted')
      else
           OPEN(1,FILE='path.info',STATUS='OLD')
      endif

      TSISOLD=.TRUE.
      J1=0
      DO 
         J1=J1+1
         IF (MACHINE) THEN
              READ(1) ENERGY
              READ(1) HORDER
         ELSE
              READ(1,*) ENERGY, HORDER
         ENDIF
         NEWEMIN=ENERGY
         IF (NOFRQS) THEN
!            DUMMY=4.675754133D0 ! 2 ln(2pi) +1
            DUMMY = 1.D0 
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ELSE
               READ(1,*) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ENDIF
            DUMMY=0.0D0
            DO J2=NFSTART,NFFINISH
               IF (FRQS(J2).LE.0.0D0) THEN
                  PRINT '(A,I8,A,G20.10)','getnewpath> ERROR - vibrational frequency ',J2,' is ',FRQS(J2)
                  STOP
               ENDIF
               DUMMY=DUMMY+LOG(FRQS(J2))
            ENDDO
         ENDIF
         IF ((J1.EQ.1).AND.(MINS.GT.0)) THEN  ! check that frequencies agree
            IF (ABS((DUMMY-FVIBMIN(MINS))/DUMMY).GT.1.0D-4) THEN
               WRITE(*,'(A,F15.5,A,F15.5)') 'getnewpath> WARNING - frequencies of starting minimum ',DUMMY,
     1           ' deviate significantly from database value ',FVIBMIN(MINS)
            ENDIF
            IF (HORDER.NE.HORDERMIN(MINS)) THEN
               WRITE(*,'(A,2I6)') 'getnewpath> ERROR, order of initial minimum does not agree with database: ',
     1                                         HORDER,HORDERMIN(MINS)
               HORDER=HORDERMIN(MINS)
               WRITE(*,'(A,I6)') 'getnewpath> Using the original value of ',HORDER
            ENDIF
         ENDIF
         NEWHORDERMIN=HORDER
         NEWFVIBMIN=DUMMY

         if (machine) then
              READ(1) (NEWPOINTSMIN(J2),J2=1,NOPT)  
         else
              READ(1,*) (NEWPOINTSMIN(J2),J2=1,NOPT)  
         endif
         LOCALPOINTS(1:NOPT)=NEWPOINTSMIN(1:NOPT)
         CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,ANGLEAXIS,NEWIXMIN,NEWIYMIN,NEWIZMIN)
         DO J2=1,NMIN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NOPT)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF

            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               NEWMIN=J2
               IF (DEBUG) PRINT '(2(A,I6))','getnewpath> path minimum ',J1,' is database minimum ',J2
               IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getnewpath> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
               ENDIF
               IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getnewpath> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                  NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                  WRITE(*,'(A,I6)') 'getnewpath> using maximum value: ',NEWHORDERMIN
               ENDIF
               GOTO 130
            ENDIF
         ENDDO
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
         WRITE(*,'(A,I6,A)') 'getnewpath> new minimum ',NMIN,' writing parameters to file min.data and points to points.min'
         TOPPOINTER(NMIN)=0
         IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
         CALL FLUSH(UMINDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
         WRITE(UMIN,REC=NMIN) (NEWPOINTSMIN(J2),J2=1,NOPT)

         CALL FLUSH(UMIN)
130      CONTINUE
C
C  If TSISOLD is .FALSE. the previous ts was new. Now we know the id of the MINUS
C  minimum we can finish the bookkeeping for this ts.
C
         IF (.NOT.TSISOLD) THEN
            MINUS(NTS)=NEWMIN
            IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
            IF (IMFRQT) THEN
               WRITE(UTSDATA,'(2F20.10,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     &                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
            ELSE
               WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     &                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
            ENDIF
            CALL FLUSH(UTSDATA)
            IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
C
C  Update ts pointers.
C
            POINTERP(NTS)=-1
            POINTERM(NTS)=-1
            IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
            IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
            TOPPOINTER(PLUS(NTS))=NTS
            TOPPOINTER(MINUS(NTS))=NTS
   
C            DO J2=NTS-1,1,-1
C               IF (PLUS(J2).EQ.PLUS(NTS)) THEN
C                  POINTERP(NTS)=J2
C                  GOTO 41
C               ELSE IF (MINUS(J2).EQ.PLUS(NTS)) THEN
C                  POINTERP(NTS)=J2
C                  GOTO 41
C               ENDIF
C            ENDDO
C41          CONTINUE
C
C            DO J2=NTS-1,1,-1
C               IF (PLUS(J2).EQ.MINUS(NTS)) THEN
C                  POINTERM(NTS)=J2
C                  GOTO 42
C               ELSE IF (MINUS(J2).EQ.MINUS(NTS)) THEN
C                  POINTERM(NTS)=J2
C                  GOTO 42
C               ENDIF
C            ENDDO
C42          CONTINUE
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
C  Update sum of rates out of the connected minima.
C
!           IF (KSUM(PLUS(NTS)).EQ.0.0D0) THEN
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(PLUS(NTS))=KPLUS(NTS)
!           ELSE
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(PLUS(NTS)) =LOG(EXP(KSUM(PLUS(NTS))-KMEAN) + EXP(KPLUS(NTS) -KMEAN)) + KMEAN
!           ENDIF
!           IF (KSUM(MINUS(NTS)).EQ.0.0D0) THEN
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(MINUS(NTS))=KMINUS(NTS)
!           ELSE
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(MINUS(NTS))=LOG(EXP(KSUM(MINUS(NTS))-KMEAN) + EXP(KMINUS(NTS)-KMEAN)) + KMEAN
!           ENDIF

            DO J2=1,NOPT
               LPOINTSTS(J2)=NEWPOINTSTS(J2)
               LPLUS(J2)=NEWPOINTSMINPREV(J2)
               LMINUS(J2)=NEWPOINTSMIN(J2)
            ENDDO
            IF (ADDPT) CALL ADDPERM(LPOINTSTS,LPLUS,LMINUS) 
            IF (ADDPT2) CALL ADDPERM2(NTS,LPOINTSTS,LPLUS,LMINUS) 
            IF (ADDPT3) CALL ADDPERM3(NTS,LPOINTSTS,LPLUS,LMINUS)
         ENDIF
         NEWPOINTSMINPREV(1:NOPT)=NEWPOINTSMIN(1:NOPT)
C
C  Read TS data - if the above minimum is the last one then there won;t be a TS here
C  and we must exit.
C
         if (machine) then
              READ(1,END=110) ENERGY
              READ(1) HORDER
         else
              READ(1,*,END=110) ENERGY, HORDER
         endif
         NEWETS=ENERGY
         IF (NOFRQS) THEN
            DUMMY=1.0D0
            NEWNEGEIG=-1.0D0
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ELSE
               READ(1,*) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ENDIF
            DUMMY=0.0D0
            DO J2=NFSTART,NFFINISH-1
               IF (FRQS(J2).LT.2.0D-6) THEN 
                  WRITE(*,'(A,I6,A,G20.10)') 'getallpaths> WARNING - eigenvalue ',J2,' of this transition state is only ',FRQS(J2)
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ELSE IF (FRQS(J2).LE.0.0D0) THEN
                  PRINT '(A,I8,A,G20.10)','getnewpath> ERROR - vibrational frequency ',J2,' is ',FRQS(J2)
                  STOP
               ELSE
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ENDIF
            ENDDO
            NEWNEGEIG = FRQS(3*(NATOMS-NGLY))
         ENDIF
C
C  Now we store the transition state coordinates.
C
         IF (MACHINE) THEN
            READ(1) (NEWPOINTSTS(J2),J2=1,NOPT)  
         ELSE
            READ(1,*) (NEWPOINTSTS(J2),J2=1,NOPT)  
         ENDIF
         LOCALPOINTS(1:NOPT)=NEWPOINTSTS(1:NOPT)
         CALL INERTIAWRAPPER(LOCALPOINTS,NOPT,ANGLEAXIS,NEWIXTS,NEWIYTS,NEWIZTS)
         NEWFVIBTS=DUMMY
         NEWHORDERTS=HORDER
C
C  Now check for new transition states.
C
         TSISOLD=.FALSE. 
!
! We should also reject if both barriers are > MAXBARRIER, but we don't know the minus energy
! yet. We can allow such transition states into the database, but exclude them from
! pathway analysis.
!
         IF (DIJINITT.AND.(NEWETS.GT.TSTHRESH)) THEN ! reject this ts
            PRINT '(A,I8,A,G20.10,A)','getnewpath> Transition state ',J1,' energy ',NEWETS,' lies above the threshold'
            TSISOLD=.TRUE.
            GOTO 120
         ENDIF
         DO J2=1,NTS
            DISTANCE=1.0D100
            IF (ABS(NEWETS-ETS(J2)).LT.EDIFFTOL) THEN
               READ(UTS,REC=J2) (LOCALPOINTS2(J3),J3=1,NOPT)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
            ENDIF
            IF ((ABS(NEWETS-ETS(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               IF (DEBUG) PRINT '(2(A,I6))','getnewpath> path ts ',J1,' is database ts ',J2
               IF (ABS(NEWFVIBTS-FVIBTS(J2))/FVIBTS(J2).GT.1.0D-4) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getnewpath> WARNING, NEWFVIBTS=',NEWFVIBTS,' should be ',FVIBTS(J2)
               ENDIF
               IF (NEWHORDERTS.NE.HORDERTS(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getnewpath> ERROR, NEWHORDERTS=',NEWHORDERTS,' should be ',HORDERTS(J2)
                  NEWHORDERTS=MAX(NEWHORDERTS,HORDERTS(J2))
                  WRITE(*,'(A,I6)') 'getnewpath> using maximum value: ',NEWHORDERTS
               ENDIF
               TSISOLD=.TRUE.
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
         WRITE(*,'(A,I6,A)') 'getnewpath> new intermediate ts ',NTS,' writing parameters to file ts.data and updating rates'
         PLUS(NTS)=NEWMIN
         WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NOPT)

         CALL FLUSH(UTS)
120      CONTINUE
      ENDDO
110   CLOSE(1)
C
C  We already know the parameters for the final minimum if MINF>0, so check that they agree!
C  
      IF (MINF.GT.0) THEN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(MINF)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=MINF) (LOCALPOINTS2(J3),J3=1,NOPT)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
            ENDIF
         IF ((ABS(NEWEMIN-EMIN(MINF)).GT.EDIFFTOL).OR.(DISTANCE.GT.GEOMDIFFTOL)) THEN
            WRITE(*,'(A,2F20.10,A,2F20.10)') 'getnewpath> ERROR, energy and DISTANCE=',ENERGY,
     1               DISTANCE,' should be ',EMIN(MINF),0.0D0
            STOP
         ENDIF
         IF (ABS((DUMMY-FVIBMIN(MINF))/DUMMY).GT.1.0D-4) THEN
            WRITE(*,'(A,F15.5,A,F15.5)') 'getnewpath> WARNING - frequencies of final minimum ',DUMMY,
     1                         ' deviate significantly from database value ',FVIBMIN(MINF)
         ENDIF
         IF (HORDER.NE.HORDERMIN(MINF)) THEN
            WRITE(*,'(A,2I6)') 'getnewpath> ERROR, order of final minimum does not agree with database: ',HORDER,HORDERMIN(MINF)
         ENDIF
      ENDIF
C
C  Find pathways out of the new intermediate minimum if required. We can;t do this before
C  because we need all the data for new minima and transition states to be assigned in the
C  right places.
C
C     WRITE(*,'(A,I6,A)') 'getnewpath> NOT checking for at least ',CONNECTIONS,' connections per minimum'
      WRITE(*,'(A,I6,A)') 'getnewpath> checking for at least ',CONNECTIONS,' connections per minimum'
      IF (CONNECTIONS.GT.2) THEN
         DO J1=NMINOLD+1,NMIN
C              CALL TSSEARCH2(NEWMIN(J1),0)
               CALL TSSEARCH(J1,0)
         ENDDO
      ENDIF

      RETURN
      END
