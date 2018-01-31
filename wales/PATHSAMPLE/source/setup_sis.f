! GPL LIcense Info {{{
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
! }}}

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C>  Read in databases of A and B minima, calculate partition functions, rate constants
C>  and sums of rate constants for all the transition states in the database.
C
C>  We need pre-existing databases to specify which minima are A and which are B.
C
      SUBROUTINE SETUP_SIS
      ! Declarations: Modules and Variables {{{
      USE PORFUNCS
      USE UTILS
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2, STATUS, J3, NDUMMY, NRANDOM, NCOUNT, NMINREMOVE, NTSREMOVE, NMINRETAIN, NTSRETAIN, ISTAT, J4
      INTEGER SENDEMIC, IENDEMIC, POSA, POSB
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS), IXM, IYM, IZM, LOCALPOINTS2(3*NATOMS), DISTANCE, RMAT(3,3), DIST2, DPRAND
      DOUBLE PRECISION PFNORM1, PFNORM2
      DOUBLE PRECISION, ALLOCATABLE :: NEWPFMIN(:)
      INTEGER, ALLOCATABLE :: CANDIDATES(:), MINPREV(:), MINREMOVE(:), TSREMOVE(:), MINRETAIN(:), TSRETAIN(:), 
     &                        MINLABEL(:,:)
      DOUBLE PRECISION NEWEMIN, NEWIXMIN, NEWIYMIN, NEWIZMIN, NEWFVIBMIN, TSEGUESS, NEWMINCURVE, NEWMINFRQ2,
     &                 TSFVIBGUESS, DUMMY, FRICTIONFAC, INVR0
      INTEGER NEWHORDERMIN!}}}

      OPEN(UNIT=UMIN,FILE='points.min',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS) 
      OPEN(UNIT=UTS,FILE='points.ts',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS) 

! JMC set up the minima
      
      NMIN=0
      DO J1=0, SMAX
         DO J2=0, IMAX
            NMIN=NMIN+1
            IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
            EMIN(NMIN)=1.0D0
            FVIBMIN(NMIN)=1.0D0
            HORDERMIN(NMIN)=1
            IXMIN(NMIN)=1.0D0
            IYMIN(NMIN)=1.0D0
            IZMIN(NMIN)=1.0D0
         END DO
      END DO
      ALLOCATE(MINLABEL(NMIN,2))
      INVR0=(SISMU+SISKAPPA)/SISBETA
      SENDEMIC=NINT(POPSS*INVR0) ! ??? is this rounding OK?
      IENDEMIC=POPSS-SENDEMIC
      J4=0
      DO J1=0, SMAX
         IF(J1==POPSS)POSB=J4+1
         DO J2=0, IMAX
            J4=J4+1
            IF(J1==SENDEMIC .AND. J2==IENDEMIC)POSA=J4
            MINLABEL(J4,1)=J1 ! S
            MINLABEL(J4,2)=J2 ! I
         END DO
      END DO

      NMINA=1
      NMINB=1
      ALLOCATE(LOCATIONA(NMINA),LOCATIONB(NMINB))
      LOCATIONA(NMINA)=POSA ! the endemic state (S=N/R0,I=N(1-1/R0)) with R0=beta/(mu+kappa)
      LOCATIONB(NMINB)=POSB ! the (S=N,I=0) state
        
      IF (DEBUG) THEN
         PRINT '(A,I6,A)','setup> there are ',NMINA,' A minima at locations:'
         PRINT '(10I6)',LOCATIONA(1:NMINA)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> locations read for ',NMINA,' min of type A'
      IF (PRINTT) WRITE(*,'(A,2I6)') 'setup> A minimum is (S,I)=',MINLABEL(LOCATIONA(NMINA),1),MINLABEL(LOCATIONA(NMINA),2)
      IF (DEBUG) THEN
         PRINT '(A,I6,A)','setup> there are ',NMINB,' B minima at locations:'
         PRINT '(10I6)',LOCATIONB(1:NMINB)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> locations read for ',NMINB,' min of type B'
      IF (PRINTT) WRITE(*,'(A,2I6)') 'setup> B minimum is (S,I)=',MINLABEL(LOCATIONB(NMINB),1),MINLABEL(LOCATIONB(NMINB),2)
      IF (PRINTT) WRITE(*,'(A,I7,A)') 'setup> parameters read for ',NMIN,' min'
!     IF (DEBUG) WRITE(*,'(I6,2F17.7,I6,3F15.5)') (J1,EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),
!    1                                         IXMIN(J1),IYMIN(J1),IZMIN(J1),J1=1,NMIN)
      DO J1=1,NMINA
         IF (LOCATIONA(J1).GT.NMIN) THEN
            PRINT '(3(A,I8))','setup> ERROR - A minimum ',J1,' is number ',LOCATIONA(J1),' but total minima=',NMIN
            STOP
         ENDIF
      ENDDO
      DO J1=1,NMINB
         IF (LOCATIONB(J1).GT.NMIN) THEN
            PRINT '(3(A,I8))','setup> ERROR - B minimum ',J1,' is number ',LOCATIONB(J1),' but total minima=',NMIN
            STOP
         ENDIF
      ENDDO

!      IF (NMIN.GT.0) THEN
!         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE") ! read used in Dijkstra
!      ENDIF

      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)

C
C  Calculate partition functions for minima. Note that the total partition function
C  is not needed, just the totals for A and B. Since A and B sets are fixed here
C  we don;t need to change the totals.
C
!     PFMEAN=0.0D0
      PFMEAN=-HUGE(1.0D0)
!     PFNORM1=0.0D0 ! use this to calculate ratios without the pe factor
!     PFNORM2=0.0D0 ! use this to calculate ratios with the pe factor

      !op226 {{{

      IF (ENSEMBLE.EQ.'T') THEN
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A,G20.10)','setup> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
         DO J1 = 1,NMIN
            PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!           PFMEAN=PFMEAN+PFMIN(J1)
!           PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
!           PFNORM2=PFNORM2+EXP(PFMIN(J1))
            IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
         ENDDO
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         DO J1 = 1,NMIN
            IF (TOTALE.GT.EMIN(J1)) THEN
               PFMIN(J1) = (KAPPA-1)*LOG(TOTALE-EMIN(J1)) - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!              PFMEAN=PFMEAN+PFMIN(J1)
!              PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
!              PFNORM2=PFNORM2+EXP(PFMIN(J1))
               IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
            ELSE
               PFMIN(J1) = -1.0D250
            ENDIF
         ENDDO
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF
      !op226 }}}
!     PFMEAN=PFMEAN/NMIN ! DJW
      IF (DEBUG) THEN
         WRITE(*,'(A,G20.10)') 'setup> mean ln Z=',PFMEAN
!        WRITE(*,'(A)') '     energy        pg order     high T/E prob       Peq'
!        DO J1=1,NMIN
!           WRITE(*,'(F20.10,I6,2G20.10)') EMIN(J1),HORDERMIN(J1), 
!    &                    EXP(-FVIBMIN(J1)/2.0D0-LOG(1.0D0*HORDERMIN(J1))-LOG(PFNORM1)), 
!    &                    EXP(PFMIN(J1)-LOG(PFNORM2))
!        ENDDO
      ENDIF
      DO J1=1,NMIN
         PFMIN(J1) = PFMIN(J1) - PFMEAN
      ENDDO

      PFTOTALB=0.0D0
      DO J1=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J1))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J1=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J1))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))
C
C  Load transition states.
C
      DO J1=1,NMIN
         TOPPOINTER(J1)=-1
      ENDDO
      
C Order the TS's by event...
      NTS=0

C First, recovery/infection
      DO J1=1,NMIN
         IF(MINLABEL(J1,2)==0 .OR. MINLABEL(J1,1)==SMAX)CYCLE ! if s=anything,i=0 or s=smax,i=anything
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         ETS(NTS)=1.0D0
         FVIBTS(NTS)=1.0D0
         HORDERTS(NTS)=1
         PLUS(NTS)=J1
         MINUS(NTS)=J1+IMAX ! corresponding to the scheme used above for setting up the minima, by row of I's.
         IXTS(NTS)=1.0D0
         IYTS(NTS)=1.0D0
         IZTS(NTS)=1.0D0
         KPLUS(NTS)=LOG(SISKAPPA*DBLE(MINLABEL(J1,2))) ! recovery
         IF (MINLABEL(J1+IMAX,2)==0) THEN
            KMINUS(NTS)=-299.0D0
         ELSE
            KMINUS(NTS)=LOG(SISBETA*DBLE(MINLABEL(J1+IMAX,1)*MINLABEL(J1+IMAX,2))/DBLE(POPSS)) ! infection
         END IF
         IF ((PLUS(NTS).GT.NMIN).OR.(MINUS(NTS).GT.NMIN)) THEN
            PRINT '(A,I6,A)','setup> ERROR - minima specified for ts ',NTS,' lie beyond those specified in min.data'
            PRINT '(A,2I6)','setup> plus and minus minima are ',PLUS(NTS),MINUS(NTS)
            STOP
         ENDIF
         IF (DEBUG) WRITE(*,'(A,3I6,5G20.10)') 'setup 1> J1,PLUS,MINUS,Ets,E+,E-,k+,k-,<k>=',NTS,PLUS(NTS),MINUS(NTS),
     1                                            ETS(NTS),EMIN(PLUS(NTS)),EMIN(MINUS(NTS)),KPLUS(NTS),KMINUS(NTS)
      ENDDO

C Second, birth/death of susceptibles
      DO J1=1,NMIN
         IF(MINLABEL(J1,1)==SMAX)CYCLE ! if s=smax, i=anything
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         ETS(NTS)=1.0D0
         FVIBTS(NTS)=1.0D0
         HORDERTS(NTS)=1
         PLUS(NTS)=J1
         MINUS(NTS)=J1+IMAX+1 ! corresponding to the scheme used above for setting up the minima, by row of I's.
         IXTS(NTS)=1.0D0
         IYTS(NTS)=1.0D0
         IZTS(NTS)=1.0D0
         KPLUS(NTS)=LOG(SISMU*DBLE(POPSS)) ! birth
         KMINUS(NTS)=LOG(SISMU*DBLE(MINLABEL(J1+IMAX+1,1))) ! death
         IF ((PLUS(NTS).GT.NMIN).OR.(MINUS(NTS).GT.NMIN)) THEN
            PRINT '(A,I6,A)','setup> ERROR - minima specified for ts ',NTS,' lie beyond those specified in min.data'
            PRINT '(A,2I6)','setup> plus and minus minima are ',PLUS(NTS),MINUS(NTS)
            STOP
         ENDIF
         IF (DEBUG) WRITE(*,'(A,3I6,5G20.10)') 'setup 2> J1,PLUS,MINUS,Ets,E+,E-,k+,k-,<k>=',NTS,PLUS(NTS),MINUS(NTS),
     1                                            ETS(NTS),EMIN(PLUS(NTS)),EMIN(MINUS(NTS)),KPLUS(NTS),KMINUS(NTS)
      ENDDO

C Third, death of infecteds
      DO J1=1,NMIN
         IF(MINLABEL(J1,2)==0)CYCLE ! if s=anything, i=0
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         ETS(NTS)=1.0D0
         FVIBTS(NTS)=1.0D0
         HORDERTS(NTS)=1
         PLUS(NTS)=J1
         MINUS(NTS)=J1-1 ! corresponding to the scheme used above for setting up the minima, by row of I's.
         IXTS(NTS)=1.0D0
         IYTS(NTS)=1.0D0
         IZTS(NTS)=1.0D0
         KPLUS(NTS)=LOG(SISMU*DBLE(MINLABEL(J1,2)))
         KMINUS(NTS)=-299.0D0 ! no birth of infecteds
         IF (DEBUG) WRITE(*,'(A,3I6,5G20.10)') 'setup 3> J1,PLUS,MINUS,Ets,E+,E-,k+,k-,<k>=',NTS,PLUS(NTS),MINUS(NTS),
     1                                            ETS(NTS),EMIN(PLUS(NTS)),EMIN(MINUS(NTS)),KPLUS(NTS),KMINUS(NTS)
      ENDDO

!      OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE",FORM="FORMATTED") ! read used in Dijkstra
      IF (PRINTT) WRITE(*,'(A,I7,A)') 'setup> parameters read for ',NTS,' ts'

      IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN
         INQUIRE(FILE='ts.attempts',EXIST=YESNO)
         TSATTEMPT(1:NTS)=0
         IF (YESNO) THEN
            PRINT '(A)','setup> Reading the number of searches for existing transition states from ts.attempts'
            OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
            J2=0
            DO J1=1,NTS
               READ(1,'(I8)',END=51) TSATTEMPT(J1)
               J2=J2+1
            ENDDO
51          CLOSE(1)
            IF (J2.LT.NTS) PRINT '(A)','setup> WARNING - end of file ts.attempts, remaining attempts set to zero'
         ENDIF
      ENDIF

      DO J1=1,NTS
         POINTERP(J1)=-1
         POINTERM(J1)=-1
      ENDDO
      DO J1=NTS,1,-1
         IF (J1.GT.TOPPOINTER(PLUS(J1)))  TOPPOINTER(PLUS(J1))=J1
         IF (J1.GT.TOPPOINTER(MINUS(J1))) TOPPOINTER(MINUS(J1))=J1
         DO J2=J1-1,1,-1
            IF (PLUS(J2).EQ.PLUS(J1)) THEN
               POINTERP(J1)=J2
               GOTO 41
            ELSE IF (MINUS(J2).EQ.PLUS(J1)) THEN
               POINTERP(J1)=J2
               GOTO 41
            ENDIF
         ENDDO
41       CONTINUE
         DO J2=J1-1,1,-1
            IF (PLUS(J2).EQ.MINUS(J1)) THEN
               POINTERM(J1)=J2
               GOTO 42
            ELSE IF (MINUS(J2).EQ.MINUS(J1)) THEN
               POINTERM(J1)=J2
               GOTO 42
            ENDIF
         ENDDO
42       CONTINUE
      ENDDO
!     IF (DEBUG) THEN
!        DO J1=1,NMIN
!           WRITE(*,'(A,I6,A,I6)') 'setup> for minimum ',J1,' last occurrence is for ts number ',TOPPOINTER(J1)
!        ENDDO
!        DO J1=1,NTS
!           WRITE(*,'(A,I6,A,4I6)') 'setup> for ts ',J1,' +,-,p+,p-:',PLUS(J1),MINUS(J1),POINTERP(J1),POINTERM(J1)
!        ENDDO
!     ENDIF
    
      IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
C
C  Sums of rates out of the intermediate minima
C
!op226 {{{
!       DO J1=1,NMIN
!          KSUM(J1)=0.0D0
!       ENDDO
!       DO J1=1,NTS
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(PLUS(J1))=KSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(MINUS(J1))=KSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
!       ENDDO
!       DO J1=1,NMIN
!          IF (KSUM(J1).GT.0.0D0) THEN
!             KSUM(J1)=LOG(KSUM(J1))+KMEAN
! !           IF (DEBUG) WRITE(*,'(A,I6,2E20.10)') 'setup> J1,KSUM=',J1,KSUM(J1)
!          ENDIF
!       ENDDO
!       DO J1=1,NTS
! !        IF (DEBUG) WRITE(*,'(A,I6,2E20.10)') 'setup> J1,k+,k-=',J1,KPLUS(J1),KMINUS(J1)
!       ENDDO
!op226 }}}

      IF (NPFOLD.GT.0) THEN
         INQUIRE(FILE='commit.data',EXIST=YESNO)
         GPFOLD(1:NMIN)=0.0D0
         IF (YESNO) THEN
            PRINT '(A)','setup> Reading initial committor probabilities read from commit.data'
            OPEN(UNIT=1,FILE='commit.data',STATUS='OLD')
            J2=0
            DO J1=1,NMIN
               READ(1,*,END=110) GPFOLD(J1)
               J2=J2+1
            ENDDO 
110         CLOSE(1)
            IF (J2.LT.NMIN) PRINT '(A)','setup> WARNING - end of file commit.data, remaining probabilities set to zero'
         ELSE
            IF (DIRECTION.EQ.'AB') THEN ! GPFOLD is PFA
               DO J1=1,NMINA
                  GPFOLD(LOCATIONA(J1))=1.0D0
               ENDDO
            ELSE ! GPFOLD is PFB
               DO J1=1,NMINB
                  GPFOLD(LOCATIONB(J1))=1.0D0
               ENDDO
            ENDIF
            PRINT '(A)','setup> Initial committor probabilities set to 0 or 1'
!           PRINT '(6G20.10)',GPFOLD(1:NMIN)
         ENDIF
      ENDIF
C
C  Read in the pairs of minima previously searched in pairs.data exists.
C
!{{{
      ALLOCATE(PAIR1(MAXPAIRS),PAIR2(MAXPAIRS))
      INQUIRE(FILE='pairs.data',EXIST=YESNO)
      NPAIRDONE=0
      IF (YESNO) THEN
         OPEN(UNIT=1,FILE='pairs.data',STATUS='OLD')
         DO 
            NPAIRDONE=NPAIRDONE+1
            IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
            READ(1,*,END=120) PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched pair number ',
     &                                      NPAIRDONE,' is ',PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
         ENDDO
120      CLOSE(1)
         NPAIRDONE=NPAIRDONE-1
         PRINT '(A,I8,A)','setup> ',NPAIRDONE,' pairs of minima already searched read from pairs.data'
      ENDIF
! }}}
C
C  Read in the minima previously searched in UNTRAP runs.
C
! {{{
      ALLOCATE(MINDONE(MAXDONE))
      INQUIRE(FILE='min.done',EXIST=YESNO)
      NMINDONE=0
      IF (YESNO) THEN
         OPEN(UNIT=1,FILE='min.done',STATUS='OLD')
         DO 
            NMINDONE=NMINDONE+1
            IF (NMINDONE.GT.MAXDONE) CALL DONEDOUBLE
            READ(1,*,END=121) MINDONE(NMINDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched minimum ',
     &                                      NMINDONE,' is ',MINDONE(NMINDONE)
         ENDDO
121      CLOSE(1)
         NMINDONE=NMINDONE-1
         PRINT '(A,I8,A)','setup> ',NMINDONE,' minima already searched read from min.done'
      ENDIF
! }}}

      IF (DOST) CALL DOS
      IF (CVT) CALL CV

      RETURN 
      END
