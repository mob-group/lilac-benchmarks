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
C  Check transition state NTS and the minima that it connects to see if
C  permutational isomers can be added to those currently known.
C
C  One potential problem: the orders of the point groups are not
C  calculated for new stationary points. Of course, they are all one for
C  this system. This affects the rates as well, so we can.t just copy
C  the corresponding quantities from the original transition state.
C
      SUBROUTINE ADDPERM(LPOINTSTS,LPLUS,LMINUS)
      USE PORFUNCS
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1,J2,DOTS,ISTAT
      DOUBLE PRECISION LPOINTS(NOPT), LDUMMY,NIX,NIY,NIZ,LPOINTSTS(NOPT),LPLUS(NOPT),
     1                 LMINUS(NOPT), FRICTIONFAC

      PRINT '(A)','WARNING- routine ADDPERM has not been tested'
      IF (NTAG.LE.0) THEN
         WRITE(*,'(A,I5)') 'WARNING - addperm called but NTAG=',NTAG
         RETURN
      ENDIF
      DOTS=NTS

      DO J1=2,NATOMS
         DO J2=1,NOPT
            LPOINTS(J2)=LPOINTSTS(J2)
         ENDDO
C
C  Permute tagged atom NTAG and atom J1
C
         LDUMMY=LPOINTS(3*(NTAG-1)+1)
         LPOINTS(3*(NTAG-1)+1)=LPOINTS(3*(J1-1)+1)
         LPOINTS(3*(J1-1)+1)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+2)
         LPOINTS(3*(NTAG-1)+2)=LPOINTS(3*(J1-1)+2)
         LPOINTS(3*(J1-1)+2)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+3)
         LPOINTS(3*(NTAG-1)+3)=LPOINTS(3*(J1-1)+3)
         LPOINTS(3*(J1-1)+3)=LDUMMY

         CALL INERTIAWRAPPER(LPOINTS,NOPT,angleAxis,NIX,NIY,NIZ)
         DO J2=1,NTS
            IF ((ABS(ETS(J2)-ETS(DOTS)).LT.EDIFFTOL).AND.
     1          (ABS(IXTS(J2)-NIX).LT.IDIFFTOL).AND.
     2          (ABS(IYTS(J2)-NIY).LT.IDIFFTOL).AND.
     3          (ABS(IZTS(J2)-NIZ).LT.IDIFFTOL)) THEN
               WRITE(*,'(A,I5,A,I5,A,I5)') 
     1                 'permuting atoms 1 and ',J1,' for ts ',DOTS,' gives old transition state ',J2
               GOTO 10
            ENDIF
         ENDDO
         WRITE(*,'(A,I5,A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for ts ',DOTS,' gives new transition state ',NTS+1
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         ETS(NTS)=ETS(DOTS)
         FVIBTS(NTS)=FVIBTS(DOTS)
         IF (FRICTIONT) NEGEIG(NTS)=NEGEIG(DOTS)
         HORDERTS(NTS)=1          ! not valid in general ! see next line
         IF ((NATOMS.NE.7).OR.(.NOT.TWOD)) THEN
            PRINT*,'addperm assumes all minima have point group order 1 - quit'
            STOP
         ENDIF
         IXTS(NTS)=NIX
         IYTS(NTS)=NIY
         IZTS(NTS)=NIZ
         WRITE(UTS,REC=NTS) (LPOINTS(J2),J2=1,NOPT)
         call flush(UTS)
C
C  Identify the connected minimum.
C
         DO J2=1,NOPT
            LPOINTS(J2)=LPLUS(J2)
         ENDDO
         LDUMMY=LPOINTS(3*(NTAG-1)+1)
         LPOINTS(3*(NTAG-1)+1)=LPOINTS(3*(J1-1)+1)
         LPOINTS(3*(J1-1)+1)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+2)
         LPOINTS(3*(NTAG-1)+2)=LPOINTS(3*(J1-1)+2)
         LPOINTS(3*(J1-1)+2)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+3)
         LPOINTS(3*(NTAG-1)+3)=LPOINTS(3*(J1-1)+3)
         LPOINTS(3*(J1-1)+3)=LDUMMY

         CALL INERTIAWRAPPER(LPOINTS,NOPT,angleAxis,NIX,NIY,NIZ)
         DO J2=1,NMIN
            IF ((ABS(EMIN(PLUS(DOTS))-EMIN(J2)).LT.EDIFFTOL).AND.
     1          (ABS(NIX-IXMIN(J2)).LT.IDIFFTOL).AND.
     2          (ABS(NIY-IYMIN(J2)).LT.IDIFFTOL).AND.
     3          (ABS(NIZ-IZMIN(J2)).LT.IDIFFTOL)) THEN
               WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for plus minimum gives old minimum ',J2
               PLUS(NTS)=J2
               IF (ENSEMBLE.EQ.'T') THEN
                  KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(J2)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1                       (FVIBMIN(J2)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J2))/TEMPERATURE
                  IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(J2)  / (2*PI*HORDERTS(NTS))) +
     1               (FVIBMIN(J2)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J2)))
                  ELSE
                      KPLUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
               IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!              IF (KSUM(J2).EQ.0.0D0) THEN
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2)=KPLUS(NTS)
!              ELSE
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2) =LOG(EXP(KSUM(J2)-KMEAN) + EXP(KPLUS(NTS) -KMEAN)) + KMEAN
!              ENDIF
               GOTO 11
            ENDIF
         ENDDO
         WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for plus minimum gives new minimum ',NMIN+1
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         EMIN(NMIN)=EMIN(PLUS(DOTS))
         FVIBMIN(NMIN)=FVIBMIN(PLUS(DOTS))
         HORDERMIN(NMIN)=1          ! not valid in general !
         IXMIN(NMIN)=NIX
         IYMIN(NMIN)=NIY
         IZMIN(NMIN)=NIZ
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

         WRITE(UMIN,REC=NMIN) (LPOINTS(J2),J2=1,NOPT)
         CALL FLUSH(UMIN)
         IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
         CALL FLUSH(UMINDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)

         PLUS(NTS)=NMIN
         IF (ENSEMBLE.EQ.'T') THEN
            KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1                 (FVIBMIN(NMIN)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
            IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
         ELSE
            IF (TEMPERATURE.GT.ETS(NTS)) THEN
               KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(NMIN)  / (2*PI*HORDERTS(NTS))) +
     1            (FVIBMIN(NMIN)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
            ELSE
               KPLUS(NTS)=-1.0D250
            ENDIF
         ENDIF
         IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
         IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!        KSUM(NMIN)=0.0D0
!        IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KPLUS(NTS)
11       CONTINUE

         DO J2=1,NOPT
            LPOINTS(J2)=LMINUS(J2)
         ENDDO
         LDUMMY=LPOINTS(3*(NTAG-1)+1)
         LPOINTS(3*(NTAG-1)+1)=LPOINTS(3*(J1-1)+1)
         LPOINTS(3*(J1-1)+1)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+2)
         LPOINTS(3*(NTAG-1)+2)=LPOINTS(3*(J1-1)+2)
         LPOINTS(3*(J1-1)+2)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+3)
         LPOINTS(3*(NTAG-1)+3)=LPOINTS(3*(J1-1)+3)
         LPOINTS(3*(J1-1)+3)=LDUMMY

         CALL INERTIAWRAPPER(LPOINTS,NOPT,angleAxis,NIX,NIY,NIZ)
         DO J2=1,NMIN
            IF ((ABS(EMIN(MINUS(DOTS))-EMIN(J2)).LT.EDIFFTOL).AND.
     1          (ABS(NIX-IXMIN(J2)).LT.IDIFFTOL).AND.
     2          (ABS(NIY-IYMIN(J2)).LT.IDIFFTOL).AND.
     3          (ABS(NIZ-IZMIN(J2)).LT.IDIFFTOL)) THEN
               WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for minus minimum gives old minimum ',J2
               MINUS(NTS)=J2
               IF (ENSEMBLE.EQ.'T') THEN
                  KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(J2) / (2.0D0 * PI*HORDERTS(NTS))) +
     1                     (FVIBMIN(J2) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J2))/TEMPERATURE
                  IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1                   (FVIBMIN(J2) - FVIBTS(NTS))/2 + 
     2                   (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J2)))
                  ELSE
                     KMINUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
               IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!              IF (KSUM(J2).EQ.0.0D0) THEN
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2)=KMINUS(NTS)
!              ELSE
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2)=LOG(EXP(KSUM(J2)-KMEAN) + EXP(KMINUS(NTS)-KMEAN)) + KMEAN
!              ENDIF
               GOTO 12
            ENDIF
         ENDDO
         WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for minus minimum gives new minimum ',NMIN+1
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         EMIN(NMIN)=EMIN(MINUS(DOTS))
         FVIBMIN(NMIN)=FVIBMIN(MINUS(DOTS))
         HORDERMIN(NMIN)=1          ! not valid in general !
         IXMIN(NMIN)=NIX
         IYMIN(NMIN)=NIY
         IZMIN(NMIN)=NIZ
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

         WRITE(UMIN,REC=NMIN) (LPOINTS(J2),J2=1,NOPT)
         call flush(UMIN)
         IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN),IZMIN(NMIN)
         CALL FLUSH(UMINDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
         MINUS(NTS)=NMIN
         IF (ENSEMBLE.EQ.'T') THEN
            KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN) / (2.0D0 * PI*HORDERTS(NTS))) +
     1               (FVIBMIN(NMIN) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
            IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
         ELSE
            IF (TEMPERATURE.GT.ETS(NTS)) THEN
               KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(NMIN) / (2*PI*HORDERTS(NTS))) +
     1             (FVIBMIN(NMIN) - FVIBTS(NTS))/2 + 
     2             (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
            ELSE
               KMINUS(NTS)=-1.0D250
            ENDIF
         ENDIF
         IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
         IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!        KSUM(NMIN)=0.0D0
!        IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KMINUS(NTS)
12       CONTINUE
C
C  Only now do we know the connected minima for sure.
C
         IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
         IF (IMFRQT) THEN
            WRITE(UTSDATA,'(2F25.15,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
         ELSE
            WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
         ENDIF
         CALL FLUSH(UTSDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
C
C  Update ts pointers.
C
         POINTERP(NTS)=-1
         POINTERM(NTS)=-1
         TOPPOINTER(PLUS(NTS))=NTS
         TOPPOINTER(MINUS(NTS))=NTS

         DO J2=NTS-1,1,-1
            IF (PLUS(J2).EQ.PLUS(NTS)) THEN
               POINTERP(NTS)=J2
               GOTO 41
            ELSE IF (MINUS(J2).EQ.PLUS(NTS)) THEN
               POINTERP(NTS)=J2
               GOTO 41
            ENDIF
         ENDDO
41       CONTINUE

         DO J2=NTS-1,1,-1
            IF (PLUS(J2).EQ.MINUS(NTS)) THEN
               POINTERM(NTS)=J2
               GOTO 42
            ELSE IF (MINUS(J2).EQ.MINUS(NTS)) THEN
               POINTERM(NTS)=J2
               GOTO 42
            ENDIF
         ENDDO
42       CONTINUE
C
C  ts pointers have been updated.
C
10       CONTINUE
      ENDDO

      RETURN
      END
C
C  Check transition state NTS and the minima that it connects to see if
C  permutation-inversion isomers can be added to those currently known.
C
C  One potential problem: the orders of the point groups are not
C  calculated for new stationary points. 
C
C  No tagging assumed here.
C
      SUBROUTINE ADDPERM2(LDOTS,LPOINTSTS,LPLUS,LMINUS)
      USE PORFUNCS
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1,J2,DOTS,ISTAT,J3,J4,NDUMMY,NPERM,INVERSE, NPSIZE,LDOTS
      INTEGER, ALLOCATABLE :: PERM(:), PERM2(:), REALPERM(:)
      DOUBLE PRECISION LPOINTS(NOPT),NIX,NIY,NIZ,LPOINTSTS(NOPT),LPLUS(NOPT),
     1                 LMINUS(NOPT), FRICTIONFAC, RMAT(3,3), DISTANCE, LOCALPOINTS2(NOPT), DIST2

      NPSIZE=PTFINISH-PTSTART+1
      PRINT '(A,3I6)','PTSTART,PTFINISH,NPSIZE=',PTSTART,PTFINISH,NPSIZE
      ALLOCATE(REALPERM(NATOMS),PERM2(NPSIZE),PERM(NPSIZE))
      DOTS=LDOTS

      DO J1=1,NATOMS
         REALPERM(J1)=J1
      ENDDO
      DO J1=1,NPSIZE
         PERM(J1)=J1
      ENDDO
      NPERM=1
      INVERSE=1
      GOTO 333 ! This branch ensures that we check the identity permutation and the corresponding inverse.
               ! Otherwise we miss the enantiomer of the original transition state for C1 symetry!
!
! Find the highest index i such that s[i] < s[i+1]. If no such index exists, the permutation is the last permutation.
! Begin enumeration of permutations.
!
116   INVERSE=1
      J2=-1
      DO J1=1,NPSIZE-1
         IF (PERM(J1).LT.PERM(J1+1)) J2=J1
      ENDDO
      IF (J2.LT.0) GOTO 20
!
! Find the highest index j > i such that s[j] > s[i]. Such a j must exist, since i+1 is such an index.
!
      DO J1=1,NPSIZE
         IF (PERM(J1).GT.PERM(J2)) J3=J1
      ENDDO
!
! Swap s[i] with s[j].
!
      NDUMMY=PERM(J2)
      PERM(J2)=PERM(J3)
      PERM(J3)=NDUMMY
!
! Reverse all the order of all of the elements after index i
!
      DO J1=1,J2
         PERM2(J1)=PERM(J1)
      ENDDO
      DO J1=J2+1,NPSIZE
         PERM2(J1)=PERM(NPSIZE+J2-J1+1)
      ENDDO
      PERM(1:NPSIZE)=PERM2(1:NPSIZE)
      NPERM=NPERM+1
      PRINT '(20I6)',NPERM,PERM(1:NPSIZE)
      DO J1=1,NPSIZE
         REALPERM(J1+PTSTART-1)=PERM(J1)+PTSTART-1
      ENDDO
333   PRINT '(20I6)',NPERM,REALPERM(1:NATOMS)

C
C  Try this permutation.
C
222   CONTINUE
      DO J1=1,NATOMS
         LPOINTS(3*(REALPERM(J1)-1)+1)=INVERSE*LPOINTSTS(3*(J1-1)+1)
         LPOINTS(3*(REALPERM(J1)-1)+2)=INVERSE*LPOINTSTS(3*(J1-1)+2)
         LPOINTS(3*(REALPERM(J1)-1)+3)=INVERSE*LPOINTSTS(3*(J1-1)+3)
      ENDDO

      tsloop: DO J3=1,NTS
         IF (ABS(ETS(J3)-ETS(DOTS)).LT.EDIFFTOL) THEN
            READ(UTS,REC=J3) (LOCALPOINTS2(J4),J4=1,NOPT)
            CALL MINPERMDIST(LPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                    RMAT,.FALSE.)
            IF (DISTANCE.LT.GEOMDIFFTOL) THEN
               PRINT '(A,I6)','addperm2> permutation gives database TS ',J3
               GOTO 130
            ENDIF
         ENDIF
      ENDDO tsloop
      PRINT '(A,2I6,A,I6)','addperm2> permutation gives new TS ',NTS+1
      NTS=NTS+1
      IF (NTS.GT.MAXTS) CALL TSDOUBLE
      ETS(NTS)=ETS(DOTS)
      FVIBTS(NTS)=FVIBTS(DOTS)
      IF (FRICTIONT) NEGEIG(NTS)=NEGEIG(DOTS)
      HORDERTS(NTS)=1          ! not valid in general ! see next line
      CALL INERTIAWRAPPER(LPOINTS,NOPT,ANGLEAXIS,NIX,NIY,NIZ)
      IXTS(NTS)=NIX
      IYTS(NTS)=NIY
      IZTS(NTS)=NIZ
      WRITE(UTS,REC=NTS) (LPOINTS(J3),J3=1,NOPT)
      CALL FLUSH(UTS)
C
C  Identify the connected minimum.
C
      DO J3=1,NATOMS
         LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LPLUS(3*(J3-1)+1)
         LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LPLUS(3*(J3-1)+2)
         LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LPLUS(3*(J3-1)+3)
      ENDDO

      DO J3=1,NMIN
         IF (ABS(EMIN(J3)-EMIN(PLUS(DOTS))).LT.EDIFFTOL) THEN
            READ(UMIN,REC=J3) (LOCALPOINTS2(J4),J4=1,NOPT)
            CALL MINPERMDIST(LPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                 RMAT,.FALSE.)
            IF (DISTANCE.LT.GEOMDIFFTOL) THEN
               WRITE(*,'(A,I6)') 'permutation for plus minimum gives old minimum ',J3
               PLUS(NTS)=J3
               IF (ENSEMBLE.EQ.'T') THEN
                  KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(J3)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1                 (FVIBMIN(J3)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J3))/TEMPERATURE
                  IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
                  ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(J3)  / (2*PI*HORDERTS(NTS))) +
     1             (FVIBMIN(J3)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J3)))
                  ELSE
                      KPLUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
               GOTO 11
            ENDIF
         ENDIF
      ENDDO
      WRITE(*,'(A,I6)') 'permutation for plus minimum gives new minimum ',NMIN+1
      NMIN=NMIN+1
      IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
      EMIN(NMIN)=EMIN(PLUS(DOTS))
      FVIBMIN(NMIN)=FVIBMIN(PLUS(DOTS))
      HORDERMIN(NMIN)=1          ! not valid in general !
      CALL INERTIAWRAPPER(LPOINTS,NOPT,ANGLEAXIS,NIX,NIY,NIZ)
      IXMIN(NMIN)=NIX
      IYMIN(NMIN)=NIY
      IZMIN(NMIN)=NIZ
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

      WRITE(UMIN,REC=NMIN) (LPOINTS(J3),J3=1,NOPT)
      CALL FLUSH(UMIN)
      IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
      CALL FLUSH(UMINDATA)
      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)

      PLUS(NTS)=NMIN
      IF (ENSEMBLE.EQ.'T') THEN
         KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1           (FVIBMIN(NMIN)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
         IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
      ELSE
         IF (TEMPERATURE.GT.ETS(NTS)) THEN
            KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(NMIN)  / (2*PI*HORDERTS(NTS))) +
     1      (FVIBMIN(NMIN)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
         ELSE
            KPLUS(NTS)=-1.0D250
         ENDIF
      ENDIF
      IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
      IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!     KSUM(NMIN)=0.0D0
!     IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KPLUS(NTS)
11    CONTINUE

      DO J3=1,NATOMS
         LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LMINUS(3*(J3-1)+1)
         LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LMINUS(3*(J3-1)+2)
         LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LMINUS(3*(J3-1)+3)
      ENDDO

      DO J3=1,NMIN
         IF (ABS(EMIN(J3)-EMIN(MINUS(DOTS))).LT.EDIFFTOL) THEN
            READ(UMIN,REC=J3) (LOCALPOINTS2(J4),J4=1,NOPT)
            CALL MINPERMDIST(LPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                 RMAT,.FALSE.)
            IF (DISTANCE.LT.GEOMDIFFTOL) THEN
               WRITE(*,'(A,I6)') 'permutation for minus minimum gives old minimum ',J3
               MINUS(NTS)=J3
               IF (ENSEMBLE.EQ.'T') THEN
                  KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(J3) / (2.0D0 * PI*HORDERTS(NTS))) +
     1               (FVIBMIN(J3) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J3))/TEMPERATURE
                  IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1             (FVIBMIN(J3) - FVIBTS(NTS))/2 + 
     2             (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J3)))
                  ELSE
                     KMINUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
               GOTO 12
            ENDIF
         ENDIF
      ENDDO
      WRITE(*,'(A,I5,A,I5)') 'permutation for minus minimum gives new minimum ',NMIN+1
      NMIN=NMIN+1
      IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
      EMIN(NMIN)=EMIN(MINUS(DOTS))
      FVIBMIN(NMIN)=FVIBMIN(MINUS(DOTS))
      HORDERMIN(NMIN)=1          ! not valid in general !
      CALL INERTIAWRAPPER(LPOINTS,NOPT,ANGLEAXIS,NIX,NIY,NIZ)
      IXMIN(NMIN)=NIX
      IYMIN(NMIN)=NIY
      IZMIN(NMIN)=NIZ
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

      WRITE(UMIN,REC=NMIN) (LPOINTS(J3),J3=1,NOPT)
      CALL FLUSH(UMIN)
      IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN),IZMIN(NMIN)
      CALL FLUSH(UMINDATA)
      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
      MINUS(NTS)=NMIN
      IF (ENSEMBLE.EQ.'T') THEN
         KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN) / (2.0D0 * PI*HORDERTS(NTS))) +
     1         (FVIBMIN(NMIN) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
         IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
      ELSE
         IF (TEMPERATURE.GT.ETS(NTS)) THEN
            KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(NMIN) / (2*PI*HORDERTS(NTS))) +
     1       (FVIBMIN(NMIN) - FVIBTS(NTS))/2 + 
     2       (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
         ELSE
            KMINUS(NTS)=-1.0D250
         ENDIF
      ENDIF
!     IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KMINUS(NTS)
12    CONTINUE
C
C  Only now do we know the connected minima for sure.
C
      IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
      IF (IMFRQT) THEN
         WRITE(UTSDATA,'(2F25.15,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
      ELSE
         WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
      ENDIF
      CALL FLUSH(UTSDATA)
      IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
C
C  Update ts pointers.
C
      POINTERP(NTS)=-1
      POINTERM(NTS)=-1
      TOPPOINTER(PLUS(NTS))=NTS
      TOPPOINTER(MINUS(NTS))=NTS

      DO J3=NTS-1,1,-1
         IF (PLUS(J3).EQ.PLUS(NTS)) THEN
            POINTERP(NTS)=J3
            GOTO 41
         ELSE IF (MINUS(J3).EQ.PLUS(NTS)) THEN
            POINTERP(NTS)=J3
            GOTO 41
         ENDIF
      ENDDO
41    CONTINUE

      DO J3=NTS-1,1,-1
         IF (PLUS(J3).EQ.MINUS(NTS)) THEN
            POINTERM(NTS)=J3
            GOTO 42
         ELSE IF (MINUS(J3).EQ.MINUS(NTS)) THEN
            POINTERM(NTS)=J3
            GOTO 42
         ENDIF
      ENDDO
42    CONTINUE
C
C  ts pointers have been updated.
C
10    CONTINUE
130   CONTINUE

      IF (INVERSE.EQ.1) THEN
         INVERSE=-1
         PRINT '(A)','addperm2> Trying permutation plus inversion'
         GOTO 222
      ENDIF
      GOTO 116

20    CONTINUE
      DEALLOCATE(REALPERM,PERM,PERM2)

      RETURN
      END
C
C  Check transition state NTS and the minima that it connects to see if
C  permutation-inversion isomers can be added to those currently known.
C
C  Use permutations specified in perm.allow.
C
      SUBROUTINE ADDPERM3(LDOTS,LPOINTSTS,LPLUS,LMINUS)
      USE PORFUNCS
      USE COMMONS
      USE UTILS,ONLY : GETUNIT
      IMPLICIT NONE
      INTEGER J1,LDOTS,ISTAT,J3,J4,NDUMMY,NPERMTOTAL,INVERSE,K1,NDUMMY2,NPSIZEMAX,K2,DOTS
      LOGICAL CHANGED, LLPERMDIST
      INTEGER, ALLOCATABLE :: NPSIZE(:), NPERM(:), J2(:)
      INTEGER, ALLOCATABLE :: PERM(:,:), PERM2(:,:), REALPERM(:)
      DOUBLE PRECISION LPOINTS(NOPT),NIX,NIY,NIZ,LPOINTSTS(NOPT),LPLUS(NOPT),
     1                 LMINUS(NOPT), FRICTIONFAC, RMAT(3,3), DISTANCE, LOCALPOINTS2(NOPT), DIST2

!
!  Need to combine all possible permutations of each set with one another:
!  a product of factorials of their size.
!  NPERMGROUP is the number of groups
!  NPERMSIZE(J1) is the number of permutable atoms in this group
!
      ALLOCATE(NPSIZE(NPERMGROUP),J2(NPERMGROUP),NPERM(NPERMGROUP))
      ALLOCATE(REALPERM(NATOMS))
      DO J1=1,NATOMS
         REALPERM(J1)=J1
      ENDDO
      DOTS=LDOTS

      NPSIZEMAX=1
      DO K1=1,NPERMGROUP
         IF (NPERMSIZE(K1).GT.NPSIZEMAX) NPSIZEMAX=NPERMSIZE(K1)
      ENDDO
      PRINT '(A,I6,A,I6)','addperm3> Largest number of permutable atoms=',NPSIZEMAX,' doing TS number ',DOTS
      ALLOCATE(PERM2(NPSIZEMAX,NPERMGROUP),PERM(NPSIZEMAX,NPERMGROUP))

      DO K1=1,NPERMGROUP
         NPSIZE(K1)=NPERMSIZE(K1)
         DO J1=1,NPSIZEMAX ! trailing entries for smaller groups are just dummies
            PERM(J1,K1)=J1
         ENDDO
      ENDDO

      LLPERMDIST=PERMDIST
      PERMDIST=.FALSE.
      NPERMTOTAL=1
      NPERM(1:NPERMGROUP)=1
      INVERSE=1
      K1=1

      GOTO 333 ! This branch ensures that we check the identity permutation and the corresponding inverted structures.
               ! Otherwise we miss the enantiomer of the original transition state for C1 symetry!
!
! Find the highest index i such that s[i] < s[i+1]. If no such index exists, the permutation is the last permutation.
! Begin enumeration of permutations.
!
116   INVERSE=1
      CHANGED=.FALSE.
      PRINT '(A,I6,A,I6)','Incrementing permutation for group ',K1,' inverse=',INVERSE
!
! Generate all possible permutations for each set of permutable atoms.
!
      NDUMMY=1
!
! This block generates the next permutation for group K1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         J2(K1)=-1
         DO J1=1,NPSIZE(K1)-1
            IF (PERM(J1,K1).LT.PERM(J1+1,K1)) J2(K1)=J1
         ENDDO
!
! We have made all the possible permutations for group K1
! If this is the last group, then we have done all possible permutations.
! If not, reset all permutations up to and including group K1 to the
! identify, and increement the permutation for K1+1.
!
         IF (J2(K1).LT.0) THEN 
            IF (K1.EQ.NPERMGROUP) THEN
               PRINT '(A,I6)','addperm3> All permutations checked, total=',NPERMTOTAL
               GOTO 20 
            ELSE
               DO K2=1,K1
                  DO J1=1,NPSIZEMAX ! trailing entries for smaller groups are just dummies
                     PERM(J1,K2)=J1
                  ENDDO
               ENDDO
               NPERM(1:K1)=1
               PRINT '(A,I6)','addperm3> Resetting permutations to the identify for groups one to ',K1
               K1=K1+1
               CHANGED=.TRUE.
               GOTO 116
            ENDIF
         ENDIF 
!
! Find the highest index j > i such that s[j] > s[i]. Such a j must exist, since i+1 is such an index.
!
         DO J1=1,NPSIZE(K1)
            IF (PERM(J1,K1).GT.PERM(J2(K1),K1)) J3=J1
         ENDDO
!
! Swap s[i] with s[j].
!
         NDUMMY2=PERM(J2(K1),K1)
         PERM(J2(K1),K1)=PERM(J3,K1)
         PERM(J3,K1)=NDUMMY2
!
! Reverse all the order of all of the elements after index i
!
         DO J1=1,J2(K1)
            PERM2(J1,K1)=PERM(J1,K1)
         ENDDO
         DO J1=J2(K1)+1,NPSIZE(K1)
            PERM2(J1,K1)=PERM(NPSIZE(K1)+J2(K1)-J1+1,K1)
         ENDDO
         PERM(1:NPSIZE(K1),K1)=PERM2(1:NPSIZE(K1),K1)
         NPERM(K1)=NPERM(K1)+1
         PRINT '(A,I6,A,I6)','addperm3> permutation ',NPERM(K1),' for group ',K1
         PRINT '(20I6)',PERM(1:NPSIZE(K1),K1)
!
! At this point we have the NPERM(K1)-1 th permutation of NPERM objects
! for group K1 in PERM(1:NPSIZE(K1),K1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (CHANGED) K1=1
         NDUMMY=1
         DO J3=1,NPERMGROUP
            DO J1=1,NPSIZE(J3)
               REALPERM(PERMGROUP(NDUMMY+J1-1))=PERMGROUP(NDUMMY+PERM(J1,J3)-1)
            ENDDO
            NDUMMY=NDUMMY+NPERMSIZE(J3)
         ENDDO

      NPERMTOTAL=NPERMTOTAL+1
333   CONTINUE 
C
C  Try this permutation.
C
222   CONTINUE
      PRINT '(A,I6,A,I6)','addperm3> Overall permutation number ',NPERMTOTAL,' inverse=',INVERSE
      PRINT '(20I6)',INVERSE*REALPERM(1:NATOMS)
      DO J1=1,NATOMS
         LPOINTS(3*(REALPERM(J1)-1)+1)=INVERSE*LPOINTSTS(3*(J1-1)+1)
         LPOINTS(3*(REALPERM(J1)-1)+2)=INVERSE*LPOINTSTS(3*(J1-1)+2)
         LPOINTS(3*(REALPERM(J1)-1)+3)=INVERSE*LPOINTSTS(3*(J1-1)+3)
      ENDDO

      tsloop: DO J3=1,NTS
!        PRINT '(A,I6,2F20.10)','J3,ETS(J3),ETS(DOTS)=',J3,ETS(J3),ETS(DOTS)
         IF (ABS(ETS(J3)-ETS(DOTS)).LT.EDIFFTOL) THEN
            READ(UTS,REC=J3) (LOCALPOINTS2(J4),J4=1,NOPT)
            CALL MINPERMDIST(LPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                    RMAT,.FALSE.)
!           PRINT '(A,2F20.10)','DISTANCE,GEOMDIFFTOL=',DISTANCE,GEOMDIFFTOL
            IF (DISTANCE.LT.GEOMDIFFTOL) THEN
               PRINT '(A,I6)','addperm3> permutation gives database TS ',J3
               GOTO 130
            ENDIF
         ENDIF
      ENDDO tsloop
      PRINT '(A,2I6,A,I6)','addperm3> permutation gives new TS ',NTS+1
      NTS=NTS+1
      IF (NTS.GT.MAXTS) CALL TSDOUBLE
      ETS(NTS)=ETS(DOTS)
      FVIBTS(NTS)=FVIBTS(DOTS)
      IF (FRICTIONT) NEGEIG(NTS)=NEGEIG(DOTS)
      HORDERTS(NTS)=1          ! not valid in general ! see next line
      CALL INERTIAWRAPPER(LPOINTS,NOPT,ANGLEAXIS,NIX,NIY,NIZ)
      IXTS(NTS)=NIX
      IYTS(NTS)=NIY
      IZTS(NTS)=NIZ
      WRITE(UTS,REC=NTS) (LPOINTS(J3),J3=1,NOPT)
      CALL FLUSH(UTS)
C
C  Identify the connected minimum.
C
      DO J3=1,NATOMS
         LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LPLUS(3*(J3-1)+1)
         LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LPLUS(3*(J3-1)+2)
         LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LPLUS(3*(J3-1)+3)
      ENDDO

      DO J3=1,NMIN
!        PRINT '(A,5I10)','addperm3> J3,DOTS,NTS,PLUS(DOTS),MINUS(DOTS)=',J3,DOTS,NTS,PLUS(DOTS),MINUS(DOTS)
         IF (ABS(EMIN(J3)-EMIN(PLUS(DOTS))).LT.EDIFFTOL) THEN
            READ(UMIN,REC=J3) (LOCALPOINTS2(J4),J4=1,NOPT)
            CALL MINPERMDIST(LPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                 RMAT,.FALSE.)
            IF (DISTANCE.LT.GEOMDIFFTOL) THEN
               WRITE(*,'(A,I6)') 'permutation for plus minimum gives old minimum ',J3
               PLUS(NTS)=J3
               IF (ENSEMBLE.EQ.'T') THEN
                  KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(J3)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1                 (FVIBMIN(J3)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J3))/TEMPERATURE
                  IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
                  ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(J3)  / (2*PI*HORDERTS(NTS))) +
     1             (FVIBMIN(J3)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J3)))
                  ELSE
                      KPLUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
               GOTO 11
            ENDIF
         ENDIF
      ENDDO
      WRITE(*,'(A,I6)') 'permutation for plus minimum gives new minimum ',NMIN+1
      NMIN=NMIN+1
      IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
      EMIN(NMIN)=EMIN(PLUS(DOTS))
      FVIBMIN(NMIN)=FVIBMIN(PLUS(DOTS))
      HORDERMIN(NMIN)=1          ! not valid in general !
      CALL INERTIAWRAPPER(LPOINTS,NOPT,ANGLEAXIS,NIX,NIY,NIZ)
      IXMIN(NMIN)=NIX
      IYMIN(NMIN)=NIY
      IZMIN(NMIN)=NIZ
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

      WRITE(UMIN,REC=NMIN) (LPOINTS(J3),J3=1,NOPT)
      CALL FLUSH(UMIN)
      IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
      CALL FLUSH(UMINDATA)
      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)

      PLUS(NTS)=NMIN
      IF (ENSEMBLE.EQ.'T') THEN
         KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1           (FVIBMIN(NMIN)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
         IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
      ELSE
         IF (TEMPERATURE.GT.ETS(NTS)) THEN
            KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(NMIN)  / (2*PI*HORDERTS(NTS))) +
     1      (FVIBMIN(NMIN)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
         ELSE
            KPLUS(NTS)=-1.0D250
         ENDIF
      ENDIF
      IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
      IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!     KSUM(NMIN)=0.0D0
!     IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KPLUS(NTS)
11    CONTINUE

      DO J3=1,NATOMS
         LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LMINUS(3*(J3-1)+1)
         LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LMINUS(3*(J3-1)+2)
         LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LMINUS(3*(J3-1)+3)
      ENDDO

      DO J3=1,NMIN
         IF (ABS(EMIN(J3)-EMIN(MINUS(DOTS))).LT.EDIFFTOL) THEN
            READ(UMIN,REC=J3) (LOCALPOINTS2(J4),J4=1,NOPT)
            CALL MINPERMDIST(LPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                 RMAT,.FALSE.)
            IF (DISTANCE.LT.GEOMDIFFTOL) THEN
               WRITE(*,'(A,I6)') 'permutation for minus minimum gives old minimum ',J3
               MINUS(NTS)=J3
               IF (ENSEMBLE.EQ.'T') THEN
                  KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(J3) / (2.0D0 * PI*HORDERTS(NTS))) +
     1               (FVIBMIN(J3) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J3))/TEMPERATURE
                  IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1             (FVIBMIN(J3) - FVIBTS(NTS))/2 + 
     2             (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J3)))
                  ELSE
                     KMINUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
               GOTO 12
            ENDIF
         ENDIF
      ENDDO
      WRITE(*,'(A,I5,A,I5)') 'permutation for minus minimum gives new minimum ',NMIN+1
      NMIN=NMIN+1
      IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
      EMIN(NMIN)=EMIN(MINUS(DOTS))
      FVIBMIN(NMIN)=FVIBMIN(MINUS(DOTS))
      HORDERMIN(NMIN)=1          ! not valid in general !
      CALL INERTIAWRAPPER(LPOINTS,NOPT,ANGLEAXIS,NIX,NIY,NIZ)
      IXMIN(NMIN)=NIX
      IYMIN(NMIN)=NIY
      IZMIN(NMIN)=NIZ
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

      WRITE(UMIN,REC=NMIN) (LPOINTS(J3),J3=1,NOPT)
      CALL FLUSH(UMIN)
      IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(UMINDATA,'(2F25.15,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN),IZMIN(NMIN)
      CALL FLUSH(UMINDATA)
      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
      MINUS(NTS)=NMIN
      IF (ENSEMBLE.EQ.'T') THEN
         KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN) / (2.0D0 * PI*HORDERTS(NTS))) +
     1         (FVIBMIN(NMIN) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
         IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
      ELSE
         IF (TEMPERATURE.GT.ETS(NTS)) THEN
            KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(NMIN) / (2*PI*HORDERTS(NTS))) +
     1       (FVIBMIN(NMIN) - FVIBTS(NTS))/2 + 
     2       (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
         ELSE
            KMINUS(NTS)=-1.0D250
         ENDIF
      ENDIF
!     IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KMINUS(NTS)
12    CONTINUE
C
C debug printing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
!     dump plus(NTS) coordinates from points.min, permuted plus coordinates,
!          ts coordinates of NTS, permuted minus coordinates, minus(nts) coordinates from points.min
!          in a file path.<nts>.xyz
!
!     WRITE(TSTRING,'(I10)') NTS
!     TSTRING='path.ts' // TRIM(ADJUSTL(TSTRING)) // '.xyz'

!     DO J3=1,NATOMS
!        LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LMINUS(3*(J3-1)+1)
!        LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LMINUS(3*(J3-1)+2)
!        LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LMINUS(3*(J3-1)+3)
!     ENDDO
!     LUNIT=GETUNIT()
!     OPEN(LUNIT,FILE=TRIM(ADJUSTL(TSTRING)),STATUS='UNKNOWN')
!     READ(UMIN,REC=PLUS(NTS)) (LOCALPOINTS2(J4),J4=1,NOPT)
!     WRITE(LUNIT,'(I10)') NATOMS
!     WRITE(LUNIT,'(A,2I10)') 'coordinates for nts,plus(nts)',NTS,PLUS(NTS)
!     WRITE(LUNIT,'(A2,1X,3F20.10)') ('LA ',LOCALPOINTS2(3*(J4-1)+1),LOCALPOINTS2(3*(J4-1)+2),LOCALPOINTS2(3*(J4-1)+3),J4=1,NATOMS)
!     WRITE(LUNIT,'(I10)') NATOMS
!     WRITE(LUNIT,'(A,2I10)') 'coordinates for permuted LPLUS:'
!     DO J3=1,NATOMS
!        LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LPLUS(3*(J3-1)+1)
!        LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LPLUS(3*(J3-1)+2)
!        LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LPLUS(3*(J3-1)+3)
!     ENDDO
!     WRITE(LUNIT,'(A2,1X,3F20.10)') ('LA ',LPOINTS(3*(J4-1)+1),LPOINTS(3*(J4-1)+2),LPOINTS(3*(J4-1)+3),J4=1,NATOMS)

!     WRITE(LUNIT,'(I10)') NATOMS
!     WRITE(LUNIT,'(A,20I6)') 'coordinates for ts with permutation ',INVERSE*REALPERM(1:NATOMS)
!     DO J1=1,NATOMS
!        LPOINTS(3*(REALPERM(J1)-1)+1)=INVERSE*LPOINTSTS(3*(J1-1)+1)
!        LPOINTS(3*(REALPERM(J1)-1)+2)=INVERSE*LPOINTSTS(3*(J1-1)+2)
!        LPOINTS(3*(REALPERM(J1)-1)+3)=INVERSE*LPOINTSTS(3*(J1-1)+3)
!     ENDDO
!     WRITE(LUNIT,'(A2,1X,3F20.10)') ('LA ',LPOINTS(3*(J4-1)+1),LPOINTS(3*(J4-1)+2),LPOINTS(3*(J4-1)+3),J4=1,NATOMS)
!     WRITE(LUNIT,'(I10)') NATOMS
!     WRITE(LUNIT,'(A,2I10)') 'coordinates for permuted LMINUS:'
!     DO J3=1,NATOMS
!        LPOINTS(3*(REALPERM(J3)-1)+1)=INVERSE*LMINUS(3*(J3-1)+1)
!        LPOINTS(3*(REALPERM(J3)-1)+2)=INVERSE*LMINUS(3*(J3-1)+2)
!        LPOINTS(3*(REALPERM(J3)-1)+3)=INVERSE*LMINUS(3*(J3-1)+3)
!     ENDDO
!     WRITE(LUNIT,'(A2,1X,3F20.10)') ('LA ',LPOINTS(3*(J4-1)+1),LPOINTS(3*(J4-1)+2),LPOINTS(3*(J4-1)+3),J4=1,NATOMS)
!     READ(UMIN,REC=MINUS(NTS)) (LOCALPOINTS2(J4),J4=1,NOPT)
!     WRITE(LUNIT,'(I10)') NATOMS
!     WRITE(LUNIT,'(A,2I10)') 'coordinates for nts,minus(nts)',NTS,MINUS(NTS)
!     WRITE(LUNIT,'(A2,1X,3F20.10)') ('LA ',LOCALPOINTS2(3*(J4-1)+1),LOCALPOINTS2(3*(J4-1)+2),LOCALPOINTS2(3*(J4-1)+3),J4=1,NATOMS)
!
!     CLOSE(LUNIT)
C
C end of debug printing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C

C
C  Only now do we know the connected minima for sure.
C
      IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
      IF (IMFRQT) THEN
         WRITE(UTSDATA,'(2F25.15,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
      ELSE
         WRITE(UTSDATA,'(2F25.15,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
      ENDIF
      CALL FLUSH(UTSDATA)
      IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
C
C  Update ts pointers.
C
      POINTERP(NTS)=-1
      POINTERM(NTS)=-1
      TOPPOINTER(PLUS(NTS))=NTS
      TOPPOINTER(MINUS(NTS))=NTS

      DO J3=NTS-1,1,-1
         IF (PLUS(J3).EQ.PLUS(NTS)) THEN
            POINTERP(NTS)=J3
            GOTO 41
         ELSE IF (MINUS(J3).EQ.PLUS(NTS)) THEN
            POINTERP(NTS)=J3
            GOTO 41
         ENDIF
      ENDDO
41    CONTINUE

      DO J3=NTS-1,1,-1
         IF (PLUS(J3).EQ.MINUS(NTS)) THEN
            POINTERM(NTS)=J3
            GOTO 42
         ELSE IF (MINUS(J3).EQ.MINUS(NTS)) THEN
            POINTERM(NTS)=J3
            GOTO 42
         ENDIF
      ENDDO
42    CONTINUE
C
C  ts pointers have been updated.
C
10    CONTINUE
130   CONTINUE

      IF (INVERSE.EQ.1) THEN
         INVERSE=-1
         PRINT '(A)','addperm3> Trying permutation plus inversion'
         GOTO 222
      ENDIF
      GOTO 116

20    CONTINUE
      DEALLOCATE(REALPERM,PERM,PERM2)
      DEALLOCATE(NPSIZE,J2,NPERM)
      PERMDIST=LLPERMDIST

      RETURN
      END
