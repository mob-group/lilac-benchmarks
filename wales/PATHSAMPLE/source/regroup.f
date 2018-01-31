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

      SUBROUTINE REGROUP(MINMAP)
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2
      INTEGER MINMAP(NMIN), NEWNMINA, NEWNMINB, BASIN(NMIN), NBASIN, NCOUNTA, NCOUNTB, NCOUNT, INVMAP(NMIN), XNCONN(NMIN)
      LOGICAL ISA(NMIN), ISB(NMIN), CHANGED, BASINA(NMIN), BASINB(NMIN)
      DOUBLE PRECISION XEMIN(NMIN), XFVIBMIN(NMIN), XIXMIN(NMIN), XIYMIN(NMIN), 
     &                 XIZMIN(NMIN), XPFMIN(NMIN), XGPFOLD(NMIN), LNPROD
      INTEGER XHORDERMIN(NMIN), XPLUS(NTS), XMINUS(NTS), NCONNDUMM(NMIN), INDEX(NMIN)

      DO J1=1,NMIN
         MINMAP(J1)=J1
      ENDDO
C
C  We always want to do the reordering because it makes GT more efficient.
C  Just set REGROUPTHRESH to a huge negative number in this case.
C
      IF (.NOT.REGROUPT) REGROUPTHRESH=-HUGE(1.0D0)
C
C  Add minima to the A and B sets according to REGROUPTHRESH.
C  We always want to do this, because it also reorders the minima A, B, I,
C  which speeds things up.
C  Two ways have been tried for this: using P^fold values read in from a previous
C  run, and via an energy threshold and a disconnectivity graph analysis.
C  Reorganise the database and rates accordingly.
C
      ISA(1:NMIN)=.FALSE.
      ISB(1:NMIN)=.FALSE.
      DO J1=1,NMINA
         ISA(LOCATIONA(J1))=.TRUE.
      ENDDO
      DO J1=1,NMINB
         ISB(LOCATIONB(J1))=.TRUE.
      ENDDO
      NEWNMINA=NMINA
      NEWNMINB=NMINB
C
C  Do a superbasin analysis at energy REGROUPTHRESH to reassign A and B minima.
C  The objective is to stop the regrouping producing a structure that ruins the KMC.
C  Minima in their own superbasin at the given threshold will have BASIN(J1)=0.
C
      DO J1=1,NMIN
         BASIN(J1)=0
      ENDDO
      NBASIN=0
      DO 
         CHANGED=.FALSE.
         DO J1=1,NTS
            IF (ETS(J1).LT.REGROUPTHRESH) THEN
               IF ((BASIN(PLUS(J1)).EQ.0).AND.(BASIN(MINUS(J1)).EQ.0)) THEN
                  CHANGED=.TRUE.
                  NBASIN=NBASIN+1
                  BASIN(PLUS(J1))=NBASIN
                  BASIN(MINUS(J1))=NBASIN
               ELSEIF (BASIN(PLUS(J1)).NE.BASIN(MINUS(J1))) THEN
                  CHANGED=.TRUE.
                  IF (BASIN(PLUS(J1)).EQ.0) THEN
                     BASIN(PLUS(J1))=BASIN(MINUS(J1))
                  ELSEIF (BASIN(MINUS(J1)).EQ.0) THEN
                     BASIN(MINUS(J1))=BASIN(PLUS(J1))
                  ELSE
                     BASIN(PLUS(J1))=MIN(BASIN(PLUS(J1)),BASIN(MINUS(J1)))
                     BASIN(MINUS(J1))=BASIN(PLUS(J1))
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF (.NOT.CHANGED) EXIT
      ENDDO
      BASINA(1:NMIN)=.FALSE.
      DO J1=1,NMIN
         IF (ISA(J1).AND.(BASIN(J1).GT.0)) BASINA(BASIN(J1))=.TRUE.
      ENDDO
      DO J1=1,NMIN
         IF (BASIN(J1).GT.0) THEN
            IF (BASINA(BASIN(J1)).AND.(.NOT.ISA(J1))) THEN
               ISA(J1)=.TRUE.
               NEWNMINA=NEWNMINA+1
               IF (DEBUG) PRINT '(A,I6,A)','regroup> minimum ',J1,' added to A set'
            ENDIF
         ENDIF
      ENDDO
C
C  Now do the B set.
C
!     DO J1=1,NMIN
!        BASIN(J1)=J1
!     ENDDO
!     DO 
!        CHANGED=.FALSE.
!        DO J1=1,NTS
!           IF (ETS(J1).LT.REGROUPTHRESH) THEN
!              IF (BASIN(PLUS(J1)).NE.BASIN(MINUS(J1))) THEN
!                 CHANGED=.TRUE.
!                 BASIN(PLUS(J1))=MIN(BASIN(PLUS(J1)),BASIN(MINUS(J1)))
!                 BASIN(MINUS(J1))=BASIN(PLUS(J1))
!              ENDIF
!           ENDIF
!        ENDDO
!        IF (.NOT.CHANGED) EXIT
!     ENDDO

      BASINB(1:NMIN)=.FALSE.
      DO J1=1,NMIN
         IF (ISB(J1).AND.(BASIN(J1).GT.0)) BASINB(BASIN(J1))=.TRUE.
      ENDDO
      DO J1=1,NMIN
         IF (BASIN(J1).GT.0) THEN
            IF (BASINB(BASIN(J1)).AND.(.NOT.ISB(J1))) THEN
               ISB(J1)=.TRUE.
               NEWNMINB=NEWNMINB+1
               IF (DEBUG) PRINT '(A,I6,A)','regroup> minimum ',J1,' added to B set'
            ENDIF
         ENDIF
      ENDDO
      DO J1=1,NMIN
         IF (ISA(J1).AND.ISB(J1)) THEN
            PRINT '(A,I6,A)','regroup> ERROR - minimum ',J1,' belongs to A and B sets'
            STOP
         ENDIF
      ENDDO
C
C  End superbasin analysis.
C
      NCOUNTA=0
      NCOUNTB=0
      NCOUNT=0
C
C  Set up MINMAP(J1), the location of minimum J1 in the original scheme.
C  Set up INVMAP(J1), the location that minimum J1 in the original scheme maps to.
C  Even if the number of A and B minima did not change there may be a reordering
C  into NMINA A, then NMINB B, then the I minima. LOCATIONA and LOCATIONB should
C  not be used hereafter unless they are saved and restored.
C
C  We also sort the I minima from most to least connected. This speeds up the GT calculation
C  significantly.
C
      DO J1=1,NMIN
         IF (ISA(J1).OR.ISB(J1)) THEN
            NCONNDUMM(J1)=-100
         ELSE
            NCONNDUMM(J1)=NCONN(J1)
         ENDIF
         INDEX(J1)=J1
      ENDDO
      CALL SORT2(NMIN,NMIN,NCONNDUMM,INDEX)

      DO J1=1,NMIN
         IF (ISA(J1)) THEN
            NCOUNTA=NCOUNTA+1
            MINMAP(NCOUNTA)=J1
            INVMAP(J1)=NCOUNTA
         ELSEIF (ISB(J1)) THEN
            NCOUNTB=NCOUNTB+1
            MINMAP(NEWNMINA+NCOUNTB)=J1
            INVMAP(J1)=NEWNMINA+NCOUNTB
         ELSE
            NCOUNT=NCOUNT+1
            MINMAP(NEWNMINA+NEWNMINB+NCOUNT)=INDEX(NCOUNT)
            INVMAP(INDEX(NCOUNT))=NEWNMINA+NEWNMINB+NCOUNT
!           MINMAP(NEWNMINA+NEWNMINB+NCOUNT)=J1
!           INVMAP(J1)=NEWNMINA+NEWNMINB+NCOUNT
         ENDIF
      ENDDO
      IF (.NOT.SHANNONT) PRINT '(2(A,I8),A,G20.10)','regroup> regrouping gives ',NEWNMINA,' A min and ',NEWNMINB,
     &                        ' B min for a threshold of ',REGROUPTHRESH
      DO J1=1,NMIN
         XEMIN(J1)=EMIN(MINMAP(J1)) 
         XFVIBMIN(J1)=FVIBMIN(MINMAP(J1))
         XIXMIN(J1)=IXMIN(MINMAP(J1))
         XIYMIN(J1)=IYMIN(MINMAP(J1))
         XIZMIN(J1)=IZMIN(MINMAP(J1))
         XHORDERMIN(J1)=HORDERMIN(MINMAP(J1))
         XPFMIN(J1)=PFMIN(MINMAP(J1))
         XNCONN(J1)=NCONN(MINMAP(J1))
         XGPFOLD(J1)=GPFOLD(MINMAP(J1))
      ENDDO
      DO J1=1,NMIN
         EMIN(J1)=XEMIN(J1)
         FVIBMIN(J1)=XFVIBMIN(J1)
         IXMIN(J1)=XIXMIN(J1)
         IYMIN(J1)=XIYMIN(J1)
         IZMIN(J1)=XIZMIN(J1)
         HORDERMIN(J1)=XHORDERMIN(J1)
         PFMIN(J1)=XPFMIN(J1)
         NCONN(J1)=XNCONN(J1)
         GPFOLD(J1)=XGPFOLD(J1)
      ENDDO
      DO J1=1,NTS
         XPLUS(J1)=INVMAP(PLUS(J1))
         XMINUS(J1)=INVMAP(MINUS(J1))
      ENDDO
      DO J1=1,NTS
         PLUS(J1)=XPLUS(J1)
         MINUS(J1)=XMINUS(J1)
      ENDDO
      NMINA=NEWNMINA
      NMINB=NEWNMINB
C
C  Now we basically repeat some of setup.
C
      DO J1=1,NMIN
         TOPPOINTER(J1)=-1
      ENDDO
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
      PFTOTALA=0.0D0
      DO J1=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(J1))
      ENDDO
      PFTOTALA=LOG(PFTOTALA)
      PFTOTALB=0.0D0
      DO J1=NMINA+1,NMINA+NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(J1))
      ENDDO
      PFTOTALB=LOG(PFTOTALB)
C
C  LOCATIONA and LOCATIONB are used in GT, so we need to reset them.
C
      DEALLOCATE(LOCATIONA,LOCATIONB)
      ALLOCATE(LOCATIONA(NMINA),LOCATIONB(NMINB))
      DO J1=1,NMINA
         LOCATIONA(J1)=J1
      ENDDO
      DO J1=1,NMINB
         LOCATIONB(J1)=NMINA+J1
      ENDDO
C
C  If we are going to analyse the min.data.regrouped.resorted and ts.data.regrouped.resorted
C  files for rates subsequently, then we have to arrange for the ln products of frequencies
C  to give us a factor of (kT/h). This can be done by setting the ln product equal to zero
C  for the transition state and 2 * ln(2*Pi*k*T/h) for all the minima. We already have h in the
C  units of kT, so this is easy. The 2*Pi factor occurs because the frequencies are assumed to be 
C  angular normal mmode frequencies, and the factor of two occurs because they are assumed
C  to be squared.
C
C  For a microcanonical ensemble the factor should be 1/h.
C
      IF (ENSEMBLE.EQ.'T') THEN
         LNPROD=2.0D0*LOG(2.0D0*3.141592654D0*TEMPERATURE/PLANCK)
      ELSE
         LNPROD=2.0D0*LOG(2.0D0*3.141592654D0/PLANCK)
      ENDIF
      OPEN(UNIT=1,FILE='min.data.regrouped.resorted',STATUS='UNKNOWN')
      DO J1=1,NMIN
         WRITE(1,'(2G20.10,I6,4F20.10)') EMIN(J1),LNPROD,1,1.0,1.0,1.0,0.0
      ENDDO
      CLOSE(1)
      OPEN(UNIT=1,FILE='ts.data.regrouped.resorted',STATUS='UNKNOWN')
      DO J1=1,NTS
         WRITE(1,'(2G20.10,3I10,3F20.10)') ETS(J1),0.0,1,PLUS(J1),MINUS(J1),1.0,1.0,1.0
      ENDDO
      CLOSE(1)
      OPEN(UNIT=1,FILE='min.A.regrouped.resorted',STATUS='UNKNOWN')
      WRITE(1,'(I6)') NMINA
      DO J1=1,NMINA
         WRITE(1,'(I6)') LOCATIONA(J1)
      ENDDO
      CLOSE(1)
      OPEN(UNIT=1,FILE='min.B.regrouped.resorted',STATUS='UNKNOWN')
      WRITE(1,'(I6)') NMINB
      DO J1=1,NMINB
         WRITE(1,'(I6)') LOCATIONB(J1)
      ENDDO
      CLOSE(1)
      IF (.NOT.SHANNONT) PRINT '(A)','regroup> WARNING - regroup reorders the minima A, then B, then I'

      RETURN
      END
