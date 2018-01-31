C
C  Calculate rate constants by direct diagonalisation of the transition matrix.
C  
      SUBROUTINE DIAGONALISE2
      USE COMMONS
      USE PORFUNCS
      IMPLICIT NONE

      INTEGER J1, J2, M2
C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER IWORK(33*NMIN)
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*NMIN)
      INTEGER INFO
      DOUBLE PRECISION WORK(33*NMIN), DIAG(NMIN), DLAMCH, ZWORK(NMIN,NDIAGEIG)

      INTEGER MINMAP(NMIN), NAVAIL, NDEAD, FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1)
      DOUBLE PRECISION ABSTOL, LKSUM(NMIN), ELAPSED, TNEW
      DOUBLE PRECISION DUMMY
      DOUBLE PRECISION, ALLOCATABLE :: WP(:,:)
      LOGICAL DEADTS(NTS)

      INTEGER I, J

      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS
      CALL CPU_TIME(ELAPSED)
!
!  REGROUP and REGROUPFREE change NMINA, NMINB, LOCATIONA, LOCATIONB, so save the values and reset
!  to call GT more than once.
!  In fact, REGROUPFREE changes NMIN, NTS, etc. so we must stop after such a run or
!  do a complete reset somehow. This is done explicitly in routines like getppair
!  using the SAVESTATE module. Not needed here because we assume that diagonalise2 cannot be
!  called more than once!
!
      IF (REGROUPFREET.OR.REGROUPFREEABT) THEN
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
         CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
         TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                              ! before calling GETNCONN
         MAXBARRIER=HUGE(1.0D0)
      ENDIF

      IF (REGROUPPERSISTT) THEN
         ! CALL GETNCONN done in getbarrier2 called from persistence
         CALL PERSISTENCE
         TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                              ! before calling GETNCONN
         MAXBARRIER=HUGE(1.0D0)
      ENDIF

      IF (REGROUPRATET.OR.REGROUPPET) THEN
         CALL REGROUPFREE
         TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                              ! before calling GETNCONN
         MAXBARRIER=HUGE(1.0D0)
      ENDIF

      CALL GETNCONN
!
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!  
      CALL REGROUP(MINMAP)

      ALLOCATE(WP(NMIN,NMIN))

      CALL RATECONST_SETUP(LKSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)

!!!!!!!!!!!!!!!!!!!   set up transition matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Degenerate rearrangements are excluded.
!  
!  NCONN may have been set to zero for minima in a disconnected region. DEADTS could still
!  be false for these, so exclude them explicitly. 
!  Self-connection branching probability is zero.
!  
      WP(1:NMIN,1:NMIN)=0.0D0
      FROMLOOP: DO M2=1,NMIN
         IF (NCONN(M2).LE.NCONNMIN) CYCLE FROMLOOP
         J1=TOPPOINTER(M2)  !  sets J1 to the TS connected to minimum M2 with the highest id
         IF (J1.LE.0) CYCLE FROMLOOP
         DO WHILE (J1.GT.0) 
            IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
               IF (PLUS(J1).EQ.M2) THEN       !  M2 M1
                  WP(MINUS(J1),M2)=WP(MINUS(J1),M2)+EXP(KPLUS(J1))
               ELSEIF (MINUS(J1).EQ.M2) THEN  !  M1 M2
                  WP(PLUS(J1),M2)=WP(PLUS(J1),M2)+EXP(KMINUS(J1))
               ENDIF
            ENDIF
            IF (PLUS(J1).EQ.M2) THEN
               J1=POINTERP(J1)
            ELSE IF (MINUS(J1).EQ.M2) THEN
               J1=POINTERM(J1)
            ENDIF
         ENDDO
      ENDDO FROMLOOP

      DO J1=1,NMIN
         DO J2=1,NMIN
            IF (J2.NE.J1) WP(J1,J1)=WP(J1,J1)-WP(J2,J1)
         ENDDO
      ENDDO
      DO J1=1,NMIN
         DO J2=1,NMIN
            WP(J2,J1)=WP(J2,J1)*EXP((PFMIN(J1)-PFMIN(J2))/2.0D0)
         ENDDO
      ENDDO

      DO J1=1,NMIN
         DUMMY=0.0D0
         DO J2=1,NMIN
            DUMMY=DUMMY+WP(J2,J1)
            IF ((WP(J1,J2).LT.0.0D0).AND.(J1.NE.J2)) WRITE(*,'(A,2I6,D20.10)') 'WARNING, J1, J2, WP=',J1,J2,WP(J1,J2)
         ENDDO
         IF (ABS(DUMMY).GT.1.0D-10) WRITE(*,'(A,I4,A,D20.10)') 'diagonalise2> WARNING *** Matrix sum for row ',J1,' is ',DUMMY
C        WP(J1,J1)=WP(J1,J1)-DUMMY
      ENDDO

      DO I=1,NMIN
         DO J=1,NMIN
            WP(J,I)=WP(J,I)*DIAGSCALE
         ENDDO
      ENDDO

      ABSTOL=DLAMCH('Safe  minimum')
      NDIAGEIG=MIN(NMIN,NDIAGEIG)
      CALL DSYEVR('N','I','U',NMIN,WP,NMIN,0.0D0,1.0D0,NMIN-NDIAGEIG+1,NMIN,ABSTOL,NFOUND,DIAG,
     1            ZWORK,NMIN,ISUPPZ,WORK,LWORK,IWORK,ILWORK,INFO )
      WRITE(*,'(A,I6,A)') 'Lowest ',NDIAGEIG,' eigenvalues before rescaling:'
      WRITE(*,'(5G20.10)') (DIAG(J1),J1=1,NDIAGEIG)
      DIAG(1:NDIAGEIG)=DIAG(1:NDIAGEIG)/DIAGSCALE
      PRINT '(A)','diagonalise2> rates calculated from each eigenvalue assuming detailed balance holds:'
      DO J1=1,NDIAGEIG
         PRINT '(I6,A,G20.10,A,G20.10)',J1,' k(A<-B)=',-DIAG(J1)/(1.0D0+EXP(PFTOTALB-PFTOTALA)), 
     &                                     ' k(B<-A)=',-DIAG(J1)/(1.0D0+EXP(PFTOTALA-PFTOTALB))
      ENDDO

      CALL CPU_TIME(TNEW)
      TGT=TGT+TNEW-ELAPSED
      DEALLOCATE(WP)

      STOP

      END
