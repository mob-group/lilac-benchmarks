C LBFGS method for optimizing image points on great circle
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C EFK: Copied from mylbfgs.f, 3/11/06
C PROJECT, TWOEND, GRADSQ, INTMIN, CHARMM, DRAGT  options removed
C PV option retained
C 
      SUBROUTINE GCLBFGS(RS,RF,XI,N,M,X,DIAGCO,EPS,MFLAG,ENERGY,RMS,ITMAX,RESET,ITDONE,PTEST)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODHESS
      USE ZWK
      USE MODUNRES
      USE MODCHARMM
      use porfuncs
      IMPLICIT NONE
      INTEGER N,M,J1,ITMAX,ITDONE,NFAIL,NCOUNT,J2
C     DOUBLE PRECISION X(*),G(3*NATOMS),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,OVERLAP
      DOUBLE PRECISION X(N),G(N),SLENGTH,DDOT,OVERLAP,DISTF,RMAT(3,3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
      DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,RMS,ALPHA,GSAVE(N),DIST2
      LOGICAL DIAGCO, RESET, PTEST
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA
      DOUBLE PRECISION GLAST(N)
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN
      LOGICAL MFLAG
      DOUBLE PRECISION GSQSCALE, GSTHRESH, DOT1, DOT2
      INTEGER NSPECIAL, NALLOW, NINFO, FRAME
      LOGICAL PVFLAG
      COMMON /PVF/ PVFLAG
      COMMON /G2/ GSTHRESH, GSQSCALE, NSPECIAL, NALLOW, NINFO
      LOGICAL PATHT, DRAGT
      INTEGER NPATHFRAME, NDECREASE, ISTAT
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
      LOGICAL RADMOVED, OVERLP
      COMMON /DISCON/ RADMOVED, OVERLP
      LOGICAL PUSH, PULL
      COMMON /MORPHDATA/ PUSH, PULL
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
C
C  SGI appears to need this SAVE statement!
C
      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT
C
C EFK: specifically for great circle stuff
C
      DOUBLE PRECISION RS(3*NATOMS), RF(3*NATOMS), XI(GCIMAGE,3*NATOMS)
      DOUBLE PRECISION IMGE(GCIMAGE), STARTE, FINE, dEdX(N), RMST

      CALL POTENTIAL(RS, STARTE, dEdX, .TRUE., .FALSE., RMSt, .FALSE., .FALSE.)
      print*, 'START ENERGY:', STARTE
      CALL POTENTIAL(RF, FINE, dEdX, .TRUE., .FALSE., RMSt, .FALSE., .FALSE.)
      print*, 'FINISH ENERGY:', FINE

      KNOWE = .FALSE.; KNOWG = .FALSE.

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun
      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         PRINT '(A,I10,A,I10,A)', 'ERROR, dimension of W=',SIZE(W,1),
     $        ' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
         STOP
      ENDIF

      IF (N.NE.3*NATOMS+1) THEN
         PRINT*,'ERROR - N and 3*NATOMS are different in GCLBFGS: ',N,3*NATOMS
         STOP
      ENDIF

      WRITE(*,*) "sn402: Changed GCLBFGS so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please set FASTOVERLAP and check carefully that this part of the code is working as you expect."
      WRITE(*,*) "Then remove these messages!"
      CALL ALIGN_DECIDE(RS,RF,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DISTF,DIST2,RIGIDBODY,RMAT)

      ALPHA=1.0D0
      NFAIL=0
      FRAME=1
      IF (RESET) ITER=0
      ITDONE=0
      IF (RESET.AND.PTEST) WRITE(*,'(A)') ' Resetting LBFGS minimiser'
      IF ((.NOT.RESET).AND.PTEST) WRITE(*,'(A)')
     $     ' Not resetting LBFGS minimiser'
1     FIXIMAGE=.FALSE.
      IF (PV) THEN
         IF (.NOT.KNOWE) THEN 
            CALL GETENERG(RS, RF, X, 3*NATOMS, GCIMAGE, ENERGY,GSAVE,
     $           RMS, XI, IMGE, DEBUG)

! Check for cold fusion
            if (ENERGY.LT.coldFusionLimit) then
               ENERGY=1.0d6
               RMS=1.0d0
               WRITE(*,'(A)') ' Cold fusion diagnosed - step discarded'
               RETURN
            endif
         ENDIF
         PVFLAG=.FALSE.
         CALL PVOPT(X,ENERGY,GSAVE)
      ENDIF
      IF ((.NOT.KNOWE).OR.(.NOT.KNOWG)) THEN
         CALL GETENERG(RS, RF, X, 3*NATOMS, GCIMAGE, ENERGY, GSAVE, RMS,
     $        XI, IMGE, DEBUG)
         if (ENERGY.LT.coldFusionLimit) then
            ENERGY=1.0d6
            RMS=1.0d0
            WRITE(*,'(A)') ' Cold fusion diagnosed - step discarded'
            RETURN
         endif
      ENDIF

      G(1:N)=GSAVE(1:N)
      GLAST(1:N)=GSAVE(1:N)
C
      IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') 
     1             ' Energy per image and RMS force=',ENERGY/GCIMAGE, RMS,
     $     ' after ',ITDONE,' LBFGS steps'

10    CALL FLUSH(6)
      IF (DEBUG) THEN
         OPEN(UNIT=45,FILE='greatcircle.xyz',STATUS='UNKNOWN')
         WRITE(45,'(I6)') NATOMS
         WRITE(45,'(A)') ' '
         WRITE(45,'(A3,3G20.10)') ('LA ',RS(3*(J1-1)+1),RS(3*(J1-1)+2),RS(3*(J1-1)+3),J1=1,NATOMS)
         DO J1=1,GCIMAGE
            WRITE(45,'(I6)') NATOMS
            WRITE(45,'(A)') ' '
            WRITE(45,'(A3,3G20.10)') ('LA ',XI(J1,3*(J2-1)+1),XI(J1,3*(J2-1)+2),XI(J1,3*(J2-1)+3),J2=1,NATOMS)
         ENDDO
         WRITE(45,'(I6)') NATOMS
         WRITE(45,'(A)') ' '
         WRITE(45,'(A3,3G20.10)') ('LA ',RF(3*(J1-1)+1),RF(3*(J1-1)+2),RF(3*(J1-1)+3),J1=1,NATOMS)
         CLOSE(45)
         OPEN(UNIT=46, FILE='imgenergies.out', STATUS='UNKNOWN')
         WRITE(46, '(I3,1x,G20.10)') 0, STARTE
         DO J2=1,GCIMAGE
            WRITE(46, '(I3,1x,G20.10)') J2, IMGE(J2)
         ENDDO
         WRITE(46, '(I3,1x,G20.10)') GCIMAGE+1, FINE
         CLOSE(46)
      ENDIF
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         MFLAG=.TRUE.
C        PRINT*,'RMS,EPS,ITDONE,NSTEPMIN=',RMS,EPS,ITDONE,NSTEPMIN
         IF (ITDONE.LT.NSTEPMIN) MFLAG=.FALSE.
         IF (PV.AND.(.NOT.PVFLAG)) MFLAG=.FALSE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (PTEST) WRITE(*,'(1x,a,g25.17)') 'Final energy is ',ENERGY
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         FIXIMAGE=.FALSE.
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            PRINT*,'using estimate of the inverse diagonal elements'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
            DO I=1,N
               DIAG(I)=DGUESS
            ENDDO
         ENDIF
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M
         IYPT= ISPT+N*M
C
C  NR step for diagonal inverse Hessian
C
         
         DO I=1,N
            W(ISPT+I)= -G(I)*DIAG(I)
            W(I)= -G(I)*DIAG(I)
         ENDDO
         GNORM= DSQRT(DDOT(N,G,1,G,1))

C
C  Make the first guess for the step length cautious.
C
         IF (GNORM.EQ.0.0D0) THEN
            GNORM=1.0D0 ! exact zero is presumably wrong!
            PRINT '(A)','WARNING - GNORM was zero in mylbfgs,
     $           resetting to one'
         ENDIF
         STP=MIN(1.0D0/GNORM,GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  Update estimate of diagonal inverse Hessian elements
C  We divide by both YS and YY at different points, so
C  they had better not be zero!
C
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) YY=1.0D0
            DUMMY1=YS/YY
            DO I=1,N
               DIAG(I)=DUMMY1
            ENDDO
         ELSE
            PRINT*,'using estimate of the inverse diagonal elements'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
                  STOP
               ENDIF
            ENDDO
         ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS         
         DO I=1,N
            W(I)= -G(I)
         ENDDO
         CP= POINT
         DO I= 1,BOUND
            CP=CP-1
            IF (CP.EQ. -1)CP=M-1
            SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)= W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO I=1,N
            W(I)=DIAG(I)*W(I)
         ENDDO

         DO I=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0
      ENDIF
C
C  If this is a BFGSTST or MORPHT  run project out the uphill direction.
C     PRINT '(A,6G20.10)','W=',W(1:6)
C
C  Store the new search direction
      IF (ITER.GT.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
      ENDIF
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) THEN         
         OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)
      ENDIF
C
      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'Search direction has positive
     $        projection onto gradient - reversing step'
         DO I=1,N
            W(ISPT+POINT*N+I)= -W(I)  ! if we reverse the step it is important not to take the ABS value of YS/YY!
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF      
      DO I=1,N
         W(I)=G(I)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.GCMXSTP) STP=GCMXSTP/SLENGTH
C
C  We now have the proposed step.
C      
      GNORM= DSQRT(DDOT(N,G,1,G,1))      
C SOMETHING@S WRONG W/ STP!
      DO J1=1,N
         X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
C
      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.
C
C At this point we have new Cartesian or internal coordinates after taking a full
C or decreased step. The gradient is not known at this geometry.
C
      NDECREASE=0
      NCOUNT=0

 20   IF (PV) THEN
         CALL GETENERG(RS, RF, X, 3*NATOMS, GCIMAGE, ENEW, GSAVE, RMS, XI,
     $        IMGE, DEBUG)
C         CALL POTENTIAL(X,ENEW,GSAVE,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         PVFLAG=.FALSE.
         CALL PVOPT(X,ENEW,GSAVE)
      ENDIF
      CALL GETENERG(RS, RF, X, 3*NATOMS, GCIMAGE, ENEW, GSAVE, RMS, XI, IMGE, DEBUG)
      
      if (ENERGY.LT.coldFusionLimit) then
         ENERGY=1.0d6
         RMS=1.0d0
         WRITE(*,'(A)') ' Cold fusion diagnosed - step discarded'
         RETURN
      endif

      G(1:N)=GSAVE(1:N)

C
C  Must allow the energy to rise during a minimisation to allow for numerical noise or
C  systematic errors due to discontinuities or SCF convergence problems.
C
      IF ((ENEW-ENERGY.LE.MAXERISE).OR.PVTS.OR.DRAGT.OR.TWOENDS.OR.RADMOVED) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW

         IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,G13.5)')
     $        ' Energy per image and RMS force=',ENERGY/GCIMAGE,RMS,' after ',ITDONE,
     1           ' LBFGS steps, step:',STP*SLENGTH
C
C  Step finished so can reset OLDQ to new XINT, OLDCART to new CART,
C  as well as the Cartesian and internal gradients.
C         
         GLAST(1:N)=GSAVE(1:N) 
      ELSE 
C
C  Energy increased - try again with a smaller step size. Must cater for possible enormous
C  values of SLENGTH. Decreasing the step size doesn;t seem to help for CASTEP.
C
         IF (((ITER.GT.1).AND.(NDECREASE.GT.2)).OR.((ITER.LE.1).AND.(NDECREASE.GT.10)).OR.
     1              ((CASTEP.OR.ONETEP.OR.CP2K).AND.(NDECREASE.GT.1))) THEN 
            NFAIL=NFAIL+1
            IF (PTEST) PRINT*,' in mylbfgs LBFGS step
     $           cannot find a lower energy, NFAIL=',NFAIL
C
C  try resetting - go back to previous coordinates, ENERGY is not set to ENEW
C  we need to save the gradient corresponding to the last successful step
C              
            ITER=0            
            DO J1=1,N
               X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
               G(J1)=GLAST(J1)
            ENDDO
            IF (NFAIL.GT.NFAILMAX) THEN
               PRINT*,' Too many failures - give up'
               FIXIMAGE=.FALSE.
               RETURN
            ENDIF
            GOTO 30
         ENDIF
C
C  Try a smaller step.
C         
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         KNOWE=.FALSE.
         KNOWG=.FALSE.
         KNOWH=.FALSE.
         STP=STP/10.0D0
         NDECREASE=NDECREASE+1
         IF (PTEST) 
     1    WRITE(*,'(A,G19.10,A,G16.10,A,G15.8)')
     $        ' energy increased from ',ENERGY,' to ',ENEW,
     2            ' decreasing step to ',STP*SLENGTH
         GOTO 20
      ENDIF
C
C     Compute the new step and gradient change. Note that the step
C     length is accounted for when the step taken is saved.
C
30    NPT=POINT*N
      
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITER.GE.FIXAFTER)) FIXIMAGE=.TRUE.
      GOTO 10

      RETURN
      END

      SUBROUTINE GETENERG(RS, RF, V, NC, NI, TOTENERG, dTEdV, RMS, XI, ENERGS, VERBOSE)
C
C Subroutine to get the total energy of the image points on a particular 
C great circle and the derivatives of the energy wrt hypersphere center
C RS and RF are start and finish geometries
C V is the displacement of the hypersphere center from the linear center
C NC is the number of coordinates
C NI is the number of image points
C TOTENERG is the output total energy of the image points
C dTEdV is the gradient of the total energy wrt to hypersphere center
C RMS is the norm of dTEdV
C
      IMPLICIT NONE
      INTEGER NC, NI
      DOUBLE PRECISION RS(NC), RF(NC), V(NC), TOTENERG, dTEdV(NC+1)
      DOUBLE PRECISION ENERGS(NI)
      DOUBLE PRECISION XI(NI,NC), dXIdV(NI,NC,NC+1)
      DOUBLE PRECISION ENERG, dEdX(NC)
      DOUBLE PRECISION RMS
      INTEGER I, A, B
      LOGICAL VERBOSE
C ---- for testing ----
      DOUBLE PRECISION TOTENERG2, dTEdV2(NC+1)
      DOUBLE PRECISION XI2(NI,NC), dXIdV2(NI,NC,NC+1)
C      DOUBLE PRECISION T(NC), dTdV(NC, NC+1), T2(NC), dTdV2(NC, NC+1)
      DOUBLE PRECISION T, dTdV(NC+1), T2, dTdV2(NC+1)
      INTEGER J
      DOUBLE PRECISION TINY, testsum, RMSt
C
      TOTENERG = 0.0D0
      dTEdV = 0.0D0

C     Get each image point and derivatives
      CALL GETIMGDERVS(RS, RF, V, NC, NI, XI, dXIdV, VERBOSE)
C     Calculate energies and energy derivatives
      DO I=1,NI
         CALL POTENTIAL(XI(I,:), ENERGS(I), dEdX, .TRUE., .FALSE., RMS,
     $        .FALSE., .FALSE.)
         TOTENERG = TOTENERG+ENERGS(I)
         DO B = 1,NC+1
            DO A = 1,NC            
               dTEdV(B) = dTEdV(B) + dEdX(A)*dXIdV(I,A,B)
            ENDDO
         ENDDO
      ENDDO
      RMS = 0.0D0      
      DO A=1,NC+1
         RMS = RMS + dTEdV(A)**2
      ENDDO
      RMS = SQRT(RMS)
C ----TESTING ----
C Test image point derivatives
C      TINY = 0.000001
C      DO J = 1,NC+1
C         V(J) = V(J) + TINY
C     Get each image point and derivatives
C         CALL GETIMGDERVS(RS, RF, V, NC, NI, XI2, dXIdV2,T2, dTdV2)
C
C         DO I=1,NI
C            print*, J, (T2-T)/TINY, dTdV(J)
C            DO A=1,NC
C               print*, A, J, (T2(A)-T(A))/TINY, dTdV(A,J)
C            ENDDO
C            DO A = 1,NC
C               print*, A, J, (XI2(I,A)-XI(I,A))/TINY, dXIdV(I,A,J)
C            ENDDO
C         ENDDO
C                  
C         V(J) = V(J) - TINY
C      ENDDO
C      STOP
C ------------
C test energy derivatives wrt V
C      TINY = 0.0001
C      DO J = 1,NC
C         V(J) = V(J) + TINY
C         TOTENERG2=0.0D0
C     Get each image point and derivatives
C         CALL GETIMGDERVS(RS, RF, V, NC, NI, XI2, dXIdV2)
C         DO I=1,NI
C     Calculate energy and energy derivatives at this point
C            CALL POTENTIAL(XI(I:), ENERG, dEdX, .TRUE., .FALSE., RMS,
C     $           .FALSE., .FALSE.)
C            print*, 'ENERGY:', ENERG
C            TOTENERG2 = TOTENERG2+ENERG
C         ENDDO
C         print*, 'TOTE2, TOTE', TOTENERG2, TOTENERG
C         print*, J, (TOTENERG2-TOTENERG)/TINY, dTEdV(J)
C         V(J) = V(J) - TINY
C      ENDDO
C      STOP
C -------------
C compare to linear interpolation
C      DO I=1,NC
C         testc(I) = (rs(I)+rf(I))/2
C      ENDDO
C      
C      DO I=1,NC
C         print*, testc(I), XI(I)
C      ENDDO
C
      END SUBROUTINE GETENERG
C

      SUBROUTINE GETIMGDERVS(RS, RF, V, NC, NI, XI, dXIdV, VERBOSE)
C Subroutine to get the interpolation points Xi on a particular 
C great circle, and the derivatives dXia/dCb
C with respect to the hypersphere center.
C RS and RF are the start and finish geometries
C V is a particular shift of the hypersphere center 
C I gives the index of the interpolation point
C NC is the number of coordinates
C N is the number of interpolation points
C DXDC is the ouput matrix of derivatives

      IMPLICIT NONE
      LOGICAL VERBOSE
      INTEGER I, NC, NI
      DOUBLE PRECISION RS(NC), RF(NC), V(NC+1)
      DOUBLE PRECISION dXIdV(NI,NC,NC+1), XI(NI,NC)
      DOUBLE PRECISION TINY, S
      PARAMETER(TINY = 0.000001)
      DOUBLE PRECISION PI, TWOPI
      PARAMETER(PI=3.141592653589793D0,TWOPI=2.0D0*PI)
      DOUBLE PRECISION RSF(NC), RSF2, RSFR, VRSF, VV
      DOUBLE PRECISION dCdV(NC,NC+1), dROdV(NC+1), dV1dV(NC,NC+1)
      DOUBLE PRECISION dTFdV(NC+1), dV2dV(NC,NC+1), FV1, V22, V2R
      DOUBLE PRECISION dV2RdV(NC+1)
      DOUBLE PRECISION C(NC), F(NC), RO, RO2, ROR, V1(NC), TF, V2(NC)
      DOUBLE PRECISION dFdV(NC, NC+1)
      DOUBLE PRECISION DUMMY, CTF, STF
      DOUBLE PRECISION STI, CTI, X, Y, dX(NC+1), dY(NC+1)
      INTEGER A,B,J
      DOUBLE PRECISION COFFSET
C for testing ----
      DOUBLE PRECISION testsum, testsum2, DDOT
      
      S = V(NC+1) ! the scalar by which center offset is scaled

C Get difference vector RSF, RSF^2, V^2, and V*RSF
      RSF2 = 0.0D0
      VV = 0.0D0
      VRSF = 0.0D0
      DO J=1,NC
         RSF(J) = RS(J)-RF(J)
         RSF2 = RSF2 + RSF(J)**2
         VV = VV + V(J)**2
         VRSF = VRSF + V(J)*RSF(J)
      ENDDO
      RSFR = 1/SQRT(RSF2)

C Get hypersphere center and its derivatives
      COFFSET = 0.0D0
      DO A=1,NC         
         F(A) = V(A) - VRSF*RSF(A)/RSF2
         C(A) = (RS(A)+RF(A))/2 + S*F(A)
         COFFSET=COFFSET+(F(A))**2
      ENDDO
      IF(COFFSET.LT.1.0D-6) print*, 'SMALL CENTER OFFSET, w/o S:', COFFSET
      DO A=1,NC
         DO B= 1,NC
            dFdV(A,B) = -RSF(A)*RSF(B)/RSF2            
         ENDDO
         dFdV(A,A) = dFdV(A,A) + 1
      ENDDO
      dCdV(1:NC,1:NC) = S*dFdV(1:NC,1:NC)
      DO A = 1,NC
         dFdV(A, NC+1) = 0.0D0
         dCdV(A, NC+1) = V(A) - RSF(A)/RSF2*VRSF
      ENDDO

C Get radius of circle & its derivatives
      RO2 = 0.0D0
      DO A=1,NC
         RO2 = RO2 + (RS(A)-C(A))**2
      ENDDO
      RO = SQRT(RO2)
      ROR = 1/RO

      dROdV = 0.0D0
      DO B = 1, NC+1
         DO A=1,NC
            dROdV(B) = dROdV(B) + ROR*(C(A)-RS(A))*dCdV(A,B)
         ENDDO
      ENDDO

C Get V1 & its derivatives
      DO A=1,NC
         V1(A) = (RS(A)-C(A))*ROR
      ENDDO
      DO A=1,NC
         DUMMY = (RS(A)-C(A))/RO2
         DO B=1,NC+1
            dV1dV(A,B) = -ROR*dCdV(A,B) - DUMMY*dROdV(B)
         ENDDO
      ENDDO

C Get V2 and its derivatives
      FV1 = 0.0D0 ! F*V1
      DO A=1,NC
         FV1 = FV1 + F(A)*V1(A)
      ENDDO

      V22 = 0.0D0
      DO A=1,NC
         V2(A) = F(A) - FV1*V1(A)
         V22 = V22 + V2(A)**2
         DO B = 1,NC+1
            dV2dV(A,B) = dFdV(A,B)-FV1*dV1dV(A,B)
            DO J = 1,NC
               dV2dV(A,B) = dV2dV(A,B) -
     $              (F(J)*dV1dV(J,B)+dFdV(J,B)*V1(J))*V1(A)
            ENDDO
         ENDDO
      ENDDO

      V2R = 1/SQRT(V22)

C get derivatives of V2R
      DO B=1,NC+1
         dV2RdV(B) = 0.0D0
         DO A = 1,NC
            dV2RdV(B) = dV2RdV(B) + V2(A)*dV2dV(A,B)*V2R
         ENDDO
         dV2RdV(B) = - dV2RdV(B) * V2R * V2R
      ENDDO

C      testsum = 0.0D0
C      testsum2 = 0.0D0
C      DO A=1,NC
C         testsum = testsum + (RF(A)-C(A))**2
C         testsum2 = testsum2 + (V1(A))**2
C      ENDDO

C Get TF and its derivatives
      X = 0.0D0
      Y = 0.0D0
      DO A = 1,NC
         X = X + (RF(A)-C(A))*V1(A)
         Y = Y + (RF(A)-C(A))*V2(A)*V2R
      ENDDO
 
      CTF = X*ROR ! cos(T_N+1)
      STF = Y*ROR ! sin(T_N+1)

C Deal with numerical errors which make abs(CTF) > 1
      IF (CTF.GT.1) THEN
         CTF = 1.0D0
      ELSE IF (CTF.LT.-1) THEN
         CTF = -1.0D0
      ENDIF

      IF (STF.GE.0) THEN
         TF = ACOS(CTF)
      ELSE IF (STF.LT.0) THEN
         TF = TWOPI-ACOS(CTF)
      ENDIF
     
      IF (VERBOSE) print '(A,1x,F10.5,1x,F10.5,1x,F10.5)',' greatcircle> Angle, Radius:', TF, RO

      DO B = 1,NC+1
         dX(B) = 0.0D0
         dY(B) = 0.0D0
         DO A = 1,NC
            dX(B) = dX(B)   -dCdV(A,B)*V1(A) + (RF(A)-C(A))*dV1dV(A,B)
            dY(B) = dY(B) + (-dCdV(A,B)*V2(A) + (RF(A)-C(A))*dV2dV(A,B)) * V2R
     &           + (RF(A)-C(A))*V2(A)*dV2RdV(B)
         ENDDO

         dTFdV(B) = (X*dY(B)-Y*dX(B))/RO2
      ENDDO

      DO I=1,NI
C Get image points Xi      
         STI = SIN(I*TF / (NI+1))
         CTI = COS(I*TF / (NI+1))

         DO A = 1,NC            
            XI(I,A) = C(A) + V1(A)*CTI*RO + V2(A)*STI*RO*V2R
         ENDDO

C Get derivatives of Xi
         DO A = 1,NC
            DO B= 1, NC+1
               dXIdV(I,A,B) = dCdV(A,B) - I*V1(A)/(NI+1)*STI*dTFdV(B)*RO +
     $              I*V2(A)/(NI+1)*CTI*dTFdV(B)*RO*V2R + 
     $              dV1dV(A,B)*CTI*RO + dV2dV(A,B)*STI*RO*V2R +
     $              V1(A)*CTI*dROdV(B)+V2(A)*STI*dROdV(B)*V2R +
     $              V2(A)*RO*STI*dV2RdV(B)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE GETIMGDERVS
