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
      SUBROUTINE OLDNEB(ITEST,PTEST,EMAX,LGMAX,PERM,Q)
      USE COMMONS
      USE MODTWOEND
      USE MODNEB
      USE KEY
      USE MODHESS
      use porfuncs
      IMPLICIT NONE
      INTEGER J1, J2, ITDONE, JMAX, NMAGNIFY, NIMAGESAVE
      DOUBLE PRECISION RMS, EMAX, DVEC(NIMAGEMAX+1), TIME, TIME0, Q(3*NATOMS),
     1                 LBWORK(3*NATOMS*NIMAGEMAX*(2*MUPDATE+1)+2*MUPDATE), RMSREAL(NIMAGEMAX), EREAL(NIMAGEMAX),
     2                 POINTS(5*NATOMS*NIMAGEMAX), SEPARATION, EINITIAL, EFINAL, RMSINITIAL, RMSFINAL,
     3                 DIAG(3*NATOMS*NIMAGEMAX), GRAD(3*NATOMS), MIDDLE(3*NATOMS), 
     5                 CMX,CMY,CMZ,OVEC(3),H1VEC(3),H2VEC(3),EFSAVE,EISAVE,
     6                 EFORIG, EIORIG, FINORIG(3*NATOMS),LGMAX
      LOGICAL MFLAG, RESET, ITEST, PTEST, PERM
      COMMON /DVNEB/ DVEC
      COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
      CHARACTER FNAME*80
      INTEGER IADD, LNOPT, LNATOMS
      COMMON /INTS/ IADD, LNOPT, LNATOMS
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

!
! The declaration of POINTS provides more space than is actually needed for TIP potential.
! 9*NATOMS*NIMAGEMAX/2 is needed in reality.
!
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
         IADD=3*(NATOMS/2)
         LNOPT=NOPT/2      ! we only want to optimise the angular degrees of freedom for TIP
         LNATOMS=NATOMS/2  ! for TIP potentials both centre-of-mass and Euler angles are now in odata
      ELSE
         IADD=0
         LNOPT=NOPT
         LNATOMS=NATOMS
      ENDIF
      NMAGNIFY=0
      IF (NIMAGE.GT.NIMAGEMAX) THEN
         PRINT*,' Too many images requested - quit'
         STOP
      ENDIF
      CALL MYCPU_TIME(TIME0,.TRUE.)
!
!  End points should already be in optimal alignment
!
!     IF (ZSYM(1)(1:1).NE.'W') THEN
!        CMX=0.0D0
!        CMY=0.0D0
!        CMZ=0.0D0
!        DO J1=1,LNATOMS
!           CMX=CMX+Q(3*(J1-1)+1)
!           CMY=CMY+Q(3*(J1-1)+2)
!           CMZ=CMZ+Q(3*(J1-1)+3)
!        ENDDO
!        CMX=CMX/LNATOMS
!        CMY=CMY/LNATOMS
!        CMZ=CMZ/LNATOMS
!        DO J1=1,LNATOMS
!           Q(3*(J1-1)+1)=Q(3*(J1-1)+1)-CMX
!           Q(3*(J1-1)+2)=Q(3*(J1-1)+2)-CMY
!           Q(3*(J1-1)+3)=Q(3*(J1-1)+3)-CMZ
!        ENDDO
!        CMX=0.0D0
!        CMY=0.0D0
!        CMZ=0.0D0
!        DO J1=1,LNATOMS
!           CMX=CMX+FIN(3*(J1-1)+1)
!           CMY=CMY+FIN(3*(J1-1)+2)
!           CMZ=CMZ+FIN(3*(J1-1)+3)
!        ENDDO
!        CMX=CMX/LNATOMS
!        CMY=CMY/LNATOMS
!        CMZ=CMZ/LNATOMS
!        DO J1=1,LNATOMS
!           FIN(3*(J1-1)+1)=FIN(3*(J1-1)+1)-CMX
!           FIN(3*(J1-1)+2)=FIN(3*(J1-1)+2)-CMY
!           FIN(3*(J1-1)+3)=FIN(3*(J1-1)+3)-CMZ
!        ENDDO
!     ENDIF
      RMSINITIAL=0.0D0
      RMSFINAL=0.0D0
!     PRINT *,'Q:'
!     PRINT '(6F15.5)',Q(1:3*NATOMS)
!     PRINT *,'FIN:'
!     PRINT '(6F15.5)',FIN(1:3*NATOMS)
      IF (ITEST) CALL POTENTIAL(Q,EINITIAL,GRAD,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
      IF (ITEST) CALL POTENTIAL(FIN,EFINAL,GRAD,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
      IF (ITEST) WRITE(*,'(A,F20.10,A,F15.6,A,F20.10,A,F15.6)') 
     1           ' Initial energy=',EINITIAL,' RMS force=',RMSINITIAL,' final energy=',EFINAL,' RMS=',RMSFINAL
      WRITE(*,'(A,F20.10,A,F15.6,A,F20.10,A,F15.6)') 
     1           ' Initial energy=',EINITIAL,' RMS force=',RMSINITIAL,' final energy=',EFINAL,' RMS=',RMSFINAL
      EIORIG=EINITIAL
      EFORIG=EFINAL
      DO J1=1,LNOPT+IADD
         FINORIG(J1)=FIN(J1)
      ENDDO
C
C  Return to this point, possibly reinitialising EFINAL, EINITIAL, Q and
C  FIN, to increase the resolution.
C
10    IF (NMAGNIFY.EQ.0) THEN
         CALL MAKEIMAGE(POINTS,PERM,ITEST,0,Q)
      ELSE
C
C  Interpolate between image JMAX-1 and JMAX - we need (NIMAGE-1)/2
C  images, so NIMAGE ought to be odd!
C
         DO J1=1,LNOPT+IADD
            MIDDLE(J1)=POINTS((LNOPT+IADD)*(JMAX-1)+J1)
         ENDDO
         IF (JMAX.GT.1) THEN  ! if jmax is the first image, Q and EINITIAL stay the same
            DO J1=1,LNOPT+IADD
               Q(J1)=POINTS((LNOPT+IADD)*(JMAX-2)+J1)
            ENDDO
            EINITIAL=EREAL(JMAX-1)
         ENDIF
         EISAVE=EINITIAL
         DO J1=1,LNOPT+IADD
            FIN(J1)=MIDDLE(J1)
         ENDDO
         EFSAVE=EFINAL
         EFINAL=EREAL(JMAX)
         NIMAGESAVE=NIMAGE
         NIMAGE=(NIMAGE-1)/2
         CALL MAKEIMAGE(POINTS,PERM,ITEST,0,Q)
         EFINAL=EFSAVE
C
C  Interpolate between images JMAX and JMAX+1.
C
         IF (JMAX.LT.NIMAGESAVE) THEN  ! if jmax is the last image, FIN and EFINAL stay the same
            DO J1=1,LNOPT+IADD
               FIN(J1)=POINTS((LNOPT+IADD)*(JMAX)+J1)
            ENDDO
            EFINAL=EREAL(JMAX+1)
         ENDIF
         DO J1=1,LNOPT+IADD
            Q(J1)=MIDDLE(J1)
         ENDDO
         EINITIAL=EREAL(JMAX)
         CALL MAKEIMAGE(POINTS,PERM,ITEST,NIMAGE+1,Q)
         DO J1=1,LNOPT+IADD
            POINTS((LNOPT+IADD)*NIMAGE+J1)=MIDDLE(J1)
         ENDDO
         NIMAGE=NIMAGESAVE
         EINITIAL=EISAVE
      ENDIF
C     PRINT*,'NIMAGE after MAKEIMAGE=',NIMAGE
C     PRINT*,'EINITIAL,EFINAL=',EINITIAL,EFINAL
      NMAGNIFY=NMAGNIFY+1
C
C  Minimisation subject to spring constraints
C
      RESET=.TRUE.
      IF (FILTH.EQ.0) THEN
         WRITE(FNAME,'(A8)') 'EofS.neb'
      ELSE 
         WRITE(FNAME,'(A)') 'EofS.neb.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF
      OPEN(UNIT=4,FILE=FNAME,STATUS='UNKNOWN')
      CALL NEBBFGS((LNOPT+IADD)*NIMAGE,MUPDATE,POINTS,.FALSE.,DIAG,RMSNEB,LBWORK,MFLAG,EREAL,RMSREAL,
     1             NSTEPNEB,RESET,ITDONE,PTEST,Q,LGMAX)
         WRITE(4,'(3F20.10)') 0.0D0,EINITIAL,RMSINITIAL
         WRITE(4,'(3F20.10)') (DVEC(J1),EREAL(J1),RMSREAL(J1),J1=1,NIMAGE)
         WRITE(4,'(3F20.10)') DVEC(NIMAGE+1),EFINAL,RMSFINAL
      CLOSE(4)
C
C  Dump pathway as an xyz file
C
      IF (FILTH.EQ.0) THEN
         WRITE(FNAME,'(A12)') 'neb.path.xyz'
      ELSE 
         WRITE(FNAME,'(A)') 'neb.path.xyz.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF
      OPEN(UNIT=3,FILE=FNAME,STATUS='UNKNOWN')
      IF (ZSYM(1)(1:1).EQ.'W') THEN
         WRITE(3,'(I6)') 3*LNATOMS
         WRITE(3,'(F20.10)') EINITIAL
         DO J1=1,LNATOMS
            CALL CONVERT(Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),
     1       Q(3*(LNATOMS+J1-1)+1),Q(3*(LNATOMS+J1-1)+2),Q(3*(LNATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
            WRITE(3,'(A2,4X,3F20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
            WRITE(3,'(A2,4X,3F20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
            WRITE(3,'(A2,4X,3F20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
         ENDDO
         DO J1=1,NIMAGE
            WRITE(3,'(I6)') 3*LNATOMS
            WRITE(3,'(F20.10)') EREAL(J1)
            DO J2=1,LNATOMS
               CALL CONVERT(POINTS((LNOPT+IADD)*(J1-1)+3*(J2-1)+1),
     1                      POINTS((LNOPT+IADD)*(J1-1)+3*(J2-1)+2),
     1                      POINTS((LNOPT+IADD)*(J1-1)+3*(J2-1)+3),
     1                      POINTS((LNOPT+IADD)*(J1-1)+3*(LNATOMS+J2-1)+1),
     1                      POINTS((LNOPT+IADD)*(J1-1)+3*(LNATOMS+J2-1)+2),
     1                      POINTS((LNOPT+IADD)*(J1-1)+3*(LNATOMS+J2-1)+3),OVEC,H1VEC,H2VEC)
               WRITE(3,'(A2,4X,3F20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
               WRITE(3,'(A2,4X,3F20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
               WRITE(3,'(A2,4X,3F20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
            ENDDO
         ENDDO
         WRITE(3,'(I6)') 3*LNATOMS
         WRITE(3,'(F20.10)') EFINAL
         DO J1=1,LNATOMS
            CALL CONVERT(FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3),
     1                FIN(3*(LNATOMS+J1-1)+1),FIN(3*(LNATOMS+J1-1)+2),FIN(3*(LNATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
            WRITE(3,'(A2,4X,3F20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
            WRITE(3,'(A2,4X,3F20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
            WRITE(3,'(A2,4X,3F20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
         ENDDO
      ELSE
         WRITE(3,'(I6)') LNATOMS
         WRITE(3,'(F20.10)') EINITIAL
         WRITE(3,'(A2,4X,3F20.10)') (ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),J1=1,LNATOMS)
         DO J1=1,NIMAGE
            WRITE(3,'(I6)') LNATOMS
            WRITE(3,'(F20.10)') EREAL(J1) 
            WRITE(3,'(A2,4X,3F20.10)') 
     1    (ZSYM(J2),POINTS(LNOPT*(J1-1)+3*(J2-1)+1),POINTS(LNOPT*(J1-1)+3*(J2-1)+2),POINTS(LNOPT*(J1-1)+3*(J2-1)+3),J2=1,LNATOMS)
         ENDDO
         WRITE(3,'(I6)') LNATOMS
         WRITE(3,'(F20.10)') EFINAL
         WRITE(3,'(A2,4X,3F20.10)') (ZSYM(J1),FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3),J1=1,LNATOMS)
      ENDIF
      CLOSE(3)
C
C  Find the highest energy image.
C
      JMAX=1
      EMAX=-1.0D100
      DO J1=1,NIMAGE
         IF (EREAL(J1).GT.EMAX) THEN
C
C  Don't choose a bad geometry for DFTB.
C
            IF ((.NOT.DFTBT).OR.(DFTBT.AND.EREAL(J1).LT.1.0D6)) THEN
               EMAX=EREAL(J1)
               JMAX=J1
            ENDIF
         ENDIF
      ENDDO
C     PRINT*,'JMAX=',JMAX
      CALL MYCPU_TIME(TIME,.FALSE.)
      WRITE(*,'(A,I4,A,G15.5,A,I6,A,F12.4,A,I4,A,F11.2)') ' Image ',JMAX,' has highest energy=',EMAX,
     1                   ' steps=',ITDONE,' RMS=',RMS,' images=',NIMAGE,'    time=',TIME-TIME0
      TIME0=TIME
      IF ((.NOT.MFLAG).AND.(NEBMAG.GE.NMAGNIFY)) GOTO 10
C
C  Return highest energy point in Q as a candidate transition state
C  We should also return the true gradient for this image
C
      DO J1=1,LNOPT+IADD
         Q(J1)=POINTS((LNOPT+IADD)*(JMAX-1)+J1)
      ENDDO
C
C  The following parameters can change if NEBMAG >0, but are assumed fixed in the rest of the 
C  program.
C
      EINITIAL=EIORIG
      EFINAL=EFORIG
      DO J1=1,LNOPT+IADD
         FIN(J1)=FINORIG(J1)
      ENDDO
      RETURN
      END
C
C****************************************************************************************
C
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C        Line search removed plus small modifications, DJW 2001
C
      SUBROUTINE NEBBFGS(N,M,X,DIAGCO,DIAG,EPS,W,MFLAG,EREAL,RMSREAL,ITMAX,RESET,ITDONE,PTEST,QSTART,LGMAX)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODNEB
      USE modcharmm
      use porfuncs
      IMPLICIT NONE
      INTEGER N,M,J1,ITMAX,ITDONE,NFAIL, ITER
      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH(NIMAGEMAX),DDOT,OVERLAP,LGMAX,LLGMAX(3*NATOMS),
     1                 QSTART(3*NATOMS),DVEC(NIMAGEMAX+1),EINITIAL,EFINAL,SLENGTHTOT,EREALOLD(NIMAGEMAX),
     2                 EPS,ETOTAL,ENEW,GNEW(3*NATOMS*NIMAGEMAX),RMS,EREAL(NIMAGEMAX),RMSREAL(NIMAGEMAX),SEPARATION
      LOGICAL DIAGCO, RESET, PTEST
      DOUBLE PRECISION GNORM,STP(NIMAGEMAX),YS,YY,SQ,YR,BETA,STPMIN
      INTEGER POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,J2
      LOGICAL MFLAG
      DOUBLE PRECISION GSQSCALE, GSTHRESH, DOT1, DOT2
      INTEGER NSPECIAL, NALLOW, NINFO, FRAME
      LOGICAL PVFLAG, CHANGENEBK, CHANGENEBKP, DISCONT, BADIMAGE(NIMAGEMAX), NOBAD
      COMMON /PVF/ PVFLAG
      COMMON /G2/ GSTHRESH, GSQSCALE, NSPECIAL, NALLOW, NINFO
      COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
      COMMON /DVNEB/ DVEC
      INTEGER NDECREASE
      INTEGER IADD, LNOPT, LNATOMS, ISTAT
      COMMON /INTS/ IADD, LNOPT, LNATOMS
      SAVE

      NFAIL=0
      FRAME=1
      IF (RESET) ITER=0
      ITDONE=0
      FIXIMAGE=.FALSE.
      FAILT=.FALSE.
      IF (RESET.AND.PTEST) WRITE(*,'(A)') ' Resetting LBFGS minimiser in NEB module'
      IF ((.NOT.RESET).AND.PTEST) WRITE(*,'(A)') ' Not resetting LBFGS minimiser in NEB module'
C     IF (PV.AND.(.NOT.BFGSTST)) THEN
C        CALL POTENTIAL(X,EREAL(J1),G,.FALSE.,.FALSE.,RMSREAL(J1),.FALSE.,.FALSE.)
C        PVFLAG=.FALSE.
C        CALL PVOPT(X,ENERGY,G)
C     ENDIF

C     IF (ZSYM(1)(1:1).EQ.'W') THEN
C        CALL MAKEGRADH(X,G,ETOTAL,EREAL,RMSREAL,QSTART,DISCONT)
C     ELSE
         CALL MAKEGRADR(X,G,ETOTAL,EREAL,RMSREAL,QSTART,DISCONT,LLGMAX)
C     ENDIF
      DO J1=1,NIMAGE
         EREALOLD(J1)=EREAL(J1)
      ENDDO
      RMS=0.0D0
      DO J1=1,NIMAGE*LNOPT
         RMS=RMS+G(J1)**2
      ENDDO
      RMS=SQRT(RMS/(NIMAGE*LNOPT))

      IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') ' Energy and RMS force=',ETOTAL/NIMAGE,RMS,' after ',0,' LBFGS steps'
C     WRITE(4,'(I6,2F20.10)') 0,EINITIAL,RMSINITIAL
C     WRITE(4,'(I6,2F20.10)') (J1,EREAL(J1),RMSREAL(J1),J1=1,NIMAGE)
C     WRITE(4,'(I6,2F20.10)') NIMAGE+1,EFINAL,RMSFINAL
C     WRITE(4,'(A)') ''

10    CALL FLUSH(6)
      MFLAG=.FALSE.
      IF ((RMS.LE.EPS).AND.(ITDONE.GT.0)) THEN
         MFLAG=.TRUE.
         IF (ITDONE.LT.NSTEPMIN) MFLAG=.FALSE.
         IF (PV.AND.(.NOT.PVFLAG)) MFLAG=.FALSE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
C           WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
            NEBDGUESS=DIAG(1)
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.GE.ITMAX) THEN
         FIXIMAGE=.FALSE.
C        WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         NEBDGUESS=DIAG(1)
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
C           INQUIRE(FILE='diag',EXIST=YESNO)
C           IF (YESNO) THEN
C              OPEN(UNIT=34,FILE='diag',STATUS='OLD')
C              READ(34,*) (DIAG(I),I=1,N)
C              PRINT*,'diag read in LBFGS'
C              WRITE(*,'(6F15.5)') (DIAG(I),I=1,N)
C           ELSE
            DO I=1,N
               DIAG(I)=NEBDGUESS
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
         DO J1=1,NIMAGE
            STP(J1)=MIN(1.0D0/GNORM,GNORM)
         ENDDO
C        STP=1.0D0
      ELSE 
         BOUND=ITER
         IF (ITDONE.GT.M) BOUND=M
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
            DO I=1,N
               DIAG(I)= YS/YY
C              DIAG(I)= MAX(ABS(YS/YY),0.0001D0)
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
         DO J1=1,NIMAGE
            STP(J1)=1.0D0
         ENDDO
      ENDIF
C
C  Store the new search direction
C
      IF (ITER.GT.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
      ENDIF

C     OVERLAP=DDOT(N,G,1,W,1)/SQRT(DDOT(N,G,1,G,1)*DDOT(N,W,1,W,1))
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

C     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
C     PRINT*,'G . G=',DDOT(N,G,1,G,1)
C     PRINT*,'W . W=',DDOT(N,W,1,W,1)
      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reversing step'
         DO I=1,N
            W(ISPT+POINT*N+I)= -W(I)
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF

      DO I=1,N
         W(I)=G(I)
      ENDDO
C
C  We should apply the maximum NEB LBFGS step to each image separately.
C
C  However, using different scale factors for different images leads to huge
C  discontinuities! Now take the minimum scale factor for all images. DJW 18/12/07
C
      STPMIN=1.0D0
      DO J2=1,NIMAGE
         SLENGTH(J2)=0.0D0
         DO J1=1,LNOPT+IADD
            SLENGTH(J2)=SLENGTH(J2)+W(ISPT+POINT*N+(LNOPT+IADD)*(J2-1)+J1)**2
         ENDDO
         SLENGTH(J2)=SQRT(SLENGTH(J2))
         IF (SLENGTH(J2).GT.MAXNEBBFGS) THEN
            STP(J2)=MAXNEBBFGS/SLENGTH(J2)
            STPMIN=MIN(STPMIN,STP(J2))
         ENDIF
      ENDDO
      STP(1:NIMAGE)=STPMIN
C
C  We now have the proposed step.
C
      DO J2=1,NIMAGE
         DO J1=1,LNOPT+IADD
            X(J1+(LNOPT+IADD)*(J2-1))=X(J1+(LNOPT+IADD)*(J2-1))+STP(J2)*W(ISPT+POINT*N+(LNOPT+IADD)*(J2-1)+J1)
         ENDDO 
      ENDDO
      CHANGENEBKP=.FALSE.
      CHANGENEBK=.FALSE.
C
C  For charmm internals must transform and back-transform!
C
      NDECREASE=0
20    CONTINUE
C     IF (ZSYM(1)(1:1).EQ.'W') THEN
C        CALL MAKEGRADH(X,GNEW,ENEW,EREAL,RMSREAL,QSTART,DISCONT)
C     ELSE
         CALL MAKEGRADR(X,GNEW,ENEW,EREAL,RMSREAL,QSTART,DISCONT,LLGMAX)
C     ENDIF
!      IF (NDECREASE.EQ.0) THEN
!         CHANGENEBKP=CHANGENEBK
!         CHANGENEBK=.FALSE.
!         IF ((DVEC(1).LT.SEPARATION/10.0D0).AND.(MOD(ITDONE,1).EQ.0)) CHANGENEBK=.TRUE.
!C        PRINT*,'J1,DVEC,SEPARATION=',1,DVEC(1),SEPARATION
!         DO J1=2,NIMAGE+1
!            IF (DVEC(J1)-DVEC(J1-1).LT.SEPARATION/10.0D0.AND.(MOD(ITDONE,1).EQ.0)) CHANGENEBK=.TRUE.
!C        PRINT*,'J1,DVEC,SEPARATION=',J1,DVEC(J1)-DVEC(J1-1),SEPARATION
!         ENDDO
!         IF (CHANGENEBK) NEBK=NEBK*1.05D0
!         IF (PTEST.AND.CHANGENEBKP) WRITE(*,'(A,F15.5)') 'Increasing force constant to ',NEBK
!      ENDIF

      NOBAD=.TRUE.
      DO J1=1,NIMAGE
         BADIMAGE(J1)=.TRUE.
         IF (((EREAL(J1)-EREALOLD(J1))/ABS(EREAL(J1)).LE.1.0D-2).OR.
     1        (PVFLAG.AND.(EREAL(J1)-EREALOLD(J1).LE.1.0D-8)).OR.PVTS.OR.
     2        (CHANGENEBKP.AND.((EREAL(J1)-EREALOLD(J1))/ABS(EREAL(J1)).LE.0.10)).OR.DISCONT) BADIMAGE(J1)=.FALSE.

C    2     (CHANGENEBKP.AND.((ENEW-ETOTAL)/ABS(ETOTAL).LE.0.10)).OR.DISCONT.OR.BULKT) THEN
C     IF (.TRUE.) THEN

         IF (BADIMAGE(J1)) NOBAD=.FALSE.
      ENDDO

      IF (NOBAD) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ETOTAL=ENEW
         DO J1=1,NIMAGE
            EREALOLD(J1)=EREAL(J1)
         ENDDO
         DO J1=1,N
            G(J1)=GNEW(J1)
         ENDDO
         SLENGTHTOT=0.0D0
         DO J1=1,NIMAGE
            SLENGTHTOT=SLENGTHTOT+STP(J1)*SLENGTH(J1)
         ENDDO
         RMS=0.0D0
         DO J1=1,NIMAGE*LNOPT
            RMS=RMS+G(J1)**2
         ENDDO
         RMS=SQRT(RMS/(NIMAGE*LNOPT))
         IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,G15.10)') ' Energy and RMS force=',ETOTAL/NIMAGE,RMS,' after ',ITDONE,
     1           ' LBFGS steps, step:',SLENGTHTOT

C        WRITE(4,'(I6,2F20.10)') 0,EINITIAL,RMSINITIAL
C        WRITE(4,'(I6,2F20.10)') (J1,EREAL(J1),RMSREAL(J1),J1=1,NIMAGE)
C        WRITE(4,'(I6,2F20.10)') NIMAGE+1,EFINAL,RMSFINAL
C        WRITE(4,'(A)') ''
      ELSE 
C
C  Energy increased for one or more of the images - try again with a smaller step size
C
         IF (NDECREASE.GT.5) THEN
            PRINT*,' in nebbfgs LBFGS step cannot find a lower energy try again'
            NFAIL=NFAIL+1
            DO J1=1,N
               G(J1)=GNEW(J1)
            ENDDO
            IF (NFAIL.GT.20) THEN
               PRINT*,' Too many failures - give up'
C
C DAE make this fail the subsequent ts search automatically (only detected for CHARMM)
C we are in a region of configuration space with discontinuities in the potential.
C
               FAILT=.TRUE.
               FIXIMAGE=.FALSE.
               NEBDGUESS=DIAG(1)
               RETURN
            ENDIF
            DO J2=1,NIMAGE
               IF (BADIMAGE(J2)) THEN
                  DO J1=1,LNOPT+IADD
                     X(J1+(LNOPT+IADD)*(J2-1))=X(J1+(LNOPT+IADD)*(J2-1))-STP(J2)*W(ISPT+POINT*N+(LNOPT+IADD)*(J2-1)+J1)
                  ENDDO 
               ENDIF
            ENDDO
            GOTO 30
         ENDIF
C
C  Should change this to decrease only the step for images whewre the energy has risen.
C
         SLENGTHTOT=0.0D0
         DO J2=1,NIMAGE
            IF (BADIMAGE(J2)) THEN
               DO J1=1,LNOPT+IADD
                  X(J1+(LNOPT+IADD)*(J2-1))=X(J1+(LNOPT+IADD)*(J2-1))-0.9*STP(J2)*W(ISPT+POINT*N+(LNOPT+IADD)*(J2-1)+J1)
               ENDDO 
               STP(J2)=STP(J2)/10.0D0
               SLENGTHTOT=SLENGTHTOT+STP(J2)*SLENGTH(J2)
               IF (PTEST) 
     1             WRITE(*,'(A,I5,A,G16.10,A,G19.10,A,F15.8)') ' function for image ',J2,' increased from ',
     2                                                     EREALOLD(J2),' to ',EREAL(J2),
     2                                                    ' decreasing step to ',STP(J2)*SLENGTH(J2)
            ENDIF
         ENDDO
         NDECREASE=NDECREASE+1
C
C  There is a big problem with FIXIMAGE - the ANV array could be different for all the
C  images! Also, if we want to recover the previous geometry in the limit of zero step
C  then ANV must return to the same array.
C
C        FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
C
C     Compute the new step and gradient change
C
30    NPT=POINT*N
      DO J2=1,NIMAGE
         DO I=1,LNOPT+IADD
            W(ISPT+NPT+(LNOPT+IADD)*(J2-1)+I)= STP(J2)*W(ISPT+NPT+(LNOPT+IADD)*(J2-1)+I)
         ENDDO
      ENDDO
      DO I=1,N
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
C     PRINT*,'diag in MYLBFGS'
C     WRITE(*,'(6F15.5)') (DIAG(I),I=1,12)
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITER.GE.FIXAFTER)) FIXIMAGE=.TRUE.
      GOTO 10

      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE MAKEGRADH(X,G,SUM,EREAL,RMSREAL,QSTART,DISCONT)
      USE COMMONS
      USE MODTWOEND
      USE MODNEB
      USE MODHESS
      USE KEY,ONLY : NEBK
      use porfuncs
      IMPLICIT NONE
      INTEGER J1,J2
      DOUBLE PRECISION DUMMY1,EREAL(NIMAGEMAX),ETOTAL,DVEC(NIMAGEMAX+1),
     1                 XX(3*NATOMS),GG(3*NATOMS),RMS,G(3*NATOMS*NIMAGEMAX),X(3*NATOMS*NIMAGEMAX),
     2                 QSTART(3*NATOMS),RMSREAL(NIMAGEMAX),SUM,EINITIAL,EFINAL,SEPARATION,PI2
      PARAMETER (PI2=6.283185307D0)
      LOGICAL RADMOVED, OVERLP, DISCONT
      COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
      COMMON /DVNEB/ DVEC
      INTEGER IADD, LNOPT, LNATOMS, ISTAT
      COMMON /INTS/ IADD, LNOPT, LNATOMS
      COMMON /DISCON/ RADMOVED, OVERLP

      ETOTAL=0.0D0
      DISCONT=.FALSE.
      DO J1=1,NIMAGE
         DO J2=1,LNOPT
            XX(J2)=X(J2+(LNOPT+IADD)*(J1-1))
         ENDDO
         DO J2=LNOPT+1,IADD
C           IF (X(J2+(LNOPT+IADD)*(J1-1)).LT.0.0D0) X(J2+(LNOPT+IADD)*(J1-1))=X(J2+(LNOPT+IADD)*(J1-1))+PI2
            XX(J2)=X(J2+(LNOPT+IADD)*(J1-1))
         ENDDO
         CALL POTENTIAL(XX,EREAL(J1),GG,.TRUE.,.FALSE.,RMSREAL(J1),.FALSE.,.FALSE.)
         IF (RADMOVED) DISCONT=.TRUE.
C
C  The points can be changed for C60 in the potential.
C
         IF (OVERLP) THEN
            DISCONT=.TRUE.
            DO J2=1,LNOPT
               X(J2+(LNOPT+IADD)*(J1-1))=XX(J2)
            ENDDO
         ENDIF
         DO J2=1,LNOPT+IADD
            G(J2+(LNOPT+IADD)*(J1-1))=GG(J2)
         ENDDO
C        WRITE(*,'(I6,2G20.10)') J1,EREAL(J1),RMSREAL(J1)
         CALL FLUSH(6)
         ETOTAL=ETOTAL+EREAL(J1)
      ENDDO

      SUM=ETOTAL
      DUMMY1=0.0D0
      DO J1=1,LNOPT+IADD
         DUMMY1=DUMMY1+(X(J1)-QSTART(J1))**2              !  (0,1)
      ENDDO
      SUM=SUM+NEBK*DUMMY1/2
      DVEC(1)=SQRT(DUMMY1)

      DO J2=1,NIMAGE-1
         DUMMY1=0.0D0
         DO J1=1,LNOPT+IADD
            DUMMY1=DUMMY1+(X(J1+(LNOPT+IADD)*(J2-1))-X(J1+(LNOPT+IADD)*J2))**2 ! (J2,J2+1)
         ENDDO
         DVEC(J2+1)=SQRT(DUMMY1)+DVEC(J2)
         SUM=SUM+NEBK*DUMMY1/2
      ENDDO

      DUMMY1=0.0D0
      DO J1=1,LNOPT+IADD
         DUMMY1=DUMMY1+(X(J1+(LNOPT+IADD)*(NIMAGE-1))-FIN(J1))**2 !  (N,N+1)
      ENDDO
      SUM=SUM+NEBK*DUMMY1/2
      DVEC(NIMAGE+1)=SQRT(DUMMY1)+DVEC(NIMAGE)

      IF (NIMAGE.GT.1) THEN
         RMS=0.0D0
         DO J1=1,LNOPT+IADD
            G(J1)=G(J1)+NEBK*(2*X(J1)-QSTART(J1)-X(J1+LNOPT+IADD))
            RMS=RMS+G(J1)**2
         ENDDO
   
         DO J1=1,LNOPT+IADD
            G(J1+(LNOPT+IADD)*(NIMAGE-1))=G(J1+(LNOPT+IADD)*(NIMAGE-1))+NEBK*(2*X(J1+(LNOPT+IADD)*(NIMAGE-1))
     1                   -X(J1+(LNOPT+IADD)*(NIMAGE-2))-FIN(J1))
            RMS=RMS+G(J1)**2
         ENDDO
         DO J2=2,NIMAGE-1
            DO J1=1,LNOPT+IADD
               G(J1+(LNOPT+IADD)*(J2-1))=G(J1+(LNOPT+IADD)*(J2-1))+NEBK*(2*X(J1+(LNOPT+IADD)*(J2-1))
     1               -X(J1+(LNOPT+IADD)*(J2-2))-X(J1+(LNOPT+IADD)*J2))
               RMS=RMS+G(J1+(LNOPT+IADD)*(J2-1))**2
            ENDDO
         ENDDO
      ELSE
         DO J1=1,LNOPT+IADD
            G(J1)=G(J1)+NEBK*(2*X(J1)-QSTART(J1)-FIN(J1))
            RMS=RMS+G(J1)**2
         ENDDO
         RMS=0.0D0
         DO J1=1,LNOPT+IADD
            RMS=RMS+G(J1)**2
         ENDDO
      ENDIF
      RMS=SQRT(RMS/((LNOPT+IADD)*NIMAGE))

      RETURN
      END
C
C****************************************************************************************
C
      SUBROUTINE MAKEGRADR(X,G,ETOTAL,EREAL,RMSREAL,QSTART,DISCONT,LLGMAX)
      USE COMMONS
      USE MODTWOEND
      USE MODNEB
      USE MODHESS
      USE KEY,ONLY : NEBK, FROZEN, FREEZE, VARIABLES
      use porfuncs
      IMPLICIT NONE
      INTEGER J1,J2
      DOUBLE PRECISION DUMMY1,DUMMY2,EREAL(NIMAGEMAX),ETOTAL,TANVEC(3*NATOMS,NIMAGEMAX),WPLUS,WMINUS,EMAX,
     1                 XX(3*NATOMS),GG(3*NATOMS),RMS,G(3*NATOMS*NIMAGEMAX),X(3*NATOMS*NIMAGEMAX),LLGMAX(3*NATOMS),
     2                 QSTART(3*NATOMS),RMSREAL(NIMAGEMAX),DVEC(NIMAGEMAX+1),EINITIAL,EFINAL,SEPARATION,PI2
      PARAMETER (PI2=6.283185307D0)
      COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
      COMMON /DVNEB/ DVEC
      INTEGER IADD, LNOPT, LNATOMS, ISTAT
      COMMON /INTS/ IADD, LNOPT, LNATOMS
      LOGICAL RADMOVED, DISCONT, OVERLP
      COMMON /DISCON/ RADMOVED, OVERLP

      ETOTAL=0.0D0
      DISCONT=.FALSE.
      EMAX=-1.0D100
      DO J1=1,NIMAGE
         DO J2=1,LNOPT+IADD
            XX(J2)=X(J2+(LNOPT+IADD)*(J1-1))
         ENDDO

         CALL POTENTIAL(XX,EREAL(J1),GG,.TRUE.,.FALSE.,RMSREAL(J1),.FALSE.,.FALSE.)
         IF (EREAL(J1).GT.EMAX) LLGMAX(1:3*NATOMS)=G(1:3*NATOMS) ! save for use in bfgsts
         IF (RADMOVED) DISCONT=.TRUE.
C
C  The points can be changed for C60 in the potential subroutine.
C
         IF (OVERLP) THEN
            DISCONT=.TRUE.
            DO J2=1,LNOPT
               X(J2+(LNOPT+IADD)*(J1-1))=XX(J2)
            ENDDO
         ENDIF
C        DO J2=1,LNOPT+IADD
         DO J2=1,LNOPT
            G(J2+(LNOPT+IADD)*(J1-1))=GG(J2)
         ENDDO
         DO J2=LNOPT+1,LNOPT+IADD
            G(J2+(LNOPT+IADD)*(J1-1))=0.0D0
         ENDDO
C        WRITE(*,'(I6,2G20.10)') J1,EREAL(J1),RMSREAL(J1)
         CALL FLUSH(6)
         ETOTAL=ETOTAL+EREAL(J1)
      ENDDO
      ETOTAL=ETOTAL
C
C  Define tangent vectors. This is the scheme from Henkelman and Jonsson JCP, 113, 9978, 2000
C
      IF ((EREAL(1).GT.EINITIAL).AND.(EREAL(1).LT.EREAL(2))) THEN
         WPLUS=1.0D0
         WMINUS=0.0D0
      ELSE IF ((EREAL(1).LT.EINITIAL).AND.(EREAL(1).GT.EREAL(2))) THEN
         WPLUS=0.0D0
         WMINUS=1.0D0
      ELSE 
         IF (EREAL(2).GT.EINITIAL) THEN
            WPLUS= MAX(ABS(EREAL(1)-EREAL(2)),ABS(EREAL(1)-EINITIAL))
            WMINUS=MIN(ABS(EREAL(1)-EREAL(2)),ABS(EREAL(1)-EINITIAL))
         ELSE
            WPLUS= MIN(ABS(EREAL(1)-EREAL(2)),ABS(EREAL(1)-EINITIAL))
            WMINUS=MAX(ABS(EREAL(1)-EREAL(2)),ABS(EREAL(1)-EINITIAL))
         ENDIF
      ENDIF
      DUMMY1=0.0D0

      DO J2=1,LNOPT
         TANVEC(J2,1)=WMINUS*(X(J2)-QSTART(J2)) + WPLUS*(X(J2+(LNOPT+IADD))-X(J2))
         DUMMY1=DUMMY1+(X(J2)-QSTART(J2))**2
      ENDDO
      DVEC(1)=SQRT(DUMMY1)

      IF (NIMAGE.GT.1) THEN
         IF ((EREAL(NIMAGE).GT.EFINAL).AND.(EREAL(NIMAGE).LT.EREAL(NIMAGE-1))) THEN
            WPLUS=0.0D0
            WMINUS=1.0D0
         ELSE IF ((EREAL(NIMAGE).LT.EFINAL).AND.(EREAL(NIMAGE).GT.EREAL(NIMAGE-1))) THEN
            WPLUS=1.0D0
            WMINUS=0.0D0
         ELSE 
            IF (EREAL(NIMAGE-1).LT.EFINAL) THEN
               WPLUS= MAX(ABS(EREAL(NIMAGE)-EREAL(NIMAGE-1)),ABS(EREAL(NIMAGE)-EFINAL))
               WMINUS=MIN(ABS(EREAL(NIMAGE)-EREAL(NIMAGE-1)),ABS(EREAL(NIMAGE)-EFINAL))
            ELSE
               WPLUS= MIN(ABS(EREAL(NIMAGE)-EREAL(NIMAGE-1)),ABS(EREAL(NIMAGE)-EFINAL))
               WMINUS=MAX(ABS(EREAL(NIMAGE)-EREAL(NIMAGE-1)),ABS(EREAL(NIMAGE)-EFINAL))
            ENDIF
         ENDIF
         DUMMY1=0.0D0
         DO J2=1,LNOPT
            TANVEC(J2,NIMAGE)=WMINUS*(X(J2+(LNOPT+IADD)*(NIMAGE-1))-X(J2+(LNOPT+IADD)*(NIMAGE-2))) 
     1                       + WPLUS*(FIN(J2)-X(J2+(LNOPT+IADD)*(NIMAGE-1)))
            DUMMY1=DUMMY1+(FIN(J2)-X(J2+(LNOPT+IADD)*(NIMAGE-1)))**2
         ENDDO
      ELSE
         DUMMY1=0.0D0
         DO J2=1,LNOPT
            TANVEC(J2,NIMAGE)=WMINUS*X(J2+(LNOPT+IADD)*(NIMAGE-1))
     1                       + WPLUS*(FIN(J2)-X(J2+(LNOPT+IADD)*(NIMAGE-1)))
            DUMMY1=DUMMY1+(FIN(J2)-X(J2+(LNOPT+IADD)*(NIMAGE-1)))**2
         ENDDO
      ENDIF
      DVEC(NIMAGE+1)=SQRT(DUMMY1)

      DO J1=2,NIMAGE-1
         IF ((EREAL(J1).GT.EREAL(J1-1)).AND.(EREAL(J1).LT.EREAL(J1+1))) THEN
            WPLUS=1.0D0
            WMINUS=0.0D0
         ELSE IF ((EREAL(J1).LT.EREAL(J1-1)).AND.(EREAL(J1).GT.EREAL(J1+1))) THEN
            WPLUS=0.0D0
            WMINUS=1.0D0
         ELSE 
            IF (EREAL(J1+1).GT.EREAL(J1-1)) THEN
               WPLUS= MAX(ABS(EREAL(J1)-EREAL(J1+1)),ABS(EREAL(J1)-EREAL(J1-1)))
               WMINUS=MIN(ABS(EREAL(J1)-EREAL(J1+1)),ABS(EREAL(J1)-EREAL(J1-1)))
            ELSE
               WPLUS= MIN(ABS(EREAL(J1)-EREAL(J1+1)),ABS(EREAL(J1)-EREAL(J1-1)))
               WMINUS=MAX(ABS(EREAL(J1)-EREAL(J1+1)),ABS(EREAL(J1)-EREAL(J1-1)))
            ENDIF
         ENDIF
         DUMMY1=0.0D0
         DO J2=1,LNOPT
            TANVEC(J2,J1)=WMINUS*(X(J2+(LNOPT+IADD)*(J1-1))-X(J2+(LNOPT+IADD)*(J1-2))) 
     1               + WPLUS*(X(J2+(LNOPT+IADD)*J1)-X(J2+(LNOPT+IADD)*(J1-1)))
            DUMMY1=DUMMY1+(X(J2+(LNOPT+IADD)*(J1-1))-X(J2+(LNOPT+IADD)*(J1-2)))**2
         ENDDO
         DVEC(J1)=SQRT(DUMMY1)+DVEC(J1-1)
      ENDDO
      DUMMY1=0.0D0
      IF (NIMAGE.GT.1) THEN
         DO J2=1,LNOPT
            DUMMY1=DUMMY1+(X(J2+(LNOPT+IADD)*(NIMAGE-2))-X(J2+(LNOPT+IADD)*(NIMAGE-1)))**2
         ENDDO
      ENDIF
      IF (NIMAGE.GT.1) THEN
         DVEC(NIMAGE)=SQRT(DUMMY1)+DVEC(NIMAGE-1)
      ENDIF
      DVEC(NIMAGE+1)= DVEC(NIMAGE+1)+DVEC(NIMAGE)
C
C  Normalise tangent vectors
C
      DO J1=1,NIMAGE
         DUMMY1=0.0D0
         DO J2=1,LNOPT
            DUMMY1=DUMMY1+TANVEC(J2,J1)**2
         ENDDO
         DUMMY1=1.0D0/SQRT(DUMMY1)
         DO J2=1,LNOPT
            TANVEC(J2,J1)=TANVEC(J2,J1)*DUMMY1
         ENDDO
      ENDDO
C
C  Gradient of the potential perpendicular to the tangent vector.
C
      RMS=0.0D0
      DO J1=1,NIMAGE
         DUMMY1=0.0D0
         DO J2=1,LNOPT
            DUMMY1=DUMMY1+G(J2+(LNOPT+IADD)*(J1-1))*TANVEC(J2,J1)
         ENDDO
         DO J2=1,LNOPT
            G(J2+(LNOPT+IADD)*(J1-1))=G(J2+(LNOPT+IADD)*(J1-1))-DUMMY1*TANVEC(J2,J1)
            RMS=RMS+G(J2+(LNOPT+IADD)*(J1-1))**2
         ENDDO
      ENDDO
      RMS=SQRT(RMS/(NIMAGE*(LNOPT+IADD)*1.0D0))
C
C  Spring energy derivatives.
C
      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,LNOPT
         DUMMY1=DUMMY1+(X(J2)-X(J2+(LNOPT+IADD)))**2
         DUMMY2=DUMMY2+(X(J2)-QSTART(J2))**2
      ENDDO
      DUMMY1=SQRT(DUMMY1)*NEBK
      DUMMY2=SQRT(DUMMY2)*NEBK
      DO J2=1,LNOPT
         G(J2)=G(J2)-(DUMMY1-DUMMY2)*TANVEC(J2,1)
      ENDDO

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      IF (NIMAGE.LT.1) THEN
         PRINT '(A)',' oldneb> no images - quit'
         STOP
      ENDIF

      IF (NIMAGE.GT.1) THEN
         DO J2=1,LNOPT
            DUMMY1=DUMMY1+(X(J2+(LNOPT+IADD)*(NIMAGE-1))-FIN(J2))**2
            DUMMY2=DUMMY2+(X(J2+(LNOPT+IADD)*(NIMAGE-1))-X(J2+(LNOPT+IADD)*(NIMAGE-2)))**2
         ENDDO
C        DO J2=LNOPT+1,LNOPT+IADD
C           DUMMY1=DUMMY1+MOD(X(J2+(LNOPT+IADD)*(NIMAGE-1))-FIN(J2),PI2)**2
C           DUMMY2=DUMMY2+MOD(X(J2+(LNOPT+IADD)*(NIMAGE-1))-X(J2+(LNOPT+IADD)*(NIMAGE-2)),PI2)**2
C        ENDDO
      ELSE
         DO J2=1,LNOPT
            DUMMY1=DUMMY1+(X(J2+(LNOPT+IADD)*(NIMAGE-1))-FIN(J2))**2
            DUMMY2=DUMMY2+(X(J2+(LNOPT+IADD)*(NIMAGE-1)))**2
         ENDDO
C        DO J2=LNOPT+1,LNOPT+IADD
C           DUMMY1=DUMMY1+MOD(X(J2+(LNOPT+IADD)*(NIMAGE-1))-FIN(J2),PI2)**2
C           DUMMY2=DUMMY2+MOD(X(J2+(LNOPT+IADD)*(NIMAGE-1)),PI2)**2
C        ENDDO
      ENDIF

C     DO J2=1,LNOPT
C        DUMMY1=DUMMY1+(X(J2+(LNOPT+IADD)*(NIMAGE-1))-FIN(J2))**2
C        DUMMY2=DUMMY2+(X(J2+(LNOPT+IADD)*(NIMAGE-1))-X(J2+(LNOPT+IADD)*(NIMAGE-2)))**2
C     ENDDO

      DUMMY1=SQRT(DUMMY1)*NEBK
      DUMMY2=SQRT(DUMMY2)*NEBK
      DO J2=1,LNOPT
         G(J2+(LNOPT+IADD)*(NIMAGE-1))=G(J2+(LNOPT+IADD)*(NIMAGE-1))-(DUMMY1-DUMMY2)*TANVEC(J2,NIMAGE)
      ENDDO

      DO J1=2,NIMAGE-1
         DUMMY1=0.0D0
         DUMMY2=0.0D0
         DO J2=1,LNOPT
            DUMMY1=DUMMY1+(X(J2+(LNOPT+IADD)*(J1-1))-X(J2+(LNOPT+IADD)*J1))**2
            DUMMY2=DUMMY2+(X(J2+(LNOPT+IADD)*(J1-1))-X(J2+(LNOPT+IADD)*(J1-2)))**2
         ENDDO
         DUMMY1=SQRT(DUMMY1)*NEBK
         DUMMY2=SQRT(DUMMY2)*NEBK
         DO J2=1,LNOPT
            G(J2+(LNOPT+IADD)*(J1-1))=G(J2+(LNOPT+IADD)*(J1-1))-(DUMMY1-DUMMY2)*TANVEC(J2,J1)
         ENDDO
      ENDDO
!   
! Set gradients on frozen atoms to zero.
!   
      IF (FREEZE) THEN
         DO J1=1,NIMAGE
            DO J2=1,NATOMS
               IF (FROZEN(J2)) THEN
                  IF (VARIABLES) THEN
                     G(LNOPT*(J1-1)+J2)=0.0D0
                  ELSE
                     G(LNOPT*(J1-1)+3*(J2-1)+1)=0.0D0
                     G(LNOPT*(J1-1)+3*(J2-1)+2)=0.0D0
                     G(LNOPT*(J1-1)+3*(J2-1)+3)=0.0D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
C
C  Put a specified number of images between specified starting and
C  finishing points using suitable interpolation.
C  OFFSET tells us how many images we already have saved in POINTS, so
C  the next image goes at position OFFSET+1
C  
      SUBROUTINE MAKEIMAGE(POINTS,PERM,ITEST,OFFSET,Q)
      USE COMMONS
      USE MODTWOEND
      USE KEY
      USE MODNEB
      IMPLICIT NONE
      INTEGER OFFSET
      LOGICAL PERM, ITEST
      INTEGER IADD, LNOPT, LNATOMS
      COMMON /INTS/ IADD, LNOPT, LNATOMS
      INTEGER J1,J2,J3
      DOUBLE PRECISION VECS(3*NATOMS), STEP, GUESS(3*NATOMS),DPRAND,TEMP(3),Q(3*NATOMS),
     1                 POINTS(5*NATOMS*NIMAGEMAX), OVEC1(3),H1VEC1(3),H2VEC1(3), C2VEC(3),
     2                 OVEC2(3),H1VEC2(3),H2VEC2(3),
     3                 HSVEL(3*NATOMS),
     4                 DUMMY, EFINAL, EINITIAL,SEPARATION,RMS,DIST
      COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: FINCART,STARTCART,POINTSCART

      IF (ZSYM(1)(1:1).EQ.'W') THEN
         ALLOCATE(FINCART(9*LNATOMS),STARTCART(9*LNATOMS),POINTSCART(9*LNATOMS))
      ELSE
         ALLOCATE(FINCART(3*NATOMS),STARTCART(3*NATOMS),POINTSCART(3*NATOMS))
      ENDIF
      IF (.NOT.ALLOCATED(LANV)) ALLOCATE(LANV(3*LNATOMS+IADD))
C
C  Bulk principal image code assumes that PARAM1, PARAM2 and PARAM3 are the box lengths.
C
      DO J1=1,LNATOMS+IADD/3
         J3=3*(J1-1)
         VECS(J3+1)=FIN(J3+1)-Q(J3+1)
         VECS(J3+2)=FIN(J3+2)-Q(J3+2)
         VECS(J3+3)=FIN(J3+3)-Q(J3+3)
         LANV(J3+1)=0
         LANV(J3+2)=0
         LANV(J3+3)=0
         IF (BULKT.AND.(.NOT.NORESET)) THEN
            LANV(J3+1)=NINT(VECS(J3+1)/PARAM1)
            LANV(J3+2)=NINT(VECS(J3+2)/PARAM2)
            LANV(J3+3)=NINT(VECS(J3+3)/PARAM3)
         ENDIF
         VECS(J3+1)=VECS(J3+1)-PARAM1*LANV(J3+1)
         VECS(J3+2)=VECS(J3+2)-PARAM2*LANV(J3+2)
         VECS(J3+3)=VECS(J3+3)-PARAM3*LANV(J3+3)
         START(J3+1)=Q(J3+1)
         START(J3+2)=Q(J3+2)
         START(J3+3)=Q(J3+3)
      ENDDO
C
C  Initial distance
C
      DIST=0.0D0
      IF (ZSYM(1)(1:1).EQ.'W') THEN
         DO J1=1,LNATOMS
            CALL CONVERT(FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3),
     1    FIN(3*(LNATOMS+J1-1)+1),FIN(3*(LNATOMS+J1-1)+2),FIN(3*(LNATOMS+J1-1)+3),OVEC1,H1VEC1,H2VEC1)
            FINCART(9*(J1-1)+1)=OVEC1(1)
            FINCART(9*(J1-1)+2)=OVEC1(2)
            FINCART(9*(J1-1)+3)=OVEC1(3)
            FINCART(9*(J1-1)+4)=H1VEC1(1)
            FINCART(9*(J1-1)+5)=H1VEC1(2)
            FINCART(9*(J1-1)+6)=H1VEC1(3)
            FINCART(9*(J1-1)+7)=H2VEC1(1)
            FINCART(9*(J1-1)+8)=H2VEC1(2)
            FINCART(9*(J1-1)+9)=H2VEC1(3)
            CALL CONVERT(START(3*(J1-1)+1),START(3*(J1-1)+2),START(3*(J1-1)+3),
     1    START(3*(LNATOMS+J1-1)+1),START(3*(LNATOMS+J1-1)+2),START(3*(LNATOMS+J1-1)+3),OVEC2,H1VEC2,H2VEC2)
            STARTCART(9*(J1-1)+1)=OVEC2(1)
            STARTCART(9*(J1-1)+2)=OVEC2(2)
            STARTCART(9*(J1-1)+3)=OVEC2(3)
            STARTCART(9*(J1-1)+4)=H1VEC2(1)
            STARTCART(9*(J1-1)+5)=H1VEC2(2)
            STARTCART(9*(J1-1)+6)=H1VEC2(3)
            STARTCART(9*(J1-1)+7)=H2VEC2(1)
            STARTCART(9*(J1-1)+8)=H2VEC2(2)
            STARTCART(9*(J1-1)+9)=H2VEC2(3)
            DIST=DIST+( OVEC1(1)- OVEC2(1))**2+( OVEC1(2)- OVEC2(2))**2+( OVEC1(3)- OVEC2(3))**2+
     1                (H1VEC1(1)-H1VEC2(1))**2+(H1VEC1(2)-H1VEC2(2))**2+(H1VEC1(3)-H1VEC2(3))**2+
     2                (H2VEC1(1)-H2VEC2(1))**2+(H2VEC1(2)-H2VEC2(2))**2+(H2VEC1(3)-H2VEC2(3))**2
         ENDDO 
      ELSE
         DO J1=1,LNOPT
            DIST=DIST+(FIN(J1)-START(J1)-PARAM3*LANV(J1))**2
         ENDDO
      ENDIF
      DIST=SQRT(DIST)
C     PRINT *,'oldneb> DIST=',DIST
C
C  Assign initial image positions.
C
      SEPARATION=DIST/(NIMAGE+1)
C     NEBK=(ABS(EFINAL)+ABS(EINITIAL))/(100*(NIMAGE+1)*SEPARATION**2)
      IF (ZSYM(1)(1:1).EQ.'W') NEBK=NEBK*0.0D0
      IF (ZSYM(1)(1:1).EQ.'W') THEN
C
C  Note that interpolating the Cartesian coordinates will generally break the geometry of the individual water molecules,
C  so there won't be a perfect conversion back to centre-of-mass/Euler angle coordinates.
C
         DO J1=1,NIMAGE
            DO J2=1,LNATOMS
               POINTSCART(9*(J2-1)+1)=STARTCART(9*(J2-1)+1)+(FINCART(9*(J2-1)+1)-STARTCART(9*(J2-1)+1))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+2)=STARTCART(9*(J2-1)+2)+(FINCART(9*(J2-1)+2)-STARTCART(9*(J2-1)+2))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+3)=STARTCART(9*(J2-1)+3)+(FINCART(9*(J2-1)+3)-STARTCART(9*(J2-1)+3))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+4)=STARTCART(9*(J2-1)+4)+(FINCART(9*(J2-1)+4)-STARTCART(9*(J2-1)+4))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+5)=STARTCART(9*(J2-1)+5)+(FINCART(9*(J2-1)+5)-STARTCART(9*(J2-1)+5))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+6)=STARTCART(9*(J2-1)+6)+(FINCART(9*(J2-1)+6)-STARTCART(9*(J2-1)+6))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+7)=STARTCART(9*(J2-1)+7)+(FINCART(9*(J2-1)+7)-STARTCART(9*(J2-1)+7))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+8)=STARTCART(9*(J2-1)+8)+(FINCART(9*(J2-1)+8)-STARTCART(9*(J2-1)+8))*J1*SEPARATION/DIST
               POINTSCART(9*(J2-1)+9)=STARTCART(9*(J2-1)+9)+(FINCART(9*(J2-1)+9)-STARTCART(9*(J2-1)+9))*J1*SEPARATION/DIST
               OVEC1(1)=POINTSCART(9*(J2-1)+1)
               OVEC1(2)=POINTSCART(9*(J2-1)+2)
               OVEC1(3)=POINTSCART(9*(J2-1)+3)
               H1VEC1(1)=POINTSCART(9*(J2-1)+4)
               H1VEC1(2)=POINTSCART(9*(J2-1)+5)
               H1VEC1(3)=POINTSCART(9*(J2-1)+6)
               H2VEC1(1)=POINTSCART(9*(J2-1)+7)
               H2VEC1(2)=POINTSCART(9*(J2-1)+8)
               H2VEC1(3)=POINTSCART(9*(J2-1)+9)
               DUMMY=(H2VEC1(1)-H1VEC1(1))**2+(H2VEC1(2)-H1VEC1(2))**2+(H2VEC1(3)-H1VEC1(3))**2 
               IF ((DUMMY.LT.0.1D0.AND.PERM)) THEN
                  C2VEC(1)=(STARTCART(9*(J2-1)+4)+FINCART(9*(J2-1)+4))/2-STARTCART(9*(J2-1)+1)
                  C2VEC(2)=(STARTCART(9*(J2-1)+5)+FINCART(9*(J2-1)+5))/2-STARTCART(9*(J2-1)+2)
                  C2VEC(3)=(STARTCART(9*(J2-1)+6)+FINCART(9*(J2-1)+6))/2-STARTCART(9*(J2-1)+3)
                  DUMMY=SQRT(C2VEC(1)**2+C2VEC(2)**2+C2VEC(3)**2)
                  C2VEC(1)=C2VEC(1)/DUMMY
                  C2VEC(2)=C2VEC(2)/DUMMY
                  C2VEC(3)=C2VEC(3)/DUMMY
                  DUMMY=C2VEC(1)*(STARTCART(9*(J2-1)+4)-STARTCART(9*(J2-1)+1))
     1                 +C2VEC(2)*(STARTCART(9*(J2-1)+5)-STARTCART(9*(J2-1)+2))
     2                 +C2VEC(3)*(STARTCART(9*(J2-1)+6)-STARTCART(9*(J2-1)+3))
                  TEMP(1)=C2VEC(1)*DUMMY+(STARTCART(9*(J2-1)+5)-STARTCART(9*(J2-1)+2))*C2VEC(3)
     1                                  -(STARTCART(9*(J2-1)+6)-STARTCART(9*(J2-1)+3))*C2VEC(2)
                  TEMP(2)=C2VEC(2)*DUMMY-(STARTCART(9*(J2-1)+4)-STARTCART(9*(J2-1)+1))*C2VEC(3)
     1                                  +(STARTCART(9*(J2-1)+6)-STARTCART(9*(J2-1)+3))*C2VEC(1)
                  TEMP(3)=C2VEC(3)*DUMMY+(STARTCART(9*(J2-1)+4)-STARTCART(9*(J2-1)+1))*C2VEC(2)
     1                                  -(STARTCART(9*(J2-1)+5)-STARTCART(9*(J2-1)+2))*C2VEC(1)
                  H1VEC1(1)=STARTCART(9*(J2-1)+1)+TEMP(1)
                  H1VEC1(2)=STARTCART(9*(J2-1)+2)+TEMP(2)
                  H1VEC1(3)=STARTCART(9*(J2-1)+3)+TEMP(3)
                  DUMMY=C2VEC(1)*(STARTCART(9*(J2-1)+7)-STARTCART(9*(J2-1)+1))
     1                 +C2VEC(2)*(STARTCART(9*(J2-1)+8)-STARTCART(9*(J2-1)+2))
     2                 +C2VEC(3)*(STARTCART(9*(J2-1)+9)-STARTCART(9*(J2-1)+3))
                  TEMP(1)=C2VEC(1)*DUMMY+(STARTCART(9*(J2-1)+8)-STARTCART(9*(J2-1)+2))*C2VEC(3)
     1                                  -(STARTCART(9*(J2-1)+9)-STARTCART(9*(J2-1)+3))*C2VEC(2)
                  TEMP(2)=C2VEC(2)*DUMMY-(STARTCART(9*(J2-1)+7)-STARTCART(9*(J2-1)+1))*C2VEC(3)
     1                                  +(STARTCART(9*(J2-1)+9)-STARTCART(9*(J2-1)+3))*C2VEC(1)
                  TEMP(3)=C2VEC(3)*DUMMY+(STARTCART(9*(J2-1)+7)-STARTCART(9*(J2-1)+1))*C2VEC(2)
     1                                  -(STARTCART(9*(J2-1)+8)-STARTCART(9*(J2-1)+2))*C2VEC(1)
                  H2VEC1(1)=STARTCART(9*(J2-1)+1)+TEMP(1)
                  H2VEC1(2)=STARTCART(9*(J2-1)+2)+TEMP(2)
                  H2VEC1(3)=STARTCART(9*(J2-1)+3)+TEMP(3)
               ENDIF
               CALL CONVERT2(OVEC1,H1VEC1,H2VEC1,
     1            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+1),
     2            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+2),
     3            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+3),
     4            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*LNATOMS+3*(J2-1)+1),
     5            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*LNATOMS+3*(J2-1)+2),
     6            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*LNATOMS+3*(J2-1)+3) )
C              WRITE(*,'(A,I5,6F15.5)') 'J2,X,Y,Z,A,B,C=',J2,
C    1            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+1),
C    2            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+2),
C    3            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+3),
C    4            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*LNATOMS+3*(J2-1)+1),
C    5            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*LNATOMS+3*(J2-1)+2),
C    6            POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*LNATOMS+3*(J2-1)+3) 
C              PRINT*,'points before:'
C              WRITE(*,'(A6,3F20.10)') 'H     ',H1VEC1(1),H1VEC1(2),H1VEC1(3)
C              WRITE(*,'(A6,3F20.10)') 'H     ',H2VEC1(1),H2VEC1(2),H2VEC1(3)
C              WRITE(*,'(A6,3F20.10)') 'O     ',OVEC1(1),OVEC1(2),OVEC1(3)

            ENDDO
C
C  If any of the centre-of-mass distances are too small then move back towards the previous image.
C  This can happen for permutational isomers where whole molecules are exchanged.
C
            IF (J1.GT.1) THEN
               DO J2=1,LNATOMS
                  DO J3=J2+1,LNATOMS
                     DUMMY=(POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+1)
     1                     -POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+1))**2
     1                    +(POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+2)
     1                     -POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+2))**2
     2                    +(POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+3)
     1                     -POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+3))**2
                     IF (DUMMY.LT.0.1D0) THEN
                        POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+1)=(POINTS((J1-2+OFFSET)*(LNOPT+IADD)+3*(J2-1)+1)+
     1                                                      POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+1))/2
                        POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+2)=(POINTS((J1-2+OFFSET)*(LNOPT+IADD)+3*(J2-1)+2)+
     1                                                      POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+2))/2
                        POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+3)=(POINTS((J1-2+OFFSET)*(LNOPT+IADD)+3*(J2-1)+3)+
     1                                                      POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J2-1)+3))/2
                        POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+1)=(POINTS((J1-2+OFFSET)*(LNOPT+IADD)+3*(J3-1)+1)+
     1                                                      POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+1))/2
                        POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+2)=(POINTS((J1-2+OFFSET)*(LNOPT+IADD)+3*(J3-1)+2)+
     1                                                      POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+2))/2
                        POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+3)=(POINTS((J1-2+OFFSET)*(LNOPT+IADD)+3*(J3-1)+3)+
     1                                                      POINTS((J1-1+OFFSET)*(LNOPT+IADD)+3*(J3-1)+3))/2
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ELSE IF (PERM) THEN
C
C  Not done for bulk. A little bit of noise is added because otherwise head-on collisions don't
C  achieve much.
C
         DO J1=1,LNATOMS
            HSVEL(3*(J1-1)+1)=FIN(3*(J1-1)+1)-START(3*(J1-1)+1)
            HSVEL(3*(J1-1)+2)=FIN(3*(J1-1)+2)-START(3*(J1-1)+2)
            HSVEL(3*(J1-1)+3)=FIN(3*(J1-1)+3)-START(3*(J1-1)+3)
            DUMMY=MAX(ABS(HSVEL(3*(J1-1)+1)),ABS(HSVEL(3*(J1-1)+2)),ABS(HSVEL(3*(J1-1)+3)))
            HSVEL(3*(J1-1)+1)=HSVEL(3*(J1-1)+1) + 0.4D0*(DPRAND()-0.5D0)*DUMMY
            HSVEL(3*(J1-1)+2)=HSVEL(3*(J1-1)+2) + 0.4D0*(DPRAND()-0.5D0)*DUMMY
            HSVEL(3*(J1-1)+3)=HSVEL(3*(J1-1)+3) + 0.4D0*(DPRAND()-0.5D0)*DUMMY
         ENDDO
         T12FAC=1.1D0
         
         CALL HSMOVE(START,GUESS,HSVEL,STEP,.TRUE.)
         
         IF (MOD(NIMAGE,2).EQ.0) THEN
            DIST=0.0D0
            DO J1=1,LNOPT
               DIST=DIST+(GUESS(J1)-START(J1))**2
            ENDDO
            DIST=SQRT(DIST)
            SEPARATION=DIST/(NIMAGE/2+1)
            DO J1=1,NIMAGE/2
               DO J2=1,LNOPT
                  POINTS(J2+LNOPT*(J1-1+OFFSET))=START(J2)+(GUESS(J2)-START(J2))*J1*SEPARATION/DIST
               ENDDO
            ENDDO
            DIST=0.0D0
            DO J1=1,LNOPT
               DIST=DIST+(GUESS(J1)-FIN(J1))**2
            ENDDO
            DIST=SQRT(DIST)
            SEPARATION=DIST/(NIMAGE/2+1)
            DO J1=1,NIMAGE/2
               DO J2=1,LNOPT
                  POINTS(J2+LNOPT*(J1+NIMAGE/2-1+OFFSET))=GUESS(J2)+(FIN(J2)-GUESS(J2))*J1*SEPARATION/DIST
               ENDDO
            ENDDO
         ELSE
            DIST=0.0D0
            DO J1=1,LNOPT
               DIST=DIST+(GUESS(J1)-START(J1))**2
            ENDDO
            DIST=SQRT(DIST)
            SEPARATION=DIST/((NIMAGE-1)/2+1)
            DO J2=1,LNOPT
               POINTS(J2+LNOPT*((NIMAGE-1)/2+1-1+OFFSET))=GUESS(J2)
            ENDDO
            DO J1=1,(NIMAGE-1)/2
               DO J2=1,LNOPT
                  POINTS(J2+LNOPT*(J1-1+OFFSET))=START(J2)+(GUESS(J2)-START(J2))*J1*SEPARATION/DIST
               ENDDO
            ENDDO
            DIST=0.0D0
            DO J1=1,LNOPT
               DIST=DIST+(GUESS(J1)-FIN(J1))**2
            ENDDO
            DIST=SQRT(DIST)
            SEPARATION=DIST/((NIMAGE-1)/2+1)
            DO J1=1,(NIMAGE-1)/2
               DO J2=1,LNOPT
                  POINTS(J2+LNOPT*(J1+(NIMAGE-1)/2+1-1+OFFSET))=GUESS(J2)+(FIN(J2)-GUESS(J2))*J1*SEPARATION/DIST
               ENDDO
            ENDDO
         ENDIF
      ELSE
         DO J1=1,NIMAGE
            DO J2=1,LNOPT
               POINTS(J2+LNOPT*(J1-1+OFFSET))=START(J2)+(FIN(J2)-START(J2))*J1*SEPARATION/DIST
            ENDDO
C           PRINT*,'in neb image ',J1+OFFSET
C           IF (VARIABLES) THEN
C              PRINT'(2G20.10)',POINTS(LNOPT*(J1-1+OFFSET)+1),POINTS(LNOPT*(J1-1+OFFSET)+2)
C           ELSE
C              DO J2=1,LNATOMS
C                 WRITE(*,'(3F20.10)') POINTS(3*(J2-1)+1+LNOPT*(J1-1+OFFSET)),
C    1                                 POINTS(3*(J2-1)+2+LNOPT*(J1-1+OFFSET)),
C    1                                 POINTS(3*(J2-1)+3+LNOPT*(J1-1+OFFSET))
C              ENDDO
C           ENDIF
         ENDDO
      ENDIF
      DEALLOCATE(FINCART,STARTCART,POINTSCART)

      RETURN
      END
