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
      SUBROUTINE ODESD(ITMAX,VARS,MFLAG,NSTP,PTEST)
      USE COMMONS
      USE KEY
      USE MODHESS
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER NVAR, I, NOK, NBAD, NSTP, J1, NINFO, ITMAX
      LOGICAL MFLAG, PTEST
      DOUBLE PRECISION VARS(3*NATOMS), YSCAL(3*NATOMS), ENERGY, VEC2(3*NATOMS), DDOT,
     1                 DYDX(3*NATOMS), STRYDID, RMS, EPS, STRYNEXT, DUMMY, REALRMS, EREAL
      DOUBLE PRECISION TINY, STRYMIN, SLENGTH, MXSTPLOCAL
      PARAMETER (TINY=1.D-30)
      COMMON /BSNEW/ SLENGTH,NOK,NBAD,EPS
      INTEGER NSPECIAL, NALLOW, ISTAT
      DOUBLE PRECISION GSQSCALE, GSTHRESH
      LOGICAL PVFLAG
      COMMON /PVF/ PVFLAG
      COMMON /G2/ GSTHRESH, GSQSCALE, NSPECIAL, NALLOW, NINFO
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

      MXSTPLOCAL=MXSTP
      STRYMIN=0.0D0
      NSTP=1
      NVAR=3*NATOMS
      IF ((FIXAFTER.GT.0).AND.(NSTP.GT.FIXAFTER)) FIXIMAGE=.TRUE.
10    IF (PV.AND.(.NOT.BFGSSTEP)) THEN
         CALL POTENTIAL(VARS,ENERGY,DYDX,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         PVFLAG=.FALSE.
         CALL PVOPT(VARS,ENERGY,DYDX)
      ENDIF
!     IF (.NOT.KNOWG) CALL POTENTIAL(VARS,ENERGY,DYDX,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      CALL POTENTIAL(VARS,ENERGY,DYDX,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      PRINT '(A,G20.10)', 'EPS=',EPS
      CALL DUMPP(VARS,ENERGY) ! assumes NODUMP is always false if we are calling odesd.f
      EREAL=ENERGY
      REALRMS=RMS

      IF (GRADSQ) THEN
         ENERGY=DDOT(3*NATOMS,DYDX,1,DYDX,1)
C        ENERGY=SQRT(DDOT(3*NATOMS,DYDX,1,DYDX,1))
         CALL DSYMV('U',3*NATOMS,2.0D0,HESS,SIZE(HESS,1),DYDX,1,0.0D0,VEC2,1)
         RMS= DSQRT(DDOT(3*NATOMS,VEC2,1,VEC2,1)/(3*NATOMS))
         IF (RMS.LT.GSTHRESH) THEN
            IF ((MOD(NSTP+1,NSPECIAL).EQ.0).AND.(NSPECIAL.GT.0)) THEN
               CALL G2SPECIAL(NATOMS,VARS,DYDX,VEC2,ENERGY,RMS,EREAL,REALRMS,MFLAG)
            ELSE IF (FIXAFTER.EQ.0) THEN
               FIXAFTER=NSTP+2
               PRINT*,' changing to atom type LC and setting FIXAFTER=',NSTP+2
               DO J1=1,NATOMS
                  ZSYM(J1)='LC'
               ENDDO
               RETURN
            ENDIF
         ENDIF
C        IF (GRADSQ) WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,RMS
         IF (MFLAG) THEN
            IF (GRADSQ) WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
            RETURN
         ENDIF
         DO J1=1,3*NATOMS
            DYDX(J1)=VEC2(J1)
C           DYDX(J1)=VEC2(J1)/ENERGY
         ENDDO
      ENDIF
      IF (GRADSQ) THEN
         IF (PTEST) WRITE(*,'(A,4F20.10,A,I6,2(A,F15.5))') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS,
     1                                     ' after ',NSTP-1,' SD steps, path length=',SLENGTH,' step=',MXSTPLOCAL
      ELSE
         IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,2(A,F15.5))') ' Energy and RMS force=',ENERGY,REALRMS,
     1                                     ' after ',NSTP-1,' SD steps, path length=',SLENGTH,' step=',MXSTPLOCAL
      ENDIF
      CALL FLUSH(6)

      IF ((RMS.LT.GMAX).AND.(NSTP.GT.NSTEPMIN)) THEN
         MFLAG=.TRUE.
         NSTP=NSTP-1
         IF (NSTP.LT.NSTEPMIN) MFLAG=.FALSE.
         IF (PV.AND.(.NOT.PVFLAG)) MFLAG=.FALSE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (GRADSQ) WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
            RETURN
         ENDIF
      ENDIF

      IF (NSTP.EQ.ITMAX) THEN
         MFLAG=.FALSE.
         FIXIMAGE=.FALSE.
         RETURN
      ENDIF

      DUMMY=RMS*SQRT(1.0D0*NVAR)
      DO I=1,NVAR
         DYDX(I)=-DYDX(I)/DUMMY
         YSCAL(I)=ABS(VARS(I))+ABS(MXSTPLOCAL*DYDX(I))+TINY
C        DYDX(I)=-DYDX(I)
      ENDDO

      IF (BSMIN) CALL BSSTEP(VARS,DYDX,NVAR,SLENGTH,MXSTPLOCAL,EPS,YSCAL,STRYDID,STRYNEXT,ENERGY,RMS,EREAL,REALRMS,PTEST)
      IF (RKMIN) CALL RKQS(VARS,DYDX,NVAR,SLENGTH,MXSTPLOCAL,EPS,YSCAL,STRYDID,STRYNEXT,ENERGY,RMS,EREAL,REALRMS,PTEST)

      WRITE(*,'(A,F20.10,A,F20.10)') 'Step length=',STRYDID,' next estimated step size=',STRYNEXT

      IF (STRYDID.GE.MXSTPLOCAL) THEN
         NOK=NOK+1
      ELSE
         NBAD=NBAD+1
      ENDIF
      PRINT*,'STRYDID,MXSTPLOCAL,NOK,NBAD=',STRYDID,MXSTPLOCAL,NOK,NBAD
      IF (ABS(STRYNEXT).LT.STRYMIN) THEN
         PRINT*, ' WARNING, stepsize < 0 in odesd'
         OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
         PRINT*,' intractable discontinuity - quit '
         WRITE(96,'(A)') 'intractable discontinuity'
         CLOSE(96)
         IF (GRADSQ) THEN 
            WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
         ELSE
            WRITE(*,'(A,2F20.10,A,I6,A,F15.10)') ' Energy and RMS force=',ENERGY,REALRMS
         ENDIF
         CALL DUMPIT(VARS,'points.final')
         STOP
      ENDIF

      MXSTPLOCAL=STRYNEXT
      MFLAG=.FALSE.
      FIXIMAGE=.FALSE.
      NSTP=NSTP+1
      GOTO 10

      RETURN
      END

      SUBROUTINE BSSTEP(Y,DYDX,NV,X,HTRY,EPS,YSCAL,HDID,HNEXT,ENERGY,RMS,EREAL,REALRMS,PTEST)
      USE PORFUNCS
      USE MODHESS
      USE KEY
      IMPLICIT NONE
      INTEGER NV,KMAXX,IMAX
      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,DYDX(NV),Y(NV),YSCAL(NV)
     *,SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX,DDOT,RMS
      PARAMETER (KMAXX=8,IMAX=KMAXX+1,SAFE1=.25D0,SAFE2=.7D0,
     *REDMAX=1.D-5,REDMIN=.7D0,TINY=1.D-30,SCALMX=.1D0)
CU    USES derivs,mmid,pzextr
      INTEGER I,IQ,K,KK,KM,KMAX,KOPT,NSEQ(IMAX)
      DOUBLE PRECISION EPS1,EPSOLD,ERRMAX,FACT,H,RED,SCALE,WORK,WRKMIN
     *,XEST,XNEW,ENERGY,EOLD,VNEW(NV),
     *A(IMAX),ALF(KMAXX,KMAXX),ERR(KMAXX),YERR(NV),YSAV(NV),EREAL,REALRMS,
     *YSEQ(NV)
      LOGICAL FIRST,REDUCT,PTEST
      SAVE A,ALF,EPSOLD,FIRST,KMAX,KOPT,NSEQ,XNEW
      INTEGER IMAXSAVE
      PARAMETER (IMAXSAVE=13)
      DOUBLE PRECISION DSAVE(NV,IMAXSAVE),FXSAVE(IMAXSAVE),XSAVE(IMAXSAVE),QCOLSAVE(NV,IMAXSAVE), EPSSAVE

      DATA FIRST/.TRUE./,EPSOLD/-1.D0/
      DATA NSEQ /2,4,6,8,10,12,14,16,18/

C     PRINT*,'in bsstep'
      EOLD=ENERGY
      EPSSAVE=EPS
      IF(EPS.NE.EPSOLD)THEN
        HNEXT=-1.D29
        XNEW=-1.D29
        EPS1=SAFE1*EPS
        A(1)=NSEQ(1)+1
        DO 11 K=1,KMAXX
          A(K+1)=A(K)+NSEQ(K+1)
11      CONTINUE
        DO 13 IQ=2,KMAXX
          DO 12 K=1,IQ-1
            ALF(K,IQ)=EPS1**((A(K+1)-A(IQ+1))/((A(IQ+1)-A(1)+1.D0)*(2*K+1)))
12        CONTINUE
13      CONTINUE
        EPSOLD=EPS
        DO 14 KOPT=2,KMAXX-1
          IF(A(KOPT+1).GT.A(KOPT)*ALF(KOPT-1,KOPT))GOTO 1
14      CONTINUE
1       KMAX=KOPT
      ENDIF
      H=HTRY
      DO 15 I=1,NV
        YSAV(I)=Y(I)
15    CONTINUE
      IF(H.NE.HNEXT.OR.X.NE.XNEW)THEN
        FIRST=.TRUE.
        KOPT=KMAX
      ENDIF
      REDUCT=.FALSE.
2     DO 17 K=1,KMAX
        XNEW=X+H
        IF ((XNEW.EQ.X).OR.(EPS.LT.1.0D-20)) THEN
           PRINT*,'stepsize underflow in bsstep'
           OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
           PRINT*,' intractable discontinuity - quit '
           WRITE(96,'(A)') 'intractable discontinuity'
           CLOSE(96)
           WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
           CALL DUMPIT(Y,'points.final')
           STOP
        ENDIF
        CALL MMID(YSAV,DYDX,NV,X,H,NSEQ(K),YSEQ)
        XEST=(H/NSEQ(K))**2
        CALL RZEXTR(K,XEST,YSEQ,Y,YERR,NV,DSAVE,XSAVE,FXSAVE)
C       CALL PZEXTR(K,XEST,YSEQ,Y,YERR,NV,QCOLSAVE,XSAVE)
        IF(K.NE.1)THEN
          ERRMAX=TINY
          DO 16 I=1,NV
            ERRMAX=MAX(ERRMAX,ABS(YERR(I)/YSCAL(I)))
16        CONTINUE
          ERRMAX=ERRMAX/EPS
          KM=K-1
          ERR(KM)=(ERRMAX/SAFE1)**(1.D0/(2*KM+1))
        ENDIF
        IF(K.NE.1.AND.(K.GE.KOPT-1.OR.FIRST))THEN
          IF (ERRMAX.LT.1.0D0) THEN
!            CALL POTENTIAL(Y,ENERGY,VNEW,GRADSQ,.FALSE.,RMS,.FALSE.,.FALSE.)
             CALL POTENTIAL(Y,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
             EREAL=ENERGY
             REALRMS=RMS
             IF (GRADSQ) THEN
                ENERGY=DDOT(NV,VNEW,1,VNEW,1)
C               CALL DSYMV('U',NV,2.0D0,HESS,NV,VNEW,1,0.0D0,VEC2,1)
C               RMS= DSQRT(DDOT(NV,VEC2,1,VEC2,1)/(NV))
C               WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
             ENDIF
!            IF (ENERGY.GT.EOLD) THEN
             IF (ENERGY-EOLD.GT.MAXERISE) THEN
                ERRMAX=2.0D0 ! reject step
                EPS=MAX(EPS/2.0D0,1.0D-10)
                IF (PTEST) WRITE(*,'(3(A,G20.10))') ' bsstep> energy rises from ',EOLD,' to ',ENERGY,
     &                                               ' reduce eps to ',EPS
             ELSE
!               EPS=MIN(10.0D0,EPS*1.5D0)
                DYDX(1:NV)=VNEW(1:NV) ! pass the energy and gradient back and we don;t have to call potential again
             ENDIF
          ENDIF
          IF (ERRMAX.LT.1.D0) GOTO 4 ! step accepted
          IF(K.EQ.KMAX.OR.K.EQ.KOPT+1)THEN
            RED=SAFE2/ERR(KM)
            GOTO 3
          ELSE IF(K.EQ.KOPT)THEN
            IF(ALF(KOPT-1,KOPT).LT.ERR(KM))THEN
              RED=1.D0/ERR(KM)
              GOTO 3
            ENDIF
          ELSE IF(KOPT.EQ.KMAX)THEN
            IF(ALF(KM,KMAX-1).LT.ERR(KM))THEN
              RED=ALF(KM,KMAX-1)*SAFE2/ERR(KM)
              GOTO 3
            ENDIF
          ELSE IF(ALF(KM,KOPT).LT.ERR(KM))THEN
            RED=ALF(KM,KOPT-1)/ERR(KM)
            GOTO 3
          ENDIF
        ENDIF
17    CONTINUE
3     RED=MIN(RED,REDMIN)
      RED=MAX(RED,REDMAX)
      H=H*RED
      REDUCT=.TRUE.
      GOTO 2
4     X=XNEW
      HDID=H
      FIRST=.FALSE.
      WRKMIN=1.D35
      DO 18 KK=1,KM
        FACT=MAX(ERR(KK),SCALMX)
        WORK=FACT*A(KK+1)
        IF(WORK.LT.WRKMIN)THEN
          SCALE=FACT
          WRKMIN=WORK
          KOPT=KK+1
        ENDIF
18    CONTINUE
      HNEXT=H/SCALE
      IF(KOPT.GE.K.AND.KOPT.NE.KMAX.AND..NOT.REDUCT)THEN
        FACT=MAX(SCALE/ALF(KOPT-1,KOPT),SCALMX)
        IF(A(KOPT+1)*FACT.LE.WRKMIN)THEN
          HNEXT=H/FACT
          KOPT=KOPT+1
        ENDIF
      ENDIF
!     EPS=EPSSAVE
      RETURN
      END

      SUBROUTINE MMID(Y,DYDX,NVAR,XS,HTOT,NSTEP,YOUT)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER NSTEP,NVAR
      DOUBLE PRECISION HTOT,XS,DYDX(NVAR),Y(NVAR),YOUT(NVAR),ENERGY
      INTEGER I,N
      DOUBLE PRECISION H,H2,SWAP,X,YM(NVAR),YN(NVAR),RMS,EREAL,DDOT,VEC2(NVAR)

      H=HTOT/NSTEP
      DO 11 I=1,NVAR
        YM(I)=Y(I)
        YN(I)=Y(I)+H*DYDX(I)
11    CONTINUE
      X=XS+H
      CALL POTENTIAL(YN,ENERGY,YOUT,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      IF (GRADSQ) THEN
         EREAL=ENERGY
         ENERGY=DDOT(NVAR,YOUT,1,YOUT,1)
C        ENERGY=SQRT(DDOT(NVAR,YOUT,1,YOUT,1))
         CALL DSYMV('U',NVAR,2.0D0,HESS,NVAR,YOUT,1,0.0D0,VEC2,1)
         DO I=1,NVAR
            YOUT(I)=VEC2(I)
CC          YOUT(I)=VEC2(I)/ENERGY
         ENDDO
         RMS=DSQRT(DDOT(NVAR,VEC2,1,VEC2,1)/(NVAR))
      ENDIF
      RMS=RMS*SQRT(1.0D0*NVAR)
      DO I=1,NVAR
         YOUT(I)=-YOUT(I)/RMS
C        YOUT(I)=-YOUT(I)
      ENDDO
      H2=2.D0*H
      DO 13 N=2,NSTEP
        DO 12 I=1,NVAR
          SWAP=YM(I)+H2*YOUT(I)
          YM(I)=YN(I)
          YN(I)=SWAP
12      CONTINUE
        X=X+H
        CALL POTENTIAL(YN,ENERGY,YOUT,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
         IF (GRADSQ) THEN
            EREAL=ENERGY
            ENERGY=DDOT(NVAR,YOUT,1,YOUT,1)
C           ENERGY=SQRT(DDOT(NVAR,YOUT,1,YOUT,1))
            CALL DSYMV('U',NVAR,2.0D0,HESS,NVAR,YOUT,1,0.0D0,VEC2,1)
            DO I=1,NVAR
               YOUT(I)=VEC2(I)
C              YOUT(I)=VEC2(I)/ENERGY
            ENDDO
            RMS=DSQRT(DDOT(NVAR,VEC2,1,VEC2,1)/(NVAR))
         ENDIF
         RMS=RMS*SQRT(1.0D0*NVAR)
         DO I=1,NVAR
            YOUT(I)=-YOUT(I)/RMS
C           YOUT(I)=-VEC2(I)
         ENDDO
13    CONTINUE
      DO 14 I=1,NVAR
        YOUT(I)=0.5D0*(YM(I)+YN(I)+H*YOUT(I))
14    CONTINUE
      RETURN
      END

!
!  QCOL and X could not be saved properly - have to pass them as arguments
!
      SUBROUTINE PZEXTR(IEST,XEST,YEST,YZ,DY,NV,QCOL,X)
      IMPLICIT NONE
      INTEGER IEST,NV,IMAX
      DOUBLE PRECISION XEST,DY(NV),YEST(NV),YZ(NV)
      PARAMETER (IMAX=13)
      INTEGER J,K1
      DOUBLE PRECISION DELTA,F1,F2,Q,D(NV),QCOL(NV,IMAX),X(IMAX)
C     SAVE qcol,x

C     PRINT*,'in pzextr'
      X(IEST)=XEST
      DO 11 J=1,NV
        DY(J)=YEST(J)
        YZ(J)=YEST(J)
11    CONTINUE
      IF(IEST.EQ.1) THEN
        DO 12 J=1,NV
          QCOL(J,1)=YEST(J)
12      CONTINUE
      ELSE
        DO 13 J=1,NV
          D(J)=YEST(J)
13      CONTINUE
        DO 15 K1=1,IEST-1
          DELTA=1.D0/(X(IEST-K1)-XEST)
          F1=XEST*DELTA
          F2=X(IEST-K1)*DELTA
          DO 14 J=1,NV
            Q=QCOL(J,K1)
            QCOL(J,K1)=DY(J)
            DELTA=D(J)-Q
            DY(J)=F1*DELTA
            D(J)=F2*DELTA
            YZ(J)=YZ(J)+DY(J)
14        CONTINUE
15      CONTINUE
        DO 16 J=1,NV
          QCOL(J,IEST)=DY(J)
16      CONTINUE
      ENDIF
      RETURN
      END

!
!  SAVE statement alone does NOT save D and X !
!  SAVE X cannot work because X is an automatic variable
!  hence we must pass these arrays and declare them in the calling routine
!
      SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,DY,NV,D,X,FX)
      IMPLICIT NONE
      INTEGER IEST,NV,IMAX
      DOUBLE PRECISION XEST,DY(NV),YEST(NV),YZ(NV)
      PARAMETER (IMAX=13)
      INTEGER J,K
      DOUBLE PRECISION B,B1,C,DDY,V,YY,D(NV,IMAX),FX(IMAX),X(IMAX)
C     SAVE d,x

      X(IEST)=XEST
      IF(IEST.EQ.1) THEN
        DO 11 J=1,NV
          YZ(J)=YEST(J)
          D(J,1)=YEST(J)
          DY(J)=YEST(J)
11      CONTINUE
      ELSE
        DO 12 K=1,IEST-1
          FX(K+1)=X(IEST-K)/XEST
12      CONTINUE
        DO 14 J=1,NV
          YY=YEST(J)
          V=D(J,1)
          C=YY
          D(J,1)=YY
          DO 13 K=2,IEST
            B1=FX(K)*V
            B=B1-C
            IF(B.NE.0.D0) THEN
              B=(C-V)/B
              DDY=C*B
              C=B1*B
            ELSE
              DDY=V
            ENDIF
            IF (K.NE.IEST) V=D(J,K)
            D(J,K)=DDY
            YY=YY+DDY
13        CONTINUE
          DY(J)=DDY
          YZ(J)=YY
14      CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE RKQS(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,ENERGY,RMS,EREAL,REALRMS,PTEST)
      USE PORFUNCS
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,DYDX(N),Y(N),YSCAL(N)
CU    USES rkck
      INTEGER I
      DOUBLE PRECISION ERRMAX,H,XNEW,YERR(N),YTEMP(N),SAFETY,PGROW
     *,PSHRNK,
     *ERRCON
      PARAMETER (SAFETY=0.9D0,PGROW=-.2D0,PSHRNK=-.25D0,ERRCON=1.89D-4)
      DOUBLE PRECISION EOLD, ENERGY, RMS, DDOT, VNEW(N), EREAL, REALRMS
      LOGICAL PTEST
      
      EOLD=ENERGY
      H=HTRY
1     CALL RKCK(Y,DYDX,N,X,H,YTEMP,YERR)
      ERRMAX=0.D0
      DO 11 I=1,N
        ERRMAX=MAX(ERRMAX,ABS(YERR(I)/YSCAL(I)))
11    CONTINUE
      ERRMAX=ERRMAX/EPS

      IF (ERRMAX.LT.1.0D0) THEN
         CALL POTENTIAL(YTEMP,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         EREAL=ENERGY
         REALRMS=RMS
         IF (GRADSQ) THEN
            ENERGY=DDOT(N,VNEW,1,VNEW,1)
C           CALL DSYMV('U',N,2.0D0,HESS,N,VNEW,1,0.0D0,VEC2,1)
C           RMS= DSQRT(DDOT(N,VEC2,1,VEC2,1)/(N))
C           WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
         ENDIF
!        IF (ENERGY.GT.EOLD) THEN
         IF (ENERGY-EOLD.GT.MAXERISE) THEN
            ERRMAX=2.0D0
            EPS=EPS/2.0D0
            IF (PTEST) WRITE(*,'(3(A,G20.10))') ' rkqs> energy rises from ',EOLD,' to ',ENERGY,
     &                                               ' reduce eps to ',EPS
         ELSE
C           EPS=EPS*1.1D0
            DYDX(1:N)=VNEW(1:N) ! pass the energy and gradient back and we don;t have to call potential again
         ENDIF
      ENDIF

      IF(ERRMAX.GT.1.D0)THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        IF(H.LT.0.1D0*H)THEN
          H=.1D0*H
        ENDIF
        XNEW=X+H
        IF ((XNEW.EQ.X).OR.(EPS.LT.1.0D-20)) THEN
           PRINT*,' stepsize underflow in rkqs'
           OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
           PRINT*,' intractable discontinuity - quit '
           WRITE(96,'(A)') 'intractable discontinuity'
           CLOSE(96)
           WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
           CALL DUMPIT(Y,'points.final')
           STOP
        ENDIF
        GOTO 1
      ELSE
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=5.D0*H
        ENDIF
        HDID=H
        X=X+H
        DO 12 I=1,N
          Y(I)=YTEMP(I)
12      CONTINUE
        RETURN
      ENDIF
      END

      SUBROUTINE RKCK(Y,DYDX,N,X,H,YOUT,YERR)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION H,X,DYDX(N),Y(N),YERR(N),YOUT(N)
      INTEGER I
      DOUBLE PRECISION AK2(N),AK3(N),AK4(N),AK5(N),AK6(N)
     *,
     *YTEMP(N),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2D0,A3=.3D0,A4=.6D0,A5=1.D0,A6=.875D0,B21=.2D0,B31
     *=3.D0/40.D0,
     *B32=9.D0/40.D0,B41=.3D0,B42=-.9D0,B43=1.2D0,B51=-11.D0/54.D0,B52
     *=2.5D0,
     *B53=-70.D0/27.D0,B54=35.D0/27.D0,B61=1631.D0/55296.D0,B62=175.D0
     */512.D0,
     *B63=575.D0/13824.D0,B64=44275.D0/110592.D0,B65=253.D0/4096.D0,C1
     *=37.D0/378.D0,
     *C3=250.D0/621.D0,C4=125.D0/594.D0,C6=512.D0/1771.D0,DC1=C1-2825.D0
     */27648.D0,
     *DC3=C3-18575.D0/48384.D0,DC4=C4-13525.D0/55296.D0,DC5=-277.D0
     */14336.D0,
     *DC6=C6-.25D0)
      DOUBLE PRECISION ENERGY, RMS, VEC2(N), EREAL, DDOT

      DO 11 I=1,N
        YTEMP(I)=Y(I)+B21*H*DYDX(I)
11    CONTINUE
C     call derivs(x+A2*h,ytemp,ak2)

      CALL POTENTIAL(YTEMP,ENERGY,AK2,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      IF (GRADSQ) THEN
         EREAL=ENERGY
         ENERGY=DDOT(N,AK2,1,AK2,1)
         CALL DSYMV('U',N,2.0D0,HESS,N,AK2,1,0.0D0,VEC2,1)
         RMS=DSQRT(DDOT(N,VEC2,1,VEC2,1)/(N))
         DO I=1,N
            AK2(I)=VEC2(I)
         ENDDO
      ENDIF

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
C        AK2(I)=-AK2(I)
         AK2(I)=-AK2(I)/RMS
      ENDDO

      DO 12 I=1,N
        YTEMP(I)=Y(I)+H*(B31*DYDX(I)+B32*AK2(I))
12    CONTINUE
C     call derivs(x+A3*h,ytemp,ak3)

      CALL POTENTIAL(YTEMP,ENERGY,AK3,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      IF (GRADSQ) THEN
         EREAL=ENERGY
         ENERGY=DDOT(N,AK3,1,AK3,1)
         CALL DSYMV('U',N,2.0D0,HESS,N,AK3,1,0.0D0,VEC2,1)
         RMS=DSQRT(DDOT(N,VEC2,1,VEC2,1)/(N))
         DO I=1,N
            AK3(I)=VEC2(I)
         ENDDO
      ENDIF

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
C        AK3(I)=-AK3(I)
         AK3(I)=-AK3(I)/RMS
      ENDDO

      DO 13 I=1,N
        YTEMP(I)=Y(I)+H*(B41*DYDX(I)+B42*AK2(I)+B43*AK3(I))
13    CONTINUE
C     call derivs(x+A4*h,ytemp,ak4)

      CALL POTENTIAL(YTEMP,ENERGY,AK4,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      IF (GRADSQ) THEN
         EREAL=ENERGY
         ENERGY=DDOT(N,AK4,1,AK4,1)
         CALL DSYMV('U',N,2.0D0,HESS,N,AK4,1,0.0D0,VEC2,1)
         RMS=DSQRT(DDOT(N,VEC2,1,VEC2,1)/(N))
         DO I=1,N
            AK4(I)=VEC2(I)
         ENDDO
      ENDIF

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
C        AK4(I)=-AK4(I)
         AK4(I)=-AK4(I)/RMS
      ENDDO

      DO 14 I=1,N
        YTEMP(I)=Y(I)+H*(B51*DYDX(I)+B52*AK2(I)+B53*AK3(I)+B54*AK4(I))
14    CONTINUE
C     call derivs(x+A5*h,ytemp,ak5)

      CALL POTENTIAL(YTEMP,ENERGY,AK5,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      IF (GRADSQ) THEN
         EREAL=ENERGY
         ENERGY=DDOT(N,AK5,1,AK5,1)
         CALL DSYMV('U',N,2.0D0,HESS,N,AK5,1,0.0D0,VEC2,1)
         RMS=DSQRT(DDOT(N,VEC2,1,VEC2,1)/(N))
         DO I=1,N
            AK5(I)=VEC2(I)
         ENDDO
      ENDIF

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
C        AK5(I)=-AK5(I)
         AK5(I)=-AK5(I)/RMS
      ENDDO

      DO 15 I=1,N
        YTEMP(I)=Y(I)+H*(B61*DYDX(I)+B62*AK2(I)+B63*AK3(I)+B64*AK4(I)+
     *B65*AK5(I))
15    CONTINUE
C     call derivs(x+A6*h,ytemp,ak6)

      CALL POTENTIAL(YTEMP,ENERGY,AK6,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
      IF (GRADSQ) THEN
         EREAL=ENERGY
         ENERGY=DDOT(N,AK6,1,AK6,1)
         CALL DSYMV('U',N,2.0D0,HESS,N,AK6,1,0.0D0,VEC2,1)
         RMS=DSQRT(DDOT(N,VEC2,1,VEC2,1)/(N))
         DO I=1,N
            AK6(I)=VEC2(I)
         ENDDO
      ENDIF

      RMS=RMS*SQRT(1.0D0*N)

      DO I=1,N
C        AK6(I)=-AK6(I)
         AK6(I)=-AK6(I)/RMS
      ENDDO

      DO 16 I=1,N
        YOUT(I)=Y(I)+H*(C1*DYDX(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
16    CONTINUE
      DO 17 I=1,N
        YERR(I)=H*(DC1*DYDX(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*
     *AK6(I))
17    CONTINUE
      RETURN
      END

