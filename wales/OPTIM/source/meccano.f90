!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE MECCANO(ITEST,PTEST,EMAX,PERM,Q,FINAL,E1,E2,RMSINITIAL,RMSFINAL)
USE COMMONS
use porfuncs
USE KEY
USE MODNEB
USE MODTWOEND
USE MODMEC
use KeyNEB, only: nitermax
IMPLICIT NONE
INTEGER :: J1, J2, ITDONE, JMAX, NVARS, NFROZEN, NPASS
DOUBLE PRECISION, ALLOCATABLE :: XVAR(:), LEVEC(:), LRMSVEC(:)
DOUBLE PRECISION :: RMS, EMAX, TIME, TIME0, Q(3*NATOMS), FINAL(3*NATOMS), &
     &              RMSVEC(0:Nimage+1), EVEC(0:Nimage+1), E1, E2, &
     &              POINTS(3*NATOMS*Nimage+1), SEPARATION, EINITIAL, EFINAL, RMSINITIAL, RMSFINAL, &
     &              QSTART(3*NATOMS), QFINISH(3*NATOMS), &
     &              GRAD(3*NATOMS), &
     &              EFORIG, EIORIG, FINORIG(3*NATOMS)
LOGICAL :: MFLAG, ITEST, PTEST, PERM, REOPTIM(Nimage)
COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
CHARACTER FNAME*19
INTEGER IADD, LNOPT, LNATOMS
COMMON /INTS/ IADD, LNOPT, LNATOMS
CHARACTER(LEN=5) ZSYMSAVE
CHARACTER(LEN=80) LFNAME
COMMON /SYS/ ZSYMSAVE

IF (.NOT.ITEST) THEN
   EINITIAL=E1
   EFINAL=E2
ELSE
   CALL POTENTIAL(Q,EINITIAL,GRAD,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
   CALL POTENTIAL(FINAL,EFINAL,GRAD,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
ENDIF
IF (ITEST.AND.PTEST) WRITE(*,'(A,F20.10,A,F15.6,A,F20.10,A,F15.6)') &
     &           ' Initial energy=',EINITIAL,' RMS force=',RMSINITIAL,' final energy=',EFINAL,' RMS=',RMSFINAL


IF (ZSYMSAVE(1:1).EQ.'W') THEN
   PRINT '(A)','MECCANO not adapted for rigid bodies'
   STOP
!  IADD=3*(NATOMS/2)
!  LNOPT=NOPT/2      ! we only want to optimise the angular degrees of freedom for TIP
!  LNATOMS=NATOMS/2  ! for TIP potentials both centre-of-mass and Euler angles are now in odata
ELSE
   IADD=0
   LNOPT=NOPT
   LNATOMS=NATOMS
ENDIF
CALL MYCPU_TIME(TIME0,.TRUE.)
EIORIG=EINITIAL
EFORIG=EFINAL
FINORIG(1:LNOPT+IADD)=FINAL(1:LNOPT+IADD)
EVEC(0)=EINITIAL
EVEC(NIMAGE+1)=EFINAL
RMSVEC(0)=RMSINITIAL
RMSVEC(NIMAGE+1)=RMSFINAL
FIN(1:LNOPT+IADD)=FINAL(1:LNOPT+IADD) ! needed by makeimage in oldneb
CALL MAKEIMAGE(POINTS,PERM,ITEST,0,Q) ! linear interpolation
! final variable for Lagrangian minimisation is the distance d.
IF (VARIABLES) THEN
   POINTS(LNOPT*NIMAGE+1)=MECDIST
ELSE
   POINTS(3*NATOMS*NIMAGE+1)=MECDIST
ENDIF
!
!  Minimisation of Lagrangian
!
IF (FILTH.EQ.0) THEN
   WRITE(FNAME,'(A8)') 'EofS.neb'
ELSE 
   WRITE(FNAME,'(A)') 'EofS.neb.'//TRIM(ADJUSTL(FILTHSTR))
ENDIF
NVARS=(LNOPT+IADD)*NIMAGE+1
ALLOCATE(LEVEC(NIMAGE),LRMSVEC(NIMAGE))

! first pass - nothing frozen

CALL MECBFGS(NVARS,MECUPDATE,POINTS,.FALSE.,MECRMSTOL,MFLAG,LEVEC,LRMSVEC, &
 &             NIterMax,ITDONE,PTEST,Q,FINAL)
EVEC(1:NIMAGE)=LEVEC(1:NIMAGE)
RMSVEC(1:NIMAGE)=LRMSVEC(1:NIMAGE)
DEALLOCATE(LEVEC,LRMSVEC)
! GOTO 12
! optimise images corresponding to TS candidates further 
! while freezing all others 

NVARS=LNOPT+IADD+1
ALLOCATE(XVAR(LNOPT+IADD+1),LEVEC(1),LRMSVEC(1))

NFROZEN=NIMAGE-1
NIMAGE=1
NPASS=0
11 NPASS=NPASS+1
PRINT*,'NPASS=',NPASS
REOPTIM(1:NIMAGE+NFROZEN)=.FALSE.
DO J1=1,NIMAGE+NFROZEN
   IF ((EVEC(J1).GT.EVEC(J1-1)).AND.(EVEC(J1).GT.EVEC(J1+1))) REOPTIM(J1)=.TRUE.
ENDDO
DO J1=1,NIMAGE+NFROZEN
   IF (REOPTIM(J1)) THEN
      IF (PTEST) PRINT*,'EVEC(J1-1)=',J1-1,EVEC(J1-1)
      IF (PTEST) PRINT*,'EVEC(J1)  =',J1,EVEC(J1)
      IF (PTEST) PRINT*,'EVEC(J1+1)=',J1+1,EVEC(J1+1)
      IF (PTEST) PRINT '(A,I5)',' reoptimising image ',J1
      IF (.NOT.PTEST) PRINT '(A,I5)',' reoptimising image ',J1
      IF (J1.EQ.1) THEN
         QSTART(1:LNOPT+IADD)=Q(1:LNOPT+IADD)
      ELSE
         QSTART(1:LNOPT+IADD)=POINTS((J1-2)*(LNOPT+IADD)+1:(J1-1)*(LNOPT+IADD))
      ENDIF
      XVAR(1:LNOPT+IADD)=POINTS((J1-1)*(LNOPT+IADD)+1:J1*(LNOPT+IADD))
      XVAR(LNOPT+IADD+1)=POINTS((NIMAGE+NFROZEN)*(LNOPT+IADD)+1)
      IF (J1.EQ.NIMAGE+NFROZEN) THEN
         QFINISH(1:LNOPT+IADD)=FINAL(1:LNOPT+IADD)
      ELSE
         QFINISH(1:LNOPT+IADD)=POINTS(J1*(LNOPT+IADD)+1:(J1+1)*(LNOPT+IADD))
      ENDIF
      CALL MECBFGS(NVARS,MECUPDATE,XVAR,.FALSE.,MECRMSTOL,MFLAG,LEVEC,LRMSVEC, &
  &             NIterMax,ITDONE,PTEST,QSTART,QFINISH)
      POINTS((J1-1)*(LNOPT+IADD)+1:J1*(LNOPT+IADD))=XVAR(1:LNOPT+IADD)
      POINTS((NIMAGE+NFROZEN)*(LNOPT+IADD)+1)=XVAR(LNOPT+IADD+1)
      EVEC(J1)=LEVEC(1)
      RMSVEC(J1)=LRMSVEC(1)
   ENDIF
ENDDO 
IF (NPASS.LE.0) GOTO 11
NIMAGE=NIMAGE+NFROZEN
! 12 CONTINUE
!
!  Find the highest energy image.
!
JMAX=1
EMAX=-1.0D100
DO J1=1,NIMAGE
   IF (EVEC(J1).GT.EMAX) THEN
!
!  Don;t choose a bad geometry for DFTB.
!
      IF ((.NOT.DFTBT).OR.(DFTBT.AND.EVEC(J1).LT.1.0D6)) THEN
         EMAX=EVEC(J1)
         JMAX=J1
      ENDIF
   ENDIF
ENDDO
!
!  Dump pathway as an xyz file
!
IF (FILTH.EQ.0) THEN
   WRITE(FNAME,'(A12)') 'neb.path.xyz'
ELSE 
   WRITE(FNAME,'(A)') 'neb.path.xyz.'//TRIM(ADJUSTL(FILTHSTR))
ENDIF
OPEN(UNIT=3,FILE=FNAME,STATUS='UNKNOWN')
WRITE(3,'(I6)') LNATOMS
WRITE(3,'(F20.10)') EINITIAL
WRITE(3,'(A2,4X,3F20.10)') (ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),J1=1,LNATOMS)
DO J1=1,NIMAGE
   WRITE(3,'(I6)') LNATOMS
   WRITE(3,'(F20.10)') EVEC(J1) 
   WRITE(3,'(A2,4X,3F20.10)') &
   &  (ZSYM(J2),POINTS(LNOPT*(J1-1)+3*(J2-1)+1),POINTS(LNOPT*(J1-1)+3*(J2-1)+2),POINTS(LNOPT*(J1-1)+3*(J2-1)+3),J2=1,LNATOMS)
ENDDO
WRITE(3,'(I6)') LNATOMS
WRITE(3,'(F20.10)') EFINAL
WRITE(3,'(A2,4X,3F20.10)') (ZSYM(J1),FINAL(3*(J1-1)+1),FINAL(3*(J1-1)+2),FINAL(3*(J1-1)+3),J1=1,LNATOMS)
CLOSE(3)
CALL MYCPU_TIME(TIME,.FALSE.)
WRITE(*,'(A,I4,A,G15.5,A,I6,A,F12.4,A,I4,A,F11.2)') ' image ',JMAX,' has highest energy=',EMAX, &
     &                   ' steps=',ITDONE,' RMS=',RMS,' images=',NIMAGE,'    time=',TIME-TIME0
TIME0=TIME
!
!  Write images to guess.xyz as candidate transition states
!
IF (FILTH.EQ.0) THEN
   OPEN(UNIT=7,FILE='guess.xyz',STATUS='UNKNOWN')
ELSE
   LFNAME='guess.xyz.'//TRIM(ADJUSTL(FILTHSTR)) ! Sun compiler fails without this intermediate step
   OPEN(UNIT=7,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
ENDIF
DO J1=0,NIMAGE+1
   IF (J1.EQ.0) THEN
      WRITE(7,'(3G20.10)') (Q(J2),J2=1,LNOPT)
   ELSE IF (J1.EQ.NIMAGE+1) THEN
      WRITE(7,'(3G20.10)') (FINAL(J2),J2=1,LNOPT)
   ELSE
      WRITE(7,'(3G20.10)') (POINTS((J1-1)*(LNOPT+IADD)+J2),J2=1,LNOPT)
   ENDIF
ENDDO
CLOSE(7)
IF (VARIABLES) THEN
    OPEN(UNIT=7,FILE='variables.images',STATUS='UNKNOWN')   
    WRITE(7,'(3G20.10)') (Q(J2),J2=1,LNOPT)
    DO J1=1,NIMAGE
       WRITE(7,'(3G20.10)') (POINTS((J1-1)*(LNOPT+IADD)+J2),J2=1,LNOPT)
    ENDDO
    WRITE(7,'(3G20.10)') (FINAL(J2),J2=1,LNOPT)
    CLOSE(7)
ENDIF

RETURN
END SUBROUTINE MECCANO
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** July 1990 ***
!
!        Line search removed plus small modifications, DJW 2001
!        If INTMINT is true, then N is NINTS, 
!        otherwise N = NOPT JMC
!        changed declaration of X(N) to X(3*NATOMS) 30/4/04
!        X is passed in and out in Cartesians.
!
!  Cartesian coordinate and gradient vectors are declared 3*NATOMS - the complication
!  is that for internal coordinate optimisations we can have a number of degrees of
!  freedom that is more or less than 3*NATOMS. N should specify this dimension.
!
SUBROUTINE MECBFGS(N,M,X,DIAGCO,EPS,MFLAG,EVEC,RMSVEC,ITMAX,ITDONE,PTEST, &
     &    QSTART,QFINISH)
USE COMMONS
USE KEY
USE MODTWOEND
USE MODUNRES
USE MODCHARMM
USE MODNEB
USE MODMEC
use porfuncs
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL,NCOUNT
DOUBLE PRECISION X(*),G(N),SLENGTH,DDOT,OVERLAP,QSTART(*),QFINISH(*)
DOUBLE PRECISION, ALLOCATABLE :: W(:)
DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,RMS,EREAL,REALRMS,ALPHA,GSAVE(N)
LOGICAL DIAGCO, PTEST
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA
DOUBLE PRECISION DELTAQ(N),EPSILON,EVEC(*),RMSVEC(*)
DOUBLE PRECISION GINT(N),XINT(N) ! JMC
DOUBLE PRECISION GLAST(N)
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN
LOGICAL MFLAG
DOUBLE PRECISION DOT1, DOT2
INTEGER FRAME
DOUBLE PRECISION :: EINITIAL, EFINAL, SEPARATION
COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
LOGICAL PATHT, DRAGT
INTEGER NPATHFRAME, NDECREASE, ISTAT
COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
LOGICAL KNOWE, KNOWG, KNOWH
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80
COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING

ALLOCATE(W(N*(2*M+1)+2*M))

ALPHA=1.0D0
NFAIL=0
FRAME=1
ITER=0
ITDONE=0

FIXIMAGE=.FALSE.

CALL MAKEGRADMEC(X,GSAVE,ENERGY,EVEC,RMSVEC,QSTART,QFINISH)
G(1:N)=GSAVE(1:N)
GLAST(1:N)=GSAVE(1:N)

! DO J1=1,N
   ! X(J1)=X(J1)+1.0D-3
   ! CALL MAKEGRADMEC(X,G,EPLUS,EVEC,RMSVEC,QSTART,QFINISH)
   ! X(J1)=X(J1)-2.0D-3
   ! CALL MAKEGRADMEC(X,G,EMINUS,EVEC,RMSVEC,QSTART,QFINISH)
   ! X(J1)=X(J1)+1.0D-3
   ! PRINT '(A,I5,3G20.10)', 'J1,anal,num,%error=', &
       ! &   J1,GLAST(J1),(EPLUS-EMINUS)/(2.0D-3),ABS(100*(GLAST(J1)-(EPLUS-EMINUS)/(2.0D-3))/GLAST(J1))
! ENDDO
! STOP

EREAL=ENERGY
REALRMS=RMS

IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') &
     &             ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
      WRITE(ESTRING,16) ' Energy for last cycle=',ENERGY,' '
16 FORMAT(A,27X,F20.10,A)

10 CALL FLUSH(6)

MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
!  PRINT*,'RMS,EPS,ITDONE,NSTEPMIN=',RMS,EPS,ITDONE,NSTEPMIN
   IF (ITDONE.LT.NSTEPMIN) MFLAG=.FALSE.
   IF (MFLAG) THEN
      FIXIMAGE=.FALSE.
      DEALLOCATE(W)
      RETURN
   ENDIF
ENDIF

IF (ITDONE.EQ.ITMAX) THEN
   FIXIMAGE=.FALSE.
!  WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',MECDGUESS
   DEALLOCATE(W)
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240  FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   IF (DIAGCO) THEN
      PRINT*,'using estimate of the inverse diagonal elements'
      DO I=1,N
         IF (MECDGUESS.LE.0.0D0) THEN
            WRITE(*,235) I
235         FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/, &
     &                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
            STOP
         ENDIF
      ENDDO
   ENDIF
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR step for diagonal inverse Hessian
!
   DO I=1,N
      W(ISPT+I)= -G(I)*MECDGUESS
      W(I)= -G(I)*MECDGUESS
   ENDDO
!  PRINT*,'I,W,W=',I,W(1),W(ISPT+1)
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  Make the first guess for the step length cautious.
!
   STP=MIN(1.0D0/GNORM,GNORM)
!  STP=1.0D0
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
!  PRINT*,'before overlap W, W: ITER,M,ISPT,IYPT,NPT=',ITER,M,ISPT,IYPT,NPT
!  WRITE(*,'(I5,2E20.10)') (J1,W(ISPT+NPT+J1),W(IYPT+NPT+J1),J1=1,10)
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
!  WRITE(*,'(A,E20.10)') 'YS=',YS
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  Update estimate of diagonal inverse Hessian elements
!  We divide by both YS and YY at different points, so
!  they had better not be zero!
!
   IF (.NOT.DIAGCO) THEN
      YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
!     WRITE(*,'(A,E20.10)') 'YY=',YY
      IF (YY.EQ.0.0D0) YY=1.0D0
      DUMMY1=YS/YY
!     DUMMY1=ABS(YS/YY)
!     WRITE(*,'(A,E20.10)') 'DUMMY1=',DUMMY1
      MECDGUESS=DUMMY1
   ELSE
      PRINT*,'using estimate of the inverse diagonal elements'
      DO I=1,N
         IF (MECDGUESS.LE.0.0D0) THEN
            WRITE(*,235) I
            STOP
         ENDIF
      ENDDO
   ENDIF
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
!  PRINT*,'W(I) gets set to -G(I):'
!  WRITE(*,'(I5,2E20.10)') (J1,W(J1),G(J1),J1=1,10)
   IF (CHRMMT.AND.INTMINT) THEN
      DO I=1,N
         W(I)= -GINT(I)
      ENDDO
   ELSE
      DO I=1,N
         W(I)= -G(I)
      ENDDO
   ENDIF
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
      W(I)=MECDGUESS*W(I)
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
!
!  Store the new search direction
!
!     PRINT*,'W(I):'
!     WRITE(*,'(I5,E20.10)') (J1,W(J1),J1=1,10)
IF (ITER.GT.0) THEN
   DO I=1,N
      W(ISPT+POINT*N+I)= W(I)
   ENDDO
ENDIF

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))

OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) THEN
   IF (CHRMMT.AND.INTMINT) THEN
      OVERLAP=DDOT(N,GINT,1,W,1)/(DOT1*DOT2)
   ELSE
      OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)
  ENDIF
ENDIF
!     PRINT*,'OVERLAP,MECDGUESS=',OVERLAP,MECDGUESS
!     PRINT*,'G . G=',DDOT(N,G,1,G,1)
!     PRINT*,'W . W=',DDOT(N,W,1,W,1)
IF (OVERLAP.GT.0.0D0) THEN
   IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reversing step'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
!        ITER=0
!        GOTO 10
ENDIF

DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MECSTEP) STP=MECSTEP/SLENGTH
!
!  We now have the proposed step.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
!  PRINT '(A,I6,2G20.10)','J1,X,W,STP=',J1,X(J1),W(ISPT+POINT*N+J1)
ENDDO 
KNOWE=.FALSE.
KNOWG=.FALSE.
KNOWH=.FALSE.
!
! At this point we have new Cartesian or internal coordinates after taking a full
! or decreased step. The gradient is not known at this geometry.
! If INTMIN must transform to Cartesians here.
!
NDECREASE=0
EPSILON=1.0D-6
NCOUNT=0
20    CONTINUE

CALL MAKEGRADMEC(X,GSAVE,ENEW,EVEC,RMSVEC,QSTART,QFINISH)

G(1:N)=GSAVE(1:N)

EREAL=ENEW
REALRMS=RMS
!
!  Must allow the energy to rise during a minimisation to allow for numerical noise or
!  systematic errors due to discontinuities or SCF convergence problems.
!
IF (ENEW-ENERGY.LE.MAXERISE) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,G13.5)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE, &
     &     ' LBFGS steps, step:',STP*SLENGTH
   WRITE(ESTRING,16) ' Energy for last cycle=',ENERGY,' '
!
!  Step finished so can reset OLDQ to new XINT, 
!  as well as the Cartesian and internal gradients.
!
   GLAST(1:N)=GSAVE(1:N)
!  IF (.NOT.BFGSTST) CALL DUMPP(X,ENERGY)
ELSE 
!
!  Energy increased - try again with a smaller step size. Must cater for possible enormous
!  values of SLENGTH. Decreasing the step size doesn;t seem to help for CASTEP.
!
   IF (((ITER.GT.1).AND.(NDECREASE.GT.2)).OR.((ITER.LE.1).AND.(NDECREASE.GT.10)).OR. &
     &        ((CASTEP.OR.ONETEP.OR.CP2K).AND.(NDECREASE.GT.1))) THEN 
      NFAIL=NFAIL+1
      IF (PTEST) PRINT*,' in mecbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
!
!  try resetting - go back to previous coordinates, ENERGY is not set to ENEW
!  we need to save the gradient corresponding to the last successful step
!              
      ITER=0
      DO J1=1,N
         X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
         G(J1)=GLAST(J1)
!        G(J1)=GSAVE(J1) ! DJW 6/5/04
      ENDDO
      IF (NFAIL.GT.20) THEN
         PRINT*,' Too many failures - give up'
! check numerical derivatives
!        DO J1=1,N
!           X(J1)=X(J1)+1.0D-5
!           CALL MAKEGRADMEC(X,G,EPLUS,EVEC,RMSVEC,QSTART,QFINISH)
!           X(J1)=X(J1)-2.0D-5
!           CALL MAKEGRADMEC(X,G,EMINUS,EVEC,RMSVEC,QSTART,QFINISH)
!           X(J1)=X(J1)+1.0D-5
!           IF (ABS(100*(GSAVE(J1)-(EPLUS-EMINUS)/(2.0D-5))/GSAVE(J1)).GT.0.01D0) THEN
!              PRINT '(A,I5,4G20.10)', 'J1,X,anal,num,%error=', &
!         &      J1,X(J1),GSAVE(J1),(EPLUS-EMINUS)/(2.0D-5),ABS(100*(GLAST(J1)-(EPLUS-EMINUS)/(2.0D-5))/GLAST(J1))
!              PRINT '(A,2G20.10)', 'EPLUS,EMINUS=',EPLUS,EMINUS
!           ENDIF
!        ENDDO
!        STOP

         FIXIMAGE=.FALSE.
         DEALLOCATE(W)
         RETURN
      ENDIF
      GOTO 30
   ENDIF
!
!  Try a smaller step.
!
   IF (CHRMMT.AND.INTMINT) THEN
      DO J1=1,N
         XINT(J1)=XINT(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.1D0
      ENDDO 
   ELSE
      DO J1=1,N
         X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
      ENDDO 
   ENDIF
   KNOWE=.FALSE.
   KNOWG=.FALSE.
   KNOWH=.FALSE.
   STP=STP/10.0D0
   NDECREASE=NDECREASE+1
   IF (PTEST) &
     &    WRITE(*,'(A,G19.10,A,G16.10,A,G15.8)') ' energy increased from ',ENERGY,' to ',ENEW, &
     &            ' decreasing step to ',STP*SLENGTH
   FIXIMAGE=.TRUE.
   GOTO 20
ENDIF
!
!     Compute the new step and gradient change. Note that the step
!     length is accounted for when the step taken is saved.
!
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
END SUBROUTINE MECBFGS
!
!****************************************************************************************
!
SUBROUTINE MAKEGRADMEC(X,G,SUM,EVEC,RMSVEC,QSTART,FINAL)

USE COMMONS
USE MODMEC
USE MODNEB
USE KEY
IMPLICIT NONE
INTEGER J1,J2,J3
DOUBLE PRECISION DUMMY1,EVEC(*),ETOTAL,FINAL(*), &
     &                 XX(3*NATOMS),GG(3*NATOMS),RMS,G(*),X(*), &
     &                 QSTART(*),RMSVEC(*),SUM,EINITIAL,EFINAL,&
     &                 SEPARATION, LARG(0:NIMAGE+1,0:NIMAGE+1), LDIST, NORM1, NORM2
COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
INTEGER IADD, LNOPT, LNATOMS, LNP
COMMON /INTS/ IADD, LNOPT, LNATOMS

LNP=6
LARG(0:NIMAGE+1,0:NIMAGE+1)=0.0D0
NORM1=1.0D0/NIMAGE
NORM2=2.0D0/((NIMAGE+1)*(NIMAGE+2))
LDIST=X((LNOPT+IADD)*NIMAGE+1)
ETOTAL=0.0D0
DO J1=1,NIMAGE
   DO J2=1,LNOPT+IADD
      XX(J2)=X(J2+(LNOPT+IADD)*(J1-1))
   ENDDO
   CALL POTENTIAL(XX,EVEC(J1),GG,.TRUE.,.FALSE.,RMSVEC(J1),.FALSE.,.FALSE.)
!
!  The points can be changed for C60 in the potential. See oldneb routine.
!
!  PRINT'(A,I5,G20.10)','J1,E=',J1,EVEC(J1)
   DO J2=1,LNOPT+IADD
      G(J2+(LNOPT+IADD)*(J1-1))=NORM1*GG(J2)
   ENDDO
   ETOTAL=ETOTAL+EVEC(J1)
ENDDO

SUM=ETOTAL*NORM1

DO J2=1,NIMAGE
   DUMMY1=0.0D0
   DO J1=1,LNOPT+IADD
      DUMMY1=DUMMY1+((X(J1+(LNOPT+IADD)*(J2-1))-QSTART(J1))/J2)**2              !  (0,J2)
   ENDDO
   LARG(J2,0)=DUMMY1-LDIST**2
   SUM=SUM+NORM2*MECLAMBDA*LARG(J2,0)**2/REAL(J2)**LNP
ENDDO

DO J2=1,NIMAGE
   DO J3=J2+1,NIMAGE
      DUMMY1=0.0D0
      DO J1=1,LNOPT+IADD
         DUMMY1=DUMMY1+((X(J1+(LNOPT+IADD)*(J2-1))-X(J1+(LNOPT+IADD)*(J3-1)))/(J2-J3))**2 ! (J2,J3)
      ENDDO
      LARG(J3,J2)=DUMMY1-LDIST**2
      SUM=SUM+NORM2*MECLAMBDA*LARG(J3,J2)**2/REAL(J2-J3)**LNP
   ENDDO
ENDDO

DO J2=1,NIMAGE
   DUMMY1=0.0D0
   DO J1=1,LNOPT+IADD
      DUMMY1=DUMMY1+((X(J1+(LNOPT+IADD)*(J2-1))-FINAL(J1))/(J2-NIMAGE-1))**2 !  (J2,N+1)
   ENDDO
   LARG(NIMAGE+1,J2)=DUMMY1-LDIST**2
   SUM=SUM+NORM2*MECLAMBDA*LARG(NIMAGE+1,J2)**2/REAL(NIMAGE+1-J2)**LNP
ENDDO

RMS=0.0D0
DO J1=1,LNOPT+IADD ! k
   DO J2=1,NIMAGE ! gamma
      DUMMY1=LARG(NIMAGE+1,J2)*(X(J1+(LNOPT+IADD)*(J2-1))-FINAL(J1))/REAL(J2-NIMAGE-1)**(2+LNP)
      DUMMY1=DUMMY1+LARG(J2,0)*(X(J1+(LNOPT+IADD)*(J2-1))-QSTART(J1))/REAL(J2)**(2+LNP)
      DO J3=1,J2-1
         DUMMY1=DUMMY1+LARG(J2,J3)*(X(J1+(LNOPT+IADD)*(J2-1))-X(J1+(LNOPT+IADD)*(J3-1)))/REAL(J2-J3)**(2+LNP)
      ENDDO
      DO J3=J2+1,NIMAGE
         DUMMY1=DUMMY1+LARG(J3,J2)*(X(J1+(LNOPT+IADD)*(J2-1))-X(J1+(LNOPT+IADD)*(J3-1)))/REAL(J2-J3)**(2+LNP)
      ENDDO
      G(J1+(LNOPT+IADD)*(J2-1))=G(J1+(LNOPT+IADD)*(J2-1))+4*NORM2*MECLAMBDA*DUMMY1
      RMS=RMS+G(J1)**2
   ENDDO
ENDDO

G((LNOPT+IADD)*NIMAGE+1)=0.0D0
LARG(NIMAGE+1,0)=0.0D0
DO J1=0,NIMAGE
   DO J2=J1+1,NIMAGE+1
      G((LNOPT+IADD)*NIMAGE+1)=G((LNOPT+IADD)*NIMAGE+1)+LARG(J2,J1)/REAL(J1-J2)**LNP
   ENDDO
ENDDO
G((LNOPT+IADD)*NIMAGE+1)=-4*NORM2*MECLAMBDA*LDIST*G((LNOPT+IADD)*NIMAGE+1)
RMS=RMS+G((LNOPT+IADD)*NIMAGE+1)**2
RMS=SQRT(RMS/((LNOPT+IADD)*NIMAGE+1))

RETURN
END SUBROUTINE MAKEGRADMEC
