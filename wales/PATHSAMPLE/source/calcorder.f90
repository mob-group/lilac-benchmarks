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
!
!  Calculate distance for lowest M minima to lowest LJ75 D5h minimum number 1
!  and incomplete icosahedron, minimum 149. Also calculate bond order parameters.
!
SUBROUTINE CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
USE COMMONS, ONLY : PFMIN, PFMEAN, BOXLX, BOXLY, BOXLZ, BULKT, TWOD, RIGIDBODY, LOCATIONA, PATOM1, &
                   PATOM2, EMIN, TEMPERATURE, FVIBMIN, HORDERMIN, OSTART, OFINISH
IMPLICIT NONE
LOGICAL DEBUG
INTEGER NATOMS, J1, J2, J3, NDUMMY, NMIN, NTS, UMIN, UTS, LUNIT, GETUNIT, NBINS, MDECA, MICOS
DOUBLE PRECISION POINTSDECA(3*NATOMS), LOCALPOINTS(3*NATOMS), PFNORM, PEQ(NMIN), REFMIN(3*NATOMS), ETOE(NMIN), GDIST(NMIN)
DOUBLE PRECISION DISTANCE, RMAT(3,3), ETOEMAX, ETOEMIN, GDISTMAX, GDISTMIN, ETOEMEAN, &
  &              ETOESIG, GDISTMEAN, GDISTSIG, DUMMY, DIST2, GMIN(3*NATOMS), POINTSICOS(3*NATOMS), ETHRESH, &
  &              DISTDECA, DISTICOS
REAL(8) Q4, Q6

MDECA=13730373   ! positions of lowest funnel minima
MICOS=17695
! ETHRESH=-394.3592130906 ! lowest 1000 from merge_with_highT database
ETHRESH=1.0D0 ! lowest 1000 from merge_with_highT database
READ(UMIN,REC=MDECA) (POINTSDECA(J2),J2=1,3*NATOMS)
READ(UMIN,REC=MICOS) (POINTSICOS(J2),J2=1,3*NATOMS)
PRINT '(A,I8,A,I8,A,G20.10)','calcorder> reference minima are ',MDECA,' and ',MICOS,' energy threshold=',ETHRESH
PRINT '(A)','     min         E                 dist deca           dist icos             Q4                  Q6'
IF (OSTART.LT.1) OSTART=1
IF (OFINISH.LT.1) OFINISH=1
IF (OSTART.GT.NMIN) OSTART=NMIN
IF (OFINISH.GT.NMIN) OFINISH=NMIN
DO J3=OSTART,OFINISH
   IF (EMIN(J3).LE.ETHRESH) THEN
      READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      CALL MINPERMDIST(POINTSDECA,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTDECA,DIST2,RIGIDBODY, &
     &                          RMAT,.FALSE.)
      CALL MINPERMDIST(POINTSICOS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTICOS,DIST2,RIGIDBODY, &
     &                          RMAT,.FALSE.)
      CALL QORDER_LJ(LOCALPOINTS,Q4,Q6)
      PRINT '(I8,5G20.10)',J3,EMIN(J3),DISTDECA,DISTICOS,Q4,Q6
!     PRINT '(I8,5G20.10)',J3,EMIN(J3),Q4,Q6
   ENDIF
ENDDO

END SUBROUTINE CALCORDER

!
!  Calculate distance between PATOM1 and PATOM2 and minimum distance to
!  first A minimum (assumed global minimum). Also calculate mean and
!  standard deviation for these distances and probability distributions.
!
SUBROUTINE CALCORDEROLD(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
USE COMMONS, ONLY : PFMIN, PFMEAN, BOXLX, BOXLY, BOXLZ, BULKT, TWOD, RIGIDBODY, LOCATIONA, PATOM1, &
                   PATOM2, EMIN, TEMPERATURE, FVIBMIN, HORDERMIN
IMPLICIT NONE
LOGICAL DEBUG
INTEGER NATOMS, J1, J2, J3, NDUMMY, NMIN, NTS, UMIN, UTS, LUNIT, GETUNIT, NBINS
DOUBLE PRECISION LOCALPOINTS(3*NATOMS), PFNORM, PEQ(NMIN), REFMIN(3*NATOMS), ETOE(NMIN), GDIST(NMIN)
DOUBLE PRECISION DISTANCE, RMAT(3,3), ETOEMAX, ETOEMIN, GDISTMAX, GDISTMIN, ETOEMEAN, &
  &              ETOESIG, GDISTMEAN, GDISTSIG, DUMMY, DIST2, GMIN(3*NATOMS)
DOUBLE PRECISION, ALLOCATABLE :: BINWEIGHT(:)

NBINS=50
ALLOCATE(BINWEIGHT(NBINS))

PFNORM=0.0D0
PFMEAN=0.0D0
DO J1=1,NMIN
   PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
   PFNORM=PFNORM+EXP(PFMIN(J1)-PFMEAN)
ENDDO
PFNORM=LOG(PFNORM)+PFMEAN
DUMMY=0.0D0
DO J1=1,NMIN
   PEQ(J1)=EXP(PFMIN(J1)-PFNORM)
   DUMMY=DUMMY+PEQ(J1)
ENDDO
PRINT '(A,G20.10)','sum of equilibrium occupation probabilities=',DUMMY
   
READ(UMIN,REC=LOCATIONA(1)) (GMIN(J2),J2=1,3*NATOMS)
PRINT '(A,I6)','calcorder> reference minimum geometry is number ',LOCATIONA(1)
DO J3=1,NMIN
   READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,3*NATOMS)
   ETOE(J3)=SQRT((LOCALPOINTS(3*(PATOM1-1)+1)-LOCALPOINTS(3*(PATOM2-1)+1))**2 &
  &             +(LOCALPOINTS(3*(PATOM1-1)+2)-LOCALPOINTS(3*(PATOM2-1)+2))**2 &
  &             +(LOCALPOINTS(3*(PATOM1-1)+3)-LOCALPOINTS(3*(PATOM2-1)+3))**2)
   CALL MINPERMDIST(GMIN,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
     &                          RMAT,.FALSE.)
   GDIST(J3)=DISTANCE
ENDDO

ETOEMAX=-1.0D100
ETOEMIN=1.0D100
GDISTMAX=-1.0D100
GDISTMIN=1.0D100
ETOEMEAN=0.0D0
ETOESIG=0.0D0
GDISTMEAN=0.0D0
GDISTSIG=0.0D0

PRINT '(A)','  min         Peq                 EtoE                Gdist'
DO J1=1,NMIN
   ETOEMEAN=ETOEMEAN+PEQ(J1)*ETOE(J1)
   ETOESIG=ETOESIG+PEQ(J1)*ETOE(J1)**2
   IF (ETOE(J1).GT.ETOEMAX) ETOEMAX=ETOE(J1)
   IF (ETOE(J1).LT.ETOEMIN) ETOEMIN=ETOE(J1)

   GDISTMEAN=GDISTMEAN+PEQ(J1)*GDIST(J1)
   GDISTSIG=GDISTSIG+PEQ(J1)*GDIST(J1)**2
   IF (GDIST(J1).GT.GDISTMAX) GDISTMAX=GDIST(J1)
   IF (GDIST(J1).LT.GDISTMIN) GDISTMIN=GDIST(J1)
   PRINT '(I8,4G20.10)',J1,PEQ(J1),ETOE(J1),GDIST(J1)
ENDDO
ETOESIG=DSQRT(ETOESIG-ETOEMEAN**2)
GDISTSIG=DSQRT(GDISTSIG-GDISTMEAN**2)
WRITE(*,'(A,4F20.10)') 'mu/sig of EtoE and Gdist  are ',ETOEMEAN,ETOESIG,GDISTMEAN,GDISTSIG

!
! Probability distribution for end-to-end distance, including Boltzmann weight.
!
PRINT '(A)','Probability distribution for end-to-end distance including Boltzmann weight'
BINWEIGHT(1:NBINS)=0.0D0
DO J1=1,NMIN
   IF (ETOE(J1).GE.ETOEMAX) THEN
      J3=NBINS
   ELSE
      J3=INT(NBINS*(ETOE(J1)-ETOEMIN)/(ETOEMAX-ETOEMIN))+1
   ENDIF
   BINWEIGHT(J3)=BINWEIGHT(J3)+EXP(PFMIN(J1)-PFNORM)
ENDDO
DUMMY=0.0D0
DO J2=1,NBINS
   DUMMY=DUMMY+BINWEIGHT(J2)
   PRINT '(2G20.10)',ETOEMIN+(J2-0.5D0)*(ETOEMAX-ETOEMIN)/NBINS,BINWEIGHT(J2)
ENDDO
PRINT '(A,G20.10)',' sum of bin weights=',DUMMY
!
! Probability distribution for distance to global minimum, including Boltzmann weight.
!
PRINT '(A)','Probability distribution for distance from global minimum including Boltzmann weight'
BINWEIGHT(1:NBINS)=0.0D0
DO J1=1,NMIN
   IF (GDIST(J1).GE.GDISTMAX) THEN
      J3=NBINS
   ELSE
      J3=INT(NBINS*(GDIST(J1)-GDISTMIN)/(GDISTMAX-GDISTMIN))+1
   ENDIF
   BINWEIGHT(J3)=BINWEIGHT(J3)+EXP(PFMIN(J1)-PFNORM)
ENDDO
DUMMY=0.0D0
DO J2=1,NBINS
   DUMMY=DUMMY+BINWEIGHT(J2)
   PRINT '(2G20.10)',GDISTMIN+(J2-0.5D0)*(GDISTMAX-GDISTMIN)/NBINS,BINWEIGHT(J2)
ENDDO
PRINT '(A,G20.10)',' sum of bin weights=',DUMMY


END SUBROUTINE CALCORDEROLD

!  
! Q4 and Q6 routines from GMIN.
! 
!---======================================---
      SUBROUTINE QORDER_LJ(Q,Q4,Q6)
      USE COMMONS, ONLY: NATOMS

      IMPLICIT NONE

      REAL(8) Q(3,NATOMS), Q4, Q6
        INTEGER J, K, NB, I,m
        DOUBLE PRECISION DISTBOND
        PARAMETER (DISTBOND=1.3909)
        DOUBLE PRECISION DX,DY,DZ,DIST,phi,costheta,coef2,arg,pi
        COMPLEX Y2,Q4bar(0:4),Q6bar(0:6)


	NB=0
	pi=dacos(-1d0)
	do m=0,4
	   Q4bar(m)=(0,0d0)
	enddo
	do m=0,6
	   Q6bar(m)=(0,0d0)
	enddo
	do J=1,NATOMS-1
	   do K=J+1,NATOMS
	      DX=Q(1,J)-Q(1,K)
	      DY=Q(2,J)-Q(2,K)
	      DZ=Q(3,J)-Q(3,K)

	      DIST=DSQRT(DX*DX+DY*DY+DZ*DZ)
	      IF (DIST.lt.DISTBOND) THEN
		 NB=NB+1
		 costheta =DZ/DIST
		 phi = datan(DY/DX)
		 if(DX.lt.0) phi=phi+pi
		 do m=0,4
		    Q4bar(m)=Q4bar(m)+Y2(4,m,costheta,phi)
		 enddo
		 do m=0,6
		    Q6bar(m)=Q6bar(m)+Y2(6,m,costheta,phi)
		 enddo
	      ENDIF
	   enddo
	enddo
	Q4=0
	do m=0,4
	   Q4=Q4+coef2(4,m)*Q4bar(m)*conjg(Q4bar(m))
	enddo
	Q4=dsqrt(Q4)/(3*Nb)
	Q6=0
	do m=0,6
	   Q6=Q6+coef2(6,m)*Q6bar(m)*conjg(Q6bar(m))
	enddo
	Q6=dsqrt(Q6/13.D0)/Nb
	return
	end
	
	DOUBLE PRECISION function coef2(l,m)
	integer l,m,k
	if(m.eq.0) then
	   coef2=2*l+1d0
	else
	   coef2=4*l+2.
	   do k=l-m+1,l+m
	      coef2=coef2/k
	   enddo
	endif
	return
	end
	
	COMPLEX function Y2(l,m,costheta,phi)
!Computes the the spherical harmonic Y2(l.m)*sqrt(4*pi)  
	implicit none
	INTEGER l,m
	DOUBLE PRECISION plgndr2,costheta,phi
	if(m.lt.0) stop 'm<0 in Y2'
	if(m.eq.0) then
	   Y2=plgndr2(l,m,costheta)
	else 
	   Y2=plgndr2(l,m,costheta)*exp(cmplx(0d0,m*phi))
	endif
	return 
	END
	
	DOUBLE PRECISION FUNCTION plgndr2(l,m,x) 
	implicit none
	INTEGER l,m 
	DOUBLE PRECISION x
!Computes the associated Legendre polynomial Pml (x).
	INTEGER i,ll 
	DOUBLE PRECISION fact,pll,pmm,pmmp1,somx2 
	if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) &
        & stop 'bad arguments in plgndr2'
	pmm=1.  
	if(m.gt.0) then 
	   somx2=sqrt((1.-x)*(1.+x)) 
	   fact=1. 
	   do i=1,m 
	      pmm=-pmm*fact*somx2 
	      fact=fact+2. 
	   enddo 
	endif 
	if(l.eq.m) then 
	   plgndr2=pmm 
	else 
	   pmmp1=x*(2*m+1)*pmm 
	   if(l.eq.m+1) then 
	      plgndr2=pmmp1 
	   else 
	      do ll=m+2,l
		 pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m) 
		 pmm=pmmp1 
		 pmmp1=pll 
	      enddo  
	      plgndr2=pll 
	   endif 
	endif 
	return 
	END



