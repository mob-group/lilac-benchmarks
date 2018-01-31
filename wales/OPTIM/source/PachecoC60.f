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
C       POTENTIAL ENERGY OF A (C60)N SYSTEM
C       PACHECO-PRATES-RAMALHO POTENTIAL (PRL 1997)
C
C      X,Y,Z: vectors

      SUBROUTINE PRC60(NATOMS,LCOORDS,V,EPPR,GTEST,SECT)
      USE MODHESS
      IMPLICIT NONE
        DOUBLE PRECISION CAT
C       PARAMETER(CAT=4752000.D0)
        PARAMETER(CAT=0.0D0)
        INTEGER J1, J2, I, J, NATOMS
      DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),LCOORDS(3*NATOMS),V3B,V2B,EPPR,V(3*NATOMS)
      DOUBLE PRECISION MIJ,FIJ,WIJ,XIJ,YIJ,ZIJ,RIJ,RIJ2,VIJ,PMORSE,WTOT,RIJ23,
     1                   fx(NATOMS),fy(NATOMS),fZ(NATOMS),ermi,dwij,dfij,dmij
        DOUBLE PRECISION DMU, DELTA
      parameter(dmu=10.05D0)
      parameter(delta=1.04D0)
        DOUBLE PRECISION DM0, TAU, D0
      parameter(dM0=0.3D0)
      parameter(tau=9.75D0)
      parameter(D0=10.3D0)
        DOUBLE PRECISION C6, C8, C10, C12
      parameter(C6=75600.D0)
      parameter(C8=9122400.D0)
      parameter(C10=2.09D8)
      parameter(C12=7.78D10)
        LOGICAL GTEST,SECT

        DO J1=1,NATOMS
           J2=3*(J1-1)
           X(J1)=LCOORDS(J2+1)
           Y(J1)=LCOORDS(J2+2)
           Z(J1)=LCOORDS(J2+3)
        ENDDO

      V3B=0.D0
      V2B=0.D0

        IF (.NOT.GTEST) THEN
      do i=1,NATOMS
         do j=i+1,NATOMS
C
C       2-body interaction
C
            xij=x(j)-x(i)
            yij=y(j)-y(i)
            zij=z(j)-z(i)
            rij2=xij*xij+yij*yij+zij*zij
            rij=sqrt(rij2)
            
            fij=1.0D0/(1.D0+exp((rij-dmu)/delta))
              pmorse=exp(tau*(1.D0-rij/d0))
            mij=dM0*pmorse*(pmorse-2.D0)
            wij=-(C6+(C8+(C10+C12/rij2)/rij2)/rij2)/rij2**3
            
            vij=fij*mij+(1.D0-fij)*wij
            
            v2b=v2b+vij
         enddo
      enddo

        ELSE

      do i=1,NATOMS
           fx(i)=0.D0
           fy(i)=0.D0
           fz(i)=0.D0
        enddo
      do i=1,NATOMS
         do j=i+1,NATOMS
C
C       2-body interaction
C
            xij=x(j)-x(i)
            yij=y(j)-y(i)
            zij=z(j)-z(i)
            rij2=xij*xij+yij*yij+zij*zij
            rij=sqrt(rij2)
            
              ermi=exp((rij-dmu)/delta)
              fij=1.0D0/(1.D0+ermi)
              dfij=-ermi/(delta*(1.D0+ermi)**2)

              pmorse=exp(tau*(1.D0-rij/d0))
            mij=dM0*pmorse*(pmorse-2.D0)
              dmij=(2.D0*tau*dM0*pmorse*(1.D0-pmorse))/d0

              rij23=rij2**3
            wij=-(C6+(C8+(C10+C12/rij2)/rij2)/rij2)/rij23
              dwij=(6*C6+(8*C8+(10*C10+12*C12/rij2)/rij2)/rij2)/(rij*rij23)
            
            vij=fij*mij+(1.D0-fij)*wij
            
            v2b=v2b+vij

              wtot=mij*dfij+fij*dmij+(1.D0-fij)*dwij-dfij*wij
              wtot=-wtot

              fx(i)=fx(i)+(x(i)-x(j))*wtot/rij
              fy(i)=fy(i)+(y(i)-y(j))*wtot/rij
              fz(i)=fz(i)+(z(i)-z(j))*wtot/rij

              fx(j)=fx(j)+(x(j)-x(i))*wtot/rij
              fy(j)=fy(j)+(y(j)-y(i))*wtot/rij
              fz(j)=fz(j)+(z(j)-z(i))*wtot/rij
       
         enddo
      enddo
        ENDIF
      
      EPPR=v2b

        IF (GTEST) THEN
           DO J1=1,NATOMS
              J2=3*(J1-1)
              V(J2+1)=-FX(J1)
              V(J2+2)=-FY(J1)
              V(J2+3)=-FZ(J1)
C             WRITE(*,'(A,I3,3F20.10)') 'J1,FX,FY,FZ=',J1,FX(J1),FY(J1),FZ(J1)
           ENDDO
        ENDIF

        IF (SECT) CALL HESSIAN(x,y,z,NATOMS)
      
      return
      end

      subroutine HESSIAN(x,y,z,Nsize)
      USE MODHESS
      IMPLICIT NONE
        INTEGER nsize,i,j
      DOUBLE PRECISION X(NSIZE),Y(NSIZE),Z(NSIZE),WYZ,WXZ,WXY,WZZ,WYY,WXX,D2V,D2V2,DV2,DVDW,VDW
      DOUBLE PRECISION MIJ,DMIJ,FIJ,DFIJ,WIJ,DWIJ,XMORSE,XX,YY,ZZ,XY,XZ,YZ,XIJ,YIJ,ZIJ,RIJ2,RIJ,MYFERMI,DFERMI,DMORSE
      DOUBLE PRECISION D2MIJ,D2FIJ,D2WIJ,D2MORSE,D2FERMI,D2VDW

      do i=1,Nsize
         xx=0.0D0
         yy=0.0D0
         zz=0.0D0
         xy=0.0D0
         xz=0.0D0
         yz=0.0D0
         do j=1,Nsize
            if (i.ne.j) then
         xij=x(j)-x(i)
         yij=y(j)-y(i)
         zij=z(j)-z(i)
         rij2=xij*xij+yij*yij+zij*zij
         rij=dsqrt(rij2)
         
         fij=myfermi(rij)
         dfij=dfermi(rij)
         d2fij=d2fermi(rij)
         mij=xmorse(rij)
         dmij=dmorse(rij)
         d2mij=d2morse(rij)
         wij=vdw(rij)
         dwij=dvdw(rij)
         d2wij=d2vdw(rij)
         
         dv2=mij*dfij+fij*dmij+(1.D0-fij)*dwij-dfij*wij
         d2v2=d2fij*mij+2.d0*dfij*dmij+d2mij*fij
         d2v2=d2v2+(1.d0-fij)*d2wij-d2fij*wij-2.d0*dfij*dwij

         d2v=dv2/rij-d2v2

         wxx=-dv2/rij+xij**2*d2v/rij2
         wyy=-dv2/rij+yij**2*d2v/rij2
         wzz=-dv2/rij+zij**2*d2v/rij2
         wxy=xij*yij*d2v/rij2
         wxz=xij*zij*d2v/rij2
         wyz=yij*zij*d2v/rij2

         HESS(3*i-2,3*j-2)=wxx
         HESS(3*i-1,3*j-1)=wyy
         HESS(3*i,3*j)    =wzz
         HESS(3*i-2,3*j-1)=wxy
         HESS(3*i-2,3*j)=wxz
         HESS(3*i-1,3*j)=wyz
         HESS(3*i-1,3*j-2)=wxy
         HESS(3*i,3*j-2)=wxz
         HESS(3*i,3*j-1)=wyz

         xx=xx+wxx
         yy=yy+wyy
         zz=zz+wzz
         xy=xy+wxy
         xz=xz+wxz
         yz=yz+wyz

            endif

         enddo

         HESS(3*i-2,3*i-2)=-xx
         HESS(3*i-1,3*i-1)=-yy
         HESS(3*i,3*i)=-zz
         HESS(3*i-2,3*i-1)=-xy
         HESS(3*i-2,3*i)=-xz
         HESS(3*i-1,3*i)=-yz
         HESS(3*i-1,3*i-2)=-xy
         HESS(3*i,3*i-2)=-xz
         HESS(3*i,3*i-1)=-yz
      
      enddo

      return
      end

c__________________________________________________________________________

      function myfermi(x)
      IMPLICIT NONE
      DOUBLE PRECISION DELTA, X, DMU, MYFERMI
      parameter(dmu=10.05D0)
      parameter(delta=1.04D0)

      myfermi=1.D0+dexp((x-dmu)/delta)
      myfermi=1.D0/myfermi

      return
      end
c__________________________________________________________________________

      function dfermi(x)
      IMPLICIT NONE
      DOUBLE PRECISION DELTA, X, DMU, ERMI, FERMI, DFERMI
      parameter(dmu=10.05D0)
      parameter(delta=1.04D0)

      ermi=dexp((x-dmu)/delta)
      fermi=1.D0+ermi
      dfermi=-ermi/delta
      dfermi=dfermi/(fermi**2)

      return
      end
c__________________________________________________________________________

      function d2fermi(x)
      IMPLICIT NONE
      DOUBLE PRECISION DELTA, X, DMU, ERMI, FERMI, DFERMI, D2FERMI
      parameter(dmu=10.05D0)
      parameter(delta=1.04D0)

      ermi=dexp((x-dmu)/delta)
      fermi=1.D0+ermi
      dfermi=ermi*(ermi-1.d0)/delta/delta
      d2fermi=dfermi/(fermi**3)

      return
      end
c__________________________________________________________________________

      function vdw(x)
      IMPLICIT NONE
      DOUBLE PRECISION C6,C8,C10,C12,X,VDW
      parameter(C6=75600.D0)
      parameter(C8=9122400.D0)
      parameter(C10=2.09D8)
      parameter(C12=7.78D10)

      vdw=-C6/x**6-C8/x**8-C10/x**10-C12/x**12

      return
      end
c__________________________________________________________________________

      function dvdw(x)
      IMPLICIT NONE
      DOUBLE PRECISION C6,C8,C10,C12,X,DVDW
      parameter(C6=75600.D0)
      parameter(C8=9122400.D0)
      parameter(C10=2.09D8)
      parameter(C12=7.78D10)

      dvdw=6*C6/x**7+8*C8/x**9+10*C10/x**11+12*C12/x**13

      return
      end
c__________________________________________________________________________

      function d2vdw(x)
      IMPLICIT NONE
      DOUBLE PRECISION C6,C8,C10,C12,X,D2VDW
      parameter(C6=75600.D0)
      parameter(C8=9122400.D0)
      parameter(C10=2.09D8)
      parameter(C12=7.78D10)

      d2vdw=-42.d0*C6/x**8-72.d0*C8/x**10
      d2vdw=d2vdw-110.d0*C10/x**12-156.d0*C12/x**14

      return
      end
c__________________________________________________________________________

      function xmorse(x)
      IMPLICIT NONE
      DOUBLE PRECISION DM0,TAU,D0,X,XMORSE
      parameter(dM0=0.3D0)
      parameter(tau=9.75D0)
      parameter(D0=10.3D0)

      xmorse=dexp(tau*(1.D0-x/d0))
      xmorse=dM0*xmorse*(xmorse-2.D0)

      return
      end
c__________________________________________________________________________

      function dmorse(x)
      IMPLICIT NONE
      DOUBLE PRECISION DM0,TAU,D0,MORSE,X,DMORSE
      parameter(dM0=0.3D0)
      parameter(tau=9.75D0)
      parameter(D0=10.3D0)

      morse=dexp(tau*(1.D0-x/d0))
      dmorse=2.D0*tau*dM0*morse*(1.D0-morse)
      dmorse=dmorse/d0

      return
      end
c__________________________________________________________________________

      function d2morse(x)
      IMPLICIT NONE
      DOUBLE PRECISION DM0,TAU,D0,X,D2MORSE
      parameter(dM0=0.3D0)
      parameter(tau=9.75D0)
      parameter(D0=10.3D0)
      DOUBLE PRECISION MORSE,DMORSE

      morse=dexp(tau*(1.D0-x/d0))
      dmorse=2.D0*tau*tau*dM0*morse*(2.D0*morse-1.d0)
      d2morse=dmorse/d0/d0

      return
      end
