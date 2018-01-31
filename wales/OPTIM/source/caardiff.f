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
C
C*************************************************************************
C
C  Subroutine CAARDIFF calculates the cartesian gradient and second
C  derivative matrix analytically. Atomic units for CaAr_n clusters
C
C*************************************************************************
C
      SUBROUTINE CAARDIFF(N, X, V, ENERGY, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, R6,
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N),
     3                 R14(N,N), F(N,N)
      DOUBLE PRECISION R1(N),R1T
      DOUBLE PRECISION SIGAR6,EPSAR,CAAR0,CAAR1,CAAR2,CAAR3
C 
C  Store distance matrices.
C

      SIGAR6=58808.93841
      EPSAR=0.00045935
      CAAR0=0.85*14.6732
      CAAR1=1.09357
      CAAR2=-0.85*606.872
      CAAR3=0.85*6.4D5

c      SIGAR6=58808.93841
c      EPSAR=0.00045935
c      CAAR0=12.472220
c      CAAR1=1.09357
c      CAAR2=-515.84120
c      CAAR3=5.55475d5

      ENERGY=0.0D0
      IF (GTEST) THEN
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               if (j1.eq.1) then
                  r1(j2)=dsqrt(r2(j2,j1))
               endif
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               if (j1.eq.1) then
                  energy=energy+caar0*dexp(-caar1*r1(j2))
                  energy=energy+caar2*r6+caar3*r6*r6
               else
                  energy=energy+4.d0*epsar*sigar6*r6*(sigar6*r6-1.d0)
               endif
               R8(J2,J1)=R2(J2,J1)*r6
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            J3=3*(J1-1)
            DO J2=J1+1,N
               J4=3*(J2-1)
               r1t=0.d0
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               if (j1.eq.1) r1t=dsqrt(r2t)
               R2T=1.0D0/R2T
               R6=R2T**3
               if (j1.eq.1) then
                  energy=energy+caar0*dexp(-caar1*r1t)
                  energy=energy+caar2*r6+caar3*r6*r6
               else
                  energy=energy+4.d0*epsar*sigar6*r6*(sigar6*r6-1.d0)
               endif
            ENDDO
         ENDDO

      ENDIF

      energy=energy/epsar

      IF (.NOT.GTEST) RETURN
      CALL CAARG(G,r1,R14,R8,V,X,N)
      
      IF (.NOT.STEST) RETURN
      CALL CAARS(G,F,r1,R2,R14,R8,X,N)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE CAARG(G,r1,R14,R8,V,X,N)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 V(3*N), X(3*N), DUMMY
      DOUBLE PRECISION R1(N)
      DOUBLE PRECISION SIGAR6,EPSAR,CAAR0,CAAR1,CAAR2,CAAR3
C
C  Calculate the g tensor.
C
      SIGAR6=58808.93841
      EPSAR=0.00045935
      CAAR0=0.85*14.6732
      CAAR1=1.09357
      CAAR2=-0.85*606.872
      CAAR3=0.85*6.4D5

c      SIGAR6=58808.93841
c      EPSAR=0.00045935
c      CAAR0=12.472220
c      CAAR1=1.09357
c      CAAR2=-515.84120
c      CAAR3=5.55475d5

      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            if (j1.eq.1) then
               g(j2,j1)=-caar1*caar0*dexp(-caar1*r1(j2))/r1(j2)
               g(j2,j1)=g(j2,j1)-6.d0*caar2*r8(j2,j1)-12.d0*caar3*r14(j2,j1)
            else
               g(j2,j1)=-24.d0*sigar6*epsar*(2.d0*sigar6*r14(j2,j1)-r8(j2,j1))
            endif
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      do j1=1,n
         do j2=1,n
            g(j1,j2)=g(j1,j2)/epsar
         enddo
      enddo

C
C  From here on down the code is system-independent!
C  First calculate the gradient analytically.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE CAARS(G,F,r1,R2,R14,R8,X,N)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, I
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), 
     2                 X(3*N)
      DOUBLE PRECISION R1(N)
      DOUBLE PRECISION SIGAR6,EPSAR,CAAR0,CAAR1,CAAR2,CAAR3,SIGAR2
      DOUBLE PRECISION X0,Y0,Z0,XX,YY,ZZ,XY,YZ,XZ,XI,YI,ZI,XI2,YI2,ZI2,RI,RI2
      DOUBLE PRECISION W1,W11,WX1,WY1,WZ1,W2,WX2,WY2,WZ2,W3,WX3,WY3,WZ3
      DOUBLE PRECISION WXX,WYY,WZZ,WXY,WYZ,WXZ,WTOT,RIJ,XIJ,YIJ,ZIJ
      DOUBLE PRECISION W,S,WX,WY,WZ

      sigar2=38.887896
      SIGAR6=58808.93841
      EPSAR=0.00045935
      CAAR0=0.85*14.6732
      CAAR1=1.09357
      CAAR2=-0.85*606.872
      CAAR3=0.85*6.4D5

c      SIGAR6=58808.93841
c      EPSAR=0.00045935
c      CAAR0=12.472220
c      CAAR1=1.09357
c      CAAR2=-515.84120
c      CAAR3=5.55475d5

c     first do the Ar-Ar part

      do j1=1,n*3
         do j2=1,n*3
            HESS(j1,j2)=0.d0
         enddo
      enddo

c     Ar-Ar part
      
      do i=2,n
         xx=0.0
         yy=0.0
         zz=0.0
         xi=x(3*i-2)
         yi=x(3*i-1)
         zi=x(3*i  )
         do j1=2,N
            if (i.ne.j1) then
               xij=xi-x(3*j1-2)
               yij=yi-x(3*j1-1)
               zij=zi-x(3*j1  )
                 rij=xij**2.+yij**2.+zij**2.
                 w=2.*(sigar2/rij)**3-1.
                 s=(7.*(sigar2/rij)**3.-2.)/rij
                 wx=w-4.*s*xij**2.
                 wy=w-4.*s*yij**2.
                 wz=w-4.*s*zij**2.
                 wx=24.*epsar*wx*(sigar2**3)/(rij**4.)
                 wy=24.*epsar*wy*(sigar2**3)/(rij**4.)
                 wz=24.*epsar*wz*(sigar2**3)/(rij**4.)
                 xx=xx-wx
                 yy=yy-wy
                 zz=zz-wz
                 
                 HESS(3*i-2,3*j1-2)=wx
                 HESS(3*i-1,3*j1-1)=wy
                 HESS(3*i,3*j1)=wz
                 
                 HESS(3*j1,3*i)=HESS(3*i,3*j1)
                 HESS(3*j1-1,3*i-1)=HESS(3*i-1,3*j1-1)
                 HESS(3*j1-2,3*i-2)=HESS(3*i-2,3*j1-2)
                 
              endif
           enddo
           HESS(3*i-2,3*i-2)=xx
           HESS(3*i-1,3*i-1)=yy
           HESS(3*i,3*i)=zz
        enddo

        do i=2,N
           xy=0.0
           xz=0.0
           yz=0.0
           xi=x(3*i-2)
           yi=x(3*i-1)
           zi=x(3*i)
           do j1=2,N
              if (i.ne.j1) then
                 xij=xi-x(3*j1-2)
                 yij=yi-x(3*j1-1)
                 zij=zi-x(3*j1)
                 rij=xij**2.+yij**2.+zij**2.
                 s=sigar6*(7.*(sigar2/rij)**3.-2.)/(rij**5.)
                 wz=epsar*96.*s*xij*yij
                 wy=epsar*96.*s*xij*zij
                 wx=epsar*96.*s*yij*zij
                 xy=xy+wz
                 xz=xz+wy
                 yz=yz+wx

                 HESS(3*i-2,3*j1-1)=-wz
                 HESS(3*i-2,3*j1)=-wy
                 HESS(3*i-1,3*j1)=-wx
                 HESS(3*j1,3*i-1)=HESS(3*i-1,3*j1)
                 HESS(3*j1,3*i-2)=HESS(3*i-2,3*j1)
                 HESS(3*j1-1,3*i-2)=HESS(3*i-2,3*j1-1)
                 
                 HESS(3*j1-2,3*i-1)=HESS(3*i-2,3*j1-1)
                 HESS(3*j1-2,3*i)=HESS(3*i-2,3*j1)
                 HESS(3*i,3*j1-2)=HESS(3*j1,3*i-2)
                 HESS(3*i-1,3*j1-2)=HESS(3*j1-1,3*i-2)
                 HESS(3*j1-1,3*i)=HESS(3*i-1,3*j1)
                 HESS(3*i,3*j1-1)=HESS(3*j1,3*i-1)
                 
              endif
           enddo
           HESS(3*i-2,3*i-1)=xy
           HESS(3*i-1,3*i-2)=xy
           HESS(3*i-2,3*i)=xz
           HESS(3*i,3*i-2)=xz
           HESS(3*i-1,3*i)=yz
           HESS(3*i,3*i-1)=yz
        enddo

C     Ca-Ar part

      X0=X(1)
      Y0=X(2)
      Z0=X(3)
      XX=0.D0
      YY=0.D0
      ZZ=0.D0
      XY=0.D0
      XZ=0.D0
      YZ=0.D0

      DO I=2,N
         XI=X0-X(3*I-2)
         YI=Y0-X(3*I-1)
         ZI=Z0-X(3*I  )
         XI2=XI**2
         YI2=YI**2
         ZI2=ZI**2
         RI2=XI2+YI2+ZI2
         XI2=XI2/RI2
         YI2=YI2/RI2
         ZI2=ZI2/RI2
         RI=DSQRT(RI2)

         W1=CAAR0*CAAR1*DEXP(-CAAR1*RI)
         W1=W1/RI
         W11=1.D0+CAAR1*RI
         WX1=W1*(XI2*W11-1.D0)
         WY1=W1*(YI2*W11-1.D0)
         WZ1=W1*(ZI2*W11-1.D0)

         W2=6.D0*CAAR2/RI2**4
         WX2=W2*(1.D0-8.D0*XI2)
         WY2=W2*(1.D0-8.D0*YI2)
         WZ2=W2*(1.D0-8.D0*ZI2)

         W3=12.D0*CAAR3/RI2**7
         WX3=W3*(1.D0-14.D0*XI2)
         WY3=W3*(1.D0-14.D0*YI2)
         WZ3=W3*(1.D0-14.D0*ZI2)

         WXX=WX1-WX2-WX3
         WYY=WY1-WY2-WY3
         WZZ=WZ1-WZ2-WZ3

c         WXX=-WXX
c         WYY=-WYY
c         WZZ=-WZZ

         XX=XX+WXX
         YY=YY+WYY
         ZZ=ZZ+WZZ

         HESS(1,3*I-2)=HESS(1,3*i-2)-WXX
         HESS(3*I-2,1)=HESS(3*I-2,1)-WXX
         HESS(2,3*I-1)=HESS(2,3*I-1)-WYY
         HESS(3*I-1,2)=HESS(3*I-1,2)-WYY
         HESS(3,3*I  )=HESS(3,3*I  )-WZZ
         HESS(3*I  ,3)=HESS(3*I  ,3)-WZZ

         HESS(3*i-2,3*i-2)=HESS(3*i-2,3*i-2)+wxx
         HESS(3*i-1,3*i-1)=HESS(3*i-1,3*i-1)+wyy
         HESS(3*i  ,3*i  )=HESS(3*i  ,3*i  )+wzz

         W1=CAAR0*CAAR1*(CAAR1+1.D0/RI)*DEXP(-CAAR1*RI)/RI2
         W2=48.D0*CAAR2/RI2**5
         W3=168.D0*CAAR3/RI2**8
         WTOT=W1+W2+W3

         WXY=WTOT*XI*YI
         WXZ=WTOT*XI*ZI
         WYZ=WTOT*ZI*YI

         XY=XY+WXY
         YZ=YZ+WYZ
         XZ=XZ+WXZ

         HESS(1,3*I-1)=HESS(1,3*I-1)-WXY
         HESS(3*I-1,1)=HESS(3*I-1,1)-WXY
         HESS(2,3*I-2)=HESS(2,3*I-2)-WXY
         HESS(3*I-2,2)=HESS(3*I-2,2)-WXY

         HESS(1,3*I  )=HESS(1,3*I  )-WXZ
         HESS(3*I  ,1)=HESS(3*I  ,1)-WXZ
         HESS(3,3*I-2)=HESS(3,3*I-2)-WXZ
         HESS(3*I-2,3)=HESS(3*I-2,3)-WXZ

         HESS(2,3*I  )=HESS(2,3*I  )-WYZ
         HESS(3*I  ,2)=HESS(3*I  ,2)-WYZ
         HESS(3,3*I-1)=HESS(3,3*I-1)-WYZ
         HESS(3*I-1,3)=HESS(3*I-1,3)-WYZ

         HESS(3*i-2,3*i-1)=HESS(3*i-2,3*i-1)+wxy
         HESS(3*i-1,3*i-2)=HESS(3*i-1,3*i-2)+wxy
         HESS(3*i-1,3*i  )=HESS(3*i-1,3*i  )+wyz
         HESS(3*i  ,3*i-1)=HESS(3*i  ,3*i-1)+wyz
         HESS(3*i-2,3*i  )=HESS(3*i-2,3*i  )+wxz
         HESS(3*i  ,3*i-2)=HESS(3*i  ,3*i-2)+wxz

      ENDDO

      HESS(1,1)=XX
      HESS(2,2)=YY
      HESS(3,3)=ZZ
      HESS(1,2)=XY
      HESS(2,1)=XY
      HESS(1,3)=XZ
      HESS(3,1)=XZ
      HESS(2,3)=YZ
      HESS(3,2)=YZ

      do j1=1,3*n
         do j2=1,3*n
            HESS(j1,j2)=HESS(j1,j2)/epsar
         enddo
      enddo

      RETURN
      END
