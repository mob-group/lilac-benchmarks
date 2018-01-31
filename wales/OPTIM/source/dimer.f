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
C      subroutine fdimer (r,b1,c1,b2,c2,a,value,grad,icall)
       subroutine fdimer (coords,value,grad,icall)
       USE MODHESS
       USE COMMONS
       IMPLICIT NONE
       INTEGER J1,NPAR,ICALL,NP,I
       DOUBLE PRECISION BPAR(1000,5),GRAD(3*NATOMS),X(5),
     &           bfit(1000,3),coords(3*NATOMS),VALUE,B,BB2,EOFF,R,B1,C1,B2,C2,A
       character(LEN=10) index(1000)
       common /bs/b,bb2,x,eoff,npar
       save /bs/

       b1=coords(2)
       c1=coords(3)
       b2=coords(4)
       c2=coords(5)
       a=coords(6)
       WRITE(*,'(6F15.7)') (COORDS(J1),J1=1,6)

       if ( icall .eq. 1) then
          np=1
          open(unit=40,file='dimer.dat',status='old')
          read (40,'(6f10.5)') eoff,(x(i),i=1,5)
110         read (40,'(a10,5e25.17)',end=210) index(np),
     &                     (bpar(np,i),i=1,4)
            np=np+1
          goto 110
210       continue
          close(40)
          npar=np-1
       endif  
       PRINT*,'ERROR, R is not set'
       STOP
!      do i=1,npar
!         bfit(i,1)=bpar(i,1)+bpar(i,2)/r**6+bpar(i,3)/r**8
!    &             +bpar(i,4)/r**5
!         bfit(i,2)=-6.d0*bpar(i,2)/r**7 - 8.d0*bpar(i,3)/r**9
!    &             -5.d0*bpar(i,4)/r**6
!         bfit(i,3)=42.d0*bpar(i,2)/r**8 + 72.d0*bpar(i,3)/r**10
!    &             +30.d0*bpar(i,4)/r**7
!      end do
       call dimer (bfit,b1,c1,b2,c2,a,index,eoff,npar,value,grad)
       return
       end
C
C   Calculation of P.E. surface
C
       subroutine dimer (b,b1,c1,b2,c2,a,index,eoff,npar,value,grad)
       USE MODHESS
       USE COMMONS
       IMPLICIT NONE
       DOUBLE PRECISION  B(1000,3),GRAD(3*NATOMS), WT(16),
     1                  B1,B2,C2,A,EOFF,VALUE,GRAD2D,GRADD,DVALUE,C1,
     2                  ANG,BDDC,BDDS,BDD1,BDD3,D1BDDC,D1BDDS,D1BDD1,D1BDD3,D2BDDC
       INTEGER NPAR,I,J,L1,K1,L2,K2,M,N,N1,NOTS
       character(LEN=10) index(*),symf(16)
 
 20    format (5i2)
       value=0.
       do i=1,6
          grad(i)=0.d0
          do j=1,6
             hess(i,j)=0.d0
          end do
       end do
       do i=1,npar
          read (index(i),20) l1,k1,l2,k2,m
          do n=1,8
             wt(n)=1.d0
             wt(n+8)=(-1.d0)**(l1+l2)
          end do
          write (symf(1),20) l1,k1,l2,k2,m
          write (symf(2),20) l1,k1,l2,k2,-m
          write (symf(3),20) l1,k1,l2,-k2,m
          write (symf(4),20) l1,k1,l2,-k2,-m
          write (symf(5),20) l1,-k1,l2,k2,m
          write (symf(6),20) l1,-k1,l2,k2,-m
          write (symf(7),20) l1,-k1,l2,-k2,m
          write (symf(8),20) l1,-k1,l2,-k2,-m
          write (symf(9),20) l2,k2,l1,k1,m
          write (symf(10),20) l2,k2,l1,k1,-m
          write (symf(11),20) l2,k2,l1,-k1,m
          write (symf(12),20) l2,k2,l1,-k1,-m
          write (symf(13),20) l2,-k2,l1,k1,m
          write (symf(14),20) l2,-k2,l1,k1,-m
          write (symf(15),20) l2,-k2,l1,-k1,m
          write (symf(16),20) l2,-k2,l1,-k1,-m
          do n = 2,16
             do n1=1,n-1
                if (symf(n).eq.symf(n1)) wt(n)=0.d0
             end do
          end do
          nots=0
          do n=1,16
             if ( wt(n).ne.0.) then
                read (symf(n),20) l1,k1,l2,k2,m

c
c
          ang = k1*c1+k2*c2+m*a
          bddc =wt(n)*b(i,1)*dvalue(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)*
     &          cos(ang)
          bdds =wt(n)*b(i,1)*dvalue(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)*
     &          sin(ang)
          bdd1 =wt(n)*b(i,1)*gradd(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)
          bdd3 =wt(n)*b(i,1)*dvalue(l1,m,k1,b1)*gradd(l2,-m,k2,b2)
          d1bddc =wt(n)*b(i,2)*dvalue(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)*
     &          cos(ang)
          d1bdds =wt(n)*b(i,2)*dvalue(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)*
     &          sin(ang)
          d1bdd1 =wt(n)*b(i,2)*gradd(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)
          d1bdd3 =wt(n)*b(i,2)*dvalue(l1,m,k1,b1)*gradd(l2,-m,k2,b2)
          d2bddc =wt(n)*b(i,3)*dvalue(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)*
     &          cos(ang)
          value=value+bddc
          grad(1)=grad(1)+d1bddc
          grad(2)=grad(2)+bdd1*cos(ang)
          grad(3)=grad(3)-k1*bdds
          grad(4)=grad(4)+bdd3*cos(ang)
          grad(5)=grad(5)-k2*bdds
          grad(6)=grad(6)-m*bdds
          hess(6,6)=hess(6,6)-m*m*bddc
          hess(6,5)=hess(6,5)-m*k2*bddc
          hess(6,4)=hess(6,4)-m*bdd3*sin(ang)
          hess(6,3)=hess(6,3)-m*k1*bddc
          hess(6,2)=hess(6,2)-m*bdd1*sin(ang)
          hess(6,1)=hess(6,1)-m*d1bdds
          hess(5,5)=hess(5,5)-k2*k2*bddc
          hess(5,4)=hess(5,4)-k2*bdd3*sin(ang)
          hess(5,3)=hess(5,3)-k2*k1*bddc
          hess(5,2)=hess(5,2)-k2*bdd1*sin(ang)
          hess(5,1)=hess(5,1)-k2*d1bdds
          hess(4,4)=hess(4,4)+cos(ang)*b(i,1)*wt(n)*
     &              dvalue(l1,m,k1,b1)*grad2d(l2,-m,k2,b2)
          hess(4,3)=hess(4,3)-k1*bdd3*sin(ang)
          hess(4,2)=hess(4,2)+cos(ang)*b(i,1)*wt(n)*
     &              gradd(l1,m,k1,b1)*gradd(l2,-m,k2,b2)
          hess(4,1)=hess(4,1)+d1bdd3*cos(ang)
          hess(3,3)=hess(3,3)-k1*k1*bddc
          hess(3,2)=hess(3,2)-k1*bdd1*sin(ang)
          hess(3,1)=hess(3,1)-k1*d1bdds
          hess(2,2)=hess(2,2)+cos(ang)*b(i,1)*wt(n)*
     &              grad2d(l1,m,k1,b1)*dvalue(l2,-m,k2,b2)
          hess(2,1)=hess(2,1)+d1bdd1*cos(ang)
          hess(1,1)=hess(1,1)+d2bddc
c
c
          end if
         end do
       end do
       do i=6,1,-1
          do j=i,1,-1
             hess(j,i)=hess(i,j)
          end do
       end do
       value=value+eoff
       return
       end

       function grad2d (j,m1,m,t)
       IMPLICIT NONE
       DOUBLE PRECISION FACTLOCAL,X2N,COEFF,VALUE,DENOM,T,T2,GRAD2D
       INTEGER M,MDIFF,N2,N1,I,IA,IB,J,M1

       coeff=sqrt (factlocal(j+m)*factlocal(j-m)*factlocal(j+m1)*factlocal(j-m1))
       mdiff=m1-m
       n2=min(j-m1,j+m)
       t2=t/2.
       value=0.
       if (mdiff.lt.0) then
           n1=-mdiff
       else
           n1=0
       end if
       do i=n1,n2
       ia=2*j-mdiff-2*i
       ib=mdiff+2*i
          denom=factlocal(i)*factlocal(mdiff+i)*factlocal(j+m-i)*factlocal(j-m1-i)*4
          value=value+(-1)**(mdiff+i)*(
     &          ia*(ia-1.)*x2n(cos(t2),ia-2)*x2n(sin(t2),ib+2)+
     &          ib*(ib-1)*x2n(cos(t2),ia+2)*x2n(sin(t2),ib-2)-
     &          2*(ia*ib+j)*x2n(cos(t2),ia)*x2n(sin(t2),ib))/denom
       end do
       grad2d=coeff*value
       return
       end

       function gradd (j,m1,m,t)
       IMPLICIT NONE
       DOUBLE PRECISION FACTLOCAL,X2N,COEFF,VALUE,DENOM,T,T2,GRADD
       INTEGER J,M1,M,MDIFF,N2,N1,I

       coeff=sqrt (factlocal(j+m)*factlocal(j-m)*factlocal(j+m1)*factlocal(j-m1))
       mdiff=m1-m
       n2=min(j-m1,j+m)
       t2=t/2.
       value=0.
       if (mdiff.lt.0) then
           n1=-mdiff
       else
           n1=0
       end if
       do i=n1,n2
          denom=factlocal(i)*factlocal(mdiff+i)*factlocal(j+m-i)*factlocal(j-m1-i)*2
          value=value+(-1)**(mdiff+i)*(-1*
     &          (2*j-mdiff-2*i)*x2n(cos(t2),2*j-mdiff-2*i-1)*
     &          x2n(sin(t2),mdiff+2*i+1)+
     &          (mdiff+2*i)*x2n(cos(t2),2*j-mdiff-2*i+1)*
     &          x2n(sin(t2),mdiff+2*i-1))/denom
       end do
       gradd=coeff*value
       return
       end

       function dvalue (j,m1,m,t)
       IMPLICIT NONE
       DOUBLE PRECISION FACTLOCAL,X2N,COEFF,VALUE,DENOM,T,DVALUE
       INTEGER J, M1, M, MDIFF, N2, N1, I

       coeff=sqrt (factlocal(j+m)*factlocal(j-m)*factlocal(j+m1)*factlocal(j-m1))
       mdiff=m1-m
       n2=min(j-m1,j+m)
       value=0.
       if (mdiff.lt.0) then
           n1=-mdiff
       else
           n1=0
       end if
       do i=n1,n2
          denom=factlocal(i)*factlocal(mdiff+i)*factlocal(j+m-i)*factlocal(j-m1-i)
          value=value+(-1)**(mdiff+i)*x2n(cos(t/2.),(2*j-mdiff-2*i))*
     &          x2n(sin(t/2.),mdiff+2*i)/denom
       end do
       dvalue=coeff*value
       return
       end

       function factlocal(n)
       IMPLICIT NONE
       INTEGER N, I
       DOUBLE PRECISION :: FACTLOCAL

       if (n.eq.0) then
           factlocal=1
       else
           factlocal=1.
           do i=2,n
              factlocal=factlocal*i
           end do
       end if
cc     write (6,*) n,fact
       return
       end

       function x2n(x,n)
       IMPLICIT NONE
       DOUBLE PRECISION X, X2N
       INTEGER N

       if (x.eq.0.) then
          if (n .eq. 0) then
              x2n=1.
          else
              x2n=0.
          end if
        else
          x2n=x**n
        end if
        return
        end




