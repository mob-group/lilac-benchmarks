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
C********************************************************************
C
C Subroutine ENERGY calculates the SC energy:
C
C********************************************************************
C
      SUBROUTINE ACK(N,X,PSC,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF,PRESSURE,GTEST,STEST)
      USE KEY, ONLY : ACKLANDID
      IMPLICIT NONE
      INTEGER N, J1, J2, J, I, MZ, MY, MX 
      LOGICAL PRESSURE, YESNO
      DOUBLE PRECISION X(3*N), POTA, POTB, DIST, CUTOFF, BOXLZ, BOXLY, BOXLX, PSC, VNEW(3*N),Rc
      COMMON /param_cut_off/Rc

      DOUBLE PRECISION RHO(3*N)
      DOUBLE PRECISION VEC(N,N,3),Rneigh(N,N,3),VSITE(N),force(3,N)
      double precision rbuf,r,zero,norm,rho_temp,vpot_temp,Fembed_d_i
      double precision Vpot,Vpot_d,rho_pot,rho_pot_d,Fembed,Fembed_d
      integer ipot,icount,ja
      integer ic(N),neigh_type(N,N),ndir(N)
      LOGICAL GTEST,STEST

      
      Rc=CUTOFF
      ipot=ACKLANDID
      zero=1.0D-12
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is
C  used:
C
      icount=0
      ic(:)=0
      
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
 
         icount=ic(J1)
         
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            MZ=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(MX)
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(MY)
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(MZ)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
            norm=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 + VEC(J1,J2,3)**2
            if (norm < Rc*Rc.and.norm>zero) then
            icount=icount+1
            ic(J2)=ic(J2)+1
            Rneigh(J2,ic(J2),1)=VEC(J1,J2,1)
            Rneigh(J2,ic(J2),2)=VEC(J1,J2,2)
            Rneigh(J2,ic(J2),3)=VEC(J1,J2,3)
            
            Rneigh(J1,icount,1)=VEC(J2,J1,1)
            Rneigh(J1,icount,2)=VEC(J2,J1,2)
            Rneigh(J1,icount,3)=VEC(J2,J1,3)
            
            neigh_type(J1,icount)=J2
            neigh_type(J2,ic(J2))=J1
            
            end if
15       CONTINUE
      ndir(J1)=icount
25    CONTINUE
C
C
C Call scl.f for lattice constant optimisation if required:
C
      IF (PRESSURE) THEN
         write(*,*) 'NO PRESUURE IMPLEMENTATION FOR THIS POTENTIAL'
         stop
C         CALL SCL(N,X,EPS,C,SIG,BOXLX,BOXLY,BOXLZ,CUTOFF)
C         PRINT*,'Energy minimised with respect to lattice constants' 
C        CUTOFF=BOXLX/2.0D0
         PRINT*,'New box length and cutoff=',BOXLX,CUTOFF
      ENDIF
C
C Store density matrix: In the case of the perfect fcc lattice,
C the infinitely extended crystal implies that every RHO(J) is
C equal to RHO(1).
C
      DO 11 I=1,N
         
         rho_temp=0.d0
         vpot_temp=0.d0
         
         DO 122 J=1,ndir(I)
             
             ja=neigh_type(I,J)
             
             rbuf=dsqrt(Rneigh(I,J,1)**2+Rneigh(I,J,2)**2+Rneigh(I,J,3)**2)
             
             rho_temp  = rho_temp  +  rho_pot(ipot,rbuf)
             vpot_temp = vpot_temp +  Vpot(ipot,rbuf)
             
122      CONTINUE
       RHO(I)=rho_temp
       VSITE(I)=vpot_temp
!C        write(*,*) I, RHO(I)
11    CONTINUE
C
C calculate the potential energy:
C
      POTA=0.0D0
      POTB=0.0D0
      DO 13 I=1,N
        POTA=POTA+VSITE(I)
        POTB=POTB - Fembed(ipot,RHO(I))
13    CONTINUE
      PSC=POTA - POTB

Cdebug      write(*,*) POTA, -POTB, PSC
C      stop
      
      force(:,:)=0.d0
      
      IF (GTEST.OR.STEST) then
          if (GTEST) then
             do I=1,N
                Fembed_d_i=Fembed_d(ipot,RHO(I))
                do J=1,ndir(I)
                   ja=neigh_type(I,J)
                   rbuf=dsqrt(Rneigh(I,J,1)**2+Rneigh(I,J,2)**2+Rneigh(I,J,3)**2)
   
                   force(1,I)= force(1,I)+Rneigh(I,J,1)*(2*Vpot_d(ipot,rbuf)+      !&
     1              (Fembed_d_i+ Fembed_d(ipot,RHO(ja)) )*                 !&
     1                     rho_pot_d(ipot,rbuf)  ) /rbuf
                   force(2,I)= force(2,I)+Rneigh(I,J,2)*(2*Vpot_d(ipot,rbuf)+      !&
     1               (Fembed_d_i+ Fembed_d(ipot,RHO(ja)) )*                !&
     1                     rho_pot_d(ipot,rbuf)  ) /rbuf
                   force(3,I)= force(3,I)+Rneigh(I,J,3)*(2*Vpot_d(ipot,rbuf)+      !&
     1               (Fembed_d_i+ Fembed_d(ipot,RHO(ja)) )*                !&
     1                     rho_pot_d(ipot,rbuf)  ) /rbuf
                end do  
                Vnew(3*(I-1)+1)=-force(1,I)   
                Vnew(3*(I-1)+2)=-force(2,I)   
                Vnew(3*(I-1)+3)=-force(3,I) 
             end do    
          end if
          
          IF (STEST) then
             write(*,*) 'THERE is no impletantion for this case'
             write(*,*) 'GTEST', GTEST 
             write(*,*) 'STEST', STEST
             stop 
          ENDIF
      ENDIF 

      RETURN
      END
