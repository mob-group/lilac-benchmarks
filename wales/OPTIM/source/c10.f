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
C Program: c10.f

C Function: This subroutine calculates the energy, gradient, and second
C derivatives for a given configuration of the 46 particle polymer chain.
C Modified for a C10 configuration 17/6/03 by DJW.
C A configuration and number of particles is passed to the subroutine and
C the energy, gradient, and matrix of second derivatives is returned.

        subroutine c10(qo, n, grad, energy, gtest, stest)
        USE MODHESS
        IMPLICIT NONE
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        logical gtest, stest
        INTEGER N
        DOUBLE PRECISION QO(3*N), GRAD(3*N)
        DOUBLE PRECISION ENERGY
        INTEGER ntype(10)
        DOUBLE PRECISION A_PARAM(N,N),
     1  b_param(n,n),
     2  d_param(n),c_param(n),
     3  x(n), y(n), z(n),
     4  xr(n,n), yr(n,n), zr(n,n),
     5  dot_prod(n,3), x_prod(n),
     6  bond_angle(n), tor_angle(n), radii(n,n)

C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(10),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n), 
C    4  xr(n,n), yr(n,n), zr(n,n), 
C    5  dot_prod(n,3), x_prod(n), 
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        call c10param_array(a_param,b_param,c_param,d_param,n)
        call c10calc_int_coords(qo,n,a_param,b_param,c_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        call c10calc_energy(qo,energy,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        call c10calc_gradient(qo,grad,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)

        IF (.NOT.STEST) RETURN
        call c10calc_dyn(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)

        return
        end

C Calculate the Internal Coordinates

        subroutine c10calc_int_coords(qo,n,a_param,b_param,c_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        IMPLICIT NONE
        INTEGER N, I, J
        DOUBLE PRECISION QO(3*N), COS_THETA, COS_PHI
        INTEGER ntype(10)
        DOUBLE PRECISION A_PARAM(N,N),
     1  b_param(n,n),
     2  c_param(n),
     3  x(n), y(n), z(n),
     4  xr(n,n), yr(n,n), zr(n,n),
     5  dot_prod(n,3), x_prod(n),
     6  bond_angle(n), tor_angle(n), radii(n,n)
        
C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(10),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n),
C    4  xr(n,n), yr(n,n), zr(n,n),
C    5  dot_prod(n,3), x_prod(n),
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        do i = 1, n
        j = (i-1)*3
        x(i) = qo(j+1)
        y(i) = qo(j+2)
        z(i) = qo(j+3)
        enddo

C Inter-particle distances

        do i = 1, n-1
        do j = i+1, n
        xr(i,j) = x(j) - x(i)
        yr(i,j) = y(j) - y(i)
        zr(i,j) = z(j) - z(i)
        radii(i,j) = dsqrt(xr(i,j)*xr(i,j) + yr(i,j)*yr(i,j) + zr(i,j)*zr(i,j))
        radii(j,i) = radii(i,j)
        enddo
        enddo

C Dot products between bond vectors

        do i = 1, n-3
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + zr(i,i+1)*zr(i,i+1)

        dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+ zr(i,i+1)*zr(i+1,i+2)

        dot_prod(i,3) = xr(i,i+1)*xr(i+2,i+3)+yr(i,i+1)*yr(i+2,i+3)+ zr(i,i+1)*zr(i+2,i+3)
        enddo

        i = n-2
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + zr(i,i+1)*zr(i,i+1)

        dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+ zr(i,i+1)*zr(i+1,i+2)

        i = n-1
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + zr(i,i+1)*zr(i,i+1)

C Cross-products between adjacent bond vectors

        do i = 1, n-2
        x_prod(i) = dot_prod(i,1)*dot_prod(i+1,1) - dot_prod(i,2)*dot_prod(i,2)   
        enddo

C Bond angles

        do i = 1, n-2
        cos_theta=-dot_prod(i,2)/(dsqrt(dot_prod(i,1)
     1  *dot_prod(i+1,1)))
        bond_angle(i+1) = dacos(cos_theta)
        enddo

C Torsional angles

        do i = 1, n-3
        cos_phi = (dot_prod(i,2)*dot_prod(i+1,2) -
     1  dot_prod(i,3)*dot_prod(i+1,1))/dsqrt(x_prod(i)*x_prod(i+1))
        IF (ABS(cos_phi).GT.1.0D0) cos_phi=cos_phi/abs(cos_phi)
        tor_angle(i+1) = dacos(cos_phi)
C       WRITE(*,'(A,I4,4F20.10)') 'i,tor_angle,cos_phi,dacos=',i,tor_angle(i+1),cos_phi,dacos(cos_phi)
        enddo

        return
        end


C Calculate the Energy

        subroutine c10calc_energy(qo,energy,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        IMPLICIT NONE
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA
        DOUBLE PRECISION ENERGY, S6, E_NBOND, E_BOND, E_BANGLE, E_TANGLE, RAD6
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        INTEGER N, I, J
        DOUBLE PRECISION QO(3*N)
        INTEGER ntype(10)
        DOUBLE PRECISION A_PARAM(N,N),
     1  b_param(n,n),
     2  d_param(n),c_param(n),
     3  x(n), y(n), z(n),
     4  xr(n,n), yr(n,n), zr(n,n),
     5  dot_prod(n,3), x_prod(n),
     6  bond_angle(n), tor_angle(n), radii(n,n)

C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(10),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n),
C    4  xr(n,n), yr(n,n), zr(n,n),
C    5  dot_prod(n,3), x_prod(n),
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        s6 = sigma*sigma*sigma*sigma*sigma*sigma
        e_nbond=0.0D0
        e_bond=0.0D0
        e_bangle=0.0D0
        e_tangle=0.0D0

        do i = 1, n-2
        do j = i+2, n

        rad6 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*
     1  radii(i,j)

        e_nbond = e_nbond + 4.0*((a_param(i,j)*s6*s6/(rad6*rad6)) + 
     1  (b_param(i,j)*s6/rad6))

        enddo
        enddo

        do i = 1, n-1

        e_bond = e_bond + 0.5*rk_r*(radii(i,i+1)-sigma)*
     1  (radii(i,i+1)-sigma)

        enddo

        
        do i = 2, n-1

        e_bangle = e_bangle + 0.5*rk_theta*(bond_angle(i)-theta_0)
     1  *(bond_angle(i)-theta_0)

        enddo

        do i = 2, n-2

        e_tangle = e_tangle + c_param(i)*(1.0 + cos(tor_angle(i))) 
     1  + d_param(i)*(1.0 + cos(3.0*tor_angle(i)))

        enddo

        energy = e_nbond + e_bond + e_bangle + e_tangle
C       WRITE(*,'(A,4F20.10)') 'nbond,bond,bangle,tangle=',e_nbond,e_bond,e_bangle,e_tangle

        return
        end

C Calculate the gradiants

        subroutine c10calc_gradient(qo,fq,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        IMPLICIT NONE
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6,theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        INTEGER J, I, N, NTYPE(10)
        DOUBLE PRECISION QO(3*N),FQ(3*N),FX(N),FY(N),
     1  fz(n)
        DOUBLE PRECISION FNB_X(N),FNB_Y(N),FNB_Z(N),
     1  fb_x(n),fb_y(n)
        DOUBLE PRECISION FB_Z(N),FBA_X(N),FBA_Y(N),FBA_Z(N)
        DOUBLE PRECISION FTA_X(N),FTA_Y(N),FTA_Z(N)
        DOUBLE PRECISION A4, COEF3, COEF2, COEF1, COEF, A3, DEN2, A2, A1, DEN1, RNUM, DEN, RVAR,
     1                   FZZ, FYY, FXX, DF, RAD14, RAD7, S6
        DOUBLE PRECISION A_PARAM(N,N),
     1  b_param(n,n),
     2  d_param(n),c_param(n),
     3  x(n), y(n), z(n),
     4  xr(n,n), yr(n,n), zr(n,n),
     5  dot_prod(n,3), x_prod(n),
     6  bond_angle(n), tor_angle(n), radii(n,n)

C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(10),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n),
C    4  xr(n,n), yr(n,n), zr(n,n),
C    5  dot_prod(n,3), x_prod(n),
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        s6 = sigma*sigma*sigma*sigma*sigma*sigma

C Gradients of potential

        do i = 1,n

        fnb_x(i) = 0.0  
        fnb_y(i) = 0.0 
        fnb_z(i) = 0.0 

        fb_x(i)  = 0.0 
        fb_y(i)  = 0.0 
        fb_z(i)  = 0.0 

        fba_x(i) = 0.0 
        fba_y(i) = 0.0 
        fba_z(i) = 0.0 

        fta_x(i) = 0.0 
        fta_y(i) = 0.0 
        fta_z(i) = 0.0 

        fx(i)= 0.0 
        fy(i)= 0.0 
        fz(i)= 0.0 

        enddo

C ..... Non-bonded interaction forces ..... 

        do i = 1, n-2
        do j = i+2, n

        rad7 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*
     1        radii(i,j)*radii(i,j)*radii(i,j)   
        rad14 = rad7*rad7 

        df = -24.0*((2.0*a_param(i,j)*s6*s6/rad14) + 
     1              (b_param(i,j)*s6/(rad7*radii(i,j))))

        fxx = df*xr(i,j) 
        fyy = df*yr(i,j) 
        fzz = df*zr(i,j) 

        fnb_x(i) = fxx + fnb_x(i)
        fnb_y(i) = fyy + fnb_y(i)
        fnb_z(i) = fzz + fnb_z(i)

        fnb_x(j) = -fxx + fnb_x(j)
        fnb_y(j) = -fyy + fnb_y(j)
        fnb_z(j) = -fzz + fnb_z(j)

        enddo
        enddo

C ... Bond interaction forces ... 

        do i = 1, n-1

        rvar = sigma/radii(i,i+1) 

        df = rk_r*(1.0 - rvar) 
        fxx = df*xr(i,i+1) 
        fyy = df*yr(i,i+1) 
        fzz = df*zr(i,i+1) 

        fb_x(i) = fxx + fb_x(i)
        fb_y(i) = fyy + fb_y(i)
        fb_z(i) = fzz + fb_z(i)

        fb_x(i+1) = -fxx + fb_x(i+1)
        fb_y(i+1) = -fyy + fb_y(i+1)
        fb_z(i+1) = -fzz + fb_z(i+1)

        enddo

C bond angle forces  particle 1
C particles 1,2,n-1, and n done outside of the loop

        i = 1
        den = dsin(bond_angle(i+1))*dsqrt(dot_prod(i+1,1)*dot_prod(i,1))

        rnum = rk_theta*(bond_angle(i+1) - theta_0)

        fba_x(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) -
     1        xr(i+1,i+2))/den

        fba_y(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*yr(i,i+1) -
     1        yr(i+1,i+2))/den

        fba_z(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*zr(i,i+1) -
     1        zr(i+1,i+2))/den


C particle 2

        i = 2
        den = dsin(bond_angle(i))
     1        *dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*dsqrt(dot_prod(i+1,1)
     1         *dot_prod(i,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        fba_x(i) = a1 + a2 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        fba_y(i) = a1 + a2 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        fba_z(i) = a1 + a2 

C particles 3 thru n-2 

        do i = 3, n-2

        den = dsin(bond_angle(i))*
     1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*
     1               dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        den2 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1         *dot_prod(i-1,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

        fba_x(i) = a1 + a2 + a3 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

        fba_y(i) = a1 + a2 + a3 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

        fba_z(i) = a1 + a2 + a3 

        enddo

C particle n-1 

        i = n-1
        den = dsin(bond_angle(i))*
     1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1         *dot_prod(i-1,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den1

        fba_x(i) = a1 + a2

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den1

        fba_y(i) = a1 + a2

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den1

        fba_z(i) = a1 + a2

C particle n

        i = n
        den = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1        *dot_prod(i-1,1))

        fba_x(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) 
     1        - xr(i-2,i-1))/den

        fba_y(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*yr(i-1,i) 
     1        - yr(i-2,i-1))/den

        fba_z(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*zr(i-1,i) 
     1        - zr(i-2,i-1))/den

C Torsional angle forces
C particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
C particle 1

        i = 1
             coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1         *dcos(tor_angle(i+1))-3.0))
     1  *(1.0/dsqrt(x_prod(i+1)*x_prod(i)))  

        fta_x(i) = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1         dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        fta_y(i) = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        fta_z(i) = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 


C particle 2

        i = 2
        coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 3.0))
     1        *(1.0/dsqrt(x_prod(i+1)*x_prod(i)))  

             coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i)*x_prod(i-1)))  

        a1 =  -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        fta_x(i) = a1 + a2 

        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2)))

        fta_y(i) = a1 + a2 
        
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        fta_z(i) = a1 + a2 

C particle 3

        i = 3
        coef=(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i+1)*x_prod(i)))  

        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        fta_x(i) = a1 + a2 + a3 
 
        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 
        
        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)))

        fta_y(i) = a1 + a2 + a3 
 
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 =  -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        fta_z(i) = a1 + a2 + a3 

C particles 4 to n-3

        do i = 4, n-3

        coef = (c_param(i+1) + d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 3.0))*(1.0/dsqrt(x_prod(i+1)*x_prod(i)))

        coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) -3.0))*(1.0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2 = (c_param(i-1) + d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1  3.0))*(1.0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3 = (c_param(i-2) + d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1  3.0))*(1.0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2 + a3 + a4 

        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i))) 

        a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2 + a3 + a4 

        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 
        
        a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 + a3 + a4 

        enddo

C particle n-2

        i = n-2
        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 3.0))*(1.0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 3.0))*(1.0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 3.0))*(1.0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) + dot_prod(i,2)*xr(i+1,i+2)))


        a2 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 
        
        a3 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2 + a3 

        a1 =  -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +  
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 =  -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)))

        a3 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2 + a3 
 
        a1 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +  
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        a3 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 + a3 

C particle n-1

        i = n-1
        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - 
     1        dot_prod(i-2,2)*xr(i-1,i) +
     1        dot_prod(i-1,2)*xr(i-2,i-1) +  dot_prod(i-1,1)*xr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) - 
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2  

        a1 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) - 
     1        dot_prod(i-2,2)*yr(i-1,i) +
     1        dot_prod(i-1,2)*yr(i-2,i-1) +  dot_prod(i-1,1)*yr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1  dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - 
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2  

        a1 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - 
     1        dot_prod(i-2,2)*zr(i-1,i) +
     1        dot_prod(i-1,2)*zr(i-2,i-1) +  dot_prod(i-1,1)*zr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - 
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 
 
C particle n

        i = n
        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1        3.0))*(1.0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        fta_x(i) = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) 
     1        - dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_y(i) = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - 
     1        dot_prod(i-2,1)*yr(i-3,i-2) 
     1        - (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) 
     1        - dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_z(i) = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - 
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

C Total up the gradients

        do i = 1, n
        fx(i) = fnb_x(i) + fb_x(i) + fba_x(i) + fta_x(i) 
        fy(i) = fnb_y(i) + fb_y(i) + fba_y(i) + fta_y(i) 
        fz(i) = fnb_z(i) + fb_z(i) + fba_z(i) + fta_z(i) 
        enddo

        do i = 1, n
        j = (i-1)*3
        fq(j+1) = -fx(i)
        fq(j+2) = -fy(i)
        fq(j+3) = -fz(i)
        enddo

        return
        end


C Calculate the Second Derivative matrix (two-sided numerical approach)

        subroutine c10calc_dyn(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        USE MODHESS
        IMPLICIT NONE
        INTEGER ntype(10), N, I, J
        DOUBLE PRECISION QO(3*N), FQ1(3*N), FQ2(3*N)
        DOUBLE PRECISION A_PARAM(N,N), RK_THETA, THETA_0, SIGMA, RMASS, EPSILON, DELTA, 
     1  b_param(n,n), rk_r,
     2  c_param(n), d_param(n),
     3  x(n), y(n), z(n),
     4  xr(n,n), yr(n,n), zr(n,n),
     5  dot_prod(n,3), x_prod(n),
     6  bond_angle(n), tor_angle(n), radii(n,n)
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4 ,delta=1.0d-4, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)

C Fill in the Hessian matrix

        do j = 1, 3*n

        qo(j) = qo(j) + delta
        call c10calc_int_coords(qo,n,a_param,b_param,c_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        call c10calc_gradient(qo,fq2,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        qo(j) = qo(j) - 2.0*delta
        call c10calc_int_coords(qo,n,a_param,b_param,c_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        call c10calc_gradient(qo,fq1,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        qo(j) = qo(j) + delta

        do i = j, 3*n

        HESS(i,j) = (fq2(i) -  fq1(i))/(2.0*delta)
        HESS(j,i) = HESS(i,j)

        enddo
        enddo

        return
        end

C Fill the parameter arrays

        subroutine c10param_array(a_param,b_param,c_param,d_param,n)
        IMPLICIT NONE
        INTEGER ntype(10), N, ICOUNT, J, I
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), EPSILON
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N)
        parameter (epsilon = 0.0100570D0)

C Amino Acid types

        ntype(1) = 1
        ntype(2) = 1
        ntype(3) = 1
        ntype(4) = 1
        ntype(5) = 1
        ntype(6) = 1
        ntype(7) = 1
        ntype(8) = 1
        ntype(9) = 1
        ntype(10) = 1

C Parameters for the dihedral angle potential

        do i = 1, n-3
        icount = 0

        do j = 0,3
        if(ntype(i+j) .eq. 3)then
        icount = icount + 1
        endif
        enddo

        if(icount .ge. 2)then
        c_param(i+1) = 0.0
        d_param(i+1) = 0.2*epsilon
        else
        c_param(i+1) = 1.2*epsilon
        d_param(i+1) = 1.2*epsilon
        endif

        icount = 0

        enddo

C  Parameters for the L-J interaction between non-bonded particles

        do i = 1, n-1
        do j = i+1, n

        if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3)then
        a_param(i,j) = 1.0*epsilon 
        b_param(i,j) = 0.0 
        a_param(j,i) = 1.0*epsilon 
        b_param(j,i) = 0.0

        elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
        a_param(i,j) =  epsilon
        b_param(i,j) = -epsilon 
        a_param(j,i) =  epsilon
        b_param(j,i) = -epsilon
        
        else

        a_param(i,j) = epsilon*2.0/3.0 
        b_param(i,j) = epsilon*2.0/3.0 
        a_param(j,i) = epsilon*2.0/3.0 
        b_param(j,i) = epsilon*2.0/3.0 

        endif

        enddo
        enddo

        return
        end
