c This subroutine accepts R in Angstroms. R is the distance between 
c centers of mass.
c Six angles are also accepted.
c The interaction energy in kcal/mol is returned in the variable 'val'.
      subroutine WATERMETHANE(COORDS,val)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=150,maxp=3000)
      parameter (nsitemax=15,ntypemax=5,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=5,maxpar2=10)
      parameter (mxl=100)
c
      dimension oa(3), ob(3), values(20), valder(20), COORDS(6)
      dimension torquea(3), torqueb(3), forceab(3)
c Only for the purpose of checking.....
      dimension valuesp(20),valuesm(20)
c
c
c The following common to pass the vector of linear parameters...
c
      common/leastsq/dat(maxp),a(maxp,maxb),a0(maxp),g(maxp),c(maxb),
     1 ata(maxb*maxb),b(maxb),works(2*maxb),iworks(maxb),np,nb
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
      common/temp_sites/ siteat(3,nsitemax), sitebt(3,nsitemax),
     1                   tranA(3,3), tranB(3,3)
      common/fit/rpt(6,maxp),dat0(maxp),sigma(maxp),iweight,ntpot,
     1 idonl, iopt
      common/type_used/ itypu(ntypemax,ntypemax)
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      common/misc/ R_0,isyst,npowers,numlin
      data bohr2a /0.529177249d0/
      data mode/-1/
      save
c
c If mode = -1 -- read input
c
      R=COORDS(1)
      OA(1)=0.0D0
      OA(2)=COORDS(2)
      OA(3)=COORDS(3)
      OB(1)=COORDS(4)
      OB(2)=COORDS(5)
      OB(3)=COORDS(6)
      if(mode.eq.-1) then
       call data
       mode=4
      endif
        call ang2cart(R,oa,ob)
      
      none=0
       val=object(none)
       if(val.le.-1.2d0)then
            write(*,*)R,val
           do ik=1,8
            write(*,*)(siteat(ki,ik),ki=1,3)
           end do
           write(*,*)'*************************************'
           do ik=1,9
            write(*,*)(sitebt(ki,ik),ki=1,3)
           end do
       endif
       RETURN
      
      end
c
c maxp - maximum nuber of data points 
c nsitemax - maximum number of sites in one monomer
c ntypemax - maximum number of site types
c maxb - maximum nuber of "basis functions", must be >=
c        0.5*nsitemax*(nsitemax+1)*number_of_terms_in_the_pair_potential
c maxpar1 - maximum number of types of one-site-type nonl parameters
c maxpar2 - maximum number of types of two-site-type nonl parameters
c
c nopar - number of optimized one-site-type nonlinear parameters
c noparab - number of optimized two-site-type nonlintypesnear parameters
c
c nparm - number of nonzero one-site-type nonlinear parameters
c nparab - number of nonzero two-site-type nonlintypesnear parameters
c
      subroutine data
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=150,maxp=3000)
      parameter (nsitemax=15,ntypemax=5,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=5,maxpar2=10)
      parameter (mxl=100)
      parameter (maxn=6*nsitemax+maxpar1*ntypemax+maxpar2*ntypemax2)
c
      character*8 label(9)
      dimension com(3)
      dimension chara(nsitemax),charb(nsitemax)
      common/npar/param(maxpar1,ntypemax),
     1 parab(maxpar2,ntypemax,ntypemax),nparm,nparab
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
      common/types/ itypea(nsitemax),itypeb(nsitemax)
      common/optdrv/ iosita(2,3*nsitemax),iositb(2,3*nsitemax),
     1 iopar(2,ntypemax),ioparab(3,ntypemax2),nosa,nosb,nopar,noparab
        common /charges/chara,charb
        common /amass/amassa(nsitemax),amassb(nsitemax)
c     
      data bohr2a /0.529177249d0/
      open(98,file='dimer',status='old')

c read in the sites of monomer A
c
      read(98,*) nsitea
      do i=1,nsitea
       read(98,*)chara(i),(sitea(j,i),j=1,3),amassa(i)
c transform to Angstroms...
       do j=1,3
        sitea(j,i) = bohr2a*sitea(j,i)
       end do
      end do
c Transfrom to COM 
      amtot=0.0d0
      do i=1,3
         com(i)=0.0d0
      end do
      
      do i=1,nsitea
         ams=amassa(i)
         com(1)=com(1)+ams*sitea(1,i)
         com(2)=com(2)+ams*sitea(2,i)
         com(3)=com(3)+ams*sitea(3,i)
         amtot=amtot+ams
      end do

c        write(*,*)'Here6',(sitea(2,jk),jk=1,nsitea)       
     
      do i=1,3
         com(i)=com(i)/amtot
      end do
      
c Shift to COM
      do i=1,nsitea
         do j=1,3
            sitea(j,i)=sitea(j,i)-com(j)
         end do
      end do
      
c     
c read in the sites of monomer B
c
      read(98,*) nsiteb
      do i=1,nsiteb
       read(98,*)charb(i),(siteb(j,i),j=1,3),amassb(i)
c transform to Angstroms...
       do j=1,3
        siteb(j,i) = bohr2a*siteb(j,i)
       end do
      end do
c shld shift COM for B also but ... it's CH4 & it's in COM already...
c for some other system, pls shift! :-)      
       close(98)
c       write(*,*)'Here5',(sitea(2,jk),jk=1,nsitea),
c     1           (siteb(1,jk),jk=1,nsiteb)
       return
       end
c
c For a given set of Euler angles (in radians) and the intermonomer separation R
c compute Cartesian coordinates of all sites. Assume that the initial
c coordinates of the monomer X sites are given in the center of mass system 
c of monomer X. Output in common/temp_sites/
c
      subroutine ang2cart(R,oa,ob)
      implicit real*8 (a-h,o-z)
c
      parameter (maxb=150,maxp=3000)
      parameter (nsitemax=15,ntypemax=5,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=5,maxpar2=10)
      parameter (mxl=100)
c
      dimension trand(3,3), tranc(3,3),oa(3),ob(3)
c
      common/sites/ sitea(3,nsitemax), siteb(3,nsitemax),nsitea,nsiteb
c
      common/temp_sites/ siteat(3,nsitemax), sitebt(3,nsitemax),
     1                   trand,tranc
c
c
c Build the rotation matrix for monomer A...
c
c      ad = 0.d0
      ad = oa(1)
      bd = oa(2)
      gd = oa(3)
      sad = dsin(ad)
      sbd = dsin(bd)
      sgd = dsin(gd)
      cad = dcos(ad)
      cbd = dcos(bd)
      cgd = dcos(gd)
      trand(1,1) = cad*cbd*cgd - sad*sgd
      trand(1,2) = -cad*cbd*sgd - sad*cgd
      trand(1,3) = cad*sbd
      trand(2,1) = sad*cbd*cgd + cad*sgd
      trand(2,2) = -sad*cbd*sgd + cad*cgd
      trand(2,3) = sad*sbd
      trand(3,1) = -sbd*cgd
      trand(3,2) = sbd*sgd
      trand(3,3) = cbd
c
c Build the rotation matrix for monomer B...
c
c      ac = ob(1)-oa(1)
      ac = ob(1)
      bc = ob(2)
      gc = ob(3)
      sac = dsin(ac)
      sbc = dsin(bc)
      sgc = dsin(gc)
      cac = dcos(ac)
      cbc = dcos(bc)
      cgc = dcos(gc)
      tranc(1,1) = cac*cbc*cgc - sac*sgc
      tranc(1,2) = -cac*cbc*sgc - sac*cgc
      tranc(1,3) = cac*sbc
      tranc(2,1) = sac*cbc*cgc + cac*sgc
      tranc(2,2) = -sac*cbc*sgc + cac*cgc
      tranc(2,3) = sac*sbc
      tranc(3,1) = -sbc*cgc
      tranc(3,2) = sbc*sgc
      tranc(3,3) = cbc
c
c Transform the coordinates of monomer A sites...
c
      do i=1,nsitea
       call matvec1(trand,sitea(1,i),siteat(1,i))
      end do
c      write(*,*)'Here3',(siteat(2,kia),kia=1,nsitea),'Here4',
c     1          (sitea(2,kia),kia=1,nsitea),(trand(2,ki),ki=1,3)
c
c Transform the coordinates of monomer B sites...
c
      do i=1,nsiteb
       call matvec1(tranc,siteb(1,i),sitebt(1,i))
      end do
c
c Shift monomer B by R in the positive z direction...
c
      do i=1,nsiteb
       sitebt(3,i) = sitebt(3,i) + R
      end do
c
      return
      end
c
c ----- Multiply vetor v by the matrix a. Store in u.
c
       subroutine matvec1(a,v,u)
       implicit real*8 (a-h,o-z)
       dimension a(3,3),v(3),u(3)
c
       n = 3
       do i=1,n
        u(i) = 0.d0
        do j=1,n
         u(i) = u(i) + a(i,j)*v(j)
        end do
       end do
       return
       end
c
c ----- Multiply vetor v by the transposed matrix a. Store in v again.
c
       subroutine matvec2(a,v)
       implicit real*8 (a-h,o-z)
       dimension a(3,3),v(3),u(3)
c
       n = 3
       do i=1,n
        u(i) = 0.d0
        do j=1,n
         u(i) = u(i) + a(j,i)*v(j)
        end do
       end do
       do i=1,n
        v(i) = u(i)
       end do
       return
       end
c function to calculate the SAPT interaction energy between methane and
c water as described in the paper: J. Chem. Phys., vol 123, p 134311 (2005).
        function object(none)
        implicit real*8 (a-h,o-z)
      parameter(NQW=50)
      DOUBLE PRECISION b(NQW),rij(7,9),c(32),h(32),s(42)
       data istat /0/
        common /encoul/ecoul
        save
        if(istat.eq.0)then
           open(99,file='dispindCoeff',status='old')
           do ik=1,27
              read(99,*)c(ik)
           end do
           close(99)
           open(99,file='damp',status='old')
           do ik=1,18
              read(99,*)h(ik)
           end do
           close(99)
           open(99,file='init_params',status='old')
           do ik=1,12
              read(99,*)b(ik)
           end do
           do ik=1,36
              read(99,*)s(ik)
           end do
           close(99)           
           istat=1
          endif

        f=0.0d0
        ep=0.0d0
        ijkl=i

        call mat2(rij)
        if(rij(3,5).le.2.05d0)then
         zbb=rij(3,5)
         object=541.85744d0*dexp(-(zbb-2.05d0))/zbb
         RETURN
        endif

        do kW=6,7
           do km=1,4
         z=1.0d0/rij(kW,kM)
           ep=ep+dexp(-b(10)/z)*(s(28)+s(29)*z+s(30)/z)
           end do
        end do

        do kW=6,7
           do km=5,5
         z=1.0d0/rij(kW,kM)
           ep=ep+dexp(-b(11)/z)*(s(31)+s(32)*z+s(33)/z)
           end do
        end do

        do kW=6,7
           do km=6,9
         z=1.0d0/rij(kW,kM)
           ep=ep+dexp(-b(12)/z)*(s(34)+s(35)*z+s(36)/z)
           end do
        end do

        do kW=1,2
           do kM=1,4
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(1)*h(1)/z
              g2=h(10)*h(10)/z
          ep=ep+dexp(-b(1)/z)*(s(1)+s(2)*z+s(3)/z)
              f=f+dT(6,g)*c(1)*z6+dT(8,g2)*c(2)*z8+
     1     dT(10,g2)*c(3)*z10
           end do
        end do

        do kW=3,3
           do kM=1,4
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(2)*h(2)/z
              g2=h(11)*h(11)/z
              ep=ep+dexp(-b(2)/z)*(s(4)+s(5)*z+s(6)/z)
              f=f+dT(6,g)*c(4)*z6+dT(8,g2)*c(5)*z8+
     1     dT(10,g2)*c(6)*z10
           end do
        end do


         do kW=1,2
           do kM=5,5
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(3)*h(3)/z
              g2=h(12)*h(12)/z
              ep=ep+dexp(-b(3)/z)*(s(7)+s(8)*z+s(9)/z)
              f=f+dT(6,g)*c(7)*z6+dT(8,g2)*c(8)*z8+
     1      dT(10,g2)*c(9)*z10
           end do
        end do

        do kW=3,3
           do kM=5,5
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              z12=z10*z2
              g=h(4)*h(4)/z
              g2=h(13)*h(13)/z
              ep=ep+dexp(-b(4)/z)*(s(10)+s(11)*z+s(12)/z)
              f=f+dT(6,g)*c(10)*z6+dT(8,g2)*c(11)*z8+
     1     dT(10,g2)*c(12)*z10
           end do
        end do
        do kW=4,5
           do kM=1,4
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(5)*h(5)/z
              g2=h(14)*h(14)/z
              ep=ep+dexp(-b(5)/z)*(s(13)+s(14)*z+s(15)/z)
              f=f+dT(6,g)*c(13)*z6+dT(8,g2)*c(14)*z8+
     1      dT(10,g2)*c(15)*z10
           end do
        end do

        do kW=4,5
           do kM=5,5
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(6)*h(6)/z
              g2=h(15)*h(15)/z
              ep=ep+dexp(-b(6)/z)*(s(16)+s(17)*z+s(18)/z)
              f=f+dT(6,g)*c(16)*z6+dT(8,g2)*c(17)*z8+
     1       dT(10,g2)*c(18)*z10
           end do
        end do

        do kW=1,2
           do kM=6,9
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(7)*h(7)/z
              g2=h(16)*h(16)/z
              ep=ep+dexp(-b(7)/z)*(s(19)+s(20)*z+s(21)/z)
              f=f+dT(6,g)*c(19)*z6+dT(8,g2)*c(20)*z8+
     1     dT(10,g2)*c(21)*z10
           end do
        end do

        do kW=3,3
           do kM=6,9
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(8)*h(8)/z
              g2=h(17)*h(17)/z
              ep=ep+dexp(-b(8)/z)*(s(22)+s(23)*z+s(24)/z)
              f=f+dT(6,g)*c(22)*z6+dT(8,g2)*c(23)*z8+
     1     dT(10,g2)*c(24)*z10
           end do
        end do

        do kW=4,5
           do kM=6,9
              z=1.0d0/rij(kW,kM)
              z2=z*z
              z4=z2*z2
              z6=z4*z2
              z8=z6*z2
              z10=z8*z2
              g=h(9)*h(9)/z
              g2=h(18)*h(18)/z
              ep=ep+dexp(-b(9)/z)*(s(25)+s(26)*z+s(27)/z)
              f=f+dT(6,g)*c(25)*z6+dT(8,g2)*c(26)*z8+
     1     dT(10,g2)*c(27)*z10
           end do
        end do

          object=f+ep+ecoul
        RETURN
        END


      function dT(n,x)
      implicit real*8 (a-h,o-z)
      if(x.le.1.d-10)then
         dT=0.0d0
         RETURN
      endif

      dT=1.0d0
      xk=x
      do k=1,n
         dT=dT+xk
         xk=xk*x/dble(k+1)
      end do

      dT=1.0d0-dexp(-x)*dT
      RETURN
      end


        subroutine  mat2(rij)
        implicit real*8 (a-h,o-z)
        parameter(deg2rad=0.017453293d0)
        data istat /0/
        real*8 rij(7,9)
        parameter (nsmax=15)
        dimension chara(nsmax),charb(nsmax)
        data a0 / 0.529177249d0/, zero /0.d0/
        common /encoul/ecoul
        common/temp_sites/ siteat(3,nsmax), sitebt(3,nsmax),
     1       tranA(3,3), tranB(3,3)
        common /charges/chara,charb
        
        nsa=7
        nsb=9
      ecoul=0.0d0
      do i=1,nsa
         do j=1,nsb
            x=siteat(1,i)-sitebt(1,j)
            y=siteat(2,i)-sitebt(2,j)
            zp=siteat(3,i)-sitebt(3,j)
            rang=dsqrt(x*x+y*y+zp*zp)
            ecoul=ecoul+chara(i)*charb(j)/rang
            rij(i,j)=rang
         end do
      end do   
       ecoul=ecoul*627.51d0*a0 
      RETURN
      end

