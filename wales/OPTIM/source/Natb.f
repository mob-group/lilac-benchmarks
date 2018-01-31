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
C       TB potential for Na, Ag, Li clusters

        SUBROUTINE NATB(NSIZE, P, GRAD, E0, GTEST, STEST)
        IMPLICIT NONE
        INTEGER ION, NSIZE
        DOUBLE PRECISION, PARAMETER :: EPSILON=1D-10
        LOGICAL GTEST, STEST
        DOUBLE PRECISION P(3*NSIZE), GRAD(3*NSIZE), E0
        character(LEN=2) metal
        common/ION/ion
        common/METAL/metal

        IF (STEST) THEN
           PRINT '(A)','ERROR - second derivatives not available in NATB routine'
           STOP
        ENDIF

        metal='NA'
        ion=0

C       if (metal.eq.'NA') then
C          print *,'Parameters optimized for SODIUM'
C       elseif (metal.eq.'AG') then
C          print *,'Parameters optimized for SILVER'
C       elseif (metal.eq.'LI') then
C          print *,'Parameters optimized for LITHIUM'
C       endif

C       if (ion.eq.1) then
C          print *,'CATION'
C       elseif (ion.eq.0) then
C          print *,'NEUTRAL'
C       else
C          print *,'ANION'
C       endif

        CALL entots(P,grad,nsize,e0)

        end
c______________________________________________________________________

        SUBROUTINE entots(P,DP,nbat,e0)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        IMPLICIT INTEGER(I-N)
        PARAMETER(ANM=41901.052,RKB=3.165D-6,erepmax=100.d0)

        DIMENSION P(3*NBAT)
        DIMENSION GRAD(3*NBAT),DP(3*NBAT)

        DIMENSION HDER(3*NBAT,NBAT*(NBAT+1)/2),R(NBAT*(NBAT+1)/2)
        DIMENSION HA(NBAT*(NBAT+1)/2)
        DIMENSION VECT(NBAT*(NBAT+1)/2),RINV(NBAT*(NBAT+1)/2)
        DIMENSION H(NBAT*(NBAT+1)/2),DIF(3*NBAT*(NBAT+1)/2)
        DIMENSION DBSZ(3*NBAT*(NBAT+1)/2),UNIT(3*NBAT*(NBAT+1)/2)
        DIMENSION VALP(NBAT),VT(NBAT,NBAT),NUM(NBAT)
        DOUBLE PRECISION HMAT(NBAT,NBAT)
        INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NBAT)
        INTEGER IWORK(33*3*NBAT), INFO, ISTAT
        DOUBLE PRECISION WORK(33*3*NBAT), ABSTOL, DLAMCH
        integer ion
        common/ION/ion

        DIMENSION DI(3),DJ(3)

        LOGICAL MONOC
        character*2 metal
        common/METAL/metal

C        PARAMETERS

         LWORK=33*3*NBAT
         ILWORK=33*3*NBAT

        if (metal.eq.'NA') then
           ESP=0.077307
           RMIN=3.
           RMAX=100.
        elseif (metal.eq.'AG') then
           ESP=0.227857
           RMIN=4.
           RMAX=20.
        elseif (metal.eq.'LI') then
           ESP=0.067908d0
           RMIN=0.
           RMAX=1000.
        endif

        NALK=NBAT
        IDIAG=1
        DISTOR=0.D0
        HH=1D-6
        FTOL=1D-4
        ITMAX=100
        IPARAM=1
        ASZ=1
        NGR=0
        IPRTH=0
        IPRTG=0
        NBLEC=NBAT-ion
        apol=0.

        ndimh=nbat ! changed by DJW

        do ij=1,nbat
           h(ij)=0.D0
           do ii=1,3*nbat
              hder(ii,ij)=0.D0
           enddo
        enddo
        do i=1,3*nbat
           grad(i)=0.D0
        enddo
        DENOM=1.D0/esp 
        
        IJ=0
        IJQ=0
        DO 2680 I=1,NBAT
           DO 2682 J=1,I-1
              IJ=IJ+1
              IJX=IJQ+1
              IJY=IJQ+2
              IJZ=IJQ+3 
              DIF(IJX)=P(3*J-2)-P(3*I-2)
              DIF(IJY)=P(3*J-1)-P(3*I-1)
              DIF(IJZ)=P(3*J)-P(3*I)
              RIJ=DIF(IJX)*DIF(IJX)+DIF(IJY)*DIF(IJY)+DIF(IJZ)*DIF(IJZ)
              rij=dsqrt(rij)
              R(IJ)=RIJ
              rinv(ij)=1.d0/rij
              UNIT(IJX)=DIF(IJX)*RINV(IJ)
              UNIT(IJY)=DIF(IJY)*RINV(IJ)
              UNIT(IJZ)=DIF(IJZ)*RINV(IJ)
 2682           IJQ=IJQ+3
           ij=ij+1
           ijq=ijq+3
 2680        CONTINUE
        
c-----------------------------
C       INITIALISATION 
C------------------------------

        DO 5500 IJ=1,NBAT*(NBAT+1)/2
           DO 5500 K=1,3*NBAT
 5500   HDER(K,IJ)=0.D0

C       
C-------------------------------------------------------
C       PRECALCUL ET STOCKAGE DES TERMES DE PERTURBATION
C-------------------------------------------------------
C       
        IJQ=0
        IJ=0
        DO 5700 I=1,NALK
           DO 5720 J=1,I-1
              IJ=IJ+1
              RIJ=R(IJ)
              IF ((RIJ.LT.RMAX).and.(rij.ge.rmin)) THEN
                 call ftsz(rij,val,der)
              else if (rij.lt.rmin) then
                 if (metal.eq.'NA') then
                    val=-3.724D-10*dexp(3.32*rij)
                    der=3.32*val
                 elseif (metal.eq.'AG') then
                    val=-3.724D-10*dexp(3.32*rij)
                    der=3.32*val
                 elseif (metal.eq.'LI') then
                    val=-3.724D-10*dexp(3.32*rij)
                    der=3.32*val
                 endif
              ELSE
                 val=0.d0
                 der=0.D0
              ENDIF
 1346              vect(ij)=val
C       
C       DERIVEE DE BSZ
C       
              RINV(IJ)=1.D0/RIJ 
              do l=1,3
                 ijql=ijq+l
                 DBSZ(ijql)=DER*unit(ijql)
              enddo
              IJQ=IJQ+3
 5720           CONTINUE
           IJQ=IJQ+3
           IJ=IJ+1
 5700        CONTINUE 
C       
        
C----------------------------------------------
C       HAMILTONIEN ET gradient DE L HAMILTONIEN
C-----------------------------------------------
        IJ=0
        III=0
        ijq=0
        DO 30 I=1,NALK
           JJJ=0
C-------------------------------------------
C       ELEMENTS EXTRA-DIAGONAUX (BETA EFFECTIF)
C-------------------------------------------
           DO 32 J=1,I-1
              IJ=IJ+1
C       
C       HAMILTONIEN SS
C       
              RIJ=R(IJ)
              IF(RIJ.GE.RMAX) THEN
                 hss=0.D0
                 der=0.D0
              ELSE
                 if (rij.ge.rmin) then
                    call ftss(rij,val,der)
                 else
                    if (metal.eq.'NA') then
                       val=-3.743D-8*dexp(2.544*rij)
                       der=2.544*val
                    elseif (metal.eq.'AG') then
                       val=-3.743D-8*dexp(2.544*rij)
                       der=2.544*val
                    elseif (metal.eq.'LI') then
                       val=-3.743D-8*dexp(2.544*rij)
                       der=2.544*val
                    endif
                 endif
                 hss=VAL
              ENDIF
C
C       gradient DES TERMES SS
C
              IDER=0
              DO 2504 l=1,3 
                 ijql=ijq+l
                 derq=der*unit(ijql)
                 HDER(IDER+I,IJ)=-DERQ
                 HDER(IDER+J,IJ)=DERQ
                 IDER=IDER+NBAT
 2504              CONTINUE    
              DO 2660 L=1,3
                 DI(L)=0.D0
 2660              DJ(L)=0.D0
              KKK=0  
C       
C       HAMILTONIEN PERTURBE SP (TERMES A 3 CORPS)
C       
              Hij=0.D0
              DO 2600 K=1,NALK
                 IF((K.EQ.J).OR.(K.EQ.I)) GO TO 2605
                 IF(I.GT.K) THEN
                    isi=1
                    IK=III+K
                 ELSE
                    IK=KKK+I
                    isi=-1
                 ENDIF
                 IF(J.GT.K) THEN
                    isj=1
                    JK=JJJ+K
                 ELSE
                    JK=KKK+J
                    isj=-1
                 ENDIF
                 IKQ=ik+ik+ik-3
                 JKQ=jk+jk+jk-3
                 SCAL=UNIT(IKQ+1)*UNIT(JKQ+1)+UNIT(IKQ+2)*UNIT(JKQ+2)
     +               +UNIT(IKQ+3)*UNIT(JKQ+3)
                 scal=scal*isi*isj
                 BIJK=VECT(IK)*VECT(JK)
                 HIJK=BIJK*SCAL
                 Hij=Hij-HIJK
                 RIJK=RINV(IK)*RINV(JK)
C       
C       gradient DES TERMES SP
C       
                 IDER=0
                 DO 2610 L=1,3
                    ikql=ikq+l
                    jkql=jkq+l
                    uik=isi*unit(ikql)
                    ujk=isj*unit(jkql)
                    aijk=bijk*rijk
                    DII=-DBSZ(IKQL)*vect(jk)*isi*scal
     1                    -isj*dif(jkql)*aijk
     2                    +hijk*uik*rinv(ik)
                    Djj=-DBSZ(jKQL)*vect(ik)*isj*scal
     1                 -isi*dif(ikql)*aijk
     2           +hijk*ujk*rinv(jk)
                    dk=-dii-djj
                    DI(L)=DI(L)+DII
                    DJ(L)=DJ(L)+DJJ
                    HDER(IDER+K,IJ)=-DK*denom
 2610            IDER=IDER+NBAT
 2605                 KKK=KKK+K      
 2600              CONTINUE
              
 3030              format( 3(2x,i3),3(2x,f10.6))
              H(IJ)=HSS+hij*DENOM
              IDER=0 
              DO 2670 L=1,3
                 HDER(IDER+I,IJ)=hder(ider+i,ij)-DI(L)*DENOM
                 HDER(IDER+J,IJ)=hder(ider+j,ij)-DJ(L)*DENOM
 2670         IDER=IDER+NBAT
              JJJ=JJJ+J
              ijq=ijq+3
 32           CONTINUE
C       
C--------------------------------------------------------
C       ELEMENTS DIAGONAUX (REPULSION)
C----------------------------------------------------------
C       
           IJ=IJ+1
           kkk=0
           hii=0.D0 
           DO 223 K=1,nalk
              if(k.eq.i) go to 221
              if(i.gt.k) then
                 ik=iii+k
                 isi=1
              else
                 ik=kkk+i
                 isi=-1
              endif
              ikq=ik+ik+ik-3
C       
C       HAMILTONIEN
C       
              val=0.d0 
              der=0.D0
              RIK=R(IK)
              IF (RIK.LE.RMAX) THEN
                 if (rik.ge.rmin) then
                    call frss(rik,val,der)
                 else
                    val=1.d0
                    der=0.d0
                 endif
                 hii=hii+val
              ENDIF
C       
C       gradient
C       
              IDER=0
              DO 2584 l=1,3
                 ikql=ikq+l
                 derq=isi*der*unit(ikql)
                 HDER(K+IDER,IJ)=derq
                 HDER(I+IDER,IJ)=HDER(I+IDER,IJ)-DERQ
                 IDER=IDER+NBAT
 2584              CONTINUE
 221              kkk=kkk+k
 223           CONTINUE
           h(ij)=hii
           III=III+I
           ijq=ijq+3
              
 30        CONTINUE

C--------------------------------
C       DIAGONALISATION DE H
C---------------------------------
c       
c       nombre d orbitales occupees
c       
           NDOC=NBLEC/2
           NOC=NDOC
           MONOC=.FALSE.
           IF ((2*NDOC).NE.NBLEC) THEN
              MONOC=.TRUE. 
              NOC=NDOC+1
           ENDIF
c
           do i=1,nalk
              num(i)=i
           enddo
c       
c       sauvegarde de l hamiltonien dans ha
c       
           do ij=1,nalk*(nalk+1)/2
              ha(ij)=h(ij)
           enddo
           if(nalk.eq.1) then
              valp(1)=-h(1)
              vect(1)=1.D0
           endif
           if(nalk.gt.1) then
c       
c       appel a givens
c       
              if(idiag.eq.1) then
c       
c       ha contient la demi matrice inferieure
c       h contient la demi matrice superieure
c       h est detruit apres la diagonalisation
c       
                 ii=0
                 ij=0
                 do i=1,nalk
                    jj=ii
                    do j=i,nalk
                       ij=ij+1
                       h(ij)=ha(jj+i)
                       jj=jj+j
                    enddo
                    ii=ii+i
                 enddo

C start new DJW

                 HMAT(1:NBAT,1:NBAT)=0.0D0
                 ii=0
                 DO I=1,NALK
                    DO J=I,NALK
                       II=II+1
                       HMAT(I,J)=h(II)
                    ENDDO
                 ENDDO
                 
                 ABSTOL=DLAMCH('Safe  minimum')
                 CALL DSYEVR('V','I','U',NALK,HMAT,NBAT,0.0D0,1.0D0,1,NOC,ABSTOL,NFOUND,VALP,
     1                        VT,NBAT,ISUPPZ,WORK,LWORK, IWORK, ILWORK, INFO )
                 IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'

C end new DJW
                 
                 DO 210 I=1,noc
                    valp(i)=-valp(i)
 210             CONTINUE
              endif
           endif
c       
           DO I=1,noc
              valp(i)=-valp(i)
              ii=(num(i)-1)*nalk
C             DO J=1,nalk           ! DJW comment
C                VT(J,I)=VECT(II+J) ! DJW comment
C             END DO                ! DJW comment
           END DO
C       
C-----------------------------------------------------
C     ENERGIE TOTALE
C-----------------------------------------------------
           vgr=0.
           en=vgr
           DO 401 I=1,NDOC
              EN=EN+VALP(I)*2.0
 401           CONTINUE
           IF (monoc) en=en+valp(noc)
c       
C---------------------------------------------
C       DENSITE ELECTRONIQUE
C----------------------------------------------
           IND=0
           E0=EN
C-------------------------------------
C       CALCUL DU gradient DE L ENERGIE
C-------------------------------------
C       
C       A LA FIN grad CONTIENT LES DERIVEES DANS L ORDRE
C       (DE/DXI,I=1,NBAT),(DE/DYI,I=1,NBAT),(DE/DZI,I=1,NBAT)
C       
           
           if(nalk.gt.0) then
              DO 2570 K=1,3*NBAT
                 GRD=0.D0
                 DO 2572 I=1,NDOC
                    DO 2585 L=1,NALK
                       DO 2585 M=1,NALK
                          IF(L.GT.M) THEN
                             ML=L*(L-1)/2+M
                          ELSE
                             ML=M*(M-1)/2+L
                          ENDIF
                          GRD=GRD+2.D0*VT(L,I)*VT(M,I)*HDER(K,ML)
 2585            CONTINUE
 2572              CONTINUE
              if(monoc) then 
                 DO 2586 L=1,NALK
                    DO 2586 M=1,NALK
                       IF(L.GT.M) THEN
                          ML=L*(L-1)/2+M
                       ELSE
                          ML=M*(M-1)/2+L
                       ENDIF
 2586                 GRD=GRD+VT(L,NOC)*VT(M,NOC)*HDER(K,ML)
              endif
 2570              grad(K)=GRD
           endif
           
           ij=0
           erep=0.0
           do i=1,nbat
              do j=1,i-1
                 ij=ij+1
                 erep=erep+1./(r(ij)**12)
              enddo
              ij=ij+1
           enddo
           
           erep=erepmax*erep
           
           do i=1,nbat
              DP(3*i-2)=grad(i)
              DP(3*i-1)=grad(i+nbat)
              DP(3*i  )=grad(i+2*nbat)
           enddo
           
           do i=1,nbat
              do j=1,nbat
                 if (i.ne.j) then
                    xij=P(3*i-2)-P(3*j-2)
                    yij=P(3*i-1)-P(3*j-1)
                    zij=P(3*i)-P(3*j)
                    rij2=xij**2+yij**2+zij**2
                    rij14=rij2**7
                    DP(3*i-2)=DP(3*i-2)-erepmax*12.*xij/rij14
                    DP(3*i-1)=DP(3*i-1)-erepmax*12.*yij/rij14
                    DP(3*i  )=DP(3*i  )-erepmax*12.*zij/rij14
                 endif
              enddo
           enddo
           
           e0=e0+erep
           
        return
        end

c-----------------------------------------------------
        subroutine ftss(r,f,df)
        implicit DOUBLE PRECISION(a-h,o-z)
        IMPLICIT INTEGER(I-N)
        dimension a(10),b(10),c(10)
        character*2 metal
        common/METAL/metal
        if (metal.eq.'NA') then
           a(1)=-0.0000600345
           a(2)=0.7366633997
           a(3)=-21.2536277450
           b(1)=7.8227026687
           b(2)=5.0784652072
           b(3)=-5.7014389921
           c(1)=1.4334140314
           c(2)=2.7260385637
           c(3)=0.D0
           nt=3
           f=0.D0
           df=0.D0
           do i=1,nt
              fi=a(i)*(r**b(i))*dexp(-r*c(i))
              dfi=(-c(i)+b(i)/r)*fi
              f=f+fi
              df=df+dfi
           enddo
        elseif (metal.eq.'AG') then
           a(1)=-.0148920093/27.21
           b(1)= 9.4301889740
           c(1)=2.2599701079
           a(2)=-11.8760873508/27.21
           b(2)=-1.0347495105
           c(2)= .5509616893
           nt=2
           f=0.D0
           df=0.D0
           do i=1,nt
              fi=a(i)*(r**b(i))*dexp(-r*c(i))
              dfi=(-c(i)+b(i)/r)*fi
              f=f+fi
              df=df+dfi
           enddo
        elseif (metal.eq.'LI') then
           aa=-0.000195d0
           bb=1.65d0
           xn=8.d0
           fexp=aa*dexp(-bb*r)
           f=(r**xn)*fexp
           df=(xn*r**(xn-1.D0)-bb*(r**xn))*fexp
        endif
        return
        end
c-----------------------------------------------------
        subroutine ftsz(r,f,df)
        implicit DOUBLE PRECISION(a-h,o-z)
        IMPLICIT INTEGER(I-N)
        dimension a(10),b(10),c(10)
        character*2 metal
        common/METAL/metal
        if (metal.eq.'NA') then
           a(1)=0.0002818237
           a(2)=-0.0025316896
           a(3)=67.9894048289
           b(1)=9.2162096989
           b(2)= 4.0958902363
           b(3)= -5.2076863938
           c(1)= 2.2554055362
           c(2)=.9408906266
           c(3)=.6157151927
           nt=3
           f=0.D0
           df=0.D0
           do i=1,nt
              fi=a(i)*(r**b(i))*dexp(-r*c(i))
              dfi=(-c(i)+b(i)/r)*fi
              f=f+fi
              df=df+dfi
           enddo
        elseif (metal.eq.'AG') then
           nt=2
           a(1)= -.0058542/27.21
           b(1)= 9.5986911821
           c(1)= 2.0836237596
           a(2)= -17.9409106/27.219
           b(2)=-3.5128659875
           c(2)=-.0000000007
           a(1)=-.0043951845/27.21
           a(2)=-13.469637022/27.219
           f=0.D0
           df=0.D0
           do i=1,nt
              fi=a(i)*(r**b(i))*dexp(-r*c(i))
              dfi=(-c(i)+b(i)/r)*fi
              f=f+fi
              df=df+dfi
           enddo
        elseif (metal.eq.'LI') then
           aa=0.0000091d0
           bb=1.80d0
           xn=10.D0
           fexp=aa*dexp(-bb*r)
           f=(r**xn)*fexp
           df=(xn*r**(xn-1.D0)-bb*(r**xn))*fexp
        endif
        
        return
        end
c-----------------------------------
        subroutine frss(r,f,df)
        implicit DOUBLE PRECISION(a-h,o-z)
        character*2 metal
        common/METAL/metal
        if (metal.eq.'NA') then
           a=1.6822175236
           c= 1.3007271463
           re=5.5425217957
           c6= 10.4276970395
           x6= 5.8999900152
           frep=a*dexp(-c*r)
           drep=-c*frep
           fvdw=c6*r**(-x6)
           dvdw=-x6*fvdw/r
           arg= (re/r-1.d0)
           if(r.gt.(re)) then
              coup=1.d0
              dcoup=0.D0
           else
              coup=dexp(-arg**2)
              darg=-re/r**2
              dcoup=-2.D0*coup*arg*darg
           endif
           f=frep-coup*fvdw
           df=drep-dvdw*coup-dcoup*fvdw
        elseif (metal.eq.'AG') then
           a=8245.404377/27.21
           c=2.292537
           re=16.726372
           c6=10679.191179/27.21
           x6=6. 
           frep=a*dexp(-c*r)
           drep=-c*frep
           fvdw=c6*r**(-x6)
           dvdw=-x6*fvdw/r
           arg= (re/r-1.d0)
           if(r.gt.(re)) then
              coup=1.d0
              dcoup=0.D0
           else
              coup=dexp(-arg**2)
              darg=-re/r**2
              dcoup=-2.D0*coup*arg*darg
           endif
           f=frep-coup*fvdw
           df=drep-dvdw*coup-dcoup*fvdw
        elseif (metal.eq.'LI') then
           a=1.6822175236d0
           c= 1.3007271463d0*5.80d0/5.05d0
           re=5.5425217957d0*5.05d0/5.80d0
           x6= 5.8999900152d0
           c6= 10.4276970395d0*(5.05d0/5.80d0)**x6
           frep=a*dexp(-c*r)
           drep=-c*frep
           fvdw=c6*r**(-x6)
           dvdw=-x6*fvdw/r
           arg= (re/r-1.d0)
           if(r.gt.re) then
              coup=1.d0
              dcoup=0.d0
           else
              coup=dexp(-arg**2)
              darg=-re/r**2
              dcoup=-2.d0*coup*arg*darg
           endif
           f=1.43d0*(frep-coup*fvdw)
           df=1.43d0*(drep-dvdw*coup-dcoup*fvdw)
        endif
        
        return
        end
