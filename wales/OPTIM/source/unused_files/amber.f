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
      SUBROUTINE AMB(XVEC,GRAD,EREAL,GRADT,SECT)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      USE MODHESS
      IMPLICIT NONE
      DOUBLE PRECISION XVEC(3*NATOMS), GRAD(3*NATOMS), EREAL
      LOGICAL GRADT, SECT
      INTEGER J1, J2
      LOGICAL CONNECTIVITY_KNOWN
      SAVE


      IF (ang.EQ.0) THEN
        count=1
        CALL AMBERMASS
        CALL AMBERS
      END IF

C     IF (count.EQ.1) THEN
C       DO J1=1,NATOMS
C         XVEC(3*(J1-1)+1)=x(J1)
C         XVEC(3*(J1-1)+2)=y(J1)
C         XVEC(3*(J1-1)+3)=z(J1)
C       ENDDO
C     ELSE
        DO J1=1,NATOMS
           x(J1)=XVEC(3*(J1-1)+1)
           y(J1)=XVEC(3*(J1-1)+2)
           z(J1)=XVEC(3*(J1-1)+3)
        ENDDO
C     END IF
      count=0

      CALL amberenergy

      EREAL=totenergy

      IF (GRADT) THEN
         CALL AMBG
         DO J1=1,NATOMS
            GRAD(3*(J1-1)+1)=dEbydx(J1)
            GRAD(3*(J1-1)+2)=dEbydy(J1)
            GRAD(3*(J1-1)+3)=dEbydz(J1)
         ENDDO
      ENDIF

      IF (SECT) THEN
         CALL secondderivs
         DO J1=1,3*NATOMS
            DO J2=1,3*NATOMS
               HESS(J2,J1)=HELL(J2,J1)
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

      SUBROUTINE one4
      USE MODAMBER
      USE MODAMBER2
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2

      DO J1=1,NATOMS
         DO J2=1,NATOMS
            one_four(J2,J1)=0
         ENDDO
      ENDDO
      DO a=1,t
        one_four(da1(a),da4(a))=1
      ENDDO

C     PRINT *,"One-four relationships"
C     DO a=1,atoms
C       DO b=a+1,atoms
C         IF (one_four(a,b).EQ.1)  PRINT *,a,"-",b
C       ENDDO
C     ENDDO

      RETURN
      END
C
C  Setup stuff only needs doing once.
C
      SUBROUTINE AMBERS
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER            marvin
      INTEGER            chiratom(4),CANINE
      COMMON /CHIR/      canine
      LOGICAL CONNECTIVITY_KNOWN

      NDIHEDRALS=0
      INQUIRE (FILE='amber_connectivity',EXIST=CONNECTIVITY_KNOWN)

      CONNECTIVITY_KNOWN=.FALSE.

      IF (CONNECTIVITY_KNOWN) THEN
         PRINT *,'Connectivity previously determined - reading file'
         STOP
      ELSE
C 
C Create entry in angarray
C 
         DO b=1,atoms
            DO c=1,atoms
               DO d=b,atoms
                  IF (b.NE.c .AND. c.NE.d .AND. b.NE.d) THEN
                     IF (.NOT.(bonds(b,c).NE.1 .OR. bonds(c,d).NE.1)) THEN
                        ang=ang+1
                        aa1(ang)=b
                        aa2(ang)=c
                        aa3(ang)=d
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         DO i=1,atoms
            DO j=1,atoms
               DO k=1,atoms
                  DO l=i,atoms
                     IF (bonds(i,j).NE.1 .OR. bonds(j,i).NE.1) GOTO 113
                     IF (bonds(j,k).NE.1 .OR. bonds(k,j).NE.1) GOTO 113
                     IF (bonds(k,l).NE.1 .OR. bonds(l,k).NE.1) GOTO 113
                     IF (i.EQ.j .OR. i.EQ.k .OR. i.EQ.l .OR. j.EQ.k .OR. j.EQ.l .OR. k.EQ.l) GOTO 113
     
                     colin=0
                     ax=x(j)-x(i)
                     ay=y(j)-y(i)
                     az=z(j)-z(i)
                     bx=x(k)-x(j)
                     by=y(k)-y(j)
                     bz=z(k)-z(j)
                     cx=x(j)-x(l)
                     cy=y(j)-y(l)
                     cz=z(j)-z(l)
                     dx=x(l)-x(k)
                     dy=y(l)-y(k)
                     dz=z(l)-z(k)
C 
C Checking for colinearity,connectivity and coincidence
C 
                     IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) 
     1                  .AND. ay*bz.GT.(by*az-TINY) 
     2                  .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
                     IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) 
     1                  .AND. by*dz.GT.(dy*bz-TINY) 
     2                  .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2
                     IF (colin.EQ.1) THEN
                       PRINT *,'Three sites',i,j,k,'are colinear'
                     ELSE IF (colin.EQ.2) THEN
                       PRINT *,'Three sites',j,k,l,'are colinear'
                     ELSE
C 
C Create entry in torsarray
C
                       t=t+1
                       da1(t)=i
                       da2(t)=j
                       da3(t)=k
                       da4(t)=l
                     ENDIF
113                  CONTINUE
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL one4

         PRINT *,' Torsions assigned'

         DO a=1,t
            b=type(da1(a))
            c=type(da2(a))
            d=type(da3(a))        
            e=type(da4(a))
            match=0
            DO f=1,60
               IF (gentorsparams(f,1).LT.0.9) GOTO 5
               g=gentorsparams(f,1)
               h=gentorsparams(f,2)
               IF (c.EQ.g .AND. d.EQ.h) THEN
                  match=1
                  did(a)=gentorsparams(f,3)
                  dvn(a)=gentorsparams(f,4)
                  ddelta(a)=gentorsparams(f,5)
                  dn(a)=gentorsparams(f,6)
                  dvn2(a)=0.0
                  dvn3(a)=0.0
                  ddelta2(a)=0.0
                  ddelta3(a)=0.0
                  dn2(a)=0.0
                  dn3(a)=0.0
               ELSE IF (d.EQ.g .AND. c.EQ.h) THEN
                  match=1
                  did(a)=gentorsparams(f,3)
                  dvn(a)=gentorsparams(f,4)
                  ddelta(a)=gentorsparams(f,5)
                  dn(a)=gentorsparams(f,6)
                  dvn2(a)=0.0
                  dvn3(a)=0.0
                  ddelta2(a)=0.0
                  ddelta3(a)=0.0
                  dn2(a)=0.0
                  dn3(a)=0.0
               ENDIF
            ENDDO
5           CONTINUE

            DO f=1,50
               IF (spectorsparams(f,1).LT.0.9) GOTO 6
               g=spectorsparams(f,1)
               h=spectorsparams(f,2)
               i=spectorsparams(f,3)
               j=spectorsparams(f,4)
               IF (b.EQ.g .AND. c.EQ.h .AND. d.EQ.i .AND. e.EQ.j) THEN
                  match=1
                  did(a)=spectorsparams(f,5)
                  dvn(a)=spectorsparams(f,6)
                  ddelta(a)=spectorsparams(f,7)
                  dn(a)=spectorsparams(f,8)
                  dvn2(a)=spectorsparams(f,9)
                  ddelta2(a)=spectorsparams(f,10)
                  dn2(a)=spectorsparams(f,11)
                  dvn3(a)=spectorsparams(f,12)
                  ddelta3(a)=spectorsparams(f,13)
                  dn3(a)=spectorsparams(f,14)
               ELSE IF (e.EQ.g .AND. d.EQ.h .AND. c.EQ.i .AND. b.EQ.j) THEN
                  match=1
                  did(a)=spectorsparams(f,5)
                  dvn(a)=spectorsparams(f,6)
                  ddelta(a)=spectorsparams(f,7)
                  dn(a)=spectorsparams(f,8)
                  dvn2(a)=spectorsparams(f,9)
                  ddelta2(a)=spectorsparams(f,10)
                  dn2(a)=spectorsparams(f,11)
                  dvn3(a)=spectorsparams(f,12)
                  ddelta3(a)=spectorsparams(f,13)
                  dn3(a)=spectorsparams(f,14)
               ENDIF
            ENDDO
6           CONTINUE
         ENDDO

         IF (match.EQ.0) THEN
            did(a)=1
            dvn(a)=0.0
            dvn2(a)=0.0
            dvn3(a)=0.0
         ENDIF

         DO i=1,atoms
          DO j=i+1,atoms
           DO k=j+1,atoms
            DO l=k+1,atoms
             IF (.NOT.(i.EQ.j .OR. i.EQ.k .OR. i.EQ.l .OR. j.EQ.k .OR. j.EQ.l .OR. k.EQ.l)) THEN
C 
C This IF construct puts the atoms in the improper in the correct order 1-2-3-4
C 
             IF (bonds(i,k).EQ.1 .AND. bonds(j,k).EQ.1 .AND. bonds(k,l).EQ.1) THEN
              imp=imp+1
              IF ((type(i) .LE. type(j)) .AND. (type(j) .LE. type(l))) THEN 
               ia1(imp)=i
               ia2(imp)=j
               ia3(imp)=k
               ia4(imp)=l
              ELSE IF ((type(i) .LE. type(l)) .AND. (type(l) .LT. type(j))) THEN
               ia1(imp)=i
               ia2(imp)=l
               ia3(imp)=k
               ia4(imp)=j
              ELSE IF ((type(j) .LT. type(i)) .AND. (type(i) .LE. type(l))) THEN
               ia1(imp)=j
               ia2(imp)=i
               ia3(imp)=k
               ia4(imp)=l
              ELSE IF ((type(j) .LE. type(l)) .AND. (type(l) .LT. type(i))) THEN
               ia1(imp)=j
               ia2(imp)=l
               ia3(imp)=k
               ia4(imp)=i
              ELSE IF ((type(l) .LT. type(i)) .AND. (type(i) .LE. type(j))) THEN
               ia1(imp)=l
               ia2(imp)=i
               ia3(imp)=k
               ia4(imp)=j
              ELSE IF ((type(l) .LT. type(j)) .AND. (type(j) .LT. type(i))) THEN
               ia1(imp)=l
               ia2(imp)=j
               ia3(imp)=k
               ia4(imp)=i
              ENDIF
            
             ELSE IF (bonds(i,j).EQ.1 .AND. bonds(i,k).EQ.1 .AND. bonds(i,l).EQ.1) THEN
              imp=imp+1
              IF ((type(j) .LE. type(k)) .AND. (type(k) .LE. type(l))) THEN 
               ia1(imp)=j
               ia2(imp)=k
               ia3(imp)=i
               ia4(imp)=l
              ELSE IF ((type(j) .LE. type(l)) .AND. (type(l) .LT. type(k))) THEN
               ia1(imp)=j
               ia2(imp)=l
               ia3(imp)=i
               ia4(imp)=k
              ELSE IF ((type(k) .LT. type(j)) .AND. (type(j) .LE. type(l))) THEN
               ia1(imp)=k
               ia2(imp)=j
               ia3(imp)=i
               ia4(imp)=l
              ELSE IF ((type(k) .LE. type(l)) .AND. (type(l) .LT. type(j))) THEN
               ia1(imp)=k
               ia2(imp)=l
               ia3(imp)=i
               ia4(imp)=j
              ELSE IF ((type(l) .LT. type(j)) .AND. (type(j) .LE. type(k))) THEN
               ia1(imp)=l
               ia2(imp)=j
               ia3(imp)=i
               ia4(imp)=k
              ELSE IF ((type(l) .LT. type(k)) .AND. (type(k) .LT. type(j))) THEN
               ia1(imp)=l
               ia2(imp)=k
               ia3(imp)=i
               ia4(imp)=j
              ENDIF
   
             ELSE IF (bonds(i,j).EQ.1 .AND. bonds(j,k).EQ.1 .AND. bonds(j,l).EQ.1) THEN
              imp=imp+1
              IF ((type(i) .LE. type(k)) .AND. (type(k) .LE. type(l))) THEN 
               ia1(imp)=i
               ia2(imp)=k
               ia3(imp)=j
               ia4(imp)=l
              ELSE IF ((type(i) .LE. type(l)) .AND. (type(l) .LT. type(k))) THEN
               ia1(imp)=i
               ia2(imp)=l
               ia3(imp)=j
               ia4(imp)=k
              ELSE IF ((type(k) .LT. type(i)) .AND. (type(i) .LE. type(l))) THEN
               ia1(imp)=k
               ia2(imp)=i
               ia3(imp)=j
               ia4(imp)=l
              ELSE IF ((type(k) .LE. type(l)) .AND. (type(l) .LT. type(i))) THEN
               ia1(imp)=k
               ia2(imp)=l
               ia3(imp)=j
               ia4(imp)=i
              ELSE IF ((type(l) .LT. type(i)) .AND. (type(i) .LE. type(k))) THEN
               ia1(imp)=l
               ia2(imp)=i
               ia3(imp)=j
               ia4(imp)=k
              ELSE IF ((type(l) .LT. type(k)) .AND. (type(k) .LT. type(i))) THEN
               ia1(imp)=l
               ia2(imp)=k
               ia3(imp)=j
               ia4(imp)=i
              ENDIF
   
             ELSE IF (bonds(i,l).EQ.1 .AND. bonds(j,l).EQ.1 .AND. bonds(k,l).EQ.1) THEN
              imp=imp+1
              IF ((type(i) .LE. type(j)) .AND. (type(j) .LE. type(k))) THEN 
               ia1(imp)=i
               ia2(imp)=j
               ia3(imp)=l
               ia4(imp)=k
              ELSE IF ((type(i) .LE. type(k)) .AND. (type(k) .LT. type(j))) THEN
               ia1(imp)=i
               ia2(imp)=k
               ia3(imp)=l
               ia4(imp)=j
              ELSE IF ((type(j) .LT. type(i)) .AND. (type(i) .LE. type(k))) THEN
               ia1(imp)=j
               ia2(imp)=i
               ia3(imp)=l
               ia4(imp)=k
              ELSE IF ((type(j) .LE. type(k)) .AND. (type(k) .LT. type(i))) THEN
               ia1(imp)=j
               ia2(imp)=k
               ia3(imp)=l
               ia4(imp)=i
              ELSE IF ((type(k) .LT. type(i)) .AND. (type(i) .LE. type(j))) THEN
               ia1(imp)=k
               ia2(imp)=i
               ia3(imp)=l
               ia4(imp)=j
              ELSE IF ((type(k) .LT. type(j)) .AND. (type(j) .LT. type(i))) THEN
               ia1(imp)=k
               ia2(imp)=j
               ia3(imp)=l
               ia4(imp)=i
              ENDIF
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO

         PRINT *,' Impropers assigned'

C
C Assign parameters to each improper
C
         DO a=1,imp
            b=type(ia1(a))
            c=type(ia2(a))        
            d=type(ia3(a))
            e=type(ia4(a))
            DO f=1,15
              g=genimpparams(f,1)
              h=genimpparams(f,2)
              IF (d.EQ.g .AND. (e.EQ.h .OR. c.EQ.h .OR. b.EQ.h)) THEN
                ivn(a)=genimpparams(f,3)
                idelta(a)=genimpparams(f,4)
                in1(a)=genimpparams(f,5)
              ENDIF
            ENDDO
            DO f=1,4
              g=midimpparams(f,1)
              h=midimpparams(f,2)
              i=midimpparams(f,3)
              IF (c.EQ.g .AND. d.EQ.h .AND. e.EQ.i) THEN
                ivn(a)=midimpparams(f,4)
                idelta(a)=midimpparams(f,5)
                in1(a)=midimpparams(f,6)
              ELSE IF (b.EQ.g .AND. c.EQ.i .AND. d.EQ.h) THEN
                ivn(a)=midimpparams(f,4)
                idelta(a)=midimpparams(f,5)
                in1(a)=midimpparams(f,6)
              ENDIF
            ENDDO
            DO f=1,15
              g=specimpparams(f,1)
              h=specimpparams(f,2)
              i=specimpparams(f,3)
              j=specimpparams(f,4)
      
              IF (b.EQ.g .AND. c.EQ.h .AND. d.EQ.i .AND. e.EQ.j) THEN
                ivn(a)=specimpparams(f,5)
                idelta(a)=specimpparams(f,6)
                in1(a)=specimpparams(f,7)
              ENDIF
            ENDDO
         ENDDO
      END IF

      DO a=1,atoms
         DO b=a+1,atoms
            rstar=vdwr(type(a))+vdwr(type(b))
            epsilon=SQRT(vdwe(type(a))*vdwe(type(b)))
            vdwa(a,b)=epsilon*rstar**12
            vdwb(a,b)=2*epsilon*rstar**6  
         ENDDO
      ENDDO

      canine=1
      DO a=1,atoms
        IF (chiral(a).EQ.1) GOTO 111

        IF (typech(a).NE.'CT') GOTO 112
        marvin=1
        DO b=1,atoms
          IF (bonds(a,b).EQ.1 .OR. bonds(b,a).EQ.1) THEN
            chiratom(marvin)=b
            marvin=marvin+1
          END IF
        END DO

        i=chiratom(1)
        j=chiratom(2)
        k=chiratom(3)
        l=chiratom(4)

        IF (typech(i).NE.typech(j) .AND. typech(i).NE.typech(k) .AND. typech(i).NE.typech(l) .AND. 
     1      typech(j).NE.typech(k) .AND. typech(j).NE.typech(l) .AND. typech(k).NE.typech(l)) THEN
          chiral(a)=1
        ELSE
          GOTO 112
        END IF
111     CONTINUE

        chiralarray(canine,1)=a
        chiralarray(canine,2)=i
        chiralarray(canine,3)=j
        chiralarray(canine,4)=k
        chiralarray(canine,5)=l
        canine=canine+1
   
112   CONTINUE
      END DO

      PRINT *,'Set up routine completed'

      canine=canine-1
      DO a=1,canine
        i=chiralarray(a,2)
        j=chiralarray(a,3)
        k=chiralarray(a,4)
        l=chiralarray(a,5)
C        PRINT *,'Atom',chiralarray(a,1),' is chiral and is bonded to atoms',i,' ',j,' ',k,' ',l
      END DO

      DO a=1,atoms
        IF (label(a).EQ.'d' .AND. typech(a).NE.'N3') THEN
          DO b=a+1,atoms
            IF (label(b).EQ.'d' .AND. bonds(a,b).EQ.1) THEN
              NDIHEDRALS=NDIHEDRALS+1
              DATOM1(NDIHEDRALS)=a
              DATOM2(NDIHEDRALS)=b
C              PRINT *,'Dihedral',NDIHEDRALS,'  ',a,'-',b
C              PRINT *,'         ',typech(a),' - ',typech(b)
            END IF
          END DO
        END IF
      END DO

      CALL aconnectdump


      RETURN
      END
C
C  Distances
C  Van der Waals and charge terms
C
      SUBROUTINE AMBERD
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER NDUMMY1, NDUMMY3
      DOUBLE PRECISION DUMMY, DUMMY2, VIXEN1, VIXEN2, VIXEN3, VIXEN4, VIXEN5

      DO a=1,atoms
         vixen4=pq(a)
         DO b=a+1,atoms
            DUMMY2=(x(a)-x(b))**2+(y(a)-y(b))**2+(z(a)-z(b))**2
            vixen5=DSQRT(DUMMY2)
            r(a,b)=vixen5
            r(b,a)=vixen5

            NDUMMY1=1-bonds(a,b)
            NDUMMY3=NDUMMY1*(1-one_three(a,b))
            IF (NDUMMY3.NE.0) THEN
               vixen3=DBLE(one_four(a,b))

               DUMMY=1.0D0/DUMMY2**3
               vixen1=(vdwa(a,b)*DUMMY-vdwb(a,b))*DUMMY
               IF (FAKEWATER) THEN
                  vixen2=vixen4*pq(b)/(dielec*vixen5**2)
               ELSE
                  vixen2=vixen4*pq(b)/(vixen5*dielec)
               END IF

               vdwenergy=vdwenergy+vixen1*(1.0D0-vixen3/2.0D0)
               qenergy=    qenergy+vixen2*(1.0D0-vixen3/6.0D0)
            ENDIF

C           vixen3=FLOAT(one_four(a,b))

C           DUMMY=1.0D0/DUMMY2**3
C           vixen1=(vdwa(a,b)*DUMMY-vdwb(a,b))*DUMMY
C           vixen2=vixen4*pq(b)/r(a,b)

C           vdwenergy=vdwenergy+vixen1*DUMMY3*(1.0D0-vixen3/2.0D0)
C           qenergy=    qenergy+vixen2*DUMMY3*(1.0D0-vixen3/6.0D0)

         ENDDO
      ENDDO
      qenergy=qenergy/dielec

      RETURN
      END


      SUBROUTINE chiraltest(CTEST)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE

      INTEGER  chiratom(4)
      INTEGER  marvin
      DOUBLE PRECISION  TEMPXH,TEMPYH,TEMPZH,TEMPXI,TEMPYI,TEMPZI,TEMPXJ,TEMPYJ,TEMPZJ
      DOUBLE PRECISION  TEMPXK,TEMPYK,TEMPZK,TEMPXL,TEMPYL,TEMPZL
      DOUBLE PRECISION  NEWXJ,NEWYJ,NEWZJ,NEWXK,NEWYK,NEWZK,GAMMA,DELTA,DOT
      DOUBLE PRECISION  CROSS(3)
      LOGICAL  CTEST
      INTEGER canine
      COMMON /CHIR/ canine

      CTEST=.FALSE.

      DO a=1,atoms
        IF (chiral(a).NE.1) GOTO 17
        marvin=1
        DO b=1,atoms
          IF (bonds(a,b).EQ.1 .OR. bonds(b,a).EQ.1) THEN
            chiratom(marvin)=b
            marvin=marvin+1
          END IF
        END DO    

        i=chiratom(1)
        j=chiratom(2)
        k=chiratom(3)
        l=chiratom(4)

        IF (typech(i).EQ.'HC' .OR. typech(j).EQ.'HC' .OR. typech(k).EQ.'HC' .OR. typech(l).EQ.'HC') THEN
C
C  We are dealing with the beta carbon of isoleucine
C
        ELSE IF (typech(i).EQ.'OH'.OR.typech(j).EQ.'OH'.OR.typech(k).EQ.'OH'.OR.typech(l).EQ.'OH') THEN
C
C  We are dealing with the beta carbon of threonine
C  Note that we order the atoms backwards as this carbon should be 'R' configuration
C
          PRINT *,' Found threonine beta carbon'
          IF (typech(i).EQ.'OH') THEN
            chiratom(4)=i
          ELSE IF (typech(j).EQ.'OH') THEN
            chiratom(4)=j
          ELSE IF (typech(k).EQ.'OH') THEN
            chiratom(4)=k
          ELSE
            chiratom(4)=l
          END IF
          IF (typech(i).EQ.'H1') THEN
            chiratom(1)=i
          ELSE IF (typech(j).EQ.'H1') THEN
            chiratom(1)=j
          ELSE IF (typech(k).EQ.'H1') THEN
            chiratom(1)=k
          ELSE
            chiratom(1)=l
          END IF
          IF (label(i).EQ.'d') THEN
            chiratom(3)=i
          ELSE IF (label(j).EQ.'d') THEN
            chiratom(3)=j
          ELSE IF (label(k).EQ.'d') THEN
            chiratom(3)=k
          ELSE 
            chiratom(3)=l
          END IF
          IF ((chiratom(1).EQ.i .OR. chiratom(3).EQ.i .OR. chiratom(4).EQ.i) .AND. (chiratom(1).EQ.j
     1       .OR. chiratom(3).EQ.j .OR. chiratom(4).EQ.j) .AND. (chiratom(1).EQ.k .OR. chiratom(3).EQ.k
     2       .OR. chiratom(4).EQ.k)) THEN
            chiratom(2)=l
          ELSE IF ((chiratom(1).EQ.i .OR. chiratom(3).EQ.i .OR. chiratom(4).EQ.i) .AND. (chiratom(1).EQ.j
     1       .OR. chiratom(3).EQ.j .OR. chiratom(4).EQ.j) .AND. (chiratom(1).EQ.l .OR. chiratom(3).EQ.l
     2       .OR. chiratom(4).EQ.l)) THEN          
            chiratom(2)=k
          ELSE IF ((chiratom(1).EQ.i .OR. chiratom(3).EQ.i .OR. chiratom(4).EQ.i) .AND. (chiratom(1).EQ.k
     1       .OR. chiratom(3).EQ.k .OR. chiratom(4).EQ.k) .AND. (chiratom(1).EQ.l .OR. chiratom(3).EQ.l
     2       .OR. chiratom(4).EQ.l)) THEN          
            chiratom(2)=j
          ELSE 
            chiratom(2)=i
          END IF

          PRINT *,'Ranking for threonine beta carbon'
          PRINT *,'1',typech(chiratom(4)),' ',chiratom(4)
          PRINT *,'2',typech(chiratom(3)),' ',chiratom(3)
          PRINT *,'3',typech(chiratom(2)),' ',chiratom(2)
          PRINT *,'4',typech(chiratom(1)),' ',chiratom(1)

        ELSE
C  We have an alpha carbon
          IF (typech(j)(1:1).EQ.'H') THEN
            chiratom(1)=j
          ELSE IF (typech(k)(1:1).EQ.'H') THEN
            chiratom(1)=k
          ELSE IF (typech(l)(1:1).EQ.'H') THEN
            chiratom(1)=l
          END IF
    
          IF ((typech(i).EQ.'N ').OR.(typech(i).EQ.'N3')) THEN
            chiratom(2)=i
          ELSE IF ((typech(k).EQ.'N ').OR.(typech(k).EQ.'N3')) THEN
            chiratom(2)=k
          ELSE IF ((typech(l).EQ.'N ').OR.(typech(l).EQ.'N3')) THEN
            chiratom(2)=l
          END IF
  
          IF (typech(i).EQ.'C ') THEN
            chiratom(3)=i
          ELSE IF (typech(j).EQ.'C ') THEN
            chiratom(3)=j
          ELSE IF (typech(l).EQ.'C ') THEN
            chiratom(3)=l
          END IF
  
          IF (typech(i).EQ.'CT') THEN
            chiratom(4)=i
          ELSE IF (typech(j).EQ.'CT') THEN
            chiratom(4)=j
          ELSE IF (typech(k).EQ.'CT') THEN
            chiratom(4)=k
          END IF
        END IF

          h=a
          i=chiratom(1)
          j=chiratom(2)
          k=chiratom(3)
          l=chiratom(4)
C
C Move group to origin
C
        tempxi=x(i)-x(h)
        tempyi=y(i)-y(h)
        tempzi=z(i)-z(h)
        tempxj=x(j)-x(h)
        tempyj=y(j)-y(h)
        tempzj=z(j)-z(h)
        tempxk=x(k)-x(h)
        tempyk=y(k)-y(h)
        tempzk=z(k)-z(h)
        tempxl=x(l)-x(h)
        tempyl=y(l)-y(h)
        tempzl=z(l)-z(h)
        tempxh=0.0
        tempyh=0.0
        tempzh=0.0
  
        gamma=-(tempxj*tempxi+tempyj*tempyi+tempzj*tempzi)/(tempxi**2+tempyi**2+tempzi**2)
        delta=-(tempxk*tempxi+tempyk*tempyi+tempzk*tempzi)/(tempxi**2+tempyi**2+tempzi**2)       

        newxj=tempxj+gamma*tempxi
        newyj=tempyj+gamma*tempyi
        newzj=tempzj+gamma*tempzi
        newxk=tempxk+delta*tempxi
        newyk=tempyk+delta*tempyi
        newzk=tempzk+delta*tempzi
  
        cross(1)=newyj*newzk-newzj*newyk
        cross(2)=newzj*newxk-newxj*newzk
        cross(3)=newxj*newyk-newyj*newxk
  
        dot=cross(1)*tempxi+cross(2)*tempyi+cross(3)*tempzi
        IF (dot .LT. 0.0) THEN
C          PRINT *,'Configuration is S'
        ELSE
          CTEST=.TRUE. 
        END IF
17    CONTINUE
      END DO
  
      RETURN
      END

      SUBROUTINE ambermass
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER m1

      DO m1=1,atoms
         IF (typech(m1).EQ.'C ') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CA') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CB') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CC') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CK') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CM') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CN') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CQ') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CR') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CT') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CV') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'CW') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'C*') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'C0') THEN
            mass(m1)= 12.01
         ELSE IF (typech(m1).EQ.'F ') THEN
            mass(m1)= 19.00
         ELSE IF (typech(m1).EQ.'H ') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HC') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H1') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H2') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H3') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HA') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H4') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'H5') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HO') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HS') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HW') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'HP') THEN
            mass(m1)= 1.008
         ELSE IF (typech(m1).EQ.'N ') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'NA') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'NB') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'NC') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'N2') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'N3') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'N*') THEN
            mass(m1)= 14.01
         ELSE IF (typech(m1).EQ.'O ') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'OW') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'OH') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'OS') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'O2') THEN
            mass(m1)= 16.00
         ELSE IF (typech(m1).EQ.'P ') THEN
            mass(m1)= 30.97
         ELSE IF (typech(m1).EQ.'S ') THEN
            mass(m1)= 32.06
         ELSE IF (typech(m1).EQ.'SH') THEN
            mass(m1)= 32.06
         END IF
      END DO


      RETURN

      END 



      SUBROUTINE AMBG
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER J1
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,DUBYDX,DUBYDY,DUBYDZ,DVBYDX,DVBYDY,DVBYDZ,U,V
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3,VIXEN4,VIXEN5,VIXEN6,VIXEN7,VIXEN8
C 
C Initialise all derivatives
C 
      DO J1=1,atoms
         dbondEbydx(J1)=0.0D0
         dbondEbydy(J1)=0.0D0
         dbondEbydz(J1)=0.0D0
         dangEbydx(J1)=0.0D0
         dangEbydy(J1)=0.0D0
         dangEbydz(J1)=0.0D0
         dtorsEbydx(J1)=0.0D0
         dtorsEbydy(J1)=0.0D0
         dtorsEbydz(J1)=0.0D0
         dvdwEbydx(J1)=0.0D0
         dvdwEbydy(J1)=0.0D0
         dvdwEbydz(J1)=0.0D0
         dqEbydx(J1)=0.0D0
         dqEbydy(J1)=0.0D0
         dqEbydz(J1)=0.0D0
      ENDDO

C      PRINT *,'Co-ordinates are:'
C      WRITE (*,FMT='(A,3F11.6)') 'O',x(1),y(1),z(1)
C      WRITE (*,FMT='(A,3F11.6)') 'H',x(2),y(2),z(2)
C      WRITE (*,FMT='(A,3F11.6)') 'H',x(3),y(3),z(3)

       DO a=1,atoms
         DO b=a+1,atoms
           IF (bonds(a,b).EQ.1) THEN
  
             dbondEbydx(a)=dbondEbydx(a)+(2.0D0*kr(type(a),type(b))*(x(a)-x(b))*(1.0D0-(ro(type(a),type(b))/r(a,b))))
             dbondEbydy(a)=dbondEbydy(a)+(2.0D0*kr(type(a),type(b))*(y(a)-y(b))*(1.0D0-(ro(type(a),type(b))/r(a,b))))
             dbondEbydz(a)=dbondEbydz(a)+(2.0D0*kr(type(a),type(b))*(z(a)-z(b))*(1.0D0-(ro(type(a),type(b))/r(a,b))))
             dbondEbydx(b)=dbondEbydx(b)+(2.0D0*kr(type(a),type(b))*(x(b)-x(a))*(1.0D0-(ro(type(a),type(b))/r(a,b))))
             dbondEbydy(b)=dbondEbydy(b)+(2.0D0*kr(type(a),type(b))*(y(b)-y(a))*(1.0D0-(ro(type(a),type(b))/r(a,b))))
             dbondEbydz(b)=dbondEbydz(b)+(2.0D0*kr(type(a),type(b))*(z(b)-z(a))*(1.0D0-(ro(type(a),type(b))/r(a,b))))

           ENDIF
         ENDDO
       ENDDO

C  PRINT *,"Bond derivatives calculated"
C 
C First calculate energy gradient from central atom in angle
C 
      DO d=1,ang    
       a=aa1(d)
        b=aa2(d)
        c=aa3(d)
 
       vixen1=(x(b)**2+y(b)**2+z(b)**2)-(x(b)*x(c)+y(b)*y(c)+z(b)*z(c))
       u=((x(c)-x(b))*x(a)+(y(c)-y(b))*y(a)+(z(c)-z(b))*z(a))+vixen1
       v=r(a,b)*r(b,c)
 
       dubydx=2.0D0*x(b)-x(c)-x(a)
       dubydy=2.0D0*y(b)-y(c)-y(a)
       dubydz=2.0D0*z(b)-z(c)-z(a)
       dvbydx=((x(b)-x(a))*r(b,c)**2+(x(b)-x(c))*r(a,b)**2)/(r(a,b)*r(b,c))
       dvbydy=((y(b)-y(a))*r(b,c)**2+(y(b)-y(c))*r(a,b)**2)/(r(a,b)*r(b,c))
       dvbydz=((z(b)-z(a))*r(b,c)**2+(z(b)-z(c))*r(a,b)**2)/(r(a,b)*r(b,c))
       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       e=type(a)
       f=type(b)
       g=type(c)

       dangEbydx(b)=dangEbydx(b)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
       dangEbydy(b)=dangEbydy(b)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
       dangEbydz(b)=dangEbydz(b)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))
C 
C Now from end atoms
C 
       dubydx=x(c)-x(b)
       dubydy=y(c)-y(b)
       dubydz=z(c)-z(b)
       dvbydx=(x(a)-x(b))*r(b,c)/r(a,b)
       dvbydy=(y(a)-y(b))*r(b,c)/r(a,b)
       dvbydz=(z(a)-z(b))*r(b,c)/r(a,b)
       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       dangEbydx(a)=dangEbydx(a)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
       dangEbydy(a)=dangEbydy(a)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
       dangEbydz(a)=dangEbydz(a)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))


       dubydx=x(a)-x(b)
       dubydy=y(a)-y(b)
       dubydz=z(a)-z(b)
       dvbydx=(x(c)-x(b))*r(a,b)/r(b,c)
       dvbydy=(y(c)-y(b))*r(a,b)/r(b,c)
       dvbydz=(z(c)-z(b))*r(a,b)/r(b,c)
       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       dangEbydx(c)=dangEbydx(c)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydx/SIN(theta(d))))
       dangEbydy(c)=dangEbydy(c)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydy/SIN(theta(d))))
       dangEbydz(c)=dangEbydz(c)+2.0D0*(kt(e,f,g)*(theta(d)-to(e,f,g))*(-dpbydz/SIN(theta(d))))

      ENDDO

C  PRINT *,"Angle derivatives calculated"

       DO a=1,atoms
         DO b=a+1,atoms
           IF (bonds(a,b).NE.1) THEN
             vixen1=(-12.0D0*vdwa(a,b)/r(a,b)**13)+(6.0D0*vdwb(a,b))/r(a,b)**7
             IF (FAKEWATER) THEN
               vixen2=2.0D0*pq(a)*pq(b)/(dielec*r(a,b)**3)
             ELSE
               vixen2=pq(a)*pq(b)/(dielec*r(a,b)**2)
             END IF
             vixen3=(x(a)-x(b))/r(a,b)
             vixen4=(y(a)-y(b))/r(a,b)
             vixen5=(z(a)-z(b))/r(a,b)
             vixen6=(x(b)-x(a))/r(a,b)
             vixen7=(y(b)-y(a))/r(a,b)
             vixen8=(z(b)-z(a))/r(a,b)
 
             IF ((one_four(a,b).EQ.1).AND.(one_three(a,b).NE.1)) THEN
               dvdwEbydx(a)=dvdwEbydx(a)+0.5D0*(vixen1*vixen3)
               dvdwEbydy(a)=dvdwEbydy(a)+0.5D0*(vixen1*vixen4)
               dvdwEbydz(a)=dvdwEbydz(a)+0.5D0*(vixen1*vixen5)
               dvdwEbydx(b)=dvdwEbydx(b)+0.5D0*(vixen1*vixen6)
               dvdwEbydy(b)=dvdwEbydy(b)+0.5D0*(vixen1*vixen7)
               dvdwEbydz(b)=dvdwEbydz(b)+0.5D0*(vixen1*vixen8)
               dqEbydx(a)=dqEbydx(a)-(1.0D0/1.2D0)*vixen2*vixen3
               dqEbydy(a)=dqEbydy(a)-(1.0D0/1.2D0)*vixen2*vixen4
               dqEbydz(a)=dqEbydz(a)-(1.0D0/1.2D0)*vixen2*vixen5
               dqEbydx(b)=dqEbydx(b)-(1.0D0/1.2D0)*vixen2*vixen6
               dqEbydy(b)=dqEbydy(b)-(1.0D0/1.2D0)*vixen2*vixen7
               dqEbydz(b)=dqEbydz(b)-(1.0D0/1.2D0)*vixen2*vixen8
             ELSE IF (one_three(a,b).EQ.1) THEN
C              DO NOTHING
             ELSE 
               dvdwEbydx(a)=dvdwEbydx(a)+(vixen1*vixen3)
               dvdwEbydy(a)=dvdwEbydy(a)+(vixen1*vixen4)
               dvdwEbydz(a)=dvdwEbydz(a)+(vixen1*vixen5)
               dvdwEbydx(b)=dvdwEbydx(b)+(vixen1*vixen6)
               dvdwEbydy(b)=dvdwEbydy(b)+(vixen1*vixen7)
               dvdwEbydz(b)=dvdwEbydz(b)+(vixen1*vixen8)
               dqEbydx(a)=dqEbydx(a)-(vixen2*vixen3)
               dqEbydy(a)=dqEbydy(a)-(vixen2*vixen4)
               dqEbydz(a)=dqEbydz(a)-(vixen2*vixen5)
               dqEbydx(b)=dqEbydx(b)-(vixen2*vixen6)
               dqEbydy(b)=dqEbydy(b)-(vixen2*vixen7)
               dqEbydz(b)=dqEbydz(b)-(vixen2*vixen8)
             END IF
           END IF
         ENDDO
       ENDDO


      CALL pderivs

      CALL impderivs
 
      DO a=1,atoms
       dEbydx(a)=dbondEbydx(a)+dangEbydx(a)+dvdwEbydx(a)+dtorsEbydx(a)+dimpEbydx(a)+dqEbydx(a)
       dEbydy(a)=dbondEbydy(a)+dangEbydy(a)+dvdwEbydy(a)+dtorsEbydy(a)+dimpEbydy(a)+dqEbydy(a)
       dEbydz(a)=dbondEbydz(a)+dangEbydz(a)+dvdwEbydz(a)+dtorsEbydz(a)+dimpEbydz(a)+dqEbydz(a)
      ENDDO

      IF (count.EQ.1) PRINT *,"Derivatives calculated"
      count=0

      RETURN
      END

      SUBROUTINE pderivs
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION XE,XF,XG,YE,YF,YG,ZE,ZF,ZG,U,V,UP,VP,P
      DOUBLE PRECISION DLAMBDABYDX,DLAMBDABYDY,DLAMBDABYDZ,DMUBYDX,DMUBYDY,DMUBYDZ
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,DUBYDX,DUBYDY,DUBYDZ,DVBYDX,DVBYDY,DVBYDZ
      DOUBLE PRECISION DPHIBYDX,DPHIBYDY,DPHIBYDZ,DUPBYDX,DUPBYDY,DUPBYDZ,DVPBYDX,DVPBYDY,DVPBYDZ
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3
      INTEGER J1
C 
C Initialise derivatives
C 
      DO J1=1,atoms
         dtorsEbydx(J1)=0.0D0
         dtorsEbydy(J1)=0.0D0
         dtorsEbydz(J1)=0.0D0
      ENDDO
C 
C Calculate derivatives for all atoms in position "a"
C 
      DO j=1,t
C 
C Initialise intermediates
C 
       dubydx=0.0D0
       dubydy=0.0D0
       dubydz=0.0D0
       dvbydx=0.0D0
       dvbydy=0.0D0
       dvbydz=0.0D0
       dpbydx=0.0D0
       dpbydy=0.0D0
       dupbydx=0.0D0
       dupbydy=0.0D0
       dupbydz=0.0D0
       dvpbydx=0.0D0
       dvpbydy=0.0D0
       dvpbydz=0.0D0
       dpbydz=0.0D0

       i=da1(j)
       IF (i.EQ.0) GOTO 20
       a=da1(j)
       b=da2(j)
       c=da3(j)
       d=da4(j)
C 
C First define all functions as on paper
C 
       xe=x(c)-x(b)
       ye=y(c)-y(b)
       ze=z(c)-z(b)
       lambda=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(a)*xe)+(y(a)*ye)+(z(a)*ze)))/(xe**2+ye**2+ze**2)
       mu=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(d)*xe)+(y(d)*ye)+(z(d)*ze)))/(xe**2+ye**2+ze**2)
       xf=x(b)-x(a)-(lambda*xe)
       yf=y(b)-y(a)-(lambda*ye)
       zf=z(b)-z(a)-(lambda*ze)
       xg=x(b)-x(d)-(mu*xe)
       yg=y(b)-y(d)-(mu*ye)
       zg=z(b)-z(d)-(mu*ze)
       vixen1=(x(b)**2+y(b)**2+z(b)**2)+mu*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)-((x(b)*x(d))+(y(b)*y(d))+(z(b)*z(d)))
       vixen3=(x(d)-x(b))*x(a)+(y(d)-y(b))*y(a)+(z(d)-z(b))*z(a)
       u=vixen1+vixen2+vixen3+(lambda*mu*r(b,c)**2)
       up=SQRT(xf**2+yf**2+zf**2)
       vp=SQRT(xg**2+yg**2+zg**2)
       v=up*vp
       p=u/v

       IF (p.LT.-1) p=p+TINY
       IF (p.GT.1) p=p-TINY
C 
C Now calculate derivatives
C 
       dlambdabydx=-xe/r(b,c)**2
       dlambdabydy=-ye/r(b,c)**2
       dlambdabydz=-ze/r(b,c)**2
       dubydx=(mu*xe)+x(d)-x(b)+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*r(b,c)**2))
       dubydy=(mu*ye)+y(d)-y(b)+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*r(b,c)**2))
       dubydz=(mu*ze)+z(d)-z(b)+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*r(b,c)**2))

       dvbydx=(vp/up)*(-xf-dlambdabydx*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydy=(vp/up)*(-yf-dlambdabydy*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydz=(vp/up)*(-zf-dlambdabydz*((xf*xe)+(yf*ye)+(zf*ze)))

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
C 
C Now add terms to energy gradients
C 
       CALL hairy(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position "b"
C 
       i=da2(j)      

       vixen1=(r(b,c)**2)*(x(a)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(x(b)-x(c))
       dlambdabydx=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(y(a)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(y(b)-y(c))
       dlambdabydy=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(z(a)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(z(b)-z(c))
       dlambdabydz=vixen1/(r(b,c)**4)

       vixen1=(r(b,c)**2)*(x(d)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(b)-x(c))
       dmubydx=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(y(d)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(b)-y(c))
       dmubydy=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(z(d)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(b)-z(c))
       dmubydz=vixen1/(r(b,c)**4)

       vixen1=mu*(2.0D0*x(b)-x(c)-x(a))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*x(b)-x(c)-x(d))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydx)+(mu*dlambdabydx))*r(b,c)**2
       dubydx=2.0D0*x(b)-x(a)-x(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(x(b)-x(c))
  
       vixen1=mu*(2.0D0*y(b)-y(c)-y(a))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*y(b)-y(c)-y(d))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydy)+(mu*dlambdabydy))*r(b,c)**2
       dubydy=2.0D0*y(b)-y(a)-y(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(y(b)-y(c))

       vixen1=mu*(2.0D0*z(b)-z(c)-z(a))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*z(b)-z(c)-z(d))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydz)+(mu*dlambdabydz))*r(b,c)**2
       dubydz=2.0D0*z(b)-z(a)-z(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(z(b)-z(c))
       
       dupbydx=(xf*(lambda+1)-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(yf*(lambda+1)-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(zf*(lambda+1)-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(xg*(mu+1)-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(yg*(mu+1)-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(zg*(mu+1)-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
      
       CALL hairy(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position 3
C 
       i=da3(j)

       dlambdabydx=(r(b,c)**2*(x(b)-x(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (x(c)-x(b)))/r(b,c)**4    
       dlambdabydy=(r(b,c)**2*(y(b)-y(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (y(c)-y(b)))/r(b,c)**4
       dlambdabydz=(r(b,c)**2*(z(b)-z(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (z(c)-z(b)))/r(b,c)**4     

       dmubydx=(r(b,c)**2*(x(b)-x(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(c)-x(b)))
     1         /r(b,c)**4
       dmubydy=(r(b,c)**2*(y(b)-y(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(c)-y(b)))
     1         /r(b,c)**4
       dmubydz=(r(b,c)**2*(z(b)-z(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(c)-z(b)))
     1         /r(b,c)**4

       vixen1=mu*(x(a)-x(b))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydx)+(mu*dlambdabydx))*r(b,c)**2+2.0D0*lambda*mu*xe
       dubydx=vixen1+lambda*(x(d)-x(b))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(y(a)-y(b))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydy)+(mu*dlambdabydy))*r(b,c)**2+2.0D0*lambda*mu*ye
       dubydy=vixen1+lambda*(y(d)-y(b))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(z(a)-z(b))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydz)+(mu*dlambdabydz))*r(b,c)**2+2.0D0*lambda*mu*ze
       dubydz=vixen1+lambda*(z(d)-z(b))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       dupbydx=(-lambda*xf-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(-lambda*yf-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(-lambda*zf-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(-mu*xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(-mu*yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(-mu*zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairy(dpbydx,dpbydy,dpbydz)       
C 
C Now for all atoms in position 4
C 
       i=da4(j)

       dmubydx=-xe/r(b,c)**2
       dmubydy=-ye/r(b,c)**2
       dmubydz=-ze/r(b,c)**2

       vixen1=lambda*dmubydx*r(b,c)**2
       vixen2=lambda*dmubydy*r(b,c)**2
       vixen3=lambda*dmubydz*r(b,c)**2
       dubydx=(lambda*xe)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydx)+x(a)-x(b)+vixen1
       dubydy=(lambda*ye)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydy)+y(a)-y(b)+vixen2
       dubydz=(lambda*ze)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydz)+z(a)-z(b)+vixen3

       dvbydx=up*(-xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvbydy=up*(-yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvbydz=up*(-zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairy(dpbydx,dpbydy,dpbydz)
      ENDDO 
20    CONTINUE

      RETURN
      END

      SUBROUTINE hairy(dpbydx,dpbydy,dpbydz)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ
      DOUBLE PRECISION VIXEN1,VIXEN2,VIXEN3

       e=type(a)
       f=type(b)
       g=type(c)
       h=type(d)
       PK=dvn(j)
       PK2=dvn2(j)
       PK3=dvn3(j)
       PN=dn(j)
       PN2=dn2(j)
       PN3=dn3(j)
       PHASE=ddelta(j)
       PHASE2=ddelta2(j)
       PHASE3=ddelta3(j)
       IDIVF=did(j)

       vixen1=PK*PN/IDIVF
       vixen2=PK2*PN2/IDIVF
       vixen3=PK3*PN3/IDIVF       

       IF (PN.EQ.2) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen1*dpbydx*2.0D0*COS(dphi(j))*COS(PHASE)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen1*dpbydy*2.0D0*COS(dphi(j))*COS(PHASE)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen1*dpbydz*2.0D0*COS(dphi(j))*COS(PHASE)
       ELSE IF (PN.EQ.3) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen1*dpbydx*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen1*dpbydy*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen1*dpbydz*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE)
       ELSE IF (PN.EQ.4) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen1*dpbydx*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen1*dpbydy*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen1*dpbydz*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE)
       ENDIF

       IF (PN2.EQ.1) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*COS(PHASE2)
       ELSE IF (PN2.EQ.2) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*2.0D0*COS(dphi(j))*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*2.0D0*COS(dphi(j))*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*2.0D0*COS(dphi(j))*COS(PHASE2)
       ELSE IF (PN2.EQ.3) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE2)
       ELSE IF (PN2.EQ.4) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen2*dpbydx*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE2)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen2*dpbydy*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE2)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen2*dpbydz*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE2)
       ENDIF

       IF (PN3.EQ.1) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*COS(PHASE3)
       ELSE IF (PN3.EQ.2) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*2.0D0*COS(dphi(j))*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*2.0D0*COS(dphi(j))*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*2.0D0*COS(dphi(j))*COS(PHASE3)
       ELSE IF (PN3.EQ.3) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*(3.0D0-4.0D0*(SIN(dphi(j)))**2)*COS(PHASE3)
       ELSE IF (PN3.EQ.4) THEN
         dtorsEbydx(i)=dtorsEbydx(i)+vixen3*dpbydx*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE3)
         dtorsEbydy(i)=dtorsEbydy(i)+vixen3*dpbydy*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE3)
         dtorsEbydz(i)=dtorsEbydz(i)+vixen3*dpbydz*(8.0D0*COS(dphi(j))**3-4.0D0*COS(dphi(j)))*COS(PHASE3)
       ENDIF

       RETURN
       END

      SUBROUTINE hairyimp(dpbydx,dpbydy,dpbydz)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,VIXEN1

      e=type(a)
      f=type(b)
      g=type(c)
      h=type(d)
      IPK=ivn(j)
      IPHASE=idelta(j)
      IPN=in1(j)

      vixen1=IPK*IPN
   
      IF (IPN.EQ.2) THEN
        dimpEbydx(i)=dimpEbydx(i)+2.0D0*dpbydx*vixen1*COS(iphi(j))*COS(IPHASE)       
        dimpEbydy(i)=dimpEbydy(i)+2.0D0*dpbydy*vixen1*COS(iphi(j))*COS(IPHASE)       
        dimpEbydz(i)=dimpEbydz(i)+2.0D0*dpbydz*vixen1*COS(iphi(j))*COS(IPHASE)       
      ELSE IF (IPN.EQ.0) THEN

      ELSE 
        PRINT *,"IPN for this torsion is not catered for in derivatives"
        STOP
      ENDIF

      RETURN
      END
       
      SUBROUTINE impderivs
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION XE,XF,XG,YE,YF,YG,ZE,ZF,ZG,U,V,UP,VP,P
      DOUBLE PRECISION DLAMBDABYDX,DLAMBDABYDY,DLAMBDABYDZ,DMUBYDX,DMUBYDY,DMUBYDZ
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,DUBYDX,DUBYDY,DUBYDZ,DVBYDX,DVBYDY,DVBYDZ
      DOUBLE PRECISION DPHIBYDX,DPHIBYDY,DPHIBYDZ,DUPBYDX,DUPBYDY,DUPBYDZ,DVPBYDX,DVPBYDY,DVPBYDZ
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3
      INTEGER J1
C 
C Initialise derivatives
C 
      DO J1=1,atoms
         dimpEbydx(J1)=0.0D0
         dimpEbydy(J1)=0.0D0
         dimpEbydz(J1)=0.0D0
      ENDDO
C 
C Calculate derivatives for all atoms in position "a"
C 
      DO j=1,imp
C 
C Initialise intermediates
C 
       dubydx=0.0D0
       dubydy=0.0D0
       dubydz=0.0D0
       dvbydx=0.0D0
       dvbydy=0.0D0
       dvbydz=0.0D0
       dpbydx=0.0D0
       dpbydy=0.0D0
       dupbydx=0.0D0
       dupbydy=0.0D0
       dupbydz=0.0D0
       dvpbydx=0.0D0
       dvpbydy=0.0D0
       dvpbydz=0.0D0
       dpbydz=0.0D0

       i=ia1(j)
       IF (i.EQ.0) GOTO 10
       a=ia1(j)
       b=ia2(j)
       c=ia3(j)
       d=ia4(j)
C 
C First define all functions as on paper
C 
       xe=x(c)-x(b)
       ye=y(c)-y(b)
       ze=z(c)-z(b)
       lambda=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(a)*xe)+(y(a)*ye)+(z(a)*ze)))/(xe**2+ye**2+ze**2)
       mu=(((x(b)*xe)+(y(b)*ye)+(z(b)*ze))-((x(d)*xe)+(y(d)*ye)+(z(d)*ze)))/(xe**2+ye**2+ze**2)
       xf=x(b)-x(a)-(lambda*xe)
       yf=y(b)-y(a)-(lambda*ye)
       zf=z(b)-z(a)-(lambda*ze)
       xg=x(b)-x(d)-(mu*xe)
       yg=y(b)-y(d)-(mu*ye)
       zg=z(b)-z(d)-(mu*ze)
       vixen1=(x(b)**2+y(b)**2+z(b)**2)+mu*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)-((x(b)*x(d))+(y(b)*y(d))+(z(b)*z(d)))
       vixen3=(x(d)-x(b))*x(a)+(y(d)-y(b))*y(a)+(z(d)-z(b))*z(a)
       u=vixen1+vixen2+vixen3+(lambda*mu*r(b,c)**2)
       up=SQRT(xf**2+yf**2+zf**2)
       vp=SQRT(xg**2+yg**2+zg**2)
       v=up*vp
       p=u/v

        IF (p.LT.-1) p=p+TINY
        IF (p.GT.1) p=p-TINY
C 
C Now calculate derivatives
C 
       dlambdabydx=-xe/r(b,c)**2
       dlambdabydy=-ye/r(b,c)**2
       dlambdabydz=-ze/r(b,c)**2
       dubydx=(mu*xe)+x(d)-x(b)+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*r(b,c)**2))
       dubydy=(mu*ye)+y(d)-y(b)+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*r(b,c)**2))
       dubydz=(mu*ze)+z(d)-z(b)+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze+(mu*r(b,c)**2))

       dvbydx=(vp/up)*(-xf-dlambdabydx*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydy=(vp/up)*(-yf-dlambdabydy*((xf*xe)+(yf*ye)+(zf*ze)))
       dvbydz=(vp/up)*(-zf-dlambdabydz*((xf*xe)+(yf*ye)+(zf*ze)))

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
C 
C Now add terms to energy gradients
C 
       CALL hairyimp(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position "b"
C 
       i=ia2(j)      

       vixen1=(r(b,c)**2)*(x(a)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(x(b)-x(c))
       dlambdabydx=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(y(a)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(y(b)-y(c))
       dlambdabydy=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(z(a)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*(z(b)-z(c))
       dlambdabydz=vixen1/(r(b,c)**4)

       vixen1=(r(b,c)**2)*(x(d)+x(c)-2.0D0*x(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(b)-x(c))
       dmubydx=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(y(d)+y(c)-2.0D0*y(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(b)-y(c))
       dmubydy=vixen1/(r(b,c)**4)
       vixen1=(r(b,c)**2)*(z(d)+z(c)-2.0D0*z(b))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(b)-z(c))
       dmubydz=vixen1/(r(b,c)**4)

       vixen1=mu*(2.0D0*x(b)-x(c)-x(a))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*x(b)-x(c)-x(d))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydx)+(mu*dlambdabydx))*r(b,c)**2
       dubydx=2.0D0*x(b)-x(a)-x(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(x(b)-x(c))
  
       vixen1=mu*(2.0D0*y(b)-y(c)-y(a))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*y(b)-y(c)-y(d))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydy)+(mu*dlambdabydy))*r(b,c)**2
       dubydy=2.0D0*y(b)-y(a)-y(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(y(b)-y(c))

       vixen1=mu*(2.0D0*z(b)-z(c)-z(a))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=lambda*(2.0D0*z(b)-z(c)-z(d))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)
       vixen3=((lambda*dmubydz)+(mu*dlambdabydz))*r(b,c)**2
       dubydz=2.0D0*z(b)-z(a)-z(d)+vixen1+vixen2+vixen3+2.0D0*lambda*mu*(z(b)-z(c))
       
       dupbydx=(xf*(lambda+1)-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(yf*(lambda+1)-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(zf*(lambda+1)-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(xg*(mu+1)-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(yg*(mu+1)-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(zg*(mu+1)-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2
      
       CALL hairyimp(dpbydx,dpbydy,dpbydz)
C 
C Now for all atoms in position 3
C 
       i=ia3(j)

       dlambdabydx=(r(b,c)**2*(x(b)-x(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (x(c)-x(b)))/r(b,c)**4    
       dlambdabydy=(r(b,c)**2*(y(b)-y(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (y(c)-y(b)))/r(b,c)**4
       dlambdabydz=(r(b,c)**2*(z(b)-z(a))-2.0D0*((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)*
     1             (z(c)-z(b)))/r(b,c)**4     

       dmubydx=(r(b,c)**2*(x(b)-x(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(x(c)-x(b)))
     1         /r(b,c)**4
       dmubydy=(r(b,c)**2*(y(b)-y(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(y(c)-y(b)))
     1         /r(b,c)**4
       dmubydz=(r(b,c)**2*(z(b)-z(d))-2.0D0*((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)*(z(c)-z(b)))
     1         /r(b,c)**4

       vixen1=mu*(x(a)-x(b))+dmubydx*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydx)+(mu*dlambdabydx))*r(b,c)**2+2.0D0*lambda*mu*xe
       dubydx=vixen1+lambda*(x(d)-x(b))+dlambdabydx*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(y(a)-y(b))+dmubydy*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydy)+(mu*dlambdabydy))*r(b,c)**2+2.0D0*lambda*mu*ye
       dubydy=vixen1+lambda*(y(d)-y(b))+dlambdabydy*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       vixen1=mu*(z(a)-z(b))+dmubydz*((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)
       vixen2=((lambda*dmubydz)+(mu*dlambdabydz))*r(b,c)**2+2.0D0*lambda*mu*ze
       dubydz=vixen1+lambda*(z(d)-z(b))+dlambdabydz*((x(d)-x(b))*xe+(y(d)-y(b))*ye+(z(d)-z(b))*ze)+vixen2

       dupbydx=(-lambda*xf-dlambdabydx*(xe*xf+ye*yf+ze*zf))/up
       dupbydy=(-lambda*yf-dlambdabydy*(xe*xf+ye*yf+ze*zf))/up
       dupbydz=(-lambda*zf-dlambdabydz*(xe*xf+ye*yf+ze*zf))/up

       dvpbydx=(-mu*xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydy=(-mu*yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvpbydz=(-mu*zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dvbydx=(up*dvpbydx)+(vp*dupbydx)
       dvbydy=(up*dvpbydy)+(vp*dupbydy)
       dvbydz=(up*dvpbydz)+(vp*dupbydz)

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairyimp(dpbydx,dpbydy,dpbydz)       
C 
C Now for all atoms in position 4
C 
       i=ia4(j)

       dmubydx=-xe/r(b,c)**2
       dmubydy=-ye/r(b,c)**2
       dmubydz=-ze/r(b,c)**2

       vixen1=lambda*dmubydx*r(b,c)**2
       vixen2=lambda*dmubydy*r(b,c)**2
       vixen3=lambda*dmubydz*r(b,c)**2
       dubydx=(lambda*xe)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydx)+x(a)-x(b)+vixen1
       dubydy=(lambda*ye)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydy)+y(a)-y(b)+vixen2
       dubydz=(lambda*ze)+(((x(a)-x(b))*xe+(y(a)-y(b))*ye+(z(a)-z(b))*ze)*dmubydz)+z(a)-z(b)+vixen3

       dvbydx=up*(-xg-dmubydx*(xe*xg+ye*yg+ze*zg))/vp
       dvbydy=up*(-yg-dmubydy*(xe*xg+ye*yg+ze*zg))/vp
       dvbydz=up*(-zg-dmubydz*(xe*xg+ye*yg+ze*zg))/vp

       dpbydx=((v*dubydx)-(u*dvbydx))/v**2
       dpbydy=((v*dubydy)-(u*dvbydy))/v**2
       dpbydz=((v*dubydz)-(u*dvbydz))/v**2

       CALL hairyimp(dpbydx,dpbydy,dpbydz)
      ENDDO 
10    CONTINUE

      RETURN
      END


      SUBROUTINE numergrad
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER i1
      DOUBLE PRECISION         NUMERDERIV(3*NATOMS)
      DOUBLE PRECISION         ET

      DO i1=1,3*NATOMS
         numerderiv(i1)= 0.0
      END DO

      DO i1=1,atoms
         x(i1)= x(i1)+ 0.00001
         CALL amberenergy
         et= totenergy
         x(i1)= x(i1)- 0.00002
         CALL amberenergy
         numerderiv(3*i1-2)= (et-totenergy)/0.00002
         x(i1)= x(i1)+ 0.00001

         y(i1)= y(i1)+ 0.00001
         CALL amberenergy
         et= totenergy
         y(i1)= y(i1)- 0.00002
         CALL amberenergy
         numerderiv(3*i1-1)= (et-totenergy)/0.00002
         y(i1)= y(i1)+ 0.00001

         z(i1)= z(i1)+ 0.00001
         CALL amberenergy
         et= totenergy
         z(i1)= z(i1)- 0.00002
         CALL amberenergy
         numerderiv(3*i1)= (et-totenergy)/0.00002
         z(i1)= z(i1)+ 0.00001
      END DO

      WRITE (*,FMT='(5X,A,5X,A,5X,A,5X,A,5X,A)') 'NUMERICAL','   *   ','ANALYTIC','   *   ','DIFFERENCE'
      WRITE (*,FMT='(5X,A,5X,A,5X,A,5X,A,5X,A)') '---------','   *   ','--------','   *   ','----------'
      DO i1=1,atoms
         WRITE (*,FMT='(F20.10,A,F20.10,A,F20.10)') numerderiv(3*i1-2),'  * ',dEbydx(i1),'   * ',dEbydx(i1)-numerderiv(3*i1-2)
         WRITE (*,FMT='(F20.10,A,F20.10,A,F20.10)') numerderiv(3*i1-1),'  * ',dEbydy(i1),'   * ',dEbydy(i1)-numerderiv(3*i1-1)
         WRITE (*,FMT='(F20.10,A,F20.10,A,F20.10)') numerderiv(3*i1),'  * ',dEbydz(i1),'   * ',dEbydz(i1)-numerderiv(3*i1)
      END DO


      RETURN

      END

      SUBROUTINE AMBERA
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE

      DO a=1,ang
         b=aa1(a)
         c=aa2(a)
         d=aa3(a)

         top=(x(b)-x(c))*(x(d)-x(c))+(y(b)-y(c))*(y(d)-y(c))+(z(b)-z(c))*(z(d)-z(c))
         bottom=r(b,c)*DSQRT((x(d)-x(c))**2+(y(d)-y(c))**2+(z(d)-z(c))**2)
         theta(a)=ACOS(top/bottom)
      END DO

      DO a=1,t
         i=da1(a)
         j=da2(a)
         k=da3(a)
         l=da4(a)
   
         colin=0
         ax=x(j)-x(i)
         ay=y(j)-y(i)
         az=z(j)-z(i)
         bx=x(k)-x(j)
         by=y(k)-y(j)
         bz=z(k)-z(j)
         cx=x(j)-x(l)
         cy=y(j)-y(l)
         cz=z(j)-z(l)
         dx=x(l)-x(k)
         dy=y(l)-y(k)
         dz=z(l)-z(k)
   
         IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) .AND. ay*bz.GT.(by*az-TINY) 
     1        .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
         IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) .AND. by*dz.GT.(dy*bz-TINY) 
     1        .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2

         IF (colin.EQ.1) THEN
            PRINT *,'Three sites',i,j,k,'are colinear'
         ELSE IF (colin.EQ.2) THEN
            PRINT *,'Three sites',j,k,l,'are colinear'
         ELSE
            lambda=(ax*bx+ay*by+az*bz)/(bx**2+by**2+bz**2)
            mu=(cx*bx+cy*by+cz*bz)/(bx**2+by**2+bz**2)
            qx(1)=x(i)+(lambda*bx)
            qy(1)=y(i)+(lambda*by)
            qz(1)=z(i)+(lambda*bz)
            qx(2)=x(j)
            qy(2)=y(j)
            qz(2)=z(j)
            qx(3)=x(l)+(mu*bx)
            qy(3)=y(l)+(mu*by)
            qz(3)=z(l)+(mu*bz)
            numer=(qx(1)-qx(2))*(qx(3)-qx(2))+(qy(1)-qy(2))*(qy(3)-qy(2))+(qz(1)-qz(2))*(qz(3)-qz(2))
            denom=SQRT(((qx(1)-qx(2))**2+(qy(1)-qy(2))**2+(qz(1)-qz(2))**2)
     1          *((qx(3)-qx(2))**2+(qy(3)-qy(2))**2+(qz(3)-qz(2))**2))
   
            IF ((numer/denom).LT. -1) numer=numer+1D-12
            IF ((numer/denom).GT.1) numer=numer-1D-12
  
            dphi(a)=ACOS(numer/denom)
         ENDIF
      END DO

      DO a=1,imp 
         e=ia1(a)
         f=ia2(a)
         g=ia3(a)
         h=ia4(a)    
         colin=0
         ax=x(f)-x(e)
         ay=y(f)-y(e)
         az=z(f)-z(e)
         bx=x(g)-x(f)
         by=y(g)-y(f)
         bz=z(g)-z(f)
         cx=x(f)-x(h)
         cy=y(f)-y(h)
         cz=z(f)-z(h)
         dx=x(h)-x(g)
         dy=y(h)-y(g)
         dz=z(h)-z(g)

         IF (ax*by.LT.(bx*ay+TINY) .AND. ax*by.GT.(bx*ay-TINY) .AND. ay*bz.LT.(by*az+TINY) .AND. ay*bz.GT.(by*az-TINY) 
     1      .AND. ax*bz.LT.(bx*az+TINY) .AND. ax*bz.GT.(az*bx-TINY)) colin=1
         IF (bx*dy.LT.(dx*by+TINY) .AND. bx*dy.GT.(dx*by-TINY) .AND. by*dz.LT.(dy*bz+TINY) .AND. by*dz.GT.(dy*bz-TINY) 
     1      .AND. bx*dz.LT.(bz*dx+TINY) .AND. bx*dz.GT.(bz*dx+TINY)) colin=2
         IF (colin.EQ.1) THEN
            PRINT *,'Three sites',e,f,g,' are colinear'
         ELSE IF (colin.EQ.2) THEN
            PRINT *,'Three sites',f,g,h,' are colinear'
         ELSE

            lambda=(ax*bx+ay*by+az*bz)/(bx**2+by**2+bz**2)
            mu=(cx*bx+cy*by+cz*bz)/(bx**2+by**2+bz**2)
            qx(1)=x(e)+(lambda*bx)
            qy(1)=y(e)+(lambda*by)
            qz(1)=z(e)+(lambda*bz)
            qx(2)=x(f)
            qy(2)=y(f)
            qz(2)=z(f)
            qx(3)=x(h)+(mu*bx)
            qy(3)=y(h)+(mu*by)
            qz(3)=z(h)+(mu*bz)
              numer=(qx(1)-qx(2))*(qx(3)-qx(2))+(qy(1)-qy(2))*(qy(3)-qy(2))+(qz(1)-qz(2))*(qz(3)-qz(2))
            denom=SQRT(((qx(1)-qx(2))**2+(qy(1)-qy(2))**2+(qz(1)-qz(2))**2)
     1          *((qx(3)-qx(2))**2+(qy(3)-qy(2))**2+(qz(3)-qz(2))**2))

            IF ((numer/denom).LT. -1) numer=numer+1D-12
            IF ((numer/denom).GT.1) numer=numer-1D-12

            iphi(a)=ACOS(numer/denom)
         ENDIF
      END DO

      RETURN

      END
      SUBROUTINE amberenergy
      USE KEY
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3, VIXEN4, VIXEN5, VIXEN6, E1, E2, E3
      CHARACTER FNAMEF*80

C 
C First initialising count variables
C 
      benergy=0.0
      tenergy=0.0
      penergy=0.0
      vdwenergy=0.0
      totenergy=0.0
      impenergy=0.0
      qenergy=0.0
C 
C Calculating r,benergy,vdwenergy,qenergy - note that we calculate ALL r as we need them for vdW energy
C 
      CALL AMBERD

      DO i=1,bondnumber
         a=bondarray(i,1)
         b=bondarray(i,2)
         benergy=benergy+kr(type(a),type(b))*((r(a,b)-ro(type(a),type(b)))**2)
      ENDDO

      CALL AMBERA

      DO a=1,ang
         b=aa1(a)
         c=aa2(a)
         d=aa3(a)

         vixen1=kt(type(b),type(c),type(d))*((theta(a)-to(type(b),type(c),type(d)))**2)
         tenergy=tenergy+vixen1
C        PRINT *,'Angle ',b,c,d,typech(b),'-',typech(c),'-',typech(d),' theta=',theta(a)*57.29577951,' energy=',vixen1
      ENDDO

      DO a=1,t
         i=da1(a)
         j=da2(a)
         k=da3(a)
         l=da4(a)
         vixen1=dvn(a)/did(a)
         vixen2=dn(a)*dphi(a)-ddelta(a)
 
         vixen3=dvn2(a)/did(a)
         vixen4=dn2(a)*dphi(a)-ddelta2(a)

         vixen5=dvn3(a)/did(a)
         vixen6=dn3(a)*dphi(a)-ddelta3(a)

         penergy=penergy+(vixen1*(1+cos(vixen2)))+(vixen3*(1+cos(vixen4)))+(vixen5*(1+cos(vixen6)))
         E1=vixen1*(1+cos(vixen2))
         E2=(vixen3*(1+cos(vixen4)))
         E3=(vixen5*(1+cos(vixen6)))

C         PRINT *,'Torsion',i,j,k,l,typech(i),'-',typech(j),'-',typech(k),'-',typech(l),' Phi=',dphi(a)*57.29577951,
C     1           ' Energy=',E1+E2+E3
C         PRINT *,'PK=',dvn(a),' IDIVF=',did(a),' PN=',dn(a),' PHASE=',ddelta(a)

      ENDDO

      DO a=1,imp
         impenergy=impenergy+(ivn(a)*(1+COS(in1(a)*iphi(a)-idelta(a))))
     
C        PRINT *,'Improper',e,f,g,h,' energy= ',ivn(a)*(1+COS(in1(a)*iphi(a)-idelta(a)))

      ENDDO

      IF (count.EQ.1) THEN 
         PRINT *,'ang= ',ang
         PRINT *,'t= ',t
         PRINT *,'imp= ',imp
      ENDIF
C     DO a=1,imp
C        PRINT *,ia1(a),ia2(a),ia3(a),ia4(a)
C        PRINT *,'phi= ',impropers(a)%phi*57.29577951
C     ENDDO

C     PRINT *,'Bond lengths'
C       DO a=1,atoms
C        DO b=a,atoms
C        IF (bonds(a,b).EQ.1) PRINT *,a,'-',b,'    ',r(a,b),'  ',typech(a),' ',typech(b)
C        ENDDO
C       ENDDO

C       PRINT *,'Bond angles'
C         DO a=1,ang
C          PRINT *,angles(a)%a1,'-',angles(a)%a2,'-',angles(a)%a3,' = ',(angles(a)%theta)*57.29577951
C         ENDDO

C       PRINT *,'Torsion angles'
C       DO a=1,t
C        PRINT *,da1(a),torsions(a)%a2,torsions(a)%a3,torsions(a)%a4,'    ',&
C                57.29577951*torsions(a)%phi
C       ENDDO
    
       PRINT *,'Total bond strain energy= ',benergy
       PRINT *,'Total angle strain energy= ',tenergy
       PRINT *,'Torsion angle energy= ',penergy
       PRINT *,'Improper torsion energy= ',impenergy
       PRINT *,'vdW energy= ',vdwenergy
       PRINT *,'qenergy= ',qenergy
       totenergy=benergy+tenergy+penergy+vdwenergy+impenergy+qenergy
       PRINT *,'Total energy= ',totenergy
C      IF (count.EQ.1) PRINT *,'Energy calculated'

       IF (TOTENERGY.LT.-1.0D4) THEN  
          WRITE(*,'(A)') 'Cold fusion diagnosed - quit'
          IF (FILTH.EQ.0) THEN
             FNAMEF='points.final'
          ELSE
             WRITE(FNAMEF,'(A)') 'points.final.'//TRIM(ADJUSTL(FILTHSTR))
          ENDIF
C         CALL amberdump(Q,FNAMEF)
          STOP
       ENDIF

      RETURN
      END

      SUBROUTINE amoviedump(frame)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER frame, moviecount

      WRITE (27,FMT=*) atoms
      WRITE (27,FMT='(A,I4)') 'Frame no. ',frame
      DO moviecount=1,atoms
         WRITE (27,FMT='(A,3F20.10)') typech(moviecount)(1:1), x(moviecount), y(moviecount), z(moviecount)
      END DO

      frame=frame+1

      RETURN

      END
      SUBROUTINE aread
      USE COMMONS
      USE KEY
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER canine
      COMMON /CHIR/ canine
      DOUBLE PRECISION VIXEN1
      DOUBLE PRECISION CRAP
      PARAMETER (crap=18.2223D0)
      CHARACTER(LEN=17) filthfile

      ALLOCATE(DATOM1(NATOMS),DATOM2(NATOMS))
      ALLOCATE(bondarray(NATOMS*4,2))
      ALLOCATE(x(NATOMS),y(NATOMS),z(NATOMS),pq(NATOMS),mass(NATOMS))
      ALLOCATE(r(NATOMS,NATOMS),vdwa(NATOMS,NATOMS),vdwb(NATOMS,NATOMS))
      ALLOCATE(aa1(NATOMS*5),aa2(NATOMS*5),aa3(NATOMS*5))
      ALLOCATE(theta(NATOMS*5))
      ALLOCATE(da1(NATOMS*5),da2(NATOMS*5),da3(NATOMS*5),da4(NATOMS*5))
      ALLOCATE(dphi(NATOMS*5),did(NATOMS*5),dvn(NATOMS*5),dvn2(NATOMS*5))
      ALLOCATE(dvn3(NATOMS*5),ddelta(NATOMS*5),ddelta2(NATOMS*5),ddelta3(NATOMS*5))
      ALLOCATE(dn(NATOMS*5),dn2(NATOMS*5),dn3(NATOMS*5))
      ALLOCATE(ia1(NATOMS*5),ia2(NATOMS*5),ia3(NATOMS*5),ia4(NATOMS*5))
      ALLOCATE(iphi(NATOMS*5),ivn(NATOMS*5),idelta(NATOMS*5),in1(NATOMS*5))
      ALLOCATE(atnum(NATOMS),bondedto(NATOMS),type(NATOMS))
      ALLOCATE(bonds(NATOMS,NATOMS))
      ALLOCATE(one_four(NATOMS,NATOMS),one_three(NATOMS,NATOMS))
      ALLOCATE(label(NATOMS))
      ALLOCATE(typech(NATOMS))
      ALLOCATE(dbondEbydx(NATOMS),dbondEbydy(NATOMS),dbondEbydz(NATOMS),dangEbydx(NATOMS), 
     &                 dangEbydy(NATOMS),dangEbydz(NATOMS), 
     &                 dtorsEbydx(NATOMS),dtorsEbydy(NATOMS),dtorsEbydz(NATOMS),dvdwEbydx(NATOMS), 
     &                 dvdwEbydy(NATOMS),dvdwEbydz(NATOMS), 
     &                 dEbydx(NATOMS),dEbydy(NATOMS),dEbydz(NATOMS),dimpEbydx(NATOMS), 
     &                 dimpEbydy(NATOMS),dimpEbydz(NATOMS), 
     &                 dqEbydx(NATOMS),dqEbydy(NATOMS),dqEbydz(NATOMS))
      ALLOCATE(bondhell(3*NATOMS,3*NATOMS), anglehell(3*NATOMS,3*NATOMS), torshell(3*NATOMS,3*NATOMS),
     1         imphell(3*NATOMS,3*NATOMS), qhell(3*NATOMS,3*NATOMS), vdwhell(3*NATOMS,3*NATOMS),
     2         hell(3*NATOMS,3*NATOMS))
      ALLOCATE(CHIRAL(NATOMS),CHIRALARRAY(NATOMS,5))

      IF (FILTH.EQ.0) THEN
         OPEN (UNIT=9,FILE='coords.amber',STATUS='OLD')
      ELSE 
         WRITE (filthfile,'(A)') 'coords.amber.'//TRIM(ADJUSTL(FILTHSTR))
         OPEN (UNIT=9,FILE=filthfile,STATUS='OLD')
      ENDIF
C 
C Format for info file   label  typech  number_in_list  bondedto  x  y  z
C 
      PRINT *,'Max. no. of atoms = ',NATOMS
      DO a=1,NATOMS
        READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
        IF (ios.LT.0) THEN
          PRINT *,'End of file before all information specified',a
          STOP
        END IF
        IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') GOTO 10
        BACKSPACE 9
        READ (UNIT=9,FMT=*) label(a),typech(a),atnum(a),bondedto(a),x(a),y(a),z(a)
        ZSYM(a)=typech(a)(1:1)
        IF (label(a).EQ.'*') THEN
          chiral(a)=1
C         PRINT *,'Found chiral atom',a
        END IF
      END DO
10    CONTINUE

      PRINT *,'Co-ordinates read successfully - bondedto(5)= ',bondedto(5)
      atoms=a-1
      PRINT *,'atoms= ',atoms
C 
C Check for any rings and close them!
C 
      DO a=1,10
       READ (UNIT=9,IOSTAT=ios,FMT='(A4)') check
       IF (check.NE.'loop' .AND. a.EQ.10) GOTO 20
       IF (.NOT.(check.NE.'loop')) THEN
          BACKSPACE 9
          READ (UNIT=9,FMT='(A4,7X,I2)') check,rings
          GOTO 20
      ENDIF
      END DO
20    CONTINUE

      bean=1
      IF (a.LT.10) THEN
       DO i=1,rings
C        READ (UNIT=9,FMT='(I3,2X,I3)') loopatom(2*i-1),loopatom(2*i)
        READ (UNIT=9,FMT=*) loopatom(2*i-1),loopatom(2*i)
         b=loopatom(2*i-1)
         c=loopatom(2*i)
         IF (b.EQ.0 .OR. c.EQ.0) THEN
          PRINT *,'No atoms specified for loop bond'
         END IF
        bonds(b,c)=1
       bonds(c,b)=1
       bondarray(bean,1)=b
       bondarray(bean,2)=c
       bean=bean+1
       END DO
      ELSE
       PRINT *,'No loop line in co-ordinate file'
       STOP
      END IF
C 
C Read charges
C 
      DO a=1,10
       READ (UNIT=9,IOSTAT=ios,FMT='(A7)') check
       IF (check.NE.'charges' .AND. a.EQ.10) GOTO 30
       IF (.NOT.(check.NE.'charges')) GOTO 30
      END DO
30    CONTINUE

      IF (a.EQ.10) THEN
       PRINT *,'No charges specified'
       STOP
      ELSE 
      DO c=1,atoms
C       READ (UNIT=9,IOSTAT=ios,FMT='(I3,2X,F7.4)') b,vixen1
       READ (UNIT=9,IOSTAT=ios,FMT=*) b,vixen1
       pq(b)=vixen1*crap
       IF (ios.LT.0) GOTO 40
      END DO
40    CONTINUE
      END IF
C 
C Converting typech to type number
C 
      DO a=1,atoms
       CALL typenumber(typech(a))
       type(a)=ans
       IF (type(a).EQ.-1) THEN
        PRINT *,'Atom ',a,'type specified incorrectly - aborting'
      STOP
       END IF
      END DO
C 
C Connecting backbone atoms
C 
      DO a=2,atoms
       IF (label(a).NE.'b') GOTO 50
       bonds(a,a-1)=1
       bonds(a-1,a)=1
      END DO       
50    CONTINUE
C 
C Connecting side-chain atoms
C 
      DO a=2,atoms
       IF (.NOT.(label(a).EQ.'b')) THEN
       bonds(a,bondedto(a))=1
       bonds(bondedto(a),a)=1
       bondarray(bean,1)=a
       bondarray(bean,2)=bondedto(a)
       bean=bean+1
       ENDIF
      END DO

      bondnumber=bean-1

      CALL one3

      RETURN
      END

      SUBROUTINE one3
      USE MODAMBER
      USE MODAMBER2
       IMPLICIT NONE
C      PRINT *,'SUBROUTINE one3'
       DO a=1,atoms
         DO b=a+1,atoms
           IF (bonds(a,b) .NE. 1) THEN
             DO c=1,atoms
               IF (bonds(a,c).EQ.1 .AND. bonds(b,c).EQ.1) THEN
                 one_three(a,b)=1
                 one_three(b,a)=1
C                PRINT *,a,b
               END IF
             END DO
           END IF
         END DO
       END DO

      RETURN

      END
      SUBROUTINE aparams
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3, VIXEN4
      DOUBLE PRECISION ANGFAC
      PARAMETER (angfac=57.29577951D0)
     
      PRINT *,'Reading parameters'
      hellcount=0

      OPEN (UNIT=9,FILE='amber.dat',STATUS='OLD')
C 
C Reading parameters
C 
C Format of params.dat for bonds  typechb typechc kr(typechb,typechc) ro(typechb,typechc)
C 
      DO a=1,200
         READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
         IF (ios.LT.0) GOTO 10
         IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') GOTO 10
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2)') typechb,typechc

         CALL typenumber(typechb)
         b=ans
         CALL typenumber(typechc)
         c=ans 
         IF (b.EQ.-1 .OR. c.EQ.-1) THEN
          PRINT *,'Bond parameters ',typechb,'-',typechc,' specified incorrectly - aborting'
          STOP
         ENDIF

         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2,2X,F5.0,4X,F6.0)') typechb,typechc,kr(b,c),ro(b,c)
C        WRITE(*,'(3I3,X,A,2X,A,2X,2F15.5)') a,b,c,typechb,typechc,kr(b,c),ro(b,c)
         kr(c,b)=kr(b,c)
         ro(c,b)=ro(b,c)
      ENDDO
10    CONTINUE
C 
C Format of params.dat for angles  type1 type2 type3 kt(1,2,3) degto(1,2,3)
C 
      DO a=1,300
         READ (UNIT=9,IOSTAT=ios,FMT='(A5)') check 
         IF (ios.LT.0) GOTO 20
         IF (.NOT.(check.EQ.'start' .OR. check.EQ.'Start' .OR. check.EQ.'START')) THEN
         IF (.NOT.(check.NE.'start' .AND. check.NE.'Start' .AND. check.NE.'START' .AND. a.EQ.1)) THEN
         IF (check.EQ.'end' .OR. check.EQ.'End' .OR. check.EQ.'END') GOTO 20
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2,X,A2)') typechb,typechc,typechd

         CALL typenumber(typechb)
         b=ans
         CALL typenumber(typechc)
         c=ans
         CALL typenumber(typechd)
         d=ans

         IF (b.EQ.-1 .OR. c.EQ.-1 .OR. d.EQ.-1) THEN
            PRINT *,'Angle ',typechb,'-',typechc,'-',typechd,' specified incorrectly - aborting'
            STOP
         ENDIF
   
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2,X,A2,3X,F5.0,6X,F6.0)') typechb,typechc,typechd,kt(b,c,d),degto(b,c,d)
         to(b,c,d)=degto(b,c,d)/angfac
         kt(d,c,b)=kt(b,c,d)
         to(d,c,b)=to(b,c,d)
         degto(d,c,b)=degto(b,c,d)
         ENDIF
         ENDIF
      ENDDO
20    CONTINUE

      IF (a.EQ.1) THEN
       PRINT *,'Error reading angle parameters - no start'
       STOP
      ENDIF 
C 
C Format of params.dat for torsions  type1 type2 type3 type4 PK IDIVF PN PHASE
C 
      DO a=1,100
         READ (UNIT=9,IOSTAT=ios,FMT='(A5)') check
         IF (ios.LT.0) GOTO 30
         IF (.NOT.(check.EQ.'start' .OR. check.EQ.'Start' .OR. check.EQ.'START')) THEN
         IF (check.NE.'start' .AND. check.NE.'Start' .AND. check.NE.'START' .AND. a.EQ.1) GOTO 30
         IF (check.EQ.'end' .OR. check.EQ.'End' .OR. check.EQ.'END')  GOTO 30
         IF (a.GT.1) BACKSPACE 9
         
         READ (UNIT=9,IOSTAT=ios,FMT='(A2)') typechb
         IF (typechb.NE.'X ') GOTO 30
         BACKSPACE 9

         READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,13X,F2.0)') typechb,typechc,typechd,
     1        typeche,vixen1,vixen2,vixen3,vixen4
         CALL typenumber(typechc)
           c=ans
         CALL typenumber(typechd)
           d=ans

         gentorsparams(a-1,1)=c
         gentorsparams(a-1,2)=d
         gentorsparams(a-1,3)=vixen1 
         gentorsparams(a-1,4)=vixen2
         gentorsparams(a-1,5)=vixen3/angfac
         gentorsparams(a-1,6)=ABS(vixen4)

         ENDIF
      ENDDO
30    CONTINUE
 
      IF (a.EQ.1) THEN
       PRINT *,'Error reading torsion parameters - no start'
      STOP
      ENDIF

      BACKSPACE 9
      DO a=1,50
        READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
C        PRINT *,"check= ",check
        IF (check.EQ."end") THEN
C        PRINT *,"End of torsion parameters reached - a=",a
        GOTO 32
        END IF
        BACKSPACE 9

        READ (UNIT=9,IOSTAT=ios,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F3.0)') typechb,typechc,
     1       typechd,typeche,vixen1,vixen2,vixen3,vixen4
        CALL typenumber(typechb)
          b=ans
        CALL typenumber(typechc)
          c=ans
        CALL typenumber(typechd)
          d=ans
        CALL typenumber(typeche)
          e=ans

        spectorsparams(a,1)=b
        spectorsparams(a,2)=c
        spectorsparams(a,3)=d
        spectorsparams(a,4)=e
        spectorsparams(a,5)=vixen1
        spectorsparams(a,6)=vixen2
        spectorsparams(a,7)=vixen3/angfac
        spectorsparams(a,8)=ABS(vixen4)

        IF (vixen4.GT.0) THEN
         spectorsparams(a,9)=0.0
         spectorsparams(a,10)=0.0
         spectorsparams(a,11)=0.0
         spectorsparams(a,12)=0.0
         spectorsparams(a,13)=0.0
         spectorsparams(a,14)=0.0
         GOTO 31
        END IF

        READ (UNIT=9,IOSTAT=ios,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F2.0)') typechb,typechc,
     1      typechd,typeche,vixen1,vixen2,vixen3,vixen4

        spectorsparams(a,9)=vixen2
        spectorsparams(a,10)=vixen3/angfac
        spectorsparams(a,11)=ABS(vixen4)
       
        IF (vixen4.GT.0) THEN
          spectorsparams(a,12)=0.0
          spectorsparams(a,13)=0.0
          spectorsparams(a,14)=0.0
          GOTO 31
        END IF

        READ (UNIT=9,IOSTAT=ios,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F2.0)') typechb,typechc,
     1    typechd,typeche,vixen1,vixen2,vixen3,vixen4

        spectorsparams(a,12)=vixen2
        spectorsparams(a,13)=vixen3/angfac
        spectorsparams(a,14)=vixen4

31     CONTINUE
       END DO
32     CONTINUE
C 
C Format for improper torsion parameters    typecha  typechb  typechc  typechd  IPK  IPHASE  IPN
C 
      DO a=1,100
       READ (UNIT=9,IOSTAT=ios,FMT='(A9)') check
       IF (ios.LT.0) THEN
        PRINT *,'End of parameters file without improper torsion parameters'
        STOP
       ENDIF
       IF (.NOT.(check.NE.'impropers' .AND. check.NE.'IMPROPERS' .AND. check.NE.'Impropers')) GOTO 40
      ENDDO
40    CONTINUE

      IF (a.GE.100) THEN 
       PRINT *,'No improper torsion parameters specified'
      STOP
      ENDIF

      bean=0
      DO a=1,100
       READ (UNIT=9,IOSTAT=ios,FMT='(A2,X,A2,X,A2,X,A2)') typechb,typechc,typechd,typeche
       IF (ios.LT.0) THEN
        PRINT *,'End of file during improper torsion parameters'
      STOP
       ENDIF
      IF (typechb.EQ.'X ' .AND. typechc.EQ.'X ') THEN
        BACKSPACE 9 
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2)') typechb,typechc,typechd,typeche
        CALL typenumber(typechd)
          d=ans
        CALL typenumber(typeche)
          e=ans
 
        genimpparams(a-bean,1)=d
        genimpparams(a-bean,2)=e

        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') typechb,typechc,typechd,typeche,
     1       vixen1,vixen2,vixen3
        genimpparams(a-bean,3)=vixen1
        genimpparams(a-bean,4)=vixen2/angfac
        genimpparams(a-bean,5)=vixen3

       ELSE IF (typechb.EQ.'X ' .AND. typechc.NE.'X ') THEN
        bean=bean+1
        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2)') typechb,typechc,typechd,typeche
         CALL typenumber(typechc)
           c=ans
         CALL typenumber(typechd)
           d=ans
         CALL typenumber(typeche)
           e=ans

        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') typechb,typechc,typechd,typeche,
     1       vixen1,vixen2,vixen3

        midimpparams(bean,1)=c
        midimpparams(bean,2)=d        
        midimpparams(bean,3)=e
        midimpparams(bean,4)=vixen1
        midimpparams(bean,5)=vixen2/angfac
        midimpparams(bean,6)=vixen3

       ELSE IF (typechb.EQ.'  ') THEN
          GOTO 50

       ELSE IF (typechb.EQ.'en' .OR. typechb.EQ.'En' .OR. typechb.EQ.'EN') THEN
          GOTO 60
       ELSE 
          GOTO 60

       ENDIF
50     CONTINUE
      ENDDO   
60    CONTINUE

      BACKSPACE 9
      DO a=1,20
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') typechb,typechc,typechd,typeche,
     1       vixen1,vixen2,vixen3
        IF (typechb.EQ.'en') GOTO 65
        CALL typenumber(typechb)
          b=ans
        CALL typenumber(typechc)
          c=ans
        CALL typenumber(typechd)
          d=ans
        CALL typenumber(typeche)
          e=ans

        specimpparams(a,1)=b
        specimpparams(a,2)=c
        specimpparams(a,3)=d
        specimpparams(a,4)=e
        specimpparams(a,5)=vixen1
        specimpparams(a,6)=vixen2/angfac
        specimpparams(a,7)=vixen3

      END DO
65    CONTINUE

      DO a=1,100
       READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
       IF (ios.LT.0) THEN
        PRINT *,'End of parameter file without vdW parameters'
        STOP
       ENDIF
       IF (.NOT.(check.NE.'vdw' .AND. check.NE.'vdW')) GOTO 70
      ENDDO
70    CONTINUE
  
      IF (a.LT.100) THEN
       DO a=1,42
        READ (UNIT=9,IOSTAT=ios,FMT='(A2)') typechb
        IF (typechb.EQ.'en') GOTO 80
        IF (ios.LT.0) GOTO 80
        CALL typenumber(typechb)
        b=ans
        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,10X,F6.4,2X,F6.4)') typechb,vdwr(b),vdwe(b)
       ENDDO
80     CONTINUE
      ELSE
       PRINT *,'Unable to find vdW parameters'
      ENDIF

      PRINT *,'Parameters read'
    
      RETURN
      END

      SUBROUTINE typenumber(typechar)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      CHARACTER(LEN=2)  typechar

      IF (typechar.EQ.'C ') THEN
         ans=1
      ELSE IF (typechar.EQ.'CA') THEN
         ans=2
      ELSE IF (typechar.EQ.'CB') THEN
         ans=3
      ELSE IF (typechar.EQ.'CC') THEN
         ans=4
      ELSE IF (typechar.EQ.'CK') THEN
         ans=5
      ELSE IF (typechar.EQ.'CM') THEN
         ans=6
      ELSE IF (typechar.EQ.'CN') THEN
         ans=7
      ELSE IF (typechar.EQ.'CQ') THEN
         ans=8
      ELSE IF (typechar.EQ.'CR') THEN
         ans=9
      ELSE IF (typechar.EQ.'CT') THEN
         ans=10
      ELSE IF (typechar.EQ.'CV') THEN
         ans=11
      ELSE IF (typechar.EQ.'CW') THEN
         ans=12
      ELSE IF (typechar.EQ.'C*') THEN
         ans=13
      ELSE IF (typechar.EQ.'C0') THEN
         ans=14
      ELSE IF (typechar.EQ.'F ') THEN
         ans=15
      ELSE IF (typechar.EQ.'H ') THEN
         ans=16
      ELSE IF (typechar.EQ.'HC') THEN
         ans=17
      ELSE IF (typechar.EQ.'H1') THEN
         ans=18
      ELSE IF (typechar.EQ.'H2') THEN
         ans=19
      ELSE IF (typechar.EQ.'H3') THEN
         ans=20
      ELSE IF (typechar.EQ.'HA') THEN
         ans=21
      ELSE IF (typechar.EQ.'H4') THEN
         ans=22
      ELSE IF (typechar.EQ.'H5') THEN
         ans=23
      ELSE IF (typechar.EQ.'HO') THEN
         ans=24
      ELSE IF (typechar.EQ.'HS') THEN
         ans=25
      ELSE IF (typechar.EQ.'HW') THEN
         ans=26
      ELSE IF (typechar.EQ.'HP') THEN
         ans=27
      ELSE IF (typechar.EQ.'N ') THEN
         ans=28
      ELSE IF (typechar.EQ.'NA') THEN
         ans=29
      ELSE IF (typechar.EQ.'NB') THEN
         ans=30
      ELSE IF (typechar.EQ.'NC') THEN
         ans=31
      ELSE IF (typechar.EQ.'N2') THEN
         ans=32
      ELSE IF (typechar.EQ.'N3') THEN
         ans=33
      ELSE IF (typechar.EQ.'N*') THEN
         ans=34
      ELSE IF (typechar.EQ.'O ') THEN
         ans=35
      ELSE IF (typechar.EQ.'OW') THEN
         ans=36
      ELSE IF (typechar.EQ.'OH') THEN
         ans=37
      ELSE IF (typechar.EQ.'OS') THEN
         ans=38
      ELSE IF (typechar.EQ.'O2') THEN
         ans=39
      ELSE IF (typechar.EQ.'P ') THEN
         ans=40
      ELSE IF (typechar.EQ.'S ') THEN
         ans=41
      ELSE IF (typechar.EQ.'SH') THEN
         ans=42
      ELSE 
         ans=-1
      ENDIF

      RETURN
      END
      SUBROUTINE secondderivs
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION      VIXEN1,VIXEN2,VIXEN3,VIXEN4,VIXEN5,VIXEN6,VIXEN7,VIXEN8,VIXEN9,JKL
      DOUBLE PRECISION      TORSCRAP(3*NATOMS)
      DOUBLE PRECISION      MUCRAP, LAMBDACRAP, RBC, RBC2, RBC4, RBC6, XE, YE, ZE, RAB, RAB5, RAB2, RAB3, RAB6
      DOUBLE PRECISION      KBOND, RBOND, RAB4, ST, CT, XAB2, YAB2, ZAB2, XAB, YAB, ZAB

      PRINT *,'In secondderivs' 
      hellcount= hellcount+1
C
C Initialise!!
C
      DO a=1,3*atoms
         torscrap(a)=0.0D0
         DO b=1,3*atoms
            bondhell(a,b)=0.0D0
            anglehell(a,b)=0.0D0
            torshell(a,b)=0.0D0
            imphell(a,b)=0.0D0
            qhell(a,b)=0.0D0
            vdwhell(a,b)=0.0D0
         END DO
      END DO

      DO a=1,atoms
         DO b=a+1,atoms
            rab2=r(a,b)**2
            rab3=r(a,b)**3
            rab4=rab2**2
            rab5=rab2*rab3
            rab6=rab3**2
            xab=x(a)-x(b)
            yab=y(a)-y(b)
            zab=z(a)-z(b)
            xab2=xab**2
            yab2=yab**2
            zab2=zab**2

            IF (bonds(a,b).NE.1) GOTO 201

            kbond=kr(type(a),type(b))
            rbond=ro(type(a),type(b))
            vixen1=rbond/r(a,b)
            vixen2=((x(a)-x(b))**2)*rbond/(rab3)
            jkl= 2.0D0*kbond*(1.0D0-vixen1+vixen2)
C
C dxadxa,dxbdxb,dxadxb,dxbdxa
C
            bondhell(3*a-2,3*a-2)= bondhell(3*a-2,3*a-2)+ jkl
            bondhell(3*b-2,3*b-2)= bondhell(3*b-2,3*b-2)+ jkl
            bondhell(3*a-2,3*b-2)= bondhell(3*a-2,3*b-2)- jkl 
            bondhell(3*b-2,3*a-2)= bondhell(3*b-2,3*a-2)- jkl 

            vixen2=((y(a)-y(b))**2)*rbond/(rab3)
            jkl= 2.0D0*kbond*(1.0D0-vixen1+vixen2)
C
C dyadya,dybdyb,dyadyb,dybdya
C
            bondhell(3*a-1,3*a-1)= bondhell(3*a-1,3*a-1)+ jkl
            bondhell(3*b-1,3*b-1)= bondhell(3*b-1,3*b-1)+ jkl
            bondhell(3*a-1,3*b-1)= bondhell(3*a-1,3*b-1)- jkl
            bondhell(3*b-1,3*a-1)= bondhell(3*b-1,3*a-1)- jkl

            vixen2=((z(a)-z(b))**2)*rbond/(rab3)
            jkl= 2.0D0*kbond*(1.0D0-vixen1+vixen2)
C
C dzadza,dzbdzb,dzadzb,dzbdza
C
            bondhell(3*a,3*a)= bondhell(3*a,3*a)+ jkl
            bondhell(3*b,3*b)= bondhell(3*b,3*b)+ jkl
            bondhell(3*a,3*b)= bondhell(3*a,3*b)- jkl
            bondhell(3*b,3*a)= bondhell(3*b,3*a)- jkl

            vixen1=2.0D0*kbond*(x(a)-x(b))*(y(a)-y(b))*rbond/rab3
C
C dxadya,dyadxa,dxbdyb,dybdxb
C
            bondhell(3*a-2,3*a-1)= bondhell(3*a-2,3*a-1)+ vixen1
            bondhell(3*a-1,3*a-2)= bondhell(3*a-1,3*a-2)+ vixen1
            bondhell(3*b-2,3*b-1)= bondhell(3*b-2,3*b-1)+ vixen1
            bondhell(3*b-1,3*b-2)= bondhell(3*b-1,3*b-2)+ vixen1

            vixen1=2.0D0*kbond*(x(a)-x(b))*(z(a)-z(b))*rbond/rab3
C
C dxadza,dzadxa,dxbdzb,dzbdxb
C
            bondhell(3*a-2,3*a)= bondhell(3*a-2,3*a)+ vixen1
            bondhell(3*a,3*a-2)= bondhell(3*a,3*a-2)+ vixen1
            bondhell(3*b-2,3*b)= bondhell(3*b-2,3*b)+ vixen1
            bondhell(3*b,3*b-2)= bondhell(3*b,3*b-2)+ vixen1

            vixen1=2.0D0*kbond*(z(a)-z(b))*(y(a)-y(b))*rbond/rab3
C
C dyadza,dzadya,dybdzb,dzbdyb
C
            bondhell(3*a-1,3*a)= bondhell(3*a-1,3*a)+ vixen1
            bondhell(3*a,3*a-1)= bondhell(3*a,3*a-1)+ vixen1
            bondhell(3*b-1,3*b)= bondhell(3*b-1,3*b)+ vixen1
            bondhell(3*b,3*b-1)= bondhell(3*b,3*b-1)+ vixen1

            vixen1=-2.0D0*kbond*(x(a)-x(b))*(y(a)-y(b))*rbond/rab3
C
C dxadyb,dybdxa,dxbdya,dyadxb
C
            bondhell(3*a-2,3*b-1)= bondhell(3*a-2,3*b-1)+ vixen1
            bondhell(3*b-1,3*a-2)= bondhell(3*b-1,3*a-2)+ vixen1
            bondhell(3*b-2,3*a-1)= bondhell(3*b-2,3*a-1)+ vixen1
            bondhell(3*a-1,3*b-2)= bondhell(3*a-1,3*b-2)+ vixen1

            vixen1=-2.0D0*kbond*(x(a)-x(b))*(z(a)-z(b))*rbond/rab3
C
C dxadzb,dzbdxa,dxbdza,dzadxb
C
            bondhell(3*a-2,3*b)= bondhell(3*a-2,3*b)+ vixen1
            bondhell(3*b,3*a-2)= bondhell(3*b,3*a-2)+ vixen1
            bondhell(3*b-2,3*a)= bondhell(3*b-2,3*a)+ vixen1
            bondhell(3*a,3*b-2)= bondhell(3*a,3*b-2)+ vixen1

            vixen1=-2.0D0*kbond*(y(a)-y(b))*(z(a)-z(b))*rbond/rab3
C
C dyadzb,dzbdya,dybdza,dzadyb
C
            bondhell(3*a-1,3*b)= bondhell(3*a-1,3*b)+ vixen1
            bondhell(3*b,3*a-1)= bondhell(3*b,3*a-1)+ vixen1
            bondhell(3*b-1,3*a)= bondhell(3*b-1,3*a)+ vixen1
            bondhell(3*a,3*b-1)= bondhell(3*a,3*b-1)+ vixen1

201      CONTINUE
C
C Electrostatic Derivatives
C
            IF (FAKEWATER) THEN
               vixen2=(1.0D0-DBLE(bonds(a,b)+one_three(a,b)))/(1.0D0+0.2D0*one_four(a,b))
C
C dxadxa,dxbdxb,dxadxb,dxbdxa
C
               qhell(3*a-2,3*a-2)= qhell(3*a-2,3*a-2)+ vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*xab2/rab2-1.0D0)/rab4
               qhell(3*b-2,3*b-2)= qhell(3*b-2,3*b-2)+ vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*xab2/rab2-1.0D0)/rab4
               qhell(3*a-2,3*b-2)= qhell(3*a-2,3*b-2)- vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*xab2/rab2-1.0D0)/rab4
               qhell(3*b-2,3*a-2)= qhell(3*b-2,3*a-2)- vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*xab2/rab2-1.0D0)/rab4
C
C dyadya,dybdyb,dyadyb,dybdya
C
               qhell(3*a-1,3*a-1)= qhell(3*a-1,3*a-1)+ vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*yab2/rab2-1.0D0)/rab4
               qhell(3*b-1,3*b-1)= qhell(3*b-1,3*b-1)+ vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*yab2/rab2-1.0D0)/rab4
               qhell(3*a-1,3*b-1)= qhell(3*a-1,3*b-1)- vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*yab2/rab2-1.0D0)/rab4
               qhell(3*b-1,3*a-1)= qhell(3*b-1,3*a-1)- vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*yab2/rab2-1.0D0)/rab4
C
C dzadza,dzbdzb,dzadzb,dzbdza
C
               qhell(3*a,3*a)= qhell(3*a,3*a)+ vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*zab2/rab2-1.0D0)/rab4
               qhell(3*b,3*b)= qhell(3*b,3*b)+ vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*zab2/rab2-1.0D0)/rab4
               qhell(3*a,3*b)= qhell(3*a,3*b)- vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*zab2/rab2-1.0D0)/rab4
               qhell(3*b,3*a)= qhell(3*b,3*a)- vixen2*2.0D0*pq(a)*pq(b)*(4.0D0*zab2/rab2-1.0D0)/rab4
C
C dxadya,dyadxa,dxbdyb,dybdxb
C
               qhell(3*a-2,3*a-1)= qhell(3*a-2,3*a-1)+ vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
               qhell(3*a-1,3*a-2)= qhell(3*a-1,3*a-2)+ vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
               qhell(3*b-2,3*b-1)= qhell(3*b-2,3*b-1)+ vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
               qhell(3*b-1,3*b-2)= qhell(3*b-1,3*b-2)+ vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
C
C dxadyb,dybdxa,dxbdya,dyadxb
C
               qhell(3*a-2,3*b-1)= qhell(3*a-2,3*b-1)- vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
               qhell(3*b-1,3*a-2)= qhell(3*b-1,3*a-2)- vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
               qhell(3*b-2,3*a-1)= qhell(3*b-2,3*a-1)- vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
               qhell(3*a-1,3*b-2)= qhell(3*a-1,3*b-2)- vixen2*8.0D0*pq(a)*pq(b)*xab*yab/rab6
C
C dxadza,dzadxa,dxbdzb,dzbdxb
C
               qhell(3*a-2,3*a)= qhell(3*a-2,3*a)+ vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
               qhell(3*a,3*a-2)= qhell(3*a,3*a-2)+ vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
               qhell(3*b-2,3*b)= qhell(3*b-2,3*b)+ vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
               qhell(3*b,3*b-2)= qhell(3*b,3*b-2)+ vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
C
C dxadzb,dzbdxa,dxbdza,dzadxb
C
               qhell(3*a-2,3*b)= qhell(3*a-2,3*b)- vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
               qhell(3*b,3*a-2)= qhell(3*b,3*a-2)- vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
               qhell(3*b-2,3*a)= qhell(3*b-2,3*a)- vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
               qhell(3*a,3*b-2)= qhell(3*a,3*b-2)- vixen2*8.0D0*pq(a)*pq(b)*xab*zab/rab6
C
C dzadya,dyadza,dzbdyb,dybdzb
C
               qhell(3*a,3*a-1)= qhell(3*a,3*a-1)+ vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
               qhell(3*a-1,3*a)= qhell(3*a-1,3*a)+ vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
               qhell(3*b,3*b-1)= qhell(3*b,3*b-1)+ vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
               qhell(3*b-1,3*b)= qhell(3*b-1,3*b)+ vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
C
C dzadyb,dybdza,dzbdya,dyadzb
C
               qhell(3*a,3*b-1)= qhell(3*a,3*b-1)- vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
               qhell(3*b-1,3*a)= qhell(3*b-1,3*a)- vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
               qhell(3*b,3*a-1)= qhell(3*b,3*a-1)- vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6
               qhell(3*a-1,3*b)= qhell(3*a-1,3*b)- vixen2*8.0D0*pq(a)*pq(b)*zab*yab/rab6


            ELSE
               vixen2= (1.0D0-DBLE(bonds(a,b)+one_three(a,b)))*pq(a)*pq(b)/
     1                 ((1.0D0+0.2D0*one_four(a,b))*dielec*rab5)
C
C dxadxa,dxbdxb,dxadxb,dxbdxa
C
               vixen1= (3.0D0*xab**2 - rab2)*vixen2
               qhell(3*a-2,3*a-2)= qhell(3*a-2,3*a-2)+ vixen1
               qhell(3*b-2,3*b-2)= qhell(3*b-2,3*b-2)+ vixen1
               qhell(3*a-2,3*b-2)= qhell(3*a-2,3*b-2)- vixen1
               qhell(3*b-2,3*a-2)= qhell(3*b-2,3*a-2)- vixen1
C
C dyadya,dybdyb,dyadyb,dybdya
C
               vixen1= (3.0D0*yab**2 - rab2)*vixen2
               qhell(3*a-1,3*a-1)= qhell(3*a-1,3*a-1)+ vixen1
               qhell(3*b-1,3*b-1)= qhell(3*b-1,3*b-1)+ vixen1
               qhell(3*a-1,3*b-1)= qhell(3*a-1,3*b-1)- vixen1
               qhell(3*b-1,3*a-1)= qhell(3*b-1,3*a-1)- vixen1
C
C dzadza,dzbdzb,dzadzb,dzbdza
C
               vixen1= (3.0D0*zab**2 - rab2)*vixen2
               qhell(3*a,3*a)= qhell(3*a,3*a)+ vixen1
               qhell(3*b,3*b)= qhell(3*b,3*b)+ vixen1
               qhell(3*a,3*b)= qhell(3*a,3*b)- vixen1
               qhell(3*b,3*a)= qhell(3*b,3*a)- vixen1

               vixen7= vixen2*3.0D0*xab*yab
               vixen8= vixen2*3.0D0*xab*zab
               vixen9= vixen2*3.0D0*yab*zab
C
C dxadya,dyadxa,dxbdyb,dybdxb
C
               qhell(3*a-2,3*a-1)= qhell(3*a-2,3*a-1)+ vixen7
               qhell(3*a-1,3*a-2)= qhell(3*a-1,3*a-2)+ vixen7
               qhell(3*b-2,3*b-1)= qhell(3*b-2,3*b-1)+ vixen7
               qhell(3*b-1,3*b-2)= qhell(3*b-1,3*b-2)+ vixen7
C
C dxadza,dzadxa,dxbdzb,dzbdxb
C
               qhell(3*a-2,3*a)= qhell(3*a-2,3*a)+ vixen8
               qhell(3*a,3*a-2)= qhell(3*a,3*a-2)+ vixen8
               qhell(3*b-2,3*b)= qhell(3*b-2,3*b)+ vixen8
               qhell(3*b,3*b-2)= qhell(3*b,3*b-2)+ vixen8
C
C dyadza,dzadya,dybdzb,dzbdyb
C
               qhell(3*a-1,3*a)= qhell(3*a-1,3*a)+ vixen9
               qhell(3*a,3*a-1)= qhell(3*a,3*a-1)+ vixen9
               qhell(3*b-1,3*b)= qhell(3*b-1,3*b)+ vixen9
               qhell(3*b,3*b-1)= qhell(3*b,3*b-1)+ vixen9
C
C dxadyb,dybdxa,dxbdya,dyadxb
C
               qhell(3*a-2,3*b-1)= qhell(3*a-2,3*b-1)- vixen7
               qhell(3*b-1,3*a-2)= qhell(3*b-1,3*a-2)- vixen7
               qhell(3*b-2,3*a-1)= qhell(3*b-2,3*a-1)- vixen7
               qhell(3*a-1,3*b-2)= qhell(3*a-1,3*b-2)- vixen7
C
C dxadzb,dzbdxa,dxbdza,dzadxb
C
               qhell(3*a-2,3*b)= qhell(3*a-2,3*b)- vixen8
               qhell(3*b,3*a-2)= qhell(3*b,3*a-2)- vixen8
               qhell(3*b-2,3*a)= qhell(3*b-2,3*a)- vixen8
               qhell(3*a,3*b-2)= qhell(3*a,3*b-2)- vixen8
C
C dyadzb,dzbdya,dybdza,dzadyb
C            
               qhell(3*a-1,3*b)= qhell(3*a-1,3*b)- vixen9
               qhell(3*b,3*a-1)= qhell(3*b,3*a-1)- vixen9
               qhell(3*b-1,3*a)= qhell(3*b-1,3*a)- vixen9
               qhell(3*a,3*b-1)= qhell(3*a,3*b-1)- vixen9
            END IF

C
C VdW derivatives
C
            vixen2= (1.0D0-DBLE(bonds(a,b)+one_three(a,b)))/(1.0D0+DBLE(one_four(a,b)))
C
C dxadxa,dxbdxb,dxadxb,dxbdxa
C
            vixen1=24.0D0*(7.0D0*(vdwa(a,b)/rab4**4)-2.0D0*(vdwb(a,b)/(rab4**2*rab2)))
            vixen3=6.0D0*(vdwb(a,b)/rab4**2-2.0D0*(vdwa(a,b)/(rab4**3*rab2)))
            jkl= vixen2*(vixen1*xab2+vixen3)
            vdwhell(3*a-2,3*a-2)= vdwhell(3*a-2,3*a-2)+ jkl
            vdwhell(3*b-2,3*b-2)= vdwhell(3*b-2,3*b-2)+ jkl
            vdwhell(3*a-2,3*b-2)= vdwhell(3*a-2,3*b-2)- jkl
            vdwhell(3*b-2,3*a-2)= vdwhell(3*b-2,3*a-2)- jkl
C
C dyadya,dybdyb,dyadyb,dybdya
C
            jkl= vixen2*(vixen1*yab2+vixen3)
            vdwhell(3*a-1,3*a-1)= vdwhell(3*a-1,3*a-1)+ jkl
            vdwhell(3*b-1,3*b-1)= vdwhell(3*b-1,3*b-1)+ jkl
            vdwhell(3*a-1,3*b-1)= vdwhell(3*a-1,3*b-1)- jkl
            vdwhell(3*b-1,3*a-1)= vdwhell(3*b-1,3*a-1)- jkl
C
C dzadza,dzbdzb,dzadzb,dzbdza
C
            jkl= vixen2*(vixen1*zab2+vixen3)
            vdwhell(3*a,3*a)= vdwhell(3*a,3*a)+ jkl
            vdwhell(3*b,3*b)= vdwhell(3*b,3*b)+ jkl
            vdwhell(3*a,3*b)= vdwhell(3*a,3*b)- jkl
            vdwhell(3*b,3*a)= vdwhell(3*b,3*a)- jkl
C
C dxadya,dyadxa,dxbdyb,dybdxb
C
            jkl= vixen2*vixen1*yab*xab
            vdwhell(3*a-2,3*a-1)= vdwhell(3*a-2,3*a-1)+ jkl
            vdwhell(3*a-1,3*a-2)= vdwhell(3*a-1,3*a-2)+ jkl
            vdwhell(3*b-2,3*b-1)= vdwhell(3*b-2,3*b-1)+ jkl
            vdwhell(3*b-1,3*b-2)= vdwhell(3*b-1,3*b-2)+ jkl
C
C dxadyb,dybdxa,dxbdya,dyadxb
C
            vdwhell(3*a-2,3*b-1)= vdwhell(3*a-2,3*b-1)- jkl
            vdwhell(3*b-1,3*a-2)= vdwhell(3*b-1,3*a-2)- jkl
            vdwhell(3*b-2,3*a-1)= vdwhell(3*b-2,3*a-1)- jkl
            vdwhell(3*a-1,3*b-2)= vdwhell(3*a-1,3*b-2)- jkl
C
C dxadza,dzadxa,dxbdzb,dzbdxb
C
            jkl= vixen2*vixen1*xab*zab
            vdwhell(3*a-2,3*a)= vdwhell(3*a-2,3*a)+ jkl
            vdwhell(3*a,3*a-2)= vdwhell(3*a,3*a-2)+ jkl
            vdwhell(3*b-2,3*b)= vdwhell(3*b-2,3*b)+ jkl
            vdwhell(3*b,3*b-2)= vdwhell(3*b,3*b-2)+ jkl
C
C dxadzb,dzbdxa,dxbdza,dzadxb
C
            vdwhell(3*a-2,3*b)= vdwhell(3*a-2,3*b)- jkl
            vdwhell(3*b,3*a-2)= vdwhell(3*b,3*a-2)- jkl
            vdwhell(3*b-2,3*a)= vdwhell(3*b-2,3*a)- jkl
            vdwhell(3*a,3*b-2)= vdwhell(3*a,3*b-2)- jkl
C
C dyadza,dzadya,dybdzb,dzbdyb
C
            jkl= vixen2*vixen1*yab*zab
            vdwhell(3*a-1,3*a)= vdwhell(3*a-1,3*a)+ jkl
            vdwhell(3*a,3*a-1)= vdwhell(3*a,3*a-1)+ jkl
            vdwhell(3*b-1,3*b)= vdwhell(3*b-1,3*b)+ jkl
            vdwhell(3*b,3*b-1)= vdwhell(3*b,3*b-1)+ jkl
C
C dyadzb,dzbdya,dybdza,dzadyb
C
            vdwhell(3*a-1,3*b)= vdwhell(3*a-1,3*b)- jkl
            vdwhell(3*b,3*a-1)= vdwhell(3*b,3*a-1)- jkl
            vdwhell(3*b-1,3*a)= vdwhell(3*b-1,3*a)- jkl
            vdwhell(3*a,3*b-1)= vdwhell(3*a,3*b-1)- jkl
         END DO
      END DO

      CALL anglesecderivs

      RETURN

      END

   
      SUBROUTINE anglesecderivs
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE

      DOUBLE PRECISION      VIXEN1,VIXEN2,VIXEN3,VIXEN4,VIXEN5,VIXEN6,VIXEN7,VIXEN8,VIXEN9
      DOUBLE PRECISION      U,DTHETAX1,DTHETAY1,DTHETAZ1,DTHETAX2,DTHETAY2,DTHETAZ2,DTHETAX3,DTHETAY3,DTHETAZ3,KTHETA,THETA0
      DOUBLE PRECISION      V,XP,YP,ZP,XPP,YPP,ZPP,DLAMBDAX1,DLAMBDAY1,DLAMBDAZ1,DLAMBDAX2,DLAMBDAY2,DLAMBDAZ2
      DOUBLE PRECISION      DLAMBDAX3,DLAMBDAY3,DLAMBDAZ3,DMUX2,DMUY2,DMUZ2,DMUX3,DMUY3,DMUZ3,DMUX4,DMUY4,DMUZ4
      DOUBLE PRECISION      TORSCRAP(3*NATOMS)
      DOUBLE PRECISION      MUCRAP, LAMBDACRAP, RBC, RBC2, RBC4, RBC6, XE, YE, ZE, RAB, RAB5, RAB2, RAB3
      DOUBLE PRECISION      KBOND, RBOND, RAB4, ST, CT
      DOUBLE PRECISION      DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      COMMON /DU/           dux1,dux2,dux3,dux4,duy1,duy2,duy3,duy4,duz1,duz2,duz3,duz4
      DOUBLE PRECISION      DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      COMMON /DV/           dvx1,dvx2,dvx3,dvx4,dvy1,dvy2,dvy3,dvy4,dvz1,dvz2,dvz3,dvz4
      DOUBLE PRECISION      DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      COMMON /DM/           dmx1,dmx2,dmx3,dmx4,dmy1,dmy2,dmy3,dmy4,dmz1,dmz2,dmz3,dmz4
      DOUBLE PRECISION      DPX1,DPX2,DPX3,DPY1,DPY2,DPY3,DPZ1,DPZ2,DPZ3
      COMMON /DP/           dpx1,dpx2,dpx3,dpy1,dpy2,dpy3,dpz1,dpz2,dpz3
      DOUBLE PRECISION      DNX1,DNX2,DNX3,DNX4,DNY1,DNY2,DNY3,DNY4,DNZ1,DNZ2,DNZ3,DNZ4
      COMMON /DN/           dnx1,dnx2,dnx3,dnx4,dny1,dny2,dny3,dny4,dnz1,dnz2,dnz3,dnz4
      DOUBLE PRECISION      DQX2,DQX3,DQX4,DQY2,DQY3,DQY4,DQZ2,DQZ3,DQZ4
      COMMON /DQ/           dqx2,dqx3,dqx4,dqy2,dqy3,dqy4,dqz2,dqz3,dqz4
      DOUBLE PRECISION      M,N
      COMMON /MN/           m,n

C
C ANGLE DERIVATIVES
C *****************
C
      DO i=1,ang
         a=aa1(i)
         b=aa2(i)
         c=aa3(i) 
C
C For comparison with notes, vixen1=m, vixen2=n, vixen3=dm, vixen4=dn
C
         ktheta=kt(type(a),type(b),type(c))
         theta0=to(type(a),type(b),type(c))
         st=SIN(theta(i))
         ct=COS(theta(i))
         u=x(b)*(x(b)-x(a)-x(c))+y(b)*(y(b)-y(a)-y(c))+z(b)*(z(b)-z(a)-z(c))+x(a)*x(c)+y(a)*y(c)+z(a)*z(c)
         vixen1=2.0D0*ktheta*(theta(i)-theta0)/st

         dthetax1=(u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c)))/st
         dthetay1=(u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c)))/st
         dthetaz1=(u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c)))/st

         vixen5= u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a))) / (r(a,b)**3*r(b,c)**3)
         vixen6= (2.0D0*x(b)-x(a)-x(c))/(r(a,b)*r(b,c))
         dthetax2=(vixen5-vixen6)/st

         vixen5=u*(r(a,b)**2*(y(b)-y(c))+r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**3)
         vixen6=(2.0D0*y(b)-y(a)-y(c))/(r(a,b)*r(b,c))
         dthetay2=(vixen5-vixen6)/st

         vixen5=u*(r(a,b)**2*(z(b)-z(c))+r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**3)
         vixen6=(2.0D0*z(b)-z(a)-z(c))/(r(a,b)*r(b,c))
         dthetaz2=(vixen5-vixen6)/st

         dthetax3=(u*(x(c)-x(b))/(r(a,b)*r(b,c)**3)-(x(a)-x(b))/(r(a,b)*r(b,c)))/st
         dthetay3=(u*(y(c)-y(b))/(r(a,b)*r(b,c)**3)-(y(a)-y(b))/(r(a,b)*r(b,c)))/st
         dthetaz3=(u*(z(c)-z(b))/(r(a,b)*r(b,c)**3)-(z(a)-z(b))/(r(a,b)*r(b,c)))/st
C
C dxadxa,dyadya,dzadza
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax1*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=(u+2.0D0*(x(a)-x(b))*(x(c)-x(b))-3.0D0*u*(x(a)-x(b))**2/r(a,b)**2)/(r(a,b)**3*r(b,c))

         anglehell(3*a-2,3*a-2)= anglehell(3*a-2,3*a-2)+ vixen1*vixen4+vixen2*vixen3 

         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay1*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=(u+2.0D0*(y(a)-y(b))*(y(c)-y(b))-3.0D0*u*(y(a)-y(b))**2/r(a,b)**2)/(r(a,b)**3*r(b,c))

         anglehell(3*a-1,3*a-1)= anglehell(3*a-1,3*a-1)+ vixen1*vixen4+vixen2*vixen3

         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz1*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=(u+2.0D0*(z(a)-z(b))*(z(c)-z(b))-3.0D0*u*(z(a)-z(b))**2/r(a,b)**2)/(r(a,b)**3*r(b,c))

         anglehell(3*a,3*a)= anglehell(3*a,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dxadya,dyadxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay1*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=((x(a)-x(b))*(y(c)-y(b))+(x(c)-x(b))*(y(a)-y(b))-3.0D0*u*(x(a)-x(b))*(y(a)-y(b))/r(a,b)**2)/(r(a,b)**3*r(b,c))
 
         anglehell(3*a-2,3*a-1)= anglehell(3*a-2,3*a-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*a-1,3*a-2)= anglehell(3*a-1,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dxadza,dzadxa
C
         vixen3=2.0D0*ktheta*dthetaz1*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=((x(a)-x(b))*(z(c)-z(b))+(x(c)-x(b))*(z(a)-z(b))-3.0D0*u*(x(a)-x(b))*(z(a)-z(b))/r(a,b)**2)/(r(a,b)**3*r(b,c))

         anglehell(3*a-2,3*a)= anglehell(3*a-2,3*a)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*a,3*a-2)= anglehell(3*a,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dyadza,dzadya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz1*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=((y(a)-y(b))*(z(c)-z(b))+(y(c)-y(b))*(z(a)-z(b))-3.0D0*u*(y(a)-y(b))*(z(a)-z(b))/r(a,b)**2)/(r(a,b)**3*r(b,c))

         anglehell(3*a-1,3*a)= anglehell(3*a-1,3*a)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*a,3*a-1)= anglehell(3*a,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dxadxb,dxbdxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=((x(a)-x(b))*(2.0D0*x(b)-x(a)-x(c))-u)/(r(a,b)**3*r(b,c)) + 1.0D0/(r(a,b)*r(b,c))
         vixen6=(u*(x(a)-x(b))*(r(a,b)**3*(x(b)-x(c)) + 3.0D0*r(b,c)**2*r(a,b)*(x(b)-x(a))))/(r(a,b)**6*r(b,c)**3)
         vixen7=(x(c)-x(b))*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**3)
         vixen4= vixen5 - vixen6 + vixen7

         anglehell(3*a-2,3*b-2)= anglehell(3*a-2,3*b-2)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-2,3*a-2)= anglehell(3*b-2,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dyadyb,dybdya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=((y(a)-y(b))*(2.0D0*y(b)-y(a)-y(c))-u)/(r(a,b)**3*r(b,c)) + 1.0D0/(r(a,b)*r(b,c))
         vixen6=(u*(y(a)-y(b))*(r(a,b)**3*(y(b)-y(c)) + 3.0D0*r(b,c)**2*r(a,b)*(y(b)-y(a))))/(r(a,b)**6*r(b,c)**3)
         vixen7=(y(c)-y(b))*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**3)
         vixen4= vixen5 - vixen6 + vixen7

         anglehell(3*a-1,3*b-1)= anglehell(3*a-1,3*b-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*a-1)= anglehell(3*b-1,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dzadzb,dzbdza
C
         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=((z(a)-z(b))*(2.0D0*z(b)-z(a)-z(c))-u)/(r(a,b)**3*r(b,c)) + 1.0D0/(r(a,b)*r(b,c))
         vixen6=(u*(z(a)-z(b))*(r(a,b)**3*(z(b)-z(c)) + 3.0D0*r(b,c)**2*r(a,b)*(z(b)-z(a))))/(r(a,b)**6*r(b,c)**3)
         vixen7=(z(c)-z(b))*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**3)
         vixen4= vixen5 - vixen6 + vixen7

         anglehell(3*a,3*b)= anglehell(3*a,3*b)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*a)= anglehell(3*b,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dxadyb,dybdxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(x(a)-x(b))*(2.0D0*y(b)-y(a)-y(c))/(r(a,b)**3*r(b,c)) 
         vixen6=u*(x(a)-x(b))*(r(a,b)**3*(y(b)-y(c)) + 3.0D0*r(a,b)*r(b,c)**2*(y(b)-y(a)))/(r(a,b)**6*r(b,c)**3)
         vixen4= vixen5 - vixen6 + (x(c)-x(b))*((y(b)-y(c))/(r(a,b)*r(b,c)**3) + (y(b)-y(a))/(r(a,b)**3*r(b,c)))

         anglehell(3*a-2,3*b-1)= anglehell(3*a-2,3*b-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*a-2)= anglehell(3*b-1,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dxadzb,dzbdxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(x(a)-x(b))*(2.0D0*z(b)-z(a)-z(c))/(r(a,b)**3*r(b,c)) 
         vixen6=u*(x(a)-x(b))*(r(a,b)**3*(z(b)-z(c))+3.0D0*r(a,b)*r(b,c)**2*(z(b)-z(a)))/(r(a,b)**6*r(b,c)**3)
         vixen4= vixen5 - vixen6 + (x(c)-x(b))*((z(b)-z(c))/(r(a,b)*r(b,c)**3) + (z(b)-z(a))/(r(a,b)**3*r(b,c)))

         anglehell(3*a-2,3*b)= anglehell(3*a-2,3*b)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*a-2)= anglehell(3*b,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dyadxb,dxbdya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(a)-y(b))*(2.0D0*x(b)-x(a)-x(c))/(r(a,b)**3*r(b,c))
         vixen6=u*(y(a)-y(b))*(r(a,b)**3*(x(b)-x(c))+3.0D0*r(a,b)*r(b,c)**2*(x(b)-x(a)))/(r(a,b)**6*r(b,c)**3)
         vixen4= vixen5 - vixen6 + (y(c)-y(b))*((x(b)-x(c))/(r(a,b)*r(b,c)**3) + (x(b)-x(a))/(r(a,b)**3*r(b,c)))

         anglehell(3*a-1,3*b-2)= anglehell(3*a-1,3*b-2)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-2,3*a-1)= anglehell(3*b-2,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dyadzb,dzbdya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(a)-y(b))*(2.0D0*z(b)-z(a)-z(c))/(r(a,b)**3*r(b,c))
         vixen6=u*(y(a)-y(b))*(r(a,b)**3*(z(b)-z(c))+3.0D0*r(a,b)*r(b,c)**2*(z(b)-z(a)))/(r(a,b)**6*r(b,c)**3)
         vixen4= vixen5 - vixen6 + (y(c)-y(b))*((z(b)-z(c))/(r(a,b)*r(b,c)**3) + (z(b)-z(a))/(r(a,b)**3*r(b,c)))

         anglehell(3*a-1,3*b)= anglehell(3*a-1,3*b)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*a-1)= anglehell(3*b,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dzadxb,dxbdza
C
         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(a)-z(b))*(2.0D0*x(b)-x(a)-x(c))/(r(a,b)**3*r(b,c))
         vixen6=u*(z(a)-z(b))*(r(a,b)**3*(x(b)-x(c))+3.0D0*r(a,b)*r(b,c)**2*(x(b)-x(a)))/(r(a,b)**6*r(b,c)**3)
         vixen4= vixen5 - vixen6 + (z(c)-z(b))*((x(b)-x(c))/(r(a,b)*r(b,c)**3) + (x(b)-x(a))/(r(a,b)**3*r(b,c)))

         anglehell(3*a,3*b-2)= anglehell(3*a,3*b-2)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-2,3*a)= anglehell(3*b-2,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dzadyb,dybdza
C
         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(a)-z(b))*(2.0D0*y(b)-y(a)-y(c))/(r(a,b)**3*r(b,c))
         vixen6=u*(z(a)-z(b))*(r(a,b)**3*(y(b)-y(c))+3.0D0*r(a,b)*r(b,c)**2*(y(b)-y(a)))/(r(a,b)**6*r(b,c)**3)
         vixen4= vixen5 - vixen6 + (z(c)-z(b))*((y(b)-y(c))/(r(a,b)*r(b,c)**3) + (y(b)-y(a))/(r(a,b)**3*r(b,c)))

         anglehell(3*a,3*b-1)= anglehell(3*a,3*b-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*a)= anglehell(3*b-1,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dxadxc,dxcdxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(x(a)-x(b))**2/(r(a,b)**3*r(b,c)) + (x(c)-x(b))**2/(r(a,b)*r(b,c)**3) - 1.0D0/(r(a,b)*r(b,c))
         vixen4= vixen5 - u*(x(a)-x(b))*(x(c)-x(b))/(r(a,b)**3*r(b,c)**3)

         anglehell(3*a-2,3*c-2)= anglehell(3*a-2,3*c-2)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-2,3*a-2)= anglehell(3*c-2,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dyadyc,dycdya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(a)-y(b))**2/(r(a,b)**3*r(b,c)) + (y(c)-y(b))**2/(r(a,b)*r(b,c)**3) - 1.0D0/(r(a,b)*r(b,c))
         vixen4= vixen5 - u*(y(a)-y(b))*(y(c)-y(b))/(r(a,b)**3*r(b,c)**3)

         anglehell(3*a-1,3*c-1)= anglehell(3*a-1,3*c-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-1,3*a-1)= anglehell(3*c-1,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dzadzc,dzcdza
C
         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(a)-z(b))**2/(r(a,b)**3*r(b,c)) + (z(c)-z(b))**2/(r(a,b)*r(b,c)**3) - 1.0D0/(r(a,b)*r(b,c))
         vixen4= vixen5 - u*(z(a)-z(b))*(z(c)-z(b))/(r(a,b)**3*r(b,c)**3)

         anglehell(3*a,3*c)= anglehell(3*a,3*c)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c,3*a)= anglehell(3*c,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dxadyc,dycdxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(x(a)-x(b))*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - u*(x(a)-x(b))*(y(c)-y(b))/(r(a,b)**3*r(b,c)**3) 
         vixen4= vixen5 + (x(c)-x(b))*(y(c)-y(b))/(r(a,b)*r(b,c)**3)

         anglehell(3*a-2,3*c-1)= anglehell(3*a-2,3*c-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-1,3*a-2)= anglehell(3*c-1,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dxadzc,dzcdxa
C
         vixen2=u*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - (x(c)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(x(a)-x(b))*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - u*(x(a)-x(b))*(z(c)-z(b))/(r(a,b)**3*r(b,c)**3) 
         vixen4= vixen5 + (x(c)-x(b))*(z(c)-z(b))/(r(a,b)*r(b,c)**3)

         anglehell(3*a-2,3*c)= anglehell(3*a-2,3*c)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c,3*a-2)= anglehell(3*c,3*a-2)+ vixen1*vixen4+vixen2*vixen3
C
C dyadzc,dzcdya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(a)-y(b))*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - u*(y(a)-y(b))*(z(c)-z(b))/(r(a,b)**3*r(b,c)**3) 
         vixen4= vixen5 + (y(c)-y(b))*(z(c)-z(b))/(r(a,b)*r(b,c)**3)

         anglehell(3*a-1,3*c)= anglehell(3*a-1,3*c)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c,3*a-1)= anglehell(3*c,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dyadxc,dxcdya
C
         vixen2=u*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - (y(c)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(a)-y(b))*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - u*(y(a)-y(b))*(x(c)-x(b))/(r(a,b)**3*r(b,c)**3) 
         vixen4= vixen5 + (y(c)-y(b))*(x(c)-x(b))/(r(a,b)*r(b,c)**3)

         anglehell(3*a-1,3*c-2)= anglehell(3*a-1,3*c-2)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-2,3*a-1)= anglehell(3*c-2,3*a-1)+ vixen1*vixen4+vixen2*vixen3
C
C dzadxc,dxcdza
C
         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(a)-z(b))*(x(a)-x(b))/(r(a,b)**3*r(b,c)) - u*(z(a)-z(b))*(x(c)-x(b))/(r(a,b)**3*r(b,c)**3) 
         vixen4= vixen5 + (z(c)-z(b))*(x(c)-x(b))/(r(a,b)*r(b,c)**3)

         anglehell(3*a,3*c-2)= anglehell(3*a,3*c-2)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-2,3*a)= anglehell(3*c-2,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dzadyc,dycdza
C
         vixen2=u*(z(a)-z(b))/(r(a,b)**3*r(b,c)) - (z(c)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay3*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(a)-z(b))*(y(a)-y(b))/(r(a,b)**3*r(b,c)) - u*(z(a)-z(b))*(y(c)-y(b))/(r(a,b)**3*r(b,c)**3) 
         vixen4= vixen5 + (z(c)-z(b))*(y(c)-y(b))/(r(a,b)*r(b,c)**3)

         anglehell(3*a,3*c-1)= anglehell(3*a,3*c-1)+ vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-1,3*a)= anglehell(3*c-1,3*a)+ vixen1*vixen4+vixen2*vixen3
C
C dxcdxc
C
         vixen2=u*(x(c)-x(b))/(r(a,b)*r(b,c)**3) - (x(a)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax3*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=(2.0D0*(x(a)-x(b))*(x(c)-x(b)) + u - 3.0D0*u*(x(c)-x(b))**2/r(b,c)**2)/(r(a,b)*r(b,c)**3)

         anglehell(3*c-2,3*c-2)= anglehell(3*c-2,3*c-2) + vixen1*vixen4+vixen2*vixen3
C
C dycdyc
C
         vixen2=u*(y(c)-y(b))/(r(a,b)*r(b,c)**3) - (y(a)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay3*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=(2.0D0*(y(a)-y(b))*(y(c)-y(b)) + u - 3.0D0*u*(y(c)-y(b))**2/r(b,c)**2)/(r(a,b)*r(b,c)**3)

         anglehell(3*c-1,3*c-1)= anglehell(3*c-1,3*c-1) + vixen1*vixen4+vixen2*vixen3
C
C dzcdzc
C
         vixen2=u*(z(c)-z(b))/(r(a,b)*r(b,c)**3) - (z(a)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz3*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=(2.0D0*(z(a)-z(b))*(z(c)-z(b)) + u - 3.0D0*u*(z(c)-z(b))**2/r(b,c)**2)/(r(a,b)*r(b,c)**3)

         anglehell(3*c,3*c)= anglehell(3*c,3*c) + vixen1*vixen4+vixen2*vixen3
C
C dxcdyc,dycdxc
C
         vixen2=u*(x(c)-x(b))/(r(a,b)*r(b,c)**3) - (x(a)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay3*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=((x(c)-x(b))*(y(a)-y(b))+(x(a)-x(b))*(y(c)-y(b)) - 3.0D0*u*(x(c)-x(b))*(y(c)-y(b))/r(b,c)**2)/(r(a,b)*r(b,c)**3)

         anglehell(3*c-2,3*c-1)=anglehell(3*c-2,3*c-1) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*c-1,3*c-2)=anglehell(3*c-1,3*c-2) + vixen1*vixen4+vixen2*vixen3 
C
C dxcdzc,dzcdxc
C
         vixen2=u*(x(c)-x(b))/(r(a,b)*r(b,c)**3) - (x(a)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz3*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=((x(c)-x(b))*(z(a)-z(b))+(x(a)-x(b))*(z(c)-z(b)) - 3.0D0*u*(x(c)-x(b))*(z(c)-z(b))/r(b,c)**2)/(r(a,b)*r(b,c)**3)

         anglehell(3*c-2,3*c)=anglehell(3*c-2,3*c) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*c,3*c-2)=anglehell(3*c,3*c-2) + vixen1*vixen4+vixen2*vixen3
C
C dycdzc,dzcdyc
C
         vixen2=u*(y(c)-y(b))/(r(a,b)*r(b,c)**3) - (y(a)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz3*(st-(theta(i)-theta0)*ct)/st**2
         vixen4=((y(c)-y(b))*(z(a)-z(b))+(y(a)-y(b))*(z(c)-z(b)) - 3.0D0*u*(y(c)-y(b))*(z(c)-z(b))/r(b,c)**2)/(r(a,b)*r(b,c)**3)

         anglehell(3*c-1,3*c)=anglehell(3*c-1,3*c) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*c,3*c-1)=anglehell(3*c,3*c-1) + vixen1*vixen4+vixen2*vixen3
C
C dxcdxb,dxbdxc
C
         vixen2=u*(x(c)-x(b))/(r(a,b)*r(b,c)**3) - (x(a)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=((x(c)-x(b))*(2.0D0*x(b)-x(c)-x(a))-u)/(r(a,b)*r(b,c)**3) + 1.0D0/(r(a,b)*r(b,c))
         vixen6=u*(x(c)-x(b))*(3.0D0*r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (x(a)-x(b))*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**3)

         anglehell(3*c-2,3*b-2)= anglehell(3*c-2,3*b-2) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-2,3*c-2)= anglehell(3*b-2,3*c-2) + vixen1*vixen4+vixen2*vixen3
C
C dycdyb,dybdyc
C
         vixen2=u*(y(c)-y(b))/(r(a,b)*r(b,c)**3) - (y(a)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=((y(c)-y(b))*(2.0D0*y(b)-y(c)-y(a))-u)/(r(a,b)*r(b,c)**3) + 1.0D0/(r(a,b)*r(b,c))
         vixen6=u*(y(c)-y(b))*(3.0D0*r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (y(a)-y(b))*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**3)

         anglehell(3*c-1,3*b-1)= anglehell(3*c-1,3*b-1) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*c-1)= anglehell(3*b-1,3*c-1) + vixen1*vixen4+vixen2*vixen3
C
C dzcdzb,dzbdzc
C
         vixen2=u*(z(c)-z(b))/(r(a,b)*r(b,c)**3) - (z(a)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=((z(c)-z(b))*(2.0D0*z(b)-z(c)-z(a))-u)/(r(a,b)*r(b,c)**3) + 1.0D0/(r(a,b)*r(b,c))
         vixen6=u*(z(c)-z(b))*(3.0D0*r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (z(a)-z(b))*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**3)

         anglehell(3*c,3*b)= anglehell(3*c,3*b) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*c)= anglehell(3*b,3*c) + vixen1*vixen4+vixen2*vixen3
C
C dxcdyb,dybdxc
C
         vixen2=u*(x(c)-x(b))/(r(a,b)*r(b,c)**3) - (x(a)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2 

         vixen5=(x(c)-x(b))*(2.0D0*y(b)-y(c)-y(a))/(r(a,b)*r(b,c)**3) 
         vixen6=u*(x(c)-x(b))*(3.0D0*r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (x(a)-x(b))*((y(b)-y(c))/r(b,c)**2 + (y(b)-y(a))/r(a,b)**2)/(r(a,b)*r(b,c))

         anglehell(3*c-2,3*b-1)= anglehell(3*c-2,3*b-1) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*c-2)= anglehell(3*b-1,3*c-2) + vixen1*vixen4+vixen2*vixen3
C
C dxcdzb,dzbdxc
C
         vixen2=u*(x(c)-x(b))/(r(a,b)*r(b,c)**3) - (x(a)-x(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(x(c)-x(b))*(2.0D0*z(b)-z(c)-z(a))/(r(a,b)*r(b,c)**3) 
         vixen6=u*(x(c)-x(b))*(3.0D0*r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (x(a)-x(b))*((z(b)-z(c))/r(b,c)**2 + (z(b)-z(a))/r(a,b)**2)/(r(a,b)*r(b,c))

         anglehell(3*c-2,3*b)= anglehell(3*c-2,3*b) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*c-2)= anglehell(3*b,3*c-2) + vixen1*vixen4+vixen2*vixen3
C
C dycdzb,dzbdyc
C
         vixen2=u*(y(c)-y(b))/(r(a,b)*r(b,c)**3) - (y(a)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(c)-y(b))*(2.0D0*z(b)-z(c)-z(a))/(r(a,b)*r(b,c)**3) 
         vixen6=u*(y(c)-y(b))*(3.0D0*r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (y(a)-y(b))*((z(b)-z(c))/r(b,c)**2 + (z(b)-z(a))/r(a,b)**2)/(r(a,b)*r(b,c))

         anglehell(3*c-1,3*b)= anglehell(3*c-1,3*b) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*c-1)= anglehell(3*b,3*c-1) + vixen1*vixen4+vixen2*vixen3
C
C dycdxb,dxbdyc
C
         vixen2=u*(y(c)-y(b))/(r(a,b)*r(b,c)**3) - (y(a)-y(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(y(c)-y(b))*(2.0D0*x(b)-x(c)-x(a))/(r(a,b)*r(b,c)**3)
         vixen6=u*(y(c)-y(b))*(3.0D0*r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (y(a)-y(b))*((x(b)-x(c))/r(b,c)**2 + (x(b)-x(a))/r(a,b)**2)/(r(a,b)*r(b,c))

         anglehell(3*c-1,3*b-2)= anglehell(3*c-1,3*b-2) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-2,3*c-1)= anglehell(3*b-2,3*c-1) + vixen1*vixen4+vixen2*vixen3 
C
C dzcdxb,dxbdzc 
C
         vixen2=u*(z(c)-z(b))/(r(a,b)*r(b,c)**3) - (z(a)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(c)-z(b))*(2.0D0*x(b)-x(c)-x(a))/(r(a,b)*r(b,c)**3)
         vixen6=u*(z(c)-z(b))*(3.0D0*r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (z(a)-z(b))*((x(b)-x(c))/r(b,c)**2 + (x(b)-x(a))/r(a,b)**2)/(r(a,b)*r(b,c))

         anglehell(3*c,3*b-2)= anglehell(3*c,3*b-2) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-2,3*c)= anglehell(3*b-2,3*c) + vixen1*vixen4+vixen2*vixen3
C
C dzcdyb,dybdzc 
C
         vixen2=u*(z(c)-z(b))/(r(a,b)*r(b,c)**3) - (z(a)-z(b))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=(z(c)-z(b))*(2.0D0*y(b)-y(c)-y(a))/(r(a,b)*r(b,c)**3)
         vixen6=u*(z(c)-z(b))*(3.0D0*r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**5)
         vixen4= vixen5 - vixen6 + (z(a)-z(b))*((y(b)-y(c))/r(b,c)**2 + (y(b)-y(a))/r(a,b)**2)/(r(a,b)*r(b,c))

         anglehell(3*c,3*b-1)= anglehell(3*c,3*b-1) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*c)= anglehell(3*b-1,3*c) + vixen1*vixen4+vixen2*vixen3
C
C dxbdxb
C
         vixen2=u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**3) - (2.0D0*x(b)-x(c)-x(a))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetax2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=u*(4.0D0*(x(b)-x(c))*(x(b)-x(a))+r(a,b)**2+r(b,c)**2) + 
     1          (2.0D0*x(b)-x(c)-x(a))*((x(b)-x(c))*r(a,b)**2+(x(b)-x(a))*r(b,c)**2)
         vixen6=3.0D0*r(a,b)*r(b,c)*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))
         vixen7=u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))
         vixen8= (vixen5*r(a,b)**3*r(b,c)**3 - vixen7*vixen6)/(r(a,b)**6*r(b,c)**6) 
         vixen9= (2.0D0*x(b)-x(c)-x(a))*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**3)
         vixen4= vixen8 - 2.0D0/(r(a,b)*r(b,c)) + vixen9

         anglehell(3*b-2,3*b-2)= anglehell(3*b-2,3*b-2) + vixen1*vixen4+vixen2*vixen3
C
C dybdyb
C
         vixen2=u*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**3) - (2.0D0*y(b)-y(c)-y(a))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=u*(4.0D0*(y(b)-y(c))*(y(b)-y(a))+r(a,b)**2+r(b,c)**2) + 
     1          (2.0D0*y(b)-y(c)-y(a))*((y(b)-y(c))*r(a,b)**2+(y(b)-y(a))*r(b,c)**2)
         vixen6=3.0D0*r(a,b)*r(b,c)*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))
         vixen7=u*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))
         vixen8= (vixen5*r(a,b)**3*r(b,c)**3 - vixen7*vixen6)/(r(a,b)**6*r(b,c)**6)
         vixen9= (2.0D0*y(b)-y(c)-y(a))*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**3)
         vixen4= vixen8 - 2.0D0/(r(a,b)*r(b,c)) + vixen9

         anglehell(3*b-1,3*b-1)= anglehell(3*b-1,3*b-1) + vixen1*vixen4+vixen2*vixen3
C
C dzbdzb
C
         vixen2=u*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**3) - (2.0D0*z(b)-z(c)-z(a))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=u*(4.0D0*(z(b)-z(c))*(z(b)-z(a))+r(a,b)**2+r(b,c)**2) + 
     1          (2.0D0*z(b)-z(c)-z(a))*((z(b)-z(c))*r(a,b)**2+(z(b)-z(a))*r(b,c)**2)
         vixen6=3.0D0*r(a,b)*r(b,c)*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))
         vixen7=u*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))
         vixen8= (vixen5*r(a,b)**3*r(b,c)**3 - vixen7*vixen6)/(r(a,b)**6*r(b,c)**6)
         vixen9= (2.0D0*z(b)-z(c)-z(a))*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))/(r(a,b)**3*r(b,c)**3)
         vixen4= vixen8 - 2.0D0/(r(a,b)*r(b,c)) + vixen9

         anglehell(3*b,3*b)= anglehell(3*b,3*b) + vixen1*vixen4+vixen2*vixen3
C
C dxbdyb,dybdxb
C
         vixen2=u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**3) - (2.0D0*x(b)-x(c)-x(a))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetay2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=u*(2.0D0*(y(b)-y(a))*(x(b)-x(c))+2.0D0*(y(b)-y(c))*(x(b)-x(a)))
         vixen5= vixen5 + (2.0D0*y(b)-y(c)-y(a))*(r(a,b)**2*(x(b)-x(c))+r(b,c)**2*(x(b)-x(a)))
C
C Above is just a trick as otherwise get 'unbalanced parentheses' error message...
C
         vixen6=3.0D0*r(a,b)*r(b,c)*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))
         vixen7=u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))
         vixen8= (vixen5*r(a,b)**3*r(b,c)**3 - vixen6*vixen7)/(r(a,b)**6*r(b,c)**6)
         vixen9=(2.0D0*x(b)-x(c)-x(a))*((y(b)-y(c))/r(b,c)**2 + (y(b)-y(a))/r(a,b)**2)/(r(a,b)*r(b,c))
         vixen4= vixen8 + vixen9

         anglehell(3*b-2,3*b-1)= anglehell(3*b-2,3*b-1) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b-1,3*b-2)= anglehell(3*b-1,3*b-2) + vixen1*vixen4+vixen2*vixen3
C
C dxbdzb,dzbdxb
C
         vixen2=u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))/(r(a,b)**3*r(b,c)**3) - (2.0D0*x(b)-x(c)-x(a))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=u*(2.0D0*(z(b)-z(a))*(x(b)-x(c))+2.0D0*(z(b)-z(c))*(x(b)-x(a)))
         vixen5= vixen5 + (2.0D0*z(b)-z(c)-z(a))*(r(a,b)**2*(x(b)-x(c))+r(b,c)**2*(x(b)-x(a)))
C
C Above is just a trick as otherwise get 'unbalanced parentheses' error message...
C
         vixen6=3.0D0*r(a,b)*r(b,c)*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))
         vixen7=u*(r(a,b)**2*(x(b)-x(c)) + r(b,c)**2*(x(b)-x(a)))
         vixen8= (vixen5*r(a,b)**3*r(b,c)**3 - vixen6*vixen7)/(r(a,b)**6*r(b,c)**6)
         vixen9=(2.0D0*x(b)-x(c)-x(a))*((z(b)-z(c))/r(b,c)**2 + (z(b)-z(a))/r(a,b)**2)/(r(a,b)*r(b,c))
         vixen4= vixen8 + vixen9

         anglehell(3*b-2,3*b)= anglehell(3*b-2,3*b) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*b-2)= anglehell(3*b,3*b-2) + vixen1*vixen4+vixen2*vixen3
C
C dybdzb,dzbdyb
C
         vixen2=u*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))/(r(a,b)**3*r(b,c)**3) - (2.0D0*y(b)-y(c)-y(a))/(r(a,b)*r(b,c))
         vixen3=2.0D0*ktheta*dthetaz2*(st-(theta(i)-theta0)*ct)/st**2

         vixen5=u*(2.0D0*(z(b)-z(a))*(y(b)-y(c))+2*(z(b)-z(c))*(y(b)-y(a)))
         vixen5= vixen5 + (2.0D0*z(b)-z(c)-z(a))*(r(a,b)**2*(y(b)-y(c))+r(b,c)**2*(y(b)-y(a)))
C
C Above is just a trick as otherwise get 'unbalanced parentheses' error message...
C
         vixen6=3.0D0*r(a,b)*r(b,c)*(r(a,b)**2*(z(b)-z(c)) + r(b,c)**2*(z(b)-z(a)))
         vixen7=u*(r(a,b)**2*(y(b)-y(c)) + r(b,c)**2*(y(b)-y(a)))
         vixen8= (vixen5*r(a,b)**3*r(b,c)**3 - vixen6*vixen7)/(r(a,b)**6*r(b,c)**6)
         vixen9=(2.0D0*y(b)-y(c)-y(a))*((z(b)-z(c))/r(b,c)**2 + (z(b)-z(a))/r(a,b)**2)/(r(a,b)*r(b,c))
         vixen4= vixen8 + vixen9

         anglehell(3*b-1,3*b)= anglehell(3*b-1,3*b) + vixen1*vixen4+vixen2*vixen3
         anglehell(3*b,3*b-1)= anglehell(3*b,3*b-1) + vixen1*vixen4+vixen2*vixen3

      END DO
C
C***************************
C TORSION ANGLE DERIVATIVES*
C***************************
C
      DO i=1,t
         a=da1(i)
         b=da2(i)
         c=da3(i)
         d=da4(i)
C         PRINT *,a,b,c,d,' : ',typech(a),'-',typech(b),'-',typech(c),'-',typech(d)

         rbc=r(b,c)
         rbc2=rbc**2
         rbc4=rbc**4
         rbc6=rbc**6
         xe=(x(c)-x(b))
         ye=(y(c)-y(b))
         ze=(z(c)-z(b))
         lambda=((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)/rbc2
         mu=    ((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)/rbc2
         xp= x(a)+lambda*x(c)-(1+lambda)*x(b)
         yp= y(a)+lambda*y(c)-(1+lambda)*y(b)
         zp= z(a)+lambda*z(c)-(1+lambda)*z(b)
         xpp= x(d)+mu*x(c)-(1+mu)*x(b)
         ypp= y(d)+mu*y(c)-(1+mu)*y(b)
         zpp= z(d)+mu*z(c)-(1+mu)*z(b)
         u= xp*xpp + yp*ypp + zp*zpp
         v= SQRT(xp**2 + yp**2 + zp**2)*SQRT(xpp**2 + ypp**2 + zpp**2)

         lambdacrap= (x(b)-x(a))*xe + (y(b)-y(a))*ye + (z(b)-z(a))*ze
         mucrap= (x(b)-x(d))*xe + (y(b)-y(d))*ye + (z(b)-z(d))*ze

         dlambdax1= (x(b)-x(c))/rbc2
         dlambday1= (y(b)-y(c))/rbc2
         dlambdaz1= (z(b)-z(c))/rbc2
         dlambdax2= (x(c)+x(a)-2.0D0*x(b))/rbc2 - 2.0D0*(x(b)-x(c))*lambdacrap/rbc4 
         dlambday2= (y(c)+y(a)-2.0D0*y(b))/rbc2 - 2.0D0*(y(b)-y(c))*lambdacrap/rbc4
         dlambdaz2= (z(c)+z(a)-2.0D0*z(b))/rbc2 - 2.0D0*(z(b)-z(c))*lambdacrap/rbc4
         dlambdax3= (x(b)-x(a))/rbc2 - 2.0D0*xe*lambdacrap/rbc4
         dlambday3= (y(b)-y(a))/rbc2 - 2.0D0*ye*lambdacrap/rbc4
         dlambdaz3= (z(b)-z(a))/rbc2 - 2.0D0*ze*lambdacrap/rbc4
         dmux2= (x(c)+x(d)-2.0D0*x(b))/rbc2 - 2.0D0*(x(b)-x(c))*mucrap/rbc4
         dmuy2= (y(c)+y(d)-2.0D0*y(b))/rbc2 - 2.0D0*(y(b)-y(c))*mucrap/rbc4
         dmuz2= (z(c)+z(d)-2.0D0*z(b))/rbc2 - 2.0D0*(z(b)-z(c))*mucrap/rbc4
         dmux3= (x(b)-x(d))/rbc2 - 2.0D0*xe*mucrap/rbc4
         dmuy3= (y(b)-y(d))/rbc2 - 2.0D0*ye*mucrap/rbc4
         dmuz3= (z(b)-z(d))/rbc2 - 2.0D0*ze*mucrap/rbc4
         dmux4= (x(b)-x(c))/rbc2
         dmuy4= (y(b)-y(c))/rbc2
         dmuz4= (z(b)-z(c))/rbc2
C
C Testing derivative maths! First the dx1 terms
C
         vixen1=SQRT(xp**2 + yp**2 + zp**2)
         vixen2=SQRT(xpp**2 + ypp**2 + zpp**2)
         dpx1=(1.0D0/vixen1)*(xp*(1.0D0+xe*dlambdax1) + yp*ye*dlambdax1 + zp*ze*dlambdax1)
         dvx1=vixen2*dpx1
         dux1=xpp*(1.0D0+xe*dlambdax1) + ypp*ye*dlambdax1 + zpp*ze*dlambdax1

         vixen5=(dvn(i)/did(i))*4.0D0*COS(dphi(i))*COS(ddelta(i))
         torscrap(3*a-2)= torscrap(3*a-2)+ vixen5*(v*dux1 - u*dvx1)/v**2

         dpy1=(1.0D0/vixen1)*(yp*(1.0D0+ye*dlambday1) + xp*xe*dlambday1 + zp*ze*dlambday1)
         dvy1=vixen2*dpy1
         duy1=ypp*(1.0D0+ye*dlambday1) + xpp*xe*dlambday1 + zpp*ze*dlambday1

         torscrap(3*a-1)= torscrap(3*a-1)+ vixen5*(v*duy1 - u*dvy1)/v**2

         dpz1=(1.0D0/vixen1)*(zp*(1.0D0+ze*dlambdaz1) + xp*xe*dlambdaz1 + yp*ye*dlambdaz1)
         dvz1=vixen2*dpz1
         duz1=zpp*(1.0D0+ze*dlambdaz1) + xpp*xe*dlambdaz1 + ypp*ye*dlambdaz1

         torscrap(3*a)= torscrap(3*a)+ vixen5*(v*duz1 - u*dvz1)/v**2
C
C Now dx2 terms
C
         dux2=xp*(xe*dmux2-1.0D0-mu) + xpp*(xe*dlambdax2-1.0D0-lambda) + ye*(yp*dmux2 + ypp*dlambdax2) +
     1          ze*(zp*dmux2 + zpp*dlambdax2)
         dpx2=(1.0D0/vixen1)*(xp*(xe*dlambdax2-1.0D0-lambda) + yp*ye*dlambdax2 + zp*ze*dlambdax2)
         vixen6=vixen2*dpx2
         dqx2=(1.0D0/vixen2)*(xpp*(xe*dmux2-1.0D0-mu) + ypp*ye*dmux2 + zpp*ze*dmux2)
         vixen7=vixen1*dqx2
         dvx2=vixen6+vixen7

         torscrap(3*b-2)= torscrap(3*b-2)+ vixen5*(v*dux2 - u*dvx2)/v**2

         duy2=yp*(ye*dmuy2-1.0D0-mu) + ypp*(ye*dlambday2-1.0D0-lambda) + xe*(xp*dmuy2 + xpp*dlambday2) +
     1          ze*(zp*dmuy2 + zpp*dlambday2)
         dpy2=(1.0D0/vixen1)*(yp*(ye*dlambday2-1.0D0-lambda) + xp*xe*dlambday2 + zp*ze*dlambday2)
         vixen6=vixen2*dpy2
         dqy2=(1.0D0/vixen2)*(ypp*(ye*dmuy2-1.0D0-mu) + xpp*xe*dmuy2 + zpp*ze*dmuy2)
         vixen7=vixen1*dqy2
         dvy2=vixen6+vixen7

         torscrap(3*b-1)= torscrap(3*b-1)+ vixen5*(v*duy2 - u*dvy2)/v**2

         duz2=zp*(ze*dmuz2-1.0D0-mu) + zpp*(ze*dlambdaz2-1.0D0-lambda) + xe*(xp*dmuz2 + xpp*dlambdaz2) +
     1          ye*(yp*dmuz2 + ypp*dlambdaz2)
         dpz2=(1.0D0/vixen1)*(zp*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*dlambdaz2 + yp*ye*dlambdaz2)
         vixen6=vixen2*dpz2
         dqz2=(1.0D0/vixen2)*(zpp*(ze*dmuz2-1.0D0-mu) + xpp*xe*dmuz2 + ypp*ye*dmuz2)
         vixen7=vixen1*dqz2
         dvz2=vixen6+vixen7

         torscrap(3*b)= torscrap(3*b)+ vixen5*(v*duz2 - u*dvz2)/v**2
C
C Now dx3 terms
C
         dux3=xp*(mu + xe*dmux3) + xpp*(lambda + xe*dlambdax3) + ye*(yp*dmux3 + ypp*dlambdax3) +
     1          ze*(zp*dmux3 + zpp*dlambdax3)
         dpx3=(1.0D0/vixen1)*(xp*(lambda + xe*dlambdax3) + yp*ye*dlambdax3 + zp*ze*dlambdax3)
         vixen6=vixen2*dpx3
         dqx3=(1.0D0/vixen2)*(xpp*(mu + xe*dmux3) + ypp*ye*dmux3 + zpp*ze*dmux3)
         vixen7=vixen1*dqx3
         dvx3=vixen6+vixen7

         torscrap(3*c-2)= torscrap(3*c-2)+ vixen5*(v*dux3 - u*dvx3)/v**2

         duy3=yp*(mu + ye*dmuy3) + ypp*(lambda + ye*dlambday3) + xe*(xp*dmuy3 + xpp*dlambday3) +
     1          ze*(zp*dmuy3 + zpp*dlambday3)
         dpy3=(1.0D0/vixen1)*(yp*(lambda + ye*dlambday3) + xp*xe*dlambday3 + zp*ze*dlambday3)
         vixen6=vixen2*dpy3
         dqy3=(1.0D0/vixen2)*(ypp*(mu + ye*dmuy3) + xpp*xe*dmuy3 + zpp*ze*dmuy3)
         vixen7=vixen1*dqy3
         dvy3=vixen6+vixen7

         torscrap(3*c-1)= torscrap(3*c-1)+ vixen5*(v*duy3 - u*dvy3)/v**2

         duz3=zp*(mu + ze*dmuz3) + zpp*(lambda + ze*dlambdaz3) + xe*(xp*dmuz3 + xpp*dlambdaz3) +
     1          ye*(yp*dmuz3 + ypp*dlambdaz3)
         dpz3=(1.0D0/vixen1)*(zp*(lambda + ze*dlambdaz3) + xp*xe*dlambdaz3 + yp*ye*dlambdaz3)
         vixen6=vixen2*dpz3
         dqz3=(1.0D0/vixen2)*(zpp*(mu + ze*dmuz3) + xpp*xe*dmuz3 + ypp*ye*dmuz3)
         vixen7=vixen1*dqz3
         dvz3=vixen6+vixen7

         torscrap(3*c)= torscrap(3*c)+ vixen5*(v*duz3 - u*dvz3)/v**2
C
C Now dx4 terms
C
         dux4=xp*(1.0D0 + xe*dmux4) + yp*ye*dmux4 + zp*ze*dmux4
         dqx4=(1.0D0/vixen2)*(xpp*(1.0D0 + xe*dmux4) + ypp*ye*dmux4 + zpp*ze*dmux4)
         dvx4=vixen1*dqx4
  
         torscrap(3*d-2)= torscrap(3*d-2)+ vixen5*(v*dux4 - u*dvx4)/v**2

         duy4=yp*(1.0D0 + ye*dmuy4) + xp*xe*dmuy4 + zp*ze*dmuy4
         dqy4=(1.0D0/vixen2)*(ypp*(1.0D0 + ye*dmuy4) + xpp*xe*dmuy4 + zpp*ze*dmuy4)
         dvy4=vixen1*dqy4

         torscrap(3*d-1)= torscrap(3*d-1)+ vixen5*(v*duy4 - u*dvy4)/v**2

         duz4=zp*(1.0D0 + ze*dmuz4) + xp*xe*dmuz4 + yp*ye*dmuz4
         dqz4=(1.0D0/vixen2)*(zpp*(1.0D0 + ze*dmuz4) + xpp*xe*dmuz4 + ypp*ye*dmuz4)
         dvz4=vixen1*dqz4

         torscrap(3*d)= torscrap(3*d)+ vixen5*(v*duz4 - u*dvz4)/v**2
C
C Now second derivatives
C
         CALL hairyhell(u,v)
C
C dx1dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= ((1.0D0+xe*dlambdax1)**2 + ye**2*dlambdax1**2 + ze**2*dlambdax1**2 - dpx1**2)/vixen1
         dnx1= -(u*vixen2/v**2)*vixen3 -2.0D0*dvx1*(v*dux1 - u*dvx1)/v**3

         torshell(3*a-2,3*a-2) =torshell(3*a-2,3*a-2) + m*dnx1 + n*dmx1
C
C dy1dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= ((1.0D0+ye*dlambday1)**2 + xe**2*dlambday1**2 + ze**2*dlambday1**2 - dpy1**2)/vixen1
         dny1= -(u*vixen2/v**2)*vixen3 -2.0D0*dvy1*(v*duy1 - u*dvy1)/v**3

         torshell(3*a-1,3*a-1) =torshell(3*a-1,3*a-1) + m*dny1 + n*dmy1
C
C dz1dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= ((1.0D0+ze*dlambdaz1)**2 + xe**2*dlambdaz1**2 + ye**2*dlambdaz1**2 - dpz1**2)/vixen1
         dnz1= -(u*vixen2/v**2)*vixen3 -2.0D0*dvz1*(v*duz1 - u*dvz1)/v**3

         torshell(3*a,3*a) =torshell(3*a,3*a) + m*dnz1 + n*dmz1
C
C dx1dy1,dy1dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= ((1.0D0+xe*dlambdax1)*xe*dlambday1 + ye*dlambdax1*(1.0D0+ye*dlambday1) 
     1          +ze**2*dlambdax1*dlambday1 - dpx1*dpy1)/vixen1
         dny1= (dux1*dvy1 - duy1*dvx1 - u*vixen2*vixen3)/v**2 - 2.0D0*dvy1*(v*dux1 - u*dvx1)/v**3
         
         torshell(3*a-2,3*a-1)= torshell(3*a-2,3*a-1) + m*dny1 + n*dmy1
         torshell(3*a-1,3*a-2)= torshell(3*a-1,3*a-2) + m*dny1 + n*dmy1
C
C dx1dz1,dz1dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= ((1.0D0+xe*dlambdax1)*xe*dlambdaz1 + ze*dlambdax1*(1.0D0+ze*dlambdaz1)
     1          +ye**2*dlambdax1*dlambdaz1 - dpx1*dpz1)/vixen1
         dnz1= (dux1*dvz1 - duz1*dvx1 - u*vixen2*vixen3)/v**2 - 2.0D0*dvz1*(v*dux1 - u*dvx1)/v**3

         torshell(3*a-2,3*a)= torshell(3*a-2,3*a) + m*dnz1 + n*dmz1
         torshell(3*a,3*a-2)= torshell(3*a,3*a-2) + m*dnz1 + n*dmz1
C
C dy1dz1,dz1dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= ((1.0D0+ye*dlambday1)*ye*dlambdaz1 + ze*dlambday1*(1.0D0+ze*dlambdaz1)
     1          +xe**2*dlambday1*dlambdaz1 - dpy1*dpz1)/vixen1
         dnz1= (duy1*dvz1 - duz1*dvy1 - u*vixen2*vixen3)/v**2 - 2.0D0*dvz1*(v*duy1 - u*dvy1)/v**3

         torshell(3*a-1,3*a)= torshell(3*a-1,3*a) + m*dnz1 + n*dmz1
         torshell(3*a,3*a-1)= torshell(3*a,3*a-1) + m*dnz1 + n*dmz1
C
C dx1dx4,dx4dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= (1.0D0+xe*dlambdax1)*(1.0D0+xe*dmux4) + (ye**2 + ze**2)*dlambdax1*dmux4
         dnx4= (v*vixen3 + dvx4*dux1 - u*dqx4*dpx1 - dvx1*dux4)/v**2 - 2.0D0*dvx4*n/v

         torshell(3*a-2,3*d-2)= torshell(3*a-2,3*d-2) + m*dnx4 + n*dmx4
         torshell(3*d-2,3*a-2)= torshell(3*d-2,3*a-2) + m*dnx4 + n*dmx4
C
C dy1dy4,dy4dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= (1.0D0+ye*dlambday1)*(1.0D0+ye*dmuy4) + (xe**2 + ze**2)*dlambday1*dmuy4
         dny4= (v*vixen3 + dvy4*duy1 - u*dqy4*dpy1 - dvy1*duy4)/v**2 - 2.0D0*dvy4*n/v

         torshell(3*a-1,3*d-1)= torshell(3*a-1,3*d-1) + m*dny4 + n*dmy4
         torshell(3*d-1,3*a-1)= torshell(3*d-1,3*a-1) + m*dny4 + n*dmy4
C
C dz1dz4,dz4dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= (1.0D0+ze*dlambdaz1)*(1.0D0+ze*dmuz4) + (ye**2 + xe**2)*dlambdaz1*dmuz4
         dnz4= (v*vixen3 + dvz4*duz1 - u*dqz4*dpz1 - dvz1*duz4)/v**2 - 2.0D0*dvz4*n/v

         torshell(3*a,3*d)= torshell(3*a,3*d) + m*dnz4 + n*dmz4
         torshell(3*d,3*a)= torshell(3*d,3*a) + m*dnz4 + n*dmz4
C
C dx1dy4,dy4dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= (1.0D0+xe*dlambdax1)*xe*dmuy4 + ye*dlambdax1*(1.0D0+ye*dmuy4) 
     1          + ze**2*dlambdax1*dmuy4
         dny4= (v*vixen3 + dvy4*dux1 - u*dqy4*dpx1 - dvx1*duy4)/v**2 - 2.0D0*dvy4*n/v

         torshell(3*a-2,3*d-1)= torshell(3*a-2,3*d-1) + m*dny4 + n*dmy4
         torshell(3*d-1,3*a-2)= torshell(3*d-1,3*a-2) + m*dny4 + n*dmy4
C
C dx1dz4,dz4dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= (1.0D0+xe*dlambdax1)*xe*dmuz4 + ze*dlambdax1*(1.0D0+ze*dmuz4)
     1          + ye**2*dlambdax1*dmuz4
         dnz4= (v*vixen3 + dvz4*dux1 - u*dqz4*dpx1 - dvx1*duz4)/v**2 - 2.0D0*dvz4*n/v

         torshell(3*a-2,3*d)= torshell(3*a-2,3*d) + m*dnz4 + n*dmz4
         torshell(3*d,3*a-2)= torshell(3*d,3*a-2) + m*dnz4 + n*dmz4
C
C dy1dz4,dz4dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= (1.0D0+ye*dlambday1)*ye*dmuz4 + ze*dlambday1*(1.0D0+ze*dmuz4)
     1          + xe**2*dlambday1*dmuz4
         dnz4= (v*vixen3 + dvz4*duy1 - u*dqz4*dpy1 - dvy1*duz4)/v**2 - 2.0D0*dvz4*n/v

         torshell(3*a-1,3*d)= torshell(3*a-1,3*d) + m*dnz4 + n*dmz4
         torshell(3*d,3*a-1)= torshell(3*d,3*a-1) + m*dnz4 + n*dmz4
C
C dy1dx4,dx4dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= (1.0D0+ye*dlambday1)*ye*dmux4 + xe*dlambday1*(1.0D0+xe*dmux4)
     1          + ze**2*dlambday1*dmux4
         dnx4= (v*vixen3 + dvx4*duy1 - u*dqx4*dpy1 - dvy1*dux4)/v**2 - 2.0D0*dvx4*n/v

         torshell(3*a-1,3*d-2)= torshell(3*a-1,3*d-2) + m*dnx4 + n*dmx4
         torshell(3*d-2,3*a-1)= torshell(3*d-2,3*a-1) + m*dnx4 + n*dmx4
C
C dz1dx4,dx4dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= (1.0D0+ze*dlambdaz1)*ze*dmux4 + xe*dlambdaz1*(1.0D0+xe*dmux4)
     1          + ye**2*dlambdaz1*dmux4
         dnx4= (v*vixen3 + dvx4*duz1 - u*dqx4*dpz1 - dvz1*dux4)/v**2 - 2.0D0*dvx4*n/v

         torshell(3*a,3*d-2)= torshell(3*a,3*d-2) + m*dnx4 + n*dmx4
         torshell(3*d-2,3*a)= torshell(3*d-2,3*a) + m*dnx4 + n*dmx4
C
C dz1dy4,dy4dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= (1.0D0+ze*dlambdaz1)*ze*dmuy4 + ye*dlambdaz1*(1.0D0+ye*dmuy4)
     1          + xe**2*dlambdaz1*dmuy4
         dny4= (v*vixen3 + dvy4*duz1 - u*dqy4*dpz1 - dvz1*duy4)/v**2 - 2.0D0*dvy4*n/v

         torshell(3*a,3*d-1)= torshell(3*a,3*d-1) + m*dny4 + n*dmy4
         torshell(3*d-1,3*a)= torshell(3*d-1,3*a) + m*dny4 + n*dmy4
C
C dx4dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+xe*dmux4)**2 + (ye**2 + ze**2)*dmux4**2 - dqx4**2)
         dnx4= -u*vixen3/v**2 - 2.0D0*dvx4*n/v

         torshell(3*d-2,3*d-2)= torshell(3*d-2,3*d-2) + m*dnx4 + n*dmx4
C
C dy4dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+ye*dmuy4)**2 + (xe**2 + ze**2)*dmuy4**2 - dqy4**2)
         dny4= -u*vixen3/v**2 - 2.0D0*dvy4*n/v

         torshell(3*d-1,3*d-1)= torshell(3*d-1,3*d-1) + m*dny4 + n*dmy4
C
C dz4dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+ze*dmuz4)**2 + (ye**2 + xe**2)*dmuz4**2 - dqz4**2)
         dnz4= -u*vixen3/v**2 - 2.0D0*dvz4*n/v

         torshell(3*d,3*d)= torshell(3*d,3*d) + m*dnz4 + n*dmz4
C
C dx4dy4,dy4dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+xe*dmux4)*xe*dmuy4 + ye*dmux4*(1.0D0+ye*dmuy4)
     1          + ze**2*dmux4*dmuy4 - dqx4*dqy4)
         dny4= (dvy4*dux4 - u*vixen3 - duy4*dvx4)/v**2 - 2.0D0*dvy4*n/v

         torshell(3*d-2,3*d-1)= torshell(3*d-2,3*d-1) + m*dny4 + n*dmy4
         torshell(3*d-1,3*d-2)= torshell(3*d-1,3*d-2) + m*dny4 + n*dmy4
C
C dx4dz4,dz4dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+xe*dmux4)*xe*dmuz4 + ze*dmux4*(1.0D0+ze*dmuz4)
     1          + ye**2*dmux4*dmuz4 - dqx4*dqz4)
         dnz4= (dvz4*dux4 - u*vixen3 - duz4*dvx4)/v**2 - 2.0D0*dvz4*n/v

         torshell(3*d-2,3*d)= torshell(3*d-2,3*d) + m*dnz4 + n*dmz4
         torshell(3*d,3*d-2)= torshell(3*d,3*d-2) + m*dnz4 + n*dmz4
C
C dy4dz4,dz4dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+ye*dmuy4)*ye*dmuz4 + ze*dmuy4*(1.0D0+ze*dmuz4)
     1          + xe**2*dmuy4*dmuz4 - dqy4*dqz4)
         dnz4= (dvz4*duy4 - u*vixen3 - duz4*dvy4)/v**2 - 2.0D0*dvz4*n/v

         torshell(3*d-1,3*d)= torshell(3*d-1,3*d) + m*dnz4 + n*dmz4
         torshell(3*d,3*d-1)= torshell(3*d,3*d-1) + m*dnz4 + n*dmz4
C
C dx1dx2,dx2dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(x(b)-x(c))**2/rbc4
         vixen5= xp*(vixen4*xe-dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dlambdax2-1.0D0-lambda)
     1           + yp*ye*vixen4 + ye**2*dlambdax1*dlambdax2 + zp*ze*vixen4 
     2           + ze**2*dlambdax1*dlambdax2 - dpx1*dpx2
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqx2 
         vixen6= xpp*(vixen4*xe-dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dmux2-1.0D0-mu) +ypp*ye*vixen4
     1           + ye**2*dlambdax1*dmux2 + zpp*ze*vixen4 + ze**2*dlambdax1*dmux2
         dnx2= (v*vixen6 + dvx2*dux1 - u*vixen3 - dux2*dvx1)/v**2 - 2.0D0*dvx2*n/v

         torshell(3*a-2,3*b-2)= torshell(3*a-2,3*b-2)+ m*dnx2 + n*dmx2
         torshell(3*b-2,3*a-2)= torshell(3*b-2,3*a-2)+ m*dnx2 + n*dmx2
C
C dy1dy2,dy2dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(y(b)-y(c))**2/rbc4
         vixen5= yp*(vixen4*ye-dlambday1) + (1.0D0+ye*dlambday1)*(ye*dlambday2-1.0D0-lambda)
     1           + xp*xe*vixen4 + xe**2*dlambday1*dlambday2 + zp*ze*vixen4
     2           + ze**2*dlambday1*dlambday2 - dpy1*dpy2
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqy2
         vixen6= ypp*(vixen4*ye-dlambday1) + (1.0D0+ye*dlambday1)*(ye*dmuy2-1.0D0-mu) +xpp*xe*vixen4
     1           + xe**2*dlambday1*dmuy2 + zpp*ze*vixen4 + ze**2*dlambday1*dmuy2
         dny2= (v*vixen6 + dvy2*duy1 - u*vixen3 - duy2*dvy1)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*a-1,3*b-1)= torshell(3*a-1,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*a-1)= torshell(3*b-1,3*a-1)+ m*dny2 + n*dmy2
C
C dz1dz2,dz2dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(z(b)-z(c))**2/rbc4
         vixen5= zp*(vixen4*ze-dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dlambdaz2-1.0D0-lambda)
     1           + xp*xe*vixen4 + xe**2*dlambdaz1*dlambdaz2 + yp*ye*vixen4
     2           + ye**2*dlambdaz1*dlambdaz2 - dpz1*dpz2
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqz2
         vixen6= zpp*(vixen4*ze-dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dmuz2-1.0D0-mu) +xpp*xe*vixen4
     1           + xe**2*dlambdaz1*dmuz2 + ypp*ye*vixen4 + ye**2*dlambdaz1*dmuz2
         dnz2= (v*vixen6 + dvz2*duz1 - u*vixen3 - duz2*dvz1)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*a,3*b)= torshell(3*a,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*a)= torshell(3*b,3*a)+ m*dnz2 + n*dmz2
C
C dx1dy2,dy2dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambday2 + yp*(ye*vixen4-dlambdax1)
     1          + ye*dlambdax1*(ye*dlambday2-1.0D0-lambda) + zp*ze*vixen4 
     2          + ze**2*dlambdax1*dlambday2 - dpx1*dpy2
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqy2
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuy2 + ypp*(ye*vixen4-dlambdax1)
     1          + ye*dlambdax1*(ye*dmuy2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dlambdax1*dmuy2
         dny2= (v*vixen6 + dvy2*dux1 - u*vixen3 - duy2*dvx1)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*a-2,3*b-1)= torshell(3*a-2,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*a-2)= torshell(3*b-1,3*a-2)+ m*dny2 + n*dmy2
C
C dx1dz2,dz2dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambdaz2 + zp*(ze*vixen4-dlambdax1)
     1          + ze*dlambdax1*(ze*dlambdaz2-1.0D0-lambda) + yp*ye*vixen4
     2          + ye**2*dlambdax1*dlambdaz2 - dpx1*dpz2
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqz2
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuz2 + zpp*(ze*vixen4-dlambdax1)
     1          + ze*dlambdax1*(ze*dmuz2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dlambdax1*dmuz2
         dnz2= (v*vixen6 + dvz2*dux1 - u*vixen3 - duz2*dvx1)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*a-2,3*b)= torshell(3*a-2,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*a-2)= torshell(3*b,3*a-2)+ m*dnz2 + n*dmz2
C
C dy1dz2,dz2dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdaz2 + zp*(ze*vixen4-dlambday1)
     1          + ze*dlambday1*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dlambday1*dlambdaz2 - dpy1*dpz2
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqz2
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmuz2 + zpp*(ze*vixen4-dlambday1)
     1          + ze*dlambday1*(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dlambday1*dmuz2
         dnz2= (v*vixen6 + dvz2*duy1 - u*vixen3 - duz2*dvy1)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*a-1,3*b)= torshell(3*a-1,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*a-1)= torshell(3*b,3*a-1)+ m*dnz2 + n*dmz2
C
C dy1dx2,dx2dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdax2 + xp*(xe*vixen4-dlambday1)
     1          + xe*dlambday1*(xe*dlambdax2-1.0D0-lambda) + zp*ze*vixen4
     2          + ze**2*dlambday1*dlambdax2 - dpy1*dpx2
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqx2
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmux2 + xpp*(xe*vixen4-dlambday1)
     1          + xe*dlambday1*(xe*dmux2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dlambday1*dmux2
         dnx2= (v*vixen6 + dvx2*duy1 - u*vixen3 - dux2*dvy1)/v**2 - 2.0D0*dvx2*n/v

         torshell(3*a-1,3*b-2)= torshell(3*a-1,3*b-2)+ m*dnx2 + n*dmx2
         torshell(3*b-2,3*a-1)= torshell(3*b-2,3*a-1)+ m*dnx2 + n*dmx2
C
C dz1dx2,dx2dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambdax2 + xp*(xe*vixen4-dlambdaz1)
     1          + xe*dlambdaz1*(xe*dlambdax2-1.0D0-lambda) + yp*ye*vixen4
     2          + ye**2*dlambdaz1*dlambdax2 - dpz1*dpx2
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqx2
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmux2 + xpp*(xe*vixen4-dlambdaz1)
     1          + xe*dlambdaz1*(xe*dmux2-1-mu) + ypp*ye*vixen4 + ye**2*dlambdaz1*dmux2
         dnx2= (v*vixen6 + dvx2*duz1 - u*vixen3 - dux2*dvz1)/v**2 - 2.0D0*dvx2*n/v

         torshell(3*a,3*b-2)= torshell(3*a,3*b-2)+ m*dnx2 + n*dmx2
         torshell(3*b-2,3*a)= torshell(3*b-2,3*a)+ m*dnx2 + n*dmx2
C
C dz1dy2,dy2dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambday2 + yp*(ye*vixen4-dlambdaz1)
     1          + ye*dlambdaz1*(ye*dlambday2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dlambdaz1*dlambday2 - dpz1*dpy2
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqy2
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmuy2 + ypp*(ye*vixen4-dlambdaz1)
     1          + ye*dlambdaz1*(ye*dmuy2-1-mu) + xpp*xe*vixen4 + xe**2*dlambdaz1*dmuy2
         dny2= (v*vixen6 + dvy2*duz1 - u*vixen3 - duy2*dvz1)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*a,3*b-1)= torshell(3*a,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*a)= torshell(3*b-1,3*a)+ m*dny2 + n*dmy2
C
C dx1dx3,dx3dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 2.0D0*(x(b)-x(c))**2/rbc4 - 1/rbc2
         vixen5= xp*(xe*vixen4+dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dlambdax3+lambda)
     1          + (yp*ye+zp*ze)*vixen4 + (ye**2+ze**2)*dlambdax1*dlambdax3 - dpx1*dpx3
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqx3
         vixen6= xpp*(xe*vixen4+dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dmux3+mu)
     1          + (ypp*ye+zpp*ze)*vixen4 + (ye**2+ze**2)*dlambdax1*dmux3
         dnx3= (v*vixen6 + dvx3*dux1 - u*vixen3 - dux3*dvx1)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*a-2,3*c-2)= torshell(3*a-2,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*a-2)= torshell(3*c-2,3*a-2)+ m*dnx3 + n*dmx3
C
C dy1dy3,dy3dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 2.0D0*(y(b)-y(c))**2/rbc4 - 1/rbc2
         vixen5= yp*(ye*vixen4+dlambday1) + (1.0D0+ye*dlambday1)*(ye*dlambday3+lambda)
     1          + (xp*xe+zp*ze)*vixen4 + (xe**2+ze**2)*dlambday1*dlambday3 - dpy1*dpy3
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqy3
         vixen6= ypp*(ye*vixen4+dlambday1) + (1.0D0+ye*dlambday1)*(ye*dmuy3+mu)
     1          + (xpp*xe+zpp*ze)*vixen4 + (xe**2+ze**2)*dlambday1*dmuy3
         dny3= (v*vixen6 + dvy3*duy1 - u*vixen3 - duy3*dvy1)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*a-1,3*c-1)= torshell(3*a-1,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*a-1)= torshell(3*c-1,3*a-1)+ m*dny3 + n*dmy3
C
C dz1dz3,dz3dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 2.0D0*(z(b)-z(c))**2/rbc4 - 1/rbc2
         vixen5= zp*(ze*vixen4+dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dlambdaz3+lambda)
     1          + (xp*xe+yp*ye)*vixen4 + (xe**2+ye**2)*dlambdaz1*dlambdaz3 - dpz1*dpz3
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqz3
         vixen6= zpp*(ze*vixen4+dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dmuz3+mu)
     1          + (xpp*xe+ypp*ye)*vixen4 + (xe**2+ye**2)*dlambdaz1*dmuz3
         dnz3= (v*vixen6 + dvz3*duz1 - u*vixen3 - duz3*dvz1)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*a,3*c)= torshell(3*a,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*a)= torshell(3*c,3*a)+ m*dnz3 + n*dmz3
C
C dx1dy3
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambday3 
     1          + yp*(ye*vixen4+dlambdax1) + ye*dlambdax1*(ye*dlambday3+lambda) 
     2          + zp*ze*vixen4 + ze**2*dlambdax1*dlambday3 - dpx1*dpy3
         vixen3=(vixen2/vixen1)*vixen5 + dpx1*dqy3
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuy3 + ypp*(ye*vixen4+dlambdax1)
     1          + ye*dlambdax1*(ye*dmuy3+mu) + zpp*ze*vixen4 + ze**2*dlambdax1*dmuy3
         dny3= (v*vixen6 + dvy3*dux1 - u*vixen3 - duy3*dvx1)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*a-2,3*c-1)= torshell(3*a-2,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*a-2)= torshell(3*c-1,3*a-2)+ m*dny3 + n*dmy3
C
C dx1dz3
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambdaz3
     1          + zp*(ze*vixen4+dlambdax1) + ze*dlambdax1*(ze*dlambdaz3+lambda)
     2          + yp*ye*vixen4 + ye**2*dlambdax1*dlambdaz3 - dpx1*dpz3
         vixen3=(vixen2/vixen1)*vixen5 + dpx1*dqz3
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuz3 + zpp*(ze*vixen4+dlambdax1)
     1          + ze*dlambdax1*(ze*dmuz3+mu) + ypp*ye*vixen4 + ye**2*dlambdax1*dmuz3
         dnz3= (v*vixen6 + dvz3*dux1 - u*vixen3 - duz3*dvx1)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*a-2,3*c)= torshell(3*a-2,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*a-2)= torshell(3*c,3*a-2)+ m*dnz3 + n*dmz3
C
C dy1dz3
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdaz3
     1          + zp*(ze*vixen4+dlambday1) + ze*dlambday1*(ze*dlambdaz3+lambda)
     2          + xp*xe*vixen4 + xe**2*dlambday1*dlambdaz3 - dpy1*dpz3
         vixen3=(vixen2/vixen1)*vixen5 + dpy1*dqz3
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmuz3 + zpp*(ze*vixen4+dlambday1)
     1          + ze*dlambday1*(ze*dmuz3+mu) + xpp*xe*vixen4 + xe**2*dlambday1*dmuz3
         dnz3= (v*vixen6 + dvz3*duy1 - u*vixen3 - duz3*dvy1)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*a-1,3*c)= torshell(3*a-1,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*a-1)= torshell(3*c,3*a-1)+ m*dnz3 + n*dmz3
C
C dy1dx3
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdax3
     1          + xp*(xe*vixen4+dlambday1) + xe*dlambday1*(xe*dlambdax3+lambda)
     2          + zp*ze*vixen4 + ze**2*dlambday1*dlambdax3 - dpy1*dpx3
         vixen3=(vixen2/vixen1)*vixen5 + dpy1*dqx3
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmux3 + xpp*(xe*vixen4+dlambday1)
     1          + xe*dlambday1*(xe*dmux3+mu) + zpp*ze*vixen4 + ze**2*dlambday1*dmux3
         dnx3= (v*vixen6 + dvx3*duy1 - u*vixen3 - dux3*dvy1)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*a-1,3*c-2)= torshell(3*a-1,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*a-1)= torshell(3*c-2,3*a-1)+ m*dnx3 + n*dmx3
C
C dz1dx3
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambdax3
     1          + xp*(xe*vixen4+dlambdaz1) + xe*dlambdaz1*(xe*dlambdax3+lambda)
     2          + yp*ye*vixen4 + ye**2*dlambdaz1*dlambdax3 - dpz1*dpx3
         vixen3=(vixen2/vixen1)*vixen5 + dpz1*dqx3
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmux3 + xpp*(xe*vixen4+dlambdaz1)
     1          + xe*dlambdaz1*(xe*dmux3+mu) + ypp*ye*vixen4 + ye**2*dlambdaz1*dmux3
         dnx3= (v*vixen6 + dvx3*duz1 - u*vixen3 - dux3*dvz1)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*a,3*c-2)= torshell(3*a,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*a)= torshell(3*c-2,3*a)+ m*dnx3 + n*dmx3
C
C dz1dy3
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambday3
     1          + yp*(ye*vixen4+dlambdaz1) + ye*dlambdaz1*(ye*dlambday3+lambda)
     2          + xp*xe*vixen4 + xe**2*dlambdaz1*dlambday3 - dpz1*dpy3
         vixen3=(vixen2/vixen1)*vixen5 + dpz1*dqy3
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmuy3 + ypp*(ye*vixen4+dlambdaz1)
     1          + ye*dlambdaz1*(ye*dmuy3+mu) + xpp*xe*vixen4 + xe**2*dlambdaz1*dmuy3
         dny3= (v*vixen6 + dvy3*duz1 - u*vixen3 - duy3*dvz1)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*a,3*c-1)= torshell(3*a,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*a)= torshell(3*c-1,3*a)+ m*dny3 + n*dmy3
C
C dx4dx2,dx2dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(x(b)-x(c))**2/rbc4
         vixen5= xpp*(xe*vixen4-dmux4) + (1.0D0+xe*dmux4)
     1          *(xe*dmux2-1.0D0-mu) + ypp*ye*vixen4
     2          + ye**2*dmux4*dmux2 + zpp*ze*vixen4
     3          + ze**2*dmux4*dmux2 - dqx4*dqx2
         vixen3=(vixen1/vixen2)*vixen5 + dpx2*dqx4
         vixen6= xp*(xe*vixen4-dmux4) + (1.0D0+xe*dmux4)
     1          *(xe*dlambdax2-1.0D0-lambda) + yp*ye*vixen4
     2          + ye**2*dmux4*dlambdax2 + zp*ze*vixen4
     3          + ze**2*dmux4*dlambdax2
         dnx2= (v*vixen6 + dvx2*dux4 - u*vixen3 - dux2*dvx4)/v**2
     1        - 2.0D0*dvx2*(v*dux4 - u*dvx4)/v**3

         torshell(3*d-2,3*b-2)= torshell(3*d-2,3*b-2)+ m*dnx2 + n*dmx2
         torshell(3*b-2,3*d-2)= torshell(3*b-2,3*d-2)+ m*dnx2 + n*dmx2
C
C dy4dy2,dy2dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(y(b)-y(c))**2/rbc4
         vixen5= ypp*(ye*vixen4-dmuy4) + (1.0D0+ye*dmuy4)
     1          *(ye*dmuy2-1.0D0-mu) + xpp*xe*vixen4
     2          + xe**2*dmuy4*dmuy2 + zpp*ze*vixen4
     3          + ze**2*dmuy4*dmuy2 - dqy4*dqy2
         vixen3=(vixen1/vixen2)*vixen5 + dpy2*dqy4
         vixen6= yp*(ye*vixen4-dmuy4) + (1.0D0+ye*dmuy4)
     1          *(ye*dlambday2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dmuy4*dlambday2 + zp*ze*vixen4
     3          + ze**2*dmuy4*dlambday2
         dny2= (v*vixen6 + dvy2*duy4 - u*vixen3 - duy2*dvy4)/v**2
     1        - 2.0D0*dvy2*(v*duy4 - u*dvy4)/v**3

         torshell(3*d-1,3*b-1)= torshell(3*d-1,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*d-1)= torshell(3*b-1,3*d-1)+ m*dny2 + n*dmy2
C
C dz4dz2,dz2dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(z(b)-z(c))**2/rbc4
         vixen5= zpp*(ze*vixen4-dmuz4) + (1.0D0+ze*dmuz4)
     1          *(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4
     2          + xe**2*dmuz4*dmuz2 + ypp*ye*vixen4
     3          + ye**2*dmuz4*dmuz2 - dqz4*dqz2
         vixen3=(vixen1/vixen2)*vixen5 + dpz2*dqz4
         vixen6= zp*(ze*vixen4-dmuz4) + (1.0D0+ze*dmuz4)
     1          *(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dmuz4*dlambdaz2 + yp*ye*vixen4
     3          + ye**2*dmuz4*dlambdaz2
         dnz2= (v*vixen6 + dvz2*duz4 - u*vixen3 - duz2*dvz4)/v**2
     1        - 2.0D0*dvz2*(v*duz4 - u*dvz4)/v**3

         torshell(3*d,3*b)= torshell(3*d,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*d)= torshell(3*b,3*d)+ m*dnz2 + n*dmz2
C
C dx4dy2,dy2dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuy2 + ypp*(ye*vixen4-dmux4)
     1          + ye*dmux4*(ye*dmuy2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dmux4*dmuy2 - dqx4*dqy2
         vixen3=(vixen1/vixen2)*vixen5 + dpy2*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambday2 + yp*(ye*vixen4-dmux4)
     1          + ye*dmux4*(ye*dlambday2-1.0D0-lambda) + zp*ze*vixen4 + ze**2*dmux4*dlambday2
         dny2= (v*vixen6 + dvy2*dux4 - u*vixen3 - duy2*dvx4)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*d-2,3*b-1)= torshell(3*d-2,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*d-2)= torshell(3*b-1,3*d-2)+ m*dny2 + n*dmy2
C
C dx4dz2,dz2dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuz2 + zpp*(ze*vixen4-dmux4)
     1          + ze*dmux4*(ze*dmuz2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dmux4*dmuz2 - dqx4*dqz2
         vixen3=(vixen1/vixen2)*vixen5 + dpz2*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambdaz2 + zp*(ze*vixen4-dmux4)
     1          + ze*dmux4*(ze*dlambdaz2-1.0D0-lambda) + yp*ye*vixen4 + ye**2*dmux4*dlambdaz2
         dnz2= (v*vixen6 + dvz2*dux4 - u*vixen3 - duz2*dvx4)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*d-2,3*b)= torshell(3*d-2,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*d-2)= torshell(3*b,3*d-2)+ m*dnz2 + n*dmz2
C
C dy4dz2,dz2dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmuz2 + zpp*(ze*vixen4-dmuy4)
     1          + ze*dmuy4*(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dmuy4*dmuz2 - dqy4*dqz2
         vixen3=(vixen1/vixen2)*vixen5 + dpz2*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdaz2 + zp*(ze*vixen4-dmuy4)
     1          + ze*dmuy4*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen4 + xe**2*dmuy4*dlambdaz2
         dnz2= (v*vixen6 + dvz2*duy4 - u*vixen3 - duz2*dvy4)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*d-1,3*b)= torshell(3*d-1,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*d-1)= torshell(3*b,3*d-1)+ m*dnz2 + n*dmz2
C
C dy4dx2,dx2dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmux2 + xpp*(xe*vixen4-dmuy4)
     1          + xe*dmuy4*(xe*dmux2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dmuy4*dmux2 - dqy4*dqx2
         vixen3=(vixen1/vixen2)*vixen5 + dpx2*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdax2 + xp*(xe*vixen4-dmuy4)
     1          + xe*dmuy4*(xe*dlambdax2-1.0D0-lambda) + zp*ze*vixen4 + ze**2*dmuy4*dlambdax2
         dnx2= (v*vixen6 + dvx2*duy4 - u*vixen3 - dux2*dvy4)/v**2 - 2.0D0*dvx2*n/v

         torshell(3*d-1,3*b-2)= torshell(3*d-1,3*b-2)+ m*dnx2 + n*dmx2
         torshell(3*b-2,3*d-1)= torshell(3*b-2,3*d-1)+ m*dnx2 + n*dmx2
C
C dz4dx2,dx2dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmux2 + xpp*(xe*vixen4-dmuz4)
     1          + xe*dmuz4*(xe*dmux2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dmuz4*dmux2 - dqz4*dqx2
         vixen3=(vixen1/vixen2)*vixen5 + dpx2*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambdax2 + xp*(xe*vixen4-dmuz4)
     1          + xe*dmuz4*(xe*dlambdax2-1.0D0-lambda) + yp*ye*vixen4 + ye**2*dmuz4*dlambdax2
         dnx2= (v*vixen6 + dvx2*duz4 - u*vixen3 - dux2*dvz4)/v**2 - 2.0D0*dvx2*n/v

         torshell(3*d,3*b-2)= torshell(3*d,3*b-2)+ m*dnx2 + n*dmx2
         torshell(3*b-2,3*d)= torshell(3*b-2,3*d)+ m*dnx2 + n*dmx2
C
C dz4dy2,dy2dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmuy2 + ypp*(ye*vixen4-dmuz4)
     1          + ye*dmuz4*(ye*dmuy2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dmuz4*dmuy2 - dqz4*dqy2
         vixen3=(vixen1/vixen2)*vixen5 + dpy2*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambday2 + yp*(ye*vixen4-dmuz4)
     1          + ye*dmuz4*(ye*dlambday2-1.0D0-lambda) + xp*xe*vixen4 + xe**2*dmuz4*dlambday2
         dny2= (v*vixen6 + dvy2*duz4 - u*vixen3 - duy2*dvz4)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*d,3*b-1)= torshell(3*d,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*d)= torshell(3*b-1,3*d)+ m*dny2 + n*dmy2
C
C dx4dx3,dx3dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 2.0D0*(x(b)-x(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= xpp*(xe*vixen4+dmux4) + (1.0D0+xe*dmux4)*(xe*dmux3+mu) + ypp*ye*vixen4
     1          + ye**2*dmux4*dmux3 + zpp*ze*vixen4 + ze**2*dmux4*dmux3 - dqx4*dqx3
         vixen3=(vixen1/vixen2)*vixen5 + dpx3*dqx4
         vixen6= xp*(xe*vixen4+dmux4) + (1.0D0+xe*dmux4)*(xe*dlambdax3+lambda) + vixen4*
     1          (yp*ye+zp*ze) + (ye**2+ze**2)*dmux4*dlambdax3
         dnx3= (v*vixen6 + dvx3*dux4 - u*vixen3 - dux3*dvx4)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*d-2,3*c-2)= torshell(3*d-2,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*d-2)= torshell(3*c-2,3*d-2)+ m*dnx3 + n*dmx3
C
C dy4dy3,dy3dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 2.0D0*(y(b)-y(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= ypp*(ye*vixen4+dmuy4) + (1.0D0+ye*dmuy4)*(ye*dmuy3+mu) + xpp*xe*vixen4
     1          + xe**2*dmuy4*dmuy3 + zpp*ze*vixen4 + ze**2*dmuy4*dmuy3 - dqy4*dqy3
         vixen3=(vixen1/vixen2)*vixen5 + dpy3*dqy4
         vixen6= yp*(ye*vixen4+dmuy4) + (1.0D0+ye*dmuy4)*(ye*dlambday3+lambda) + vixen4*
     1          (xp*xe+zp*ze) + (xe**2+ze**2)*dmuy4*dlambday3
         dny3= (v*vixen6 + dvy3*duy4 - u*vixen3 - duy3*dvy4)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*d-1,3*c-1)= torshell(3*d-1,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*d-1)= torshell(3*c-1,3*d-1)+ m*dny3 + n*dmy3
C
C dz4dz3,dz3dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 2.0D0*(z(b)-z(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= zpp*(ze*vixen4+dmuz4) + (1.0D0+ze*dmuz4)*(ze*dmuz3+mu) + xpp*xe*vixen4
     1          + xe**2*dmuz4*dmuz3 + ypp*ye*vixen4 + ye**2*dmuz4*dmuz3 - dqz4*dqz3
         vixen3=(vixen1/vixen2)*vixen5 + dpz3*dqz4
         vixen6= zp*(ze*vixen4+dmuz4) + (1.0D0+ze*dmuz4)*(ze*dlambdaz3+lambda) + vixen4*
     1          (xp*xe+yp*ye) + (xe**2+ye**2)*dmuz4*dlambdaz3
         dnz3= (v*vixen6 + dvz3*duz4 - u*vixen3 - duz3*dvz4)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*d,3*c)= torshell(3*d,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*d)= torshell(3*c,3*d)+ m*dnz3 + n*dmz3
C
C dx4dy3,dy3dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuy3 + ypp*(ye*vixen4+dmux4)
     1          + ye*dmux4*(ye*dmuy3+mu) + zpp*ze*vixen4 + ze**2*dmux4*dmuy3 - dqx4*dqy3
         vixen3=(vixen1/vixen2)*vixen5 + dpy3*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambday3 + yp*(ye*vixen4+dmux4)
     1          + ye*dmux4*(ye*dlambday3+lambda) + zp*ze*vixen4 + ze**2*dmux4*dlambday3
         dny3= (v*vixen6 + dvy3*dux4 - u*vixen3 - duy3*dvx4)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*d-2,3*c-1)= torshell(3*d-2,3*c-1) + m*dny3 + n*dmy3
         torshell(3*c-1,3*d-2)= torshell(3*c-1,3*d-2) + m*dny3 + n*dmy3
C
C dx4dz3,dz3dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuz3 + zpp*(ze*vixen4+dmux4)
     1          + ze*dmux4*(ze*dmuz3+mu) + ypp*ye*vixen4 + ye**2*dmux4*dmuz3 - dqx4*dqz3
         vixen3=(vixen1/vixen2)*vixen5 + dpz3*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambdaz3 + zp*(ze*vixen4+dmux4)
     1          + ze*dmux4*(ze*dlambdaz3+lambda) + yp*ye*vixen4 + ye**2*dmux4*dlambdaz3
         dnz3= (v*vixen6 + dvz3*dux4 - u*vixen3 - duz3*dvx4)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*d-2,3*c)= torshell(3*d-2,3*c) + m*dnz3 + n*dmz3
         torshell(3*c,3*d-2)= torshell(3*c,3*d-2) + m*dnz3 + n*dmz3
C
C dy4dz3,dz3dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmuz3 + zpp*(ze*vixen4+dmuy4)
     1          + ze*dmuy4*(ze*dmuz3+mu) + xpp*xe*vixen4 + xe**2*dmuy4*dmuz3 - dqy4*dqz3
         vixen3=(vixen1/vixen2)*vixen5 + dpz3*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdaz3 + zp*(ze*vixen4+dmuy4)
     1          + ze*dmuy4*(ze*dlambdaz3+lambda) + xp*xe*vixen4 + xe**2*dmuy4*dlambdaz3
         dnz3= (v*vixen6 + dvz3*duy4 - u*vixen3 - duz3*dvy4)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*d-1,3*c)= torshell(3*d-1,3*c) + m*dnz3 + n*dmz3
         torshell(3*c,3*d-1)= torshell(3*c,3*d-1) + m*dnz3 + n*dmz3
C
C dy4dx3,dx3dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmux3 + xpp*(xe*vixen4+dmuy4)
     1          + xe*dmuy4*(xe*dmux3+mu) + zpp*ze*vixen4 + ze**2*dmuy4*dmux3 - dqy4*dqx3
         vixen3=(vixen1/vixen2)*vixen5 + dpx3*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdax3 + xp*(xe*vixen4+dmuy4)
     1          + xe*dmuy4*(xe*dlambdax3+lambda) + zp*ze*vixen4 + ze**2*dmuy4*dlambdax3
         dnx3= (v*vixen6 + dvx3*duy4 - u*vixen3 - dux3*dvy4)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*d-1,3*c-2)= torshell(3*d-1,3*c-2) + m*dnx3 + n*dmx3
         torshell(3*c-2,3*d-1)= torshell(3*c-2,3*d-1) + m*dnx3 + n*dmx3
C
C dz4dx3,dx3dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmux3 + xpp*(xe*vixen4+dmuz4)
     1          + xe*dmuz4*(xe*dmux3+mu) + ypp*ye*vixen4 + ye**2*dmuz4*dmux3 - dqz4*dqx3
         vixen3=(vixen1/vixen2)*vixen5 + dpx3*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambdax3 + xp*(xe*vixen4+dmuz4)
     1          + xe*dmuz4*(xe*dlambdax3+lambda) + yp*ye*vixen4 + ye**2*dmuz4*dlambdax3
         dnx3= (v*vixen6 + dvx3*duz4 - u*vixen3 - dux3*dvz4)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*d,3*c-2)= torshell(3*d,3*c-2) + m*dnx3 + n*dmx3
         torshell(3*c-2,3*d)= torshell(3*c-2,3*d) + m*dnx3 + n*dmx3
C
C dz4dy3,dy3dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmuy3 + ypp*(ye*vixen4+dmuz4)
     1          + ye*dmuz4*(ye*dmuy3+mu) + xpp*xe*vixen4 + xe**2*dmuz4*dmuy3 - dqz4*dqy3
         vixen3=(vixen1/vixen2)*vixen5 + dpy3*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambday3 + yp*(ye*vixen4+dmuz4)
     1          + ye*dmuz4*(ye*dlambday3+lambda) + xp*xe*vixen4 + xe**2*dmuz4*dlambday3
         dny3= (v*vixen6 + dvy3*duz4 - u*vixen3 - duy3*dvz4)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*d,3*c-1)= torshell(3*d,3*c-1) + m*dny3 + n*dmy3
         torshell(3*c-1,3*d)= torshell(3*c-1,3*d) + m*dny3 + n*dmy3
C
C dx3dx3
C
         n=(v*dux3 - u*dvx3)/v**2
         vixen4= 8.0D0*(xe**2-rbc2)*mucrap/rbc6 - 4.0D0*xe*(x(b)-x(d))/rbc4
         vixen5= 8.0D0*(xe**2-rbc2)*lambdacrap/rbc6 - 4.0D0*xe*(x(b)-x(a))/rbc4
         vixen7= xpp*(2.0D0*dmux3+xe*vixen4) + (mu+xe*dmux3)**2 + vixen4*(ypp*ye+zpp*ze)
     1          + dmux3**2*(ye**2+ze**2) - dqx3**2
         vixen8= xp*(2.0D0*dlambdax3+xe*vixen5) + (lambda+xe*dlambdax3)**2 + vixen5*(yp*ye+zp*ze)
     1          + dlambdax3**2*(ye**2+ze**2) - dpx3**2
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + 2.0D0*dpx3*dqx3
         vixen6= xp*(2.0D0*dmux3+xe*vixen4) + 2.0D0*(mu+xe*dmux3)*(lambda+xe*dlambdax3) 
     1          + xpp*(2.0D0*dlambdax3+xe*vixen5) + vixen4*(yp*ye+zp*ze) 
     2          + vixen5*(ypp*ye+zpp*ze) + 2.0D0*dlambdax3*dmux3*(ye**2+ze**2)
         dnx3= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*c-2,3*c-2)= torshell(3*c-2,3*c-2)+ m*dnx3 + n*dmx3
C
C dy3dy3
C
         n=(v*duy3 - u*dvy3)/v**2
         vixen4= 8.0D0*(xe**2-rbc2)*mucrap/rbc6 - 4.0D0*xe*(x(b)-x(d))/rbc4
         vixen5= 8.0D0*(xe**2-rbc2)*lambdacrap/rbc6 - 4.0D0*xe*(x(b)-x(a))/rbc4
         vixen7= ypp*(2.0D0*dmuy3+ye*vixen4) + (mu+ye*dmuy3)**2 + vixen4*(xpp*xe+zpp*ze)
     1          + dmuy3**2*(xe**2+ze**2) - dqy3**2
         vixen8= yp*(2.0D0*dlambday3+ye*vixen5) + (lambda+ye*dlambday3)**2 + vixen5*(xp*xe+zp*ze)
     1          + dlambday3**2*(xe**2+ze**2) - dpy3**2
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + 2.0D0*dpy3*dqy3
         vixen6= yp*(2.0D0*dmuy3+ye*vixen4) + 2.0D0*(mu+ye*dmuy3)*(lambda+ye*dlambday3)
     1          + ypp*(2.0D0*dlambday3+ye*vixen5) + vixen4*(xp*xe+zp*ze)
     2          + vixen5*(xpp*xe+zpp*ze) + 2.0D0*dlambday3*dmuy3*(xe**2+ze**2)
         dny3= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*c-1,3*c-1)= torshell(3*c-1,3*c-1)+ m*dny3 + n*dmy3
C
C dz3dz3
C
         n=(v*duz3 - u*dvz3)/v**2
         vixen4= 8.0D0*(xe**2-rbc2)*mucrap/rbc6 - 4.0D0*xe*(x(b)-x(d))/rbc4
         vixen5= 8.0D0*(xe**2-rbc2)*lambdacrap/rbc6 - 4.0D0*xe*(x(b)-x(a))/rbc4
         vixen7= zpp*(2.0D0*dmuz3+ze*vixen4) + (mu+ze*dmuz3)**2 + vixen4*(xpp*xe+ypp*ye)
     1          + dmuz3**2*(xe**2+ye**2) - dqz3**2
         vixen8= zp*(2.0D0*dlambdaz3+ze*vixen5) + (lambda+ze*dlambdaz3)**2 + vixen5*(xp*xe+yp*ye)
     1          + dlambdaz3**2*(xe**2+ye**2) - dpz3**2
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + 2.0D0*dpz3*dqz3
         vixen6= zp*(2.0D0*dmuz3+ze*vixen4) + 2.0D0*(mu+ze*dmuz3)*(lambda+ze*dlambdaz3)
     1          + zpp*(2.0D0*dlambdaz3+ze*vixen5) + vixen4*(xp*xe+yp*ye)
     2          + vixen5*(xpp*xe+ypp*ye) + 2.0D0*dlambdaz3*dmuz3*(xe**2+ye**2)
         dnz3= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*c,3*c)= torshell(3*c,3*c)+ m*dnz3 + n*dmz3
C
C dx3dy3,dy3dx3
C
         n=(v*dux3 - u*dvx3)/v**2
         vixen4= 8.0D0*mucrap*xe*ye/rbc6 - 2.0D0*(ye*(x(b)-x(d))+xe*(y(b)-y(d)))/rbc4
         vixen5= 8.0D0*mucrap*xe*ye/rbc6 - 2.0D0*(ye*(x(b)-x(a))+xe*(y(b)-y(a)))/rbc4
         vixen7= xpp*(dmuy3+xe*vixen4) + (mu+xe*dmux3)*xe*dmuy3 + ypp*(ye*vixen4+dmux3)
     1          + ye*dmux3*(mu+ye*dmuy3) + zpp*ze*vixen4 + ze**2*dmux3*dmuy3 - dqx3*dqy3
         vixen8= xp*(dlambday3+xe*vixen5) + (lambda+xe*dlambdax3)*xe*dlambday3 
     1          + yp*(ye*vixen5 + dlambdax3) + ye*dlambdax3*(lambda+ye*dlambday3) 
     2          + zp*ze*vixen5 + ze**2*dlambdax3*dlambday3 - dpx3*dpy3
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + dpx3*dqy3 + dpy3*dqx3
         vixen6= xp*(dmuy3+xe*vixen4) + (mu+xe*dmux3)*xe*dlambday3 + xpp*(dlambday3+xe*vixen5)
     1          + (lambda+xe*dlambdax3)*xe*dmuy3 + yp*(ye*vixen4+dmux3) 
     2          + ye*dmux3*(ye*dlambday3+lambda) + ypp*(ye*vixen5 + dlambdax3) 
     3          + ye*dlambdax3*(mu+ye*dmuy3) + zp*ze*vixen4 + zpp*ze*vixen5
     4          + ze**2*(dlambdax3*dmuy3+dmux3*dlambday3)
         dny3= (v*vixen6 + dvy3*dux3 - u*vixen3 - duy3*dvx3)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*c-2,3*c-1)= torshell(3*c-2,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*c-2)= torshell(3*c-1,3*c-2)+ m*dny3 + n*dmy3
C
C dx3dz3,dz3dx3
C
         n=(v*dux3 - u*dvx3)/v**2
         vixen4= 8.0D0*mucrap*xe*ze/rbc6 - 2.0D0*(ze*(x(b)-x(d))+xe*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*mucrap*xe*ze/rbc6 - 2.0D0*(ze*(x(b)-x(a))+xe*(z(b)-z(a)))/rbc4
         vixen7= xpp*(dmuz3+xe*vixen4) + (mu+xe*dmux3)*xe*dmuz3 + zpp*(ze*vixen4+dmux3)
     1          + ze*dmux3*(mu+ze*dmuz3) + ypp*ye*vixen4 + ye**2*dmux3*dmuz3 - dqx3*dqz3
         vixen8= xp*(dlambdaz3+xe*vixen5) + (lambda+xe*dlambdax3)*xe*dlambdaz3
     1          + zp*(ze*vixen5 + dlambdax3) + ze*dlambdax3*(lambda+ze*dlambdaz3)
     2          + yp*ye*vixen5 + ye**2*dlambdax3*dlambdaz3 - dpx3*dpz3
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + dpx3*dqz3 + dpz3*dqx3
         vixen6= xp*(dmuz3+xe*vixen4) + (mu+xe*dmux3)*xe*dlambdaz3 + xpp*(dlambdaz3+xe*vixen5)
     1          + (lambda+xe*dlambdax3)*xe*dmuz3 + zp*(ze*vixen4+dmux3)
     2          + ze*dmux3*(ze*dlambdaz3+lambda) + zpp*(ze*vixen5 + dlambdax3)
     3          + ze*dlambdax3*(mu+ze*dmuz3) + yp*ye*vixen4 + ypp*ye*vixen5
     4          + ye**2*(dlambdax3*dmuz3+dmux3*dlambdaz3)
         dnz3= (v*vixen6 + dvz3*dux3 - u*vixen3 - duz3*dvx3)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*c-2,3*c)= torshell(3*c-2,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*c-2)= torshell(3*c,3*c-2)+ m*dnz3 + n*dmz3
C
C dy3dz3,dz3dy3
C
         n=(v*duy3 - u*dvy3)/v**2
         vixen4= 8.0D0*mucrap*ye*ze/rbc6 - 2.0D0*(ze*(y(b)-y(d))+ye*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*mucrap*ye*ze/rbc6 - 2.0D0*(ze*(y(b)-y(a))+ye*(z(b)-z(a)))/rbc4
         vixen7= ypp*(dmuz3+ye*vixen4) + (mu+ye*dmuy3)*ye*dmuz3 + zpp*(ze*vixen4+dmuy3)
     1          + ze*dmuy3*(mu+ze*dmuz3) + xpp*xe*vixen4 + xe**2*dmuy3*dmuz3 - dqy3*dqz3
         vixen8= yp*(dlambdaz3+ye*vixen5) + (lambda+ye*dlambday3)*ye*dlambdaz3
     1          + zp*(ze*vixen5 + dlambday3) + ze*dlambday3*(lambda+ze*dlambdaz3)
     2          + xp*xe*vixen5 + xe**2*dlambday3*dlambdaz3 - dpy3*dpz3
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + dpy3*dqz3 + dpz3*dqy3
         vixen6= yp*(dmuz3+ye*vixen4) + (mu+ye*dmuy3)*ye*dlambdaz3 + ypp*(dlambdaz3+ye*vixen5)
     1          + (lambda+ye*dlambday3)*ye*dmuz3 + zp*(ze*vixen4+dmuy3)
     2          + ze*dmuy3*(ze*dlambdaz3+lambda) + zpp*(ze*vixen5 + dlambday3)
     3          + ze*dlambday3*(mu+ze*dmuz3) + xp*xe*vixen4 + xpp*xe*vixen5
     4          + xe**2*(dlambday3*dmuz3+dmuy3*dlambdaz3)
         dnz3= (v*vixen6 + dvz3*duy3 - u*vixen3 - duz3*dvy3)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*c-1,3*c)= torshell(3*c-1,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*c-1)= torshell(3*c,3*c-1)+ m*dnz3 + n*dmz3
C
C dx2dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= -2.0D0/rbc2 - (4.0D0*(x(b)-x(c))*(x(c)+x(d)+2.0D0*x(b))-
     1          2.0D0*mucrap)/rbc4 + 8.0D0*(x(b)-x(c))**2*mucrap/rbc6
         vixen5= -2.0D0/rbc2 - (4.0D0*(x(b)-x(c))*(x(c)+x(a)+2.0D0*x(b))-
     1    2.0D0*lambdacrap)/rbc4 + 8.0D0*(x(b)-x(c))**2*lambdacrap/rbc6
         vixen7= xp*(xe*vixen5-2.0D0*dlambdax2) + (xe*dlambdax2
     1    -1.0D0-lambda)**2 + (yp*ye+zp*ze)*vixen5 - dpx2**2
     2    + (ye**2+ze**2)*dlambdax2**2
         vixen8= xpp*(xe*vixen4-2.0D0*dmux2) + (xe*dmux2-1.0D0-mu)
     1    **2 + (ypp*ye+zpp*ze)*vixen4 - dqx2**2
     2    + (ye**2+ze**2)*dmux2**2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8+2.0D0*dpx2*dqx2
         vixen6= xp*(xe*vixen4-2.0D0*dmux2) + 2.0D0*(xe*dmux2-1-mu)
     1    *(xe*dlambdax2-1.0D0-lambda) + xpp*(xe*vixen5-2.0D0*
     2    dlambdax2) + vixen4*(yp*ye+zp*ze) + vixen5*
     3    (ypp*ye+zpp*ze) + 2.0D0*dmux2*dlambdax2*(ye
     4    **2+ze**2)
         dnx2= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvx2*n/v

         torshell(3*b-2,3*b-2)= torshell(3*b-2,3*b-2)+ m*dnx2 + n*dmx2
C
C dy2dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= -2.0D0/rbc2 - (4.0D0*(y(b)-y(c))*(y(c)+y(d)+2.0D0*y(b))-
     1          2.0D0*mucrap)/rbc4 + 8.0D0*(y(b)-y(c))**2*mucrap/rbc6
         vixen5= -2.0D0/rbc2 - (4.0D0*(y(b)-y(c))*(y(c)+y(a)+2.0D0*y(b))-
     1    2.0D0*lambdacrap)/rbc4 + 8.0D0*(y(b)-y(c))**2*lambdacrap/rbc6
         vixen7= yp*(ye*vixen5-2.0D0*dlambday2) + (ye*dlambday2
     1    -1.0D0-lambda)**2 + (xp*xe+zp*ze)*vixen5 - dpy2**2
     2    + (xe**2+ze**2)*dlambday2**2
         vixen8= ypp*(ye*vixen4-2.0D0*dmuy2) + (ye*dmuy2-1.0D0-mu)
     1    **2 + (xpp*xe+zpp*ze)*vixen4 - dqy2**2
     2    + (xe**2+ze**2)*dmuy2**2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8+2.0D0*dpy2*dqy2
         vixen6= yp*(ye*vixen4-2.0D0*dmuy2) + 2.0D0*(ye*dmuy2-1.0D0-mu)
     1    *(ye*dlambday2-1.0D0-lambda) + ypp*(ye*vixen5-2.0D0*
     2    dlambday2) + vixen4*(xp*xe+zp*ze) + vixen5*
     3    (xpp*xe+zpp*ze) + 2.0D0*dmuy2*dlambday2*(xe
     4    **2+ze**2)
         dny2= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*b-1,3*b-1)= torshell(3*b-1,3*b-1)+ m*dny2 + n*dmy2
C
C dz2dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= -2.0D0/rbc2 - (4.0D0*(z(b)-z(c))*(z(c)+z(d)+2.0D0*z(b))-
     1          2.0D0*mucrap)/rbc4 + 8.0D0*(z(b)-z(c))**2*mucrap/rbc6
         vixen5= -2.0D0/rbc2 - (4.0D0*(z(b)-z(c))*(z(c)+z(a)+2.0D0*z(b))-
     1    2.0D0*lambdacrap)/rbc4 + 8.0D0*(z(b)-z(c))**2*lambdacrap/rbc6
         vixen7= zp*(ze*vixen5-2.0D0*dlambdaz2) + (ze*dlambdaz2
     1    -1.0D0-lambda)**2 + (xp*xe+yp*ye)*vixen5 - dpz2**2
     2    + (xe**2+ye**2)*dlambdaz2**2
         vixen8= zpp*(ze*vixen4-2.0D0*dmuz2) + (ze*dmuz2-1.0D0-mu)
     1    **2 + (xpp*xe+ypp*ye)*vixen4 - dqz2**2
     2    + (xe**2+ye**2)*dmuz2**2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8+2.0D0*dpz2*dqz2
         vixen6= zp*(ze*vixen4-2.0D0*dmuz2) + 2.0D0*(ze*dmuz2-1.0D0-mu)
     1    *(ze*dlambdaz2-1.0D0-lambda) + zpp*(ze*vixen5-2.0D0*
     2    dlambdaz2) + vixen4*(xp*xe+yp*ye) + vixen5*
     3    (xpp*xe+ypp*ye) + 2.0D0*dmuz2*dlambdaz2*(xe
     4    **2+ye**2)
         dnz2= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*b,3*b)= torshell(3*b,3*b)+ m*dnz2 + n*dmz2
C
C dx2dy2,dy2dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= (-2.0D0*(y(b)-y(c))*(x(c)+x(d)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(y(c)+y(d)-2.0D0*y(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(y(b)-y(c))*mucrap
     2          /rbc6
         vixen5= (-2.0D0*(y(b)-y(c))*(x(c)+x(a)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(y(c)+y(a)-2.0D0*y(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(y(b)-y(c))
     2          *lambdacrap/rbc6
         vixen7= xp*(xe*vixen5-dlambday2) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambday2 + yp*(ye
     1          *vixen5-dlambdax2) + ye*dlambdax2*(ye*dlambday2-1.0D0-lambda) + zp*ze*vixen5
     2          + ze**2*dlambdax2*dlambday2 - dpx2*dpy2
         vixen8= xpp*(xe*vixen4-dmuy2) + (xe*dmux2-1.0D0-mu)*xe*dmuy2 + ypp*(ye
     1          *vixen4-dmux2) + ye*dmux2*(ye*dmuy2-1.0D0-mu) + zpp*ze*vixen4
     2          + ze**2*dmux2*dmuy2 - dqx2*dqy2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqy2 + dqx2*dpy2
         vixen6= xp*(xe*vixen4-dmuy2) + (xe*dmux2-1.0D0-mu)*xe*dlambday2 + xpp*(xe*vixen5-dlambday2)
     1          + (xe*dlambdax2-1.0D0-lambda)*xe*dmuy2 + yp*(ye*vixen4-dmux2) + ye*dmux2
     2          *(ye*dlambday2-1.0D0-lambda) + ypp*(ye*vixen5-dlambdax2) + ye*dlambdax2
     3          *(ye*dmuy2-1.0D0-mu) + zp*ze*vixen4 + zpp*ze*vixen5 + ze**2*(dmux2*dlambday2
     4          +dmuy2*dlambdax2)
         dny2= (v*vixen6 + dvy2*dux2 - u*vixen3 - dvx2*duy2)/v**2 - 2.0D0*dvy2*n/v

         torshell(3*b-2,3*b-1)= torshell(3*b-2,3*b-1)+ m*dny2 + n*dmy2
         torshell(3*b-1,3*b-2)= torshell(3*b-1,3*b-2)+ m*dny2 + n*dmy2
C
C dx2dz2,dz2dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= (-2.0D0*(z(b)-z(c))*(x(c)+x(d)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(z(c)+z(d)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(z(b)-z(c))*mucrap
     2          /rbc6
         vixen5= (-2.0D0*(z(b)-z(c))*(x(c)+x(a)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(z(c)+z(a)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(z(b)-z(c))
     2          *lambdacrap/rbc6
         vixen7= xp*(xe*vixen5-dlambdaz2) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambdaz2 + zp*(ze
     1          *vixen5-dlambdax2) + ze*dlambdax2*(ze*dlambdaz2-1.0D0-lambda) + yp*ye*vixen5
     2          + ye**2*dlambdax2*dlambdaz2 - dpx2*dpz2
         vixen8= xpp*(xe*vixen4-dmuz2) + (xe*dmux2-1.0D0-mu)*xe*dmuz2 + zpp*(ze
     1          *vixen4-dmux2) + ze*dmux2*(ze*dmuz2-1.0D0-mu) + ypp*ye*vixen4
     2          + ye**2*dmux2*dmuz2 - dqx2*dqz2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqz2 + dqx2*dpz2
         vixen6= xp*(xe*vixen4-dmuz2) + (xe*dmux2-1.0D0-mu)*xe*dlambdaz2 + xpp*(xe*vixen5-dlambdaz2)
     1          + (xe*dlambdax2-1.0D0-lambda)*xe*dmuz2 + zp*(ze*vixen4-dmux2) + ze*dmux2
     2          *(ze*dlambdaz2-1.0D0-lambda) + zpp*(ze*vixen5-dlambdax2) + ze*dlambdax2
     3          *(ze*dmuz2-1.0D0-mu) + yp*ye*vixen4 + ypp*ye*vixen5 + ye**2*(dmux2*dlambdaz2
     4          +dmuz2*dlambdax2)
         dnz2= (v*vixen6 + dvz2*dux2 - u*vixen3 - dvx2*duz2)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*b-2,3*b)= torshell(3*b-2,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*b-2)= torshell(3*b,3*b-2)+ m*dnz2 + n*dmz2
C
C dy2dz2,dz2dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= (-2.0D0*(z(b)-z(c))*(y(c)+y(d)-2.0D0*y(b))-2.0D0*(y(b)-y(c))*(z(c)+z(d)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(y(b)-y(c))*(z(b)-z(c))*mucrap
     2          /rbc6
         vixen5= (-2.0D0*(z(b)-z(c))*(y(c)+y(a)-2.0D0*y(b))-2.0D0*(y(b)-y(c))*(z(c)+z(a)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(y(b)-y(c))*(z(b)-z(c))
     2          *lambdacrap/rbc6
         vixen7= yp*(ye*vixen5-dlambdaz2) + (ye*dlambday2-1.0D0-lambda)*ye*dlambdaz2 + zp*(ze
     1          *vixen5-dlambday2) + ze*dlambday2*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen5
     2          + xe**2*dlambday2*dlambdaz2 - dpy2*dpz2
         vixen8= ypp*(ye*vixen4-dmuz2) + (ye*dmuy2-1.0D0-mu)*ye*dmuz2 + zpp*(ze
     1          *vixen4-dmuy2) + ze*dmuy2*(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4
     2          + xe**2*dmuy2*dmuz2 - dqy2*dqz2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqz2 + dqy2*dpz2
         vixen6= yp*(ye*vixen4-dmuz2) + (ye*dmuy2-1.0D0-mu)*ye*dlambdaz2 + ypp*(ye*vixen5-dlambdaz2)
     1          + (ye*dlambday2-1.0D0-lambda)*ye*dmuz2 + zp*(ze*vixen4-dmuy2) + ze*dmuy2
     2          *(ze*dlambdaz2-1.0D0-lambda) + zpp*(ze*vixen5-dlambday2) + ze*dlambday2
     3          *(ze*dmuz2-1.0D0-mu) + xp*xe*vixen4 + xpp*xe*vixen5 + xe**2*(dmuy2*dlambdaz2
     4          +dmuz2*dlambday2)
         dnz2= (v*vixen6 + dvz2*duy2 - u*vixen3 - dvy2*duz2)/v**2 - 2.0D0*dvz2*n/v

         torshell(3*b-1,3*b)= torshell(3*b-1,3*b)+ m*dnz2 + n*dmz2
         torshell(3*b,3*b-1)= torshell(3*b,3*b-1)+ m*dnz2 + n*dmz2
C
C dx2dx3,dx3dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= 1.0D0/rbc2 - (2.0D0*xe*(x(c)+x(d)-2.0D0*x(b))+2.0D0*((x(b)-x(c))*(x(b)-x(d))-mucrap))/rbc4 
     1          - 8.0D0*(x(b)-x(c))**2*mucrap/rbc6
         vixen5= 1.0D0/rbc2 - (2.0D0*xe*(x(c)+x(a)-2.0D0*x(b))+2.0D0*((x(b)-x(c))*(x(b)-x(a))-lambdacrap))/rbc4 
     1          - 8.0D0*(x(b)-x(c))**2*lambdacrap/rbc6
         vixen7= xp*(xe*vixen5+dlambdax2-dlambdax3) + (xe*dlambdax2-1.0D0-lambda)*(xe*dlambdax3+lambda)
     1          + (yp*ye+zp*ze)*vixen5 + (ye**2+ze**2)*dlambdax2*dlambdax3 - dpx2*dpx3
         vixen8= xpp*(xe*vixen4+dmux2-dmux3) + (xe*dmux2-1.0D0-mu)*(xe*dmux3+mu)
     1          + (ypp*ye+zpp*ze)*vixen4 + (ye**2+ze**2)*dmux2*dmux3 - dqx2*dqx3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqx3 + dqx2*dpx3
         vixen6= xp*(xe*vixen4+dmux2-dmux3) + (xe*dmux2-1.0D0-mu)*(xe*dlambdax3+lambda)
     1          + xpp*(xe*vixen5+dlambdax2-dlambdax3) + (xe*dlambdax2-1.0D0-lambda)*(xe*dmux3+mu)
     2          + vixen4*(yp*ye+zp*ze) + vixen5*(ypp*ye+zpp*ze) 
     3          + (dmux2*dlambdax3+dmux3*dlambdax2)*(ye**2+ze**2)
         dnx3= (v*vixen6 + dvx3*dux2 - u*vixen3 - dux3*dvx2)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*b-2,3*c-2)= torshell(3*b-2,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*b-2)= torshell(3*c-2,3*b-2)+ m*dnx3 + n*dmx3
C
C dy2dy3,dy3dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= 1.0D0/rbc2 - (2.0D0*ye*(y(c)+y(d)-2.0D0*y(b))+2.0D0*((y(b)-y(c))*(y(b)-y(d))-mucrap))/rbc4 
     1          - 8.0D0*(y(b)-y(c))**2*mucrap/rbc6
         vixen5= 1.0D0/rbc2 - (2.0D0*ye*(y(c)+y(a)-2.0D0*y(b))+2.0D0*((y(b)-y(c))*(y(b)-y(a))-lambdacrap))/rbc4 
     1          - 8.0D0*(y(b)-y(c))**2*lambdacrap/rbc6
         vixen7= yp*(ye*vixen5+dlambday2-dlambday3) + (ye*dlambday2-1.0D0-lambda)*(ye*dlambday3+lambda)
     1          + (xp*xe+zp*ze)*vixen5 + (xe**2+ze**2)*dlambday2*dlambday3 - dpy2*dpy3
         vixen8= ypp*(ye*vixen4+dmuy2-dmuy3) + (ye*dmuy2-1.0D0-mu)*(ye*dmuy3+mu)
     1          + (xpp*xe+zpp*ze)*vixen4 + (xe**2+ze**2)*dmuy2*dmuy3 - dqy2*dqy3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqy3 + dqy2*dpy3
         vixen6= yp*(ye*vixen4+dmuy2-dmuy3) + (ye*dmuy2-1.0D0-mu)*(ye*dlambday3+lambda)
     1          + ypp*(ye*vixen5+dlambday2-dlambday3) + (ye*dlambday2-1.0D0-lambda)*(ye*dmuy3+mu)
     2          + vixen4*(xp*xe+zp*ze) + vixen5*(xpp*xe+zpp*ze) 
     3          + (dmuy2*dlambday3+dmuy3*dlambday2)*(xe**2+ze**2)
         dny3= (v*vixen6 + dvy3*duy2 - u*vixen3 - duy3*dvy2)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*b-1,3*c-1)= torshell(3*b-1,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*b-1)= torshell(3*c-1,3*b-1)+ m*dny3 + n*dmy3
C
C dz2dz3,dz3dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= 1.0D0/rbc2 - (2.0D0*ze*(z(c)+z(d)-2.0D0*z(b))+2.0D0*((z(b)-z(c))*(z(b)-z(d))-mucrap))/rbc4 
     1          - 8.0D0*(z(b)-z(c))**2*mucrap/rbc6
         vixen5= 1.0D0/rbc2 - (2.0D0*ze*(z(c)+z(a)-2.0D0*z(b))+2.0D0*((z(b)-z(c))*(z(b)-z(a))-lambdacrap))/rbc4 
     1          - 8.0D0*(z(b)-z(c))**2*lambdacrap/rbc6
         vixen7= zp*(ze*vixen5+dlambdaz2-dlambdaz3) + (ze*dlambdaz2-1.0D0-lambda)*(ze*dlambdaz3+lambda)
     1          + (xp*xe+yp*ye)*vixen5 + (xe**2+ye**2)*dlambdaz2*dlambdaz3 - dpz2*dpz3
         vixen8= zpp*(ze*vixen4+dmuz2-dmuz3) + (ze*dmuz2-1.0D0-mu)*(ze*dmuz3+mu)
     1          + (xpp*xe+ypp*ye)*vixen4 + (xe**2+ye**2)*dmuz2*dmuz3 - dqz2*dqz3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpz2*dqz3 + dqz2*dpz3
         vixen6= zp*(ze*vixen4+dmuz2-dmuz3) + (ze*dmuz2-1.0D0-mu)*(ze*dlambdaz3+lambda)
     1          + zpp*(ze*vixen5+dlambdaz2-dlambdaz3) + (ze*dlambdaz2-1-lambda)*(ze*dmuz3+mu)
     2          + vixen4*(xp*xe+yp*ye) + vixen5*(xpp*xe+ypp*ye) 
     3          + (dmuz2*dlambdaz3+dmuz3*dlambdaz2)*(xe**2+ye**2)
         dnz3= (v*vixen6 + dvz3*duz2 - u*vixen3 - duz3*dvz2)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*b,3*c)= torshell(3*b,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*b)= torshell(3*c,3*b)+ m*dnz3 + n*dmz3
C
C dx2dy3,dy3dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= 8.0D0*(x(b)-x(c))*ye*mucrap/rbc6 - 2.0D0*((x(c)+x(d)-2*x(b))*ye
     1          +(x(b)-x(c))*(y(b)-y(d)))/rbc4
         vixen5= 8.0D0*(x(b)-x(c))*ye*lambdacrap/rbc6 - 2.0D0*((x(c)+x(a)-2*x(b))*ye
     1          +(x(b)-x(c))*(y(b)-y(a)))/rbc4
         vixen7= xp*(xe*vixen5-dlambday3) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambday3 + 
     1          yp*(ye*vixen5+dlambdax2) + ye*dlambdax2*(ye*dlambday3+lambda) 
     2          + zp*ze*vixen5 + ze**2*dlambdax2*dlambday3 - dpx2*dpy3
         vixen8= xpp*(xe*vixen4-dmuy3) + (xe*dmux2-1.0D0-mu)*xe*dmuy3 + 
     1          ypp*(ye*vixen4+dmux2) + ye*dmux2*(ye*dmuy3+mu) + 
     2          zpp*ze*vixen4 + ze**2*dmux2*dmuy3 - dqx2*dqy3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqy3 + dqx2*dpy3
         vixen6= xp*(xe*vixen4-dmuy3) + (xe*dmux2-1.0D0-mu)*xe*dlambday3
     1          + xpp*(xe*vixen5-dlambday3) + (xe*dlambdax2-1.0D0-lambda)*xe*dmuy3
     2          + yp*(ye*vixen4+dmux2) + ye*dmux2*(ye*dlambday3+lambda)
     3          + ypp*(ye*vixen5+dlambdax2) + ye*dlambdax2*(ye*dmuy3+mu)
     4          + zp*ze*vixen4 + ze**2*(dmux2*dlambday3+dlambdax2*dmuy3) + zpp*ze*vixen5
         dny3= (v*vixen6 + dvy3*dux2 - u*vixen3 - duy3*dvx2)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*b-2,3*c-1)= torshell(3*b-2,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*b-2)= torshell(3*c-1,3*b-2)+ m*dny3 + n*dmy3
C
C dx2dz3,dz3dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= 8.0D0*(x(b)-x(c))*ze*mucrap/rbc6 - 2.0D0*((x(c)+x(d)-2*x(b))*ze
     1          +(x(b)-x(c))*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*(x(b)-x(c))*ze*lambdacrap/rbc6 - 2.0D0*((x(c)+x(a)-2*x(b))*ze
     1          +(x(b)-x(c))*(z(b)-z(a)))/rbc4
         vixen7= xp*(xe*vixen5-dlambdaz3) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambdaz3 + 
     1          zp*(ze*vixen5+dlambdax2) + ze*dlambdax2*(ze*dlambdaz3+lambda) 
     2          + yp*ye*vixen5 + ye**2*dlambdax2*dlambdaz3 - dpx2*dpz3
         vixen8= xpp*(xe*vixen4-dmuz3) + (xe*dmux2-1.0D0-mu)*xe*dmuz3 + 
     1          zpp*(ze*vixen4+dmux2) + ze*dmux2*(ze*dmuz3+mu) + 
     2          ypp*ye*vixen4 + ye**2*dmux2*dmuz3 - dqx2*dqz3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqz3 + dqx2*dpz3
         vixen6= xp*(xe*vixen4-dmuz3) + (xe*dmux2-1.0D0-mu)*xe*dlambdaz3
     1          + xpp*(xe*vixen5-dlambdaz3) + (xe*dlambdax2-1.0D0-lambda)*xe*dmuz3
     2          + zp*(ze*vixen4+dmux2) + ze*dmux2*(ze*dlambdaz3+lambda)
     3          + zpp*(ze*vixen5+dlambdax2) + ze*dlambdax2*(ze*dmuz3+mu)
     4          + yp*ye*vixen4 + ye**2*(dmux2*dlambdaz3+dlambdax2*dmuz3) + ypp*ye*vixen5
         dnz3= (v*vixen6 + dvz3*dux2 - u*vixen3 - duz3*dvx2)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*b-2,3*c)= torshell(3*b-2,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*b-2)= torshell(3*c,3*b-2)+ m*dnz3 + n*dmz3
C
C dy2dz3,dz3dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= 8.0D0*(y(b)-y(c))*ze*mucrap/rbc6 - 2.0D0*((y(c)+y(d)-2*y(b))*ze
     1          +(y(b)-y(c))*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*(y(b)-y(c))*ze*lambdacrap/rbc6 - 2.0D0*((y(c)+y(a)-2*y(b))*ze
     1          +(y(b)-y(c))*(z(b)-z(a)))/rbc4
         vixen7= yp*(ye*vixen5-dlambdaz3) + (ye*dlambday2-1.0D0-lambda)*ye*dlambdaz3 +
     1          zp*(ze*vixen5+dlambday2) + ze*dlambday2*(ze*dlambdaz3+lambda)
     2          + xp*xe*vixen5 + xe**2*dlambday2*dlambdaz3 - dpy2*dpz3
         vixen8= ypp*(ye*vixen4-dmuz3) + (ye*dmuy2-1.0D0-mu)*ye*dmuz3 +
     1          zpp*(ze*vixen4+dmuy2) + ze*dmuy2*(ze*dmuz3+mu) +
     2          xpp*xe*vixen4 + xe**2*dmuy2*dmuz3 - dqy2*dqz3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqz3 + dqy2*dpz3
         vixen6= yp*(ye*vixen4-dmuz3) + (ye*dmuy2-1.0D0-mu)*ye*dlambdaz3
     1          + ypp*(ye*vixen5-dlambdaz3) + (ye*dlambday2-1.0D0-lambda)*ye*dmuz3
     2          + zp*(ze*vixen4+dmuy2) + ze*dmuy2*(ze*dlambdaz3+lambda)
     3          + zpp*(ze*vixen5+dlambday2) + ze*dlambday2*(ze*dmuz3+mu)
     4          + xp*xe*vixen4 + xe**2*(dmuy2*dlambdaz3+dlambday2*dmuz3) + xpp*xe*vixen5
         dnz3= (v*vixen6 + dvz3*duy2 - u*vixen3 - duz3*dvy2)/v**2 - 2.0D0*dvz3*n/v

         torshell(3*b-1,3*c)= torshell(3*b-1,3*c)+ m*dnz3 + n*dmz3
         torshell(3*c,3*b-1)= torshell(3*c,3*b-1)+ m*dnz3 + n*dmz3
C
C dy2dx3,dx3dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= 8.0D0*(y(b)-y(c))*xe*mucrap/rbc6 - 2.0D0*((y(c)+y(d)-2*y(b))*xe
     1          +(y(b)-y(c))*(x(b)-x(d)))/rbc4
         vixen5= 8.0D0*(y(b)-y(c))*xe*lambdacrap/rbc6 - 2.0D0*((y(c)+y(a)-2*y(b))*xe
     1          +(y(b)-y(c))*(x(b)-x(a)))/rbc4
         vixen7= yp*(ye*vixen5-dlambdax3) + (ye*dlambday2-1.0D0-lambda)*ye*dlambdax3 +
     1          xp*(xe*vixen5+dlambday2) + xe*dlambday2*(xe*dlambdax3+lambda)
     2          + zp*ze*vixen5 + ze**2*dlambday2*dlambdax3 - dpy2*dpx3
         vixen8= ypp*(ye*vixen4-dmux3) + (ye*dmuy2-1.0D0-mu)*ye*dmux3 +
     1          xpp*(xe*vixen4+dmuy2) + xe*dmuy2*(xe*dmux3+mu) +
     2          zpp*ze*vixen4 + ze**2*dmuy2*dmux3 - dqy2*dqx3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqx3 + dqy2*dpx3
         vixen6= yp*(ye*vixen4-dmux3) + (ye*dmuy2-1.0D0-mu)*ye*dlambdax3
     1          + ypp*(ye*vixen5-dlambdax3) + (ye*dlambday2-1.0D0-lambda)*ye*dmux3
     2          + xp*(xe*vixen4+dmuy2) + xe*dmuy2*(xe*dlambdax3+lambda)
     3          + xpp*(xe*vixen5+dlambday2) + xe*dlambday2*(xe*dmux3+mu)
     4          + zp*ze*vixen4 + ze**2*(dmuy2*dlambdax3+dlambday2*dmux3) + zpp*ze*vixen5
         dnx3= (v*vixen6 + dvx3*duy2 - u*vixen3 - dux3*dvy2)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*b-1,3*c-2)= torshell(3*b-1,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*b-1)= torshell(3*c-2,3*b-1)+ m*dnx3 + n*dmx3
C
C dz2dx3,dx3dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= 8.0D0*(z(b)-z(c))*xe*mucrap/rbc6 - 2.0D0*((z(c)+z(d)-2*z(b))*xe
     1          +(z(b)-z(c))*(x(b)-x(d)))/rbc4
         vixen5= 8.0D0*(z(b)-z(c))*xe*lambdacrap/rbc6 - 2.0D0*((z(c)+z(a)-2*z(b))*xe
     1          +(z(b)-z(c))*(x(b)-x(a)))/rbc4
         vixen7= zp*(ze*vixen5-dlambdax3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dlambdax3 +
     1          xp*(xe*vixen5+dlambdaz2) + xe*dlambdaz2*(xe*dlambdax3+lambda)
     2          + yp*ye*vixen5 + ye**2*dlambdaz2*dlambdax3 - dpz2*dpx3
         vixen8= zpp*(ze*vixen4-dmux3) + (ze*dmuz2-1.0D0-mu)*ze*dmux3 +
     1          xpp*(xe*vixen4+dmuz2) + xe*dmuz2*(xe*dmux3+mu) +
     2          ypp*ye*vixen4 + ye**2*dmuz2*dmux3 - dqz2*dqx3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpz2*dqx3 + dqz2*dpx3
         vixen6= zp*(ze*vixen4-dmux3) + (ze*dmuz2-1.0D0-mu)*ze*dlambdax3
     1          + zpp*(ze*vixen5-dlambdax3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dmux3
     2          + xp*(xe*vixen4+dmuz2) + xe*dmuz2*(xe*dlambdax3+lambda)
     3          + xpp*(xe*vixen5+dlambdaz2) + xe*dlambdaz2*(xe*dmux3+mu)
     4          + yp*ye*vixen4 + ye**2*(dmuz2*dlambdax3+dlambdaz2*dmux3) + ypp*ye*vixen5
         dnx3= (v*vixen6 + dvx3*duz2 - u*vixen3 - dux3*dvz2)/v**2 - 2.0D0*dvx3*n/v

         torshell(3*b,3*c-2)= torshell(3*b,3*c-2)+ m*dnx3 + n*dmx3
         torshell(3*c-2,3*b)= torshell(3*c-2,3*b)+ m*dnx3 + n*dmx3
C
C dz2dy3,dy3dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= 8.0D0*(z(b)-z(c))*ye*mucrap/rbc6 - 2.0D0*((z(c)+z(d)-2*z(b))*ye
     1          +(z(b)-z(c))*(y(b)-y(d)))/rbc4
         vixen5= 8.0D0*(z(b)-z(c))*ye*lambdacrap/rbc6 - 2.0D0*((z(c)+z(a)-2*z(b))*ye
     1          +(z(b)-z(c))*(y(b)-y(a)))/rbc4
         vixen7= zp*(ze*vixen5-dlambday3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dlambday3 +
     1          yp*(ye*vixen5+dlambdaz2) + ye*dlambdaz2*(ye*dlambday3+lambda)
     2          + xp*xe*vixen5 + xe**2*dlambdaz2*dlambday3 - dpz2*dpy3
         vixen8= zpp*(ze*vixen4-dmuy3) + (ze*dmuz2-1.0D0-mu)*ze*dmuy3 +
     1          ypp*(ye*vixen4+dmuz2) + ye*dmuz2*(ye*dmuy3+mu) +
     2          xpp*xe*vixen4 + xe**2*dmuz2*dmuy3 - dqz2*dqy3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpz2*dqy3 + dqz2*dpy3
         vixen6= zp*(ze*vixen4-dmuy3) + (ze*dmuz2-1.0D0-mu)*ze*dlambday3
     1          + zpp*(ze*vixen5-dlambday3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dmuy3
     2          + yp*(ye*vixen4+dmuz2) + ye*dmuz2*(ye*dlambday3+lambda)
     3          + ypp*(ye*vixen5+dlambdaz2) + ye*dlambdaz2*(ye*dmuy3+mu)
     4          + xp*xe*vixen4 + xe**2*(dmuz2*dlambday3+dlambdaz2*dmuy3) + xpp*xe*vixen5
         dny3= (v*vixen6 + dvy3*duz2 - u*vixen3 - duy3*dvz2)/v**2 - 2.0D0*dvy3*n/v

         torshell(3*b,3*c-1)= torshell(3*b,3*c-1)+ m*dny3 + n*dmy3
         torshell(3*c-1,3*b)= torshell(3*c-1,3*b)+ m*dny3 + n*dmy3



      END DO

C      PRINT *,' '
C      PRINT *,'TORSION TEST....'
C      PRINT *,' '
C      DO i=1,atoms
C         vixen1=dtorsEbydx(i)-torscrap(3*i-2)
C         WRITE(*,FMT='(A,I2,A,F20.10,5X,F20.10,A,F20.10)') 'x(',i,')',dtorsEbydx(i),torscrap(3*i-2),' Diff ',vixen1
C         vixen1=dtorsEbydy(i)-torscrap(3*i-1)
C         WRITE(*,FMT='(A,I2,A,F20.10,5X,F20.10,A,F20.10)') 'y(',i,')',dtorsEbydy(i),torscrap(3*i-1),' Diff ',vixen1
C         vixen1=dtorsEbydz(i)-torscrap(3*i)
C         WRITE(*,FMT='(A,I2,A,F20.10,5X,F20.10,A,F20.10)') 'z(',i,')',dtorsEbydz(i),torscrap(3*i),' Diff ',vixen1
C      END DO
C
C************************************
C IMPROPER TORSION ANGLE DERIVATIVES*
C************************************
C
      DO i=1,imp
         a=ia1(i)
         b=ia2(i)
         c=ia3(i)
         d=ia4(i)
C         PRINT *,a,b,c,d,' : ',typech(a),'-',typech(b),'-',typech(c),'-',typech(d)

         rbc=r(b,c)
         rbc2=rbc**2
         rbc4=rbc**4
         rbc6=rbc**6
         xe=x(c)-x(b)
         ye=y(c)-y(b)
         ze=z(c)-z(b)
         lambda=((x(b)-x(a))*xe+(y(b)-y(a))*ye+(z(b)-z(a))*ze)/rbc2
         mu=    ((x(b)-x(d))*xe+(y(b)-y(d))*ye+(z(b)-z(d))*ze)/rbc2
         xp= x(a)+lambda*x(c)-(1.0D0+lambda)*x(b)
         yp= y(a)+lambda*y(c)-(1.0D0+lambda)*y(b)
         zp= z(a)+lambda*z(c)-(1.0D0+lambda)*z(b)
         xpp= x(d)+mu*x(c)-(1.0D0+mu)*x(b)
         ypp= y(d)+mu*y(c)-(1.0D0+mu)*y(b)
         zpp= z(d)+mu*z(c)-(1.0D0+mu)*z(b)
         u= xp*xpp + yp*ypp + zp*zpp
         v= SQRT(xp**2 + yp**2 + zp**2)*SQRT(xpp**2 + ypp**2 + zpp**2)

         lambdacrap= (x(b)-x(a))*xe + (y(b)-y(a))*ye + (z(b)-z(a))*ze
         mucrap= (x(b)-x(d))*xe + (y(b)-y(d))*ye + (z(b)-z(d))*ze

         dlambdax1= (x(b)-x(c))/rbc2
         dlambday1= (y(b)-y(c))/rbc2
         dlambdaz1= (z(b)-z(c))/rbc2
         dlambdax2= (x(c)+x(a)-2.0D0*x(b))/rbc2 - 2.0D0*(x(b)-x(c))*lambdacrap/rbc4 
         dlambday2= (y(c)+y(a)-2.0D0*y(b))/rbc2 - 2.0D0*(y(b)-y(c))*lambdacrap/rbc4
         dlambdaz2= (z(c)+z(a)-2.0D0*z(b))/rbc2 - 2.0D0*(z(b)-z(c))*lambdacrap/rbc4
         dlambdax3= (x(b)-x(a))/rbc2 - 2.0D0*xe*lambdacrap/rbc4
         dlambday3= (y(b)-y(a))/rbc2 - 2.0D0*ye*lambdacrap/rbc4
         dlambdaz3= (z(b)-z(a))/rbc2 - 2.0D0*ze*lambdacrap/rbc4
         dmux2= (x(c)+x(d)-2.0D0*x(b))/rbc2 - 2.0D0*(x(b)-x(c))*mucrap/rbc4
         dmuy2= (y(c)+y(d)-2.0D0*y(b))/rbc2 - 2.0D0*(y(b)-y(c))*mucrap/rbc4
         dmuz2= (z(c)+z(d)-2.0D0*z(b))/rbc2 - 2.0D0*(z(b)-z(c))*mucrap/rbc4
         dmux3= (x(b)-x(d))/rbc2 - 2.0D0*xe*mucrap/rbc4
         dmuy3= (y(b)-y(d))/rbc2 - 2.0D0*ye*mucrap/rbc4
         dmuz3= (z(b)-z(d))/rbc2 - 2.0D0*ze*mucrap/rbc4
         dmux4= (x(b)-x(c))/rbc2
         dmuy4= (y(b)-y(c))/rbc2
         dmuz4= (z(b)-z(c))/rbc2

         vixen1=SQRT(xp**2 + yp**2 + zp**2)
         vixen2=SQRT(xpp**2 + ypp**2 + zpp**2)
         dpx1=(1.0D0/vixen1)*(xp*(1.0D0+xe*dlambdax1) + yp*ye*dlambdax1 + zp*ze*dlambdax1)
         dvx1=vixen2*dpx1
         dux1=xpp*(1.0D0+xe*dlambdax1) + ypp*ye*dlambdax1 + zpp*ze*dlambdax1

         dpy1=(1.0D0/vixen1)*(yp*(1.0D0+ye*dlambday1) + xp*xe*dlambday1 + zp*ze*dlambday1)
         dvy1=vixen2*dpy1
         duy1=ypp*(1.0D0+ye*dlambday1) + xpp*xe*dlambday1 + zpp*ze*dlambday1

         dpz1=(1.0D0/vixen1)*(zp*(1.0D0+ze*dlambdaz1) + xp*xe*dlambdaz1 + yp*ye*dlambdaz1)
         dvz1=vixen2*dpz1
         duz1=zpp*(1.0D0+ze*dlambdaz1) + xpp*xe*dlambdaz1 + ypp*ye*dlambdaz1

         dux2=xp*(xe*dmux2-1.0D0-mu) + xpp*(xe*dlambdax2-1.0D0-lambda) + ye*(yp*dmux2 + ypp*dlambdax2) +
     1          ze*(zp*dmux2 + zpp*dlambdax2)
         dpx2=(1.0D0/vixen1)*(xp*(xe*dlambdax2-1.0D0-lambda) + yp*ye*dlambdax2 + zp*ze*dlambdax2)
         vixen6=vixen2*dpx2
         dqx2=(1.0D0/vixen2)*(xpp*(xe*dmux2-1.0D0-mu) + ypp*ye*dmux2 + zpp*ze*dmux2)
         vixen7=vixen1*dqx2
         dvx2=vixen6+vixen7

         duy2=yp*(ye*dmuy2-1.0D0-mu) + ypp*(ye*dlambday2-1.0D0-lambda) + xe*(xp*dmuy2 + xpp*dlambday2) +
     1          ze*(zp*dmuy2 + zpp*dlambday2)
         dpy2=(1.0D0/vixen1)*(yp*(ye*dlambday2-1.0D0-lambda) + xp*xe*dlambday2 + zp*ze*dlambday2)
         vixen6=vixen2*dpy2
         dqy2=(1.0D0/vixen2)*(ypp*(ye*dmuy2-1.0D0-mu) + xpp*xe*dmuy2 + zpp*ze*dmuy2)
         vixen7=vixen1*dqy2
         dvy2=vixen6+vixen7

         duz2=zp*(ze*dmuz2-1.0D0-mu) + zpp*(ze*dlambdaz2-1.0D0-lambda) + xe*(xp*dmuz2 + xpp*dlambdaz2) +
     1          ye*(yp*dmuz2 + ypp*dlambdaz2)
         dpz2=(1.0D0/vixen1)*(zp*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*dlambdaz2 + yp*ye*dlambdaz2)
         vixen6=vixen2*dpz2
         dqz2=(1.0D0/vixen2)*(zpp*(ze*dmuz2-1.0D0-mu) + xpp*xe*dmuz2 + ypp*ye*dmuz2)
         vixen7=vixen1*dqz2
         dvz2=vixen6+vixen7

         dux3=xp*(mu + xe*dmux3) + xpp*(lambda + xe*dlambdax3) + ye*(yp*dmux3 + ypp*dlambdax3) +
     1          ze*(zp*dmux3 + zpp*dlambdax3)
         dpx3=(1.0D0/vixen1)*(xp*(lambda + xe*dlambdax3) + yp*ye*dlambdax3 + zp*ze*dlambdax3)
         vixen6=vixen2*dpx3
         dqx3=(1.0D0/vixen2)*(xpp*(mu + xe*dmux3) + ypp*ye*dmux3 + zpp*ze*dmux3)
         vixen7=vixen1*dqx3
         dvx3=vixen6+vixen7

         duy3=yp*(mu + ye*dmuy3) + ypp*(lambda + ye*dlambday3) + xe*(xp*dmuy3 + xpp*dlambday3) +
     1          ze*(zp*dmuy3 + zpp*dlambday3)
         dpy3=(1.0D0/vixen1)*(yp*(lambda + ye*dlambday3) + xp*xe*dlambday3 + zp*ze*dlambday3)
         vixen6=vixen2*dpy3
         dqy3=(1.0D0/vixen2)*(ypp*(mu + ye*dmuy3) + xpp*xe*dmuy3 + zpp*ze*dmuy3)
         vixen7=vixen1*dqy3
         dvy3=vixen6+vixen7

         duz3=zp*(mu + ze*dmuz3) + zpp*(lambda + ze*dlambdaz3) + xe*(xp*dmuz3 + xpp*dlambdaz3) +
     1          ye*(yp*dmuz3 + ypp*dlambdaz3)
         dpz3=(1.0D0/vixen1)*(zp*(lambda + ze*dlambdaz3) + xp*xe*dlambdaz3 + yp*ye*dlambdaz3)
         vixen6=vixen2*dpz3
         dqz3=(1.0D0/vixen2)*(zpp*(mu + ze*dmuz3) + xpp*xe*dmuz3 + ypp*ye*dmuz3)
         vixen7=vixen1*dqz3
         dvz3=vixen6+vixen7

         dux4=xp*(1.0D0 + xe*dmux4) + yp*ye*dmux4 + zp*ze*dmux4
         dqx4=(1.0D0/vixen2)*(xpp*(1.0D0 + xe*dmux4) + ypp*ye*dmux4 + zpp*ze*dmux4)
         dvx4=vixen1*dqx4
  
         duy4=yp*(1.0D0 + ye*dmuy4) + xp*xe*dmuy4 + zp*ze*dmuy4
         dqy4=(1.0D0/vixen2)*(ypp*(1.0D0 + ye*dmuy4) + xpp*xe*dmuy4 + zpp*ze*dmuy4)
         dvy4=vixen1*dqy4

         duz4=zp*(1.0D0 + ze*dmuz4) + xp*xe*dmuz4 + yp*ye*dmuz4
         dqz4=(1.0D0/vixen2)*(zpp*(1.0D0 + ze*dmuz4) + xpp*xe*dmuz4 + ypp*ye*dmuz4)
         dvz4=vixen1*dqz4
C
C Now second derivatives
C
         CALL hairyimphell(u,v)
C
C dx1dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= ((1.0D0+xe*dlambdax1)**2 + ye**2*dlambdax1**2 + ze**2*dlambdax1**2 - dpx1**2)/vixen1
         dnx1= -(u*vixen2/v**2)*vixen3 -2.0D0*dvx1*(v*dux1 - u*dvx1)/v**3

         imphell(3*a-2,3*a-2) =imphell(3*a-2,3*a-2) + m*dnx1 + n*dmx1
C
C dy1dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= ((1.0D0+ye*dlambday1)**2 + xe**2*dlambday1**2 + ze**2*dlambday1**2 - dpy1**2)/vixen1
         dny1= -(u*vixen2/v**2)*vixen3 -2.0D0*dvy1*(v*duy1 - u*dvy1)/v**3

         imphell(3*a-1,3*a-1) =imphell(3*a-1,3*a-1) + m*dny1 + n*dmy1
C
C dz1dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= ((1.0D0+ze*dlambdaz1)**2 + xe**2*dlambdaz1**2 + ye**2*dlambdaz1**2 - dpz1**2)/vixen1
         dnz1= -(u*vixen2/v**2)*vixen3 -2.0D0*dvz1*(v*duz1 - u*dvz1)/v**3

         imphell(3*a,3*a) =imphell(3*a,3*a) + m*dnz1 + n*dmz1
C
C dx1dy1,dy1dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= ((1.0D0+xe*dlambdax1)*xe*dlambday1 + ye*dlambdax1*(1.0D0+ye*dlambday1) 
     1          +ze**2*dlambdax1*dlambday1 - dpx1*dpy1)/vixen1
         dny1= (dux1*dvy1 - duy1*dvx1 - u*vixen2*vixen3)/v**2 - 2.0D0*dvy1*(v*dux1 - u*dvx1)/v**3
         
         imphell(3*a-2,3*a-1)= imphell(3*a-2,3*a-1) + m*dny1 + n*dmy1
         imphell(3*a-1,3*a-2)= imphell(3*a-1,3*a-2) + m*dny1 + n*dmy1
C
C dx1dz1,dz1dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= ((1.0D0+xe*dlambdax1)*xe*dlambdaz1 + ze*dlambdax1*(1.0D0+ze*dlambdaz1)
     1          +ye**2*dlambdax1*dlambdaz1 - dpx1*dpz1)/vixen1
         dnz1= (dux1*dvz1 - duz1*dvx1 - u*vixen2*vixen3)/v**2 - 2.0D0*dvz1*(v*dux1 - u*dvx1)/v**3

         imphell(3*a-2,3*a)= imphell(3*a-2,3*a) + m*dnz1 + n*dmz1
         imphell(3*a,3*a-2)= imphell(3*a,3*a-2) + m*dnz1 + n*dmz1
C
C dy1dz1,dz1dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= ((1.0D0+ye*dlambday1)*ye*dlambdaz1 + ze*dlambday1*(1.0D0+ze*dlambdaz1)
     1          +xe**2*dlambday1*dlambdaz1 - dpy1*dpz1)/vixen1
         dnz1= (duy1*dvz1 - duz1*dvy1 - u*vixen2*vixen3)/v**2 - 2.0D0*dvz1*(v*duy1 - u*dvy1)/v**3

         imphell(3*a-1,3*a)= imphell(3*a-1,3*a) + m*dnz1 + n*dmz1
         imphell(3*a,3*a-1)= imphell(3*a,3*a-1) + m*dnz1 + n*dmz1
C
C dx1dx4,dx4dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= (1.0D0+xe*dlambdax1)*(1.0D0+xe*dmux4) + (ye**2 + ze**2)*dlambdax1*dmux4
         dnx4= (v*vixen3 + dvx4*dux1 - u*dqx4*dpx1 - dvx1*dux4)/v**2 - 2.0D0*dvx4*n/v

         imphell(3*a-2,3*d-2)= imphell(3*a-2,3*d-2) + m*dnx4 + n*dmx4
         imphell(3*d-2,3*a-2)= imphell(3*d-2,3*a-2) + m*dnx4 + n*dmx4
C
C dy1dy4,dy4dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= (1.0D0+ye*dlambday1)*(1.0D0+ye*dmuy4) + (xe**2 + ze**2)*dlambday1*dmuy4
         dny4= (v*vixen3 + dvy4*duy1 - u*dqy4*dpy1 - dvy1*duy4)/v**2 - 2.0D0*dvy4*n/v

         imphell(3*a-1,3*d-1)= imphell(3*a-1,3*d-1) + m*dny4 + n*dmy4
         imphell(3*d-1,3*a-1)= imphell(3*d-1,3*a-1) + m*dny4 + n*dmy4
C
C dz1dz4,dz4dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= (1.0D0+ze*dlambdaz1)*(1.0D0+ze*dmuz4) + (ye**2 + xe**2)*dlambdaz1*dmuz4
         dnz4= (v*vixen3 + dvz4*duz1 - u*dqz4*dpz1 - dvz1*duz4)/v**2 - 2.0D0*dvz4*n/v

         imphell(3*a,3*d)= imphell(3*a,3*d) + m*dnz4 + n*dmz4
         imphell(3*d,3*a)= imphell(3*d,3*a) + m*dnz4 + n*dmz4
C
C dx1dy4,dy4dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= (1.0D0+xe*dlambdax1)*xe*dmuy4 + ye*dlambdax1*(1.0D0+ye*dmuy4) 
     1          + ze**2*dlambdax1*dmuy4
         dny4= (v*vixen3 + dvy4*dux1 - u*dqy4*dpx1 - dvx1*duy4)/v**2 - 2.0D0*dvy4*n/v

         imphell(3*a-2,3*d-1)= imphell(3*a-2,3*d-1) + m*dny4 + n*dmy4
         imphell(3*d-1,3*a-2)= imphell(3*d-1,3*a-2) + m*dny4 + n*dmy4
C
C dx1dz4,dz4dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen3= (1.0D0+xe*dlambdax1)*xe*dmuz4 + ze*dlambdax1*(1.0D0+ze*dmuz4)
     1          + ye**2*dlambdax1*dmuz4
         dnz4= (v*vixen3 + dvz4*dux1 - u*dqz4*dpx1 - dvx1*duz4)/v**2 - 2.0D0*dvz4*n/v

         imphell(3*a-2,3*d)= imphell(3*a-2,3*d) + m*dnz4 + n*dmz4
         imphell(3*d,3*a-2)= imphell(3*d,3*a-2) + m*dnz4 + n*dmz4
C
C dy1dz4,dz4dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= (1.0D0+ye*dlambday1)*ye*dmuz4 + ze*dlambday1*(1.0D0+ze*dmuz4)
     1          + xe**2*dlambday1*dmuz4
         dnz4= (v*vixen3 + dvz4*duy1 - u*dqz4*dpy1 - dvy1*duz4)/v**2 - 2.0D0*dvz4*n/v

         imphell(3*a-1,3*d)= imphell(3*a-1,3*d) + m*dnz4 + n*dmz4
         imphell(3*d,3*a-1)= imphell(3*d,3*a-1) + m*dnz4 + n*dmz4
C
C dy1dx4,dx4dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen3= (1.0D0+ye*dlambday1)*ye*dmux4 + xe*dlambday1*(1.0D0+xe*dmux4)
     1          + ze**2*dlambday1*dmux4
         dnx4= (v*vixen3 + dvx4*duy1 - u*dqx4*dpy1 - dvy1*dux4)/v**2 - 2.0D0*dvx4*n/v

         imphell(3*a-1,3*d-2)= imphell(3*a-1,3*d-2) + m*dnx4 + n*dmx4
         imphell(3*d-2,3*a-1)= imphell(3*d-2,3*a-1) + m*dnx4 + n*dmx4
C
C dz1dx4,dx4dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= (1.0D0+ze*dlambdaz1)*ze*dmux4 + xe*dlambdaz1*(1.0D0+xe*dmux4)
     1          + ye**2*dlambdaz1*dmux4
         dnx4= (v*vixen3 + dvx4*duz1 - u*dqx4*dpz1 - dvz1*dux4)/v**2 - 2.0D0*dvx4*n/v

         imphell(3*a,3*d-2)= imphell(3*a,3*d-2) + m*dnx4 + n*dmx4
         imphell(3*d-2,3*a)= imphell(3*d-2,3*a) + m*dnx4 + n*dmx4
C
C dz1dy4,dy4dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen3= (1.0D0+ze*dlambdaz1)*ze*dmuy4 + ye*dlambdaz1*(1.0D0+ye*dmuy4)
     1          + xe**2*dlambdaz1*dmuy4
         dny4= (v*vixen3 + dvy4*duz1 - u*dqy4*dpz1 - dvz1*duy4)/v**2 - 2.0D0*dvy4*n/v

         imphell(3*a,3*d-1)= imphell(3*a,3*d-1) + m*dny4 + n*dmy4
         imphell(3*d-1,3*a)= imphell(3*d-1,3*a) + m*dny4 + n*dmy4
C
C dx4dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+xe*dmux4)**2 + (ye**2 + ze**2)*dmux4**2 - dqx4**2)
         dnx4= -u*vixen3/v**2 - 2.0D0*dvx4*n/v

         imphell(3*d-2,3*d-2)= imphell(3*d-2,3*d-2) + m*dnx4 + n*dmx4
C
C dy4dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+ye*dmuy4)**2 + (xe**2 + ze**2)*dmuy4**2 - dqy4**2)
         dny4= -u*vixen3/v**2 - 2.0D0*dvy4*n/v

         imphell(3*d-1,3*d-1)= imphell(3*d-1,3*d-1) + m*dny4 + n*dmy4
C
C dz4dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+ze*dmuz4)**2 + (ye**2 + xe**2)*dmuz4**2 - dqz4**2)
         dnz4= -u*vixen3/v**2 - 2.0D0*dvz4*n/v

         imphell(3*d,3*d)= imphell(3*d,3*d) + m*dnz4 + n*dmz4
C
C dx4dy4,dy4dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+xe*dmux4)*xe*dmuy4 + ye*dmux4*(1.0D0+ye*dmuy4)
     1          + ze**2*dmux4*dmuy4 - dqx4*dqy4)
         dny4= (dvy4*dux4 - u*vixen3 - duy4*dvx4)/v**2 - 2.0D0*dvy4*n/v

         imphell(3*d-2,3*d-1)= imphell(3*d-2,3*d-1) + m*dny4 + n*dmy4
         imphell(3*d-1,3*d-2)= imphell(3*d-1,3*d-2) + m*dny4 + n*dmy4
C
C dx4dz4,dz4dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+xe*dmux4)*xe*dmuz4 + ze*dmux4*(1.0D0+ze*dmuz4)
     1          + ye**2*dmux4*dmuz4 - dqx4*dqz4)
         dnz4= (dvz4*dux4 - u*vixen3 - duz4*dvx4)/v**2 - 2.0D0*dvz4*n/v

         imphell(3*d-2,3*d)= imphell(3*d-2,3*d) + m*dnz4 + n*dmz4
         imphell(3*d,3*d-2)= imphell(3*d,3*d-2) + m*dnz4 + n*dmz4
C
C dy4dz4,dz4dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen3= (vixen1/vixen2)*((1.0D0+ye*dmuy4)*ye*dmuz4 + ze*dmuy4*(1.0D0+ze*dmuz4)
     1          + xe**2*dmuy4*dmuz4 - dqy4*dqz4)
         dnz4= (dvz4*duy4 - u*vixen3 - duz4*dvy4)/v**2 - 2.0D0*dvz4*n/v

         imphell(3*d-1,3*d)= imphell(3*d-1,3*d) + m*dnz4 + n*dmz4
         imphell(3*d,3*d-1)= imphell(3*d,3*d-1) + m*dnz4 + n*dmz4
C
C dx1dx2,dx2dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(x(b)-x(c))**2/rbc4
         vixen5= xp*(vixen4*xe-dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dlambdax2-1.0D0-lambda)
     1           + yp*ye*vixen4 + ye**2*dlambdax1*dlambdax2 + zp*ze*vixen4 
     2           + ze**2*dlambdax1*dlambdax2 - dpx1*dpx2
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqx2 
         vixen6= xpp*(vixen4*xe-dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dmux2-1.0D0-mu) +ypp*ye*vixen4
     1           + ye**2*dlambdax1*dmux2 + zpp*ze*vixen4 + ze**2*dlambdax1*dmux2
         dnx2= (v*vixen6 + dvx2*dux1 - u*vixen3 - dux2*dvx1)/v**2 - 2.0D0*dvx2*n/v

         imphell(3*a-2,3*b-2)= imphell(3*a-2,3*b-2)+ m*dnx2 + n*dmx2
         imphell(3*b-2,3*a-2)= imphell(3*b-2,3*a-2)+ m*dnx2 + n*dmx2
C
C dy1dy2,dy2dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(y(b)-y(c))**2/rbc4
         vixen5= yp*(vixen4*ye-dlambday1) + (1.0D0+ye*dlambday1)*(ye*dlambday2-1.0D0-lambda)
     1           + xp*xe*vixen4 + xe**2*dlambday1*dlambday2 + zp*ze*vixen4
     2           + ze**2*dlambday1*dlambday2 - dpy1*dpy2
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqy2
         vixen6= ypp*(vixen4*ye-dlambday1) + (1.0D0+ye*dlambday1)*(ye*dmuy2-1.0D0-mu) +xpp*xe*vixen4
     1           + xe**2*dlambday1*dmuy2 + zpp*ze*vixen4 + ze**2*dlambday1*dmuy2
         dny2= (v*vixen6 + dvy2*duy1 - u*vixen3 - duy2*dvy1)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*a-1,3*b-1)= imphell(3*a-1,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*a-1)= imphell(3*b-1,3*a-1)+ m*dny2 + n*dmy2
C
C dz1dz2,dz2dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(z(b)-z(c))**2/rbc4
         vixen5= zp*(vixen4*ze-dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dlambdaz2-1.0D0-lambda)
     1           + xp*xe*vixen4 + xe**2*dlambdaz1*dlambdaz2 + yp*ye*vixen4
     2           + ye**2*dlambdaz1*dlambdaz2 - dpz1*dpz2
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqz2
         vixen6= zpp*(vixen4*ze-dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dmuz2-1.0D0-mu) +xpp*xe*vixen4
     1           + xe**2*dlambdaz1*dmuz2 + ypp*ye*vixen4 + ye**2*dlambdaz1*dmuz2
         dnz2= (v*vixen6 + dvz2*duz1 - u*vixen3 - duz2*dvz1)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*a,3*b)= imphell(3*a,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*a)= imphell(3*b,3*a)+ m*dnz2 + n*dmz2
C
C dx1dy2,dy2dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambday2 + yp*(ye*vixen4-dlambdax1)
     1          + ye*dlambdax1*(ye*dlambday2-1.0D0-lambda) + zp*ze*vixen4 
     2          + ze**2*dlambdax1*dlambday2 - dpx1*dpy2
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqy2
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuy2 + ypp*(ye*vixen4-dlambdax1)
     1          + ye*dlambdax1*(ye*dmuy2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dlambdax1*dmuy2
         dny2= (v*vixen6 + dvy2*dux1 - u*vixen3 - duy2*dvx1)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*a-2,3*b-1)= imphell(3*a-2,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*a-2)= imphell(3*b-1,3*a-2)+ m*dny2 + n*dmy2
C
C dx1dz2,dz2dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambdaz2 + zp*(ze*vixen4-dlambdax1)
     1          + ze*dlambdax1*(ze*dlambdaz2-1.0D0-lambda) + yp*ye*vixen4
     2          + ye**2*dlambdax1*dlambdaz2 - dpx1*dpz2
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqz2
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuz2 + zpp*(ze*vixen4-dlambdax1)
     1          + ze*dlambdax1*(ze*dmuz2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dlambdax1*dmuz2
         dnz2= (v*vixen6 + dvz2*dux1 - u*vixen3 - duz2*dvx1)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*a-2,3*b)= imphell(3*a-2,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*a-2)= imphell(3*b,3*a-2)+ m*dnz2 + n*dmz2
C
C dy1dz2,dz2dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdaz2 + zp*(ze*vixen4-dlambday1)
     1          + ze*dlambday1*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dlambday1*dlambdaz2 - dpy1*dpz2
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqz2
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmuz2 + zpp*(ze*vixen4-dlambday1)
     1          + ze*dlambday1*(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dlambday1*dmuz2
         dnz2= (v*vixen6 + dvz2*duy1 - u*vixen3 - duz2*dvy1)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*a-1,3*b)= imphell(3*a-1,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*a-1)= imphell(3*b,3*a-1)+ m*dnz2 + n*dmz2
C
C dy1dx2,dx2dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdax2 + xp*(xe*vixen4-dlambday1)
     1          + xe*dlambday1*(xe*dlambdax2-1.0D0-lambda) + zp*ze*vixen4
     2          + ze**2*dlambday1*dlambdax2 - dpy1*dpx2
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqx2
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmux2 + xpp*(xe*vixen4-dlambday1)
     1          + xe*dlambday1*(xe*dmux2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dlambday1*dmux2
         dnx2= (v*vixen6 + dvx2*duy1 - u*vixen3 - dux2*dvy1)/v**2 - 2.0D0*dvx2*n/v

         imphell(3*a-1,3*b-2)= imphell(3*a-1,3*b-2)+ m*dnx2 + n*dmx2
         imphell(3*b-2,3*a-1)= imphell(3*b-2,3*a-1)+ m*dnx2 + n*dmx2
C
C dz1dx2,dx2dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambdax2 + xp*(xe*vixen4-dlambdaz1)
     1          + xe*dlambdaz1*(xe*dlambdax2-1.0D0-lambda) + yp*ye*vixen4
     2          + ye**2*dlambdaz1*dlambdax2 - dpz1*dpx2
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqx2
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmux2 + xpp*(xe*vixen4-dlambdaz1)
     1          + xe*dlambdaz1*(xe*dmux2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dlambdaz1*dmux2
         dnx2= (v*vixen6 + dvx2*duz1 - u*vixen3 - dux2*dvz1)/v**2 - 2.0D0*dvx2*n/v

         imphell(3*a,3*b-2)= imphell(3*a,3*b-2)+ m*dnx2 + n*dmx2
         imphell(3*b-2,3*a)= imphell(3*b-2,3*a)+ m*dnx2 + n*dmx2
C
C dz1dy2,dy2dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambday2 + yp*(ye*vixen4-dlambdaz1)
     1          + ye*dlambdaz1*(ye*dlambday2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dlambdaz1*dlambday2 - dpz1*dpy2
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqy2
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmuy2 + ypp*(ye*vixen4-dlambdaz1)
     1          + ye*dlambdaz1*(ye*dmuy2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dlambdaz1*dmuy2
         dny2= (v*vixen6 + dvy2*duz1 - u*vixen3 - duy2*dvz1)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*a,3*b-1)= imphell(3*a,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*a)= imphell(3*b-1,3*a)+ m*dny2 + n*dmy2
C
C dx1dx3,dx3dx1
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 2.0D0*(x(b)-x(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= xp*(xe*vixen4+dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dlambdax3+lambda)
     1          + (yp*ye+zp*ze)*vixen4 + (ye**2+ze**2)*dlambdax1*dlambdax3 - dpx1*dpx3
         vixen3= (vixen2/vixen1)*vixen5 + dpx1*dqx3
         vixen6= xpp*(xe*vixen4+dlambdax1) + (1.0D0+xe*dlambdax1)*(xe*dmux3+mu)
     1          + (ypp*ye+zpp*ze)*vixen4 + (ye**2+ze**2)*dlambdax1*dmux3
         dnx3= (v*vixen6 + dvx3*dux1 - u*vixen3 - dux3*dvx1)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*a-2,3*c-2)= imphell(3*a-2,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*a-2)= imphell(3*c-2,3*a-2)+ m*dnx3 + n*dmx3
C
C dy1dy3,dy3dy1
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 2.0D0*(y(b)-y(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= yp*(ye*vixen4+dlambday1) + (1.0D0+ye*dlambday1)*(ye*dlambday3+lambda)
     1          + (xp*xe+zp*ze)*vixen4 + (xe**2+ze**2)*dlambday1*dlambday3 - dpy1*dpy3
         vixen3= (vixen2/vixen1)*vixen5 + dpy1*dqy3
         vixen6= ypp*(ye*vixen4+dlambday1) + (1.0D0+ye*dlambday1)*(ye*dmuy3+mu)
     1          + (xpp*xe+zpp*ze)*vixen4 + (xe**2+ze**2)*dlambday1*dmuy3
         dny3= (v*vixen6 + dvy3*duy1 - u*vixen3 - duy3*dvy1)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*a-1,3*c-1)= imphell(3*a-1,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*a-1)= imphell(3*c-1,3*a-1)+ m*dny3 + n*dmy3
C
C dz1dz3,dz3dz1
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 2.0D0*(z(b)-z(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= zp*(ze*vixen4+dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dlambdaz3+lambda)
     1          + (xp*xe+yp*ye)*vixen4 + (xe**2+ye**2)*dlambdaz1*dlambdaz3 - dpz1*dpz3
         vixen3= (vixen2/vixen1)*vixen5 + dpz1*dqz3
         vixen6= zpp*(ze*vixen4+dlambdaz1) + (1.0D0+ze*dlambdaz1)*(ze*dmuz3+mu)
     1          + (xpp*xe+ypp*ye)*vixen4 + (xe**2+ye**2)*dlambdaz1*dmuz3
         dnz3= (v*vixen6 + dvz3*duz1 - u*vixen3 - duz3*dvz1)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*a,3*c)= imphell(3*a,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*a)= imphell(3*c,3*a)+ m*dnz3 + n*dmz3
C
C dx1dy3
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambday3 
     1          + yp*(ye*vixen4+dlambdax1) + ye*dlambdax1*(ye*dlambday3+lambda) 
     2          + zp*ze*vixen4 + ze**2*dlambdax1*dlambday3 - dpx1*dpy3
         vixen3=(vixen2/vixen1)*vixen5 + dpx1*dqy3
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuy3 + ypp*(ye*vixen4+dlambdax1)
     1          + ye*dlambdax1*(ye*dmuy3+mu) + zpp*ze*vixen4 + ze**2*dlambdax1*dmuy3
         dny3= (v*vixen6 + dvy3*dux1 - u*vixen3 - duy3*dvx1)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*a-2,3*c-1)= imphell(3*a-2,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*a-2)= imphell(3*c-1,3*a-2)+ m*dny3 + n*dmy3
C
C dx1dz3
C
         n=(v*dux1 - u*dvx1)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dlambdaz3
     1          + zp*(ze*vixen4+dlambdax1) + ze*dlambdax1*(ze*dlambdaz3+lambda)
     2          + yp*ye*vixen4 + ye**2*dlambdax1*dlambdaz3 - dpx1*dpz3
         vixen3=(vixen2/vixen1)*vixen5 + dpx1*dqz3
         vixen6= xpp*xe*vixen4 + (1.0D0+xe*dlambdax1)*xe*dmuz3 + zpp*(ze*vixen4+dlambdax1)
     1          + ze*dlambdax1*(ze*dmuz3+mu) + ypp*ye*vixen4 + ye**2*dlambdax1*dmuz3
         dnz3= (v*vixen6 + dvz3*dux1 - u*vixen3 - duz3*dvx1)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*a-2,3*c)= imphell(3*a-2,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*a-2)= imphell(3*c,3*a-2)+ m*dnz3 + n*dmz3
C
C dy1dz3
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdaz3
     1          + zp*(ze*vixen4+dlambday1) + ze*dlambday1*(ze*dlambdaz3+lambda)
     2          + xp*xe*vixen4 + xe**2*dlambday1*dlambdaz3 - dpy1*dpz3
         vixen3=(vixen2/vixen1)*vixen5 + dpy1*dqz3
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmuz3 + zpp*(ze*vixen4+dlambday1)
     1          + ze*dlambday1*(ze*dmuz3+mu) + xpp*xe*vixen4 + xe**2*dlambday1*dmuz3
         dnz3= (v*vixen6 + dvz3*duy1 - u*vixen3 - duz3*dvy1)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*a-1,3*c)= imphell(3*a-1,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*a-1)= imphell(3*c,3*a-1)+ m*dnz3 + n*dmz3
C
C dy1dx3
C
         n=(v*duy1 - u*dvy1)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= yp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dlambdax3
     1          + xp*(xe*vixen4+dlambday1) + xe*dlambday1*(xe*dlambdax3+lambda)
     2          + zp*ze*vixen4 + ze**2*dlambday1*dlambdax3 - dpy1*dpx3
         vixen3=(vixen2/vixen1)*vixen5 + dpy1*dqx3
         vixen6= ypp*ye*vixen4 + (1.0D0+ye*dlambday1)*ye*dmux3 + xpp*(xe*vixen4+dlambday1)
     1          + xe*dlambday1*(xe*dmux3+mu) + zpp*ze*vixen4 + ze**2*dlambday1*dmux3
         dnx3= (v*vixen6 + dvx3*duy1 - u*vixen3 - dux3*dvy1)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*a-1,3*c-2)= imphell(3*a-1,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*a-1)= imphell(3*c-2,3*a-1)+ m*dnx3 + n*dmx3
C
C dz1dx3
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambdax3
     1          + xp*(xe*vixen4+dlambdaz1) + xe*dlambdaz1*(xe*dlambdax3+lambda)
     2          + yp*ye*vixen4 + ye**2*dlambdaz1*dlambdax3 - dpz1*dpx3
         vixen3=(vixen2/vixen1)*vixen5 + dpz1*dqx3
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmux3 + xpp*(xe*vixen4+dlambdaz1)
     1          + xe*dlambdaz1*(xe*dmux3+mu) + ypp*ye*vixen4 + ye**2*dlambdaz1*dmux3
         dnx3= (v*vixen6 + dvx3*duz1 - u*vixen3 - dux3*dvz1)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*a,3*c-2)= imphell(3*a,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*a)= imphell(3*c-2,3*a)+ m*dnx3 + n*dmx3
C
C dz1dy3
C
         n=(v*duz1 - u*dvz1)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dlambday3
     1          + yp*(ye*vixen4+dlambdaz1) + ye*dlambdaz1*(ye*dlambday3+lambda)
     2          + xp*xe*vixen4 + xe**2*dlambdaz1*dlambday3 - dpz1*dpy3
         vixen3=(vixen2/vixen1)*vixen5 + dpz1*dqy3
         vixen6= zpp*ze*vixen4 + (1.0D0+ze*dlambdaz1)*ze*dmuy3 + ypp*(ye*vixen4+dlambdaz1)
     1          + ye*dlambdaz1*(ye*dmuy3+mu) + xpp*xe*vixen4 + xe**2*dlambdaz1*dmuy3
         dny3= (v*vixen6 + dvy3*duz1 - u*vixen3 - duy3*dvz1)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*a,3*c-1)= imphell(3*a,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*a)= imphell(3*c-1,3*a)+ m*dny3 + n*dmy3
C
C dx4dx2,dx2dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(x(b)-x(c))**2/rbc4
         vixen5= xpp*(xe*vixen4-dmux4) + (1.0D0+xe*dmux4)
     1          *(xe*dmux2-1.0D0-mu) + ypp*ye*vixen4
     2          + ye**2*dmux4*dmux2 + zpp*ze*vixen4
     3          + ze**2*dmux4*dmux2 - dqx4*dqx2
         vixen3=(vixen1/vixen2)*vixen5 + dpx2*dqx4
         vixen6= xp*(xe*vixen4-dmux4) + (1.0D0+xe*dmux4)
     1          *(xe*dlambdax2-1.0D0-lambda) + yp*ye*vixen4
     2          + ye**2*dmux4*dlambdax2 + zp*ze*vixen4
     3          + ze**2*dmux4*dlambdax2
         dnx2= (v*vixen6 + dvx2*dux4 - u*vixen3 - dux2*dvx4)/v**2
     1        - 2.0D0*dvx2*(v*dux4 - u*dvx4)/v**3

         imphell(3*d-2,3*b-2)= imphell(3*d-2,3*b-2)+ m*dnx2 + n*dmx2
         imphell(3*b-2,3*d-2)= imphell(3*b-2,3*d-2)+ m*dnx2 + n*dmx2
C
C dy4dy2,dy2dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(y(b)-y(c))**2/rbc4
         vixen5= ypp*(ye*vixen4-dmuy4) + (1.0D0+ye*dmuy4)
     1          *(ye*dmuy2-1.0D0-mu) + xpp*xe*vixen4
     2          + xe**2*dmuy4*dmuy2 + zpp*ze*vixen4
     3          + ze**2*dmuy4*dmuy2 - dqy4*dqy2
         vixen3=(vixen1/vixen2)*vixen5 + dpy2*dqy4
         vixen6= yp*(ye*vixen4-dmuy4) + (1.0D0+ye*dmuy4)
     1          *(ye*dlambday2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dmuy4*dlambday2 + zp*ze*vixen4
     3          + ze**2*dmuy4*dlambday2
         dny2= (v*vixen6 + dvy2*duy4 - u*vixen3 - duy2*dvy4)/v**2
     1        - 2.0D0*dvy2*(v*duy4 - u*dvy4)/v**3

         imphell(3*d-1,3*b-1)= imphell(3*d-1,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*d-1)= imphell(3*b-1,3*d-1)+ m*dny2 + n*dmy2
C
C dz4dz2,dz2dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 1.0D0/rbc2 - 2.0D0*(z(b)-z(c))**2/rbc4
         vixen5= zpp*(ze*vixen4-dmuz4) + (1.0D0+ze*dmuz4)
     1          *(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4
     2          + xe**2*dmuz4*dmuz2 + ypp*ye*vixen4
     3          + ye**2*dmuz4*dmuz2 - dqz4*dqz2
         vixen3=(vixen1/vixen2)*vixen5 + dpz2*dqz4
         vixen6= zp*(ze*vixen4-dmuz4) + (1.0D0+ze*dmuz4)
     1          *(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen4
     2          + xe**2*dmuz4*dlambdaz2 + yp*ye*vixen4
     3          + ye**2*dmuz4*dlambdaz2
         dnz2= (v*vixen6 + dvz2*duz4 - u*vixen3 - duz2*dvz4)/v**2
     1        - 2.0D0*dvz2*(v*duz4 - u*dvz4)/v**3

         imphell(3*d,3*b)= imphell(3*d,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*d)= imphell(3*b,3*d)+ m*dnz2 + n*dmz2
C
C dx4dy2,dy2dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuy2 + ypp*(ye*vixen4-dmux4)
     1          + ye*dmux4*(ye*dmuy2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dmux4*dmuy2 - dqx4*dqy2
         vixen3=(vixen1/vixen2)*vixen5 + dpy2*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambday2 + yp*(ye*vixen4-dmux4)
     1          + ye*dmux4*(ye*dlambday2-1.0D0-lambda) + zp*ze*vixen4 + ze**2*dmux4*dlambday2
         dny2= (v*vixen6 + dvy2*dux4 - u*vixen3 - duy2*dvx4)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*d-2,3*b-1)= imphell(3*d-2,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*d-2)= imphell(3*b-1,3*d-2)+ m*dny2 + n*dmy2
C
C dx4dz2,dz2dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= -2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuz2 + zpp*(ze*vixen4-dmux4)
     1          + ze*dmux4*(ze*dmuz2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dmux4*dmuz2 - dqx4*dqz2
         vixen3=(vixen1/vixen2)*vixen5 + dpz2*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambdaz2 + zp*(ze*vixen4-dmux4)
     1          + ze*dmux4*(ze*dlambdaz2-1.0D0-lambda) + yp*ye*vixen4 + ye**2*dmux4*dlambdaz2
         dnz2= (v*vixen6 + dvz2*dux4 - u*vixen3 - duz2*dvx4)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*d-2,3*b)= imphell(3*d-2,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*d-2)= imphell(3*b,3*d-2)+ m*dnz2 + n*dmz2
C
C dy4dz2,dz2dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmuz2 + zpp*(ze*vixen4-dmuy4)
     1          + ze*dmuy4*(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dmuy4*dmuz2 - dqy4*dqz2
         vixen3=(vixen1/vixen2)*vixen5 + dpz2*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdaz2 + zp*(ze*vixen4-dmuy4)
     1          + ze*dmuy4*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen4 + xe**2*dmuy4*dlambdaz2
         dnz2= (v*vixen6 + dvz2*duy4 - u*vixen3 - duz2*dvy4)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*d-1,3*b)= imphell(3*d-1,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*d-1)= imphell(3*b,3*d-1)+ m*dnz2 + n*dmz2
C
C dy4dx2,dx2dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= -2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmux2 + xpp*(xe*vixen4-dmuy4)
     1          + xe*dmuy4*(xe*dmux2-1.0D0-mu) + zpp*ze*vixen4 + ze**2*dmuy4*dmux2 - dqy4*dqx2
         vixen3=(vixen1/vixen2)*vixen5 + dpx2*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdax2 + xp*(xe*vixen4-dmuy4)
     1          + xe*dmuy4*(xe*dlambdax2-1.0D0-lambda) + zp*ze*vixen4 + ze**2*dmuy4*dlambdax2
         dnx2= (v*vixen6 + dvx2*duy4 - u*vixen3 - dux2*dvy4)/v**2 - 2.0D0*dvx2*n/v

         imphell(3*d-1,3*b-2)= imphell(3*d-1,3*b-2)+ m*dnx2 + n*dmx2
         imphell(3*b-2,3*d-1)= imphell(3*b-2,3*d-1)+ m*dnx2 + n*dmx2
C
C dz4dx2,dx2dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmux2 + xpp*(xe*vixen4-dmuz4)
     1          + xe*dmuz4*(xe*dmux2-1.0D0-mu) + ypp*ye*vixen4 + ye**2*dmuz4*dmux2 - dqz4*dqx2
         vixen3=(vixen1/vixen2)*vixen5 + dpx2*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambdax2 + xp*(xe*vixen4-dmuz4)
     1          + xe*dmuz4*(xe*dlambdax2-1.0D0-lambda) + yp*ye*vixen4 + ye**2*dmuz4*dlambdax2
         dnx2= (v*vixen6 + dvx2*duz4 - u*vixen3 - dux2*dvz4)/v**2 - 2.0D0*dvx2*n/v

         imphell(3*d,3*b-2)= imphell(3*d,3*b-2)+ m*dnx2 + n*dmx2
         imphell(3*b-2,3*d)= imphell(3*b-2,3*d)+ m*dnx2 + n*dmx2
C
C dz4dy2,dy2dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= -2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmuy2 + ypp*(ye*vixen4-dmuz4)
     1          + ye*dmuz4*(ye*dmuy2-1.0D0-mu) + xpp*xe*vixen4 + xe**2*dmuz4*dmuy2 - dqz4*dqy2
         vixen3=(vixen1/vixen2)*vixen5 + dpy2*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambday2 + yp*(ye*vixen4-dmuz4)
     1          + ye*dmuz4*(ye*dlambday2-1.0D0-lambda) + xp*xe*vixen4 + xe**2*dmuz4*dlambday2
         dny2= (v*vixen6 + dvy2*duz4 - u*vixen3 - duy2*dvz4)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*d,3*b-1)= imphell(3*d,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*d)= imphell(3*b-1,3*d)+ m*dny2 + n*dmy2
C
C dx4dx3,dx3dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 2.0D0*(x(b)-x(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= xpp*(xe*vixen4+dmux4) + (1.0D0+xe*dmux4)*(xe*dmux3+mu) + ypp*ye*vixen4
     1          + ye**2*dmux4*dmux3 + zpp*ze*vixen4 + ze**2*dmux4*dmux3 - dqx4*dqx3
         vixen3=(vixen1/vixen2)*vixen5 + dpx3*dqx4
         vixen6= xp*(xe*vixen4+dmux4) + (1.0D0+xe*dmux4)*(xe*dlambdax3+lambda) + vixen4*
     1          (yp*ye+zp*ze) + (ye**2+ze**2)*dmux4*dlambdax3
         dnx3= (v*vixen6 + dvx3*dux4 - u*vixen3 - dux3*dvx4)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*d-2,3*c-2)= imphell(3*d-2,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*d-2)= imphell(3*c-2,3*d-2)+ m*dnx3 + n*dmx3
C
C dy4dy3,dy3dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 2.0D0*(y(b)-y(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= ypp*(ye*vixen4+dmuy4) + (1.0D0+ye*dmuy4)*(ye*dmuy3+mu) + xpp*xe*vixen4
     1          + xe**2*dmuy4*dmuy3 + zpp*ze*vixen4 + ze**2*dmuy4*dmuy3 - dqy4*dqy3
         vixen3=(vixen1/vixen2)*vixen5 + dpy3*dqy4
         vixen6= yp*(ye*vixen4+dmuy4) + (1.0D0+ye*dmuy4)*(ye*dlambday3+lambda) + vixen4*
     1          (xp*xe+zp*ze) + (xe**2+ze**2)*dmuy4*dlambday3
         dny3= (v*vixen6 + dvy3*duy4 - u*vixen3 - duy3*dvy4)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*d-1,3*c-1)= imphell(3*d-1,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*d-1)= imphell(3*c-1,3*d-1)+ m*dny3 + n*dmy3
C
C dz4dz3,dz3dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 2.0D0*(z(b)-z(c))**2/rbc4 - 1.0D0/rbc2
         vixen5= zpp*(ze*vixen4+dmuz4) + (1.0D0+ze*dmuz4)*(ze*dmuz3+mu) + xpp*xe*vixen4
     1          + xe**2*dmuz4*dmuz3 + ypp*ye*vixen4 + ye**2*dmuz4*dmuz3 - dqz4*dqz3
         vixen3=(vixen1/vixen2)*vixen5 + dpz3*dqz4
         vixen6= zp*(ze*vixen4+dmuz4) + (1.0D0+ze*dmuz4)*(ze*dlambdaz3+lambda) + vixen4*
     1          (xp*xe+yp*ye) + (xe**2+ye**2)*dmuz4*dlambdaz3
         dnz3= (v*vixen6 + dvz3*duz4 - u*vixen3 - duz3*dvz4)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*d,3*c)= imphell(3*d,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*d)= imphell(3*c,3*d)+ m*dnz3 + n*dmz3
C
C dx4dy3,dy3dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(y(b)-y(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuy3 + ypp*(ye*vixen4+dmux4)
     1          + ye*dmux4*(ye*dmuy3+mu) + zpp*ze*vixen4 + ze**2*dmux4*dmuy3 - dqx4*dqy3
         vixen3=(vixen1/vixen2)*vixen5 + dpy3*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambday3 + yp*(ye*vixen4+dmux4)
     1          + ye*dmux4*(ye*dlambday3+lambda) + zp*ze*vixen4 + ze**2*dmux4*dlambday3
         dny3= (v*vixen6 + dvy3*dux4 - u*vixen3 - duy3*dvx4)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*d-2,3*c-1)= imphell(3*d-2,3*c-1) + m*dny3 + n*dmy3
         imphell(3*c-1,3*d-2)= imphell(3*c-1,3*d-2) + m*dny3 + n*dmy3
C
C dx4dz3,dz3dx4
C
         n=(v*dux4 - u*dvx4)/v**2
         vixen4= 2.0D0*(x(b)-x(c))*(z(b)-z(c))/rbc4
         vixen5= xpp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dmuz3 + zpp*(ze*vixen4+dmux4)
     1          + ze*dmux4*(ze*dmuz3+mu) + ypp*ye*vixen4 + ye**2*dmux4*dmuz3 - dqx4*dqz3
         vixen3=(vixen1/vixen2)*vixen5 + dpz3*dqx4
         vixen6= xp*xe*vixen4 + (1.0D0+xe*dmux4)*xe*dlambdaz3 + zp*(ze*vixen4+dmux4)
     1          + ze*dmux4*(ze*dlambdaz3+lambda) + yp*ye*vixen4 + ye**2*dmux4*dlambdaz3
         dnz3= (v*vixen6 + dvz3*dux4 - u*vixen3 - duz3*dvx4)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*d-2,3*c)= imphell(3*d-2,3*c) + m*dnz3 + n*dmz3
         imphell(3*c,3*d-2)= imphell(3*c,3*d-2) + m*dnz3 + n*dmz3
C
C dy4dz3,dz3dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(z(b)-z(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmuz3 + zpp*(ze*vixen4+dmuy4)
     1          + ze*dmuy4*(ze*dmuz3+mu) + xpp*xe*vixen4 + xe**2*dmuy4*dmuz3 - dqy4*dqz3
         vixen3=(vixen1/vixen2)*vixen5 + dpz3*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdaz3 + zp*(ze*vixen4+dmuy4)
     1          + ze*dmuy4*(ze*dlambdaz3+lambda) + xp*xe*vixen4 + xe**2*dmuy4*dlambdaz3
         dnz3= (v*vixen6 + dvz3*duy4 - u*vixen3 - duz3*dvy4)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*d-1,3*c)= imphell(3*d-1,3*c) + m*dnz3 + n*dmz3
         imphell(3*c,3*d-1)= imphell(3*c,3*d-1) + m*dnz3 + n*dmz3
C
C dy4dx3,dx3dy4
C
         n=(v*duy4 - u*dvy4)/v**2
         vixen4= 2.0D0*(y(b)-y(c))*(x(b)-x(c))/rbc4
         vixen5= ypp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dmux3 + xpp*(xe*vixen4+dmuy4)
     1          + xe*dmuy4*(xe*dmux3+mu) + zpp*ze*vixen4 + ze**2*dmuy4*dmux3 - dqy4*dqx3
         vixen3=(vixen1/vixen2)*vixen5 + dpx3*dqy4
         vixen6= yp*ye*vixen4 + (1.0D0+ye*dmuy4)*ye*dlambdax3 + xp*(xe*vixen4+dmuy4)
     1          + xe*dmuy4*(xe*dlambdax3+lambda) + zp*ze*vixen4 + ze**2*dmuy4*dlambdax3
         dnx3= (v*vixen6 + dvx3*duy4 - u*vixen3 - dux3*dvy4)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*d-1,3*c-2)= imphell(3*d-1,3*c-2) + m*dnx3 + n*dmx3
         imphell(3*c-2,3*d-1)= imphell(3*c-2,3*d-1) + m*dnx3 + n*dmx3
C
C dz4dx3,dx3dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(x(b)-x(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmux3 + xpp*(xe*vixen4+dmuz4)
     1          + xe*dmuz4*(xe*dmux3+mu) + ypp*ye*vixen4 + ye**2*dmuz4*dmux3 - dqz4*dqx3
         vixen3=(vixen1/vixen2)*vixen5 + dpx3*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambdax3 + xp*(xe*vixen4+dmuz4)
     1          + xe*dmuz4*(xe*dlambdax3+lambda) + yp*ye*vixen4 + ye**2*dmuz4*dlambdax3
         dnx3= (v*vixen6 + dvx3*duz4 - u*vixen3 - dux3*dvz4)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*d,3*c-2)= imphell(3*d,3*c-2) + m*dnx3 + n*dmx3
         imphell(3*c-2,3*d)= imphell(3*c-2,3*d) + m*dnx3 + n*dmx3
C
C dz4dy3,dy3dz4
C
         n=(v*duz4 - u*dvz4)/v**2
         vixen4= 2.0D0*(z(b)-z(c))*(y(b)-y(c))/rbc4
         vixen5= zpp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dmuy3 + ypp*(ye*vixen4+dmuz4)
     1          + ye*dmuz4*(ye*dmuy3+mu) + xpp*xe*vixen4 + xe**2*dmuz4*dmuy3 - dqz4*dqy3
         vixen3=(vixen1/vixen2)*vixen5 + dpy3*dqz4
         vixen6= zp*ze*vixen4 + (1.0D0+ze*dmuz4)*ze*dlambday3 + yp*(ye*vixen4+dmuz4)
     1          + ye*dmuz4*(ye*dlambday3+lambda) + xp*xe*vixen4 + xe**2*dmuz4*dlambday3
         dny3= (v*vixen6 + dvy3*duz4 - u*vixen3 - duy3*dvz4)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*d,3*c-1)= imphell(3*d,3*c-1) + m*dny3 + n*dmy3
         imphell(3*c-1,3*d)= imphell(3*c-1,3*d) + m*dny3 + n*dmy3
C
C dx3dx3
C
         n=(v*dux3 - u*dvx3)/v**2
         vixen4= 8.0D0*(xe**2-rbc2)*mucrap/rbc6 - 4.0D0*xe*(x(b)-x(d))/rbc4
         vixen5= 8.0D0*(xe**2-rbc2)*lambdacrap/rbc6 - 4.0D0*xe*(x(b)-x(a))/rbc4
         vixen7= xpp*(2.0D0*dmux3+xe*vixen4) + (mu+xe*dmux3)**2 + vixen4*(ypp*ye+zpp*ze)
     1          + dmux3**2*(ye**2+ze**2) - dqx3**2
         vixen8= xp*(2.0D0*dlambdax3+xe*vixen5) + (lambda+xe*dlambdax3)**2 + vixen5*(yp*ye+zp*ze)
     1          + dlambdax3**2*(ye**2+ze**2) - dpx3**2
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + 2.0D0*dpx3*dqx3
         vixen6= xp*(2.0D0*dmux3+xe*vixen4) + 2.0D0*(mu+xe*dmux3)*(lambda+xe*dlambdax3) 
     1          + xpp*(2.0D0*dlambdax3+xe*vixen5) + vixen4*(yp*ye+zp*ze) 
     2          + vixen5*(ypp*ye+zpp*ze) + 2.0D0*dlambdax3*dmux3*(ye**2+ze**2)
         dnx3= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*c-2,3*c-2)= imphell(3*c-2,3*c-2)+ m*dnx3 + n*dmx3
C
C dy3dy3
C
         n=(v*duy3 - u*dvy3)/v**2
         vixen4= 8.0D0*(xe**2-rbc2)*mucrap/rbc6 - 4.0D0*xe*(x(b)-x(d))/rbc4
         vixen5= 8.0D0*(xe**2-rbc2)*lambdacrap/rbc6 - 4.0D0*xe*(x(b)-x(a))/rbc4
         vixen7= ypp*(2.0D0*dmuy3+ye*vixen4) + (mu+ye*dmuy3)**2 + vixen4*(xpp*xe+zpp*ze)
     1          + dmuy3**2*(xe**2+ze**2) - dqy3**2
         vixen8= yp*(2.0D0*dlambday3+ye*vixen5) + (lambda+ye*dlambday3)**2 + vixen5*(xp*xe+zp*ze)
     1          + dlambday3**2*(xe**2+ze**2) - dpy3**2
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + 2.0D0*dpy3*dqy3
         vixen6= yp*(2.0D0*dmuy3+ye*vixen4) + 2.0D0*(mu+ye*dmuy3)*(lambda+ye*dlambday3)
     1          + ypp*(2.0D0*dlambday3+ye*vixen5) + vixen4*(xp*xe+zp*ze)
     2          + vixen5*(xpp*xe+zpp*ze) + 2.0D0*dlambday3*dmuy3*(xe**2+ze**2)
         dny3= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*c-1,3*c-1)= imphell(3*c-1,3*c-1)+ m*dny3 + n*dmy3
C
C dz3dz3
C
         n=(v*duz3 - u*dvz3)/v**2
         vixen4= 8.0D0*(xe**2-rbc2)*mucrap/rbc6 - 4.0D0*xe*(x(b)-x(d))/rbc4
         vixen5= 8.0D0*(xe**2-rbc2)*lambdacrap/rbc6 - 4.0D0*xe*(x(b)-x(a))/rbc4
         vixen7= zpp*(2.0D0*dmuz3+ze*vixen4) + (mu+ze*dmuz3)**2 + vixen4*(xpp*xe+ypp*ye)
     1          + dmuz3**2*(xe**2+ye**2) - dqz3**2
         vixen8= zp*(2.0D0*dlambdaz3+ze*vixen5) + (lambda+ze*dlambdaz3)**2 + vixen5*(xp*xe+yp*ye)
     1          + dlambdaz3**2*(xe**2+ye**2) - dpz3**2
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + 2.0D0*dpz3*dqz3
         vixen6= zp*(2.0D0*dmuz3+ze*vixen4) + 2.0D0*(mu+ze*dmuz3)*(lambda+ze*dlambdaz3)
     1          + zpp*(2.0D0*dlambdaz3+ze*vixen5) + vixen4*(xp*xe+yp*ye)
     2          + vixen5*(xpp*xe+ypp*ye) + 2.0D0*dlambdaz3*dmuz3*(xe**2+ye**2)
         dnz3= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*c,3*c)= imphell(3*c,3*c)+ m*dnz3 + n*dmz3
C
C dx3dy3,dy3dx3
C
         n=(v*dux3 - u*dvx3)/v**2
         vixen4= 8.0D0*mucrap*xe*ye/rbc6 - 2.0D0*(ye*(x(b)-x(d))+xe*(y(b)-y(d)))/rbc4
         vixen5= 8.0D0*mucrap*xe*ye/rbc6 - 2.0D0*(ye*(x(b)-x(a))+xe*(y(b)-y(a)))/rbc4
         vixen7= xpp*(dmuy3+xe*vixen4) + (mu+xe*dmux3)*xe*dmuy3 + ypp*(ye*vixen4+dmux3)
     1          + ye*dmux3*(mu+ye*dmuy3) + zpp*ze*vixen4 + ze**2*dmux3*dmuy3 - dqx3*dqy3
         vixen8= xp*(dlambday3+xe*vixen5) + (lambda+xe*dlambdax3)*xe*dlambday3 
     1          + yp*(ye*vixen5 + dlambdax3) + ye*dlambdax3*(lambda+ye*dlambday3) 
     2          + zp*ze*vixen5 + ze**2*dlambdax3*dlambday3 - dpx3*dpy3
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + dpx3*dqy3 + dpy3*dqx3
         vixen6= xp*(dmuy3+xe*vixen4) + (mu+xe*dmux3)*xe*dlambday3 + xpp*(dlambday3+xe*vixen5)
     1          + (lambda+xe*dlambdax3)*xe*dmuy3 + yp*(ye*vixen4+dmux3) 
     2          + ye*dmux3*(ye*dlambday3+lambda) + ypp*(ye*vixen5 + dlambdax3) 
     3          + ye*dlambdax3*(mu+ye*dmuy3) + zp*ze*vixen4 + zpp*ze*vixen5
     4          + ze**2*(dlambdax3*dmuy3+dmux3*dlambday3)
         dny3= (v*vixen6 + dvy3*dux3 - u*vixen3 - duy3*dvx3)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*c-2,3*c-1)= imphell(3*c-2,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*c-2)= imphell(3*c-1,3*c-2)+ m*dny3 + n*dmy3
C
C dx3dz3,dz3dx3
C
         n=(v*dux3 - u*dvx3)/v**2
         vixen4= 8.0D0*mucrap*xe*ze/rbc6 - 2.0D0*(ze*(x(b)-x(d))+xe*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*mucrap*xe*ze/rbc6 - 2.0D0*(ze*(x(b)-x(a))+xe*(z(b)-z(a)))/rbc4
         vixen7= xpp*(dmuz3+xe*vixen4) + (mu+xe*dmux3)*xe*dmuz3 + zpp*(ze*vixen4+dmux3)
     1          + ze*dmux3*(mu+ze*dmuz3) + ypp*ye*vixen4 + ye**2*dmux3*dmuz3 - dqx3*dqz3
         vixen8= xp*(dlambdaz3+xe*vixen5) + (lambda+xe*dlambdax3)*xe*dlambdaz3
     1          + zp*(ze*vixen5 + dlambdax3) + ze*dlambdax3*(lambda+ze*dlambdaz3)
     2          + yp*ye*vixen5 + ye**2*dlambdax3*dlambdaz3 - dpx3*dpz3
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + dpx3*dqz3 + dpz3*dqx3
         vixen6= xp*(dmuz3+xe*vixen4) + (mu+xe*dmux3)*xe*dlambdaz3 + xpp*(dlambdaz3+xe*vixen5)
     1          + (lambda+xe*dlambdax3)*xe*dmuz3 + zp*(ze*vixen4+dmux3)
     2          + ze*dmux3*(ze*dlambdaz3+lambda) + zpp*(ze*vixen5 + dlambdax3)
     3          + ze*dlambdax3*(mu+ze*dmuz3) + yp*ye*vixen4 + ypp*ye*vixen5
     4          + ye**2*(dlambdax3*dmuz3+dmux3*dlambdaz3)
         dnz3= (v*vixen6 + dvz3*dux3 - u*vixen3 - duz3*dvx3)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*c-2,3*c)= imphell(3*c-2,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*c-2)= imphell(3*c,3*c-2)+ m*dnz3 + n*dmz3
C
C dy3dz3,dz3dy3
C
         n=(v*duy3 - u*dvy3)/v**2
         vixen4= 8.0D0*mucrap*ye*ze/rbc6 - 2.0D0*(ze*(y(b)-y(d))+ye*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*mucrap*ye*ze/rbc6 - 2.0D0*(ze*(y(b)-y(a))+ye*(z(b)-z(a)))/rbc4
         vixen7= ypp*(dmuz3+ye*vixen4) + (mu+ye*dmuy3)*ye*dmuz3 + zpp*(ze*vixen4+dmuy3)
     1          + ze*dmuy3*(mu+ze*dmuz3) + xpp*xe*vixen4 + xe**2*dmuy3*dmuz3 - dqy3*dqz3
         vixen8= yp*(dlambdaz3+ye*vixen5) + (lambda+ye*dlambday3)*ye*dlambdaz3
     1          + zp*(ze*vixen5 + dlambday3) + ze*dlambday3*(lambda+ze*dlambdaz3)
     2          + xp*xe*vixen5 + xe**2*dlambday3*dlambdaz3 - dpy3*dpz3
         vixen3= (vixen1/vixen2)*vixen7 + (vixen2/vixen1)*vixen8 + dpy3*dqz3 + dpz3*dqy3
         vixen6= yp*(dmuz3+ye*vixen4) + (mu+ye*dmuy3)*ye*dlambdaz3 + ypp*(dlambdaz3+ye*vixen5)
     1          + (lambda+ye*dlambday3)*ye*dmuz3 + zp*(ze*vixen4+dmuy3)
     2          + ze*dmuy3*(ze*dlambdaz3+lambda) + zpp*(ze*vixen5 + dlambday3)
     3          + ze*dlambday3*(mu+ze*dmuz3) + xp*xe*vixen4 + xpp*xe*vixen5
     4          + xe**2*(dlambday3*dmuz3+dmuy3*dlambdaz3)
         dnz3= (v*vixen6 + dvz3*duy3 - u*vixen3 - duz3*dvy3)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*c-1,3*c)= imphell(3*c-1,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*c-1)= imphell(3*c,3*c-1)+ m*dnz3 + n*dmz3
C
C dx2dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= -2.0D0/rbc2 - (4.0D0*(x(b)-x(c))*(x(c)+x(d)+2.0D0*x(b))-
     1          2.0D0*mucrap)/rbc4 + 8.0D0*(x(b)-x(c))**2*mucrap/rbc6
         vixen5= -2.0D0/rbc2 - (4.0D0*(x(b)-x(c))*(x(c)+x(a)+2.0D0*x(b))-
     1    2.0D0*lambdacrap)/rbc4 + 8.0D0*(x(b)-x(c))**2*lambdacrap/rbc6
         vixen7= xp*(xe*vixen5-2.0D0*dlambdax2) + (xe*dlambdax2
     1    -1.0D0-lambda)**2 + (yp*ye+zp*ze)*vixen5 - dpx2**2
     2    + (ye**2+ze**2)*dlambdax2**2
         vixen8= xpp*(xe*vixen4-2.0D0*dmux2) + (xe*dmux2-1.0D0-mu)
     1    **2 + (ypp*ye+zpp*ze)*vixen4 - dqx2**2
     2    + (ye**2+ze**2)*dmux2**2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8+2.0D0*dpx2*dqx2
         vixen6= xp*(xe*vixen4-2.0D0*dmux2) + 2.0D0*(xe*dmux2-1.0D0-mu)
     1    *(xe*dlambdax2-1.0D0-lambda) + xpp*(xe*vixen5-2.0D0*
     2    dlambdax2) + vixen4*(yp*ye+zp*ze) + vixen5*
     3    (ypp*ye+zpp*ze) + 2.0D0*dmux2*dlambdax2*(ye
     4    **2+ze**2)
         dnx2= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvx2*n/v

         imphell(3*b-2,3*b-2)= imphell(3*b-2,3*b-2)+ m*dnx2 + n*dmx2
C
C dy2dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= -2.0D0/rbc2 - (4.0D0*(y(b)-y(c))*(y(c)+y(d)+2.0D0*y(b))-
     1          2.0D0*mucrap)/rbc4 + 8.0D0*(y(b)-y(c))**2*mucrap/rbc6
         vixen5= -2.0D0/rbc2 - (4.0D0*(y(b)-y(c))*(y(c)+y(a)+2.0D0*y(b))-
     1    2.0D0*lambdacrap)/rbc4 + 8.0D0*(y(b)-y(c))**2*lambdacrap/rbc6
         vixen7= yp*(ye*vixen5-2.0D0*dlambday2) + (ye*dlambday2
     1    -1.0D0-lambda)**2 + (xp*xe+zp*ze)*vixen5 - dpy2**2
     2    + (xe**2+ze**2)*dlambday2**2
         vixen8= ypp*(ye*vixen4-2.0D0*dmuy2) + (ye*dmuy2-1.0D0-mu)
     1    **2 + (xpp*xe+zpp*ze)*vixen4 - dqy2**2
     2    + (xe**2+ze**2)*dmuy2**2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8+2.0D0*dpy2*dqy2
         vixen6= yp*(ye*vixen4-2.0D0*dmuy2) + 2.0D0*(ye*dmuy2-1.0D0-mu)
     1    *(ye*dlambday2-1.0D0-lambda) + ypp*(ye*vixen5-2.0D0*
     2    dlambday2) + vixen4*(xp*xe+zp*ze) + vixen5*
     3    (xpp*xe+zpp*ze) + 2.0D0*dmuy2*dlambday2*(xe
     4    **2+ze**2)
         dny2= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*b-1,3*b-1)= imphell(3*b-1,3*b-1)+ m*dny2 + n*dmy2
C
C dz2dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= -2.0D0/rbc2 - (4.0D0*(z(b)-z(c))*(z(c)+z(d)+2.0D0*z(b))-
     1          2.0D0*mucrap)/rbc4 + 8.0D0*(z(b)-z(c))**2*mucrap/rbc6
         vixen5= -2.0D0/rbc2 - (4.0D0*(z(b)-z(c))*(z(c)+z(a)+2.0D0*z(b))-
     1    2.0D0*lambdacrap)/rbc4 + 8.0D0*(z(b)-z(c))**2*lambdacrap/rbc6
         vixen7= zp*(ze*vixen5-2.0D0*dlambdaz2) + (ze*dlambdaz2
     1    -1.0D0-lambda)**2 + (xp*xe+yp*ye)*vixen5 - dpz2**2
     2    + (xe**2+ye**2)*dlambdaz2**2
         vixen8= zpp*(ze*vixen4-2.0D0*dmuz2) + (ze*dmuz2-1.0D0-mu)
     1    **2 + (xpp*xe+ypp*ye)*vixen4 - dqz2**2
     2    + (xe**2+ye**2)*dmuz2**2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8+2.0D0*dpz2*dqz2
         vixen6= zp*(ze*vixen4-2*dmuz2) + 2.0D0*(ze*dmuz2-1.0D0-mu)
     1    *(ze*dlambdaz2-1.0D0-lambda) + zpp*(ze*vixen5-2.0D0*
     2    dlambdaz2) + vixen4*(xp*xe+yp*ye) + vixen5*
     3    (xpp*xe+ypp*ye) + 2.0D0*dmuz2*dlambdaz2*(xe
     4    **2+ye**2)
         dnz2= (v*vixen6 - u*vixen3)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*b,3*b)= imphell(3*b,3*b)+ m*dnz2 + n*dmz2
C
C dx2dy2,dy2dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= (-2.0D0*(y(b)-y(c))*(x(c)+x(d)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(y(c)+y(d)-2.0D0*y(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(y(b)-y(c))*mucrap
     2          /rbc6
         vixen5= (-2.0D0*(y(b)-y(c))*(x(c)+x(a)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(y(c)+y(a)-2.0D0*y(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(y(b)-y(c))
     2          *lambdacrap/rbc6
         vixen7= xp*(xe*vixen5-dlambday2) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambday2 + yp*(ye
     1          *vixen5-dlambdax2) + ye*dlambdax2*(ye*dlambday2-1.0D0-lambda) + zp*ze*vixen5
     2          + ze**2*dlambdax2*dlambday2 - dpx2*dpy2
         vixen8= xpp*(xe*vixen4-dmuy2) + (xe*dmux2-1.0D0-mu)*xe*dmuy2 + ypp*(ye
     1          *vixen4-dmux2) + ye*dmux2*(ye*dmuy2-1.0D0-mu) + zpp*ze*vixen4
     2          + ze**2*dmux2*dmuy2 - dqx2*dqy2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqy2 + dqx2*dpy2
         vixen6= xp*(xe*vixen4-dmuy2) + (xe*dmux2-1.0D0-mu)*xe*dlambday2 + xpp*(xe*vixen5-dlambday2)
     1          + (xe*dlambdax2-1.0D0-lambda)*xe*dmuy2 + yp*(ye*vixen4-dmux2) + ye*dmux2
     2          *(ye*dlambday2-1.0D0-lambda) + ypp*(ye*vixen5-dlambdax2) + ye*dlambdax2
     3          *(ye*dmuy2-1.0D0-mu) + zp*ze*vixen4 + zpp*ze*vixen5 + ze**2*(dmux2*dlambday2
     4          +dmuy2*dlambdax2)
         dny2= (v*vixen6 + dvy2*dux2 - u*vixen3 - dvx2*duy2)/v**2 - 2.0D0*dvy2*n/v

         imphell(3*b-2,3*b-1)= imphell(3*b-2,3*b-1)+ m*dny2 + n*dmy2
         imphell(3*b-1,3*b-2)= imphell(3*b-1,3*b-2)+ m*dny2 + n*dmy2
C
C dx2dz2,dz2dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= (-2.0D0*(z(b)-z(c))*(x(c)+x(d)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(z(c)+z(d)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(z(b)-z(c))*mucrap
     2          /rbc6
         vixen5= (-2.0D0*(z(b)-z(c))*(x(c)+x(a)-2.0D0*x(b))-2.0D0*(x(b)-x(c))*(z(c)+z(a)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(x(b)-x(c))*(z(b)-z(c))
     2          *lambdacrap/rbc6
         vixen7= xp*(xe*vixen5-dlambdaz2) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambdaz2 + zp*(ze
     1          *vixen5-dlambdax2) + ze*dlambdax2*(ze*dlambdaz2-1.0D0-lambda) + yp*ye*vixen5
     2          + ye**2*dlambdax2*dlambdaz2 - dpx2*dpz2
         vixen8= xpp*(xe*vixen4-dmuz2) + (xe*dmux2-1.0D0-mu)*xe*dmuz2 + zpp*(ze
     1          *vixen4-dmux2) + ze*dmux2*(ze*dmuz2-1.0D0-mu) + ypp*ye*vixen4
     2          + ye**2*dmux2*dmuz2 - dqx2*dqz2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqz2 + dqx2*dpz2
         vixen6= xp*(xe*vixen4-dmuz2) + (xe*dmux2-1.0D0-mu)*xe*dlambdaz2 + xpp*(xe*vixen5-dlambdaz2)
     1          + (xe*dlambdax2-1.0D0-lambda)*xe*dmuz2 + zp*(ze*vixen4-dmux2) + ze*dmux2
     2          *(ze*dlambdaz2-1.0D0-lambda) + zpp*(ze*vixen5-dlambdax2) + ze*dlambdax2
     3          *(ze*dmuz2-1.0D0-mu) + yp*ye*vixen4 + ypp*ye*vixen5 + ye**2*(dmux2*dlambdaz2
     4          +dmuz2*dlambdax2)
         dnz2= (v*vixen6 + dvz2*dux2 - u*vixen3 - dvx2*duz2)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*b-2,3*b)= imphell(3*b-2,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*b-2)= imphell(3*b,3*b-2)+ m*dnz2 + n*dmz2
C
C dy2dz2,dz2dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= (-2.0D0*(z(b)-z(c))*(y(c)+y(d)-2.0D0*y(b))-2.0D0*(y(b)-y(c))*(z(c)+z(d)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(y(b)-y(c))*(z(b)-z(c))*mucrap
     2          /rbc6
         vixen5= (-2.0D0*(z(b)-z(c))*(y(c)+y(a)-2.0D0*y(b))-2.0D0*(y(b)-y(c))*(z(c)+z(a)-2.0D0*z(b)))/rbc4 
     1          + 8.0D0*(y(b)-y(c))*(z(b)-z(c))
     2          *lambdacrap/rbc6
         vixen7= yp*(ye*vixen5-dlambdaz2) + (ye*dlambday2-1.0D0-lambda)*ye*dlambdaz2 + zp*(ze
     1          *vixen5-dlambday2) + ze*dlambday2*(ze*dlambdaz2-1.0D0-lambda) + xp*xe*vixen5
     2          + xe**2*dlambday2*dlambdaz2 - dpy2*dpz2
         vixen8= ypp*(ye*vixen4-dmuz2) + (ye*dmuy2-1.0D0-mu)*ye*dmuz2 + zpp*(ze
     1          *vixen4-dmuy2) + ze*dmuy2*(ze*dmuz2-1.0D0-mu) + xpp*xe*vixen4
     2          + xe**2*dmuy2*dmuz2 - dqy2*dqz2
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqz2 + dqy2*dpz2
         vixen6= yp*(ye*vixen4-dmuz2) + (ye*dmuy2-1.0D0-mu)*ye*dlambdaz2 + ypp*(ye*vixen5-dlambdaz2)
     1          + (ye*dlambday2-1.0D0-lambda)*ye*dmuz2 + zp*(ze*vixen4-dmuy2) + ze*dmuy2
     2          *(ze*dlambdaz2-1.0D0-lambda) + zpp*(ze*vixen5-dlambday2) + ze*dlambday2
     3          *(ze*dmuz2-1.0D0-mu) + xp*xe*vixen4 + xpp*xe*vixen5 + xe**2*(dmuy2*dlambdaz2
     4          +dmuz2*dlambday2)
         dnz2= (v*vixen6 + dvz2*duy2 - u*vixen3 - dvy2*duz2)/v**2 - 2.0D0*dvz2*n/v

         imphell(3*b-1,3*b)= imphell(3*b-1,3*b)+ m*dnz2 + n*dmz2
         imphell(3*b,3*b-1)= imphell(3*b,3*b-1)+ m*dnz2 + n*dmz2
C
C dx2dx3,dx3dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= 1.0D0/rbc2 - (2.0D0*xe*(x(c)+x(d)-2.0D0*x(b))+2.0D0*((x(b)-x(c))*(x(b)-x(d))-mucrap))/rbc4 
     1          - 8.0D0*(x(b)-x(c))**2*mucrap/rbc6
         vixen5= 1.0D0/rbc2 - (2.0D0*xe*(x(c)+x(a)-2.0D0*x(b))+2.0D0*((x(b)-x(c))*(x(b)-x(a))-lambdacrap))/rbc4 
     1          - 8.0D0*(x(b)-x(c))**2*lambdacrap/rbc6
         vixen7= xp*(xe*vixen5+dlambdax2-dlambdax3) + (xe*dlambdax2-1.0D0-lambda)*(xe*dlambdax3+lambda)
     1          + (yp*ye+zp*ze)*vixen5 + (ye**2+ze**2)*dlambdax2*dlambdax3 - dpx2*dpx3
         vixen8= xpp*(xe*vixen4+dmux2-dmux3) + (xe*dmux2-1.0D0-mu)*(xe*dmux3+mu)
     1          + (ypp*ye+zpp*ze)*vixen4 + (ye**2+ze**2)*dmux2*dmux3 - dqx2*dqx3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqx3 + dqx2*dpx3
         vixen6= xp*(xe*vixen4+dmux2-dmux3) + (xe*dmux2-1.0D0-mu)*(xe*dlambdax3+lambda)
     1          + xpp*(xe*vixen5+dlambdax2-dlambdax3) + (xe*dlambdax2-1.0D0-lambda)*(xe*dmux3+mu)
     2          + vixen4*(yp*ye+zp*ze) + vixen5*(ypp*ye+zpp*ze) 
     3          + (dmux2*dlambdax3+dmux3*dlambdax2)*(ye**2+ze**2)
         dnx3= (v*vixen6 + dvx3*dux2 - u*vixen3 - dux3*dvx2)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*b-2,3*c-2)= imphell(3*b-2,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*b-2)= imphell(3*c-2,3*b-2)+ m*dnx3 + n*dmx3
C
C dy2dy3,dy3dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= 1.0D0/rbc2 - (2.0D0*ye*(y(c)+y(d)-2.0D0*y(b))+2.0D0*((y(b)-y(c))*(y(b)-y(d))-mucrap))/rbc4 
     1          - 8.0D0*(y(b)-y(c))**2*mucrap/rbc6
         vixen5= 1.0D0/rbc2 - (2.0D0*ye*(y(c)+y(a)-2.0D0*y(b))+2.0D0*((y(b)-y(c))*(y(b)-y(a))-lambdacrap))/rbc4 
     1          - 8.0D0*(y(b)-y(c))**2*lambdacrap/rbc6
         vixen7= yp*(ye*vixen5+dlambday2-dlambday3) + (ye*dlambday2-1.0D0-lambda)*(ye*dlambday3+lambda)
     1          + (xp*xe+zp*ze)*vixen5 + (xe**2+ze**2)*dlambday2*dlambday3 - dpy2*dpy3
         vixen8= ypp*(ye*vixen4+dmuy2-dmuy3) + (ye*dmuy2-1.0D0-mu)*(ye*dmuy3+mu)
     1          + (xpp*xe+zpp*ze)*vixen4 + (xe**2+ze**2)*dmuy2*dmuy3 - dqy2*dqy3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqy3 + dqy2*dpy3
         vixen6= yp*(ye*vixen4+dmuy2-dmuy3) + (ye*dmuy2-1.0D0-mu)*(ye*dlambday3+lambda)
     1          + ypp*(ye*vixen5+dlambday2-dlambday3) + (ye*dlambday2-1.0D0-lambda)*(ye*dmuy3+mu)
     2          + vixen4*(xp*xe+zp*ze) + vixen5*(xpp*xe+zpp*ze) 
     3          + (dmuy2*dlambday3+dmuy3*dlambday2)*(xe**2+ze**2)
         dny3= (v*vixen6 + dvy3*duy2 - u*vixen3 - duy3*dvy2)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*b-1,3*c-1)= imphell(3*b-1,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*b-1)= imphell(3*c-1,3*b-1)+ m*dny3 + n*dmy3
C
C dz2dz3,dz3dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= 1.0D0/rbc2 - (2.0D0*ze*(z(c)+z(d)-2.0D0*z(b))+2.0D0*((z(b)-z(c))*(z(b)-z(d))-mucrap))/rbc4 
     1          - 8.0D0*(z(b)-z(c))**2*mucrap/rbc6
         vixen5= 1.0D0/rbc2 - (2.0D0*ze*(z(c)+z(a)-2.0D0*z(b))+2.0D0*((z(b)-z(c))*(z(b)-z(a))-lambdacrap))/rbc4 
     1          - 8.0D0*(z(b)-z(c))**2*lambdacrap/rbc6
         vixen7= zp*(ze*vixen5+dlambdaz2-dlambdaz3) + (ze*dlambdaz2-1.0D0-lambda)*(ze*dlambdaz3+lambda)
     1          + (xp*xe+yp*ye)*vixen5 + (xe**2+ye**2)*dlambdaz2*dlambdaz3 - dpz2*dpz3
         vixen8= zpp*(ze*vixen4+dmuz2-dmuz3) + (ze*dmuz2-1.0D0-mu)*(ze*dmuz3+mu)
     1          + (xpp*xe+ypp*ye)*vixen4 + (xe**2+ye**2)*dmuz2*dmuz3 - dqz2*dqz3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpz2*dqz3 + dqz2*dpz3
         vixen6= zp*(ze*vixen4+dmuz2-dmuz3) + (ze*dmuz2-1.0D0-mu)*(ze*dlambdaz3+lambda)
     1          + zpp*(ze*vixen5+dlambdaz2-dlambdaz3) + (ze*dlambdaz2-1.0D0-lambda)*(ze*dmuz3+mu)
     2          + vixen4*(xp*xe+yp*ye) + vixen5*(xpp*xe+ypp*ye) 
     3          + (dmuz2*dlambdaz3+dmuz3*dlambdaz2)*(xe**2+ye**2)
         dnz3= (v*vixen6 + dvz3*duz2 - u*vixen3 - duz3*dvz2)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*b,3*c)= imphell(3*b,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*b)= imphell(3*c,3*b)+ m*dnz3 + n*dmz3
C
C dx2dy3,dy3dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= 8.0D0*(x(b)-x(c))*ye*mucrap/rbc6 - 2.0D0*((x(c)+x(d)-2.0D0*x(b))*ye
     1          +(x(b)-x(c))*(y(b)-y(d)))/rbc4
         vixen5= 8.0D0*(x(b)-x(c))*ye*lambdacrap/rbc6 - 2.0D0*((x(c)+x(a)-2.0D0*x(b))*ye
     1          +(x(b)-x(c))*(y(b)-y(a)))/rbc4
         vixen7= xp*(xe*vixen5-dlambday3) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambday3 + 
     1          yp*(ye*vixen5+dlambdax2) + ye*dlambdax2*(ye*dlambday3+lambda) 
     2          + zp*ze*vixen5 + ze**2*dlambdax2*dlambday3 - dpx2*dpy3
         vixen8= xpp*(xe*vixen4-dmuy3) + (xe*dmux2-1.0D0-mu)*xe*dmuy3 + 
     1          ypp*(ye*vixen4+dmux2) + ye*dmux2*(ye*dmuy3+mu) + 
     2          zpp*ze*vixen4 + ze**2*dmux2*dmuy3 - dqx2*dqy3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqy3 + dqx2*dpy3
         vixen6= xp*(xe*vixen4-dmuy3) + (xe*dmux2-1.0D0-mu)*xe*dlambday3
     1          + xpp*(xe*vixen5-dlambday3) + (xe*dlambdax2-1.0D0-lambda)*xe*dmuy3
     2          + yp*(ye*vixen4+dmux2) + ye*dmux2*(ye*dlambday3+lambda)
     3          + ypp*(ye*vixen5+dlambdax2) + ye*dlambdax2*(ye*dmuy3+mu)
     4          + zp*ze*vixen4 + ze**2*(dmux2*dlambday3+dlambdax2*dmuy3) + zpp*ze*vixen5
         dny3= (v*vixen6 + dvy3*dux2 - u*vixen3 - duy3*dvx2)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*b-2,3*c-1)= imphell(3*b-2,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*b-2)= imphell(3*c-1,3*b-2)+ m*dny3 + n*dmy3
C
C dx2dz3,dz3dx2
C
         n=(v*dux2 - u*dvx2)/v**2
         vixen4= 8.0D0*(x(b)-x(c))*ze*mucrap/rbc6 - 2.0D0*((x(c)+x(d)-2.0D0*x(b))*ze
     1          +(x(b)-x(c))*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*(x(b)-x(c))*ze*lambdacrap/rbc6 - 2.0D0*((x(c)+x(a)-2.0D0*x(b))*ze
     1          +(x(b)-x(c))*(z(b)-z(a)))/rbc4
         vixen7= xp*(xe*vixen5-dlambdaz3) + (xe*dlambdax2-1.0D0-lambda)*xe*dlambdaz3 + 
     1          zp*(ze*vixen5+dlambdax2) + ze*dlambdax2*(ze*dlambdaz3+lambda) 
     2          + yp*ye*vixen5 + ye**2*dlambdax2*dlambdaz3 - dpx2*dpz3
         vixen8= xpp*(xe*vixen4-dmuz3) + (xe*dmux2-1.0D0-mu)*xe*dmuz3 + 
     1          zpp*(ze*vixen4+dmux2) + ze*dmux2*(ze*dmuz3+mu) + 
     2          ypp*ye*vixen4 + ye**2*dmux2*dmuz3 - dqx2*dqz3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpx2*dqz3 + dqx2*dpz3
         vixen6= xp*(xe*vixen4-dmuz3) + (xe*dmux2-1.0D0-mu)*xe*dlambdaz3
     1          + xpp*(xe*vixen5-dlambdaz3) + (xe*dlambdax2-1.0D0-lambda)*xe*dmuz3
     2          + zp*(ze*vixen4+dmux2) + ze*dmux2*(ze*dlambdaz3+lambda)
     3          + zpp*(ze*vixen5+dlambdax2) + ze*dlambdax2*(ze*dmuz3+mu)
     4          + yp*ye*vixen4 + ye**2*(dmux2*dlambdaz3+dlambdax2*dmuz3) + ypp*ye*vixen5
         dnz3= (v*vixen6 + dvz3*dux2 - u*vixen3 - duz3*dvx2)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*b-2,3*c)= imphell(3*b-2,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*b-2)= imphell(3*c,3*b-2)+ m*dnz3 + n*dmz3
C
C dy2dz3,dz3dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= 8.0D0*(y(b)-y(c))*ze*mucrap/rbc6 - 2.0D0*((y(c)+y(d)-2.0D0*y(b))*ze
     1          +(y(b)-y(c))*(z(b)-z(d)))/rbc4
         vixen5= 8.0D0*(y(b)-y(c))*ze*lambdacrap/rbc6 - 2.0D0*((y(c)+y(a)-2.0D0*y(b))*ze
     1          +(y(b)-y(c))*(z(b)-z(a)))/rbc4
         vixen7= yp*(ye*vixen5-dlambdaz3) + (ye*dlambday2-1.0D0-lambda)*ye*dlambdaz3 +
     1          zp*(ze*vixen5+dlambday2) + ze*dlambday2*(ze*dlambdaz3+lambda)
     2          + xp*xe*vixen5 + xe**2*dlambday2*dlambdaz3 - dpy2*dpz3
         vixen8= ypp*(ye*vixen4-dmuz3) + (ye*dmuy2-1.0D0-mu)*ye*dmuz3 +
     1          zpp*(ze*vixen4+dmuy2) + ze*dmuy2*(ze*dmuz3+mu) +
     2          xpp*xe*vixen4 + xe**2*dmuy2*dmuz3 - dqy2*dqz3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqz3 + dqy2*dpz3
         vixen6= yp*(ye*vixen4-dmuz3) + (ye*dmuy2-1.0D0-mu)*ye*dlambdaz3
     1          + ypp*(ye*vixen5-dlambdaz3) + (ye*dlambday2-1.0D0-lambda)*ye*dmuz3
     2          + zp*(ze*vixen4+dmuy2) + ze*dmuy2*(ze*dlambdaz3+lambda)
     3          + zpp*(ze*vixen5+dlambday2) + ze*dlambday2*(ze*dmuz3+mu)
     4          + xp*xe*vixen4 + xe**2*(dmuy2*dlambdaz3+dlambday2*dmuz3) + xpp*xe*vixen5
         dnz3= (v*vixen6 + dvz3*duy2 - u*vixen3 - duz3*dvy2)/v**2 - 2.0D0*dvz3*n/v

         imphell(3*b-1,3*c)= imphell(3*b-1,3*c)+ m*dnz3 + n*dmz3
         imphell(3*c,3*b-1)= imphell(3*c,3*b-1)+ m*dnz3 + n*dmz3
C
C dy2dx3,dx3dy2
C
         n=(v*duy2 - u*dvy2)/v**2
         vixen4= 8.0D0*(y(b)-y(c))*xe*mucrap/rbc6 - 2.0D0*((y(c)+y(d)-2.0D0*y(b))*xe
     1          +(y(b)-y(c))*(x(b)-x(d)))/rbc4
         vixen5= 8.0D0*(y(b)-y(c))*xe*lambdacrap/rbc6 - 2.0D0*((y(c)+y(a)-2.0D0*y(b))*xe
     1          +(y(b)-y(c))*(x(b)-x(a)))/rbc4
         vixen7= yp*(ye*vixen5-dlambdax3) + (ye*dlambday2-1.0D0-lambda)*ye*dlambdax3 +
     1          xp*(xe*vixen5+dlambday2) + xe*dlambday2*(xe*dlambdax3+lambda)
     2          + zp*ze*vixen5 + ze**2*dlambday2*dlambdax3 - dpy2*dpx3
         vixen8= ypp*(ye*vixen4-dmux3) + (ye*dmuy2-1.0D0-mu)*ye*dmux3 +
     1          xpp*(xe*vixen4+dmuy2) + xe*dmuy2*(xe*dmux3+mu) +
     2          zpp*ze*vixen4 + ze**2*dmuy2*dmux3 - dqy2*dqx3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpy2*dqx3 + dqy2*dpx3
         vixen6= yp*(ye*vixen4-dmux3) + (ye*dmuy2-1.0D0-mu)*ye*dlambdax3
     1          + ypp*(ye*vixen5-dlambdax3) + (ye*dlambday2-1.0D0-lambda)*ye*dmux3
     2          + xp*(xe*vixen4+dmuy2) + xe*dmuy2*(xe*dlambdax3+lambda)
     3          + xpp*(xe*vixen5+dlambday2) + xe*dlambday2*(xe*dmux3+mu)
     4          + zp*ze*vixen4 + ze**2*(dmuy2*dlambdax3+dlambday2*dmux3) + zpp*ze*vixen5
         dnx3= (v*vixen6 + dvx3*duy2 - u*vixen3 - dux3*dvy2)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*b-1,3*c-2)= imphell(3*b-1,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*b-1)= imphell(3*c-2,3*b-1)+ m*dnx3 + n*dmx3
C
C dz2dx3,dx3dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= 8.0D0*(z(b)-z(c))*xe*mucrap/rbc6 - 2.0D0*((z(c)+z(d)-2.0D0*z(b))*xe
     1          +(z(b)-z(c))*(x(b)-x(d)))/rbc4
         vixen5= 8.0D0*(z(b)-z(c))*xe*lambdacrap/rbc6 - 2.0D0*((z(c)+z(a)-2.0D0*z(b))*xe
     1          +(z(b)-z(c))*(x(b)-x(a)))/rbc4
         vixen7= zp*(ze*vixen5-dlambdax3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dlambdax3 +
     1          xp*(xe*vixen5+dlambdaz2) + xe*dlambdaz2*(xe*dlambdax3+lambda)
     2          + yp*ye*vixen5 + ye**2*dlambdaz2*dlambdax3 - dpz2*dpx3
         vixen8= zpp*(ze*vixen4-dmux3) + (ze*dmuz2-1.0D0-mu)*ze*dmux3 +
     1          xpp*(xe*vixen4+dmuz2) + xe*dmuz2*(xe*dmux3+mu) +
     2          ypp*ye*vixen4 + ye**2*dmuz2*dmux3 - dqz2*dqx3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpz2*dqx3 + dqz2*dpx3
         vixen6= zp*(ze*vixen4-dmux3) + (ze*dmuz2-1.0D0-mu)*ze*dlambdax3
     1          + zpp*(ze*vixen5-dlambdax3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dmux3
     2          + xp*(xe*vixen4+dmuz2) + xe*dmuz2*(xe*dlambdax3+lambda)
     3          + xpp*(xe*vixen5+dlambdaz2) + xe*dlambdaz2*(xe*dmux3+mu)
     4          + yp*ye*vixen4 + ye**2*(dmuz2*dlambdax3+dlambdaz2*dmux3) + ypp*ye*vixen5
         dnx3= (v*vixen6 + dvx3*duz2 - u*vixen3 - dux3*dvz2)/v**2 - 2.0D0*dvx3*n/v

         imphell(3*b,3*c-2)= imphell(3*b,3*c-2)+ m*dnx3 + n*dmx3
         imphell(3*c-2,3*b)= imphell(3*c-2,3*b)+ m*dnx3 + n*dmx3
C
C dz2dy3,dy3dz2
C
         n=(v*duz2 - u*dvz2)/v**2
         vixen4= 8.0D0*(z(b)-z(c))*ye*mucrap/rbc6 - 2.0D0*((z(c)+z(d)-2.0D0*z(b))*ye
     1          +(z(b)-z(c))*(y(b)-y(d)))/rbc4
         vixen5= 8.0D0*(z(b)-z(c))*ye*lambdacrap/rbc6 - 2.0D0*((z(c)+z(a)-2.0D0*z(b))*ye
     1          +(z(b)-z(c))*(y(b)-y(a)))/rbc4
         vixen7= zp*(ze*vixen5-dlambday3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dlambday3 +
     1          yp*(ye*vixen5+dlambdaz2) + ye*dlambdaz2*(ye*dlambday3+lambda)
     2          + xp*xe*vixen5 + xe**2*dlambdaz2*dlambday3 - dpz2*dpy3
         vixen8= zpp*(ze*vixen4-dmuy3) + (ze*dmuz2-1.0D0-mu)*ze*dmuy3 +
     1          ypp*(ye*vixen4+dmuz2) + ye*dmuz2*(ye*dmuy3+mu) +
     2          xpp*xe*vixen4 + xe**2*dmuz2*dmuy3 - dqz2*dqy3
         vixen3= (vixen2/vixen1)*vixen7 + (vixen1/vixen2)*vixen8 + dpz2*dqy3 + dqz2*dpy3
         vixen6= zp*(ze*vixen4-dmuy3) + (ze*dmuz2-1.0D0-mu)*ze*dlambday3
     1          + zpp*(ze*vixen5-dlambday3) + (ze*dlambdaz2-1.0D0-lambda)*ze*dmuy3
     2          + yp*(ye*vixen4+dmuz2) + ye*dmuz2*(ye*dlambday3+lambda)
     3          + ypp*(ye*vixen5+dlambdaz2) + ye*dlambdaz2*(ye*dmuy3+mu)
     4          + xp*xe*vixen4 + xe**2*(dmuz2*dlambday3+dlambdaz2*dmuy3) + xpp*xe*vixen5
         dny3= (v*vixen6 + dvy3*duz2 - u*vixen3 - duy3*dvz2)/v**2 - 2.0D0*dvy3*n/v

         imphell(3*b,3*c-1)= imphell(3*b,3*c-1)+ m*dny3 + n*dmy3
         imphell(3*c-1,3*b)= imphell(3*c-1,3*b)+ m*dny3 + n*dmy3
      END DO

      DO i=1,3*atoms
         DO j=1,3*atoms
            hell(i,j)=bondhell(i,j)+anglehell(i,j)+torshell(i,j)+imphell(i,j)+qhell(i,j)+vdwhell(i,j)
         END DO
      END DO

      RETURN

      END


      SUBROUTINE hairyhell(u,v)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION      DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      COMMON /DU/           dux1,dux2,dux3,dux4,duy1,duy2,duy3,duy4,duz1,duz2,duz3,duz4
      DOUBLE PRECISION      DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      COMMON /DV/           dvx1,dvx2,dvx3,dvx4,dvy1,dvy2,dvy3,dvy4,dvz1,dvz2,dvz3,dvz4
      DOUBLE PRECISION      DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      COMMON /DM/           dmx1,dmx2,dmx3,dmx4,dmy1,dmy2,dmy3,dmy4,dmz1,dmz2,dmz3,dmz4
      DOUBLE PRECISION      U,V
      DOUBLE PRECISION      M,N
      COMMON /MN/           m,n

C      PRINT *,'In Hairyhell',i
      m=0.0D0
      dmx1=0.0D0
      dmx2=0.0D0
      dmx3=0.0D0
      dmx4=0.0D0
      dmy1=0.0D0
      dmy2=0.0D0
      dmy3=0.0D0
      dmy4=0.0D0
      dmz1=0.0D0
      dmz2=0.0D0
      dmz3=0.0D0
      dmz4=0.0D0

      IF (dn(i).EQ.2) THEN
C         PRINT *,'PN1=2'
         m= (dvn(i)/did(i))*4.0D0*COS(dphi(i))*COS(ddelta(i))
         dmx1= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= (dvn(i)/did(i))*4.0D0*COS(ddelta(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE IF (dn(i).EQ.3) THEN
C         PRINT *,'PN1=3'
         m= (dvn(i)/did(i))*(9.0D0-12.0D0*(SIN(dphi(i)))**2)*COS(ddelta(i))
         dmx1= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= (dvn(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE IF (dn(i).EQ.4) THEN
C         PRINT *,'PN1=4'
         m= (dvn(i)/did(i))*(32.0D0*(COS(dphi(i)))**3 - 16.0D0*COS(dphi(i)))*COS(ddelta(i))
         dmx1= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= (dvn(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta(i))*((v*duz4 - u*dvz4)/v**2)
      END IF

      IF (dn2(i).EQ.1) THEN
C         PRINT *,'PN2=1'
         m= m+(dvn2(i)/did(i))*COS(ddelta2(i))
C
C Do Nothing
C
      ELSE IF (dn2(i).EQ.2) THEN
C         PRINT *,'PN2=2'
         m= m+ (dvn2(i)/did(i))*4.0D0*COS(dphi(i))*COS(ddelta(i))
         dmx1= dmx1+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= dmx2+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= dmx3+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= dmx4+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= dmy1+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= dmy2+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= dmy3+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= dmy4+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= dmz1+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= dmz2+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= dmz3+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= dmz4+ (dvn2(i)/did(i))*4.0D0*COS(ddelta2(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE IF (dn2(i).EQ.3) THEN
C         PRINT *,'PN2=3'
         m= m+ (dvn2(i)/did(i))*(9.0D0-12.0D0*(SIN(dphi(i)))**2)*COS(ddelta2(i))
         dmx1= dmx1+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= dmx2+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= dmx3+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= dmx4+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= dmy1+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= dmy2+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= dmy3+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= dmy4+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= dmz1+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= dmz2+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= dmz3+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= dmz4+ (dvn2(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta2(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE IF (dn2(i).EQ.4) THEN
C         PRINT *,'PN2=4'
         m= m+ (dvn2(i)/did(i))*(32.0D0*(COS(dphi(i)))**3 - 16.0D0*COS(dphi(i)))*COS(ddelta2(i))
         dmx1= dmx1+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= dmx2+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= dmx3+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= dmx4+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= dmy1+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= dmy2+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= dmy3+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= dmy4+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= dmz1+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= dmz2+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= dmz3+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= dmz4+ (dvn2(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta2(i))*((v*duz4 - u*dvz4)/v**2)
      END IF

      IF (dn3(i).EQ.1) THEN
C         PRINT *,'PN3=1'
         m= m+(dvn3(i)/did(i))*COS(ddelta3(i))
C
C Do Nothing
C
      ELSE IF (dn3(i).EQ.2) THEN
C         PRINT *,'PN3=2'
         m= m+ (dvn3(i)/did(i))*4.0D0*COS(dphi(i))*COS(ddelta3(i))
         dmx1= dmx1+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= dmx2+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= dmx3+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= dmx4+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= dmy1+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= dmy2+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= dmy3+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= dmy4+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= dmz1+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= dmz2+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= dmz3+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= dmz4+ (dvn3(i)/did(i))*4.0D0*COS(ddelta3(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE IF (dn3(i).EQ.3) THEN
C         PRINT *,'PN3=3'
         m= m+ (dvn3(i)/did(i))*(9.0D0-12.0D0*(SIN(dphi(i)))**2)*COS(ddelta3(i))
         dmx1= dmx1+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= dmx2+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= dmx3+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= dmx4+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= dmy1+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= dmy2+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= dmy3+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= dmy4+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= dmz1+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= dmz2+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= dmz3+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= dmz4+ (dvn3(i)/did(i))*24.0D0*COS(dphi(i))*COS(ddelta3(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE IF (dn3(i).EQ.4) THEN
C         PRINT *,'PN3=4'
         m= m+ (dvn3(i)/did(i))*(32.0D0*(COS(dphi(i)))**3 - 16.0D0*COS(dphi(i)))*COS(ddelta3(i))
         dmx1= dmx1+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= dmx2+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= dmx3+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= dmx4+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= dmy1+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= dmy2+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= dmy3+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= dmy4+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= dmz1+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= dmz2+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= dmz3+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= dmz4+ (dvn3(i)/did(i))*(96.0D0*(COS(dphi(i))**2)-16.0D0)*COS(ddelta3(i))*((v*duz4 - u*dvz4)/v**2)
      END IF


      RETURN

      END


      SUBROUTINE hairyimphell(u,v)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION      DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      COMMON /DU/           dux1,dux2,dux3,dux4,duy1,duy2,duy3,duy4,duz1,duz2,duz3,duz4
      DOUBLE PRECISION      DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      COMMON /DV/           dvx1,dvx2,dvx3,dvx4,dvy1,dvy2,dvy3,dvy4,dvz1,dvz2,dvz3,dvz4
      DOUBLE PRECISION      DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      COMMON /DM/           dmx1,dmx2,dmx3,dmx4,dmy1,dmy2,dmy3,dmy4,dmz1,dmz2,dmz3,dmz4
      DOUBLE PRECISION      U,V
      DOUBLE PRECISION      M,N
      COMMON /MN/           m,n

C      PRINT *,'In Hairyimphell',i
      m=0.0D0
      dmx1=0.0D0
      dmx2=0.0D0
      dmx3=0.0D0
      dmx4=0.0D0
      dmy1=0.0D0
      dmy2=0.0D0
      dmy3=0.0D0
      dmy4=0.0D0
      dmz1=0.0D0
      dmz2=0.0D0
      dmz3=0.0D0
      dmz4=0.0D0

      IF (in1(i).EQ.2) THEN
C         PRINT *,'PN1=2'
         m= ivn(i)*4.0D0*COS(iphi(i))*COS(idelta(i))
         dmx1= ivn(i)*4.0D0*COS(idelta(i))*((v*dux1 - u*dvx1)/v**2)
         dmx2= ivn(i)*4.0D0*COS(idelta(i))*((v*dux2 - u*dvx2)/v**2)
         dmx3= ivn(i)*4.0D0*COS(idelta(i))*((v*dux3 - u*dvx3)/v**2)
         dmx4= ivn(i)*4.0D0*COS(idelta(i))*((v*dux4 - u*dvx4)/v**2)
         dmy1= ivn(i)*4.0D0*COS(idelta(i))*((v*duy1 - u*dvy1)/v**2)
         dmy2= ivn(i)*4.0D0*COS(idelta(i))*((v*duy2 - u*dvy2)/v**2)
         dmy3= ivn(i)*4.0D0*COS(idelta(i))*((v*duy3 - u*dvy3)/v**2)
         dmy4= ivn(i)*4.0D0*COS(idelta(i))*((v*duy4 - u*dvy4)/v**2)
         dmz1= ivn(i)*4.0D0*COS(idelta(i))*((v*duz1 - u*dvz1)/v**2)
         dmz2= ivn(i)*4.0D0*COS(idelta(i))*((v*duz2 - u*dvz2)/v**2)
         dmz3= ivn(i)*4.0D0*COS(idelta(i))*((v*duz3 - u*dvz3)/v**2)
         dmz4= ivn(i)*4.0D0*COS(idelta(i))*((v*duz4 - u*dvz4)/v**2)
      ELSE 
C         PRINT *,'Eeeek! Cock up in Hairyimphell!!!'
      END IF


      RETURN

      END

      SUBROUTINE amberdump(acoords,afilename)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION  ACOORDS(3*NATOMS)
      CHARACTER(LEN=*)      afilename
      INTEGER           dumpint

      OPEN(UNIT=19,FILE=afilename,STATUS='UNKNOWN')
      DO dumpint=1,atoms
         WRITE (19,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F20.10,X,F20.10,X,F20.10)')
     1      label(dumpint),typech(dumpint),dumpint,bondedto(dumpint),
     2      acoords(3*dumpint-2),acoords(3*dumpint-1),acoords(3*dumpint)
      END DO
      WRITE (19,FMT='(A)') 'end'
      WRITE (19,FMT='(A)') ' '
      WRITE (19,FMT='(A4,7X,I2)') 'loop',rings
      DO dumpint=1,rings
         WRITE (19,FMT='(I3,2X,I3)') loopatom(2*dumpint-1),loopatom(2*dumpint)
      END DO
      WRITE (19,FMT='(A)') ' '
      WRITE (19,FMT='(A)') 'charges'
      DO dumpint=1,atoms
         WRITE (19,FMT='(I3,2X,F7.4)') dumpint,pq(dumpint)/18.2223
      END DO
      WRITE (19,FMT='(A)') 'end'
      CLOSE (19)


      RETURN

      END


      SUBROUTINE aconnectdump
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER       ci,cj

      OPEN (UNIT=17,FILE='amber_connectivity',STATUS='UNKNOWN')
      WRITE (17,FMT='(A)') 'ANGLES'
      DO ci=1,ang
         WRITE (17,FMT='(I3,A,I3,A,I3)') aa1(ci),'-',aa2(ci),'-',aa3(ci)
      END DO
      WRITE (17,FMT='(A)') 'TORSIONS'
      DO ci=1,t
         WRITE (17,FMT='(I3,A,I3,A,I3,A,I3,A,4F10.5)') da1(ci),'-',da2(ci),'-',da3(ci),'-',da4(ci),'....',dvn(ci),did(ci),ddelta(ci)
     1         ,dn(ci)
      END DO
      WRITE (17,FMT='(A)') 'IMPROPERS'
      DO ci=1,imp
         WRITE (17,FMT='(I3,A,I3,A,I3,A,I3)') ia1(ci),'-',ia2(ci),'-',ia3(ci),'-',ia4(ci)
      END DO
      WRITE (17,FMT='(A)') 'PAIRWISE PARAMETERS'
      DO ci=1,atoms
        DO cj=1,atoms
          WRITE (17,FMT='(I3,2X,I3,2X,I1,2X,I1,2X,I1,2X,F20.7,2X,F20.7)') 
     1          ci,cj,bonds(ci,cj),
     2          one_three(ci,cj),one_four(ci,cj),vdwa(ci,cj),vdwb(ci,cj)
        END DO
      END DO
      CLOSE (17)


      RETURN

      END
