      SUBROUTINE MKTRAP(N,X,VNEW,ENERGY,SSTEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N,I,J,J1,J2,J3,J4,J5,J6
      DOUBLE PRECISION V0,DELTA,OMEGA,M,B,A,VW,C1,SUM,DIST,
     & SUM2_11,SUM2_12,SUM2_13,SUM2_22,SUM2_23,SUM2_33,VNEW(3*N),ENERGY,SUM1_1,SUM1_2,SUM1_3

      PARAMETER(OMEGA=0.37D-6)
      PARAMETER(M=1.6D4)
      PARAMETER(V0=0.33D-6)
!     PARAMETER(DELTA=0.001575D0) ! weak wall
      PARAMETER(DELTA=0.004725D0) ! strong wall
      PARAMETER(B=0.0139D-4)
      PARAMETER(A=5.29D-6) ! 5.29 micrometres

      DOUBLE PRECISION X(3*N)
      DOUBLE PRECISION R(N,N),R3(N,N),R5(N,N),V(N)
      DOUBLE PRECISION V1(N,3)
      DOUBLE PRECISION V2(N,9)
      LOGICAL SSTEST

      VW=V0*DELTA
      VNEW(1:3*N)=0.0D0
      IF (SSTEST) HESS(1:3*N,1:3*N)=0.0D0

C     N IS THE TOTAL NUMBER OF ATOMS.
C     V0,VW,OMEGA,M,B ARE PARAMETERS CHARACTERIZING THE POTENTIAL.
C     A IS THE LENGTH UNIT.
C     X(3*N) CONTAIN THE COORDINATES OF THE ATOMS IN THE 3 SPATIAL DIRECTIONS X,Y,Z.
C     R(N,N) IS THE MATRIX WHOSE IJ-TH ELEMENT IS 1/(DISTANCE BETWEEN ITH AND JTH ATOMS) 
C     R3(N,N) IS THE MATRIX WHOSE IJ-TH ELEMENT IS 1/(DISTANCE BETWEEN ITH AND JTH ATOMS)**3 
C     R5(N,N) IS THE MATRIX WHOSE IJ-TH ELEMENT IS 1/(DISTANCE BETWEEN ITH AND JTH ATOMS)**5 
C     V(I) IS THE POTENTIAL ENERGY OF THE ITH ATOM.
C     V1(I,1) IS \PARTIAL_X V(I).
C     V1(I,2) IS \PARTIAL_Y V(I).
C     V1(I,3) IS \PARTIAL_Z V(I).
C     V2(I,1) IS \PARTIAL_XX V(I).
C     V2(I,2) IS \PARTIAL_XY V(I).
C     V2(I,3) IS \PARTIAL_XZ V(I).
C     V2(I,4) IS \PARTIAL_YX V(I).
C     V2(I,5) IS \PARTIAL_YY V(I).
C     V2(I,6) IS \PARTIAL_YZ V(I).
C     V2(I,7) IS \PARTIAL_ZX V(I).
C     V2(I,8) IS \PARTIAL_ZY V(I).
C     V2(I,9) IS \PARTIAL_ZZ V(I).

      C1=0.5D0*(OMEGA-M*OMEGA**2-V0)

C     COMPUTING THE DISTANCE MATRIX.
      CALL RMATMK(N,X,R,R3,R5)

      ENERGY=0.0D0
      DO I=1,N

c     THE POTENTIAL. 
       SUM=0.0D0
       DO J=I+1,I+N-1
        IF(J.GT.N)THEN
         J1=MOD(J,N)
        ELSE
         J1=J
        ENDIF
        SUM=SUM+R(I,J1)
       ENDDO 
       V(I)=V0*X(3*(I-1)+3)**2
     & +C1*(X(3*(I-1)+1)**2+X(3*(I-1)+2)**2)
     & +VW*(X(3*(I-1)+1)**2-X(3*(I-1)+2)**2)
     & +B*SUM
        ENERGY=ENERGY+V(I)

c     THE FIRST DERIVATIVES OF THE POTENTIAL. 
       SUM1_1=0.0D0
       SUM1_2=0.0D0
       SUM1_3=0.0D0
       DO J=I+1,I+N-1
        IF(J.GT.N)THEN
         J1=MOD(J,N)
        ELSE
         J1=J
        ENDIF
        SUM1_1=SUM1_1+(X(3*(I-1)+1)-X(3*(J1-1)+1))*R3(I,J1)
        SUM1_2=SUM1_2+(X(3*(I-1)+2)-X(3*(J1-1)+2))*R3(I,J1)
        SUM1_3=SUM1_3+(X(3*(I-1)+3)-X(3*(J1-1)+3))*R3(I,J1)
       ENDDO 
       V1(I,1)=2*(C1+VW)*X(3*(I-1)+1)-B*SUM1_1
       V1(I,2)=2*(C1-VW)*X(3*(I-1)+2)-B*SUM1_2
       V1(I,3)=2*V0*X(3*(I-1)+3)-B*SUM1_3
       VNEW(3*(I-1)+1)=V1(I,1)
       VNEW(3*(I-1)+2)=V1(I,2)
       VNEW(3*(I-1)+3)=V1(I,3)

c     THE SECOND DERIVATIVES OF THE POTENTIAL.
      IF (SSTEST) THEN
       SUM=0.D0
       SUM2_11=0.0D0
       SUM2_12=0.0D0
       SUM2_13=0.0D0
       SUM2_22=0.0D0
       SUM2_23=0.0D0
       SUM2_33=0.0D0
       DO J=I+1,I+N-1
        IF(J.GT.N)THEN
         J1=MOD(J,N)
        ELSE
         J1=J
        ENDIF
        SUM=SUM+R3(I,J1)
        SUM2_11=SUM2_11+((X(3*(I-1)+1)-X(3*(J1-1)+1))**2)*R5(I,J1)
        SUM2_12=SUM2_12+
     &  ((X(3*(I-1)+1)-X(3*(J1-1)+1))*(X(3*(I-1)+2)-X(3*(J1-1)+2)))
     &  *R5(I,J1)
        SUM2_13=SUM2_13+
     &  ((X(3*(I-1)+1)-X(3*(J1-1)+1))*(X(3*(I-1)+3)-X(3*(J1-1)+3)))
     &  *R5(I,J1)
        SUM2_22=SUM2_22+((X(3*(I-1)+2)-X(3*(J1-1)+2))**2)*R5(I,J1)
        SUM2_23=SUM2_23+
     &  ((X(3*(I-1)+2)-X(3*(J1-1)+2))*(X(3*(I-1)+3)-X(3*(J1-1)+3)))
     &  *R5(I,J1)
        SUM2_33=SUM2_33+((X(3*(I-1)+3)-X(3*(J1-1)+3))**2)*R5(I,J1)
       ENDDO 
       V2(I,1)=2*(C1+VW)-B*SUM+3*B*SUM2_11
       V2(I,2)=3*B*SUM2_12
       V2(I,3)=3*B*SUM2_13
       V2(I,5)=2*(C1-VW)-B*SUM+3*B*SUM2_22
       V2(I,6)=3*B*SUM2_23
       V2(I,9)=2*V0-B*SUM+3*B*SUM2_33
       HESS(3*(I-1)+1,3*(I-1)+1)=V2(I,1)
       HESS(3*(I-1)+1,3*(I-1)+2)=V2(I,2)
       HESS(3*(I-1)+1,3*(I-1)+3)=V2(I,3)
       HESS(3*(I-1)+2,3*(I-1)+1)=V2(I,2)
       HESS(3*(I-1)+2,3*(I-1)+2)=V2(I,5)
       HESS(3*(I-1)+2,3*(I-1)+3)=V2(I,6)
       HESS(3*(I-1)+3,3*(I-1)+1)=V2(I,3)
       HESS(3*(I-1)+3,3*(I-1)+2)=V2(I,6)
       HESS(3*(I-1)+3,3*(I-1)+3)=V2(I,9)
      ENDIF
      ENDDO 

      ENERGY=ENERGY/2.0D0 ! energy was out by a factor of two DJW 25/1/16

      IF (.NOT.SSTEST) RETURN

!
!  Include missing second derivatives. DJW 25/1/16
!

C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               J5=3*(J4-1)+J2
               HESS(J3,3*(J4-1)+J2)=-B*3.0D0*R5(J4,J1)*(X(J3)-X(3*(J4-1)+J2))**2+B*R3(J4,J1)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
           ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=-B*R5(J4,J1)*3.0D0
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(J6,J3)=HESS(J3,J6)
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J3,J6)
                  HESS(3*(J1-1)+J5,3*(J4-1)+J2)=HESS(J3,J6)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
      
      SUBROUTINE RMATMK(N,X,R,R3,R5)
      IMPLICIT NONE
      INTEGER N,J1,J2
      DOUBLE PRECISION X(3*N),DISTSQ
      DOUBLE PRECISION R(N,N),R3(N,N),R5(N,N) 
      DO J1=1,N
         R(J1,J1)=0.0D0
         DO J2=J1+1,N
          DISTSQ=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     &        +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     &        +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
          R(J1,J2)=1.0D0/DISTSQ**(0.5D0)
          R3(J1,J2)=1.0D0/DISTSQ**(1.5D0)
          R5(J1,J2)=1.0D0/DISTSQ**(2.5D0)
          R(J2,J1)=R(J1,J2)
          R3(J2,J1)=R3(J1,J2)
          R5(J2,J1)=R5(J1,J2)
         ENDDO
      ENDDO
      RETURN
      END

