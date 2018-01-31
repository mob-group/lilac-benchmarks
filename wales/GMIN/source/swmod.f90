MODULE SWMOD
USE COMMONS, ONLY : NATOMS

IMPLICIT NONE
INTEGER                       :: I, K
DOUBLE PRECISION              :: BOXLX, BOXLY, BOXLZ, BOX(3), LAT(3,3), RI(3)
DOUBLE PRECISION              :: V
DOUBLE PRECISION              :: VSW
DOUBLE PRECISION              :: VFLM
LOGICAL                       :: GTEST, PERIODIC(3)
INTEGER                       :: P, Q
DOUBLE PRECISION              :: AA, BB, A, CSTHTA, LMBDA, GMA, SGMASW, EPSLNSW, SGMAFLM, EPSLNFLM

PUBLIC :: MWINIT, SWINIT, SWTYPE
PRIVATE
CONTAINS

SUBROUTINE MWINIT

AA=7.049556277D0
BB=0.6022245584D0
P=4
Q=0
A=1.8D0
CSTHTA=-1.D0/3.D0
GMA=1.2D0
LMBDA=23.15D0
SGMASW=2.3925D0 ! AA
EPSLNSW=25.89D0 ! kJ mol-1
PERIODIC(:)=.FALSE.

END SUBROUTINE MWINIT

SUBROUTINE SWINIT

AA=7.049556277D0
BB=0.6022245584D0
P=4
Q=0
A=1.8D0
CSTHTA=-1.D0/3.D0
GMA=1.2D0
LMBDA=21.00D0
SGMASW=2.0951D0 ! AA
EPSLNSW=209.2D0 ! kJ mol-1
PERIODIC(:)=.FALSE.

END SUBROUTINE SWINIT

SUBROUTINE SWTYPE (X, G, V, GTEST)
USE COMMONS, ONLY: VT
! calculates energy and gradient for Stillinger-Weber type potentials in reduced units
IMPLICIT NONE
INTEGER                       :: I, J, K
DOUBLE PRECISION              :: X(NATOMS*3), COORDS(NATOMS*3), DMAT(NATOMS,NATOMS), RMAT(NATOMS,NATOMS,3), RI(3), RIJ(3), RET(10)
DOUBLE PRECISION              :: G(NATOMS*3), V
LOGICAL                       :: GTEST
!DOUBLE PRECISION              :: MYROUND

COORDS(:)=X(:)/SGMASW

V   = 0.D0
VT  = 0.D0
G   = 0.D0
DMAT= 0.D0
RMAT= 0.D0

! get minimum vectors and distances
DO I=1,NATOMS
   DO J=I+1,NATOMS
      RIJ=COORDS(J*3-2:J*3)-COORDS(I*3-2:I*3)
      DO K=1,3
         IF (PERIODIC(K)) RIJ=RIJ-LAT(:,K)*MYROUND(RIJ(K)/LAT(K,K))
      ENDDO
      IF (ANY(RIJ(:)>A)) THEN
         DMAT(I,J)=1.0D10
         DMAT(J,I)=DMAT(I,J)
      ELSE
         DMAT(I,J)=SQRT(DOT_PRODUCT(RIJ,RIJ))
         DMAT(J,I)=DMAT(I,J)
         RMAT(I,J,:)=RIJ
         RMAT(J,I,:)=-1.D0*RIJ
      ENDIF
   ENDDO
ENDDO

! calculate SW energy and gradient
DO I=1,NATOMS
   DO J=1,NATOMS
      DO K=I,NATOMS
         IF (((K==I).AND.(I<J)).AND.(DMAT(I,J)<A)) THEN !two-body
            CALL SWPHI2(DMAT(I,J), RET(1), RET(2:), GTEST)
            V=V+RET(1)
            VT(I)=VT(I)+RET(1)
            VT(J)=VT(J)+RET(1)
            IF (GTEST) THEN
               G(I*3-2:I*3)=G(I*3-2:I*3)+RET(2)*RMAT(J,I,:)/DMAT(J,I)
               G(J*3-2:J*3)=G(J*3-2:J*3)+RET(2)*RMAT(I,J,:)/DMAT(I,J)
            ENDIF
         ELSEIF (((I/=J).AND.(J/=K).AND.(I/=K)).AND.((DMAT(J,I)<A).AND.(DMAT(J,K)<A))) THEN !three-body
            CALL SWPHI3(RMAT(J,I,:), RMAT(J,K,:), DMAT(J,I), DMAT(J,K), &
                      RET(1), RET(2:), GTEST)
            V=V+RET(1)
            VT(J)=VT(J)+RET(1)
            IF (GTEST) THEN
               G(I*3-2:I*3)=G(I*3-2:I*3)+RET(2:4)
               G(J*3-2:J*3)=G(J*3-2:J*3)+RET(5:7)
               G(K*3-2:K*3)=G(K*3-2:K*3)+RET(8:10)
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

V=V*EPSLNSW
G=G*EPSLNSW/SGMASW

RETURN
END SUBROUTINE SWTYPE

SUBROUTINE SWPHI2 (DIJ, V, G, GTEST)
! calculates two-body contributions to the energy and gradient
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)  :: DIJ
LOGICAL, INTENT(IN)           :: GTEST
DOUBLE PRECISION, INTENT(OUT) :: V, G(9)
DOUBLE PRECISION              :: DAS, DIJ1

IF (DIJ>A) RETURN
V=0.D0
G=0.D0
DAS=1.D0/(DIJ-A)
DIJ1=1.D0/DIJ
V=AA*(BB*DIJ1**P-DIJ1**Q)*EXP(DAS)

IF (.NOT.GTEST) RETURN

G(1)=AA*DIJ1*(Q*DIJ1**Q-BB*P*DIJ1**P)*EXP(DAS)-V*DAS**2

RETURN
END SUBROUTINE SWPHI2

SUBROUTINE SWPHI3 (RJI, RJK, DJI, DJK, V, G, GTEST)
! calculates three-body contributions to the energy and gradient
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)  :: RJI(3), RJK(3), DJI, DJK
LOGICAL, INTENT(IN)           :: GTEST
DOUBLE PRECISION, INTENT(OUT) :: V, G(9)
DOUBLE PRECISION              :: CSTHTAIJK, DJIAS, DJKAS

IF ((DJI>A).OR.(DJK>A)) RETURN
V=0.D0
G=0.D0
CSTHTAIJK=DOT_PRODUCT(RJI,RJK)/(DJI*DJK)
DJIAS=1.D0/(DJI-A)
DJKAS=1.D0/(DJK-A)
V=LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA)**2

IF (.NOT.GTEST) RETURN
G(1:3)=-1.D0*V*GMA*RJI/DJI*DJIAS**2 &
       +2.D0*LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA) &
       *(RJK/(DJI*DJK)-CSTHTAIJK*RJI/DJI**2)

G(4:6)=+1.D0*V*GMA*(RJI/DJI*DJIAS**2+RJK/DJK*DJKAS**2) &
       -2.D0*LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA) &
       *((RJI+RJK)/(DJI*DJK)-CSTHTAIJK*(RJI/DJI**2+RJK/DJK**2))

G(7:9)=-1.D0*V*GMA*RJK/DJK*DJKAS**2 &
       +2.D0*LMBDA*EXP(GMA*DJIAS+GMA*DJKAS)*(CSTHTAIJK-CSTHTA) &
       *(RJI/(DJI*DJK)-CSTHTAIJK*RJK/DJK**2)

RETURN
END SUBROUTINE SWPHI3

SUBROUTINE HYDROPHOBICPLATES (COORDS, NATOMS, V, G, GTEST, HEIGHT)
! calculates component of energy and gradient corresponding to water-plate interaction in reduced units
IMPLICIT NONE
INTEGER                       :: NATOMS
DOUBLE PRECISION              :: HEIGHT, COORDS(NATOMS*3)
LOGICAL                       :: GTEST
DOUBLE PRECISION              :: G(NATOMS*3), V
INTEGER                       :: I
DOUBLE PRECISION              :: RIZ1, RIZ2

V=0.D0
G=0.D0

DO I=1,NATOMS
   RIZ1=1.D0/MAX(COORDS(I*3),TINY(V))
   RIZ2=1.D0/MAX((HEIGHT-COORDS(I*3)),TINY(V))
   V=V+(((2.D0/15.D0)*RIZ1**9-RIZ1**3)+((2.D0/15.D0)*RIZ2**9-RIZ2**3))
   IF (GTEST) THEN
      G(I*3)=((3.D0*RIZ1**4-9.D0*(2.D0/15.D0)*RIZ1**10)-(3.D0*RIZ2**4-9.D0*(2.D0/15.D0)*RIZ2**10))
   ENDIF
ENDDO

RETURN

END SUBROUTINE HYDROPHOBICPLATES
 
FUNCTION MYROUND(X)
! rounds to nearest integer (no universal intrinsic for this)
! no facility for edge cases (probably very unlikely)
IMPLICIT NONE
DOUBLE PRECISION :: X,MYROUND
MYROUND=SIGN(1.D0*INT(INT(2.D0*SIGN(X,1.D0)+1)/2),X)
END FUNCTION MYROUND

END MODULE SWMOD
