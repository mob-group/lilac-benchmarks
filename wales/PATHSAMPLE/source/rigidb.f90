!
      SUBROUTINE RBMINDIST(RA,RB,NATOMS,DIST,Q2,DEBUG)

!     Follows the prescription of Kearsley, Acta Cryst. A, 45, 208-210, 1989, making necessary changes
!     to conform to right-handed rotation in the right-handed coordinate system.

!     Returns DIST as the actual distance, rather than the squared distance

      USE COMMONS, ONLY: NTSITES, RBSITE, DBPT, DBPTDT, DMBLPYT, MSSTOCKT, STOCKAAT, EFIELDT, MULTISITEPYT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG

      IF ((DBPT .OR. DBPTDT .OR. DMBLPYT .OR. MSSTOCKT .OR. STOCKAAT .OR. MULTISITEPYT) .AND. EFIELDT) THEN
     

         CALL FLDMINDIST(RA,RB,NATOMS,DIST,DEBUG,Q2)
         RETURN

      ENDIF 

      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     MOVE CENTRE OF COORDINATES TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXA = CMXA + RA(J2+1)
         CMYA = CMYA + RA(J2+2)
         CMZA = CMZA + RA(J2+3)
      ENDDO
      CMXA = 2.D0*CMXA/NATOMS; CMYA = 2.D0*CMYA/NATOMS; CMZA = 2.D0*CMZA/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RA(J2+1) = RA(J2+1) - CMXA
         RA(J2+2) = RA(J2+2) - CMYA
         RA(J2+3) = RA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXB = CMXB + RB(J2+1)
         CMYB = CMYB + RB(J2+2)
         CMZB = CMZB + RB(J2+3)
      ENDDO
      CMXB = 2.D0*CMXB/NATOMS; CMYB = 2.D0*CMYB/NATOMS; CMZB = 2.D0*CMZB/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      CALL SITEPOS(RA,XA)
      CALL SITEPOS(RB,XB)



      QMAT(1:4,1:4) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1) 
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - YP*ZM + YM*ZP
         QMAT(1,3) = QMAT(1,3) - XM*ZP + XP*ZM
         QMAT(1,4) = QMAT(1,4) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + YP**2 + ZP**2 + XM**2
         QMAT(2,3) = QMAT(2,3) + XM*YM - XP*YP
         QMAT(2,4) = QMAT(2,4) + XM*ZM - XP*ZP
         QMAT(3,3) = QMAT(3,3) + XP**2 + ZP**2 + YM**2
         QMAT(3,4) = QMAT(3,4) + YM*ZM - YP*ZP
         QMAT(4,4) = QMAT(4,4) + XP**2 + YP**2 + ZM**2
      ENDDO

      QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4) 
      QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)

      CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

      IF (INFO /= 0) PRINT '(A,I6,A)','newmindist> WARNING - INFO=',INFO,' in DSYEV'

      MINV = 1.0D100
      DO J1 = 1,4
         IF (DIAG(J1).LT.MINV) THEN
            JMIN = J1
            MINV = DIAG(J1)
         ENDIF
      ENDDO
      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)

!     IF (DEBUG) PRINT '(A,G20.10,A,I6)',' rbmindist2> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN

      Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE RBMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBNEWROTGEOM(NATOMS,COORDS,Q2,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: COORDS(3*NATOMS), RM(3,3), CX, CY, CZ, R(3), P(3), Q1(4), Q2(4), Q(4)
      DOUBLE PRECISION :: THETA, THETAH, ST, FCT

!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS/2
  
         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) + CZ
      
!     CONVERT THE ANGLE-AXIS COORDINATES

         J      = 3*NATOMS/2 + J
         P(:)   = COORDS(J+1:J+3)

         CALL QROTAA(Q2,P)

         COORDS(J+1:J+3) = P(1:3)

      ENDDO

      END SUBROUTINE RBNEWROTGEOM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE FLDMINDIST(RA,RB,NATOMS,DIST,DEBUG,Q2)

!     returns DIST as the actual distance, rather than the squared distance

      USE COMMONS, ONLY: NTSITES, RBSITE, STOCKAAT, MULTISITEPYT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NTSITES), RB(3*NTSITES), DIST, QMAT(2,2), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: DEBUG

      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     MOVE CENTRE OF COORDINATES TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXA = CMXA + RA(J2+1)
         CMYA = CMYA + RA(J2+2)
         CMZA = CMZA + RA(J2+3)
      ENDDO
      CMXA = 2.D0*CMXA/NATOMS; CMYA = 2.D0*CMYA/NATOMS; CMZA = 2.D0*CMZA/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RA(J2+1) = RA(J2+1) - CMXA
         RA(J2+2) = RA(J2+2) - CMYA
! do not move along the z axis in case of a gravity field (currently implemented for MULTISITEPY)
         IF(MULTISITEPYT) CMZA=0.0D0
         RA(J2+3) = RA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXB = CMXB + RB(J2+1)
         CMYB = CMYB + RB(J2+2)
         CMZB = CMZB + RB(J2+3)
      ENDDO
      CMXB = 2.D0*CMXB/NATOMS; CMYB = 2.D0*CMYB/NATOMS; CMZB = 2.D0*CMZB/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
! do not move along the z axis in case of a gravity field (currently implemented for MULTISITEPY)
         IF(MULTISITEPYT) CMZB=0.0D0
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      IF (STOCKAAT) THEN
         XA(:) = RA(:)
         XB(:) = RB(:)
      ELSE
         CALL SITEPOS(RA,XA)
         CALL SITEPOS(RB,XB)
      ENDIF

      QMAT(1:2,1:2) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1) 
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + XP**2 + YP**2 + ZM**2
      ENDDO

      IF (DEBUG) print *, QMAT
!     QMAT IS SYMMETRIC; QMAT(2,1) = QMAT(1,2)

      MINV = 0.5D0*(QMAT(1,1) + QMAT(2,2) - SQRT(4.D0*QMAT(1,2)*QMAT(1,2) + (QMAT(1,1) - QMAT(2,2))**2.D0))

      Q2(2) = 0.D0; Q2(3) = 0.D0
      
      IF (((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)) < 1.D-12 .OR. QMAT(1,2) < 1.D-12) THEN
         Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
         Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))
      ELSE
         Q2(1) = SQRT(QMAT(1,2)*QMAT(1,2)/((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)))
         Q2(4) = (MINV - QMAT(1,1))*Q2(1)/QMAT(1,2)
      ENDIF

      IF ((QMAT(1,1)<1.D-12).AND.(QMAT(1,2)<1.D-12).AND.(QMAT(2,2)<1.D-12)) THEN
         Q2(1) = 1.D0
         Q2(4) = 0.D0
      ELSEIF (((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)) < 1.D-12 .OR. QMAT(1,2) < 1.D-12) THEN
         Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
         Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))
      ELSE
         Q2(1) = SQRT(QMAT(1,2)*QMAT(1,2)/((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)))
         Q2(4) = (MINV - QMAT(1,1))*Q2(1)/QMAT(1,2)
      ENDIF

      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)
  
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
! do not move along the z axis in case of a gravity field (currently implemented for MULTISITEPY)
         IF(MULTISITEPYT) CMZB=0.0D0
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      IF (DEBUG) print *, 'Q2', Q2
      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE FLDMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBCOMMINDIST(RA,RB,NATOMS,DIST,RM,DEBUG)

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG 

      NSIZE = NATOMS/2
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))
      XA(1:3*NSIZE) = RA(1:3*NSIZE); XB(1:3*NSIZE) = RB(1:3*NSIZE)

!     MOVE CENTRE OF COORDINATES OF XA AND XB TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXA = CMXA + XA(J2+1)
         CMYA = CMYA + XA(J2+2)
         CMZA = CMZA + XA(J2+3)
      ENDDO
      CMXA = CMXA/NSIZE; CMYA = CMYA/NSIZE; CMZA = CMZA/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XA(J2+1) = XA(J2+1) - CMXA
         XA(J2+2) = XA(J2+2) - CMYA
         XA(J2+3) = XA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXB = CMXB + XB(J2+1)
         CMYB = CMYB + XB(J2+2)
         CMZB = CMZB + XB(J2+3)
      ENDDO
      CMXB = CMXB/NSIZE; CMYB = CMYB/NSIZE; CMZB = CMZB/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XB(J2+1) = XB(J2+1) - CMXB
         XB(J2+2) = XB(J2+2) - CMYB
         XB(J2+3) = XB(J2+3) - CMZB
      ENDDO

      QMAT(1:4,1:4) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - YP*ZM + YM*ZP 
         QMAT(1,3) = QMAT(1,3) - XM*ZP + XP*ZM
         QMAT(1,4) = QMAT(1,4) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + YP**2 + ZP**2 + XM**2
         QMAT(2,3) = QMAT(2,3) + XM*YM - XP*YP 
         QMAT(2,4) = QMAT(2,4) + XM*ZM - XP*ZP
         QMAT(3,3) = QMAT(3,3) + XP**2 + ZP**2 + YM**2
         QMAT(3,4) = QMAT(3,4) + YM*ZM - YP*ZP 
         QMAT(4,4) = QMAT(4,4) + XP**2 + YP**2 + ZM**2
      ENDDO

      QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4)
      QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)
      CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

      IF (INFO /= 0) PRINT '(A,I6,A)','newmindist> WARNING - INFO=',INFO,' in DSYEV'

      MINV = 1.0D100
      DO J1 = 1,4
         IF (DIAG(J1).LT.MINV) THEN
            JMIN = J1
            MINV = DIAG(J1)
         ENDIF
      ENDDO
      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)

      IF (DEBUG) PRINT '(A,G20.10,A,I6)',' rbmindist2> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN

      Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE RBCOMMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTMAT(P, RM)

      INTEGER          :: K
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, CT, ST, I3(3,3), E(3,3), RM(3,3)

      I3(:,:) = 0.D0
      DO K = 1, 3
         I3(K,K) = 1.D0
      ENDDO

      THETA2 = DOT_PRODUCT(P,P)

      IF (THETA2 == 0.D0) THEN 

         RM(:,:) = I3(:,:)

      ELSE

         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)

         RM(:,:) = I3(:,:) + (1.D0-CT)*MATMUL(E(:,:),E(:,:)) + ST*E(:,:)

      ENDIF

      END SUBROUTINE ROTMAT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDUM(DBSIGBB)

      USE COMMONS, ONLY: RBSITE

      IMPLICIT NONE

      DOUBLE PRECISION :: DBSIGAA, DBSIGBB

      DBSIGAA = 1.D0

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.5D0*DBSIGAA

      RBSITE(2,1) = 0.D0
      RBSITE(2,2) = 0.D0
      RBSITE(2,3) =-0.5D0*DBSIGBB

      RBSITE(3,1) = 0.D0
      RBSITE(3,2) = 0.D0
      RBSITE(3,3) = 0.D0

      END SUBROUTINE DEFDUM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFDMBLPY(DBSIGBB) 

      USE COMMONS, ONLY: RBSITE

      IMPLICIT NONE

      DOUBLE PRECISION :: DBSIGAA, DBSIGBB

      DBSIGAA = 1.D0

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.5D0*DBSIGAA

      RBSITE(2,1) = 0.D0
      RBSITE(2,2) = 0.D0 
      RBSITE(2,3) =-0.5D0*DBSIGBB

      RBSITE(3,1) = 0.D0
      RBSITE(3,2) = 0.D0
      RBSITE(3,3) = 0.D0

      END SUBROUTINE DEFDMBLPY

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFMULTSTOCK()

      USE COMMONS, ONLY: NRBSITES, RBSITE

      IMPLICIT NONE

      IF (NRBSITES == 4) THEN

         RBSITE(1,:) = (/0.D0, 1.D0/DSQRT(3.D0), 0.D0/)
         RBSITE(2,:) = (/0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
         RBSITE(3,:) = (/-0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
         RBSITE(4,:) = (/ 0.D0, 0.D0, 0.D0/)

      ELSEIF (NRBSITES == 3) THEN

         RBSITE(1,:) = (/0.D0, 1.D0/DSQRT(3.D0), 0.D0/)
         RBSITE(2,:) = (/0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)
         RBSITE(3,:) = (/-0.5D0, -0.5D0/DSQRT(3.D0), 0.D0/)

      ENDIF

      END SUBROUTINE DEFMULTSTOCK

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFBENZENE()

      USE COMMONS, ONLY: RBSITE

      IMPLICIT NONE

!     C6H6

      RBSITE(1,:)  = (/ 2.63923430843701,   0.00000000000000,   0.00000000000000/)        ! C1
      RBSITE(2,:)  = (/ 1.31961715421850,  -2.28564395764590,   0.00000000000000/)        ! C2
      RBSITE(3,:)  = (/-1.31961715421850,  -2.28564395764590,   0.00000000000000/)        ! C3
      RBSITE(4,:)  = (/-2.63923430843701,   0.00000000000000,   0.00000000000000/)        ! C4
      RBSITE(5,:)  = (/-1.31961715421850,   2.28564395764590,   0.00000000000000/)        ! C5
      RBSITE(6,:)  = (/ 1.31961715421850,   2.28564395764590,   0.00000000000000/)        ! C6
      RBSITE(7,:)  = (/ 4.69338981379532,   0.00000000000000,   0.00000000000000/)        ! H1
      RBSITE(8,:)  = (/ 2.34669490689766,  -4.06459480860986,   0.00000000000000/)        ! H2
      RBSITE(9,:)  = (/-2.34669490689766,  -4.06459480860986,   0.00000000000000/)        ! H3
      RBSITE(10,:) = (/-4.69338981379532,   0.00000000000000,   0.00000000000000/)        ! H4
      RBSITE(11,:) = (/-2.34669490689766,   4.06459480860986,   0.00000000000000/)        ! H5
      RBSITE(12,:) = (/ 2.34669490689766,   4.06459480860986,   0.00000000000000/)        ! H6

      END SUBROUTINE DEFBENZENE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTD(RBSITETD)

      IMPLICIT NONE

      DOUBLE PRECISION :: RBSITETD(4,3), FCTR


      FCTR        = 1.D0/DSQRT(8.D0)
!      RBSITETD(1,:) = FCTR*(/ 0.D0, 0.D0, 0.D0/)
!      RBSITETD(2,:) = FCTR*(/ 0.D0, 0.D0, 0.D0/)
!      RBSITETD(3,:) = FCTR*(/ 0.D0, 0.D0, 0.D0/)
!      RBSITETD(4,:) = FCTR*(/ 0.D0, 0.D0, 0.D0/)
!      SIGTD(:)    = (/0.8D0, 0.8D0, 0.8D0, 0.8D0/)
! hk286
      RBSITETD(1,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      RBSITETD(2,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      RBSITETD(3,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)
      RBSITETD(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)
!      SIGTD(:)    = (/1.D0, 0.9D0, 0.8D0, 0.7D0/)


      END SUBROUTINE DEFTD

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPATCHYD()

      USE COMMONS, ONLY: RBSITE

      IMPLICIT NONE

      RBSITE(1,:)= (5.D-1/SQRT(3.D0))*(/  1.D0,  1.D0,  1.D0/)
      RBSITE(2,:)= (5.D-1/SQRT(3.D0))*(/ -1.D0, -1.D0,  1.D0/)
      RBSITE(3,:)= (5.D-1/SQRT(3.D0))*(/ -1.D0,  1.D0, -1.D0/)
      RBSITE(4,:)= (5.D-1/SQRT(3.D0))*(/  1.D0, -1.D0, -1.D0/)

      END SUBROUTINE DEFPATCHYD

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP4(M)
!     TIP4P water

      USE COMMONS, ONLY: NRBSITES, RBSITE

      IMPLICIT NONE

      INTEGER          :: J1
      DOUBLE PRECISION :: M(NRBSITES), MASS, CM(3)
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.D0

      RBSITE(2,1) = 0.D0
      RBSITE(2,2) = SIN(0.5D0*THETA)*ROH
      RBSITE(2,3) = COS(0.5D0*THETA)*ROH

      RBSITE(3,1) = 0.D0
      RBSITE(3,2) = -SIN(0.5D0*THETA)*ROH
      RBSITE(3,3) = COS(0.5D0*THETA)*ROH

      RBSITE(4,1) = 0.D0
      RBSITE(4,2) = 0.D0
      RBSITE(4,3) = ROM

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO J1 = 1, NRBSITES
         CM(:) = CM(:) + M(J1)*RBSITE(J1,:)
         MASS = MASS + M(J1)
      ENDDO
      CM(:) = CM(:)/MASS
      DO J1 = 1, NRBSITES
         RBSITE(J1,:) = RBSITE(J1,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP4

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFPAP()
! jwrm2> Patch anti-patch potential

      USE COMMONS, ONLY: NRBSITES, PAPALP, RBSITE, PAPID

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: LFCT

! Set the distance of sites from the CoM as one third of the equilibrium body-body distance
      LFCT = SQRT(1.D0 + (2.D0/PAPALP)**(1.D0/3.D0))/3.D0

      RBSITE(1,:) = (/0.D0, 0.D0, 0.D0/)

      IF (PAPID == 4) THEN
! Create a tetrahedral arrangement of patches and antipatches
        RBSITE(2,:)= 1.D0/SQRT(3.D0)*(/  1.D0,  1.D0,  1.D0/)
        RBSITE(3,:)= 1.D0/SQRT(3.D0)*(/ -1.D0, -1.D0,  1.D0/)
        RBSITE(4,:)= 1.D0/SQRT(3.D0)*(/ -1.D0,  1.D0, -1.D0/)
        RBSITE(5,:)= 1.D0/SQRT(3.D0)*(/  1.D0, -1.D0, -1.D0/)
      ELSE IF (PAPID == 3) THEN
! Create a patch along the positive z axis and an antipatch along the negative z axis
!        RBSITE(2,:) = (/ 0.D0, 0.D0, 1.D0/)
!        RBSITE(3,:) = (/ 0.D0, 0.D0,-1.D0/)
! Bernal patch 1 only
!        RBSITE(2,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
!        RBSITE(3,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
! Bernal patch 2 only
!        RBSITE(2,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
!        RBSITE(3,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
! Bernal patch 3 only
        RBSITE(2,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
        RBSITE(3,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)

      ELSE IF (PAPID == 2) THEN
! Bernal patches 1 and 2 only
        RBSITE(2,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
        RBSITE(3,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
        RBSITE(4,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
        RBSITE(5,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
! Bernal patches 2 and 3 only
!        RBSITE(2,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
!        RBSITE(3,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
!        RBSITE(4,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
!        RBSITE(5,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)
! Bernal patches 1 and 3 only
!        RBSITE(2,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
!        RBSITE(3,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
!        RBSITE(4,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
!        RBSITE(5,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)

      ELSE IF (PAPID == 1) THEN
! Bernal spiral patches
        RBSITE(2,:)= (/-0.866026D0, 0.387297D0, 0.316227D0/)
        RBSITE(3,:)= (/-0.577348D0,-0.516396D0, 0.632454D0/)
        RBSITE(4,:)= (/-0.096225D0, 0.301232D0, 0.948682D0/)
        RBSITE(5,:)= (/-0.866026D0,-0.387297D0,-0.316227D0/)
        RBSITE(6,:)= (/-0.577348D0, 0.516396D0,-0.632454D0/)
        RBSITE(7,:)= (/-0.096225D0,-0.301232D0,-0.948682D0/)
      ENDIF

! Adjust the distance from the CoM
      DO I = 2, NRBSITES
        RBSITE(I,:) = RBSITE(I,:)*LFCT
      ENDDO

      END SUBROUTINE DEFPAP

!     ---------------------------------------------------------------------------------------------

      SUBROUTINE DEFSILANE()

      USE COMMONS, ONLY: NRBSITES, RBSITE

      IMPLICIT NONE
      DOUBLE PRECISION :: FCTR, M(NRBSITES)
!     SiH4 : Si - 1, H - 2

      FCTR      = 0.85391048D0
      RBSITE(1,:) = (/ 0.0D0, 0.0D0, 0.0D0/)
      RBSITE(2,:) = FCTR*(/ 1.D0, 1.D0, 1.D0/)
      RBSITE(3,:) = FCTR*(/-1.D0,-1.D0, 1.D0/)
      RBSITE(4,:) = FCTR*(/-1.D0, 1.D0,-1.D0/)
      RBSITE(5,:) = FCTR*(/ 1.D0,-1.D0,-1.D0/)

      M(:)  = (/28.1D0, 1.D0, 1.D0, 1.D0, 1.D0/)

      END SUBROUTINE DEFSILANE
