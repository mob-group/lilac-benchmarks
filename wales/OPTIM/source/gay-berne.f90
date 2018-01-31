MODULE PYMODULE

 INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, OFFSET, REALNATOMS
 DOUBLE PRECISION, ALLOCATABLE :: RMIvec(:,:,:), DPI1RMvec(:,:,:), DPI2RMvec(:,:,:), DPI3RMvec(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: PSCALEFAC1vec(:),PSCALEFAC2vec(:),epsilon1(:,:,:),AEZR1(:,:,:), AEZR2(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: pyrepp(:), pyattrp(:)

 DOUBLE PRECISION :: angle,angle2,sigma1(4),cut,ron,ron2,range2inv3,attr(4),vecsbf(3)


 DOUBLE PRECISION :: BOXLX,BOXLY,BOXLZ
 DOUBLE PRECISION ::I3(3,3) 
!
! sf344> multisite PY additions
!
DOUBLE PRECISION, ALLOCATABLE :: PST(:,:),OST(:,:),ELLST1(:,:),ELLST2(:,:),ELLMAT(:,:,:),SITECOORDS(:,:) 
DOUBLE PRECISION, ALLOCATABLE :: SMRBF(:,:,:), SMABF(:,:,:),SMRDET(:),SMADET(:)
LOGICAL, ALLOCATABLE :: RADIFTARRAY(:),RADIFTMATRIX(:,:)
DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793, twopi=6.283185307179586

      ! LJGSITE additions
      INTEGER :: NLJSITE
      DOUBLE PRECISION, ALLOCATABLE :: LJGSITECOORDS(:,:), ljattrp(:), ljrepp(:)

      ! MULTISITEPY2 with periodic boundary conditions
      DOUBLE PRECISION, ALLOCATABLE :: LONGRSAX(:), LONGASAX(:)

END MODULE PYMODULE

MODULE LJCAPSIDMODULE

 INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, OFFSET, REALNATOMS
 DOUBLE PRECISION, ALLOCATABLE :: RMIvec(:,:,:), DPI1RMvec(:,:,:), DPI2RMvec(:,:,:), DPI3RMvec(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: PSCALEFAC1vec(:),PSCALEFAC2vec(:),epsilon1(:,:,:),AEZR1(:,:,:), AEZR2(:,:,:)

 DOUBLE PRECISION :: angle,angle2,sigma1,cut,ron,ron2,range2inv3

 DOUBLE PRECISION ::I3(3,3) 
DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793, twopi=6.283185307179586

END MODULE LJCAPSIDMODULE

SUBROUTINE INITIALISEPYGPERIODIC

use key, only: PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,VT,&
                &       PEPSILONATTR, PSIGMAATTR, LJSITEATTR, LJSITECOORDST, LJSITECOORDS, REALIGNXYZ
use commons, only: natoms, rbsite, nrbsites
use pymodule

implicit none
   WRITE(*,*) ' initialising variables for PY ',NATOMS,LJSITE,BLJSITE,PARAMONOVCUTOFF 
! allocate arrays

    ALLOCATE(RMIvec(natoms/2,3,3),DPI1RMvec(natoms/2,3,3), DPI2RMvec(natoms/2,3,3), DPI3RMvec(natoms/2,3,3))
    ALLOCATE(PSCALEFAC1vec(natoms/2),PSCALEFAC2vec(natoms/2),epsilon1(4,natoms/2,natoms/2))
    ALLOCATE(AEZR1(NATOMS/2,3,3), AEZR2(NATOMS/2,3,3))
    IF(.NOT.ALLOCATED(VT)) ALLOCATE(VT(NATOMS/2))
    IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
    IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))


          vecsbf(1)=1.0D0
          vecsbf(2)=0.0D0
          vecsbf(3)=0.0D0

    IF(LJSITECOORDST) THEN
        vecsbf(:)=LJSITECOORDS(:)
        WRITE(*,'(A,3F8.3)') ' repulsive LJ site coordinates will be ', LJSITECOORDS(:)
    END IF
         NRBSITES = 2
         IF(.NOT.ALLOCATED(RBSITE)) ALLOCATE(RBSITE(NRBSITES,3))

         RBSITE(1,:)=0.0D0
         RBSITE(2,1)=LJSITECOORDS(1)
         RBSITE(2,2)=LJSITECOORDS(2)
         RBSITE(2,3)=LJSITECOORDS(3)

      I3(:,:)    = 0.D0
      AEZR1(:,:,:) = 0.D0
      AEZR2(:,:,:) = 0.D0

!     pi=ATAN(1.0D0)*4.0
!     twopi=2.0D0*pi

     IF(LJSITE) THEN
       epsilon1(:,:,:)=PEPSILON1(1)   ! we are going to use epsilon1 for the extra LJ sites
       sigma1(:)=1.0D0 ! PSIGMA1 is nonexistent from now on, except for the attractive secondary apex sites
       PSCALEFAC1vec(:)=PSCALEFAC1(1)
       PSCALEFAC2vec(:)=PSCALEFAC2(1)
       IF(LJSITEATTR) THEN !attractive secondary apex site is turned on
        attr(:)=1.0D0
        epsilon1(1,:,:)=PEPSILONATTR(1)
        epsilon1(2,:,:)=0.0D0!SQRT(PEPSILONATTR(1)*PEPSILONATTR(2))
        epsilon1(3,:,:)=epsilon1(2,:,:)
        epsilon1(4,:,:)=PEPSILONATTR(2)

        sigma1(1)=PSIGMAATTR(1)
        sigma1(2)=0.5D0*(PSIGMAATTR(1)+PSIGMAATTR(2))
        sigma1(3)=sigma1(2)
        sigma1(4)=PSIGMAATTR(2)
       END IF
     ELSE
       epsilon1(:,:,:)=0.0D0
       sigma1(:)=0.0D0
     END IF

! sanity checks
       IF(PYBINARYT.AND.LJSITE.AND..NOT.BLJSITE) THEN
        WRITE(*,*) 'ERROR --- for binary PY systems with extra LJ sites '// &
                        & 'you have to specify the parameters for both types separately (>3 arguments after EXTRALJSITE)! '
        STOP
       END IF
       IF(BLJSITE.AND..NOT.PYBINARYT) THEN
        WRITE(*,*) 'ERROR --- binary LJ sites specified, but no binary PY particles. '// &
                        & 'EXTRALJSITE should have only 3 arguments!'
        STOP
       END IF


       IF(PYBINARYT.AND.BLJSITE) THEN
        DO J1=1,NATOMS/2
          IF(J1<=PYBINARYTYPE1) THEN
                PSCALEFAC1vec(J1)=PSCALEFAC1(1)
                PSCALEFAC2vec(J1)=PSCALEFAC2(1)
          ELSE
                PSCALEFAC1vec(J1)=PSCALEFAC1(2)
                PSCALEFAC2vec(J1)=PSCALEFAC2(2)
          END IF
        END DO

        DO J1=1,NATOMS/2-1
          DO J2=J1+1,NATOMS/2
           IF(J1<=PYBINARYTYPE1.AND.J2<=PYBINARYTYPE1) THEN
                epsilon1(:,J1,J2)=PEPSILON1(1)
           ELSE IF(J1<=PYBINARYTYPE1.AND.J2>PYBINARYTYPE1) THEN
                epsilon1(:,J1,J2)=PEPSILON1(3)
           ELSE
                epsilon1(:,J1,J2)=PEPSILON1(2)
           END IF
          END DO
        END DO
!        WRITE(*,*) 'epsilon1: ', epsilon1(:,:,:)
       END IF
! cutoff stuff for the extra LJ sites
       cut = PCUTOFF*PCUTOFF ! cutoff squared
       ron = PCUTOFF*0.9D0
       ron2 = ron*ron
       range2inv3=1.0D0/(cut-ron2)**3
!     FROM IUT PARAMETERS

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
      
      DO K1 = 1, 3
        I3(K1,K1) = 1.0D0
      ENDDO
      DO J1=1,REALNATOMS
       DO K1 = 1, 3
         AEZR1(J1,K1,K1) = 1.D0/(PYA1bin(J1,K1)*PYA1bin(J1,K1))
         AEZR2(J1,K1,K1) = 1.D0/(PYA2bin(J1,K1)*PYA2bin(J1,K1))
       ENDDO
      END DO
      IF(REALIGNXYZ) THEN
        CALL PYREALIGNXYZ
        STOP
      END IF
END SUBROUTINE INITIALISEPYGPERIODIC

! PY potential, dc430's implementation
! with PBC and continuous cutoff added
! plus extra LJ site
      SUBROUTINE PYGPERIODIC (X, G, ENERGY, GTEST, STEST)

       use key, only: PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,VT,PEPSILONATTR,PSIGMAATTR
       use commons, only: natoms
       use pymodule
      IMPLICIT NONE

      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), RM(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, SRTFI(2), FMIN, LAMDAC1, LAMDAC2, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DVDF1, DVDF2 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3), ESQ(3,3), DE1(3,3), DE2(3,3), DE3(3,3)
      DOUBLE PRECISION :: DPI1RM(3,3), DPI2RM(3,3), DPI3RM(3,3), DPJ1RM(3,3), DPJ2RM(3,3), DPJ3RM(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      double precision :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST


! sf344 additions

      DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2),CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2), &
                          & CLJ13(2),CLJ14(2),DFDR(2,3),DFP(2,6),dr(6),dCLJ1(2,12),dVDUMMY(12),LJ1(2),DUMMY,DUMMY1,DUMMY2,&
                          & dLJ1(2,12),VDUMMY
      DOUBLE PRECISION :: term2(maxinteractions),term3(maxinteractions), &
                          & xlj(maxinteractions,2,3),rljvec(maxinteractions,3),rljunitvec(maxinteractions,3),&
                          & drlj(maxinteractions,12),rlj(maxinteractions),rlj2(maxinteractions)
      DOUBLE PRECISION :: LLJ(12,maxinteractions), dLLJ1(maxinteractions,12)
      INTEGER          :: k
      LOGICAL          :: STEST
      DOUBLE PRECISION :: cosangle, sinangle, IMAT(3,3), tildematrix(3,3), rotmat(3,3), crossvector(3),tempcrd(6)
     VT(1:NATOMS/2)=0.0D0
     term2(:)=1.0D0
     term3(:)=0.0D0
!       IF(PYBINARYT) THEN
!        DO J1=1,NATOMS/2-1
!          DO J2=J1+1,NATOMS/2
!           IF(J1<=PYBINARYTYPE1.AND.J2>PYBINARYTYPE1) THEN
!                epsilon1(J1,J2)=PEPSILON1/10.0D0
!           ELSE
!                epsilon1(J1,J2)=PEPSILON1
!          END DO
!        END DO
!       END IF

!      DO J1=1,REALNATOMS
!        J2=3*J1
!
!        IF (PARAMONOVPBCX) THEN
!!ensure x component of particle 1 vector is within BoxLx/2 of zero. 
!!If it isn't then subtract integer number of boxlx's such that it is.
!                X(J2-2)=X(J2-2)-BOXLX*NINT(X(J2-2)/BOXLX)
!        ENDIF
!
!        IF (PARAMONOVPBCY) THEN
!!ensure y component of particle 1 vector is within BoxLy/2 of zero. 
!!If it isn't then subtract integer number of boxly's such that it is.
!                X(J2-1)=X(J2-1)-BOXLY*NINT(X(J2-1)/BOXLY)
!        END IF
!
!        IF (PARAMONOVPBCZ) THEN
!!ensure z component of particle 1 vector is within BoxLz/2 of zero. 
!!If it isn't then subtract integer number of boxlz's such that it is.
!                X(J2)=X(J2)-BOXLZ*NINT(X(J2)/BOXLZ)
!        ENDIF
!         
!      END DO

         ENERGY = 0.D0
         G(:)   = 0.D0
        
        DO J1=1, REALNATOMS

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
            angle=sqrt(dot_product(P,P))
 
! temporary stuff to create spheroidal starting geometry for global optimisation
! starting coordinates are taken from global minima of Thomson clusters of the appropriate size
! - we need to figure out the rotation matrix that rotates the rigid body into the vector
! defined by the position of the particle (since that is a unit vector)
! cos(angle)=dot_product(vecsbf,coords(3*J-2:3*J))
!           crossvector(:)=0.0D0
!           cosangle=dot_product(vecsbf,X(J3-2:J3))
!           call  CROSSOPT(vecsbf,X(J3-2:J3),crossvector,0)
!           WRITE(*,*) 'cross vector: ', crossvector(:) 
!           sinangle=sqrt(dot_product(crossvector,crossvector))/(sqrt(dot_product(vecsbf,vecsbf))*sqrt(dot_product(X(J3-2:J3),X(J3-2:J3))))
!           IMAT(:,:)=0.0D0
!           WRITE(*,*) 'cosangle', cosangle, sinangle
!           WRITE(*,*) 'cosangle2+sinangle2=', cosangle*cosangle+sinangle*sinangle
!           DO i=1,3
!              IMAT(i,i)=1.0D0
!           END DO
!           tildematrix(1,1)=0.0D0
!           tildematrix(1,2)=-crossvector(3)
!           tildematrix(1,3)= crossvector(2)
!           tildematrix(2,1)=-tildematrix(1,2)
!           tildematrix(2,2)=0.0D0
!           tildematrix(2,3)=-crossvector(1)
!           tildematrix(3,1)=-tildematrix(1,3)
!           tildematrix(3,2)=-tildematrix(3,2)
!           tildematrix(3,3)=0.0D0
!
!           rotmat=IMAT+sinangle*tildematrix+(1.0D0-cosangle)*MATMUL(tildematrix,tildematrix)
!
!           WRITE(*,'(A,3F8.3)') 'original vector: ', vecsbf(:)
!           WRITE(*,'(A,3F8.3)') 'original vector rotated ', MATMUL(vecsbf,rotmat)
!           WRITE(*,'(A,3F8.3)') 'should have rotated into this: ', X(J3-2:J3)
!           tempcrd(1:6)=0.0D0
!           tempcrd(4)=1.0D0
!           CALL RBNEWROTGEOMMYORIENT(2,tempcrd,rotmat,0.0,0.0,0.0)
!           X(J5-2:J5)=tempcrd(4:6)
!           WRITE(*,'(A,3F8.3)') 'angle-axis after rotation: ',  X(J5-2:J5)
!           IF(J1==264) STOP



            if(angle>twopi) then
! normalise angle-axis coordinates
                X(J5-2:J5)=X(J5-2:J5)/angle
                do
                  angle=angle-twopi
                  if(angle<2*pi) exit
                end do
! multiply with new angle
                X(J5-2:J5)=X(J5-2:J5)*angle
            end if

            CALL RMDRVT(P, RMIvec(J1,:,:), DPI1RMvec(J1,:,:), DPI2RMvec(J1,:,:), DPI3RMvec(J1,:,:), .TRUE.)

        END DO

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
!     ROTATION MATRIX

!            CALL RMDRVT(P, RMI, DPI1RM, DPI2RM, DPI3RM, GTEST)
            RMI(:,:)=RMIvec(J1,:,:)
            DPI1RM(:,:)=DPI1RMvec(J1,:,:)
            DPI2RM(:,:)=DPI2RMvec(J1,:,:)
            DPI3RM(:,:)=DPI3RMvec(J1,:,:)

            AE1 = MATMUL(RMI,(MATMUL(AEZR1(J1,:,:),(TRANSPOSE(RMI)))))

            IF (RADIFT) THEN

               AE2 = MATMUL(RMI,(MATMUL(AEZR2(J1,:,:),(TRANSPOSE(RMI)))))         

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

!               CALL RMDRVT(P, RMJ, DPJ1RM, DPJ2RM, DPJ3RM, GTEST)               
               RMJ(:,:)=RMIvec(J2,:,:)
               DPJ1RM(:,:)=DPI1RMvec(J2,:,:)
               DPJ2RM(:,:)=DPI2RMvec(J2,:,:)
               DPJ3RM(:,:)=DPI3RMvec(J2,:,:)
     
               BE1 = MATMUL(RMJ,(MATMUL(AEZR1(J2,:,:),(TRANSPOSE(RMJ)))))

               IF (RADIFT) THEN
   
                  BE2 = MATMUL(RMJ,(MATMUL(AEZR2(J2,:,:),(TRANSPOSE(RMJ)))))

               ENDIF

!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NR     = RIJ / ABSRIJ

                IF(PARAMONOVCUTOFF.AND.RIJSQ>cut) GOTO 124

!     CALCULATE ECF

               CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC1, FMIN)

               FCNT1   = - FMIN
               SRTFI1  = 1.D0 / DSQRT(FCNT1)
               APB     = LAMDAC1 * AE1 + (1.D0 - LAMDAC1) * BE1
               
               CALL MTRXIN (APB, APBINV)

               ARIBRJ =  LAMDAC1 * MATMUL(AE1,RI) + (1.D0 - LAMDAC1) * MATMUL(BE1,RJ)
               XC     =  MATMUL(APBINV, ARIBRJ)
               XCMRI  = XC - RI
               XCMRJ  = XC - RJ
               DF1DR  = - 2.D0 * LAMDAC1 * MATMUL(AE1,XCMRI)

               D1ABEZ = MATMUL(DPI1RM,AEZR1(J1,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPI2RM,AEZR1(J1,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPI3RM,AEZR1(J1,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

               DF1PI1 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI))
               DF1PI2 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI))
               DF1PI3 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI))

               D1ABEZ = MATMUL(DPJ1RM,AEZR1(J2,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPJ2RM,AEZR1(J2,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPJ3RM,AEZR1(J2,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
               DF1PJ1 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ))
               DF1PJ2 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ))
               DF1PJ3 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ))

               RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
               RHO1SQ = RHO1*RHO1
               RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
               RHO112 = RHO16 * RHO16

               FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
               DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
               DVDF1  = -2.D0*RHO112*RHO1*FCTR1

               IF (RADIFT) THEN

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC2, FMIN)

                  FCNT2   = - FMIN
                  SRTFI2  = 1.D0 / DSQRT(FCNT2)
                  APB     = LAMDAC2 * AE2 + (1.D0 - LAMDAC2) * BE2

                  CALL MTRXIN (APB, APBINV)

                  ARIBRJ =  LAMDAC2 * MATMUL(AE2,RI) + (1.D0 - LAMDAC2) * MATMUL(BE2,RJ)
                  XC     =  MATMUL(APBINV, ARIBRJ)
                  XCMRI  = XC - RI
                  XCMRJ  = XC - RJ
                  DF2DR  = - 2.D0 * LAMDAC2 * MATMUL(AE2,XCMRI)

                  RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                  RHO2SQ = RHO2*RHO2
                  RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                  FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                  DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                  DVDF2  = RHO26*RHO2*FCTR2

                  D1ABEZ = MATMUL(DPI1RM,AEZR2(J1,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPI2RM,AEZR2(J1,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPI3RM,AEZR2(J1,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                  DF2PI1 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D1ABE),XCMRI)
                  DF2PI2 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D2ABE),XCMRI)
                  DF2PI3 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D3ABE),XCMRI)

                  D1ABEZ = MATMUL(DPJ1RM,AEZR2(J2,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPJ2RM,AEZR2(J2,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPJ3RM,AEZR2(J2,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                  DF2PJ1 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D1ABE),XCMRJ)
                  DF2PJ2 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D2ABE),XCMRJ)
                  DF2PJ3 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D3ABE),XCMRJ)

               ELSE

                  RHO2   = RHO1
                  RHO26  = RHO16
                  DG2DR  = DG1DR
                  DVDF2  = RHO26*RHO2*FCTR1
                  DF2PI1 = DF1PI1
                  DF2PI2 = DF1PI2
                  DF2PI3 = DF1PI3
                  DF2PJ1 = DF1PJ1
                  DF2PJ2 = DF1PJ2
                  DF2PJ3 = DF1PJ3

                  SRTFI2 = SRTFI1
                  DF2DR  = DF1DR
                  RHO2SQ = RHO1SQ

               ENDIF


! Correction terms to the potential if we require a cutoff at rc. 
                 SRTFI(1)=SRTFI1
                 SRTFI(2)=SRTFI2
                 DFP(1,1)=DF1PI1
                 DFP(1,2)=DF1PI2
                 DFP(1,3)=DF1PI3
                 DFP(1,4)=DF1PJ1
                 DFP(1,5)=DF1PJ2
                 DFP(1,6)=DF1PJ3
                 DFP(2,1)=DF2PI1
                 DFP(2,2)=DF2PI2
                 DFP(2,3)=DF2PI3
                 DFP(2,4)=DF2PJ1
                 DFP(2,5)=DF2PJ2
                 DFP(2,6)=DF2PJ3
                 DFDR(1,:)=DF1DR(:)
                 DFDR(2,:)=DF2DR(:)
                 LJ1(1)=RHO1
                 LJ1(2)=RHO2
                   IF (PARAMONOVCUTOFF) THEN
! r/SQRT(F1(k)) = f(A,B), a parameter ONLY dependent on orientation!
                         dr(1)=1.0D0/ABSRIJ*(RI(1)-RJ(1))
                         dr(2)=1.0D0/ABSRIJ*(RI(2)-RJ(2))
                         dr(3)=1.0D0/ABSRIJ*(RI(3)-RJ(3))
                         dr(4)=-dr(1)
                         dr(5)=-dr(2)
                         dr(6)=-dr(3)
                     do k=1,2
                          CLJ1(k)=PYSIGNOT/(PCUTOFF-ABSRIJ*SRTFI(k)+PYSIGNOT)
                          CLJ2(k)=CLJ1(k)**2
                          CLJ3(k)=CLJ2(k)*CLJ1(k)
                          CLJ4(k)=CLJ2(k)**2
                          CLJ5(k)=CLJ4(k)*CLJ1(k)
                          CLJ6(k)=CLJ4(k)*CLJ2(k)
                          CLJ7(k)=CLJ6(k)*CLJ1(k)
                          CLJ8(k)=CLJ6(k)*CLJ2(k)
                          CLJ11(k)=CLJ5(k)*CLJ6(k)
                          CLJ12(k)=CLJ6(k)**2
                          CLJ13(k)=CLJ12(k)*CLJ1(k)
                          CLJ14(k)=CLJ7(k)**2

                         DUMMY=CLJ1(k)/PYSIGNOT
                         DUMMY=DUMMY**2
                         DUMMY1=SRTFI(k)
                         DUMMY2=DUMMY1**3

                         DO j=1,3
                            dCLJ1(k,j) = 1.0D0*PYSIGNOT*DUMMY*(-1.0D0*DUMMY1*dr(j)+0.5D0*ABSRIJ*DUMMY2*DFDR(k,j))
                            dCLJ1(k,j+3)= -dCLJ1(k,j)
                         END DO
                         DO j=7,12
                            dCLJ1(k,j) =-1.0D0*PYSIGNOT*DUMMY*(0.5D0*ABSRIJ*DUMMY2*DFP(k,j-6)) !derivatives wrt to orientation
                         END DO
                     end do
                   END IF
 !     CALCULATE PY POTENTIAL ENERGY
123 CONTINUE 
        VDUMMY = 0.0D0
                 
               IF(LJSITE) THEN  
                do k=1,maxinteractions ! k=1 -- interaction between repulsive primary 'apex' sites
                         ! k=2 and k=3 -- interaction between secondary and primary 'apex' sites 
                         ! k=4 -- interaction between secondary 'apex' sites (normal LJ interaction) 
                        ! trying to modify code to allow for binary systems. 
                        ! apex site heights will be defined in absolute units,
                        ! hence PYA1bin(J1,1) etc. will be removed from below

                   IF(k==1) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==2) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==3) THEN
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE 
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   END IF
                        ! first particle
                        xlj(k,1,:)=RI+DUMMY1*MATMUL(RMI,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid

                        ! second particle
                        xlj(k,2,:)=RJ+DUMMY2*MATMUL(RMJ,vecsbf)

                        ! separation between the LJ sites
                        rlj2(k)=(xlj(k,2,1)-xlj(k,1,1))**2+(xlj(k,2,2)-xlj(k,1,2))**2+(xlj(k,2,3)-xlj(k,1,3))**2
                        rlj(k)=sqrt(rlj2(k))
                        rljvec(k,1)=xlj(k,2,1)-xlj(k,1,1)
                        rljvec(k,2)=xlj(k,2,2)-xlj(k,1,2)
                        rljvec(k,3)=xlj(k,2,3)-xlj(k,1,3)

                        DUMMY=1.0D0/rlj(k)
                        rljunitvec(k,:)=rljvec(k,:)*DUMMY !/rlj(k)

                        drlj(k,1)=DUMMY*(xlj(k,2,1)-xlj(k,1,1))         !drlj/dx1
                        drlj(k,2)=DUMMY*(xlj(k,2,2)-xlj(k,1,2))         !drlj/dy1
                        drlj(k,3)=DUMMY*(xlj(k,2,3)-xlj(k,1,3))         !drlj/dz1
                        drlj(k,4)=-drlj(k,1)                               !drlj/dx2
                        drlj(k,5)=-drlj(k,2)                               !drlj/dy2
                        drlj(k,6)=-drlj(k,3)                               !drlj/dz2
                        drlj(k,7) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI1RM,vecsbf)) !drlj/dpx1
                        drlj(k,8) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI2RM,vecsbf)) !drlj/dpy1
                        drlj(k,9) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI3RM,vecsbf)) !drlj/dpz1
                        drlj(k,10) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ1RM,vecsbf)) !drlj/dpx2
                        drlj(k,11) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ2RM,vecsbf)) !drlj/dpy2
                        drlj(k,12) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ3RM,vecsbf)) !drlj/dpz2

              ! interaction between the extra LJ sites:
                        LLJ(1,k)=sigma1(k)*DUMMY !/rlj(k)
                        LLJ(2,k)=LLJ(1,k)**2
                        LLJ(3,k)=LLJ(2,k)*LLJ(1,k)
                        LLJ(4,k)=LLJ(2,k)**2
                        LLJ(5,k)=LLJ(4,k)*LLJ(1,k)
                        LLJ(6,k)=LLJ(4,k)*LLJ(2,k)
                        LLJ(7,k)=LLJ(6,k)*LLJ(1,k)
                        LLJ(11,k)=LLJ(5,k)*LLJ(6,k)
                        LLJ(12,k)=LLJ(6,k)*LLJ(6,k)

!                            DUMMY=1.0D0/rlj(k)
!                            DUMMY=DUMMY**2
                            DO j=1,12
                                dLLJ1(k,j) =-sigma1(k)*DUMMY*DUMMY*drlj(k,j)
                            END DO

!                add corrections to
!                derivatives are zero at rc, and vanish smoothly with no discontinuities. 
!                VDUMMY=epsilon1*(LLJ12-0.0D0*LLJ6) !inner hard core
                 IF (PARAMONOVCUTOFF) THEN
                   IF (rlj(k)>=PCUTOFF) THEN
!                           VDUMMY = 0.0D0
                           term2(k)=1.0D0
                           term3(k)=0.0D0
                   ELSE

                        ! work out the spline terms for smooth cutoff of extra LJ sites
                        if(rlj2(k)<ron2) then
                           term2(k)=1.0D0
                           term3(k)=0.0D0
                        else if(rlj2(k)>ron2) then
                           term2(k)=(cut-rlj(k)**2)*(cut-rlj(k)**2)*(cut+2.0D0*rlj(k)**2-3.0D0*ron2)*range2inv3
                           term3(k)=rlj(k)*12.0D0*(cut-rlj(k)**2)*(ron2-rlj(k)**2)*range2inv3 ! d(term2)/dr
                        end if
                   END IF
!                   VDUMMY = 0.0D0
                 END IF ! IF (PARAMONOVCUTOFF)
!                 VDUMMY=4.0D0*epsilon1*(LLJ12-LLJ6)*term2 !extra LJ site
                 VDUMMY=VDUMMY+4.0D0*epsilon1(k,J1,J2)*term2(k)*(LLJ(12,k) - attr(k)*LLJ(6,k)) ! extra LJ sites
              end do ! k=1,4
             END IF ! IF(LJSITE)

        IF(PARAMONOVCUTOFF) THEN
                !repulsive potential and periodic cutoff corrections
                VDUMMY = VDUMMY + 4.0D0 * ( RHO112 + 6.0D0*CLJ12(1)*CLJ2(1)/RHO1SQ-7.0D0*CLJ12(1))
                    !                      LJ12(1)+(6.0D0*CLJ12(1))*CLJ2(1)/LJ2(1)-7.0D0*CLJ12(1)
                !attractive potential and periodic cutoff corrections
                VDUMMY = VDUMMY + 4.0D0 * (- RHO26 - 3.0D0* CLJ6(2)*CLJ2(2)/RHO2SQ+4.0D0* CLJ6(2))
                    !                      (-LJ6(2)+(-3.0D0*CLJ6(2))*CLJ2(2)/LJ2(2) +4.0D0*CLJ6(2))
        ELSE
                VDUMMY = VDUMMY + 4.0D0 * (RHO112 - RHO26)
        END IF

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! pair potentials

            dVDUMMY(:) = 0.0D0

!     CALCULATE GRADIENT
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0

          IF(PARAMONOVCUTOFF) THEN
            dVDUMMY(:) = 0.0D0
            ! with respect to cartesians
            do j=1,3
              IF(LJSITE) THEN
               do k=1,maxinteractions
!               dVDUMMY(j) = 24.0D0*epsilon1*(2.0D0*LLJ11*dLLJ1(j)-LLJ5*dLLJ1(j))*1.0D0*term2 + &
!                          & 4.0D0*epsilon1*(LLJ12-LLJ6)*term3*drlj(j)! extra LJ site derivatives
                dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, now only repulsive
                dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
               end do
              ELSE
               dVDUMMY(j) = 0.0D0
              END IF
!                dVDUMMY(j) = dLLJ1(j)
               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * ( 2.0D0*RHO1SQ*DG1DR(j)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
                   &+14.0D0*dCLJ1(1,j)*(CLJ13(1)/RHO1SQ-CLJ11(1))) !repulsive derivative

               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * (-1.0D0*RHO2SQ*DG2DR(j)*(RHO2SQ*RHO2SQ*RHO2-CLJ8(2)/(RHO2SQ*RHO2))&
                   &-4.0D0*dCLJ1(2,j)*(CLJ7(2)/RHO2SQ-CLJ5(2))) !attractive derivative
!               dVDUMMY(j) = dCLJ1(1,j)
               FIJ(j) = dVDUMMY(j)
            end do

            do k=1,2
                DUMMY=LJ1(k)/PYSIGNOT
                DUMMY=DUMMY**2
                DUMMY1=SRTFI(k)
                DUMMY2=DUMMY1**3
               do j=7,12
                  dLJ1(k,j) =-1.0D0*PYSIGNOT*DUMMY*(0.5D0*ABSRIJ*DUMMY2*DFP(k,j-6))
               end do
            end do
            ! with respect to orientation vectors
            do j=7,12
             IF(LJSITE) THEN
              do k=1,maxinteractions
!               dVDUMMY(j) = 24.0D0*epsilon1*(2.0D0*LLJ1**11*dLLJ1(j)-LLJ1**5*dLLJ1(j))*term2 + &
!                               & 4.0D0*epsilon1*(LLJ12-LLJ6)*term3*drlj(j)! extra LJ site derivatives
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                               & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, now only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
              end do
             ELSE
               dVDUMMY(j) = 0.0
             END IF
!                dVDUMMY(j) = dLLJ1(j)
               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * ( 2.0D0*dLJ1(1,j)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
                   &+14.0D0*dCLJ1(1,j)*(CLJ13(1)/RHO1SQ-CLJ11(1))) !repulsive derivative

               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * (-1.0D0*dLJ1(2,j)*(RHO2SQ*RHO2SQ*RHO2- CLJ8(2)/(RHO2SQ*RHO2))&
                   &- 4.0D0*dCLJ1(2,j)*( CLJ7(2)/RHO2SQ- CLJ5(2))) !attractive derivative
!               write(*,*) 'dLJ1(2,j)', DFP(2,j-6),j-6
!               dVDUMMY(j) = dCLJ1(1,j)
            end do
            do j=1,3
               TIJ(j) = dVDUMMY(6+j)
               TJI(j) = dVDUMMY(9+j)
            end do
!        ELSE
                 
!        DO j=1,12
!          dvdummy(j)=drlj(j)
!        END DO

          ELSE  !no cutoff
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0
             dVDUMMY(:)=0.0D0

           IF(LJSITE) THEN
            do k=1,maxinteractions
             do j=1,3
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
               FIJ(j) = dVDUMMY(j)
             end do
             do j=7,12
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
             end do
            end do
           END IF
!4.0D0*epsilon0*(12.0D0*LJ11(1)*dLJ1(1,j)-6.0D0*LJ5(2)*dLJ1(2,j))
               FIJ    = FIJ + 24.0D0*(2.D0*RHO112*RHO1*DG1DR - RHO26*RHO2*DG2DR)
               TIJ(1) = DVDF1*DF1PI1 + DVDF2*DF2PI1
               TIJ(2) = DVDF1*DF1PI2 + DVDF2*DF2PI2
               TIJ(3) = DVDF1*DF1PI3 + DVDF2*DF2PI3
               TJI(1) = DVDF1*DF1PJ1 + DVDF2*DF2PJ1
               TJI(2) = DVDF1*DF1PJ2 + DVDF2*DF2PJ2
               TJI(3) = DVDF1*DF1PJ3 + DVDF2*DF2PJ3
               TIJ = dVDUMMY(7:9) + 24.0D0 * TIJ
               TJI = dVDUMMY(10:12) + 24.0D0 * TJI
          END IF

!                G(J3-2)=G(J3-2)+dVDUMMY(1)
!                G(j4-2)=G(j4-2)+dVDUMMY(4)
!                G(J3-1)=G(J3-1)+dVDUMMY(2)
!                G(j4-1)=G(j4-1)+dVDUMMY(5)
!                G(J3  )=G(J3  )+dVDUMMY(3)
!                G(j4  )=G(j4  )+dVDUMMY(6)



!              G(J3-2:J3) = G(J3-2:J3) + dvdummy(1:3)
!              G(J4-2:J4) = G(J4-2:J4) + dvdummy(4:6)
!              G(J5-2:J5) = G(J5-2:J5) + dvdummy(7:9)
!              G(J6-2:J6) = G(J6-2:J6) + dvdummy(10:12)

               G(J3-2:J3) = G(J3-2:J3) - FIJ
               G(J4-2:J4) = G(J4-2:J4) + FIJ
               G(J5-2:J5) = G(J5-2:J5) + TIJ
               G(J6-2:J6) = G(J6-2:J6) + TJI

!     END INNER LOOP OVER PARTICLES
124 CONTINUE
            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO
!        IF(ENERGY<-1.0D5) THEN
!          G=1.0D0
!          ENERGY=0.0D0
!        ELSE
         ENERGY = PYEPSNOT*ENERGY
         G      = PYEPSNOT*G

!        END IF
!998     CONTINUE
        RETURN
      END SUBROUTINE PYGPERIODIC 



! PY potential, dc430's implementation
! with PBC and continuous cutoff added
! plus extra LJ site
!
! This version is addressable with interactions scaled by the input epsilon matrix PYADDEPS
!
      SUBROUTINE PYGPERIODICADD (X, G, ENERGY, GTEST, STEST)

       use key, only: PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,VT,PEPSILONATTR,PSIGMAATTR,PYADDEPS, NADDTARGET
       use commons, only: natoms
       use pymodule
      IMPLICIT NONE

      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), RM(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, SRTFI(2), FMIN, LAMDAC1, LAMDAC2, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DVDF1, DVDF2 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3), ESQ(3,3), DE1(3,3), DE2(3,3), DE3(3,3)
      DOUBLE PRECISION :: DPI1RM(3,3), DPI2RM(3,3), DPI3RM(3,3), DPJ1RM(3,3), DPJ2RM(3,3), DPJ3RM(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      double precision :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST


! sf344 additions

      DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2),CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2), &
                          & CLJ13(2),CLJ14(2),DFDR(2,3),DFP(2,6),dr(6),dCLJ1(2,12),dVDUMMY(12),LJ1(2),DUMMY,DUMMY1,DUMMY2,&
                          & dLJ1(2,12),VDUMMY
      DOUBLE PRECISION :: term2(maxinteractions),term3(maxinteractions), &
                          & xlj(maxinteractions,2,3),rljvec(maxinteractions,3),rljunitvec(maxinteractions,3),&
                          & drlj(maxinteractions,12),rlj(maxinteractions),rlj2(maxinteractions)
      DOUBLE PRECISION :: LLJ(12,maxinteractions), dLLJ1(maxinteractions,12)
      INTEGER          :: k, MJ1, MJ2
      LOGICAL          :: STEST
      DOUBLE PRECISION :: cosangle, sinangle, IMAT(3,3), tildematrix(3,3), rotmat(3,3), crossvector(3),tempcrd(6)
     VT(1:NATOMS/2)=0.0D0
     term2(:)=1.0D0
     term3(:)=0.0D0

         ENERGY = 0.D0
         G(:)   = 0.D0
        
        DO J1=1, REALNATOMS

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
            angle=sqrt(dot_product(P,P))
            if(angle>twopi) then
! normalise angle-axis coordinates
                X(J5-2:J5)=X(J5-2:J5)/angle
                do
                  angle=angle-twopi
                  if(angle<2*pi) exit
                end do
! multiply with new angle
                X(J5-2:J5)=X(J5-2:J5)*angle
            end if

            CALL RMDRVT(P, RMIvec(J1,:,:), DPI1RMvec(J1,:,:), DPI2RMvec(J1,:,:), DPI3RMvec(J1,:,:), .TRUE.)

        END DO

         DO J1 = 1, REALNATOMS - 1

            MJ1=MOD(J1-1,NADDTARGET)+1
            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
!     ROTATION MATRIX

!            CALL RMDRVT(P, RMI, DPI1RM, DPI2RM, DPI3RM, GTEST)
            RMI(:,:)=RMIvec(J1,:,:)
            DPI1RM(:,:)=DPI1RMvec(J1,:,:)
            DPI2RM(:,:)=DPI2RMvec(J1,:,:)
            DPI3RM(:,:)=DPI3RMvec(J1,:,:)

            AE1 = MATMUL(RMI,(MATMUL(AEZR1(J1,:,:),(TRANSPOSE(RMI)))))

            IF (RADIFT) THEN

               AE2 = MATMUL(RMI,(MATMUL(AEZR2(J1,:,:),(TRANSPOSE(RMI)))))         

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               MJ2=MOD(J2-1,NADDTARGET)+1
               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

!               CALL RMDRVT(P, RMJ, DPJ1RM, DPJ2RM, DPJ3RM, GTEST)               
               RMJ(:,:)=RMIvec(J2,:,:)
               DPJ1RM(:,:)=DPI1RMvec(J2,:,:)
               DPJ2RM(:,:)=DPI2RMvec(J2,:,:)
               DPJ3RM(:,:)=DPI3RMvec(J2,:,:)
     
               BE1 = MATMUL(RMJ,(MATMUL(AEZR1(J2,:,:),(TRANSPOSE(RMJ)))))

               IF (RADIFT) THEN
   
                  BE2 = MATMUL(RMJ,(MATMUL(AEZR2(J2,:,:),(TRANSPOSE(RMJ)))))

               ENDIF

!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NR     = RIJ / ABSRIJ

                IF(PARAMONOVCUTOFF.AND.RIJSQ>cut) GOTO 124

!     CALCULATE ECF

               CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC1, FMIN)

               FCNT1   = - FMIN
               SRTFI1  = 1.D0 / DSQRT(FCNT1)
               APB     = LAMDAC1 * AE1 + (1.D0 - LAMDAC1) * BE1
               
               CALL MTRXIN (APB, APBINV)

               ARIBRJ =  LAMDAC1 * MATMUL(AE1,RI) + (1.D0 - LAMDAC1) * MATMUL(BE1,RJ)
               XC     =  MATMUL(APBINV, ARIBRJ)
               XCMRI  = XC - RI
               XCMRJ  = XC - RJ
               DF1DR  = - 2.D0 * LAMDAC1 * MATMUL(AE1,XCMRI)

               D1ABEZ = MATMUL(DPI1RM,AEZR1(J1,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPI2RM,AEZR1(J1,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPI3RM,AEZR1(J1,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

               DF1PI1 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI))
               DF1PI2 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI))
               DF1PI3 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI))

               D1ABEZ = MATMUL(DPJ1RM,AEZR1(J2,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPJ2RM,AEZR1(J2,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPJ3RM,AEZR1(J2,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
               DF1PJ1 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ))
               DF1PJ2 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ))
               DF1PJ3 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ))

               RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
               RHO1SQ = RHO1*RHO1
               RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
               RHO112 = RHO16 * RHO16

               FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
               DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
               DVDF1  = -2.D0*RHO112*RHO1*FCTR1

               IF (RADIFT) THEN

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC2, FMIN)

                  FCNT2   = - FMIN
                  SRTFI2  = 1.D0 / DSQRT(FCNT2)
                  APB     = LAMDAC2 * AE2 + (1.D0 - LAMDAC2) * BE2

                  CALL MTRXIN (APB, APBINV)

                  ARIBRJ =  LAMDAC2 * MATMUL(AE2,RI) + (1.D0 - LAMDAC2) * MATMUL(BE2,RJ)
                  XC     =  MATMUL(APBINV, ARIBRJ)
                  XCMRI  = XC - RI
                  XCMRJ  = XC - RJ
                  DF2DR  = - 2.D0 * LAMDAC2 * MATMUL(AE2,XCMRI)

                  RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                  RHO2SQ = RHO2*RHO2
                  RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                  FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                  DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                  DVDF2  = RHO26*RHO2*FCTR2

                  D1ABEZ = MATMUL(DPI1RM,AEZR2(J1,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPI2RM,AEZR2(J1,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPI3RM,AEZR2(J1,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                  DF2PI1 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D1ABE),XCMRI)
                  DF2PI2 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D2ABE),XCMRI)
                  DF2PI3 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D3ABE),XCMRI)

                  D1ABEZ = MATMUL(DPJ1RM,AEZR2(J2,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPJ2RM,AEZR2(J2,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPJ3RM,AEZR2(J2,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                  DF2PJ1 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D1ABE),XCMRJ)
                  DF2PJ2 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D2ABE),XCMRJ)
                  DF2PJ3 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D3ABE),XCMRJ)

               ELSE

                  RHO2   = RHO1
                  RHO26  = RHO16
                  DG2DR  = DG1DR
                  DVDF2  = RHO26*RHO2*FCTR1
                  DF2PI1 = DF1PI1
                  DF2PI2 = DF1PI2
                  DF2PI3 = DF1PI3
                  DF2PJ1 = DF1PJ1
                  DF2PJ2 = DF1PJ2
                  DF2PJ3 = DF1PJ3

                  SRTFI2 = SRTFI1
                  DF2DR  = DF1DR
                  RHO2SQ = RHO1SQ

               ENDIF


! Correction terms to the potential if we require a cutoff at rc. 
                 SRTFI(1)=SRTFI1
                 SRTFI(2)=SRTFI2
                 DFP(1,1)=DF1PI1
                 DFP(1,2)=DF1PI2
                 DFP(1,3)=DF1PI3
                 DFP(1,4)=DF1PJ1
                 DFP(1,5)=DF1PJ2
                 DFP(1,6)=DF1PJ3
                 DFP(2,1)=DF2PI1
                 DFP(2,2)=DF2PI2
                 DFP(2,3)=DF2PI3
                 DFP(2,4)=DF2PJ1
                 DFP(2,5)=DF2PJ2
                 DFP(2,6)=DF2PJ3
                 DFDR(1,:)=DF1DR(:)
                 DFDR(2,:)=DF2DR(:)
                 LJ1(1)=RHO1
                 LJ1(2)=RHO2
 !     CALCULATE PY POTENTIAL ENERGY
123 CONTINUE 
        VDUMMY = 0.0D0
                 
               IF(LJSITE) THEN  
                do k=1,maxinteractions ! k=1 -- interaction between repulsive primary 'apex' sites
                         ! k=2 and k=3 -- interaction between secondary and primary 'apex' sites 
                         ! k=4 -- interaction between secondary 'apex' sites (normal LJ interaction) 
                        ! trying to modify code to allow for binary systems. 
                        ! apex site heights will be defined in absolute units,
                        ! hence PYA1bin(J1,1) etc. will be removed from below

                   IF(k==1) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==2) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==3) THEN
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE 
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   END IF
                        ! first particle
                        xlj(k,1,:)=RI+DUMMY1*MATMUL(RMI,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid

                        ! second particle
                        xlj(k,2,:)=RJ+DUMMY2*MATMUL(RMJ,vecsbf)

                        ! separation between the LJ sites
                        rlj2(k)=(xlj(k,2,1)-xlj(k,1,1))**2+(xlj(k,2,2)-xlj(k,1,2))**2+(xlj(k,2,3)-xlj(k,1,3))**2
                        rlj(k)=sqrt(rlj2(k))
                        rljvec(k,1)=xlj(k,2,1)-xlj(k,1,1)
                        rljvec(k,2)=xlj(k,2,2)-xlj(k,1,2)
                        rljvec(k,3)=xlj(k,2,3)-xlj(k,1,3)

                        DUMMY=1.0D0/rlj(k)
                        rljunitvec(k,:)=rljvec(k,:)*DUMMY !/rlj(k)

                        drlj(k,1)=DUMMY*(xlj(k,2,1)-xlj(k,1,1))         !drlj/dx1
                        drlj(k,2)=DUMMY*(xlj(k,2,2)-xlj(k,1,2))         !drlj/dy1
                        drlj(k,3)=DUMMY*(xlj(k,2,3)-xlj(k,1,3))         !drlj/dz1
                        drlj(k,4)=-drlj(k,1)                               !drlj/dx2
                        drlj(k,5)=-drlj(k,2)                               !drlj/dy2
                        drlj(k,6)=-drlj(k,3)                               !drlj/dz2
                        drlj(k,7) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI1RM,vecsbf)) !drlj/dpx1
                        drlj(k,8) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI2RM,vecsbf)) !drlj/dpy1
                        drlj(k,9) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI3RM,vecsbf)) !drlj/dpz1
                        drlj(k,10) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ1RM,vecsbf)) !drlj/dpx2
                        drlj(k,11) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ2RM,vecsbf)) !drlj/dpy2
                        drlj(k,12) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ3RM,vecsbf)) !drlj/dpz2

              ! interaction between the extra LJ sites:
                        LLJ(1,k)=sigma1(k)*DUMMY !/rlj(k)
                        LLJ(2,k)=LLJ(1,k)**2
                        LLJ(3,k)=LLJ(2,k)*LLJ(1,k)
                        LLJ(4,k)=LLJ(2,k)**2
                        LLJ(5,k)=LLJ(4,k)*LLJ(1,k)
                        LLJ(6,k)=LLJ(4,k)*LLJ(2,k)
                        LLJ(7,k)=LLJ(6,k)*LLJ(1,k)
                        LLJ(11,k)=LLJ(5,k)*LLJ(6,k)
                        LLJ(12,k)=LLJ(6,k)*LLJ(6,k)

                            DO j=1,12
                                dLLJ1(k,j) =-sigma1(k)*DUMMY*DUMMY*drlj(k,j)
                            END DO

                 VDUMMY=VDUMMY+4.0D0*epsilon1(k,J1,J2)*term2(k)*(LLJ(12,k) - attr(k)*LLJ(6,k)) ! extra LJ sites
              end do ! k=1,4
             END IF ! IF(LJSITE)
                IF (MJ1.EQ.MJ2) THEN ! same address - take just repulsion?
                   VDUMMY = VDUMMY + 4.0D0 * RHO112*PYADDEPS(MJ2,MJ1)
               ELSE
                   VDUMMY = VDUMMY + 4.0D0 * (RHO112 - 1.0D0 * RHO26)*PYADDEPS(MJ2,MJ1)
               ENDIF

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! pair potentials

            dVDUMMY(:) = 0.0D0

!     CALCULATE GRADIENT
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0
             dVDUMMY(:)=0.0D0

           IF(LJSITE) THEN
            do k=1,maxinteractions
             do j=1,3
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
               FIJ(j) = dVDUMMY(j)
             end do
             do j=7,12
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
             end do
            end do
           END IF
              IF (MJ1.EQ.MJ2) THEN ! same address - take just repulsion?
                   FIJ    = FIJ + 24.0D0*2.D0*RHO112*RHO1*DG1DR*PYADDEPS(MJ2,MJ1)
                   TIJ(1) = DVDF1*DF1PI1
                   TIJ(2) = DVDF1*DF1PI2
                   TIJ(3) = DVDF1*DF1PI3
                   TJI(1) = DVDF1*DF1PJ1
                   TJI(2) = DVDF1*DF1PJ2
                   TJI(3) = DVDF1*DF1PJ3
                   TIJ = dVDUMMY(7:9) + 24.0D0 * TIJ*PYADDEPS(MJ2,MJ1)
                   TJI = dVDUMMY(10:12) + 24.0D0 * TJI*PYADDEPS(MJ2,MJ1)
               ELSE
                   FIJ    = FIJ + 24.0D0*(2.D0*RHO112*RHO1*DG1DR - RHO26*RHO2*DG2DR)*PYADDEPS(MJ2,MJ1)
                   TIJ(1) = DVDF1*DF1PI1 + 1.0D0* DVDF2*DF2PI1
                   TIJ(2) = DVDF1*DF1PI2 + 1.0D0*  DVDF2*DF2PI2
                   TIJ(3) = DVDF1*DF1PI3 + 1.0D0*  DVDF2*DF2PI3
                   TJI(1) = DVDF1*DF1PJ1 + 1.0D0*  DVDF2*DF2PJ1
                   TJI(2) = DVDF1*DF1PJ2 + 1.0D0*  DVDF2*DF2PJ2
                   TJI(3) = DVDF1*DF1PJ3 + 1.0D0*  DVDF2*DF2PJ3
                   TIJ = dVDUMMY(7:9) + 24.0D0 * TIJ*PYADDEPS(MJ2,MJ1)
                   TJI = dVDUMMY(10:12) + 24.0D0 * TJI*PYADDEPS(MJ2,MJ1)
               ENDIF

               G(J3-2:J3) = G(J3-2:J3) - FIJ
               G(J4-2:J4) = G(J4-2:J4) + FIJ
               G(J5-2:J5) = G(J5-2:J5) + TIJ
               G(J6-2:J6) = G(J6-2:J6) + TJI

!     END INNER LOOP OVER PARTICLES
124 CONTINUE
            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO
         ENERGY = PYEPSNOT*ENERGY
         G      = PYEPSNOT*G

        RETURN
      END SUBROUTINE PYGPERIODICADD


! PY potential, dc430's implementation
! with PBC and continuous cutoff added
! plus extra LJ site
!
! This version is addressable with interactions scaled by the input epsilon matrix PYADDEPS 
! with separate rep and att entries
!
      SUBROUTINE PYGPERIODICADD2 (X, G, ENERGY, GTEST, STEST)

       use key, only: PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,VT,PEPSILONATTR,PSIGMAATTR,PYADDREP, PYADDATT, NADDTARGET
       use commons, only: natoms
       use pymodule
      IMPLICIT NONE

      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), RM(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, SRTFI(2), FMIN, LAMDAC1, LAMDAC2, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DVDF1, DVDF2 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3), ESQ(3,3), DE1(3,3), DE2(3,3), DE3(3,3)
      DOUBLE PRECISION :: DPI1RM(3,3), DPI2RM(3,3), DPI3RM(3,3), DPJ1RM(3,3), DPJ2RM(3,3), DPJ3RM(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      double precision :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST


! sf344 additions

      DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2),CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2), &
                          & CLJ13(2),CLJ14(2),DFDR(2,3),DFP(2,6),dr(6),dCLJ1(2,12),dVDUMMY(12),LJ1(2),DUMMY,DUMMY1,DUMMY2,&
                          & dLJ1(2,12),VDUMMY
      DOUBLE PRECISION :: term2(maxinteractions),term3(maxinteractions), &
                          & xlj(maxinteractions,2,3),rljvec(maxinteractions,3),rljunitvec(maxinteractions,3),&
                          & drlj(maxinteractions,12),rlj(maxinteractions),rlj2(maxinteractions)
      DOUBLE PRECISION :: LLJ(12,maxinteractions), dLLJ1(maxinteractions,12)
      INTEGER          :: k, MJ1, MJ2
      LOGICAL          :: STEST
      DOUBLE PRECISION :: cosangle, sinangle, IMAT(3,3), tildematrix(3,3), rotmat(3,3), crossvector(3),tempcrd(6)
     VT(1:NATOMS/2)=0.0D0
     term2(:)=1.0D0
     term3(:)=0.0D0

         ENERGY = 0.D0
         G(:)   = 0.D0
        
        DO J1=1, REALNATOMS

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
            angle=sqrt(dot_product(P,P))
            if(angle>twopi) then
! normalise angle-axis coordinates
                X(J5-2:J5)=X(J5-2:J5)/angle
                do
                  angle=angle-twopi
                  if(angle<2*pi) exit
                end do
! multiply with new angle
                X(J5-2:J5)=X(J5-2:J5)*angle
            end if

            CALL RMDRVT(P, RMIvec(J1,:,:), DPI1RMvec(J1,:,:), DPI2RMvec(J1,:,:), DPI3RMvec(J1,:,:), .TRUE.)

        END DO

         DO J1 = 1, REALNATOMS - 1

            MJ1=MOD(J1-1,NADDTARGET)+1
            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
!     ROTATION MATRIX

!            CALL RMDRVT(P, RMI, DPI1RM, DPI2RM, DPI3RM, GTEST)
            RMI(:,:)=RMIvec(J1,:,:)
            DPI1RM(:,:)=DPI1RMvec(J1,:,:)
            DPI2RM(:,:)=DPI2RMvec(J1,:,:)
            DPI3RM(:,:)=DPI3RMvec(J1,:,:)

            AE1 = MATMUL(RMI,(MATMUL(AEZR1(J1,:,:),(TRANSPOSE(RMI)))))

            IF (RADIFT) THEN

               AE2 = MATMUL(RMI,(MATMUL(AEZR2(J1,:,:),(TRANSPOSE(RMI)))))         

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               MJ2=MOD(J2-1,NADDTARGET)+1
               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

!               CALL RMDRVT(P, RMJ, DPJ1RM, DPJ2RM, DPJ3RM, GTEST)               
               RMJ(:,:)=RMIvec(J2,:,:)
               DPJ1RM(:,:)=DPI1RMvec(J2,:,:)
               DPJ2RM(:,:)=DPI2RMvec(J2,:,:)
               DPJ3RM(:,:)=DPI3RMvec(J2,:,:)
     
               BE1 = MATMUL(RMJ,(MATMUL(AEZR1(J2,:,:),(TRANSPOSE(RMJ)))))

               IF (RADIFT) THEN
   
                  BE2 = MATMUL(RMJ,(MATMUL(AEZR2(J2,:,:),(TRANSPOSE(RMJ)))))

               ENDIF

!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NR     = RIJ / ABSRIJ

                IF(PARAMONOVCUTOFF.AND.RIJSQ>cut) GOTO 124

!     CALCULATE ECF

               CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC1, FMIN)

               FCNT1   = - FMIN
               SRTFI1  = 1.D0 / DSQRT(FCNT1)
               APB     = LAMDAC1 * AE1 + (1.D0 - LAMDAC1) * BE1
               
               CALL MTRXIN (APB, APBINV)

               ARIBRJ =  LAMDAC1 * MATMUL(AE1,RI) + (1.D0 - LAMDAC1) * MATMUL(BE1,RJ)
               XC     =  MATMUL(APBINV, ARIBRJ)
               XCMRI  = XC - RI
               XCMRJ  = XC - RJ
               DF1DR  = - 2.D0 * LAMDAC1 * MATMUL(AE1,XCMRI)

               D1ABEZ = MATMUL(DPI1RM,AEZR1(J1,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPI2RM,AEZR1(J1,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPI3RM,AEZR1(J1,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

               DF1PI1 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI))
               DF1PI2 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI))
               DF1PI3 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI))

               D1ABEZ = MATMUL(DPJ1RM,AEZR1(J2,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPJ2RM,AEZR1(J2,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPJ3RM,AEZR1(J2,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
               DF1PJ1 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ))
               DF1PJ2 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ))
               DF1PJ3 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ))

               RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
               RHO1SQ = RHO1*RHO1
               RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
               RHO112 = RHO16 * RHO16

               FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
               DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
               DVDF1  = -2.D0*RHO112*RHO1*FCTR1

               IF (RADIFT) THEN

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC2, FMIN)

                  FCNT2   = - FMIN
                  SRTFI2  = 1.D0 / DSQRT(FCNT2)
                  APB     = LAMDAC2 * AE2 + (1.D0 - LAMDAC2) * BE2

                  CALL MTRXIN (APB, APBINV)

                  ARIBRJ =  LAMDAC2 * MATMUL(AE2,RI) + (1.D0 - LAMDAC2) * MATMUL(BE2,RJ)
                  XC     =  MATMUL(APBINV, ARIBRJ)
                  XCMRI  = XC - RI
                  XCMRJ  = XC - RJ
                  DF2DR  = - 2.D0 * LAMDAC2 * MATMUL(AE2,XCMRI)

                  RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                  RHO2SQ = RHO2*RHO2
                  RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                  FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                  DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                  DVDF2  = RHO26*RHO2*FCTR2

                  D1ABEZ = MATMUL(DPI1RM,AEZR2(J1,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPI2RM,AEZR2(J1,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPI3RM,AEZR2(J1,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                  DF2PI1 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D1ABE),XCMRI)
                  DF2PI2 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D2ABE),XCMRI)
                  DF2PI3 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D3ABE),XCMRI)

                  D1ABEZ = MATMUL(DPJ1RM,AEZR2(J2,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPJ2RM,AEZR2(J2,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPJ3RM,AEZR2(J2,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                  DF2PJ1 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D1ABE),XCMRJ)
                  DF2PJ2 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D2ABE),XCMRJ)
                  DF2PJ3 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D3ABE),XCMRJ)

               ELSE

                  RHO2   = RHO1
                  RHO26  = RHO16
                  DG2DR  = DG1DR
                  DVDF2  = RHO26*RHO2*FCTR1
                  DF2PI1 = DF1PI1
                  DF2PI2 = DF1PI2
                  DF2PI3 = DF1PI3
                  DF2PJ1 = DF1PJ1
                  DF2PJ2 = DF1PJ2
                  DF2PJ3 = DF1PJ3

                  SRTFI2 = SRTFI1
                  DF2DR  = DF1DR
                  RHO2SQ = RHO1SQ

               ENDIF


! Correction terms to the potential if we require a cutoff at rc. 
                 SRTFI(1)=SRTFI1
                 SRTFI(2)=SRTFI2
                 DFP(1,1)=DF1PI1
                 DFP(1,2)=DF1PI2
                 DFP(1,3)=DF1PI3
                 DFP(1,4)=DF1PJ1
                 DFP(1,5)=DF1PJ2
                 DFP(1,6)=DF1PJ3
                 DFP(2,1)=DF2PI1
                 DFP(2,2)=DF2PI2
                 DFP(2,3)=DF2PI3
                 DFP(2,4)=DF2PJ1
                 DFP(2,5)=DF2PJ2
                 DFP(2,6)=DF2PJ3
                 DFDR(1,:)=DF1DR(:)
                 DFDR(2,:)=DF2DR(:)
                 LJ1(1)=RHO1
                 LJ1(2)=RHO2
 !     CALCULATE PY POTENTIAL ENERGY
123 CONTINUE 
        VDUMMY = 0.0D0
                 
               IF(LJSITE) THEN  
                do k=1,maxinteractions ! k=1 -- interaction between repulsive primary 'apex' sites
                         ! k=2 and k=3 -- interaction between secondary and primary 'apex' sites 
                         ! k=4 -- interaction between secondary 'apex' sites (normal LJ interaction) 
                        ! trying to modify code to allow for binary systems. 
                        ! apex site heights will be defined in absolute units,
                        ! hence PYA1bin(J1,1) etc. will be removed from below

                   IF(k==1) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==2) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==3) THEN
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE 
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   END IF
                        ! first particle
                        xlj(k,1,:)=RI+DUMMY1*MATMUL(RMI,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid

                        ! second particle
                        xlj(k,2,:)=RJ+DUMMY2*MATMUL(RMJ,vecsbf)

                        ! separation between the LJ sites
                        rlj2(k)=(xlj(k,2,1)-xlj(k,1,1))**2+(xlj(k,2,2)-xlj(k,1,2))**2+(xlj(k,2,3)-xlj(k,1,3))**2
                        rlj(k)=sqrt(rlj2(k))
                        rljvec(k,1)=xlj(k,2,1)-xlj(k,1,1)
                        rljvec(k,2)=xlj(k,2,2)-xlj(k,1,2)
                        rljvec(k,3)=xlj(k,2,3)-xlj(k,1,3)

                        DUMMY=1.0D0/rlj(k)
                        rljunitvec(k,:)=rljvec(k,:)*DUMMY !/rlj(k)

                        drlj(k,1)=DUMMY*(xlj(k,2,1)-xlj(k,1,1))         !drlj/dx1
                        drlj(k,2)=DUMMY*(xlj(k,2,2)-xlj(k,1,2))         !drlj/dy1
                        drlj(k,3)=DUMMY*(xlj(k,2,3)-xlj(k,1,3))         !drlj/dz1
                        drlj(k,4)=-drlj(k,1)                               !drlj/dx2
                        drlj(k,5)=-drlj(k,2)                               !drlj/dy2
                        drlj(k,6)=-drlj(k,3)                               !drlj/dz2
                        drlj(k,7) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI1RM,vecsbf)) !drlj/dpx1
                        drlj(k,8) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI2RM,vecsbf)) !drlj/dpy1
                        drlj(k,9) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI3RM,vecsbf)) !drlj/dpz1
                        drlj(k,10) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ1RM,vecsbf)) !drlj/dpx2
                        drlj(k,11) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ2RM,vecsbf)) !drlj/dpy2
                        drlj(k,12) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ3RM,vecsbf)) !drlj/dpz2

              ! interaction between the extra LJ sites:
                        LLJ(1,k)=sigma1(k)*DUMMY !/rlj(k)
                        LLJ(2,k)=LLJ(1,k)**2
                        LLJ(3,k)=LLJ(2,k)*LLJ(1,k)
                        LLJ(4,k)=LLJ(2,k)**2
                        LLJ(5,k)=LLJ(4,k)*LLJ(1,k)
                        LLJ(6,k)=LLJ(4,k)*LLJ(2,k)
                        LLJ(7,k)=LLJ(6,k)*LLJ(1,k)
                        LLJ(11,k)=LLJ(5,k)*LLJ(6,k)
                        LLJ(12,k)=LLJ(6,k)*LLJ(6,k)

                            DO j=1,12
                                dLLJ1(k,j) =-sigma1(k)*DUMMY*DUMMY*drlj(k,j)
                            END DO

                 VDUMMY=VDUMMY+4.0D0*epsilon1(k,J1,J2)*term2(k)*(LLJ(12,k) - attr(k)*LLJ(6,k)) ! extra LJ sites
              end do ! k=1,4
             END IF ! IF(LJSITE)
                VDUMMY = VDUMMY + 4.0D0 * (RHO112*PYADDREP(MJ2,MJ1) - RHO26*PYADDATT(MJ2,MJ1))

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! pair potentials

            dVDUMMY(:) = 0.0D0

!     CALCULATE GRADIENT
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0
             dVDUMMY(:)=0.0D0

           IF(LJSITE) THEN
            do k=1,maxinteractions
             do j=1,3
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
               FIJ(j) = dVDUMMY(j)
             end do
             do j=7,12
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
             end do
            end do
           END IF

                   FIJ    = FIJ + 24.0D0*(2.D0*RHO112*RHO1*DG1DR*PYADDREP(MJ2,MJ1) - 1.0D0 * RHO26*RHO2*DG2DR*PYADDATT(MJ2,MJ1))
                   TIJ(1) = DVDF1*DF1PI1*PYADDREP(MJ2,MJ1) + 1.0D0* DVDF2*DF2PI1*PYADDATT(MJ2,MJ1)
                   TIJ(2) = DVDF1*DF1PI2*PYADDREP(MJ2,MJ1) + 1.0D0*  DVDF2*DF2PI2*PYADDATT(MJ2,MJ1)
                   TIJ(3) = DVDF1*DF1PI3*PYADDREP(MJ2,MJ1) + 1.0D0*  DVDF2*DF2PI3*PYADDATT(MJ2,MJ1)
                   TJI(1) = DVDF1*DF1PJ1*PYADDREP(MJ2,MJ1) + 1.0D0*  DVDF2*DF2PJ1*PYADDATT(MJ2,MJ1)
                   TJI(2) = DVDF1*DF1PJ2*PYADDREP(MJ2,MJ1) + 1.0D0*  DVDF2*DF2PJ2*PYADDATT(MJ2,MJ1)
                   TJI(3) = DVDF1*DF1PJ3*PYADDREP(MJ2,MJ1) + 1.0D0*  DVDF2*DF2PJ3*PYADDATT(MJ2,MJ1)
                   TIJ = dVDUMMY(7:9) + 24.0D0 * TIJ
                   TJI = dVDUMMY(10:12) + 24.0D0 * TJI

               G(J3-2:J3) = G(J3-2:J3) - FIJ
               G(J4-2:J4) = G(J4-2:J4) + FIJ
               G(J5-2:J5) = G(J5-2:J5) + TIJ
               G(J6-2:J6) = G(J6-2:J6) + TJI

!     END INNER LOOP OVER PARTICLES
124 CONTINUE
            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO
         ENERGY = PYEPSNOT*ENERGY
         G      = PYEPSNOT*G

        RETURN
      END SUBROUTINE PYGPERIODICADD2


SUBROUTINE PARAMONOVNUMFIRSTDER(OLDX,STEST)
use commons, only : natoms
!use key
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),ENERGY,X(3*NATOMS),NUMGRAD(3*NATOMS),TRUEGRAD(3*NATOMS),&
        &            OLDX(3*NATOMS),ETEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST


ksi=0.00001D0
X(:)=OLDX(:)

!         CALL OLDPARAMONOV(X,TRUEGRAD,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,TRUEGRAD,ENERGY,.TRUE.,.FALSE.)


!WRITE(*,*) 'sf344> in Paramonovsecder'
DO i=1,3*NATOMS

         X(i)=X(i)-ksi

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.,.FALSE.)

         ETEMP(1,i)=ENERGY

         X(i)=X(i)+2.0D0*ksi

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.,.FALSE.)
    
         ETEMP(2,i)=ENERGY
         NUMGRAD(i)=(ETEMP(2,i)-ETEMP(1,i))/(2.0D0*ksi)

         X(i)=X(i)-ksi

WRITE(*,'(3F20.10)') TRUEGRAD(i),NUMGRAD(i), 100.0D0*(TRUEGRAD(i)-NUMGRAD(i))/TRUEGRAD(i)

END DO

STOP
!WRITE(*,*) 'sf344> exiting AMBERsecder'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)
END SUBROUTINE PARAMONOVNUMFIRSTDER

SUBROUTINE PYGSECDER(OLDX,STEST)
use commons
!use keyword
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST


ksi=0.00001D0
X(:)=OLDX(:)

!WRITE(*,*) 'sf344> in Paramonovsecder', NATOMS, X(1)
!ksi=0.0001D0
!X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
!DO i=1,3*NATOMS
!
!         X(i)=X(i)-ksi
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(1,:)=V(:)
!
!         X(i)=X(i)+2.0D0*ksi
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(2,:)=V(:)
!
!                DO j=i,3*NATOM
!                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
!                        HESS(j,i)=HESS(i,j)
!                END DO
!END DO
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO i=1,3*NATOMS

!       IF ((i.GT.3*NATOMS/2).and.(MOD(i,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(i)=X(i)-ksi
 
!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYG(X,V,EGB,.TRUE.,.FALSE.)
!         WRITE(*,*) 'EGB1',EGB
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYG(X,V,EGB,.TRUE.,.FALSE.)
!         WRITE(*,*) 'EGB2',EGB

         VTEMP(2,:)=V(:)

                DO j=i,3*NATOMS
                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                        HESS(j,i)=HESS(i,j)
                END DO
!        END IF
END DO
!WRITE(*,*) 'sf344> exiting Paramonovsecder', NATOMS, X(1)
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'(12F10.3)') HESS(:,:)

END SUBROUTINE PYGSECDER


SUBROUTINE PYGPERIODICSECDER(OLDX,STEST)
use commons
!use keyword
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST


ksi=0.00001D0
X(:)=OLDX(:)

!WRITE(*,*) 'sf344> in Paramonovsecder', NATOMS, X(1)
!ksi=0.0001D0
!X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
!DO i=1,3*NATOMS
!
!         X(i)=X(i)-ksi
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(1,:)=V(:)
!
!         X(i)=X(i)+2.0D0*ksi
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(2,:)=V(:)
!
!                DO j=i,3*NATOM
!                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
!                        HESS(j,i)=HESS(i,j)
!                END DO
!END DO
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO i=1,3*NATOMS

!       IF ((i.GT.3*NATOMS/2).and.(MOD(i,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(i)=X(i)-ksi
 
!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,EGB,GTEST,.FALSE.)
!         WRITE(*,*) 'EGB1',EGB
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,EGB,GTEST,.FALSE.)
!         WRITE(*,*) 'EGB2',EGB

         VTEMP(2,:)=V(:)

                DO j=i,3*NATOMS
                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                        HESS(j,i)=HESS(i,j)
                END DO
!        END IF
END DO
!WRITE(*,*) 'sf344> exiting Paramonovsecder', NATOMS, X(1)
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'(12F10.3)') HESS(:,:)

END SUBROUTINE PYGPERIODICSECDER

SUBROUTINE PYGPERIODICSECDERADD(OLDX,STEST)
use commons
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST

ksi=0.00001D0
X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO i=1,3*NATOMS

         X(i)=X(i)-ksi
 
         CALL PYGPERIODICADD(X,V,EGB,.TRUE.,.FALSE.)
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

         CALL PYGPERIODICADD(X,V,EGB,.TRUE.,.FALSE.)

         VTEMP(2,:)=V(:)

                DO j=i,3*NATOMS
                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                        HESS(j,i)=HESS(i,j)
                END DO
END DO

END SUBROUTINE PYGPERIODICSECDERADD


SUBROUTINE PYGPERIODICSECDERADD2(OLDX,STEST)
use commons
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST

ksi=0.00001D0
X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO i=1,3*NATOMS

         X(i)=X(i)-ksi
 
         CALL PYGPERIODICADD2(X,V,EGB,.TRUE.,.FALSE.)
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

         CALL PYGPERIODICADD2(X,V,EGB,.TRUE.,.FALSE.)

         VTEMP(2,:)=V(:)

                DO j=i,3*NATOMS
                        HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                        HESS(j,i)=HESS(i,j)
                END DO
END DO

END SUBROUTINE PYGPERIODICSECDERADD2

SUBROUTINE EllipsoidsAAtoPolar(px1,py1,pz1,alpha,beta,gamma,alphadeg,betadeg,gammadeg)
! converts angle-axis coordinates px, py, pz to "polar-like" angles alpha and beta
! px=cos(alpha)*cos(alpha)
! py=cos(alpha)*sin(beta)
! pz=sin(alpha)
!use commons
implicit none

DOUBLE PRECISION        :: px1, py1, pz1,px,py,pz

DOUBLE PRECISION,intent(out)       :: alpha, beta, gamma, alphadeg, betadeg, gammadeg
DOUBLE PRECISION                   :: PI,twopi

! alpha: angle of the vector with the xy plane
! beta: angle of the vector with the x axis
! gamma: angle of rotation about the z axis (angle-axis convention),
!        provided by the magnitude of the vector (sqrt(px**2+py**2+pz**2)-1)
!        so that a vector of magnitude 1 means no rotation is performed

!IF(pz*py*pz<0) THEN
!       px = -px
!       py = -py
!       pz = -pz
!END IF

PI = atan(1.0)*4.0D0
twopi = 2.0D0*PI

     gamma = sqrt(px1**2+py1**2+pz1**2)

     px = px1/gamma
     py = py1/gamma
     pz = pz1/gamma

     IF(px.eq.0.0d0) THEN
       IF(py>=0.0D0) THEN
        alpha = PI/2           ! Euler angle alpha 
       ELSE
        alpha = -1.0D0*PI/2
       END IF
     ELSE
        IF(py>=0.0D0) THEN
          IF(px>0.0D0) THEN
            alpha = 1.0D0*atan(py/px)   ! first quadrant
          ELSE    ! px<0
            alpha = 1.0D0*atan(py/px)+PI       ! should be in the second quadrant
          END IF
        ELSE IF(py<0.0D0) THEN
          IF(px>0.0D0) THEN
            alpha = 1.0D0*atan(py/px)    ! fourth quadrant
          ELSE    ! px<0
            alpha = 1.0D0*atan(py/px)-PI            ! third quadrant
          END IF
        END IF 
     END IF
        beta = 1.0D0*acos(pz)

alphadeg=alpha*180/Pi
betadeg=beta*180/Pi
gammadeg=gamma*180/Pi

!write(*,*) 'exiting EllipsoidsAAtoPolar'
END SUBROUTINE EllipsoidsAAtoPolar

      SUBROUTINE AAtoEuler(px,py,pz,phi,theta,chi)

      USE COMMONS
      USE key

      IMPLICIT NONE
      
      INTEGER          :: I, J, K
      DOUBLE PRECISION :: R(3), P(3), A(3), RM(3,3),PI
      DOUBLE PRECISION, intent(out) :: PHI, THETA, CHI, px, py, pz
 
      PI = 4.D0 * DATAN(1.D0)

            P(1) = px
            P(2) = py
            P(3) = pz

            CALL ELLIPSOIDROTATION(P, RM)

            PHI   = DATAN2(RM(2,3),RM(1,3)) 
            IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI
            
            THETA = DACOS(RM(3,3))

            CHI   = - DATAN2(RM(3,2),RM(3,1))     
            IF (CHI <= 0.D0) CHI = CHI + 2.D0*PI

            PHI    = PHI*180.D0/PI
            THETA  = THETA*180.D0/PI
            CHI    = CHI*180.D0/PI
           
!            WRITE(3,'(a5,2x,3f20.10,2x,a8,6f20.10)')                 & 
!                 'O', R(1), R(2), R(3),                              &
!                 'ellipse', 2.D0*A(1), 2.D0*A(2), 2.D0*A(3), PHI, THETA, CHI 


      END SUBROUTINE AAtoEuler

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ELLIPSOIDROTATION (P, ROTMAT)
! gives back the three Euler angles needed for displaying an ellipsoid in Xmakemol
      IMPLICIT NONE

      INTEGER          :: K1, K2
      DOUBLE PRECISION :: E(3,3), I3(3,3), ROTMAT(3,3), P(3), THETA, THETA2

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0
      I3(2,2) = 1.D0
      I3(3,3) = 1.D0

      THETA   = DSQRT(DOT_PRODUCT(P,P))

      IF (THETA == 0.D0) THEN

         ROTMAT = I3

      ELSE

         THETA2  = THETA * THETA
         E(:,:)  = 0.D0
         E(1,2)  = -P(3)
         E(1,3)  =  P(2)
         E(2,3)  = -P(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         E       = E/THETA

         ROTMAT = I3 + (1.D0-COS(THETA))*MATMUL(E,E) + E*SIN(THETA)

      ENDIF

      END SUBROUTINE ELLIPSOIDROTATION

SUBROUTINE PYREALIGNXYZ
!use key
use commons, only : natoms
!use pymodule
implicit none

integer nlines,nframes,J1,J2,J3,newnatoms,realnatoms
character(len=4) atomlabels(natoms/2)
character(len=30) dummyline
double precision COORDSREF(3*NATOMS),COORDSNEW(3*NATOMS),d,dist2,RMAT(3,3),COORDSDUMMY(3*NATOMS)

 realnatoms=natoms/2
 OPEN(unit=133,file='path.xyz',status='old') 
 OPEN(unit=134,file='newpath.xyz',status='unknown')

 CALL DETERMINELINES(133,nlines)
 WRITE(*,*) 'number of lines in the path.xyz file determined as ', nlines
 nframes=nlines/(natoms/2 + 2)
 WRITE(*,*) 'number of frames in the path.xyz file determined as ', nframes

! first frame
   READ(133,*) newnatoms
   READ(133,*) dummyline
   IF(newnatoms/=realnatoms) THEN
        WRITE(*,'(A)') 'ERROR - number of atoms determined from OPTIM does not match '// &
&                        ' with the number of atoms determined from the path.xyz file!!!'
        STOP
   END IF
   WRITE(134,*) REALNATOMS
   WRITE(134,*) dummyline
   DO J2=1,REALNATOMS
        READ(133,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSREF(3*J2-2),COORDSREF(3*J2-1),COORDSREF(3*J2), &
        & COORDSREF(3*(J2+REALNATOMS)-2),COORDSREF(3*(J2+REALNATOMS)-1),COORDSREF(3*(J2+REALNATOMS))
        WRITE(134,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSREF(3*J2-2),COORDSREF(3*J2-1),COORDSREF(3*J2), &
        & COORDSREF(3*(J2+REALNATOMS)-2),COORDSREF(3*(J2+REALNATOMS)-1),COORDSREF(3*(J2+REALNATOMS))
   END DO


 DO J1=1,nframes-1
   READ(133,*) newnatoms
   READ(133,*) dummyline
   IF(newnatoms/=realnatoms) THEN
        WRITE(*,'(A)') 'ERROR - number of atoms determined from OPTIM does not match '// &
&                        ' with the number of atoms determined from the path.xyz file!!!'
        STOP
   END IF
   DO J2=1,REALNATOMS
        READ(133,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSNEW(3*J2-2),COORDSNEW(3*J2-1),COORDSNEW(3*J2), &
        & COORDSNEW(3*(J2+REALNATOMS)-2),COORDSNEW(3*(J2+REALNATOMS)-1),COORDSNEW(3*(J2+REALNATOMS))
   END DO
   COORDSDUMMY(:)=COORDSNEW(:)
!   COORDSDUMMY(50)=COORDSNEW(50)*10.0D0
WRITE(*,*) "sn402: Changed PYREALIGNXYZ so that it calls ALIGN_DECIDE instead of MINPERMDIST"
WRITE(*,*) "I haven't tested this change and suspect that it may not work (especially if the coordinates are non-cartesian)." 
WRITE(*,*) "Please set FASTOVERLAP and check carefully that this part of the code is working as you expect, then remove these messages!"
   CALL ALIGN_DECIDE(COORDSREF,COORDSDUMMY,NATOMS,.false.,1.0,1.0,1.0,.false.,.false.,D,DIST2,.false.,RMAT)
   D=D**2 ! since ALIGN_DECIDE now returns the distance
   COORDSNEW(:)=COORDSDUMMY(:)
   WRITE(*,*) 'writing frame ', J1+1
   WRITE(134,*) REALNATOMS
   WRITE(134,*) dummyline
   DO J2=1,REALNATOMS
     WRITE(134,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSNEW(3*J2-2),COORDSNEW(3*J2-1),COORDSNEW(3*J2), &
        & COORDSNEW(3*(J2+REALNATOMS)-2),COORDSNEW(3*(J2+REALNATOMS)-1),COORDSNEW(3*(J2+REALNATOMS))
   END DO
   COORDSREF(:)=COORDSNEW(:)
 END DO


END SUBROUTINE PYREALIGNXYZ

SUBROUTINE PYRANDOMSWAP(COORDSREF,COORDSDUMMY,D)
use commons, only : natoms
use key, only : pybinarytype1
implicit none

integer :: NP, RANDOM1, RANDOM2, J1, J2, J3
double precision :: COORDSSTORE(3,4),DPRAND,COORDSREF(3*NATOMS),COORDSDUMMY(3*NATOMS),D
double precision :: rmat(3,3),BESTCOORDS(3*NATOMS),DBEST

DBEST=1.0D10
BESTCOORDS(:)=COORDSDUMMY(:)
DO J3=1,50

DO J2=1,1000000
COORDSDUMMY(:)=BESTCOORDS(:)
! first swap 
        DO
          RANDOM1=INT(DPRAND()*100)
          IF(RANDOM1<=PYBINARYTYPE1.AND.RANDOM1>=1) EXIT
        END DO
        DO
          RANDOM2=INT(DPRAND()*100)
          IF(RANDOM2<=PYBINARYTYPE1.AND.RANDOM2>=1) EXIT
        END DO

        !RANDOM1=PYSWAP(1)
        !RANDOM2=PYSWAP(2)
        COORDSSTORE(:,:)=0.0D0
        ! select a body from the two types
        DO J1=1,3
         COORDSSTORE(J1,1)=COORDSDUMMY(3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,2)=COORDSDUMMY(3*(RANDOM2-1)+J1)
         COORDSSTORE(J1,3)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,4)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)
        END DO

        ! swap coordinates
        !WRITE(MYUNIT,*) COORDSSTORE(:,1)
        !WRITE(MYUNIT,*) COORDSSTORE(:,3)
        !WRITE(MYUNIT,*) COORDSSTORE(:,2)
        !WRITE(MYUNIT,*) COORDSSTORE(:,4)
        DO J1=1,3
         COORDSDUMMY(3*(RANDOM1-1)+J1)=COORDSSTORE(J1,2)
         COORDSDUMMY(3*(RANDOM2-1)+J1)=COORDSSTORE(J1,1)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)=COORDSSTORE(J1,4)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)=COORDSSTORE(J1,3)
        END DO

        DO
          RANDOM1=INT(DPRAND()*100)
          IF(RANDOM1>PYBINARYTYPE1.AND.RANDOM1<=NATOMS/2) EXIT
        END DO
        DO
          RANDOM2=INT(DPRAND()*100)
          IF(RANDOM2>PYBINARYTYPE1.AND.RANDOM2<=NATOMS/2) EXIT
        END DO

        !RANDOM1=PYSWAP(1)
        !RANDOM2=PYSWAP(2)
        COORDSSTORE(:,:)=0.0D0
        ! select a body from the two types
        DO J1=1,3
         COORDSSTORE(J1,1)=COORDSDUMMY(3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,2)=COORDSDUMMY(3*(RANDOM2-1)+J1)
         COORDSSTORE(J1,3)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,4)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)
        END DO

        ! swap coordinates
        !WRITE(MYUNIT,*) COORDSSTORE(:,1)
        !WRITE(MYUNIT,*) COORDSSTORE(:,3)
        !WRITE(MYUNIT,*) COORDSSTORE(:,2)
        !WRITE(MYUNIT,*) COORDSSTORE(:,4)
        DO J1=1,3
         COORDSDUMMY(3*(RANDOM1-1)+J1)=COORDSSTORE(J1,2)
         COORDSDUMMY(3*(RANDOM2-1)+J1)=COORDSSTORE(J1,1)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)=COORDSSTORE(J1,4)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)=COORDSSTORE(J1,3)
        END DO

   IF(MOD(J2,1000)==0)   then
WRITE(*,*) "sn402: Changed PYRANDOMSWAP so that it calls ALIGN_DECIDE instead of MINPERMDIST"
WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
WRITE(*,*) "Please set FASTOVERLAP and check carefully that this part of the code is working as you expect, then remove these messages!" 
   CALL ALIGN_DECIDE(COORDSREF,COORDSDUMMY,NATOMS,.false.,1.0,1.0,1.0,.false.,.false.,D,D,.false.,RMAT)
   D=D**2
        IF(D<DBEST) THEN
          DBEST=D
          BESTCOORDS(:)=COORDSDUMMY(:)
!        ELSE IF (D-DBEST<0.1D0.AND.DPRAND()<0.5D0) THEN
!          DBEST=D
!          BESTCOORDS(:)=COORDSDUMMY(:)
        END IF
   END IF
END DO
WRITE(*,*) 'after 10000 swaps, DBEST= ', DBEST

END DO

END SUBROUTINE PYRANDOMSWAP

SUBROUTINE DETERMINELINES(nunit,nlines)
implicit none
integer nunit, nlines, iostatus
character(len=10) check

REWIND(nunit)

nlines=0
do
  IF(iostatus<0) EXIT
  nlines = nlines + 1
  READ(nunit,*,iostat=iostatus) check
!  write(*,*) check,nunit
end do
  nlines = nlines - 1
  REWIND(nunit)
RETURN


END SUBROUTINE DETERMINELINES

      SUBROUTINE TAKESTEPELPSD (COORDS)

!     THIS ROUTINE TAKES STEP FOR SINGLE-SITE ELLIPSOIDAL BODIES ENSURING NO OVERLAP
      USE commons, only : natoms
      USE key
      use pymodule, only : angle,angle2,twopi

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET
      INTEGER          :: J1, J2, J3, J4, J5, J6
!      INTEGER          :: PTINDX, I, J, K
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2, SAVECOORDS(3*NATOMS), COORDS(3*NATOMS), COORDSO(3*NATOMS)
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2, VAT(3*NATOMS),VATO(3*NATOMS)
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP, STEP, ASTEP, OSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI, BOXLX, BOXLY, BOXLZ
      DOUBLE PRECISION :: AE(3,3), BE(3,3) 
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3)
      DOUBLE PRECISION :: P1(3), P2(3), ABSRIJ, RCUT, XMIN, FMIN, ECFVAL
! sf344 additions
      LOGICAL          :: DISSOCIATED(NATOMS/2),DISSOC,TMOVE,OMOVE,LJCAPSIDT
      INTEGER          :: CLOSESTATOMINDEX, K1, K2
      DOUBLE PRECISION :: MINDISTANCE, x1, x2, y1, y2, z1, z2, RCUTARRAY(NATOMS/2), RVEC(NATOMS/2,NATOMS/2)
!      DOUBLE PRECISION, ALLOCATABLE :: RVEC(:)

!      WRITE(*,*) 'NATOMS in takestepelpsd ', NATOMS
!      IF(ALLOCATED(RVEC)) DEALLOCATE(RVEC)

      PI         = 4.D0*ATAN(1.0D0)
      IF (GBT .OR. GBDT) THEN

         RCUT       = GBEPSNOT*MAX(GBKAPPA,1.D0)

      ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
        DO J1=1,NATOMS/2
         RCUTARRAY(J1) = 2.0D0 * MAXVAL(PYA1bin(J1,:)) 
        END DO
         RCUT = MAXVAL(RCUTARRAY)
      ELSE IF (LJCAPSIDT) THEN
        PYA1bin(:,:)=0.5D0
        DO J1=1,NATOMS/2
         RCUTARRAY(J1) = 2.0D0
        END DO
         RCUT = 2.0D0
      ENDIF

      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS

                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)

      DO J1 = 1,REALNATOMS

         J2     = 3*J1
        IF (PARAMONOVPBCX) THEN
!ensure x component of particle 1 vector is within BoxLx/2 of zero.
!if it isn't then subtract integer number of boxlx's such that it is.
                COORDS(J2-2)=COORDS(J2-2)-BOXLX*NINT(COORDS(J2-2)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
!ensure y component of particle 1 vector is within BoxLy/2 of zero.
!if it isn't then subtract integer number of boxly's such that it is.
                COORDS(J2-1)=COORDS(J2-1)-BOXLY*NINT(COORDS(J2-1)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
!ensure z component of particle 1 vector is within BoxLz/2 of zero.
!if it isn't then subtract integer number of boxlz's such that it is.
                COORDS(J2)=COORDS(J2)-BOXLZ*NINT(COORDS(J2)/BOXLZ)
        ENDIF

         DUMMY2 = COORDS(J2-2)**2 + COORDS(J2-1)**2 + COORDS(J2)**2

         IF (DUMMY2 .GT. RADIUS) THEN
! bring back the molecule within the radius
                COORDS(J2-2)=COORDS(J2-2)-SQRT(RADIUS)*NINT(COORDS(J2-2)/SQRT(RADIUS))
                COORDS(J2-1)=COORDS(J2-1)-SQRT(RADIUS)*NINT(COORDS(J2-1)/SQRT(RADIUS))
                COORDS(J2)=COORDS(J2)-SQRT(RADIUS)*NINT(COORDS(J2)/SQRT(RADIUS))
    WRITE(MYUNIT,'(A,2F20.10)') 'initial coordinate outside container -- bringing molecule back within the container radius'

!            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2), &
!                                       COORDS(J2-1), COORDS(J2)
!            PRINT*, 'initial coordinate outside container -- increase container radius'
!            STOP
         END IF

      END DO

      DO J1 = 1,3*NATOMS
         COORDSO(J1) = COORDS(J1)
      END DO

      DO J1 = 1,NATOMS/2
         VATO(J1) = VAT(J1)
      END DO
!     FIND THE CENTRE OF MASS

      XMASS = 0.0D0
      YMASS = 0.0D0
      ZMASS = 0.0D0

      DO J1 = 1,NATOMS/2

         XMASS = XMASS + COORDS(3*(J1-1)+1)
         YMASS = YMASS + COORDS(3*(J1-1)+2)
         ZMASS = ZMASS + COORDS(3*(J1-1)+3)

      ENDDO

      XMASS = XMASS/(REALNATOMS)
      YMASS = YMASS/(REALNATOMS)
      ZMASS = ZMASS/(REALNATOMS)

!     Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!     and the pair energy of the most tightly bound atom, VMIN. An angular step is
!     taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!     DMAX (or CMMAX from CM of the cluster).
      JMAX  = 1
      DMAX  =  1.0D0
      VMAX  = -1.0D3
      VMAX2 = -1.0D3
      VMIN  =  0.0D0
      CMMAX =  1.0D0

      DO J1 = 1, REALNATOMS

         J2 = 3*J1
         DIST(J1)   = DSQRT( COORDS(J2-2)**2 + COORDS(J2-1)**2 + COORDS(J2)**2)
         CMDIST(J1) = DSQRT((COORDS(J2-2)-XMASS)**2+(COORDS(J2-1)-YMASS)**2+(COORDS(J2)-ZMASS)**2)
         IF (CMDIST(J1) .GT. CMMAX) CMMAX = CMDIST(J1)
         IF (DIST(J1) .GT. DMAX) DMAX = DIST(J1)
         IF (VAT(J1) .GT. VMAX) THEN
            VMAX = VAT(J1)
            JMAX = J1
         ELSE IF ((VAT(J1).LT. VMAX) .AND. (VAT(J1) .GT. VMAX2)) THEN
              VMAX2 = VAT(J1)
              JMAX2 = J1
         ENDIF
         IF (VAT(J1) .LT. VMIN) VMIN = VAT(J1)

      ENDDO

      IF (VAT(JMAX) > (ASTEP*VMIN) .AND. (.NOT.NORESET)) THEN

         J2 = 3*JMAX
         THETA           = DPRAND()*PI
         PHI             = DPRAND()*PI*2.0D0
         COORDS(J2-2) = XMASS + (CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
         COORDS(J2-1) = YMASS + (CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
         COORDS(J2)   = ZMASS + (CMMAX+1.0D0)*DCOS(THETA)
         DUMMY           = COORDS(J2-2)**2 + COORDS(J2-1)**2 + COORDS(J2)**2

         IF (DUMMY > RADIUS) THEN
!          RADIUS=1.0
            DUMMY           = DSQRT(RADIUS*0.99D0/DUMMY)
            COORDS(J2-2) = COORDS(J2-2)*DUMMY
            COORDS(J2-1) = COORDS(J2-1)*DUMMY
            COORDS(J2)   = COORDS(J2)*DUMMY
         END IF

      ENDIF

      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3

!     CHECK FOR OVERLAP

         OVRLPT = .TRUE.
!        WRITE(*,*) 'in here'



95       DO WHILE (OVRLPT)
            LOCALSTEP = 0.0D0
            IF (TMOVE) LOCALSTEP = STEP

!            IF(FROZEN(J1))  LOCALSTEP = 0.0D0

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-2) = COORDS(J3-2) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-1) = COORDS(J3-1) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3)   = COORDS(J3) + LOCALSTEP*RANDOM

            LOCALSTEP = 0.0D0
            IF (OMOVE) LOCALSTEP = OSTEP 

            IF(PYBINARYT) LOCALSTEP = 0.0D0

            IF(FROZEN(J1)) LOCALSTEP = 0.0D0

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-2) = COORDS(J5-2) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-1) = COORDS(J5-1) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5)   = COORDS(J5) + LOCALSTEP*RANDOM
!   
            OVRLPT = .FALSE.

            RI(:)  = COORDS(J3-2:J3)
            P1(:)  = COORDS(J5-2:J5)

!     ROTATION MATRIX

            CALL ELPSDRTMT (P1, PYA1bin(J1,:), AE)

            DO J2 = 1, REALNATOMS

               IF (J2 == J1) CYCLE
 
               J4    = 3*J2
               J6    = OFFSET + J4
               RJ(:) = COORDS(J4-2:J4)
               P2(:) = COORDS(J6-2:J6)

!            IF(FROZEN(J2))  LOCALSTEP = 0.0D0

!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J4-2) = COORDS(J4-2) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J4-1) = COORDS(J4-1) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J4)   = COORDS(J4) + LOCALSTEP*RANDOM

!            LOCALSTEP = 0.0D0
!            IF (OMOVE) LOCALSTEP = OSTEP 
!            IF(PYBINARYT) LOCALSTEP = 0.0D0

!            IF(FROZEN(J2)) LOCALSTEP = 0.0D0

!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J6-2) = COORDS(J6-2) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J6-1) = COORDS(J6-1) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J6)   = COORDS(J6) + LOCALSTEP*RANDOM

!     ROTATION MATRIX

               CALL ELPSDRTMT (P2, PYA1bin(J2,:), BE)

               RIJ    = RI - RJ
               ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        
               IF (ABSRIJ < RCUT) THEN 

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE, BE, RIJ, XMIN, FMIN)

                  ECFVAL = - FMIN
! allow for a slight overlap 
                  IF (ECFVAL < PYOVERLAPTHRESH) THEN
!                        WRITE(*,*) 'atoms overlapping', J1, J2, ECFVAL, ABSRIJ
                     OVRLPT = .TRUE.
                     GO TO 95
                  ENDIF

               ENDIF

            ENDDO  ! END LOOP OVER J2
         ENDDO  ! END WHILE
      ENDDO  ! END LOOP OVER J1   

                     IF (FREEZE) THEN
                        DO J2=1,NATOMS
                           IF (FROZEN(J2)) THEN
                              COORDS(3*(J2-1)+1:3*(J2-1)+3)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                           ENDIF
                        ENDDO
                     ENDIF


! now determine whether any particle has dissociated
 IF(PYGPERIODICT) THEN
  DO
  DISSOC=.FALSE.
   DO J1=1,natoms/2-1
    J2=3*J1
    x1=COORDS(J2-2)
    y1=COORDS(J2-1)
    z1=COORDS(J2)
   
    DISSOCIATED(J1)=.TRUE.
    DO K1=1,natoms/2
      K2=3*K1
      x2=COORDS(K2-2)
      y2=COORDS(K2-1)
      z2=COORDS(K2)
      RVEC(J1,K1)=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
!      WRITE(*,*) J1,K1,RVEC(J1,K1), 4*RCUT
      IF(RVEC(J1,K1)<4*RCUT) DISSOCIATED(J1)=.FALSE.
!      WRITE(*,*) 'in here, cutoff=',pcutoff,J1
    END DO
    IF(DISSOCIATED(J1)) THEN
        ! determine the closest neighbor
        MINDISTANCE = HUGE(1.0D0)
        DO K1=1,natoms/2
         IF(K1==J1) CYCLE
         IF(RVEC(J1,K1)<MINDISTANCE) THEN
           MINDISTANCE=RVEC(J1,K1)
           CLOSESTATOMINDEX=K1
         END IF
        END DO
        WRITE(MYUNIT,'(A,I6,A,F10.6,I6)') 'Atom ', J1, ' is dissociated, distance to the closest atom: ', &
        & MINDISTANCE, CLOSESTATOMINDEX
        DISSOC=.TRUE.
        ! move the dissociated atom closer to the closest atom
        COORDS(J2-2)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-2)+COORDS(J2-2))
        COORDS(J2-1)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-1)+COORDS(J2-1))
        COORDS(J2  )=0.5D0*(COORDS(3*CLOSESTATOMINDEX  )+COORDS(J2  ))
    END IF
   END DO
   IF(.NOT.DISSOC) EXIT
  END DO
        DO J1=NATOMS/2+1,NATOMS
            J2=3*J1
            angle2=dot_product(COORDS(J2-2:J2),COORDS(J2-2:J2))

            if(angle2>twopi**2) then
                angle2=sqrt(angle2)
                angle = angle2 - int(angle2/twopi)*twopi
                COORDS(J2-2:J2)=COORDS(J2-2:J2)/angle2 * angle
            end if
        END DO
  END IF
END SUBROUTINE TAKESTEPELPSD

      SUBROUTINE ELPSDRTMT (P, ESA, A)

      IMPLICIT NONE

      INTEGER          :: K
      DOUBLE PRECISION :: E(3,3), I3(3,3), ROT2(3,3), ROTMAT(3,3), P(3), THETA, THETA2 
      DOUBLE PRECISION :: AEZR1(3,3), A(3,3), ESA(3)

      I3(:,:)    = 0.D0
      AEZR1(:,:) = 0.D0
      I3(1,1)    = 1.D0
      I3(2,2)    = 1.D0
      I3(3,3)    = 1.D0
       
      THETA   = DSQRT(DOT_PRODUCT(P,P))

      IF (THETA == 0.D0) THEN

         ROTMAT = I3

      ELSE

         THETA2  = THETA * THETA
         E(:,:)  = 0.D0
         E(1,2)  = -P(3)
         E(1,3)  =  P(2)
         E(2,3)  = -P(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         E       = E/THETA

         ROTMAT = I3 + (1.D0-COS(THETA))*MATMUL(E,E) + E*SIN(THETA)

      ENDIF

      DO K = 1, 3
         AEZR1(K,K) = 1.D0/(ESA(K)*ESA(K))
      ENDDO

      A = MATMUL(ROTMAT,(MATMUL(AEZR1,(TRANSPOSE(ROTMAT)))))

      END SUBROUTINE ELPSDRTMT 
