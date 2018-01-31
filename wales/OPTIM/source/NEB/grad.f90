!   NEB module is an implementati on of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
MODULE GRADIENTS
     IMPLICIT NONE
     CONTAINS

     ! This work was described in detail in: S. A. Trygubenko and D. J. Wales, `A Doubly Nudged Elastic Band Method for Finding
     ! Transition States', J. Chem. Phys., 120, 2082-2094 (2004). Summary is available online at
     ! http://www-wales.ch.cam.ac.uk/~sat39/DNEBtests/ 
     SUBROUTINE NEBGRADIENT
          USE KEYGRAD
          USE KEYTAU
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          USE TANGENT
          USE SPFUNCTS
          USE NEBUTILS
          USE MODUNRES, ONLY: C,NRES
          USE GENRIGID, ONLY: RIGIDINIT, DEGFREEDOMS, RB_DISTANCE
          USE KEY, ONLY: UNRST, FROZEN, FREEZE, NEBK, BULKT, NEBRESEEDT, NEBRESEEDEMAX, NEBRESEEDBMAX, &
  &                      RBAAT, NREPULSIVE, NEBRESEEDINT, INTCONFRAC, INTCONSTRAINTT, &
  &                      ORDERI,ORDERJ,EPSALPHA, REPPOW,PERMDIST, TWOD, RIGIDBODY, NEBRESEEDDEL1, &
  &                      DISTREF, ADDREPT, NEBRESEEDDEL2, NEBRESEEDPOW1, NEBRESEEDPOW2, CHECKCONINT, VARIABLES
          USE COMMONS, ONLY : PARAM1, PARAM2, PARAM3, DEBUG
          IMPLICIT NONE
           
          INTEGER :: J1,J2,K,J3,J4,J5,JOUT,JIN,NSHRINK,NSHRINKATOM(NATOMS),J6,NATTACH,NATTACHS,NMAXINT,NMININT,LUNIT
          DOUBLE PRECISION :: PERPCOMP=0.1D0
          DOUBLE PRECISION :: DIHEDIST,TMPRMS,QINT(NINTS*(NIMAGE+2)),DIFFM(NINTS),DIFFP(NINTS) ! JMC
          DOUBLE PRECISION DIST
          DOUBLE PRECISION :: PI=3.141592653589793D0
          DOUBLE PRECISION TEMP1(NOPT), TEMP2(NOPT), GGGSAVE(NOPT*(NIMAGE+2)), DUMMY, DUMMY2, SPRING(NOPT*(NIMAGE+2)), &
  &                        TRUEPERP(NOPT*(NIMAGE+2)), IMCOORDS(NOPT), RSITE(NOPT), RMAT(3,3), D, DIST2, LX(3), LV(3), &
  &                        REPGRAD(NOPT), MAGREP, MEAND, VS(3), VF(3), MP1S(3), MP2S(3), MP1F(3), MP2F(3), DS, DF, D1, D2, DI, DJ, &
  &                        MP2SS(3), MP2FF(3), MP1SS(3), MP1FF(3), DSS, DFF, DISTA, DISTB

          ! sn402: rigid body variables
          DOUBLE PRECISION :: DISTL, DISTR, GRADL1(DEGFREEDOMS), GRADR1(DEGFREEDOMS), GRADL2(DEGFREEDOMS), GRADR2(DEGFREEDOMS)

          INTEGER, ALLOCATABLE :: IREPTEMP(:)
          LOGICAL SHRINKT(NATOMS)
          DOUBLE PRECISION EEETMP(NIMAGE+2), MYGTMP(NOPT*(NIMAGE+2)), ETMP, RMSTMP, GETUNIT

          GSPR(:) = 0.0D0  ! sn402: added
          GGG(:) = 0.0D0   ! sn402: added
          SSS(:) = 0.0D0   ! sn402: added
          
          CALL TRUEPOTEG

!         IF (INTCONSTRAINTT.AND.(INTCONFRAC.LT.1.0D0)) THEN
!            EEETMP(1:NIMAGE+2)=EEE(1:NIMAGE+2)
!            MYGTMP(1:NOPT*(NIMAGE+2))=GGG(1:NOPT*(NIMAGE+2))
!            IF (CHECKCONINT) THEN
!               CALL CONGRAD2(NMAXINT,NMININT,ETMP,XYZ,GGG,EEE,IMGFREEZE,RMSTMP)
!            ELSE
!               CALL CONGRAD(NMAXINT,NMININT,ETMP,XYZ,GGG,EEE,IMGFREEZE,RMSTMP)
!            ENDIF
!            LUNIT=GETUNIT()
!            OPEN(LUNIT,FILE='intcon_energies',STATUS='UNKNOWN')
!            WRITE(LUNIT,'I6,2G20.10') (J1,EEE(J1),EEETMP(J1),J1=1,NIMAGE+2)
!            CLOSE(LUNIT)
!            GGG(1:NOPT*(NIMAGE+2))=INTCONFRAC*MYGTMP(1:NOPT*(NIMAGE+2))+(1.0D0-INTCONFRAC)*GGG(1:NOPT*(NIMAGE+2))
!            EEE(1:NIMAGE+2)=INTCONFRAC*EEETMP(1:NIMAGE+2)+(1.0D0-INTCONFRAC)*EEE(1:NIMAGE+2)
!         ENDIF

          CALL MEPTANGENT

          IF (BADTAU) RETURN

!
! RBAAT actually works best with standard DNEB procedure. Hmm.
!

!
! For RBAAT rigid bodies freeze one set of coordinates by setting gradient
! components to zero. We are then just minimising with respect to the other set
! for each image.
!
!         IF (RBAAT) THEN
!            DO J1=2,NIMAGE+1
!               GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+(NOPT/2))=0.0D0 ! freeze centres of mass
!               GGG(NOPT*(J1-1)+(NOPT/2)+1:NOPT*J1)=0.0D0 ! freeze angles
!            ENDDO
!            GOTO 555
!         ENDIF
!          
! GGGSAVE contains the true gradient on images
!
          GGGSAVE(1:NOPT*(NIMAGE+2))=GGG(1:NOPT*(NIMAGE+2))
          ! GGG contains the perpendicular component of the true gradient after this block
          DO J1=2,NIMAGE+1  !  MEP-perpendicular component of true gradient only: G = G - G|| = G - (G,TAU)*TAU
               IF (UNRST) THEN
                  GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS) = GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS) - &
               &    DOT_PRODUCT(GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS),TANVEC(:NINTS,J1-1))*TANVEC(:NINTS,J1-1)
               ELSEIF (RIGIDINIT) THEN
                  GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) - &
               &    DOT_PRODUCT(GGG(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS),TANVEC(:DEGFREEDOMS,J1-1))*TANVEC(:DEGFREEDOMS,J1-1)
               ELSE
                  GGG(NOPT*(J1-1)+1:NOPT*J1) = GGG(NOPT*(J1-1)+1:NOPT*J1) - &
               &    DOT_PRODUCT(GGG(NOPT*(J1-1)+1:NOPT*J1),TANVEC(:,J1-1))*TANVEC(:,J1-1)
               ENDIF
               TRUEPERP(1:NOPT*(NIMAGE+2))=GGG(1:NOPT*(NIMAGE+2))
          ENDDO
          IF (UNRST) THEN
             TMPRMS=0.0D0
             DO J1=1,NIMAGE
                TMPRMS = TMPRMS + DOT_PRODUCT(G(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS),G(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS))
             ENDDO
             RMS=SQRT( TMPRMS/(NIMAGE*NINTS) ) 
          ELSE
             RMS=SQRT( DOT_PRODUCT(G,G)/(NIMAGE*NOPT) ) 
          ENDIF

          IF (BULKT.AND.(GRADTYPE.NE.'dneb')) THEN
             PRINT '(A)',' nebgradient> ERROR - minimum image distances only coded for GRADTYPE dneb'
             STOP
          ENDIF

          ! spring gradient on images
          SELECT CASE (GRADTYPE) ! THIS CALCULATES SPRING GRADIENT

             CASE("jnew") ! springs in new implementation of neb by Jonsson
                IF(RIGIDINIT) THEN  ! Warning: not tested
                   WRITE(*,*) "Warning: using the new rigid body distance measure with NEB gradient type ", GRADTYPE
                   WRITE(*,*) "This combination has not been tested"
                   DO J1=1,NIMAGE
                      CALL RB_DISTANCE(DISTR, XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), &
                      &                   GRADL1, GRADR1, .FALSE.)
                      CALL RB_DISTANCE(DISTL, XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), &
                      &                   GRADL2, GRADR2, .FALSE.)
                      GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = - NEBK*(DISTR - DISTL)*TANVEC(:DEGFREEDOMS,J1)  ! Not sure what to do about TANVEC
                   ENDDO
                ELSE
                   DO J1=1,NIMAGE
                      GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                      & - NEBK*( SQRT(SUM( ( XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) &
                      &         - SQRT(SUM( ( XYZ(NOPT*(J1-1)+1:NOPT*J1)     - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) )*TANVEC(:,J1)
                   ENDDO
                ENDIF

             CASE("spr") ! full springs:  Gspr = k [2x(j) -x(j-1) -x(j+1)]
                  IF(RIGIDINIT) THEN   ! Not tested
                   WRITE(*,*) "Warning: using the new rigid body distance measure with NEB gradient type ", GRADTYPE
                   WRITE(*,*) "This combination has not been tested"
                     DO J1=1,NIMAGE
                        ! GRADL1 is the distance derivative with respect to the first coordinate argument (X_j1)
                        CALL RB_DISTANCE(DISTR, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), &
                      &                   GRADL1, GRADR1, .TRUE.)
                        CALL RB_DISTANCE(DISTL, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), &
                      &                   GRADL2, GRADR2, .TRUE.)
                        GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = 0.5D0*NEBK*(GRADL1(:)+GRADL2(:))
                     ENDDO
                  ELSE
                     DO J1=1,NIMAGE
                        GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                        & NEBK*( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )
                     ENDDO
                  ENDIF

             CASE("jold") ! springs MEP-parallel component only (original formulation): G = (Gspr,Tau)*Tau
                IF(RIGIDINIT) THEN  ! Not tested
                   WRITE(*,*) "Warning: using the new rigid body distance measure with NEB gradient type ", GRADTYPE
                   WRITE(*,*) "This combination has not been tested"
                   DO J1=1,NIMAGE
                      ! GRADL1 is the distance derivative with respect to the first coordinate argument (X_j1)
                      CALL RB_DISTANCE(DISTR, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), &
                    &                   GRADL1, GRADR1, .TRUE.)
                      CALL RB_DISTANCE(DISTL, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), &
                    &                   GRADL2, GRADR2, .TRUE.)
                      GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = 0.5D0*NEBK*DOT_PRODUCT((GRADL1(:)+GRADL2(:)), &
                    &         TANVEC(:DEGFREEDOMS,J1))*TANVEC(:DEGFREEDOMS,J1)
                   ENDDO
                ELSE
                   DO J1=1,NIMAGE
                      GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                   &  NEBK*DOT_PRODUCT( (2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2))), &
                   &  TANVEC(:,J1) )*TANVEC(:,J1)
                   ENDDO
                ENDIF

             CASE("nper")! springs MEP-parallel component (original formulation) + n% of perpendicular component:
                IF(RIGIDINIT) THEN  ! Not tested
                   WRITE(*,*) "Warning: using the new rigid body distance measure with NEB gradient type ", GRADTYPE
                   WRITE(*,*) "This combination has not been tested"
                   DO J1=1,NIMAGE
                      ! GRADL1 is the distance derivative with respect to the first coordinate argument (X_j1)
                      CALL RB_DISTANCE(DISTR, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), &
                    &                   GRADL1, GRADR1, .TRUE.)
                      CALL RB_DISTANCE(DISTL, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), &
                    &                   GRADL2, GRADR2, .TRUE.)
                      GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = 0.5D0*NEBK*(GRADL1(:)+GRADL2(:)) - &
                    &      0.5D0*NEBK*DOT_PRODUCT((GRADL1(:)+GRADL2(:)),TANVEC(:DEGFREEDOMS,J1))*TANVEC(:DEGFREEDOMS,J1)*PERPCOMP

                      GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) + &
                    &      0.5*NEBK*DOT_PRODUCT((GRADL1(:)+GRADL2(:)),TANVEC(:DEGFREEDOMS,J1))*TANVEC(:DEGFREEDOMS,J1)
                   ENDDO
                ELSE
                   DO J1=1,NIMAGE   !GJ = GJ + G|_*0.1 + GSPR|| == GJ + [ GSPR - (GSPR,TAU)*TAU ]*N + ...
                      GSPR(NOPT*(J1-1)+1:NOPT*J1) = ( &
                   &   NEBK*           ( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )&
                   & - NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                   &   TANVEC(:,J1) )*TANVEC(:,J1) )*PERPCOMP
                 
                      GSPR(NOPT*(J1-1)+1:NOPT*J1) = GSPR(NOPT*(J1-1)+1:NOPT*J1) &  !... (GSPR,TAU)*TAU
                   &  + NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                   &    TANVEC(:,J1) )*TANVEC(:,J1)
                   ENDDO
                ENDIF

             CASE("nper2")! same as above but optimised for efficiency?
                IF(RIGIDINIT) THEN  ! Not tested
                   WRITE(*,*) "Warning: using the new rigid body distance measure with NEB gradient type ", GRADTYPE
                   WRITE(*,*) "This combination has not been tested"
                   DO J1=1,NIMAGE
                      ! GRADL1 is the distance derivative with respect to the first coordinate argument (X_j1)
                      CALL RB_DISTANCE(DISTR, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), &
                    &                   GRADL1, GRADR1, .TRUE.)
                      CALL RB_DISTANCE(DISTL, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), &
                    &                   GRADL2, GRADR2, .TRUE.)
                      GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS) = 0.5D0*NEBK*(GRADL1(:)+GRADL2(:))*PERPCOMP + &
                    &     (1-PERPCOMP)*DOT_PRODUCT(0.5D0*NEBK*(GRADL1(:)+GRADL2(:)),TANVEC(:DEGFREEDOMS,J1))*TANVEC(:DEGFREEDOMS,J1)
                   ENDDO
                ELSE
                   DO J1=1,NIMAGE
                      ! calculate gs (spring gradient on image j1)
                      GSPR(NOPT*(J1-1)+1:NOPT*J1) = NEBK*(&
                            2*XYZ(NOPT*J1+1:NOPT*(J1+1))-XYZ(NOPT*(J1-1)+1:NOPT*J1)-XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) &
                                                       &)
                      ! calculate [ f*gs + (1-f)(gs,tau)tau ]
                      GSPR(NOPT*(J1-1)+1:NOPT*J1) = &
                 &    PERPCOMP*GSPR(NOPT*(J1-1)+1:NOPT*J1) &
                 &        +(1-PERPCOMP)*DOT_PRODUCT(GSPR(NOPT*(J1-1)+1:NOPT*J1),TANVEC(:,J1))*TANVEC(:,J1)
                   ENDDO
                ENDIF
              CASE("dneb")  ! this appears to be the default
                 DO J1=1,NIMAGE
                 ! calculate SSS = ~g perp in DNEB paper, the perpendicular part of the spring gradient vector
!
!  Here we are calculating spring gradient vector - (spring gradient vector . tangent vector) tangent vector
!  This gives the spring gradient perpendicular component, which is put in SSS.
!
!  Changed to use dynamically adjusted NEBK DJW 21/10/08
!  Only done for default CASE("dneb")
!
                    IF(RIGIDINIT) THEN ! sn402: added
                       ! calculate SSS = ~g perp in DNEB paper, the perpendicular part of the spring gradient vector
                       ! Note, the RB_DISTANCE call takes account of PBCs, if necessary.
                       CALL RB_DISTANCE(DISTL, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), &
                    &                   GRADL1, GRADR1, .TRUE.)
                       CALL RB_DISTANCE(DISTR, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), &
                    &                   GRADL2, GRADR2, .TRUE.)

                       SSS(NOPT*J1+1:NOPT*J1+DEGFREEDOMS) = 0.5D0*NEWNEBK(J1+1)*GRADL1(:)+0.5D0*NEWNEBK(J1)*GRADL2(:) - &
                    &      0.5D0*NEWNEBK(J1+1)*DOT_PRODUCT( GRADL1(:),TANVEC(1:DEGFREEDOMS,J1) )*TANVEC(1:DEGFREEDOMS,J1) - &
                    &      0.5D0*NEWNEBK(J1)*DOT_PRODUCT( GRADL2(:),TANVEC(1:DEGFREEDOMS,J1) )*TANVEC(1:DEGFREEDOMS,J1)

                    ELSEIF (BULKT) THEN ! minimum image convention for distances!
                       DO K=1,NATOMS
                          ! TEMP1 holds the minimum-image displacement of image J1 from image J1-1
                          TEMP1(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                         -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                          TEMP1(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                         -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                          IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                         -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)
                          ! TEMP2 holds the minimum-image displacement of image J1 from image J1+1
                          TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1 ) &
   &                          -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1))/PARAM1)
                          TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2 ) &
   &                         -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2))/PARAM2)
                          IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3 ) &
   &                         -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3))/PARAM3)
                       ENDDO
                       SSS(NOPT*J1+1:NOPT*(J1+1)) = NEWNEBK(J1)*TEMP1(1:NOPT)+NEWNEBK(J1+1)*TEMP2(1:NOPT) &
   &                      -NEWNEBK(J1)  *DOT_PRODUCT(TEMP1(1:NOPT),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1) &
   &                      -NEWNEBK(J1+1)*DOT_PRODUCT(TEMP2(1:NOPT),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1)
                    ELSE
                       ! sn402: So far as I can tell, this variable (SPRING) is not used for anything. I will experiment with
                       ! commenting it out.
                       SPRING(NOPT*J1+1:NOPT*(J1+1))=  &
   &                        NEWNEBK(J1)  *( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) ) &
   &                        + NEWNEBK(J1+1)*( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )

                       SSS(NOPT*J1+1:NOPT*(J1+1)) = &
   &                       NEWNEBK(J1)  *( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) ) &
   &                       + NEWNEBK(J1+1)*( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ) &
   &      - NEWNEBK(J1)  *DOT_PRODUCT( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1)    , TANVEC(:,J1) )*TANVEC(:,J1) &
   &      - NEWNEBK(J1+1)*DOT_PRODUCT( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)), TANVEC(:,J1) )*TANVEC(:,J1)
                    ENDIF
!                 PRINT '(A,I6)','grad> SSS spring gradient perpendicular for image ',J1
!                 PRINT '(6G20.10)',SSS(1+NOPT*J1:(J1+1)*NOPT)

!
! SSS now contains the spring gradient perpendicular part.
! SPRING contains the complete spring gradient.
!
                       CALL NUDGE(J1,J1)
!                 PRINT '(A,I6)','grad> SSS  g~* for image ',J1
!                 PRINT '(6G20.10)',SSS(1+NOPT*J1:(J1+1)*NOPT)
!
! On exit from NUDGE SSS contains g~* from equation (13) of the DNEB paper.
! GGG still contains g perp, the perpendicular part of the true gradient.
! These are two of the three terms in equation (12) of the DNEB paper.
!

                 ENDDO  ! End of loop over images (J1)
!
! The call to NUDGE calculated a new SSS containing the nudged spring gradient g~*.
! We now add this to GGG. We are still missing ~g parallel.
! GGG now contains those bits of equation (12) that we have calculated thus far.
!
                 GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))
!                 PRINT '(A,I6)','grad> GGG g perp + g~* for image ',J1
!                 PRINT '(6G20.10)',GGG(1+NOPT*J1:(J1+1)*NOPT)
!
! Calculate ~g parallel and store in SSS. This is not added to GGG until after the whole SELECT block.
! In fact ~g parallel is not used as the formula expected from projection, but uses an
! alternative formulation from equation (5) of the DNEB paper.
!
                 IF(RIGIDINIT) THEN
                    DO J1=1,NIMAGE
                       CALL RB_DISTANCE(DISTR, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1+1)+1:NOPT*(J1+1)+DEGFREEDOMS), &
                    &                   GRADL1, GRADR1, .FALSE.)
                       CALL RB_DISTANCE(DISTL, XYZ(NOPT*J1+1:NOPT*J1+DEGFREEDOMS), XYZ(NOPT*(J1-1)+1:NOPT*(J1-1)+DEGFREEDOMS), &
                    &                   GRADL2, GRADR2, .FALSE.)
                       SSS(NOPT*J1+1:NOPT*J1+DEGFREEDOMS)= -0.5D0*(NEWNEBK(J1+1)*DISTR - NEWNEBK(J1)*DISTL)*TANVEC(1:DEGFREEDOMS,J1)
                    ENDDO
                 ELSE
                    DO J1=1,NIMAGE
                       IF (BULKT) THEN ! minimum image convention for distances!

                          DISTA=0.0D0; DISTB=0.0D0
                          DO J2=1,NATOMS
                             DISTA=DISTA + (XYZ(NOPT*(J1+1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1) - &
  &                            PARAM1*NINT((XYZ(NOPT*(J1+1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1))/PARAM1))**2+ &
  &                                        (XYZ(NOPT*(J1+1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2) - &
  &                            PARAM2*NINT((XYZ(NOPT*(J1+1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2))/PARAM2))**2
                IF (.NOT.TWOD) DISTA=DISTA+(XYZ(NOPT*(J1+1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3) - &
  &                            PARAM3*NINT((XYZ(NOPT*(J1+1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3))/PARAM3))**2

                             DISTB=DISTB + (XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1) - &
  &                            PARAM1*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1))/PARAM1))**2+ &
  &                                        (XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2) - &
  &                            PARAM2*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2))/PARAM2))**2
                IF (.NOT.TWOD) DISTB=DISTB+(XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3) - &
  &                            PARAM3*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3))/PARAM3))**2
                          ENDDO
                          SSS(NOPT*J1+1:NOPT*(J1+1))= - ( NEWNEBK(J1+1)*SQRT(DISTA) - NEWNEBK(J1)*SQRT(DISTB) )*TANVEC(1:NOPT,J1)
                       ELSE
                          SSS(NOPT*J1+1:NOPT*(J1+1))= &
  &               - ( NEWNEBK(J1+1)*SQRT(SUM( ( XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) &
  &                 - NEWNEBK(J1)  *SQRT(SUM( ( XYZ(NOPT*(J1-1)+1:NOPT*J1)     - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) )*TANVEC(:,J1)
                       ENDIF
                    ENDDO
                 ENDIF

              CASE("dneb2") ! same as "dneb", except spring gradient usage is more consistent;
               DO J1=1,NIMAGE ! THANKS TO DR.~DOMINIC R. ALFONSO FOR POINTING THAT OUT
                 ! calculate Gspr perp
                 SSS(NOPT*J1+1:NOPT*(J1+1)) = &
                 &   NEBK*           ( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )&
                 & - NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                 & TANVEC(:,J1) )*TANVEC(:,J1)
                 CALL NUDGE(J1,J1)
               ENDDO
!
! GGG contains the true gradient component perpendicular and SSS now contains g~*,
! so the next line gives the perpendicular component of the DNEB gradient.
!
               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))

               DO J1=1,NIMAGE ! NOW CALCULATE GSPR PARALLEL (USING OLD FORMULA)
                  IF (BULKT) THEN ! minimum image convention for distances!
                     DO K=1,NATOMS
                        TEMP1(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                        TEMP1(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)

                        TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1))/PARAM1)
                        TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3))/PARAM3)
                     ENDDO
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
   &                    NEBK*DOT_PRODUCT(TEMP1(1:NOPT)+TEMP2(1:NOPT),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1)
                  ELSE
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
                 & NEBK*DOT_PRODUCT( 2*XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) ,&
                 &   TANVEC(:,J1) )*TANVEC(:,J1)
                  ENDIF
               ENDDO
              CASE("dihen") ! DAE Distance between images calculated in dihedral angle space; for use with CHARMM and UNRES.
                            ! Use new formulation of spring force to simplify
               SSS=0.0D0    ! INITIALISE SPRING GRADIENT
               DO J1=1,NIMAGE
                 IF (UNRST) THEN
                    CALL UNRESGETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ELSE
                    CALL GETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ENDIF
                 GSPR(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS) = - NEBK * DIHEDIST * TANVEC(:NINTS,J1)
               ENDDO
              CASE("dnebu") ! DNEB for unres
               SSS=0.0D0    ! INITIALISE SPRING GRADIENT
               DO J1=1,NIMAGE+2 ! GET INTERNALS FOR ENDPOINTS AND IMAGES
                  DO J2=1,NRES
                     C(1,J2)=XYZ(6*(J2-1)+1+NOPT*(J1-1))
                     C(2,J2)=XYZ(6*(J2-1)+2+NOPT*(J1-1))
                     C(3,J2)=XYZ(6*(J2-1)+3+NOPT*(J1-1))
                     C(1,J2+NRES)=XYZ(6*(J2-1)+4+NOPT*(J1-1))
                     C(2,J2+NRES)=XYZ(6*(J2-1)+5+NOPT*(J1-1))
                     C(3,J2+NRES)=XYZ(6*(J2-1)+6+NOPT*(J1-1))
                  ENDDO
                  CALL UPDATEDC
                  CALL INT_FROM_CART(.TRUE.,.FALSE.)
                  CALL GEOM_TO_VAR(NINTS,QINT(NINTS*(J1-1)+1:NINTS*J1))
               ENDDO

               DO J1=1,NIMAGE
                 ! calculate Gspr perp
                 DIFFM(1:NINTS) = QINT(NINTS*J1+1:NINTS*(J1+1)) - QINT(NINTS*(J1-1)+1:NINTS*J1)
                 DIFFP(1:NINTS) = QINT(NINTS*J1+1:NINTS*(J1+1)) - QINT(NINTS*(J1+1)+1:NINTS*(J1+2))
                 DO J2=1,NINTS
                    IF (DIFFM(J2).LT.-PI) DIFFM(J2) = DIFFM(J2)+2.0D0*PI
                    IF (DIFFP(J2).GT.PI) DIFFP(J2) = DIFFP(J2)-2.0D0*PI
                 ENDDO

                 SSS(NOPT*J1+1:NOPT*J1+NINTS) = NEBK*( DIFFM(1:NINTS) + DIFFP(1:NINTS) ) &
                 & - NEBK*DOT_PRODUCT( ( DIFFM(1:NINTS) + DIFFP(1:NINTS) ) , TANVEC(1:NINTS,J1) )*TANVEC(1:NINTS,J1)
                 CALL NUDGE(J1,J1) ! ZERO-ORDER NUDGING
               ENDDO

               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1)) ! G PERP + G TILDE STAR

               DO J1=1,NIMAGE ! NOW CALCULATE GSPR PARALLEL AS FOR DIHEN
                 IF (UNRST) THEN
                    CALL UNRESGETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ELSE
                    CALL GETDIHEDIST(DIHEDIST,XYZ(NOPT*(J1-1)+1:NOPT*J1),XYZ(NOPT*J1+1:NOPT*(J1+1)),&
                    &XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))
                 ENDIF
                 SSS(NOPT*J1+1:NOPT*J1+NINTS) = - NEBK * DIHEDIST * TANVEC(:NINTS,J1)
                 !print *,'gspr ',sss(nopt*j1+1:nopt*j1+NINTS)
               ENDDO
              CASE("dneb3")  ! uses deviation from the average separation in the spring term instead of absolute distance
               CALL DISTANCES ! sets up image distances in vector DVEC from 1 to NIMAGE+1
               MEAND=0.0D0
               DO J1=1,NIMAGE+1
                  MEAND=MEAND+DVEC(J1)
               ENDDO
               MEAND=MEAND/(NIMAGE+1)
               PRINT '(A,F15.3)','grad> Mean image separation=',MEAND
               
               DO J1=1,NIMAGE
                 ! calculate SSS = Gspr perp  i.e. ~g perp in DNEB paper
!
!  Here we are calculating spring gradient vector - (spring gradient vector . tangent vector) tangent vector
!  This gives the spring gradient perpendicular component
!
                  IF (BULKT) THEN ! minimum image convention for distances!
                     DO K=1,NATOMS
                        TEMP1(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                        TEMP1(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)

                        TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1+1)+3*(K-1)+1))/PARAM1)
                        TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1+1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1+1)+3*(K-1)+3))/PARAM3)
                     ENDDO
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = NEWNEBK(J1)*TEMP1(1:NOPT)*(DVEC(J1)-MEAND)/DVEC(J1) &
   &                                             +NEWNEBK(J1+1)*TEMP2(1:NOPT)*(DVEC(J1+1)-MEAND)/DVEC(J1+1) &
   &                   -NEWNEBK(J1)  *DOT_PRODUCT(TEMP1(1:NOPT)*(DVEC(J1)-MEAND)/DVEC(J1),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1) &
   &                   -NEWNEBK(J1+1)*DOT_PRODUCT(TEMP2(1:NOPT)*(DVEC(J1+1)-MEAND)/DVEC(J1+1),TANVEC(1:NOPT,J1) )*TANVEC(1:NOPT,J1)
                  ELSE
                     SSS(NOPT*J1+1:NOPT*(J1+1)) = &
      &   NEWNEBK(J1)  *( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1) )*(DVEC(J1)-MEAND)/DVEC(J1) &
      & + NEWNEBK(J1+1)*( XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) )*(DVEC(J1+1)-MEAND)/DVEC(J1+1) &
      & - NEWNEBK(J1)  *DOT_PRODUCT( (XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1-1)+1:NOPT*J1))*(DVEC(J1)-MEAND)/DVEC(J1)    &
                                    , TANVEC(:,J1) )*TANVEC(:,J1) &
      & - NEWNEBK(J1+1)*DOT_PRODUCT( (XYZ(NOPT*J1+1:NOPT*(J1+1)) - XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)))*(DVEC(J1+1)-MEAND)/DVEC(J1+1) &
      &                             , TANVEC(:,J1) )*TANVEC(:,J1)
                  ENDIF
                  SPRING(NOPT*(J1-1)+1:NOPT*J1)=SSS(NOPT*(J1-1)+1:NOPT*J1)
                  CALL NUDGE(J1,J1)
               ENDDO
!
! Here we complete the DNEB gradient vector as in equation (12) of the DNEB paper.
! The call to NUDGE calculated a new SSS containing the nudged spring gradient.
!
               GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))

               DO J1=1,NIMAGE ! NOW CALCULATE GSPR PARALLEL (USING FORMULA THAT KEEPS IMAGES EQUISPACED WITH NEW TANGENT)
                 SSS(NOPT*J1+1:NOPT*(J1+1)) = &
                 & - ( NEWNEBK(J1+1)*SQRT(SUM( ( XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) &
                 &   - NEWNEBK(J1)  *SQRT(SUM( ( XYZ(NOPT*(J1-1)+1:NOPT*J1)     - XYZ(NOPT*J1+1:NOPT*(J1+1)) )**2 )) )*TANVEC(:,J1)
               ENDDO
              CASE DEFAULT
               PRINT *,'ERROR: unknown gradient type "'//trim(adjustl(GradType))//'"'; stop
          END SELECT

!----------------------------------------------------------------------------------------------------------------------------------
!
! End of large select block. Now we add on SSS to GGG, which completes equation (12) of the DNEB paper.
! At this point GG contains g perp + ~g* and we need to add on ~g parallel.
!

          GGG(NOPT+1:NOPT*(NIMAGE+1)) = GGG(NOPT+1:NOPT*(NIMAGE+1)) + SSS(NOPT+1:NOPT*(NIMAGE+1))

          IF (NEBRESEEDT.AND.(NREPULSIVE.GT.0).AND.ADDREPT) THEN
             DO J2=1,NREPULSIVE
                DO J1=2,NIMAGE+1
                   IF ((BADIMAGE(J1).AND.(REPPOW(J2).GT.0)).OR.(BADPEPTIDE(J1).AND.(REPPOW(J2).LT.0))) THEN
!
!  Add repulsive/constraint terms defined as interatomic distances.
!
                      DIST=SQRT((XYZ(NOPT*(J1-1)+3*(ORDERI(J2)-1)+1)-XYZ(NOPT*(J1-1)+3*(ORDERJ(J2)-1)+1))**2 &
  &                            +(XYZ(NOPT*(J1-1)+3*(ORDERI(J2)-1)+2)-XYZ(NOPT*(J1-1)+3*(ORDERJ(J2)-1)+2))**2 &
  &                            +(XYZ(NOPT*(J1-1)+3*(ORDERI(J2)-1)+3)-XYZ(NOPT*(J1-1)+3*(ORDERJ(J2)-1)+3))**2)

                      REPGRAD(1:NOPT)=0.0D0

                      DUMMY=DIST-DISTREF(J2)  !  DJW BUG this J2 was J1!
                      IF (DUMMY.EQ.0.0D0) DUMMY=1.0D-10

                      REPGRAD(3*(ORDERI(J2)-1)+1:3*(ORDERI(J2)-1)+3)= &
     &                 -REPPOW(J2)*EPSALPHA(J2) &
     &                  *(XYZ((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3) &
     &                   -XYZ((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3)) &
     &                         /(DIST*DUMMY**(REPPOW(J2)+1))

                      REPGRAD(3*(ORDERJ(J2)-1)+1:3*(ORDERJ(J2)-1)+3)= &
     &                 -REPPOW(J2)*EPSALPHA(J2) &
     &                  *(XYZ((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3) &
     &                   -XYZ((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3)) &
     &                         /(DIST*DUMMY**(REPPOW(J2)+1))

                      GGG((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3)= &
    &                 GGG((J1-1)*NOPT+3*(ORDERI(J2)-1)+1:(J1-1)*NOPT+3*(ORDERI(J2)-1)+3)  &
    &                   +REPGRAD(3*(ORDERI(J2)-1)+1:3*(ORDERI(J2)-1)+3)
  
                      GGG((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3)= &
    &                 GGG((J1-1)*NOPT+3*(ORDERJ(J2)-1)+1:(J1-1)*NOPT+3*(ORDERJ(J2)-1)+3)  &
    &                   +REPGRAD(3*(ORDERJ(J2)-1)+1:3*(ORDERJ(J2)-1)+3)

                   ENDIF
                ENDDO
             ENDDO
          ENDIF
!
! Set gradients on frozen atoms to zero.
!
555       IF (FREEZE) THEN
             DO J1=2,NIMAGE+1  
                DO J2=1,NATOMS
                   IF (FROZEN(J2)) THEN
                      IF (VARIABLES) THEN
                         GGG(NOPT*(J1-1)+J2)=0.0D0
                      ELSE
                         GGG(NOPT*(J1-1)+3*(J2-1)+1)=0.0D0
                         GGG(NOPT*(J1-1)+3*(J2-1)+2)=0.0D0
                         GGG(NOPT*(J1-1)+3*(J2-1)+3)=0.0D0
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          RMS=0.0D0
          DO J1=2,NIMAGE+1
             DO J2=1,NOPT
                RMS=RMS+GGG(NOPT*(J1-1)+J2)**2
             ENDDO
          ENDDO
          RMS=SQRT(RMS/(NOPT*NIMAGE))

          ETOTAL=SUM(EIMAGE)

          IF (ORT) THEN ! REMOVE OVERALL ROTATION/TRANSLATION BY FREEZING 6 DEGREES OF FREEDOM
               G(1::NOPT)=0.0D0
               G(2::NOPT)=0.0D0
               G(3::NOPT)=0.0D0
               G(4::NOPT)=0.0D0
               G(5::NOPT)=0.0D0
               G(7::NOPT)=0.0D0
          ENDIF

     END SUBROUTINE NEBGRADIENT
          
     SUBROUTINE NUDGE(JS,JT) ! SUBROUTINE DOES HIGHER-ORDER NUDGINGS; IT MODIFIES GSPR, EATS BOTH GSPR AND G; PREREQ.
          ! - Gspr and G are calculated; js - image on which gsper is nudged by gtper on image jt
          USE NEBDATA
          USE KEY, ONLY: UNRST
          USE GENRIGID, ONLY: RIGIDINIT, DEGFREEDOMS
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: JS,JT
          DOUBLE PRECISION DUMMY
!
! On entry, SSS contains the perpendicular spring gradient component
! GGG must contain a non-unit vector with the perpendicular component of the true potential
! On exit, SSS contains g~* from equation (13) of the DNEB paper
!
          IF (UNRST) THEN
             IF (DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*JS+NINTS),GGG(NOPT*JT+1:NOPT*JT+NINTS)) < 0.0D0) THEN
                  ! calculate the beast Gspr perp - (Gspr perp, Gt perp)* Gt perp/(Gt perp, Gt perp)
                  SSS(NOPT*JS+1:NOPT*JS+NINTS) = SSS(NOPT*JS+1:NOPT*JS+NINTS) - &
                  & DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*JS+NINTS),GGG(NOPT*JT+1:NOPT*JT+NINTS))*GGG(NOPT*JT+1:NOPT*JT+NINTS)/ &
                  & DOT_PRODUCT(GGG(NOPT*JT+1:NOPT*JT+NINTS),GGG(NOPT*JT+1:NOPT*JT+NINTS))
             ENDIF
          ELSE IF(RIGIDINIT) THEN
!             sn402: I don't know why the next line is here. So I've commented it out
!             IF (DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*JS+DEGFREEDOMS),GGG(NOPT*JT+1:NOPT*JT+DEGFREEDOMS)) < 0.0D0) THEN
                  ! calculate the beast Gspr perp - (Gspr perp, Gt perp)* Gt perp/(Gt perp, Gt perp)
                DUMMY=DOT_PRODUCT(GGG(NOPT*JT+1:NOPT*JT+DEGFREEDOMS),GGG(NOPT*JT+1:NOPT*JT+DEGFREEDOMS))
                IF (DUMMY.EQ.0.0D0) DUMMY=1.0D-100
                  SSS(NOPT*JS+1:NOPT*JS+DEGFREEDOMS) = SSS(NOPT*JS+1:NOPT*JS+DEGFREEDOMS) - &
                  & DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*JS+DEGFREEDOMS),GGG(NOPT*JT+1:NOPT*JT+DEGFREEDOMS))*  &
                  & GGG(NOPT*JT+1:NOPT*JT+DEGFREEDOMS) / DUMMY
!             ENDIF
          ELSE
             IF (DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*(JS+1)),GGG(NOPT*JT+1:NOPT*(JT+1))) < 0.0D0) THEN
                  ! calculate the beast Gspr perp - (Gspr perp, Gt perp)* Gt perp/(Gt perp, Gt perp)
                DUMMY=DOT_PRODUCT(GGG(NOPT*JT+1:NOPT*(JT+1)),GGG(NOPT*JT+1:NOPT*(JT+1)))
                IF (DUMMY.EQ.0.0D0) DUMMY=1.0D-100
                  SSS(NOPT*JS+1:NOPT*(JS+1)) = SSS(NOPT*JS+1:NOPT*(JS+1)) - &
                  & DOT_PRODUCT(SSS(NOPT*JS+1:NOPT*(JS+1)),GGG(NOPT*JT+1:NOPT*(JT+1)))*GGG(NOPT*JT+1:NOPT*(JT+1))/ &
                  & DUMMY
             ENDIF
          ENDIF
     END SUBROUTINE NUDGE

     SUBROUTINE TRUEPOTEG
       USE PORFUNCS
       USE NEBDATA
       USE KEYNEB,ONLY: NIMAGE
       USE MODUNRES, ONLY: C,NRES
       USE KEY, ONLY: UNRST, FREEZENODEST, INTEPSILON
       USE INTCOMMONS, ONLY : DESMINT, NNZ, KD, NINTC, PREVDIH, DIHINFO
       USE MODCHARMM, ONLY: CHRMMT, ACESOLV, NCHENCALLS, ACEUPSTEP
       IMPLICIT NONE

       INTEGER :: I,J
       DOUBLE PRECISION :: RMSTMP,HARVEST,DPRAND

       DOUBLE PRECISION :: DUMINT(NINTC)

       ! energy and gradient for images
!      PRINT '(A)',' '
!      PRINT '(A,I8,A,G20.10,A)','image ',1,' energy ',EEE(1),' points:'
!      PRINT '(3G20.10)',XYZ(NOPT*(1-1)+1:NOPT*(1-1)+3)
       DO I=2,NIMAGE+1

          IF (FREEZENODEST) THEN
             IF (DESMINT) THEN
                print*, 'nebgrad>> DESMINT not implemented with freezenodes.'
                STOP
             ENDIF
             ! for nodes that didn't move, set true gradient to previous value
             IF (IMGFREEZE(I-1)) THEN
                GGG(NOPT*(I-1)+1:NOPT*I) = TRUEGRAD(NOPT*(I-1)+1:NOPT*I)
                CYCLE
             ENDIF
          ENDIF

          IF (UNRST) THEN
             IF (DESMINT) THEN
                print*, 'nebgrad>> UNRST not implemented with DESMINT!'
                STOP
             ENDIF
             DO J=1,NRES
                C(1,J)=XYZ(6*(J-1)+1+NOPT*(I-1))
                C(2,J)=XYZ(6*(J-1)+2+NOPT*(I-1))
                C(3,J)=XYZ(6*(J-1)+3+NOPT*(I-1))
                C(1,J+NRES)=XYZ(6*(J-1)+4+NOPT*(I-1))
                C(2,J+NRES)=XYZ(6*(J-1)+5+NOPT*(I-1))
                C(3,J+NRES)=XYZ(6*(J-1)+6+NOPT*(I-1))
             ENDDO
             CALL UPDATEDC
             CALL INT_FROM_CART(.TRUE.,.FALSE.)
             CALL CHAINBUILD
          ENDIF

!bs360: update ACE Born radii for each image
          IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
          IF (DESMINT) THEN

             CALL POTENTIAL(XYZCART(NOPT*(I-1)+1:NOPT*I),EEE(I), &
  &                   GGGCART(NOPT*(I-1)+1:NOPT*I),.TRUE.,.FALSE.,RMSTMP,.FALSE.,.FALSE.)
             PREVDIH => DIHINFO(I,:)
             CALL TRANSFORM(XYZCART(NOPT*(I-1)+1:NOPT*I),GGGCART(NOPT*(I-1)+1:NOPT*I),&
                  & DUMINT,GGG(NOPT*(I-1)+1:NOPT*I), NINTC,NOPT,NNZ,.TRUE.,.FALSE.,KD,INTEPSILON)
             TRUEGRAD(NOPT*(I-1)+1:NOPT*I)=GGGCART(NOPT*(I-1)+1:NOPT*I) ! SAVE FOR PASSING TO BFGSTS LATER
          ELSE
             CALL POTENTIAL(XYZ(NOPT*(I-1)+1:NOPT*I),EEE(I),GGG(NOPT*(I-1)+1:NOPT*I),.TRUE.,.FALSE.,RMSTMP,.FALSE.,.FALSE.)
             TRUEGRAD(NOPT*(I-1)+1:NOPT*I)=GGG(NOPT*(I-1)+1:NOPT*I) ! SAVE FOR PASSING TO BFGSTS LATER
          ENDIF
!         PRINT '(A,I8,A,G20.10,A)','image ',I,' energy ',EEE(I),' points:'
!         PRINT '(3G20.10)',XYZ(NOPT*(I-1)+1:NOPT*(I-1)+3)

          IF ( EEE(I) > HUGE(EEE(I)) ) THEN ! BAD GUESS - HIGH-ENERGY IMAGE
             PRINT *, "IMAGE",I," IS BAD! - TRYING TO LOWER IT's energy..."
             IF (DESMINT) THEN
                print*, 'nebgrad>> bad energy perturbation not implemented with DESMINT yet.'
                STOP
             ENDIF
             DO J=1,NOPT ! CHANGING GEOMETRY RANDOMLY
                !                        call random_number(harvest)
                HARVEST=DPRAND()
                XYZ(NOPT*(I-1)+J) = XYZ(NOPT*(I-1)+J) + HARVEST*0.01
             ENDDO


             IF (UNRST) THEN
                DO J=1,NRES
                   C(1,J)=XYZ(6*(J-1)+1+NOPT*(I-1))
                   C(2,J)=XYZ(6*(J-1)+2+NOPT*(I-1))
                   C(3,J)=XYZ(6*(J-1)+3+NOPT*(I-1))
                   C(1,J+NRES)=XYZ(6*(J-1)+4+NOPT*(I-1))
                   C(2,J+NRES)=XYZ(6*(J-1)+5+NOPT*(I-1))
                   C(3,J+NRES)=XYZ(6*(J-1)+6+NOPT*(I-1))
                ENDDO
                CALL UPDATEDC
                CALL INT_FROM_CART(.TRUE.,.FALSE.)
                CALL CHAINBUILD
             ENDIF

             CALL POTENTIAL(XYZ(NOPT*(I-1)+1:NOPT*I),EEE(I),GGG(NOPT*(I-1)+1:NOPT*I),.TRUE.,.FALSE. &
                  & ,RMSTMP,.FALSE.,.FALSE.)
             TRUEGRAD(NOPT*(I-1)+1:NOPT*I)=GGG(NOPT*(I-1)+1:NOPT*I) ! SAVE FOR PASSING TO BFGSTS LATER

             IF ( EEE(I) > HUGE(EEE(I)) ) THEN
                PRINT *, "Failed."
                CALL TSUMMARY
                STOP
             ENDIF
          ENDIF
       ENDDO
!      PRINT '(A,I8,A,G20.10,A)','image ',NIMAGE+2,' energy ',EEE(NIMAGE+2),' points:'
!      PRINT '(3G20.10)',XYZ(NOPT*(NIMAGE+2-1)+1:NOPT*(NIMAGE+2-1)+3)
     END SUBROUTINE TRUEPOTEG

END MODULE GRADIENTS
