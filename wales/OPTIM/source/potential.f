! 
! OPTIM: A program for optimizing geometries and calculating reaction pathways
! Copyright (C) 1999-2006 David J. Wales
! This file is part of OPTIM.
! 
! OPTIM is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! at your option) any later version.
! 
! OPTIM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
      SUBROUTINE POTENTIAL(COORDS,ENERGY,VNEW,GTEST,STEST,RMS,PTEST,BOXTEST)
         USE COMMONS
         USE KEY
         USE MODNEB
         USE MODHESS
         USE MODCHARMM
         use PORFUNCS
         use SDWATER, ONLY : SDPOTENTIAL, SDGRAD, SDHESS
         USE MCY, ONLY : MCYPOT=>POTENTIAL
         USE BOWMANWATER, ONLY : BOWMANPOT
         USE MBPOLMOD, ONLY: MBPOL
         USE FINITE_DIFFERENCES
         USE MODAMBER9,only : ifswitch,goodstructure1,irespa,cisarray1,checkcistransalways,checkcistransalwaysdna,
     1   checkcistransalwaysrna
         ! hk286
         USE GENRIGID
         USE MULTIPOT, ONLY: MULTIPOT_CALL
         USE CHIRALITY, ONLY: CIS_TRANS_CHECK, CHIRALITY_CHECK
         USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_ENERGY_AND_GRADIENT, POT_ENE_REC_C, AMBER12_NUM_HESS
         USE MODCUDALBFGS, ONLY: CUDA_ENEGRAD_WRAPPER
         USE OPEP_INTERFACE_MOD, ONLY: OPEP_ENERGY_AND_GRADIENT,OPEP_NUM_HESS
         USE ORBITALS_MOD, ONLY: GET_ORBITAL_LOCALITY
         ! use AMHGLOBALS

         IMPLICIT NONE

         DOUBLE PRECISION COORDS(NOPT)
         DOUBLE PRECISION ENERGY
         DOUBLE PRECISION VNEW(NOPT)
         LOGICAL GTEST, STEST
         DOUBLE PRECISION RMS
         LOGICAL PTEST, BOXTEST, file_exists

         INTEGER J1, J2, J3, NN, MM, IPOT, NELEMENTS, NTYPE(105), NCOUNT, ISTART, NDUM, NPCALL, ECALL, FCALL, SCALL,
     1   J4, J5, J6, JSTART, ISTAT, NDUMMY,IDUM1,IDUM2, NATOMSSAVE
         DOUBLE PRECISION P2, P3, BA, XLAMBDA, AMAT(3,3), AINV(3,3), TEMPX, BOXLX,
     1   TEMPXX, TEMPH(3,3), TEMPV(3), ENERGY1, ENERGY2,DOUBLEDUM1,
     1   EDOUBLE, XTEMP, YTEMP, ZTEMP, GAMESR(3,3), GAMEST(3),
     2   TEMPA, TEMP(6), COORDSO(3*NATOMS), GRADO(3*NATOMS), GEMAX,
     4   ETIME, FTIME, STIME, TIME, TIME0, DUMMY1, DUMMY2, DUMMY3, EPLUS, EMINUS, DUMMY, DOTOPT, DIST, DIST1, DIST2, RMAT(3,3)
         INTEGER NSTART, NFINISH, NSTEP
         INTEGER LUNIT, LUNIT2, GETUNIT
         DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::  HGAUSS
         CHARACTER(LEN=87) ESTRING
         CHARACTER(LEN=80) GPSTRING, NSTRING, FSTRING, FNAME, FNAME2, GSTRING,FROMFILE
         CHARACTER(LEN=132) STRING
         COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
         DOUBLE PRECISION C1, C2, C3, IZ, ROT, EDISP, EIND
         LOGICAL ETEST, SSTEST, YESNO, AMIDEFAIL,GAUSSIANTEST
         COMMON /CAS/ AMAT, AINV, NELEMENTS, NTYPE
         COMMON /PCALL/ NPCALL, ECALL, FCALL, SCALL, ETIME, FTIME, STIME
         LOGICAL KNOWE, KNOWG, KNOWH
         COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

         DOUBLE PRECISION VPLUS(3*NATOMS), VMINUS(3*NATOMS), DIFF, EPLUS1, EPLUS2, EMINUS1, EMINUS2, EDUMMY

         ! double precision upperE, lowerE, deltaCoord, numericalGrad(3*NATOMS), RMSdiff
         ! double precision dummyGrad(3*NATOMS), upperGrad(3*NATOMS), lowerGrad(3*NATOMS)
         ! double precision numericalSD, tempHess(3*NATOMS,3*NATOMS)
         double precision HESSDUM(3*NATOMS,3*NATOMS)

         ! sf344> NAB & AMBER additions
         DOUBLE PRECISION,dimension(:),allocatable  ::  temphess
         DOUBLE PRECISION  :: GRAD1(3*NATOMS)
         integer i,j,k
         ! hk286 - local rigid body
         DOUBLE PRECISION :: XCOORDS(3*NATOMS), XRIGIDCOORDS(DEGFREEDOMS), XRIGIDGRAD(DEGFREEDOMS)
         DOUBLE PRECISION :: XRIGIDHESS(DEGFREEDOMS, DEGFREEDOMS)
         DOUBLE PRECISION,ALLOCATABLE :: FOCK(:,:), MOCOEFF(:,:)
         DOUBLE PRECISION :: GRADATOMS(3*NATOMS)
         TYPE(POT_ENE_REC_C) :: ENERGY_DECOMP
         SAVE

! #ifndef _OPTIMLIBRARY
        IF (RIGIDMOLECULEST) THEN   ! Should this be changed to IF(RIGIDINIT)?
            ! sn402: addition to allow generalised rigid bodies to be used for arbitrary potential.
            IF (ATOMRIGIDCOORDT .EQV. .FALSE.) THEN
                IF (DEBUG) THEN
                    IF(ANY(ABS(COORDS(DEGFREEDOMS+1:)) .GT. 1.0E-12)) THEN
                        WRITE(*,*) "Warning: called POTENTIAL with ATOMRIGIDCOORDT = FALSE but coords appear to be atomistic"
!                        WRITE(*,*) "Largest absolute value:", MAXVAL(ABS(COORDS(DEGFREEDOMS+1:)))
                        WRITE(*,*) COORDS(DEGFREEDOMS+1:)
!                         WRITE(*,*) COORDS
                    ENDIF
                ENDIF
                XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
                CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS, XRIGIDCOORDS)
            ELSE IF (DEBUG) THEN
                IF(ANY(ABS(COORDS(DEGFREEDOMS+1:)) .LT. 1.0D-12)) THEN
                    WRITE(*,*) "Warning: called POTENTIAL with ATOMRIGIDCOORDT = TRUE but coords appear to be angleaxis"
!                    WRITE(*,*) "Smallest absolute value:", MINVAL(ABS(COORDS(DEGFREEDOMS+1:)))
                    WRITE(*,*) COORDS(DEGFREEDOMS+1:)
                ENDIF
            ENDIF
        ENDIF
! #endif
         ! 
         ! SSTEST needs to be true if we really want an analytic Hessian this step.
         ! 
         CALL MYCPU_TIME(TIME0,.FALSE.)
         NPCALL=NPCALL+1
         CALL CHANGEP
         SSTEST=.FALSE.
         IF (STEST.AND.(.NOT.HUPDATE)) SSTEST=.TRUE.
         IF (STEST.AND.(NHUP+1.EQ.NSTHUP)) SSTEST=.TRUE.
         IF (STEST.AND.((INTHUP.GT.0).AND.(NHUP+1.GT.NSTHUP))) THEN
            IF (MOD(NHUP+1,INTHUP).EQ.0) SSTEST=.TRUE.
         ENDIF
         KNOWE=.TRUE.
         IF (GTEST) KNOWG=.TRUE.
         IF (SSTEST) KNOWH=.TRUE.

         SHIFTED=.FALSE.
         IF (.NOT.ALLOCATED(HESS)) THEN
            IF (STEST.OR.READHESS.OR.HUPDATE) THEN
               IF (VARIABLES) THEN
                  IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(NATOMS,NATOMS))
                  IF (DEBUG) PRINT '(A,I10)', ' potential> allocating hessian with dimension ',NATOMS
                  ! ELSE IF (RINGPOLYMERT) THEN
               ELSE IF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'AECK') THEN ! asymmetric Eckart barrier
                  IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(NATOMS,NATOMS))
                  IF (DEBUG) PRINT '(A,I10)', ' potential> allocating hessian with dimension ',NATOMS
               ELSE
                  IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
                  IF (DEBUG) PRINT '(A,I10)', ' potential> allocating hessian with dimension ',3*NATOMS
               ENDIF
            ENDIF
         ENDIF
         IF (READHESS) THEN
            WRITE(*,'(A)') ' potential> Reading Hessian from file derivs'
            OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
            IF (CADPAC) THEN
               DO J1=1,3+NATOMS
                  READ(15,*,ERR=666)
               ENDDO
               READ(15,*,ERR=666) (VNEW(J1),J1=1,3*NATOMS)
               READ(15,*,ERR=666)
               READ(15,*,ERR=666) (((HESS(3*(J1-1)+J3,J2),J3=1,3),J2=1,3*NATOMS),J1=1,NATOMS)
            ELSE IF (GAMESSUS) THEN
               OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
11             READ(15,'(A80)') GSTRING
               IF (GSTRING(1:6).NE.' $GRAD') GOTO 11
               READ(15,'(A80)') GSTRING
               DO J1=1,NATOMS
                  READ(15,'(15X,3E20.10)',ERR=666) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ENDDO
               READ(15,*,ERR=666) GSTRING
               READ(15,*,ERR=666) GSTRING
               READ(15,*,ERR=666) GSTRING
               DO J1=1,3*NATOMS
                  JSTART=1
124               READ(15,'(5X,5E15.8)',ERR=666) (HESS(J1,J3),J3=1+5*(JSTART-1),MIN(5+5*(JSTART-1),3*NATOMS))
                  IF (5*JSTART.LT.3*NATOMS) THEN
                     JSTART=JSTART+1
                     GOTO 124
                  ENDIF
               ENDDO
            ELSE IF (GAMESSUK) THEN
               ! 
               ! Since the molecule may well have been reoriented by GAMESS-UK as part of the point group symmetry
               ! analysis, these gradients will have to be transformed if they are required to correspond to the orientation of the
               ! molecule as it was input to GAMESS-UK, and account must also be taken of the reordering of atoms. The
               ! TRANSFORM keyword requests a block of type tr_matrix, containing a rotation and translation matrices
               ! denoted here as R and T respectively) which may be used to construct the gradient in the original, input frame.
               ! There are 3 records, record number i containing R(i,1), R(i,2), R(i,3), T(i) in format (2x,4f15.7). Taking the
               ! coordinates c in the GAMESS-UK coordinate system (i.e. after symmetry adaption, as found in the
               ! coordinates block) the transformation R c-T yields the coordinates of the atom as it was input to
               ! GAMESS-UK, and R g, (where g is the gradient from the gradients block) is the gradient in the initial
               ! coordinate system. Note that the record structure written out in the tr_matrix block is suitable for inclusion
               ! without the block header) after the ORIENT directive.
               ! 
               OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
121            READ(15,'(A80)') GSTRING
               IF (GSTRING(1:17).NE.'block = gradients') GOTO 121
               DO J1=1,NATOMS
                  READ(15,*) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ENDDO
               READ(15,*) GSTRING
               READ(15,*) GSTRING
               READ(15,*) GSTRING
               READ(15,*) ((HESS(J1,J3),J3=1,3*NATOMS),J1=1,3*NATOMS)
               REWIND(15)
161            READ(15,'(A80)') GSTRING
               IF (GSTRING(1:17).NE.'block = tr_matrix') GOTO 161
               READ(15,*) GAMESR(1,1), GAMESR(1,2), GAMESR(1,3), GAMEST(1)
               READ(15,*) GAMESR(2,1), GAMESR(2,2), GAMESR(2,3), GAMEST(2)
               READ(15,*) GAMESR(3,1), GAMESR(3,2), GAMESR(3,3), GAMEST(3)
               DO J1=1,NATOMS
                  DO J3=1,3
                     TEMPX=0.0D0
                     DO J5=1,3
                        TEMPX=TEMPX+GAMESR(J3,J5)*VNEW(3*(J1-1)+J5)
                     ENDDO
                     TEMPV(J3)=TEMPX
                  ENDDO
                  DO J3=1,3
                     VNEW(3*(J1-1)+J3)=TEMPV(J3)
                  ENDDO

                  DO J2=1,NATOMS
                     DO J3=1,3
                        DO J4=1,3
                           TEMPXX=0.0D0
                           DO J5=1,3
                              DO J6=1,3
                                 TEMPXX=TEMPXX+GAMESR(J3,J5)*HESS(3*(J1-1)+J5,3*(J2-1)+J6)*GAMESR(J4,J6)
                              ENDDO
                           ENDDO
                           TEMPH(J3,J4)=TEMPXX
                        ENDDO
                     ENDDO
                     DO J3=1,3
                        DO J4=1,3
                           HESS(3*(J1-1)+J3,3*(J2-1)+J4)=TEMPH(J3,J4)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE
               WRITE(*,'(A)') ' potential> READHESS currently only works for CADPAC, GAMESS-US and GAMESS-UK run types'
            ENDIF
            CLOSE(15)
         ENDIF
         ! 
         ! Spherical container
         ! 
         IF (CONTAINER) CALL RAD(COORDS)
         ! IF (CONTAINER) CALL RAD(COORDS,ENERGY,VNEW,GTEST)

! #ifndef _OPTIMLIBRARY
         IF (RINGPOLYMERT) THEN
            ! 
            ! Get the energy and derivatives corresponding to the true potential for each bead and
            ! construct the total energy, gradient and (if necessary) Hessian.
            ! 
            IF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'AECK') THEN ! asymmetric Eckart barrier
               NZERO=0
               ENERGY=0.0D0
               IF (GTEST) VNEW(1:NOPT)=0.0D0
               IF (STEST) HESS(1:NOPT,1:NOPT)=0.0D0
               DUMMY2=RPMASSES(1)/(RPBETA/RPIMAGES)**2
               DO J1=1,RPIMAGES
                  DUMMY1=COORDS(J1)
                  ENERGY=ENERGY-18.0D0/(3.141592653589793D0 + 3.141592653589793D0*EXP(-0.7674950309598664D0*DUMMY1)) +
     &            4.297183463481175D0/COSH(0.3837475154799332*DUMMY1)**2
                  IF (GTEST) VNEW(J1)=(-18.0D0- 54.0D0*TANH(0.3837475154799332D0*DUMMY1))/
     &            (8.186613663571908D0 + 8.186613663571908D0*COSH(0.7674950309598664D0*DUMMY1))
                  IF (STEST) HESS(J1,J1)=(-2.53125D0 + 1.265625D0*COSH(0.7674950309598664D0*DUMMY1) +
     &            0.421875D0*SINH(0.7674950309598664D0*DUMMY1))/COSH(0.3837475154799332*DUMMY1)**4
               ENDDO
               ! 
               ! Extra terms due to springs between images.
               ! 
               ENERGY=ENERGY+DUMMY2*(COORDS(1)-COORDS(RPIMAGES))**2/2.0D0
               IF (RPIMAGES.GT.2) THEN
                  DO J1=2,RPIMAGES
                     ENERGY=ENERGY+DUMMY2*(COORDS(J1)-COORDS(J1-1))**2/2.0D0
                  ENDDO
               ENDIF
               ! PRINT '(A)','coordinates:'
               ! PRINT '(6G20.10)',COORDS(1:NOPT)
               ! PRINT '(A,G20.10)','energy=',ENERGY
               IF (GTEST) THEN
                  IF (RPIMAGES.GT.2) THEN
                     VNEW(1)=VNEW(1)+(2.0D0*COORDS(1)-COORDS(RPIMAGES)-COORDS(2))*DUMMY2
                     DO J1=2,RPIMAGES-1
                        VNEW(J1)=VNEW(J1)+(2.0D0*COORDS(J1)-COORDS(J1-1)-COORDS(J1+1))*DUMMY2
                     ENDDO
                     VNEW(RPIMAGES)=VNEW(RPIMAGES)+(2*COORDS(RPIMAGES)-COORDS(RPIMAGES-1)-COORDS(1))*DUMMY2
                  ELSEIF (RPIMAGES.EQ.2) THEN
                     VNEW(1)=VNEW(1)+(COORDS(1)-COORDS(2))*DUMMY2
                     VNEW(2)=VNEW(2)+(COORDS(2)-COORDS(1))*DUMMY2
                  ENDIF
                  ! PRINT '(A)','gradient:'
                  ! PRINT '(6G20.10)',VNEW(1:NOPT)
               ENDIF
               IF (STEST) THEN
                  IF (RPIMAGES.GT.2) THEN
                     HESS(1,1)=HESS(1,1)+2.0D0*DUMMY2
                     HESS(1,RPIMAGES)=HESS(1,RPIMAGES)-DUMMY2
                     HESS(1,2)=HESS(1,2)-DUMMY2
                     DO J1=2,RPIMAGES-1
                        HESS(J1,J1)=HESS(J1,J1)+2.0D0*DUMMY2
                        HESS(J1,J1-1)=HESS(J1,J1-1)-DUMMY2
                        HESS(J1,J1+1)=HESS(J1,J1+1)-DUMMY2
                     ENDDO
                     HESS(RPIMAGES,RPIMAGES)=HESS(RPIMAGES,RPIMAGES)+2*DUMMY2
                     HESS(RPIMAGES,RPIMAGES-1)=HESS(RPIMAGES,RPIMAGES-1)-DUMMY2
                     HESS(RPIMAGES,1)=HESS(RPIMAGES,1)-DUMMY2
                  ELSEIF (RPIMAGES.EQ.2) THEN
                     HESS(1,1)=HESS(1,1)+2.0D0*DUMMY2
                     HESS(1,2)=HESS(1,2)-DUMMY2
                     HESS(2,1)=HESS(2,1)-DUMMY2
                     HESS(2,2)=HESS(2,2)+2*DUMMY2
                  ENDIF
                  ! PRINT '(A)','Hessian:'
                  ! PRINT '(6G20.10)',HESS(1:NOPT,1:NOPT)
               ENDIF
            ELSE
               NZERO=0
               ENERGY=0.0D0
               IF (GTEST) VNEW(1:NOPT)=0.0D0
               IF (STEST) HESS(1:NOPT,1:NOPT)=0.0D0
               IF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'NIMET') THEN ! nickel(100) metal + hydrogen
                  DUMMY2=1.0D0/(RPBETA/RPIMAGES)**2 ! atomic units habr =1
                  DO J1=1,RPIMAGES
                     CALL NIMETP(COORDS(RPDOF*(J1-1)+1:RPDOF*J1), VNEW(RPDOF*(J1-1)+1:RPDOF*J1), DUMMY1, GTEST, SSTEST)
                     print*, ENERGY, 'EneTest'
                     ! ENERGY
                     ! IF (GTEST.OR.STEST)
                  ENDDO

               ELSEIF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'SD') THEN ! Stillinger-David flexible water potential
                  DUMMY2=1.0D0/(RPBETA*15.1787D0/RPIMAGES)**2 ! hbar = 15.1787 fs kcal / mol
                  DO J1=1,RPIMAGES
                     ENERGY=ENERGY+SDPOTENTIAL(COORDS(RPDOF*(J1-1)+1:RPDOF*J1))
                     IF (GTEST.OR.STEST) VNEW(RPDOF*(J1-1)+1:RPDOF*J1)=SDGRAD(COORDS(RPDOF*(J1-1)+1:RPDOF*J1))
                     IF (STEST) HESS(RPDOF*(J1-1)+1:RPDOF*J1,RPDOF*(J1-1)+1:RPDOF*J1) =
     &               SDHESS(COORDS(RPDOF*(J1-1)+1:RPDOF*J1),VNEW(RPDOF*(J1-1)+1:RPDOF*J1))
                  ENDDO
               ELSEIF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'TT') THEN ! Xantheas' TTM3-F water potential
                  DUMMY2=1.0D0/(RPBETA*15.1787D0/RPIMAGES)**2 ! hbar = 15.1787 fs kcal / mol
                  DO J1=1,RPIMAGES
                     CALL TTM3FCALL(RPDOF/9,COORDS(RPDOF*(J1-1)+1:RPDOF*J1),ENERGY,VNEW(RPDOF*(J1-1)+1:RPDOF*J1))
                     IF (STEST) THEN
                        PRINT *, 'no hessian for TTM3-F'
                        STOP
                     ENDIF
                  ENDDO
               ELSEIF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'MCY') THEN !VRT(MCY-5f) water potential
                  DUMMY2=1.0D0/(RPBETA/RPIMAGES)**2 ! atomic units
                  DO J1=1,RPIMAGES
                     ENERGY=ENERGY+MCYPOT(COORDS(RPDOF*(J1-1)+1:RPDOF*J1))
                     IF (GTEST) VNEW(RPDOF*(J1-1)+1:RPDOF*J1)=FINDIFGRAD(COORDS(RPDOF*(J1-1)+1:RPDOF*J1), MCYPOT, 1.0D-3, GRAD4T)
                     IF (STEST) HESS(RPDOF*(J1-1)+1:RPDOF*J1,RPDOF*(J1-1)+1:RPDOF*J1) =
     &               FINDIFHESS_POT(COORDS(RPDOF*(J1-1)+1:RPDOF*J1),MCYPOT,1.0D-3)
                  ENDDO
               ELSEIF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'JB') THEN ! Joel Bowman's water potential
                  DUMMY2=1.0D0/(RPBETA/RPIMAGES)**2 ! atomic units (hbar = 1)
                  DO J1=1,RPIMAGES
                     ENERGY=ENERGY+BOWMANPOT(COORDS(RPDOF*(J1-1)+1:RPDOF*J1)) ! Hartrees
                     IF (GTEST) VNEW(RPDOF*(J1-1)+1:RPDOF*J1)=FINDIFGRAD(COORDS(RPDOF*(J1-1)+1:RPDOF*J1), BOWMANPOT, 1.0D-3, GRAD4T)
                     IF (STEST) HESS(RPDOF*(J1-1)+1:RPDOF*J1,RPDOF*(J1-1)+1:RPDOF*J1) =
     &               FINDIFHESS_POT(COORDS(RPDOF*(J1-1)+1:RPDOF*J1),BOWMANPOT,1.0D-3)
                  ENDDO
               ELSEIF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'MB') THEN ! MBPOL
                  DUMMY2=1.0D0/(RPBETA/RPIMAGES)**2 ! atomic units (hbar = 1)
                  DO J1=1,RPIMAGES
                     CALL MBPOL(COORDS(RPDOF*(J1-1)+1:RPDOF*J1),EDUMMY,VNEW(RPDOF*(J1-1)+1:RPDOF*J1),GTEST)
                     ENERGY=ENERGY+(EDUMMY*0.00159360144367D0)
                     IF (STEST) THEN
                        PRINT *, 'no hessian for MBPOL'
                        STOP
                     ENDIF
                  ENDDO
                  VNEW(:)= VNEW*0.00159360144367D0*0.52917721067D0
               ELSEIF (TRIM(ADJUSTL(RPSYSTEM)).EQ.'MAL') THEN ! malonaldehyde
                  DUMMY2=1.0D0/(RPBETA/RPIMAGES)**2 ! atomic units (hbar = 1)
                  DO J1=1,RPIMAGES
                     CALL MALPES(9, COORDS(RPDOF*(J1-1)+1:RPDOF*J1),  VNEW(RPDOF*(J1-1)+1:RPDOF*J1), DUMMY, GTEST, SSTEST)
                     ENERGY=ENERGY+DUMMY ! Hartrees
                     IF (SSTEST) THEN
                        PRINT *, 'no hessian for TTM3-F'
                        STOP
                     ENDIF
                  ENDDO
               ELSE
                  PRINT '(A)',' potential> ERROR *** unrecognised RP system type ',TRIM(ADJUSTL(RPSYSTEM))
                  STOP
               ENDIF
               ! 
               ! Extra terms due to springs between images.
               ! 
               IF (RPCYCLICT) THEN
                  DO J2=1,RPDOF
                     ENERGY=ENERGY+RPMASSES(J2)*DUMMY2*(COORDS(J2)-COORDS(RPDOF*(RPIMAGES-1)+J2))**2/2.0D0
                  ENDDO
               ELSEIF (RPFIXT) THEN
                  DO J2=1,RPDOF
                     ENERGY=ENERGY+RPMASSES(J2)*DUMMY2*(COORDS(J2)-XMINA(J2))**2/2.0D0
                     ENERGY=ENERGY+RPMASSES(J2)*DUMMY2*(COORDS(RPDOF*(RPIMAGES-1)+J2)-XMINB(J2))**2/2.0D0
                  ENDDO
               ENDIF
               IF (RPIMAGES.GT.1) THEN
                  DO J1=1,RPIMAGES-1
                     DO J2=1,RPDOF
                        ENERGY=ENERGY+RPMASSES(J2)*DUMMY2*(COORDS(RPDOF*(J1-1)+J2)-COORDS(RPDOF*J1+J2))**2/2.0D0
                     ENDDO
                  ENDDO
               ENDIF
               ! PRINT '(A)','coordinates:'
               ! PRINT '(6G20.10)',COORDS(1:NOPT)
               ! PRINT '(A,G20.10)','energy=',ENERGY
               IF (GTEST) THEN
                  IF (RPIMAGES.GT.1) THEN
                     IF (RPCYCLICT) THEN
                        DO J2=1,RPDOF
                           VNEW(J2)=VNEW(J2)+RPMASSES(J2)*DUMMY2*(2*COORDS(J2)-COORDS(RPDOF*(RPIMAGES-1)+J2)-COORDS(RPDOF+J2))
                        ENDDO
                        DO J2=1,RPDOF
                           VNEW(RPDOF*(RPIMAGES-1)+J2)=VNEW(RPDOF*(RPIMAGES-1)+J2)+
     &                     RPMASSES(J2)*DUMMY2*(2*COORDS(RPDOF*(RPIMAGES-1)+J2)-COORDS(RPDOF*(RPIMAGES-2)+J2)-COORDS(J2))
                        ENDDO
                     ELSEIF (RPFIXT) THEN
                        DO J2=1,RPDOF
                           VNEW(J2)=VNEW(J2)+RPMASSES(J2)*DUMMY2*(2*COORDS(J2)-XMINA(J2)-COORDS(RPDOF+J2))
                        ENDDO
                        DO J2=1,RPDOF
                           VNEW(RPDOF*(RPIMAGES-1)+J2)=VNEW(RPDOF*(RPIMAGES-1)+J2)+
     &                     RPMASSES(J2)*DUMMY2*(2*COORDS(RPDOF*(RPIMAGES-1)+J2)-COORDS(RPDOF*(RPIMAGES-2)+J2)-XMINB(J2))
                        ENDDO
                     ELSE
                        DO J2=1,RPDOF
                           VNEW(J2)=VNEW(J2)+RPMASSES(J2)*DUMMY2*(COORDS(J2)-COORDS(RPDOF+J2))
                        ENDDO
                        DO J2=1,RPDOF
                           VNEW(RPDOF*(RPIMAGES-1)+J2)=VNEW(RPDOF*(RPIMAGES-1)+J2)+
     &                     RPMASSES(J2)*DUMMY2*(COORDS(RPDOF*(RPIMAGES-1)+J2)-COORDS(RPDOF*(RPIMAGES-2)+J2))
                        ENDDO
                     ENDIF
                     DO J1=2,RPIMAGES-1
                        DO J2=1,RPDOF
                           VNEW(RPDOF*(J1-1)+J2)=VNEW(RPDOF*(J1-1)+J2)+
     &                     RPMASSES(J2)*DUMMY2*(2*COORDS(RPDOF*(J1-1)+J2)-COORDS(RPDOF*(J1-2)+J2)-COORDS(RPDOF*J1+J2))
                        ENDDO
                     ENDDO
                  ENDIF
                  ! PRINT '(A)','gradient:'
                  ! PRINT '(6G20.10)',VNEW(1:NOPT)
               ENDIF
               IF (STEST) THEN
                  IF (RPIMAGES.GT.1) THEN
                     ! J1 = 1
                     DO J2=1,RPDOF
                        HESS(J2,J2) = HESS(J2,J2) + 2.0D0*RPMASSES(J2)*DUMMY2
                     END DO
                     DO J1=2,RPIMAGES-1
                        DO J2=1,RPDOF
                           HESS(RPDOF*(J1-1)+J2,RPDOF*(J1-1)+J2) = HESS(RPDOF*(J1-1)+J2,RPDOF*(J1-1)+J2) + 2.0D0*RPMASSES(J2)*DUMMY2
                           HESS(RPDOF*(J1-1)+J2,RPDOF*(J1-2)+J2) = HESS(RPDOF*(J1-1)+J2,RPDOF*(J1-2)+J2) - RPMASSES(J2)*DUMMY2
                           HESS(RPDOF*(J1-1)+J2,RPDOF*(J1)+J2) = HESS(RPDOF*(J1-1)+J2,RPDOF*(J1)+J2) - RPMASSES(J2)*DUMMY2
                        END DO
                     END DO
                     ! J1 = RPIMAGES
                     DO J2=1,RPDOF
                        HESS(RPDOF*(RPIMAGES-1)+J2,RPDOF*(RPIMAGES-1)+J2) = HESS(RPDOF*(RPIMAGES-1)+J2,RPDOF*(RPIMAGES-1)+J2)
     &                  + 2.0D0*RPMASSES(J2)*DUMMY2
                        HESS(RPDOF*(RPIMAGES-1)+J2,RPDOF*(RPIMAGES-2)+J2) = HESS(RPDOF*(RPIMAGES-1)+J2,RPDOF*(RPIMAGES-2)+J2)
     &                  - RPMASSES(J2)*DUMMY2
                        HESS(RPDOF*(RPIMAGES-1)+J2,J2) = HESS(RPDOF*(RPIMAGES-1)+J2,J2) - RPMASSES(J2)*DUMMY2
                     END DO
                     IF (RPCYCLICT) THEN
                        ! J1 = 1
                        DO J2=1,RPDOF
                           HESS(J2,RPDOF*(RPIMAGES-1)+J2) = HESS(J2,RPDOF*(RPIMAGES-1)+J2) - RPMASSES(J2)*DUMMY2
                           HESS(J2,RPDOF+J2) = HESS(J2,RPDOF+J2) - RPMASSES(J2)*DUMMY2
                        END DO
                        ! J1 = RPIMAGES
                        DO J2=1,RPDOF
                           HESS(RPDOF*(RPIMAGES-1)+J2,RPDOF*(RPIMAGES-2)+J2) = HESS(RPDOF*(RPIMAGES-1)+J2,RPDOF*(RPIMAGES-2)+J2)
     &                     - RPMASSES(J2)*DUMMY2
                           HESS(RPDOF*(RPIMAGES-1)+J2,J2) = HESS(RPDOF*(RPIMAGES-1)+J2,J2) - RPMASSES(J2)*DUMMY2
                        END DO
                     ENDIF
                  ENDIF
                  IF (DEBUG) PRINT *, ' potential> hessian created'
                  ! PRINT '(A)','Hessian:'
                  ! PRINT '(6G20.10)',HESS(1:NOPT,1:NOPT)
               ENDIF
            ENDIF
            IF (PTEST) THEN
               WRITE(*,10) ' Energy for last cycle=',ENERGY
               WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY
            ENDIF
            
         ELSE IF(MULTIPOTT) THEN
            CALL MULTIPOT_CALL(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
            ENDIF

         ELSE IF(EX1DT) THEN
            
            CALL example1Dpotential(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
            ENDIF


         ELSE IF (VARIABLES) THEN

            IF (MLPVB3T) THEN
               CALL MLPVB3(COORDS, VNEW, ENERGY, GTEST, SSTEST)

!               DIFF=1.0D-4
!               PRINT*,'analytic and numerical gradients:'
!               IF (.NOT.(ALLOCATED(HESS))) ALLOCATE(HESS(NMLP,NMLP))
!               CALL MLPVB3(COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
!               PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!               HESSDUM(1:NMLP,1:NMLP)=HESS(1:NMLP,1:NMLP)
!               DO J1=1,NATOMS
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  CALL MLPVB3(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                  CALL MLPVB3(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  IF ((ABS(VNEW(J1)).NE.0.0D0).AND. 
!    &               (ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D0)) THEN
!                     WRITE(*,'(A,I5,2F20.10,A)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),'   X'
!                  ELSE
!                     WRITE(*,'(A,I5,2F20.10,A)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!                  ENDIF
!               ENDDO
!               PRINT*,'analytic and numerical second derivatives:'
!               DO J1=1,NATOMS
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  CALL MLPVB3(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                  CALL MLPVB3(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  DO J2=1,NATOMS
!                     IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND. 
!    &                  (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
!                     WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
!                     ELSE
!                        WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!                     ENDIF
!                  ENDDO
!                ENDDO
!                STOP
            ELSEIF (MLPB3T) THEN
               CALL MLPB3(COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ELSEIF (MLP3T) THEN
               CALL MLP3(COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ELSEIF (MLQT) THEN
               CALL MLQ(COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ELSEIF (ORBITALS) THEN
!              IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
!              IF (.NOT.(ALLOCATED(HESS))) THEN
!                   PRINT*,'allocating hessian'
!                   ALLOCATE(HESS(1:NROTS,1:NROTS))
!              END IF
               CALL GET_ORBITAL_LOCALITY(COORDS, VNEW, ENERGY, GTEST, SSTEST)
!                DIFF=1.0D-6
!                PRINT*,'analytic and numerical gradients:'
!                IF (.NOT.(ALLOCATED(HESS))) THEN
!                     PRINT*,'allocating hessian'
!                     ALLOCATE(HESS(NROTS,NROTS))
!                END IF
!                !CALL MLPVB3(COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
!                CALL GET_ORBITAL_LOCALITY(COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
!                !PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!                !HESSDUM(1:NROTS,1:NROTS)=HESS(1:NROTS,1:NROTS)
!                !DO J1=1,NATOMS
!                !   COORDS(J1)=COORDS(J1)+DIFF
! !                  CALL MLPVB3(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
!                !   CALL GET_ORBITAL_LOCALITY(COORDS, VPLUS, EPLUS, .FALSE., .FALSE.)
!                !   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
! !                  CALL MLPVB3(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
!                !   CALL GET_ORBITAL_LOCALITY(COORDS, VMINUS, EMINUS, .FALSE., .FALSE.)
!                !   COORDS(J1)=COORDS(J1)+DIFF
!                !   IF ((ABS(VNEW(J1)).NE.0.0D0).AND. 
!      &         !      (ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D0)) THEN
!                !      WRITE(*,'(A,I5,2F20.10,A)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),'   X'
!                !   ELSE
!                !      WRITE(*,'(A,I5,2F20.10)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!                !   ENDIF
!                !ENDDO
!                !PRINT*,'analytic and numerical second derivatives:'
!                !DO J1=1,NATOMS
!                   !COORDS(J1)=COORDS(J1)+DIFF
! !                  CALL MLPVB3(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
!                   !CALL GET_ORBITAL_LOCALITY(COORDS, VPLUS, ENERGY, .TRUE., .FALSE.)
!                   !COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                   !CALL MLPVB3(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
!                   !CALL GET_ORBITAL_LOCALITY(COORDS, VMINUS, ENERGY, .TRUE., .FALSE.)
!                   !COORDS(J1)=COORDS(J1)+DIFF
!                   !DO J2=1,NATOMS
!                   !   IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND. 
!      &            !      (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
!                   !   WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
!                   !   ELSE
!                   !      WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!                   !   ENDIF
!                   !ENDDO
!                 !ENDDO
!                !STOP
                IF (PTEST) THEN
                   WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                   WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
                ENDIF
            ELSE
               CALL FUNCTIONAL( COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ENDIF
            IF (PTEST.OR.MLPPROB.OR.MLQPROB) THEN
                WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
                WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            ENDIF
!           IF (MLPPROB) STOP
            ! CALL CTEST(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ! CALL TWODFUNC(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! CALL MB(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! CALL P4DIFF(NATOMS,COORDS,VNEW,ENERGY,PARAM1,GTEST,SSTEST)
            ! CALL P4DIFF(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,GTEST,SSTEST)
!           CALL WATERMETHANE(COORDS,ENERGY)
!           DIFF=1.0D-4
!           DO J1=1,6
!              COORDS(J1)=COORDS(J1)+DIFF
!              CALL WATERMETHANE(COORDS,EPLUS)
!              COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!              CALL WATERMETHANE(COORDS,EMINUS)
!              COORDS(J1)=COORDS(J1)+DIFF
!              VNEW(J1)=(EPLUS-EMINUS)/(2*DIFF)
!           ENDDO
!           CALL WATERMETHANE(COORDS,ENERGY)

            IF (PHI4MODT) THEN
                CALL PHI4MODEL(COORDS,VNEW,ENERGY,GTEST,SSTEST)
                IF (PTEST) THEN
                    WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
                    WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
                ENDIF

         ELSE   IF (ONEDAPBCT) THEN

               CALL ENERGY_1D_APBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)
               ! Debug tools
               ! 
               ! DIFF=1.0D-2
               ! PRINT*,'analytic and numerical gradients: NATOMS=',NATOMS
               ! DO J1=1,NATOMS
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! CALL ENERGY_1D_APBC(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
               ! CALL ENERGY_1D_APBC(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*ABS((VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.0.0D0)) THEN
               ! WRITE(*,'(I5,2G20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
               ! ENDIF
               ! ENDDO
               ! PRINT*,'analytic and numerical second derivatives:'
               ! DO J1=1,NATOMS
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! CALL ENERGY_1D_APBC(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
               ! CALL ENERGY_1D_APBC(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! DO J2=1,NATOMS
               ! IF (ABS(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)).GT.1.0D-1) THEN
               ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
               ! ELSE
               ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
               ! ENDIF
               ! ENDDO
               ! ENDDO

               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF

            ELSE IF(ONEDPBCT) THEN

               CALL ENERGY_1D_PBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)

               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF


            ELSE IF(INVTONEDPBCT) THEN

!              CALL ENERGY_1D_PBC_INVT(COORDS,VNEW,ENERGY,GTEST,SSTEST)

               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF

            ELSE IF(INVTTWODPBCT) THEN

               CALL ENERGY_2D_PBC_INVT(COORDS,VNEW,ENERGY,GTEST,SSTEST)

               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF


            ELSE IF(TWODAPBCT) THEN

               IF(NATOMS.NE.(NONEDAPBC**2)) THEN
                  PRINT *, 'potential.f > Number of lattice sites specified does not equal number uunder variables in odata'
               ENDIF

               CALL ENERGY_2D_APBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)

               ! Debug tools
               ! DIFF=1.0D-2
               ! PRINT*,'analytic and numerical gradients: N_lattice, NATOMS=',(NONEDAPBC), NATOMS
               ! DO J1=1,NATOMS
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! CALL ENERGY_2D_APBC(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
               ! CALL ENERGY_2D_APBC(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*ABS((VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.0.0D0)) THEN
               ! WRITE(*,'(I5,2G20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
               ! ENDIF
               ! ENDDO
               ! PRINT*,'analytic and numerical second derivatives:'
               ! DO J1=1,NATOMS
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! CALL ENERGY_2D_APBC(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
               ! CALL ENERGY_2D_APBC(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
               ! COORDS(J1)=COORDS(J1)+DIFF
               ! DO J2=1,NATOMS
               ! IF (ABS(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)).GT.1.0D-1) THEN
               ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
               ! ELSE
               ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
               ! ENDIF
               ! ENDDO
               ! ENDDO

               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF

            ELSE IF(TWODPBCT) THEN

               IF(NATOMS.NE.(NONEDAPBC**2)) THEN
                  PRINT *, 'potential.f > Number of lattice sites specified does not equal number uunder variables in odata'
               ENDIF

               CALL ENERGY_2D_PBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)

               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF




            ELSE IF(THREEDAPBCT) THEN

               IF(NATOMS.NE.(NONEDAPBC**3)) THEN
                  PRINT *, 'potential.f > Number of lattice sites specified does not equal number uunder variables in odata'
               ENDIF

               CALL ENERGY_3D_APBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)


               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF

            ELSE IF( FOURDAPBCT) THEN

               IF(NATOMS.NE.(NONEDAPBC**4)) THEN
                  PRINT *, 'potential.f > Number of lattice sites specified does not equal number uunder variables in odata'
               ENDIF

               CALL ENERGY_4D_APBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)


               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF

            ELSE IF(FOURDPBCT) THEN

               IF(NATOMS.NE.(NONEDAPBC**4)) THEN
                  PRINT *, 'potential.f > Number of lattice sites specified does not equal number uunder variables in odata'
               ENDIF

               CALL ENERGY_4D_PBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)


               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF


            ELSE IF(THREEDPBCT) THEN

               IF(NATOMS.NE.(NONEDAPBC**3)) THEN
                  PRINT *, 'potential.f > Number of lattice sites specified does not equal number under variables in odata'
               ENDIF

               CALL ENERGY_3D_PBC(COORDS,VNEW,ENERGY,GTEST,SSTEST)


               IF (PTEST) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' '
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' '
               ENDIF

            END IF

            ! IF (RESTART) THEN
            ! C
            ! C  This is an inefficient dirty fix
            ! C
            ! CALL FDIMER(COORDS,ENERGY,VNEW,1)
            ! OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
            ! IF (PTEST) THEN
            ! PRINT*,'Reading external derivative information'
            ! ENDIF
            ! READ(15,*) ENERGY
            ! READ(15,*) (VNEW(J1),J1=1,NOPT)
            ! READ(15,*) ((HESS(J2,J3),J3=1,NOPT),J2=1,NOPT)
            ! CLOSE(15)
            ! ELSE
            ! CALL FDIMER(COORDS,ENERGY,VNEW,ITER)
            ! PRINT*,'Analytic Hessian:'
            ! WRITE(*,'(3F20.10)') ((HESS(J1,J2),J1=1,6),J2=1,6)
            ! CALL DIFF(COORDS, 6, VNEW, HESS)
            ! PRINT*,'Numerical Hessian:'
            ! WRITE(*,'(3F20.10)') ((HESS(J1,J2),J1=1,6),J2=1,6)
            ! STOP
            ! PRINT*,'Reading external derivative information'
            ! OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
            ! READ(15,*) ENERGY
            ! READ(15,*) (VNEW(J1),J1=1,NOPT)
            ! READ(15,*) ((HESS(J2,J3),J3=1,NOPT),J2=1,NOPT)
            ! CLOSE(15)
            ! ENDIF
            ! WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' hartree'
            ! WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ! 
            ! The Noel-Onuchic-Whitford Structure-based model: P. C. Whitford, J.
            ! K. Noel, S. Gosavi, A. Schug, K. Y. Sanbonmatsu, and J. N.
            ! Onuchic. An all-atom structure-based potential for proteins: Bridging
            ! minimal models
            ! with all-atom empirical forcefields. Prot. Struct. Func. Bioinfo. 75,
            ! 430-441, 2009.
            ! 
         ELSE IF (ZSYM(NATOMS).EQ.'SB')then


          IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
               XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
               CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS,XRIGIDCOORDS)
            ENDIF

         IF (SSTEST) PRINT '(A)',' potential> ERROR - calling SBM with SSTEST true'
         CALL SBM(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST)
         IF (PTEST) THEN
            WRITE(*,10) ' Energy for last cycle=',ENERGY
            WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY
         ENDIF

            IF ( RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.) ) THEN
               IF (SSTEST) THEN
                  CALL TRANSFORMHESSIAN(HESS, VNEW, XRIGIDCOORDS,XRIGIDHESS, RBAANORMALMODET)
                  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
               ENDIF
               CALL TRANSFORMGRAD(VNEW, XRIGIDCOORDS, XRIGIDGRAD)
               COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
               VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
            ENDIF


         ! 
         ! qSPCFw  flexible water model introduced by Paesani et al. (JCP 125, 184507 (2006))
         ! Coded by Javier.
         ! 
      ELSE IF (QSPCFWT) THEN
         CALL QSPCFW((NATOMS/3),COORDS,VNEW,ENERGY,GTEST)
         IF (PTEST) THEN
            WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
            WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY,' kcal/mol'
         ENDIF

      ELSE IF (MBPOLT) THEN
         CALL MBPOL (COORDS, ENERGY, VNEW, GTEST)
         IF (STEST) THEN
            WRITE(*,*) 'No hessian for MBPol'
            STOP
            END IF

         ! 
         ! qTIP4PF flexible water model introduced by Habershon et al. (JCP 131, 024501 (2009))
         ! Coded by Javier.
         ! 
      ELSE  IF (QTIP4PFT) THEN
         CALL QTIP4PF((NATOMS/3),COORDS,VNEW,ENERGY,GTEST)
         if (STEST) THEN
            print *, 'no hessian for QTIP4PF'
            stop
            end if
            IF (PTEST) THEN
               WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            ! 
            ! Jeremy Richardson's Stillinger-David model
            ! 
         ELSE IF (SDT) THEN
            ENERGY=SDPOTENTIAL(COORDS(1:NOPT))
            IF (GTEST.OR.STEST) VNEW(1:NOPT)=SDGRAD(COORDS(1:NOPT))
            IF (STEST) HESS(1:NOPT,1:NOPT)=SDHESS(COORDS(1:NOPT),VNEW(1:NOPT))
            IF (PTEST) THEN
               WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            ! 
            ! Yimin Wang and Joel Bowman's water potential
            ! 
         ELSE IF (BOWMANT) THEN
            ENERGY=BOWMANPOT(COORDS(1:NOPT)) ! Hartrees
            IF (GTEST) VNEW(1:NOPT)=FINDIFGRAD(COORDS(1:NOPT), BOWMANPOT, 1.0D-3, GRAD4T)
            IF (STEST) THEN
               HESS(1:NOPT,1:NOPT)=FINDIFHESS_POT(COORDS(1:NOPT),BOWMANPOT,1.0D-3)
            END IF
            IF (PTEST) THEN
               WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' hartree'
               WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY,' hartree'
            ENDIF
         ELSE IF (MALONALDEHYDE) THEN
            CALL MALPES(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ENDIF
         ELSE IF (NATBT) THEN
            CALL NATB(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (GLJT) THEN
            ! generalised LJ with no cutoff
            CALL GLJ(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'SV') THEN
            CALL DRVMSEVB(3*NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' MSEVB Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) ' MSEVB Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (WELCH) THEN
            CALL WEL(NATOMS, COORDS, VNEW, ENERGY, APP, AMM, APM, RHO, XQP, XQM, ALPHAP, ALPHAM, ZSYM, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' hartree'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ENDIF
         ELSE IF (TOSI) THEN
            CALL TOSIFUMI(NATOMS, COORDS, VNEW, ENERGY, PARAM1, PARAM2, PARAM3, PARAM4, ZSYM, GTEST, SSTEST)
            IF (TOSIC6) THEN
               CALL TOSIFUMIC6(NATOMS, COORDS, VNEW, EDISP, C6PP, C6MM, C6PM, ZSYM, GTEST, SSTEST)
               IF (PTEST) WRITE(*,'(A,F20.10,A)') ' Dispersion energy=',EDISP,' hartree'
               ENERGY=ENERGY+EDISP
            ENDIF
            IF (TOSIPOL) THEN
               CALL TOSIFUMIPOL(NATOMS, COORDS, VNEW, EIND, ALPHAP, ALPHAM, ZSYM, DAMP, GTEST, SSTEST)
               IF (PTEST) WRITE(*,'(A,F20.10,A)') ' First order induction energy=',EIND,' hartree'
               ENERGY=ENERGY+EIND
            ENDIF
            IF (PTEST) WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
            IF (PTEST) WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
         ELSE IF (SIO2T) THEN
            CALL SIO2(NATOMS, COORDS, VNEW, ENERGY, PARAM1, PARAM2, PARAM3, PARAM4, ZSYM, GTEST, SSTEST)
            IF (SIO2C6T) THEN
               CALL SIO2C6(NATOMS, COORDS, VNEW, EDISP, C6PP, C6MM, C6PM, ZSYM, GTEST, SSTEST)
               IF (PTEST) WRITE(*,'(A,F20.10,A)') ' Dispersion energy=',EDISP,' hartree'
               ENERGY=ENERGY+EDISP
            ENDIF
            ! 
            ! Check numerical first and second derivatives
            ! 
            ! DIFF=1.0D-4
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL SIO2(NATOMS,COORDS,VPLUS,EPLUS,PARAM1, PARAM2, PARAM3, PARAM4,ZSYM,.FALSE.,.FALSE.)
            ! IF (SIO2C6T) THEN
            ! CALL SIO2C6(NATOMS, COORDS, VPLUS, EDISP, C6PP, C6MM, C6PM, ZSYM, .FALSE.,.FALSE.)
            ! EPLUS=EPLUS+EDISP
            ! ENDIF
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL SIO2(NATOMS,COORDS,VMINUS,EMINUS,PARAM1, PARAM2, PARAM3, PARAM4,ZSYM,.FALSE.,.FALSE.)
            ! IF (SIO2C6T) THEN
            ! CALL SIO2C6(NATOMS, COORDS, VMINUS, EDISP, C6PP, C6MM, C6PM, ZSYM, .FALSE.,.FALSE.)
            ! EMINUS=EMINUS+EDISP
            ! ENDIF
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
            ! WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! PRINT*,'analytic and numerical second derivatives:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL SIO2(NATOMS,COORDS,VPLUS,EPLUS,PARAM1, PARAM2, PARAM3, PARAM4,ZSYM,.TRUE.,.FALSE.)
            ! IF (SIO2C6T) CALL SIO2C6(NATOMS, COORDS, VPLUS, EDISP, C6PP, C6MM, C6PM, ZSYM, .TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL SIO2(NATOMS,COORDS,VMINUS,EMINUS,PARAM1, PARAM2, PARAM3, PARAM4,ZSYM,.TRUE.,.FALSE.)
            ! IF (SIO2C6T) CALL SIO2C6(NATOMS, COORDS, VMINUS, EDISP, C6PP, C6MM, C6PM, ZSYM, .TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! DO J2=1,3*NATOMS
            ! IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
            ! 1             (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
            ! ELSE
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! ENDDO
            ! STOP

            IF (PTEST) WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
            IF (PTEST) WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ! ELSE IF (ZSYM(NATOMS).EQ.'CL') THEN
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            ! ELSE IF (ZSYM(NATOMS).EQ.'CL') THEN
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            ! CALL KDIFF(NATOMS, COORDS, VNEW, ENERGY)
            ! CALL KPAIRS(NATOMS, 3*NATOMS, COORDS, ENERGY)
            ! WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
            ! WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
         ELSE IF (SIO2PT) THEN
             CALL SIO2PSHIFT(NATOMS, COORDS, VNEW, ENERGY, PARAM1,PARAM2, PARAM3, PARAM4, GTEST, SSTEST, PTEST, BOXTEST)
         ELSE IF (ZSYM(NATOMS).EQ.'AZ') THEN
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL AZIZ(NATOMS,COORDS,VNEW,ENERGY,1)
            IF (PARAM1.NE.0.0D0) THEN
               ZSTAR=PARAM1
               CALL AXDIFF(NATOMS, COORDS, VNEW, ZSTAR, GTEST, SSTEST)
               CALL AXPAIRS (NATOMS, COORDS, P2, P3, TEMPA, ZSTAR)
               P2=ENERGY
               ENERGY=ENERGY+P3
            ENDIF
            WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
20          FORMAT(A,27X,F20.10)
            WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            IF (PARAM1.NE.0.0D0) THEN
               WRITE(*,20) ' Two-body contribution=   ',P2
               WRITE(*,20) ' Three-body contribution= ',P3
               WRITE(*,20) ' Z parameter=             ',ZSTAR
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'Z1') THEN
            CALL Z1(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST,PARAM1,PARAM2,PARAM3)
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'ZF') THEN
            CALL Z2FASTER(NATOMS,COORDS,PARAM1,PARAM2,PARAM3,ENERGY,VNEW,SSTEST)
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'Z2') THEN
            CALL Z2(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST,PARAM1,PARAM2,PARAM3)
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'DZ') THEN
            CALL DZUGUTOV(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (LJADD4T) THEN
            CALL LJADD4(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
!                 DIFF=1.0D-4
!                 PRINT*,'analytic and numerical gradients:'
!                 IF (.NOT.(ALLOCATED(HESS))) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
!                 CALL LJADD4(NATOMS, COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
!                 PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!                 HESSDUM(1:3*NATOMS,1:3*NATOMS)=HESS(1:3*NATOMS,1:3*NATOMS)
!                 DO J1=1,3*NATOMS
!                    COORDS(J1)=COORDS(J1)+DIFF
!                    CALL LJADD4(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
!                    COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                    CALL LJADD4(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
!                    COORDS(J1)=COORDS(J1)+DIFF
!                    IF ((ABS(VNEW(J1)).NE.0.0D0).AND.
!     &                (ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D0)) THEN
!                       WRITE(*,'(A,I5,3G20.10,A)') 'anal, num, ratio ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),
!     &                            VNEW(J1)*2.0D0*DIFF/(EPLUS-EMINUS),'   X'
!                    ELSE
!                       WRITE(*,'(A,I5,2F20.10)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!                    ENDIF
!                 ENDDO
!!                PRINT*,'analytic and numerical second derivatives:'
!!                DO J1=1,NATOMS
!!                   COORDS(J1)=COORDS(J1)+DIFF
!!                   CALL LJADD4(NATOMS,COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
!!                   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!!                   CALL LJADD4(NATOMS,COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
!!                   COORDS(J1)=COORDS(J1)+DIFF
!!                   DO J2=1,NATOMS
!!                      IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
!!    &                   (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
!!                      WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
!!                      ELSE
!!                         WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!!                      ENDIF
!!                   ENDDO
!!                ENDDO
!                 STOP

            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (LJADD3T) THEN
            IF (CUDAT) THEN
               ! This call copies CPU coordinates to GPU, calculates energy/gradient and copies energy/gradient back to CPU
               CALL CUDA_ENEGRAD_WRAPPER(NATOMS, COORDS, ENERGY, GRADATOMS)
               VNEW(1:3*NATOMS) = GRADATOMS(:)
            ELSE
               CALL LJADD3(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            END IF
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (LJADD2T) THEN
            CALL LJADD2(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
!                DIFF=1.0D-4
!                PRINT*,'analytic and numerical gradients:'
!                IF (.NOT.(ALLOCATED(HESS))) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
!                CALL LJADD2(NATOMS, COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
!                PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!                HESSDUM(1:3*NATOMS,1:3*NATOMS)=HESS(1:3*NATOMS,1:3*NATOMS)
!                DO J1=1,NATOMS
!                   COORDS(J1)=COORDS(J1)+DIFF
!                   CALL LJADD2(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
!                   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                   CALL LJADD2(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
!                   COORDS(J1)=COORDS(J1)+DIFF
!                   IF ((ABS(VNEW(J1)).NE.0.0D0).AND.
!     &                (ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D0)) THEN
!                      WRITE(*,'(A,I5,2F20.10,A)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),'   X'
!                   ELSE
!                      WRITE(*,'(A,I5,2F20.10)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!                   ENDIF
!                ENDDO
!                PRINT*,'analytic and numerical second derivatives:'
!                DO J1=1,NATOMS
!                   COORDS(J1)=COORDS(J1)+DIFF
!                   CALL LJADD2(NATOMS,COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
!                   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                   CALL LJADD2(NATOMS,COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
!                   COORDS(J1)=COORDS(J1)+DIFF
!                   DO J2=1,NATOMS
!                      IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
!     &                   (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
!                      WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
!                      ELSE
!                         WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!                      ENDIF
!                   ENDDO
!                ENDDO
!                STOP

            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (LJADDT) THEN
            CALL LJADD(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'AX') THEN
            ZSTAR=PARAM1
            IF (ZSTAR.EQ.0.0D0) THEN
               CALL LJDIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
               P2=ENERGY; P3=0.0D0
            ELSE
               CALL LJDIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST) ! 2-body derivatives
               CALL AXDIFF(NATOMS, COORDS, VNEW, ZSTAR, GTEST, SSTEST)  ! 3-body derivatives added
               CALL AXPAIRS (NATOMS, COORDS, P2, P3, ENERGY, ZSTAR)     ! AxTell energy
            ENDIF
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
               IF (PARAM1.NE.0.0D0) THEN
                  WRITE(*,20) ' Two-body contribution=   ',P2
                  WRITE(*,20) ' Three-body contribution= ',P3
                  WRITE(*,20) ' Z parameter=             ',ZSTAR
               ENDIF
            ENDIF
!
!  polymer in a pore
!
         ELSE IF (ZSYM(NATOMS).EQ.'PP') THEN
            CALL LJDIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            CALL POLY2(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            CALL PORE8(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
         ELSE IF ((ZSYM(NATOMS).EQ.'SW').OR.(ZSYM(NATOMS).EQ.'SM')) THEN
            IF (ZSYM(NATOMS).EQ.'SW') XLAMBDA=21.0D0
            IF (ZSYM(NATOMS).EQ.'SM') XLAMBDA=21.0D0*1.5D0
            IF (PRESSURE) THEN
               CALL SWLATMIN(NATOMS,COORDS,PARAM1,PARAM2,PARAM3,VNEW,XLAMBDA)
               PRINT*,'Lattice constant optimised'
               PRINT*,'New Box length=',PARAM1
            ENDIF
            CALL SWTWO(NATOMS, COORDS, VNEW, P2, P3, PARAM1, PARAM2, PARAM3, GTEST, SSTEST, XLAMBDA)
            ENERGY=P2+P3
            ! WRITE(*,'(A,3F20.10)') 'In potential, U,PV,H=',ENERGY,PRESS*PARAM1*PARAM2*PARAM3,ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3

            ! PRINT*,'Analytic derivatives:'
            ! WRITE(*,'(3E20.10)') (VNEW(J1),J1=1,3*NATOMS)
            ! WRITE(*,'(3F20.10)') ((HESS(J1,J2),J1=1,3*NATOMS),J2=1,3*NATOMS)
            ! CALL DIFF(COORDS,NATOMS,VNEW,HDUM)
            ! PRINT*,'Numerical derivatives:'
            ! WRITE(*,'(3E20.10)') (VNEW(J1),J1=1,3*NATOMS)
            ! WRITE(*,'(3F20.10)') ((HESS(J1,J2),J1=1,3*NATOMS),J2=1,3*NATOMS)
            ! DO J1=1,3*NATOMS
            ! DO J2=1,3*NATOMS
            ! IF (HESS(J1,J2).NE.0.0D0) THEN
            ! IF (ABS((HESS(J1,J2)-HDUM(J1,J2))/HESS(J1,J2)).GT.0.01D0) THEN
            ! WRITE(*,'(2I4,3F20.10)') J1,J2,HESS(J1,J2),HDUM(J1,J2)
            ! ENDIF
            ! ENDIF
            ! ENDDO
            ! ENDDO

            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(*,20) ' Two-body contribution=   ',P2
               WRITE(*,20) ' Three-body contribution= ',P3
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'QI') THEN
            CALL OPTIM_QUIP_WRAPPER(NATOMS,COORDS,VNEW,ENERGY)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
!                DIFF=1.0D-3
!                PRINT*,'analytic and numerical gradients:'
!                IF (.NOT.(ALLOCATED(HESS))) ALLOCATE(HESS(NOPT,NOPT))
!            CALL OPTIM_QUIP_WRAPPER(NATOMS,COORDS,VNEW,ENERGY)
!                PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!                HESSDUM(1:NOPT,1:NOPT)=HESS(1:NOPT,1:NOPT)
!                DO J1=1,NATOMS
!                   COORDS(J1)=COORDS(J1)+DIFF
!            CALL OPTIM_QUIP_WRAPPER(NATOMS,COORDS,VPLUS,EPLUS)
!                   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                   CALL MLPVB3(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
!            CALL OPTIM_QUIP_WRAPPER(NATOMS,COORDS,VMINUS,EMINUS)
!                   COORDS(J1)=COORDS(J1)+DIFF
!                   IF ((ABS(VNEW(J1)).NE.0.0D0).AND. 
!     &               (ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D0)) THEN
!                      WRITE(*,'(A,I5,3F20.10,A)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF), 
!     &                   ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)),'   X'
!                   ELSE
!                      WRITE(*,'(A,I5,3F20.10,A)') 'anal and num ',J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF), 
!     &                   ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1))
!                   ENDIF
!                ENDDO
         ELSE IF (ZSYM(NATOMS).EQ.'JC') THEN
            CALL JMEC(NATOMS, COORDS, P2, P3, VNEW,ENERGY, PARAM4,GTEST,SSTEST)
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            WRITE(*,10) 'Two-body contribution=',P2,' eV'
            WRITE(*,10) 'Three-body contribution=',P3,' eV'
         ELSE IF (ZSYM(NATOMS).EQ.'CC') THEN
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL JMECC(NATOMS, COORDS, P2, P3, VNEW,ENERGY)
            ! CALL JM2CC(NATOMS, COORDS, VNEW)
            ! CALL JM3CC(NATOMS, COORDS, VNEW)
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            WRITE(*,10) 'Two-body contribution=',P2,' eV'
            WRITE(*,10) 'Three-body contribution=',P3,' eV'
         ELSE IF (ZSYM(NATOMS).EQ.'JM') THEN
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL JMEP(NATOMS,COORDS,P2,P3,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4)
            ! CALL JM2P(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,PARAM4)
            ! CALL JM3P(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,PARAM4)
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            WRITE(*,10) 'Two-body contribution=',P2,' eV'
            WRITE(*,10) 'Three-body contribution=',P3,' eV'
         ELSE IF (ZSYM(NATOMS).EQ.'M') THEN
            IF (BULKT) THEN
               CALL MORSE_BULK_WRAPPER(COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ELSE
               CALL MORSE(NATOMS,COORDS,ENERGY,VNEW,PARAM1,GTEST,SSTEST)
            ENDIF
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(*,20) ' RHO=',PARAM1
            ENDIF
            ! gd351 ->
         ELSE IF (ZSYM(NATOMS).EQ.'AO') THEN
            CALL ASAOOSPOT(NATOMS,COORDS,ENERGY,VNEW,PARAM1,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(*,20) ' SIGMAP=',PARAM1
            ENDIF
            ! gd351
            ! ELSE IF (ZSYM(NATOMS).EQ.'M2') THEN
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            ! CALL M2(NATOMS,COORDS,VNEW,ENERGY,MALPHA2,PARAM2,1)
            ! WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
            ! WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ! WRITE(*,'(A23,7X,2F20.10)') ' RHO and DELTA=',PARAM1, PARAM2
            ! ELSE IF (ZSYM(NATOMS).EQ.'MV') THEN
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            ! CALL MAV(NATOMS,COORDS,VNEW,ENERGY,MALPHA1,PARAM2,1)
            ! WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
            ! WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ! WRITE(*,'(A23,7X,2F20.10)') ' RHO and DELTA=',PARAM1, PARAM2
            ! ELSE IF (ZSYM(NATOMS).EQ.'GV') THEN
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            ! CALL GAV(NATOMS,COORDS,VNEW,ENERGY,GALPHA,PARAM2,1)
            ! WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
            ! WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ! WRITE(*,'(A23,7X,2F20.10)') ' RHO and DELTA=',PARAM1, PARAM2

         ELSEIF (TTM3T) THEN ! Xantheas' TTM3-F water potential
            ENERGY=0.0D0
            VNEW(1:3*NATOMS)=0.0D0
            CALL TTM3FCALL(NATOMS/3,COORDS,ENERGY,VNEW)
            IF (PTEST) THEN
               WRITE(*,'(A,F20.10,A)') ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY
            ENDIF
            IF (STEST) THEN
               PRINT *, 'potential> No hessian for TTM3-F - making a numerical Hessian'
               DIFF=1.0D-4
               DO J1=1,3*NATOMS
                  COORDS(J1)=COORDS(J1)+DIFF
                  VPLUS(1:3*NATOMS)=0.0D0
                  CALL TTM3FCALL(NATOMS/3,COORDS,EPLUS,VPLUS)
                  COORDS(J1)=COORDS(J1)-2.0D0*DIFF
                  VMINUS(1:3*NATOMS)=0.0D0
                  CALL TTM3FCALL(NATOMS/3,COORDS,EMINUS,VMINUS)
                  COORDS(J1)=COORDS(J1)+DIFF
                  DO J2=1,3*NATOMS
                     HESS(J1,J2)=(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
                  ENDDO
               ENDDO
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'TT') THEN
            ! 
            ! PARAM1, PARAM2 and PARAM3 are the boxlengths and PARAM4 is the cutoff fraction.
            ! 
            CALL TIGHTE(NATOMS,COORDS,VNEW,ENERGY,GTEST,PARAM1,PARAM2,PARAM3)
            IF (SSTEST) THEN
               FIXIMAGE=.TRUE.
               CALL SECSI(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3)
               FIXIMAGE=.FALSE.
            ENDIF
            IF (PTEST) THEN
               WRITE(*,'(A,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'LP') THEN
            ! 
            ! PARAM1, PARAM2 and PARAM3 are the boxlengths, PARAM4 is the cutoff.
            ! 
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            IF (BINARY) THEN
               CALL  LJPBIN(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST,PTEST)
            ELSE
               CALL LJPDIFF(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST,PTEST)
            ENDIF
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            ! FIXIMAGE=.TRUE.
            ! CALL DIFF(NATOMS,COORDS,HDUM,VDUM,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4)
            ! PRINT*,'Analytic and numerical derivatives:'
            ! WRITE(*,'(2I4,3F20.10)') ((J1,J2,HESS(J1,J2),HDUM(J1,J2),HESS(J1,J2)/HDUM(J1,J2),J1=1,NOPT),J2=1,NOPT)
            ! STOP
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'LC') THEN
            ! 
            ! PARAM1, PARAM2 and PARAM3 are the boxlengths, PARAM4 is the cutoff.
            ! 
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            IF (.NOT.BINARY) THEN
               PRINT*,'error atom type LC only works for binary LJ'
               STOP
            ENDIF
            CALL  LJPSHIFTBIN2(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST,PTEST)
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'BC') THEN
            ! 
            ! PARAM1 is the cutoff.
            ! 
            CALL LJPSHIFTBINC(NATOMS,COORDS,VNEW,ENERGY,PARAM1,GTEST,SSTEST,PTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF

         ELSE IF (ZSYM(NATOMS).EQ.'LS') THEN
            ! 
            ! PARAM1, PARAM2 and PARAM3 are the boxlengths, PARAM4 is the cutoff.
            ! 
            IF (BINARY) THEN
               CALL  LJPSHIFTBIN(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST,PTEST,BOXTEST)
            ELSE
               IF (PRESSURE) THEN
                  CALL LJPSLATMIN(NATOMS,COORDS,PARAM1,PARAM2,PARAM3,PARAM4,VNEW)
                  PRINT*,'Lattice constant optimised'
                  PRINT*,'New Box length=',PARAM1
               ENDIF
               CALL  LJPSHIFT(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST,PTEST)
            ENDIF
            ! FIXIMAGE=.TRUE.
            ! CALL DIFF(NATOMS,COORDS,HDUM,VDUM,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4)
            ! PRINT*,'Analytic and numerical derivatives:'
            ! DO J1=1,NOPT
            ! DO J2=1,NOPT
            ! IF (DABS(HDUM(J1,J2)).GT.1.0D-10) THEN
            ! IF (DABS(DABS(HESS(J1,J2)/HDUM(J1,J2))-1.0D0).GT.0.01D0)
            ! 1               WRITE(*,'(2I4,3F20.10)') J1,J2,HESS(J1,J2),HDUM(J1,J2),HESS(J1,J2)/HDUM(J1,J2)
            ! ENDIF
            ! ENDDO
            ! ENDDO
            ! DO J1=1,NOPT
            ! IF (DABS(VDUM(J1)).GT.1.0D-10) THEN
            ! IF (DABS(DABS(VNEW(J1)/VDUM(J1))-1.0D0).GT.0.01D0)
            ! 1            WRITE(*,'(I4,3F20.10)') J1,VNEW(J1),VDUM(J1),VNEW(J1)/VDUM(J1)
            ! ENDIF
            ! ENDDO
            ! STOP
            ! WRITE(*,'(A,3F20.10)') 'ENERGY,PV,H=',ENERGY,PRESS*PARAM1*PARAM2*PARAM3,ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'LK') THEN
            ! 
            ! PARAM1, PARAM2 and PARAM3 are the boxlengths, PARAM4 is the cutoff.
            ! 
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            IF (BINARY) THEN
               CALL  LJPKOB(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST,PTEST)
            ELSE
               PRINT*,'BINARY keyword expected for atom type LK - quit'
               STOP
            ENDIF
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'LM') THEN
            ! 
            ! Parameters are epsilon, rm and gamma.
            ! 
            CALL LJMS(NATOMS, PARAM1, PARAM2, PARAM3, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'MP') THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            IF (PRESSURE) THEN
               CALL MLATMIN(NATOMS,COORDS,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5)
               PRINT*,'Lattice constant optimised'
               PRINT*,'New Box length in x =',PARAM2
               PRINT*,'New Box length in y =',PARAM3
               PRINT*,'New Box length in z =',PARAM4
            ENDIF
            CALL MPDIFF(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,GTEST,SSTEST)
            IF (PTEST) WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
         ELSE IF (ZSYM(NATOMS).EQ.'GP') THEN
            CALL  GUPTA(NATOMS,COORDS,VNEW,ENERGY,GTEST,GUPTATYPE)
            IF (PTEST) WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
         ELSE IF ((ZSYM(NATOMS).EQ.'G1') .OR. (ZSYM(NATOMS).EQ.'G2')) THEN
            CALL  BGUPTA(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'

            ! DIFF=1.0D-5
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL BGUPTA(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL BGUPTA(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
            ! WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! PRINT*,'analytic and numerical second derivatives:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL BGUPTA(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL BGUPTA(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! DO J2=1,3*NATOMS
            ! IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
            ! 1             (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
            ! ELSE
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! ENDDO
            ! STOP
            ! ds656>
         ELSE IF (ZSYM(NATOMS) .EQ. 'GY') THEN
            ! CALL  GLJYPOT(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
            ! <ds656
         ELSE IF (ZSYM(NATOMS).EQ.'DS') THEN
            CALL MPDIFFDS(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,GTEST,SSTEST)

            ! DIFF=1.0D-5
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! DO J1=18,18
            ! IF (FROZEN((J1-1)/3+1)) CYCLE
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL MPDIFFDS(NATOMS,COORDS,VPLUS,EPLUS,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL MPDIFFDS(NATOMS,COORDS,VMINUS,EMINUS,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
            ! WRITE(*,'(I5,2G20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! STOP

            IF (PTEST) WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
         ELSE IF (ZSYM(NATOMS).EQ.'MS') THEN
            IF (PRESSURE) THEN
               CALL MSLATMIN(NATOMS,COORDS,PARAM1,PARAM2,PARAM3,PARAM4)
               PRINT*,'Lattice constant optimised'
               PRINT*,'New Box length in x=',PARAM2
               PRINT*,'New Box length in y=',PARAM3
            ENDIF
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL MSDIFF(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4)
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
         ELSE IF (ZSYM(NATOMS).EQ.'CK') THEN
            CALL ECTRAP(NATOMS,COORDS,ENERGY,C1,C2,C3)
            CALL DCTRAP(NATOMS,COORDS,VNEW,C1,C2,C3)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' units'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' units'
            ENDIF
         ELSE IF (EYTRAPT) THEN
            CALL EYETRAP(NATOMS,COORDS,VNEW,ENERGY,C1,C2,C3)
            IF (SSTEST) CALL EYDTRAP(NATOMS,COORDS,VNEW,C1,C2,C3,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ENDIF
         ELSE IF (MKTRAPT) THEN
            CALL MKTRAP(NATOMS,COORDS,VNEW,ENERGY,SSTEST)

!              DIFF=1.0D-4
!              PRINT*,'analytic and numerical gradients:'
!              IF (.NOT.(ALLOCATED(HESS))) ALLOCATE(HESS(NOPT,NOPT))
!              CALL MKTRAP(NATOMS, COORDS, VNEW, ENERGY,.TRUE.)
!              PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!              HESSDUM(1:NOPT,1:NOPT)=HESS(1:NOPT,1:NOPT)
!              DO J1=1,NATOMS
!                 COORDS(J1)=COORDS(J1)+DIFF
!                 CALL MKTRAP(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.)
!                 COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                 CALL MKTRAP(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.)
!                 COORDS(J1)=COORDS(J1)+DIFF
!                 IF ((ABS(VNEW(J1)).NE.0.0D0).AND.ABS(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.-1.0D0) THEN
!                    WRITE(*,'(I5,2G25.15)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!                 ENDIF
!              ENDDO
!              PRINT*,'analytic and numerical second derivatives:'
!              DO J1=1,NATOMS
!                 COORDS(J1)=COORDS(J1)+DIFF
!                 CALL MKTRAP(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.)
!                 COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                 CALL MKTRAP(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.)
!                 COORDS(J1)=COORDS(J1)+DIFF
!                 DO J2=1,NATOMS
!                    IF ((ABS(HESS(J1,J2)).GT.1.0D-100)) THEN
!                       IF ((ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D-3)) THEN
!                          WRITE(*,'(2I5,2G25.15,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
!                       ELSE
!                          WRITE(*,'(2I5,2G25.15,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!                       ENDIF
!                    ELSE
!                       WRITE(*,'(2I5,2G25.15,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!                    ENDIF
!                 ENDDO
!              ENDDO
!          PRINT*,'coords in potential:'
!          WRITE(*,'(3G25.15)') (COORDS(J1),J1=1,NOPT)
!              STOP

            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY ! ,' hartree'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY ! ,' hartree'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'C1') THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            CALL C10(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' arbs'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' arbs'
            ENDIF
         ELSE IF (DF1T) THEN
            CALL DF1GRAD(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST,PARAM1,PARAM2)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (BLNT) THEN
            CALL BLN(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST.AND.(.NOT.PULLT).AND.(.NOT.TWISTT)) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'PL') THEN
            CALL P46MERDIFF(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' arbs'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' arbs'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'GL') THEN
            CALL G46MERDIFF(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' arbs'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' arbs'
            ENDIF
            ! ELSE IF (ZSYM(NATOMS).EQ.'AK') THEN
         ELSE IF (ACK1.OR.ACK2) THEN
            CALL ACK(NATOMS,COORDS,ENERGY,VNEW,PARAM1,PARAM2,PARAM3,PARAM4,PRESSURE,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'SC') THEN
            CALL OESCP(NATOMS,COORDS,ENERGY,VNEW,PARAM1,PARAM2,PARAM3,PARAM4,PRESSURE,GTEST,SSTEST)
            ! IF (GTEST.OR.SSTEST) CALL DSCP(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,PARAM4,GTEST,SSTEST)
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'PR') THEN
            CALL PRC60(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            ZSTAR=PARAM1
            IF (ZSTAR.NE.0.0D0) THEN
               CALL AXDIFF(NATOMS, COORDS, VNEW, ZSTAR, GTEST, SSTEST)
               CALL AXPAIRS (NATOMS, COORDS, P2, P3, TEMPA, ZSTAR)
               P2=ENERGY
               ENERGY=ENERGY+P3
            ENDIF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'C6') THEN
            CALL C60DIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'P6') THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            ! 
            ! Lattice constant optimisation if required: (zero external pressure)
            ! The coordinates may get changed for lattice optimisation, so use Q.
            ! 
            IF (PRESSURE) THEN
               ! ENERGY=GOLDEN(NATOMS,COORDS,PARAM1+0.001D0,PARAM1,PARAM1-0.001D0,
               ! 1           1.0D-10,XMIN,PARAM4)
               CALL LATMIN(NATOMS,COORDS,PARAM1,PARAM4)
               PRINT*,'Lattice constant optimised'
               PARAM2=PARAM1
               PARAM3=PARAM1
               PRINT*,'New Box length=',PARAM1
               PRINT*,'Cutoff (changed in proportion) is now ',PARAM4
            ENDIF
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL C60P(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4)
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ! ELSE IF (ZSYM(NATOMS).EQ.'TB') THEN
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            ! CALL TBE(NATOMS,COORDS,ENERGY)
            ! WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            ! WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ! CALL TBDIFF(NATOMS,COORDS,VNEW)
            ! CALL DIFF1(NATOMS,COORDS,HESS)
            ! CALL DIFF2(NATOMS,COORDS,HESS)
         ELSE IF (ZSYM(NATOMS).EQ.'FH') THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            BA=0.5291772D0
            PRINT*,'NATOMS,ZSYM(NATOMS)=',NATOMS,ZSYM(NATOMS)
            OPEN(UNIT=15,FILE='fhderivs',STATUS='OLD')
            READ(15,*) (VNEW(J1),J1=1,3*NATOMS)
            READ(15,*) (((HESS(3*(J1-1)+J3,J2),J3=1,3),J2=1,3*NATOMS),J1=1,NATOMS)
            DO J1=1,3*NATOMS
               VNEW(J1)=VNEW(J1)*BA
               PRINT*,'J2=',J2
               DO J2=1,NOPT
                  HESS(J2,J1)=HESS(J2,J1)*BA*BA
               ENDDO
            ENDDO
         ELSE IF (ZSYM(NATOMS).EQ.'ME') THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL EMIE(NATOMS,COORDS,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PRESSURE)
            CALL MIED(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6)
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
         ELSE IF (ZSYM(NATOMS).EQ.'SY') THEN
            CALL STOCK(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! DIFF=1.0D-5
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL STOCK(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL STOCK(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
            ! WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! PRINT*,'analytic and numerical second derivatives:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL STOCK(NATOMS,COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL STOCK(NATOMS,COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! DO J2=1,3*NATOMS
            ! IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
            ! 1             (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
            ! ELSE
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! ENDDO
            ! STOP

            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'T4') THEN
            CALL GENRIGIDTIP(COORDS,VNEW,ENERGY,GTEST,SSTEST)
         ELSE IF (DBPT) THEN
            CALL DUMBBELLP(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (DBPTDT) THEN

            CALL DMBLTD(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (DMBLPYT) THEN
            CALL DUMBBELLPOLARYUKAWA(COORDS,VNEW,ENERGY,GTEST,SSTEST)


         ELSE IF (NIMET) THEN
            CALL NIMETP(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (NIHEAM7T) THEN
            CALL NIHEAM7(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (NIH2LEPST) THEN
            CALL NIH2LEPS(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (LWOTPT) THEN
            IF (BULKT) THEN
               CALL LWOTPBOX(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ELSE
               CALL LWOTPGH(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ENDIF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF
         ELSE IF (MSSTOCKT) THEN
            CALL MSSTOCKGH(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF
         ELSE IF (MULTISITEPYT) THEN
            CALL MULTISITEPY2 (COORDS, VNEW, ENERGY, GTEST)
            IF (SSTEST) THEN
               CALL MULTISITEPYSECDER(COORDS,SSTEST)
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF
         ELSE IF (DJWRBT) THEN
            IF (DJWRBID.EQ.1) THEN
               CALL DJWGR1(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)

!               DIFF=1.0D-4
!               PRINT*,'analytic and numerical gradients:'
!               IF (.NOT.(ALLOCATED(HESS))) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
!               CALL DJWGR1(NATOMS,COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
!               PRINT '(A,I8)','SIZE(HESS)=',SIZE(HESS)
!               HESSDUM(1:3*NATOMS,1:3*NATOMS)=HESS(1:3*NATOMS,1:3*NATOMS)
!               DO J1=1,3*NATOMS
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  CALL DJWGR1(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                  CALL DJWGR1(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!               ENDDO
!               PRINT*,'analytic and numerical second derivatives:'
!               DO J1=1,3*NATOMS
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  CALL DJWGR1(NATOMS,COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)-2.0D0*DIFF
!                  CALL DJWGR1(NATOMS,COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
!                  COORDS(J1)=COORDS(J1)+DIFF
!                  DO J2=1,3*NATOMS
!                     DUMMY1=(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
!                     IF ((ABS(DUMMY1).GT.1.0D-10).AND. 
!     &                   (ABS(100.0D0*(HESS(J1,J2)-DUMMY1)/DUMMY1).GT.1.0D0)) THEN
!                     WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),DUMMY1,'   X'
!                     ELSE
!                        WRITE(*,'(2I5,2G20.10,A)') J1,J2,HESS(J1,J2),DUMMY1
!                     ENDIF
!                  ENDDO
!               ENDDO
!               STOP

             ENDIF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            ENDIF
         ELSE IF (NCAPT) THEN
            CALL NEWCAPSID (COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF
         ELSE IF (NTIPT) THEN
            CALL NEWTIP(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (GBT) THEN
            IF (BULKT) THEN
               CALL GBBOX(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ELSE
               CALL GB(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ENDIF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            ENDIF

         ELSE IF (GBDT) THEN

            CALL GBD(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (PAHAT) THEN
            CALL PAHAGH(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (PAPT) THEN
            CALL PAP(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (PATCHYDT) THEN
            CALL PATCHYD(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (PYGT) THEN
            CALL PYG (COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (SSTEST) THEN
               CALL PYGSECDER(COORDS,SSTEST)
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF
         ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
            IF (PYADDT) THEN
               CALL PYGPERIODICADD(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ELSEIF (PYADD2T) THEN
               CALL PYGPERIODICADD2(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ELSE
               CALL PYGPERIODIC (COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ENDIF
            ! CALL PARAMONOVNUMFIRSTDER (COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! STOP
            ! WRITE(*,*) 'VNEW='
            ! DO J1=1,3*NATOMS
            ! WRITE(*,*) VNEW(J1)
            ! END DO
            ! STOP
            IF (SSTEST) THEN
               IF (PYADDT) THEN
                  CALL PYGPERIODICSECDERADD(COORDS,SSTEST)
               ELSE IF (PYADD2T) THEN 
                  CALL PYGPERIODICSECDERADD2(COORDS,SSTEST)
               ELSE
                  CALL PYGPERIODICSECDER(COORDS,SSTEST)
               ENDIF
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF
         ELSE IF (PYT) THEN
          call py(COORDS,VNEW,ENERGY,.true.)
            IF (SSTEST) THEN
               CALL py_secder(COORDS,SSTEST)
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (PTSTSTT) THEN
            CALL PTSTST(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF


         ELSE IF (STOCKAAT) THEN
            CALL STOCKGHAA(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            ENDIF
         ELSE IF (SANDBOXT) THEN
            CALL SANDBOX (COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (SSTEST) THEN
               CALL SANDBOXSECDER(COORDS,SSTEST)
            END IF

            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            ENDIF
       
         ELSE IF (SILANET) THEN
            CALL SILANE(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'         '
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'         '
            END IF

         ELSE IF (VASP) THEN
            FNAME='POSCAR'
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            OPEN(UNIT=8,FILE='POSCAR.new',STATUS='UNKNOWN')

            DO J1=1, 8                      !Writes out first 8 lines of file
               READ(7,'(A80)') FNAME
               WRITE(8,'(A80)') FNAME
            ENDDO
            DO J1=1, NATOMS                                !Updates coordinates from array COORDS
               READ(7,*) FNAME
               XTEMP=AINV(1,1)*COORDS(3*(J1-1)+1)+AINV(1,2)*COORDS(3*(J1-1)+2)+AINV(1,3)*COORDS(3*(J1-1)+3)
               YTEMP=AINV(2,1)*COORDS(3*(J1-1)+1)+AINV(2,2)*COORDS(3*(J1-1)+2)+AINV(2,3)*COORDS(3*(J1-1)+3)
               ZTEMP=AINV(3,1)*COORDS(3*(J1-1)+1)+AINV(3,2)*COORDS(3*(J1-1)+2)+AINV(3,3)*COORDS(3*(J1-1)+3)
               WRITE(8,'(3F20.10)') XTEMP,YTEMP,ZTEMP
               ! WRITE(*,'(3F20.10)') XTEMP,YTEMP,ZTEMP
            ENDDO

            DO
               READ(7,'(A80)') FNAME           !Reads initial velocites (0. 0. 0.) at end of file)
               IF (FNAME(1:80) .EQ. " ") THEN
                  WRITE(8,'(A80)')
               ELSE
                  WRITE(8,'(A80)') FNAME
                  DO J1=1, NATOMS-1
                     READ(7,'(A80)') FNAME
                     WRITE(8,'(A80)') FNAME
                  ENDDO
                  EXIT
               END IF
            ENDDO
            CLOSE(7)
            CLOSE(8)

            CALL SYSTEM(' cp POSCAR POSCAR.old' )
            CALL SYSTEM(' mv OUTCAR OUTCAR.old')
            CALL SYSTEM(' mv POSCAR.new POSCAR')
            IF (DEBUG) WRITE(*,'(A)') ' potential> Calling VASP energy and gradient'
            CALL SYSTEM(VASPJOB)


            ! Searches OUTCAR for 'energy without entropy' string to get energies
            FNAME='OUTCAR'
            CALL SYSTEM (' grep "energy without entropy" OUTCAR | grep energy | tail -1 > temp_eng.outcar1')
            CALL SYSTEM(' sed -e "s/^\([^=]*=\)\{2\}//" temp_eng.outcar1 | sed -e "s/^[ \t]*//" > temp_eng.outcar2 ')
            OPEN(UNIT=7,FILE='temp_eng.outcar2',STATUS='OLD')
            READ(7,*) ENERGY1
            ENERGY=ENERGY1
            CLOSE(7)
            CALL SYSTEM ('rm temp_eng.outcar1')
            CALL SYSTEM ('rm temp_eng.outcar2')
            ! CALL SYSTEM('grep "Elapsed time" OUTCAR | sed -e "s/^\([^:]*:\)\{1\}//" | sed -e "s/^[ \t]*//" > temp_time.outcar1 ')
            ! OPEN (UNIT=7,FILE='temp_time.outcar1',STATUS='OLD')
            ! READ(7,'(A)') !FNAME
            ! CLOSE(7)
            WRITE(*,'(A,A,A,F20.10)') ' potential> VASP ',TRIM(FNAME),' Energies=',ENERGY1

            ! Searches POSCAR for 'POSITION' string to get forces, updates into vector VNEW
            CALL SYSTEM ('sed -e "1,/POSITION/d" OUTCAR | sed "1,1d" > test_force.outcar1 ')
            OPEN(UNIT=7,FILE='test_force.outcar1',STATUS='OLD')
            DO J1=1,NATOMS
               READ(7,*) DUMMY1, DUMMY2, DUMMY3, VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ! May not need this bit, copied from castep!
               VNEW(3*(J1-1)+1)=-VNEW(3*(J1-1)+1)
               VNEW(3*(J1-1)+2)=-VNEW(3*(J1-1)+2)
               VNEW(3*(J1-1)+3)=-VNEW(3*(J1-1)+3)
               ! WRITE (*,'(3F15.8)') VNEW(3*(J1-1)+1), VNEW(3*(J1-1)+2), VNEW(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            CALL SYSTEM(' rm test_force.outcar1 ')


            ! 
            ! Implementation of REAXFF interface here:
         ELSE IF (REAXFFT) THEN
            ! copy input file here but replace coordinates with optimised ones
            FNAME=SYS(1:LSYS)
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            OPEN(UNIT=8,FILE='reax_geo.new',STATUS='UNKNOWN')
            DO
               READ(7,'(A80)',END=835) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:6).EQ.'HETATM') THEN
                  BACKSPACE(7)
                  DO J1=1,NATOMS
                     READ(7,*) FROMFILE(1:6),MM,FROMFILE(7:8),
     1               DUMMY1, DUMMY2, DUMMY3,FROMFILE(9:10),
     1               IDUM1,IDUM2,DOUBLEDUM1
                     ! PRINT*,FROMFILE(1:6),MM,FROMFILE(7:8),
                     ! 1               COORDS(3*J1-2),COORDS(3*J1-1),COORDS(3*J1),FROMFILE(9:10),
                     ! 1              IDUM1,IDUM2,DOUBLEDUM1
                     WRITE(8,'(a6,1x,i5,1x,a2,3x,1x,3x,1x,1x,1x,5x,3f20.15,1x,a5,i3,i2,1x,f8.5)')
     1               FROMFILE(1:6),MM,FROMFILE(7:8),
     1               COORDS(3*J1-2),COORDS(3*J1-1),COORDS(3*J1),FROMFILE(9:10),
     1               IDUM1,IDUM2,DOUBLEDUM1
                  ENDDO
               ELSE
                  WRITE(8,'(A80)') FNAME
               ENDIF
            ENDDO

835         CONTINUE
            CLOSE(7)
            CLOSE(8)

            ! Fiddle with input and output files and execute ReaxFF

            CALL SYSTEM(' cp ' // SYS(1:LSYS) // ' ' // SYS(1:LSYS) // '.old' )
            ! CALL SYSTEM(' cp ' // SYS(1:(LSYS-4)) // '.out' // ' ' // SYS(1:(LSYS-4)) // '.old' ) Not sure if this works!
            CALL SYSTEM(' rm ' // SYS(1:LSYS))
            CALL SYSTEM(' mv reax_geo.new ' // SYS(1:LSYS) )
            IF (DEBUG) WRITE(*,'(A)') ' potential> Calling ReaxFF energy and gradient'
            WRITE(STRING,*) TRIM(ADJUSTL(REAXFFJOB)) !REAXFFJOB is executable
            IF (DEBUG) WRITE(*,'(A,A)') ' potential> ',STRING
            CALL SYSTEM(STRING)

            ! read out energy from ReaxFF output energies aus fort.57 und grad geom nach xmolout
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='fort.57',STATUS='OLD')
            DO
               READ(LUNIT,'(A80)',END=813) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:6).EQ.'     0') THEN
                  FNAME2=FNAME2(12:29)
                  READ(FNAME2,'(F20.10)') ENERGY1
               ENDIF
            ENDDO
813         CONTINUE
            CLOSE(LUNIT)
            ENERGY1=ENERGY1/23.06054D0 !conversion from kcal to eV (taken from NIST)
            WRITE(*,'(A,2F20.10)')' Energy=',ENERGY1
            ENERGY=ENERGY1

            ! read gradients and freeze them
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='xmolout',STATUS='OLD') !assuming that the output is called like this!

            READ(LUNIT,'(A80)',END=812) FNAME
            READ(LUNIT,'(A80)',END=812) FNAME
            DO J1=1,NATOMS
               READ(LUNIT,*) STRING,XTEMP,YTEMP,ZTEMP,
     1         VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3),NDUMMY
               IF (FROZEN(J1)) THEN
                  IF (VARIABLES) THEN
                     VNEW(J1)=0.D0
                  ELSE
                     VNEW(3*(J1-1)+1)=0.D0
                     VNEW(3*(J1-1)+2)=0.D0
                     VNEW(3*(J1-1)+3)=0.D0
                  ENDIF
               ENDIF
            ENDDO

812         CONTINUE
            CLOSE(LUNIT)
            VNEW(1:3*NATOMS)=VNEW(1:3*NATOMS)/23.06054D0   !kcal/Angstrom to eV/Angstrom.
            ! WRITE(*,*) 'VNEW: ',VNEW



            ! 
            ! Implementation of MOLPRO interface here:
         ELSE IF (MOLPRO) THEN
            ! copy input file here but replace coordinates with optimised ones
            FNAME=SYS(1:LSYS)
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            OPEN(UNIT=8,FILE='molpro_input.new',STATUS='UNKNOWN')
            DO
               READ(7,'(A80)',END=815) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:9).EQ.' GEOMETRY') THEN
                  WRITE(8,'(A)') FNAME
                  DO J1=1,NATOMS
                     READ(7,*) FNAME
                     WRITE(8,'(A2,3F20.10)') FNAME(1:2),COORDS(3*(J1-1)+1),
     &               COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3)
                     WRITE(*,'(A2,3F20.10)') FNAME(1:2),COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3)
                  ENDDO
                  READ(7,'(A24)') FNAME
                  IF (FNAME(1:1).NE.'}') THEN
                     PRINT '(3A)',' potential> ERROR - string ',FNAME(1:1),' should be }'
                  ELSE
                     WRITE(8,'(A24)') FNAME
                  ENDIF
               ELSE
                  WRITE(8,'(A80)') FNAME
               ENDIF
            ENDDO

            ! Fiddle with input and output files and execute MOLPRO
815         CONTINUE
            CLOSE(7)
            CLOSE(8)
            CALL SYSTEM(' cp ' // SYS(1:LSYS) // ' ' // SYS(1:LSYS) // '.old' )
            ! CALL SYSTEM(' cp ' // SYS(1:(LSYS-4)) // '.out' // ' ' // SYS(1:(LSYS-4)) // '.old' ) Not sure if this works!
            CALL SYSTEM(' rm ' // SYS(1:LSYS))
            CALL SYSTEM(' mv molpro_input.new ' // SYS(1:LSYS) )
            IF (DEBUG) WRITE(*,'(A)') ' potential> Calling MOLPRO energy and gradient'
            CALL SYSTEM('rm input.out* input.xml* input.log*') !removing last run's temporary/output files
            WRITE(STRING,*) TRIM(ADJUSTL(MOLPROJOB)) // ' ' // TRIM(ADJUSTL(MOLPROJOBPARAMS)) !MOLPROJOB is executable, ~PARAMS varies
            IF (DEBUG) WRITE(*,'(A,A)') ' potential> ',STRING
            CALL SYSTEM(STRING)

            ! read out energy from MOLPRO output
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='input.out',STATUS='OLD')
            DO
               READ(LUNIT,'(A80)',END=816) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:13).EQ.' TOTAL ENERGY') THEN
                  FNAME2=FNAME2(26:36)
                  READ(FNAME2,'(F20.10)') ENERGY1
               ENDIF
            ENDDO
816         CONTINUE
            CLOSE(LUNIT)
            ENERGY1=ENERGY1/3.67493237981D-2 !conversion from E_H to eV
            WRITE(*,'(A,2F20.10)')' Energy=',ENERGY1
            ENERGY=ENERGY1

            ! read gradients
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='input.out',STATUS='OLD') !assuming that the output is called like this!
            DO
               READ(LUNIT,'(A80)',END=817) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:13).EQ.' TOTAL ENERGY') THEN
                  READ(LUNIT,'(A80)',END=817) FNAME
                  READ(LUNIT,'(A80)',END=817) FNAME
                  DO J1=1,NATOMS
                     READ(LUNIT,*) NDUMMY,VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
                  ENDDO
               ENDIF
            ENDDO
            WRITE(*,*) 'VNEW: ',VNEW
817         CONTINUE
            VNEW(1:3*NATOMS)=VNEW(1:3*NATOMS)/3.67493237981D-2 !conversion from E_H to eV
            ! VNEW(1:3*NATOMS)=VNEW(1:3*NATOMS)/MOLPROSCALE ! conversion from a0 to AA, should not be required here


            ! 



         ELSE IF (QCHEM) THEN
            IF (.NOT.QCHEMES) THEN
               FNAME= SYS(1:LSYS)
               IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
               FNAME= SYS(1:LSYS)
               OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
               OPEN(UNIT=8,FILE='qchem.new',STATUS='UNKNOWN')
               DO
                  READ(7,'(A80)',END=891) FNAME
                  FNAME2=FNAME
                  CALL UPPERCASE(FNAME2)
                  IF (FNAME2(1:4).EQ.'$MOL') THEN
                     WRITE(8,'(A21)') FNAME
                     READ(7,'(A80)',END=891) FNAME
                     WRITE(8,'(A80)') FNAME
                     DO J1=1,NATOMS
                        READ(7,*) FNAME
                        WRITE(8,'(A2,3F20.10)') FNAME(1:2),COORDS(3*(J1-1)+1),
     &                  COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3)
                     ENDDO
                     READ(7,'(A24)') FNAME
                     FNAME2=FNAME
                     CALL UPPERCASE(FNAME2)
                     IF (FNAME2(1:4).NE.'$END') THEN
                        PRINT '(3A)',' potential> ERROR - string ',FNAME(1:4),' should be $END'
                     ELSE
                        WRITE(8,'(A24)') FNAME
                     ENDIF
                  ELSE
                     WRITE(8,'(A80)') FNAME
                  ENDIF
               ENDDO
891            CONTINUE
               CLOSE(7)
               CLOSE(8)
               CALL SYSTEM(' cp ' // SYS(1:LSYS) // ' ' // SYS(1:LSYS) // '.old' )
               CALL SYSTEM(' cp qchem.out qchem.out.old')
               CALL SYSTEM(' rm ' // SYS(1:LSYS))
               CALL SYSTEM(' mv qchem.new ' // SYS(1:LSYS) )
               IF (DEBUG) WRITE(*,'(A)') ' potential> Calling QCHEM energy and gradient'
               ! 
               ! This assignment didn't work - the WRITE statement does!
               ! 
               ! FNAME=QCHEMJOB // ' qchem.out'
               ! 
               WRITE(STRING,*) TRIM(ADJUSTL(QCHEMJOB)) // ' qchem.out '//TRIM(ADJUSTL(QCHEMJOBPARAMS))
               IF (DEBUG) WRITE(*,'(A,A)') ' potential> ',STRING
               CALL SYSTEM(STRING)

               CALL SYSTEM(' grep "Convergence criterion met" qchem.out > temp.qchem')
               CALL SYSTEM(' sed -e "s/.....//" temp.qchem > temp2.qchem ')
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='temp2.qchem',STATUS='OLD')
               READ(LUNIT,*) ENERGY1
               CLOSE(LUNIT)
               CALL SYSTEM('grep "Total job time:" qchem.out > temp')
               OPEN (UNIT=LUNIT,FILE='temp',STATUS='OLD')
               READ(LUNIT,'(A)') FNAME
               CLOSE(LUNIT)
               WRITE(*,'(A,A,A,2F20.10)') ' potential> QCHEM ',TRIM(FNAME),' Energy=',ENERGY1
               ENERGY=ENERGY1
               CALL SYSTEM(' sed -e "1,/Gradient of/d" qchem.out > edited.0')
               CALL SYSTEM(' sed -e "1,/Gradient of/d" edited.0 > edited.1')
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='edited.1',STATUS='OLD')
               ! IF (DEBUG) PRINT '(A)',' potential> QCHEM forces:'
               J4=0
               DO J1=1,NATOMS/6
                  READ(LUNIT,*) FNAME
                  READ(LUNIT,*) NDUMMY,(VNEW(3*(J4+J2-1)+1),J2=1,6)
                  READ(LUNIT,*) NDUMMY,(VNEW(3*(J4+J2-1)+2),J2=1,6)
                  READ(LUNIT,*) NDUMMY,(VNEW(3*(J4+J2-1)+3),J2=1,6)
                  J4=J4+6
               ENDDO
               IF (J4.LT.NATOMS) THEN
                  READ(LUNIT,*) FNAME
                  READ(LUNIT,*) NDUMMY,(VNEW(3*(J4+J2-1)+1),J2=1,NATOMS-J4)
                  READ(LUNIT,*) NDUMMY,(VNEW(3*(J4+J2-1)+2),J2=1,NATOMS-J4)
                  READ(LUNIT,*) NDUMMY,(VNEW(3*(J4+J2-1)+3),J2=1,NATOMS-J4)
               ENDIF
               VNEW(1:3*NATOMS)=VNEW(1:3*NATOMS)/QCHEMSCALE ! unit conversion if required
               ! 
               ! QCHEMES - a rather different sort of job!
               ! 
            ELSE
               ! 
               ! First need to run QCHEM. Not quite the same as above!
               ! Need to use the same geometry, but different MO coefficients
               ! and do a single point energy. Write the new coefficients,
               ! which are the OPTIM variables, and run QCHEM.
               ! Can copy some of this from above.
               ! 
               ! This block will read the new output. Formats as in fetchz.f
               ! 
               CALL SYSTEM(' sed -e "1,/Calculated Alpha Fock Matrix/d" qchem.out > edited.1')
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='edited.1',STATUS='OLD')
               READ(LUNIT,*) FNAME
               IF (DEBUG) PRINT '(A)',' potential> Reading Fock matrix:'
               IF (DEBUG) PRINT '(A)', FNAME
               ALLOCATE(FOCK(QCHEMESNAO,QCHEMESNAO))
               J4=0
               DO J1=1,QCHEMESNAO/6
                  DO J2=1,QCHEMESNAO
                     READ(LUNIT,*) NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
                     WRITE(*,'(I5,6F12.7)') NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
                  ENDDO
                  J4=J4+6
               ENDDO
               IF (J4.LT.QCHEMESNAO) THEN
                  DO J2=1,QCHEMESNAO
                     READ(LUNIT,*) NDUMMY,(FOCK(J3,J2),J3=J4+1,QCHEMESNAO)
                     WRITE(*,'(I5,6F12.7)') NDUMMY,(FOCK(J3,J2),J3=J4+1,QCHEMESNAO)
                  ENDDO
               ENDIF
               CLOSE(LUNIT)

               CALL SYSTEM(' sed -e "1,/Alpha MO Coefficients/d" qchem.out > edited.1')
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='edited.1',STATUS='OLD')
               IF (DEBUG) PRINT '(A)',' potential> Reading MO coefficients:'
               READ(LUNIT,*) FNAME
               ALLOCATE(MOCOEFF(QCHEMESNELEC,QCHEMESNAO))
               J4=0
               DO J1=1,QCHEMESNELEC/6
                  DO J2=1,QCHEMESNAO
                     READ(LUNIT,*) NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
                     WRITE(*,'(I5,6F12.7)') NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
                  ENDDO
                  J4=J4+6
               ENDDO
               IF (J4.LT.QCHEMESNELEC) THEN
                  DO J2=1,QCHEMESNAO
                     READ(LUNIT,*) NDUMMY,(MOCOEFF(J3,J2),J3=J4+1,QCHEMESNELEC)
                     WRITE(*,'(I5,6F12.7)') NDUMMY,(MOCOEFF(J3,J2),J3=J4+1,QCHEMESNELEC)
                  ENDDO
               ENDIF
               CLOSE(LUNIT)
               ! 
               ! Now we need a matrix multiply to get the actual gradient?
               ! 
               DEALLOCATE(FOCK,MOCOEFF)
               STOP
            ENDIF
            ! IF (DEBUG) WRITE(*,'(3G20.10)') VNEW(1:3*NATOMS)
            CLOSE(LUNIT)
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV/atom'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV/atom'
         ELSE IF (CASTEP) THEN
            IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
            FNAME= SYS(1:LSYS) // '.cell'
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            OPEN(UNIT=8,FILE='new.cell',STATUS='UNKNOWN')
            ! WRITE(*,'(A)') 'The absolute coordinates in potentential.f (before CASTEP)'
            ! WRITE(*,'(6F15.5)') (COORDS(J1),J1=1,NOPT)
            DO
               READ(7,'(A80)',END=831) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:16).EQ.'%BLOCK POSITIONS') THEN
                  WRITE(8,'(A21)') FNAME
                  ! IF ((FNAME2(18:20).EQ.'ABS').AND.DEBUG) PRINT '(A)',' potential> CASTEP absolute coordinates'
                  ! IF ((FNAME2(18:20).EQ.'FRA').AND.DEBUG) PRINT '(A)',' potential> CASTEP fractional coordinates'
                  DO J1=1,NATOMS
                     READ(7,*) FNAME
                     XTEMP=AINV(1,1)*COORDS(3*(J1-1)+1)+AINV(1,2)*COORDS(3*(J1-1)+2)+AINV(1,3)*COORDS(3*(J1-1)+3)
                     YTEMP=AINV(2,1)*COORDS(3*(J1-1)+1)+AINV(2,2)*COORDS(3*(J1-1)+2)+AINV(2,3)*COORDS(3*(J1-1)+3)
                     ZTEMP=AINV(3,1)*COORDS(3*(J1-1)+1)+AINV(3,2)*COORDS(3*(J1-1)+2)+AINV(3,3)*COORDS(3*(J1-1)+3)
                     WRITE(8,'(A2,3F20.10)') FNAME(1:2),XTEMP,YTEMP,ZTEMP
                     ! WRITE(*,'(A2,3F20.10)') FNAME(1:2),XTEMP,YTEMP,ZTEMP
                  ENDDO
                  READ(7,'(A24)') FNAME
                  FNAME2=FNAME
                  CALL UPPERCASE(FNAME2)
                  IF (FNAME2(1:19).NE.'%ENDBLOCK POSITIONS') THEN
                     PRINT '(3A)',' potential> ERROR - string ',FNAME(1:19),' should be %ENDBLOCK POSITIONS'
                  ELSE
                     WRITE(8,'(A24)') FNAME
                  ENDIF
               ELSE
                  WRITE(8,'(A80)') FNAME
               ENDIF
            ENDDO
831         CONTINUE
            CLOSE(7)
            CLOSE(8)
            ! CALL SYSTEM(' cat ' // SYS(1:LSYS) // '.cell >>' // SYS(1:LSYS) // '.cell.old' )
            ! CALL SYSTEM(' cat ' // SYS(1:LSYS) // '.castep >>' // SYS(1:LSYS) // '.castep.old' )
            CALL SYSTEM(' cat ' // SYS(1:LSYS) // '.cell >' // SYS(1:LSYS) // '.cell.old' )
            CALL SYSTEM(' cat ' // SYS(1:LSYS) // '.castep >' // SYS(1:LSYS) // '.castep.old' )
            CALL SYSTEM(' rm ' // SYS(1:LSYS) // '.castep ' )
            CALL SYSTEM(' mv new.cell ' // SYS(1:LSYS) // '.cell' )
            ! WRITE(*,'(A)') ' potential> Calling CASTEP energy and gradient'

            ! 
            ! Change job submission to use a string following the CASTEP keyword in odata
            ! The PARALLEL keyword will no longer be needed if the number of processors
            ! is specified here instead.
            ! The examples below are retained for reference.
            ! 
            CALL SYSTEM(CASTEPJOB)

            ! IF (PARALLEL) THEN
            ! CALL SYSTEM(' ( mpirun -np ' // NPROC // ' /export/home/wales/bin/castepexe.new ' // SYS(1:LSYS) // ' ) ')
            ! CALL SYSTEM(' ( mpirun -np ' // NPROC // ' /export/home/wales/bin/castepexe.new ' // SYS(1:LSYS) // ' ) >& /dev/null')
            ! CALL SYSTEM(' ( mprun -n -np ' // NPROC // ' /export/home/wales/bin/castepexe.new ' // SYS(1:LSYS) // ' ) ')
            ! CALL SYSTEM(' ( lamwrapper  /home/wales/bin/castep4.1.mpi ' // SYS(1:LSYS) // ' ) ')
            ! CALL SYSTEM(' ( scrunwrapper  /home/wales/bin/castepexe.new ' // SYS(1:LSYS) // ' ) ')
            ! CALL SYSTEM(' ( mpichwrapper  /home/wales/bin/castep4.1.mpi ' // SYS(1:LSYS) // ' ) ')
            ! CALL SYSTEM(' ( mpichwrapper  /home/wales/bin/castep.mpi ' // SYS(1:LSYS) // ' ) ')
            ! CALL SYSTEM(' ( mpirun  /home/wales/bin/castep ' // SYS(1:LSYS) // ' ) ')
            ! 
            ! Next version for SiCortex
            ! 
            ! CALL SYSTEM(' srun -p sca -n ' // NPROC // ' /home/wales/bin/castep ' // SYS(1:LSYS) )
            ! 
            ! Next version for darwin.
            ! The sleep 10 line seems to be needed on darwin where one job can start
            ! before the previous one has been cleaned up.
            ! 
            ! CALL SYSTEM(' ( mpirun  -np '//NPROC//' -machinefile machine.file /home/dw34/bin/castep.mpi '//SYS(1:LSYS)//' ) ')
            ! CALL SYSTEM('sleep 10') ! new DJW
            ! CALL SYSTEM('mpiexec -comm none killall -9 castep.mpi') ! new DJW
            ! ELSE
            ! CALL SYSTEM(' ( castep4.1 ' // SYS(1:LSYS) // ' ) >& /dev/null')
            ! CALL SYSTEM(' ( castep4.1 ' // SYS(1:LSYS) // ' ) ')
            ! ENDIF

            CALL SYSTEM(' grep "Final" ' // SYS(1:LSYS) // '.castep | grep energy | tail -1 > temp.castep1')
            CALL SYSTEM(' grep "Final energy" ' // SYS(1:LSYS) // '.castep | grep energy | tail -1 > temp.castep2')
            CALL SYSTEM(' sed -e "s/[a-zA-Z]//g" -e "s/=//" -e "s/,//" -e "s/(.*)//" temp.castep1 > temp2.castep1 ')
            CALL SYSTEM(' sed -e "s/[a-zA-Z]//g" -e "s/=//" -e "s/,//" -e "s/(.*)//" temp.castep2 > temp2.castep2 ')
            OPEN(UNIT=7,FILE='temp2.castep1',STATUS='OLD')
            READ(7,*) ENERGY1
            CLOSE(7)
            OPEN(UNIT=7,FILE='temp2.castep2',STATUS='OLD')
            READ(7,*) ENERGY2
            CLOSE(7)
            CALL SYSTEM('grep "Total time" ' // SYS(1:LSYS) // '.castep | tail -1 | sed -e "s/  */ /g" > temp')
            OPEN (UNIT=7,FILE='temp',STATUS='OLD')
            READ(7,'(A)') FNAME
            CLOSE(7)
            WRITE(*,'(A,A,A,2F20.10)') ' potential> CASTEP ',TRIM(FNAME),' Energies=',ENERGY1, ENERGY2
            ENERGY=ENERGY1
            ! 
            ! Note that CASTEP4.1 reorders the atoms! Looks like they have to be in order of
            ! increasing atomic number.
            ! 
            CALL SYSTEM(' sed -e "1,/Force/d" ' // SYS(1:LSYS) // '.castep > edited.1')
            ! CALL SYSTEM(' sed -e "1,/x/d" -e "s/.........//" -e "s/\*//g" edited.1 > edited.castep')
            ! 
            ! Changed to deal with new output format for Force printing in castep 4.4 DJW 27/7/10
            ! 
            CALL SYSTEM(' sed -e "1,/x/d" -e "s/.....//" -e "s/\*//g" edited.1 > edited.castep')

            OPEN(UNIT=7,FILE='edited.castep',STATUS='OLD')
            ! IF (DEBUG) PRINT '(A)',' potential> CASTEP forces:'
            DO J1=1,NATOMS
               ! 
               ! Changed to deal with new output format for Force printing in castep 4.4 DJW 27/7/10
               ! 
               ! READ(7,*) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               READ(7,*) NDUMMY,VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ! IF (DEBUG) WRITE(*,'(3F20.10)') VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               VNEW(3*(J1-1)+1)=-VNEW(3*(J1-1)+1)
               VNEW(3*(J1-1)+2)=-VNEW(3*(J1-1)+2)
               VNEW(3*(J1-1)+3)=-VNEW(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            ! 
            ! Do we need to project out overall rotation and translation? Only if it s a
            ! cluster.
            ! 

            ! el316
            ! IF (CASTEPC) CALL ORTHOGOPT(VNEW,COORDS,.FALSE.)
            ! el316
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV/atom'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV/atom'
            ENDIF
            CALL SYSTEM(' grep -c "elec_restore_file" ' // SYS(1:LSYS) // '.param > temp')
            OPEN(UNIT=7,FILE='temp',STATUS='OLD')
            READ(7,*) ISTART
            CLOSE(7)
            INQUIRE(FILE=SYS(1:LSYS) // '.wvfn.1',EXIST=YESNO)
            FNAME='echo elec_restore_file  :  ' // SYS(1:LSYS) // '.wvfn >> ' // SYS(1:LSYS) // '.param' !the .1 should be excluded
            IF (YESNO.AND.(ISTART.EQ.0)) CALL SYSTEM(FNAME)
         ELSE IF (ONETEP) THEN
            IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
            FNAME= SYS(1:LSYS) // '.dat'
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            OPEN(UNIT=8,FILE='new.dat',STATUS='UNKNOWN')
            DO
               READ(7,'(A80)',END=889) FNAME
               FNAME2=FNAME
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:16).EQ.'%BLOCK POSITIONS') THEN
                  WRITE(8,'(A21)') FNAME
                  ! IF ((FNAME2(18:20).EQ.'ABS').AND.DEBUG) PRINT '(A)',' potential> ONETEP absolute coordinates'
                  ! IF ((FNAME2(18:20).EQ.'FRA').AND.DEBUG) PRINT '(A)',' potential> ONETEP fractional coordinates'
                  DO J1=1,NATOMS
                     READ(7,*) FNAME
                     XTEMP=AINV(1,1)*COORDS(3*(J1-1)+1)+AINV(1,2)*COORDS(3*(J1-1)+2)+AINV(1,3)*COORDS(3*(J1-1)+3)
                     YTEMP=AINV(2,1)*COORDS(3*(J1-1)+1)+AINV(2,2)*COORDS(3*(J1-1)+2)+AINV(2,3)*COORDS(3*(J1-1)+3)
                     ZTEMP=AINV(3,1)*COORDS(3*(J1-1)+1)+AINV(3,2)*COORDS(3*(J1-1)+2)+AINV(3,3)*COORDS(3*(J1-1)+3)
                     WRITE(8,'(A2,3F20.10)') FNAME(1:2),XTEMP,YTEMP,ZTEMP
                     ! WRITE(*,'(A2,3F20.10)') FNAME(1:2),XTEMP,YTEMP,ZTEMP
                  ENDDO
                  READ(7,'(A24)') FNAME
                  FNAME2=FNAME
                  CALL UPPERCASE(FNAME2)
                  IF (FNAME2(1:19).NE.'%ENDBLOCK POSITIONS') THEN
                     PRINT '(3A)',' potential> ERROR - string ',FNAME(1:19),' should be %ENDBLOCK POSITIONS'
                  ELSE
                     WRITE(8,'(A24)') FNAME
                  ENDIF
               ELSE
                  WRITE(8,'(A80)') FNAME
               ENDIF
            ENDDO
889         CONTINUE
            CLOSE(7)
            CLOSE(8)
            CALL SYSTEM(' cat ' // SYS(1:LSYS) // '.dat >' // SYS(1:LSYS) // '.dat.old' )
            CALL SYSTEM(' cat ' // SYS(1:LSYS) // '.onetep >' // SYS(1:LSYS) // '.onetep.old' )
            CALL SYSTEM(' rm ' // SYS(1:LSYS) // '.onetep ' )
            CALL SYSTEM(' mv new.dat ' // SYS(1:LSYS) // '.dat' )
            ! WRITE(*,'(A)') ' potential> Calling ONETEP energy and gradient'

            CALL SYSTEM(ONETEPJOB)

            ! IF (PARALLEL) THEN
            ! CALL SYSTEM(' ( mpirun  /home/wales/bin/onetep ' // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.onetep ) ')
            ! 
            ! Next version for SiCortex
            ! 
            ! CALL SYSTEM(' srun -p sca -n ' // NPROC // ' /home/wales/bin/onetep ' // SYS(1:LSYS) )
            ! 
            ! ELSE
            ! CALL SYSTEM(' ( mpirun  /home/wales/bin/onetep ' // SYS(1:LSYS) // ' >& ' // SYS(1:LSYS) // '.onetep ) ')
            ! ENDIF
            ! CALL SYSTEM(' grep "total_energy" ' // SYS(1:LSYS) // '.onetep > temp.onetep')
            CALL SYSTEM(' grep "Total" ' // SYS(1:LSYS) // '.onetep | tail -1 > temp.onetep')
            CALL SYSTEM(' sed -e "s/.*://" temp.onetep > temp2.onetep ')
            OPEN(UNIT=7,FILE='temp2.onetep',STATUS='OLD')
            READ(7,*) ENERGY1
            CLOSE(7)
            CALL SYSTEM('grep "TOTAL TIME" ' // SYS(1:LSYS) // '.onetep | sed -e "s/  */ /g" > temp')
            OPEN (UNIT=7,FILE='temp',STATUS='OLD')
            READ(7,'(A)') FNAME
            CLOSE(7)
            WRITE(*,'(A,A,A,2F20.10)') ' potential> ONETEP ',TRIM(FNAME),' Energy=',ENERGY1
            ENERGY=ENERGY1
            ! 
            ! ONETEP does not seem to reorder the atoms.
            ! 
            CALL SYSTEM(' sed -e "1,/Forces/d" ' // SYS(1:LSYS) // '.onetep > edited.1')
            CALL SYSTEM(' sed -e "1,/x/d" -e "s/.........//" -e "s/\*//g" edited.1 > edited.onetep')

            OPEN(UNIT=7,FILE='edited.onetep',STATUS='OLD')
            ! IF (DEBUG) PRINT '(A)',' potential> ONETEP coords:'
            ! 
            ! It seems impossible to change force output from hartree/angstrom units, so need
            ! to change to hartree/bohr here.
            ! 
            DO J1=1,NATOMS
               READ(7,*) NDUMMY,VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ! IF (DEBUG) WRITE(*,'(3F20.10)') VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ! IF (DEBUG) WRITE(*,'(3F20.10)') COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3)
               VNEW(3*(J1-1)+1)=-VNEW(3*(J1-1)+1)*0.529177249D0
               VNEW(3*(J1-1)+2)=-VNEW(3*(J1-1)+2)*0.529177249D0
               VNEW(3*(J1-1)+3)=-VNEW(3*(J1-1)+3)*0.529177249D0
            ENDDO
            CLOSE(7)
            ! 
            ! Do we need to project out overall rotation and translation?
            ! 
            ! IF (ONETEPC) CALL ORTHOGOPT(VNEW,COORDS,.FALSE.)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV/atom'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV/atom'
            ENDIF
            INQUIRE(FILE=SYS(1:LSYS) // '.denskern',EXIST=YESNO)
            CALL SYSTEM(' grep -c "read_denskern" ' // SYS(1:LSYS) // '.dat > temp')
            OPEN(UNIT=7,FILE='temp',STATUS='OLD')
            READ(7,*) ISTART
            CLOSE(7)
            FNAME='echo "read_denskern       : TRUE" >> ' // SYS(1:LSYS) // '.dat'
            IF (YESNO.AND.(ISTART.EQ.0)) CALL SYSTEM(FNAME)

            ! INQUIRE(FILE=SYS(1:LSYS) // '.tightbox_ngwfs',EXIST=YESNO)
            ! CALL SYSTEM(' grep -c "read_tightbox_ngwfs" ' // SYS(1:LSYS) // '.dat > temp')
            ! OPEN(UNIT=7,FILE='temp',STATUS='OLD')
            ! READ(7,*) ISTART
            ! CLOSE(7)
            ! FNAME='echo "read_tightbox_ngwfs       : TRUE" >> ' // SYS(1:LSYS) // '.dat'
            ! IF (YESNO.AND.(ISTART.EQ.0)) CALL SYSTEM(FNAME)

         ELSE IF (CP2K) THEN
            IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
            FNAME= SYS(1:LSYS) // '.inp'
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            OPEN(UNIT=8,FILE='driver.inp',STATUS='UNKNOWN')
            OPEN(UNIT=9,FILE='trajectory.xyz',POSITION='APPEND')
            DO
               READ(7,'(A80)',END=890) FNAME
               FNAME2=FNAME
               FNAME2=ADJUSTL(FNAME2)
               CALL UPPERCASE(FNAME2)
               IF (FNAME2(1:6).EQ.'&COORD') THEN
                  WRITE(8,'(A80)') FNAME
                  WRITE(9,*) NATOMS
                  READ(7,'(A80)') FNAME ! Reading and writing the setting of scaled away
                  WRITE(8,'(A80)') FNAME
                  WRITE(9,*) "    "
                  DO J1=1,NATOMS
                     READ(7,*) FNAME
                     XTEMP=AINV(1,1)*COORDS(3*(J1-1)+1)+AINV(1,2)*COORDS(3*(J1-1)+2)+AINV(1,3)*COORDS(3*(J1-1)+3)
                     YTEMP=AINV(2,1)*COORDS(3*(J1-1)+1)+AINV(2,2)*COORDS(3*(J1-1)+2)+AINV(2,3)*COORDS(3*(J1-1)+3)
                     ZTEMP=AINV(3,1)*COORDS(3*(J1-1)+1)+AINV(3,2)*COORDS(3*(J1-1)+2)+AINV(3,3)*COORDS(3*(J1-1)+3)
                     WRITE(8,'(A2,3F20.10)') FNAME(1:2),XTEMP,YTEMP,ZTEMP
                     WRITE(9,'(A2,3F20.10)') FNAME(1:2),COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3)
                  ENDDO
                  READ(7,'(A24)') FNAME
                  FNAME2=FNAME
                  FNAME2=ADJUSTL(FNAME2)
                  CALL UPPERCASE(FNAME2)
                  IF (FNAME2(1:10).NE.'&END COORD') THEN
                     PRINT '(3A)',' potential> ERROR - string ',FNAME(1:10),' should be &END COORD'
                  ELSE
                     WRITE(8,'(A24)') FNAME
                  ENDIF
               ELSE
                  WRITE(8,'(A80)') FNAME
               ENDIF
            ENDDO
890         CONTINUE
            CLOSE(7)
            CLOSE(8)
            CLOSE(9)
            CALL SYSTEM(' cat ' // 'driver.inp >>' // SYS(1:LSYS) // '.inp.old' )
            WRITE(*,'(A)') ' potential> Calling CP2K energy and gradient'

            CALL SYSTEM(' ( ' // CP2KJOB // ' > ' // SYS(1:LSYS) // '.out ) ' )

            ! 
            ! Mek-quake
            ! 
            ! IF (PARALLEL) THEN
            ! FNAME='/usr/local/openmpi/1.3/intel101/bin/'
            ! FNAME2=' /home/el316/cp2koptim2/exe/Linux-x86-64-intel/'
            ! CALL SYSTEM(' ( ' // FNAME(1:36) //'mpirun ' // FNAME2(1:47) // 'cp2kdriver' // ' > ' // SYS(1:LSYS) // '.out ) ' )
            ! ELSE
            ! FNAME2=' /home/el316/cp2koptim/exe/Linux-x86-64-intel/'
            ! CALL SYSTEM(' ( ' // FNAME2(1:46) // 'cp2kdriver ' // ' > ' // SYS(1:LSYS) // '.out ) ' )
            ! ENDIF

            ! 
            ! Clust
            ! 
            ! IF (PARALLEL) THEN
            ! FNAME='/usr/local/openmpi-1.2.6-intel-10.0/bin/'
            ! FNAME2=' /home/el316/cp2koptim2/exe/Linux-x86-64-intel/'
            ! CALL SYSTEM(FNAME(1:40)//'mpirun -np '//NPROC// FNAME2(1:47)//'cp2kdriver'//' > '//SYS(1:LSYS)//'.out')
            ! ELSE
            ! CALL SYSTEM(' ( /home/el316/cp2koptim/exe/Linux-x86-64-intel/cp2kdriver ' // ' > ' // SYS(1:LSYS) // '.out ) ')
            ! ENDIF

            CALL SYSTEM(' cat ' // 'driver.out >>' // SYS(1:LSYS) // '.cp2k.old' )
            CALL SYSTEM(' rm ' // 'driver.inp' )
            CALL SYSTEM(' rm ' // 'driver.out' )
            CALL SYSTEM(' rm ' // 'dump.out' )
            FNAME='coordenergrad.out'
            OPEN(UNIT=7,FILE=FNAME,STATUS='UNKNOWN')
            READ(7,'(A)') FNAME
            DO J1=1,NATOMS
               READ(7,*) FNAME,FNAME,FNAME
            ENDDO
            READ(7,'(A)') FNAME
            READ(7,*) ENERGY
            READ(7,'(A)') FNAME
            DO J1=1,NATOMS
               READ(7,*) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            CALL SYSTEM(' cat coordenergrad.out >> ' // SYS(1:LSYS) // '.out.old' )
            CALL SYSTEM(' rm coordenergrad.out' )

            DO J1=1,3*NATOMS
               VNEW(J1)=-VNEW(J1)/0.529177249D0
            ENDDO

         ELSE IF (CPMD) THEN
            IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
            INQUIRE(FILE='RESTART.1',EXIST=YESNO)
            ! IF (NPCALL.GT.1) THEN
            OPEN(UNIT=8,FILE='newgeom',STATUS='UNKNOWN')
            DO J1=1,NATOMS
               WRITE(8,'(6F20.10)') COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3),0.0D0,0.0D0,0.0D0
            ENDDO
            CLOSE(8)
            CALL SYSTEM(' mv newgeom GEOMETRY ')
            IF ((NPCALL.EQ.1).OR.((NPCALL.EQ.0).AND.YESNO)) THEN
               IF (PRESSURE) THEN
                  NCOUNT=0
                  FNAME=SYS(1:LSYS)
                  OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
111               READ(7,'(A)',ERR=561) STRING
                  NCOUNT=NCOUNT+1
                  NDUM=INDEX(STRING,'CELL')
                  IF (NDUM.EQ.0) GOTO 111
                  READ(7,*) BOXLX
                  CLOSE(7)
                  ! CALL SYSTEM(' grep -c ANGSTROM ' // SYS(1:LSYS) // ' > temp')
                  ! OPEN(UNIT=7,FILE='temp',STATUS='OLD')
                  ! READ(7,*) J1
                  ! CLOSE(7)
                  ! IF (J1.EQ.1) THEN
                  ! WRITE(*,'(A)') ' Converting cell size from Angstrom to Bohr'
                  ! BOXLX=BOXLX*1.889726164D0
                  ! ENDIF
                  GOTO 567
561               WRITE(*,'(A)') 'CELL not found in input data set - quit'
                  STOP
567               CONTINUE
               ENDIF
               CALL SYSTEM(' sed -e "s/DUMMY/RESTART WAVEFUNCTION GEOFILE COORDINATES LATEST/" ' //  SYS(1:LSYS) // ' > temp ')
               IF (PRESSURE) THEN
                  OPEN(UNIT=7,FILE='temp',STATUS='OLD')
                  OPEN(UNIT=8,FILE='temp2',STATUS='UNKNOWN')
                  DO J1=1,NCOUNT
                     READ(7,'(A)') STRING
                     WRITE(8,'(A)') STRING
                  ENDDO
                  WRITE(8,'(A)') '  CELLSIZE'
                  READ(7,'(A)') STRING
664               READ(7,'(A)',ERR=665) STRING
                  WRITE(8,'(A)') STRING
                  GOTO 664
665               CONTINUE
                  ! CALL SYSTEM(' sed -e "s/ANGSTROM/noangstrom/" temp2 > temp ')
                  CALL SYSTEM(' cp temp2 temp ')
               ENDIF
               CALL SYSTEM(' mv temp ' // SYS(1:LSYS) // '.restart')
            ENDIF

            IF (PRESSURE) THEN
               WRITE(*,'(A)') ' Box length optimization'
               CALL CPMDLATMIN(NATOMS,COORDS,ENERGY,VNEW,BOXLX)
               YESNO=.TRUE.
            ELSE
               ! WRITE(*,'(A)') ' Calling CPMD energy and gradient'
               CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out >& /dev/null ')
               IF (.NOT.YESNO) THEN
                  IF (SCORE_QUEUE) THEN
                     CALL SYSTEM(' ( /usr/local/bin/scrunwrapper ' // TRIM(CPMD_COMMAND) // ' '
     1               // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
                  ELSE IF (PARALLEL) THEN
                     CALL SYSTEM(' ( mprun -n -np '// NPROC // ' ' // TRIM(CPMD_COMMAND) // ' '
     1               // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
                  ELSE
                     CALL SYSTEM(' ( ' // TRIM(CPMD_COMMAND) // ' ' // SYS(1:LSYS) // ' > '
     &                                // SYS(1:LSYS) // '.out ) >& /dev/null')
                  ENDIF
               ELSE
                  IF (SCORE_QUEUE) THEN
                     CALL SYSTEM(' ( /usr/local/bin/scrunwrapper ' // TRIM(CPMD_COMMAND) // ' '
     1               // SYS(1:LSYS) // '.restart > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
                  ELSE IF (PARALLEL) THEN
                     CALL SYSTEM(' ( mprun -n -np ' // NPROC // ' ' // TRIM(CPMD_COMMAND) // ' '
     1               // SYS(1:LSYS) // '.restart > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
                  ELSE
                     CALL SYSTEM(' ( ' // TRIM(CPMD_COMMAND) // ' ' // SYS(1:LSYS) // '.restart > '
     1               // SYS(1:LSYS) // '.out ) >& /dev/null')
                  ENDIF
               ENDIF
               CLOSE(7)
               OPEN (UNIT=7,FILE='ENERGY',STATUS='OLD')
               READ(7,*) ENERGY, GEMAX
               CLOSE(7)
               CALL SYSTEM('grep "CPU TIME" ' // SYS(1:LSYS) //
     1         '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD/" -e "s/  */ /g" > temp')
               OPEN (UNIT=7,FILE='temp',STATUS='OLD')
               READ(7,'(A)') FNAME
               WRITE(*,'(A,A,F20.10,A,F20.10)') TRIM(FNAME),' potential> Energy=',ENERGY,' GEMAX=',GEMAX
               ! IF (GEMAX.GT.1.0D-5) THEN
               ! WRITE(*,'(A,G15.5,A)') 'WARNING, GEMAX=',GEMAX,' CPMD wavefunction convergence suspect'
               ! ENDIF
               OPEN(UNIT=7,FILE='GEOMETRY',STATUS='OLD')
               DO J1=1,NATOMS
                  ! READ(7,*) GEMAX,GEMAX,GEMAX,VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
                  READ(7,*) COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3),
     1            VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
                  ! WRITE(*,'(6F20.10)') COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3),
                  ! 1                              VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
                  VNEW(3*(J1-1)+1)=-VNEW(3*(J1-1)+1)
                  VNEW(3*(J1-1)+2)=-VNEW(3*(J1-1)+2)
                  VNEW(3*(J1-1)+3)=-VNEW(3*(J1-1)+3)
               ENDDO
               CLOSE(7)
            ENDIF
            ! 
            ! Do we need to project out overall rotation and translation? Only if it s a
            ! cluster.
            ! 
            IF (CPMDC) CALL ORTHOGOPT(VNEW,COORDS,.FALSE.)
            IF (PV) ENERGY=ENERGY+PRESS*PARAM1*PARAM2*PARAM3
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ENDIF
         ELSE IF (CADPAC) THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
            OPEN(UNIT=15,FILE='temppoints',STATUS='UNKNOWN')
            WRITE(15,'(3F20.10)') (COORDS(J1),J1=1,NOPT)
            CLOSE(15)
            ! INQUIRE(FILE='derivs',EXIST=YESNO)
            IF (SSTEST) THEN
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling CADPAC energy and first and second derivatives'
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints > ' // SYS(1:LSYS) // '.dat ) >& /dev/null')
                  CALL SYSTEM('echo SECDER >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo PUNCH FCM >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo START >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo FINISH >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('runcadpac ' // SYS(1:LSYS) // ' >& /dev/null')
               ELSE
                  WRITE(*,'(A)') ' Using existing CADPAC derivs file for first and second derivatives'
               ENDIF

               OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
               DO J1=1,3+NATOMS
                  READ(15,*,ERR=666)
               ENDDO
               READ(15,*,ERR=666) (VNEW(J1),J1=1,3*NATOMS)
               READ(15,*,ERR=666)
               READ(15,*,ERR=666) (((HESS(3*(J1-1)+J3,J2),J3=1,3),J2=1,3*NATOMS),J1=1,NATOMS)
               CLOSE(15)
            ELSE IF (GTEST) THEN
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling CADPAC energy and first derivatives'
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints > ' // SYS(1:LSYS) // '.dat ) >& /dev/null')
                  CALL SYSTEM('echo GRADIENT >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo PUNCH GRADIENT >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo START >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo FINISH >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('runcadpac ' // SYS(1:LSYS) // ' >& /dev/null')
               ELSE
                  WRITE(*,'(A)') ' Using existing CADPAC derivs file for energy and first derivatives'
               ENDIF
               OPEN(UNIT=15,FILE='derivs',STATUS='OLD',ERR=666)
               DO J1=1,4+NATOMS
                  READ(15,*,ERR=666)
               ENDDO
               READ(15,*,ERR=666) (VNEW(J1),J1=1,3*NATOMS)
            ELSE
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling CADPAC energy'
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints > ' // SYS(1:LSYS) // '.dat ) >& /dev/null')
                  CALL SYSTEM('echo START >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('echo FINISH >> ' // SYS(1:LSYS) // '.dat')
                  CALL SYSTEM('runcadpac ' // SYS(1:LSYS) // ' >& /dev/null')
               ELSE
                  WRITE(*,'(A)') ' Using existing CADPAC derivs file for energy'
               ENDIF
            ENDIF
            ! 
            ! There is no provision for punching the energy consistently in CADPAC.
            ! 
            CALL SYSTEM('grep "Final SCF" ' // SYS(1:LSYS) // '.out > temp')
            CALL SYSTEM('grep "Energy (RMP2)" ' // SYS(1:LSYS) // '.out >> temp')
            CALL SYSTEM('grep "Final DFT" ' // SYS(1:LSYS) // '.out | head -1 >> temp')
            CALL SYSTEM('tail -1 temp | sed -e "s/y.*-/-/"  -e "s/[a-zA-Z]//g" > abenergy')
            CALL SYSTEM('mv ' // SYS(1:LSYS) // '.dat ' // SYS(1:LSYS) // '.old.dat')
            CALL SYSTEM('mv ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out')
            IF (GTEST.OR.SSTEST) CALL SYSTEM('mv derivs derivs.old')
         ELSE IF (GAMESSUS) THEN
            OPEN(UNIT=15,FILE='temppoints',STATUS='UNKNOWN')
            WRITE(15,'(3F20.10)') (COORDS(J1),J1=1,NOPT)
            CLOSE(15)
            ! INQUIRE(FILE='derivs',EXIST=YESNO)
            IF (SSTEST) THEN
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' potential> Calling GAMESS-US energy and first and second derivatives'
                  CALL SYSTEM("echo ' $CONTRL RUNTYP=HESSIAN $END' > " // SYS(1:LSYS) // ".inp")
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints >> ' // SYS(1:LSYS) // '.inp ) >& /dev/null')
                  CALL SYSTEM('rungms ' // SYS(1:LSYS) // ' 01 ' // NPROC // ' >& ' // SYS(1:LSYS) // '.out ')
                  CALL SYSTEM(' cp ~/scr/' // TRIM(ADJUSTL(SYS(1:LSYS))) // '.dat derivs ')
               ELSE
                  WRITE(*,'(A)') ' Using existing GAMESS-US derivs file for first and second derivatives'
               ENDIF

               OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
122            READ(15,'(A80)') GSTRING
               IF (GSTRING(1:6).NE.' $GRAD') GOTO 122
               READ(15,'(A80)') GSTRING
               DO J1=1,NATOMS
                  READ(15,'(15X,3E20.10)',ERR=666) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ENDDO
               READ(15,*,ERR=666)
               READ(15,*,ERR=666)
               READ(15,*,ERR=666)
               DO J1=1,3*NATOMS
                  JSTART=1
123               READ(15,'(5X,5E15.8)',ERR=666) (HESS(J1,J3),J3=1+5*(JSTART-1),MIN(5+5*(JSTART-1),3*NATOMS))
                  IF (5*JSTART.LT.3*NATOMS) THEN
                     JSTART=JSTART+1
                     GOTO 123
                  ENDIF
                  ! READ(15,'(5X,5E15.8)',ERR=666) (HESS(J1,J3),J3=6,10)
                  ! READ(15,'(5X,5E15.8)',ERR=666) (HESS(J1,J3),J3=11,15)
                  ! READ(15,'(5X,3E15.8)',ERR=666) (HESS(J1,J3),J3=16,18)
               ENDDO
               CLOSE(15)
            ELSE IF (GTEST) THEN
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' potential> Calling GAMESS-US energy and first derivatives'
                  CALL SYSTEM("echo  ' $CONTRL RUNTYP=GRADIENT $END' > " // SYS(1:LSYS) // ".inp")
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints >> ' // SYS(1:LSYS) // '.inp ) ')
                  PRINT '(A)',' running: ' // 'rungms ' // SYS(1:LSYS) // ' 01 ' // NPROC // ' >& ' // SYS(1:LSYS) // '.out'
                  CALL SYSTEM('rungms ' // SYS(1:LSYS) // ' 01 ' // NPROC // ' >& ' // SYS(1:LSYS) // '.out')
                  CALL SYSTEM(' cp ~/scr/' // TRIM(ADJUSTL(SYS(1:LSYS))) // '.dat derivs ')
               ELSE
                  WRITE(*,'(A)') ' Using existing GAMESS-US derivs file for energy and first derivatives'
               ENDIF
               OPEN(UNIT=15,FILE='derivs',STATUS='OLD',ERR=666)
133            READ(15,'(A80)') GSTRING
               IF (GSTRING(1:6).NE.' $GRAD') GOTO 133
               READ(15,'(A80)') GSTRING
               DO J1=1,NATOMS
                  READ(15,'(15X,3E20.10)',ERR=666) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
                  ! WRITE(*,'(15X,3E20.10)') VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ENDDO
               CLOSE(15)
            ELSE
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling GAMESS-US energy'
                  CALL SYSTEM("echo  ' $CONTRL RUNTYP=ENERGY $END' > " // SYS(1:LSYS) // ".inp")
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints >> ' // SYS(1:LSYS) // '.inp ) >& /dev/null')
                  CALL SYSTEM('rungms ' // SYS(1:LSYS) // ' 01 ' // NPROC // ' >& ' // SYS(1:LSYS) // '.out')
                  CALL SYSTEM(' cp ~/scr/' // TRIM(ADJUSTL(SYS(1:LSYS))) // '.dat derivs ')
               ELSE
                  WRITE(*,'(A)') ' Using existing GAMESS-US derivs file for energy'
               ENDIF
            ENDIF
            ! 
            ! There is no provision for punching the energy consistently in GAMESS-US.
            ! 
            CALL SYSTEM('grep "FINAL" ' // SYS(1:LSYS) // '.out | grep ENERGY > temp')
            CALL SYSTEM('grep "E(MP2)" ' // SYS(1:LSYS) // '.out >> temp')
            CALL SYSTEM('grep "FINAL MCSCF ENERGY" ' // SYS(1:LSYS) // '.out >> temp')
            CALL SYSTEM('tail -1 temp | sed -e "s/.*-/-/"  -e "s/[a-zA-Z]//g" > abenergy')
            CALL SYSTEM('mv ' // SYS(1:LSYS) // '.inp ' // SYS(1:LSYS) // '.old.inp')
            CALL SYSTEM('mv ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out')
            IF (GTEST.OR.SSTEST) CALL SYSTEM('mv derivs derivs.old')
         ELSE IF (GAMESSUK) THEN
            ! 
            ! The system.in file will need the line:
            ! PUNCH SCFENERGY TRANSFORM GRADIENT SECDER
            ! 
            IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in potential'
            OPEN(UNIT=15,FILE='temppoints',STATUS='UNKNOWN')
            WRITE(15,'(3F20.10)') (COORDS(J1),J1=1,NOPT)
            CLOSE(15)
            ! INQUIRE(FILE='derivs',EXIST=YESNO)
            IF (SSTEST) THEN
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling GAMESS-UK energy and first and second derivatives'
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints > ' // SYS(1:LSYS) // '.in ) >& /dev/null')
                  CALL SYSTEM("echo 'runtype hessian' >> " // SYS(1:LSYS) // ".in")
                  CALL SYSTEM("echo 'enter' >> " // SYS(1:LSYS) // ".in")
                  CALL SYSTEM('rungamess ' // SYS(1:LSYS) // '>&/dev/null' )
                  ! CALL SYSTEM('rungamess ' // SYS(1:LSYS) // ' ' // NPROC // ' >& ' // SYS(1:LSYS) // '.out ')
                  CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.pun derivs ')
               ELSE
                  WRITE(*,'(A)') ' Using existing GAMESS derivs file for first and second derivatives'
               ENDIF

               OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
12             READ(15,'(A80)') GSTRING
               IF (GSTRING(1:17).NE.'block = gradients') GOTO 12
               DO J1=1,NATOMS
                  READ(15,*,ERR=666) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ENDDO
               READ(15,*) GSTRING
               READ(15,*) GSTRING
               READ(15,*) GSTRING
               READ(15,*,ERR=666) ((HESS(J1,J3),J3=1,3*NATOMS),J1=1,3*NATOMS)

               REWIND(15)
162            READ(15,'(A80)') GSTRING
               IF (GSTRING(1:17).NE.'block = tr_matrix') GOTO 162
               READ(15,*) GAMESR(1,1), GAMESR(1,2), GAMESR(1,3), GAMEST(1)
               READ(15,*) GAMESR(2,1), GAMESR(2,2), GAMESR(2,3), GAMEST(2)
               READ(15,*) GAMESR(3,1), GAMESR(3,2), GAMESR(3,3), GAMEST(3)
               DO J1=1,NATOMS
                  DO J3=1,3
                     TEMPX=0.0D0
                     DO J5=1,3
                        TEMPX=TEMPX+GAMESR(J3,J5)*VNEW(3*(J1-1)+J5)
                     ENDDO
                     TEMPV(J3)=TEMPX
                  ENDDO
                  DO J3=1,3
                     VNEW(3*(J1-1)+J3)=TEMPV(J3)
                  ENDDO
                  DO J2=1,NATOMS
                     DO J3=1,3
                        DO J4=1,3
                           TEMPXX=0.0D0
                           DO J5=1,3
                              DO J6=1,3
                                 TEMPXX=TEMPXX+GAMESR(J3,J5)*HESS(3*(J1-1)+J5,3*(J2-1)+J6)*GAMESR(J4,J6)
                              ENDDO
                           ENDDO
                           TEMPH(J3,J4)=TEMPXX
                        ENDDO
                     ENDDO
                     DO J3=1,3
                        DO J4=1,3
                           HESS(3*(J1-1)+J3,3*(J2-1)+J4)=TEMPH(J3,J4)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ELSE IF (GTEST) THEN
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling GAMESS-UK energy and first derivatives'
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints > ' // SYS(1:LSYS) // '.in ) >& /dev/null')
                  CALL SYSTEM("echo 'runtype gradient' >> " // SYS(1:LSYS) // ".in")
                  CALL SYSTEM("echo 'enter' >> " // SYS(1:LSYS) // ".in")
                  CALL SYSTEM('rungamess ' // SYS(1:LSYS) // '>&/dev/null' )
                  CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.pun derivs ')
               ELSE
                  WRITE(*,'(A)') ' Using existing GAMESS-UK derivs file for energy and first derivatives'
               ENDIF

               OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
13             READ(15,'(A80)') GSTRING
               IF (GSTRING(1:17).NE.'block = gradients') GOTO 13
               DO J1=1,NATOMS
                  READ(15,'(1X,3F16.6)',ERR=666) VNEW(3*(J1-1)+1),VNEW(3*(J1-1)+2),VNEW(3*(J1-1)+3)
               ENDDO
               REWIND(15)
163            READ(15,'(A80)') GSTRING
               IF (GSTRING(1:17).NE.'block = tr_matrix') GOTO 163
               READ(15,'(1X,4F16.6)') GAMESR(1,1), GAMESR(1,2), GAMESR(1,3), GAMEST(1)
               READ(15,'(1X,4F16.6)') GAMESR(2,1), GAMESR(2,2), GAMESR(2,3), GAMEST(2)
               READ(15,'(1X,4F16.6)') GAMESR(3,1), GAMESR(3,2), GAMESR(3,3), GAMEST(3)
               DO J1=1,NATOMS
                  DO J3=1,3
                     TEMPX=0.0D0
                     DO J5=1,3
                        TEMPX=TEMPX+GAMESR(J3,J5)*VNEW(3*(J1-1)+J5)
                     ENDDO
                     TEMPV(J3)=TEMPX
                  ENDDO
                  DO J3=1,3
                     VNEW(3*(J1-1)+J3)=TEMPV(J3)
                  ENDDO
               ENDDO
            ELSE
               IF (.NOT.YESNO) THEN
                  WRITE(*,'(A)') ' Calling GAMESS-UK energy'
                  CALL SYSTEM(' ( ' // EDITIT // ' < temppoints > ' // SYS(1:LSYS) // '.in ) >& /dev/null')
                  CALL SYSTEM("echo 'runtype energy' >> " // SYS(1:LSYS) // ".in")
                  CALL SYSTEM("echo 'enter' >> " // SYS(1:LSYS) // ".in")
                  CALL SYSTEM('rungamess ' // SYS(1:LSYS) // '>&/dev/null' )
                  CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.pun derivs ')
               ELSE
                  WRITE(*,'(A)') ' Using existing GAMESS derivs file for energy'
               ENDIF
            ENDIF
            ! 
            ! There is provision for punching the energy consistently in GAMESS-UK.
            ! 
            REWIND(15)
14          READ(15,'(A80)') GSTRING
            IF (GSTRING(1:20).NE.'block = total_energy') GOTO 14
            READ(15,*) ENERGY
            CLOSE(15)
            OPEN(UNIT=15,FILE='abenergy',STATUS='UNKNOWN')
            WRITE(15,'(F20.10)') ENERGY
            CLOSE(15)
            CALL SYSTEM('mv ' // SYS(1:LSYS) // '.in ' // SYS(1:LSYS) // '.old.in')
            CALL SYSTEM('mv ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out')
            IF (GTEST.OR.SSTEST) CALL SYSTEM('mv derivs derivs.old')
         ELSE IF (GAUSSIAN) THEN
            IF ((.NOT.ALLOCATED(HGAUSS)).AND.SSTEST) ALLOCATE(HGAUSS(3*NATOMS,3*NATOMS))
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            PRINT*,'Reading GAUSSIAN derivative information'
            OPEN(UNIT=15,FILE='derivs',STATUS='OLD')
            READ(15,*) (VNEW(J1),J1=1,3*NATOMS)
            ! PRINT*,'Gaussian forces:'
            ! WRITE(*,30) (VNEW(J1),J1=1,3*NATOMS)
C30         FORMAT(3F15.6)
            ! 
            ! Are they forces or first derivatives of the energy??
            ! 
            DO J1=1,3*NATOMS
               VNEW(J1)=-VNEW(J1)
            ENDDO
            NSTEP=1
            NSTART=1
            NFINISH=5
            ! PRINT*,'Gaussian second derivatives:'
40          DO J1=NSTEP,3*NATOMS
               READ(15,*) (HESS(J1,J2),J2=NSTART,MIN(NFINISH,J1))
               DO J2=NSTART,MIN(NFINISH,J1)
                  HESS(J2,J1)=HESS(J1,J2)
               ENDDO
               ! WRITE(*,116) (HESS(J1,J2),J2=NSTART,MIN(NFINISH,J1))
            ENDDO
            IF (NFINISH.EQ.3*NATOMS) GOTO 50
            NSTART=NSTART+5
            NFINISH=MIN(NFINISH+5,3*NATOMS)
            NSTEP=NSTEP+5
            GOTO 40

50          CONTINUE
            ! 
            ! Permute these matrix elements from the daft Gaussian convention
            ! to something sensible.
            ! 
            DO J1=1,3*NATOMS
               J3=MOD(J1,9)
               IF (J3.EQ.0) J3=9
               J3=3*(J3-1)+(J1-1)/9+1
               DO J2=1,3*NATOMS
                  J4=MOD(J2,9)
                  IF (J4.EQ.0) J4=9
                  J4=3*(J4-1)+(J2-1)/9+1
                  HGAUSS(J4,J3)=HESS(J2,J1)
               ENDDO
            ENDDO
            DO J1=1,3*NATOMS
               DO J2=1,3*NATOMS
                  HESS(J2,J1)=HGAUSS(J2,J1)
               ENDDO
            ENDDO
            CLOSE(15)
         ELSE IF (GAUSSIAN03.OR.GAUSSIAN09.OR.GAUSSIAN16) THEN
            IF (GAUSSIAN03) WRITE(*,'(A)') ' potential> Calling GAUSSIAN03'
            IF (GAUSSIAN09) WRITE(*,'(A)') ' potential> Calling GAUSSIAN09'
            IF (GAUSSIAN16) WRITE(*,'(A)') ' potential> Calling GAUSSIAN16'
            CALL SYSTEM('cp gaussian.inp.pre gaussian.inp' )
            LUNIT=GETUNIT()
            OPEN(LUNIT,file="gaussian.inp",STATUS="OLD",POSITION="APPEND")
            INQUIRE(FILE='gaussian.chk',EXIST=GAUSSIANTEST) 
            IF(GAUSSIANTEST)THEN
               WRITE(LUNIT,"(A10)") "GUESS=READ"
            ENDIF
            WRITE(LUNIT,*) ""
            WRITE(LUNIT,*) "OPTIM INTERFACE"
            WRITE(LUNIT,*) ""
            WRITE(LUNIT,"(2I4)") GAUSSIANCHARGE,GAUSSIANMULTI
            LUNIT2=GETUNIT()
            OPEN(LUNIT2,file="gaussian.name")
            DO J1=1,NATOMS
               READ(LUNIT2,*) STRING
               WRITE(LUNIT,"(A2,3F20.10)") ADJUSTL(TRIM(STRING)),COORDS(3*J1-2:3*J1)
            ENDDO
            CLOSE(LUNIT2)
            CLOSE(LUNIT)
            INQUIRE(FILE='gaussian.inp.after',EXIST=GAUSSIANTEST)
            IF(.NOT.GAUSSIANTEST) THEN
               CALL SYSTEM('echo "" > gaussian.inp.after' )
               CALL SYSTEM('echo "" >> gaussian.inp.after' )
               CALL SYSTEM('echo "" >> gaussian.inp.after' )
            ENDIF
            CALL SYSTEM('cat gaussian.inp.after >> gaussian.inp' )
            IF (SSTEST) THEN
               CALL SYSTEM('sed "s/Force/Freq=Noraman/Ig" gaussian.inp >gaussian.inp.temp')
               CALL SYSTEM('mv gaussian.inp.temp gaussian.inp')
            ENDIF
            INQUIRE(FILE='gaussian.out',EXIST=GAUSSIANTEST) 
            IF(GAUSSIANTEST) THEN
               CALL SYSTEM('mv gaussian.out gaussian.out.old')
            ENDIF
            IF (GAUSSIAN03) CALL SYSTEM('g03 <gaussian.inp> gaussian.out')
            IF (GAUSSIAN09) CALL SYSTEM('g09 <gaussian.inp> gaussian.out')
            IF (GAUSSIAN16) CALL SYSTEM('g16 <gaussian.inp> gaussian.out')
            IF (SSTEST) THEN
               !WRITE(*,'(A)') ' potential> extract GAUSSIAN second derivatives'
               J1=(1+3*NATOMS)*3*NATOMS/2
               IF (.NOT.ALLOCATED(HGAUSS)) ALLOCATE(HGAUSS(J1,1))
               IF(MOD(J1,5)==0) then
                  J1=J1/5
               ELSE
                  J1=J1/5+1
               ENDIF
               WRITE (STRING, *) J1
               STRING='grep -A' //
     &         ADJUSTL(TRIM(STRING)) // 
     &         ' "Cartesian Force Constants" Test.FChk|tail -' //
     &         ADJUSTL(TRIM(STRING)) //
     &         ' >gaussian.hess'
               CALL SYSTEM(ADJUSTL(TRIM(STRING)))
               OPEN(15,file="gaussian.hess")
               J1=(1+3*NATOMS)*3*NATOMS/2
               READ(15,*) HGAUSS(1:J1,1)
               CLOSE(15)
               J3=1
               DO J1=1,3*NATOMS
                  DO J2=1,J1
                     HESS(J1,J2)=HGAUSS(J3,1)
                     J3=J3+1
                     HESS(J2,J1)=HESS(J1,J2)
                  ENDDO
               ENDDO
            ENDIF
            IF (GTEST) THEN
               !WRITE(*,'(A)') ' potential> extract GAUSSIAN second derivatives'
               J1=3*NATOMS
               IF(MOD(J1,5)==0) then
                  J1=J1/5
               ELSE
                  J1=J1/5+1
               ENDIF
               WRITE (STRING, *) J1
               STRING='grep -A' //
     &         ADJUSTL(TRIM(STRING)) //
     &         ' "Cartesian Gradient" Test.FChk|tail -' //
     &         ADJUSTL(TRIM(STRING)) //
     &         ' >gaussian.grad'
               CALL SYSTEM(STRING)
               J1=3*NATOMS
               OPEN(15,file="gaussian.grad")
               READ(15,*) VNEW(1:3*NATOMS)
               CLOSE(15)
            ENDIF
            CALL SYSTEM('grep "Total Energy" Test.FChk|sed "s/Total Energy//g"|sed "s/R//g" >gaussian.energy')
            OPEN(15,file="gaussian.energy")
            READ(15,*) ENERGY
            IF (PTEST) WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' hartree' 
            CLOSE(15)
            OPEN(15,file="gaussian.xyz",STATUS="OLD",POSITION="APPEND")
            WRITE(15,*) NATOMS
            WRITE(15,*) "AFTER STEP",ITG03,"Energy=",ENERGY
            OPEN(7,file="gaussian.name")
            DO J1=1,NATOMS
               READ(7,*) STRING
               WRITE(15,"(A2,3F20.10)") ADJUSTL(TRIM(STRING)),COORDS(3*J1-2:3*J1)
            ENDDO
            CLOSE(7)
            CLOSE(15)
            ELSE IF (DFTBT) THEN
            CALL DFTB(NATOMS,COORDS,VNEW,ENERGY,GTEST,PARAM1,PARAM2,PARAM3,PARAM4)
            IF (SSTEST) CALL SECDFTB(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,PARAM4)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (AMBER12T) THEN
            VNEW = 0.0D0
            IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
               XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
               CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS, XRIGIDCOORDS)
            ENDIF
! Check cis/trans isomerisation - actually, chirality of sugar ring carbons for nucleic acids
            IF (CHECKCISTRANSALWAYSDNA .OR. CHECKCISTRANSALWAYSRNA) THEN ! CHECKCHIRALT removed from IF statement
               CALL CHIRALITY_CHECK(COORDS, GOODSTRUCTURE1)
            END IF
! Check cis/trans isomerisation for proteins
            IF (CHECKCISTRANSALWAYS) THEN
               CALL CIS_TRANS_CHECK(COORDS, GOODSTRUCTURE1)
            END IF
            IF (CUDAT) THEN
               ! This call copies CPU coordinates to GPU, calculates energy/gradient and copies energy/gradient back to CPU
               CALL CUDA_ENEGRAD_WRAPPER(NATOMS, COORDS, ENERGY, GRADATOMS)
            ELSE
               CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS,
     &                                          COORDS,
     &                                          ENERGY,
     &                                          GRADATOMS,
     &                                          ENERGY_DECOMP)
            END IF
            VNEW(1:3*NATOMS) = GRADATOMS(:)
! Calculate the numerical hessian
            IF (STEST) THEN
               IF (.NOT. ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS, 3*NATOMS))
               CALL AMBER12_NUM_HESS(NATOMS, COORDS, DELTA=1.0D-4, HESSIAN=HESS(:, :))
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT) ) THEN
               IF (STEST) THEN
                  CALL TRANSFORMHESSIAN(HESS, GRADATOMS, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
                  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
               ENDIF
               CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
               COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
               COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
               VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
            ENDIF
! AMBER 9 energy calculation
         ELSE IF (AMBERT) THEN
            VNEW = 0.0D0
            IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
               XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
               CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS, XRIGIDCOORDS)
            ENDIF
            IF (CHECKCISTRANSALWAYS) CALL CHECK_CISTRANS_PROTEIN(COORDS,NATOMS,goodstructure1,MINOMEGA,CISARRAY1)
            IF (CHECKCISTRANSALWAYSDNA) CALL CHECK_CISTRANS_DNA(COORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
            IF (CHECKCISTRANSALWAYSRNA) CALL CHECK_CISTRANS_RNA(COORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
            CALL AMBERENERGIES(COORDS,GRADATOMS,ENERGY,GTEST,SSTEST)
            VNEW(1:3*NATOMS) = GRADATOMS(:)
            IF (STEST) CALL AMBERSECDER(COORDS,.TRUE.)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT) ) THEN
               IF (STEST) THEN
                  CALL TRANSFORMHESSIAN(HESS, GRADATOMS, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
                  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
               ENDIF
               CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
               COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
               COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
               VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
            ENDIF
            ! DO i=1,3*NATOMS
            ! WRITE(*,*) VNEW(i)
            ! END DO
         ELSE IF (NABT) THEN
            ! CALL CHECK_CISTRANS_RNA(COORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
            ! WRITE(*,*) 'GOODSTRUCTURE1=', GOODSTRUCTURE1
            ! STOP
            ! NAB structures initialised in keywords.f
            VNEW = 0.0D0
            IF(STEST) THEN
               ! Analytical second derivatives
               ! hk286
               IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
                  XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
                  CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS, XRIGIDCOORDS)
               ENDIF
               ! hk286 - numerical second derivatives
               ! CALL NABSECDER(COORDS,SSTEST)

               IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
               IF (.NOT.ALLOCATED(TEMPHESS)) ALLOCATE(TEMPHESS(9*NATOMS*NATOMS))
               TEMPHESS(:) = 0.0D0
               CALL MME2WRAPPER(COORDS,ENERGY,GRADATOMS,TEMPHESS,ATMASS,GRAD1)
               IF (SQRT(ABS(ENERGY)) > 1.0D40) THEN
                  WRITE(*, *) "Linear dihedral detected in NAB routines (sff2.c)."
               END IF
               VNEW(1:3*NATOMS) = GRADATOMS(:)
               ALLOCATE(HESS(3*NATOMS,3*NATOMS))
               k=1
               DO i=1,3*NATOMS
                  DO j=1,3*NATOMS
                     HESS(i,j) = TEMPHESS(k)
                     k=k+1
                  END DO
               END DO
               DEALLOCATE(TEMPHESS)
               ! hk286
               IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
                  CALL TRANSFORMHESSIAN(HESS, GRADATOMS, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
                  CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
                  COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
                  VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
                  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
               ENDIF
            ELSE
               ! hk286
               IF (.NOT. ATOMRIGIDCOORDT) THEN
                  XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
                  CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS, XRIGIDCOORDS)
               ENDIF
               ! CALL AMBERNUMFIRSTDER(COORDS,GTEST)

               ! check cis-trans isomerisation for DNA or RNA
               IF (CHECKCISTRANSALWAYSDNA) CALL CHECK_CISTRANS_DNA(COORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
               IF (CHECKCISTRANSALWAYSRNA) CALL CHECK_CISTRANS_RNA(COORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
               IF (CHECKCISTRANSALWAYS) CALL CHECK_CISTRANS_PROTEIN(COORDS,NATOMS,goodstructure1,MINOMEGA,CISARRAY1)
               ! CALL CHECK_CISTRANS_DNA(COORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
               ! STOP
               CALL AMBERENERGIES(COORDS,GRADATOMS,ENERGY,GTEST,SSTEST)
               VNEW(1:3*NATOMS) = GRADATOMS(:)
               ! CALL MME(ENERGY,COORDS,VNEW,1)
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            ! hk286
            IF (.NOT. ATOMRIGIDCOORDT) THEN
               CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
               COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
               COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
               VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
            ENDIF

         ELSE IF (AMHT) THEN
            ! IF (DEBUG) WRITE(6,*)'Entering WALESAMH_INTERFACE'
            ! DO J1=1,NATOMS
            ! WRITE(6,*)'COORDS WALESAMH_INTERFACE',COORDS(J1)
            ! ENDDO

            CALl WALESAMH_INTERFACE(COORDS,VNEW,ENERGY)

            IF (PTEST) THEN
               WRITE(*,10) ' Energy for last cycle=',ENERGY,' Wolynes Units'
               WRITE(ESTRING,10) ' Energy for last cycle=',ENERGY,' Wolynes Units'
            ENDIF

            ! DIFF=1.0D-3
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! WRITE(*,'(F20.10,2x,F20.10,2xI5)')X(J1),GRAD(J1),J1
            ! X(J1)=X(J1)+DIFF
            ! CALL WALESAMH_INTERFACE(X,GRADDUM,EPLUS)
            ! X(J1)=X(J1)-2.0D0*DIFF
            ! CALL WALESAMH_INTERFACE(X,GRADDUM,EMINUS)
            ! X(J1)=X(J1)+DIFF
            ! IF (GRAD(J1).NE.0.0D0) WRITE(*,'(A5,I5,3F20.10)') fff,J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),
            ! 1                          100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
            ! GRAD(J1)=(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDDO

         ELSE IF (CHRMMT) THEN
            ! hk286
            IF (RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.)) THEN
               XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
               CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS, XRIGIDCOORDS)
            ENDIF
            ! IF (CHECKOMEGAT.AND.DEBUG) THEN ! DJW
            ! AMIDEFAIL=.FALSE.
            ! CALL CHECKOMEGA(COORDS,AMIDEFAIL)
            ! IF (AMIDEFAIL) PRINT '(A,L5)',' potential> WARNING *** cis peptide bond detected'
            ! ENDIF

            CALL OCHARMM(COORDS,VNEW,ENERGY,GTEST,SSTEST)

            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            ! hk286
            IF ( RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.) ) THEN
               IF (SSTEST) THEN
                  CALL TRANSFORMHESSIAN(HESS, VNEW, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
                  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
               ENDIF
               CALL TRANSFORMGRAD(VNEW, XRIGIDCOORDS, XRIGIDGRAD)
               COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
               VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
            ENDIF

         ELSE IF (OPEPT) THEN
            !used AMBER12 as example and simply replaced potential and
            !gradient calls accordingly
            VNEW = 0.0D0
            IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT)) THEN
               XRIGIDCOORDS(1:DEGFREEDOMS) = COORDS(1:DEGFREEDOMS)
               CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, COORDS,XRIGIDCOORDS)
            ENDIF
            CALL OPEP_ENERGY_AND_GRADIENT(NATOMS,COORDS,GRADATOMS,ENERGY,.TRUE.)
            VNEW(1:3*NATOMS) = GRADATOMS(:)
            ! Calculate the numerical hessian
            IF (STEST) THEN
               IF (.NOT. ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
               CALL OPEP_NUM_HESS(NATOMS, COORDS, DELTA=1.0D-5,HESSIAN=HESS(:,:))
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,'kcal/mol'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,'kcal/mol'
            ENDIF
            IF (RIGIDINIT .AND. (.NOT. ATOMRIGIDCOORDT) ) THEN
               IF (STEST) THEN
                  CALL TRANSFORMHESSIAN(HESS, GRADATOMS, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
                  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) =XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
               ENDIF
               CALL TRANSFORMGRAD(GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD)
               COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
               COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
               VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
               VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
            ENDIF

         ELSE IF (UNRST) THEN
            CALL UENERGY(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' kcal/mol'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' kcal/mol'
            ENDIF
            ! 
            ! Potentials identified by atom types corresponding to real atomic symbols must
            ! come at the end in case someone wants to use this atom in one of the electronic
            ! strcuture codes identified by keywords above.
            ! 
         ELSE IF (ZSYM(NATOMS).EQ.'BE') THEN
            CALL ETRAP(NATOMS,COORDS,ENERGY,C1,C2,C3)
            CALL DTRAP(NATOMS,COORDS,VNEW,C1,C2,C3,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
            ENDIF
         ELSE IF ((ZSYM(NATOMS).EQ.'AU').OR.(ZSYM(NATOMS).EQ.'AG')
     &   .OR.(ZSYM(NATOMS).EQ.'NI').OR.(ZSYM(NATOMS).EQ.'PT')) THEN
            IF ((ZSYM(NATOMS).EQ.'AU').OR.(ZSYM(NATOMS).EQ.'PT')) THEN
               NN=10
               MM=8
            ELSE IF (ZSYM(NATOMS).EQ.'AG') THEN
               NN=12
               MM=6
            ELSE IF (ZSYM(NATOMS).EQ.'NI') THEN
               NN=9
               MM=6
            ENDIF
            PARAM4=FLOAT(NN)
            PARAM5=FLOAT(MM)
            ! 
            ! PARAM1 is really EPS
            ! PARAM2 is really C
            ! PARAM3 is really SIG
            ! 
            ! PRINT*,' WARNING - GTEST and SSTEST ignored'
            IF (.not.BULKT) then
            CALL SCDIFF(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,NN,MM,ENERGY)
            else
            CALL SCDIFF_periodic(NATOMS,COORDS,VNEW,PARAM1,PARAM2,PARAM3,NN,MM,ENERGY)
            endif
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
               WRITE(ESTRING,10)' potential> Energy for last cycle=',ENERGY,' eV'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'IN') THEN
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL IONS(NATOMS,COORDS,VNEW,ENERGY,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,1)
            WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' hartree'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
10          FORMAT(A,27X,F20.10,A)

         ELSE IF (ZSYM(1).EQ.'CA') THEN
            CALL CAARDIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'AR') THEN
            CALL LJDIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'RM') THEN
            CALL LJDIFF(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'TH') THEN
            CALL THOMSON(NATOMS, COORDS, VNEW, ENERGY, GTEST, SSTEST)

            ! hk286
         ELSE IF (GTHOMSONT) THEN
            CALL GTHOMSON(COORDS, VNEW, ENERGY, GTEST, SSTEST)

            ! 
            ! Check numerical first and second derivatives
            ! 
            ! DIFF=1.0D-4
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL THOMSON(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL THOMSON(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
            ! WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! PRINT*,'analytic and numerical second derivatives:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL THOMSON(NATOMS,COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL THOMSON(NATOMS,COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! DO J2=1,3*NATOMS
            ! IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
            ! 1             (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
            ! ELSE
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! ENDDO
            IF (PTEST) THEN
               WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'SI') THEN
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL JM2(NATOMS, COORDS, VNEW)
            CALL JM3(NATOMS, COORDS, VNEW)
            CALL JME(NATOMS, COORDS, P2, P3, ENERGY)
            WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' eV'
            WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' eV'
            WRITE(*,10) 'Two-body contribution=',P2,' eV'
            WRITE(*,10) 'Three-body contribution=',P3,' eV'
         ELSE IF (ZSYM(NATOMS)(1:1).EQ.'W') THEN
            IF (ZSYM(NATOMS).EQ.'W5') IPOT=5
            IF (ZSYM(NATOMS).EQ.'W4') IPOT=4
            IF (ZSYM(NATOMS).EQ.'W3') IPOT=3
            IF (ZSYM(NATOMS).EQ.'W2') IPOT=2
            IF (ZSYM(NATOMS).EQ.'W1') IPOT=1
            IF (ANGLEAXIS) THEN
               CALL TIPNP(IPOT,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ELSE
               IF (IPOT.EQ.5) THEN
                  PRINT*,' TIP5P not yet coded in C of M/Euler coordinates'
                  STOP
               ENDIF
               CALL H2O(NATOMS/2,IPOT,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ENDIF

            ! Test derivatives

            ! deltaCoord = 1.0d-5
            ! RMSdiff = 0.0d0

            ! do j1 = 1, 3*NATOMS
            ! coords(j1) = coords(j1) + deltaCoord
            ! call TIPnP(IPOT,COORDS,dummyGrad,upperE,.FALSE.,.FALSE.)
            ! coords(j1) = coords(j1) - 2.0d0*deltaCoord
            ! call TIPnP(IPOT,COORDS,dummyGrad,lowerE,.FALSE.,.FALSE.)
            ! coords(j1) = coords(j1) + deltaCoord

            ! numericalGrad(j1) = (upperE-lowerE)/(2.0d0*deltaCoord)
            ! RMSdiff = RMSdiff + (numericalGrad(j1)-VNEW(j1))**2
            ! print *, j1, VNEW(j1), numericalGrad(j1)
            ! enddo

            ! RMSdiff = DSQRT(RMSdiff/(3*NATOMS))
            ! write(11,*) 'RMS difference:', RMSdiff

            ! if (SSTEST) then
            ! tempHess = HESS

            ! do j1 = 1, 3*NATOMS
            ! coords(j1) = coords(j1) + deltaCoord
            ! call TIPnP(IPOT,COORDS,upperGrad,upperE,.TRUE.,.FALSE.)
            ! coords(j1) = coords(j1) - 2.0d0*deltaCoord
            ! call TIPnP(IPOT,COORDS,lowerGrad,lowerE,.TRUE.,.FALSE.)
            ! coords(j1) = coords(j1) + deltaCoord

            ! do j2 = 1, j1
            ! numericalSD = (upperGrad(j2)-lowerGrad(j2))/(2.0d0*deltaCoord)
            ! write(12, *) j2, j1, tempHess(j2,j1), numericalSD, (tempHess(j2,j1)-numericalSD)
            ! RMSdiff = RMSdiff + (tempHess(j2,j1)-numericalSD)**2
            ! enddo
            ! end do

            ! HESS = tempHess

            ! RMSdiff = DSQRT(RMSdiff/((9*NATOMS*NATOMS+3*NATOMS)/2))
            ! write(11,*) 'RMS difference:', RMSdiff
            ! endif

            IF (PTEST) THEN
               IF (IPOT.LE.4) THEN
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' kJ/mol'
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' kJ/mol'
               ELSE
                  WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' hartree'
                  WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' hartree'
               ENDIF
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'CD') THEN
            CALL FCAPSID(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'SS') THEN
            ! param1 is the soft sphere power
            CALL SOFT_SPHERE_POTENTIAL(NATOMS,COORDS,VNEW,HESS,ENERGY, 
     &              BULK_BOXVEC, SOFT_SPHERE_RADII, PARAM1, GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY
            ENDIF
            ! davidg: introduced userpot here
!
! USERPOT is the only potential required for the Biovia build.
!
         ELSE IF (USERPOTT) THEN
! #endif
! #ifdef _OPTIMLIBRARY
!        IF (USERPOTT) THEN
! #endif
            CALL USERPOT_POTENTIALHESS(3*NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST,HESS)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' meV'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' meV'
            ENDIF
! #ifdef _OPTIMLIBRARY
!        ENDIF
! #else
            ELSE IF (STEALTHYT) THEN
               CALL STEALTHY(COORDS, VNEW, ENERGY, GTEST, SSTEST)
 
               !DIFF=1.0D-6
               !OPEN(UNIT=1500, FILE="grad-test")
               !WRITE(1500,*) 'ANALYTIC AND NUMERICAL GRADIENTS:NATOMS=',NATOMS
               !DO J1=1, 3*NATOMS
               !   COORDS(J1)=COORDS(J1)+DIFF
               !   CALL STEALTHY(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
               !   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
               !   CALL STEALTHY(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
               !   COORDS(J1)=COORDS(J1)+DIFF               
               !   IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*ABS((VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D-3)) THEN
               !      WRITE(1500,'(I5,2G20.10,A)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF), '  X'
               !   ELSE
               !      WRITE(1500,'(I5,2G20.10,A)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
               !   END IF
               !END DO
  
               !CALL STEALTHY(COORDS, VNEW, ENERGY, .TRUE., .TRUE.)
               !OPEN(UNIT=1499, FILE="hessian-matrix")
               !WRITE(1499,*) 'ANALYTIC AND NUMERICAL SECOND DERIVATIVES:'
               !DO J1=1, 3*NATOMS
               !   COORDS(J1)=COORDS(J1)+DIFF
               !   CALL STEALTHY(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
               !   COORDS(J1)=COORDS(J1)-2.0D0*DIFF
               !   CALL STEALTHY(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
               !   COORDS(J1)=COORDS(J1)+DIFF
               !   DO J2=1, 3*NATOMS
               !      IF (ABS(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)).GT.(1.0D-5*HESS(J1,J2))) THEN
               !         WRITE(1499,'(2I5,2F20.15,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'  X'
               !      ELSE
               !         WRITE(1499,'(2I5,2F20.15)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
               !      END IF
               !   END DO
               !END DO




            ! 
            ! Check numerical first and second derivatives
            ! 
            ! DIFF=1.0D-4
            ! PRINT*,'analytic and numerical gradients:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL FCAPSID(NATOMS,COORDS,HDUMM,VPLUS,EPLUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL FCAPSID(NATOMS,COORDS,HDUMM,VPLUS,EMINUS,.FALSE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
            ! WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! PRINT*,'analytic and numerical second derivatives:'
            ! DO J1=1,3*NATOMS
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! CALL FCAPSID(NATOMS,COORDS,HDUMM,VPLUS,EPLUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! CALL FCAPSID(NATOMS,COORDS,HDUMM,VMINUS,EMINUS,.TRUE.,.FALSE.)
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! DO J2=1,3*NATOMS
            ! IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
            ! 1             (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
            ! ELSE
            ! WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
            ! ENDIF
            ! ENDDO
            ! ENDDO
            ! DC430 >
         ENDIF
!   #endif
         ! 
         ! End of possible potentials - now add extra terms required, if any ------------------------------
         !
         !ds656> Apply a substrate field for keyword MIE_FIELD
         
! #ifndef _OPTIMLIBRARY         
!ds656> Test gradient and Hessian with MIEF
         IF(MIEFT) CALL MIEF(NATOMS,COORDS,VNEW,ENERGY,GTEST,SSTEST)
         IF (ALLOCATED(HESS).AND.MIEFT.AND..FALSE.) THEN
         DIFF=1.0D-4
         PRINT*,'analytic and numerical gradients:'
         DO J1=1,3*NATOMS
            COORDS(J1)=COORDS(J1)+DIFF
            CALL BGUPTA(COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
            CALL MIEF(NATOMS,COORDS,VPLUS,EPLUS,.FALSE.,.FALSE.)
            COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            CALL BGUPTA(COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
            CALL MIEF(NATOMS,COORDS,VMINUS,EMINUS,.FALSE.,.FALSE.)
            COORDS(J1)=COORDS(J1)+DIFF
            IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*(VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1).GT.1.0D0)) THEN
               WRITE(*,'(I5,2F20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
            ENDIF
         ENDDO
         PRINT*,'analytic and numerical second derivatives:'
         DO J1=1,3*NATOMS
            COORDS(J1)=COORDS(J1)+DIFF
            CALL BGUPTA(COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
            CALL MIEF(NATOMS,COORDS,VPLUS,EPLUS,.TRUE.,.FALSE.)
            COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            CALL BGUPTA(COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
            CALL MIEF(NATOMS,COORDS,VMINUS,EMINUS,.TRUE.,.FALSE.)
            COORDS(J1)=COORDS(J1)+DIFF
            DO J2=1,3*NATOMS
               IF ((ABS(HESS(J1,J2)).NE.0.0D0).AND.
     1              (ABS(100.0D0*(HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D0)) THEN
                  WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),'   X'
               ELSE
                  WRITE(*,'(2I5,2F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF)
               ENDIF
            ENDDO
         ENDDO
         STOP
         ENDIF
         !<ds656 End of test for gradient and Hessian with MIEF
         
         IF (PULLT) THEN
            ENERGY=ENERGY-PFORCE*(COORDS(3*(PATOM1-1)+3)-COORDS(3*(PATOM2-1)+3))
            VNEW(3*(PATOM1-1)+3)=VNEW(3*(PATOM1-1)+3)-PFORCE
            VNEW(3*(PATOM2-1)+3)=VNEW(3*(PATOM2-1)+3)+PFORCE
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ENDIF
         IF (TWISTT) THEN
            CALL TWIST(COORDS,NATOMS,VNEW,ENERGY,GTEST,SSTEST)
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle=',ENERGY,' epsilon'
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY,' epsilon'
            ENDIF
         ENDIF
         IF (AAORIENTT) THEN
            ! IF (DEBUG) PRINT '(A)',' potential> Adding additional angle-axis potential'
            ! ENERGY=0.0D0  !!!!!!!!!!!!!!!!!!!!! debug
            ! VNEW(1:NOPT)=0.0D0  !!!!!!!!!!!!!!!!!!!!! debug
            ! HESS(1:NOPT,1:NOPT)=0.0D0  !!!!!!!!!!!!!!!!!!!!! debug
            IF (SIGMAAA .EQ. 0.0D0) THEN
               CALL AAORIENT(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ELSE
               CALL AAORIENTSR(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            END IF
            IF (PTEST) THEN
               WRITE(*,10) ' potential> Energy for last cycle with additional angle-axis potential=',ENERGY
               WRITE(ESTRING,10) 'Energy for last cycle=',ENERGY
            ENDIF
            ! Debug tools
            ! 
            ! DIFF=1.0D-4
            ! PRINT*,'analytic and numerical gradients: NATOMS=',NATOMS
            ! DO J1=1,3*NATOMS
            ! EPLUS=0.0D0
            ! EMINUS=0.0D0
            ! VPLUS(1:NOPT)=0.0D0
            ! VMINUS(1:NOPT)=0.0D0
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF (SIGMAAA .EQ. 0.0D0) THEN
            ! CALL AAORIENT(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! ELSE
            ! CALL AAORIENTSR(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! END IF
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! IF (SIGMAAA .EQ. 0.0D0) THEN
            ! CALL AAORIENT(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! ELSE
            ! CALL AAORIENTSR(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! END IF
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF ((ABS(VNEW(J1)).NE.0.0D0).AND.(100.0D0*ABS((VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.0.0D0)) THEN
            ! IF (ABS((VNEW(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/VNEW(J1)).GT.1.0D-2) THEN
            ! WRITE(*,'(I5,3G20.10,A)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),(EPLUS-EMINUS)/(2.0D0*DIFF*VNEW(J1)),'   X'
            ! ELSE
            ! WRITE(*,'(I5,3G20.10)') J1,VNEW(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),(EPLUS-EMINUS)/(2.0D0*DIFF*VNEW(J1))
            ! ENDIF
            ! ENDIF
            ! ENDDO
            ! PRINT*,'analytic and numerical second derivatives:'
            ! DO J1=1,3*NATOMS
            ! EPLUS=0.0D0
            ! EMINUS=0.0D0
            ! VPLUS(1:NOPT)=0.0D0
            ! VMINUS(1:NOPT)=0.0D0
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! IF (SIGMAAA .EQ. 0.0D0) THEN
            ! CALL AAORIENT(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! ELSE
            ! CALL AAORIENTSR(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! END IF
            ! COORDS(J1)=COORDS(J1)-2.0D0*DIFF
            ! IF (SIGMAAA .EQ. 0.0D0) THEN
            ! CALL AAORIENT(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! ELSE
            ! CALL AAORIENTSR(COORDS,VNEW,ENERGY,GTEST,SSTEST)
            ! END IF
            ! COORDS(J1)=COORDS(J1)+DIFF
            ! DO J2=1,3*NATOMS
            ! IF (ABS(HESS(J1,J2)).GT.0.0D0) THEN
            ! IF (ABS((HESS(J1,J2)-(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF))/HESS(J1,J2)).GT.1.0D-2) THEN
            ! WRITE(*,'(2I5,3F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),
            ! &                                                        ABS((VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF*HESS(J1,J2))),'   X'
            ! ELSE
            ! WRITE(*,'(2I5,3F20.10,A)') J1,J2,HESS(J1,J2),(VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF),
            ! &                                                        ABS((VPLUS(J2)-VMINUS(J2))/(2.0D0*DIFF*HESS(J1,J2)))
            ! ENDIF
            ! ENDIF
            ! ENDDO
            ! ENDDO
            ! STOP
         ENDIF
         ! 
         ! Add on attractive term for "closest" minimum if this is a REDOPATH run
         ! Here we have to define closest in terms of which minimum we have moved
         ! towards, not in terms of absolute distance, because of asymmetric pathways.
         ! Avoid division by zero for D1INIT and D2INIT!
         ! 
! #endif
!          IF (REDOKADD.AND.REDOPATH.AND.(.NOT.REDOPATHXYZ).AND.(REDOK.NE.0.0D0).AND.
!      &   (ALLOCATED(MIN1REDO)).AND.(ALLOCATED(MIN2REDO)).AND.(D1INIT*D2INIT.NE.0.0D0)) THEN
! 
!             CALL NEWMINDIST(COORDS,MIN1REDO,NATOMS,DIST1,BULKT,TWOD,'AX   ',.FALSE.,RIGIDBODY,DEBUG,RMAT)
!             CALL NEWMINDIST(COORDS,MIN2REDO,NATOMS,DIST2,BULKT,TWOD,'AX   ',.FALSE.,RIGIDBODY,DEBUG,RMAT)
! 
!             IF     ((DIST1/D1INIT.LT.REDOFRAC).AND.REDOPATH1) THEN
!                DUMMY1=0.0D0
!                DUMMY2=0.0D0
!                DUMMY3=0.0D0
!                DO J1=1,NOPT
!                   ! VNEW(J1)=VNEW(J1)+REDOK*(COORDS(J1)-MIN1REDO(J1))
!                   DUMMY1=DUMMY1+VNEW(J1)*(MIN1REDO(J1)-COORDS(J1))
!                   DUMMY2=DUMMY2+VNEW(J1)**2
!                   DUMMY3=DUMMY3+(COORDS(J1)-MIN1REDO(J1))**2
!                ENDDO
!                DUMMY1=DUMMY1/SQRT(DUMMY3)
!                IF (DEBUG) PRINT '(A,2G15.5,A,G20.10)',' potential> distances/initial distance are ',DIST1/D1INIT,DIST2/D2INIT,
!      &         ' grad % towards first minimum=',DUMMY1*100/SQRT(DUMMY2)
!                ! IF (DUMMY1.GT.0.0D0) DUMMY1=-DUMMY1
!                DO J1=1,NOPT
!                   VNEW(J1)=VNEW(J1)+REDOK*(COORDS(J1)-MIN1REDO(J1))
!                   ! VNEW(J1)=VNEW(J1)+REDOK*(COORDS(J1)-MIN1REDO(J1))*SQRT(DUMMY2/DUMMY3)
!                ENDDO
!             ELSEIF ((DIST2/D2INIT.LT.REDOFRAC).AND.REDOPATH2) THEN
!                DUMMY1=0.0D0
!                DUMMY2=0.0D0
!                DUMMY3=0.0D0
!                DO J1=1,NOPT
!                   ! VNEW(J1)=VNEW(J1)+REDOK*(COORDS(J1)-MIN2REDO(J1))
!                   DUMMY1=DUMMY1+VNEW(J1)*(MIN2REDO(J1)-COORDS(J1))
!                   DUMMY2=DUMMY2+VNEW(J1)**2
!                   DUMMY3=DUMMY3+(COORDS(J1)-MIN2REDO(J1))**2
!                ENDDO
!                DUMMY1=DUMMY1/SQRT(DUMMY3)
!                IF (DEBUG) PRINT '(A,2G15.5,A,G20.10)',' potential> distances/initial distance are ',DIST1/D1INIT,DIST2/D2INIT,
!      &         ' grad % for second minimum=',DUMMY1*100/SQRT(DUMMY2)
!                ! IF (DUMMY1.GT.0.0D0) DUMMY1=-DUMMY1
!                DO J1=1,NOPT
!                   VNEW(J1)=VNEW(J1)+REDOK*(COORDS(J1)-MIN2REDO(J1))
!                   ! VNEW(J1)=VNEW(J1)+REDOK*(COORDS(J1)-MIN2REDO(J1))*SQRT(DUMMY2/DUMMY3)
!                ENDDO
!             ENDIF
!          ENDIF
! #ifndef _OPTIMLIBRARY
         ! 
         ! Add on terms for rotation about the z axis
         ! 
         IF (RTEST) THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            IF (JZ.NE.0.0D0) THEN
               CALL ROTD(NATOMS, COORDS, VNEW, 1.0D0, JZ, .FALSE., ROT)
            ELSE
               CALL ROTENERGY(NATOMS, COORDS, OMEGA, 1.0D0, IZ, ROT)
               CALL ROTDERIV(NATOMS, COORDS, VNEW, 1.0D0, OMEGA, IZ)
            ENDIF
            ENERGY=ENERGY+ROT
            WRITE(*,10) ' potential> Energy for last cycle including rotation=',ENERGY,' epsilon'
            WRITE(ESTRING,'(A,9X,F20.10,A)') ' potential> Energy for last cycle including rotation=',ENERGY,' epsilon'
         ENDIF

         IF (FIELDT) THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            PRINT*,' WARNING - GTEST and SSTEST ignored'
            CALL FD(NATOMS,COORDS,VNEW,ENERGY)
            WRITE(*,20) ' potential> Energy for last cycle including field=',ENERGY
            WRITE(ESTRING,20) ' potential> Energy for last cycle including field=',ENERGY
         ENDIF

         ! IF (FTEST) THEN
         ! IF (GFRACTION.NE.0.0D0) THEN
         ! PRINT*,' WARNING - GTEST and SSTEST ignored'
         ! CALL GAV(NATOMS,COORDS,VNEW,ENERGY,GALPHA,PARAM2,1)
         ! WRITE(*,'(A,F20.10)') ' Fraction of non-local Gaussian potential used=',GFRACTION
         ! WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
         ! WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
         ! WRITE(*,'(A23,7X,2F20.10)') ' RHO and DELTA=',PARAM1, PARAM2
         ! ENDIF
         ! IF (MFRACTION1.NE.0.0D0) THEN
         ! PRINT*,' WARNING - GTEST and SSTEST ignored'
         ! CALL MAV(NATOMS,COORDS,VNEW,ENERGY,GALPHA,PARAM2,1)
         ! WRITE(*,'(A,F20.10)') ' Fraction of non-local Morse1 potential used=',MFRACTION1
         ! WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
         ! WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
         ! WRITE(*,'(A23,7X,2F20.10)') ' RHO and DELTA=',PARAM1, PARAM2
         ! ENDIF
         ! IF (MFRACTION2.NE.0.0D0) THEN
         ! PRINT*,' WARNING - GTEST and SSTEST ignored'
         ! CALL M2(NATOMS,COORDS,VNEW,ENERGY,GALPHA,PARAM2,1)
         ! WRITE(*,'(A,F20.10)') ' Fraction of non-local Morse2 potential used=',MFRACTION2
         ! WRITE(*,20) ' potential> Energy for last cycle=',ENERGY
         ! WRITE(ESTRING,20) ' potential> Energy for last cycle=',ENERGY
         ! WRITE(*,'(A23,7X,2F20.10)') ' RHO and DELTA=',PARAM1, PARAM2
         ! ENDIF
         ! ENDIF
         ! 
         ! Double well potential between the first two atoms.
         ! 
         IF (DOUBLET) THEN
            PRINT*,'WARNING - this potential has not been tested in OPTIM.3.0'
            CALL DOUBLE(NATOMS,COORDS,VNEW,EDOUBLE,GTEST,SSTEST,PARAM4,PARAM5,PARAM6)
            ENERGY=ENERGY+EDOUBLE
            IF (PTEST) THEN
               WRITE(*,'(A,F20.10,A)') ' potential> Energy for last cycle including double well=     ',ENERGY,' epsilon'
               WRITE(ESTRING,'(A,F20.10,A)') ' potential> Energy for last cycle including double well=     ',ENERGY,' epsilon'
            ENDIF
         ENDIF

         IF (GAUSSIAN.OR.CADPAC.OR.GAMESSUK.OR.GAMESSUS) THEN
            INQUIRE(FILE='abenergy',EXIST=ETEST)
            IF (ETEST) THEN
               OPEN(UNIT=91,FILE='abenergy',STATUS='OLD')
               READ(91,*) ENERGY
               IF (PTEST) WRITE(*,'(A,27X,F20.10,A)') ' potential> Energy for last cycle=',ENERGY,' hartree'
               CLOSE(91)
            ELSE
               WRITE(*,'(A)') ' potential> Error - abenergy file not found'
               STOP
            ENDIF
         ENDIF
! #endif

         IF (GTEST) THEN
            ! PRINT '(A,L5)',' potential> FREEZE,VNEW=',FREEZE
            ! PRINT '(3G20.10)',VNEW(1:3*NATOMS)
            IF (FREEZE) THEN
               DO J1=1,NATOMS
                  IF (FROZEN(J1)) THEN
                     IF (VARIABLES) THEN
                        VNEW(J1)=0.0D0
                     ELSE
                        VNEW(3*(J1-1)+1)=0.0D0
                        VNEW(3*(J1-1)+2)=0.0D0
                        VNEW(3*(J1-1)+3)=0.0D0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
! #ifndef _OPTIMLIBRARY
            ! sn402: other half of hack to allow genrigid to work with arbitrary potential
            IF (RIGIDMOLECULEST .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.)) THEN
                IF (STEST) THEN
                     CALL TRANSFORMHESSIAN(HESS, VNEW, XRIGIDCOORDS, XRIGIDHESS, RBAANORMALMODET)
                         HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
                         HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0
                         HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = XRIGIDHESS(1:DEGFREEDOMS,1:DEGFREEDOMS)
                ENDIF
                CALL TRANSFORMGRAD(VNEW, XRIGIDCOORDS, XRIGIDGRAD)
                VNEW(1:DEGFREEDOMS) = XRIGIDGRAD(1:DEGFREEDOMS)
                VNEW(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
            ENDIF
! #endif

            IF (UNRST) THEN
               CALL VSTAT(VNEW,TEMP,NINTS,NOPT)
               IF (PTEST) WRITE(*,'(A,43X,F15.10,2X,A,G15.10)') ' potential> RMS force: ',TEMP(5),' |gradient|=',
     &         TEMP(5)*SQRT(1.0D0*(NINTS))
            ELSE IF ( RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.) ) THEN !hk286
               CALL VSTAT(VNEW(1:DEGFREEDOMS),TEMP,DEGFREEDOMS,DEGFREEDOMS)
               IF (PTEST) WRITE(*,'(A,43X,F15.10,2X,A,G15.10)') ' potential> RMS force: ',TEMP(5),' |gradient|=',
     &         TEMP(5)*SQRT(1.0D0*(DEGFREEDOMS))
               IF (AACONVERGENCET .AND. (TEMP(5) < 5.0D0 * GMAX)) THEN
                  CALL AACONVERGENCE (GRADATOMS, XRIGIDCOORDS, XRIGIDGRAD, TEMP(5))
               END IF
            ELSE
               CALL VSTAT(VNEW,TEMP,NOPT,NOPT)
               IF (PTEST) WRITE(*,'(A,43X,G15.10,2X,A,G15.10)') ' potential> RMS force: ',TEMP(5),' |gradient|=',
     &         TEMP(5)*SQRT(1.0D0*(NOPT))
            ENDIF
            RMS=TEMP(5)
            IF(DEBUG.AND.(RMS.NE.RMS)) THEN
               WRITE(*,'(A)' ) 'potential> WARNING - RMS force is NaN - if using AMBER igb=1, can be due to negative Born radii'
            END IF
            ! PRINT '(A,G20.10)',' potential> RMS=',RMS
            IF (CPMD.AND.(RMS.EQ.0.0D0)) RMS=1.0D0  !  to prevent convergence when CPMD SCF fails

            IF (REDOKADD.AND.REDOPATH.AND.ALLOCATED(MIN1REDO).AND.ALLOCATED(MIN2REDO).AND.(PLUSSIDET.OR.MINUSSIDET)) THEN
!!!!!!!!!!!!!!!!!!!! DEBUG DJW
!           ZSTAR=PARAM1
!           CALL LJDIFF(NATOMS, MIN1REDO, VPLUS, EPLUS, .FALSE., .FALSE.)
!           CALL LJDIFF(NATOMS, MIN2REDO, VMINUS, EMINUS, .FALSE., .FALSE.)
!           PRINT '(A,2G20.10)','potential> plus and minus minima energies=',EPLUS,EMINUS
!!!!!!!!!!!!!!!!!!!! DEBUG DJW

               IF (PLUSSIDET) THEN
!              IF (DEBUG) PRINT '(A,G20.10)','Adding bias term to gradient for plus minimum with force constant ',REDOK
                  DO J1=1,NOPT
                     VNEW(J1)=VNEW(J1)+2.0D0*REDOK*RMS*(COORDS(J1)-MIN1REDO(J1))
!                    VNEW(J1)=VNEW(J1)+2.0D0*REDOK*(COORDS(J1)-MIN1REDO(J1))
                  ENDDO
               ELSEIF (MINUSSIDET) THEN
!              IF (DEBUG) PRINT '(A,G20.10)','Adding bias term to gradient for minus minimum with force constant ',REDOK
                  DO J1=1,NOPT
                     VNEW(J1)=VNEW(J1)+2.0D0*REDOK*RMS*(COORDS(J1)-MIN2REDO(J1))
!                    VNEW(J1)=VNEW(J1)+2.0D0*REDOK*(COORDS(J1)-MIN2REDO(J1))
                  ENDDO
               ENDIF
            ENDIF

         ENDIF
         ! 
         ! If the Hessian gets overwritten by diagonalisation we must read it back in before updating!
         ! 
         IF (GTEST.AND.(STEST.AND.HUPDATE)) THEN
            NHUP=NHUP+1
            IF (INTHUP.EQ.-1) THEN
               WRITE(*,'(A)') ' potential> Not updating Hessian'
               IF (NHUP.GT.1) THEN
                  OPEN(UNIT=34,FILE='hessdump',FORM='UNFORMATTED',STATUS='UNKNOWN')
                  READ(34) ((HESS(J2,J1),J2=1,NOPT),J1=1,NOPT)
                  CLOSE(34)
                  ! DO J1=1,NOPT
                  ! DO J2=1,NOPT
                  ! HESS(J2,J1)=HSAVE(J2,J1)
                  ! ENDDO
                  ! ENDDO
               ENDIF
            ELSE IF ((NHUP.GT.1).AND.(.NOT.SSTEST)) THEN
               WRITE(*,'(A)') ' potential> Updating Hessian'
               OPEN(UNIT=34,FILE='hessdump',FORM='UNFORMATTED',STATUS='UNKNOWN')
               READ(34) ((HESS(J2,J1),J2=1,NOPT),J1=1,NOPT)
               CLOSE(34)
               ! CALL HUPD(HSAVE,COORDS,COORDSO,VNEW,GRADO,PHI)
               CALL HUPD(COORDS,COORDSO,VNEW,GRADO,PHIG)
            ELSE IF ((NSTHUP.NE.1).AND.(.NOT.READHESS)) THEN
               DO J1=1,NOPT
                  DO J2=J1+1,NOPT
                     HESS(J2,J1)=0.0D0
                     HESS(J1,J2)=0.0D0
                  ENDDO
                  HESS(J1,J1)=1.0D0
               ENDDO
            ENDIF
            OPEN(UNIT=34,FILE='hessdump',FORM='UNFORMATTED',STATUS='UNKNOWN')
            WRITE(34) ((HESS(J2,J1),J2=1,NOPT),J1=1,NOPT)
            CLOSE(34)
            DO J1=1,NOPT
               COORDSO(J1)=COORDS(J1)
               GRADO(J1)=VNEW(J1)
               ! DO J2=1,NOPT
               ! HSAVE(J2,J1)=HESS(J2,J1)
               ! ENDDO
            ENDDO
         ENDIF
         READHESS=.FALSE.

         IF (TWOD) THEN
            DO J1=1,NATOMS
               J2=3*J1
               COORDS(J2)=0.0D0
               IF (GTEST) VNEW(J2)=0.0D0
               IF (STEST.AND.(.NOT.BFGSTST)) THEN
                  DO J3=1,NATOMS
                     HESS(J2,3*(J3-1)+1)=0.0D0
                     HESS(J2,3*(J3-1)+2)=0.0D0
                     HESS(J2,3*(J3-1)+3)=0.0D0
                     HESS(3*(J3-1)+1,J2)=0.0D0
                     HESS(3*(J3-1)+2,J2)=0.0D0
                     HESS(3*(J3-1)+3,J2)=0.0D0
                  ENDDO
                  ! 
                  ! We must not shift here! This needs to be done in shifth.f.
                  ! We need to have correct zero eigenvalues for path.info, for example.
                  ! 
                  ! HESS(J2,J2)=SHIFTV
               ENDIF
            ENDDO
         ENDIF
         IF (STEST) THEN
            IF (FREEZE) THEN
               DO J1=1,NATOMS
                  IF (FROZEN(J1)) THEN
                     IF (VARIABLES) THEN
                        DO J2=1,NOPT
                           HESS(J1,J2)=0.0D0
                           HESS(J2,J1)=0.0D0
                        ENDDO
                     ELSE
                        DO J2=1,NOPT
                           HESS(3*(J1-1)+1,J2)=0.0D0
                           HESS(3*(J1-1)+2,J2)=0.0D0
                           HESS(3*(J1-1)+3,J2)=0.0D0
                           HESS(J2,3*(J1-1)+1)=0.0D0
                           HESS(J2,3*(J1-1)+2)=0.0D0
                           HESS(J2,3*(J1-1)+3)=0.0D0
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         IF (INVERTPT) THEN
            ENERGY=-ENERGY
            IF (GTEST) VNEW(1:NOPT)=-VNEW(1:NOPT)
            IF (SSTEST) HESS(1:NOPT,1:NOPT)=-HESS(1:NOPT,1:NOPT)
         ENDIF

      IF (RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.) ) THEN
        ! sn402: transform back to rigid coords if required
          COORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
          COORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
      ENDIF

!          WRITE(*,'(A,G30.20)') 'Energy in potential:',ENERGY
!          PRINT*,'GTEST,STEST,SSTEST=',GTEST,STEST,SSTEST
!          WRITE(*,'(A,G20.10)') 'RMS in potential=',RMS
!          PRINT*,'coords in potential:'
!          WRITE(*,'(6G25.15)') (COORDS(J1),J1=1,NOPT)
         ! CALL FLUSH(6)
         ! PRINT*,'PARAMS'
         ! WRITE(*,'(3F20.10)') PARAM1,PARAM2,PARAM3
!          IF (GTEST) PRINT*,'grad:'
!          IF (GTEST) WRITE(*,'(6F15.5)') (VNEW(J1),J1=1,NOPT)
!          IF (SSTEST) PRINT*,'hess:'
!          IF (SSTEST) WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)

         CALL MYCPU_TIME(TIME,.FALSE.)
         IF (SSTEST) THEN
            SCALL=SCALL+1
            STIME=STIME+TIME-TIME0
         ELSE IF (GTEST) THEN
            FCALL=FCALL+1
            FTIME=FTIME+TIME-TIME0
         ELSE
            ECALL=ECALL+1
            ETIME=ETIME+TIME-TIME0
         ENDIF
         RETURN

666      WRITE(*,'(A)') ' potential> Error reading CADPAC or GAMES output'
         STOP

      END
      ! 
      ! See Bofill and Comajuan, J. Comp. Chem., 11, 1326, 1995.
      ! PHI=1 is Powell update and PHI=0 is Murtagh-Sargent.
      ! 
      ! SUBROUTINE HUPD(HSAVE,COORDS,COORDSO,GRAD,GRADO,PHI)
      SUBROUTINE HUPD(COORDS,COORDSO,GRAD,GRADO,PHIG)
         USE COMMONS
         USE MODHESS
         IMPLICIT NONE
         INTEGER J1, J2
         DOUBLE PRECISION COORDS(3*NATOMS), COORDSO(3*NATOMS),
     1   GRAD(3*NATOMS), GRADO(3*NATOMS), VECJ(3*NATOMS), DUMMY1,
     2   DTJ, DTD, DUMMY2, PHIG, VECD(3*NATOMS), VECZ(3*NATOMS)
         ! 3                 HSAVE(3*NATOMS,3*NATOMS)

         PRINT*,'WARNING - Hessian updating has not been tested'
         DTD=0.0D0
         DO J1=1,NOPT
            VECD(J1)=COORDS(J1)-COORDSO(J1)
            DTD=DTD+VECD(J1)*VECD(J1)
         ENDDO

         DO J1=1,NOPT
            DUMMY1=0.0D0
            DO J2=1,NOPT
               DUMMY1=DUMMY1+HESS(J1,J2)*VECD(J2)
            ENDDO
            VECJ(J1)=-DUMMY1
         ENDDO

         DTJ=0.0D0
         DO J1=1,NOPT
            VECJ(J1)=GRAD(J1)-GRADO(J1)+VECJ(J1)
            DTJ=DTJ+VECD(J1)*VECJ(J1)
         ENDDO

         DO J1=1,NOPT
            VECZ(J1)=VECD(J1)/DTD - VECJ(J1)/DTJ
         ENDDO

         DO J1=1,NOPT
            DUMMY1=VECJ(J1)/DTJ
            DUMMY2=PHIG*DTJ*VECZ(J1)
            DO J2=J1,NOPT
               ! HESS(J2,J1)=HSAVE(J2,J1)+VECJ(J2)*DUMMY1 - VECZ(J2)*DUMMY2
               ! HESS(J1,J2)=HESS(J2,J1)
               HESS(J2,J1)=HESS(J2,J1)+VECJ(J2)*DUMMY1 - VECZ(J2)*DUMMY2
               HESS(J1,J2)=HESS(J2,J1)
            ENDDO
         ENDDO

         RETURN
      END
