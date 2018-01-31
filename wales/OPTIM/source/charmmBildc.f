
C**************************************************************


      SUBROUTINE CH_SEED(I,J,K,Q,BOND1,BOND2, ANGLE1,
     &   ANGLE2, PHI)
      use modamber9
      use commons, only : NATOMS
      IMPLICIT NONE
 
      DOUBLE PRECISION, INTENT(OUT) :: Q(*)
      DOUBLE PRECISION, INTENT(IN)  :: BOND1(*),BOND2(*),ANGLE1(*),ANGLE2(*),PHI(*)
!      INTEGER, INTENT(IN) :: IC_COORDS(nres*95+5*nres,4)
!      CHARACTER(LEN=3)    :: IC_IMPROP(nres*95+5*nres)
      INTEGER, INTENT(IN) :: I,J,K
      DOUBLE PRECISION              :: RIJ, RJK, THETA
      INTEGER             :: iic
      DOUBLE PRECISION              :: PI
      data PI /3.1415926535897931D0/

      RIJ = 0.0
      RJK = 0.0
      THETA = 0.0
      DO iic = 1,LENIC
         IF(BOND1(iic).GT.0.0) THEN
            IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,2).EQ.J.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,2).EQ.I.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,3).EQ.J.AND.
     &          IC_IMPROP(iic).EQ."IMP")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,3).EQ.I.AND.
     &          IC_IMPROP(iic).EQ."IMP")
     &            RIJ=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,2).EQ.J.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RJK=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,2).EQ.K.AND.
     &          IC_IMPROP(iic).EQ."NOT")
     &            RJK=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,3).EQ.J.AND. 
     &          IC_IMPROP(iic).EQ."IMP")
     &            RJK=BOND1(iic)
            IF(IC_COORDS(iic,1).EQ.J.AND.IC_COORDS(iic,3).EQ.K.AND.
     &          IC_IMPROP(iic).EQ."IMP")
     &            RJK=BOND1(iic)
         ENDIF
         IF(BOND2(iic).GT.0.0) THEN
            IF(IC_COORDS(iic,4).EQ.I.AND.IC_COORDS(iic,3).EQ.J) RIJ=BOND2(iic)
            IF(IC_COORDS(iic,4).EQ.J.AND.IC_COORDS(iic,3).EQ.I) RIJ=BOND2(iic)
            IF(IC_COORDS(iic,4).EQ.K.AND.IC_COORDS(iic,3).EQ.J) RJK=BOND2(iic)
            IF(IC_COORDS(iic,4).EQ.J.AND.IC_COORDS(iic,3).EQ.K) RJK=BOND2(iic)
         ENDIF
C
         IF(IC_COORDS(iic,3).EQ.J) THEN
            IF(ANGLE2(iic).GT.0.0) THEN
               IF(IC_COORDS(iic,4).EQ.I.AND.IC_COORDS(iic,2).EQ.K) 
     &                 THETA=ANGLE2(iic)
               IF(IC_COORDS(iic,4).EQ.K.AND.IC_COORDS(iic,2).EQ.I) 
     &                 THETA=ANGLE2(iic)
            ENDIF
            IF(ANGLE1(iic).GT.0.0) THEN
               IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,2).EQ.K.AND.
     &                IC_IMPROP(iic).EQ."IMP")
     &                  THETA=ANGLE1(iic)
               IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,2).EQ.I.AND.
     &                IC_IMPROP(iic).EQ."IMP")
     &                  THETA=ANGLE1(iic)
            ENDIF
         ELSE
            IF(IC_COORDS(iic,2).EQ.J.AND.IC_IMPROP(iic).EQ."NOT"
     &           .AND.ANGLE1(iic).GT.0.0)THEN
               IF(IC_COORDS(iic,1).EQ.I.AND.IC_COORDS(iic,3).EQ.K)
     &                THETA=ANGLE1(iic)
               IF(IC_COORDS(iic,1).EQ.K.AND.IC_COORDS(iic,3).EQ.I) 
     &                THETA=ANGLE1(iic)
            ENDIF
         ENDIF
      ENDDO
C
      IF(RIJ.EQ.0.0.OR.RJK.EQ.0.0.OR.THETA.EQ.0.0) THEN
         PRINT*, "Error in SEED"
         PRINT '(A,3(1X,I5),3(1X,F7.2))',
     &    ' IC SEED> I,J,K,RIJ,RJK,THETA=',I,J,K,RIJ,RJK,THETA
         RETURN
      ENDIF
C
      Q(3*(I-1)+1) = 0.0
      Q(3*(I-1)+2) = 0.0
      Q(3*(I-1)+3) = 0.0 
      Q(3*(J-1)+1) = RIJ
      Q(3*(J-1)+2) = 0.0
      Q(3*(J-1)+3) = 0.0
      THETA=THETA*(PI/180.0)
      Q(3*(K-1)+1) = RIJ-RJK*COS(THETA)
      Q(3*(K-1)+2) = RJK*SIN(THETA)
      Q(3*(K-1)+3) = 0.0

      RETURN 
      END SUBROUTINE CH_SEED

C****************************************************************************
C
C from PERTDIHE in twist.src
C
C      SUBROUTINE PERTDIHAM(Q,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
C      use modamber9
C      use commons
C      use KEY, ONLY: BHDEBUG, BHSTEPSIZE
C      IMPLICIT NONE
C      DOUBLE PRECISION::  Q(3*NATOMS)
C      DOUBLE PRECISION::  P,ANGLE,DPRAND,MYRANDOM
C      DOUBLE PRECISION::  CHPMIN,CHPMAX
C      INTEGER::           ATOT,A,B,C,I1,J1,IICD, III
!      INTEGER::           RESNUMseg,J3,ISEG,IRES,TOTPHIPSI,TOTSIDECHAIN
!      INTEGER::           CHNMIN,CHNMAX,ISEED
!     DOUBLE PRECISION::    X(NATOMS),Y(NATOMS),Z(NATMS)
!      LOGICAL::   TPP(NATOMS),TS(NATOMS),TO(NATOMS),TT
!      INTEGER::   lenic2
!      DOUBLE PRECISION :: BOND1(nbonh+nbona+nbper), BOND2(nbonh+nbona+nbper),
!     1       THET1(ntheth+ntheta+ngper), THET2(ntheth+ntheta+ngper),
!     2       PHI(nphia+nphih)
!      DOUBLE PRECISION :: ANGLE_SCALE(lenic)
!      INTEGER:: count, P_RESPAR, P_RESMAX
!      DOUBLE PRECISION:: ANGMAX,ANGMIN
!      DOUBLE PRECISION :: PROBABIL
!      DOUBLE PRECISION :: slope
!
!      count=0
!C      write(*,*)'CHPMIN,CHPMAX,CHNMIN,CHNMAX= ',CHPMIN,CHPMAX,CHNMIN,CHNMAX
!
!      !msb50 - preliminary -  how many residues have a lin de/increasing probability
!       P_RESPAR = 6 !number of residues for which twisting prob decreases
!       P_RESMAX = 15
!      !ANGMAX, ANGMIN for angles rescaling - twist more when at chain end
!       ANGMAX= 1.2
!       ANGMIN= 0.6
!
!C  fill IC table with actual Cartesians
!      CALL CHGETICVAL(Q,BOND1, BOND2, THET1, THET2, PHI, .FALSE.)
!      
!C initialise random nmber generator with input seed
!      IF(ISEED.GT.-1) CALL SDPRND(ISEED)
!C
!C will be sent back to 192 if too many or too few dihedrals are altered
!C as determined by CHNMIN and CHNMAX
!192   CONTINUE
!
!      count = count +1 !too avoid endless loop
!C
!C phi/psi angles, omega and sidechain dihedrals are stored in separate lists
!C
!C phi/psi angles
!      B=0
!      ATOT=0
!      DO ISEG=1,NSEG !not at the moment - check twist.src for amendments
!         DO A=1,NPHIPSI(ISEG)
!            TPP(ATOT+A)=.FALSE.
!C
!C  Calculate P, the probability of twisting
!C
!            IF (AMBOLDPERTT) THEN
!                IF (REAL(A).LE.(0.5*(NPHIPSI(ISEG)+1))) THEN
!                   !P=CHPMAX-A*((CHPMAX-CHPMIN)/(NPHIPSI(ISEG)*0.5))
!                   slope = (CHPMAX-CHPMIN)/(0.5*(NPHIPSI(ISEG)-1)-1)
!                   P=-slope*A+(CHPMAX+slope)
!                ELSE
!                   !P=CHPMIN+(A-0.5*NPHIPSI(ISEG))*((CHPMAX-CHPMIN)/(NPHIPSI(ISEG)*0.5))
!                   slope = (CHPMAX-CHPMIN)/(0.5*(NPHIPSI(ISEG)+1))
!                   P = slope*A +(CHPMAX-slope*NPHIPSI(ISEG))
!                END IF
!            ELSE
!            !msb50
!            !take 2* P_RESMAX as NPHIPSI(ISEG) = 2*RESNUMseg(ISEG) 
!            !  one phi, one psi angle per chain
!                P = PROBABIL(A,NPHIPSI(ISEG), CHPMAX,CHPMIN,2*P_RESMAX,2*P_RESPAR)
!               !msb50 - calculation of P depends on position of angle in chain only
!               ! can't rely on TW_ANGLES being space out equally along the chain -->
!               ! would not have P=CHNMIN in the middle of the chain if P= P(NTW_ANGLES)
!                IF (AMBPERTT) THEN
!                   !PRINT '(a15, 2f7.3)', "AMBERPERT probs",P,  TW_DIFFP(PHIPSI(A))
!                   P = (P + TW_DIFFP(PHIPSI(A)))/2!rescaled prob+difference prob-then rescale
!                ENDIF
!
!                ! Calculate scaling factor for angle
!                ANGLE_SCALE(PHIPSI(ATOT+A)) = PROBABIL(A,NPHIPSI(ISEG),
!     &                        ANGMAX,ANGMIN,2*P_RESMAX,2*P_RESPAR)
!            ENDIF
!C            IF(BHDEBUG) WRITE(*,*)'P phipsi = ',P,ATOT+A
!            MYRANDOM=DPRAND()
!            !PRINT '(i4,a24,2f7.3)',PHIPSI(A), "nphipsi probabilities", MYRANDOM,P
!            MYRANDOM=DPRAND()
!            IF (MYRANDOM.LT.P) THEN
!               TPP(ATOT+A)=.TRUE.
!               B=B+1
!            ENDIF
!         ENDDO
!         ATOT=ATOT+NPHIPSI(ISEG)
!      ENDDO
!      TOTPHIPSI=ATOT
!C      IF(BHDEBUG) WRITE(*,*)'PHIPSI: TOT=',ATOT
!
!C amide bond
!      IF (TOMEGAC) THEN
!         ATOT=0
!         DO ISEG=1,NSEG
!            DO A=1,NOMEGAC(ISEG)
!               TO(ATOT+A)=.FALSE.
!C
!C  Calculate P, the probability of twisting
!C
!               IF (REAL(A).LE.(0.5*NOMEGAC(ISEG))) THEN
!                  P=CHPMAX-A*((CHPMAX-CHPMIN)/(NOMEGAC(ISEG)*0.5))
!               ELSE
!                  P=CHPMIN+(A-0.5*NOMEGAC(ISEG))*((CHPMAX-CHPMIN)/(NOMEGAC(ISEG)*0.5))
!               END IF
!C               IF(BHDEBUG) WRITE(*,*)'P omega= ',P
!               MYRANDOM=DPRAND()
!               IF (MYRANDOM.LT.P) THEN
!                  TO(ATOT+A)=.TRUE.
!                  B=B+1
!               ENDIF
!            ENDDO
!            ATOT=ATOT+NOMEGAC(ISEG)
!         ENDDO
!C         IF(BHDEBUG) WRITE(*,*)'OMEGA: TOT=',ATOT
!      ENDIF
!
!C  sidechains
!      ATOT=0
!      DO ISEG=1,NSEG
!         DO A=1,NSIDECHAIN(ISEG)
!            TS(ATOT+A)=.FALSE.
!C
!C  Calculate P, the probability of twisting
!C  Make the probability dependent on the residue number
!C  this means that all dihedrals in the same sidechain have the same P
!C
!            IICD=TW_SIDECHAIN(ATOT+A)
!            j3 = IC_COORDS(IICD,2)
!            DO III = 1,NRES
!               IF (ix(i02+III-1).LE.j3 .AND.(ix(i02+III).GT.j3)) THEN
!               ires = III
!               EXIT
!            ENDIF
!            ENDDO
!            IF (ISEG .GT. 1) THEN !NICTOT different from CHARMM!!
!               RESNUMseg=NICTOT(ISEG)-NICTOT(ISEG-1)     ! number of residues in this segment
!               IRES=IRES-NICTOT(ISEG-1) !no of residue in segment
!            ELSE 
!               RESNUMseg=NICTOT(ISEG)
!            ENDIF
!C            IF(BHDEBUG) WRITE(*,*)'IICD, JAR, IRES, RESNUM : ',IICD,JAR,IRES,RESNUM
!            IF (AMBOLDPERTT) THEN
!                IF (REAL(IRES).LE.(0.5*(RESNUMseg+1))) THEN
!                   P=CHPMAX-IRES*(CHPMAX-CHPMIN)/(RESNUMseg*0.5)
!                   slope = (CHPMAX-CHPMIN)/(0.5*(RESNUMseg-1)-1) !msb50
!                   P=-slope*IRES+(CHPMAX+slope)
!                ELSE
!                   P=CHPMIN+(IRES-0.5*RESNUMseg)*(CHPMAX-CHPMIN)/(RESNUMseg*0.5)
!                   slope = (CHPMAX-CHPMIN)/(0.5*(RESNUMseg+1))
!                   P = slope*IRES +(CHPMAX-slope*RESNUMseg)
!                END IF
!            ELSE
!                P = PROBABIL(IRES,RESNUMseg,CHPMAX,CHPMIN,P_RESMAX, P_RESPAR)
!                IF (AMBPERTT) THEN
!                   !PRINT '(a20, 2f7.3)', "AMBERPERT probs", P, TW_DIFFP(IICD)
!                   P = (P + TW_DIFFP((IICD)))/2 !rescaled prob + difference prob - then rescale
!                ENDIF
!                !  Calculate the angle scaling factor
!                ANGLE_SCALE(IICD) = PROBABIL(IRES,RESNUMseg,
!     &                        ANGMAX,ANGMIN,P_RESMAX,P_RESPAR)
!            ENDIF
!            !PRINT '(i5,a20,i3,a,2f10.7)',IICD, "side probabs for", ires,":", MYRANDOM,P
!
!            MYRANDOM=DPRAND()
!
!C            IF(BHDEBUG) WRITE(*,*)'P sidechain =',P,ATOT+A
!            MYRANDOM=DPRAND()
!            IF (MYRANDOM.LT.P) THEN
!               TS(ATOT+A)=.TRUE.
!               B=B+1
!            ENDIF
!         ENDDO
!         ATOT=ATOT+NSIDECHAIN(ISEG)
!      ENDDO
!      TOTSIDECHAIN=ATOT
!C      IF(BHDEBUG) WRITE(*,*)'SIDECHAIN: TOT=',ATOT
!C
!C      shifting b dihedrals, should be NTEST1 < B < NTEST2
!       IF (B.LT.CHNMIN .OR. B.GT.CHNMAX) THEN
!         WRITE (*,'(A)') 'Too many dihedrals shifted - retrying'
!         IF (count<11) THEN
!            GOTO 192
!         ELSEIF (11.LE.count .AND. count.LE.10000) THEN
!            CHNMAX = CHNMAX +1 !msb50 - quite often tries to twist too many     
!            !if AMBPERTONLY as CHNMAX = NTW_ANGLES << NPHIPSI+NSIDE
!            GOTO 192
!         ELSE
!             PRINT*, "loop in perdiham failed"
!             STOP
!         ENDIF
!       ENDIF
!
!C  twisting phi/psi angles
!       ATOT=0
!       DO ISEG=1,NSEG
!          DO A=1,NPHIPSI(ISEG)
!             IF (TPP(ATOT+A)) THEN
!                IICD=PHIPSI(ATOT+A)
!C                IF(BHDEBUG) PRINT *,'PERTDIHE> changing phipsi ',IICD
!                IF (AMBOLDPERTT) THEN
!                    ANGLE=(DPRAND()-0.5)*2.0*BHSTEPSIZE
!                    PRINT*, "pertdih, changing phipsi",IICD
!                    PRINT '(a20,i4,2f10.3)', "phipsi change", IICD,PHI(IICD),PHI(IICD)+ANGLE
!                ELSE
!                    PRINT*, "pertdih, changing phipsi", IICD
!                    ANGLE=((DPRAND()-0.5)*2.0*BHSTEPSIZE)*ANGLE_SCALE(IICD)
!                    PRINT '(i4,a6,f7.3,a6,f10.3)',IICD,"SCALE",ANGLE_SCALE(IICD),"ANGLE", ANGLE
!                ENDIF
!                PHI(IICD) = PHI(IICD) + ANGLE
!             ENDIF
!          ENDDO
!          ATOT=ATOT+A
!       ENDDO
!C
!C  twisting amide bond`
!       IF (TOMEGAC) THEN
!          ATOT=0
!          DO ISEG=1,NSEG
!             DO A=1,NOMEGAC(ISEG)
!               IF (TO(ATOT+A)) THEN
!                  IICD=OMEGAC(ATOT+A)
!C                  IF(BHDEBUG) PRINT *,'PERTDIHE> changing omega ',IICD
!                  ANGLE=(DPRAND()-0.5)*2.0*BHSTEPSIZE
!                  PHI(IICD) = PHI(IICD) +ANGLE
!                  !PRINT '(a15,i4,2f10.3)', "omega change",IICD, PHI(IICD)-ANGLE,PHI(IICD)
!               ENDIF
!             ENDDO
!             ATOT=ATOT+A
!          ENDDO
!       ENDIF
!C
!C  twisting sidechains
!       ATOT=0
!       DO ISEG=1,NSEG
!          DO A=1,NSIDECHAIN(ISEG)
!             IF (TS(ATOT+A)) THEN
!                IICD=TW_SIDECHAIN(ATOT+A)
!                IF(BHDEBUG) PRINT *,'PERTDIHE> changing sidechain ',IICD
!                IF (AMBOLDPERTT) THEN
!                   ANGLE=(DPRAND()-0.5)*2.0*BHSTEPSIZE
!                   PRINT '(a20,i4,2f10.3)', "side change", IICD,PHI(IICD),PHI(IICD)+ANGLE
!                ELSE
!                   ANGLE=((DPRAND()-0.5)*2.0*BHSTEPSIZE)*ANGLE_SCALE(IICD)
!                   PRINT '(i4,a6,f7.3,a6,f10.3)',IICD,"SCALE",ANGLE_SCALE(IICD),"ANGLE", ANGLE
!                ENDIF
!                PHI(IICD) = PHI(IICD) +ANGLE
!             ENDIF
!          ENDDO
!          ATOT=ATOT+A
!       ENDDO
!
!       IF(PHI(IICD).LT.-180.0) PHI(IICD)=PHI(IICD)+360.0
!       IF(PHI(IICD).GT.180.0) PHI(IICD)=PHI(IICD)-360.0
! 
!       CALL CHREBUILD(Q, BOND1, BOND2,THET1,
!     &        THET2, PHI)
!C      DO II =1, NATOMS
!C         PRINT '(a4,4f11.5)',ih(m04+II-1),Q(3*(II-1)+1),
!C     &          Q(3*(II-1)+2),Q(3*(II-1)+3)
!C      ENDDO
!
!
!     END SUBROUTINE PERTDIHAM

! **********************************************************************


C      SUBROUTINE SETDIHEAM()
C      use modamber9
C      use commons
C      IMPLICIT NONE
C      INTEGER IICD,ISLCT,OLDUSD,I3, J3,K3, L3,IRES,JRES,ISEG,A, NSEGLAST
C      INTEGER AP,AO,AS,AC
C      INTEGER i,j
C      INTEGER SEG_START(NATOMS)
C      LOGICAL LPHIPSI,LOMEGAC,LSIDECHAIN,LCHIRAL,LASTSIDECHAIN
C      INTEGER NPHIPSITOT, NOMEGACTOT,NSIDECHAINTOT, NCHIRALTOT
C      CHARACTER*4 TYPEI,TYPEJ,TYPEK
C      CHARACTER*8 RESLAB
C
C! msb50 - Warning: NICTOT defined differently from CHARMM!
C!         CHARMM - residues per segment: NICTOT(ISEG+1) - NICTOT(ISEG) 
C!              i.e. NICTOT(ISEG) is at which residue new segment starts
C!         AMBER  - residues per segment: NICTOT(ISEG) - NICTOT(ISEG-1) 
C!              i.e. NICTOT(ISEG) is at which residue segment finishes
C
C      IF (.NOT. ALLOCATED(NICTOT)) ALLOCATE(NICTOT(NATOMS))
C
C      IF (.NOT. ALLOCATED(IC_COORDS)) THEN
C         CALL GETICCOORDS()
C         PRINT*, "lenic", lenic
C      ENDIF
C
C      IF (.NOT. ALLOCATED(IS_SIDECHAIN)) THEN
C          ALLOCATE(IS_SIDECHAIN(lenic))
C          DO IICD=1, lenic
C             i3 = IC_COORDS(IICD,1); j3 = IC_COORDS(IICD,2)
C             k3 = IC_COORDS(IICD,3); l3 = IC_COORDS(IICD,4)
C             CALL CHECK_SIDECHAIN(i3,j3,k3,l3,IICD,IS_SIDECHAIN)
C          ENDDO
C      ENDIF
C
C      nseg = 0
C      i = 1
C      !SEG_START _ first atom of new segment
C      SEG_START(1) =1
C      DO WHILE (i .LE. NATOMS)
C      IF (ih(m04+i-1).EQ.'OXT ') THEN !CTerm
C         nseg = nseg +1
C         DO j = 1, nres  !ix(i02+1) = no of atoms in 1st residue +1
C            IF (ix(i02+j-1).LT.i .AND. ix(i02+j).GT.i) THEN
C            NICTOT(nseg) = j !number of residues prior to segment + in segment
C            EXIT
C            ENDIF
C         ENDDO
C
C         i=i+1
C         IF (i .LT. NATOMS) THEN
C            IF (ih(m04+i-1).EQ.'N') THEN
C               SEG_START(nseg+1) = i
C            ELSEIF (ih(m04+i).EQ.'CH3 ') THEN !starts at i+1+1
C               SEG_START(nseg+1) = i-1
C            ELSE
C               PRINT*, "setdiham: segment unknown start", ih(m04+i-1), ih(m04+i-1)
C            ENDIF
C         ENDIF
C
C      ELSEIF (ih(m04+i-1).EQ.'CH3 '.AND.ih(m04+i).EQ.'HH31') THEN !C2Term
C         nseg = nseg +1
C         DO j = 1, nres  !ix(i02+1) 0 no of atoms in 1st residue +1
C            IF (ix(i02+j-1).LT.i .AND. ix(i02+j).GT.i) THEN
C            NICTOT(nseg) = j !number of residues per segment 
C            EXIT
C            ENDIF
C         ENDDO
C         
C         i=i+4
C         IF (i .LT. NATOMS) THEN
C            IF (ih(m04+i-1).EQ.'N') THEN
C               SEG_START(nseg+1) = i
C            ELSEIF (ih(m04+i).EQ.'CH3 ') THEN !starts at i+1+1
C               SEG_START(nseg+1) = i-1
C            ELSE
C               PRINT*, "setdiham: segment unknown start", ih(m04+i-1), ih(m04+i-1)
C            ENDIF
C         ENDIF
C
CC     ELSEIF other terminals
C      
C      ELSE !no terminal
C         i=i+1
C      ENDIF   
C      ENDDO
C
C      IF (nseg .EQ. 0) THEN
C         PRINT*, "in setdiham: no terminal found - unknown type"
C         STOP
C      ENDIF
C
C       write(*,*)'NSEG=',NSEG
C       IF (.NOT. ALLOCATED(NPHIPSI)) ALLOCATE(NPHIPSI(NSEG))
C       IF (.NOT. ALLOCATED(NOMEGAC)) ALLOCATE(NOMEGAC(NSEG))
C       IF (.NOT. ALLOCATED(NSIDECHAIN)) ALLOCATE(NSIDECHAIN(NSEG))
C       IF (.NOT. ALLOCATED(NCHIRAL)) ALLOCATE(NCHIRAL(NSEG))
C       IF (.NOT. ALLOCATED(PHIPSI)) ALLOCATE(PHIPSI(NATOMS))
C       IF (.NOT. ALLOCATED(OMEGAC)) ALLOCATE(OMEGAC(NATOMS))
C       IF (.NOT. ALLOCATED(TW_SIDECHAIN)) ALLOCATE(TW_SIDECHAIN(NATOMS))
C       IF (.NOT. ALLOCATED(CHIRAL)) ALLOCATE(CHIRAL(NATOMS))
C
C       NPHIPSI(1:NSEG)=0 
C       NOMEGAC(1:NSEG)=0
C       NSIDECHAIN(1:NSEG)=0
C       NCHIRAL(1:NSEG)=0
C       PHIPSI(1:NATOMS)=0
C       OMEGAC(1:NATOMS)=0
C       TW_SIDECHAIN(1:NATOMS)=0
C       CHIRAL(1:NATOMS)=0
C
C       AP=0
C       AO=0
C       AS=0
C       AC=0
C
C       NSEGLAST=1
C       ISEG = 1
C       DO IICD=1,LENIC
C          LPHIPSI=.FALSE.
C          LOMEGAC=.FALSE.
C          LSIDECHAIN=.FALSE.
C          LCHIRAL=.FALSE.
C          i3 = IC_COORDS(IICD,1)
C          DO i = 1, NSEG
C             IF (i3 .GE. SEG_START(i)) THEN
C             ISEG =i
C             ENDIF
C          ENDDO
C          IF (ISEG .GT.NSEGLAST) THEN
C             AP=AP+NPHIPSI(NSEGLAST)
C             AO=AO+NOMEGAC(NSEGLAST)
C             AS=AS+NSIDECHAIN(NSEGLAST)
C             AC=AC+NCHIRAL(NSEGLAST)
C             NSEGLAST=ISEG
C          ENDIF
C          CALL ICTYPECHECKAM(LPHIPSI,LOMEGAC,LSIDECHAIN, LCHIRAL,IICD)
C          IF (LPHIPSI) THEN
C             NPHIPSI(ISEG)=NPHIPSI(ISEG)+1
C             PHIPSI(AP+NPHIPSI(ISEG))=IICD
CC            print *,'LP: IICD=', IICD,AP+NPHIPSI(ISEG)
C          ENDIF
C          IF (LOMEGAC) THEN
C             NOMEGAC(ISEG)=NOMEGAC(ISEG)+1
C             OMEGAC(AO+NOMEGAC(ISEG))=IICD
CC            print *,'LO: IICD=', IICD,AO+NOMEGAC(ISEG)
C          ENDIF
C          IF (LSIDECHAIN) THEN
C             NSIDECHAIN(ISEG)=NSIDECHAIN(ISEG)+1
C             TW_SIDECHAIN(AS+NSIDECHAIN(ISEG))=IICD
C            print *,'LS: IICD=', IICD, AS+NSIDECHAIN(ISEG)
C          ENDIF
C          IF (LCHIRAL) THEN
C             NCHIRAL(ISEG)=NCHIRAL(ISEG)+1
C             CHIRAL(AC+NCHIRAL(ISEG))=IICD
CC            print *,'LC: IICD=', IICD, AC+NCHIRAL(ISEG)
C          ENDIF
C       ENDDO
C
C       NPHIPSITOT=0
C       NOMEGACTOT=0
C       NSIDECHAINTOT=0
C       NCHIRALTOT=0
C       DO ISEG=1,NSEG
C          NPHIPSITOT=NPHIPSITOT+NPHIPSI(ISEG)
C          NOMEGACTOT=NOMEGACTOT+NOMEGAC(ISEG)
C          NSIDECHAINTOT=NSIDECHAINTOT+NSIDECHAIN(ISEG)
C          NCHIRALTOT=NCHIRALTOT+NCHIRAL(ISEG)
C          WRITE(*,'(A,I4)')'number of internal coordinates for segment ',ISEG
C          print *,'setdiheam NPHIPSI',NPHIPSI(ISEG)
C          print *,'setdiheam NOMEGAC',NOMEGAC(ISEG)
C          print *,'setdiheam NSIDECHAIN',NSIDECHAIN(ISEG)
C          print *,'setdiheam NCHIRAL',NCHIRAL(ISEG)
C       ENDDO
C       WRITE(*,'(A)')'total number of internal coordinates'
C       print *,'setdiheam> NPHIPSITOT= ',NPHIPSITOT
C       print *,'setdiheam> NOMEGACTOT= ',NOMEGACTOT
C       print *,'setdiheam> NSIDECHAINTOT= ',NSIDECHAINTOT
C       print *,'setdiheam> NCHIRALTOT= ',NCHIRALTOT
C set NSEGATOMS, ie. number of atoms per segment
C       ALLOCATE(NSEGATOMS(NSEG))
C       DO ISEG=2,NSEG+1
C          NSEGATOMS(ISEG-1)=IBASE(NICTOT(ISEG)+1)-IBASE(NICTOT(ISEG-1)+1)
C       ENDDO
C       DO ISEG=1,NSEG
C          WRITE(*,'(A,I4,I6)')'SEG-NR,NSEGATOMS : ',ISEG,NSEGATOMS(ISEG)
C       ENDDO

C      END SUBROUTINE SETDIHEAM




      SUBROUTINE AMB_PATOM(ATOM_NO, ATRES, ATOM)
!msb50: give it a residue number and e.g. "CD1" and in gives you the number of CD1
      use modamber9
      IMPLICIT NONE
      CHARACTER(LEN=4),INTENT(IN) :: ATOM 
      INTEGER,INTENT(IN) :: ATRES
      INTEGER, INTENT(OUT) :: ATOM_NO
      INTEGER ::k
      
!     FOR GOD'S SAKE REMEMBER
!     ix(i02) STARTS WITH 1
!
      ATOM_NO = -9999
!     nat_res =ix(i02+ATRES)- ix(i02+ATRES-1)  !number of atoms in current residue
      DO k=ix(i02+ATRES-1),ix(i02+ATRES)
          IF (ATOM==ih(m04+k-1)) THEN
             ATOM_NO = k
             EXIT
          ENDIF
      ENDDO 
      IF (ATOM_NO.LT.0) THEN
         PRINT*, "ERROR: atom not found", ATOM
      ENDIF
      END SUBROUTINE AMB_PATOM
