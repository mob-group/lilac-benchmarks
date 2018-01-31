! all this is concerning the alignment of start and finish for AMBER (all atom)
! aim is to fix "chirality" around every X - CH_2 - Y carbon so that 
! interpolation is possible
! msb50

      SUBROUTINE FINDCHIRALH(PTEST)
      USE COMMONS,ONLY: NATOMS
      USE MODAMBER9
      USE INTCOMMONS, ONLY: MBBONDNUM, MBADJACENT 

      ! prochiralh - is this atom next to a prochiral one - len natoms, but can 
      ! be true only for hydrogens
      ! prochiralcnt - centreatom which is prochiral
      ! prochiralcnt(H1) = prochiralcnt(H2) = C
      ! prochiralcnt(C) - 4th atom to calculate improper to fix prochirality
      ! prochiralcnt(prochiral_centre) = 4th atom (X or Y, whichever has
      !                                 lowest atom index)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: PTEST
      INTEGER :: ATM1, ATM2, RNUM,I, CNTR, B, A
      CHARACTER(LEN=4):: PHERING(6), PRORING(5), TRPRING5(5),  &
     &     TRPRING6(6), TRPFRNG(6), HISRING(5), RES_RNUM
      LOGICAL :: RINGATM(NATOMS)
      INTEGER :: proch(2), hyd, other, nonhyd(4)

      NBOND = nbona + nbonh
      IF(.NOT.ALLOCATED(IBIB)) THEN
         ALLOCATE(IBIB(nbona+nbonh))
         DO I =1,nbona
            IBIB(I) = ix(iiba+I-1)/3+1
         ENDDO
         DO I=1, nbonh
            IBIB(nbona+I) = ix(iibh+I-1)/3+1
         ENDDO
      ENDIF

      IF(.NOT.ALLOCATED(JBJB)) THEN
          ALLOCATE(JBJB(nbona+nbonh))
          DO I =1,nbona
             JBJB(I) = ix(ijba+I-1)/3+1
          ENDDO
          DO I=1, nbonh
             JBJB(nbona+I) = ix(ijbh+I-1)/3+1
          ENDDO
      ENDIF

      IF (.NOT.ALLOCATED(MBBONDNUM)) ALLOCATE(MBBONDNUM(NATOMS))
      IF (.NOT.ALLOCATED(MBADJACENT)) ALLOCATE(MBADJACENT(NATOMS,6))
      MBBONDNUM(:)=0

      ! ring definitions for known ring residues
      PHERING = (/'CG  ', 'CD1 ', 'CE1 ', 'CZ  ', 'CE2 ', 'CD2 '/)
      PRORING = (/'N   ', 'CA  ', 'CB  ', 'CG  ', 'CD  '/)
      TRPRING5 = (/'CG  ', 'CD2 ', 'CE2 ', 'NE1 ', 'CD1 '/)
      TRPRING6 = (/'CD2 ', 'CE3 ', 'CZ3 ', 'CH2 ', 'CZ2 ', 'CE2 '/)
      TRPFRNG = (/'CE2 ', 'CD2 ', 'NE1 ', 'CG  ', 'CZ2 ', 'CE3 '/)
      HISRING = (/'CG  ', 'ND1 ', 'CE1 ', 'NE2 ', 'CD2 '/)

      DO B=1,NBOND
         ATM1 = IBIB(B); ATM2 = JBJB(B) 
         MBBONDNUM(ATM1) = MBBONDNUM(ATM1) + 1
         MBBONDNUM(ATM2) = MBBONDNUM(ATM2) + 1

         IF (MBBONDNUM(ATM1).GT.6.OR.MBBONDNUM(ATM2).GT.6) THEN
            print*, 'getnatint>> Error! Atom has more than 6 neighbors!'
            print*, ATM1, MBBONDNUM(ATM1), ATM2, MBBONDNUM(ATM2)
            STOP
         ENDIF

         MBADJACENT(ATM1,MBBONDNUM(ATM1)) = ATM2
         MBADJACENT(ATM2,MBBONDNUM(ATM2)) = ATM1
      ENDDO

      RINGATM(:) = .FALSE.
      DO RNUM = 1, NRES
         ! pull out rings from the ring residues
         ! label ring atoms in RINGATM list
         RES_RNUM = ih(m02+RNUM-1)
         IF (RES_RNUM.EQ.'PHE'.OR.RES_RNUM.EQ.'TYR' &
     &      .OR.RES_RNUM.EQ.'NPHE'.OR.RES_RNUM.EQ.'NTYR' &
     &      .OR.RES_RNUM.EQ.'CPHE'.OR.RES_RNUM.EQ.'CTYR') THEN
            DO A = 1,6
               CALL AMB_PATOM(ATM1, RNUM, PHERING(A))
               IF (ATM1.LT.0) THEN
                  print*, 'Error! could not find atom: ', PHERING(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
         ELSE IF (RES_RNUM.EQ.'HIE'.OR.RES_RNUM.EQ.'HID'.OR.  &
     &          RES_RNUM.EQ.'HSP'    &
     &          .OR.RES_RNUM.EQ.'NHIE'.OR.RES_RNUM.EQ.'NHID'.OR. &
     &          RES_RNUM.EQ.'NHSP'.OR.RES_RNUM.EQ.'CHIE'.OR.RES_RNUM &
     &          .EQ.'CHID'.OR. RES_RNUM.EQ.'CHSP') THEN
            DO A = 1,5
               CALL AMB_PATOM(ATM1, RNUM, HISRING(A))
               IF (ATM1.LT.0) THEN
                  print*, 'Error! could not find atom: ', HISRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
         ELSE IF ((RES_RNUM.EQ.'PRO'.OR.RES_RNUM.EQ.'NPRO' &
     &        .OR. RES_RNUM.EQ.'CPRO')) THEN
            DO A = 1,5
               CALL AMB_PATOM(ATM1, RNUM, PRORING(A))
               IF (ATM1.LT.0) THEN
                  print*, 'Error! could not find atom: ', PRORING(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
         ELSE IF (RES_RNUM.EQ.'TRP'.OR.RES_RNUM.EQ.'NTRP'&
     &      .OR.RES_RNUM.EQ.'CTRP') THEN
            DO A = 1,5
               CALL AMB_PATOM(ATM1, RNUM, TRPRING5(A))
               IF (ATM1.LT.0) THEN
                  print*, 'Error! could not find atom: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
            DO A = 1,6
               CALL AMB_PATOM(ATM1, RNUM, TRPRING6(A))
               IF (ATM1.LT.0) THEN
                  print*, 'Error! could not find atom: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
        ENDIF
      ENDDO

      !find prochiral centres: 4 bonds, exactly 2 H's 
      IF (.NOT.ALLOCATED(PROCHIRALH)) ALLOCATE(PROCHIRALH(NATOMS))
      IF (.NOT.ALLOCATED(PROCHIRALCNT)) ALLOCATE(PROCHIRALCNT(NATOMS))
      PROCHIRALH(:)=.FALSE.
      DO CNTR=1, NATOMS
         IF (RINGATM(CNTR)) CYCLE
         IF (MBBONDNUM(CNTR).EQ.4) THEN
            hyd = 0; other = 0
            DO I = 1, MBBONDNUM(CNTR)
               ATM1=MBADJACENT(CNTR,I)
               IF (ih(m04+ATM1-1)(:1).EQ."H") THEN
                  hyd = hyd + 1
                  proch(hyd)=ATM1
               ELSE
                  other = other + 1
                  nonhyd(other) = ATM1
               ENDIF
            ENDDO       
            IF (hyd == 2) THEN
               PROCHIRALH(proch(1)) = .TRUE.
               PROCHIRALH(proch(2)) = .TRUE.
               PROCHIRALCNT(proch(1)) = CNTR
               PROCHIRALCNT(proch(2)) = CNTR
               IF (nonhyd(1) .LT. nonhyd(2)) THEN
                 PROCHIRALCNT(CNTR) = nonhyd(1)
               ELSE
                 PROCHIRALCNT(CNTR) = nonhyd(2) 
                ENDIF 
            ENDIF
         ENDIF
      ENDDO

        !DO I=1,NATOMS
        !  PRINT*, I, "prochiral", PROCHIRALH(I), PROCHIRALCNT(I)
        !ENDDO

      END SUBROUTINE FINDCHIRALH 


!***********************************************************************
      SUBROUTINE ALLCHIRALH_ALIGN(RS,RF,ATM1,ATM2,PTEST)
      use commons, only: NATOMS
      use modamber9      
      use intcommons, only : INTMINPERMT
      IMPLICIT NONE
      
      double precision, intent(in)    :: RS(3*NATOMS)
      double precision, intent(inout) :: RF(3*NATOMS)
      logical, intent(in)             :: ptest
      ! atoms which have to be permuted maintaining prochirality 
      integer, intent(in)             :: atm1, atm2
      double precision                :: dihed1, dihed2
      double precision                :: tmp(3)
      integer                         :: i1, j1
      double precision                :: distance1, distance2, intdist
      double precision                :: COORDS_copy(3*NATOMS), RSTART(3*NATOMS)
      double precision                ::  TMP2(3), TMP3(3)
 
      i1 = PROCHIRALCNT(ATM1)
      j1 = PROCHIRALCNT(i1)
      CALL AMBERDIHEDR(RS, NATOMS,ATM1,ATM2,i1,j1, dihed1)
      CALL AMBERDIHEDR(RF, NATOMS,ATM1,ATM2,i1,j1, dihed2)
      IF (PTEST) PRINT*, "allchiralH, dihed1, dihed2", dihed1, dihed2
      IF (dihed1*dihed2.LT.0) THEN !swap RF hydrogens
         IF (PTEST) PRINT*, "allchiralh> swapp", atm1, atm2
         tmp(:) = RF(3*(ATM1-1)+1:3*ATM1)
         RF(3*(ATM1-1)+1:3*ATM1) = RF(3*(ATM2-1)+1:3*ATM2)
         RF(3*(ATM2-1)+1:3*ATM2) = tmp(:)
      ENDIF

 
      END SUBROUTINE ALLCHIRALH_ALIGN

!*************************************************************************   

! ***********************************************************************
      SUBROUTINE GLYINTPERM(RS, RF, PTEST)
      USE intcommons, only: NDIH, PERMNEIGHBOURS, PERMCHAIN, NGLYDIH, GLYDIH !msb50
      USE commons, only: NATOMS
      IMPLICIT NONE
    
      DOUBLE PRECISION, INTENT(IN)    :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS)
      LOGICAL, INTENT(IN)             :: PTEST
       
      INTEGER          :: III, JJJ
      DOUBLE PRECISION :: GLYSTART(4), GLYFIN(4), GLYDIFF(4), GLYDIST1, GLYDIST
      DOUBLE PRECISION :: dihed1, RFTMP(3*NATOMS)
      INTEGER          :: i1, i2, i3, i4, PLUS

      RFTMP(:) = RF(:)
      DO III=1,NGLYDIH
         DO JJJ =2,3
            i1=GLYDIH(III,JJJ+3)
            i2=GLYDIH(III,JJJ)
            i3=GLYDIH(III,1)
            CALL AMBERDIHEDR(RS,NATOMS,i1,i2,i3,GLYDIH(III,7), dihed1)
            GLYSTART(2*(JJJ-1)-1)=dihed1
            CALL AMBERDIHEDR(RS,NATOMS,i1,i2,i3,GLYDIH(III,8), dihed1)
            GLYSTART(2*(JJJ-1))=dihed1
            CALL AMBERDIHEDR(RF,NATOMS,i1,i2,i3,GLYDIH(III,7), GLYFIN(2*(JJJ-1)-1))
            CALL AMBERDIHEDR(RF,NATOMS,i1,i2,i3,GLYDIH(III,8), GLYFIN(2*(JJJ-1)))
         ENDDO
         GLYDIFF(:)=GLYFIN(:)-GLYSTART(:)
         DO JJJ=1,4
            IF (GLYDIFF(JJJ).GT.180.0D0) THEN
                GLYDIFF(JJJ)=360.0D0-GLYDIFF(JJJ)
            ELSEIF (GLYDIFF(JJJ).LT.-180.0D0) THEN
                GLYDIFF(JJJ)=360+GLYDIFF(JJJ)
            ENDIF
         ENDDO
         GLYDIST1 = DOT_PRODUCT(GLYDIFF,GLYDIFF)
         IF (PTEST) PRINT*, "glydist", GLYDIST1
         RF(3*(GLYDIH(III,7)-1)+1:3*GLYDIH(III,7)) = RFTMP(3*(GLYDIH(III,8)-1)+1:3*GLYDIH(III,8))
         RF(3*(GLYDIH(III,8)-1)+1:3*GLYDIH(III,8)) = RFTMP(3*(GLYDIH(III,7)-1)+1:3*GLYDIH(III,7))
         DO JJJ =2,3
            i1=GLYDIH(III,JJJ+3)
            i2=GLYDIH(III,JJJ)
            i3=GLYDIH(III,1)
            CALL AMBERDIHEDR(RF,NATOMS,i1,i2,i3,GLYDIH(III,7), GLYFIN(2*(JJJ-1)-1))
            CALL AMBERDIHEDR(RF,NATOMS,i1,i2,i3,GLYDIH(III,8), GLYFIN(2*(JJJ-1)))
         ENDDO
         GLYDIFF(:)=GLYFIN(:)-GLYSTART(:)
         DO JJJ=1,4
            IF (GLYDIFF(JJJ).GT.180.0D0) THEN
                GLYDIFF(JJJ)=360.0D0-GLYDIFF(JJJ)
            ELSEIF (GLYDIFF(JJJ).LT.-180.0D0) THEN
                GLYDIFF(JJJ)=360+GLYDIFF(JJJ)
            ENDIF
         ENDDO
         GLYDIST = DOT_PRODUCT(GLYDIFF,GLYDIFF)
         IF (PTEST) PRINT*, "glydist", GLYDIST
         IF (GLYDIST.LT.GLYDIST1) THEN
            IF (PTEST) print*, "keep permutation"
            !RFMP(:) = RF(:)
         ELSE
            RF(:)=RFTMP(:)
         ENDIF
      ENDDO

      END SUBROUTINE GLYINTPERM

! **************************************************************************
      SUBROUTINE MINPERM_CHIRAL(RS, RF,DISTANCE,RMAT,PTEST)
      !minimise Cartesian distance without permuting then
      !permute and calculate Cart distance using a dumb dot product
      !MINPERM would be better, but here prochiral methyl carbons are
      !fixed to have same chirality
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS
      USE modamber9, only:nres, ih, m02,  ix, i02, m04 , PROCHIRALH
      USE commons, only : NATOMS
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS), RMAT(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: DISTANCE
      LOGICAL, INTENT(IN) :: PTEST

      INTEGER :: PG, SW, I, J, K, A1, A2, A3, B1, B2, B3, START,III
      DOUBLE PRECISION :: CURDIST, DIST 
      DOUBLE PRECISION :: RFTMP(3*NATOMS)
      DOUBLE PRECISION :: TMP(3), TMP2(3), TMP3(3)

      CURDIST = SQRT(DOT_PRODUCT(RF-RS, RF-RS))
      START = 1
      DO PG = 1,NPERMGROUP

        IF (NPERMSIZE(PG).EQ.2) THEN
           RFTMP(:) = RF(:)
           ! A1 and A2 are the atom numbers to permute
           A1 = PERMGROUP(START)
           A2 = PERMGROUP(START+1) 
 
           IF (.NOT.PROCHIRALH(A1)) THEN 

              IF (PTEST) print*, 'Swapping atoms: ', A1, A2
              RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
              RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
            
              !move any other groups that must be dragged along
              DO SW = 1,NSETS(PG)
!                A1 = SWAP1(PG,SW)
!                A2 = SWAP2(PG,SW)
                 A1= SETS(PERMGROUP(START),SW)
                 A2= SETS(PERMGROUP(START+1),SW)
            
                 IF (PTEST) print*, 'Dragging swap: ', A1, A2
            
                 RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
                 RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
              ENDDO
           ELSE
               CALL ALLCHIRALH_ALIGN(RS,RF,A1,A2,PTEST)
           ENDIF

           DIST = SQRT(DOT_PRODUCT(RF-RS, RF-RS))
           IF (PTEST) PRINT*, "distances", DIST, "current", CURDIST

           IF (DIST.GT.CURDIST.AND..NOT.PROCHIRALH(A1)) THEN
              IF (PTEST) PRINT*, "undo permutation"
              RF(:) = RFTMP(:)
           ELSEIF (DIST.LT.CURDIST) THEN
              IF (PTEST) PRINT*, "keep permutation"
              CURDIST = DIST 
           ELSEIF (PROCHIRALH(A1)) THEN
              CURDIST = DIST
           ENDIF

        ELSE IF (NPERMSIZE(PG).EQ.3) THEN
          A1 = PERMGROUP(START)
          A2 = PERMGROUP(START+1)
          A3 = PERMGROUP(START+2)

          TMP(:) = RF(3*(A1-1)+1:3*A1)
          TMP2(:) = RF(3*(A2-1)+1:3*A2)
          TMP3(:) = RF(3*(A3-1)+1:3*A3)

          RFTMP(:) = RF(:)
          DO I = 0,2
             DO J = 0,2
                IF (J.EQ.I.OR.(I.EQ.0.AND.J.EQ.1)) CYCLE
                DO K = 0,2
                   IF (K.EQ.I.OR.K.EQ.J) CYCLE

                   B1 = PERMGROUP(START+I)
                   B2 = PERMGROUP(START+J)
                   B3 = PERMGROUP(START+K)

                   IF (PTEST) print*, 'permuting: ', B1, B2, B3

                   RF(3*(B1-1)+1:3*B1) = TMP(:)
                   RF(3*(B2-1)+1:3*B2) = TMP2(:)
                   RF(3*(B3-1)+1:3*B3) = TMP3(:)

                   DIST = SQRT(DOT_PRODUCT(RF-RS, RF-RS))

                   IF (PTEST) print*, 'new distance: ', DIST

                   IF (DIST.LT.CURDIST) THEN
                      RFTMP(:) = RF(:)
                      CURDIST = DIST
                      IF (PTEST) print*, 'keep permutation'
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          RF(:) = RFTMP(:)
       ELSE
          print*, 'Error! INTMINPERM is only set up for permutation groups &
               &          of size 2 & 3 for now', PG, NPERMSIZE(PG)
          STOP
       ENDIF

       START = START + NPERMSIZE(PG)
       ENDDO

       DISTANCE = DOT_PRODUCT(RF-RS,RF-RS)
       IF (PTEST) PRINT*, "minperm_chiral> final distance", distance


      END SUBROUTINE MINPERM_CHIRAL
