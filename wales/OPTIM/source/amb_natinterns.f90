! msb50 --> efk_ CHARMM/natinterns.src


      SUBROUTINE AMBGETNATINTERN()
! generate the lists specifying the natural internal coordinates for the system
      USE COMMONS, ONLY : NATOMS
      USE MODAMBER9
      USE INTCOMMONS
      
      IMPLICIT NONE

      INTEGER A, B, I, R
      CHARACTER*4 PHERING(6), PRORING(5), TRPRING5(5),  &
     &     TRPRING6(6), TRPFRNG(6), HISRING(5), TMPTYPE, &
     &     CRING(6), ARING5(5), ARING6(6), AFRNG(6), RRING(5)
      INTEGER RINGID(NATOMS)   !which ring is this atoms part of
      INTEGER ATM1, ATM2, RNUM
      INTEGER PATATM, PATRES
      CHARACTER PTYP
      INTEGER NTERM
      CHARACTER(LEN=4) :: RES_RNUM
      LOGICAL ALLRING, RINGATM(NATOMS)
      LOGICAL BBATM1, BBATM2, BBONLY
      INTEGER GLYCA(NATOMS), nindex !msb50
      LOGICAL not_nterm, ringneighb

      ! testing only
      INTEGER TEST

      GLYCA(:)=0
      !this array is to show which CA's are glycine i.e. not backbone

      !msb50 - try to fake CHARMM
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

!      !msb50
!      IF(.NOT.ALLOCATED(IBIB)) THEN
!      ALLOCATE(IBIB(nbona+nbonh))
!      ALLOCATE(JBJB(nbona+nbonh))
!      IBIB(1) = 1; JBJB(1) = 2
!      IBIB(2) = 2; JBJB(2) = 4
!      IBIB(3) = 2; JBJB(3) = 3
!      IBIB(4) = 4; JBJB(4) = 6
!      IBIB(5) = 4; JBJB(5) = 5
!      ENDIF

      CALL AMB_GETBACKBONE
      IF (.NOT.ALLOCATED(MBBONDNUM)) ALLOCATE(MBBONDNUM(NATOMS))
      IF (.NOT.ALLOCATED(MBADJACENT)) ALLOCATE(MBADJACENT(NATOMS,6)) 
! msb50 deallocate in FINDPERMDIH
      MBBONDNUM(:)=0
      NCNT2 = 0
      NGLYDIH=0     
 
      ! ring definitions for known ring residues
      RRING = (/"C4' ", "O4' ", "C1' ", "C2' ", "C3' "/)
      ARING5 = (/"N9  ", "C8  ", "N7  ", "C5  ", "C4  "/)
      ARING6 = (/"C5  ", "C6  ", "N1  ", "C2  ", "N3  ", "C4  "/)
      AFRNG = (/"C4  ", "C5  ", "N7  ", "N9  ", "N3  ", "C6  "/)
      CRING = (/'N1  ', 'C6  ', 'C5  ', 'C4  ', 'N3  ', 'C2  '/)
      !GRING5 = (/'N9 ', 'C8  ', 'N7  ', 'C5  ', 'C4  '/)
      !GRING6 = (/'C4  ', 'C5  ', 'C6  ', 'N1  ', 'C2  ', 'N3  '/)
      !GFRNG = (/'C4  ', 'C5  ', 'N9  ', 'N7  ', 'N3  ', 'C6  '/)
      !URING = (/'N1  ', 'C6  ', 'C5  ', 'C4  ', 'N3  ', 'C2  '/)
      PHERING = (/'CG  ', 'CD1 ', 'CE1 ', 'CZ  ', 'CE2 ', 'CD2 '/)
      PRORING = (/'N   ', 'CA  ', 'CB  ', 'CG  ', 'CD  '/)
      TRPRING5 = (/'CG  ', 'CD2 ', 'CE2 ', 'NE1 ', 'CD1 '/)
      TRPRING6 = (/'CD2 ', 'CE3 ', 'CZ3 ', 'CH2 ', 'CZ2 ', 'CE2 '/)
      TRPFRNG = (/'CE2 ', 'CD2 ', 'NE1 ', 'CG  ', 'CZ2 ', 'CE3 '/)
      HISRING = (/'CG  ', 'ND1 ', 'CE1 ', 'NE2 ', 'CD2 '/) 

      IF (NUBONDS.GT.0) THEN
         ! add on user-specified bonds
         DO B = 1,NUBONDS
            IBIB(NBOND+B) = UBONDS(B,1)
            JBJB(NBOND+B) = UBONDS(B,2)
         ENDDO
         NBOND = NBOND + NUBONDS
      ENDIF

      ! list number of bonds and adjacent atoms for each atom
      NBDS = 0
      DO B = 1,NBOND
         ATM1 = IBIB(B); ATM2 = JBJB(B)

         IF (.NOT.USECART(ATM1).OR..NOT.USECART(ATM2)) NBDS = NBDS + 1

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

      NINTC = NBDS + 3*NCRT

      NRNG = 0; NFRG = 0; NDIH = 0
      RINGATM(:) = .FALSE.
      RINGS(:,:) = -1
      RINGID(:) =-1
      DO RNUM = 1,CARTRESSTART - 1
         ! pull out rings from the ring residues
         ! label ring atoms in RINGATM list
         RES_RNUM = ih(m02+RNUM-1)
         IF (RES_RNUM.EQ.'PHE'.OR.RES_RNUM.EQ.'TYR' &
     &      .OR.RES_RNUM.EQ.'NPHE'.OR.RES_RNUM.EQ.'NTYR' &
     &      .OR.RES_RNUM.EQ.'CPHE'.OR.RES_RNUM.EQ.'CTYR') THEN
            NINTC = NINTC + 6
            NDIH = NDIH + 6 ! keep track of individual dihedrals
            NRNG = NRNG+1
            DO A = 1,6               
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, PHERING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', PHERING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF (RES_RNUM.EQ.'HIE'.OR.RES_RNUM.EQ.'HID'.OR.  &
     &          RES_RNUM.EQ.'HSP'    &
     &          .OR.RES_RNUM.EQ.'NHIE'.OR.RES_RNUM.EQ.'NHID'.OR. &
     &          RES_RNUM.EQ.'NHSP'.OR.RES_RNUM.EQ.'CHIE'.OR.RES_RNUM &
     &          .EQ.'CHID'.OR. RES_RNUM.EQ.'CHSP') THEN
            NINTC = NINTC + 4
            NRNG = NRNG + 1
            NDIH = NDIH + 5
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, HISRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', HISRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF ((RES_RNUM.EQ.'PRO'.OR.RES_RNUM.EQ.'NPRO' &
     &        .OR. RES_RNUM.EQ.'CPRO').AND..NOT.BBCART) THEN
            NINTC = NINTC + 4
            NRNG = NRNG + 1 
            NDIH = NDIH + 5
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, PRORING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', PRORING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
!CJF41 added a section to deal with hydroxyproline for collagen; treating the same as proline
         ELSE IF ((RES_RNUM.EQ.'HYP'.OR.RES_RNUM.EQ.'NHYP' &
     &        .OR. RES_RNUM.EQ.'CHYP').AND..NOT.BBCART) THEN
            NINTC = NINTC + 4
            NRNG = NRNG + 1 
            NDIH = NDIH + 5
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, PRORING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', PRORING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF (RES_RNUM.EQ.'TRP'.OR.RES_RNUM.EQ.'NTRP'&
     &      .OR.RES_RNUM.EQ.'CTRP') THEN
            NINTC = NINTC + 11
            NRNG = NRNG + 1
            NDIH = NDIH + 13
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, TRPRING5(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
            NRNG = NRNG + 1 
            DO A = 1,6
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, TRPRING6(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG -1 !have same ringid!!!
            ENDDO
            NFRG = NFRG + 1
            DO A = 1,6
               CALL AMB_PATOM(FRINGS(NFRG,A), RNUM, TRPFRNG(A))
               IF (FRINGS(NFRG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
            ENDDO
         !msb50 glycine
         ELSE IF (RES_RNUM.EQ.'GLY'.OR.RES_RNUM.EQ.'NGLY'&
     &       .OR.RES_RNUM.EQ.'CGLY') THEN
             DO A = ix(i02+RNUM-1), ix(i02+RNUM)-1
                IF (ih(m04+A-1).EQ.'CA  ') THEN
                  GLYCA(NGLYDIH+1)=A
                  NGLYDIH=NGLYDIH+1
                ENDIF
             ENDDO
!msb50 rna
         ELSE IF (RES_RNUM.EQ.'RA  '.OR.RES_RNUM.EQ.'RA3 '&
             .OR.RES_RNUM.EQ.'RA5 '.OR.RES_RNUM.EQ.'RAN ' &
             .OR. RES_RNUM.EQ.'RG  '.OR.RES_RNUM.EQ.'RG3 ' &
             .OR.RES_RNUM.EQ.'RG5 '.OR.RES_RNUM.EQ.'RGN ') THEN
             NINTC = NINTC +11
             NRNG = NRNG + 1
             NDIH = NDIH + 13
             DO A = 1,5
                CALL AMB_PATOM(RINGS(NRNG,A), RNUM, ARING5(A))
                IF (RINGS(NRNG, A).LT.0)  THEN
                   print*, 'Error! could not find atom: ', ARING5(A), RNUM
                   STOP
                ENDIF
                RINGATM(RINGS(NRNG,A)) =.TRUE.
                RINGID(RINGS(NRNG,A)) = NRNG
             ENDDO
             NRNG = NRNG + 1
             DO A =1,6
                CALL AMB_PATOM(RINGS(NRNG, A), RNUM, ARING6(A))
                IF (RINGS(NRNG, A).LT.0) THEN
                   print*, 'Error! could not find atom: ', ARING6(A), RNUM
                   STOP
                ENDIF
                RINGATM(RINGS(NRNG,A)) =.TRUE.
                RINGID(RINGS(NRNG,A)) = NRNG - 1 ! have same ringID as same ring
             ENDDO
             NFRG = NFRG + 1
             DO A = 1,6
               CALL AMB_PATOM(FRINGS(NFRG,A), RNUM,  AFRNG(A))
               IF (FRINGS(NFRG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', AFRNG(A),RNUM
                  STOP
               ENDIF
             ENDDO
!ribose
             NINTC = NINTC + 4
             NRNG = NRNG + 1
             NDIH = NDIH + 5
             DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, RRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', RRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF (RES_RNUM.EQ.'RC  '.OR.RES_RNUM.EQ.'RC3 '&
             .OR.RES_RNUM.EQ.'RC5 '.OR.RES_RNUM.EQ.'RCN ' &
             .OR. RES_RNUM.EQ.'RU  '.OR.RES_RNUM.EQ.'RU3 ' &
             .OR.RES_RNUM.EQ.'RU5 '.OR.RES_RNUM.EQ.'RUN ') THEN
            NINTC = NINTC + 6
            NDIH = NDIH + 6 ! keep track of individual dihedrals
            NRNG = NRNG+1
            DO A = 1,6
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, CRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', CRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
!ribose
             NINTC = NINTC + 4
             NRNG = NRNG + 1
             NDIH = NDIH + 5
             DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, RRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  print*, 'Error! could not find atom: ', RRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
             ENDDO

         ENDIF
      ENDDO
 
      !store angles for glycine
      IF (.NOT.ALLOCATED(GLYDIH)) ALLOCATE(GLYDIH(NGLYDIH,8))
      DO I=1,NGLYDIH
        nindex = 0
        not_nterm = .FALSE.
        GLYDIH(I,1) = GLYCA(I)
        A = GLYCA(I)  
        DO B=1,MBBONDNUM(GLYCA(I)) 
           IF (MBBONDNUM(MBADJACENT(A,B)).EQ.1) THEN
              GLYDIH(I,8-nindex) = MBADJACENT(A,B)
              nindex = nindex +1
           ELSEIF (ih(m04+MBADJACENT(A,B)-1).EQ.'C   ') THEN
              GLYDIH(I,2) = MBADJACENT(A,B)
              GLYDIH(I,4) = MBADJACENT(MBADJACENT(A,B),1)
              GLYDIH(I,5) = MBADJACENT(MBADJACENT(A,B),2)
           ELSE !N
              GLYDIH(I,3) = MBADJACENT(A,B)
              DO R = 1,MBBONDNUM(MBADJACENT(A,B))
                IF (MBBONDNUM(MBADJACENT(MBADJACENT(A,B),R)).GT.1&
     &           .AND.MBADJACENT(MBADJACENT(A,B),R).NE.A) THEN 
                   GLYDIH(I,6) = MBADJACENT(MBADJACENT(A,B),R)
                   not_nterm = .TRUE.
                ENDIF
              ENDDO
              IF (.NOT.not_nterm) THEN
                 IF (MBADJACENT(MBADJACENT(A,B),1).NE.A) THEN
                    GLYDIH(I,6) = MBADJACENT(MBADJACENT(A,B),1)
                 ELSE
                    GLYDIH(I,6) = MBADJACENT(MBADJACENT(A,B),2)
                 ENDIF
              ENDIF
            ENDIF
         ENDDO
      ENDDO
      

      ! put in the extra (user-specified) rings
      DO R = 1,NURINGS
         NRNG = NRNG + 1
         NDIH = NDIH + URINGS(R,0)
         IF (URINGS(R,0).EQ.5) THEN
            NINTC = NINTC + 4
         ELSE IF (URINGS(R,0).EQ.6) THEN
            NINTC = NINTC + 6
         ELSE
            print*, 'getnatinterns>> Error!&
              & User-specified rings must have 5 or 6 atoms'
            STOP
         ENDIF
            
         DO I = 1, URINGS(R,0)
            RINGS(NRNG,I) = URINGS(R,I)
            RINGATM(URINGS(R,I)) = .TRUE.
         ENDDO
      ENDDO

      ! put in the non-terminal centers
      RNUM = 0
      NCNT = 0; NIMP = 0
      DO A = 1,CARTATMSTART - 1
         IF (RNUM.LT.NRES) THEN
            IF (A.EQ.ix(i02+RNUM)) THEN
                RNUM = RNUM + 1 !new residue started
             ENDIF
         ENDIF
         RES_RNUM = ih(m02+RNUM-1)        

         IF (MBBONDNUM(A).EQ.0) THEN
            print*, 'getnatint>> Error! Disconnected atom.', A
            STOP
         ELSE IF (MBBONDNUM(A).EQ.1) THEN
            CYCLE ! ignore terminal atoms
         ENDIF

         IF (RINGATM(A)) THEN
            ! ignore centers where all substituents are ring atms
            ALLRING = .TRUE.
            DO I = 1,MBBONDNUM(A)
               IF (.NOT.RINGATM(MBADJACENT(A,I)) .OR. &
     &            (RINGID(A).NE.RINGID(MBADJACENT(A,I)))) THEN
                  ALLRING = .FALSE.              
                  EXIT
               ENDIF
            ENDDO
            IF (ALLRING) THEN
               CYCLE 
            ENDIF
         ENDIF

         NCNT = NCNT + 1
         CENTERS(NCNT,1) = A

         NTERM = 0 ! number of terminal or quasi-terminal neighbors        
         RINGNEIGHB = .FALSE.
         DO I=1, MBBONDNUM(A)
            IF (RINGATM(A).AND.(RINGID(MBADJACENT(A,I)).LT.RINGID(A))) THEN
               RINGNEIGHB=.TRUE. 
               CYCLE
            ENDIF
         ENDDO
         !previous loop only determined whether two rings are linked

         DO I = 1,MBBONDNUM(A)
            ! treat all ring substituents as terminal
            ! CB is quasiterminal for CA center
            ! the NH's of ARG are quasi-terminal
            !msb50 make o5 quasiterminal for p to make p methine sp3
            ! linked rings: linkatom in following ring quasi-terminal
            ! its neighbours also quasiterminal to fix following ring position
     !&           .OR.(RINGATM(A).AND.RINGNEIGHB.AND.                    &
     !&              (RINGID(MBADJACENT(A,I)).EQ.RINGID(A))).OR.         &
            IF (MBBONDNUM(MBADJACENT(A,I)).EQ.1.OR.                     & 
     &           (RINGATM(A).AND..NOT.RINGATM(MBADJACENT(A,I))).OR.     &
     &           (RINGATM(A).AND.(RINGID(MBADJACENT(A,I)).NE.RINGID(A))).OR.&
     &           ih(m04+A-1).EQ.'CA'.AND.ih(m04+(MBADJACENT(A,I))-1).EQ.&
     &           'CB'.OR.                                               &
     &           (RES_RNUM.EQ.'ARG'.AND.                                &
     &           (ih(m04+(MBADJACENT(A,I))-1).EQ.'NH1'.OR.              &
     &           ih(m04+(MBADJACENT(A,I))-1).EQ.'NH2')).OR.             &
     &           (ih(m04+A-1).EQ.'P   '.AND.                            &
     &           ih(m04+(MBADJACENT(A,I))-1).EQ."O5' ")) THEN
               
               NTERM = NTERM + 1
               ! put terminal atms at end
               CENTERS(NCNT,MBBONDNUM(A)-NTERM+2) = MBADJACENT(A,I)
               !msb50 - proline exception
               IF ((RES_RNUM.EQ.'PRO'.OR.RES_RNUM.EQ.'NPRO'.OR.         &
     &         RES_RNUM.EQ.'CPRO').AND.ih(m04+MBADJACENT(A,I)-1).EQ.'HA' &
     &             .AND.(MBBONDNUM(A)-NTERM+2).EQ.3) THEN
                    CENTERS(NCNT, 3)= CENTERS(NCNT, 4)
                    CENTERS(NCNT, 4) = MBADJACENT(A,I)
               ENDIF
               !cjf41 added hydroxyproline exception - treat the same as proline
               IF ((RES_RNUM.EQ.'HYP'.OR.RES_RNUM.EQ.'NHYP'.OR.         &
     &         RES_RNUM.EQ.'CHYP').AND.ih(m04+MBADJACENT(A,I)-1).EQ.'HA' &
     &             .AND.(MBBONDNUM(A)-NTERM+2).EQ.3) THEN
                    CENTERS(NCNT, 3)= CENTERS(NCNT, 4)
                    CENTERS(NCNT, 4) = MBADJACENT(A,I)
               ENDIF
 
            ELSE
               CENTERS(NCNT,I-NTERM+1) = MBADJACENT(A,I)
            ENDIF
         ENDDO
          

         ! set center type
         IF (RINGATM(A)) THEN
            IF (MBBONDNUM(A).EQ.4) THEN
              ! IF (RINGNEIGHB.AND.NTERM == 2 ) THEN
              !    CENTERS(NCNT,0) = 2
              !    NINTC = NINTC + 5 
              ! ELSE IF (RINGNEIGHB.AND.NTERM == 3) THEN
              !    CENTERS(NCNT,0) = 1
              !    NINTC = NINTC + 5 
              ! ELSE !normal - no linked rings
                  CENTERS(NCNT,0) = 8
                  NINTC = NINTC + 4
              ! ENDIF
            ELSE IF (MBBONDNUM(A).EQ.3) THEN
               CENTERS(NCNT,0) = 7
               NINTC = NINTC + 2
               NIMP = NIMP + 1
               IMPDIH(NIMP,1) = CENTERS(NCNT,1)
               IMPDIH(NIMP,2) = CENTERS(NCNT,4)
               IMPDIH(NIMP,3:4) = CENTERS(NCNT,2:3)
            ELSE
               print*, 'ERROR! can have at most two ring substituents'
               print*, A, MBBONDNUM(A)
               STOP
            ENDIF
         ELSE IF (MBBONDNUM(A).EQ.4) THEN
            NINTC = NINTC + 5
            SELECT CASE(NTERM)
            CASE (3)
               CENTERS(NCNT,0) = 1
            CASE (2)
               CENTERS(NCNT,0) = 2
               IF (INTMINPERMT) THEN !msb50
                   IF (.NOT.ALLOCATED(CENTER2)) THEN
                       ALLOCATE(CENTER2(NATOMS))
                       CENTER2(:) = 0
                   ENDIF
                   IF (ih(m04+A-1).NE.'CA  ') THEN!no CAs, only in am acid sidechains
                      NCNT2 = NCNT2+1
                      CENTER2(NCNT2) = NCNT 
                   ENDIF
               ENDIF
            CASE (1)
               CENTERS(NCNT,0) = 3
            CASE DEFAULT
               print*, 'Error! unknown center type. Atm,bondnum,Nterm: '&
     &             ,A, MBBONDNUM(A), NTERM
               STOP
            END SELECT
         ELSE IF (MBBONDNUM(A).EQ.3) THEN           
            NIMP = NIMP + 1
            NINTC = NINTC + 3
            SELECT CASE(NTERM)
            CASE (2)
               CENTERS(NCNT,0) = 4
               IMPDIH(NIMP,1:4) = CENTERS(NCNT,1:4)
            CASE (1)
               CENTERS(NCNT,0) = 5
               IMPDIH(NIMP,1) = CENTERS(NCNT,1)
               IMPDIH(NIMP,2) = CENTERS(NCNT,4)
               IMPDIH(NIMP,3:4) = CENTERS(NCNT,2:3)
            CASE DEFAULT
               print*, 'Error! unknown center type. Atm,bondnum,Nterm: ', &
     &              A, MBBONDNUM(A), NTERM
               STOP
            END SELECT
         ELSE IF (MBBONDNUM(A).EQ.2) THEN
            NINTC = NINTC + 1
            CENTERS(NCNT,0) = 6
         ENDIF

         ! counts as backbone center only if first 3 atoms are in backbone
         IF (BBCART) THEN
            BBONLY = .TRUE.
            DO I = 1,3
               IF (.NOT.AM_BACKBONE(CENTERS(NCNT,I))) BBONLY = .FALSE.
            ENDDO
           
            IF (BBONLY) THEN
               SELECT CASE (CENTERS(NCNT,0))
               CASE(5)          !treat like ring substituent
                  CENTERS(NCNT,0) = 7
                  NINTC = NINTC - 1
               CASE(2)          ! treat like ring substituent
                  CENTERS(NCNT,0) = 8
                  NINTC = NINTC - 1
               CASE(6)          !skip coord
                  NINTC = NINTC -1
                  NCNT = NCNT - 1
               CASE(7)          !skip coord
                  NINTC = NINTC - 2
                  NCNT = NCNT - 1
               CASE DEFAULT
                  print*, 'ERROR! If BBCART set,'
                  print*, 'expect all backbone centers to be'
                  print*, 'of type 2,5,6,or 7', RNUM, RES_RNUM
                  print*, CENTERS(NCNT,:)
                  print*, "msb50 - and 8 - fix this if it's 8"
                  STOP
               END SELECT
            ENDIF
         ENDIF
      ENDDO

      ! get linear dihedrals
      NLDH = 0
      DO B = 1,NBOND
         ATM1 = IBIB(B); ATM2 = JBJB(B)
 
         ! ignore bonds involving a terminal atom, bonds in rings
         IF ((MBBONDNUM(ATM1).EQ.1.OR.MBBONDNUM(ATM2).EQ.1).OR.&
     &        (RINGATM(ATM1).AND.RINGATM(ATM2).AND.&
     &         RINGID(ATM1).EQ.RINGID(ATM2))) CYCLE        
             !lindh required if linked rings 

         NLDH = NLDH + 1
         LINDIH(NLDH,1) = ATM1; LINDIH(NLDH,2) = ATM2


         BBATM1 = ih(m04+ATM1-1).EQ.'N'.OR.ih(m04+ATM1-1).EQ.'CA'.OR.   & 
     &       ih(m04+ATM1-1).EQ.'C'.OR.ih(m04+ATM1-1).EQ.'CH3 '.OR.      &
     &       ih(m04+ATM1-1).EQ."C5' ".OR.ih(m04+ATM1-1).EQ."C3' ".OR.   &
     &       ih(m04+ATM1-1).EQ."O3' ".OR.ih(m04+ATM1-1).EQ."P"
         BBATM2 = ih(m04+ATM2-1).EQ.'N'.OR.ih(m04+ATM2-1).EQ.'CA'.OR.   &
     &        ih(m04+ATM2-1).EQ.'C'.OR.ih(m04+ATM2-1).EQ.'CH3 '.OR.     &
     &       ih(m04+ATM2-1).EQ."C5' ".OR.ih(m04+ATM2-1).EQ."C3' ".OR.   &
     &       ih(m04+ATM2-1).EQ."O3' ".OR.ih(m04+ATM2-1).EQ."P"
         BBONLY = .FALSE.


         ! when applicable, use backbone bonds only
         IF (BBATM1.AND.BBATM2) THEN
            BBONLY = .FALSE.
            DO I = 1,MBBONDNUM(ATM1)
               TMPTYPE = ih(m04+(MBADJACENT(ATM1,I))-1)
               IF (TMPTYPE.NE.ih(m04+ATM2-1).AND.                       &
     &              (TMPTYPE.EQ.'N'.OR.TMPTYPE.EQ.'CA'.OR.              &
     &              TMPTYPE.EQ.'C'.OR.TMPTYPE.EQ.'CH3 '.OR.             &
     &              TMPTYPE.EQ."C5' ".OR.TMPTYPE.EQ."C3' ".OR.           &
     &              TMPTYPE.EQ."O3' ".OR.TMPTYPE.EQ."P   ".OR.          &
     &              TMPTYPE.EQ."O5' ")) THEN
                  LINDIH(NLDH,-1) = 1
                  LINDIH(NLDH,3) = MBADJACENT(ATM1,I)
                  BBONLY = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         IF (.NOT.BBONLY) THEN
            LINDIH(NLDH,-1)= MBBONDNUM(ATM1) - 1
            A = 2
            DO I = 1,MBBONDNUM(ATM1)
               IF (MBADJACENT(ATM1,I).NE.ATM2) THEN
                  A = A + 1
                  LINDIH(NLDH,A) = MBADJACENT(ATM1,I)
               ENDIF
            ENDDO
         ENDIF

         IF (BBATM1.AND.BBATM2) THEN
            BBONLY = .FALSE.
            DO I = 1,MBBONDNUM(ATM2)
               TMPTYPE = ih(m04+MBADJACENT(ATM2,I)-1)
               IF (TMPTYPE.NE.ih(m04+ATM1-1).AND.                       &
     &              (TMPTYPE.EQ.'N'.OR.TMPTYPE.EQ.'CA'.OR.              &
     &                TMPTYPE.EQ.'C'.OR.TMPTYPE.EQ.'CH3 '.OR.           &
     &              TMPTYPE.EQ."C5' ".OR.TMPTYPE.EQ."C3' ".OR.          &
     &              TMPTYPE.EQ."O3' ".OR.TMPTYPE.EQ."P   ".OR.          &
     &              TMPTYPE.EQ."O5' ")) THEN 
                  LINDIH(NLDH,0) = 1
                  LINDIH(NLDH,LINDIH(NLDH,-1)+3) = MBADJACENT(ATM2,I)
                  BBONLY = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         IF (.NOT.BBONLY) THEN
            LINDIH(NLDH,0)= MBBONDNUM(ATM2) - 1
            A = LINDIH(NLDH,-1)+2
            DO I = 1,MBBONDNUM(ATM2)
               IF (MBADJACENT(ATM2,I).NE.ATM1) THEN
                  A = A + 1
                  LINDIH(NLDH,A) = MBADJACENT(ATM2,I)
               ENDIF
            ENDDO
         ENDIF

         ! ignore if using cartesians for all atoms involved
         BBONLY = .TRUE.
         DO I = 1,LINDIH(NLDH,-1)+LINDIH(NLDH,0)+2
            IF (.NOT.USECART(LINDIH(NLDH,I))) THEN
               BBONLY = .FALSE.
            ENDIF
         ENDDO
      
         IF (BBONLY) THEN
            NLDH = NLDH - 1
         ELSE
            NDIH = NDIH + LINDIH(NLDH,-1)*LINDIH(NLDH,0)
            NINTC = NINTC + 1
         ENDIF
                  
      ENDDO

      IF (printcoords) THEN
         print*, 'TOTAL COORDS:', NINTC
         print*, 'CARTESIAN ATOMS:', NCRT

         DO I = 1, NCRT
            print*, CARTATMS(I)
         ENDDO

         print*, 'BONDS:', NBDS
         DO I=1,NBOND
            IF (.NOT.(USECART(IBIB(I)).AND.USECART(JBJB(I))))            &
     &            print*, IBIB(I), JBJB(I)
         ENDDO
         print*, 'CENTERS:', NCNT
         DO I=1,NCNT
            print '(I6,A,6I6)', I, '|', CENTERS(I,0:5)
         ENDDO
         print*, 'RINGS:', NRNG
         DO I=1,NRNG
            print*, RINGS(I,1:6)
            print*, '---'
         ENDDO
         print*, 'FUSED RINGS:', NFRG
         DO I=1,NFRG
            print*, FRINGS(I,1:6)
            print*, '---'
         ENDDO
         print*, 'LINEAR DIHEDRALS', NLDH, NDIH
         DO I=1,NLDH
            print*, I, '|', LINDIH(I,-1:8)
            print*, '---'
         ENDDO
         print*, 'IMPROPER DIHEDRALS', NIMP
         DO I=1,NIMP
            print*, IMPDIH(I,:)
         ENDDO
      ENDIF

      IF ((BBCART.AND.NINTC.NE.3*NATOMS).OR.                            &
     &          (.NOT.BBCART.AND.NINTC.NE.3*NATOMS-6)) THEN
         print*, 'ERROR!'
         print*, 'number of natural internals should be 3*natoms-6,     &
     &        or 3*natoms if bbcart', BBCART, NINTC, NATOMS, 3*NATOMS-6
      ENDIF

      RETURN
      END SUBROUTINE


!msb50 - for AMBER
      SUBROUTINE AMB_GETBACKBONE
       ! which of the atoms are part of the backbone?
      USE MODAMBER9 
      USE COMMONS, ONLY: NATOMS
      USE INTCOMMONS, ONLY: INTMINPERMT
 
         IMPLICIT NONE
 
         INTEGER :: A, R
       
         IF(.NOT.ALLOCATED(AM_BACKBONE)) ALLOCATE(AM_BACKBONE(NATOMS))
         AM_BACKBONE(1:NATOMS) = .FALSE.
       DO A = 1,NATOMS
          IF (ih(m04+A-1).EQ.'C   '.OR.ih(m04+A-1).EQ.'CA  '.OR.    &
     &       ih(m04+A-1).EQ.'N   '.OR.ih(m04+A-1).EQ.'CH3 '.OR. &
     &       ih(m04+A-1).EQ."C5' ".OR.ih(m04+A-1).EQ."C3' ".OR. &
     &       ih(m04+A-1).EQ."O3' ".OR.ih(m04+A-1).EQ."P") AM_BACKBONE(A) = .TRUE.
             !last case - terminal - need to double check this
       ENDDO

       IF ((ih(m04+A-1).EQ."C5' ".OR.ih(m04+A-1).EQ."C3' ".OR. &
     &       ih(m04+A-1).EQ."O3' ".OR.ih(m04+A-1).EQ."P").AND. &
     &       INTMINPERMT) THEN
           PRINT*, "*************ERROR*******************"
           PRINT*,"INTMINPERM does NOT WORK with RNA BACKBONE settings!"
           PRINT*,"use NOPERMPROCHIRAL in addition!!"
        ENDIF 

       ! label all proline atoms (except O) as backbone
       DO R = 1,NRES
          IF ((ih(m02+R-1)).EQ.'PRO'.OR. (ih(m02+R-1)).EQ.'NPRO'.OR.     &
     &        (ih(m02+R-1)).EQ.'CPRO') THEN
             PRINT*, "YAY PROLINE"
             DO A = ix(i02+R-1),ix(i02+R)-1
                IF(ih(m04+A-1)(1:1).NE.'O') AM_BACKBONE(A) = .TRUE.
             ENDDO
          ENDIF
       ENDDO
 
       ! label all hydroxyproline atoms (except O) as backbone
       DO R = 1,NRES
          IF ((ih(m02+R-1)).EQ.'HYP'.OR. (ih(m02+R-1)).EQ.'NHYP'.OR.     &
     &        (ih(m02+R-1)).EQ.'CHYP') THEN
             PRINT*, "YAY HYDROXYPROLINE"
             DO A = ix(i02+R-1),ix(i02+R)-1
                IF(ih(m04+A-1)(1:1).NE.'O') AM_BACKBONE(A) = .TRUE.
             ENDDO
          ENDIF
       ENDDO
 

      
       RETURN
 
      END SUBROUTINE AMB_GETBACKBONE

! ***********************************************************************
!yay for amber - again - see efk
      SUBROUTINE AMB_NATINTSETUP
        ! setup up various arrays necessary for natural internals
        ! this is the routine that defines the natural internal coordinates in terms of angle and torsion coefficients
         USE MODAMBER9
         USE INTCOMMONS
         USE COMMONS, ONLY:NATOMS   
  
         IMPLICIT NONE
         INTEGER, PARAMETER :: LARGENEG = -99999999
         DOUBLE PRECISION :: S6R, S2R, S26R, S18R, S12R, A, B, SABR, SABR2
         
         NBDS = 0; NCNT = 0; NRNG = 0; NFRG = 0; NLDH = 0; NIMP = 0
         NCRT = 0; NDIH = 0
         
         !msb50 RINGS(3*NRES,6) as can have ribose and 2 joined rings per residue
         ALLOCATE(CENTERS(NATOMS,0:5), RINGS(3*NRES,6), FRINGS(NRES,6))
         ALLOCATE(LINDIH(NATOMS,-1:8), IMPDIH(NATOMS,4), CARTATMS(NATOMS))
         
         IF (CARTRESSTART.EQ.0) THEN
             PRINT*, "amb_natintsetup> start from residue 1, not 0!"
             CARTRESSTART=1
         ENDIF

         IF (CARTRESSTART.GT.NRES+1.OR.CARTRESSTART.LT.0) THEN
            CARTRESSTART = NRES + 1
            CARTATMSTART = NATOMS + 1
         ELSE
            CARTATMSTART = ix(i02+CARTRESSTART-1)
         ENDIF
         
         !     Set up COEFF table
         S6R = 1/SQRT(6.0)
         S2R = 1/SQRT(2.0)
         S26R = 1/SQRT(26.0)
         S18R = 1/SQRT(18.0)
         S12R = 1/SQRT(12.0)

         COEFF(1,1,1:6) = (/-1,-1,-1,1,1,1/)*S6R
         COEFF(1,2,1:6) = (/0,0,0,2,-1,-1/)*S6R
         COEFF(1,3,1:6) = (/0,0,0,0,1,-1/)*S2R
         COEFF(1,4,1:6) = (/2,-1,-1,0,0,0/)*S6R
         COEFF(1,5,1:6) = (/0,1,-1,0,0,0/)*S2R     
       
         COEFF(2,1,1:6) = (/0,1,-1,1,-1,0/)*0.5
         COEFF(2,2,1:6) = (/0,-1,-1,1,1,0/)*0.5
         COEFF(2,3,1:6) = (/0,-1,1,1,-1,0/)*0.5
         COEFF(2,4,1:6) = (/1,0,0,0,0,5/)*S26R
         COEFF(2,5,1:6) = (/5,0,0,0,0,1/)*S26R
       
         COEFF(3,1,1:6) = (/0,0,2,0,-1,-1/)*S6R
         COEFF(3,2,1:6) = (/0,0,0,0,1,-1/)*S2R
         COEFF(3,3,1:6) = (/1,4,0,1,0,0/)*S18R
         COEFF(3,4,1:6) = (/1,1,0,4,0,0/)*S18R
         COEFF(3,5,1:6) = (/4,1,0,1,0,0/)*S18R
       
         COEFF(4,1,1:3)=(/-1,-1,2/)*S6R
         COEFF(4,2,1:3)=(/1,-1,0/)*S2R
       
         COEFF(5,1,1:3)=(/0,1,-1/)*S2R
         COEFF(5,2,1:3)=(/2,-1,-1/)*S6R
       
         COEFF(6,1,1:1) = (/1/)
       
         COEFF(7,1,1:2) = (/1,-1/)*S2R
       
         COEFF(8,1,1:5) = (/1,-1,1,-1,0/)*0.5
         COEFF(8,2,1:5) = (/-1,-1,1,1,0/)*0.5
         COEFF(8,3,1:5) = (/-1,1,1,-1,0/)*0.5
         COEFF(8,4,1:5) = (/0,0,0,0,1/)
       
         A = COS(144.0*PI/180.0)
         B = COS(72.0*PI/180.0)
         SABR = 1/SQRT(1+2*A*A+2*B*B)
         SABR2 = 1/SQRT(2*(A-B)**2+2*(1-A)**2)
         COEFF(25,1,1:5) = (/1.0D0,A,B,B,A/)*SABR
         COEFF(25,2,1:5) = (/0.0D0,A-B,1.0D0-A,A-1.0D0,B-A/)*SABR2
         COEFF(25,3,1:5) = (/B,A,1.0D0,A,B/)*SABR
         COEFF(25,4,1:5) = (/A-1.0D0,B-A,0.0D0,A-B,1.0D0-A/)*SABR2
       
         COEFF(26,1,1:6) = (/1,-1,1,-1,1,-1/)*S6R
         COEFF(26,2,1:6) = (/2,-1,-1,2,-1,-1/)*S6R
         COEFF(26,3,1:6) = (/0,1,-1,0,1,-1/)*0.5D0
         COEFF(26,4,1:6) = (/1,-1,1,-1,1,-1/)*S6R
         COEFF(26,5,1:6) = (/1,0,-1,1,0,-1/)*0.5D0
         COEFF(26,6,1:6) = (/-1,2,-1,-1,2,-1/)*S12R
       
         !     Set up ATOMxPOSINC table
         !     large negative number listed for those atoms not involved
         ATOMxPOSinC(1,1,1:5) = (/0,3,6,9,12/)
         ATOMxPOSinC(1,2,1:5) = (/15,LARGENEG,18,21,24/)
         ATOMxPOSinC(1,3,1:5) = (/27,LARGENEG,30,33,36/)
         ATOMxPOSinC(1,4,1:5) = (/39,42,45,48,51/)
         ATOMxPOSinC(1,5,1:5) = (/54,57,LARGENEG,60,63/)

         ATOMxPOSinC(2,1,1:5) = (/0,3,6,9,12/)
         ATOMxPOSinC(2,2,1:5) = (/15,18,21,24,27/)
         ATOMxPOSinC(2,3,1:5) = (/30,33,36,39,42/)
         ATOMxPOSinC(2,4,1:5) = (/45,48,51,54,57/)
         ATOMxPOSinC(2,5,1:5) = (/60,63,66,69,72/)
       
         ATOMxPOSinC(3,1,1:5) = (/0,3,6,9,12/)
         ATOMxPOSinC(3,2,1:5) = (/15,LARGENEG,18,21,24/)
         ATOMxPOSinC(3,3,1:5) = (/27,30,33,36,LARGENEG/)
         ATOMxPOSinC(3,4,1:5) = (/39,42,45,48,LARGENEG/)
         ATOmxPOSinC(3,5,1:5) = (/51,54,57,60,LARGENEG/)
       
         ATOMxPOSinC(4,1,1:4)=(/0,3,6,9/)
         ATOMxPOSinC(4,2,1:4)=(/12,15,18,21/)
       
         ATOMxPOSinC(5,1,1:4)=(/0,3,6,9/)
         ATOMxPOSinC(5,2,1:4)=(/12,15,18,21/)
       
         ATOMxPOSinC(6,1,1:3) = (/0,3,6/)
       
         ATOMxPOSinC(7,1,1:4) = (/0,3,6,9/)      
      
         ATOMxPOSinC(8,1,1:5)=(/0,3,6,9,12/)
         ATOMxPOSinC(8,2,1:5)=(/15,18,21,24,27/)
         ATOMxPOSinC(8,3,1:5)=(/30,33,36,39,42/)
         ATOMxPOSinC(8,4,1:5)=(/45, LARGENEG, LARGENEG,48,51/)
       END SUBROUTINE AMB_NATINTSETUP
       


! ************************************************************************************       
      SUBROUTINE AMB_GETBEENAT(XCART,XINT,VAL,INDX,JNDX,NINTC,NCART,NNZ,GCS,NOCOOR,NODERV,KD)
      USE COMMONS, ONLY : NATOMS
      USE MODAMBER9, ONLY: IBIB,JBJB, NBONA,NBONH,ih,m04
      USE INTCOMMONS, ONLY : CARTATMSTART, PRINTCOORDS, COEFF, CENTERS, RINGS,&
     &     FRINGS, LINDIH, IMPDIH, CARTATMS, ATOMxPOSINC, NBDS, NCNT, NRNG,&
     &     NFRG, NLDH, NIMP, NCRT, USECART, PREVDIH
      
      IMPLICIT NONE
! Get internal coordinates and derivative matrix for natural internal coords


!-----------------------------------------------------------------------
!
! NNZ is number of non-zero elements in B-matrix
! NINTC is number of internal coordinates
! NCART is number of cartesians
!
       INTEGER NNZ,NINTC,NCART
!
       INTEGER INDX(NNZ),JNDX(NNZ)
       INTEGER MMM,ITHH,IPHII,II,JJ,KK,LL,IIKD,JJKD,KKKD,LLKD,KD1,KD
       DOUBLE PRECISION XCART(NCART),XINT(NINTC)
       DOUBLE PRECISION X(NCART),Y(NCART),Z(NCART)
       DOUBLE PRECISION VAL(NNZ),GCS(2*KD+1,NCART)
       DOUBLE PRECISION RX,RY,RZ,S2,S,SR,RXSR,RYSR,RZSR
       DOUBLE PRECISION DXI,DYI,DZI,DXJ,DYJ,DZJ,RI2,RJ2,RIR,RJR,RI,RJ
       DOUBLE PRECISION DXIR,DYIR,DZIR,DXJR,DYJR,DZJR,CST,ISNT
       DOUBLE PRECISION FX,FY,FZ,GX,GY,GZ,HX,HY,HZ,AX,AY,AZ,BX,BY,BZ
       DOUBLE PRECISION CSTTWO,SNTTWO2,SNTTWO2R,CSTTHREE,SNTTHREE2,SNTTHREE2R,DUMMY,DUMMY2
       DOUBLE PRECISION VXI,VYI,VZI,VXJ,VYJ,VZJ,VXK,VYK,VZK,VXL,VYL,VZL
       INTEGER MM,ITH,IPHI,I,J,K,L,IC,I1,J1,SUMPHI,DIFF,ATOMWIDTH
       INTEGER II1,JJ1,KK1,LL1,IIKD1,JJKD1,KKKD1,LLKD1
       INTEGER II2,JJ2,KK2,LL2,IIKD2,JJKD2,KKKD2,LLKD2
       LOGICAL NOCOOR, NODERV
!jmc
       DOUBLE PRECISION MYSCALAR,MYTX,MYTY,MYTZ

      INTEGER CURCOORD, CURBind, DIHC
      INTEGER XIND(8)
      INTEGER Ix, Jx, Kx, Lx
      INTEGER K1, L1, X1, X2, A1, C1, TCOUNT, NANGC, NDIHC
      INTEGER CNT, LDH, IMP, RNG, FRG, CNTATMS, RATMS
      DOUBLE PRECISION CURCOEFF
      INTEGER INFO(-1:8), CNTTYPE, RTYPE
      INTEGER ATOMxSHIFT, ISHIFT, JSHIFT
      DOUBLE PRECISION PHI, THETA, dTdI(0:2), dTdJ(0:2), dTdK(0:2)
      DOUBLE PRECISION ANGLE, dI(0:2), dJ(0:2), dK(0:2), dL(0:2)
      DOUBLE PRECISION dPdI(0:2), dPdJ(0:2), dPdK(0:2), dPdL(0:2)
      DOUBLE PRECISION FNANGC, SNANGC, SNANGCR, NANGCR
      DOUBLE PRECISION INDLDH(NLDH, 9)
      INTEGER CURPHI
      INTEGER STARTLINDH
      DOUBLE PRECISION PHI2, dPdI2(0:2), dPdJ2(0:2), dPdK2(0:2), dPdL2(0:2)
      DOUBLE PRECISION TEST(3*NATOMS,3*NATOMS)
 
      DOUBLE PRECISION PHI1, PREVPHI1

      INTEGER NBOND

      TEST(:,:) = 0.0

      NBOND = nbona +nbonh
      IF (.NOT.ALLOCATED(IBIB).OR..NOT.ALLOCATED(JBJB)) THEN
      PRINT*, "AMBERT: BONDS NOT ALLOCATED. CHECK AMBGETNATINTERN. STOP"
      STOP
      ENDIF

      IF (NCART.NE.3*NATOMS) THEN
         print*, 'getbeenat>> ERROR! NCART not 3*NATOMS', NCART, NATOMS
         STOP
      ENDIF
!---------------------------------------------------

!Zero the internal coord vector if not nocoord
       IF (.NOT.NOCOOR) XINT(1:NINTC) = 0.0D0

!Zero the B matrix values
       IF (.NOT.NODERV) THEN
          VAL(1:NNZ) = 0.0D0
          INDX(1:NNZ) = 0
          JNDX(1:NNZ) = 0
       ENDIF
!put cartesians from XCART into X,Y,Z
!
      DO I1=1,NATOMS
        X(I1)=XCART(3*(I1-1)+1)
        Y(I1)=XCART(3*(I1-1)+2)
        Z(I1)=XCART(3*(I1-1)+3)
      ENDDO
!
!KD1 is dimension of GCS sparse GC matrix
!in fact using matrix double this size so can calculate upper and lower triangle separately
!will combine the two at the end

      KD1=KD+1
!
!zero GCS here
!
      IF (.NOT.NODERV) GCS(:,:) = 0.0D0

!
!the following loops will put internals in XINT and pass them back
!if the flag NOCOOR is false
!and put derivatives in VAL if NODERV is false

      CURBIND = 1 ! current index in B-matrix lists
      DIHC = 0

      MM = 0
      DO 10 TCOUNT=1,NBOND                 
        I=IBIB(TCOUNT)
        J=JBJB(TCOUNT)
        ! ignore bond if using cartesians for both atoms
        IF ((USECART(I).AND.USECART(J)) ) CYCLE

        MM = MM + 1
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
 
!msb50 - no idea what the next three lines would correspond to in amber, or why they are there... please fix this if you can.
!      IC=ICB(TCOUNT)
!      IF(IC.EQ.0) GOTO 10
!      IF(CBC(IC).EQ.0) GOTO 10

        RX=X(I)-X(J)
        RY=Y(I)-Y(J)
        RZ=Z(I)-Z(J)
        S2=RX*RX + RY*RY + RZ*RZ
        S=SQRT(S2)
        IF (.NOT.NOCOOR) THEN
!dae put S into XINT: it is internal coordinate
          XINT(MM)=S
        ENDIF
        IF (PRINTCOORDS) print '(a14,i4,2(i4,a4),f7.3)', 'GETBEE> bond', MM, I, ih(m04+I-1),J,ih(m04+J-1), XINT(MM)

        IF (.NOT.NODERV) THEN
!dae calculate B-matrix elements for this bond
!dae VAL is 1-D array of non-zero B elements
!INDX is array of the ic(row) indices of VAL
!JNDX is array of the cart(column) indices of VAL
        SR=1/S
        RXSR=RX*SR
        RYSR=RY*SR
        RZSR=RZ*SR

        VAL(CURBIND)=RXSR
        INDX(CURBIND)=MM
        JNDX(CURBIND)=3*(I-1)+1

        TEST(MM, 3*(I-1)+1) = RXSR
!
        VAL(CURBIND+1)=RYSR
        INDX(CURBIND+1)=MM
        JNDX(CURBIND+1)=3*(I-1)+2

        TEST(MM,3*(I-1)+2) = RYSR
!
        VAL(CURBIND+2)=RZSR
        INDX(CURBIND+2)=MM
        JNDX(CURBIND+2)=3*(I-1)+3

        TEST(MM,3*(I-1)+3) = RZSR
!
        VAL(CURBIND+3)=-RXSR
        INDX(CURBIND+3)=MM
        JNDX(CURBIND+3)=3*(J-1)+1

        TEST(MM,3*(J-1)+1) = -RXSR
!
        VAL(CURBIND+4)=-RYSR
        INDX(CURBIND+4)=MM
        JNDX(CURBIND+4)=3*(J-1)+2

        TEST(MM,3*(J-1)+2) = -RYSR
!
        VAL(CURBIND+5)=-RZSR
        INDX(CURBIND+5)=MM
        JNDX(CURBIND+5)=3*(J-1)+3

        TEST(MM,3*(J-1)+3) = -RZSR
!
!update GEE for this ic
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1

        VXI=VAL(CURBIND)
        VYI=VAL(CURBIND+1)
        VZI=VAL(CURBIND+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(IIKD-II1,II1)=GCS(IIKD-II1,II1)+VXI*VYI
        GCS(IIKD-II2,II2)=GCS(IIKD-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(IIKD1-II2,II2)=GCS(IIKD1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
!
        VXJ=VAL(CURBIND+3)
        VYJ=VAL(CURBIND+4)
        VZJ=VAL(CURBIND+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(JJKD-JJ1,JJ1)=GCS(JJKD-JJ1,JJ1)+VXJ*VYJ
        GCS(JJKD-JJ2,JJ2)=GCS(JJKD-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(JJKD1-JJ2,JJ2)=GCS(JJKD1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
!
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
!
        CURBIND = CURBIND + 6
        ENDIF
   10  CONTINUE

      CURCOORD = NBDS+1
!     print*, 'AFTER BONDS, CURBIND=', CURBIND

      DO 20 CNT=1,NCNT

        INFO(0:5) = CENTERS(CNT,0:5)
        CNTTYPE = INFO(0)

        IF (CNTTYPE.LE.3) THEN
           CNTATMS = 4 ! number of atoms bound to center
           NANGC = 5 ! number of angle-sum natural coords
        ELSE IF (CNTTYPE.EQ.6) THEN
           CNTATMS = 2
           NANGC = 1
        ELSE IF (CNTTYPE.EQ.7) THEN
           CNTATMS = 3
           NANGC = 1
       ! ELSE IF (CNTTYPE.EQ.8) THEN
       !    CNTATMS = 4
       !    NANGC = 4
       !    CNTTYPE = 3
       ! ELSE IF (CNTTYPE.EQ.9) THEN
       !    CNTATMS = 4
       !    NANGC = 4
       !    CNTTYPE = 2
        ELSE IF (CNTTYPE.EQ.8) THEN
           CNTATMS = 4 
           NANGC =4 
        ELSE IF (CNTTYPE.EQ.4.OR.CNTTYPE.EQ.5) THEN
           CNTATMS = 3
           NANGC = 2
        ENDIF

!Get X coordinate indices in XCART for each atom involved
        DO I1=1,CNTATMS+1
           XIND(I1) = 3*(INFO(I1) - 1) + 1
        ENDDO

!dae calculate B-matrix elements for these atoms (Wilson Decius + Cross ch.4)
!EFK: add together the B-matrix elements from the relevant primitive internals        
!get angles and derivatives for all angles with the given central atom
        TCOUNT=1
        DO I1=2,CNTATMS
           DO 110 K1=I1+1,CNTATMS+1
              ! if center is part of ring, skip angle within ring
              IF((CNTTYPE.EQ.7.OR.CNTTYPE.EQ.8).AND.I1.EQ.2.AND.K1.EQ.3) GOTO 110
!
              CALL AMBGETANGLE(XCART(XIND(I1):XIND(I1)+2),                 &
     &             XCART(XIND(1):XIND(1)+2), XCART(XIND(K1):XIND(K1)+2), &
     &             THETA, dTdI, dTdJ, dTdK, NOCOOR,NODERV)              
              IF (PRINTCOORDS) print '(a15,3(i4,a4),f7.2)','GETBEE> &
          &angle', INFO(I1), ih(m04+ INFO(I1)-1), INFO(1), ih(m04+INFO(1)-1),&
     &      INFO(K1),ih(m04+INFO(K1)-1), THETA

               DO 100 C1=1,NANGC ! counts coordinates for this center
                 CURCOEFF = COEFF(CNTTYPE, C1, TCOUNT)

                 IF (CURCOEFF.EQ.0) GOTO 100

                 ! update internal coords list
                 IF (.NOT.NOCOOR) THEN
                    XINT(CURCOORD+C1-1) = XINT(CURCOORD+C1-1)           & 
     &                   + CURCOEFF * THETA
                 ENDIF

                 IF (.NOT.NODERV) THEN
                    ! update derivatives in B matrix
                    DO X1=0,2   ! counts x,y,z coordinates
                       VAL(CURBIND + ATOMxPOSINC(CNTTYPE,C1,I1) + X1) = &
     &                      VAL(CURBIND + ATOMxPOSINC(CNTTYPE,C1,I1)+X1)&
     &                      + CURCOEFF * dTdI(X1)
                       VAL(CURBIND + ATOMxPOSINC(CNTTYPE,C1,1) + X1) =  &
     &                      VAL(CURBIND + ATOMxPOSINC(CNTTYPE,C1,1)+X1) &
     &                      + CURCOEFF * dTdJ(X1)
                       VAL(CURBIND + ATOMxPOSINC(CNTTYPE,C1,K1) + X1) = &
     &                      VAL(CURBIND + ATOMxPOSINC(CNTTYPE,C1,K1)+X1)&
     &                      + CURCOEFF * dTdK(X1)

                    ENDDO
                 ENDIF
 100          CONTINUE
              TCOUNT = TCOUNT + 1
 110       CONTINUE
        ENDDO

        IF (.NOT.NODERV) THEN
!    Set B matrix row and column values
           DO C1=1,NANGC
              DO A1 = 1,CNTATMS+1
                 ATOMxSHIFT = ATOMxPOSINC(CNTTYPE,C1,A1)
                 IF (.NOT.(ATOMxSHIFT.LT.0)) THEN                 
                    DO X1=0,2
                       INDX(CURBIND + ATOMxSHIFT+X1) = CURCOORD+C1-1
                       JNDX(CURBIND + ATOMxSHIFT+X1) = XIND(A1)+X1
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
       
 
!    Update GEE matrix elements
           DO I1=1,CNTATMS+1
              DO J1=I1,CNTATMS+1
                 DO C1=1,NANGC
                    IF(ATOMxPOSinC(CNTTYPE,C1,I1).GE.0 .AND.            & 
     &                   ATOMxPOSinC(CNTTYPE,C1,J1).GE.0) THEN
                                ! both atoms are involved in this coordinate
                       IX = XIND(I1) ! index of x coordinate of atom i
                       JX = XIND(J1)
                                ! find how far offset from current b index are
                                ! the derivatives wrt each atom
                       ISHIFT = ATOMxPOSINC(CNTTYPE,C1,I1) 
                       JSHIFT = ATOMxPOSINC(CNTTYPE,C1,J1)
                       
                       IF (I1.EQ.J1) THEN
                          DO X1=0,2
                             DO X2=X1,2
                                GCS(KD1+IX+X1-(JX+X2), JX+X2) =         &
     &                               GCS(KD1+IX+X1-(JX+X2),JX+X2) +     &
     &                               VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                             ENDDO
                          ENDDO
                       ELSE
                          DO X1=0,2
                             DO X2=0,2
                                GCS(KD1+IX+X1-(JX+X2), JX+X2) =         &
     &                               GCS(KD1+IX+X1-(JX+X2),JX+X2) +     &
     &                               VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                             ENDDO
                          ENDDO
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO

        IF(CNTTYPE.NE.3) THEN
           CURBIND = CURBIND + ATOMxPOSinC(CNTTYPE, NANGC, CNTATMS+1)+3
        ELSE
           CURBIND = CURBIND + ATOMxPOSinC(CNTTYPE, NANGC, CNTATMS)+3
        ENDIF
      ENDIF

        CURCOORD = CURCOORD + NANGC
 20   CONTINUE

!    print*, 'AFTER CNT CURBIND=', CURBIND

!!
!Deal with rings
!!
      DO 25 RNG=1,NRNG

         IF(RINGS(RNG,6).EQ.-1) THEN
            RATMS = 5 ! 5 atoms in ring
            NANGC = 2; NDIHC = 2 ! 2 angle & 2 dihedral coords
            RTYPE = 25
         ELSE 
            RATMS = 6 ! 6 atoms in ring
            NANGC = 3; NDIHC = 3            
            RTYPE = 26
         ENDIF
         
         DO I1=1,RATMS
            XIND(I1) = 3*(RINGS(RNG,I1) - 1) + 1
         ENDDO
        
         TCOUNT = 1
         DO I1 = 1,RATMS
            J1 = MOD(I1,RATMS)+1
            K1 = MOD(I1+1,RATMS)+1
            L1 = MOD(I1+2,RATMS)+1

            CALL AMBGETANGLE(XCART(XIND(I1):XIND(I1)+2),                   &
     &           XCART(XIND(J1):XIND(J1)+2), XCART(XIND(K1):XIND(K1)+2),&
     &           THETA, dTdI, dTdJ, dTdK, NOCOOR,NODERV)              

            CALL AMBGETTORSION(XCART(XIND(I1):XIND(I1)+2),                 &
     &           XCART(XIND(J1):XIND(J1)+2), XCART(XIND(K1):XIND(K1)+2),&
     &           XCART(XIND(L1):XIND(L1)+2), PHI, dPdI, dPdJ, dPdK, dPdL,&
     &           NOCOOR,NODERV)


            DIHC = DIHC + 1

            IF (TCOUNT.EQ.1) PREVPHI1 = PREVDIH(DIHC)
          
 
            IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,TCOUNT)
            IF (TCOUNT.EQ.1) PHI1 = PHI

            DO C1 = 1,NANGC+NDIHC
               IF (C1.LE.NANGC) THEN
                  ANGLE = THETA; dI(:) = dTdI(:)
                  dJ(:)=dTdJ(:); dK(:)=dTdK(:)
               ELSE
                  ANGLE = PHI
                  dI(:)=dPdI(:); dJ(:)=dPdJ(:)
                  dK(:)=dPdK(:); dL(:)=dPdL(:)
               ENDIF

               CURCOEFF = COEFF(RTYPE, C1, TCOUNT)
               
               IF(.NOT.NOCOOR) XINT(CURCOORD+C1-1) = &
     &              XINT(CURCOORD+C1-1)+CURCOEFF*ANGLE

               IF (.NOT.NODERV) THEN
!    Update B-matrix elements
               DO X1=0,2        ! counts x,y,z coordinates
                  VAL(CURBIND + RATMS*3*(C1-1) + 3*(I1-1) + X1) =&
     &                 VAL(CURBIND + RATMS*3*(C1-1)+ 3*(I1-1) + X1)&
     &                 + CURCOEFF * dI(X1)
                  VAL(CURBIND + RATMS*3*(C1-1) + 3*(J1-1) + X1) =&
     &                 VAL(CURBIND + RATMS*3*(C1-1) + 3*(J1-1) + X1)&
     &                 + CURCOEFF * dJ(X1)
                  VAL(CURBIND + RATMS*3*(C1-1) + 3*(K1-1) + X1) =&
     &                 VAL(CURBIND + RATMS*3*(C1-1)+ 3*(K1-1) + X1)&
     &                 + CURCOEFF * dK(X1)                  
               ENDDO
               IF (C1.GT.NANGC) THEN
                  DO X1=0,2
                     VAL(CURBIND + RATMS*3*(C1-1) + 3*(L1-1) + X1) =&
     &                    VAL(CURBIND + RATMS*3*(C1-1)+ 3*(L1-1) + X1)&
     &                    + CURCOEFF * dL(X1)
                  ENDDO
               ENDIF
               ENDIF
            ENDDO
         TCOUNT = TCOUNT + 1
         ENDDO

!Set B matrix row and column values
         IF (.NOT.NODERV) THEN
         DO C1=0,NANGC+NDIHC-1  ! count coords
            DO A1 = 0,RATMS-1   ! count atoms
               DO X1=0,2        ! count x,y,z
                  INDX(CURBIND+RATMS*3*C1+3*A1+X1) = CURCOORD+C1
                  JNDX(CURBIND+RATMS*3*C1+3*A1+X1) = XIND(A1+1)+X1
               ENDDO
            ENDDO
         ENDDO
         
!Update G matrix elements
         DO I1=1,RATMS
            DO J1=I1,RATMS
               DO C1=1,NANGC+NDIHC

                  IX = XIND(I1)
                  JX = XIND(J1)
                  ISHIFT = RATMS*3*(C1-1)+3*(I1-1)
                  JSHIFT = RATMS*3*(C1-1)+3*(J1-1)

                  IF (I1.EQ.J1) THEN
                     DO X1=0,2
                        DO X2 = X1,2
                           GCS(KD1+IX+X1-(JX+X2), JX+X2) = &
     &                          GCS(KD1+IX+X1-(JX+X2), JX+X2) &
     &                          + VAL(CURBIND+ISHIFT+X1) * VAL(CURBIND+JSHIFT+X2)
                        ENDDO
                     ENDDO
                  ELSE
                     DO X1=0,2
                        DO X2 = 0,2
                           GCS(KD1+IX+X1-(JX+X2), JX+X2) = &
     &                          GCS(KD1+IX+X1-(JX+X2), JX+X2) &
     &                          + VAL(CURBIND+ISHIFT+X1) * VAL(CURBIND+JSHIFT+X2)                           
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CURBIND = CURBIND + RATMS*3*(NANGC + NDIHC)
         ENDIF
         CURCOORD = CURCOORD + NANGC + NDIHC

 25   CONTINUE
!     print*, 'AFTER RNG, CURBIND:', CURBIND

!
!Deal with fused rings
!Get difference of dihedrals 3-1-2-6 & 4-2-1-5
!
      DO FRG = 1,NFRG
         DO I1=1,6
            XIND(I1) = 3*(FRINGS(FRG,I1) - 1)+1
         ENDDO

         CALL AMBGETTORSION(XCART(XIND(3):XIND(3)+2),&
     &           XCART(XIND(1):XIND(1)+2), XCART(XIND(2):XIND(2)+2), &
     &           XCART(XIND(6):XIND(6)+2), PHI, dPdI, dPdJ, dPdK, dPdL, &
     &           NOCOOR,NODERV)   
         DIHC = DIHC + 1
         PREVPHI1 = PREVDIH(DIHC)
         IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, 0.0D0, 0.0D0, 1)         

         CALL AMBGETTORSION(XCART(XIND(5):XIND(5)+2),&
     &           XCART(XIND(1):XIND(1)+2), XCART(XIND(2):XIND(2)+2),   &
     &           XCART(XIND(4):XIND(4)+2), PHI2, dPdI2, dPdJ2, dPdK2, dPdL2,&
     &           NOCOOR,NODERV)
         DIHC = DIHC + 1
         
         IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI2, DIHC, PHI2-PHI, PREVDIH(DIHC)-PREVPHI1, 2)         
         
         IF (.NOT.NOCOOR) XINT(CURCOORD) = XINT(CURCOORD) + (PHI - PHI2)/SQRT(2.0D0)


         IF (.NOT.NODERV) THEN

!Update B matrix
         DO X1 = 0,2
            VAL(CURBIND+X1) = VAL(CURBIND+X1) + (dPdJ(X1) - dPdJ2(X1))/SQRT(2.0D0)
            VAL(CURBIND+X1+3) = VAL(CURBIND+X1+3) + (dPdK(X1)-dPdK2(X1))/SQRT(2.0D0)
            VAL(CURBIND+X1+6) = VAL(CURBIND+X1+6) + dPdI(X1)/SQRT(2.0D0)
            VAL(CURBIND+X1+9) = VAL(CURBIND+X1+9) - dPdL2(X1)/SQRT(2.0D0)
            VAL(CURBIND+X1+12) = VAL(CURBIND+X1+12) - dPdI2(X1)/SQRT(2.0D0)
            VAL(CURBIND+X1+15) = VAL(CURBIND+X1+15) + dPdL(X1)/SQRT(2.0D0)
            DO A1=1,6
               INDX(CURBIND+3*(A1-1)+X1) = CURCOORD
               JNDX(CURBIND+3*(A1-1)+X1) = XIND(A1)+X1
            ENDDO
         ENDDO
!Update GCS
         DO I1=1,6
            DO J1=I1,6
               IX = XIND(I1)
               JX = XIND(J1)

               ISHIFT = 3*(I1-1)
               JSHIFT = 3*(J1-1)

               IF (I1.EQ.J1) THEN
                  DO X1=0,2
                     DO X2=X1,2
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =        &
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2) +   &
     &                       VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                     ENDDO
                  ENDDO
               ELSE
                  DO X1=0,2
                     DO X2=0,2
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =&
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2) +&
     &                       VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CURBIND = CURBIND + 18
         ENDIF
         CURCOORD = CURCOORD + 1
      ENDDO
!     print*, 'AFTER FRG CURBIND:', CURBIND

!
!Deal with linear dihedrals
!
      STARTLINDH = DIHC +1 !line 1272 --> for AMBALIGN-CALL
      DO 30 LDH=1,NLDH
         INFO(-1:8) = LINDIH(LDH,-1:8)
         !INFO(-1) and INFO(0) denote the number of possible I and L atoms
         DO I1=1,INFO(0)+INFO(-1)+2
            XIND(I1) = 3*(INFO(I1) - 1) + 1
         ENDDO
         NANGC = INFO(-1)*INFO(0) ! number of torsions involved

         NANGCR = 1.0D0/NANGC

         CURPHI = 0         
         DO I1=3,INFO(-1)+2
            DO L1=INFO(-1)+3,INFO(-1)+INFO(0)+2
               CURPHI = CURPHI + 1
               
               Ix = XIND(I1)
               Jx = XIND(1)
               Kx = XIND(2)
               Lx = XIND(L1)

               CALL AMBGETTORSION(XCART(Ix:Ix+2), XCART(Jx:Jx+2),  &
     &              XCART(Kx:Kx+2), XCART(Lx:Lx+2), PHI, dPdI, dPdJ,   &
     &              dPdK, dPdL, NOCOOR,NODERV)                              

               DIHC = DIHC + 1
               IF (CURPHI.EQ.1) PREVPHI1 = PREVDIH(DIHC)
 
               IF (PRINTCOORDS) print '(a15,3i4,4(i4,a4),f8.2)','GETBEE>& 
            &dih',LDH, CURCOORD,DIHC, INFO(I1), ih(m04+INFO(I1)-1),INFO(1),ih(m04+INFO(1)-1),&
              INFO(2), ih(m04+INFO(2)-1),INFO(L1),ih(m04+INFO(L1)-1), PHI
               
               IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,CURPHI)

               IF (CURPHI.EQ.1) PHI1 = PHI                  

               IF (.NOT.NOCOOR) XINT(CURCOORD) = XINT(CURCOORD) + PHI*NANGCR
               IF (.NOT.NODERV) THEN
               ! update B-matrix elements with this torsion
               DO X1=0,2
                  VAL(CURBIND+3*(I1-1)+X1)=VAL(CURBIND+3*(I1-1)+X1) + dPdI(X1)*NANGCR
                  VAL(CURBIND+X1) = VAL(CURBIND+X1) + dPdJ(X1)*NANGCR
                  VAL(CURBIND+X1+3) = VAL(CURBIND+X1+3) + dPdK(X1)*NANGCR
                  VAL(CURBIND+3*(L1-1)+X1)=VAL(CURBIND+3*(L1-1)+X1) + dPdL(X1)*NANGCR
               ENDDO
               ENDIF
            ENDDO
         ENDDO
          
         IF (.NOT.NODERV) THEN
         ! set B-matrix row and column indices
         DO A1=1,INFO(0)+INFO(-1)+2
            DO X1=0,2
               INDX(CURBIND+3*(A1-1)+X1) = CURCOORD
               JNDX(CURBIND+3*(A1-1)+X1) = XIND(A1) + X1
            ENDDO
         ENDDO

         ! update the Gee matrix elements
         DO I1=1,INFO(0)+INFO(-1)+2
            DO J1=I1,INFO(0)+INFO(-1)+2
               Ix = XIND(I1)
               Jx = XIND(J1)
               IF(I1.EQ.J1) THEN
                  DO X1=0,2
                     DO X2=X1,2
                        GCS(KD1+Ix+X1-(Jx+X2), Jx+X2) = &
     &                       GCS(KD1+Ix+X1-(Jx+X2), Jx+X2) &
     &                      + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)
                     ENDDO
                  ENDDO
               ELSE
                  DO X1=0,2
                     DO X2=0,2
                        GCS(KD1+Ix+X1-(Jx+X2), Jx+X2) =&
     &                       GCS(KD1+Ix+X1-(Jx+X2), Jx+X2)&
     &                      + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CURBIND = CURBIND + 3*(INFO(0) + INFO(-1) + 2)
         ENDIF
         CURCOORD = CURCOORD + 1
         

   30 CONTINUE
      !print*, 'AFTER LDH CURBIND=', CURBIND

      DO 40 IMP=1,NIMP
!improper dihedral ('wagging') coordinates
         INFO(1:4) = IMPDIH(IMP,1:4)
         Ix = 3*(INFO(1)-1)+1
         Jx = 3*(INFO(2)-1)+1
         Kx = 3*(INFO(3)-1)+1
         Lx = 3*(INFO(4)-1)+1

         CALL AMBGETOUTOFPLANE(XCART(Ix:Ix+2), XCART(Jx:Jx+2),&
     &              XCART(Kx:Kx+2), XCART(Lx:Lx+2), PHI, dPdI, dPdJ,&
     &              dPdK, dPdL, NOCOOR,NODERV)


         IF(.NOT.NOCOOR) XINT(CURCOORD) = PHI

         IF (.NOT.NODERV) THEN
!Get B-matrix elements
         DO X1=0,2
            VAL(CURBIND+X1) = dPdI(X1)
            INDX(CURBIND+X1) = CURCOORD
            JNDX(CURBIND+X1) = Ix+X1

            VAL(CURBIND+X1+3) = dPdJ(X1)
            INDX(CURBIND+X1+3) = CURCOORD
            JNDX(CURBIND+X1+3) = Jx+X1
            
            VAL(CURBIND+X1+6) = dPdK(X1)
            INDX(CURBIND+X1+6) = CURCOORD
            JNDX(CURBIND+X1+6) = Kx+X1

            VAL(CURBIND+X1+9) = dPdL(X1)
            INDX(CURBIND+X1+9) = CURCOORD
            JNDX(CURBIND+X1+9) = Lx+X1

         ENDDO

!Update GEE matrix
         DO I1=1,4
            DO J1=I1,4
               Ix = 3*(INFO(I1)-1)+1
               Jx = 3*(INFO(J1)-1)+1
               IF (I1.EQ.J1) THEN
                  DO X1=0,2
                     DO X2=X1,2                     
                        GCS(KD1+Ix+X1-(Jx+X2), Jx+X2) =                 &    
     &                       GCS(KD1+Ix+X1-(Jx+X2), Jx+X2)              &
     &                       + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)
                     ENDDO
                  ENDDO
               ELSE
                  DO X1=0,2
                     DO X2=0,2                     
                        GCS(KD1+Ix+X1-(Jx+X2), Jx+X2) =                 &
     &                       GCS(KD1+Ix+X1-(Jx+X2), Jx+X2)              &
     &                       + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)                        
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CURBIND = CURBIND + 12
         ENDIF
         CURCOORD = CURCOORD + 1
         
   40 CONTINUE
!     print*, 'AFTER IMP CURBIND=', CURBIND

      DO 50 MM = 1,NCRT
         I = CARTATMS(MM)
         IF (.NOT.NOCOOR) XINT(CURCOORD:CURCOORD+2) = XCART(3*(I-1)+1:3*I)

         IF (.NOT.NODERV) THEN
            DO X1 = 0,2
               VAL(CURBIND+X1) = 1.0D0
               INDX(CURBIND+X1) = CURCOORD+X1
               JNDX(CURBIND+X1) = 3*(I-1)+1+X1
               GCS(KD1,3*(I-1)+1+X1) = GCS(KD1,3*(I-1)+1+X1) + 1.0D0
            ENDDO
            CURBIND = CURBIND + 3
         ENDIF

         CURCOORD = CURCOORD + 3         
 50   CONTINUE
!     print*, 'AFTER CART CURBIND=', CURBIND

      IF (.NOT.NODERV) THEN
!add lower diagonal elements to upper diagonal. Avoid doubling diagonal
!by looping only to J1-1
       DO J1=1,KD
         DO I1=1,J1-1
           GCS(KD+1+I1-J1,J1)=GCS(KD+1+I1-J1,J1)+GCS(KD+1+J1-I1,I1)
         ENDDO
       ENDDO
       DO J1=KD+1,NCART
         DO I1=J1-KD,J1-1
           GCS(KD+1+I1-J1,J1)=GCS(KD+1+I1-J1,J1)+GCS(KD+1+J1-I1,I1)
         ENDDO
       ENDDO
      ENDIF
      
      IF (PRINTCOORDS) THEN
         print*, 'done printing coords'
         PRINT*, XINT
!         STOP
      ENDIF

      RETURN
      END SUBROUTINE


