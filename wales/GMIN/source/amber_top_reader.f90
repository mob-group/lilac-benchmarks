SUBROUTINE TOPOLOGY_READER(NBOND)
USE COMMONS,ONLY : BONDS,MYUNIT
USE PORFUNCS 
   IMPLICIT NONE
   CHARACTER(100) ENTRY
   INTEGER :: MYUNIT2,GETUNIT
   INTEGER,INTENT(OUT) :: NBOND
   INTEGER :: J1,START_IND,END_IND,NBONDH,NBONDA,HENTRIES,J3,J4,J5,NDUMMY,INTDUM,J6
   INTEGER , PARAMETER :: NWORDS=20
   CHARACTER(25) :: ENTRIES(NWORDS)=''
   !INTEGER, DIMENSION(:,:), ALLOCATABLE :: BONDS

   MYUNIT2=GETUNIT()
   OPEN(MYUNIT2,FILE='coords.prmtop',STATUS='OLD')
reading:DO

98    READ(MYUNIT2,'(A)',END=99) ENTRY
      CALL READ_LINE(ENTRY,NWORDS,ENTRIES)      !get all words in line
      IF (ENTRIES(2).EQ.'POINTERS') THEN        !get number of bonds
         READ(MYUNIT2,*)                             !ignore format identifier after flag
         READ(MYUNIT2,'(A)',END=99) ENTRY
         CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
         READ(ENTRIES(3),'(I8)') NBONDH
         READ(ENTRIES(4),'(I8)') NBONDA
         NBOND = NBONDH + NBONDA
         WRITE(MYUNIT,'(A,I8)') 'readtopology> Number of bonds:',NBOND
         ALLOCATE(BONDS(NBOND,2))
      ENDIF
      IF (ENTRIES(2).EQ. 'BONDS_INC_HYDROGEN') THEN
         READ(MYUNIT2,*)                             !ignore format identifier after flag
         HENTRIES=(NBONDH*3)/10
         HENTRIES=HENTRIES+((NBONDH*3)-(HENTRIES*10)) !number of lines of entries
         NDUMMY=1
         J5=1
         DO J3=1,HENTRIES                             !go through all lines
            READ(MYUNIT2,'(A)',END=99) ENTRY               !read line
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
            J4=1
            DO WHILE(J4.LE.10)
               IF (NDUMMY.LE.NBONDH) THEN
                  IF (J5.EQ.1) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,1) = INTDUM/3+1        !atom1
                     J5=2
                  ELSE IF (J5.EQ.2) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,2) = INTDUM/3+1        !atom2
                     J5=3
                  ELSE
                     J5=1
                     NDUMMY=NDUMMY+1
                  ENDIF
               ELSE
                  GOTO 98
               ENDIF
               J4=J4+1
            ENDDO
         ENDDO
      ENDIF
      IF (ENTRIES(2).EQ. 'BONDS_WITHOUT_HYDROGEN') THEN
         READ(MYUNIT2,*)                             !ignore format identifier after flag
         HENTRIES=(NBONDA*3)/10
         HENTRIES=HENTRIES+((NBONDA*3)-(HENTRIES*10)) !number of lines of entries
         NDUMMY=NBONDH+1
         J5=1
         DO J3=1,HENTRIES                             !go through all lines
            READ(MYUNIT2,'(A)',END=99) ENTRY               !read line
            CALL READ_LINE(ENTRY,NWORDS,ENTRIES)
            J4=1
            DO WHILE(J4.LE.10)
               IF (NDUMMY.LE.(NBONDH+NBONDA)) THEN
                  IF (J5.EQ.1) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,1) = INTDUM/3+1
                     J5=2
                  ELSE IF (J5.EQ.2) THEN
                     READ(ENTRIES(J4),'(I8)') INTDUM
                     BONDS(NDUMMY,2) = INTDUM/3+1
                     J5=3
                  ELSE
                     J5=1
                     NDUMMY=NDUMMY+1
                  ENDIF
               ELSE
                  GOTO 98
               ENDIF
               J4=J4+1
            ENDDO
         ENDDO
      ENDIF

   ENDDO reading
99 CLOSE(MYUNIT2)
!  DO J6=1,NBOND
!     WRITE(MYUNIT,'(A,I8,A,I8)') 'readtopology> Bond between',BONDS(J6,1),' and',BONDS(J6,2)
!  ENDDO

END SUBROUTINE
  
SUBROUTINE READ_LINE(LINE,NWORDS,WORDSOUT)
      CHARACTER(*), INTENT(IN) :: LINE
      INTEGER, INTENT(IN) :: NWORDS
      CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
      INTEGER:: J1,START_IND,END_IND,J2
      CHARACTER(25) :: WORD
      START_IND=0
      END_IND=0
      J1=1
      J2=0
      DO WHILE(J1.LE.LEN(LINE))
          IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
             START_IND=J1
          ENDIF
          IF (START_IND.GT.0) THEN
             IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
             IF (J1.EQ.LEN(LINE)) END_IND=J1
             IF (END_IND.GT.0) THEN
                J2=J2+1
                WORD=LINE(START_IND:END_IND)
                WORDSOUT(J2)=TRIM(WORD)
                START_IND=0
                END_IND=0
             ENDIF
          ENDIF
          J1=J1+1
      ENDDO
END SUBROUTINE


