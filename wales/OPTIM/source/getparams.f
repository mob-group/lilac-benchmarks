C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      PROGRAM OPTIM3
      USE COMMONS
      USE PORFUNCS
      USE OPTIMHEADER
      USE KEY, ONLY: FILTHSTR,SEQ,NUMGLY,TARFL,CASTEPJOB,CP2KJOB,ONETEPJOB,QCHEMJOB,QCHEMJOBPARAMS,QCHEMSCALE,VASPJOB, 
     &               QCHEMES,QCHEMESNAO,QCHEMESNELEC,QCHEMESNMO, MOLPROJOB, MOLPROSCALE, MOLPROJOBPARAMS,
     &               REAXFFJOB, AMBER12T,OPEPT
      USE MODAMBER9, ONLY: AMBERSTR,AMBERSTR1,INPCRD,ATMASS1
      USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_SETUP
      USE OPEP_INTERFACE_MOD, ONLY: OPEP_GET_NATOMS

      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, J1
      LOGICAL VARIABLES,CASTEP,ONETEP,CP2K,DFTP,CPMD,END,CAT,SKIPBL,CLEAR
      LOGICAL ECHO,AMBERT,NABT,RINGPOLYMERT,USERPOTT,QCHEM,VASP, MOLPRO,REAXFF
      COMMON /BUFINF/ ITEM, NITEMS, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION DUMMY
      CHARACTER ZDUM*5
      INTEGER J2, NELEMENTS, LSYS, NTYPE(105), IOS, NARG, FILTH, FILTH2, NTYPES
      INTEGER, ALLOCATABLE :: ELEMENT_NUMS(:)
      CHARACTER FNAME*80, TSTRING*80
      CHARACTER(LEN=80) :: SYS
      CHARACTER WORD*16
      CHARACTER(LEN=10)  check
      CHARACTER(LEN=20) OTEMP, OSTRING, CSTRING, FIRSTAMBERSTR
      CHARACTER(LEN=21) DSTRING1, DSTRING2
      CHARACTER(LEN=80) ARGSTRING, MYLINE
      LOGICAL AMH
      INTEGER :: NRES,I_RES,NOGLY
      INTEGER LUNIT,GETUNIT, EPOCH

! TIME LIMITING BINARIES - uncomment the below
! csw34> Workshop time limiting
! find the number of seconds since 0000 on 1/1/1970 (unix epoch)
!      EPOCH=TIME()
! compare this to the number that will have passed at 1am on 5/11/15
!      IF (EPOCH.LT.1449277200) THEN
!         WRITE(*, '(A)') "NOTE: Binary only licensed for use at AlgoSB2015"
!      ELSE
! if after the workshop has ended, stop execution with error
!         WRITE(*, '(A)') "ERROR: license missing, contact dw34@cam.ac.uk. Stopping."
!         STOP   
!      END IF

      CALL PRINTHEADER
      CASTEP=.FALSE.
      QCHEM=.FALSE.
      MOLPRO=.FALSE.
      REAXFF=.FALSE.
      VASP=.FALSE.
      ONETEP=.FALSE.
      CP2K=.FALSE. 
      CPMD=.FALSE.
      VARIABLES=.FALSE.
      RINGPOLYMERT=.FALSE.
      AMBER12T=.FALSE.
      AMBERT=.FALSE.
      NABT=.FALSE.
      USERPOTT=.FALSE.

C
C  The function iargc returns the argument count
C  The statement call getarg( k , arg ) gets the  kth  command-
C  line argument and puts it into arg.
C  The 0th argument is the command name.
C
C  This provides a new way to filthify OPTIM.
C
      CALL iargc_subr(NARG) ! SAT: subroutine interface to iargc
      FILTH=0
      FILTH2=0
      IF (NARG.GT.0) THEN
         CALL GETARG(1,ARGSTRING)
C
C  Both methods to convert the character to an integer seem to work.
C
C  It would be better to use ARGSTRING as a string, so that non-numerical
C  extensions would work. This would require corresponding changes to Filthy_Phyllis.
C  We could also remove the odata.read file if Filthy_Phyllis were changed to work
C  with fork and wait!
C
C        FILTH2=ATOI(ARGSTRING)
         READ(ARGSTRING,*) FILTH2
         WRITE(*,'(A,I6,A)') ' Command line argument read, using extension number ',FILTH2,' for all I/O'
         FILTH=FILTH2
      ELSE
         ARGSTRING=''
      ENDIF

      IF (FILTH2.EQ.0) THEN
         OPEN(5,FILE='odata',STATUS='OLD')
      ELSE
         WRITE(OTEMP,*) FILTH2
         WRITE(OSTRING,'(A)') 'odata.' // TRIM(ADJUSTL(OTEMP))
         OPEN(5,FILE=OSTRING,STATUS='OLD')
      ENDIF

190   CALL INPUT(END)
      IF (.NOT. END) CALL READU(WORD)
      IF (END.OR.WORD.EQ.'STOP'.OR.WORD.EQ.'POINTS') GOTO 200
      IF (WORD.EQ.'    ' .OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'.OR. WORD .EQ. '\\') THEN
         GOTO 190
      ELSE IF ((WORD.EQ.'CPMD').OR.(WORD.EQ.'CPMDC')) THEN
         CPMD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') ' ERROR - no CPMD system specified'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 112
            ENDIF
         ENDDO
112      CONTINUE
      ELSE IF ((WORD.EQ.'ONETEP').OR.(WORD.EQ.'ONETEPC')) THEN
         ONETEP=.TRUE. 
         IF (WORD.EQ.'ONETEP') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(ONETEPJOB)
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') 'getparams> ERROR - ONETEP input mangled'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 231
            ENDIF
         ENDDO
231      CONTINUE
      ELSE IF ((WORD.EQ.'CP2K').OR.(WORD.EQ.'CP2KC')) THEN 
         CP2K=.TRUE.
         IF (WORD.EQ.'CP2K') DFTP=.TRUE. 
         IF (NITEMS.GT.2) THEN
            CALL READA(CP2KJOB)
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') 'getparams> ERROR - CP2K input mangled'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 271
            ENDIF
         ENDDO
271      CONTINUE

      ELSE IF (WORD.EQ.'REAXFF') THEN
         REAXFF=.TRUE.
         CALL READA(REAXFFJOB)
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 814
            ENDIF
         ENDDO
814      CONTINUE


      ELSE IF (WORD.EQ.'MOLPRO') THEN
         MOLPRO=.TRUE.
         CALL READA(MOLPROJOB)
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 815
            ENDIF
         ENDDO
815      CONTINUE
         IF(NITEMS.GT.3) THEN
            CALL READA(MOLPROJOBPARAMS)
         ELSE
            MOLPROJOBPARAMS=''
         ENDIF


      ELSE IF (WORD.EQ.'QCHEM') THEN
         QCHEM=.TRUE.
         CALL READA(QCHEMJOB)
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 212
            ENDIF
         ENDDO
212      CONTINUE
         IF(NITEMS.GT.3) THEN
            CALL READA(QCHEMJOBPARAMS)
         ELSE
            QCHEMJOBPARAMS=''
         ENDIF


      ELSE IF (WORD.EQ.'QCHEMES') THEN
         QCHEM=.TRUE.
         CALL READA(QCHEMJOB)
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 912
            ENDIF
         ENDDO
912      CONTINUE
         QCHEMES=.TRUE.
         CALL READI(QCHEMESNAO)
         CALL READI(QCHEMESNMO)
         CALL READI(QCHEMESNELEC)
      ELSE IF (WORD.EQ.'VASP') THEN
         VASP=.TRUE.
         CALL READA(VASPJOB)
      ELSE IF ((WORD.EQ.'CASTEP').OR.(WORD.EQ.'CASTEPC')) THEN
         CASTEP=.TRUE.
         IF (WORD.EQ.'CASTEP') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(CASTEPJOB)
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') 'getparams> ERROR - CASTEP input mangled'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 211
            ENDIF
         ENDDO
211      CONTINUE
      ELSE IF (WORD.EQ.'RINGPOLYMER') THEN
         RINGPOLYMERT=.TRUE.
         GOTO 200
      ELSE IF (WORD.EQ.'VARIABLES') THEN
         VARIABLES=.TRUE.
         GOTO 200
      ELSE IF (WORD .EQ. 'AMBER12') THEN
         AMBER12T = .TRUE.
         WRITE(OSTRING,'(A)') 'coords.inpcrd'
         IF (NITEMS == 3) THEN
            CALL READA(FIRSTAMBERSTR)
            IF (FILTH2 .NE. 0) THEN
               WRITE(OTEMP, *) FILTH2
               WRITE(OSTRING,'(A)') TRIM(ADJUSTL(FIRSTAMBERSTR))//'.'//TRIM(ADJUSTL(OTEMP))
               WRITE(*,*) 'ostring=', OSTRING
            ELSE
               WRITE(OSTRING,'(A)') TRIM(ADJUSTL(FIRSTAMBERSTR))
            END IF
            CALL READA(FIRSTAMBERSTR)
            IF (TRIM(ADJUSTL(FIRSTAMBERSTR)) /= 'inpcrd') THEN
               WRITE(OSTRING,'(A)') 'coords.inpcrd'
            END IF
         END IF
         CALL AMBER12_SETUP(NATOMS, OSTRING, LEN(OSTRING))
! Go to the end of the routine once we have the correct number of atoms.
         GOTO 400
      ELSE IF (WORD.EQ.'AMBER9'.OR.(WORD.EQ.'NAB')) THEN
         IF(WORD.EQ.'AMBER9') AMBERT=.TRUE.
         IF(WORD.EQ.'NAB') NABT=.TRUE.
!         NATOMS=0
        INPCRD='coords.inpcrd'
       IF(NITEMS<3) then
         CALL READA(amberstr)
          IF (FILTH2.NE.0) THEN
            WRITE(OTEMP,*) FILTH2
            write(ostring,'(A)') trim(adjustl(amberstr))//'.'//trim(adjustl(otemp))
            WRITE(*,*) 'ostring=', ostring
          ELSE
            write(ostring,'(A)') trim(adjustl(amberstr))
          END IF
          WRITE(*,'(A)') ' getparams> input coordinates for AMBER9 system will be read from ',trim(adjustl(amberstr)),ostring
         CALL AMBERINTERFACE(NATOMS,2,INPCRD,6)
         CALL amber_readcoords(ostring)
       ELSE IF(NITEMS==3) then
         CALL READA(amberstr)
         CALL READA(amberstr1)
          WRITE(*,'(A)') ' getparams> input coordinates for AMBER9 system will be read from ', trim(adjustl(amberstr)),
     &                         'type: ', trim(adjustl(amberstr1))
          IF(TRIM(ADJUSTL(AMBERSTR1)).EQ.'inpcrd') then
           INPCRD=AMBERSTR
           WRITE(*,'(A)') ' getparams> reading AMBER inpcrd coordinate format'
          ELSE
           WRITE(*,'(A)') ' getparams> ERROR - no other types defined currently than inpcrd'
          END IF
           CALL AMBERINTERFACE(NATOMS,2,INPCRD,6)
       END IF
        GOTO 400           ! finished counting atoms, go to the end of the subroutine
      ELSE IF (WORD.EQ.'AMH') THEN

         AMH=.TRUE.
         NATOMS=0
         WRITE(6,*)'Entering GETPARAMS'

         OPEN(UNIT=30,FILE='pro.list',STATUS='OLD',FORM='FORMATTED')
         READ (30,1000)TARFL
1000     FORMAT(A5)
         CLOSE(30)

         OPEN(30,FILE='proteins/'//TARFL,STATUS='OLD')
            READ(30,*)
            READ(30,*)NRES
            IF (NRES.GT.500) THEN
                WRITE(6,*) 'FAILURE NRES GR THAN 500 COUNTATOMS'
                STOP
            ENDIF
            READ (30,25)(SEQ(I_RES),I_RES=1,NRES)
25         FORMAT(25(I2,1X))
          CLOSE(30)
               
          NOGLY = 0
          NUMGLY = 0
           DO I_RES=1,NRES
             IF (SEQ(I_RES).NE.8) NOGLY = NOGLY +1
             IF (SEQ(I_RES).EQ.8) NUMGLY = NUMGLY +1
           ENDDO
            NATOMS = NOGLY*3 + NUMGLY*2
            WRITE(6,*)'NATOMS NOGLY  NUMGLY  ',NATOMS,NOGLY,NUMGLY
           GOTO 400
!
!OPEP potential - use conf_initiale.pdb as we need it anyway 
!
       ELSE IF (WORD.EQ.'OPEP') THEN
           OPEPT=.TRUE.
           !leave the initialisation for keywords, only get natoms here
           CALL OPEP_GET_NATOMS(NATOMS)
           GOTO 400

C davidg: introduced userpot here:
       ELSE IF (WORD.EQ.'USERPOT') THEN
C           PRINT *,'IN THIS IF!'
           USERPOTT=.TRUE.
           CALL USERPOT_INIT
           CALL USERPOT_GET_NATOMS(NATOMS)
           GOTO 400
           
       ELSE IF (WORD.EQ.'CHARMM') THEN
! DAE We are going to assume that there is a charmm format file in the directory called input.crd, and the first
! line of it will tell us the number of atoms. In 99% of old (OPTIM<3) runs this is what I did, but the old versions
! were more flexible in that any filename could be specified in the CHARMM bit of the odata file

         IF (FILTH2.NE.0) THEN
            WRITE(OTEMP,*) FILTH2
            WRITE(OSTRING,'(A)') 'input.crd.' // TRIM(ADJUSTL(OTEMP))
         ELSE
            WRITE(OSTRING,'(A)') 'input.crd'
         ENDIF
         OPEN (UNIT=9,FILE=OSTRING,STATUS='OLD',IOSTAT=ios)

         IF (ios /= 0) THEN
            WRITE(OSTRING,'(A)') 'input.crd'
            OPEN (UNIT=9,FILE=OSTRING,STATUS='OLD',IOSTAT=ios)
            if (ios == 0) THEN
         else
            WRITE(*,'(2A)') 'Thanks to our new dynamic memory allocation overlords, there must be a charmm-format file called ',
     &    '"input.crd" for CHARMM to find out the number of atoms. Feel free to recode to enable any filename to work properly.'
            STOP
         ENDIF
         ENDIF
         do
              read(9,*) myline
              if (myline(1:1)=='*') then ! SAT This is the goddamn CHARMM comment line
                    cycle
              else
                    read(myline,*) NATOMS
                    exit
              endif
         enddo
         CLOSE(9)

! DAE We also need to find out what MAXAIM is in CHARMM, and set MXATMS in OPTIM to be the same, so that those arrays which
! are passed between the two can be declared correctly. MXATMS is now stored in modmxatms.
       
         CALL GETMAXAIM

         GOTO 400

      ELSE IF (WORD.EQ.'UNRES') THEN
! jmc We are going to assume that there is a coords file. The first line of it will tell us the number of atoms.

         IF (FILTH2.EQ.0) THEN
            OPEN (UNIT=9,FILE='coords',STATUS='OLD',IOSTAT=ios)
         ELSE
            WRITE(CSTRING,'(A)') 'coords.'//TRIM(ADJUSTL(OTEMP))
            OPEN (UNIT=9,FILE=CSTRING,STATUS='OLD',IOSTAT=ios)
         ENDIF
         IF (ios /= 0) THEN
            WRITE(*,'(2A)') 'Thanks to our new dynamic memory allocation overlords, there must be a coords file present ',
     &    ' for OPTIM3 to find out the number of atoms for UNRES. Please make it so!'
            STOP
         ENDIF
         READ(9,*) NATOMS
         CLOSE(9)

         GOTO 400

      ENDIF
      GOTO 190

200   CONTINUE

      IF (VASP) THEN
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file POSCAR'
         CALL SYSTEM('grep "TITEL" POTCAR | wc > temp_NTYPE') !Searches POTCAR for 'TITEL', to calculate how many element types (NTYPES) are present)
         FNAME='temp_NTYPE'
         OPEN(UNIT=7, FILE=FNAME, STATUS='OLD')
         READ (7,*) NTYPES
         CLOSE (7)
   
         ALLOCATE(ELEMENT_NUMS(NTYPES))

         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE='POSCAR',STATUS='OLD')
         DO J1=1,6                                                  !First 2 lines are always the same
            READ(LUNIT,*)
         ENDDO
         READ(LUNIT,*) (ELEMENT_NUMS(J1), J1=1,NTYPES)                  
         NATOMS=SUM(ELEMENT_NUMS)                                   !Will be taken from seperate routine

! Now we need to count the atoms in POSCAR by reading through the file
! and set NATOMS
!
         CLOSE(LUNIT)
         WRITE(*,'(A,I4,A)') 'getparams> VASP run for ',NATOMS,' atoms'
         CONTINUE

      ELSEIF (MOLPRO) THEN
         FNAME=SYS(1:LSYS)
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE=FNAME,STATUS='OLD')         
         MOLPROSCALE=0.5291772D0
         NATOMS=0
         coordsq2: DO
            READ(LUNIT,'(A19)') DSTRING1
            CALL UPPERCASE(DSTRING1)
            IF (DSTRING1(1:9).EQ.' GEOMETRY') THEN
               DO
                  READ(LUNIT,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
                  IF (DSTRING1(1:1).EQ.'}') EXIT coordsq2
                  NATOMS=NATOMS+1
               ENDDO
               WRITE(*,'(A,I4,A)') 'getparams> MOLPRO run for ',NATOMS,' atoms'
               EXIT
            ENDIF
         ENDDO coordsq2
         CLOSE(LUNIT)         

      ELSEIF (QCHEMES) THEN
         NATOMS=QCHEMESNAO*QCHEMESNELEC ! this mimics the behaviour for VARIABLES keyword, and should get declarations right
         WRITE(*,'(A,I0)') ' getparams> Setting the number of variables to ',NATOMS
      ELSEIF (QCHEM) THEN
         FNAME=SYS(1:LSYS) 
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         LUNIT=GETUNIT()
         OPEN(UNIT=LUNIT,FILE=FNAME,STATUS='OLD')

         QCHEMSCALE=0.5291772D0 ! coordinates in Angstrom, gradients in a.u.
         bohr: DO
            READ(LUNIT,'(A19)',END=864) DSTRING1
            CALL UPPERCASE(DSTRING1)
            IF (DSTRING1(1:15).EQ.'INPUT_BOHR TRUE') THEN
               QCHEMSCALE=1.0D0 ! coordinates and gradient in atomic units
               PRINT '(A)',' getparams> coordinates in atomic units - no gradient/second derivative scaling required'
               EXIT bohr
            ENDIF
         ENDDO bohr
864      REWIND(LUNIT)

         NATOMS=0
         coordsq: DO
            READ(LUNIT,'(A19)') DSTRING1
            CALL UPPERCASE(DSTRING1)
            IF (DSTRING1(1:4).EQ.'$MOL') THEN
               READ(LUNIT,'(A21)') DSTRING1
               DO
                  READ(LUNIT,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
                  IF (DSTRING1(1:4).EQ.'$END') EXIT coordsq
                  NATOMS=NATOMS+1
               ENDDO
               WRITE(*,'(A,I4,A)') 'getparams> QCHEM run for ',NATOMS,' atoms'
               EXIT
            ENDIF
         ENDDO coordsq
         CLOSE(LUNIT)
      ELSEIF (CASTEP) THEN
         FNAME=SYS(1:LSYS) // '.cell'
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop: DO
            READ(7,'(A21)') DSTRING1
            CALL UPPERCASE(DSTRING1)
C           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:16).EQ.'%BLOCK POSITIONS') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
C                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:19).EQ.'%ENDBLOCK POSITIONS') EXIT coordsloop
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop
C        WRITE(*,'(A,I4,A)') 'getparams> CASTEP run for ',NATOMS,' atoms'
         CLOSE(7)
      ELSEIF (ONETEP) THEN
         FNAME=SYS(1:LSYS) // '.dat'
         WRITE(*,'(A,A)') 'getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop2: DO
            READ(7,'(A21)') DSTRING1
            CALL UPPERCASE(DSTRING1)
C           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:16).EQ.'%BLOCK POSITIONS') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
C                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:19).EQ.'%ENDBLOCK POSITIONS') EXIT coordsloop2
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop2
         WRITE(*,'(A,I4,A)') 'getparams> ONETEP run for ',NATOMS,' atoms'
         CLOSE(7)
      ELSEIF (CP2K) THEN 
         FNAME=SYS(1:LSYS) // '.inp'
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop3: DO
            READ(7,'(A21)') DSTRING1
            DSTRING1=ADJUSTL(DSTRING1)
            CALL UPPERCASE(DSTRING1)
C           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:6).EQ.'&COORD') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  DSTRING1=ADJUSTL(DSTRING1)
                  CALL UPPERCASE(DSTRING1)
C                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:10).EQ.'&END COORD') EXIT coordsloop3
                  IF (DSTRING1(1:6).EQ.'SCALED') NATOMS=NATOMS-1
                  IF (DSTRING1(1:4).EQ.'UNIT') NATOMS=NATOMS-1
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop3
         WRITE(*,'(A,I4,A)') ' getparams> CP2K run for ',NATOMS,' atoms'
         CLOSE(7) 
      ELSE IF (CPMD) THEN
         FNAME=SYS(1:LSYS)
         WRITE(*,'(A,A)') 'getparams> Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
11       READ(7,'(A)') FNAME
         IF (FNAME(1:6).EQ.'&ATOMS') THEN
            J1=0
12          READ(7,'(A)') TSTRING
            IF (TSTRING(1:1).EQ.'*') THEN
               J1=J1+1
               READ(7,'(A)') FNAME
               READ(7,*) NTYPE(J1)
               DO J2=1,NTYPE(J1)
                  READ(7,*) DUMMY
               ENDDO
               NATOMS=NATOMS+NTYPE(J1)
               GOTO 12
            ELSE IF (TSTRING(1:1).EQ.' ') THEN
               GOTO 12
            ELSE IF (TSTRING(1:4).EQ.'&END') THEN
               GOTO 13
            ENDIF
         ELSE
            GOTO 11
         ENDIF
13       CONTINUE
      ELSE IF (VARIABLES) THEN
         NATOMS=0
110      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READF(DUMMY)
            NATOMS=NATOMS+1
            GOTO 110
         ENDIF
      ELSE IF (RINGPOLYMERT) THEN
         NATOMS=0
111      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READF(DUMMY)
            NATOMS=NATOMS+NITEMS
            GOTO 111
         ENDIF
      ELSE 
         NATOMS=0
300      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READA(ZDUM)
            IF (ZDUM(1:2).EQ.'  ') GOTO 300
            CALL READF(DUMMY)
            CALL READF(DUMMY)
            CALL READF(DUMMY)
            NATOMS=NATOMS+1
            GOTO 300
         ENDIF
      ENDIF
400   CLOSE(5)
  
      WRITE(*,'(A,I6)') ' getparams> Number of atoms (or variables)  determined as ',NATOMS
      NOPT=3*NATOMS
      IF (VARIABLES) NOPT=NATOMS
      WRITE(*,'(A,I6)') ' getparams> Number of optimisation degrees of freedom ',NOPT

      CALL OPTIM(FILTH,FILTH2,ARGSTRING)


      STOP      
      END PROGRAM OPTIM3
