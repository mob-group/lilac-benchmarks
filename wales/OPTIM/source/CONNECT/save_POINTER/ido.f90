!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
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
MODULE IDOMODULE
     USE CONNECTDATA
     USE CONNECTUTILS
     USE KEYCONNECT
     IMPLICIT NONE
     CONTAINS

SUBROUTINE INITIALISE(NA,EI,Q,EF,FIN,ENDPOINTSEP)
     USE KEY, ONLY: READSP
     IMPLICIT NONE
     
     INTEGER,INTENT(IN)           :: NA
     DOUBLE PRECISION,INTENT(IN)           :: ENDPOINTSEP
     DOUBLE PRECISION,POINTER              :: EI,EF
     DOUBLE PRECISION,POINTER,DIMENSION(:) :: Q,FIN

     INTEGER J2, NPLUS, NMINUS, MINPOS, NMINA, NMINB, NTSDUMP, NDUMMY, NMINDUMP, IOERROR
     LOGICAL MINNEW
     DOUBLE PRECISION, POINTER :: ETEMP, LOCALPOINTS(:)
     DOUBLE PRECISION, TARGET :: ZEROTARGET
     DOUBLE PRECISION  DJWDUMMY
     INTEGER, ALLOCATABLE :: ACTUALMIN(:)
     LOGICAL PERMUTE
     INTEGER INVERT,INDEX(NATOMS),IMATCH

     NULLIFY(ETEMP, LOCALPOINTS)
     NATOMS=NA
     NOPT=3*NATOMS
     ALLOCATE(G(NOPT),MI(MINRACKSIZE),TS(TSRACKSIZE))

     ! endpoints initialisation
     NMIN=2
     MI(1)%DATA%E => EI
     MI(2)%DATA%E => EF
     MI(1)%DATA%X => Q
     MI(2)%DATA%X => FIN
     ALLOCATE(MI(2)%DATA%D(1),MI(2)%DATA%NTRIES(1))
     MI(2)%DATA%D(1) = ENDPOINTSEP
     IF (INTERPCOSTFUNCTION) THEN
        ALLOCATE(MI(2)%DATA%INTERP(1))
        MI(2)%DATA%INTERP(1)=INTERPVALUE(Q,FIN,ENDPOINTSEP)
     ENDIF
     MI(2)%DATA%NTRIES(1)=0
     MI(1)%DATA%S=.TRUE.
     MI(1)%DATA%F=.FALSE.
     MI(2)%DATA%S=.FALSE.
     MI(2)%DATA%F=.TRUE.
     NULLIFY( MI(1)%DATA%CTS,MI(1)%DATA%CMIN,MI(2)%DATA%CTS,MI(2)%DATA%CMIN )

     ! S and F are not connected to anything yet.
     ALLOCATE(MI(1)%DATA%C(1),MI(2)%DATA%C(1))
     MI(1)%DATA%C(1) = 0
     MI(2)%DATA%C(1) = 0

     NTS=0
     DEPTH=0 !RECURSIVE SUBROUTINES LEVEL INITIALIZATION

     NULLIFY(START,NEW,DUMMY,START2,NEW2,DUMMY2,ONEUP)

     ! path subroutine arrays allocations
     ALLOCATE(FRQSTS(3*NATOMS), EVEC(3*NATOMS), FRQSPLUS(3*NATOMS), FRQSMINUS(3*NATOMS), GDUMMY(3*NATOMS))
     FRQSTS(1:3*NATOMS)=0.0D0 ! BUGS FIXED 21/12/05 DJW
     EVEC(1:3*NATOMS)=0.0D0
     FRQSPLUS(1:3*NATOMS)=0.0D0
     FRQSMINUS(1:3*NATOMS)=0.0D0
     GDUMMY(1:3*NATOMS)=0.0D0

     IF (READSP) THEN

! DJW new provision to read minima and ts dumped by previous OPTIM or PATHSAMPLE run
! First we need to know how many minima are in min.A, min.B and min.data

        PRINT'(A)',' Reading previous stationary point information'
        NMINA=0
        OPEN(UNIT=82,FILE='min.A',STATUS='OLD')
        DO
          READ(82,*,END=11) DJWDUMMY
          NMINA=NMINA+1
        ENDDO    
11      CLOSE(UNIT=82)

        NMINB=0
        OPEN(UNIT=82,FILE='min.B',STATUS='OLD')
        DO
          READ(82,*,END=12) DJWDUMMY
          NMINB=NMINB+1
        ENDDO    
12      CLOSE(UNIT=82)

        NMINDUMP=0
        OPEN(UNIT=82,FILE='min.data',STATUS='OLD')
        DO
          READ(82,*,END=13) DJWDUMMY
          NMINDUMP=NMINDUMP+1
        ENDDO    
13      CLOSE(UNIT=82)

        DO ! ALLOCATE ENOUGH SPACE FOR ALL THE POSSIBLE MINIMA
           ! addnewmin would actually have done this anyway!
          IF (MINRACKSIZE.GE.2+NMINA+NMINB+NMINDUMP) EXIT
          CALL REALLOCATEMINRACK
        ENDDO

! Now we know how big ACTUALMIN needs to be

        ALLOCATE(ACTUALMIN(NMINA+NMINB+NMINDUMP))

! Read in all minima from min.A, min.B and min.data, assigning
! energy, coordinates, distances and number of tries using addnewmin

        INQUIRE(IOLENGTH=NDUMMY) Q(1:NOPT)
        OPEN(13,FILE='points.min',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=NDUMMY)

        NDUMMY=0 ! NDUMMY COUNTS MINIMA IN SAVED DATABASE
        OPEN(UNIT=82,FILE='min.A',STATUS='OLD')
        DO
          ALLOCATE(ETEMP,LOCALPOINTS(1:NOPT))
          READ(82,*,END=111) ETEMP
          NDUMMY=NDUMMY+1
          READ(13,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
          CALL ISNEWMIN(ETEMP,LOCALPOINTS,MINPOS,MINNEW,.FALSE.,PERMUTE,INVERT,INDEX,IMATCH)
          IF (MINNEW) THEN
             ACTUALMIN(NDUMMY)=NMIN+1 ! ACTUALMIN(DATABASE INDEX)=INDEX IN OPTIM
             CALL ADDNEWMIN(ETEMP,LOCALPOINTS)
          ELSE
             ACTUALMIN(NDUMMY)=MINPOS 
          ENDIF
          NULLIFY(ETEMP,LOCALPOINTS)
        ENDDO    
111     CLOSE(UNIT=82)

        OPEN(UNIT=82,FILE='min.B',STATUS='OLD')
        DO
          ALLOCATE(ETEMP,LOCALPOINTS(1:NOPT))
          READ(82,*,END=121) ETEMP
          NDUMMY=NDUMMY+1
          READ(13,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
          CALL ISNEWMIN(ETEMP,LOCALPOINTS,MINPOS,MINNEW,.FALSE.,PERMUTE,INVERT,INDEX,IMATCH)
          IF (MINNEW) THEN
             ACTUALMIN(NDUMMY)=NMIN+1 ! ACTUALMIN(DATABASE INDEX)=INDEX IN OPTIM
             CALL ADDNEWMIN(ETEMP,LOCALPOINTS)
          ELSE
             ACTUALMIN(NDUMMY)=MINPOS
          ENDIF
          NULLIFY(ETEMP,LOCALPOINTS)
        ENDDO    
121     CLOSE(UNIT=82)

        OPEN(UNIT=82,FILE='min.data',STATUS='OLD')
        DO
          ALLOCATE(ETEMP,LOCALPOINTS(1:NOPT))
          READ(82,*,END=131) ETEMP ! AS SOON AS ETEMP IS RESET, SO ARE ALL THE ENERGIES THAT POINT AT IT!!!!
          NDUMMY=NDUMMY+1
          READ(13,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,NOPT)
          CALL ISNEWMIN(ETEMP,LOCALPOINTS,MINPOS,MINNEW,.FALSE.,PERMUTE,INVERT,INDEX,IMATCH)
          IF (MINNEW) THEN
             CALL ADDNEWMIN(ETEMP,LOCALPOINTS)
             ACTUALMIN(NDUMMY)=NMIN ! ACTUALMIN(DATABASE INDEX)=INDEX IN OPTIM
          ELSE
             ACTUALMIN(NDUMMY)=MINPOS
          ENDIF
          NULLIFY(ETEMP,LOCALPOINTS) ! THIS LINE ALLOWS THE CURRENT VALUES IN POINTERS THAT POINT TO
                                     ! these variables to be preserved: they are not overwritten when
                                     ! ETEMP,LOCALPOINTS are reset for the next new minimum
                                     ! The association with the original memory locations is broken
                                     ! for ETEMP and LOCALPOINTS, but mi%data%E and mi%data%X remain
                                     ! pointing to the correct locations.
        ENDDO    
131     CLOSE(UNIT=82)

        CLOSE(13)
        INQUIRE(IOLENGTH=NDUMMY) LOCALPOINTS(1:NOPT)
        OPEN(15,FILE='points.ts',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=NDUMMY,IOSTAT=IOERROR)

        IF (IOERROR /= 0) THEN
            PRINT *,' Please be advised that the points.ts file is absent. Is this deliberate?'
            GOTO 141
        ENDIF

        NTSDUMP=0
        OPEN(UNIT=82,FILE='ts.data',STATUS='OLD')
        ZEROTARGET=0.0D0 
        DO
          ALLOCATE(ETEMP,LOCALPOINTS(1:NOPT))
          READ(82,*,END=14) ETEMP,DJWDUMMY,NDUMMY,NPLUS,NMINUS
          NTSDUMP=NTSDUMP+1
          READ(15,REC=NTSDUMP) (LOCALPOINTS(J2),J2=1,NOPT)
          IF (NTS==TSRACKSIZE) CALL REALLOCATETSRACK
          NTS=NTS+1
          TS(NTS)%DATA%E => ETEMP
          TS(NTS)%DATA%X => LOCALPOINTS
          TS(NTS)%DATA%EVALMIN => ZEROTARGET ! WE DON;T ASSUME THE E/VALUE AND E/VECTOR WERE SAVED
                                             ! must NOT use pointer = thing - pointer is dereferenced, not set!
                                             ! Same thing would be needed for pointer vecs, but we shouldn't use it!
          TS(NTS)%DATA%P = ACTUALMIN(NPLUS)  ! THE FOLLOWING COMPONENTS ARE NOT POINTERS!
          TS(NTS)%DATA%M = ACTUALMIN(NMINUS)
          ! set zero distance between connected minima for dijkstra weight
          CALL SETDISTANCE(ACTUALMIN(NPLUS),ACTUALMIN(NMINUS),0.0D0)
          IF (INTERPCOSTFUNCTION) CALL SETINTERP(ACTUALMIN(NPLUS),ACTUALMIN(NMINUS),0.0D0)
          TS(NTS)%DATA%SLENGTH = 0.0D0
          TS(NTS)%DATA%SLP = 0.0D0
          TS(NTS)%DATA%DISP = 0.0D0
          TS(NTS)%DATA%GAMMA = 0.0D0
          TS(NTS)%DATA%NTILDE = 0.0D0
          CALL NEWCONNECTION(ACTUALMIN(NPLUS),ACTUALMIN(NMINUS),NTS)
          NULLIFY(ETEMP,LOCALPOINTS)
        ENDDO    
14      CLOSE(UNIT=82)
        CLOSE(15)

        DEALLOCATE(ETEMP,LOCALPOINTS)
        DEALLOCATE(ACTUALMIN)

! DJW end of old minima and transition states

141  CONTINUE

     ENDIF
END SUBROUTINE INITIALISE

SUBROUTINE DEINITIALISE
     IMPLICIT NONE

     !DOUBLE PRECISION,POINTER :: P(:)

     DEALLOCATE(G,MI,TS)
     DEALLOCATE(FRQSTS,EVEC,FRQSPLUS,FRQSMINUS,GDUMMY)
     !do i=1,Nts
     !     p => ts(i)%data%X
     !     nullify(ts(i)%data%X)
     !     deallocate(p)
     !enddo
END SUBROUTINE DEINITIALISE

SUBROUTINE OUTPUT
     USE PORFUNCS
     USE KEYCONNECT
     USE CONNECTUTILS
     USE KEY, ONLY: DUMPSP, MACHINE, DUMPALLPATHS
     IMPLICIT NONE

     INTEGER :: I
     INTEGER,POINTER :: FINALPATHTS(:)
     DOUBLE PRECISION :: ESAVE,NTILDEAV
     DOUBLE PRECISION :: ETSMAX,SMAX,DMAX,NTILDEMAX ! WHERE E IS ENERGY, S IS THE INTEGRATED PATHLENGTH AND D IS THE DISTANCE
     DOUBLE PRECISION :: ETSMIN,SMIN,DMIN,NTILDEMIN ! BETWEEN MINIMA, AND NTILDE IS THE COOPERATIVITY INDEX
     DOUBLE PRECISION :: BARPLUSMAX,BARPLUSMIN,BARMINUSMAX,BARMINUSMIN
     DOUBLE PRECISION :: BAIMAX,BAIMIN,PAIMAX,PAIMIN ! BARRIER AND PATHLENGTH ASYMMETRY INDECES: {0..1}; 0 - ABSOLUTE SYMMETRY

     LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST
     INTEGER NCONNECT
     DOUBLE PRECISION TEMPERATURE, HRED
     COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
     
     NULLIFY(FINALPATHTS)
     ETSMIN = HUGE(ETSMAX)
     SMIN   = HUGE(SMAX)
     DMIN   = HUGE(DMAX)
     NTILDEMIN = HUGE(NTILDEMAX)
     ETSMAX = -ETSMIN
     SMAX   = -SMIN
     DMAX   = -DMIN
     NTILDEMAX = -NTILDEMIN
     NTILDEAV = 0.0D0
     BARPLUSMIN = HUGE(BARPLUSMIN)
     BARMINUSMIN = HUGE(BARMINUSMIN)
     BARPLUSMAX = -BARPLUSMIN
     BARMINUSMAX = -BARMINUSMIN
     BAIMIN = HUGE(BAIMIN)
     PAIMIN = HUGE(PAIMIN)
     BAIMAX = -BAIMIN
     PAIMAX = -PAIMIN
     IF (DUMPALLPATHS) CLOSE(88)

     IF (FINISHED) THEN
          PRINT *, 'Connected path found'
          IF (ASSOCIATED(START)) THEN
               DUMMY=>START
          ELSE
               PRINT *, 'Start is _not_ present!!!'
               RETURN
          ENDIF

          WRITE(*,'(A)')  &
          & '  ts        E+         Ets - E+          Ets       Ets - E-          E-          S       D' // &
          & '      gamma   ~N'

          DUMMY=>START
          I=1
          DO
               IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
               ESAVE = TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%E
               WRITE(*,'(I4,1X,F16.7,G12.5,F16.7,G12.5,F16.7,F8.3,F8.3,F8.3,F8.3)') &
               & MI(DUMMY%I)%DATA%CTS(DUMMY%J), &
               & MI(DUMMY%I)%DATA%E, ESAVE - MI(DUMMY%I)%DATA%E, ESAVE, &
               & ESAVE-MI(DUMMY%NEXT%I)%DATA%E,MI(DUMMY%NEXT%I)%DATA%E, &
               & TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%SLENGTH, &
               & TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%DISP, &
               & TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%GAMMA, &
               & TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%NTILDE
               DUMMY=>DUMMY%NEXT

               I=I+1
          ENDDO
 
          WRITE(*,'(/1x,a,i7)') "Number of TS in the path       =",i-1
          ALLOCATE(FINALPATHTS(I-1))
          DUMMY=>START
          I=1
          DO
               IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
               FINALPATHTS(I)=MI(DUMMY%I)%DATA%CTS(DUMMY%J)
               DUMMY=>DUMMY%NEXT
               I=I+1
          ENDDO

          WRITE(*,'(1x,a,i7)') "Number of cycles               =",NConDone
          IF (DUMPSP) CALL DODUMPSP
          CALL MERGEXYZEOFS
          IF (DUMPPATH) CALL MAKEPATHINFO
     ELSE
          PRINT *,'Number of allowed connection steps reached - exit newconnect'
          IF (DUMPSP) CALL DODUMPSP
          CALL TSUMMARY
          CALL EXIT(10)
     ENDIF

     IF (MACHINE) CALL DUMPDB(FINISHED,FINALPATHTS) ! DON't pour the water of wisdom on the dry sand
END SUBROUTINE OUTPUT

END MODULE IDOMODULE
