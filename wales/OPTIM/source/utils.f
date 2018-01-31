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

      SUBROUTINE NORMAL(X,N)
C
C NORMALIZES VECTOR X
C
      IMPLICIT NONE
      INTEGER I,N
      DOUBLE PRECISION X(N),Q,P

      Q=0.D0
      DO I=1,N
         Q=Q+X(I)**2
      ENDDO
      P=DSQRT(Q)
      IF(P.LT.1D-14)THEN
      WRITE(6,*)' null vector returned from NORMAL'
      RETURN
      ENDIF
      DO I=1,N
         X(I)=X(I)/P    
      ENDDO
      RETURN
      END 
C
C  Normalize vector VEC1
C
      SUBROUTINE VECNORM(VEC1,NOPT)
      IMPLICIT NONE
      INTEGER J2, NOPT
      DOUBLE PRECISION DUMMY2, VEC1(NOPT)

      DUMMY2=0.0D0
      DO J2=1,NOPT
         DUMMY2=DUMMY2+VEC1(J2)**2
      ENDDO
!     PRINT '(A,G20.10)','vecnorm> DUMMY2=',DUMMY2
      IF (DUMMY2.GT.0.0D0) THEN
         DUMMY2=1.0D0/DSQRT(DUMMY2)
         DO J2=1,NOPT
            VEC1(J2)=VEC1(J2)*DUMMY2
         ENDDO
      ENDIF

      RETURN
      END
    
C      SUBROUTINE CALCVEC(A,B,V,IX)
CC
CC CALCULATES THE VECTOR V BETWEEN CARTESIAN POINTS A AND B
CC
C      IMPLICIT NONE
C      DOUBLE PRECISION A(3),B(3),V(3)
C      INTEGER IX,I
C
C      DO I=1,3
C         V(I)=B(I)-A(I)
C      ENDDO
C      IF(IX.EQ.1)CALL NORMAL(V,3)
C
C      RETURN
C      END 

      SUBROUTINE VADD(A,B,C,N,IP)
      IMPLICIT NONE
      INTEGER N,IP,I
      DOUBLE PRECISION A(N),B(N),C(N)

      DO I=1,N
         A(I)=B(I)+C(I)*IP
      ENDDO
      RETURN
      END         

      SUBROUTINE VSTAT(V,ZQ,LENGTH,N)
C
C RETURNS STATISTICAL INFO ABOUT VECTOR V in ZQ
C     ZQ(1)  Largest absolute magnitude
C     ZQ(2)  Smallest absolute magnitude
C     ZQ(3)  Largest value
C     ZQ(4)  Smallest value
C     ZQ(5)  2-norm
C     ZQ(6)  Dynamic range of the vector (abs. min. - abs. max.)
C   

      IMPLICIT NONE
      INTEGER N,I,LENGTH
      DOUBLE PRECISION V(N),ZQ(6),U

      U=0.D0
      ZQ(1:6)=0.0D0
      ZQ(2)=DABS(V(1))
      ZQ(4)=V(1)
C     PRINT*,'Statistics reported for vector:'
C     WRITE(*,10) (J,V(J),J=1,LENGTH)
C10    FORMAT(I4,F20.10)
      DO I=1,LENGTH
         ZQ(1)=MAX(ZQ(1),DABS(V(I)))
         ZQ(2)=MIN(ZQ(2),DABS(V(I)))
         ZQ(3)=MAX(ZQ(3),V(I))
         ZQ(4)=MIN(ZQ(4),V(I))
         U=U+V(I)*V(I)
      ENDDO
      If (Length .ne. 0.0) ZQ(5)=DSQRT(U/LENGTH)
      ZQ(6)=ZQ(2)-ZQ(1)
      RETURN
      END

       INTEGER FUNCTION ATOI (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Convert STRING to an integer value
C
C Arguments:    STRING   character string (input only)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C 
C Revision 4.0  89/03/14  01:15:20  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C 
C Revision 3.0  89/01/29  23:09:23  bernhold
C First working release for VAX
C 
C Revision 2.1  89/01/02  20:35:02  bernhold
C To keep consistent with .u file just checked in.
C 
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       CHARACTER(LEN=80) STRING
       CHARACTER(LEN=1) C
       INTEGER I
       LOGICAL NEG
C
       ATOI = 0
       NEG = .FALSE.
       I = 1
       C = STRING(I:I)
C
C      Pass over any leading spaces
C
 100   IF (C .EQ. ' ') THEN
          I = I + 1
          C = STRING(I:I)
          GOTO 100
       ENDIF
C
C      See if first character makes the number negative
C      Accept + as valid character before the digits start
C
       IF (C .EQ. '-') THEN
          NEG = .TRUE.
          I = I + 1
          C = STRING(I:I)
       ELSEIF (C .EQ. '+') THEN
          I = I + 1
          C = STRING(I:I)
       ENDIF
C
C      Continue as long as its a digit...
C
 200   IF (LGE(C, '0') .AND. LLE(C,'9')) THEN
C            Shift number over & add new digit
          ATOI = 10*ATOI + ICHAR(C)-48
          I = I + 1
          C = STRING(I:I)
          GOTO 200
       ENDIF
C
C      Negate the result if necessary
C
       IF (NEG) ATOI = -ATOI
       RETURN
       END
C*****************************************************
!
! multiply matrices: A^T=B^T * C
!
      SUBROUTINE MATMULV(A,B,C,NA,NB,NC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NA, NB, NC
      DOUBLE PRECISION, INTENT(IN) :: B(NB,NA), C(NB,NC)
      DOUBLE PRECISION, INTENT(OUT) :: A(NC,NA)
      INTEGER I,K,J
      DOUBLE PRECISION Z

      DO I=1,NA
         DO K=1,NC
            Z=0.D0
            DO J=1,NB
               Z=Z+B(J,I)*C(J,K)
C              PRINT*,'I,K,J,B,C,Z=',I,K,J,B(J,I),C(J,K),Z
            ENDDO
            A(K,I)=Z
         ENDDO
      ENDDO

      RETURN 
      END
       CHARACTER*(*) FUNCTION ITOA(NR, FRCPLS)
! Doxygen {{{
!> \name ITOA 
!> \brief Convert number to a left justified string
!> \param[in] NR number to be converted
!> \param[in] FRCPLS Force leading '+' if NR positive
! }}}
C Description {{{
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:      Convert NR to a left justified string
C
C Arguments:
C     NR       number to be converted (input only)
C     FRCPLS   Force leading '+' if NR positive (input only)
C
C Limitations:
C     May return with incomplete conversion if length of ITOA is too
C     short.  Puts '*' in last position of ITOA to indicate overlow.
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C Revision 4.0  89/03/14  01:15:45  bernhold
C Baseline for Sun & VAX prior to porting everywhere
C 
C Revision 3.0  89/01/29  23:10:22  bernhold
C First working release for VAX
C 
C Revision 2.1  89/01/02  20:36:12  bernhold
C To keep consistent with .u file just checked in.
C 
C     Revision 1.1  88/12/07  13:38:51  bernhold
C     Initial revision
C     
C
C System:       Standard FORTRAN 77
C
C Copyright 1988 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C }}}
      INTEGER NR, FRCPLS, NRABS
      INTEGER I, J, N, NDIG
C
C     Clear out the string
C
      DO 10 I = 1, LEN(ITOA)
         ITOA = ' '
 10   CONTINUE
C
C     Start counting position in string
C
      J = 1
      NRABS = ABS (NR)
C
C     Put in sign as appropriate
C
      IF (NR .LT. 0) THEN
         ITOA(J:J) = '-'
         J = J + 1
      ENDIF
      IF (FRCPLS .NE. 0 .AND. NR .GT. 0) THEN
         ITOA(J:J) = '+'
         J = J + 1
      ENDIF
C
C     Check if we are about to overflow the string
C
      IF (J .GT. LEN(ITOA)) THEN
         ITOA(J-1:J-1) = '*'
         RETURN
      ENDIF
C
C     Loop over nr of digits in number
C
      NDIG = INT( LOG10( FLOAT(NRABS)) ) + 1
      DO 100 I = NDIG, 1, -1
         N = MOD ( ( NRABS / (10**(I-1) ) ), 10)
         ITOA(J:J) = CHAR(N + 48)
         J = J + 1
C
C        Check for overflow of the string, but if this is last digit
C        then its okay.
C
         IF (J .GT. LEN(ITOA) .AND. I .GT. 1) THEN
            ITOA(J-1:J-1) = '*'
            RETURN
         ENDIF
 100  CONTINUE
C
      RETURN
      END
      CHARACTER(LEN=4) FUNCTION STRING(NAME,Iord)
      CHARACTER(LEN=4) NAME, ItoA*3
      INTEGER IORD
      If (Name(2:2) .eq. 'N') then
         String(1:1) = Name(1:1)
         String(2: ) = ItoA(Iord, 0)
         String(LNBlnk(String)+1:) = Name(3:3)
      Else
         String = Name
      EndIf
      RETURN
      END
      function dotopt(a,b,n)
      IMPLICIT NONE
      INTEGER I,N
      DOUBLE PRECISION A(N),B(N), DOTOPT

      DOTOPT=0.D0
      DO I=1,N
         DOTOPT=DOTOPT+A(I)*B(I)
      ENDDO
      RETURN
      END

      SUBROUTINE CROSSOPT(A,B,C,IX)
C
C CALCULATES THE (OPTIONALLY) NORMALIZED VECTOR CROSS PRODUCT C=A x B
C
      IMPLICIT NONE
      INTEGER IX
      DOUBLE PRECISION A(3),B(3),C(3)

      C(3)=A(1)*B(2)-B(1)*A(2)
      C(2)=-A(1)*B(3)+A(3)*B(1)
      C(1)=A(2)*B(3)-A(3)*B(2)
      IF(IX.EQ.1)CALL NORMAL(C,3)
      RETURN
      END
C
C CALCULATES THE DISTANCE BETWEEN TWO POINTS IN CARTESIAN SPACE
C
      function dist(a,b)
      IMPLICIT NONE
      INTEGER I
      DOUBLE PRECISION A(3),B(3),Z,DIST

      Z=0.D0
      DO I=1,3
         Z=Z+(A(I)-B(I))**2
      ENDDO
      DIST=DSQRT(Z)
      RETURN
      END

C
C     ROBUST EQUIVALENCE CHECK - DO WELL DEFINED SORT ON COORDINATE
C     MATRIX AND COMPARE ELEMENT BY ELEMENT.  SHOULD BE FOOLPROOF.
C
C     VEC      coordinate vector to be checked (modified)
C     VECR     sorted reference coordinate vector (input only)
C     NORD     ???
C     ICOMP    number of coordinates outside of TOL (output only)
C     TOL      tolerance for comparison of coords (input only)
C
      SUBROUTINE COMPARE2(VEC,VECR,NORD,ICOMP,TOL)
      USE COMMONS
      IMPLICIT NONE
C     Maximum number of atoms currently allowed
      INTEGER NORD(2*NATOMS+1),I,JAP,ICOMP
      DOUBLE PRECISION VEC(3*NATOMS),VECR(3*NATOMS),Z,TOL
C
      ICOMP=0
      IF (IPRNT.GT.130) THEN
         WRITE(6,*)'Input vector before sorting'
         WRITE(6,80)(VEC(JAP),JAP=1,3*NATOMS)
 80   FORMAT(3(1X,F10.5))
      ENDIF
      CALL SORTXYZ(VEC,VEC,NORD(NATOMS+1),TOL)
      IF (IPRNT.GT.130) THEN
         WRITE(6,*)'Comparison vector'
         WRITE(6,80)(VECR(JAP),JAP=1,3*NATOMS)
         WRITE(6,*)'Sorted input vector'
         WRITE(6,80)(VEC(JAP),JAP=1,3*NATOMS)
      ENDIF
      DO 30 I=1,NATOMS*3
         Z = DABS( VECR(I)-VEC(I) )
C        PRINT*,'I,Z,VECR,VEC=',I,Z,VECR(I),VEC(I)
         IF (Z.GT.TOL) THEN
            ICOMP = ICOMP + 1
            RETURN
         ENDIF
30    CONTINUE
      RETURN
      END

       INTEGER FUNCTION LNBLNK (STRING)
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Purpose:   Returns the position of the last non-blank character
C
C Arguments: STRING   character string (input only)
C
C Remarks:   All FORTRAN 77 character variables are blank padded on the
C            right.  The intrinsic function LEN returns the dimension
C            of the character object, not the length of the contents.
C            
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    Revision 0.0  87/07/24  bernholdt (VAX)
C
C    Revision 1.1  88/01/11  22:08:15  bernhold
C    Initial revision
C    
C
C System:     Standard FORTRAN 77
C
C Copyright 1987 David E. Bernholdt
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       INTEGER I
       CHARACTER*(*) STRING
       CHARACTER(LEN=1) BLANK
       PARAMETER (BLANK = ' ')
C
C      Start at the end and work backwards
C
       DO 100 I = LEN(STRING), 1, -1
C         Look for first non-whitespace character
          IF (STRING(I:I) .NE. BLANK) THEN
             LNBLNK = I
             RETURN
          ENDIF
  100  CONTINUE
C
C      If we get this far, the string is empty
       LNBLNK = 0
       RETURN
       END

      SUBROUTINE MTRANSOPT(A,AT,NR,NC,NTR,NTC)
      IMPLICIT NONE
      INTEGER I,J,NTC,NTR,NR,NC
      DOUBLE PRECISION A(NR,NC),AT(NC,NR)

      DO I=1,NTR
         DO J=1,NTC
            AT(J,I)=A(I,J)
         ENDDO
      ENDDO
      RETURN
      END
        function angmag(orbit,d,iorder,ireturn)
        IMPLICIT NONE
        INTEGER IRETURN, IORDER
        DOUBLE PRECISION DTOR, ORDER, A1, X, TOP, BOT, ORBIT, D, ANGMAG
        IRETURN=0
        ANGMAG=0.D0
        DTOR=DACOS(-1.D0)/180.D0
        ORDER=FLOAT(IORDER)
        A1=180.D0/ORDER
        IF(DABS(A1-90.D0).LT.1.D-3)THEN
          ANGMAG=0.D0
        RETURN
        ENDIF
        X=ORBIT**2-0.25D0*D**2/DSIN(DTOR*A1)**2
        TOP=D/DTAN(DTOR*A1)
        IF(X.LT.-1.D-12)THEN
C
C THE ROTATION IS IMPOSSIBLE
C
         IRETURN=1
         RETURN
        ENDIF
        IF(DABS(X).LT.1.D-6)THEN
          ANGMAG=90.D0
        ELSE
        BOT=2.D0*DSQRT(X)
          ANGMAG=DATAN(TOP/BOT)/DTOR
        ENDIF
        RETURN
        END         
C
C REFLECTS POINTS IN PLANE
C
      SUBROUTINE REFLECT(A,SCRATCH,NATOM,IPLANE)
      IMPLICIT NONE
      INTEGER I,J,IX,NATOM,IPLANE
      DOUBLE PRECISION A(3,NATOM),SCRATCH(3,NATOM)
C
      IX(I,J)=MIN(I,J)/MAX(I,J)
C
      DO I=1,NATOM
         DO J=1,3
            SCRATCH(J,I)=A(J,I)*(1.D0-2.D0*IX(IPLANE,J))
         ENDDO
      ENDDO
      RETURN
      END

        SUBROUTINE PIKSR2(N,ARR,NLIST)
        IMPLICIT NONE
        INTEGER N,NLIST(N),I,J,NLST
        DOUBLE PRECISION ARR(N),A

        DO 12 J=2,N
        A=ARR(J)
        NLST=NLIST(J)
        DO 11 I=J-1,1,-1
         IF(ARR(I).LE.A)GOTO 10
         ARR(I+1)=ARR(I)
         NLIST(I+1)=NLIST(I)
   11   CONTINUE
        I=0
   10   ARR(I+1)=A
        NLIST(I+1)=NLST
   12   CONTINUE
        RETURN
        END

      SUBROUTINE SCDOT(X,V,N)
      IMPLICIT NONE
      INTEGER I,N
      DOUBLE PRECISION V(N), X

      DO I=1,N
         V(I)=V(I)*X
      ENDDO
      RETURN
      END

      SUBROUTINE FILE_OPEN(FILE_NAME, FILE_UNIT, APPEND)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         LOGICAL, INTENT(IN)           :: APPEND
         INTEGER, INTENT(OUT)          :: FILE_UNIT
         INTEGER                       :: GETUNIT
         
         FILE_UNIT = GETUNIT()
         IF (APPEND) THEN
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN', POSITION='APPEND')
         ELSE
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN')
         END IF
      
      END SUBROUTINE FILE_OPEN

      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
      !
      ! start checking for available units > 103, to avoid system default units
      ! 100, 101 and 102 are stdin, stdout and stderr respectively.
      ! 
         INTEGER :: UNITNUM

         INUSE=.TRUE.
         UNITNUM=103

         DO WHILE (INUSE)
            INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
            IF (.NOT.INUSE) THEN
               GETUNIT=UNITNUM 
            ELSE     
               UNITNUM=UNITNUM+1
            ENDIF
         ENDDO
      END FUNCTION GETUNIT

      INTEGER FUNCTION FILE_LENGTH(FILE_NAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         INTEGER                       :: FILE_UNIT
         CHARACTER(LEN=80)             :: DUMMY_LINE
         INTEGER                       :: IO_STATUS

         CALL FILE_OPEN(FILE_NAME, FILE_UNIT, .FALSE.)
         FILE_LENGTH = 0
         DO
            READ(FILE_UNIT, *, IOSTAT=IO_STATUS) DUMMY_LINE
            IF (IO_STATUS /= 0) EXIT
            FILE_LENGTH = FILE_LENGTH + 1
         ENDDO

         CLOSE(FILE_UNIT)
      END FUNCTION FILE_LENGTH 

