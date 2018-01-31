!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

C
C     This subprogram performs a sort on the input data and
C     arranges it from biggest to smallest. The exchange-sort
C     algorithm is used.
C
      SUBROUTINE SORT(N,J3,A,NA)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2, NA(J3), NTEMP
      DOUBLE PRECISION TEMP, A(J3)
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).LT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
20    CONTINUE
!     WRITE(*,'(A,G20.10,A,I6)') 'In sort the minimum saved rate is ',A(N-1),' saved at position ',NA(N-1)
!     WRITE(*,'(I6,G20.10,I6)') (J1,A(J1),NA(J1),J1=1,N)
      RETURN
      END
C
C     This subprogram performs a sort on the input data and
C     arranges it from biggest to smallest. The exchange-sort
C     algorithm is used.
C
      SUBROUTINE SORT2(N,J3,A,NA)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2, NA(J3), NTEMP, A(J3)
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).LT.A(J2)) L=J2
10       CONTINUE
         NTEMP=A(L)
         A(L)=A(J1)
         A(J1)=NTEMP
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
20    CONTINUE
      RETURN
      END
C
C     This subprogram performs a sort on the input data and
C     arranges it from smallest to largest. The exchange-sort
C     algorithm is used. In this case we sort on NA
C
      SUBROUTINE SORT4(N,J3,A,NA)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2, NA(J3), NTEMP
      DOUBLE PRECISION TEMP, A(J3)
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (NA(L).GT.NA(J2)) L=J2
10       CONTINUE
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
20    CONTINUE
C     WRITE(*,'(A,G20.10,A,I6)') 'In sort the minimum saved rate is ',A(N-1),' saved at position ',NA(N-1)
C     WRITE(*,'(I6,G20.10,I6)') (J1,A(J1),NA(J1),J1=1,N)
      RETURN
      END

C
C     This subprogram performs a sort on the input data and
C     arranges it from smallest to largest. The exchange-sort
C     algorithm is used. In this case we sort on A, rather than NA
C
      SUBROUTINE SORT5(N,J3,A,NA)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2, NA(J3), NTEMP
      DOUBLE PRECISION TEMP, A(J3)
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
10       CONTINUE
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
20    CONTINUE
C     WRITE(*,'(A,G20.10,A,I6)') 'In sort the minimum saved rate is ',A(N-1),' saved at position ',NA(N-1)
C     WRITE(*,'(I6,G20.10,I6)') (J1,A(J1),NA(J1),J1=1,N)
      RETURN
      END

