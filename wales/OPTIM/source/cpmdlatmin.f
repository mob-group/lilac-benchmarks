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
      SUBROUTINE CPMDLATMIN(N,XC,F1,V,BOXLX)
      USE KEY
      USE COMMONS,ONLY : DEBUG
      use porfuncs
      IMPLICIT NONE
      INTEGER N, NCOUNT, J1
      DOUBLE PRECISION XC(3*N),F1,F2,F3,GRAD,SECOND,XSAVE(3*N), V(3*N), BSAVE, BOXLX,
     1                 DIF, TEMP1, GEMAX, THRESH, EPREV, THRESHMAX
      LOGICAL YESNO
      CHARACTER BOXSTRING*80
C
C  Value of DIF is the order of magnitude to which the lattice
C  constant can be optimised. Setting it smaller than 10^(-7)
C  causes numerical problems on the DEC.
C

      THRESH=0.08D0
      THRESHMAX=0.08D0
      EPREV=0.0D0
      INQUIRE(FILE='RESTART.1',EXIST=YESNO)

      DIF=1.0D-3
      NCOUNT=1
      DO J1=1,3*N
         XSAVE(J1)=XC(J1)
      ENDDO

      BSAVE=BOXLX

10    TEMP1=BOXLX

      CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out >& /dev/null ')
      IF (.NOT.YESNO) THEN
         IF (PARALLEL) THEN
            CALL SYSTEM(' ( mpirun -np ' // NPROC // ' /export/home/wales/bin/cpmd.x.mpi ' 
     1                    // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
         ELSE
            CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
         ENDIF
      ELSE
         OPEN(UNIT=8,FILE='newgeom',STATUS='UNKNOWN')
         DO J1=1,N
            WRITE(8,'(6F20.10)') XC(3*(J1-1)+1),XC(3*(J1-1)+2),XC(3*(J1-1)+3),0.0D0,0.0D0,0.0D0
         ENDDO
         CLOSE(8)
         CALL SYSTEM(' mv newgeom GEOMETRY ')
         WRITE(BOXSTRING,'(F20.10,5F12.4)') BOXLX, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0
         CALL SYSTEM(' sed -e "s/CELLSIZE/' // BOXSTRING // '/" ' // SYS(1:LSYS) 
     &                    // '.restart > '  // SYS(1:LSYS) // '.latmin ')
         IF (PARALLEL) THEN
            CALL SYSTEM(' ( mpirun -np ' // NPROC // ' /export/home/wales/bin/cpmd.x.mpi ' 
     1                  // SYS(1:LSYS) // '.latmin > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
         ELSE
            CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // '.latmin > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
         ENDIF
      ENDIF

      YESNO=.TRUE.

      OPEN (UNIT=7,FILE='ENERGY',STATUS='OLD')
      READ(7,*) F1, GEMAX
      IF (DEBUG) WRITE(*,'(A,3F20.10)') 'ENERGY,GEMAX,BOXLX=',F1,GEMAX,BOXLX
C     IF (GEMAX.GT.1.0D-5) THEN
C        WRITE(*,'(A,G15.5,A)') 'WARNING, GEMAX=',GEMAX,' CPMD wavefunction convergence suspect'
C     ENDIF
      CLOSE(7)
      CALL SYSTEM('grep "CPU TIME" ' // SYS(1:LSYS) // 
     1             '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" > temp')
C    1             '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" -e "/s/  / /g" > temp')
      OPEN (UNIT=7,FILE='temp',STATUS='OLD')
      READ(7,'(A)') BOXSTRING
      WRITE(*,'(A,A,F20.10,A,F20.10)') BOXSTRING,' Energy=',F1,' GEMAX=',GEMAX
      CLOSE(7)
      OPEN(UNIT=7,FILE='GEOMETRY',STATUS='OLD')
      DO J1=1,N
         READ(7,*) GEMAX,GEMAX,GEMAX,V(3*(J1-1)+1),V(3*(J1-1)+2),V(3*(J1-1)+3)
         V(3*(J1-1)+1)=-V(3*(J1-1)+1)
         V(3*(J1-1)+2)=-V(3*(J1-1)+2)
         V(3*(J1-1)+3)=-V(3*(J1-1)+3)
      ENDDO
      CLOSE(7)

      IF (EPREV.NE.0.0D0) THEN
         IF (F1.GT.EPREV) THEN
            THRESH=THRESH/1.2D0
         ELSE
            THRESH=MIN(THRESHMAX,THRESH*1.1D0)
         ENDIF
      ENDIF

      BOXLX=TEMP1+DIF
      DO J1=1,3*N
         XC(J1)=XSAVE(J1)*BOXLX/BSAVE
      ENDDO
      OPEN(UNIT=8,FILE='newgeom',STATUS='UNKNOWN')
      DO J1=1,N
         WRITE(8,'(6F20.10)') XC(3*(J1-1)+1),XC(3*(J1-1)+2),XC(3*(J1-1)+3),0.0D0,0.0D0,0.0D0
      ENDDO
      CLOSE(8)
      CALL SYSTEM(' mv newgeom GEOMETRY ')
      WRITE(BOXSTRING,'(F20.10,5F12.4)') BOXLX, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0
      CALL SYSTEM(' sed -e "s/CELLSIZE/' // BOXSTRING // '/" ' // SYS(1:LSYS) // '.restart > '  // SYS(1:LSYS) // '.latmin ')
      IF (PARALLEL) THEN
         CALL SYSTEM(' ( mpirun -np ' // NPROC // ' /export/home/wales/bin/cpmd.x.mpi ' 
     1                // SYS(1:LSYS) // '.latmin > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
      ELSE
         CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // '.latmin > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
      ENDIF
      OPEN (UNIT=7,FILE='ENERGY',STATUS='OLD')
      READ(7,*) F2, GEMAX
      IF (DEBUG) WRITE(*,'(A,3F20.10)') 'ENERGY,GEMAX,BOXLX=',F2,GEMAX,BOXLX
C     IF (GEMAX.GT.1.0D-5) THEN
C        WRITE(*,'(A,G15.5,A)') 'WARNING, GEMAX=',GEMAX,' CPMD wavefunction convergence suspect'
C     ENDIF
      CLOSE(7)
      CALL SYSTEM('grep "CPU TIME" ' // SYS(1:LSYS) // 
     1             '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" > temp')
C    1             '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" -e "/s/  / /g" > temp')
      OPEN (UNIT=7,FILE='temp',STATUS='OLD')
      READ(7,'(A)') BOXSTRING
      WRITE(*,'(A,A,F20.10,A,F20.10)') BOXSTRING,' Energy=',F2,' GEMAX=',GEMAX
      CLOSE(7)

      BOXLX=TEMP1-DIF
      DO J1=1,3*N
         XC(J1)=XSAVE(J1)*BOXLX/BSAVE
      ENDDO
      OPEN(UNIT=8,FILE='newgeom',STATUS='UNKNOWN')
      DO J1=1,N
         WRITE(8,'(6F20.10)') XC(3*(J1-1)+1),XC(3*(J1-1)+2),XC(3*(J1-1)+3),0.0D0,0.0D0,0.0D0
      ENDDO
      CLOSE(8)
      CALL SYSTEM(' mv newgeom GEOMETRY ')
      WRITE(BOXSTRING,'(F20.10,5F12.4)') BOXLX, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0
      CALL SYSTEM(' sed -e "s/CELLSIZE/' // BOXSTRING // '/" ' // SYS(1:LSYS) // '.restart > '  // SYS(1:LSYS) // '.latmin ')
      IF (PARALLEL) THEN
         CALL SYSTEM(' ( mpirun -np ' // NPROC // ' /export/home/wales/bin/cpmd.x.mpi ' 
     1               // SYS(1:LSYS) // '.latmin > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
      ELSE
         CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // '.latmin > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
      ENDIF
      OPEN (UNIT=7,FILE='ENERGY',STATUS='OLD')
      READ(7,*) F3, GEMAX
      IF (DEBUG) WRITE(*,'(A,3F20.10)') 'ENERGY,GEMAX,BOXLX=',F3,GEMAX,BOXLX
C     IF (GEMAX.GT.1.0D-5) THEN
C        WRITE(*,'(A,G15.5,A)') 'WARNING, GEMAX=',GEMAX,' CPMD wavefunction convergence suspect'
C     ENDIF
      CLOSE(7)
      CALL SYSTEM('grep "CPU TIME" ' // SYS(1:LSYS) // 
     1             '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" > temp')
C    1             '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" -e "/s/  / /g" > temp')
      OPEN (UNIT=7,FILE='temp',STATUS='OLD')
      READ(7,'(A)') BOXSTRING
      WRITE(*,'(A,A,F20.10,A,F20.10)') BOXSTRING,' Energy=',F3,' GEMAX=',GEMAX
      CLOSE(7)

      GRAD=(F2-F3)/(2.0D0*DIF)
      
      IF (GRAD.EQ.0.0D0) STOP

      SECOND=(F3+F2-2.0D0*F1)/(DIF*DIF)
      WRITE(*,'(A,I5,A,T50,F20.10)') 'Energy for lattice cycle ',NCOUNT,' is ',F1
      WRITE(*,'(A,T50,F20.10)') 'Gradient wrt box length=',GRAD
      WRITE(*,'(A,T50,F20.10)') 'Second derivative wrt box length=',SECOND
      WRITE(*,'(A,T50,F20.10)') 'Full step size=',-GRAD/SECOND
      EPREV=F1
      NCOUNT=NCOUNT+1
      IF (DABS(GRAD/SECOND).GT.THRESH) THEN
         BOXLX=TEMP1-THRESH*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         WRITE(*,'(A,T50,F20.10,A,F15.5)') 'Scaled step=',-THRESH*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD)),' new box length=',BOXLX
         DO J1=1,3*N
            XC(J1)=XSAVE(J1)*BOXLX/BSAVE
         ENDDO
         GOTO 10
      ELSE
         BOXLX=TEMP1-GRAD/SECOND
         WRITE(*,'(A,T50,F20.10,A,F15.5)') 'Scaled step=',-GRAD/SECOND,' new box length=',BOXLX
         DO J1=1,3*N
            XC(J1)=XSAVE(J1)*BOXLX/BSAVE
         ENDDO
         IF (DABS(GRAD/SECOND).GT.DIF/1.0D1) GOTO 10
      ENDIF

      WRITE(*,'(A,I6,A)') 'Box length optimisation converged in ',NCOUNT,' steps. Coordinates written to points.latmin'
      DO J1=1,3*N
         XC(J1)=XSAVE(J1)*BOXLX/BSAVE
      ENDDO
      IF (FILTH.NE.0) PRINT*,'*** WARNING FILTH not equal to zero in cpmdlatmin'
      OPEN(UNIT=7,FILE='points.latmin',STATUS='UNKNOWN')
      WRITE(7,'(3F20.10)') (XC(J1),J1=1,3*N)
      STOP

      RETURN
      END
