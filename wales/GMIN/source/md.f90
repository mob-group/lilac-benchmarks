!  GMIN: A program for finding global minima
!  Copyright (C) 1999-2006 David J. Wales
!  This file is part of GMIN.
!
!  GMIN is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  GMIN is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
!
MODULE MOLECULAR_DYNAMICS
  !
  USE COMMONS, ONLY : NATOMS, MYUNIT, ATMASS, PRTFRQ, DUMPXYZUNIT, &
       TSTART
  ! 
  IMPLICIT NONE
  !
  ! Parameters set in the 'keywords' routine are to be saved!
  LOGICAL, SAVE :: MDT
  INTEGER, SAVE :: MD_NWAIT, MD_NFREQ, MD_UNIT, MD_NSTEPS
  DOUBLE PRECISION, SAVE :: MD_GAMMA, MD_TSTEP
  !
  ! MD-specific variables
  DOUBLE PRECISION :: EPOT, EKIN
  DOUBLE PRECISION, ALLOCATABLE :: ACC(:), VEL(:) 
  !
  DOUBLE PRECISION, PARAMETER :: PI=3.141592652D0
  !
CONTAINS
  !
  !============================================================
  !
  SUBROUTINE MDRUN(POS,NSTEPS,TEMP,DUMP)
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DUMP
    INTEGER, INTENT(IN) :: NSTEPS
    DOUBLE PRECISION, INTENT(IN) :: TEMP ! k_B*T in units of energy
    DOUBLE PRECISION, INTENT(INOUT) :: POS(3*NATOMS)
    !
    INTEGER :: I, J, K
    DOUBLE PRECISION :: TIME
    !
    IF (DUMP) WRITE(MYUNIT, '(A)') "mdrun> Start of MD run."
       
    !
    ALLOCATE(ACC(3*NATOMS), VEL(3*NATOMS))
    !
    CALL INIVEL(POS, TEMP) 
    CALL POTENTIAL(POS,ACC,EPOT,.TRUE.,.FALSE.) ! get gradient
    K=0
    DO I=1,NATOMS
       DO J=1,3
          K=K+1
          ACC(K) = -ACC(K)/ATMASS(I) ! get ACC
       ENDDO
    ENDDO
    !
    I=0
    IF (DUMP) WRITE(MYUNIT, '(A,I8,A,F18.8,A,F14.8)') &
         'mdrun> Step ', I,' EPOT= ',EPOT,' EKIN= ',EKIN
    !
    DO I=1,NSTEPS
       !
       CALL MDSTEP(POS, TEMP)
       !
       IF(DUMP) THEN ! Dump data to files
          !
          IF(MOD(I,PRTFRQ).EQ.0) THEN ! Basic data
             WRITE(MYUNIT, '(A,I8,A,F18.8,A,F14.8)') &
                  'mdrun> Step ', I,' EPOT= ',EPOT,' EKIN= ',EKIN
          ENDIF
          !
          IF(I == MD_NWAIT) THEN
             CALL MYCPU_TIME(TIME)
             WRITE(MYUNIT, '(A,F11.1)') &
                  "mdrun> Finished MD equilibration t= ",TIME-TSTART
          ELSEIF(I > MD_NWAIT .AND. &
               MOD(I-MD_NWAIT,MD_NFREQ)==0) THEN ! XYZ dump
             WRITE(DUMPXYZUNIT(1),'(I4)') NATOMS
             WRITE(DUMPXYZUNIT(1),'(A,I9,A,F20.10)')  &
                  'MD step ',I,'  EPOT= ', EPOT
             DO J=1,NATOMS
                WRITE(DUMPXYZUNIT(1),'(A,3(1X,F10.5))') &
                     'X ', (POS(3*(J-1)+K), K=1,3)
             ENDDO             
          ENDIF
          !
       ENDIF
       !
    ENDDO
    !
    DEALLOCATE(ACC, VEL)
    !
    CALL MYCPU_TIME(TIME)
    IF (DUMP) WRITE(MYUNIT, '(A,F11.1)') &
         "mdrun> Finished MD run t= ", TIME-TSTART
    !
    RETURN
    !
  END SUBROUTINE MDRUN
  !
  !============================================================
  !
  SUBROUTINE MDSTEP(POS, TEMP)
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT(IN) :: TEMP
    DOUBLE PRECISION, INTENT(INOUT) :: POS(3*NATOMS)
    
    INTEGER :: I,J,K
    DOUBLE PRECISION :: HDT, GFRIC, NOISE(NATOMS), &
         DPRAND,R1,R2,NR,NR2(3*NATOMS)
    !
    HDT=0.5D0*MD_TSTEP
    GFRIC=1.0D0-MD_GAMMA*HDT
    NOISE(:) = DSQRT(TEMP*MD_GAMMA*MD_TSTEP/ATMASS(:))
    !
    K=0
    DO I=1,NATOMS
       DO J=1,3
          K=K+1
          ! Generate normal random numbers
          R1=DSQRT(-2.0D0*LOG(DPRAND()))
          R2=2.0*PI*DPRAND()
          NR=R1*COS(R2)
          NR2(K)=R1*SIN(R2)          
          ! Update mid-step velocity
          VEL(K) = GFRIC*VEL(K) + ACC(K)*HDT + NR*NOISE(I)
          ! Update full-step position
          POS(K) = POS(K) + VEL(K)*MD_TSTEP
       ENDDO
    ENDDO
    !
    CALL POTENTIAL(POS,ACC,EPOT,.TRUE.,.FALSE.) ! GET GRADIENT
    !
    K=0
    EKIN=0.0D0
    DO I=1,NATOMS
       DO J=1,3
          K=K+1
          ACC(K)= -ACC(K)/ATMASS(I) ! get ACC
          ! Update full-step velocity
          VEL(K) = GFRIC*VEL(K) + ACC(K)*HDT + NR2(K)*NOISE(I)
          EKIN = EKIN + ATMASS(I)*VEL(K)*VEL(K)
       ENDDO
    ENDDO
    EKIN=0.5D0*EKIN
    !
    RETURN
    !
  END SUBROUTINE MDSTEP
  !
  !============================================================
  !
  SUBROUTINE INIVEL(POS, TEMP)
    !
    IMPLICIT NONE
    !
    DOUBLE PRECISION, INTENT (IN) :: TEMP, POS(3*NATOMS)
    !
    INTEGER :: I, J, K, L, IPIVOT(3)
    DOUBLE PRECISION :: COM(3), COM_MOM(3), ANG_MOM(3), &
         MOI(3,3), MOI_INV(3,3), DPRAND, SCALE, TMASS, &
         P(3), R(3), W(3), R2, X(NATOMS,3)
    !
    COM(:) = 0.0D0
    COM_MOM(:) = 0.0D0
    TMASS = 0.0D0
    !
    ! On first pass randomly generate velocities
    K=0    
    DO I=1,NATOMS
       TMASS = TMASS + ATMASS(I) 
       DO J=1,3
          K=K+1
          ! Pick random number from normal distribution with mean of
          ! of zero and unit variance (See Allen and Tildesley p347)
          VEL(K) = DSQRT(-2.0D0*LOG(DBLE(DPRAND())))* &
               COS(2.0D0*PI*DBLE(DPRAND()))
          !
          ! Will need to subtract net translation and rotation!
          COM(J) = COM(J) + POS(K)*ATMASS(I)
          COM_MOM(J) = COM_MOM(J) + VEL(K)*ATMASS(I)
          !
       ENDDO
    ENDDO
    COM(:) = COM(:)/TMASS
    !
    WRITE(MYUNIT,'(A,3(1x,F12.6))') &
         'inivel> Centre of mass:', COM
    WRITE(MYUNIT,'(A,3(1x,F12.6))') &
         'inivel> Initial COM_MOM:', COM_MOM
    ! On second pass substract net translation and determine rotation
    ANG_MOM(:) = 0.0D0 ! L = Iw
    MOI(:,:) = 0.0D0
    K=0    
    DO I=1,NATOMS
       R2=0.0D0
       DO J=1,3
          K=K+1
          VEL(K) = VEL(K) - COM_MOM(J)/TMASS
          R(J) = POS(K) - COM(J)  ! pos rel. to CoM
          X(I,J) = R(J)
          P(J) = VEL(K)*ATMASS(I) ! momentum component
          R2 = R2 + R(J)*R(J)
       ENDDO
       ANG_MOM = ANG_MOM + CROSS(R,P)
       DO J=1,3
          MOI(J,J) = MOI(J,J) + (R2 - R(J)*R(J))*ATMASS(I) ! diagonal
          DO L=J+1,3 ! off-diagonal
             MOI(J,L) = MOI(J,L) - R(J)*R(L)*ATMASS(I)
             MOI(L,J) = MOI(J,L) ! impose symmetry
          ENDDO
       ENDDO
    ENDDO
    !
    WRITE(MYUNIT,'(A,3(1X,F12.6))') &
         'inivel> Initial ANG_MOM=', ANG_MOM
    !
    ! Invert MOI using routines from LAPACK 
    CALL DGETRF(3,3,MOI,3,IPIVOT,I) ! MOI is modified!!!   
    CALL DGETRI(3,MOI,3,IPIVOT,R,3,I)
    !
    IF(I /= 0) THEN
       WRITE(MYUNIT,'(A,I1)') &
            'md> Failed to invert inertia tensor! Error:', I
       STOP
    ENDIF
    !
    ! Compute net angular velocity and store in W
    DO I=1,3
       W(I)=DOT_PRODUCT(MOI(I,1:3), ANG_MOM(1:3)) 
    ENDDO
    !
    ! On third pass eliminate rotation and compute kinetic energy
    EKIN=0
    K=0
    ANG_MOM(:) = 0.0D0 ! L = Iw
    COM_MOM(:) = 0.0D0
    DO I=1,NATOMS
       R=CROSS(W,X(I,1:3))
       DO J=1,3
          K=K+1
          VEL(K) = VEL(K) - R(J)
          EKIN = EKIN + ATMASS(I)*VEL(K)*VEL(K)
          P(J) = VEL(K)*ATMASS(I) ! momentum component
          COM_MOM(J) = COM_MOM(J) + VEL(K)*ATMASS(I)
       ENDDO
       ANG_MOM = ANG_MOM + CROSS(X(I,1:3),P)
       !write(myunit, *) ang_mom
    ENDDO
    !
    WRITE(MYUNIT,'(A,3(1x,F12.6))') &
         'inivel> Final COM_MOM=', COM_MOM
    WRITE(MYUNIT,'(A,3(1x,F12.6))') &
         'inivel> Final ANG_MOM=', ANG_MOM
    !
    SCALE=DSQRT(DBLE(K)*TEMP/EKIN)
    !
    ! On fourth pass rescale velocities to correct temperature
    K=0 
    EKIN=0.D0
    DO I=1,NATOMS
       DO J=1,3
          K=K+1
          VEL(K) = SCALE*VEL(K)
          EKIN = EKIN + ATMASS(I)*VEL(K)*VEL(K)
       ENDDO
    ENDDO
    EKIN=0.5D0*EKIN
    !
    RETURN
    !
  END SUBROUTINE INIVEL
  !
  !============================================================
  !
  FUNCTION CROSS(A,B)
    !
    DOUBLE PRECISION, DIMENSION(3) :: CROSS
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: A,B
    !
    CROSS(1) = A(2) * B(3) - A(3) * B(2)
    CROSS(2) = A(3) * B(1) - A(1) * B(3)
    CROSS(3) = A(1) * B(2) - A(2) * B(1)
    !
  END FUNCTION CROSS
  !
END MODULE MOLECULAR_DYNAMICS
