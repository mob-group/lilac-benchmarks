!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
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
MODULE MINIMISER2
     IMPLICIT NONE
     CONTAINS

     ! This work was described in detail in: S. A. Trygubenko and D. J. Wales, `A Doubly Nudged Elastic Band Method for Finding
     ! Transition States', J. Chem. Phys., 120, 2082-2094 (2004). Summary is available online at
     ! http://www-wales.ch.cam.ac.uk/~sat39/DNEBtests/
     SUBROUTINE NEBSQVV(N)
          USE KEYMINIMIZER
          USE NEBDATA
          USE KEYSQVV
          USE KEYNEB,ONLY: NIMAGE,MOREPRINTING,DUMPNEBXYZ,DUMPNEBXYZFREQ,DUMPNEBEOS,DUMPNEBEOSFREQ,DEBUG,&
          &DUMPNEBPTS,DUMPNEBPTSFREQ
          USE GRADIENTS
          USE NEBUTILS
          USE CHARUTILS
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: N

          INTEGER(4) :: I,Q,IMAX
          DOUBLE PRECISION :: PROJ,GNORM
          DOUBLE PRECISION,PARAMETER   :: MASS=1.0D0
          DOUBLE PRECISION,ALLOCATABLE :: V(:),DX(:)
          CHARACTER(LEN=1) :: QUENCHING,QVVTYPE

          ! where to quench velocity?
          QVVTYPE = "2" ! 2 == SQVV

          ALLOCATE(V(N),DX(N))
          IF (MOREPRINTING) THEN
               PRINT '(/1x,a)', repeat(">",21)//" ENTERING SQVV MINIMISER "//repeat(">",21)
               WRITE(*,'(1x,a)') repeat("_",91)
               WRITE(*,'(a6,a20,a20,a9,2a8,2a9)') 'Iter','Energy','RMS Force','Av.Dev','Path','Step','Vel'
               WRITE(*,'(1x,a)') repeat("_",91)
          ENDIF

          IF (DEBUG) CALL DUMPFILES("b")
          
          IF (SQVVGUESS) THEN
               IMAX = NITERSQVVGUESSMAX
          ELSE
               IMAX = NITERMAX
          ENDIF

          DO NITERDONE=1,IMAX
               CALL NEBGRADIENT
               IF (QVVTYPE == "4") call quench
               IF (BADTAU) EXIT
               
               ! second velocity update
               IF (READVEL) THEN
                    PRINT *, "Velocity read from file vel."
                    ! read velocity vector from previous iteration, assume
                    ! geometry was read as well already
                    OPEN(UNIT=333,FILE="vel",status="old",form="unformatted")
                    READ(UNIT=333) V
                    CLOSE(333)
                    V = V - (0.5D0*DT/MASS)*G
               ELSE
                    IF (NITERDONE==1) THEN
                         V = 0.0D0
                    ELSE
                         V = V - (0.5D0*DT/MASS)*G
                    ENDIF
               ENDIF

                    IF (QVVTYPE == "1") call quench
                    
               ! calculate the step and limit it if appropriate
               DX = V*DT - (0.5D0*DT**2/MASS)*G ! ZSUNYLYSYA PO KOORDYNATAH; V CHASI - NA DT
               IF (LIMITSQVVSTEP) WHERE (ABS(DX) > STEPDOFMAX) DX = STEPDOFMAX*DX/ABS(DX)
               
               ! calculate the step per image using dX
               DO I=1,NIMAGE
                    STEPIMAGE(I) = SQRT(DOT_PRODUCT(DX(NOPT*(I-1)+1:NOPT*I),DX(NOPT*(I-1)+1:NOPT*I)))
               ENDDO
               STEPTOT = SUM(STEPIMAGE)/NIMAGE
               
               X = X + DX ! COORDINATE UPDATE
                    IF (QVVTYPE == "2") call quench
               V =     V    - (0.5D0*DT  / MASS)*G ! FIRST VELOCITY UPDATE
                    IF (QVVTYPE == "3") call quench
               
               CALL IMAGEDISTRIBUTION
               IF (MOREPRINTING) WRITE(*,'(i6,2f20.10,f8.2,a,f8.3,2f9.3)') &
               & NITERDONE,ETOTAL/NIMAGE,RMS,AVDEV,'%',Separation,StepTot,Sqrt(Dot_Product(V,V)/(Nimage*Nopt))

               IF (SQVVGUESS) THEN
                    IF (RMS < SQVVGUESSRMSTOL .OR. NITERDONE > NITERSQVVGUESSMAX) EXIT
               ELSEIF ( RMS<=RMSTOL.AND.NITERDONE>1.AND.NITERDONE>NITERMIN ) THEN
                    EXITSTATUS=1
                    EXIT
               ENDIF

               IF (DEBUG) THEN
                    IF (MOD(NITERDONE,10)==0) CALL DUMPFILES("m")
                    !if (TwoD) call writexyztmp
               ENDIF
               IF (DUMPNEBXYZ.AND.MOD(NITERDONE,DUMPNEBXYZFREQ)==0) CALL RWG("w",.False.,NITERDONE)
               IF (DUMPNEBPTS.AND.MOD(NITERDONE,DUMPNEBPTSFREQ)==0) CALL SAVEBANDCOORD
               IF (DUMPNEBEOS.AND.MOD(NITERDONE,DUMPNEBEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE)

               ! we need to dump velocity vector
               OPEN(UNIT=333,FILE="vel",status="replace",form="unformatted")
               WRITE(UNIT=333) V
               CLOSE(333)
          ENDDO
          NITERDONE = NITERDONE - 1

          IF (NITERDONE == NITERMAX) EXITSTATUS=2
          IF (SQVVGUESS) THEN
               NITERDONESAVE = NITERDONE
               INTSTR = WI(NITERDONE)
               PRINT *, 'Initial relaxation: '//trim(IntStr)//' SQVV steps.'
          ENDIF
          IF (DEBUG.AND..NOT.SQVVGUESS) CALL DUMPFILES("e")
          DEALLOCATE(V,DX)

          CONTAINS

          SUBROUTINE QUENCH
               QUENCHING = "j"
               SELECT CASE(QUENCHING)
               CASE ("j") ! Quenching by Jonsson et al.
                    DO I=1,NIMAGE
                         PROJ = DOT_PRODUCT(V(NOPT*(I-1)+1:NOPT*I),G(NOPT*(I-1)+1:NOPT*I)) ! NE ZOVISIM PROJ BO G NE NORMOVANY
                         IF (PROJ > 0.0D0) THEN
                              !print '(1x,a,i3,a)', 'For image',i,' velocity was zeroed'
                              V(NOPT*(I-1)+1:NOPT*I) = 0.0D0
                         ELSE
                              GNORM = DOT_PRODUCT(G(NOPT*(I-1)+1:NOPT*I),G(NOPT*(I-1)+1:NOPT*I)) ! NORMA V KVADRATI (DYV *)! 
                              V(NOPT*(I-1)+1:NOPT*I) = -(PROJ/GNORM)*G(NOPT*(I-1)+1:NOPT*I)
                         ENDIF
                    ENDDO
               CASE ("c") ! Crehuet and Field quenching
                    Q = 3 ! ROZMIR KOMIRKY :)
                    DO I=1,NIMAGE*NOPT-(Q-1),Q
                         PROJ = DOT_PRODUCT(V(I:I+Q-1),G(I:I+Q-1))
                         IF (PROJ > 0.0D0) THEN
                              V(I:I+Q-1) = 0.0D0
                         ELSE
                              GNORM = DOT_PRODUCT(G(I:I+Q-1),G(I:I+Q-1))
                              V(I:I+Q-1) = -(PROJ/GNORM)*G(I:I+Q-1)
                         ENDIF
                    ENDDO
               CASE ("g") ! neighbours gradients are used in projections; this is for 2D only
                    Q = 2
                    DO I=1,NIMAGE*NOPT-(2*Q-1),Q
                         PROJ = DOT_PRODUCT(V(I:I+Q-1),G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1))
                         IF (PROJ > 0.0D0) THEN
                              V(I:I+Q-1) = 0.0D0
                         ELSE
                              GNORM = DOT_PRODUCT(G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1),G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1))
                              V(I:I+Q-1) = -(PROJ/GNORM)*(G(I:I+Q-1)+G(I-Q:I-1)+G(I+Q:I+2*Q-1))
                         ENDIF
                    ENDDO
               END SELECT
          END SUBROUTINE QUENCH
     END SUBROUTINE NEBSQVV
     
END MODULE MINIMISER2
