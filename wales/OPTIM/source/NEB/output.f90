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
MODULE NEBOUTPUT
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE TSLOCATOR(TSRESET)
          USE KEY,ONLY: BFGSTST,UNRST,NSTEPS,MACHINE, GROWSTRINGT, INTEPSILON, REDOTSIM, EDIFFTOL, &
  &                     NONEBMAX, INTCONSTRAINTT, CUDAT, FLATPATHT, FLATTESTT, FLATEDIFF, MAXIMFACTOR
          USE GSDATA, ONLY: EVOLVESTRINGT
          USE KEYOUTPUT
          USE MODCHARMM
          USE NEBDATA
          USE KEYNEB,ONLY:NIMAGE,DEBUG
          USE NEBTOCONNECT
          USE CHARUTILS
          USE MODGUESS
          USE MODMEC
          USE LINKEDLIST
          USE MODEFOL
          USE INTCOMMONS, ONLY : DESMINT, NINTC, NNZ, KD, INTNEWT
          USE COMMONS, ONLY : REDOPATH, REDOPATHNEB
          USE MODCUDABFGSTS, ONLY : CUDA_BFGSTS_WRAPPER
          USE GENRIGID, ONLY: DEGFREEDOMS, RIGIDINIT

          IMPLICIT NONE
          
          INTEGER :: I,J,NT,NM,ITDONE=0,J1,RECLEN
          INTEGER,PARAMETER :: MAXPRINTOUT = 50, ITMAX  = 30
          DOUBLE PRECISION :: EDUMMY,EVALMIN,EVALMAX,MAXE,VECSNORM
          LOGICAL :: TSCONVERGED,T,TSRESET
          DOUBLE PRECISION,DIMENSION(NOPT) :: LGDUMMY, VECS, DIAG, THISXYZ
          INTEGER :: MLOC
          DOUBLE PRECISION :: TIME, TIME0
          DOUBLE PRECISION :: DPRAND
          LOGICAL :: KNOWE, KNOWG, KNOWH ! JMC
          COMMON /KNOWN/ KNOWE, KNOWG, KNOWH ! JMC
          CHARACTER(LEN=256) :: FILENAME, METHSTR
          INTEGER TSPOS(NIMAGE+2)!, MINPOS(NIMAGE+2)
          INTEGER BIGN
          !DOUBLE PRECISION :: SCA
          LOGICAL :: TMPINTNEWT, FAILED
          NT = 0
          VECS(:) = 0 ! sn402: to avoid uninitialised value problems

          IF (REDOPATHNEB) THEN
             NT=1
             MAXE=-1.0D100
             MLOC=REDOTSIM+1
             PRINT '(A,F20.10)',' tslocator> transition state has energy ',EEE(REDOTSIM+1)
             TSPOS(NT)=MLOC
          ELSE
             ! IDENTIFY TS CANDIDATES
             SELECT CASE(CANDIDATES) ! Unless the keyword 'CANDIDATES' was used, this should
             ! always be "maxim"
             CASE('high')
                  NT=1
                  MAXE=-1.0D100
                  DO J1=2,NIMAGE+1
                     IF (EEE(J1).GT.MAXE) THEN
                        MLOC=J1
                        MAXE=EEE(J1)
                     ENDIF
                  ENDDO
                  TSPOS(1)=MLOC
             CASE('all','maxim')
                  DO I=2,NIMAGE+1 
                       T=.FALSE.
                       IF (CANDIDATES=='maxim') then
                            IF ( EEE(I-1)+EDIFFTOL*MAXIMFACTOR < EEE(I) .AND. EEE(I) > EEE(I+1)+EDIFFTOL*MAXIMFACTOR ) THEN
                                 T=.TRUE.
                            ELSE
                                 T=.FALSE.
                            ENDIF
                       ENDIF
!                      PRINT '(A,I6,3F20.10,L5)','I,EEE(I-1),EEE(I),EEE(I+1),T=',I,EEE(I-1),EEE(I),EEE(I+1),T
                       IF (T) THEN
                            NT=NT+1
                            TSPOS(NT)=I ! IS A POSITION OF A MAXIMUM IN ARRAY XYZ
                       ENDIF
                  ENDDO
             END SELECT
          ENDIF

          NONEBMAX=.FALSE.
          IF (NT.EQ.0) THEN
!
! This should cope with highly asymmetric profiles, which otherwise looks monotonic
! until we try a huge number of images. 
!
             PRINT '(1x,a)', 'No maximum in profile - using highest image'
             NT=1
             IF (EEE(2).GT.EEE(NIMAGE+1)) THEN
                TSPOS(1)=2
             ELSE
                TSPOS(1)=NIMAGE+1
             ENDIF
          ENDIF

          WRITE(*,'(1X,A,I4,A)',advance='No') 'Following ',NT,' images are candidates for TS:'
          !open(UNIT=323, FILE="energydiff")
          DO J=1,NT
             WRITE(*,'(i5)',advance='No') TSPOS(J)-1
          !   write(323,*) EEE(TSPOS(J))-EEE(1), EEE(TSPOS(J))-EEE(NIMAGE+2)
          ENDDO
          PRINT *,' '

!
!sy349: This block is to test whether it is a flat path after the dneb calculation. If
!the energy difference between highest-energy image and two end points is
!smaller than FLATEDIFF. If it is a flat path, it will regard the highest-energy image directly
!as transition state, and these two end points will be regarded connected.
!
          IF (FLATTESTT) THEN

             FLATPATHT=.FALSE.
             BIGN=TSPOS(1)
             DO J=2, NT                
                IF (EEE(BIGN)<EEE(TSPOS(J))) BIGN=TSPOS(J)
             ENDDO

             IF (ABS(EEE(1)-EEE(NIMAGE+2)).LT.EDIFFTOL) THEN
                IF (ABS(EEE(BIGN)-EEE(1)).LT.FLATEDIFF &
                   &.OR.ABS(EEE(BIGN)-EEE(NIMAGE+2)).LT.FLATEDIFF) THEN
                   FLATPATHT=.TRUE.
                   IF (TSRESET) NTSFOUND=0
                   NTSFOUND=NTSFOUND+1
                   ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),&
                           &TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
                   TSFOUND(NTSFOUND)%VECS=VECS(1:NOPT)
                   TSFOUND(NTSFOUND)%COORD=XYZ(NOPT*(BIGN-1)+1:NOPT*BIGN)
                   TSFOUND(NTSFOUND)%E=EEE(BIGN)
                   TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
                ENDIF
             ENDIF

             IF (DEBUG .AND. FLATPATHT) THEN
                 write(*,*) "Index of first candidate is ", TSPOS(1)  ! sn402
                 write(*,*) "Image energies are", EEE(:) !sn402
             ENDIF
          ENDIF

          IF (OPTIMIZETS.AND. .NOT.FLATPATHT) THEN
             IF (DEBUG) THEN
                 write(*,*) "Index of first candidate is ", TSPOS(1)  ! sn402
                 write(*,*) "Image energies are", EEE(:) !sn402
             ENDIF
             IF (TSRESET) NTSFOUND=0
             CALL MYCPU_TIME(STARTTIME,.FALSE.)
             DO J=1,NT
                !IF (FLATPATHT(J)) CYCLE
                CALL MYCPU_TIME(TIME0,.FALSE.)
                EDUMMY=EEE(TSPOS(J))
                LGDUMMY(1:NOPT)=TRUEGRAD((TSPOS(J)-1)*NOPT+1:TSPOS(J)*NOPT)
                KNOWE=.TRUE.
                KNOWG=.TRUE.
                IF (REDOPATH) THEN
                   KNOWG = .FALSE.
                   KNOWE = .FALSE.
                ENDIF
                IF (BFGSTST) THEN
                   IF (UNRST) THEN ! JMC
                      KNOWG=.FALSE. ! Is this needed now that gdummy is set? DJW
                      VECS(1:NINTS)=TANVEC(1:NINTS,TSPOS(J)-1)
                      VECSNORM=SUM(VECS(1:NINTS)**2)
                      IF (VECSNORM.EQ.0.0D0) THEN  ! Just in case TANVEC is somehow not set? e.g. for redopath !
                         IF (DEBUG) PRINT '(A)', ' output> setting random initial vector for eigenvector'
                         DO J1=1,NINTS
                            VECS(J1)=DPRAND()*2-1.0D0
                         ENDDO
                         CALL VECNORM(VECS,NINTS)
                      ENDIF
                      CALL INTBFGSTS(NSTEPS,XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)),  &
                 &     EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,DEBUG)
                   ELSE
                      IF (DESMINT) THEN
                         TMPINTNEWT = INTNEWT
                         INTNEWT = .FALSE. ! linear transformation only
                         ! convert internal tangents to cartesians
                         CALL TRANSBACKDELTA(TANVEC(1:NOPT,TSPOS(J)-1),VECS,XYZCART(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)), &
                              & NINTC,NOPT,NNZ,KD,FAILED,DEBUG,INTEPSILON)                             
                         INTNEWT = TMPINTNEWT
                         VECSNORM=SUM(VECS(1:NOPT)**2)
                         IF (VECSNORM.EQ.0.0D0) THEN  ! TANVEC ISN't set for GUESSPATH, MECCANO, UNMECCANO
                            IF (DEBUG) PRINT '(A)', ' output> setting random initial vector for eigenvector'
                            DO J1=1,NOPT
                               VECS(J1)=DPRAND()*2-1.0D0
                            ENDDO
                            CALL VECNORM(VECS,NOPT)
                         ENDIF
                      ELSE
                         VECS(1:NOPT)=TANVEC(1:NOPT,TSPOS(J)-1)
                         VECSNORM=SUM(VECS(1:NOPT)**2)
                         IF (VECSNORM.EQ.0.0D0) THEN  ! TANVEC ISN't set for GUESSPATH, MECCANO, UNMECCANO
                            IF (DEBUG) PRINT '(A)', ' output> setting random initial vector for eigenvector'
                            ! This IF block probably doesn't make any difference, but it stops complaints
                            ! about coordinate transformations further down the line.
                            IF (RIGIDINIT) THEN
                               DO J1=1,DEGFREEDOMS
                                   VECS(J1)=DPRAND()*2-1.0D0
                               ENDDO
                               VECS(DEGFREEDOMS+1:) = 0.0D0
                            ELSE
                               DO J1=1,NOPT
                                  VECS(J1)=DPRAND()*2-1.0D0
                               ENDDO
                            ENDIF
                            CALL VECNORM(VECS,NOPT)
                         ENDIF
                      ENDIF
                      IF (GROWSTRINGT.OR.REDOPATH) THEN
                         KNOWG = .FALSE.
                         KNOWE = .FALSE.
                      ENDIF

                      IF (DESMINT) THEN
                         CALL BFGSTS(NSTEPS,XYZCART(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)),  &
                              &       EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
                      ELSE
                         IF (CUDAT) THEN
                            CALL CUDA_BFGSTS_WRAPPER(NSTEPS,XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)),  &
                                 &       EDUMMY,TSCONVERGED,RMS,EVALMIN,VECS,ITDONE)
                         ELSE
                            THISXYZ = XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J))  ! Irritating hack to stop ifort complaining
                            CALL BFGSTS(NSTEPS,THISXYZ,  &
                                 &       EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
                            XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)) = THISXYZ
                         END IF
                      ENDIF
                   ENDIF
                ELSE
                   IF (DESMINT) THEN
                      CALL EFOL(XYZCART(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)),TSCONVERGED, &
                           &   NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,DIAG,2)
                   ELSE
                      CALL EFOL(XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)),TSCONVERGED, &
                           &   NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,DIAG,2)
                   ENDIF
                ENDIF
                CALL MYCPU_TIME(TIME,.FALSE.)

                IF (TSCONVERGED) THEN
                   NTSFOUND=NTSFOUND+1
                   IF (DESMINT) THEN
                      ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),&
                           &TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
                      TSFOUND(NTSFOUND)%VECS=VECS(1:NOPT)
                      TSFOUND(NTSFOUND)%COORD=XYZCART(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J))
                   ELSE
                      ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),&
                           &TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
                      TSFOUND(NTSFOUND)%VECS=VECS(1:NOPT)
                      TSFOUND(NTSFOUND)%COORD=XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J))
                   ENDIF
                   TSFOUND(NTSFOUND)%E=EDUMMY
                   TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
                ENDIF
                IF (TSCONVERGED) THEN
                   WRITE(*,'(1X,A,I6)') 'Converged to TS (number of iterations):     ',ITDONE
                ELSE
                   WRITE(*,'(1X,A,I6)') 'Failed to converge to TS (number of iterations):     ',ITDONE
                ENDIF
             ENDDO
             CALL MYCPU_TIME(ENDTIME,.FALSE.)

             WRITE(INTSTR,'(I10)') NTSFOUND

             IF (MECCANOT) THEN                  
                WRITE(METHSTR,'(A)') 'MECCANO'
             ELSE IF (GROWSTRINGT) THEN
                IF (EVOLVESTRINGT) THEN
                   WRITE(METHSTR,'(A)') 'ES'
                ELSE
                   WRITE(METHSTR,'(A)') 'GS'
                ENDIF
             ELSE
                WRITE(METHSTR,'(A)') 'DNEB'
             ENDIF              

             WRITE(*, '(1x,a,f10.2)',advance='yes') trim(METHSTR)//' run yielded '//trim(adjustl(IntStr))// &
                          &' true transition state(s) time=',EndTime-StartTime
             WRITE(*,*) "Energies:"
             DO J=1,NTSFOUND
                 write(*,*) TSFOUND(J)%E
             ENDDO
          ENDIF
          
          IF (SAVECANDIDATES) THEN
             DO J=1,NTSFOUND
                IF (DESMINT) THEN
                   INQUIRE(IOLENGTH=RECLEN) (XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)))
                ELSE                       
                   INQUIRE(IOLENGTH=RECLEN) (XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)))
                ENDIF
                WRITE(FILENAME,'(i10)') J
                FILENAME='points'//trim(adjustl(filename))//'.out'
                OPEN(UNIT=40,FILE=FILENAME,STATUS='unknown',form='unformatted',access='direct',recl=reclen)

                IF (DESMINT) THEN
                   WRITE(40,REC=1) ( XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)) )
                ELSE
                   WRITE(40,REC=1) ( XYZ(NOPT*(TSPOS(J)-1)+1:NOPT*TSPOS(J)) )
                ENDIF

                CLOSE(40)
             ENDDO
          ENDIF

          RETURN

      END SUBROUTINE TSLOCATOR

SUBROUTINE CONTSLOCATOR
USE KEY,ONLY: BFGSTST,UNRST,NSTEPS,MACHINE, GROWSTRINGT, INTEPSILON, REDOTSIM
USE GSDATA, ONLY: EVOLVESTRINGT
USE KEYOUTPUT
USE MODCHARMM
USE NEBDATA
USE KEYNEB,ONLY:NIMAGE,DEBUG
USE NEBTOCONNECT
USE CHARUTILS
USE MODGUESS
USE MODMEC
USE LINKEDLIST
USE MODEFOL
USE INTCOMMONS, ONLY : DESMINT, NINTC, NNZ, KD, INTNEWT
USE COMMONS, ONLY : REDOPATH, REDOPATHNEB
IMPLICIT NONE
          
INTEGER :: I,J,ITDONE=0,J1,RECLEN,J2,MYTSMAX,NTS
INTEGER,PARAMETER :: MAXPRINTOUT = 50, ITMAX  = 30
DOUBLE PRECISION :: EDUMMY,EVALMIN,EVALMAX,MAXE,VECSNORM
LOGICAL :: TSCONVERGED
DOUBLE PRECISION,DIMENSION(NOPT) :: LGDUMMY, VECS, DIAG, XLOCAL
DOUBLE PRECISION ELOCAL(NIMAGE+2)
INTEGER :: MLOC
DOUBLE PRECISION :: TIME, TIME0
DOUBLE PRECISION :: DPRAND
LOGICAL :: KNOWE, KNOWG, KNOWH ! JMC
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH ! JMC
CHARACTER(LEN=256) :: FILENAME, METHSTR
LOGICAL :: TMPINTNEWT, FAILED
DOUBLE PRECISION, ALLOCATABLE :: TSGUESS(:,:), TSTEMP(:,:), LTANVEC(:,:)

MYTSMAX=10
 PRINT *,' A ALLOCATED(TSGUESS)=',ALLOCATED(TSGUESS)
 PRINT *,' A ALLOCATED(LTANVEC)=',ALLOCATED(LTANVEC)
IF (ALLOCATED(TSGUESS)) DEALLOCATE(TSGUESS)
IF (ALLOCATED(LTANVEC)) DEALLOCATE(LTANVEC)
 PRINT *,' A2 ALLOCATED(TSGUESS)=',ALLOCATED(TSGUESS)
 PRINT *,' A2 ALLOCATED(LTANVEC)=',ALLOCATED(LTANVEC)
ALLOCATE(TSGUESS(MYTSMAX,NOPT),LTANVEC(MYTSMAX,NOPT))
NTS=0
LGDUMMY = 0 ! sn402 added
IF (REDOPATHNEB) THEN
   PRINT '(A,F20.10)',' contslocator> ERROR *** REDOPATH cannot be set with NEBCONSTRAINT'
   STOP
ELSE
   DO I=1,NIMAGE+1
      DO J2=1,NIMAGE+2 ! extra interpolation using the same number of images
         XLOCAL(1:NOPT)=( (NIMAGE+2-J2)*XYZ((I-1)*NOPT+1:I*NOPT)+(J2-1)*XYZ(I*NOPT+1:(I+1)*NOPT) )/(NIMAGE+1)
         CALL POTENTIAL(XLOCAL,ELOCAL(J2),LGDUMMY,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         PRINT '(3(A,I6),A,G20.10)',' contslocator> energy at position ',J2,' between images ',I,' and ',I+1, &
  &                                  ' E=',ELOCAL(J2)
      ENDDO
      IF (ELOCAL(2).LT.ELOCAL(1)) THEN
         NTS=NTS+1
         IF (NTS.GT.MYTSMAX) THEN ! increase storage as required for TS candidates
            ALLOCATE(TSTEMP(MYTSMAX,NOPT))
             PRINT *,' B ALLOCATED(TSGUESS)=',ALLOCATED(TSGUESS)
             PRINT *,' B ALLOCATED(LTANVEC)=',ALLOCATED(LTANVEC)
            TSTEMP(1:MYTSMAX,1:NOPT)=TSGUESS(1:MYTSMAX,1:NOPT)
            DEALLOCATE(TSGUESS)
            ALLOCATE(TSGUESS(2*MYTSMAX,NOPT))
            TSGUESS(1:MYTSMAX,1:NOPT)=TSTEMP(1:MYTSMAX,1:NOPT)
            TSTEMP(1:MYTSMAX,1:NOPT)=LTANVEC(1:MYTSMAX,1:NOPT)
            DEALLOCATE(LTANVEC)
            ALLOCATE(LTANVEC(2*MYTSMAX,NOPT))
            LTANVEC(1:MYTSMAX,1:NOPT)=TSTEMP(1:MYTSMAX,1:NOPT)
            DEALLOCATE(TSTEMP)
            MYTSMAX=2*MYTSMAX
         ENDIF
         PRINT '(3(A,I6),A,G20.10)',' contslocator> adding ts candidate at position ',1,' between images ',I,' and ',I+1, &
  &                               ' E=',ELOCAL(1)
         TSGUESS(NTS,1:NOPT)=XYZ((I-1)*NOPT+1:I*NOPT)
         LTANVEC(NTS,1:NOPT)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
      ENDIF
      DO J2=2,NIMAGE+1 
         IF ( (ELOCAL(J2-1).LT.ELOCAL(J2)) .AND. (ELOCAL(J2).GT.ELOCAL(J2+1)) ) THEN
            NTS=NTS+1
            IF (NTS.GT.MYTSMAX) THEN ! increase storage as required for TS candidates
              PRINT *,' C ALLOCATED(TSGUESS)=',ALLOCATED(TSGUESS)
              PRINT *,' C ALLOCATED(LTANVEC)=',ALLOCATED(LTANVEC)
               ALLOCATE(TSTEMP(MYTSMAX,NOPT))
               TSTEMP(1:MYTSMAX,1:NOPT)=TSGUESS(1:MYTSMAX,1:NOPT)
               DEALLOCATE(TSGUESS)
               ALLOCATE(TSGUESS(2*MYTSMAX,NOPT))
               TSGUESS(1:MYTSMAX,1:NOPT)=TSTEMP(1:MYTSMAX,1:NOPT)
               TSTEMP(1:MYTSMAX,1:NOPT)=LTANVEC(1:MYTSMAX,1:NOPT)
               DEALLOCATE(LTANVEC)
               ALLOCATE(LTANVEC(2*MYTSMAX,NOPT))
               LTANVEC(1:MYTSMAX,1:NOPT)=TSTEMP(1:MYTSMAX,1:NOPT)
               DEALLOCATE(TSTEMP)
               MYTSMAX=2*MYTSMAX
            ENDIF
            PRINT '(3(A,I6),A,G20.10)',' contslocator> adding ts candidate at position ',J2,' between images ',I,' and ',I+1, &
  &                                  ' E=',ELOCAL(J2)
            TSGUESS(NTS,1:NOPT)=( (NIMAGE+2-J2)*XYZ((I-1)*NOPT+1:I*NOPT)+(J2-1)*XYZ(I*NOPT+1:(I+1)*NOPT) )/(NIMAGE+1)
            LTANVEC(NTS,1:NOPT)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
         ENDIF
      ENDDO
      IF (ELOCAL(NIMAGE+1).LT.ELOCAL(NIMAGE+2)) THEN
         NTS=NTS+1
         IF (NTS.GT.MYTSMAX) THEN ! increase storage as required for TS candidates
            ALLOCATE(TSTEMP(MYTSMAX,NOPT))
             PRINT *,' D ALLOCATED(TSGUESS)=',ALLOCATED(TSGUESS)
             PRINT *,' D ALLOCATED(LTANVEC)=',ALLOCATED(LTANVEC)
            TSTEMP(1:MYTSMAX,1:NOPT)=TSGUESS(1:MYTSMAX,1:NOPT)
            DEALLOCATE(TSGUESS)
            ALLOCATE(TSGUESS(2*MYTSMAX,NOPT))
            TSGUESS(1:MYTSMAX,1:NOPT)=TSTEMP(1:MYTSMAX,1:NOPT)
            TSTEMP(1:MYTSMAX,1:NOPT)=LTANVEC(1:MYTSMAX,1:NOPT)
            DEALLOCATE(LTANVEC)
            ALLOCATE(LTANVEC(2*MYTSMAX,NOPT))
            LTANVEC(1:MYTSMAX,1:NOPT)=TSTEMP(1:MYTSMAX,1:NOPT)
            DEALLOCATE(TSTEMP)
            MYTSMAX=2*MYTSMAX
         ENDIF
         PRINT '(3(A,I6),A,G20.10)',' contslocator> adding ts candidate at position ',NIMAGE+2,' between images ',I,' and ',I+1, &
  &                               ' E=',ELOCAL(NIMAGE+2)
         TSGUESS(NTS,1:NOPT)=XYZ(I*NOPT+1:(I+1)*NOPT) 
         LTANVEC(NTS,1:NOPT)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
      ENDIF
   ENDDO
ENDIF

IF (NTS.EQ.0) THEN
   PRINT '(A)',' contslocator> No ts candidates to optimise'
   STOP
ENDIF

! WRITE(*,'(1x,a)',advance='No') 'Converged to TS (number of iterations):     '

NTSFOUND=0
CALL MYCPU_TIME(STARTTIME,.FALSE.)
DO J=1,NTS
   CALL MYCPU_TIME(TIME0,.FALSE.)
   KNOWE=.FALSE.
   KNOWG=.FALSE.
   IF (BFGSTST) THEN
      IF (UNRST) THEN 
         PRINT '(A)',' contslocator> ERROR *** not coded for UNRES'
         STOP
      ELSE
         VECS(1:NOPT)=LTANVEC(J,1:NOPT)
                            PRINT *,'in output before bfgsts B'
         CALL BFGSTS(NSTEPS,TSGUESS(J,1:NOPT),EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
                            PRINT *,'in output after bfgsts B'
      ENDIF
   ELSE
      CALL EFOL(TSGUESS(J,1:NOPT),TSCONVERGED,NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,DIAG,2)
   ENDIF
   CALL MYCPU_TIME(TIME,.FALSE.)

   IF (TSCONVERGED) THEN
      NTSFOUND=NTSFOUND+1
      ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
      TSFOUND(NTSFOUND)%VECS=VECS(1:NOPT)
      TSFOUND(NTSFOUND)%COORD=TSGUESS(J,1:NOPT)
      TSFOUND(NTSFOUND)%E=EDUMMY
      TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
      WRITE(*,'(1X,A,I6)') 'Converged to TS (number of iterations):     ',ITDONE
!     WRITE(*,'(i5)',advance='No') itdone
   ELSE
!     WRITE(*,'(a5)',advance='No') '   :('
!     WRITE(*,'(A)',advance='No') 'Failed to converge to TS'
      WRITE(*,'(1X,A,I6)') 'Failed to converge to TS (number of iterations):     ',ITDONE
   ENDIF
ENDDO
CALL MYCPU_TIME(ENDTIME,.FALSE.)

WRITE(*,'(a)') '.'

WRITE(INTSTR,'(i10)') NTSfound

WRITE(*, '(A,F7.2)',advance='yes') ' Constrained potential run yielded '//trim(adjustl(IntStr))// &
                  &' true transition state(s) time=',EndTime-StartTime

 PRINT *,' end ALLOCATED(TSGUESS)=',ALLOCATED(TSGUESS)
 PRINT *,' end ALLOCATED(LTANVEC)=',ALLOCATED(LTANVEC)
IF (ALLOCATED(TSGUESS)) DEALLOCATE(TSGUESS)
RETURN

END SUBROUTINE CONTSLOCATOR

      SUBROUTINE CHECKTS(DUMMY,EVALMIN,TSCONVERGED)
          USE NEBDATA
          USE MODCHARMM
          USE LINKEDLIST
          IMPLICIT NONE

          LOGICAL :: FAILCHECK,TSCONVERGED
          TYPE(CHAIN),POINTER :: DUMMY
          DOUBLE PRECISION :: EVALMIN

          ! DAE If EVALMIN large in magnitude, this TS is likely to be bogus, and cause problems
          ! when then the connected minima have to be found
!          IF (CHRMMT.AND.(EVALMIN.LT.-100.D0)) THEN
!               TSConverged=.FALSE.
!               WRITE(*,'(A,F20.10,A)') 'checkts> Eigenvalue ',EVALMIN,' too negative, TS search failed'
!               ! DAE for CHARMM check this transition state to see if its geometry has become unfeasible
!!               CALL CHECKPOINT(xyz(nopt*(dummy%i-1):nopt*dummy%i),FAILCHECK)
!!               CALL CHECKPOINT(xyz(nopt*(dummy%i-1)+1:nopt*dummy%i),FAILCHECK) ! bs360
!!               IF (FAILCHECK) THEN
!!                    WRITE(*,'(A)') 'checkts> Transition state has unphysical geometry, TS search failed'
!!                    TSConverged=.FALSE.
!!               ENDIF
!          ENDIF
      END SUBROUTINE CHECKTS
END MODULE NEBOUTPUT
