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
MODULE KEYNEB
     USE KEY,ONLY: NEBK,NEBKINITIAL,NEBKFINAL,NEBFACTOR
     USe COMMONS,ONLY: DEBUG
     USE KEYGRAD
     USE KEYOUTPUT
     USE KEYMINIMIZER
     USE KEYSQVV
     USE KEYLBFGS
     USE KEYTAU
     IMPLICIT NONE
     SAVE
     INTEGER :: NIMAGE              = 3
     LOGICAL :: MOREPRINTING        = .FALSE.
     LOGICAL :: DUMPNEBXYZ          = .FALSE.
     INTEGER :: DUMPNEBXYZFREQ      = 100
     LOGICAL :: DUMPNEBPTS          = .FALSE.
     INTEGER :: DUMPNEBPTSFREQ      = 100
     LOGICAL :: DUMPNEBEOS          = .FALSE.
     INTEGER :: DUMPNEBEOSFREQ      = 100
     LOGICAL :: READGUESS           = .FALSE.
     CHARACTER(LEN=10) :: XYZFILE   = 'neb.xyz   '
     CHARACTER(LEN=12) :: RBXYZFILE = 'rbneb.xyz   '
     CHARACTER(LEN=10) :: PTSFILE   = 'neb.out   '
     CHARACTER(LEN=10) :: EOFSFILE  = 'neb.EofS  '
     CHARACTER(LEN=80) :: GUESSFILE = 'guess.xyz '
     LOGICAL :: OLDCONNECT          = .FALSE.
     !logical :: debug               = .False.
     CONTAINS

     SUBROUTINE KEYNEBPRINT(VARIABLEIN)
          USE CHARUTILS
          USE KEY, ONLY : DESMAXAVGE, DESMAXEJUMP, DNEBSWITCH
          IMPLICIT NONE
          LOGICAL,INTENT(IN),OPTIONAL :: VARIABLEIN
          
          LOGICAL :: VARIABLE
          
          WRITE(*,'(1X,A,2G19.10,A,G19.10)') 'KeyNEB> Initial and final NEB force constants ', &
  &                                           NEBKINITIAL,NEBKFINAL,' factor=',NEBFACTOR
          VARIABLE=.FALSE.
          IF (PRESENT(VARIABLEIN)) VARIABLE=VARIABLEIN
          IF (DNEBSWITCH.GT.0.0D0) WRITE(*,'(A,G20.10)') &
  &               ' KeyNEB> Gradient will switch to single-nudged formulation for RMS<',DNEBSWITCH
          IF (VARIABLE) THEN
               WRITE(*,'(1x,a)') 'KeyNEB> Number of images will vary depending on the separation of the endpoints'
          ELSE
               INTSTR=WI(NIMAGE)
               WRITE(*,'(1x,a)') 'KeyNEB> Using '//trim(IntStr)//' images'
          ENDIF
          
          IF (DUMPNEBXYZ) THEN
               INTSTR=WI(DUMPNEBXYZFREQ)
               WRITE(*,'(1x,a)') 'KeyNEB> NEB coordinates will be saved to xyz file "'//trim(XyzFile)//'" every '&
               &//TRIM(INTSTR)//' iterations'
          ENDIF
          IF (DUMPNEBPTS) THEN
               INTSTR=WI(DUMPNEBPTSFREQ)
               WRITE(*,'(1x,a)') 'KeyNEB> NEB coordinates will be saved to binary file "'//trim(PtsFile)//'" every '&
               &//TRIM(INTSTR)//' iterations'
          ENDIF
          IF (DUMPNEBEOS) THEN
               INTSTR=WI(DUMPNEBEOSFREQ)
               WRITE(*,'(1x,a)') 'KeyNEB> Energy profile will be saved to file "'//trim(EofSFile)//'" every '&
               &//TRIM(INTSTR)//' iterations'
          ENDIF

          IF (READGUESS) WRITE(*,'(1x,a)') 'KeyNEB> Guess will be read from file "'//trim(GuessFile)//'"'
          IF (OLDCONNECT) WRITE(*,'(1x,a)') "KeyNEB> Transition state candidate will be returned in first endpoint's coord array"
          IF (DEBUG) THEN
               WRITE(*,'(1x,a)') "KeyNEB> Verbose printing is on"
               WRITE(*,'(1x,a)') 'KeyNEB> Evolution of AvDev, E, Rms and S will be saved to '&
               &//'files "AvDevofI", "EofI", "RmsofI" and "SofI"'
          ENDIF
          IF (DESMAXAVGE.LT.0.99*HUGE(1.0D0)) WRITE(*,'(1x,a,G20.10)') 'KeyGS> Max average energy: ', DESMAXAVGE
          IF (DESMAXEJUMP.LT.0.99*HUGE(1.0D0)) WRITE(*,'(1x,a,G20.10)') 'KeyGS> Max energy jump per image: ', DESMAXEJUMP
          
     END SUBROUTINE KEYNEBPRINT

     SUBROUTINE ALLKEYNEBPRINT(VARIABLE)
          IMPLICIT NONE
          LOGICAL,INTENT(IN),OPTIONAL :: VARIABLE
          IF (PRESENT(VARIABLE)) THEN
               CALL KEYNEBPRINT(VARIABLE)
          ELSE
               CALL KEYNEBPRINT
          ENDIF
          CALL KEYGRADPRINT
          CALL KEYOUTPUTPRINT
          IF (PRESENT(VARIABLE)) THEN
               CALL KEYMINIMIZERPRINT(VARIABLE)
          ELSE
               CALL KEYMINIMIZERPRINT
          ENDIF
          CALL KEYTAUPRINT
     END SUBROUTINE ALLKEYNEBPRINT
END MODULE KEYNEB
