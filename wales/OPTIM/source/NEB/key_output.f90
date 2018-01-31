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
MODULE KEYOUTPUT
     IMPLICIT NONE
     SAVE
     LOGICAL :: OPTIMIZETS          = .FALSE.
     LOGICAL :: PRINTOPTIMIZETS     = .FALSE.
     LOGICAL :: SAVECANDIDATES      = .FALSE.
     CHARACTER(LEN=5) :: CANDIDATES = 'maxim'
     CONTAINS

     SUBROUTINE KEYOUTPUTPRINT
          IMPLICIT NONE

          IF (OPTIMIZETS) THEN
               WRITE(*,'(1x,a)') 'KeyOutput> Transition state candidates will be optimized'
               IF (PRINTOPTIMIZETS) THEN
                    WRITE(*,'(1x,a)') 'KeyOutput> Verbose printing during transition states optimization'
               ELSE
                    WRITE(*,'(1x,a)') 'KeyOutput> Concise printing during transition states optimization'
               ENDIF
          ELSE
               WRITE(*,'(1x,a)') 'KeyOutput> Transition state candidates will NOT be optimized'
          ENDIF

          IF (SAVECANDIDATES) &
          &WRITE(*,'(1x,a)') 'KeyOutput> Coordinates for transition states candidates will be saved to files "points*.out"'

          SELECT CASE (CANDIDATES)
          CASE ('maxim')
               WRITE(*,'(1x,a)') 'KeyOutput> Transition state candidates are maxima along NEB'
          CASE ('all')
               WRITE(*,'(1x,a)') 'KeyOutput> Transition state candidates are all images'
          CASE ('high')
               WRITE(*,'(1x,a)') 'KeyOutput> Transition state candidate is the highest energy image'
          END SELECT
     END SUBROUTINE KEYOUTPUTPRINT
END MODULE KEYOUTPUT
