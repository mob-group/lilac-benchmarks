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
      SUBROUTINE UNRSETDIHE
      USE MODUNRES
      IMPLICIT NONE

      INTEGER STEPNUMBER,RUNNUMBER,REMAINDER,J1

      DOUBLE PRECISION XU(NVARU)
      INTEGER FIRSTSIDE,LASTSIDE

      STEPNUMBER=1
      RUNNUMBER=1
      REMAINDER=1

      PRINT *,'unrsetdihe here0'

C     CALL ASSIGNRAMACH(STEPNUMBER,REMAINDER,RUNNUMBER)

      CALL geom_to_var(nvaru,XU)

      PRINT *,'unrsetdihe: number of residues ',nres
      PRINT *,'unrsetdihe: number of virtual-bond dihedral angles ',nphi
      PRINT *,'unrsetdihe: number of side chain dihedral angles ',nside
      IF (nres.ne.nside) THEN
         PRINT *,nres-nside,' residue(s) are glycine - no side chain dihedral or polar angles.'
C jmc need 'if itype() .ne. 10 then...' itype is in COMMON.INTERACT.
      END IF
      DO J1=1,nphi
         PRINT *,'PHI ',XU(J1)
      END DO

      FIRSTSIDE=2*nres-4+nside
      LASTSIDE=2*nres-5+2*nside

      DO J1=FIRSTSIDE,LASTSIDE
         PRINT *,'BETA ',XU(J1)
      END DO
 
      RETURN
      END
