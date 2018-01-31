C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C*************************************************************************
C
C  Here we calculate the Gupta potential and gradient
C                                        
C*************************************************************************
C
      SUBROUTINE GUPTA(NATOMS,X,V,PG,GRADT,ATOMTYPE)
      IMPLICIT NONE 
      INTEGER J1, J2, atomtype,NATOMS
      DOUBLE PRECISION X(3*NATOMS), PG, DIST, PTEMP, V(3*NATOMS),
     1       GA, GRHO(NATOMS), VTEMP, VTEMP1, VTEMP2, VTEMP3, DUMMY, 
     2       RTEMP, p, q, r0, Gxi, twoq, RytoeV
      LOGICAL GRADT

      PG=0.0D0

C Cleri and Rosato Parameters
C PRB 48, 22 (1993)
C
C Fcc metals
C 

      if (atomtype.eq.1) then

C Nickel
        p=16.999d0
        q=1.189d0
        GA=0.0376d0
        Gxi=1.070d0

      else if (atomtype.eq.2) then

C Copper
        p=10.960d0
        q=2.278d0
        GA=0.0855d0
        Gxi=1.224d0

      else if (atomtype.eq.3) then

C Rhodium
        p=18.450d0
        q=1.867d0
        GA=0.0629d0
        Gxi=1.660d0

      else if (atomtype.eq.4) then

C Palladium
        p=10.867d0
        q=3.742d0
        GA=0.1746d0
        Gxi=1.718d0

      else if (atomtype.eq.5) then

C Silver
        p=10.928d0
        q=3.139d0
        GA=0.1028d0
        Gxi=1.178d0

      else if (atomtype.eq.6) then

C Iridium
        p=16.980d0
        q=2.691d0
        GA=0.1156d0
        Gxi=2.289d0

      else if (atomtype.eq.7) then

C Platinum
        p=10.612d0
        q=4.004d0
        GA=0.2975d0
        Gxi=2.695d0

      else if (atomtype.eq.8) then

C Gold
        p=10.229d0
        q=4.036d0
        GA=0.2061d0
        Gxi=1.790d0

      else if (atomtype.eq.9) then

C Aluminium
        p=8.612d0
        q=2.516d0
        GA=0.1221d0
        Gxi=1.316d0

      else if (atomtype.eq.10) then

C Lead
        p=9.576d0
        q=3.648d0
        GA=0.0980d0
        Gxi=0.914d0
C
C Hcp metals
C 

      else if (atomtype.eq.11) then

C Titaniium I
        p=8.620d0
        q=2.390d0
        GA=0.1519d0
        Gxi=1.8112d0

      else if (atomtype.eq.12) then

C Titaniium II
        p=11.418d0
        q=1.643d0
        GA=0.0741d0
        Gxi=1.4163d0

      else if (atomtype.eq.13) then

C Zirconium I
        p=8.250d0
        q=2.249d0
        GA=0.1934d0
        Gxi=2.2792d0

      else if (atomtype.eq.14) then

C Zirconium II
        p=13.940d0
        q=1.071d0
        GA=0.0523d0
        Gxi=1.4489d0

      else if (atomtype.eq.15) then

C Cobalt
        p=11.604d0
        q=2.286d0
        GA=0.0950d0
        Gxi=1.4880d0

      else if (atomtype.eq.16) then

C Cadmium I
        p=10.612d0
        q=5.206d0
        GA=0.1420d0
        Gxi=0.8117d0

      else if (atomtype.eq.17) then

C Cadmium II
        p=13.639d0
        q=3.908d0
        GA=0.0416d0
        Gxi=0.4720d0

      else if (atomtype.eq.18) then

C Zinc
        p=9.689d0
        q=4.602d0
        GA=0.1477d0
        Gxi=0.8900d0

      else if (atomtype.eq.19) then

C Magnesium
        p=12.820d0
        q=2.257d0
        GA=0.0290d0
        Gxi=0.4992d0
C
C bcc
C

      else if (atomtype.eq.20) then

C Vanadium
        p=5.206d0
        q=1.22d0
        GA=0.6124d0
        Gxi=2.441d0
      

      else if (atomtype.eq.21) then

C Parameters for Sodium
  
        p=10.13d0
        q=1.3d0
C r0 in bohr
C        r0=6.99d0
C energies in eV
C        GA=0.015955d0
C        Gxi=0.29113d0
C energies in Rydberg
        GA=1.1727d-3
        Gxi=21.398d-3
C convert to eV
        RytoeV=4.359748d-18/(2*1.6021773d-19)
        GA=GA*RytoeV
        Gxi=Gxi*RytoeV

      else if (atomtype.eq.22) then

C Strontium - Wang + Blaisten-Barojas JCP 115 3640 (2001)

        p=17.649d0
        q=(6.0d0-1.1132d0)*(2.0d0/natoms)**0.65d0+1.1132d0
        GA=(0.018089d0-0.015761d0)*(2.0d0/natoms)**0.52588d0+0.015761d0
        Gxi=(0.13977d0-0.3832d0)*(2.0d0/natoms)**0.8d0+0.3832d0
C       distance unit
C        r0=(5.6d0-4.301d0)*(2.0d0/natoms)**3.6d0+4.301d0

      else if (atomtype.eq.23) then

C Au - as used by Garzon et al

        p=10.15d0
        q=4.13d0
        GA=0.118438d0/2.0d0
        Gxi=0.5d0

      endif

C r0 in reduced units
      r0=1.0d0

      twoq=2*q

      DO 22 J1=1,NATOMS
         RTEMP=0.0d0
         PTEMP=0.0d0
         DO 23 J2=1,NATOMS
            if (j1.ne.j2) then
               DIST=DSQRT(( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1           ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2           ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2)
               dist=dist/r0
               PTEMP=PTEMP+dexp(p*(1-dist))
               RTEMP=RTEMP+dexp(twoq*(1-dist))
            endif
23       CONTINUE
         GRHO(J1)=DSQRT(RTEMP)
         PG=PG+GA*PTEMP-Gxi*GRHO(J1)
22    CONTINUE

C
C Now calculate the gradient analytically.
C

      if (GRADT) THEN
      DO J1=1,NATOMS
         VTEMP1=0.0D0
         VTEMP2=0.0D0
         VTEMP3=0.0D0
         DUMMY=1.0D0/GRHO(J1)
         DO J2=1,NATOMS
            if (j1.ne.j2) then
               DIST=DSQRT(( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1           ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2           ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2)
               dist=dist/r0
               VTEMP=(q*Gxi*(DUMMY+1.0D0/GRHO(J2))*dexp(twoq*(1-dist))
     1                -2.0d0*GA*p*dexp(p*(1-dist)))/(r0*dist)
               VTEMP1=VTEMP1+VTEMP*(X(3*(J1-1)+1)-X(3*(J2-1)+1))
               VTEMP2=VTEMP2+VTEMP*(X(3*(J1-1)+2)-X(3*(J2-1)+2))
               VTEMP3=VTEMP3+VTEMP*(X(3*(J1-1)+3)-X(3*(J2-1)+3))
            endif
         ENDDO
         V(3*(J1-1)+1)=VTEMP1
         V(3*(J1-1)+2)=VTEMP2
         V(3*(J1-1)+3)=VTEMP3
      ENDDO
      endif

      RETURN
      END

