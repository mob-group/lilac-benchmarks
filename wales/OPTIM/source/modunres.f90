!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
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
MODULE modunres
   IMPLICIT NONE
   SAVE

! jmc unres related keywords

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: UREFCOORD,UREFPPSANGLE ! STORE REFERENCE CARTESIANS AND INTERNALS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: INTSTEP ! INTERNAL COORDINATE STEP VECTOR OFF TS
! these 2 lines are for use with CONSEC (surprisingly enough) keyword to interpolate internal coords over a section of the 
! molecule only (selected a/c to residue numbers and stored in [start/end]res() as appropriate.)
      LOGICAL :: CONSECT
      INTEGER :: STARTRES(10),ENDRES(10),NUMSEC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: MYQMINSAVE
      INTEGER COUNTER

! jmc unres.h
! added COMMON /interact/ to know itype()

      INTEGER :: UNmaxres, UNmaxres2, maxdim, maxvar, ntyp, ntyp1
      PARAMETER (UNmaxres=100,UNmaxres2=2*UNmaxres,maxdim=(UNmaxres-1)*(UNmaxres-2)/2,maxvar=4*UNmaxres)
      PARAMETER (ntyp=20,ntyp1=ntyp+1)

      INTEGER :: nres,nres0,ialph,ivar,ntheta,nphi,nside,nvaru,nfl,icg
      DOUBLE PRECISION :: THETA,PHI,ALPH,OMEG
      DOUBLE PRECISION :: C,DC,XLOC,XROT,DC_NORM, &
     &                    dcdv,dxdv,dxds,gradx,gradc,gvdwc,gelc,gradx_scp, &
     &                    gvdwc_scp,ghpbx,ghpbc,gloc,gradcorr,gradxorr,gvdwx

      COMMON /chain/ c(3,UNmaxres2+2),dc(3,UNmaxres2),xloc(3,UNmaxres), &
     & xrot(3,UNmaxres),dc_norm(3,UNmaxres2),nres,nres0

      DOUBLE PRECISION :: DSC,VBL,VBLINV,VBLINV2,VBL_CIS,VBL0
      COMMON /sclocal/ dsc(ntyp1)
      COMMON /peptbond/ vbl,vblinv,vblinv2,vbl_cis,vbl0

      COMMON /var/ theta(UNmaxres),phi(UNmaxres),alph(UNmaxres),omeg(UNmaxres), &
     &          ialph(UNmaxres,2),ivar(4*UNmaxres2),ntheta,nphi,nside,nvaru

      COMMON /derivat/ dcdv(6,maxdim),dxdv(6,maxdim),dxds(6,UNmaxres), &
     &                 gradx(3,UNmaxres,2),gradc(3,UNmaxres,2),gvdwx(3,UNmaxres), &
     &                 gvdwc(3,UNmaxres),gelc(3,UNmaxres),gradx_scp(3,UNmaxres), &
     &                 gvdwc_scp(3,UNmaxres),ghpbx(3,UNmaxres),ghpbc(3,UNmaxres), &
     &                 gloc(maxvar,2),gradcorr(3,UNmaxres),gradxorr(3,UNmaxres),nfl,icg

      CHARACTER(LEN=3) :: restyp
      CHARACTER(LEN=1) :: onelet
      COMMON /names/ restyp(ntyp+1),onelet(ntyp+1)

      INTEGER :: nnt,nct,nint_gr,istart,iend,itype,maxint_gr,expon,expon2, &
     &           itel,itypro,ielstart,ielend,nscp_gr,iscpstart,iscpend,iatsc_s, &
     &           iatsc_e,iatel_s,iatel_e,iatscp_s,iatscp_e,ispp,iscp
      PARAMETER (maxint_gr=2)
      DOUBLE PRECISION :: AA,BB,AUGM,AAD,BAD,APPU,BPP,AEL6,AEL3

      COMMON /interact/aa(ntyp,ntyp),bb(ntyp,ntyp),augm(ntyp,ntyp), &
     & aad(ntyp,2),bad(ntyp,2),appu(2,2),bpp(2,2),ael6(2,2),ael3(2,2), &
     & expon,expon2,nnt,nct,nint_gr(UNmaxres),istart(UNmaxres,maxint_gr), &
     & iend(UNmaxres,maxint_gr),itype(UNmaxres),itel(UNmaxres),itypro, &
     & ielstart(UNmaxres),ielend(UNmaxres),nscp_gr(UNmaxres), &
     & iscpstart(UNmaxres,maxint_gr),iscpend(UNmaxres,maxint_gr), &
     & iatsc_s,iatsc_e,iatel_s,iatel_e,iatscp_s,iatscp_e,ispp,iscp

!     INTEGER LTM(UNMAXRES,200)
!     COMMON /RESIDINFO/ LTM

! mass info for Calpha, 'D' capping group (gly) and for the side chains.
! Side chain masses are the sums of the atomic masses.
      DOUBLE PRECISION :: MASSES(NTYP+2)
      COMMON /unresmass/ MASSES

! Values set in initialize_p.F via block data.
! order is (following unres' restyp (in initialize_p.F):
! C, CYS, MET , PHE , ILE , LEU , VAL , TRP , TYR , ALA , GLY , THR ,
! SER , GLN , ASN , GLU , ASP , HIS , ARG , LYS , PRO , D 
! Data from http://www.pharmacy.purdue.edu/~mcmp304/AATutorial/
! all at physiological pH, therefore arg,asp,glu,lys are charged.
!     MASSES(1)=12.01
!     MASSES(2)=47.1
!     MASSES(3)=75.1
!     MASSES(4)=91.1
!     MASSES(5)=57.1
!     MASSES(6)=57.1
!     MASSES(7)=43.1
!     MASSES(8)=130.2
!     MASSES(9)=107.1
!     MASSES(10)=15.0
!     MASSES(11)=1.0
!     MASSES(12)=45.1
!     MASSES(13)=31.0
!     MASSES(14)=72.1
!     MASSES(15)=58.1
!     MASSES(16)=72.1
!     MASSES(17)=58.1
!     MASSES(18)=81.1
!     MASSES(19)=101.2
!     MASSES(20)=73.1
!     MASSES(21)=41.1
!     MASSES(22)=0.0

END MODULE
