SUBROUTINE SETUPAMB(filename)
implicit none
        CHARACTER(*) :: filename
        WRITE(*,*) 'sf344> shouldn`t be in this routine!  Compiled correctly?'
END SUBROUTINE SETUPAMB

SUBROUTINE AMBERINTERFACE(DUMMY1,DUMMY2,DUMMY3,DUMMY4)
implicit none
        INTEGER         :: DUMMY1, DUMMY2, DUMMY4
CHARACTER(LEN=*) DUMMY3

END SUBROUTINE AMBERINTERFACE

subroutine amber_readcoords(ostring)
implicit none
character(len=20) ostring

end subroutine amber_readcoords

SUBROUTINE AMBERFINALIO(nsave,nunit,node,ostring,filth,coords2)
implicit none

integer nsave,nunit,filth,node
character(len=80) ostring
double precision coords2(*)

END SUBROUTINE AMBERFINALIO

SUBROUTINE AMBERSECDER(OLDX,STEST)
implicit none
DOUBLE PRECISION  :: OLDX(*)
LOGICAL    :: STEST
END SUBROUTINE AMBERSECDER

SUBROUTINE AMBERENERGIES(y,grad,ereal,gradt,stest)
implicit none
double precision :: y(*),grad(*),ereal
logical          :: gradt,stest

END SUBROUTINE AMBERENERGIES

subroutine check_cistrans_rna(coords,natoms,atomlabels,goodstructure)
implicit none

integer,intent(in) :: natoms
character(len=5),intent(in) :: atomlabels(natoms)
double precision,intent(in) :: coords(3*natoms)
logical,intent(out) :: goodstructure

end subroutine check_cistrans_rna

subroutine check_cistrans_dna(coords,natoms,atomlabels,goodstructure)
implicit none

integer,intent(in) :: natoms
character(len=5),intent(in) :: atomlabels(natoms)
double precision,intent(in) :: coords(3*natoms)
logical,intent(out) :: goodstructure

end subroutine check_cistrans_dna

subroutine set_check_chiral(coords,natoms,goodstructure,chiarray)
! alternative subroutine for checking chirality that can be used 
! for preserving arbitrary starting chiralities in peptide chain
! (should detect inversion of the chiral centre in the same way 
! as we are diagnosing cis/trans isomerisation)
implicit none

logical :: goodstructure
integer,intent(in) :: natoms
double precision,intent(in) :: coords(3*natoms)
integer :: chiarray(natoms)
end subroutine set_check_chiral

SUBROUTINE MME2WRAPPER(COORDS,ENERGY,VNEW,TEMPHESS,ATMASS,GRAD1)
implicit none

double precision coords(*),energy,vnew(*),temphess(*),atmass(*),grad1(*)

END SUBROUTINE MME2WRAPPER

SUBROUTINE MMEINITWRAPPER(prmtop,igb,saltcon,rgbmax,cut)
implicit none

character       :: prmtop
integer         :: igb
double precision :: saltcon,rgbmax,cut

END SUBROUTINE MMEINITWRAPPER

SUBROUTINE CHECK_CISTRANS_PROTEIN(coords,natoms,goodstructure,minomega,cisarray)
implicit none

! stuff for detecting cis-trans isomerisation of omega angles in peptide bonds
logical,intent(out) :: goodstructure
integer,intent(in) :: natoms
double precision,intent(in) :: coords(3*natoms),minomega
integer :: cisarray(natoms)

END SUBROUTINE CHECK_CISTRANS_PROTEIN

subroutine check_chirality(coords,natoms,goodstructure)
implicit none

integer,intent(in) :: natoms
double precision,intent(in) :: coords(3*natoms)
logical :: goodstructure
end subroutine check_chirality

subroutine amber_deallocate_stacks()

end subroutine amber_deallocate_stacks

subroutine takestepamber(np)
implicit none
integer np
end subroutine takestepamber


SUBROUTINE  AMBERDIHEDR(Q,NATOMS,I,J,K,L,ORDER)
IMPLICIT NONE
DOUBLE PRECISION Q(*),ORDER
INTEGER NATOMS,I,J,K,L
PRINT*,'ERROR in AMBERDIHEDR - this dummy routine should never be called - wrong executable?'
RETURN
END


subroutine  GETSEEDATM(dummy1,dummy2,dummy3,at1,at2,at3,IRES)

character :: dummy1, dummy2, dummy3
integer :: at1, at2, at3, IRES

end subroutine GETSEEDATM

subroutine GETICCOORDS()

end subroutine GETICCOORDS

subroutine CHECK_SIDECHAIN(i3,j3,k3,l3,IICD,IS_SIDECHAIN)
integer :: i3,j3,k3,l3,IICD
logical :: IS_SIDECHAIN(100)

end subroutine CHECK_SIDECHAIN

subroutine CHGETICVAL(Q,BOND1, BOND2, THET1, THET2, PHI, DIHED_ONLY)
double precision :: Q(3*5000), BOND1(5000), BOND2(5000), THET1(5000)
double precision :: THET2(5000), PHI(5000)
logical :: DIHED_ONLY
end subroutine CHGETICVAL

double precision function PHI_INTERP(ST_PHI,FI_PHI,SFRAC)
double precision:: ST_PHI,FI_PHI,PHI,SFRAC
end function PHI_INTERP

double precision function PROBABIL(x, xmax, CHPMAX,CHPMIN, max_x_func, slope_par)
double precision :: CHPMAX, CHPMIN
integer :: x, xmax, max_x_func, slope_par
end function PROBABIL

subroutine  TAKESTEPAMDIHED(COORDS, CSTART, CFINISH,SFRAC)
double precision :: COORDS(3*5000), CSTART(3*5000), CFINISH(3*5000), SFRAC
end subroutine TAKESTEPAMDIHED

subroutine SETDIHEAM()
end subroutine SETDIHEAM

subroutine ICTYPECHECKAM(LPHIPSI,LOMEGAC,LSIDECHAIN, LCHIRAL,IICD)
logical:: LPHIPSI,LOMEGAC,LSIDECHAIN,LCHIRAL
integer :: IICD
end subroutine ICTYPECHECKAM

subroutine GET_TWISTABLE(ST_PHI, FI_PHI)
double precision:: ST_PHI(3000), FI_PHI(3000)
end subroutine GET_TWISTABLE

subroutine CHARMMBILDC(NSTART,NSTOP,Q, BOND1,BOND2,THET1,THET2,PIC)
integer::NSTART,NSTOP
double precision ::Q(3*5000), BOND1(5000), BOND2(5000), THET1(5000)
double precision :: THET2(5000), PIC(5000)
end subroutine CHARMMBILDC

subroutine  CHREBUILD(Q, BOND1, BOND2, THET1, THET2, PIC)
double precision :: Q(3*5000), BOND1(5000), BOND2(5000), THET1(5000)
double precision :: THET2(5000), PIC(5000)
end subroutine CHREBUILD

subroutine TAKESTEPAMM(Q, DEBUG, BHSTEPSIZE)
double precision :: Q(3*5000)
logical :: debug
double precision :: BHSTEPSIZE
end subroutine TAKESTEPAMM

subroutine PERTDIHAM(Q,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED, debug, BHSTEPSIZE)
double precision :: Q(3*5000)
integer::  CHNMIN,CHNMAX,ISEED, CHPMIN,CHPMAX
logical :: debug
double precision :: BHSTEPSIZE
end subroutine  PERTDIHAM

!SUBROUTINE AMBALIGNDIH(PHI,DIHC,RELPHI,PREVRELPHI,N)
!double precision::PHI, RELPHI, PREVRELPHI, DRP, PHIORIG
!INTEGER:: DIHC, N
!end subroutine AMBALIGNDIH

!SUBROUTINE AMBGETOUTOFPLANE(I, J, K, L, PHI, dPdI, dPdJ, dPdK, dPdL, NOCOOR,NODERV)
!double precision:: I(0:2), J(0:2), K(0:2), L(0:2), PHI
!double precision:: dPdI(0:2), dPdJ(0:2), dPdK(0:2), dPdL(0:2)
!logical:: NOCOOR, NODERV
!end subroutine AMBGETOUTOFPLANE

!SUBROUTINE AMBGETANGLE(I, J, K, THETA, dTdI, dTdJ, dTdK, NOCOOR,NODERV)
!LOGICAL:: NOCOOR,NODERV
!DOUBLE PRECISION:: I(3), J(3), K(3), THETA
!DOUBLE PRECISION:: dTdI(3), dTdJ(3), dTdK(3)
!end subroutine AMBGETANGLE

!SUBROUTINE AMBGETTORSION(I, J, K, L, PHI, dPdI, dPdJ, dPdK, dPdL, NOCOOR, NODERV)
!LOGICAL:: NOCOOR, NODERV
!DOUBLE PRECISION:: I(3), J(3), K(3), L(3), PHI
!DOUBLE PRECISION:: dPdI(3), dPdJ(3), dPdK(3), dPdL(3)
!end subroutine AMBGETTORSION

!SUBROUTINE AMB_GETKDNAT(KD)
!INTEGER::KD
!end subroutine AMB_GETKDNAT

!subroutine AMBGETNNZNAT(NNZ)
!INTEGER::NNZ
!end subroutine AMBGETNNZNAT

!subroutine AMBTRANSFORM(XCART,GCART,XINT,GINT,NINTC,NCART,NNZ,NOCOOR,NODERV,KD,INTEPSILON)
!INTEGER:: NNZ,NINTC,NCART,KD
!DOUBLE PRECISION:: XCART(3*MXATMS),GCART(3*MXATMS),XINT(NINTC),GINT(NINTC)
!DOUBLE PRECISION:: INTEPSILON
!LOGICAL:: NOCOOR, NODERV
!end subroutine AMBTRANSFORM

!subroutine AMB_TRANSBACKDELTA(X,CARTX,COORDS,NINTC,NCART,NNZ,KD,FAILED,PTEST2,INTEPSILON)
!double precision INTEPSILON,CARTX(NCART),X(NINTC),COORDS(3*NATOM)
!integer NNZ,NINTC,NCART
!logical PTEST2,FAILEd
!end subroutine AMB_TRANSBACKDELTA

subroutine check_valleu_chirality(COORDSB, COORDSA,DEBUG)
double precision :: COORDSB(*), COORDSA(*)
logical :: debug
end subroutine check_valleu_chirality

SUBROUTINE valleu_swap(COORDSB,COORDSA,i1,i4, DEBUG)
double precision, intent(in):: COORDSB(:)
double precision, intent(inout) :: COORDSA(:)
integer, intent(in) ::i1,i4
logical :: DEBUG
END SUBROUTINE valleu_swap

subroutine A9RESTOATOM(frozenres,frozen,nfreeze)
integer i1, resnumber, currentresidue,nfreeze
logical frozenres(:), frozen(:)

end subroutine A9RESTOATOM

subroutine A9ATOMTORES(atom,residue)
integer, intent(in) :: atom
integer, intent(out) :: residue
integer :: i1, resnumber,currentresidue

end subroutine A9ATOMTORES

subroutine A9DUMPMODES(DIAG,ZT,N,M)
INTEGER :: N,M
DOUBLE PRECISION :: DIAG(M)
LOGICAL :: ZT(M)

end subroutine A9DUMPMODES


!      SUBROUTINE GET_TWISTABLE(ST_PHI,FI_PHI)
!      IMPLICIT NONE
!      DOUBLE PRECISION ST_PHI(*), FI_PHI(*)
!
!      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
!      STOP
!
!      RETURN
!      END

 !     SUBROUTINE TAKESTEPAMM(Q,DEBUG, BHSTEPSIZE)
 !     IMPLICIT NONE
 !     DOUBLE PRECISION Q(*), BHSTEPSIZE
 !     LOGICAL             DEBUG
 !
 !     PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
 !     STOP
 !
 !     RETURN
 !     END

