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
MODULE msevb_common
  IMPLICIT NONE
  SAVE
  
! Parameters

  DOUBLE PRECISION, PARAMETER :: TRIGUNITY = 0.99999999999999D0
  DOUBLE PRECISION, PARAMETER :: MINUSTRIGUNITY = -0.99999999999999D0
! Arbitrary large force for discontinuities in bond angle terms etc
  DOUBLE PRECISION, PARAMETER :: LARGE_FORCE = 1.0D5 
  DOUBLE PRECISION, PARAMETER :: NUMERICALZEROLIMIT = 1.0D-15

! Store interactions of various kinds

! Coulomb interactions

! each_coulomb(i,j) - atom i is part of exchange species, atom j is a normal water atom
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: EACH_COULOMB
! both atoms are normal water atoms
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: WATER_INTER_COULOMB
! inter_coulomb(i,j) - atom i is H3O+ atom, atom j is a normal H2O atom
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: INTER_COULOMB
! Sum of each_coulomb for each atom
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ATOM_COULOMB

! H2O-H2O LJ interaction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: LJR
! H2O-H3O+ LJ interaction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: LJ_INTER

! H2O-H3O+ O-O repulsion term
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: REPULSE_INTER

! Exchange species terms

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ZUNDEL_F
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ZUNDEL_G

!$omp threadprivate (each_coulomb,water_inter_coulomb,ljr,inter_coulomb,lj_inter,repulse_inter,atom_coulomb)
!$omp threadprivate (zundel_f,zundel_g)

!----important parameters needed throughout the program

  INTEGER :: num_hyd,num_eig,ground_state,reduced_num_eig,maxNumVBstates

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XX, YY, ZZ
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PSIX, PSIY, PSIZ  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: INTERATOMICR

!$omp threadprivate (ground_state,reduced_num_eig,maxNumVBstates,xx,yy,zz,psix,psiy,psiz,interAtomicR)

! Parameters for the assignment of the VB states

  INTEGER :: shellsToCount
  DOUBLE PRECISION :: MAXHBONDLENGTH, MINHBONDANGLE

! Geometry check parameter

  DOUBLE PRECISION :: OOCLASH_SQ

! Exchange interactions between different VB states

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VIJEXCH

!$omp threadprivate (vijexch)

! Array to hold the identities of the exchange hydronium species

  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: zundel_species
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: statesInteract  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HAM, AIJ
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: atmpl  
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GRSTWFU

!$omp threadprivate (zundel_species,statesInteract,ham,aij,atmpl,grstwfu)

! Flag to print out coefficients in ground state at the end of minimisation

  logical :: printCoefficients

! MSEVB potential parameters

  DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979D0
  DOUBLE PRECISION, PARAMETER :: RADS = PI/1.8D02
!----Newish Terms:
  DOUBLE PRECISION, PARAMETER :: H2OINTEREPSILON = 0.1522D0
  DOUBLE PRECISION, PARAMETER :: H2OINTERSIGMA = 3.1506D0
!----avogadro's constant
  DOUBLE PRECISION, PARAMETER :: NA = 6.02214D23
!----diss. energy of o-h bond in h3o (kJ/mol)
  DOUBLE PRECISION, PARAMETER :: DOH = 266.3D0
!----exponent prefactor eq o-h bond length in h3o (m)
  DOUBLE PRECISION, PARAMETER :: AOHEQ = 1.285D0
!----eq o-h bond length in h3o (m)
  DOUBLE PRECISION, PARAMETER :: ROHEQ = 0.980D0
!----force constant for h3o (kJ/mol.deg^2)
  DOUBLE PRECISION, PARAMETER :: KALPHA = 73.27D0
!----eq h3o angle (deg)
  DOUBLE PRECISION, PARAMETER :: ALPHAEQ = 116.0D0*RADS
!----LJ epsilon for h3o-h2o intermolecular term (m)
  DOUBLE PRECISION, PARAMETER :: EPSILON_MSEVB = 0.155D0
!  DOUBLE PRECISION, PARAMETER :: EPSILON_MIX = 4.00D0*DSQRT(EPSILON_MSEVB*H2OINTEREPSILON)
  DOUBLE PRECISION, PARAMETER :: EPSILON_MIX = 4.00D0*0.15359362D0
!----LJ sigma for h3o-h2o intermolecular term (m)
  DOUBLE PRECISION, PARAMETER :: SIGMA = 3.164D0 
  DOUBLE PRECISION, PARAMETER :: SIGMA_MIX = (H2OINTERSIGMA+SIGMA)/2.0D0
!----coulombic oxygen charge for h3o-h2o intermolecular term (C)
  DOUBLE PRECISION, PARAMETER :: QINTERO =  -0.5D0
!----coulombic oxygen charge for h3o-h2o intermolecular term (C)
  DOUBLE PRECISION, PARAMETER :: QINTERH = 0.5D0
!----o-o repulsive coefficient
  DOUBLE PRECISION, PARAMETER :: BIG_B = 2.591D0
!----o-o repulsive coefficient exponent coefficient
  DOUBLE PRECISION, PARAMETER :: SMALL_B = 3.50D0
!----repulsive eq separation
  DOUBLE PRECISION, PARAMETER :: DOOEQ = 2.50D0
!----alpha term from off-diagonal (f3) fitting
  DOUBLE PRECISION, PARAMETER :: ALPHA_BLEH = 15.0D0
!----off-diagonal eq roo term
  DOUBLE PRECISION, PARAMETER :: ROOEQ = 1.90D0
!----beta term in repulsive thing (m)
  DOUBLE PRECISION, PARAMETER :: BETA = 4.500D0
!----who cares
  DOUBLE PRECISION, PARAMETER :: BIG_ROOEQ = 3.14D0
!----another ff-diagonal eq doo term
  DOUBLE PRECISION, PARAMETER :: DOO = 2.875D0
!----P term in off-diagonal elements
  DOUBLE PRECISION, PARAMETER :: P = 0.270D0
!----k term in off-diagonal element (m^-2)
  DOUBLE PRECISION, PARAMETER :: KEXP = 11.50D0
!----gamma for off-diagonal term (m^-2)
  DOUBLE PRECISION, PARAMETER :: GAMMA_MSEVB = 1.850D0
!----zeta in off-diagonal terms (rad^-2)
  DOUBLE PRECISION, PARAMETER :: ETA = 1.500D0
!----omega in off-diagonal terms (deg)
  DOUBLE PRECISION, PARAMETER :: OMEGAEQ = 174.3D0*RADS
!----the coupling terms (kJ/mol)
!----Dell
  DOUBLE PRECISION, PARAMETER :: VIJ = -32.925
!----exchange charge for oxygen
  DOUBLE PRECISION, PARAMETER :: QEXCHO = -1.0005D-01
!----exchange charge for hydrogen
  DOUBLE PRECISION, PARAMETER :: QEXCHH = 4.002D-02
!----Dang Parameterized TIP3P
  DOUBLE PRECISION, PARAMETER :: ALJ = 5.80D05
  DOUBLE PRECISION, PARAMETER :: CLJ = 5.25D02
  DOUBLE PRECISION, PARAMETER :: H2OROHEQ = 0.960D0
  DOUBLE PRECISION, PARAMETER :: H2OKTHETA = 68.087D0
!----Dell
  DOUBLE PRECISION, PARAMETER :: H2OKB = 105.9162D01
  DOUBLE PRECISION, PARAMETER :: THETAEQ = 104.5D0*RADS
!----Dell
  DOUBLE PRECISION, PARAMETER :: QTIP3P_O = -0.834D0
  DOUBLE PRECISION, PARAMETER :: QTIP3P_H = 0.417D0
  DOUBLE PRECISION, PARAMETER :: DAMPCHGE = 7.6D-01

!----coulomb unit term
  DOUBLE PRECISION, PARAMETER :: FOURPIEO = 332.0638186797117D0

! Debugging variables

!	real*8 LJ_values(mxeig), coulomb_values(mxeig)
!	real*8 rep_values(mxeig)
!	real*8 h3o_intra_energy(mxeig), h2o_intra_energy(mxeig)
!	real*8 h3o_h2o_energy(mxeig), h2o_h2o_energy(mxeig)
!	real*8 LJ_forces(mxeig,3*MXATOM), coulomb_forces(mxeig,3*MXATOM)
!	real*8 rep_forces(mxeig,3*MXATOM)

!  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: H3O_INTRA_FORCES
!  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: H2O_INTRA_FORCES
!  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: H3O_H2O_FORCES
!  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: H2O_H2O_FORCES

!	real*8 off_diag_forces(mxeig,mxeig,3*MXATOM)
!	real*8 off_diag_energies(mxeig,mxeig)

END MODULE msevb_common













