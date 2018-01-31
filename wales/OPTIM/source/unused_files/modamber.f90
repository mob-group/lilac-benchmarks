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

MODULE MODAMBER
      IMPLICIT NONE
      SAVE

      DOUBLE PRECISION  TOP,BOTTOM,AX,BX,CX,AY,BY,CY,AZ,BZ,CZ,MU,LAMBDA,NUMER,DENOM,DX,DY,DZ
      DOUBLE PRECISION  TOTENERGY,BENERGY,TENERGY,PENERGY,ANSWER,IMPENERGY,QENERGY,VDWENERGY
      
      INTEGER ATOMS,ANG,IMP,COUNT,T
      INTEGER  a,b,ios,d,e,f,g,h,i,j,k,l,colin,fish,ans,rings
      INTEGER  match,bean,bondnumber
      INTEGER NDIHEDRALS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DATOM1, DATOM2
      INTEGER loopatom(20)
      INTEGER,ALLOCATABLE, DIMENSION (:,:) :: bondarray
      DOUBLE PRECISION TINY
      PARAMETER  (TINY=1.0D-12)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::  X,Y,Z,PQ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) ::  R,VDWA,VDWB
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::  MASS

      DOUBLE PRECISION GENTORSPARAMS(100,6),SPECTORSPARAMS(100,14),SPECIMPPARAMS(15,7)
      DOUBLE PRECISION MIDIMPPARAMS(4,6),GENIMPPARAMS(15,5)
  
      INTEGER, ALLOCATABLE, DIMENSION (:) :: aa1,aa2,aa3
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) ::  THETA

      INTEGER, ALLOCATABLE, DIMENSION (:) :: da1,da2,da3,da4
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: DPHI,DID,DVN,DVN2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: DVN3,DDELTA,DDELTA2,DDELTA3
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: DN,DN2,DN3

      INTEGER, ALLOCATABLE, DIMENSION (:) :: ia1,ia2,ia3,ia4
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: IPHI,IVN,IDELTA,IN1

      DOUBLE PRECISION  QX(3),QY(3),QZ(3)
      INTEGER, ALLOCATABLE, DIMENSION (:) :: atnum,bondedto,type
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: bonds
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: one_four,one_three
      CHARACTER(LEN=2)  typechb,typechc,typechd,typeche
      CHARACTER(LEN=20)  filename
      CHARACTER, ALLOCATABLE, DIMENSION(:) ::  label
      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:) :: typech
!
! Variables for energy calculation
!
      CHARACTER(LEN=10)  check
      DOUBLE PRECISION DIELEC
      PARAMETER  (dielec=1)
      LOGICAL FAKEWATER
      DOUBLE PRECISION RSTAR,EPSILON
      DOUBLE PRECISION VDWR(42),VDWE(42)
      DOUBLE PRECISION KR(42,42),RO(42,42)

      DOUBLE PRECISION KT(42,42,42),TO(42,42,42),DEGTO(42,42,42)
      DOUBLE PRECISION PK,PHASE,PK2,PK3,PHASE2,PHASE3,IPK,IPHASE

      INTEGER PN,IDIVF,PN2,PN3,IPN
!
! Variables for derivative calculation
!
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:) :: DBONDEBYDX,DBONDEBYDY,DBONDEBYDZ,DANGEBYDX, &
     &                 dangEbydy,dangEbydz, &
     &                 dtorsEbydx,dtorsEbydy,dtorsEbydz,dvdwEbydx, &
     &                 dvdwEbydy,dvdwEbydz, &
     &                 dEbydx,dEbydy,dEbydz,dimpEbydx, &
     &                 dimpEbydy,dimpEbydz, &
     &                 dqEbydx,dqEbydy,dqEbydz
!
! Variables for second derivative calculation
!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BONDHELL,ANGLEHELL,TORSHELL,  &
     &                 imphell,qhell,vdwhell,hell
      INTEGER hellcount
      INTEGER,ALLOCATABLE,DIMENSION(:) :: chiral
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: chiralarray

END MODULE MODAMBER
