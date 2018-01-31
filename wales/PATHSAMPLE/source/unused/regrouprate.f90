!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Regroup/reorder A, B and I !!!!!!!!!!!!!!!!!!!!!!!!!
!
!  In this case we want a many-to-one mapping, which changes the number of "minima" and
!  the effective "transition states".
!
!
SUBROUTINE REGROUPRATE(MINGROUP)
USE COMMONS
IMPLICIT NONE
INTEGER J1, J2, NGROUPS
INTEGER MINGROUP(NMIN), NEWNMINA, NEWNMINB, NCOUNT, NDUMMY, NEMPTY, GROUPMAP(NMIN)
LOGICAL ISA(NMIN), ISB(NMIN), CHANGED, GROUPA(NMIN), GROUPB(NMIN)
DOUBLE PRECISION LOGREGROUPRATETHRESH, NEWEMIN(NMIN), NEWETS(NTS), NEWPFMIN(NMIN), NEWKPLUS(NTS), NEWKMINUS(NTS)
INTEGER NEWNMIN, NEWNTS, NEWPLUS(NTS), NEWMINUS(NTS), NMINGROUP(NMIN), NTSGROUP(NTS)

IF (ENSEMBLE.EQ.'E') THEN
   PRINT '(A)','regrouprate> Regrouped entropies not yet coded for microcanonical ensemble'
   STOP
ENDIF
!
!  All minima belonging to the same group as an A minimum become A, and
!  similarly for B.
!
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO

NGROUPS=0
!
!  Assign minima to new groups.
!
LOGREGROUPRATETHRESH=LOG(REGROUPRATETHRESH)
DO J1=1,NMIN
   MINGROUP(J1)=0 ! MINGROUP(J1) is the index of the group containing minimum J1
ENDDO
DO
   CHANGED=.FALSE.
   DO J1=1,NTS
      IF ((KPLUS(J1).GT.LOGREGROUPRATETHRESH).AND.(KMINUS(J1).GT.LOGREGROUPRATETHRESH)) THEN
         IF ((MINGROUP(PLUS(J1)).EQ.0).AND.(MINGROUP(MINUS(J1)).EQ.0)) THEN
            CHANGED=.TRUE.
            NGROUPS=NGROUPS+1
            MINGROUP(PLUS(J1))=NGROUPS
            MINGROUP(MINUS(J1))=NGROUPS
         ELSEIF (MINGROUP(PLUS(J1)).NE.MINGROUP(MINUS(J1))) THEN
            CHANGED=.TRUE.
            IF (MINGROUP(PLUS(J1)).EQ.0) THEN
               MINGROUP(PLUS(J1))=MINGROUP(MINUS(J1))
            ELSEIF (MINGROUP(MINUS(J1)).EQ.0) THEN
               MINGROUP(MINUS(J1))=MINGROUP(PLUS(J1))
            ELSE
               MINGROUP(PLUS(J1))=MIN(MINGROUP(PLUS(J1)),MINGROUP(MINUS(J1)))
               MINGROUP(MINUS(J1))=MINGROUP(PLUS(J1))
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   IF (.NOT.CHANGED) EXIT
ENDDO

NDUMMY=NGROUPS
DO J1=1,NMIN
   IF (MINGROUP(J1).EQ.0) THEN
      NGROUPS=NGROUPS+1
      MINGROUP(J1)=NGROUPS
   ENDIF
ENDDO
NEMPTY=0
DO J1=1,NGROUPS
   NCOUNT=0
   DO J2=1,NMIN
      IF (MINGROUP(J2).EQ.J1) THEN
         NCOUNT=NCOUNT+1
      ENDIF
   ENDDO
   NMINGROUP(J1)=NCOUNT
   IF (NCOUNT.EQ.0) THEN
      IF (DEBUG) PRINT '(A,I6,A)','regrouprate> Group ',J1,' is empty'
      NEMPTY=NEMPTY+1
   ENDIF
ENDDO
PRINT '(3(A,I6))','regrouprate> Number of groups=',NGROUPS,' number of size one=',NGROUPS-NDUMMY,' number of empty groups=',NEMPTY
!
!  A set.
!
NEWNMINA=0
GROUPA(1:NGROUPS)=.FALSE.
DO J1=1,NMIN
   IF (ISA(J1)) GROUPA(MINGROUP(J1))=.TRUE.
ENDDO
DO J1=1,NMIN
   IF (GROUPA(MINGROUP(J1))) THEN
      ISA(J1)=.TRUE.
      NEWNMINA=NEWNMINA+1
      IF (DEBUG) PRINT '(A,I6,A)','regroup> minimum ',J1,' is in the A set'
   ENDIF
ENDDO
!
!  B set.
!
NEWNMINB=0
GROUPB(1:NGROUPS)=.FALSE.
DO J1=1,NMIN
   IF (ISB(J1)) GROUPB(MINGROUP(J1))=.TRUE.
ENDDO
DO J1=1,NMIN
   IF (GROUPB(MINGROUP(J1))) THEN
      ISB(J1)=.TRUE.
      NEWNMINB=NEWNMINB+1
      IF (DEBUG) PRINT '(A,I6,A)','regroup> minimum ',J1,' is in the B set'
   ENDIF
ENDDO
PRINT '(2(A,I6))','regrouprate> After regrouping number of A minima=',NEWNMINA,' number of B minima=',NEWNMINB
DO J1=1,NMIN
   IF (ISA(J1).AND.ISB(J1)) THEN
      PRINT '(A,I6,A)','regroup> ERROR - minimum ',J1,' belongs to A and B sets'
      STOP
   ENDIF
ENDDO
!
!  Need to reset NMIN, NTS, EMIN and ETS to the corresponding grouped
!  quantities and free energies.
!  Note that some of the NGROUPS of minima are actually empty (NEMPTY).
!  LOCATIONA and LOCATIONB are used in GT, so we need to reset them.
!  Probably also need to redo the TOPPOINTER and POINTER stuff for the
!  regrouped database.
!  We can't renumber everything until we have calculated the free energy
!  of the grouped transition states.
!
!  Only the odd factor of Plank's constant that shifts transition state
!  free energies from minima is included.
!
NEWEMIN(1:NMIN)=HUGE(1.0D0)
NEWPFMIN(1:NMIN)=0.0D0
NEWNMIN=NGROUPS
DO J1=1,NMIN
   NEWPFMIN(MINGROUP(J1))=NEWPFMIN(MINGROUP(J1))+EXP(PFMIN(J1))
ENDDO
DO J1=1,NGROUPS
   IF (NMINGROUP(J1).EQ.0) CYCLE
   NEWEMIN(J1)=-TEMPERATURE*LOG(NEWPFMIN(J1))
   IF (DEBUG) PRINT '(A,I6,A,G20.10,A,I6,A,G20.10)','regrouprate> For group ',J1,' Z(T)=',NEWPFMIN(J1),' size ',NMINGROUP(J1), &
  &                           ' free energy=',NEWEMIN(J1)
ENDDO

NEWKPLUS(1:NTS)=0.0D0
NEWKMINUS(1:NTS)=0.0D0
NEWNTS=0
NTSGROUP(1:NTS)=0
tsloop: DO J1=1,NTS
   IF (MINGROUP(PLUS(J1)).EQ.MINGROUP(MINUS(J1))) CYCLE ! Ignore intragroup rates
   DO J2=1,NEWNTS
      IF ((MINGROUP(PLUS(J1)).EQ.NEWPLUS(J2)).AND.(MINGROUP(MINUS(J1)).EQ.NEWMINUS(J2))) THEN
         NEWKPLUS(J2)=NEWKPLUS(J2)+EXP(PFMIN(PLUS(J1))+KPLUS(J1))/NEWPFMIN(MINGROUP(PLUS(J1)))
         NEWKMINUS(J2)=NEWKMINUS(J2)+EXP(PFMIN(MINUS(J1))+KMINUS(J1))/NEWPFMIN(MINGROUP(MINUS(J1)))
         NTSGROUP(J2)=NTSGROUP(J2)+1 
         CYCLE tsloop
      ELSEIF ((MINGROUP(PLUS(J1)).EQ.NEWMINUS(J2)).AND.(MINGROUP(MINUS(J1)).EQ.NEWPLUS(J2))) THEN
         NEWKPLUS(J2)=NEWKPLUS(J2)+EXP(PFMIN(MINUS(J1))+KMINUS(J1))/NEWPFMIN(MINGROUP(MINUS(J1)))
         NEWKMINUS(J2)=NEWKMINUS(J2)+EXP(PFMIN(PLUS(J1))+KPLUS(J1))/NEWPFMIN(MINGROUP(PLUS(J1)))
         NTSGROUP(J2)=NTSGROUP(J2)+1 
         CYCLE tsloop
      ENDIF
   ENDDO
   NEWNTS=NEWNTS+1
   NEWPLUS(NEWNTS)=MINGROUP(PLUS(J1))
   NEWMINUS(NEWNTS)=MINGROUP(MINUS(J1))
   NEWKPLUS(NEWNTS)=EXP(PFMIN(PLUS(J1))+KPLUS(J1))/NEWPFMIN(MINGROUP(PLUS(J1)))
   NEWKMINUS(NEWNTS)=EXP(PFMIN(MINUS(J1))+KMINUS(J1))/NEWPFMIN(MINGROUP(MINUS(J1)))
   NTSGROUP(NEWNTS)=1
ENDDO tsloop

PRINT '(A,I6)','regrouprate> Number of intergroup transition states=',NEWNTS
DO J1=1,NEWNTS
   IF (DEBUG) PRINT '(3(A,I6),2(A,G20.10),A,I6)','regrouprate> Grouped ts ',J1,' between minima groups ',NEWPLUS(J1), &
  &    ' and ',NEWMINUS(J1), &
  &    ' k+=',NEWKPLUS(J1),' k-=',NEWKMINUS(J1),' members=',NTSGROUP(J1)
!  PRINT '(A,2G20.10)','regrouprate> detailed balance - these numbers should be equal: ',NEWKPLUS(J1)*NEWPFMIN(NEWPLUS (J1)), &
! &                                                                                      NEWKMINUS(J1)*NEWPFMIN(NEWMINUS(J1))

   IF ((NEWKPLUS(J1).EQ.0.0D0).OR.(NEWKMINUS(J1).EQ.0.0D0)) THEN
      IF (DEBUG) PRINT '(A,I6,2G20.10)','regrouprate> WARNING - J1,NEWKPLUS,NEWKMINUS=',J1,NEWKPLUS(J1),NEWKMINUS(J1)
      NEWETS(J1)=HUGE(1.0D0)
   ELSE
      NEWETS(J1)=NEWEMIN(NEWPLUS(J1))-TEMPERATURE*LOG(NEWKPLUS(J1)*PLANCK/TEMPERATURE)
      IF (DEBUG) PRINT '(2(A,G20.10))','regrouprate> Grouped ts free energy=', &
  &                     NEWETS(J1), &
  &              ' or ',NEWEMIN(NEWMINUS(J1))-TEMPERATURE*LOG(NEWKMINUS(J1)*PLANCK/TEMPERATURE)
      IF (NEWETS(J1).NE.0.0D0) THEN ! Check for consistency
         IF (ABS((NEWETS(J1)-NEWEMIN(NEWMINUS(J1))+TEMPERATURE*LOG(NEWKMINUS(J1)*PLANCK/TEMPERATURE))/NEWETS(J1)).GT.0.01D0) &
  &      PRINT '(A,I6,A,2G20.10)','regrouprate> WARNING - free energies for ts group ',J1,' are ',  &
  &                  NEWETS(J1),NEWEMIN(NEWMINUS(J1))-TEMPERATURE*LOG(NEWKMINUS(J1)*PLANCK/TEMPERATURE)
      ENDIF
   ENDIF
ENDDO
!
! We now have all the free energies for minima and transition state groups. 
! Some of the original groups for the minima will generally be empty, so 
! now we renumber.
!
! POINTERS have not been renumbered (yet).
!
NCOUNT=0
NDUMMY=0
OPEN(UNIT=1,FILE='min.data.regrouped',STATUS='UNKNOWN')
NMINA=0
NMINB=0
!
!  If we are going to analyse the min.data.regrouped.resorted and ts.data.regrouped.resorted
!  files for rates subsequently, then we have to arrange for the ln products of frequencies
!  to give us a factor of (kT/h). This can be done by setting the ln product equal to zero
!  for the transition state and 2 * ln(2*Pi*k*T/h) for all the minima. We already have h in the
!  units of kT, so this is easy. The 2*Pi factor occurs because the frequencies are assumed to be
!  angular normal mmode frequencies, and the factor of two occurs because they are assumed
!  to be squared.
!
LNPROD=2.0D0*LOG(2.0D0*3.141592654D0*TEMPERATURE/PLANCK))

DO J1=1,NGROUPS
   IF (NMINGROUP(J1).GT.0) THEN
      NCOUNT=NCOUNT+1
      GROUPMAP(J1)=NCOUNT
      NMINGROUP(NCOUNT)=NMINGROUP(J1)
      PFMIN(NCOUNT)=LOG(NEWPFMIN(J1))
      EMIN(NCOUNT)=NEWEMIN(J1)
      GROUPA(NCOUNT)=GROUPA(J1)
      GROUPB(NCOUNT)=GROUPB(J1)
      IF (GROUPA(NCOUNT)) THEN
         NMINA=NMINA+1
         LOCATIONA(NMINA)=NCOUNT
      ENDIF
      IF (GROUPB(NCOUNT)) THEN
         NMINB=NMINB+1
         LOCATIONB(NMINB)=NCOUNT
      ENDIF
      NDUMMY=NDUMMY+NMINGROUP(NCOUNT)
      WRITE(1,'(2G20.10,I6,4F20.10)') EMIN(NCOUNT),LNPROD,1,1.0,1.0,1.0
   ENDIF
ENDDO
CLOSE(1)
NGROUPS=NCOUNT
PRINT '(4(A,I6))','regrouprate> Number of groups after removing empty sets=',NGROUPS,' total minima=',NDUMMY, &
  &               ' # A: ',NMINA,' # B: ',NMINB
IF (NDUMMY.NE.NMIN) THEN
   PRINT '(A,I6)','regrouprate> ERROR - number of minima in groups should be ',NMIN
   STOP
ENDIF
DO J1=1,NMIN
   MINGROUP(J1)=GROUPMAP(MINGROUP(J1))
ENDDO
NTS=NEWNTS
OPEN(UNIT=1,FILE='ts.data.regrouped',STATUS='UNKNOWN')
DO J1=1,NTS
   PLUS(J1)=GROUPMAP(NEWPLUS(J1))
   MINUS(J1)=GROUPMAP(NEWMINUS(J1))
   IF (NEWKPLUS(J1).GT.0.0D0) THEN
      KPLUS(J1)=LOG(NEWKPLUS(J1))
   ELSE
      KPLUS(J1)=-HUGE(1.0D0)
   ENDIF
   IF (NEWKMINUS(J1).GT.0.0D0) THEN
      KMINUS(J1)=LOG(NEWKMINUS(J1))
   ELSE
      KMINUS(J1)=-HUGE(1.0D0)
   ENDIF
   ETS(J1)=NEWETS(J1)
   WRITE(1,'(2G20.10,3I10,3F20.10)') ETS(J1),0.0,1,PLUS(J1),MINUS(J1),1.0,1.0,1.0
ENDDO
CLOSE(1)
OPEN(UNIT=1,FILE='min.A.regrouped',STATUS='UNKNOWN')
WRITE(1,'(I6)') NMINA
DO J1=1,NMIN
   IF (GROUPA(J1)) WRITE(1,'(I6)') J1
ENDDO
CLOSE(1)
OPEN(UNIT=1,FILE='min.B.regrouped',STATUS='UNKNOWN')
WRITE(1,'(I6)') NMINB
DO J1=1,NMIN
   IF (GROUPB(J1)) WRITE(1,'(I6)') J1
ENDDO
CLOSE(1)
PFTOTALA=0.0D0
PFTOTALB=0.0D0
DO J1=1,NMINB
   PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J1)))
ENDDO
PFTOTALB=LOG(PFTOTALB)
 DO J1=1,NMINA
   PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J1)))
ENDDO
PFTOTALA=LOG(PFTOTALA)
NMIN=NGROUPS

PRINT '(A)','regrouprate>  NOTE: from here on down min and ts refer to the new groups!'

RETURN

END
