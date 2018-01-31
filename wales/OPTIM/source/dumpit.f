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
C
C  Produce a new input deck at the end of the run.
C
      SUBROUTINE DUMPIT(Q,FNAMEF)
      USE COMMONS
      USE KEY
      USE AMBER12_INTERFACE_MOD
      USE OPEP_INTERFACE_MOD, ONLY: OPEP_WRITE_PDB
      IMPLICIT NONE
      INTEGER J1, NTYPEA, J2, J3
      DOUBLE PRECISION S, EPSAB, EPSBB, SIGAB, SIGBB, CAPSCOORDS(18), Q(3*NATOMS), DSHIFT
      CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80, FNAME*16
      CHARACTER(LEN=*) FNAMEF
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      DOUBLE PRECISION CAPSRHO, EPS2, RAD, HEIGHT, XS(3*NTSITES)
      COMMON /CAPS/ CAPSRHO, EPS2, RAD, HEIGHT

      IF (FILTH.EQ.0) THEN
         FNAME='odata.new'
      ELSE
         WRITE(FNAME,'(A)') 'odata.new.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

      OPEN (3,FILE=FNAME,STATUS='UNKNOWN')

      WRITE(3,'(A,A87)') 'COMMENT  ',ESTRING
      IF (GPSTRING(1:5).EQ.' The ') WRITE(3,'(A,T15,A80)') 'COMMENT',GPSTRING
C     WRITE(3,'(A,T15,A80)') 'COMMENT',NSTRING
C     WRITE(3,'(A,T15,A80)') 'COMMENT',FSTRING
      IF ((.NOT.BFGSTST).AND.(.NOT.BFGSMINT).AND.(.NOT.(BSMIN.OR.RKMIN))) WRITE(3,'(A,T15,I5)') 'SEARCH',INR
      IF (.NOT.BFGSMINT) WRITE(3,'(A,T15,F20.10)') 'MAXSTEP',STPMAX(NOPT)
      IF (MAXMAX.NE.0.5D0) WRITE(3,'(A,T15,F20.10)') 'MAXMAX',MAXMAX
      IF (MINMAX.NE.0.01D0) WRITE(3,'(A,T15,F20.10)') 'MINMAX',MINMAX
      IF (GRADSQ) WRITE(3,'(A)') 'GRADSQ'
      IF (CONVU.NE.1.0D-5) THEN
         IF (CONVR.EQ.1.0D-5) THEN
            WRITE(3,'(A,T15,F20.10)') 'CONVERGE',CONVU
         ELSE
            WRITE(3,'(A,T15,2F20.10)') 'CONVERGE', CONVU, CONVR
         ENDIF
      ENDIF
      IF ((GFRACTION.NE.0.0D0).OR.(MFRACTION1.NE.0.0D0).OR.
     1    (MFRACTION2.NE.0.0D0)) WRITE(3,'(A,T15,3F20.10)') 'NONLOCAL', GFRACTION, MFRACTION1, 
     2     MFRACTION2
      IF (SYMCUT.NE.0.001D0) WRITE(3,'(A,T15,F20.10)') 'SYMCUT',SYMCUT
      IF ((TRAD.NE.0.4D0).AND.(.NOT.(BFGSMINT.OR.BSMIN.OR.RKMIN))) WRITE(3,'(A,T15,F20.10)') 'TRAD',TRAD
C     IF (RESIZE.NE.1.0D0) WRITE(3,'(A,T15,F20.10)') 'RESIZE',RESIZE
      IF (RTEST) THEN
         IF (JZ.NE.0.0D0) THEN
            WRITE(3,'(A,T15,F20.10)') 'ROT JZ',JZ
         ELSE
            WRITE(3,'(A,T15,F20.10)') 'ROT OMEGA',OMEGA
         ENDIF
      ENDIF
      IF (PV) WRITE(3,'(A,3F20.10)') 'PV ',PRESS,PVCONV,PVTOL
      IF (PUSHOFF.NE.0.01D0) WRITE(3,'(A,T15,F20.10)') 'PUSHOFF',PUSHOFF
      IF (DFTBT) WRITE(3,'(A)') 'DFTB'
      IF (CADPAC) WRITE(3,'(A,A,A)') 'CADPAC ',SYS,EDITIT
      IF (GAMESSUS) WRITE(3,'(A,A,A)') 'GAMESS-US ',SYS,EDITIT
      IF (GAMESSUK) WRITE(3,'(A,A,A)') 'GAMESS-UK ',SYS,EDITIT
      IF (GAUSSIAN) WRITE(3,'(A)') 'GAUSSIAN'
      IF (GAUSSIAN03) WRITE(3,'(A)') 'GAUSSIAN03'
      IF (GAUSSIAN09) WRITE(3,'(A)') 'GAUSSIAN09'
      IF (GAUSSIAN09) WRITE(3,'(A)') 'GAUSSIAN16'
      IF (ONETEP) THEN
         IF (DFTP) THEN
            WRITE(3,'(A,A)') 'ONETEP ',ONETEPJOB
         ELSE
            WRITE(3,'(A,A)') 'ONETEPC ',ONETEPJOB
         ENDIF
      ENDIF
      IF (CASTEP) THEN
         IF (DFTP) THEN
            WRITE(3,'(A,A)') 'CASTEP ',CASTEPJOB
         ELSE
            WRITE(3,'(A,A)') 'CASTEPC ',CASTEPJOB
         ENDIF
      ENDIF
      IF (CP2K) THEN 
         IF (DFTP) THEN
            WRITE(3,'(A,A)') 'CP2K ',CP2KJOB
         ELSE
            WRITE(3,'(A,A)') 'CP2KC ',CP2KJOB
         ENDIF
      ENDIF
      IF (CPMD) THEN
         if (TRIM(CPMD_COMMAND).ne.'/home/trj25/bin/cpmd.x') then
            write(3, '(A,A)') 'CPMD_COMMAND ', TRIM(CPMD_COMMAND)
         endif
         IF (CPMDC) THEN
            WRITE(3,'(A,A,A)') 'CPMDC ',SYS
         ELSE
            WRITE(3,'(A,A,A)') 'CPMD ',SYS
         ENDIF
      ENDIF
      IF (ISTCRT.NE.10) WRITE(3,'(A,T15,I5)') 'SCALE',ISTCRT
      IF (IPRNT.NE.0) WRITE(3,'(A,T15,I5)') 'PRINT',IPRNT
      IF (DTEST) WRITE(3,'(A,T15,F15.6,2I5)') 'GDIIS',PCUT,NDIIA,NINTV
      IF (MASST) WRITE(3,'(A)') 'MASS  ON'
      IF (NVALUES.NE.20) WRITE(3,'(A,T15,I5)') 'VALUES',NVALUES
      IF (EFSTEPST) WRITE(3,190) EFSTEPS
190   FORMAT('EFSTEPS',T15,I5)
      IF (NSTEPS.NE.1) WRITE(3,200) NSTEPS
200   FORMAT('STEPS',T15,I5)
      IF (DUMPV) WRITE(3,210)
210   FORMAT('DUMPVECTOR')
      IF (PGRAD) WRITE(3,220) NGRADIENTS
220   FORMAT('GRADIENTS',T15,I5)
      IF (VECTORST) WRITE(3,230) NVECTORS
230   FORMAT('VECTORS',T15,I5)
      IF (NSUMMARY.NE.20) WRITE(3,240) NSUMMARY
240   FORMAT('SUMMARY',T15,I5)
      IF (SHIFTV.NE.1.0D6) WRITE(3,253) SHIFTV
253   FORMAT('SHIFTV',T15,G20.10)
      IF (ADMT) WRITE(3,250) NADM
250   FORMAT('ADM',T15,I5)
      IF (TOSI) THEN
         WRITE(3,245) PARAM1,PARAM2,PARAM3,PARAM4
245      FORMAT('TOSI',1X,4F15.8)
         IF (TOSIC6) WRITE(3,'(A,3F20.10)') 'TOSIC6 ',C6PP,C6MM,C6PM
         IF (TOSIPOL) WRITE(3,'(A,3F20.10)') 'TOSIPOL ',ALPHAP,ALPHAM,DAMP
      ELSE IF (WELCH) THEN
         WRITE(3,247) APP,AMM,APM,RHO,XQP,XQM,ALPHAP,ALPHAM
247      FORMAT('WELCH',1X,4F15.8,' +++'/'      ',4F15.8)
      ELSE IF (SIO2T) THEN
         WRITE(3,'(A,4G20.10)') 'SIO2 ', PARAM1,PARAM2,PARAM3,PARAM4
      ELSE IF (PARAM1.NE.0.0D0) THEN
         WRITE(3,260) PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7
260      FORMAT('PARAMS',1X,3F20.10,3X,'+++',/,7X,4F20.10)
      ENDIF
      IF (FTEST) WRITE(3,265) GALPHA, MALPHA1, MALPHA2
265   FORMAT('ALPHA',T15,3F20.10)
      IF (D5HT) WRITE(3,241) FD5H
241   FORMAT('D5H',T15,F15.5)
      IF (OHT) WRITE(3,242) FOH
242   FORMAT('OH',T15,F15.5)
      IF (IHT) WRITE(3,243) FIH
243   FORMAT('IH',T15,F15.5)
      IF (TDT) WRITE(3,244) FTD
244   FORMAT('TD',T15,F15.5)
      IF ((MUPDATE.NE.4).OR.(XMUPDATE.NE.4)) WRITE(3,'(A,T15,2I4)') 'UPDATES',MUPDATE,XMUPDATE
      IF ((DGUESS.NE.0.1D0).OR.(XDGUESS.NE.0.1D0)) WRITE(3,'(A,T15,2F20.10)') 'DGUESS',DGUESS,XDGUESS
      IF (TWOD) WRITE(3,248) 
248   FORMAT('2D')
      IF (TTM3T) WRITE(3,'(A)') 'TTM3'
      IF (BOWMANT) WRITE(3,'(A,I6,A)') 'BOWMANT ',BOWMANPES,TRIM(ADJUSTL(BOWMANDIR))
      IF (BSMIN) WRITE(3,'(A,F20.10)') 'BSMIN    ',GMAX
      IF (RKMIN) WRITE(3,'(A,F20.10)') 'RKMIN    ',GMAX
      IF (BFGSMINT) THEN
         WRITE(3,'(A)') 'BFGSMIN    '
         IF (GMAX.NE.0.001) WRITE(3,'(A,2F20.10)') 'BFGSCONV    ',GMAX
      ENDIF
      IF (NOIT) WRITE(3,'(A)') 'NOIT'
      IF (NORESET) WRITE(3,'(A)') 'NORESET'
      IF (HINDEX.GT.1) WRITE(3,'(A,I6)') 'INDEX ',HINDEX
      IF (BFGSTST) WRITE(3,'(A,2F20.10)') 'BFGSCONV    ',GMAX
      IF (BFGSTST.AND.NOHESS) WRITE(3,'(A,3I4,F20.10)') 'BFGSTS     ',NEVS,NBFGSMAX1,NBFGSMAX2,CEIG
      IF (BFGSTST.AND.(.NOT.NOHESS)) WRITE(3,'(A,3I6,F20.10,I6)') 'BFGSTS     ',NEVS,NBFGSMAX1,NBFGSMAX2,CEIG,NEVL
      IF (NOHESS) WRITE(3,'(A)') 'NOHESS'
      IF (CHECKINDEX) WRITE(3,'(A)') 'CHECKINDEX'
      IF (CHECKCONT) WRITE(3,'(A)') 'CHECKCONT'
      IF (BINARY) WRITE(3,'(A,I6,4F15.5)') 'BINARY ',NTYPEA,EPSAB,EPSBB,SIGAB,SIGBB
      IF (NHCHECK.NE.6) WRITE(3,'(A,I4)') 'AXIS      ',NHCHECK
      S=1.0D0
      IF ((ZSYM(NATOMS).EQ.'AR').AND.(ZSYM(1).NE.'CA')) S=3.4D0
      IF ((ZSYM(NATOMS).EQ.'CD'))  then
         WRITE(3,'(A,3G20.10)') 'CAPSID  ', CAPSRHO, EPS2, RAD
      else if (ANGLEAXIS) then
         write(3, '(A)') 'ANGLEAXIS'    ! Superfluous for CAPSID potential
      endif
      if (efield.ne.0.0d0) WRITE(3,'(A,G20.10)') 'EFIELD ', efield
      IF (NFAILMAX.NE.2) WRITE(3, '(A,I6)') 'MAXFAIL ', NFAILMAX
      if (SCORE_QUEUE) write(3,'(A)') 'SCORE_QUEUE'
      if (coldFusionLimit.ne.-1.0D6) WRITE(3,'(A,G20.10)') 'COLDFUSION ', coldFusionLimit
      IF (STOCKT) WRITE(3,'(A,2G20.10)') 'STOCK ', STOCKMU, STOCKLAMBDA
      IF (SIO2C6T) WRITE(3,'(A,3G20.10)') 'SIO2C6 ', C6PP,C6MM,C6PM
      IF (VARIABLES) THEN
         WRITE(3,275)
275      FORMAT('VARIABLES')
         WRITE(3,'(G20.10)') (Q(J1),J1=1,NATOMS)
C       
C       sf344> AMBER 9 interface
C
      ELSE IF ((AMBERT.OR.NABT).AND.DUMPSTRUCTURES) THEN
        CALL AMBERFINALIO(1,20,0,filthstr,filth,Q)
      ELSE IF (AMBER12T .AND. DUMPSTRUCTURES) THEN
         IF (FILTH == 0) THEN
! points.final.xyz
            OPEN(UNIT=20, FILE='points.final.xyz', STATUS='unknown')
            WRITE(20, '(I6)') NATOMS
            WRITE(20, '(A)') 'points.final'
            DO J1 = 1, NATOMS
               WRITE(20, '(A4,3F20.10)') AMBER12_ATOMS(J1) % NAME, Q(3*J1-2:3*J1)
            END DO
! points.final.rst
            CALL AMBER12_WRITE_RESTART(Q(:), 'points.final.rst', LEN('points.final.rst'))
! points.final.pdb
            CALL AMBER12_WRITE_PDB(Q(:), 'points.final.pdb', LEN('points.final.pdb'))
         ELSE
! points.final.xyz.<run_number>
            OPEN(UNIT=20, FILE='points.final.xyz.'//TRIM(ADJUSTL(FILTHSTR)), STATUS='unknown')
            WRITE(20, '(I6)') NATOMS
            WRITE(20, '(A)') 'points.final.'//TRIM(ADJUSTL(FILTHSTR))
            DO J1 = 1, NATOMS
               WRITE(20, '(A4,3F20.10)') AMBER12_ATOMS(J1) % NAME, Q(3*J1-2:3*J1)
            END DO
! points.final.rst.<run_number>
            CALL AMBER12_WRITE_RESTART(Q(:),
     &                                 'points.final.rst.'//TRIM(ADJUSTL(FILTHSTR)),
     &                                 LEN('points.final.rst.'//TRIM(ADJUSTL(FILTHSTR))))
! points.final.pdb.<run_number>
            CALL AMBER12_WRITE_PDB(Q(:),
     &                             'points.final.pdb.'//TRIM(ADJUSTL(FILTHSTR)),
     &                             LEN('points.final.pdb.'//TRIM(ADJUSTL(FILTHSTR))))
         END IF
         CLOSE(20)

!OPEP dump
      ELSE IF (OPEPT.AND.DUMPSTRUCTURES) THEN
         IF (FILTH == 0) THEN
            CALL OPEP_WRITE_PDB(NATOMS,Q(:),'points.final.pdb')
         ELSE
            CALL OPEP_WRITE_PDB(NATOMS,Q(:),'points.final.pdb.'//TRIM(ADJUSTL(FILTHSTR)))     
         END IF
         CLOSE(20)

      ELSE IF (VARIABLES) THEN
         WRITE(3,'(A)') 'VARIABLES'
         WRITE(3,'(G20.10)') (Q(J1),J1=1,NATOMS)
      ELSE IF (RINGPOLYMERT) THEN
         WRITE(3,'(A,I10,G20.10)') 'RINGPOLYMER  '//RPSYSTEM,RPIMAGES,RPBETA
         WRITE(3,'(G20.10)') (Q(J1),J1=1,NOPT)
      ELSE
         WRITE(3,270)
270      FORMAT('POINTS')
         WRITE(3,'(A2,3X,3G20.10)') (ZSYM(J1),Q(3*J1-2)*S,Q(3*J1-1)*S,Q(3*J1)*S,J1=1,NATOMS)
      ENDIF
      CLOSE(3)
      OPEN(UNIT=3,FILE=FNAMEF,STATUS='UNKNOWN')
      IF (VARIABLES) THEN
         WRITE(3,'(G25.15)') (Q(J1),J1=1,NATOMS)
      ELSE
         WRITE(3,'(3G25.15)') (Q(3*J1-2)*S,Q(3*J1-1)*S,Q(3*J1)*S,J1=1,NATOMS)
      ENDIF
      CLOSE(3)
      IF (ZSYM(NATOMS).EQ.'CD') THEN
         OPEN(UNIT=26,FILE='capsid.xyz',STATUS='UNKNOWN')
         WRITE(26,'(I6)') (NATOMS/2)*6
         WRITE(26,'(I5)') J1
         DO J2=1,NATOMS/2 
            CALL CAPSIDIO(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
     1                    Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),CAPSCOORDS,RAD,HEIGHT)
            DO J3=1,5
               WRITE(26,'(A4,3F20.10)') 'C1 ',CAPSCOORDS(3*(J3-1)+1),CAPSCOORDS(3*(J3-1)+2),
     1                                        CAPSCOORDS(3*(J3-1)+3)
            ENDDO
                 WRITE(26,'(A4,3F20.10)') 'C4 ',CAPSCOORDS(16),CAPSCOORDS(17),CAPSCOORDS(18)
         ENDDO
         CLOSE(26)
      ENDIF
      IF (STOCKT) THEN
         DSHIFT=0.5D0
         OPEN(UNIT=26,FILE='Stock.xyz',STATUS='UNKNOWN')
         WRITE(26,'(I6)') (NATOMS/2)
         WRITE(26,'(A)') ' '
         DO J2=1,(NATOMS/2) 
            WRITE(26,'(A,3G20.10,A,3G20.10)') 'LA  ',Q(3*J2-2),Q(3*J2-1),Q(3*J2), ' atom_vector ',
     &                                     SIN(Q(3*((NATOMS/2)+J2-1)+1))*COS(Q(3*((NATOMS/2)+J2-1)+2)),
     &                                     SIN(Q(3*((NATOMS/2)+J2-1)+1))*SIN(Q(3*((NATOMS/2)+J2-1)+2)),
     &                                     COS(Q(3*((NATOMS/2)+J2-1)+1))
         ENDDO
         CLOSE(26)
      ENDIF
      IF (RBAAT) THEN
         OPEN(UNIT=26,FILE='rbsites.xyz',STATUS='UNKNOWN')
         WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
         WRITE(26,'(A)') ' '
         CALL SITEPOS(Q,XS)
         IF (PTSTSTT) THEN
            DO J1 = 1, NATOMS/2
               DO J2 = 1, NRBSITES
                  J3 = ((J1-1)*NRBSITES + (J2-1))*3
                  IF (J2 == NRBSITES) THEN
                     WRITE(26,'(A5,1X,3F20.10)') 'O ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSEIF (J2 <= (NRBSITES-1)/2) THEN
                     WRITE(26,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSE
                     WRITE(26,'(A5,1X,3F20.10)') 'F ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF (PAHAT) THEN
            DO J1 = 1, NATOMS/2
               DO J2 = 1, NRBSITES
                  J3 = ((J1-1)*NRBSITES + (J2-1))*3
                  IF (J2 <= NCARBON) THEN
                     WRITE(26,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSE
                     WRITE(26,'(A5,1X,3F20.10)') 'H ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         CLOSE(26)
      ENDIF

      RETURN
      END
C
C  SUBROUTINE to convert capsid CofM and DV coordinates to penatgons.
C
      SUBROUTINE CAPSIDIO(X1, Y1, Z1, L1, M1, N1,COORDS,RAD,HEIGHT)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(18), HEIGHT, C2A1,
     2                 M1, L1, N1, ALPHA1, RAD, CA1, S1, C3A1,
     3                 NUM1, NUM2, NUM3, NUM4, NUM5, L12, M12, N12

C     HEIGHT=RAD/2.0D0
      NUM1=-(1.0D0+SQRT(5.0D0))/4.0D0
      NUM2=SQRT((5.0D0-SQRT(5.0D0))/2.0D0)/2.0D0
      NUM3=SQRT((5.0D0+SQRT(5.0D0))/2.0D0)/2.0D0
      NUM4=(SQRT(5.0D0)-1.0D0)/4.0D0
      NUM5=-(1.0D0+SQRT(5.0D0))/4.0D0

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=RAD*CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=RAD*(-ALPHA1/2+ALPHA1**3/24) ! bug
         C3A1=RAD*(-0.5D0+ALPHA1**2/24.0D0)
         S1=RAD*(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=RAD*(CA1-1.0D0)/ALPHA1**2
         S1=RAD*SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) =     c2a1 - c3a1*l12 + x1
      COORDS(2) =     -(c3a1*l1*m1) - n1*s1 + y1
      COORDS(3) =     -(c3a1*l1*n1) + m1*s1 + z1
      COORDS(4) =     c2a1*num4 - c3a1*l1*(m1*num3 + l1*num4) + n1*num3*s1 + x1
      COORDS(5) =     c2a1*num3 - c3a1*m1*(m1*num3 + l1*num4) - n1*num4*s1 + y1
      COORDS(6) =     -(c3a1*n1*(m1*num3 + l1*num4)) - l1*num3*s1 + m1*num4*s1 + z1
      COORDS(7) =     c2a1*num1 - c3a1*l1*(l1*num1 + m1*num2) + n1*num2*s1 + x1
      COORDS(8) = c2a1*num2 - c3a1*m1*(l1*num1 + m1*num2) - n1*num5*s1 + y1
      COORDS(9) = -(c3a1*n1*(l1*num1 + m1*num2)) + m1*num1*s1 - l1*num2*s1 + z1
      COORDS(10) = c2a1*num1 + c3a1*l1*(-(l1*num1) + m1*num2) - n1*num2*s1 + x1
      COORDS(11) = -(c2a1*num2) + c3a1*m1*(-(l1*num1) + m1*num2) - n1*num5*s1 + y1
      COORDS(12) = -(c3a1*l1*n1*num1) + c3a1*m1*n1*num2 + m1*num1*s1 + l1*num2*s1 + z1
      COORDS(13) = c2a1*num4 + c3a1*l1*(m1*num3 - l1*num4) - n1*num3*s1 + x1
      COORDS(14) = -(c2a1*num3) + c3a1*m1*(m1*num3 - l1*num4) - n1*num4*s1 + y1
      COORDS(15) = c3a1*n1*(m1*num3 - l1*num4) + l1*num3*s1 + m1*num4*s1 + z1
      COORDS(16)= (-(c3a1*l1*n1) - m1*s1 + 2*x1)/2.
      COORDS(17)= -(c3a1*m1*n1)/2. + (l1*s1)/2. + y1
      COORDS(18)= (c2a1 - c3a1*n12 + 2*z1)/2.
                
      RETURN
      END
