!
! Lump a minima.pdf file output by calc.Cv.f90 to merge quench bins
!
PROGRAM MINIMALUMP
IMPLICIT NONE
INTEGER NLUMP, NENTRIES, MINQBIN, MAXQBIN, NQBIN, NDOING, J1
CHARACTER(LEN=80) FILENAME, SDUMMY
DOUBLE PRECISION MAXEQ, MINEQ, DELTAE, XXPI, XXST, XPI, XST, DUMMY, NEWE

!
! Get the number of the bins to combine and the minima.pdf file name 
!

OPEN(UNIT=1,FILE='minima.lump.data',STATUS='OLD')
READ(1,*) NLUMP, FILENAME
CLOSE(1)
!
! This is the minima.pdf file format:
!
!#  bin  quench energy    ln(PI isomers/P)  PI isomers/P   PI isomers/total  ln(structures)    structures     structures/total
!     1    -133.6000       -0.6931472        0.5000000        0.3347172E-14     0.000000         1.000000        0.6687736E-14
!    13    -133.3168       -0.6931472        0.5000000        0.3347172E-14     0.000000         1.000000        0.6687736E-14
!    18    -133.1988        0.1887379E-14     1.000000        0.6694343E-14    0.1887379E-14     1.000000        0.6687736E-14

OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FILENAME)),STATUS='OLD')
READ(1,*) SDUMMY
READ(1,*) MINQBIN, MINEQ
NENTRIES=1
DO
   READ(1,*,END=321) MAXQBIN, MAXEQ
   NENTRIES=NENTRIES+1
ENDDO
321 CONTINUE
REWIND(1)

PRINT '(A,I10,A,2I10)','Number of quench bin entries=',NENTRIES,' quench bin index range: ',MINQBIN,MAXQBIN
DELTAE=(MAXEQ-MINEQ)/(MAXQBIN-MINQBIN)
PRINT '(A,G20.10,A,I10,A,G20.10)','Quench bin width=',DELTAE,' lumping groups of ',NLUMP,' to give width ',NLUMP*DELTAE

PRINT '(A)',' bin quench energy    ln(PI isomers/P)  PI isomers/P   ln(structures)    structures     '
READ(1,*) SDUMMY
XXPI=0.0D0; XXST=0.0D0
NDOING=MINQBIN
DO J1=1,NENTRIES
   READ(1,*,END=123) NQBIN, DUMMY, DUMMY, XPI, DUMMY, DUMMY, XST
   IF (NQBIN.GE.NDOING+NLUMP) THEN
      NEWE=MINEQ+(NDOING-MINQBIN)*DELTAE+NLUMP*DELTAE/2.0D0
      PRINT '(A,I10,A,G20.10)','move to new lumped bin for quench bin ',NQBIN,' lumped bin centre energy=',NEWE
      PRINT '(5G20.10)',NEWE,LOG(XXPI),XXPI,LOG(XXST),XXST
!
! Print entry for finished lumped bin and reinitialise
! Start next lumped bin from next non-zero entry - don't worry about empty ranges
!
      NDOING=NQBIN
      XXPI=0.0D0; XXST=0.0D0
   ENDIF
   XXPI=XXPI+XPI
   XXST=XXST+XST
ENDDO

123 CONTINUE
NEWE=MINEQ+(NDOING-MINQBIN)*DELTAE+NLUMP*DELTAE/2.0D0
PRINT '(A,I10,A,G20.10)','last lumped bin for quench bin ',NQBIN,' lumped bin centre energy=',NEWE
PRINT '(5G20.10)',NEWE,LOG(XXPI),XXPI,LOG(XXST),XXST

END PROGRAM MINIMALUMP
