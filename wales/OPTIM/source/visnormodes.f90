SUBROUTINE VISNORMODES(COORDS)
!----------------------------------------------------------------------!
!Program to apply the non mass-weighted hessian eigenvectors from      !
!a vector.dump file (OPTIM) to the coordinates in points.final.xyz     !
!                                                                      !
!The output is a set of nmode.X.xyz files where X is the vector number !
!and the numbering goes from soft->hard modes.                         !
!                                                                      !
!Each nmode.X file contains 11 'frames' where the vector has been      !
!applied -5,-4,-3,-2,-1,0,1,2,3,4 and 5 times. The first column is the !
!atom types, then the x,y and z for that frame followed by the         !
!components of the vector itself which is the same for each frame. The !
!seventh column contains the scaled atomic displacement (0->1) that can!
!be used to shade the atoms in VMD.
!                                                                      !
!To obtain the vector.dump file, you need to do an OPTIM run with an   !
!odata file containing at least the following keywords: SEARCH 0,      !
!STEPS <enough to converge>, DUMPVECTOR ALLVECTORS, DUMPSTRUCTURES and !
!NAB coords.inpcrd inpcrd. The cutoff in min.in should also be 999.    !
!NAB is used instead of AMBER9 to call analytical 2nd derivatives.     !
!----------------------------------------------------------------------!
USE commons, only : NRBSITES, NATOMS
USE key, only : MULTISITEPYT, NTSITES


IMPLICIT NONE

! Variable declarations
INTEGER :: J1, J2, J3, FILENUM, NUMPRINT, PRINTCOUNTER, NBINS, BINDEX
INTEGER, ALLOCATABLE :: BINS(:)
DOUBLE PRECISION :: SCALING, VECSCALE, MAXMAG, LOWB, UPB, BININC, COEF2, COEF4, COORDS(3*NATOMS), X(3*NATOMS)
CHARACTER (LEN=20) :: DUMMY, OUTNUM, OUTNAME, MODENUM, MODENAME
DOUBLE PRECISION, ALLOCATABLE :: EVALUE(:),MODES(:,:),MAG(:),PARATIO(:),RBCOORDS(:)
CHARACTER (LEN=4), ALLOCATABLE :: TYPES(:)
LOGICAL :: RESTRICTPRINT

! SETUP 
! First, each eiganvector is divided by SCALING to prevent too much motion
SCALING=2.0
! The scaling for the VMD plotted vectors defined in modefile is set next
VECSCALE=2.0
! How many modes do you want to print? Set below to true to enable
RESTRICTPRINT=.FALSE.
NUMPRINT=20
! How many bins do you want in the displacement histogram?
NBINS=50
! STEP 1: Read in the points.final.xyz file
!OPEN(UNIT=123,FILE='points.final.xyz',STATUS='OLD')
! How many atoms are we dealing with
!READ(123,*) NATOMS
! Allocate the size of the atom types, coordinates, normal mode energy and
! vector components arrays
ALLOCATE(TYPES(NATOMS))
ALLOCATE(MAG(NATOMS))
ALLOCATE(BINS(NBINS))
!ALLOCATE(COORDS(3*NATOMS))
ALLOCATE(EVALUE(3*NATOMS-6))
ALLOCATE(PARATIO(3*NATOMS-6))
ALLOCATE(MODES((3*NATOMS-6),3*NATOMS))

ALLOCATE(RBCOORDS(3*NTSITES))

! Read the file name into the dummy variable cause we don't care about it!
!READ(123,*) DUMMY
! Read in coordinates
!DO J1=1,NATOMS
!        READ(123,*) TYPES(J1), COORDS((3*J1)-2),COORDS((3*J1)-1),COORDS((3*J1))
!ENDDO
!CLOSE(UNIT=123)
TYPES(:)='O'
! STEP 2: Read in the vector.dump file
! Note. there are (3*NATOMS-6) hessian eigenvectors for non-linear molecules
!OPEN(UNIT=1233,FILE='vector.dump',STATUS='UNKNOWN')
! vector.dump is unit 44 in OPTIM
REWIND(44)
DO J1=1,(3*NATOMS-6)
        READ(44,*) EVALUE(J1)
! For each mode, read in the vector components per atom
        DO J2=1,NATOMS
                READ(44,*) MODES(J1,(3*J2)-2),MODES(J1,(3*J2)-1),MODES(J1,(3*J2))
        ENDDO
ENDDO
!CLOSE(UNIT=1233)
! STEP 3: Produce the output
! Produce plot of mode index vs eiganvalue
OPEN(UNIT=123,FILE='evaluevsmode.dat',STATUS='UNKNOWN')
DO J1=1,(3*NATOMS-6)
        WRITE(123,'(I5,F20.10)') J1, EVALUE(J1)
ENDDO
CLOSE(123)
! Open file to load all modes into VMD at once
OPEN(UNIT=124,FILE='modefile',STATUS='UNKNOWN')
! How many are we going to print? Depends if RESTRICTPRINT is set!
PRINTCOUNTER = (3*NATOMS-6)
IF (RESTRICTPRINT) PRINTCOUNTER = NUMPRINT
! The outer loop processes each hessian eigenvector
DO J1=1,PRINTCOUNTER
! This reverses the order of the modes to soft->hard (by changing file numbers)
! It seems to not be needed any more! Strange :/
!        FILENUM=((3*NATOMS-6)-(J1-1))
         FILENUM=J1
! Trick to produce the output file name as a string
        WRITE(OUTNUM,*) FILENUM
        OUTNAME='nmode.'//TRIM(ADJUSTL(OUTNUM))//'.xyz'
        OPEN(UNIT=33,FILE=OUTNAME,STATUS='UNKNOWN')
! Setup for entry in modefile (we want the entries from 1->(3*NATOMS-6))
        WRITE(MODENUM,*) J1 
        MODENAME='nmode.'//TRIM(ADJUSTL(MODENUM))//'.xyz'
! Write an entry into modefile to load this eigenvector into VMD after
        WRITE(124,'(2A,F3.1)') 'initxyz ',MODENAME//' ',VECSCALE
! Before I apply the vector, I need to work out the magnitude of the
! displacement for each atom and scale it between 0->1. This will allow me to
! colour the atoms by their displacement in the mode/vector. First, set the 
! current maximum displacement to zero, empty all the bins and set the coefficients 
! used in the participation ratio to zero.
        MAXMAG=0.0D0
        BINS(:)=0.0D0
        COEF2=0.0D0
        COEF4=0.0D0
! Now, loop over each atom, calculating the displacement and comparing it to the
! current max (MAXMAG), replacing as necessary.
        DO J3=1,NATOMS
                MAG(J3)=DSQRT(MODES(J1,3*J3-2)**2+MODES(J1,3*J3-1)**2+MODES(J1,3*J3)**2)
                IF (MAG(J3).GT.MAXMAG) MAXMAG=MAG(J3)
! While we're looping over NATOMS, why not calculate the participation ratio
! (measure of mode localisation) for the current mode
        COEF2=COEF2+(MODES(J1,3*J3-2)**2+MODES(J1,3*J3-1)**2+MODES(J1,3*J3)**2)
        COEF4=COEF4+(MODES(J1,3*J3-2)**4+MODES(J1,3*J3-1)**4+MODES(J1,3*J3)**4)
        ENDDO
! Calculate the participation ratio for mode J1
        PARATIO(J1)=(COEF2**2/(COEF4*NATOMS))
! Construct the magnitude histogram for mode J1 (cannot combine with above loop
! as we need MAXMAG to do the binning)
        DO J3=1,NATOMS
                LOWB=0.0D0
                UPB=1.0D0/NBINS
                BININC=1.0D0/NBINS
                DO BINDEX=1,NBINS
! The IF statement does the 'binning' is the value falls inside the bin limits
                        IF ((MAG(J3)/MAXMAG.GT.LOWB).AND.(MAG(J3)/MAXMAG.LE.UPB)) THEN 
                                BINS(BINDEX)=BINS(BINDEX)+1
                        ENDIF
                        LOWB=LOWB+BININC
                        UPB=UPB+BININC
                ENDDO
        ENDDO
! Output magnitude histogram for mode J1
OPEN(UNIT=123,FILE='mag_histo.'//TRIM(ADJUSTL(MODENUM)),STATUS='UNKNOWN')
LOWB=0.0D0
DO J3=1,NBINS
WRITE(123,'(F20.10,I5)') LOWB,BINS(J3)
LOWB=LOWB+BININC
ENDDO
CLOSE(123)
! Now I have the scaled magnitudes, I can apply the vector and print some output
! The inner loop applies the vector multiple times from -5 to +5
    IF(MULTISITEPYT) THEN

        DO J2=-5,5
           X(:)=COORDS(:)+J2*(MODES(J1,:))/SCALING
           CALL SITEPOS(X,RBCOORDS)
                WRITE(33,'(I6)') NTSITES 
                WRITE(33,'(F20.10)') EVALUE(J1) 
!                DO J3=1,NATOMS/2
!                        WRITE(33,'(A4,7F20.10)') TYPES(J3), (COORDS(3*J3-2)+(J2*(MODES(J1,3*J3-2))/SCALING)),(COORDS(3*J3-1)+(J2*(MODES(J1,3*J3-1))/SCALING)),&
!                                               &(COORDS(3*J3)+(J2*(MODES(J1,3*J3))/SCALING)),MODES(J1,3*J3-2),MODES(J1,3*J3-1),MODES(J1,3*J3),MAG(J3)/MAXMAG
!                ENDDO
                DO J3=1,NTSITES
                        WRITE(33,'(A4,7F20.10)') TYPES(1), RBCOORDS(3*J3-2), RBCOORDS(3*J3-1), RBCOORDS(3*J3)
                ENDDO
        ENDDO

    ELSE

        DO J2=-5,5
                WRITE(33,'(I6)') NATOMS
                WRITE(33,'(F20.10)') EVALUE(J1) 
! The inner-inner loop prints the line for each atom
                DO J3=1,NATOMS
                   WRITE(33,'(A4,7F20.10)') TYPES(J3), (COORDS(3*J3-2)+(J2*(MODES(J1,3*J3-2))/SCALING)), &
  &                                                    (COORDS(3*J3-1)+(J2*(MODES(J1,3*J3-1))/SCALING)),&
  &                                                    (COORDS(3*J3)+(J2*(MODES(J1,3*J3))/SCALING)), &
  &                     MODES(J1,3*J3-2),MODES(J1,3*J3-1),MODES(J1,3*J3),MAG(J3)/MAXMAG
                ENDDO
        ENDDO
    END IF  
      CLOSE(UNIT=33)
ENDDO
CLOSE(UNIT=124)
! Output participation ratio vs mode index and vs eiganvalue files
OPEN(UNIT=124,FILE='partratiovsmode.dat',STATUS='UNKNOWN')
OPEN(UNIT=125,FILE='partratiovsevalue.dat',STATUS='UNKNOWN')
DO J1=1,3*NATOMS-6
        WRITE(124,'(I5,F20.10)') J1, PARATIO(J1)
        WRITE(125,'(2F20.10)') EVALUE(J1), PARATIO(J1)
ENDDO
CLOSE(124)
CLOSE(125)
END SUBROUTINE VISNORMODES
