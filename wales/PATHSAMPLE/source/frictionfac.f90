!
! Calculate the correction factor for the TST canonical rate constant accoring to the
! formulation of Berezhkovsii, Pollak, Zitserman, JCP, 97, 2422, 1992
!
! GAMMAFRICTION is the phenomenological friction coefficient, which must be read in the
! correct frequency units (reduced or otherwise).
!
! NEGEVALUE is the unique negative Hessian eigenvalue. The square root of the absolute
! value is the modulus of the corresponding angular frequency, omega, and we have to divide 
! this by 2*pi to get frequency units. This cancels a factor of two in the denominators and
! puts a pi on top.
!
DOUBLE PRECISION FUNCTION FRICTIONFAC(NEGEVALUE)
USE COMMONS,ONLY : GAMMAFRICTION
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
DOUBLE PRECISION DUMMY, NEGEVALUE

!
! DUMMY is actually the friction factor divided by 2 times the modulus of the imaginary frequency
!
DUMMY=GAMMAFRICTION*PI/SQRT(ABS(NEGEVALUE))

!
! For large values of DUMMY FRICTIONFAC tends to 1/(2*DUMMY)
! The error is 1.25D-7 for DUMMY=1.0D2
!
IF (DUMMY.GT.1.0D2) THEN
   FRICTIONFAC=1.0D0/(2*DUMMY)
ELSE
   FRICTIONFAC=SQRT(1+DUMMY**2)-DUMMY
ENDIF

END FUNCTION FRICTIONFAC
