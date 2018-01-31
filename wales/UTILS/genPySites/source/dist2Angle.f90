      SUBROUTINE DIST2ANGLE(ld,ls,lu1,lu2,lb1,lb2,la)
      !for a given d, computes the angles b1, b2 and a which are consistent with the three input lengths: u1, u2 and s. Returns the result in radians.
      IMPLICIT NONE

      DOUBLE PRECISION :: ld,ls,lu1,lu2,lb1,lb2,la
      DOUBLE PRECISION :: pi,g,h,k,gSq,hSq,kSq
      DOUBLE PRECISION :: tanAlphaSq1, tanAlphaSq2, tanAlpha1a, tanAlpha1b, tanAlpha2a, tanAlpha2b
      DOUBLE PRECISION :: a1a, a1b, a2a, a2b, b21a, b21b, b22a, b22b
      INTEGER :: alphaSolutionExists, alphaSolution1Exists, alphaSolution2Exists

      pi=3.141592654

      !assume that all lengths are given a zero or positive.

      !If s is zero (lines are not skew) then the two ellipsoids are co-planar. Therefore la1 must be zero.
      IF (ls < DBLE(1e-16)) THEN
          la=DBLE(0.0)

          !With alpha zero, the two equations for beta and delta beta are quite simple. Check for some impossible situations.
          IF (ABS(lu2**2 + ls**2 -lu1**2 -ld**2)>DBLE(2.0)*lu1*ld) THEN
            Write(6,*) 'No solution found for angle b1 where s is zero. The magnitude of an argument for ASIN > 1. Increase d or u1 and/or decrease u2. I assume you want s=0 in this case, or else you could increase s.'
            STOP
          ELSE
            !safe to compute lb1.
            lb1=ASIN((lu2**2  - lu1**2 - ld**2)/(2.0*lu1*ld))

            !perform check for lb2 with alpha is zero.
          IF (ABS(lu2**2 + ld**2 - lu1**2 - ls**2) > DBLE(2.0) * lu2 * ld)  THEN
            Write(6,*) 'No solution found for angle b2 in a situation where s is zero. The magnitude of an argument for ASIN > 1. Increase d or u2 and/or decrease u1. I assume you want s=0 in this case, or else you could increase s.'
              STOP
            ELSE
              lb2=ASIN(lu2**2 + ld**2 - lu1**2)/(2 * lu2 * ld)
            END IF !check for lb2 with la zero
          END IF !check for lb1 with la zero
      ELSE !check for la zero.

          !alpha not zero. This is where the fun begins.

          !Compute some intermediate values that make the equations and restrictions on solutions more simple to understand, allegedly
          g=(-ld**2 + ls**2 - lu1**2 + lu2**2)/(DBLE(2.0)* ld* lu1)
          gSq=g*g

          h=(ld**2 - ls**2 - lu1**2 + lu2**2)/(DBLE(2.0)*ld*lu2)
          hSq=h*h

          kSq=(ld**2 * gSq + ls**2 - ld**2)/ls**2
          !The following check arises from making sure that SQRT(kSq) is not complex, which prohibits a sensible value for K and hence alpha.
          !IF ABS(ld*ld - ls*ls - lu1 *lu1 - lu2*lu2)<2*SQRT[u1*u2] THEN
            !Write(6,*) 'No solution found for alpha. decrease d, or increase s, u1 or u2'
            !STOP
          !ELSE
            !k=SQRT(kSq)

            !The following monster equation for (tan(a))**2 is found by eliminating db 
            !from the following two equations which arise in turn from eliminating b1 
            !from the 3 equations (used in ANGLE2DIST) which arise from considering the 
            !vector constraints on the point of closest approach.

            !Cot[a]*Sin[db] == k 
            !g Cos[db] - Sqrt(1-g^2) Sin [db] ==  h Sec[a]

            !where the constants defined above simply make the equation readable. They 
            !may have physical significance based on which spheres and solutions are allowed. 
            !Eliminating db from these two equations yields a solution for alpha which 
            !can be expressed as follows which is taken straight from mathematica. I'm 
            !sure it can be simplified but what's the point. we're computing numerically 
            !anyway. This program runs v quickly and is not iterative. So can make it 
            !expensive if we want to for now while we test.  There are two solutions for 
            !Tan(a)^2 and therefore four solutions for Tan Alpha
          tanAlphaSq1= - (gSq - hSq)**2/ ( (hSq*hSq) - (hSq*kSq) + gSq*(- kSq + hSq*(DBLE(2.0)*kSq - DBLE(1.0))) +  DBLE(2.0)*SQRT( gSq*hSq*kSq*(gSq - DBLE(1.0))*(- gSq + hSq + (hSq - DBLE(1.0))*kSq)) )

          tanAlphaSq2= - (gSq - hSq)**2/ ( (hSq*hSq) - (hSq*kSq) + gSq*(- kSq + hSq*(DBLE(2.0)*kSq - DBLE(1.0))) -  DBLE(2.0)*SQRT( gSq*hSq*kSq*(gSq - DBLE(1.0))*(- gSq + hSq + (hSq - DBLE(1.0))*kSq)) )

          alphaSolutionExists=0
          alphaSolution1Exists=0
          alphaSolution2Exists=0
          IF (tanAlphaSq1>0) THEN
            alphaSolutionExists=1
            alphaSolution1Exists=1
            tanAlpha1a=SQRT(tanAlphaSq1)
            a1a=ATAN(tanAlpha1a)
            tanAlpha1b=-SQRT(tanAlphaSq1)
            a1b=ATAN(tanAlpha1b)
          END IF

          IF (tanAlphaSq1>0) THEN
            tanAlpha2a=SQRT(tanAlphaSq2)
            a2a=ATAN(tanAlpha2a)
            tanAlpha2b=-SQRT(tanAlphaSq2)
            a2b=ATAN(tanAlpha2b)
            alphaSolutionExists=1
            alphaSolution2Exists=1
          END IF

          IF (alphaSolutionExists==0) THEN
            Write(6,*) 'No solution found for alpha. Not sure what to suggest. Eqn is too crazy. Play about. Maybe increase d, decrease s, u1 or u2 based on k being complex. Might not be that reason though.'
            STOP
          ELSE
            !With alpha not zero compute b1 and b2
            IF (ABS(lu2**2 + ls**2 -lu1**2 -ld**2)>DBLE(2.0)*lu1*ld) THEN
              Write(6,*) 'No solution found for angle b1 alpha is non-zero. The magnitude of an argument for ASIN > 1. Increase d or u1 and/or decrease u2 or s.'
              STOP
            ELSE
              !safe to compute lb1.
              lb1=ASIN((lu2**2  - lu1**2 - ld**2)/(2.0*lu1*ld))

              !perform check for lb2 with alpha .
              IF (ABS(lu2**2 + ld**2 -lu1**2 - ls**2)>DBLE(2.0)*lu2*ld)  THEN
                Write(6,*) 'No solution found for angle b2 in a situation where s is zero. The magnitude of an argument for ASIN > 1. Increase d or u2 and/or decrease u1. I assume you want s=0 in this case, or else you could increase s.'
                STOP
              ELSE
                !Four solutions for b2 corresponding to four solutions for alpha. The solution that is picked will the first solution from the last set that exists.
                IF (alphaSolution1Exists) THEN
                 b21a=ASIN(lu2**2 + ld**2 - lu1**2 -ls**2)/(2 * lu2 * ld*Cos(a1a))
                 b21b=ASIN(lu2**2 + ld**2 - lu1**2 -ls**2)/(2 * lu2 * ld*Cos(a1b))
                 lb2=b21a
                 la=a1a
                 write (6,*) 'Solution set 1: (',a1a,',',b21a,') (',a1b,',',b21b,') '
                END IF
                IF (alphaSolution2Exists) THEN
                 b22a=ASIN(lu2**2 + ld**2 - lu1**2 -ls**2)/(2 * lu2 * ld*Cos(a2a))
                 b22b=ASIN(lu2**2 + ld**2 - lu1**2 -ls**2)/(2 * lu2 * ld*Cos(a2b))
                 lb2=b22a
                 la=a2a
                 write (6,*) 'Solution set 2: (',a2a,',',b22a,') (',a2b,',',b22b,')'
                END IF
              END IF !check for lb2 with la non zero and solutions exists
            END IF !check for lb1 with la non zero and solutions exists
          END IF !check for alpha solution exists
      END IF !check for non zero alpha

      !Convert solutions to degrees. Always get to this point because on errors program crashes out.
      !la=180*la/pi
      !lb1=180*lb1/pi
      !lb2=180*lb2/pi

      return
      END SUBROUTINE DIST2ANGLE
