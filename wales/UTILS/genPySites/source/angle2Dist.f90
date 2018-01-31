      SUBROUTINE ANGLE2DIST(ld,ls,lu1,lu2,lb1,lb2,la)
      !converts the three angles lb1, lb2 and la into three distances u1, u2 and s. the units will be whatever the units of ld are.
      !assume angles are input as radians

      IMPLICIT NONE

      DOUBLE PRECISION :: ld,ls,lu1,lu2,lb1,lb2,la,pi,ldb
      DOUBLE PRECISION :: WRAP
      DOUBLE PRECISION :: dummyVal
      pi=3.141592654

      write(6,*) 'Converting Angles to Distances'

      !compute delta_B before wrapping
      ldb=lb1-lb2

10 FORMAT('la=',f5.2,' lb1=',f5.2,' lb2=',f5.2,' ldb=',f5.2,' d= ',f5.2)
      !wrap angles so they are between -pi/2 and pi/2
      la=WRAP(la,pi,- pi/DBLE(2.0))
      lb1=WRAP(lb1,pi,-pi/DBLE(2.0))
      lb2=WRAP(lb2,pi,-pi/DBLE(2.0))
      ldb=WRAP(ldb,pi,-pi/DBLE(2.0))
      write(6,10)la,lb1,lb2,ldb,ld

      !check for zero condition for alpha
      IF (ABS(la)<DBLE(1e-6)) THEN
       write (6,*) 'Coplanar: Zero alpha'
       ls=0
       !Check for parallel condition
       IF (ABS(ldb)<DBLE(1e-6)) THEN
          write(6,*) 'Parallel ellipsoids. So u1 and u2 are infinite.'
          dummyVal=0;
          lu1=1/dummyVal
          lu2=1/dummyVal
       ELSE
          write (6,*) 'not parallel'
          !Compute u1 from angle data
          lu1 = ld*(SIN(lb1) -  (COS(la)**2)*COS(ldb)*SIN(lb1 - ldb))/(-1 + (COS(la)**2)*COS(ldb)**2)
          lu2 = ld*(COS(la)*COS(lb1)*SIN(ldb))/(-1 + (COS(la)**2)*COS(ldb)**2)
       END IF
      ELSE
        write (6,*) 'not coplanar: alpha nonzero'

       !Check for perfectly parallel condition
       IF (ABS(ldb)<DBLE(1e-16)) THEN
          write (6,*) 'db=0'
          lu1=DBLE(0.0)
          lu2=DBLE(0.0)
          ls=ld*COS(lb1)
        ELSE
          write (6,*) 'not parallel'
          !All special cases accounted for. (I think!)
          !Compute u1 from angle data
          lu1 = ld*(SIN(lb1) -  (COS(la)**2)*COS(ldb)*SIN(lb1 - ldb))/(- 1 + (COS(la)**2)*COS(ldb)**2)
          lu2 = ld*(COS(la)*COS(lb1)*SIN(ldb))/ (- 1 + (COS(la)**2)*COS(ldb)**2 )
          ls  = ld*COS(lb1)/(  SQRT( DBLE(1.0) + ( SIN(ldb)**2 * COS(la)**2 /SIN(la)**2 ) ) )
        END IF
      END IF

      return
      END SUBROUTINE ANGLE2DIST
