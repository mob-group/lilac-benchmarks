SUBROUTINE DISTANCE(iErr,nP,nC,nE,lPARAMS,lESIZE,lEPOS,lERMAT,lCPOS,lCRMAT)
!This subroutine takes a pair of ellipsoids, places one of the pair at the origin and the other at a point d along the x axis.

!The ellipsoids are then oriented so that points of closest approach of their b-vectors (points p1 and p2 respectively) are a distance s apart, and also a distance u1 from the first ellipsoid and u2 from the second ellipsoid. 

!Each ellipsoid is then rotated about its b vector by angles g1 and g2 respectively.

!each cluster is then stacked along the z axis a distance V apart.

implicit none


      INTEGER :: nP, nC, nE
      INTEGER :: iErr
      DOUBLE PRECISION :: lPARAMS(1:nP)
      DOUBLE PRECISION :: lESIZE(1:nE,1:6), lEPOS(1:nE,1:3)
      DOUBLE PRECISION :: lCPOS(1:nC,1:3)
      DOUBLE PRECISION :: lERMAT(1:nE,1:3,1:3),lCRMAT(1:nC,1:3,1:3)
      DOUBLE PRECISION :: V, d, s, u1,u2, a, b1, b2, g1, g2, dg, aa, ab, ac, ra, rb, rc
      DOUBLE PRECISION, parameter :: pi=3.1415926
      INTEGER :: curE, curP, curC, numPairs, outI
  

! Initialise output error, zero good.
      iErr=0

!Copy params across
      aa=lPARAMS(1)
      ab=lPARAMS(2)
      ac=lPARAMS(3)
      ra=lPARAMS(4)
      rb=lPARAMS(5)
      rc=lPARAMS(6)
      d=lPARAMS(7)
      s=lPARAMS(8)
      u1=lPARAMS(9)
      u2=lPARAMS(10)
      g1=lPARAMS(11)*pi/DBLE(180)
      dg=lPARAMS(12)*pi/DBLE(180)
      V=lPARAMS(13)

      g2=g1-dg

      !make sure there are an even number of ellipses in the cluster if not, then exit with an error code
      IF (mod(nE,2).EQ.1) THEN
       write(6,*) 'nE is not a multiple of two'
       iErr=1
       RETURN
      ENDIF

      !compute angle parameters from distance parameters.        
      CALL DIST2ANGLE(d,s,u1,u2,b1,b2,a)

      DO curE=1,nE,2
 
         !gs are the same in both angle and distance definitions; use the computed angles to orient the vectors.
         CALL ANGLECLUSTER(DBLE(0.0),b1,g1,DBLE(0.0),lEPOS(curE,1:3),lERMAT(curE,1:3,1:3))
         CALL ANGLECLUSTER(a,b2,g2,d,lEPOS(curE+1,1:3),lERMAT(curE+1,1:3,1:3))
        
         !sort out sizes
         lESIZE(curE, 1) =aa
         lESIZE(curE, 2) =ab
         lESIZE(curE, 3) =ac
         lESIZE(curE, 4) =ra
         lESIZE(curE, 5) =rb
         lESIZE(curE, 6) =rc
         lESIZE(curE+1, 1) =aa
         lESIZE(curE+1, 2) =ab
         lESIZE(curE+1, 3) =ac
         lESIZE(curE+1, 4) =ra
         lESIZE(curE+1, 5) =rb
         lESIZE(curE+1, 6) =rc
      END DO

      numPairs=nE/2

      !Loop through each cluster
      DO curC=1,nC,1
         !compute cluster position
         lCPOS(curC,1)=0
         lCPOS(curC,2)=0
         lCPOS(curC,3)=DBLE((curC-1)*numPairs)*V
         !No rotation
         lCRMAT(curC,1,1)=1
         lCRMAT(curC,1,2)=0
         lCRMAT(curC,1,3)=0
         lCRMAT(curC,2,1)=0
         lCRMAT(curC,2,2)=1
         lCRMAT(curC,2,3)=0
         lCRMAT(curC,3,1)=0
         lCRMAT(curC,3,2)=0
         lCRMAT(curC,3,3)=1
     END DO

     RETURN
END SUBROUTINE DISTANCE
