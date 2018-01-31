SUBROUTINE DISPLACEMENT(iErr,nP,nC,nE,lPARAMS,lESIZE,lEPOS,lERMAT,lCPOS,lCRMAT)
!This subroutine takes a pair of ellipsoids, places one of the pair at the origin aligned with xyz axes.
!It then displaces one of them by a vector r, defined in the body axes of the first ellipsoid. 
!BOdy axes: a->x, b->y, c->Z
!It then rotates the second of the ellipsoids by three euler angles (a,b and gamma)(ZYZ convention)
!The resulting pair of ellipsoids are then stacked vertically above each other (z axis) spaced by an amount v
!until there are nE/2 pairs of ellipsoids in the cluster.  
!nC clusters are then stacked vertically with a displacement of  nE/2 * V

implicit none


      INTEGER :: nP, nC, nE
      INTEGER :: iErr
      DOUBLE PRECISION :: lPARAMS(1:nP)
      DOUBLE PRECISION :: lESIZE(1:nE,1:6), lEPOS(1:nE,1:3)
      DOUBLE PRECISION :: lCPOS(1:nC,1:3)
      DOUBLE PRECISION :: lERMAT(1:nE,1:3,1:3),lCRMAT(1:nC,1:3,1:3)
      DOUBLE PRECISION :: D(1:3)
      DOUBLE PRECISION :: V, a, b, g, aa, ab, ac, ra, rb, rc, ca, sa, cb, sb, cg, sg
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
      D(1)=lPARAMS(7)
      D(2)=lPARAMS(8)
      D(3)=lPARAMS(9)
      a=lPARAMS(10)*pi/DBLE(180)
      b=lPARAMS(11)*pi/DBLE(180)
      g=lPARAMS(12)*pi/DBLE(180)
      V=lPARAMS(13)
 
      ca=cos(a)
      sa=sin(a)
      cb=cos(b)
      sb=sin(b)
      cg=cos(g)
      sg=sin(g)

      !make sure there are an even number of ellipses in the cluster if not, then exit with an error code
      IF (mod(nE,2).EQ.1) THEN
       write(6,*) 'nE is not a multiple of two'
       iErr=1
       RETURN
      ENDIF

      DO curE=1,nE,2
         !first ellipsoid at origin
         lEPOS(curE,1)=0
         lEPOS(curE,2)=0
         lEPOS(curE,3)=0

         CALL EULER2RMAT(DBLE(0),DBLE(0),DBLE(0),lERMAT(curE,:,:))

         !second ellipsoid at D
         lEPOS(curE+1,:)=D

         CALL EULER2RMAT(a,b,g,lERMAT(curE+1,:,:))

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
END SUBROUTINE DISPLACEMENT
