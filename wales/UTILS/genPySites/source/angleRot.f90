SUBROUTINE ANGLEROT(iErr,nP,nC,nE,lPARAMS,lESIZE,lEPOS,lERMAT,lCPOS,lCRMAT)
!This subroutine takes a pair of ellipsoids, places one of the pair at the origin and rotates it about a vertical axis by angle beta1 and then about y' axis by angle gamma1 where the ' denotes the new y-axis after the first rotation.

!The second ellipsoid is placed a distance d along the x-axis and rotated by an angle beta2 about the z-axis and then about angle alpha about the x' axis and then about the y'' axis by angle gamma2.

!each cluster is then stacked along the z axis a distance V apart.

implicit none


      INTEGER :: nP, nC, nE
      INTEGER :: iErr
      DOUBLE PRECISION :: lPARAMS(1:nP)
      DOUBLE PRECISION :: lESIZE(1:nE,1:6), lEPOS(1:nE,1:3)
      DOUBLE PRECISION :: lCPOS(1:nC,1:3)
      DOUBLE PRECISION :: lERMAT(1:nE,1:3,1:3),lCRMAT(1:nC,1:3,1:3)
      DOUBLE PRECISION :: spacer,V, d, a, b1, b2, db, g1, g2, dg,  aa, ab, ac, ra, rb, rc, s, u1, u2, ca, cb, cg
      DOUBLE PRECISION :: ob1, ob2, odb, oa
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
      a=lPARAMS(8)*pi/DBLE(180)
      b1=lPARAMS(9)*pi/DBLE(180)
      db=lPARAMS(10)*pi/DBLE(180)
      g1=lPARAMS(11)*pi/DBLE(180)
      dg=lPARAMS(12)*pi/DBLE(180)
      spacer=lPARAMS(13)
      V=spacer*ac
      ca=lPARAMS(14)*pi/DBLE(180)
      cb=lPARAMS(15)*pi/DBLE(180)
      cg=lPARAMS(16)*pi/DBLE(180)

      !compute angles from delta angles
      b2=b1+db
      g2=g1+dg

      !make sure there are an even number of ellipses in the cluster if not, then exit with an error code
      IF (mod(nE,2).EQ.1) THEN
       write(6,*) 'nE is not a multiple of two'
       iErr=1
       RETURN
      ENDIF

      OPEN(9, FILE='ellipsoid.model.distances',STATUS='replace')
      write(9,*) 'd s u1 u2 b1 b2 db a'
      DO curE=1,nE,2
 
         lEPOS(curE,1)=-d/2.0D0
         lEPOS(curE,2)=DBLE(0)
         lEPOS(curE,3)=DBLE(0)
         lEPOS(curE+1,1)=d/2.0D0
         lEPOS(curE+1,2)=DBLE(0)
         lEPOS(curE+1,3)=DBLE(0)

         CALL EULER2RMAT(DBLE(0),b1,g1,lERMAT(curE,1:3,1:3))
         CALL EULER2RMAT(a,b2,g2,lERMAT(curE+1,1:3,1:3))
 
         CALL ANGLE2DIST(d,s,u1,u2,b1,b2,a) 

         !convert angles to degrees for output to file   
         ob1=b1*DBLE(180)/pi
         ob2=b2*DBLE(180)/pi
         odb=db*DBLE(180)/pi
         oa=a*DBLE(180)/pi

         write(6,*) 'd s u1 u2 b1 b2 db a'
         write(6,10),d,s,u1,u2,ob1,ob2,odb,oa
         write(9,10),d,s,u1,u2,ob1,ob2,odb,oa
10 FORMAT(8(1xf10.4))
 
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
      CLOSE (9)
      numPairs=nE/2

      !Loop through each cluster - create positions and body axes for the cluster frame.
      DO curC=1,nC,1
         !compute cluster position define a space curve as a function of parameter curC
         lCPOS(curC,1)=0
         lCPOS(curC,2)=0
         lCPOS(curC,3)=DBLE((curC-1)*numPairs)*V
 
        !Set up a sequence of the three euler angles over the parameter curC - defines the orientation of the cluster at each point on the space curve
         !ca=!DBLE(0)!(-pi/DBLE(2.0)) !(DBLE(curC)- DBLE(1))*DBLE(2.0)*pi/(DBLE(nC) - DBLE(1.0))
         !cb=!DBLE(0)!(DBLE(curC)- DBLE(1))*DBLE(2.0)*pi/(DBLE(nC) - DBLE(1.0))
         !cg=!DBLE(0)!-cb !DBLE(pi/DBLE(2.0)) !(DBLE(curC)- DBLE(1))*DBLE(2.0)*pi/(DBLE(nC) - DBLE(1.0))

         CALL EULER2RMAT(ca,cb,cg,lCRMAT(curC,:,:))

     END DO

     RETURN
END SUBROUTINE ANGLEROT
