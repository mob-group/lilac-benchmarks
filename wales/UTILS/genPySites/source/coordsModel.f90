SUBROUTINE COORDSMODEL(iErr,nP,nC,nE,lPARAMS,lESIZE,lEPOS,lERMAT,lCPOS,lCRMAT)
!This subroutine takes a pair of ellipsoids, places one of the pair at the origin and rotates it about a vertical axis by angle beta1 and then about y' axis by angle gamma1 where the ' denotes the new y-axis after the first rotation.

!The second ellipsoid is placed a distance d along the x-axis and rotated by an angle beta2 about the z-axis and then about angle alpha about the x' axis and then about the y'' axis by angle gamma2.

!each cluster is then stacked along the z axis a distance V apart.

      implicit none


      INTEGER :: nP, nC, nE
      INTEGER :: iErr
      DOUBLE PRECISION :: lPARAMS(1:nP)
      DOUBLE PRECISION :: lESIZE(1:nE,1:6), lEPOS(1:nE,1:3)
      DOUBLE PRECISION :: lCPOS(1:nC,1:3), POS(1:3), P(1:3)
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

      !compute angles from delta angles
      b2=b1+db
      g2=g1+dg

      !make sure there are an even number of ellipses in the cluster if not, then exit with an error code
      IF (mod(nE,2).EQ.1) THEN
         write(6,*) 'nE is not a multiple of two'
         iErr=1
         RETURN
      ENDIF

      DO curE=1,nE,2
 
         lEPOS(curE,1)=-d/2.0D0
         lEPOS(curE,2)=DBLE(0)
         lEPOS(curE,3)=DBLE(0)
         lEPOS(curE+1,1)=d/2.0D0
         lEPOS(curE+1,2)=DBLE(0)
         lEPOS(curE+1,3)=DBLE(0)

         CALL EULER2RMAT(DBLE(0),b1,g1,lERMAT(curE,1:3,1:3))
         CALL EULER2RMAT(a,b2,g2,lERMAT(curE+1,1:3,1:3))
 
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

      OPEN(9, FILE='coordsInp',STATUS='old')
      !Loop through each cluster - create positions and body axes for the cluster frame.
      DO curC=1,nC,1
         read(9,*) POS
         !compute cluster position define a space curve as a function of parameter curC
         lCPOS(curC,:)=POS
      END DO
      DO curC=1,nC,1
         read(9,*) P
         CALL AA2RMAT(P,lCRMAT(curC,:,:))
      END DO
      CLOSE(9)
      RETURN
END SUBROUTINE COORDSMODEL
