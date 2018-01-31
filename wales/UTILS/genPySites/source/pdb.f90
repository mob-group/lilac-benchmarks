SUBROUTINE PDB(iErr,nP,nC,nE,lPARAMS,lESIZE,lEPOS,lERMAT,lCPOS,lCRMAT)
!This subroutine takes a pair of ellipsoids, places one of the pair at the origin and rotates it about a vertical axis by angle beta1 and then about y' axis by angle gamma1 where the ' denotes the new y-axis after the first rotation.

!The second ellipsoid is placed a distance d along the x-axis and rotated by an angle beta2 about the z-axis and then about angle alpha about the x' axis and then about the y'' axis by angle gamma2.

!each cluster is then placed at the position and orientation specified in the file coords.in

implicit none


      INTEGER :: nP, nC, nE
      INTEGER :: iErr
      DOUBLE PRECISION :: lPARAMS(1:nP)
      DOUBLE PRECISION :: lESIZE(1:nE,1:6), lEPOS(1:nE,1:3)
      DOUBLE PRECISION :: lCPOS(1:nC,1:3), lP(1:nC,1:3)
      DOUBLE PRECISION :: lERMAT(1:nE,1:3,1:3),lCRMAT(1:nC,1:3,1:3)
      DOUBLE PRECISION :: spacer,V, d, a, b1, b2, db, g1, g2, dg,  aa, ab, ac, ra, rb, rc, s, u1, u2, ca, cb, cg
      DOUBLE PRECISION :: ob1, ob2, odb, oa
      double precision :: q(1:4)
      double precision :: thetah
      DOUBLE PRECISION, parameter :: pi=3.1415926
      INTEGER :: curE, curP, curC, numPairs, outI
      double precision :: qN(1:4), COG(1:3)
      double precision :: sq(1:4), tmp1, tmp2
      integer i
      double precision :: epsilon
      epsilon=0.0100570

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

      OPEN(9, FILE='ellipsoid.model.distances',STATUS='replace')
      write(9,*) 'd s u1 u2 b1 b2 db a'
      DO curE=1,nE,2
 
         lEPOS(curE,:)=0
         lEPOS(curE+1,1)=d
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


      !open the coords file 
      open(9,FILE='coords.in',STATUS='old')

      !read in the cluster positions
      DO curC=1,nC,1
         !read in cluster position 
         READ(9,156) lCPOS(curC,1),lCPOS(curC,2),lCPOS(curC,3)
156      FORMAT(F8.3,1X,F8.3,1X,F8.3)
      END DO 

      DO curC=1,nC,1
         !read in cluster orientations from file
         READ(9,156) lP(curc,1),lP(curC,2),lP(curC,3)

         ! convert angle axis to quaternion
         thetah = 0.5d0 * SQRT(lP(curC,1)**2 + lP(curC,2)**2 + lP(curC,3)**2)
         q(1) = cos(thetah)

         ! do linear expansion for small epsilon
         if(thetah < epsilon) then
            q(2:4) = 0.5d0 * lP(curC,:)
         else
            q(2:4) = 0.5d0 * sin(thetah) * lP(curC,:) / thetah
         endif
         ! make sure to have normal form
         if(q(1) < 0d0) q = -q

         !Now convert quaternion to matrix
         qN = q / sqrt(dot_product(q,q))

         do i=1,4
            sq(i) = qN(i)*qN(i)
         enddo

         lCRMAT(curC,1,1) = ( sq(2) - sq(3) - sq(4) + sq(1))
         lCRMAT(curC,2,2) = (-sq(2) + sq(3) - sq(4) + sq(1))
         lCRMAT(curC,3,3) = (-sq(2) - sq(3) + sq(4) + sq(1))

         tmp1 = qN(2)*qN(3)
         tmp2 = qN(1)*qN(4)
         lCRMAT(curC,2,1) = 2.0d0 * (tmp1 + tmp2)
         lCRMAT(curC,1,2) = 2.0d0 * (tmp1 - tmp2)

         tmp1 = qN(2)*qN(4)
         tmp2 = qN(3)*qN(1)
         lCRMAT(curC,3,1) = 2.0d0 * (tmp1 - tmp2)
         lCRMAT(curC,1,3) = 2.0d0 * (tmp1 + tmp2)
         tmp1 = qN(3)*qN(4)
         tmp2 = qN(1)*qN(2)
         lCRMAT(curC,3,2) = 2.0d0 * (tmp1 + tmp2)
         lCRMAT(curC,2,3) = 2.0d0 * (tmp1 - tmp2)

      END DO 
      CLOSE (9)
     RETURN
END SUBROUTINE PDB
