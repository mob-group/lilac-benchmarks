      SUBROUTINE VERTICALHELIX(iErr,nP,nC,nE,lPARAMS,lESIZE,lEPOS,lERMAT,lCPOS,lCRMAT)
!This routine arranges nE ellipsoids in a cluster with nS sheets and nL layers per cluster. 
!An error is thrown if nS*nL.ne.nE.
!Each ellipsoid belongs to one of nS different helices. 
!Each helix is defined by the same pitchlength, L, and radius, R, but the helices have a phase difference of Phi 
!with their neighbour and are corkscrewed by a distance sReg along the helical space curve with respect to their neighbours. 
!All the helices share the same principal axis. (taken to be the z-axis).
!Each ellipsoid is placed at position parameterised by the standard helical argument, t, on its own helix.
!Each cluster has nL layers of ellipsoids, each layer containing nS ellipsoids (one for each helix).
!The t is defined by which layer an ellipsoid is in within the cluster.
!The layers are spaced by dT where dT is computed to give a sensible spacing based on the 
!parameter spacer*aa, where aa is the size of the ellipsoid in the 'a' direction and spacer is just a multiplier.
!Thus all the ellipsoids are equally spaced along the helix, and tMax is as large as it needs to be, with however many twists are required.

!The cluster has a zero point, which is at the point where t=0, sReg=0 and phi=0 and the positions 
!of all ellipsoids within the cluster are given relative to this point. the helices coincide at this point.

!In the final analysis each ellipsoid is oriented in the Frenet frame of each helix so that 
!c axis of ellipsoid is parallel with the T vector of Frenet Frame of its helix
!b axis of ellipsoid is parallel with the N vector of Frenet Frame of its helix
!a axis is ellipsoid is parallel with the B vector of Frenet Frame of its helix

!The positions and orientations of a sequence of clusters is then computed such that the clusters are spread uniformly over a vertical line with no rotation

!The params in the order that they appear in the .model file, and therefore the PARAMS array, are:
!space is the multiple of aa which the ellipses along one helix are spaced by.
!aa ab ac ra rb rc - the sizes of the ellipses are all the same in this model
!L is the pitch length of the helix in the same units as the sizes of the ellipsoids
!R is the radius of the helix in the same units as the size of the ellipsoids
!sReg is the registration shift between two adjacent filaments in the model measured as an arc length displacement in the same units as the sizes of the ellipsoids
!Phi is the phase rotation between the two helices.
!Chirality was an experimental parameter to change the chirality of the cluster ordering helix, but not the chirality of the helices within each cluster.
!nS is the number of helices
!nL is thenumber of layers within each cluster
      implicit none


      INTEGER :: nP, nC, nE, nS, nL
      INTEGER :: iErr
      DOUBLE PRECISION :: lPARAMS(1:nP) 
      DOUBLE PRECISION :: lESIZE(1:nE,1:6), lEPOS(1:nE,1:3) 
      DOUBLE PRECISION :: lCPOS(1:nC,1:3) 
      DOUBLE PRECISION :: lCRMAT(1:nC,1:3,1:3),lERMAT(1:nE,1:3,1:3)
      DOUBLE PRECISION :: T(1:3),N(1:3),B(1:3)
      DOUBLE PRECISION :: Chirality, spacer, tMin, tMax, dT, tP, tC, C, aa, ab, ac, ra, rb, rc, R, L, sReg, Phi
      DOUBLE PRECISION, parameter :: pi=3.1415926
      INTEGER :: curE, curL, curS, curC, outI
  
! Initialise output error, zero good.
      iErr=0

!Copy params across
      spacer=lParams(1)
      aa=lPARAMS(2)
      ab=lPARAMS(3)
      ac=lPARAMS(4)
      ra=lPARAMS(5)
      rb=lPARAMS(6)
      rc=lPARAMS(7)
      L=lPARAMS(8)
      R=lPARAMS(9)
      sReg=lPARAMS(10)
      Phi=lPARAMS(11)
      Chirality=lParams(12)
      nS=lParams(13)
      nL=lParams(14)

!convert the pitch length parameter into length gained per radian.
      C=L/(2*DBLE(pi))

!convert the angular separation between ellipses into radians. 
      Phi=Phi*DBLE(pi)/180.0

!make sure there are a sensible number of ellipses in the cluster if not, then exit with an error code
      IF (nS*nL.eq.nE) THEN
       write(6,*)"Consistent number of ellipsoids in cluster."
      ELSE
       write(6,*)"Ooops. Inconsistent number of ellipsoids in cluster."
       iErr=1
       RETURN
      ENDIF

!compute the range of t. 
      tMin=0.0

!range of t finishes such that the arc length of the helix equals the total number of 
!ellipsoids in one filament times the size of each ellipsoid plus a 
!margin defined by the variable spacer. avoids overlaps and excessive forces.
      tMax=DBLE(nC*nL-1)*spacer*aa/SQRT(R*R+C*C)+tMin

!compute twist angle per layer of ellipsoid: totalNumLayers=nL*nC
      IF ((nC*nL).eq.1) THEN
        dT=0
      ELSE
        dT=((tMax-tMin)/DBLE(nC*nL-1))  
      END IF

!Design the cluster. The ellipsoids position and orientation is given relative to the lab frame origin and basis.
      DO curL=1,nL,1

        !compute the helical argument of the ellipsoids in the current layer
        tP = dT*DBLE(curL-1)

        !loop through the nS sheets and compute the relative offsets of each helix sReg and Phi.
        DO curS=1,nS,1

           !Compute out index
           outI=nS*(curL-1)+curS

           !compute frenet frame for helix with axis aligned to z axis in lab frame and x axis aligned to the helix at the point t=0, sreg=0, phi=0.
           CALL FRENETHELIX(R,Chirality*C,tP,DBLE(curS-1)*sReg,DBLE(curS-1)*phi, T(1:3), N(1:3), B(1:3),lEPOS(outI,1:3))

           !The first row vector in the rotation matrix is the orientation of the a-axis, which is aligned with B vector of TNB frame.
           lERMAT(outI,1,:)=B(1:3)
           !second row vector is the orientation of b-axis aligned with the Normal vector of TNB frame.
           lERMAT(outI,2,:)=N(1:3)
           !third vector is the orientation of the c-axis which is aligned with the T vector of TNB frame.
           lERMAT(outI,3,:)=T(1:3)


          !output size information
           lESIZE(outI, 1) =aa
           lESIZE(outI, 2) =ab
           lESIZE(outI, 3) =ac
           lESIZE(outI, 4) =ra
           lESIZE(outI, 5) =rb
           lESIZE(outI, 6) =rc
        END DO
      END DO


      !Loop through each cluster
      DO curC=1,nC,1

         !compute the argument of the cluster position
         IF (nC.eq.1) THEN 
           tC=0
         ELSE
           tC=DBLE(curC-1)*(tMax-DBLE(nL-1)*dT-tMin)/DBLE(nC-1) + tMin
         END IF

         !compute the direction cosines and position for the cluster, output to the right matrix
         CALL FLATHELIX(R, Chirality*C, tC, DBLE(0.0), DBLE(0.0), lCRMAT(curC,1:3,1:3),lCPOS(curC, 1:3))

         !Using FLATHELIX sets the z-scaling in terms of L and R so it's the same as helixModel, but the spacing can be adjusted with spacer.
         lCPOS(curC,1)=0
         lCPOS(curC,2)=0
         !remove the orientations of the cluster
         lCRMAT(curC,1,1)=1
         lCRMAT(curC,1,2)=0
         lCRMAT(curC,1,3)=0
         lCRMAT(curC,2,1)=0
         lCRMAT(curC,2,2)=1
         lCRMAT(curC,2,3)=0
         lCRMAT(curC,3,1)=0
         lCRMAT(curC,3,2)=0
         lCRMAT(curC,3,3)=1


    !     write(6,*) lCPOS(curC,1), lCPOS(curC,2), lCPOS(curC,3)
   !        write(6,*) lCRMAT(curC,1,1), lCRMAT(curC,2,1), lCRMAT(curC,3,1)
   !        write(6,*) lCRMAT(curC,1,2), lCRMAT(curC,2,2), lCRMAT(curC,3,2)
   !        write(6,*) lCRMAT(curC,1,3), lCRMAT(curC,2,3), lCRMAT(curC,3,3)
      END DO

      return
      END SUBROUTINE VERTICALHELIX
