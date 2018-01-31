!
! This program calculates the Voronoi polygons for the Thomson problem
!
PROGRAM TWODVORONOI
IMPLICIT NONE
INTEGER NPOINTS, J1, J2, J3, J4, NVPT, NEXTONE, NCURR, NPREV, N4, N5, N6, N7, N8, N9
INTEGER, PARAMETER :: NNEIGH=9
DOUBLE PRECISION DISTANCE, NDIST(NNEIGH), DIST1, DIST2, DUMMY
DOUBLE PRECISION,ALLOCATABLE :: X(:), Y(:), Z(:), NNX(:), NNY(:), NNZ(:)
DOUBLE PRECISION x0,y0,z0,xi,yi,zi,xj,yj,zj,X1,Y1,Z1,X2,Y2,Z2,VPX,VPY,VPZ,VX(NNEIGH),VY(NNEIGH),VZ(NNEIGH)
INTEGER NID(NNEIGH), NDONE
CHARACTER(LEN=2) ZSYM
CHARACTER(LEN=80) FNAME
LOGICAL IGNORE, USED(NNEIGH), DEBUG

DEBUG=.TRUE.
DEBUG=.FALSE.
NDONE=0
10 CONTINUE
NDONE=NDONE+1
READ(*,*,END=20) NPOINTS
ALLOCATE(X(NPOINTS),Y(NPOINTS),Z(NPOINTS),NNX(NNEIGH),NNY(NNEIGH),NNZ(NNEIGH))
READ(*,*) 
READ(*,*) (ZSYM,X(J1),Y(J1),Z(J1),J1=1,NPOINTS)
IF (NDONE.GT.99) THEN
   WRITE(FNAME,'(A10,I3)'), 'mathinput.',NDONE
ELSE IF (NDONE.GT.9) THEN
   WRITE(FNAME,'(A10,I2)'), 'mathinput.',NDONE
ELSE 
   WRITE(FNAME,'(A10,I1)'), 'mathinput.',NDONE
ENDIF

OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='UNKNOWN')
PRINT '(A,I8,A)','NDONE,FNAME=',NDONE,FNAME
WRITE(1,'(A)') 'Graphics3D[{EdgeForm[{GrayLevel[0.5], Thickness[0.005]}],'

! Find the NNEIGH nearest neighbours of each point and construct the
! corresponding Voronoi polygon.

N4=0; N5=0; N6=0; N7=0; N8=0; N9=0
DO J1=1,NPOINTS
   NDIST(1:NNEIGH)=1.0D100
   DO J2=1,NPOINTS
      IF (J2.EQ.J1) CYCLE
      DISTANCE=(X(J1)-X(J2))**2+(Y(J1)-Y(J2))**2+(Z(J1)-Z(J2))**2
      neighloop: DO J3=1,NNEIGH
         IF (DISTANCE.LE.NDIST(J3)) THEN
            DO J4=NNEIGH,J3+1,-1 ! move the others down to make room
               NDIST(J4)=NDIST(J4-1)
               NID(J4)=NID(J4-1)
               NNX(J4)=NNX(J4-1)
               NNY(J4)=NNY(J4-1)
               NNZ(J4)=NNZ(J4-1)
            ENDDO
            NDIST(J3)=DISTANCE
            NID(J3)=J2
            NNX(J3)=X(J2)
            NNY(J3)=Y(J2)
            NNZ(J3)=Z(J2)
            EXIT neighloop
         ENDIF
      ENDDO neighloop
   ENDDO
   IF (DEBUG) THEN
      PRINT '(A,I5,A,9I5)','neighbours of ion ',J1,' are ',NID(1:NNEIGH)
      PRINT '(A,9G20.10)','distances: ',NDIST(1:NNEIGH)
      PRINT '(A)','points:'
      PRINT '(3G20.10)',(NNX(J2),NNY(J2),NNZ(J2),J2=1,NNEIGH)
   ENDIF

   NVPT=0
   x0=X(J1);y0=Y(J1);z0=Z(J1)
   DO J2=1,NNEIGH 
      neighloop2: DO J3=J2+1,NNEIGH ! find the points of intersection of the great circles associated with each 
                        ! pair of neighbours

         xj=NNX(J3);yj=NNY(J3);zj=NNZ(J3)
         xi=NNX(J2);yi=NNY(J2);zi=NNZ(J2)
         IF (DEBUG) PRINT '(A,2I8)',' trying neighbours: ',NID(J2),NID(J3)
         IF (DEBUG) PRINT '(A,3G20.10)',' x0: ',x0,y0,z0
         IF (DEBUG) PRINT '(A,3G20.10)',' xi: ',xi,yi,zi
         IF (DEBUG) PRINT '(A,3G20.10)',' xj: ',xj,yj,zj
         DUMMY=Sqrt(x0**2*yi**2-2*x0**2*yi*yj+x0**2*yj**2+yi**2*z0**2-2*yi*yj*z0**2+yj**2*z0**2+xj**2* &
     &        ((y0-yi)**2+(z0-zi)**2)-2*y0*yi*z0*zi+ &
     &        2*y0*yj*z0*zi+2*yi*yj*z0*zi-2*yj**2*z0*zi+x0**2*zi**2+y0**2*zi**2-2*y0*yj*zi**2+yj**2*zi**2- &
     &        2*xi*(xj*((y0-yi)*(y0-yj)+(z0-zi)*(z0-zj))+x0*((y0-yj)*(yi-yj)+(z0-zj)*(zi-zj)))+xi**2*((y0-yj)**2 +(z0 - zj)**2)+&
     &    2*x0*xj*((y0 - yi)*(yi - yj) + (z0 - zi)*(zi - zj)) - 2*(x0**2*zi + (y0 - yi)*(-(yi*z0) + yj*z0 + y0*zi - yj*zi))*zj +&
     &    (x0**2 + (y0 - yi)**2)*zj**2)

         X1 = (yj*(-z0 + zi) + yi*(z0 - zj) + y0*(-zi + zj))/DUMMY
         Y1 = (xj*(z0 - zi) + x0*(zi - zj) + xi*(-z0 + zj))/DUMMY
         Z1 = (xj*(-y0 + yi) + xi*(y0 - yj) + x0*(-yi + yj))/DUMMY
         X2 = -X1
         Y2 = -Y1
         Z2 = -Z1
         IF (DUMMY.EQ.0.0D0) THEN
            IF (DEBUG) PRINT '(A,G20.10)','denom=',DUMMY
            IF (DEBUG) PRINT '(6G20.10)',X1,Y1,Z1,X2,Y2,Z2
            IF (DEBUG) PRINT '(A)','avoiding division by zero'
            CYCLE neighloop2
         ENDIF

        DIST1=(X(J1)-X1)**2+(Y(J1)-Y1)**2+(Z(J1)-Z1)**2
        DIST2=(X(J1)-X2)**2+(Y(J1)-Y2)**2+(Z(J1)-Z2)**2
        IF (DIST1.LT.DIST2) THEN
           VPX=X1; VPY=Y1; VPZ=Z1
        ELSE
           VPX=X2; VPY=Y2; VPZ=Z2
           DIST1=DIST2
        ENDIF
        IF (DEBUG) PRINT '(A,3I8,3G20.10)','proposed point for J1,NID(J2),NID(J3): ',J1,NID(J2),NID(J3),VPX,VPY,VPZ
        IF (VPX**2+VPY**2+VPZ**2.LT.0.1D0) THEN
           IF (DEBUG) PRINT '(A)','ignoring proposed point at the origin'
           CYCLE neighloop2
        ENDIF
        neighloop3: DO J4=1,NNEIGH ! ignore the proposed Voronoi point if it is nearer to one of the nearest neighbours
                                   ! other than the two for which we have just found the intersection of the great circles.
                                   ! These should be equidistant, but the numerical errors can confuse the program.
           IF ((J4.EQ.J2).OR.(J4.EQ.J3)) CYCLE neighloop3
           DIST2=(NNX(J4)-VPX)**2+(NNY(J4)-VPY)**2+(NNZ(J4)-VPZ)**2
           IF (DEBUG) PRINT '(A,I8,2G20.10)','NID(J4),DIST2,DIST1=',NID(J4),DIST2,DIST1
           IF (DIST2.LT.DIST1-1.0D-6) THEN  
              IF (DEBUG) PRINT '(A,I8,2G20.10)','rejecting point: NID(J4),DIST2,DIST1=',NID(J4),DIST2,DIST1
              CYCLE neighloop2
           ENDIF
        ENDDO neighloop3
        DO J4=1,NVPT ! ignore the proposed Voronoi point if it is not new
           DIST2=(VX(J4)-VPX)**2+(VY(J4)-VPY)**2+(VZ(J4)-VPZ)**2
           IF (DIST2.LT.1.0D-6) THEN
              CYCLE neighloop2
           ENDIF
        ENDDO 
        NVPT=NVPT+1
        VX(NVPT)=VPX; VY(NVPT)=VPY; VZ(NVPT)=VPZ
     ENDDO neighloop2
   ENDDO 
   IF (DEBUG) PRINT '(I6,A,I6)',NVPT,' Voronoi points for ion ',J1
   IF (NVPT.EQ.4) THEN
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[1,1,0],RGBColor[0.5,0.5,0.5],100],'
      N4=N4+1
   ELSE IF (NVPT.EQ.5) THEN
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[1,0,0],RGBColor[0.5,0.5,0.5],100],'
      N5=N5+1
   ELSE IF (NVPT.EQ.6) THEN
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[0,1,0],RGBColor[0.5,0.5,0.5],100],'
      N6=N6+1
   ELSE IF (NVPT.EQ.7) THEN
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[0,0,1],RGBColor[0.5,0.5,0.5],100],'
      N7=N7+1
   ELSE IF (NVPT.EQ.8) THEN
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[1,0,1],RGBColor[0.5,0.5,0.5],100],'
      N8=N8+1
   ELSE IF (NVPT.EQ.9) THEN
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[1,0,1],RGBColor[0.5,0.5,0.5],100],'
      N9=N9+1
   ELSE 
      WRITE(1,'(A)') 'SurfaceColor[RGBColor[0,1,1],RGBColor[0.5,0.5,0.5],100],'
   ENDIF
   WRITE(1,'(A)') 'Polygon[{'
   WRITE(1,'(A1,F20.10,A1,F20.10,A1,F20.10,A1)') '{',VX(1),',',VY(1),',',VZ(1),'}'
   NCURR=1
   USED(2:NVPT)=.FALSE.
   USED(1)=.TRUE.
   DO J2=1,NVPT-1 ! find the nearest neighbour different from the previous Voronoi point
      WRITE(1,'(A1)') ','
      DIST1=1.0D100
      DO J3=1,NVPT
         IF (USED(J3)) CYCLE
         DIST2=(VX(NCURR)-VX(J3))**2+(VY(NCURR)-VY(J3))**2+(VZ(NCURR)-VZ(J3))**2
         IF (DIST2.LT.DIST1) THEN
            NEXTONE=J3
            DIST1=DIST2
         ENDIF
      ENDDO
      NCURR=NEXTONE
      USED(NEXTONE)=.TRUE.
      WRITE(1,'(A1,F20.10,A1,F20.10,A1,F20.10,A1)') '{',VX(NEXTONE),',',VY(NEXTONE),',',VZ(NEXTONE),'}'
   ENDDO
   IF (J1.EQ.NPOINTS) THEN
      WRITE(1,'(A)') '}]'
   ELSE
      WRITE(1,'(A)') '}],'
   ENDIF
ENDDO
WRITE(1,'(A)') '}]'
CLOSE(1)
IF (N4+N5+N6+N7+N8+N9.NE.NPOINTS) PRINT '(2(A,I6))','WARNING - total is ',N4+N5+N6+N7+N8+N9,' should be ',NPOINTS
PRINT '(A,6I6)','Number of 4, 5, 6, 7, 8, 9 vertices is ',N4,N5,N6,N7,N8,N9
GOTO 10
20 CONTINUE

END PROGRAM TWODVORONOI
