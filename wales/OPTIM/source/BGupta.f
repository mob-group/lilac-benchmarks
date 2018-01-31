C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C   Loop structure recoded by J.A. Elliott 2009
C
C   GMIN is free software; you can redistribute it and/or modIFy
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; IF not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C*************************************************************************
C
C   BINARY GUPTA POTENTIAL
C   CJH 03/2011


c    TO RUN WITH OPTIM, SET BGUPTAT, BGUPTATAB AND BGUPTATBB AS LINES IN THE 
c    ODATA FILE (FOR KEYWORDS.F TO PICK UP). FOR EACH OF THESE LINES GIVE A,P,Q,Z 
C    AND R0, BUT FOR THE BGUPTAT LINE ONLY, PRECEED THESE PARAMETERS WITH NTYPEA
c    (THE NUMBER OF A TYPE ATOMS)
c    PUT COORDINATES AT THE BOTTOM OF THE ODATA FILE, LABELLED BY G1 AND G2.
c
c    SO HAVING BGUPTAT,BGUPTATAB AND BGUPTBB TELLS OPTIM TO GET THE POTENTIAL
c    PARAMETERS.  HAVING G1 AND G2 TELLS OPTIM TO CALCULATE WITH THE BINARY GUPTA
c    POTENTIAL, WITH NATOMS-NTYPEA "B" ATOMS

C
C*************************************************************************
C

        MODULE BGUPMOD
 
       DOUBLE PRECISION AAA,AAB,ABB,PAA,PAB,PBB,QAA,QAB,QBB,ZAA,ZAB,ZBB,
     1 R0AA,R0AB,R0BB
       
        END MODULE BGUPMOD

      SUBROUTINE BGUPTA(X,V,PG,GRADT,STEST)
      USE commons
      USE BGUPMOD
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, ATOMTYPE, NTYPEA, J3
      DOUBLE PRECISION X(3*NATOMS), PG, DIST, PTEMP, V(3*NATOMS),
     1       GRHO(NATOMS), VTEMP, DUMMY, DISTANCE_MATRIX(NATOMS,NATOMS),
     2       VTEMP1(NATOMS), VTEMP2(NATOMS), VTEMP3(NATOMS),
     3       RTEMP, DX, DY, DZ, WW2, WW3,
     4       GA, P,Q,R0,EPSAB,EPSBB,SIGAB,SIGBB
C             R0AA, AAB, ABB, 
c     5       P, PAA, PAB, PBB,
c     6       Q, QAA, QAB, QBB, 
c     7       GXI, ZAA, ZAB, ZBB, 
c     8       R0, R0AA, R0AB, R0BB 
      LOGICAL GRADT,STEST
 
C      COMMON /BIN/ NTYPEA, AAA, AAB, ABB, PAA, PAB, PBB, QAA, QAB,
C     &   QBB,  ZAA,  ZAB, ZBB, R0AA, R0AB, R0BB
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB,NTYPEA       

! MAKE ARRAYS FOR THE VALUES OF THE 5 GUPTA PARAMETERS. EASY TO SCALE
! TO HIGHER MULTIMETALLICS (AND HOPEFULLY QUICKER THAN USING IF LOOPS

       DOUBLE PRECISION AARRAY(NATOMS,NATOMS) 
       DOUBLE PRECISION ZARRAY(NATOMS,NATOMS)
       DOUBLE PRECISION PARRAY(NATOMS,NATOMS)
       DOUBLE PRECISION QARRAY(NATOMS,NATOMS)
       REAL R0ARRAY(NATOMS,NATOMS)

       DOUBLE PRECISION XX,YY,ZZ,XY,YZ,XZ,WW0,WW1,W0,WX,WY,WZ,W1,W2,W3,WXM,WYM, 
     &                  WZM,DIST13,DIST23,W12,W13,W23,WWX1,WWY1,WWZ1,WWX2,WWY2,WWZ2,WWX3,WWY3,WWZ3,
     &                  WX3,WY3,WZ3,WX3M,WY3M,WZ3M
 
c        WRITE(*,*) "******************************"
c        WRITE(*,*) "******************************"
c        WRITE(*,*) "NATOMS IS " ,NATOMS
c        WRITE(*,*) "******************************"
c        WRITE(*,*) "******************************"
!       R0ARRAY = 1.0d0

! FOR TESTING OF PARAMETERS (COMMONS BLOCKS ARE EVIL)
!       WRITE(*,*) 'AT START OF BGUPTA.F ROUTINE'
!       WRITE(*,*) 'NTYPEA IS', NTYPEA
!       WRITE(*,*) 'AAA AAB ABB ARE    ',AAA,AAB,ABB
!       WRITE(*,*) 'PAA PAB PBB ARE    ',PAA,PAB,PBB
!       WRITE(*,*) 'QAA QAB QBB ARE    ',QAA,QAB,QBB
!       WRITE(*,*) 'ZAA ZAB ZBB ARE    ',ZAA,ZAB,ZBB
!       WRITE(*,*) 'ROAA ROAB ROBB ARE ',R0AA,R0AB,R0BB
       DO J1=1,NTYPEA
          DO J2=1,NTYPEA
              AARRAY(J1,J2)=2.0d0*AAA
              PARRAY(J1,J2)=PAA
              QARRAY(J1,J2)=QAA
              ZARRAY(J1,J2)=ZAA
              R0ARRAY(J1,J2)=R0AA
c              WRITE(*,*) 'aaaarray is ' ,AARRAY(J1,J2)
c              WRITE(*,*) 'paaarray is ' ,PARRAY(J1,J2)
c              WRITE(*,*) 'qaaarray is ' ,QARRAY(J1,J2)
c              WRITE(*,*) 'zaaarray is ' ,ZARRAY(J1,J2)
c              WRITE(*,*) 'raaarray is ' ,R0ARRAY(J1,J2)
          ENDDO
       ENDDO
c       WRITE(*,*) "***********FINISHED AAARRAY LOOP**************"
       DO J1=1,NTYPEA
          DO J2=NTYPEA+1,NATOMS
              AARRAY(J1,J2)=2.0d0*AAB
              PARRAY(J1,J2)=PAB
              QARRAY(J1,J2)=QAB
              ZARRAY(J1,J2)=ZAB
              R0ARRAY(J1,J2)=R0AB
c              WRITE(*,*) 'aabarray is ' ,AARRAY(J1,J2)
c              WRITE(*,*) 'pabarray is ' ,PARRAY(J1,J2)
c              WRITE(*,*) 'qabarray is ' ,QARRAY(J1,J2)
c              WRITE(*,*) 'zabarray is ' ,ZARRAY(J1,J2)
c              WRITE(*,*) 'rabarray is ' ,R0ARRAY(J1,J2)
          END DO
       ENDDO
       DO J1=NTYPEA+1,NATOMS
          DO J2=1,NTYPEA
              AARRAY(J1,J2)=2.0d0*AAB
              PARRAY(J1,J2)=PAB
              QARRAY(J1,J2)=QAB
              ZARRAY(J1,J2)=ZAB
              R0ARRAY(J1,J2)=R0AB
c              WRITE(*,*) 'aabarray is ' ,AARRAY(J1,J2)
c              WRITE(*,*) 'pabarray is ' ,PARRAY(J1,J2)
c              WRITE(*,*) 'qabarray is ' ,QARRAY(J1,J2)
c              WRITE(*,*) 'zabarray is ' ,ZARRAY(J1,J2)
c              WRITE(*,*) 'rabarray is ' ,R0ARRAY(J1,J2)
          END DO
       ENDDO
c       WRITE(*,*) "***********FINISHED AABRRAY LOOP**************"
       DO J1=NTYPEA+1,NATOMS
          DO J2=NTYPEA+1,NATOMS
              AARRAY(J1,J2)=2.0d0*ABB
              PARRAY(J1,J2)=PBB
              QARRAY(J1,J2)=QBB
              ZARRAY(J1,J2)=ZBB
              R0ARRAY(J1,J2)=R0BB
c              WRITE(*,*) 'abbarray is ' ,AARRAY(J1,J2)
c              WRITE(*,*) 'pbbarray is ' ,PARRAY(J1,J2)
c              WRITE(*,*) 'qbbarray is ' ,QARRAY(J1,J2)
c              WRITE(*,*) 'zbbarray is ' ,ZARRAY(J1,J2)
c              WRITE(*,*) 'rbbarray is ' ,R0ARRAY(J1,J2)
          END DO
       ENDDO
c       WRITE(*,*) "***********FINISHED ABBRRAY LOOP**************"


!C Put R0 in reduced units, and adjust some constant
       R0=1.0D0 !makes bond lengths correct
       GA=1.0D0     !double A inside loops and not here for better energies.  ! Note that in orignal Cleri & Rosato paper, they define
                                                   ! 2-boDY term as a sum over all atom pairs, but it is much quicker to
                                                   ! evaluate sum over just J2>J1 and DOuble the constant A
!      twoq=2.0D0*q

      DO J2=1,NATOMS                               ! initialise distance and density arrays
        DO J1=1,NATOMS
          DISTANCE_MATRIX(J1,J2)=0.0D0
        ENDDO
        GRHO(J2)=0.0D0
      ENDDO

      RTEMP = 0.0D0                               ! initialise energy accumulators
      PTEMP = 0.0D0

      DO J1=1,NATOMS-1                              ! begin outer loop over all atoms except last
         DO J2=J1+1,NATOMS                           ! begin first inner loop over all J2>J1
               DISTANCE_MATRIX(J2,J1)=                           ! calc and store distance between atoms J1,J2
     1         DSQRT(( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2         ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3         ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2)

               DIST=DISTANCE_MATRIX(J2,J1)                       ! store the actual distance used in inner loop
               DISTANCE_MATRIX(J1,J2)=DIST                       ! impose symmetry on distance matrix
 
               DIST=DIST/R0ARRAY(J1,J2)
               PTEMP=PTEMP+(AARRAY(J1,J2)*DEXP(PARRAY(J1,J2) * (1-DIST)))                ! accumulate two-body potential term
               RTEMP=DEXP(2.0D0 * QARRAY(J1,J2) * (1-DIST))                   ! calculate many-body potential term
               GRHO(J1)=GRHO(J1)+((ZARRAY(J1,J2)**2)*RTEMP)
               GRHO(J2)=GRHO(J2)+((ZARRAY(J1,J2)**2)*RTEMP)
         ENDDO                                       ! END inner loop over all J2>J1
      ENDDO                                         ! END outer loop over all atoms except last

!c     Now, sum the potential energy over all atoms

      PG=PTEMP                                   ! initialise potential energy
c      WRITE(*,*) 'PG = ' ,PG

      DO J1=1,NATOMS
        GRHO(J1)=DSQRT(GRHO(J1))                    ! square root density
        PG=PG-(GRHO(J1) )                         ! accumulate potential energy
      ENDDO

c      WRITE(*,*) 'PG -GRHO = ' ,PG
!c     Calculate gradient terms, if required

      IF (GRADT) THEN
         DO J1=1,NATOMS                            ! initialise total gradient terms
             V(3*(J1-1)+1)=0.d0
             V(3*(J1-1)+2)=0.d0
             V(3*(J1-1)+3)=0.d0
             VTEMP1(J1)=0.0D0
             VTEMP2(J1)=0.0D0
             VTEMP3(J1)=0.0D0
         ENDDO

      DO J1=1,NATOMS-1                                ! begin outer loop over all atoms except last
          DUMMY=1.0D0/GRHO(J1)                          ! store reciprocal of density element for atom J1
          DO J2=J1+1,NATOMS                             ! begin inner loop over all J2>J1
              DIST=DISTANCE_MATRIX(J1,J2)                                        ! recall distance from earlier loop
              DIST=DIST/R0ARRAY(J1,J2)
              VTEMP=(QARRAY(J1,J2)*(ZARRAY(J1,J2)**2)*(DUMMY+1.0D0/GRHO(J2))*DEXP(2.0D0 * QARRAY(J1,J2) *(1-DIST))      ! calculate gradient term
     1        -AARRAY(J1,J2)*PARRAY(J1,J2)*DEXP(PARRAY(J1,J2)*(1-DIST)))/(DISTANCE_MATRIX(J1,J2)*R0ARRAY(J1,J2))

              DX=(X(3*(J1-1)+1)-X(3*(J2-1)+1))         ! calculate Cartesian components of distance
              DY=(X(3*(J1-1)+2)-X(3*(J2-1)+2))
              DZ=(X(3*(J1-1)+3)-X(3*(J2-1)+3))

              VTEMP1(J1)=VTEMP1(J1)+VTEMP*DX           ! accumulate primary gradient components
              VTEMP2(J1)=VTEMP2(J1)+VTEMP*DY
              VTEMP3(J1)=VTEMP3(J1)+VTEMP*DZ

              VTEMP1(J2)=VTEMP1(J2)-VTEMP*DX           ! accumulate symmetric gradient components
              VTEMP2(J2)=VTEMP2(J2)-VTEMP*DY
              VTEMP3(J2)=VTEMP3(J2)-VTEMP*DZ
          ENDDO
      ENDDO

!c       Finally, sum the gradient terms over all atoms

      DO J1=1,NATOMS
            V(3*(J1-1)+1)=V(3*(J1-1)+1)+VTEMP1(J1)
            V(3*(J1-1)+2)=V(3*(J1-1)+2)+VTEMP2(J1)
            V(3*(J1-1)+3)=V(3*(J1-1)+3)+VTEMP3(J1)
c      write(*,*)  'VX AT END OF BGUPTA IS  ' ,V(3*(J1-1)+1)
c      write(*,*)  'VY AT END OF BGUPTA IS  ' ,V(3*(J1-1)+2)
c      write(*,*)  'VZ AT END OF BGUPTA IS  ' ,V(3*(J1-1)+3)
      ENDDO
c      WRITE(*,*)"**************END OF VX,VY,VZ LOOP*************"

      ENDIF

!      WRITE(*,*) 'AT END OF BGUPTA.F ROUTINE'

      IF (STEST) THEN

         DO J1=1,3*NATOMS       ! initialise hessian matrix
            DO J2=1,3*NATOMS
               HESS(J1,J2)=0.D0
            ENDDO
         ENDDO
         
         DO J1=1,NATOMS
            XX=0.D0
            YY=0.D0
            ZZ=0.D0
            
            XY=0.D0
            YZ=0.D0
            XZ=0.D0
            
            DO J2=1,NATOMS
               
               IF (J1.NE.J2) THEN
                  DX=(X(3*(J1-1)+1)-X(3*(J2-1)+1))
                  DY=(X(3*(J1-1)+2)-X(3*(J2-1)+2))
                  DZ=(X(3*(J1-1)+3)-X(3*(J2-1)+3))
                  
                  DIST=DISTANCE_MATRIX(J1,J2)
                  DIST=DIST/R0ARRAY(J1,J2)
                  
                  WW0=QARRAY(J1,J2)*(1.D0/GRHO(J1)+1.D0/GRHO(J2))
     $                 *(ZARRAY(J1,J2)**2)*DEXP(2.0D0 * QARRAY(J1,J2)
     $                 *(1-DIST))
                  WW1=PARRAY(J1,J2)*AARRAY(J1,J2)*DEXP(PARRAY(J1,J2)*(1-DIST))
                  W0=(WW0-WW1)/(DIST*R0ARRAY(J1,J2)**2)
                  
                  WX=W0
                  WY=W0
                  WZ=W0
                  
                  W1=-WW0*(1.D0/DISTANCE_MATRIX(J1,J2)+2.D0*QARRAY(J1,J2
     $                 )/R0ARRAY(J1,J2))
                  W1=W1/(R0ARRAY(J1,J2)*DISTANCE_MATRIX(J1,J2)**2)
                  
                  W2=WW1*(1.D0/DISTANCE_MATRIX(J1,J2)+PARRAY(J1,J2)
     $                 /R0ARRAY(J1,J2))
                  W2=W2/(R0ARRAY(J1,J2)*DISTANCE_MATRIX(J1,J2)**2)
                  
                  W3=QARRAY(J1,J2)*(ZARRAY(J1,J2)**2)*DEXP(2.0D0
     $                 * QARRAY(J1,J2) *(1-DIST))
                  W3=W3/(R0ARRAY(J1,J2)*DISTANCE_MATRIX(J1,J2))
                  W3=W3*W3*(1.D0/GRHO(J1)**3+1.D0/GRHO(J2)**3)
                  
                  WX=WX+(W1+W2+W3)*DX*DX
                  WY=WY+(W1+W2+W3)*DY*DY
                  WZ=WZ+(W1+W2+W3)*DZ*DZ
                  
                  WZM=(W1+W2+W3)*(X(3*(J1-1)+1)-X(3*(J2-1)+1))*(X(3*(J1
     $                 -1)+2)-X(3*(J2-1)+2))
                  WYM=(W1+W2+W3)*(X(3*(J1-1)+1)-X(3*(J2-1)+1))*(X(3*(J1
     $                 -1)+3)-X(3*(J2-1)+3))
                  WXM=(W1+W2+W3)*(X(3*(J1-1)+3)-X(3*(J2-1)+3))*(X(3*(J1
     $                 -1)+2)-X(3*(J2-1)+2))

                  WX3=0.D0
                  WY3=0.D0
                  WZ3=0.D0
                  
                  WX3M=0.D0
                  WY3M=0.D0
                  WZ3M=0.D0
                  
                  DO J3=1,NATOMS
                     
                     IF (J1.NE.J3.AND.J2.NE.J3) THEN
                        DIST13=DISTANCE_MATRIX(J1,J3)
                        DIST13=DIST13/R0ARRAY(J1,J3)
                        DIST23=DISTANCE_MATRIX(J2,J3)
                        DIST23=DIST23/R0ARRAY(J2,J3)
                        
                        W13=(ZARRAY(J1,J3)**2)*DEXP(2.0D0 * QARRAY(J1,J3
     $                       ) *(1-DIST13))
                        W13=W13/DISTANCE_MATRIX(J1,J3)
                        W12=(ZARRAY(J1,J2)**2)*DEXP(2.0D0 * QARRAY(J1,J2
     $                       ) *(1-DIST))
                        W12=W12/DISTANCE_MATRIX(J1,J2)
                        W23=(ZARRAY(J2,J3)**2)*DEXP(2.0D0 * QARRAY(J2,J3
     $                       ) *(1-DIST23))
                        W23=W23/DISTANCE_MATRIX(J2,J3)

                        WW1=QARRAY(J1,J2)*QARRAY(J1,J3)*W12*W13*(1.D0
     $                       /GRHO(J1)**3)/R0ARRAY(J1,J2)/R0ARRAY(J1,J3)
                        WW2=QARRAY(J1,J2)*QARRAY(J2,J3)*W12*W23*(1.D0
     $                       /GRHO(J2)**3)/R0ARRAY(J1,J2)/R0ARRAY(J2,J3)
                        WW3=QARRAY(J1,J3)*QARRAY(J2,J3)*W13*W23*(1.D0
     $                       /GRHO(J3)**3)/R0ARRAY(J1,J3)/R0ARRAY(J2,J3)

                        WWX1=-DX*(X(3*(J1-1)+1)-X(3*(J3-1)+1))*WW1
                        WWX2= DX*(X(3*(J2-1)+1)-X(3*(J3-1)+1))*WW2
                        WWX3= (X(3*(J1-1)+1)-X(3*(J3-1)+1))*(X(3*(J2-1)
     $                       +1)-X(3*(J3-1)+1))*WW3
                        
                        WWY1=-DY*(X(3*(J1-1)+2)-X(3*(J3-1)+2))*WW1
                        WWY2= DY*(X(3*(J2-1)+2)-X(3*(J3-1)+2))*WW2
                        WWY3= (X(3*(J1-1)+2)-X(3*(J3-1)+2))*(X(3*(J2-1)
     $                       +2)-X(3*(J3-1)+2))*WW3
                        
                        WWZ1=-DZ*(X(3*(J1-1)+3)-X(3*(J3-1)+3))*WW1
                        WWZ2= DZ*(X(3*(J2-1)+3)-X(3*(J3-1)+3))*WW2
                        WWZ3= (X(3*(J1-1)+3)-X(3*(J3-1)+3))*(X(3*(J2-1)
     $                       +3)-X(3*(J3-1)+3))*WW3
                        
                        WX3=WX3+WWX1+WWX2+WWX3
                        WY3=WY3+WWY1+WWY2+WWY3
                        WZ3=WZ3+WWZ1+WWZ2+WWZ3
                        
                        WWZ1=(X(3*(J1-1)+1)-X(3*(J3-1)+1))*(X(3*(J2-1)+2
     $                       )-X(3*(J1-1)+2))*WW1
                        WWZ2=(X(3*(J1-1)+1)-X(3*(J2-1)+1))*(X(3*(J2-1)+2
     $                       )-X(3*(J3-1)+2))*WW2
                        WWZ3=(X(3*(J1-1)+1)-X(3*(J3-1)+1))*(X(3*(J2-1)+2
     $                       )-X(3*(J3-1)+2))*WW3
                        
                        WWY1=(X(3*(J1-1)+1)-X(3*(J3-1)+1))*(X(3*(J2-1)+3
     $                       )-X(3*(J1-1)+3))*WW1
                        WWY2=(X(3*(J1-1)+1)-X(3*(J2-1)+1))*(X(3*(J2-1)+3
     $                       )-X(3*(J3-1)+3))*WW2
                        WWY3=(X(3*(J1-1)+1)-X(3*(J3-1)+1))*(X(3*(J2-1)+3
     $                       )-X(3*(J3-1)+3))*WW3
                        
                        WWX1=(X(3*(J1-1)+2)-X(3*(J3-1)+2))*(X(3*(J2-1)+3
     $                       )-X(3*(J1-1)+3))*WW1
                        WWX2=(X(3*(J1-1)+2)-X(3*(J2-1)+2))*(X(3*(J2-1)+3
     $                       )-X(3*(J3-1)+3))*WW2
                        WWX3=(X(3*(J1-1)+2)-X(3*(J3-1)+2))*(X(3*(J2-1)+3
     $                       )-X(3*(J3-1)+3))*WW3
                        
                        WX3M=WX3M+WWX1+WWX2+WWX3
                        WY3M=WY3M+WWY1+WWY2+WWY3
                        WZ3M=WZ3M+WWZ1+WWZ2+WWZ3
                        
                     ENDIF
                     
                  ENDDO

                  WX3=-WX3
                  WY3=-WY3
                  WZ3=-WZ3
                  
                  WX=WX+WX3
                  WY=WY+WY3
                  WZ=WZ+WZ3
                  
                  XX=XX-WX
                  YY=YY-WY
                  ZZ=ZZ-WZ
                  
                  HESS(3*J1-2,3*J2-2)=-WX
                  HESS(3*J1-1,3*J2-1)=-WY
                  HESS(3*J1  ,3*J2  )=-WZ
                  
                  WX3M=-WX3M
                  WY3M=-WY3M
                  WZ3M=-WZ3M
                  
                  WXM=WXM+WX3M
                  WYM=WYM+WY3M
                  WZM=WZM+WZ3M
                  
                  XY=XY-WZM
                  XZ=XZ-WYM
                  YZ=YZ-WXM

                  HESS(3*J1-2,3*J2-1)=-WZM
                  HESS(3*J1-2,3*J2  )=-WYM
                  HESS(3*J1-1,3*J2  )=-WXM
                  
                  HESS(3*J2  ,3*J1-1)=HESS(3*J1-1,3*J2  )
                  HESS(3*J2  ,3*J1-2)=HESS(3*J1-2,3*J2  )
                  HESS(3*J2-1,3*J1-2)=HESS(3*J1-2,3*J2-1)

               ENDIF
               
            ENDDO

            HESS(3*J1-2,3*J1-2)=-XX
            HESS(3*J1-1,3*J1-1)=-YY
            HESS(3*J1  ,3*J1  )=-ZZ

            HESS(3*J1-2,3*J1-1)=-XY
            HESS(3*J1-1,3*J1-2)=-XY
            HESS(3*J1-2,3*J1  )=-XZ
            HESS(3*J1  ,3*J1-2)=-XZ
            HESS(3*J1-1,3*J1  )=-YZ
            HESS(3*J1  ,3*J1-1)=-YZ

         ENDDO
         
      ENDIF


      RETURN
      END


c------------------------------------------------------------

c   Attempt swaps between A and B atoms in the BGupta potential
c   TAKES "COULN" KEYWORD FROM COMMONS (USED IN LJCOULOMB ORIGINALLY)
     
!!      SUBROUTINE TAKESTEPBGUP(NP)
!
!      USE COMMONS ,  ONLY : NATOMS, COORDS, COULN
!      IMPLICIT NONE
!      INTEGER NP
!      DOUBLE PRECISION DPRAND
!      DOUBLE PRECISION X, Y, Z
!      INTEGER J1, J2
!
!
!      J1 = 3 * (INT(DPRAND()*COULN))
!      J2 = 3 * (INT(DPRAND()*(NATOMS-COULN)) + COULN)
!
!      X = COORDS(J1+1, NP)
!      Y = COORDS(J1+2, NP)
!      Z = COORDS(J1+3, NP)
!
!      COORDS(J1+1, NP) = COORDS(J2+1, NP) 
!      COORDS(J1+2, NP) = COORDS(J2+2, NP) 
!      COORDS(J1+3, NP) = COORDS(J2+3, NP) 
!
!      COORDS(J2+1, NP) = X
!      COORDS(J2+2, NP) = Y
!      COORDS(J2+3, NP) = Z
!
!      RETURN
!      END
