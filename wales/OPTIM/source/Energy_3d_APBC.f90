! 
! Three-dimensional anti periodic XY model. ch558
!
SUBROUTINE Energy_3d_APBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE KEY, ONLY : NONEDAPBC, XYPHI ! XYPHI has dim 2N*N, OPTIM reads NATOMS in as N**2 so we have to change it 
  IMPLICIT NONE
  INTEGER N, i1, i2, i3
  DOUBLE PRECISION, dimension((NONEDAPBC)**3) :: theta, GRAD
  DOUBLE PRECISION :: Energy, APBC1, APBC2, APBC3
  LOGICAL GTEST,SECT
    
  N=NONEDAPBC
!  theta(N*N)=0
  Energy=0.0D0

  i3 = 0
  Do WHILE( i3.LT.(N*N*N))
     i2=0
     DO WHILE ( i2.LT.(N*N))
        i1=0
        DO WHILE(i1.LT.(N))
           !1D part
           Energy = Energy + ( cos(xyphi(i1+i2+i3+1) + APBC1(i1+1,N)*theta( modulo(i1+1,N) + i2 +i3+ 1 ) - theta(i1+i2+i3+1)))

           !2D part
           Energy = Energy + (cos(xyphi(i1+i2+i3+(N*N)+1) + APBC2(i2+N,N)*theta(i1+ modulo(i2+N, N*N) +i3+1)- theta(i1+i2+i3+1)))

           !3D part
           Energy = Energy + (cos(xyphi(i1+i2+i3+(N*N*N)+1) + APBC3(i3+N*N,N)*theta(i1+i2+ modulo(i3+N*N, N*N*N) +1)- theta(i1+i2+i3+1)))
           
           i1=i1+1
        ENDDO
        i2 = i2 + N
     END DO
     i3 = i3+ N*N
  ENDDO
  
  Energy = 3.0D0 - (Energy/(N*N*N))
  !Energy = N*N*N*3.0D0 - Energy




  IF (.NOT.(GTEST.OR.SECT)) RETURN
  


  i3 = 0
  Do WHILE( i3.LT.(N*N*N))
     i2=0
     DO WHILE ( i2.LT.(N*N))
        i1=0
        DO WHILE(i1.LT.(N))           
           
           grad(i1+i2+i3+1)= &
                !1D part
                -sin(xyphi(i1+i2+i3+1) + APBC1(i1+1,N)*theta( modulo(i1+1,N) + i2+i3+ 1)-theta(i1+i2+i3+1)) &
                + APBC1(i1-1,N)*sin(xyphi(modulo(i1-1,N) + i2 +i3+ 1) + APBC1(i1-1,N)*theta(i1+i2+i3+1)-theta(modulo(i1-1,N) + i2 +i3+1)) & 
                !2D part
                - sin( xyphi(i1+i2+i3+N*N+1) + APBC2(i2+N,N)*theta(i1+ modulo(i2+N,N*N)+i3 +1)-theta(i1+i2+i3+1)) &
                + APBC2(i2-N,N)*sin(xyphi(i1+ modulo(i2-N, N*N) +i3+ N*N+1)+ APBC2(i2-N,N)*theta(i1+i2+i3+1)-theta(i1+ modulo(i2-N, N*N) +i3+1)) &
                !3D part
                - sin( xyphi(i1+i2+i3+N*N*N+1) + APBC3(i3+N*N,N)*theta(i1+ i2+modulo(i3+N*N,N*N*N) +1)-theta(i1+i2+i3+1)) &
                + APBC3(i3-N*N, N)*sin(xyphi(i1+ i2+modulo(i3-N*N, N*N*N) + N*N*N+1)+ APBC3(i3-N*N,N)*theta(i1+i2+i3+1)-theta(i1+ i2+modulo(i3-N*N, N*N*N) +1)) 
           
           
           i1=i1+1
        ENDDO
        i2=i2+N
     ENDDO
     i3 = i3+ N*N
  ENDDO
  
 ! grad(N*N)=0
  
  IF (.NOT.SECT) RETURN
  
  CALL Hessian_3d_APBC(THETA)

END SUBROUTINE Energy_3d_APBC

SUBROUTINE Hessian_3d_APBC(THETA)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  USE MODHESS
  IMPLICIT NONE
  INTEGER :: i1, i2,i3, N
  DOUBLE PRECISION, DIMENSION(NONEDAPBC**3) :: theta
  DOUBLE PRECISION :: APBC1, APBC2, APBC3

  N=NONEDAPBC
  ! THETA(N*N)=0
  HESS(:,:)=0.0D0


  i3 = 0
  Do WHILE( i3.LT.(N*N*N))
     i2=0
     DO WHILE ( i2.LT.(N*N))
        i1=0
        DO WHILE(i1.LT.(N))
           
           HESS(1+i1+i2+i3,1+i1+i2+i3)= &
                !1D part
                cos(xyphi(1+i1+i2+i3)+APBC1(i1+1,N)*theta( modulo(1+i1,N) +1+i2+i3)-theta(1+i1+i2+i3)) &
                +cos(xyphi(1+ modulo(i1-1,N) +i2+i3)+APBC1(i1-1,N)*theta(1+i1+i2+i3)-theta(1+ modulo(i1-1,N) +i2+i3)) &
                !2d part
                +cos(xyphi(1+i1+i2+i3+(N*N))+ APBC2(i2+N,N)*theta(1+i1+ modulo(i2+N,N*N) +i3)-theta(1+i1+i2+i3)) &
                +cos(xyphi(1+i1+modulo(i2-N,N*N) +i3+(N*N))+APBC2(i2-N,N)*theta(1+i1+i2+i3)-theta(1+i1+ modulo(i2-N,N*N)+i3  )) &
                !3D part
                +cos(xyphi(1+i1+i2+i3+(N*N*N))+ APBC3(i3+N*N,N)*theta(1+i1+i2+ modulo(i3+N*N,N*N*N) )-theta(1+i1+i2+i3)) &
                +cos(xyphi(1+i1+i2+modulo(i3-N*N,N*N*N) +(N*N*N))+APBC3(i3-N*N,N)*theta(1+i1+i2+i3)-theta(1+i1+ i2+ modulo(i3-N*N,N*N*N)  ))


           !1D part
           HESS(1+i1+i2+i3,modulo(1+i1,N) +1+i2+i3)=  - APBC1(i1+1,N)*cos(xyphi(1+i1+i2+i3)+ APBC1(i1+1,N)*theta(modulo(1+i1,N) +1+i2+i3)-theta(1+i1+i2+i3))
           
           !2D part
           HESS(1+i1+i2+i3,1+i1+ modulo(i2+N,N*N) +i3)= -APBC2(i2+N,N)*cos(xyphi(1+i1+i2+i3+(N*N))+APBC2(i2+N,N)*theta(1+i1+ modulo(i2+N,N*N) +i3)-theta(1+i1+i2+i3))

           !3D part
           HESS(1+i1+i2+i3,1+i1+ i2+ modulo(i3+N*N,N*N*N) )= -APBC3(i3+N*N,N)*cos(xyphi(1+i1+i2+i3+(N*N*N))+APBC3(i3+N*N,N)*theta(1+i1+ i2+ modulo(i3+N*N,N*N*N) )-theta(1+i1+i2+i3))
           
           !Symmetric terms
           HESS(modulo(1+i1,N) +1+i2+i3, 1+i1+i2+i3 ) =      HESS(1+i1+i2+i3,modulo(1+i1,N) +1+i2+i3 ) 
           HESS(1+i1+ modulo(i2+N,N*N)+i3, 1+i1+i2 +i3)=     HESS(1+i1+i2+i3,1+i1+ modulo(i2+N,N*N) +i3) 
           HESS(1+i1+ i2+modulo(i3+N*N,N*N*N), 1+i1+i2+i3 )=     HESS(1+i1+i2+i3,1+i1+ i2+modulo(i3+N*N,N*N*N) ) 

           
           i1=i1+1
        END DO
        i2=i2+N
     END DO
     i3 = i3+ N*N
  ENDDO
  
ENDSUBROUTINE Hessian_3d_APBC


!APBC3 ( for thord direction) is a function to get minus sign for boundary conditions if needed
! theta[i+N] = -theta[i]

FUNCTION APBC3(i,N) 
  IMPLICIT NONE
  DOUBLE PRECISION :: APBC3
  INTEGER, INTENT(IN) :: i, N
  
  if ( (i.EQ.-N*N).OR.(i.EQ.N*N*N) ) THEN
     APBC3 = -1.0D0
  else
     APBC3 = 1.0D0
  ENDIF
END FUNCTION APBC3
