! 
! Two-dimensional periodic XY model. ch558
!
SUBROUTINE Energy_2d_APBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE KEY, ONLY : NONEDAPBC, XYPHI ! XYPHI has dim 2N*N, OPTIM reads NATOMS in as N**2 so we have to change it 
  IMPLICIT NONE
  INTEGER N, i1, i2
  DOUBLE PRECISION, dimension((NONEDAPBC)*(NONEDAPBC)) :: theta, GRAD
  DOUBLE PRECISION :: Energy, APBC1, APBC2
  LOGICAL GTEST,SECT
    
  N=NONEDAPBC
!  theta(N*N)=0
  Energy=0.0D0
  i2=0
  ! indices go from 0 to N-1 for i1
  ! indices go from 0:N:N**2 for i2

  DO WHILE ( i2.LT.(N*N))
     i1=0
     DO WHILE(i1.LT.(N))

        Energy = Energy + ( cos(xyphi(i1+i2+1) + APBC1(i1+1,N)*theta( modulo(i1+1,N) + i2 + 1 ) - theta(i1+i2+1)));

        Energy = Energy + (cos(xyphi(i1+i2+(N*N)+1) + APBC2(i2+N,N)*theta(i1+ modulo(i2+N, N*N) +1)- theta(i1+i2+1)));
        
        i1=i1+1
     ENDDO
     i2 = i2 + N
  END DO

  Energy = 2.0D0 - (Energy/(N*N))
  !Energy = N*N*2.0D0 - Energy

  IF (.NOT.(GTEST.OR.SECT)) RETURN
  
  i2=0
  DO WHILE(i2.LT.(N*N))
     i1=0
     DO WHILE(i1.LT.N)      
        
        grad(i1+i2+1)=  -sin(xyphi(i1+i2+1) + APBC1(i1+1,N)*theta( modulo(i1+1,N) + i2+ 1)-theta(i1+i2+1)) &
             + APBC1(i1-1,N)*sin(xyphi(modulo(i1-1,N) + i2 + 1) + APBC1(i1-1,N)*theta(i1+i2+1)-theta(modulo(i1-1,N) + i2 +1)) & 
             ! we need a -1 boundary condition here, where -1
             ! theta[-1] = theta[N-1] because of indices
             - sin( xyphi(i1+i2+N*N+1) + APBC2(i2+N,N)*theta(i1+ modulo(i2+N,N*N) +1)-theta(i1+i2+1)) &
             + APBC2(i2-N,N)*sin(xyphi(i1+ modulo(i2-N, N*N) + N*N+1)+ APBC2(i2-N,N)*theta(i1+i2+1)-theta(i1+ modulo(i2-N, N*N) +1))
        
        i1=i1+1
     ENDDO
     i2=i2+N
  ENDDO
  
 ! grad(N*N)=0
  
  IF (.NOT.SECT) RETURN
  
  CALL Hessian_2d_APBC(THETA)

END SUBROUTINE Energy_2d_APBC

SUBROUTINE Hessian_2d_APBC(THETA)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  USE MODHESS
  IMPLICIT NONE
  INTEGER :: i1, i2, N
  DOUBLE PRECISION, DIMENSION(NONEDAPBC**2) :: theta
  DOUBLE PRECISION :: APBC1, APBC2

  N=NONEDAPBC
  ! THETA(N*N)=0
  HESS(:,:)=0.0D0

  i2=0
  DO WHILE (i2.LT.(N*N)) 
     i1=0
     DO WHILE(i1.LT.(N))
        
        HESS(1+i1+i2,1+i1+i2)=  cos(xyphi(1+i1+i2)+APBC1(i1+1,N)*theta( modulo(1+i1,N) +1+i2)-theta(1+i1+i2)) &
             +cos(xyphi(1+ modulo(i1-1,N) +i2)+APBC1(i1-1,N)*theta(1+i1+i2)-theta(1+ modulo(i1-1,N) +i2)) &
             +cos(xyphi(1+i1+i2+(N*N))+ APBC2(i2+N,N)*theta(1+i1+ modulo(i2+N,N*N) )-theta(1+i1+i2)) &
             +cos(xyphi(1+i1+modulo(i2-N,N*N) +(N*N))+APBC2(i2-N,N)*theta(1+i1+i2)-theta(1+i1+ modulo(i2-N,N*N)  ))
        

        HESS(1+i1+i2,modulo(1+i1,N) +1+i2)=  - APBC1(i1+1,N)*cos(xyphi(1+i1+i2)+ APBC1(i1+1,N)*theta(modulo(1+i1,N) +1+i2)-theta(1+i1+i2))
        HESS(1+i1+i2,1+i1+ modulo(i2+N,N*N) )= -APBC2(i2+N,N)*cos(xyphi(1+i1+i2+(N*N))+APBC2(i2+N,N)*theta(1+i1+ modulo(i2+N,N*N) )-theta(1+i1+i2))

     HESS(modulo(1+i1,N) +1+i2, 1+i1+i2 ) =      HESS(1+i1+i2,modulo(1+i1,N) +1+i2 ) 
     HESS(1+i1+ modulo(i2+N,N*N), 1+i1+i2 )=     HESS(1+i1+i2,1+i1+ modulo(i2+N,N*N) ) 

        i1=i1+1
   END DO
   i2=i2+N
END DO

ENDSUBROUTINE Hessian_2d_APBC

!APBC1 ( for first direction) is a function to get minus sign for boundary conditions if needed
! theta[i+N] = -theta[i]

FUNCTION APBC1(i,N) 
  IMPLICIT NONE
  DOUBLE PRECISION :: APBC1
  INTEGER, INTENT(IN) :: i, N
  
  if ( (i.EQ.-1).OR.(i.EQ.(N)) ) THEN
     APBC1 = -1.0D0
  else
     APBC1 = 1.0D0
  ENDIF
END FUNCTION APBC1

!APBC2 ( for second direction) is a function to get minus sign for boundary conditions if needed
! theta[i+N] = -theta[i]

FUNCTION APBC2(i,N) 
  IMPLICIT NONE
  DOUBLE PRECISION :: APBC2
  INTEGER, INTENT(IN) :: i, N
  
  if ( (i.EQ.-N).OR.(i.EQ.N*N) ) THEN
     APBC2 = -1.0D0
  else
     APBC2 = 1.0D0
  ENDIF
END FUNCTION APBC2
