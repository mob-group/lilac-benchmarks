! 
! Four-dimensional periodic XY model. ch558
!
SUBROUTINE Energy_4d_PBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE KEY, ONLY : NONEDAPBC, XYPHI ! XYPHI has dim 2N*N, OPTIM reads NATOMS in as N**2 so we have to change it 
  IMPLICIT NONE
  INTEGER N, i1, i2, i3, i4
  DOUBLE PRECISION, dimension((NONEDAPBC)**4) :: theta, GRAD
  DOUBLE PRECISION :: Energy
  LOGICAL GTEST,SECT
    
  N=NONEDAPBC

  Energy=0.0D0

  i4 = 0
  Do WHILE( i4.LT.(N**4))
     i3 = 0
     Do WHILE( i3.LT.(N**3))
        i2=0
        DO WHILE ( i2.LT.(N**2))
           i1=0
           DO WHILE(i1.LT.(N))
              ! 1D part
              Energy = Energy + ( cos(xyphi(i1+i2 + i3+i4+1) + theta( modulo(i1+1,N) + i2 + i3+i4 +1 ) - theta(i1+i2+i3+i4+1)));
              
              ! 2D part          
              Energy = Energy + (cos(xyphi(i1+i2+i3+i4+(N*N)+1) + theta(i1+ modulo(i2+N, N*N) + i3+i4 +1)- theta(i1+i2+i3+i4+1)));
              
              ! 3D part
              Energy = Energy + (cos(xyphi(i1+i2+i3+i4+(N*N*N)+1) + theta(i1+ i2 + modulo(i3 + N*N, N*N*N)+i4 + 1)- theta(i1+i2+i3+i4+1)));

              ! 4D part
              Energy = Energy + (cos(xyphi(i1+i2+i3+i4+(N*N*N*N)+1) + theta(i1+ i2 + i3+ modulo(i4 + N*N*N, N*N*N*N) + 1)- theta(i1+i2+i3+i4+1)));
              
              i1=i1+1
           ENDDO
           i2 = i2 + N
        ENDDO
        i3 = i3 + N**2
     ENDDO
     i4 = i4 + N**3
  ENDDO
     
  Energy = 4.0D0 - (Energy/(N**4))
  !Energy = N*N*N*N*4.0D0 - Energy

  IF (.NOT.(GTEST.OR.SECT)) RETURN
  

  i4 = 0
  Do WHILE( i4.LT.(N**4))
     i3 = 0
     Do WHILE( i3.LT.(N**3))
        i2=0
        DO WHILE ( i2.LT.(N**2))
           i1=0
           DO WHILE(i1.LT.N)      
              
              grad(i1+i2+i3+i4+1)=  &
                   !1D part
                   -sin(xyphi(i1+i2+i3+i4+1)+theta( modulo(i1+1,N) + i2+ i3+i4+1)-theta(i1+i2+i3+i4+1)) &
                   + sin(xyphi(modulo(i1-1,N) + i2 + i3+i4+1)+theta(i1+i2+i3+i4+1)-theta(modulo(i1-1,N) + i2 +i3+i4+1)) & 
                   !2D part
                   - sin( xyphi(i1+i2+i3+i4+ N*N+1) + theta(i1+ modulo(i2+N,N*N) +i3+i4+1)-theta(i1+i2+i3+i4+1)) &
                   + sin(xyphi(i1+ modulo(i2-N, N*N) + i3+i4+N*N+1)+theta(i1+i2+i3+i4+1)-theta(i1+ modulo(i2-N, N*N) +i3+i4+1 )) &                
                   !3D part
                   - sin( xyphi(i1+i2+i3+i4+N*N*N+1) + theta(i1+ i2+modulo(i3+N*N,N*N*N) +i4+1)-theta(i1+i2+i3+i4+1)) &
                   + sin(xyphi(i1+ i2+modulo(i3-N*N, N*N*N) + i4+ N*N*N+1)+theta(i1+i2+i3+i4+1)-theta(i1+ i2+modulo(i3-N*N, N*N*N) +i4+1 ))  &
                   !4D part
                   - sin( xyphi(i1+i2+i3+i4+N*N*N*N+1) + theta(i1+ i2+i3+ modulo(i4+N*N*N,N*N*N*N) +1)-theta(i1+i2+i3+i4+1)) &
                   + sin(xyphi(i1+ i2+i3+ modulo(i4-N*N*N, N*N*N*N) + N*N*N*N+1)+theta(i1+i2+i3+i4+1)-theta(i1+ i2+i3+modulo(i4-N*N*N, N*N*N*N) +1 )) 
              
              i1=i1+1
           ENDDO
           i2=i2+N
        ENDDO
        i3 = i3 + N*N
     ENDDO
     i4 = i4 + N**3
  ENDDO
  
 ! grad(N*N)=0
  
  IF (.NOT.SECT) RETURN
  
  CALL Hessian_4d_PBC(THETA)

END SUBROUTINE Energy_4d_PBC

SUBROUTINE Hessian_4d_PBC(THETA)
  USE KEY, ONLY : NONEDAPBC, XYPHI
  USE MODHESS
  IMPLICIT NONE
  INTEGER :: i1, i2,i3,i4, N
  DOUBLE PRECISION, DIMENSION(NONEDAPBC**4) :: theta
  
  N=NONEDAPBC
  ! THETA(N*N)=0
  HESS(:,:)=0.0D0


  i4 = 0
  Do WHILE( i4.LT.(N**4))
     i3 = 0
     Do WHILE( i3.LT.(N**3))
        i2=0
        DO WHILE ( i2.LT.(N**2))
           i1=0
           DO WHILE(i1.LT.N)      
              
              HESS(1+i1+i2+i3+i4,1+i1+i2+i3+i4)= &
                   !1D Part
                   cos(xyphi(1+i1+i2+i3+i4)+theta( modulo(1+i1,N) +1+i2+i3+i4)-theta(1+i1+i2+i3+i4)) &
                   +cos(xyphi(1+ modulo(i1-1,N) +i2+i3+i4)+theta(1+i1+i2+i3+i4)-theta(1+ modulo(i1-1,N) +i2+i3+i4)) &                
                   !2D part
                   +cos(xyphi(1+i1+i2+i3+i4+(N*N))+theta(1+i1+ modulo(i2+N,N*N)+i3 +i4)-theta(1+i1+i2+i3+i4)) &
                   +cos(xyphi(1+i1+modulo(i2-N,N*N) +(N*N)+i3+i4)+theta(1+i1+i2+i3+i4)-theta(1+i1+ modulo(i2-N,N*N)  +i3+i4))&
                   !3D part
                   +cos(xyphi(1+i1+i2+i3+i4+(N*N*N))+theta(1+i1+i2+ modulo(i3+N*N,N*N*N) +i4)-theta(1+i1+i2+i3+i4)) &
                   +cos(xyphi(1+i1+i2+modulo(i3-N*N,N*N*N) +i4+(N*N*N))+theta(1+i1+i2+i3+i4)-theta(1+i1+ i2+ modulo(i3-N*N,N*N*N)   +i4)) &
                   !4D part
                   +cos(xyphi(1+i1+i2+i3+i4+(N*N*N*N))+theta(1+i1+i2+ i3+modulo(i4+N*N*N, N*N*N*N) )-theta(1+i1+i2+i3+i4)) &
                   +cos(xyphi(1+i1+i2+i3+modulo(i4-N*N*N,N*N*N*N) +(N*N*N*N))+theta(1+i1+i2+i3+i4)-theta(1+i1+i2+i3 +modulo(i4-N*N*N,N*N*N*N)))

              
              !1D part
              HESS(1+i1+i2+i3+i4,modulo(1+i1,N) +1+i2+i3+i4)=  -cos(xyphi(1+i1+i2+i3+i4)+theta(modulo(1+i1,N) +1+i2+i3+i4)-theta(1+i1+i2+i3+i4))
              !2D Part
              HESS(1+i1+i2+i3+i4,1+i1+ modulo(i2+N,N*N)+i3 +i4)=   - cos(xyphi(1+i1+i2+(N*N)+i3+i4)+theta(1+i1+ modulo(i2+N,N*N) +i3+i4)-theta(1+i1+i2+i3+i4))
              !3D Part
              HESS(1+i1+i2+i3+i4,1+i1+ i2+modulo(i3+N*N,N*N*N) +i4)=   - cos(xyphi(1+i1+i2+(N*N*N)+i3+i4)+theta(1+i1+ i2+modulo(i3+N*N,N*N*N)+i4)-theta(1+i1+i2+i3+i4))
              !4D Part
              HESS(1+i1+i2+i3+i4,1+i1+ i2+i3+modulo(i4+N*N*N,N*N*N*N) )=   - cos(xyphi(1+i1+i2+(N*N*N*N)+i3+i4)+theta(1+i1+ i2+i3+modulo(i4+N*N*N,N*N*N*N))-theta(1+i1+i2+i3+i4))
              
              !Symmetric terms
              HESS(modulo(1+i1,N) +1+i2+i3+i4, 1+i1+i2 +i3+i4) =    HESS(1+i1+i2+i3+i4,modulo(1+i1,N) +1+i2+i3+i4) 
              HESS(1+i1+ modulo(i2+N,N*N)+i3+i4, 1+i1+i2 +i3+i4)=  HESS(1+i1+i2+i3+i4,1+i1+ modulo(i2+N,N*N)+i3+i4 )
              HESS(1+i1+ i2+modulo(i3+N*N,N*N*N)+i4, 1+i1+i2 +i3+i4) = HESS(1+i1+i2+i3+i4,1+i1+ i2+modulo(i3+N*N,N*N*N) +i4)
              HESS(1+i1+ i2+i3+modulo(i4+N*N*N,N*N*N*N), 1+i1+i2 +i3+i4)=HESS(1+i1+i2+i3+i4,1+i1+ i2+i3+modulo(i4+N*N*N,N*N*N*N) )
             
              
              i1=i1+1
           END DO
           i2=i2+N
        END DO
        i3 = i3 + N*N
     ENDDO
     i4 = i4 + N**3
  ENDDO
  
ENDSUBROUTINE Hessian_4d_PBC

!                                                             
! Shift zero eigenvalue to the value specified by SHIFTL(1)
!                                                                           
SUBROUTINE SHIFTFOURDPBCT
  USE KEY, ONLY : NONEDAPBC, SHIFTV
  USE MODHESS
  IMPLICIT NONE
  INTEGER :: N,J1,J2
  
  N=NONEDAPBC
  DO J1=1,N*N*N*N
     DO J2=1,N*N*N*N
        HESS(J2,J1)=HESS(J2,J1)+SHIFTV
     ENDDO
  ENDDO
  
END SUBROUTINE SHIFTFOURDPBCT
