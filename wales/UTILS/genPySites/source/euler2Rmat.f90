      SUBROUTINE EULER2RMAT(a,b,g,lRMAT)
      DOUBLE PRECISION :: a,b,g,lRMAT(1:3,1:3)

!      lRMAT(1,1)=COS(b)*COS(g)-COS(a)*SIN(b)*SIN(g)
!      lRMAT(1,2)=-COS(a)*COS(g)*SIN(b) - COS(b)*SIN(g) 
!      lRMAT(1,3)=SIN(a)*SIN(b)
!      lRMAT(2,1)=COS(g)*SIN(b)+COS(a)*COS(b)*SIN(g)
!      lRMAT(2,2)=COS(a)*COS(b)*COS(g)-SIN(b)*SIN(g)
!      lRMAT(2,3)=-COS(b)*SIN(a)
!      lRMAT(3,1)=SIN(a)*SIN(g)
!      lRMAT(3,2)=COS(g)*SIN(a)
!      lRMAT(3,3)=COS(a) 


      lRMAT(1,1)=COS(b)*COS(g) + SIN(a)*SIN(b)*SIN(g)
      lRMAT(1,2)=-COS(a)*SIN(b) 
      lRMAT(1,3)=COS(g)*SIN(a)*SIN(b)-COS(b)*SIN(g)
      lRMAT(2,1)=COS(g)*SIN(b)-COS(b)*SIN(a)*SIN(g)
      lRMAT(2,2)=COS(a)*COS(b)
      lRMAT(2,3)=-COS(b)*COS(g)*SIN(a)-SIN(b)*SIN(g)
      lRMAT(3,1)=COS(a)*SIN(g)
      lRMAT(3,2)=SIN(a)
      lRMAT(3,3)=COS(a)*COS(g) 

      RETURN
      END SUBROUTINE EULER2RMAT
