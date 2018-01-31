!     This function takes a list of xyz coords and a list of atomic indices 
!     provided by the COLLINDICES array read in after specifying the COLLAGENOP
!     keyword.
!     This computes the dihedral between the four atoms specified
!     in the list of 1-based indices read in after the COLLAGENOP keyword
!     in odata. 
      SUBROUTINE COLLAGENOP_CALC(c,alpha,NATOMS)
         USE KEY
         IMPLICIT NONE

         ! calculates dihedral angle alpha defined by four atoms a0,a1,a2,a3 for
         ! configuration c
         INTEGER NATOMS
         DOUBLE PRECISION v1(3),v2(3),v3(3),r(3),n1(3),n2(3),c(3*NATOMS)

         DOUBLE PRECISION v2norm(3),rdotv2norm,alpha,norm,n1dotn2,normr
         INTEGER a0,a1,a2,a3

         !extract the four indices defined in the odata file with the keyword COLLAGENOP
         a0=COLLINDICES(1)
         a1=COLLINDICES(2)
         a2=COLLINDICES(3)
         a3=COLLINDICES(4)


         !write(*,*) "Collagen Order Parameter from atoms: ",a0,a1,a2,a3
         
         ! calculate current dihedral angle corresponding to particular rotation
         v1(1)=c(3*a1-2)-c(3*a0-2)
         v1(2)=c(3*a1-1)-c(3*a0-1)
         v1(3)=c(3*a1)-c(3*a0  )
         v2(1)=c(3*a2-2)-c(3*a1-2)
         v2(2)=c(3*a2-1)-c(3*a1-1)
         v2(3)=c(3*a2)-c(3*a1)
         v3(1)=c(3*a3-2)-c(3*a2-2)
         v3(2)=c(3*a3-1)-c(3*a2-1)
         v3(3)=c(3*a3  )-c(3*a2)
 
         !WRITE(*,*), "my_dihed initial vectors> ", v1(1),v1(2),v1(3),v2(1),v2(2),v2(3),v3(1),v3(2),v3(3)
 
     
         CALL cross_product(n1,v1,v2)
         CALL cross_product(n2,v2,v3)
      
         CALL cross_product(r,n1,n2)
      
         n1dotn2 = DOT_PRODUCT(n1,n2) 
         norm = DOT_PRODUCT(v2,v2) 
         norm = SQRT(norm)
      
         normr = DOT_PRODUCT(r,r)
         normr = SQRT(normr)
         !v2norm(:) = v2norm(:) / norm
         v2norm(1) = v2(1) / norm
         v2norm(2) = v2(2) / norm
         v2norm(3) = v2(3) / norm
         rdotv2norm = DOT_PRODUCT(r,v2norm)
         alpha = 180.0*ATAN2(rdotv2norm,n1dotn2)/3.14159265358979323846264
         !WRITE(*,*), "my_dihed> ", r(1),r(2),r(3),rdotv2norm,normr
      
      END SUBROUTINE COLLAGENOP_CALC
      
      SUBROUTINE cross_product(c,u,v)
      
         ! calculates cross product c of vectors u and v (3 dimensions)
         IMPLICIT NONE
      
         DOUBLE PRECISION u(3),v(3), c(3)
      
         c(1) = u(2)*v(3)-u(3)*v(2)
         c(2) = u(3)*v(1)-u(1)*v(3)
         c(3) = u(1)*v(2)-u(2)*v(1)
      
      
      END SUBROUTINE cross_product



