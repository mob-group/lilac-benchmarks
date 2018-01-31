        SUBROUTINE GENERATEBASIS(x,y,z,newBasis)

        !This subroutine take a vector defined by coords x,y,z.
        !Normalises that vector to make it a z-axis 
        !Then produces an x-vector and y-vector which are both ortho normal to that.
        !The rotation about the z axis is hard coded by the choice of the new x-axis
 
        DOUBLE PRECISION :: newBasis(1:3,1:3),x,y,z

        !generate a global rotation (input from file)
        newBasis(3,1)=x
        newBasis(3,2)=y
        newBasis(3,3)=z
        newBasis(3,:)=newBasis(3,:)/SQRT(x**2+y**2+z**2)
      
        !generate a vector which is definitely not the z-axis.
        newBasis(1,1)=x
        newBasis(1,2)=y+10
        newBasis(1,3)=z
        
        !project the x-vector and normalise to create an orthonormal vector to z.
        newBasis(1,:) = newBasis(1,:) - (newBasis(3,1)*newBasis(1,1)+newBasis(3,2)*newBasis(1,2)+newBasis(3,3)*newBasis(1,3))*newBasis(3,:)
        newBasis(1,:) = newBasis(1,:)/SQRT(newBasis(1,1)**2+newBasis(1,2)**2+newBasis(1,3)**2)
      
        !Cross product for third vector and make sure its normalised
        newBasis(2,1)=newBasis(3,2)*newBasis(1,3)-newBasis(3,3)*newBasis(1,2)
        newBasis(2,2)=newBasis(3,3)*newBasis(1,1)-newBasis(3,1)*newBasis(1,3)
        newBasis(2,3)=newBasis(3,1)*newBasis(1,2)-newBasis(3,2)*newBasis(1,1)
        newBasis(2,:) = newBasis(2,:)/SQRT(newBasis(2,1)**2+newBasis(2,2)**2+newBasis(2,3)**2)

        write (6,*) 'New basis x:',newBasis(1,1),newBasis(1,2),newBasis(1,3)
        write (6,*) 'New basis y:',newBasis(2,1),newBasis(2,2),newBasis(2,3)
        write (6,*) 'New basis z:',newBasis(3,1),newBasis(3,2),newBasis(3,3)

        a=newBasis(1,1)**2+newBasis(1,2)**2+newBasis(1,3)**2
        b=newBasis(2,1)**2+newBasis(2,2)**2+newBasis(2,3)**2
        c=newBasis(3,1)**2+newBasis(3,2)**2+newBasis(3,3)**2
        d=newBasis(1,1)*newBasis(2,1)+newBasis(1,2)*newBasis(2,2)+newBasis(1,3)*newBasis(2,3)
        e=newBasis(1,1)*newBasis(3,1)+newBasis(1,2)*newBasis(3,2)+newBasis(1,3)*newBasis(3,3)
        f=newBasis(3,1)*newBasis(2,1)+newBasis(3,2)*newBasis(2,2)+newBasis(3,3)*newBasis(2,3)
        
        write (6,*) 'newBasis orthonormality: ',a,b,c,d,e,f

        RETURN
        END SUBROUTINE GENERATEBASIS
