        SUBROUTINE GLOBALROTATE(vectors,newBasis,n)

        !This subroutine takes a list of vectors expressed with respect to cartesian reference frame 
        !and returns the co-ordinates of those vectors in a new basis.
        !The new basis is supplied as the row vectors of a 3 by 3 matrix. Those vectors are expressed 
        !in the same reference frame as the original vectors.
        !the new vectors overwrite the old vectors

        INTEGER :: n, curN
        DOUBLE PRECISION :: vectors(n,1:3), vectorTemp(1:3), newBasis(1:3,1:3)

        DO curN=1,n,1
        !     vectorTemp(1)=vectors(curN,1)*newBasis(1,1)+vectors(curN,2)*newBasis(1,2)+vectors(curN,3)*newBasis(1,3)
        !     vectorTemp(2)=vectors(curN,1)*newBasis(2,1)+vectors(curN,2)*newBasis(2,2)+vectors(curN,3)*newBasis(2,3)
        !     vectorTemp(3)=vectors(curN,1)*newBasis(3,1)+vectors(curN,2)*newBasis(3,2)+vectors(curN,3)*newBasis(3,3)
 
        !     write(6,10) vectorTemp(1), vectorTemp(2), vectorTemp(3)
             vectorTemp(:)=matmul(vectors(curN,:),newBasis)
        
             write(6,20) vectorTemp(1), vectorTemp(2), vectorTemp(3)

             vectors(curN,:)=vectorTemp(:)
       END DO
10 FORMAT('mine: ', f10.3, f10.3, f10.3)
20 FORMAT('matmul: ', f10.3, f10.3, f10.3)

        RETURN
        END SUBROUTINE GLOBALROTATE
