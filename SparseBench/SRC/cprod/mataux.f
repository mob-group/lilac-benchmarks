      subroutine crs_find_diagonal(idx,val,nnz,ptr,dia,size)
      implicit none
C Arguments
      integer size,nnz
      integer idx(nnz),ptr(size+1),dia(size)
      real*8 val(nnz)
C Local
      integer row,col

      do row=1,size
         do col=ptr(row),ptr(row+1)-1
            if (idx(col).eq.row) then
               if (val(col).le.0.d0) then
                  print *,'nonpos diagonal at row',row
                  stop
               end if
               dia(row) = col
            end if
         end do
      end do

      return
      end
