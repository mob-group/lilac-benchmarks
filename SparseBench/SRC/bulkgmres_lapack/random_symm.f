C----------------------------------------------------------------
C
C This file is part of the sparse benchmark suite.
C Copyright 2000
C Jack Dongarra, Victor Eijkhout, Henk van der Vorst
C
C version 0.9.7
C
C This file last generated:
C Tue Jan 23 13:30:15 EST 2001
C
C----------------------------------------------------------------
C
C Matrix vector product
C The symmetric storage requires us to do the lower triangle
C by transposing the upper triangle, and working by columns.
C
      subroutine random_crs_matprod(val,idx,nnz,ptr,x,y,size)
      implicit none
C Arguments
      integer nnz,size
      integer idx(nnz),ptr(size+1)
      real*8 val(nnz),x(size),y(size)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer row,col
      real*8 sum
      real*8 t

      t = starttimer()
      do row=1,size
         sum = 0.d0
         do col=ptr(row),ptr(row+1)-1
            sum = sum + val(col)*x(idx(col))
         end do
         y(row) = sum
      end do
      do row=1,size
         do col=ptr(row)+1,ptr(row+1)-1
            y(idx(col)) = y(idx(col)) + val(col)*x(row)
         end do
      end do
      call add_mult_flops(4*(ptr(size+1)-1)-2*size)
      t = stoptimer()-t
      call add_mult_time(t)

      return
      end
C
      subroutine dump_random_crs_matrix(val,idx,nnz,ptr,size)
      implicit none
C Arguments
      integer nnz,size
      real*8 val(nnz)
      integer idx(nnz),ptr(size+1)
C Dump channels
      integer vecdump_unit,matdump_unit
      common /dump/vecdump_unit,matdump_unit
C Local
      integer row,col

      write(matdump_unit,*) 'A=sparse(',size,',',size,');'
      do row=1,size
         do col=ptr(row),ptr(row+1)-1
            write(matdump_unit,*)
     >           'A(',row,',',idx(col),') = ',val(col),';'
            if (idx(col).ne.row) write(matdump_unit,*)
     >           'A(',idx(col),',',row,') = ',val(col),';'
         end do
      end do

      return
      end
C
      subroutine random_crs_jacobi(val,prec,dia,nnz,size)
      implicit none
C Arguments
      integer nnz,size
      integer dia(size)
      real*8 val(nnz),prec(size)
C Local
      integer row

      do row=1,size
         prec(row) = 1.d0/val(dia(row))
      end do

      return
      end
C
      subroutine random_crs_ilufact(val,prec,dia,nnz,size)
      implicit none
C Arguments
      integer nnz,size
      integer dia(size)
      real*8 val(nnz),prec(size)
C Local
      integer row

      do row=1,size
         prec(row) = 1.d0/val(dia(row))
      end do
      call add_fac_flops(size)

      return
      end
C
C Solve a system with an ILU-D factorisation of a CRS matrix.
C The forward solve with symmetric storage is somewhat tricky:
C we can not use the ordinary row algorithm, so we solve
C the upper triangle transposed by columns.
C
      subroutine random_crs_ilusolve(mat,prec,x,y,tmp,
     >     ptr,idx,nnz,dia, size)
      implicit none
C Arguments
      integer nnz
      integer ptr(*),idx(nnz),dia(*),size
      real*8 mat(*),prec(*),x(size),y(size)
     >     ,tmp(size)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer row,col
      real*8 sum
      real*8 t

      t = starttimer()
C Forward solve (D+L)tmp=x
C by columns
      do row=1,size
         tmp(row) = x(row)
      end do
      do row=1,size
         tmp(row) = prec(row)*tmp(row)
         do col=dia(row)+1,ptr(row+1)-1
            tmp(idx(col)) = tmp(idx(col)) - mat(col)*tmp(row)
         end do
      end do

C Backward solve (I+DinvU)y=tmp
      y(size) = tmp(size)
      do row=size-1,1,-1
         sum = 0.d0
         do col=dia(row)+1,ptr(row+1)-1
            sum = sum + mat(col)*y(idx(col))
         end do
         y(row) = tmp(row) - prec(row) * sum
      end do
      call add_prec_flops(4*(ptr(size+1)-1)-size)
      t = stoptimer()-t
      call add_prec_time(t)

      return
      end
