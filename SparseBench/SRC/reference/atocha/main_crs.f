      program Generate_random_matrix
      implicit none
C Data
c     splits
      integer max_levels,max_blocks,max_var
      parameter (max_levels=20,max_blocks=300000,max_var=5500000/3)
      integer splits(4,max_blocks),b_per_level(max_levels+1),
     >     ordering(max_var),inv_order(max_var)
      integer dsize,nsize,ord,dump,matlab
c     actual matrix
      integer max_elements
      parameter (max_elements=5500000/8)
      integer dia(max_elements),
     >     ptr(max_elements),idx(max_elements),
     >     ptr_f(max_elements),idx_f(max_elements)
      real*8 val(max_elements),val_f(max_elements),tmp(max_var)
      integer n_elements
      real*8 decay,cutoff,d,c,unbalance
C Externals
      external hbw
      integer hbw
C Local
      integer itmp

      write(6,*) 'Domain size? (matrix size will be cube of this)'
      read(5,*) dsize
      nsize = dsize*dsize*dsize
      write(6,*) '.. size=',nsize
      if (nsize.gt.max_var) then
         print *,'Too many variables: ',nsize,' s/b ',max_var
         stop
      end if

c      write(6,*) 'Dimension? (2-3)'
c      read(5,*) d
c      write(6,*) 'Unbalance proportion of the division? (0-1)'
c      read(5,*) unbalance
c      write(6,*) 'Number of edges crossing the cut? (>1)'
c      read(5,*) c
c      write(6,*) 'Decay?'
c      read(5,*) decay
c      write(6,*) 'Cutoff?'
c      read(5,*) cutoff

      write(6,*) 'Order as 1=RecBis 2=CuthMck 3=MultiCol?'
      read(5,*) ord
      write(6,*) 'Dump? (1=yes, 0=no, -1=statistics)'
      read(5,*) dump
      write(6,*) 'Dump to matlab too?'
      read(5,*) matlab

      call generate_crs_matrix(
     >     splits,max_blocks, b_per_level,max_levels,
     >     ptr,idx,val,max_elements,tmp,
     >     d,unbalance,c,decay,cutoff,
     >     dsize,nsize,n_elements)

      if (dump.eq.-1) then
C Print statistics for different orderings
         print *,'Only printing statistics'
         call compute_ordering(ordering,inv_order,2,ptr,idx,
     >        n_elements,nsize,tmp,1)
         itmp = hbw(ptr,idx,n_elements,nsize,ordering,inv_order)
         print *,'Cuthill-McKee: hbw=',itmp
         print *,' '
         call compute_ordering(ordering,inv_order,3,ptr,idx,
     >        n_elements,nsize,tmp,1)
      else
C Permute to correct ordering, sort, and
         print *,'.. permuting to ordering ',ord
         call permute_crs_matrix(ord,ordering,inv_order,tmp,
     >        ptr,idx,val,ptr_f,idx_f,val_f, n_elements,nsize)
         call crs_find_diagonal
     >        (idx_f,val_f,n_elements,ptr_f,dia,nsize)
C Dump to file
         call dump_crs_matrix(dump,matlab,
     >        ptr_f,idx_f,val_f,dia,tmp, n_elements,dsize,nsize)
      end if

      end
c
      subroutine veccopy(x,y,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),y(*)
C Local
      integer i

      do i=1,n
         x(i) = y(i)
      end do

      return
      end
c
      subroutine iveccopy(x,y,n)
      implicit none
C Arguments
      integer n,x(*),y(*)
C Local
      integer i

      do i=1,n
         x(i) = y(i)
      end do

      return
      end
C
      subroutine irsort(iar,ar,len)
      implicit none
C Arguments
      integer len,iar(len)
      real*8 ar(len)
C Local
      integer itmp,i,j
      real*8 double
      
      do i=len-1,1,-1
         do j=1,i
            if (iar(j).gt.iar(j+1)) then
               double = ar(j)
               itmp = iar(j)
               ar(j) = ar(j+1)
               iar(j) = iar(j+1)
               ar(j+1) = double
               iar(j+1) = itmp
            end if
         end do
      end do

      return
      end
c
      subroutine veczero(x,n)
      implicit none
C Arguments
      integer n
      real*8 x(*)
C Local
      integer i

      do i=1,n
         x(i) = 0.d0
      end do

      return
      end
C
C Init common data
C
      block data unit_init
      implicit none
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
      data vecdump_unit,matlab_unit,matdump_unit/7,8,9/
      end
