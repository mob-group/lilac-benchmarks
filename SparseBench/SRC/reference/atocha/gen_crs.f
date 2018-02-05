      subroutine generate_crs_matrix(
     >     splits,max_blocks, b_per_level,max_levels,
     >     idx,jdx,val, max_elements,tmp,
     >     d_p,unbalance_p,c_p,decay_p,cutoff_p,
     >     dsize,nsize,nnz)
      implicit none
C Arguments
      integer max_levels,max_blocks, max_elements
      integer dsize,nsize,nnz
      integer splits(4,max_blocks),b_per_level(max_levels+1)
      integer idx(max_elements),jdx(max_elements)
      real*8 val(max_elements),tmp(nsize)
      real*8 unbalance_p,decay_p,cutoff_p,d_p,c_p
C Local
      integer nblocks,nlevels,n_elements
      real*8 d,unbalance,c,decay,cutoff

      call init_rand

C Preset parameters
      d = 3.d0
      unbalance = .3d0
      c = 2.d0
      decay = .1d0
      cutoff = .8d0

C Little sanity check
      if (d.le.0.d0) then
         print *,'Invalid dimension value: ',d
         stop
      end if
      if (d.lt.2.d0.or.d.gt.4.d0)
     >     print *,'Suspicious dimension value: ',d
      print *,'Emulating dimension ',d
      if (unbalance.lt.0.d0.or.unbalance.gt.1.d0) then
         print *,'Invalid unbalance value: ',unbalance
         stop
      end if
      if (c.le.0.d0) then
         print *,'Invalid cut value: ',c
         stop
      end if
      if (c.lt.1)
     >     print *,'Suspiciously low cut value: ',c
      print *,'Cross-boundary ',c
      if (decay.le.0.d0) then
         print *,'Invalid decay value: ',decay
         stop
      end if
      print *,'Decay ',decay
      print *,'Cutoff ',cutoff

C Generate the nested dissection structure
      call make_blocks(splits,max_blocks,b_per_level,max_levels,
     >     unbalance,nsize,nlevels,nblocks)

      call init_rand

C Fill in random elements
      call fill_matrix(splits,nblocks,idx,jdx,val,nsize,
     >     max_elements,decay,cutoff,d,c,n_elements,tmp)
      print *,'.. number of elements=',n_elements
      nnz = n_elements

C Convert to CRS
      print *,'.. converting to crs'
      call convert_to_crs(idx,jdx,val,n_elements,nsize)
      if (n_elements.ne.idx(nsize+1)-1) then
         print *,'Error: nnzeros is/sb',n_elements,idx(nsize+1)-1
         stop
      end if

C Kludge
      call positive_diagonal(jdx,idx,val,n_elements,nsize)

      return
      end
C================================================================
C Permute CRS matrix
C
C this computes an ordering (ord => ordering,inv_order)
C and copies the matrix in permuted form from
C idx,jdx,val into {idx,jdx,val}_f
C================================================================
      subroutine permute_crs_matrix(ord,ordering,inv_order,tmp,
     >     jdx,idx,val,jdx_f,idx_f,val_f, nnz,size)
      implicit none
C Arguments
      integer nnz,size,ord
      integer ordering(size),inv_order(size),
     >     idx(nnz),jdx(size+1),idx_f(nnz),jdx_f(size+1)
      real*8 val(nnz),val_f(nnz),tmp(size)
C Local
      integer irow,row,icol,top,s_loc,t_loc,i

      if (ord.eq.1) then
         call iveccopy(idx_f,idx,nnz)
         call iveccopy(jdx_f,jdx,size+1)
         call  veccopy(val_f,val,nnz)
      else
         call compute_ordering(ordering,inv_order,ord,
     >        idx,jdx, nnz,size,tmp,0)
         do irow=1,size
            row = ordering(irow)
            jdx_f(irow) = jdx(row+1)-jdx(row)
            if (jdx_f(irow).le.0) then
               print *,'Error: zero row or worse',
     >              irow,jdx(row),jdx(row+1)
               stop
            end if
         end do
         top = 1
         do irow=1,size
            top = top+jdx_f(irow)
            jdx_f(irow) = top-jdx_f(irow)
         end do
         jdx_f(size+1) = top
         if (top.ne.jdx(size+1).or.top.ne.nnz+1) then
            print *,'Error: reordered top',top,jdx(size+1),nnz+1
            stop
         end if
         do irow=1,size
            row = ordering(irow)
            do icol=1,jdx(row+1)-jdx(row)
               t_loc = jdx_f(irow)+icol-1
               s_loc = jdx(row)+icol-1
               idx_f(t_loc) = inv_order(idx(s_loc))
               val_f(t_loc) = val(s_loc)
            end do
            call irsort(idx_f(jdx_f(irow)),val_f(jdx_f(irow)),
     >           jdx_f(irow+1)-jdx_f(irow))
         end do 
         print *,(inv_order(i),idx_f(i),val_f(i),
     >        i=jdx_f(1),jdx_f(1+1)-1)
      end if

      return
      end
C
      subroutine dump_crs_matrix(dump,matlab,
     >     ptr,idx,val,dia,tmp, n_elements,dsize,nsize)
      implicit none
C Arguments
      integer dump,matlab,n_elements,dsize,nsize
      integer idx(n_elements),ptr(nsize+1),dia(nsize),tmp(*)
      real*8 val(n_elements)

      call dump_to_crs(ptr,idx,val,dia,n_elements,nsize,dsize)
      if (matlab.gt.0) then
         call dump_to_matlab(ptr,idx,val,n_elements,nsize)
      end if

      return
      end
C================================================================
C Recursive bisection of the matrix
C     The `splits' array contains in location
C     1: begin point 2: split point 3: end point 
C     of each block of variables. Location 4 is the direction (i=1, j=2)
C     in which the variables are numbered.
C================================================================
      subroutine make_blocks(splits,max_blocks,b_per_level,max_levels,
     >     unbalance,n,nlevels,nblocks)
      implicit none
C Arguments
      real*8 unbalance
      integer max_levels,max_blocks
      integer splits(4,max_blocks),b_per_level(max_levels+1)
      integer first,last,blocks,level,top,block,n,nlevels,nblocks
C Externals
      external bmrand
      real*8 bmrand

      splits(1,1) = 1
      splits(2,1) = (1+unbalance*(2*bmrand()-1))*(1+n)/2+1
      splits(3,1) = n
      splits(4,1) = 1
      first = 1
      b_per_level(1) = 1
      top = 1
      do level=1,max_levels
         blocks = b_per_level(level)
         b_per_level(level+1) = 0
         last = first+blocks-1
         if (last.gt.max_blocks) then
            print *,'how did you manage to write at ',last
            stop
         end if
         do block=first,last
            if (splits(2,block)-1.lt.splits(1,block)) then
               print *,'block ',block,' at level ',level,
     >              ' has negative side'
               stop
            else if (splits(2,block)-1.gt.splits(1,block)) then
               top = top+1
               if (top.gt.max_blocks) then
                  print *,'splits overflow: max=',max_blocks
                  stop
               end if
               b_per_level(level+1) = b_per_level(level+1)+1
               splits(1,top) = splits(1,block)
               splits(2,top) = (splits(1,block)+splits(2,block))/2
               splits(3,top) = splits(2,block)-1
               splits(4,top) = 3-splits(4,blocks)
            end if
            if (splits(3,block).lt.splits(2,block)) then
               print *,'block ',block,' at level ',level,
     >              ' has negative side'
               stop
            else if (splits(3,block).gt.splits(2,block)) then
               top = top+1
               if (top.gt.max_blocks) then
                  print *,'splits overflow: max=',max_blocks
                  stop
               end if
               b_per_level(level+1) = b_per_level(level+1)+1
               splits(1,top) = splits(2,block)
               splits(2,top) = (splits(2,block)+splits(3,block)+1)/2
               splits(3,top) = splits(3,block)
               splits(4,top) = 3-splits(4,blocks)
            end if
         end do
         first = last+1
         if (b_per_level(level+1).eq.0) goto 10
      end do
 10   continue
      nlevels = level
      nblocks = top

      return
      end
C================================================================
C Fill the matrix by constructing an (i,j,v) representation.
C This is a sequentialised recursion: for each block we fill
C only the s1-s2 x s2-s3 part; that is, the (1,2) off-diagonal block
C in a 2x2 representation of the submatrix. In the end
C this will fill the whole matrix.
C================================================================
      subroutine fill_matrix(splits,nblocks,idx,jdx,val,size,
     >     max_elements,decay,cutoff,d,c,n_elements,tmp)
      implicit none
C Arguments
      real*8 decay,cutoff,d,c
      integer nblocks,max_elements,size
      integer splits(4,nblocks),idx(max_elements),jdx(max_elements)
      real*8 val(max_elements),tmp(size)
c     output
      integer n_elements
C Externals
      external bmrand,reg_damp,min_i_damp,min_ij_damp
      real*8 bmrand,reg_damp,min_i_damp,min_ij_damp
C Local
      real*8 damp,allow, md,rval
      integer block,top, i1,i2,j1,j2,i,j,itmp,bdry,
     >     hinv_i,hinv_j, allow_i,allow_j, cl,ic,jc

C Create diagonal
      do i=1,size
         idx(i) = i
         jdx(i) = i
         val(i) = bmrand()/size
         tmp(i) = 1.d0
      end do

C Create off-diagonal elements
      top = size
      do block=1,nblocks
         i1 = splits(1,block)
         i2 = splits(2,block)-1
         j1 = splits(2,block)
         j2 = splits(3,block)
         itmp = top
c     special case of small blocks: just put all nonzeros in them.
         if (i1.eq.i2 .or. j1.eq.j2) then
            do i=i1,i2
               do j=j1,j2
                  call fill_elt(idx,jdx,val,block,i,j,top)
               end do
            end do
         else
            hinv_i = int(1+(i2-i1+1.001d0)**(1.d0/d))
            allow_i = hinv_i**(d-1)
            hinv_j = int(1+(j2-j1+1.001d0)**(1.d0/d))
            allow_j = hinv_j**(d-1)
            bdry = min(i2-i1+1-allow_i,j2-j1+1-allow_j)
            allow = allow_i*allow_j
            if (splits(4,block).eq.1) then
               do i=i1,i2
                  md = min_i_damp(i,i1,i2,j1,j2,bdry,decay,top)
                  if (md.gt.1.d0/cutoff) goto 10
                  do j=j1,j2
                     if (top.gt.max_elements) then
                        print *,'Element overflow in block ',block,
     >                       ' of ',nblocks
                        print *,'space available was ',max_elements
                        stop
                     end if
                     damp = reg_damp(splits(4,block),i,j,i1,j1,
     >                    hinv_i,hinv_j,bdry,decay)
                     rval = bmrand()
                     if (rval/damp.gt.cutoff) then
                        call fill_elt(idx,jdx,val,block,i,j,top)
                     end if
                  end do
 10               continue
               end do
            else
               cl = max(1,(i2-i1)/hinv_i)
               do ic=i1,i2,cl
                  do jc=j1,j2,cl
                     md = min_ij_damp(ic,jc,i1,j1,hinv_i,decay)
                     if (md.gt.1.d0/cutoff) goto 20
                     do i=ic,min(i2,ic+cl-1)
                        do j=jc,min(j2,jc+cl-1)
                        if (top.gt.max_elements) then
                           print *,'Element overflow in block ',block,
     >                          ' of ',nblocks
                           print *,'space available was ',max_elements
                           stop
                        end if
                        damp = reg_damp(splits(4,block),i,j,i1,j1,
     >                       hinv_i,hinv_j,bdry,decay)
                        rval = bmrand()
                        if (rval/damp.gt.cutoff) then
                           call fill_elt(idx,jdx,val,block,i,j,top)
                        end if
                        end do
                     end do
 20                  continue
                  end do
               end do
            end if
         end if
      end do
      n_elements = top

      return
      end
c                  if (bmrand()*allow.lt.c) then
c                     tmp(i) = tmp(i)/hinv_i
c                     tmp(j) = tmp(j)/hinv_j
C
      function min_i_damp(i,i1,i2,j1,j2,bdry,decay,t)
      implicit none
C Arguments
      integer i,i1,i2,j1,j2,bdry,t
      real*8 decay,min_i_damp
C Local
      real*8 v
      integer aim

      aim = (i-i1)-bdry
C     satisfy j-j1=aim
      if (aim.ge.0 .and. j1+aim.le.j2) then
         min_i_damp  = 1.d0
      else if (aim.lt.0) then
         v = -aim*.5d0+1.d0
         min_i_damp = v**decay
      else if (j1+aim.gt.j2) then
         v = (j1+aim-j2)*.5d0+1.d0
         min_i_damp = v**decay
      end if

      return
      end
C 
      function min_ij_damp(ic,jc,i1,j1,hinv_i,decay)
      implicit none
C Arguments
      integer ic,jc,i1,j1,hinv_i
      real*8 decay,min_ij_damp
C Externals
      external ceil
      integer ceil
C Local

      min_ij_damp = ( ceil(abs(ic-i1-jc+j1)/(1.d0*hinv_i)) + 1.d0
     >     )**decay

      return
      end
C
      function reg_damp(dir,i,j,i1,j1, hinv_i,hinv_j,bdry,decay)
      implicit none
C Arguments
      integer dir,i,j,i1,j1, hinv_i,hinv_j,bdry
      real*8 decay,reg_damp
C Externals
      external ceil
      integer ceil
C Local
      integer il,jl
      real*8 damp

      if (dir.eq.1) then
c     if the block is split in `i' direction, we put nonzeros
c     with decreasing probability along the line that connects
c     the last block of i-variables (i-i1-bdry) to the first
c     block of j-variables (j-j1).
         damp = (abs((j-j1)-(i-i1) + bdry)*.5d0+1)**decay
      else
c     if the block is split in `j' direction, we put nonzeros
c     with decreasing probability connecting the last i-variable
c     in each column (mod(i-i1,hinv) to the first j-variable
c     (mod(j-j1,hinv)
         il = mod(i-i1,hinv_i)
         jl = mod(j-j1,hinv_i)
         damp = (
     >        ( hinv_i-1-il+jl + ceil(abs(i-i1-j+j1)/(1.d0*hinv_i)) )
     >        +1.d0 )**decay
      end if
      reg_damp = damp

      return
      end
C
      function ceil(x)
      implicit none
C Arguments
      real*8 x
      integer ceil

      if (x.eq.int(x)) then
         ceil = int(x)
      else
         ceil = int(x)+1
      end if

      return
      end
C
      subroutine fill_elt(idx,jdx,val,block,i,j,top)
      implicit none
C Arguments
      integer idx(*),jdx(*),block,i,j,top
      real*8 val(*)
C Externals
      external bmrand
      real*8 bmrand
C Local
      real*8 rval

      rval = 2*bmrand()-1.d0
      top = top+1
      idx(top) = i
      jdx(top) = j
      val(top) = rval
      top = top+1
      idx(top) = j
      jdx(top) = i
      val(top) = rval
      if (i.ne.j) then
         if (block.eq.1) rval = 1.1*rval
         val(i) = val(i)+abs(rval)
         val(j) = val(j)+abs(rval)
      end if

      return
      end
C
C Sort matrix elements by `i' coordinate, and compress then to CRS format
C
      subroutine convert_to_crs(idx,jdx,val,n_elements,size)
      implicit none
C Arguments
      integer n_elements,size
      integer idx(n_elements),jdx(n_elements)
      real*8 val(n_elements)
C Local
      integer el,elt,itmp,col
      real*8 rtmp

C sort by jdx
      call qsijr(jdx,idx,val,n_elements)

      elt = 1
      itmp = 1
      do el=2,n_elements
         if (idx(el).eq.idx(el-1)) then
            itmp = itmp+1
         else
            idx(elt) = itmp
            itmp = 1
            elt = elt+1
         end if
      end do
      idx(elt) = itmp
      if (elt.ne.size) then
         print *,'Oops ',elt,' s/b ',size
         stop
      end if
      itmp = 1
      do el=1,size
         itmp = itmp+idx(el)
         idx(el) = itmp-idx(el)
      end do
      idx(size+1) = itmp

C now sort within each row
      do elt=1,size
         do el=idx(elt+1)-2,idx(elt),-1
            do col=idx(elt),el
               if (jdx(col).gt.jdx(col+1)) then
                  itmp = jdx(col)
                  jdx(col) = jdx(col+1)
                  jdx(col+1) = itmp
                  rtmp = val(col)
                  val(col) = val(col+1)
                  val(col+1) = rtmp
               end if
            end do
         end do
      end do

      return
      end
C
C Compute a suitable permutation from nested dissection to whatever
C
      subroutine compute_ordering(ordering,inv_order,ord,idx,jdx,
     >     n_elements,size,tmp,print_stat)
      implicit none
C Arguments
      integer ord,size,n_elements,print_stat
      real*8 tmp(size)
      integer ordering(size),inv_order(size),
     >     idx(n_elements),jdx(size+1)
C Externals
      external bmrand
      real*8 bmrand
C Local
      integer i,itmp,icol,col,row,top,low,prev, ncolours

      if (ord.eq.1) then
C Nested Dissection: just return what we have
         do i=1,size
            ordering(i) = i
         end do
      else if (ord.eq.2) then
C Cuthill-McKee
         do itmp=1,size
            ordering(itmp) = 0
            inv_order(itmp) = 0
         end do
         ordering(1) = 1
         inv_order(1) = 1
c     last element ordered
         top = 1
c     last element unexamined
         low = 1
c     loop until completion; maximum number of steps is size
         do itmp=1,size
c     restart code: there is no candidate in ordering(low),
c     so find the lowest col s/t inv_order(col)==0
            if (ordering(low).eq.0) then
               do icol=1,size
                  col = icol
                  if (inv_order(icol).eq.0) goto 50
               end do
               stop
 50            continue
               top = low
               ordering(low) = col
               inv_order(col) = 1
            end if
            row = ordering(low)
c     order everything this row is connected to
            do icol=jdx(row),jdx(row+1)-1
               col = idx(icol)
c     if not already ordered
               if (inv_order(col).eq.1) then
                  goto 10
               end if
               top = top+1
               ordering(top) = col
               inv_order(col) = 1
               if (top.eq.size) goto 20
 10            continue
            end do
            low = low+1
            if (low.gt.size) goto 20
         end do
 20      continue
c      print *,'Ordering:',(ordering(row),row=1,size)
c     reverse
      else if (ord.eq.3) then
C Multicolour
         do i=1,size
            tmp(i) = bmrand()
         end do
         top = 0
         ncolours = 0
         do itmp=1,size
            ncolours = ncolours+1
            prev = top
            do row=1,size
               if (tmp(row).eq.0.d0) goto 30
               do col=jdx(row),jdx(row+1)-1
                  if (tmp(idx(col)).gt.tmp(row)) goto 30
               end do
               top = top+1
               ordering(top) = row
               if (top.eq.size) goto 40
 30            continue
            end do
            do i=prev+1,top
               tmp(ordering(i)) = 0.d0
            end do
         end do
 40      continue
         if (print_stat.gt.0)
     >        print *,'Multi-Colour: #colours=',ncolours
      else
         print *,'Unknown ordering ',ord
         stop
      end if

      do row=1,size-1
         do col=row+1,size
            if (ordering(col).eq.ordering(row)) then
               print *,'element occurs twice:',ordering(col),
     >              ' at',row,col
               stop
            end if
         end do
      end do

C Also construct the inverse ordering
      do row=1,size
         inv_order(ordering(row)) = row
      end do

c      print *,'Inverse ordering:',(ordering(row),row=1,size)

      return
      end
C
      subroutine positive_diagonal(idx,jdx,val,nnz,size)
      implicit none
C Arguments
      integer nnz,size
      integer idx(nnz),jdx(size+1)
      real*8 val(nnz)
C Externals
      external bmrand
      real*8 bmrand
C Local
      integer row,col

      do row=1,size
         do col=jdx(row),jdx(row+1)-1
            if (idx(col).eq.row) then
               if (val(col).eq.0.d0) then
                  val(col) = bmrand()
               else if (val(col).lt.0.d0) then
                  print *,'Error: negative diagonal:',row,val(col)
                  stop
               end if
            end if
         end do
      end do

      return
      end
C
C Dump to Matlab. Because of silly restrictions in matlab input
C buffers we have to dump it like this, element by element.
C This takes forever to load.
C
      subroutine dump_to_matlab(jdx,idx,val,n_elements,size)
      implicit none
C Arguments
      integer n_elements,size
      integer idx(n_elements),jdx(size+1)
      real*8 val(n_elements)
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      integer row,col

      print *,'.. dumping to matlab file matdump.m'
      print *,' '
      open(matdump_unit,file='matdump.m')
      write(matdump_unit,*) 'A=sparse(',size,',',size,');'
      do row=1,size
         do col=jdx(row),jdx(row+1)-1
            write(matdump_unit,*)
     >           'A(',row,',',idx(col),')=',val(col),';'
         end do
      end do
      close(matdump_unit)

      return
      end
C
C Dump to some acceptable variant of CRS
C
      subroutine dump_to_crs(jdx,idx,val,dia,n_elements,size,dsize)
      implicit none
C Arguments
      integer size,n_elements,dsize
      integer idx(n_elements),jdx(size+1),dia(size)
      real*8 val(n_elements)
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      character*10 name
      integer row,col, top

      call crs_filename(dsize,name)
      print *,'.. dumping to file ',name
      print *,' '
      open(matdump_unit,file=name)
C write size and number of elements
      write(matdump_unit,*) size,n_elements
C write size+1 row pointers
      top = 1
      do row=1,size
         write(matdump_unit,*) top
         top = top + jdx(row+1)-jdx(row)
      end do
      write(matdump_unit,*) top
C write (idx,val) pair for each element
      do row=1,size
         do col=jdx(row),jdx(row+1)-1
            write(matdump_unit,*) idx(col),val(col)
         end do
      end do
      close(matdump_unit)

      return
      end
C
      function can_read_crs_from_file(size,nnz)
      implicit none
C Arguments
      integer size,nnz
      logical can_read_crs_from_file
C Local
      character *10 name
      integer tsize

      call crs_filename(size,name)
      open(7,file=name)
      read(7,*,end=999) tsize,nnz
      can_read_crs_from_file = tsize.eq.size*size*size
      close(7)

      return
 999  can_read_crs_from_file = .false.
      close(7)
      return
      end
C
      subroutine random_crs_from_file(val,idx,nnz,ptr,dsize,size)
      implicit none
C Arguments
      integer nnz,dsize,size
      real*8 val(nnz)
      integer idx(nnz),ptr(size+1)
C Local
      character*10 name
      integer itmp,jtmp

      call crs_filename(dsize,name)
      open(7,file=name)
      read(7,*) itmp,jtmp
      if (itmp.ne.size .or. jtmp.ne.nnz) then
         print *,'Strange crs read error: ',itmp,' s/b ',size,
     >        '; ',jtmp,' s/b ',nnz
         stop
      end if
      read(7,*) (ptr(itmp),itmp=1,size+1)
      do itmp=1,size
         do jtmp=ptr(itmp),ptr(itmp+1)-1
            read(7,*) idx(jtmp),val(jtmp)
            if (idx(jtmp).lt.1.or.idx(jtmp).gt.size) then
               print *,'Illegal matrix element ',
     >              itmp,jtmp,idx(jtmp),val(jtmp)
               stop
            end if
         end do
      end do
      close(7)
      
      return
      end
C
C Halfbandwidth calculation
C
      function hbw(idx,jdx,n_elements,nsize,ordering,inv_order)
      implicit none
C Arguments
      integer n_elements,nsize,hbw
      integer idx(n_elements),jdx(nsize+1),
     >     ordering(nsize),inv_order(nsize)
C Local
      integer mbw,irow,row,col

      mbw = 0
      do irow=1,nsize
         row = ordering(irow)
c         print *,'row',irow,'; ordered as',row
c         print *,(idx(col),ordering(idx(col)),
c     >        col=jdx(irow),jdx(irow+1)-1)
         do col=jdx(irow),jdx(irow+1)-1
            mbw = max(mbw,abs(inv_order(idx(col))-inv_order(irow)))
         end do
      end do
      hbw = mbw

      return
      end
C
C File name generation
C
      subroutine crs_filename(size,name)
      implicit none
C Arguments
      integer size
      character*10 name
C Local
      character*10 numb
      data numb/'0123456789'/
      integer i

      name = 'crsmat   u'
      i = size
      name(7:7) = numb(i/100+1:i/100+1)
      i = i-100*(i/100)
      name(8:8) = numb(i/10+1:i/10+1)
      i = i-10*(i/10)
      name(9:9) = numb(i+1:i+1)

      return
      end
