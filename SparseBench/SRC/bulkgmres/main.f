C----------------------------------------------------------------
C
C This file is part of the sparse benchmark suite.
C Copyright 2000
C Jack Dongarra, Victor Eijkhout, Henk van der Vorst
C
C version 0.9.7
C
C This file last generated:
C Tue Jan 23 13:30:12 EST 2001
C
C----------------------------------------------------------------
      program sparse_benchmark
      implicit none
C Our own memory management
      integer rsize,isize,rptr,iptr,rptr_save,iptr_save
      data rptr,iptr/1,1/
      parameter (rsize=5500000,isize=2000000)
      real*8 rmemory(rsize)
      integer imemory(isize)
C Problem parameters
      integer
     >     restart,
     >     structure,method,preconditioner,hbw,lodi,
     >     side,xside,yside,zside,vecsize,maxit
      data maxit/10/
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
      external bmrand
      real*8 bmrand
      external allocatei,allocater
      integer allocatei,allocater
      external can_read_reg_from_file,can_read_crs_from_file
      logical can_read_reg_from_file,can_read_crs_from_file
C Data structures
      integer block_width,nnz
      integer mat_ptr,x_ptr,b_ptr,prec_ptr,
     >     tmp_ptr,tmp1_ptr,tmp2_ptr,
     >     ptr_ptr,idx_ptr,dia_ptr, wav_ptr,
     >     mat_ptt,ptr_ptt,idx_ptt,
     >     hess_ptr,qhess_ptr,tri_ptr,ttri_ptr,hist_ptr
c     guard against indexing errors
      data mat_ptr,x_ptr,b_ptr,prec_ptr,
     >     tmp_ptr,tmp1_ptr,tmp2_ptr,
     >     ptr_ptr,idx_ptr,dia_ptr, wav_ptr,
     >     mat_ptt,ptr_ptt,idx_ptt,
     >     hess_ptr,qhess_ptr,tri_ptr,ttri_ptr,hist_ptr
     >     /19*1/
c     init potentially unused variables
      data nnz/0/
C Whole lot of junk for crs construction
      integer max_levels,max_elements,max_blocks,
     >     splits_ptr,bpl_ptr,ord_ptr,inv_ptr
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      real*8 fac_time,t_iter
      integer i

C Identify ourselves to the astaunished onlookers
      call banner

      imemory(1) = -4711
      rmemory(1) = -47.11
C
C Basic input data
C
      write(6,*) 'Input the numer of points along a domain side;'
      write(6,*) 'the matrix size will be the cube of this:'
      read(5,*) side
      if (side.le.0) then
         write(6,'(''Invalid data:'',i5)') side
         goto 999
      endif
      xside = side
      yside = side
      zside = side
      vecsize = xside*yside*zside
      write(6,'(''Dimension'',i9)') vecsize

C
C Matrix construction
C
      write(6,*) 'Choose data structure:'
      write(6,*) '1=regular, 2=irregular'
      read(5,*) structure
      if (structure.lt.1 .or. structure.gt.2) then
         write(6,*) 'Invalid data structure choice:',structure
         goto 999
      endif

      if (structure.eq.1) then
         lodi = -3
         if (.not.can_read_reg_from_file(xside,yside,zside)) then
            print *,' '
            print *,'I can not find a matrix on file for this size'
            print *,'so I will generate it for you;'
            print *,'If you are going to run this problem size again,'
            print *,'you can save time by generating the matrix on file'
            print *,'using the reg_gen utility.'
            print *,' '
            rptr_save = rptr
            tmp_ptr = allocater(rmemory,
     >           7*xside*yside*zside,
     >           rptr,rsize,'reg coefs')
            mat_ptr = allocater(rmemory,
     >           7*vecsize,rptr,rsize,
     >           '7point matrix')
            call seven_point_coefs(rmemory(tmp_ptr),
     >           xside,yside,zside,.1d0)
            call seven_point_matrix(rmemory(mat_ptr),rmemory(tmp_ptr),
     >           xside,yside,zside)
            ptr_ptr = mat_ptr
            rptr = rptr_save
            mat_ptr = allocater(rmemory,
     >           7*vecsize,rptr,rsize,
     >           '7point matrix')
            call rmove_left(rmemory,ptr_ptr,mat_ptr,7*vecsize)
         else
            mat_ptr = allocater(rmemory,
     >           7*vecsize,rptr,rsize,
     >           '7point matrix')
            print *,' '
            print *,'reading the matrix from file.'
            print *,' '
            call regular_matrix_from_file(rmemory(mat_ptr),
     >           xside,yside,zside)
         end if
      else if (structure.eq.2) then
         dia_ptr = allocatei(imemory,
     >        vecsize,iptr,isize,'diagonal pointers')
         if (.not.can_read_crs_from_file(xside,nnz)) then
            print *,' '
            print *,'I can not find a matrix on file for this size'
            print *,'so I will generate it for you;'
            print *,'If you are going to run this problem size again,'
            print *,'you can save time by generating the matrix on file'
            print *,'using the crs_gen utility.'
            print *,' '
            iptr_save = iptr
            rptr_save = rptr
            tmp_ptr = allocater(rmemory,
     >           vecsize,rptr,rsize,'crs tmp')
            max_elements = 15*vecsize
            mat_ptr = allocater(rmemory,
     >           max_elements,rptr,rsize,'crs val t')
            idx_ptr = allocatei(imemory,
     >           max_elements,iptr,isize,'crs idx t')
            ptr_ptr = allocatei(imemory,
     >           max_elements,iptr,isize,'crs jdx t')
            max_blocks = 3*vecsize
            splits_ptr = allocatei(imemory,5*max_blocks,
     >           iptr,isize,'splits')
            max_levels = 500
            bpl_ptr = allocatei(imemory,max_levels+1,iptr,isize,'bpl')
            call generate_crs_matrix(
     >           imemory(splits_ptr),max_blocks,
     >           imemory(bpl_ptr),max_levels,
     >           imemory(idx_ptr),imemory(ptr_ptr),rmemory(mat_ptr),
     >           max_elements, rmemory(tmp_ptr),
     >           0.d0,0.d0,0.d0,0.d0,0.d0,
     >           xside,vecsize,nnz)
            ord_ptr = allocatei(imemory,
     >           vecsize,iptr,isize,'ordering')
            inv_ptr = allocatei(imemory,
     >           vecsize,iptr,isize,'inv ordering')
            mat_ptt = allocater(rmemory,
     >           nnz,rptr,rsize,'crs val')
            idx_ptt = allocatei(imemory,
     >           nnz,iptr,isize,'crs idx')
            ptr_ptt = allocatei(imemory,
     >           vecsize+1,iptr,isize,'crs jdx')
            call permute_crs_matrix(2,
     >           imemory(ord_ptr),imemory(inv_ptr),rmemory(tmp_ptr),
     >           imemory(idx_ptr),imemory(ptr_ptr),rmemory(mat_ptr),
     >           imemory(idx_ptt),imemory(ptr_ptt),rmemory(mat_ptt),
     >           nnz,vecsize)
            iptr = iptr_save
            rptr = rptr_save
            mat_ptr = allocater(rmemory,
     >           nnz,rptr,rsize,'crs val')
            call rmove_left(rmemory,mat_ptt,mat_ptr,nnz)
            idx_ptr = allocatei(imemory,
     >           nnz,iptr,isize,'crs idx')
            call imove_left(imemory,idx_ptt,idx_ptr,nnz)
            ptr_ptr = allocatei(imemory,
     >           vecsize+1,iptr,isize,'crs jdx')
            call imove_left(imemory,ptr_ptt,ptr_ptr,vecsize+1)
         else
            idx_ptr = allocatei(imemory,
     >           nnz,iptr,isize,'crs idx')
            ptr_ptr = allocatei(imemory,
     >           vecsize+1,iptr,isize,'crs jdx')
            mat_ptr = allocater(rmemory,
     >           nnz,rptr,rsize,'crs val')
            print *,' '
            print *,'reading the matrix from file.'
            print *,' '
            call random_crs_from_file(rmemory(mat_ptr),
     >           imemory(idx_ptr),nnz,imemory(ptr_ptr),
     >           xside,vecsize)
         end if
         call crs_find_diagonal(
     >        imemory(idx_ptr),rmemory(mat_ptr),nnz,imemory(ptr_ptr),
     >        imemory(dia_ptr),vecsize)
      else
         print *,'Sorry, matrix',structure,
     >        'unknown or not yet implemented'
         goto 999
      endif

      call test_memi(imemory)
      call test_memr(rmemory)

C Flops initialisation for factorisation
      call facflops_init

C
C Preconditioner
C
      write(6,*) 'Choose the preconditioner:'
      write(6,*) '0=None, 1=Jacobi, 2=ILU'
      if (structure.eq.1)
     >     write(6,*) ' 3=Block Jacobi, 4=Block ILU'
      read(5,*) preconditioner
      if (preconditioner.lt.0
     >    .or. (structure.eq.1 .and. preconditioner.gt.3)
     >    .or. (structure.ne.1 .and. preconditioner.gt.2) ) then
         write(6,*) 'Invalid preconditioner choice:',preconditioner
         goto 999
      endif

      if (preconditioner.eq.1) then
         prec_ptr = allocater(rmemory,
     >        vecsize,rptr,rsize,'preconditioner')
         if (structure.eq.1) then
            fac_time = starttimer()
            call seven_point_jacobi(rmemory(mat_ptr),
     >           rmemory(prec_ptr), xside,yside,zside)
            fac_time = stoptimer()
            call add_fac_time(fac_time)
         else if (structure.eq.2) then
            fac_time = starttimer()
            call random_crs_jacobi(rmemory(mat_ptr),rmemory(prec_ptr),
     >           imemory(dia_ptr), nnz,vecsize)
            fac_time = stoptimer()
            call add_fac_time(fac_time)
         else
            print *,'Jacobi not implemented for structure',structure
            stop
         end if
      else if (preconditioner.eq.2) then
         prec_ptr = allocater(rmemory,
     >        vecsize,rptr,rsize,'preconditioner')
         if (structure.eq.1) then
            fac_time = starttimer()
            call seven_point_ilufact(rmemory(mat_ptr),
     >           rmemory(prec_ptr), xside,yside,zside)
            fac_time = stoptimer()
            call add_fac_time(fac_time)
         else if (structure.eq.2) then
            prec_ptr = allocater(rmemory,
     >           vecsize,rptr,rsize,'preconditioner')
            fac_time = starttimer()
            call random_crs_ilufact(rmemory(mat_ptr),rmemory(prec_ptr),
     >           imemory(dia_ptr), nnz,vecsize)
            fac_time = stoptimer()
            call add_fac_time(fac_time)
         else
            print *,'ILU not implemented for structure',structure
            stop
         end if
      else if (preconditioner.eq.3) then
         if (structure.eq.1) then
            prec_ptr = allocater(rmemory,
     >           vecsize,rptr,rsize,'preconditioner')
            call seven_point_jacobi(rmemory(mat_ptr),
     >           rmemory(prec_ptr), xside,yside,zside)
         else
            print *,'Sorry, preconditioner not available for structure'
            goto 999
         end if
      else if (preconditioner.eq.4) then
         if (structure.eq.1) then
            hbw = 3
            prec_ptr = allocater(rmemory,
     >           (2*hbw+1)*vecsize,rptr,rsize,
     >           'preconditioner')
            tmp1_ptr = allocater(rmemory,
     >           (2*hbw+1)*xside,rptr,rsize,
     >           'prec tmp space 1')
            tmp2_ptr = allocater(rmemory,
     >           (2*hbw+1)*xside,rptr,rsize,
     >           'prec tmp space 2')
            fac_time = starttimer()
            call line_ilufact(rmemory(mat_ptr),rmemory(prec_ptr),
     >           rmemory(tmp1_ptr),rmemory(tmp2_ptr),
     >           xside,yside,zside,hbw)
            fac_time = stoptimer()
            call add_fac_time(fac_time)
         else
            print *,'Sorry, Block ILU not available for structure',
     >           structure
            stop
         end if
      else if (preconditioner.ne.0) then
         print *,'Sorry, preconditioner',preconditioner,
     >        ' unknown or not yet implemented'
         goto 999
      end if

C
C Vector allocation
C
      x_ptr = allocater(rmemory,
     >     vecsize,rptr,rsize,'x')
      b_ptr = allocater(rmemory,
     >     vecsize,rptr,rsize,'rhs')

      call test_memi(imemory)
      call test_memr(rmemory)

C
C Iterative method
C
      write(6,*) 'Choose the iterative method:'
      write(6,*) '1=BiCG; 2=GMRES(r)'
      read(5,*) method
      if (method.eq.1) then
         write(6,*) 'Iterative method BiCG chosen:'
      else if (method.eq.2) then
         write(6,*) 'Iterative method GMRES chosen:'
      else
         write(6,*) 'Invalid iterative method choice:',method
         goto 999
      endif
      restart = maxit

C make initial guess and right hand side
C sophisticated, huh?
      do i=1,vecsize
         rmemory(x_ptr+i-1) = 0.d0
         rmemory(b_ptr+i-1) = 1.d0
      end do

C Flops initialisation for iterative method
      call vecflops_init
      call matflops_init
      call itflops_init

C allocate temps, and call method
      hist_ptr = allocater(rmemory,
     >     10,rptr,rsize,'convergence history')
      if (method.eq.1) then
         block_width = 9
         tmp_ptr = allocater(rmemory,
     >        block_width*vecsize,rptr,rsize,
     >        'bicg temp space')
         call bicg(structure,preconditioner,
     >        rmemory(mat_ptr),rmemory(prec_ptr),
     >        rmemory(b_ptr),rmemory(x_ptr),
     >        rmemory(tmp_ptr),block_width,
     >        maxit,1.d-6,rmemory(hist_ptr),
     >        xside,yside,zside,
     >        imemory(idx_ptr),nnz,imemory(ptr_ptr),imemory(dia_ptr),
     >        hbw, t_iter)
      else if (method.eq.2) then
         block_width = 5+restart+1
         tmp_ptr = allocater(rmemory,
     >        block_width*vecsize,rptr,rsize,
     >        'gmres temp space')
         hess_ptr = allocater(rmemory,
     >        (restart+1)*(restart+1),rptr,rsize,
     >        'hessenberg matrix')
         qhess_ptr = allocater(rmemory,
     >        (restart+1)*(restart+1),rptr,rsize,
     >        'H q factor')
         tri_ptr = allocater(rmemory,
     >        (restart+1)*(restart+1),rptr,rsize,
     >        'H tri factor')
         ttri_ptr = allocater(rmemory,
     >        (restart+1)*(restart+1),rptr,rsize,
     >        'Temp tri')
         call gmres(structure,preconditioner,
     >        rmemory(mat_ptr),rmemory(prec_ptr),
     >        rmemory(hess_ptr),rmemory(qhess_ptr),
     >        rmemory(tri_ptr),rmemory(ttri_ptr),
     >        rmemory(b_ptr),rmemory(x_ptr),
     >        rmemory(tmp_ptr),block_width, restart,maxit,1.d-6,
     >        rmemory(hist_ptr),
     >        xside,yside,zside,vecsize,hbw,
     >        imemory(idx_ptr),nnz,imemory(ptr_ptr),imemory(dia_ptr),
     >        t_iter)
      else
         print *,'Sorry, iterative method',method,
     >        ' unknown or not yet implemented'
         goto 999
      end if
      call history_print(rmemory(hist_ptr),abs(maxit))
      call report_mem_usage(rptr,iptr)

      call test_memi(imemory)
      call test_memr(rmemory)

      write(6,*) ' '
c      call factor_flops()
      call print_matrix_flops()
      call print_iterative_flops()
      call print_vector_flops()
      call print_total_flops(t_iter)

 999  continue
      end
      subroutine banner
      write(6,*) ' '
      write(6,*) '**************************************************'
      write(6,*) '**                                              **'
      write(6,*) '**  Iterative methods benchmark                 **'
      write(6,*) '**  copyright 2000                              **'
      write(6,*) '**  Jack Dongarra                               **'
      write(6,*) '**  Victor Eijkhout                             **'
      write(6,*) '**  Henk van der Vorst                          **'
      write(6,*) '**                                              **'
      write(6,*) '**  general code                                **'
      write(6,*) '**  -- bulk orthogonalisation optimisation      **'
      write(6,*) '**     in gmres                                 **'
      write(6,*) '**                                              **'
      write(6,*) '**************************************************'
      write(6,*) ' '
      return
      end
C
      subroutine history_print(hist,n)
      implicit none
C Arguments
      integer n
      real*8 hist(n)
C Local
      integer i

      do i=1,n
         if (hist(i).ne.0.d0)
     >        write(6,'(
     >        ''Iteration:'',1x,i2,1x,''residual norm:'',1x,e14.7)
     >        ') i,hist(i)
      end do

      return
      end
C
C Utility routines
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
C
C Flop stuff
C
      subroutine facflops_init
      implicit none
C Factorisation flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time
C Initialisation
      flops = 0
      time = 0.d0
      return
      end
C
      subroutine add_fac_flops(n)
      implicit none
C Arguments
      integer n
C Factor flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time

      flops = flops+n

      return
      end
C
      subroutine add_fac_time(t)
      implicit none
C Arguments
      real*8 t
C Factor flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time

      time = time+t

      return
      end
C
      subroutine print_total_flops(ttotal)
      implicit none
C Arguments
      real*8 ttotal
C Externals
      external n_matrix_flops,n_prec_flops,n_vector_flops,
     >     n_itflops
      real*8 n_matrix_flops,n_prec_flops,n_vector_flops,
     >     n_itflops
C Local
      real*8 flops,flop_rate

      flops = n_matrix_flops()+n_prec_flops()+n_vector_flops()
     >     +n_itflops()
      if (ttotal.eq.0.d0) then
         flop_rate = 0.d0
      else
         flop_rate = flops/ttotal
      end if
      write(6,'(
     >     ''Overall statistics'',/,
     >     ''Total time:'',4x,f9.3,/,
     >     ''Total Mops:'',4x,f9.3,/,
     >     ''Mflop rate:'',4x,f9.3
     >     )') ttotal,flops,flop_rate

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
