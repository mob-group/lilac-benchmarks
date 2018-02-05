C----------------------------------------------------------------
C
C This file is part of the sparse benchmark suite.
C Copyright 2000
C Jack Dongarra, Victor Eijkhout, Henk van der Vorst
C
C version 0.9.7
C
C This file last generated:
C Tue Jan 23 13:30:14 EST 2001
C
C----------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Conjugate gradients method
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cg(struct,prec,a,m,rhs,x,block,block_width,
     >     its,rtol,hist,
     >     n1,n2,n3,idx,nnz,ptr,dia,hbw,t)
      implicit none
C Arguments
      integer struct,prec,its,block_width,nnz
      integer n1,n2,n3, idx(nnz),ptr(*),dia(*), hbw
      real*8 a(*),m(*),x(n1*n2*n3),rhs(n1*n2*n3),
     >     block(n1*n2*n3,block_width), rtol,hist(its), t
C Externals
      external dotprod,vecnorm
      real*8 dotprod,vecnorm
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      real*8 rr,rrp,rn,rn0, alpha,beta
      integer p,ap,r,z,zz,tmp,len,i,it
      parameter (p=1,ap=2,r=3,z=4,tmp=5)
      data zz/z/

      if (prec.eq.0) zz = r
      len = n1*n2*n3

c      print *,'Starting Conjugate Gradient Method'
      t = starttimer()
      call matprod(struct,a,x,block(1,tmp),'n',
     >     n1,n2,n3,idx,nnz,ptr)
      do i=1,len
         block(i,r) = block(i,tmp)-rhs(i)
      end do

      do it=1,its
C rn = norm(r); test for convergence
         rn = vecnorm(block(1,r),len)
c         write(6,*) 'Iteration',it,' error =',rn
         hist(it) = rn
         if (it.eq.1) rn0 = rn
         if (rn.lt.rtol*rn0) then
            its = it
            goto 10
         end if
C z = M \ r
         if (prec.gt.0) then
            call prec_solve(struct,prec,a,m,block(1,r),block(1,z),
     >           block(1,tmp),'n',n1,n2,n3,idx,nnz,ptr,dia,hbw)
         end if
C rr_p = rr; rr = z'*r; test for breakdown
         if (it.gt.1) rrp = rr
         rr = dotprod(block(1,r),block(1,zz),len)
         if (rr.le.0) then
            write(6,*) 'Indefinite or nonsymmetric preconditioner:',rr
            stop
         end if
C p = z + beta * p
         if (it.eq.1) then
            call veccopy(block(1,p),block(1,zz),len)
         else
            beta = rr/rrp
            call x_is_ax_plus_y(block(1,p),beta,block(1,zz),len)
         end if
C ap = A*p
         call matprod(struct,a,block(1,p),block(1,ap),'n',
     >        n1,n2,n3,idx,nnz,ptr)
C alpha = r'*r / p'*A*p
         alpha = rr / dotprod(block(1,p),block(1,ap),len)
C x = x - alpha p
C r = r - alpha A*p
         call x_is_x_plus_ay(x,-alpha,block(1,p),len)
         call x_is_x_plus_ay(block(1,r),-alpha,block(1,ap),len)
      end do
      its = -its

 10   continue
      t = stoptimer()-t

      return
      end
C
C Matrix-vector product
C wrapper routine around specific products
C
      subroutine matprod(struct,a,x,y,trans,n1,n2,n3,idx,nnz,ptr)
      implicit none
C Arguments
      integer struct,nnz
      integer n1,n2,n3,idx(nnz),ptr(*)
      real*8 a(*),x(*),y(*)
      character*1 trans

      if (struct.eq.1) then
         call seven_point_matvec(a,n1,n2,n3,x,y)
      else if (struct.eq.2) then
         call random_crs_matprod(a,idx,nnz,ptr,x,y,n1*n2*n3)
      else
         print *,'matvec not implemented for structure',struct
      end if
      
      return
      end
C
C Preconditioner solve
C wrapper routine around specific solves
C
      subroutine prec_solve(struct,prec,a,m,in,out,tmp,trans,
     >     n1,n2,n3,idx,nnz,ptr,dia,hbw)
      implicit none
C Arguments
      integer nnz
      integer struct,prec, n1,n2,n3, idx(nnz),ptr(*),dia(*),hbw
      real*8 in(*),out(*),a(*),m(*),tmp(*)
      character*1 trans
C Local
      integer flops
      real*8 time

      if (prec.eq.1) then
         call vector_pointwise_multiply(out,in,m,n1*n2*n3,time,flops)
         call add_prec_flops(flops)
         call add_prec_time(time)
      else if (prec.eq.2) then
         if (struct.eq.1) then
            if (trans.eq.'t') then
               print *,'Should not happen: transpose solve of symmetric'
               stop
            else
               call seven_point_ilusolve(a,m,in,out,tmp,
     >              n1,n2,n3)
            end if
         else if (struct.eq.2) then
            if (trans.eq.'t') then
               print *,'Should not happen: transpose solve of symmetric'
               stop
            else
               call random_crs_ilusolve(a,m,in,out,tmp,
     >              ptr,idx,nnz,dia, n1*n2*n3)
            end if
         else
            print *,'ILU solve not implemented for structure',struct
            stop
         end if
      else if (prec.eq.3) then
         if (struct.eq.1) then
            if (trans.eq.'t') then
               print *,'Should not happen: transpose solve of symmetric'
               stop
            else 
               call diagonal_bjacobi_solve(a,m,in,out,tmp,
     >              n1,n2,n3)
            end if
         else
            print *,'BJacobi not implemented for structure',struct
         end if
      else if (prec.eq.4) then
         if (struct.eq.1) then
            if (trans.eq.'t') then
               print *,'Should not happen: transpose solve of',
     >                 ' symmetric line block'
               stop
            else
               call lineblock_ilusolve(a,m,in,out,tmp, n1,n2,n3,hbw)
            end if
         else
            print *,'Line ILU not implemented for structure',struct
         end if
      else
         print *,'Preconditioner not implemented',prec
         stop
      end if

      return
      end
C
C Orthogonalisation stuff
C
C Vec_hess: from Krylov sequence construct H from U
C
      subroutine vec_hess(v,wid,len,h,u,ut,
     >     qr_scale,iwork,rwork,lrwork,flops)
      implicit none
C Arguments
      integer wid,len,flops
      integer lrwork, iwork(wid+1)
      real*8 qr_scale(wid+1),rwork(lrwork)
      real*8 v(len,wid+1),h(wid+1,wid),
     >     u(wid+1,wid+1),ut(wid+1,wid+1)
C Externals
      external dotprod,vecnorm
      real*8 dotprod,vecnorm
C Local
      integer i,j,it
      integer info
      real*8 s

      call veczero(u,(wid+1)*(wid+1))
      call veczero(ut,(wid+1)*(wid+1))
      call veczero(h,(wid+1)*wid)
C QR factorisation of V into VU
      call dgeqrf(len,wid+1,v,len,qr_scale,rwork,lrwork, info)
      if (info.ne.0) then
         write(6,'(''DGEQRF error:'',1x,i4)') info
         stop
      end if
      do i=1,wid+1
         do j=i,wid+1
            u(i,j) = v(i,j)
         end do
      end do
	flops = flops + (wid*(wid+1)/2)*4*len + (wid+1)*3*len
C Invert V
      do i=1,wid+1
         do j=i,wid+1
            ut(i,j) = u(i,j)
         end do
      end do
      call dgetrf(wid+1,wid+1,ut,wid+1,iwork, info)
      if (info.ne.0) then
         write(6,'(''DGETRF error:'',1x,i4)') info
         stop
      end if
      call dgetri(wid+1,ut,wid+1,iwork,rwork,lrwork,info)
      if (info.ne.0) then
         write(6,'(''DGETRI error:'',1x,i4)') info
         stop
      end if
      flops = flops + 2*wid*wid*wid/3 * wid*wid/2
C Construct H = U J Uinv
C (the matrices are too small to bother writing this up with BLAS)
      do j=1,wid
         do i=1,wid+1
            s = 0.d0
            do it=2,wid+1
               s = s + u(i,it)*ut(it-1,j)
            end do
            h(i,j) = s
         end do
      end do
      flops = flops + 2*wid*wid*(wid+1)

      return
      end
C
C Flop stuff
C
      subroutine add_mult_flops(flops)
      implicit none
C Arguments
      integer flops
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time

      if (mult_inst.gt.max_inst) then
         print *,'Too many instances'
         stop
      end if
      mult_flops(mult_inst) = mult_flops(mult_inst)+flops

      return
      end
C
      subroutine add_mult_time(time)
      implicit none
C Arguments
      real*8 time
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time

      if (mult_inst.gt.max_inst) then
         print *,'Too many instances'
         stop
      end if
      mult_time(mult_inst) = time
      mult_inst = mult_inst+1
      if (mult_inst.le.max_inst) mult_flops(mult_inst) = 0

      return
      end
C
      subroutine add_prec_flops(flops)
      implicit none
C Arguments
      integer flops
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time

      if (prec_inst.gt.max_inst) then
         print *,'Too many instances'
         stop
      end if
      prec_flops(prec_inst) = prec_flops(prec_inst)+flops

      return
      end
C
      subroutine add_prec_time(time)
      implicit none
C Arguments
      real*8 time
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time

      if (prec_inst.gt.max_inst) then
         print *,'Too many instances'
         stop
      end if
      prec_time(prec_inst) = time
      prec_inst = prec_inst+1
      prec_flops(prec_inst) = 0

      return
      end
C
      function n_matrix_flops()
      implicit none
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time
C Local
      real*8 n_matrix_flops
      integer i

      n_matrix_flops = 0
      do i=1,mult_inst-1
         n_matrix_flops = n_matrix_flops+mult_flops(i)*1.d-6
      end do

      return
      end
C
      function n_prec_flops()
      implicit none
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time
C Local
      real*8 n_prec_flops
      integer i

      n_prec_flops = 0
      do i=1,prec_inst-1
         n_prec_flops = n_prec_flops+mult_flops(i)*1.d-6
      end do

      return
      end
C
      subroutine print_matrix_flops
      implicit none
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time
C Externals
      external n_matrix_flops,n_prec_flops
      real*8 n_matrix_flops,n_prec_flops
C Local
      real*8 time,tot_time,min_time,max_time,tot_flops,flop_rate
      integer i

      min_time = mult_time(1)
      max_time = 0.d0
      tot_time = 0.d0
      tot_flops = n_matrix_flops()
      do i=1,mult_inst-1
         time = mult_time(i)
         max_time = max(max_time,time)
         min_time = min(min_time,time)
         tot_time = tot_time+time
      end do

      if (tot_time.eq.0.d0) then
         flop_rate = 0.d0
      else
         flop_rate = tot_flops/tot_time
      end if

 30   format(
     >     a23,/,
     >     'Total time:',4x,f9.3,1x,'in',1x,i2,1x,'instances.',1x,
     >     '(range:',2(1x,f7.3),')',/,
     >     'Total Mops:',1x,f12.3,/,
     >     'Mflop rate:',4x,f9.3,/
     >     )
      write(6,30) 'Matrix multiply',
     >     tot_time,mult_inst-1,min_time,max_time,
     >     tot_flops,flop_rate

      if (prec_inst.gt.1) then
         min_time = prec_time(1)
         max_time = 0.d0
         tot_time = 0.d0
         tot_flops = n_prec_flops()
         do i=1,prec_inst-1
            time = prec_time(i)
            max_time = max(max_time,time)
            min_time = min(min_time,time)
            tot_time = tot_time+prec_time(i)
         end do
         if (tot_time.eq.0.d0) then
            flop_rate = 0.d0
         else
            flop_rate = tot_flops/tot_time
         end if
         write(6,30) 'Preconditioner solve',
     >        tot_time,prec_inst-1,min_time,max_time,
     >        tot_flops,flop_rate
      end if

      return
      end
C
      subroutine matflops_init
      implicit none
C Matrix flop counting
      integer mult_flops,prec_flops
      integer mult_inst,prec_inst
      real*8 mult_time,prec_time
      integer max_inst
      parameter (max_inst=2*10+3)
      dimension mult_flops(max_inst),prec_flops(max_inst),
     >    mult_time(max_inst),prec_time(max_inst)
      common /matinst/mult_inst,prec_inst
      common /matflops/mult_flops,prec_flops
      common /mattime/mult_time,prec_time
C Local
      integer i

C Initialisation
      mult_inst = 1
      prec_inst = 1
      do i=1,2*10+1
         mult_flops(i) = 0
         prec_flops(i) = 0
         mult_time(i) = 0.d0
         prec_time(i) = 0.d0
      end do

      return
      end
C
      subroutine print_iterative_flops
      implicit none
C Iterative method flop counting
      integer it_inst
      real*8 it_flops,it_time
      common /itinst/it_inst
      common /itflops/it_flops
      common /ittime/it_time
C Local
      real*8 flop_rate

      if (it_inst.eq.0) return
      if (it_time.eq.0.d0) then
         flop_rate = 0.d0
      else
         flop_rate = it_flops/it_time
      end if
      write(6,'(
     >     ''Iterative method'',/,
     >     ''Total time:'',4x,f9.3,1x,
     >     ''in'',1x,i4,1x,''instances.'',1x,/,
     >     ''Total Mops:'',1x,f12.3,/,
     >     ''Mflop rate:'',4x,f9.3,/
     >     )') it_time,it_inst,it_flops,flop_rate

      return
      end
C
      subroutine add_it_flops(flops,time)
      implicit none
C Arguments
      integer flops
      real*8 time
C Iterative method flop counting
      integer it_inst
      real*8 it_flops,it_time
      common /itinst/it_inst
      common /itflops/it_flops
      common /ittime/it_time

      it_inst = it_inst +1
      it_flops = it_flops+flops*1.e-6
      it_time = it_time+time

      return
      end
C
      function n_itflops()
      implicit none
C Iterative method flop counting
      integer it_inst
      real*8 it_flops,it_time
      common /itinst/it_inst
      common /itflops/it_flops
      common /ittime/it_time
C Local
      real*8 n_itflops

      n_itflops = it_flops
      
      return
      end
      subroutine itflops_init
      implicit none
C Iterative method flop counting
      integer it_inst
      real*8 it_flops,it_time
      common /itinst/it_inst
      common /itflops/it_flops
      common /ittime/it_time

C Initialisation
      it_inst = 0
      it_flops = 0.d0
      it_time = 0.d0

      return
      end
