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
C BiConjugate gradients method
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bicg(struct,prec,a,m,rhs,x,block,block_width,
     >     its,rtol,hist,
     >     n1,n2,n3,idx,nnz,ptr,dia, hbw, t)
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
      integer p,ap,pl,apl,r,rl,z,zl,zz,zzl,tmp,len,i,it
      parameter (p=1,ap=2,pl=3,apl=4,r=5,rl=6,z=7,zl=8,tmp=9)
      data zz,zzl/z,zl/

      if (prec.eq.0) then
         zz = r
         zzl = rl
      end if
      len = n1*n2*n3

c      print *,'Starting BiConjugate Gradient Method'
      t = starttimer()
      call matprod(struct,a,x,block(1,tmp),'n',
     >     n1,n2,n3,idx,nnz,ptr)
      do i=1,len
         block(i,r) = block(i,tmp)-rhs(i)
         block(i,rl) = block(i,r)
      end do

      do it=1,its
C rn = norm(r); test for convergence
         rn = vecnorm(block(1,r),len)
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
            call prec_solve(struct,prec,a,m,block(1,rl),block(1,zl),
     >           block(1,tmp),'t',n1,n2,n3,idx,nnz,ptr,dia,hbw)
         end if
C rr_p = rr; rr = zz'*r
         if (it.gt.1) rrp = rr
         rr = dotprod(block(1,r),block(1,zzl),len)
C p = z + beta * p
         if (it.eq.1) then
            call veccopy(block(1,p),block(1,zz),len)
            call veccopy(block(1,pl),block(1,zzl),len)
         else
            beta = rr/rrp
            call x_is_ax_plus_y(block(1,p),beta,block(1,zz),len)
            call x_is_ax_plus_y(block(1,pl),beta,block(1,zzl),len)
         end if
C ap = A*p
         call matprod(struct,a,block(1,p),block(1,ap),'n',
     >        n1,n2,n3,idx,nnz,ptr)
         call matprod(struct,a,block(1,pl),block(1,apl),'t',
     >        n1,n2,n3,idx,nnz,ptr)
C alpha = r'*r / p'*A*p
         alpha = rr / dotprod(block(1,pl),block(1,ap),len)
C x = x - alpha p
C r = r - alpha A*p
         call x_is_x_plus_ay(x,-alpha,block(1,p),len)
         call x_is_x_plus_ay(block(1,r),-alpha,block(1,ap),len)
         call x_is_x_plus_ay(block(1,rl),-alpha,block(1,apl),len)
      end do
      its = -its
      
 10   continue
      t = stoptimer()-t
      
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C GMRES: Generalised Minimum Residual method
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine gmres(struct,prec,
     >     a,m,h,q,u,ut,
     >     rhs,x, block,block_width,
     >     restart,its,tol,hist,
     >     n1,n2,n3,len,hbw, idx,nnz,ptr,dia,
     >     qr_scale,iwork,rwork,lrwork, t_global)
      implicit none
C Arguments
      integer nnz,restart,len
      integer lrwork,iwork(restart+1)
      real*8 qr_scale(restart+1),rwork(lrwork)
      integer struct,prec, its,block_width, hbw,
     >     n1,n2,n3, idx(nnz),ptr(len+1),dia(len)
      real*8 a(*),m(*),
     >     h(restart+1,restart+1),q(restart+1,restart+1),
     >     u(restart+1,restart+1),ut(restart+1,restart+1),
     >     rhs(*),x(*), block(len,block_width),
     >     tol,hist(its), t_global
C Externals
      external vecnorm,dotprod
      real*8 vecnorm,dotprod
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer r,tmp1,tmp2,av,mv,mvv,v
      parameter (r=1,tmp1=2,tmp2=3,mv=4,av=5,v=6)
      integer i,j,it,cycle,flops
      real*8 t_ortho,t_sol
      integer info
      real*8 err0,err1,s,est,true_error
      data mvv/mv/

c      print *,'Starting Generalised Minimum Residual Method'
c      print *,'Bulk GmRes'

      t_global = starttimer()
C Initial residual
      call matprod(struct,a,x,block(1,tmp1),'n',
     >     n1,n2,n3,idx,nnz,ptr)
      call x_is_ay_plus_bz(block(1,r),
     >     1.d0,rhs,-1.d0,block(1,tmp1),len)
      err0 = vecnorm(block(1,r),len)

C Restart loop
      do cycle=1,its/restart+1
         err1 = vecnorm(block(1,r),len)
         call vecscale(block(1,v),1.d0/err1,block(1,r),len)
         do it=1,restart
C     next Krylov vector is AM\inv 
            if (prec.gt.0)  then
               call prec_solve(struct,prec,a,m,
     >              block(1,v+it-1),block(1,mv),
     >              block(1,tmp2),'n',n1,n2,n3,idx,nnz,ptr,dia,hbw)
            else
               mvv = v+it-1
            end if
            call matprod(struct,a,block(1,mvv),block(1,v+it),'n',
     >           n1,n2,n3,idx,nnz,ptr)
         end do
C     Orthonormalise Krylov sequence; store coefficients in H
         t_ortho = starttimer()
         flops = 0
         call vec_hess(block(1,v),restart,len,H,u,ut,
     >      qr_scale,iwork,rwork,lrwork,flops)
C     QR decomposition of the Hessenberg matrix
         do it=1,restart
            do i=1,it-1
               s = 0.d0
               do j=1,i+1
                  s = s + q(j,i)*h(j,it)
               end do
               u(i,it) = s
               do j=1,i+1
                  h(j,it) = h(j,it) - u(i,it)*q(j,i)
               end do
            end do
            s = 0.d0
            do i=1,it+1
               s = s+h(i,it)*h(i,it)
            end do
            u(it,it) = sqrt(s)
            s = 1.d0/u(it,it)
            do i=1,it+1
               q(i,it) = s * h(i,it)
            end do
            flops = flops + 2*it*it + 4*(it+1)
         end do
         it = restart
C     orthonormal last column
         do i=2,it+1
            q(i,it+1) = 0.d0
         end do
         q(1,it+1) = 1
         do i=1,it
            s = 0.d0
            do j=1,i
               s = s + q(j,i)*q(j,it+1)
            end do
            do j=1,i+1
               q(j,it+1) = q(j,it+1) - s*q(j,i)
            end do
         end do
         s = 0.d0
         do i=1,it+1
            s = s+q(i,it+1)*q(i,it+1)
         end do
         q(1,it+1) = q(1,it+1)/sqrt(s)
         flops = flops + 2*it*(it+1) + 2*(it+1)
C     error estimate
         est = err1*abs(q(1,it+1))
         hist((cycle-1)*restart+it) = est
         call add_it_flops(flops,stoptimer()-t_ortho)
C     solution and residual update not fully tested.
         flops = 0
         t_sol = starttimer()
C     Solve the coefficients of the V vectors, into block(*,tmp2)
         do i=1,restart
            block(i,tmp1) = q(1,i)*err1
         end do
         call usolve(u,restart+1,
     >        block(1,tmp1),block(1,tmp2),restart,flops)
C     Now take the actual combinations
         flops = flops + 2*restart*len
         if (prec.gt.0) then
            call prec_solve(struct,prec,a,m,
     >           block(1,tmp1),block(1,mv),
     >           block(1,tmp2),'n',
     >           n1,n2,n3,idx,nnz,ptr,dia,hbw)
            mvv = mv
         else
            mvv = tmp1
         end if
         call x_is_x_plus_ay(x,-1.d0,block(1,mvv),len)
         call matprod(struct,a,block(1,mvv),block(1,av),'n',
     >        n1,n2,n3,idx,nnz,ptr)
         call x_is_x_plus_ay(block(1,r),-1.d0,block(1,av),len)
         true_error = vecnorm(block(1,r),len)
         call add_it_flops(flops,stoptimer()-t_sol)
         if (est.lt.tol*err0
     >        .or. (cycle-1)*restart+restart.eq.its) goto 10
      end do
 10   continue
      t_global = stoptimer()-t_global
      write(6,'(''GMRES ended with true error'',1x,f12.5)') true_error

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
         if (trans.eq.'t') then
            call seven_point_matvec_t(a,n1,n2,n3,x,y)
         else
            call seven_point_matvec(a,n1,n2,n3,x,y)
         end if
      else if (struct.eq.2) then
         if (trans.eq.'t') then
            call random_crs_matprod_t(a,idx,nnz,ptr,x,y,n1*n2*n3)
         else
            call random_crs_matprod(a,idx,nnz,ptr,x,y,n1*n2*n3)
         end if
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
               call seven_point_ilusolve_t(a,m,in,out,tmp,
     >              n1,n2,n3)
            else
               call seven_point_ilusolve(a,m,in,out,tmp,
     >              n1,n2,n3)
            end if
         else if (struct.eq.2) then
            if (trans.eq.'t') then
               call random_crs_ilusolve_t(a,m,in,out,
     >              ptr,idx,nnz,dia, n1*n2*n3)
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
               call diagonal_bjacobi_solve_t(a,m,in,out,tmp,
     >              n1,n2,n3)
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
               call lineblock_ilusolve_t(a,m,in,out,tmp, n1,n2,n3,hbw)
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
