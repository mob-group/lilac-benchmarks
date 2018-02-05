      function dotprod(x,y,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),y(*),dotprod
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      real*8 t
      integer i
      
      t = starttimer()
      dotprod = 0.d0
      do i=1,n
         dotprod = dotprod+x(i)*y(i)
      end do
      t = stoptimer()-t
      call add_vec_flops(2*n,t)

      return
      end
c
      function vecnorm(x,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),vecnorm
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      real*8 t
      integer i

      t = starttimer()
      vecnorm = 0.d0
      do i=1,n
         vecnorm = vecnorm+x(i)*x(i)
      end do
      vecnorm = sqrt(vecnorm)
      t = stoptimer()-t
      call add_vec_flops(2*n,0.d0)

      return
      end
C
      subroutine vector_pointwise_multiply(res,x,y,n,time,flops)
      implicit none
C Arguments
      integer n,flops
      real*8 time
      real*8 res(*),x(*),y(*)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer i
      real*8 t

      t = starttimer()
      do i=1,n
         res(i) = x(i)*y(i)
      end do
      time = stoptimer()-t
      flops = n

      return
      end
c
      subroutine x_is_x_plus_ay(x,a,y,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),a,y(*)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer i
      real*8 t

      t = starttimer()
      do i=1,n
         x(i) = x(i) + a * y(i)
      end do
      t = stoptimer()-t
      call add_vec_flops(2*n,t)
         
      return 
      end
c
      subroutine x_is_ax_plus_y(x,a,y,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),a,y(*)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer i
      real*8 t

      t = starttimer()
      do i=1,n
         x(i) = a*x(i) + y(i)
      end do
      t = stoptimer()-t
      call add_vec_flops(2*n,t)
         
      return 
      end
c
      subroutine x_is_ay_plus_bz(x,a,y,b,z,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),a,y(*),b,z(*)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer i
      real*8 t

      t = starttimer()
      do i=1,n
         x(i) = a*y(i) + b*z(i)
      end do
      t = stoptimer()-t
      call add_vec_flops(3*n,t)
         
      return 
      end
c
      subroutine veccopy(x,y,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),y(*)
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer i
      real*8 t

      t = starttimer()
      do i=1,n
         x(i) = y(i)
      end do
      t = stoptimer()-t
      call add_vec_flops(0,t)

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
c
      subroutine vecscale(x,a,y,n)
      implicit none
C Arguments
      integer n
      real*8 x(*),y(*),a
C Externals
      external starttimer,stoptimer
      real*8 starttimer,stoptimer
C Local
      integer i
      real*8 t

      t = starttimer()
      do i=1,n
         x(i) = a*y(i)
      end do
      t = stoptimer()-t
      call add_vec_flops(n,t)

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
c
      subroutine usolve(u,lda,x,y,n,flops)
      implicit none
C Arguments
      integer lda,n,flops
      real*8 u(lda,n),x(n),y(n)
C Local
      integer i,j
      real*8 s
 
      do i=n,1,-1
         s = 0.d0
         do j=i+1,n
            s = s+u(i,j)*y(j)
         end do
         y(i) = (x(i)-s)/u(i,i)
      end do
      flops = flops+ n*(n+1)/2 + 2*n

      return
      end
c
      subroutine vecdump(v,nam,num,len)
      implicit none
C Arguments
      integer len,num
      character*(*) nam
      real*8 v(*)
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      integer i

      write(vecdump_unit,*) nam,num
      do i=1,len
         write(vecdump_unit,*) 'v(',i,') = ',v(i)
      end do

      return
      end
c
      subroutine vecview(v,nam,num,len)
      implicit none
C Arguments
      integer len,num
      character*(*) nam
      real*8 v(*)
C Local
      integer i

      write(6,*) nam,num
      write(6,*) (v(i),i=1,len)

      return
      end
C
C Flop stuff
C
      subroutine add_vec_flops(n,t)
      implicit none
C Arguments
      integer n
      real*8 t
C Vector flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time

      inst = inst+1
      flops = flops+n
      time = time+t

      return
      end
C
      function n_vector_flops()
      implicit none
C Vector flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time
C Local
      real*8 n_vector_flops

      n_vector_flops = flops*1.d-6

      return
      end
C
      subroutine print_vector_flops
      implicit none
C Vector flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time
C Local
      real*8 flop_rate

      if (time.eq.0.d0) then
         flop_rate = 0.d0
      else
         flop_rate = flops/(time*1.d6)
      end if
      write(6,'(
     >     ''Vector operations'',/,
     >     ''Total time:'',4x,f9.3,1x,
     >     ''in'',1x,i4,1x,''instances.'',1x,/,
     >     ''Total Mops:'',1x,f12.3,/,
     >     ''Mflop rate:'',4x,f9.3,/
     >     )') time,inst,flops/1.d6,flop_rate

      return
      end
C
      subroutine vecflops_init
      implicit none
C Vector flop counting
      integer flops,inst
      real*8 time
      common /vecflops/inst,flops
      common /vectime/time

C Initialisation
      inst = 0
      flops = 0
      time = 0.d0

      return
      end
