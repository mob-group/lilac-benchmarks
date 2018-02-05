      function allocatei(mem,n,ptr,max,txt)
      implicit none
C Arguments
      integer mem(*),n,ptr,max
      character*(*) txt
      integer allocatei

      if (ptr+n+3.gt.max) then
         print *,'Memory overflow in ',txt,'; need at least',ptr+n
         stop
      endif
      mem(ptr) = 4711
      mem(ptr+1) = n
      allocatei = ptr+2
      ptr = ptr+n+2
      mem(ptr) = -4711

      return
      end
C
      subroutine test_memi(mem)
      implicit none
C Arguments
      integer mem(*)
C Local
      integer i

      i = 1
 10   continue
      if (abs(mem(i)).ne.4711) then
         print *,'iMemory corrupted at ',i,': ',mem(i)
         stop
      end if
      if (mem(i).lt.0) return
      i = i+2+mem(i+1)
      goto 10

      end
C
      function allocater(mem,n,ptr,max,txt)
      implicit none
C Arguments
      real*8 mem(*)
      integer n,ptr,max
      character*(*) txt
      integer allocater

      if (ptr+n+3.gt.max) then
         print *,'Memory overflow in ',txt,'; need at least',ptr+n
         stop
      endif
      mem(ptr) = 47.11
      mem(ptr+1) = n
      allocater = ptr+2
      ptr = ptr+n+2
      mem(ptr) = -47.11

      return
      end
C
      subroutine test_memr(mem)
      implicit none
C Arguments
      real*8 mem(*)
C Local
      integer i

      i = 1
 10   continue
      if (abs(mem(i)).ne.47.11) then
         print *,'rMemory corrupted at ',i,': ',mem(i)
         stop
      end if
      if (mem(i).lt.0) return
      i = i+2+mem(i+1)
      goto 10

      end
C
      subroutine imove_left(mem,old,new,n)
      implicit none
C Arguments
      integer mem(*),old,new,n
C Local
      integer i

c      print *,'moving i',n,', from to',old,new
      do i=1,n
         mem(new+i-1) = mem(old+i-1)
      end do

      return
      end
C
      subroutine rmove_left(mem,old,new,n)
      implicit none
C Arguments
      integer old,new,n
      real*8 mem(*)
C Local
      integer i

c      print *,'moving r',n,', from to',old,new
      do i=1,n
         mem(new+i-1) = mem(old+i-1)
      end do

      return
      end
c
      subroutine report_mem_usage(rptr,iptr)
      implicit none
C Arguments
      integer rptr,iptr

      write(6,'(/,''Memory used'')')
      write(6,10)
     >     'Integer',iptr,iptr*4.d0*2.d0**(-20),
     >     'Real',   rptr,rptr*8.d0*2.d0**(-20)
 10   format(4x,a7,':',1x,i10,1x,'=',1x,f5.1,1x,'Mbytes')

      return
      end
