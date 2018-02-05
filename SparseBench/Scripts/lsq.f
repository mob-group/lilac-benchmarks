c
c Least Squares fit through performance points
c
c We assume that performance will level out beyond a certain
c problem size, so we fit a y=a+b/x curve through the data points,
c and report 'a' as the asymptotic performance.
c
      implicit none
      integer ms
      parameter (ms=20)
      double precision x(ms),y(ms),mat(2,2),rhs(2),sol(2)
      data mat,rhs/4*0.d0,2*0.d0/
      integer v,m,n,i,j,ii,off,error
      data error/0/
      double precision d,t
      character*20 name,aspect*7
      logical found

c
c read the problem from stdin
c
      read(5,*,end=998) n
      if (n.lt.2) goto 996
      read(5,*,end=992) v
      off = 0
      do i=1,n
         ii = 1+mod(i-1-off,ms)
         read(5,*,end=997) x(ii),y(ii)
         if (y(ii).eq.0) off = off+1
      end do
      n = n-off
      if (n.eq.0) goto 999
      if (n.lt.2) goto 995
c
c sort
c
      do i=n,2,-1
         do j=1,i-1
            if (x(j+1).lt.x(j)) then
               t = x(j+1)
               x(j+1) = x(j)
               x(j) = t
               t = y(j+1)
               y(j+1) = y(j)
               y(j) = t
            end if
         end do
      end do
c
c find if there is a peak
c
      if (v.gt.0) then
         found = .false.
      else
         do i=1,n
            m = y(i)
            found = .true.
            do ii=i+1,n
               if (y(ii).gt.m) found = .false.
            end do
            if (found) then
               m = i
               goto 10
            end if
         end do
 10      continue
      end if
c
c if there is a peak, discard values before it, since they
c don't fit the y=a+b/x curve
c
      if (found) then
         off = m-1
         do i=m,n
            x(i-off) = x(i)
            y(i-off) = y(i)
         end do
         n = n-off
         if (n.le.3) then
            error = -1
         end if
         if (n.lt.2) goto 995
      end if
      if (error.ge.0) error = n
c
c from dy/da=0 and dy/db=0 we get two equations in
c the unknowns a and b.
c we solve those.
c
      mat(1,1) = n*1.d0
      do i=1,n
         mat(1,2) = mat(1,2)+1.d0/x(i)
         rhs(1) = rhs(1)+y(i)
         mat(2,1) = mat(2,1)+1.d0/x(i)
         mat(2,2) = mat(2,2)+1.d0/(x(i)*x(i))
         rhs(2) = rhs(2)+y(i)/x(i)
      end do

      d = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      t = mat(2,2)
      mat(2,2) = mat(1,1)/d
      mat(1,1) = t/d
      mat(1,2) = -mat(1,2)/d
      mat(2,1) = -mat(2,1)/d

      sol(1) = mat(1,1)*rhs(1)+mat(1,2)*rhs(2)
      sol(2) = mat(2,1)*rhs(1)+mat(2,2)*rhs(2)

c
c the values are written out
c
      write(6,*) y(1),sol(1),error

      goto 1000
 994  write(6,*) 'LSQ error reading first string input: name'
      goto 999
 993  write(6,*) 'LSQ error reading second string input: aspect'
      goto 999
 992  write(6,*) 'LSQ error reading vector int input'
      goto 999
 995  write(6,*) 'LSQ error analyzing ',name,aspect,':'
      write(6,*) '    too few nonzero data points: ',n
      goto 999
 996  write(6,*) 'LSQ error analyzing ',name,aspect,':'
      write(6,*) '    n should be at least 2; found ',n
      goto 999
 998  write(6,*) 'LSQ error analyzing ',name,aspect,':'
      write(6,*) '    could not read n'
      goto 999
 997  write(6,*) 'LSQ error analyzing ',name,aspect,':'
      write(6,*) '    no data for record ',i
      goto 999
 999  continue
      open(7,file='asympt.log')
      write(7,*) 0,0,1
      close(7)
 1000 continue
      end
