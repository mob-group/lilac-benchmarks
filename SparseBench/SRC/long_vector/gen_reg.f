C
C Matrix generation routines for 7-diagonal storage.
C These routines are used both in the standalone
C program 'reg_gen' and in the main benchmark program
C
      subroutine seven_point_coefs(cof,n1,n2,n3,u)
      implicit none
C Arguments
      integer n1,n2,n3
      real*8 u
      real*8 cof(n1,n2,n3,-3:3)
C Externals
      external bmrand
      real*8 bmrand
C Local
      integer os1(3),os2(3),os3(3),o1,o2,o3
      data os1,os2,os3/1,0,0, 0,1,0, 0,0,1/
      real*8 r1,r2
      integer i1,i2,i3,p

      call init_rand

      call veczero(cof,3*(n1+1)*(n2+1)*(n3+1))
      do p=1,3
         o1 = os1(p)
         o2 = os2(p)
         o3 = os3(p)
         do i1=1,n1
            do i2=1,n2
               do i3=1,n3
                  r1 = bmrand()
                  r2 = bmrand()*u
                  cof(i1,i2,i3,p) = r1*(1.d0+r2)
                  if (i1+o1.le.n1.and.i2+o2.le.n2.and.i3+o3.le.n3)
     >                 cof(i1+o1,i2+o2,i3+o3,-p) = r1*(1.d0-r2)
               end do
            end do
         end do
      end do

      return
      end
C
      subroutine seven_point_matrix(a,cof,n1,n2,n3)
      implicit none
C Arguments
      integer n1,n2,n3
      real*8 a(n1,n2,n3,-3:3),cof(n1,n2,n3,-3:3)
C Local variables
      integer i1,i2,i3,p
      real*8 c
      logical no1,no2,no3,no

      call veczero(a,7*n1*n2*n3)
      do p=-3,3
         no = .false.
         if (p.ne.0) then
            do i3=1,n3
               no3 = (p.eq.-3.and.i3.eq.1)
     >              .or. (p.eq.3.and.i3.eq.n3)
               do i2=1,n2
                  no2 = (p.eq.-2.and.i2.eq.1)
     >                 .or. (p.eq.2.and.i2.eq.n2)
                  do i1=1,n1
                     no1 = (p.eq.-1.and.i1.eq.1)
     >                    .or.(p.eq.1.and.i1.eq.n1)
                     c = cof(i1,i2,i3,p)
                     if (.not.(no1.or.no2.or.no3))
     >                    a(i1,i2,i3,p) = -c
                     a(i1,i2,i3,0) = a(i1,i2,i3,0) + c
                  end do
               end do
            end do
         end if
      end do
      
      return
      end
C
      subroutine regular_matrix_to_matlab(a,n1,n2,n3)
      implicit none
C Arguments
      integer n1,n2,n3
      real*8 a(n1,n2,n3,-3:3)
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      integer size,i1,i2,i3,row,n12

      print *,'.. dumping to matlab file matdump.m'
      print *,' '
      open(matdump_unit,file='matdump.m')
      size = n1*n2*n3
      n12 = n1*n2
      write(matdump_unit,*) 'A=sparse(',size,',',size,');'
      row = 1
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               write(matdump_unit,*)
     >              'A(',row,',',row,') = ',a(i1,i2,i3,0),';'
               if (i1.lt.n1) write(matdump_unit,*)
     >              'A(',row,',',row+1,') = ',a(i1,i2,i3,1),';'
               if (i2.lt.n2) write(matdump_unit,*)
     >              'A(',row,',',row+n1,') = ',a(i1,i2,i3,2),';'
               if (i3.lt.n3) write(matdump_unit,*)
     >              'A(',row,',',row+n12,') = ',a(i1,i2,i3,3),';'
               if (i1.gt.1) write(matdump_unit,*)
     >              'A(',row,',',row-1,') = ',a(i1,i2,i3,-1),';'
               if (i2.gt.1) write(matdump_unit,*)
     >              'A(',row,',',row-n1,') = ',a(i1,i2,i3,-2),';'
               if (i3.gt.1) write(matdump_unit,*)
     >              'A(',row,',',row-n12,') = ',a(i1,i2,i3,-3),';'
               row = row+1
            end do
         end do
      end do

      return
      end
C
      function can_read_reg_from_file(n1,n2,n3)
      implicit none
C Arguments
      integer n1,n2,n3
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      logical can_read_reg_from_file
      character*10 name
      integer nt

      call reg_filename(n1,name)
      open(matdump_unit,file=name)
      read(matdump_unit,*,end=999) nt,n2,n3
      can_read_reg_from_file = nt.eq.n1
      close(matdump_unit) 

      goto 998
 999  can_read_reg_from_file = .false.
 998  return
      end
C
      subroutine regular_matrix_from_file(mat,n1,n2,n3)
      implicit none
C Arguments
      integer n1,n2,n3
      real*8 mat(n1,n2,n3,-3:3)
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      character*10 name
      integer i1,i2,i3,p

      call reg_filename(n1,name)
      open(matdump_unit,file=name)
      read(matdump_unit,*,end=999) n1,n2,n3

      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               read(matdump_unit,*) (mat(i1,i2,i3,p),p=-3,3)
            end do
         end do
      end do
      close(matdump_unit) 

 999  return
      end
C
      subroutine regular_matrix_to_file(mat,n1,n2,n3)
      implicit none
C Arguments
      integer n1,n2,n3
      real*8 mat(n1,n2,n3,-3:3)
C Dump channels
      integer vecdump_unit,matlab_unit,matdump_unit
      common /dump/vecdump_unit,matlab_unit,matdump_unit
C Local
      character*10 name
      integer i1,i2,i3,p

      if (n1.ne.n2.or.n1.ne.n3) then
         print *,'Sorry, cube domain assumed in dump to file'
         stop
      end if
      call reg_filename(n1,name)
      print *,'.. dumping to file ',name
      print *,' '
      open(matdump_unit,file=name,err=999)

      write(matdump_unit,*) n1,n2,n3
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               write(matdump_unit,*) (mat(i1,i2,i3,p),p=-3,3)
            end do
         end do
      end do
      close(matdump_unit) 

      return
 999  print *,'Could not open file ',name
      stop
      end
C
C File name generation
C
      subroutine reg_filename(size,name)
      implicit none
C Arguments
      integer size
      character*10 name
C Local
      character*10 numb
      data numb/'0123456789'/
      integer i

      name = 'regmat   u'
      i = size
      name(7:7) = numb(i/100+1:i/100+1)
      i = i-100*(i/100)
      name(8:8) = numb(i/10+1:i/10+1)
      i = i-10*(i/10)
      name(9:9) = numb(i+1:i+1)

      return
      end
