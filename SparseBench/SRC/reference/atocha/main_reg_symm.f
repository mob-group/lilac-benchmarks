      program Generate_regular_matrix
      implicit none
C Data
      integer max_dsize
      parameter (max_dsize=90)
      real*8 a(max_dsize*max_dsize*max_dsize*4),
     >     cof((max_dsize+1)*(max_dsize+1)*(max_dsize+1)*3)
C Local
      integer xside,yside,zside,dsize, matlab

      write(6,*) 'Domain size? (matrix size will be cube of this)'
      read(5,*) dsize
      if (dsize.gt.max_dsize) then
         print *,'Sorry, max is',max_dsize
         stop
      end if
      write(6,*) 'dump to matlab too?'
      read(5,*) matlab
      xside = dsize
      yside = dsize
      zside = dsize
      call seven_point_coefs(cof, xside,yside,zside,.1d0)
      call seven_point_matrix(a,cof, xside,yside,zside)

      call regular_matrix_to_file(a,xside,yside,zside)
      if (matlab.gt.0)
     >     call regular_matrix_to_matlab(a,xside,yside,zside)

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
