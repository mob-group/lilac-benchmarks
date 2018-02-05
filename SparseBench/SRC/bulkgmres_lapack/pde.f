      function diff_coef(x1,x2,x3,dir)
      implicit none
      real*8 x1,x2,x3,diff_coef
      integer dir
      diff_coef = 1.d0
      return
      end
C
      function conv_coef(x1,x2,x3,dir)
      implicit none
      real*8 x1,x2,x3,conv_coef
      integer dir
      conv_coef = 3.d-1
      return
      end
