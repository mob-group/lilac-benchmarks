!=================================================!
! this subroutine is used to initalize the        !
!   pes from bas                                  !
!=================================================!
subroutine pes_init_3b(dname)
  use pes, only: pes0_init,pes1_init
  character (len=*) :: dname
  character (len=*), parameter:: pes_x6y3_sysall=&
       'x1 y1 x2 x1y1 y2 x3 x2y1 x1y2 y3 x4 x3y1 x2y2 x1y3 '// &
       'x5 x4y1 x3y2 x2y3 x6 x5y1 x4y2 x3y3 x6y1 x5y2 x4y3 x5y2 x4y3 '// &
       'x6y2 x5y3 x6y3'

  call pes0_init (dir=trim(dname))
  call pes1_init (pes_x6y3_sysall)

  return
end subroutine pes_init_3b
!=========================================================!
! this is the main function to get the potential          !
! energy from fitting code                                !
!=========================================================!
function fpes(x) 
  use pes,only:pes_x6y3_pot
  real,dimension(0:2,0:8),intent(in)::x
  real::fpes
  ! ::::::::::::::::::::

  fpes=pes_x6y3_pot(x)

  return
end function fpes
