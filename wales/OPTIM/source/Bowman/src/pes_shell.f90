module pes_shell
  use hb3b_coef,only:conn
  implicit none

  ! Define constants
  real,parameter::emass=1822.88848
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  real,parameter::aukcal=627.51

  ! water functions map
  integer,save::iwaterfcn
  ! 1: hbb12
  ! 2: hbb12+3b
  ! 3: hbb12+hb3b

contains
  !=================================================!
  ! this subroutine is used to initalize the pes    !
  !=================================================!
  subroutine pes_init(nw, dname)
    integer,intent(in)::nw
    character (len=*), intent(in) :: dname
    !::::::::::::::::::::
    integer::i
    integer,dimension(nw)::idx_o

    ! 3-body pot
    call pes_init_3b(dname)

    ! dimer pot
    call prepot()

    ! H-bonded 3b
    allocate(conn(nw*2,nw-1))
    idx_o=(/(i,i=1,nw)/)
    do i=1,nw
       idx_o(i)=0
       conn(i*2-1,:)=pack(idx_o,mask=idx_o.ne.0)
       conn(i*2,:)=conn(i*2-1,:)
       idx_o(i)=i
    end do

    return
  end subroutine pes_init

  !==================================================!
  ! water potentials                                 !
  !==================================================!
  function f(x) result(pot)
    real,dimension(:,:),intent(in)::x
    real::pot
    ! ::::::::::::::::::::
    real,dimension(3,1:size(x,2))::xnew,xxa
    real::p1,p2,p3
    integer::natm,nw,i

    natm=size(x,2); nw=natm/3
    xnew=x
    xxa=x*auang

    p1=0.d0; p2=0.d0; p3=0.d0
    select case(iwaterfcn)
    case(1) ! 1: hbb12
       call pot12bhbb(natm,xnew,p2)
    case(2) ! 2: hbb12+3b
       call pot12bhbb(natm,xnew,p2)
       call pot3b(natm,xnew,p3)          
    case(3) ! 3: hbb12+hb3b
       call pot12bhbb(natm,xnew,p2)
       call pot_hb3b(nw,xxa,p3)
    end select

    pot=p3+p2+p1

    return
  end function f
  
  !==================================================
  ! Calculate the three body potentials
  !==================================================
  subroutine pot3b(natm,xx,pot)
    integer,intent(in)::natm
    real,dimension(3,natm),intent(in)::xx
    real,intent(inout)::pot
    !::::::::::::::::::::
    real,dimension(3,9)::x3
    real::e3
    real,external::fpes
    integer::i,j,k,fo

    fo=natm/3*2;e3=0.d0;pot=0.d0
    do i=1,natm/3-2
       x3(:,7)=xx(:,fo+i)          ! O1
       x3(:,1)=xx(:,i*2-1)         ! H1
       x3(:,2)=xx(:,i*2)           ! H1'
       do j=i+1,natm/3-1
          x3(:,8)=xx(:,fo+j)       ! O2
          x3(:,3)=xx(:,j*2-1)      ! H2
          x3(:,4)=xx(:,j*2)        ! H2'
          do k=j+1,natm/3
             x3(:,9)=xx(:,fo+k)       ! O3
             x3(:,5)=xx(:,k*2-1)      ! H3
             x3(:,6)=xx(:,k*2)        ! H3'
             e3=fpes(x3)      ! three-body correction
             pot=pot+e3
          end do
       end do
    end do

    return
  end subroutine pot3b

end module pes_shell
