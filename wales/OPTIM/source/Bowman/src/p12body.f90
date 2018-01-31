!==================================================
! Calculate the monomer and dimer potentials, using hbb only
!==================================================
subroutine pot12bhbb(natm,xx,pot)
  integer,intent(in)::natm
  real,dimension(3,natm),intent(in)::xx
  real,intent(inout)::pot
  !::::::::::::::::::::
  real,dimension(6,3)::x1,x2
  real::e1,e2,vect(3)
  integer::i,j,fo

  fo=natm/3*2
  do i=1,natm/3-1
     x2(1,:)=xx(:,fo+i)        ! O1
     x2(2,:)=xx(:,i*2-1)       ! H1
     x2(3,:)=xx(:,i*2)         ! H2
     do j=i+1,natm/3
        x2(4,:)=xx(:,fo+j)     ! O2
        x2(5,:)=xx(:,j*2-1)    ! H1
        x2(6,:)=xx(:,j*2)      ! H2
        vect(:)=x2(4,:)-x2(1,:)
        x1=x2
        vect(:)=vect(:)*200.d0
        x1(4,:)=x2(4,:)+vect(:)
        x1(5,:)=x2(5,:)+vect(:)
        x1(6,:)=x2(6,:)+vect(:)
        call calcpot(e2,x2) 
        call calcpot(e1,x1) 
        pot=pot+e2-e1*(1.d0-1.d0/dble(natm/3-1))
     end do
  end do

  return
end subroutine pot12bhbb
