!==================================================!
! This subroutine calculates 3-body using Skinner's!
! functions.                                       !
!==================================================!  
subroutine pot_hb3b(nw,x,pot)
  use hb3b_coef
  integer,intent(in)::nw
  real,dimension(3,nw*3),intent(in)::x
  real,intent(out)::pot
  !::::::::::::::::::::
  real,dimension(nw*2,nw)::rr,ff,gg,hh
  integer::fo,i,j,k,jo,ko,ih1,ih2,jh1,jh2,kh1,kh2
  real::fa,fb,fc
  integer::na,nb,nc

  fa=0.d0;fb=0.d0;fc=0.d0
  fo=nw*2
  rr=0.d0

  ! H-bond lengthes
  do i=1,nw*2
     do j=1,nw-1
        jo=conn(i,j)
        rr(i,jo)=sqrt(sum((x(:,i)-x(:,fo+jo))**2))
     end do
  end do

  ! Morse_variables
  ff=exp(-para_a2*rr)
  gg=exp(-para_b2*rr)
  hh=exp(-para_c2*rr)

  !Type-A
  na=0;nb=0;nc=0;
  do i=1,nw
     ih1=i*2-1
     ih2=i*2
     do j=1,nw-2
        jo=conn(ih1,j)
        jh1=jo*2-1
        jh2=jo*2
        do k=j+1,nw-1
           ko=conn(ih1,k)          
           kh1=ko*2-1
           kh2=ko*2
           ! Type-A
           fa=fa+ff(ih1,jo)*ff(ih2,ko)+ff(ih1,ko)*ff(ih2,jo)

           ! Type-B
           fb=fb+gg(ih1,jo)*gg(jh1,ko)+gg(ih2,jo)*gg(jh1,ko)&
                +gg(ih1,jo)*gg(jh2,ko)+gg(ih2,jo)*gg(jh2,ko)&
                +gg(ih1,ko)*gg(kh1,jo)+gg(ih2,ko)*gg(kh1,jo)&
                +gg(ih1,ko)*gg(kh2,jo)+gg(ih2,ko)*gg(kh2,jo)

           ! Type-C
           fc=fc+hh(jh1,i)*hh(kh1,i)+hh(jh1,i)*hh(kh2,i)&
                +hh(jh2,i)*hh(kh1,i)+hh(jh2,i)*hh(kh2,i)

!!$             na=na+2
!!$             nb=nb+8
!!$             nc=nc+4
        end do
     end do
  end do

  pot=(fa*para_a1+fb*para_b1+fc*para_c1)*aukj

  return
end subroutine pot_hb3b
