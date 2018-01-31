module hb3b_coef
  implicit none

  ! global variable
  integer,dimension(:,:),allocatable::conn

  real,parameter::aukj=3.808798e-4

  ! Bowman's parameters, unpublished results.
  real,parameter::para_a1=  578.1989081946 !kJ/mol
  real,parameter::para_a2=  1.220088616116 !angstrom^-1
  real,parameter::para_b1= -827.9746906000 !kJ/mol
  real,parameter::para_b2=  1.363990629068 !angstrom^-1
  real,parameter::para_c1=  1000.066749653 !kJ/mol
  real,parameter::para_c2=  1.381253121492 !angstrom^-1

end module hb3b_coef
