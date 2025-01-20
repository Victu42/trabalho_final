

  integer :: v, v_phi, M
  real :: a. a_phi
  
  v = 5000
  v_phi = 34
  M = 2*10**30
  a = ((4*M**3)-(4*M**2*r)-(4*M**2*r**2)+(4*M*r**4)*v_phi**2-(r**5*v_phi**2)+r**2*(M-3*M*v**2))/((2*M-r)*r**3)
  a_phi = 2*(-3*M+r)*v*v_phi/((2*M-r)*r)