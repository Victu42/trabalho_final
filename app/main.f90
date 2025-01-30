program main

  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
  use calculo_orbitas
  
  implicit none
  real(dp), allocatable :: t(:), r(:), phi(:), vr(:), vphi(:), pr(:)
  
  ! Parâmetros da simulação

  M = 3e36 ! massa do objeto
  tf = 20 ! Tempo para simulação 
  r0 =  4*M! Raio inicial 
  phi0 = 0 ! phi inicial
  vr0 = 0.02*c_UA ! velocidade radial inicial
  vphi0 = 1e-30 ! velocidade angular inicial
  h = 1e-5! tamanho do passo           

  call orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)
  call ver_orbitas(h, tf, r0, phi0, c_m_s, vphi0, t, r, phi, vr, vphi)

end program main

