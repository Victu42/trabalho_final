program main

  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
  use calculo_orbitas
  
  implicit none
  real(dp), allocatable :: t(:), r(:), phi(:), vr(:), vphi(:), pr(:)
  
  ! Parâmetros da simulação
  M = 1 ! massa do objeto
  tf = 3 ! Tempo para simulação em anos
  r0 = 2.3 ! Raio inicial em U.A.
  phi0 = pi/2 ! phi inicial
  vr0 = 1e-30 !  velocidade radial inicial
  vphi0 = 3.5*pi ! velocidade angular inicial
  h = 1e-4 ! tamanho do passo           

  call orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)
  
end program main
