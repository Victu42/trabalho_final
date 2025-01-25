program main

  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
  use calculo_orbitas
  
  implicit none
  real(dp), allocatable :: t(:), r(:), phi(:), vr(:), vphi(:), pr(:)
  
  ! Parâmetros da simulação
  tf = 1000 ! Tempo para simulação em anos
  phi0 = 0.0 ! phi inicial
  vr0 = 0.0 !  velocidade radial inicial
  vphi0 = 0.0  ! velocidade angular inicial
  M = G*4e36/c**2 ! massas solares
  r0 = 5*M ! Raio inicial em U.A.
  h = 1 ! tamanho do passo           
  
  !call checar_parametros(M, tf, r0, phi0, vr0, vphi0, h, t, r, phi, vr, vphi)

  call orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)
  
end program main
