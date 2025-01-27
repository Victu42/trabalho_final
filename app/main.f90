program main

  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
  use calculo_orbitas
  
  implicit none
  real(dp), allocatable :: t(:), r(:), phi(:), vr(:), vphi(:), pr(:)
  
  ! Parâmetros da simulação

  M_kg = 1 ! massa do objeto
  tf = 64 ! Tempo para simulação em anos
  r0 = 3*M(M_kg) ! Raio inicial em U.A.
  phi0 = 0. ! phi inicial
  vr0 = 1e-3  ! velocidade radial inicial
  vphi0 = 0 ! velocidade angular inicial
  h = 1e-4 ! tamanho do passo           

  call orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)
  
end program main
