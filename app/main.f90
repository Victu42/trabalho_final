program main

  use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
  use calculo_orbitas
  
  implicit none

  ! Parâmetros da simulação
  real(dp) :: tf = 940 ! Tempo para simulação
  real(dp) :: r0 = 60 ! Raio inicial em U.A.
  real(dp) :: phi0 = 0 ! phi inicial
  real(dp) :: vr0 = 0.01*c !  velocidade radial inicial
  real(dp) :: vphi0 = 0.01*pi  ! velocidade angular inicial
  real(dp) :: h = 1e-2 ! tamanho do passo
  real(dp), allocatable :: t(:), r(:), phi(:), vr(:), vphi(:), r_teste(:)
  real(dp) :: vr_valido(20000,2)
  integer(i8) :: i,j, n_steps

  j = 1
  n_steps = int(tf/h)
 

  do while (vphi0 < 6.3)
     vphi0 = vphi0 + 0.01*pi
     vr0 = 0.1*c
     do while (vr0 < c) 
        vr0 = vr0 + 0.01*c
        r_teste = rk4(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)
        if ((isnan(r_teste(n_steps-1)) .eqv. (.false.)) .and. (r_teste(n_steps-1) < 1e10)) then
           vr_valido(j,1) = vr0
           vr_valido(j,2) = vphi0
           j = j + 1
        end if
     end do
  end do
             
  print *, j
  open(unit=10, file='parametros_validos.txt', status='replace')
  do i = 1, int(tf/h) 
     write (10,*) vr_valido
  end do
  close(10)
  
end program main
