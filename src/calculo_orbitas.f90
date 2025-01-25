module calculo_orbitas

use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
implicit none

real(dp), parameter :: pi = 16*ATAN(1./5.) - 4*ATAN(1./239.)
real(dp), parameter :: c = 299792458 ! m/s
real(dp) :: M, tf, r0, phi0, vr0, vphi0, h
real(dp), parameter :: G = 6.670430E-11 


contains
  
  subroutine checar_parametros(M, tf, r0, phi0, vr0, vphi0, h, t, r, phi, vr, vphi)

  real(dp) :: M, vphi0, vr0, h, tf, r0, phi0
  real(dp), allocatable :: pr(:), t(:), r(:), phi(:), vr(:), vphi(:)
  real(dp), allocatable :: vr_valido(:), vphi_valido(:), M_valido(:), r_valido(:), temp_vphi(:), temp_vr(:), temp_M(:), temp_r(:)
  integer(i8) :: n_steps, j

  j = 1
  allocate(vr_valido(j), vphi_valido(j), M_valido(j), r_valido(j))
  n_steps = int(tf/h)
  open(unit=11, file='parametros_validos.txt', status='replace')
  do while (r0 < 1000)
      do while (M < 4e6)
         do while (vphi0 < 4*pi)
            do while (vr0 < c)
               pr = rk4(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)
               if ((isnan(pr(n_steps-1)) .eqv. .false.)) then
                  allocate(temp_vphi(j), temp_vr(j), temp_M(j), temp_r(j))
                  vr_valido(j) = vr0
                  vphi_valido(j) = vphi0
                  M_valido(j) = M
                  r_valido(j) = r0
                  write (11,*) vr_valido(j), vphi_valido(j), M_valido(j), r_valido(j)
                  temp_vphi = vphi_valido
                  temp_vr = vr_valido 
                  temp_M = M_valido
                  temp_r = r_valido
                  deallocate(vr_valido)
                  deallocate(vphi_valido) 
                  deallocate(M_valido)
                  deallocate(r_valido)
                  j = j + 1
                  allocate(vr_valido(j), vphi_valido(j), M_valido(j), r_valido(j))
                  vphi_valido(1:j-1) = temp_vphi(1:j-1)
                  vr_valido(1:j-1) = temp_vr(1:j-1)
                  M_valido(1:j-1) = temp_M(1:j-1)
                  r_valido(1:j-1) = temp_r(1:j-1)
                  deallocate(temp_vr)
                  deallocate(temp_vphi) 
                  deallocate(temp_M)
                  deallocate(temp_r)
               end if
               vr0 = vr0 + 0.001*c
            end do
            vphi0 = vphi0 + 0.01*pi
            vr0 = 0.001*c
         end do 
         M = M + 5
         vphi = 0.01*pi
      end do
      M = 1
      r0 = r0 + 10
  end do
             
  close(11)
     

  end subroutine checar_parametros

  function ar_luz(r,dphi) result(f) 
    
    real(dp) :: r, dphi, f

    f = ((-4*M**2)+(2*M*r)+(r-5*M)*(r**3)*(dphi**2))/r**3

  end function ar_luz

  function ar(r, dr, dphi) result(f) 
    
    real(dp) :: r, dphi, f, dr

    f = ((4*(M**3))-(4*r*(M**2))-(4*(M**2)*(r**3)*(dphi**2))+(4*M*(r**4)*(dphi**2))-((r**5)*dphi**2)+(r**2)* &
    (M-3*M*(dr**2)))/((2*M-r)*(r**3))

  end function ar

  function aphi(r,dr,dphi) result(f)
     
    real(dp) :: r, dr, dphi, f
  
    f = 2*(-3*M+r)*dr*dphi/((2*M-r)*r)

  end function aphi

  function rk4(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi) result(fn)

   real(dp), intent(in) :: tf, r0, phi0, vphi0, vr0, h
   real(dp) :: k1_r, k1_phi, k1_vr, k1_vphi, k2_phi, k2_r, k2_vphi, k2_vr,k3_phi,k3_r,k3_vphi,k3_vr,k4_phi,k4_r,k4_vphi, &
               k4_vr 
   real(dp), allocatable, intent(out) :: t(:), r(:), phi(:), vr(:), vphi(:)
   real(dp), allocatable :: fn(:)
   integer(i8) :: i, steps

   steps =  int(tf/h)
   allocate(t(steps),r(steps),phi(steps),vr(steps),vphi(steps), fn(steps))
   t(1) = 0
   r(1) = r0
   phi(1) = phi0
   vr(1) = vr0
   vphi(1) = vphi0

   do i = 1, steps - 1
    
    if (i /= 1) then
       t(i) = t(i-1) + tf/steps 
    end if

    k1_r = r(i)
    k1_vr = vr(i)
    k1_phi = phi(i)
    k1_vphi = vphi(i)
    k2_r = vr(i)+h*k1_vr/2
    k2_vr = ar(r(i)+h*k1_r/2,vr(i)+h*k1_vr/2,vphi(i)+h*k1_vphi/2)
    k2_phi = vphi(i)+h*k1_vphi/2
    k2_vphi = aphi(r(i)+h*k1_r/2,vr(i)+h*k1_vr/2,vphi(i)+h*k1_vphi/2)
    k3_r = vr(i)+h*k2_vr/2
    k3_vr = ar(r(i)+h*k2_r/2,vr(i)+h*k2_vr/2,vphi(i)+h*k2_vphi/2)
    k3_phi = vphi(i)+h*k2_vphi/2
    k3_vphi = aphi(r(i)+h*k2_r/2,vr(i)+h*k2_vr/2,vphi(i)+h*k2_vphi/2)
    k4_r = vr(i)+h*k3_vr
    k4_vr = ar(r(i)+h*k3_r,vr(i)+h*k3_vr,vphi(i)+h*k3_vphi)
    k4_phi = vphi(i)+h*k3_vphi
    k4_vphi = aphi(r(i)+h*k3_r,vr(i)+h*k3_vr,vphi(i)+h*k3_vphi)

    r(i+1) = r(i) + h*(k1_r+2*k2_r+2*k3_r+k4_r)/6
    vr(i+1) = (vr(i) + h*(k1_vr+2*k2_vr+2*k3_vr+k4_vr)/6)
    phi(i+1) = phi(i) + h*(k1_phi+2*k2_phi+2*k3_phi+k4_phi)/6
    vphi(i+1) = vphi(i) + h*(k1_vphi+2*k2_vphi+2*k3_vphi+k4_vphi)/6

   end do
   
   fn = r

  end function rk4

  subroutine orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)

   real(dp), intent(in) :: tf, r0, phi0, vphi0, vr0, h
   real(dp) :: k1_r, k1_phi, k1_vr, k1_vphi, k2_phi, k2_r, k2_vphi, k2_vr,k3_phi,k3_r,k3_vphi,k3_vr,k4_phi,k4_r,k4_vphi, &
               k4_vr 
   real(dp), allocatable, intent(out) :: t(:), r(:), phi(:), vr(:), vphi(:)
   real(dp), allocatable :: fn(:)
   integer(i8) :: i, steps

   steps =  int(tf/h)
   allocate(t(steps),r(steps),phi(steps),vr(steps),vphi(steps), fn(steps))
   t(1) = 0
   r(1) = r0
   phi(1) = phi0
   vr(1) = vr0
   vphi(1) = vphi0
    
   open(unit=10, file='orbit_data.txt', status='replace')  
   write(10,*) 't(s) phi(rad) r(U.A)'

   do i = 1, steps - 1
    
    if (i /= 1) then
       t(i) = t(i-1) + tf/steps 
    end if

    k1_r = r(i)
    k1_vr = vr(i)
    k1_phi = phi(i)
    k1_vphi = vphi(i)
    k2_r = vr(i)+h*k1_vr/2
    k2_vr = ar(r(i)+h*k1_r/2,vr(i)+h*k1_vr/2,vphi(i)+h*k1_vphi/2)
    k2_phi = vphi(i)+h*k1_vphi/2
    k2_vphi = aphi(r(i)+h*k1_r/2,vr(i)+h*k1_vr/2,vphi(i)+h*k1_vphi/2)
    k3_r = vr(i)+h*k2_vr/2
    k3_vr = ar(r(i)+h*k2_r/2,vr(i)+h*k2_vr/2,vphi(i)+h*k2_vphi/2)
    k3_phi = vphi(i)+h*k2_vphi/2
    k3_vphi = aphi(r(i)+h*k2_r/2,vr(i)+h*k2_vr/2,vphi(i)+h*k2_vphi/2)
    k4_r = vr(i)+h*k3_vr
    k4_vr = ar(r(i)+h*k3_r,vr(i)+h*k3_vr,vphi(i)+h*k3_vphi)
    k4_phi = vphi(i)+h*k3_vphi
    k4_vphi = aphi(r(i)+h*k3_r,vr(i)+h*k3_vr,vphi(i)+h*k3_vphi)

    r(i+1) = r(i) + h*(k1_r+2*k2_r+2*k3_r+k4_r)/6
    vr(i+1) = (vr(i) + h*(k1_vr+2*k2_vr+2*k3_vr+k4_vr)/6)
    phi(i+1) = phi(i) + h*(k1_phi+2*k2_phi+2*k3_phi+k4_phi)/6
    vphi(i+1) = vphi(i) + h*(k1_vphi+2*k2_vphi+2*k3_vphi+k4_vphi)/6
    
    write(10,*) t(i), phi(i), r(i)

   end do
   
   fn = r
   

   close(10)

  end subroutine orbit_data

end module calculo_orbitas
