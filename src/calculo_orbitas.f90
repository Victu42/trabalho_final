module calculo_orbitas

use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
implicit none

real(dp), parameter :: pi = 16*ATAN(1./5.) - 4*ATAN(1./239.)
real(dp), parameter :: c_m_s = 299792458 
real(dp), parameter :: c_UA = 63197.7909261 
real(dp) :: M, tf, r0, phi0, vr0, vphi0, h, vr0_luz, vphi0_luz
real(dp), parameter :: G = 6.670430E-11


contains
  
  function ar_luz(r,dphi) result(f) 
    
    real(dp) :: r, dphi, f

    f = ((-4*M**2)+(2*M*r)+(r-5*M)*(r**3)*(dphi**2))/r**3

  end function ar_luz

  function ar(r, dr, dphi) result(f) 
    
    real(dp) :: r, dphi, f, dr

    f = ((4*(M**3))-(4*r*(M**2))-(4*(M**2)*(r**3)*(dphi**2))+(4*M*(r**4)*(dphi**2))-((r**5) &
    *dphi**2)+(r**2)*(M-3*M*(dr**2)))/((2*M-r)*(r**3))

  end function ar

  function aphi(r,dr,dphi) result(f)
     
    real(dp) :: r, dr, dphi, f
  
    f = 2*(-3*M+r)*dr*dphi/((2*M-r)*r)

  end function aphi

  subroutine orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)

   real(dp), intent(in) :: tf, r0, phi0, vphi0, vr0, h
   real(dp) :: k1_r, k1_phi, k1_vr, k1_vphi, k2_phi, k2_r, k2_vphi, k2_vr,k3_phi,k3_r,k3_vphi,k3_vr,k4_phi,k4_r,k4_vphi, &
               k4_vr
   real(dp), allocatable, intent(out) :: t(:), r(:), phi(:), vr(:), vphi(:)
   integer(i8) :: i, steps

   steps =  int(tf/h)
   allocate(t(steps),r(steps),phi(steps),vr(steps),vphi(steps))

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

    k1_r = h*vr(i)
    k1_vr = h*ar(r(i),vr(i),vphi(i))
    k1_phi = h*vphi(i)
    k1_vphi = h*aphi(r(i), vr(i), vphi(i))
    k2_r = h*vr(i)+h*k1_vr/2
    k2_vr = h*ar(r(i)+k1_r/2,vr(i)+k1_vr/2,vphi(i)+k1_vphi/2)
    k2_phi = h*vphi(i)+h*k1_vphi/2
    k2_vphi = h*aphi(r(i)+k1_r/2,vr(i)+k1_vr/2,vphi(i)+k1_vphi/2)
    k3_r = h*vr(i)+h*k2_vr/2
    k3_vr = h*ar(r(i)+k2_r/2,vr(i)+k2_vr/2,vphi(i)+k2_vphi/2)
    k3_phi = h*vphi(i)+h*k2_vphi/2
    k3_vphi = h*aphi(r(i)+k2_r/2,vr(i)+k2_vr/2,vphi(i)+k2_vphi/2)
    k4_r = h*vr(i)+h*k3_vr
    k4_vr = h*ar(r(i)+k3_r,vr(i)+k3_vr,vphi(i)+k3_vphi)
    k4_phi = h*vphi(i)+h*k3_vphi
    k4_vphi = h*aphi(r(i)+k3_r,vr(i)+k3_vr,vphi(i)+k3_vphi)

    r(i+1) = r(i) + (k1_r+2*k2_r+2*k3_r+k4_r)/6
    vr(i+1) = vr(i) + (k1_vr+2*k2_vr+2*k3_vr+k4_vr)/6
    phi(i+1) = phi(i) + (k1_phi+2*k2_phi+2*k3_phi+k4_phi)/6
    vphi(i+1) = vphi(i) + (k1_vphi+2*k2_vphi+2*k3_vphi+k4_vphi)/6
    
    write(10,*) t(i), phi(i), r(i)

   end do

   close(10)

  deallocate(t,r,phi,vr,vphi)

  end subroutine orbit_data



subroutine ver_orbitas(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)

   real(dp), intent(in) :: tf, r0, phi0, vphi0, vr0, h
   real(dp) :: k1_r, k1_phi, k1_vr, k1_vphi, k2_phi, k2_r, k2_vphi, k2_vr,k3_phi,k3_r,k3_vphi,k3_vr,k4_phi,k4_r,k4_vphi, &
               k4_vr
   real(dp), allocatable, intent(out) :: t(:), r(:), phi(:), vr(:), vphi(:)
   integer(i8) :: i, steps

   steps =  int(tf/h)
   allocate(t(steps),r(steps),phi(steps),vr(steps),vphi(steps))

   t(1) = 0
   r(1) = r0
   phi(1) = phi0
   vr(1) = vr0
   vphi(1) = vphi0
    
   open(unit=11, file='eye_data.txt', status='replace')  
   write(11,*) 't(s) phi(rad) r(U.A)'

   do i = 1, steps - 1
    
    if (i /= 1) then
       t(i) = t(i-1) + tf/steps 
    end if

    k1_r = h*vr(i)
    k1_vr = h*ar_luz(r(i), vphi(i))
    k1_phi = h*vphi(i)
    k1_vphi = h*aphi(r(i), vr(i), vphi(i))
    k2_r = h*vr(i)+h*k1_vr/2
    k2_vr = h*ar_luz(r(i)+k1_r/2, vphi(i)+k1_vphi/2)
    k2_phi = h*vphi(i)+h*k1_vphi/2
    k2_vphi = h*aphi(r(i)+k1_r/2,vr(i)+k1_vr/2,vphi(i)+k1_vphi/2)
    k3_r = h*vr(i)+h*k2_vr/2
    k3_vr = h*ar_luz(r(i)+k2_r/2, vphi(i)+k2_vphi/2)
    k3_phi = h*vphi(i)+h*k2_vphi/2
    k3_vphi = h*aphi(r(i)+k2_r/2,vr(i)+k2_vr/2,vphi(i)+k2_vphi/2)
    k4_r = h*vr(i)+h*k3_vr
    k4_vr = h*ar_luz(r(i)+k3_r, vphi(i)+k3_vphi)
    k4_phi = h*vphi(i)+h*k3_vphi
    k4_vphi = h*aphi(r(i)+k3_r,vr(i)+k3_vr,vphi(i)+k3_vphi)

    r(i+1) = r(i) + (k1_r+2*k2_r+2*k3_r+k4_r)/6
    vr(i+1) = vr(i) + (k1_vr+2*k2_vr+2*k3_vr+k4_vr)/6
    phi(i+1) = phi(i) + (k1_phi+2*k2_phi+2*k3_phi+k4_phi)/6
    vphi(i+1) = vphi(i) + (k1_vphi+2*k2_vphi+2*k3_vphi+k4_vphi)/6
    
    write(11,*) t(i), phi(i), r(i)

   end do

   close(11)

  deallocate(t,r,phi,vr,vphi)

  end subroutine ver_orbitas

end module calculo_orbitas
