module calculo_orbitas

use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
implicit none
  

  real(dp), parameter :: c = 2.99792458e8
  real(dp), parameter :: G = 6.67430e-11
  real(dp), parameter :: M = 2.0e36
  real(dp), parameter :: Rs = 2.0*G*M/c**2
  real(dp), parameter :: h = 0.1       ! Passo de tempo
  integer :: n_steps, i
  logical :: is_valid
    
    

  contains

    subroutine rk4_step(t0, r0, phi0, dr0, dphi0, &
                        t1, r1, phi1, dr1, dphi1, is_valid)
        real(dp), intent(in) :: t0, r0, phi0, dr0, dphi0
        real(dp), intent(out) :: t1, r1, phi1, dr1, dphi1
        logical, intent(out) :: is_valid
        real(dp) :: k1(4), k2(4), k3(4), k4(4)
        real(dp) :: r_temp, phi_temp, dr_temp, dphi_temp
        
        is_valid = .true.
        
        if (r0 <= Rs) then
            is_valid = .false.
            return
        end if
        
        ! k1
        call derivatives(r0, phi0, dr0, dphi0, k1, is_valid)
        if (.not. is_valid) return
        
        ! k2
        r_temp = r0 + h*k1(1)/2
        if (r_temp <= Rs) then
            is_valid = .false.
            return
        end if
        phi_temp = phi0 + h*k1(2)/2
        dr_temp = dr0 + h*k1(3)/2
        dphi_temp = dphi0 + h*k1(4)/2
        call derivatives(r_temp, phi_temp, dr_temp, dphi_temp, k2, is_valid)
        if (.not. is_valid) return
        
        ! k3
        r_temp = r0 + h*k2(1)/2
        if (r_temp <= Rs) then
            is_valid = .false.
            return
        end if
        phi_temp = phi0 + h*k2(2)/2
        dr_temp = dr0 + h*k2(3)/2
        dphi_temp = dphi0 + h*k2(4)/2
        call derivatives(r_temp, phi_temp, dr_temp, dphi_temp, k3, is_valid)
        if (.not. is_valid) return
        
        ! k4
        r_temp = r0 + h*k3(1)
        if (r_temp <= Rs) then
            is_valid = .false.
            return
        end if
        phi_temp = phi0 + h*k3(2)
        dr_temp = dr0 + h*k3(3)
        dphi_temp = dphi0 + h*k3(4)
        call derivatives(r_temp, phi_temp, dr_temp, dphi_temp, k4, is_valid)
        if (.not. is_valid) return
        
        ! Atualiza variáveis com verificação
        t1 = t0 + h
        r1 = r0 + h*(k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6
        
        if (r1 <= Rs .or. isnan(r1)) then
            is_valid = .false.
            return
        end if
        
        phi1 = phi0 + h*(k1(2) + 2*k2(2) + 2*k3(2) + k4(2))/6
        dr1 = dr0 + h*(k1(3) + 2*k2(3) + 2*k3(3) + k4(3))/6
        dphi1 = dphi0 + h*(k1(4) + 2*k2(4) + 2*k3(4) + k4(4))/6
        
        if (isnan(r1) .or. isnan(phi1) .or. isnan(dr1) .or. isnan(dphi1)) then
            is_valid = .false.
        end if
    end subroutine rk4_step
    
    subroutine derivatives(r, phi, dr, dphi, dxdt, is_valid)
        real(dp), intent(in) :: r, phi, dr, dphi
        real(dp), intent(out) :: dxdt(4)
        logical, intent(out) :: is_valid
        real(dp) :: f, Gamma_r_rr, Gamma_r_phiphi, Gamma_phi_rphi
        
        is_valid = .true.
        
        if (r <= Rs) then
            is_valid = .false.
            return
        end if
        
        f = 1.0_dp - 2.0_dp*G*M/(r*c**2)
        if (abs(f) < 1.0e-10) then
            is_valid = .false.
            return
        end if
        
        ! Símbolos de Christoffel relevantes
        Gamma_r_rr = -(G*M)/(r**2*c**2*f)
        Gamma_r_phiphi = -r*f
        Gamma_phi_rphi = 1.0/r
        
        ! Derivadas em relação ao tempo coordenado
        dxdt(1) = dr
        dxdt(2) = dphi
        dxdt(3) = -Gamma_r_rr*dr**2 - Gamma_r_phiphi*dphi**2
        dxdt(4) = -2.0*Gamma_phi_rphi*dr*dphi
        
        if (any(isnan(dxdt))) then
            is_valid = .false.
        end if
    end subroutine derivatives
    
    subroutine save_results(t, r, phi, dr_dt, dphi_dt, n_steps)
        integer, intent(in) :: n_steps
        real(dp), intent(in) :: t(0:n_steps), r(0:n_steps), phi(0:n_steps)
        real(dp), intent(in) :: dr_dt(0:n_steps), dphi_dt(0:n_steps)
        integer :: i
        
        open(unit=10, file='orbit_data.txt', status='replace')
        write(10,*) '# t(s)  r(m)  phi(rad)  dr/dt(m/s)  dphi/dt(rad/s)'
        do i = 0, n_steps
            write(10,*) t(i), r(i), phi(i), dr_dt(i), dphi_dt(i)
        end do
        close(10)
    end subroutine save_results


end module calculo_orbitas
