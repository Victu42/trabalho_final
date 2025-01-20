program main

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
    use calculo_orbitas
    
    implicit none

    ! Parâmetros da simulação
    
    real(dp), parameter :: t0 = 0.0
    real(dp), parameter :: tf = 10.0     ! Tempo final
    real(dp), parameter :: r0 = 8.0e20_dp*Rs ! Raio inicial 

    real(dp), allocatable :: t(:), r(:), phi(:)
  real(dp), allocatable :: dr_dt(:), dphi_dt(:)
  
  n_steps = int((tf - t0) / h)
  
  allocate(t(0:n_steps), r(0:n_steps), phi(0:n_steps))
  allocate(dr_dt(0:n_steps), dphi_dt(0:n_steps))
  

  t(0) = t0
  r(0) = r0   
  phi(0) = 0.0
  
  dr_dt(0) = 0.0
  dphi_dt(0) = c*sqrt(G*M/(r(0)**3))

    
    ! Integração com verificação de validade
    do i = 0, n_steps-1
        call rk4_step(t(i), r(i), phi(i), dr_dt(i), dphi_dt(i), &
                     t(i+1), r(i+1), phi(i+1), dr_dt(i+1), dphi_dt(i+1), &
                     is_valid)
        
        if (.not. is_valid .or. r(i+1) <= Rs) then
            print *, "Simulação interrompida no passo:", i
            print *, "Tempo coordenado (s):", t(i)
            print *, "Raio (Rs):", r(i)/Rs
            n_steps = i
            exit
        end if
    end do
    
    call save_results(t, r, phi, dr_dt, dphi_dt, n_steps)
    
    deallocate(t, r, phi, dr_dt, dphi_dt)
    


    

  
end program main
