program main

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, i4 => int32, i8 => int64
    use trabalho_final.f90
    
    implicit none

    ! Parâmetros da simulação
    real(dp), parameter :: M = 2.0e36
    real(dp), parameter :: t0 = 0.0
    real(dp), parameter :: tf = 10.0     ! Tempo final
    real(dp), parameter :: h = 0.1       ! Passo de tempo
    real(dp), parameter :: r0 = 8.0e20_dp*Rs ! Raio inicial 

    
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
