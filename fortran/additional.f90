subroutine calc_fit_all
    
    call calc_fit_CO 
    call calc_fit_CI
    
    end

subroutine calculate_chi_all(chi,f)
    
    use CO_data
    use CI_data
    use collision_data  
    use synth_param
    
    integer f
    double precision chi,x
    
    syn(1)%free_parameters = 0
    chi = 0
    
    if(num_of_level_CI >1) then 
        do i = 1, num_of_level_CI
            x = NJ_CI_fit(i,1) - syn(1)%val(4) -dlog10(NJ_CI(i))
            if (x .ge. 0) then 
                chi = chi + x*x/(NJ_CI_fit(i,2))**2                 ! err++
            else
                chi = chi + x*x/(NJ_CI_fit(i,3))**2                 ! err--
            end if
        end do    
    end if
    
    syn(1)%free_parameters = syn(1)%free_parameters + num_of_level_CI -  num_of_var
    
    
    
    if(num_of_level_CO >1) then 
        do i = 1, num_of_level_CO
            x = NJ_CO_fit(i,1) - syn(1)%val(6) -dlog10(NJ_CO(i))
            if (x .le. 0) then 
                chi = chi + x*x/(NJ_CO_fit(i,2))**2
            else
                chi = chi + x*x/(NJ_CO_fit(i,3))**2
            end if
        end do    
    end if
    if(syn(1)%chi > 1) chi = syn(1)%chi
    syn(1)%free_parameters = syn(1)%free_parameters + num_of_level_CO -  num_of_var
    
    
    if ((f == 2) .and. (syn(1)%free_parameters > 0)) then
        chi = chi/syn(1)%free_parameters
    endif

    end subroutine
    
subroutine calc_lnL_plane


    use CI_data
    use CO_data
    use MCMC
    use synth_param

    
    character*200 file_name, buf
    
    external calc_fit_CI,calculate_chi_CI
    double precision chi_uv, x, exc_data(7), n, tkin,chi,res(100,100),cmax 
    
    call SetFunction_CI
    
    file_name = TRIM(ADJUSTL('data.sss'))
    call OpenAnalysis(file_name)
    
!    open(11,file='lnL.dat')
    tcmb = 2.725*3.626
    if(1 == 1) then 
        t_kin = 100.
        do i =1,100
            print*, i
            n = 0.+5.0*(i-1)/99
            syn(1)%val(1) = n
            syn(1)%val(2) = log10(t_kin)
            syn(1)%val(3) = tcmb   ! tcmb
            syn(1)%val(5) = -10.
            call calc_fit_CI
            cmax = NJ_CI(3)
            do j=1,100
                chi_uv = -1.+4.0*(j-1)/99
                
                
                syn(1)%val(1) = n
                syn(1)%val(2) = log10(t_kin)
                syn(1)%val(3) = tcmb   ! tcmb
                syn(1)%val(5) = chi_uv                
                
                call calc_fit_CI
                call calculate_chi_CI(chi,1)
!                write(11,'(3f20.10)') n, chi_uv, -chi
                res(i,j) = NJ_CI(3)/cmax
            end do    
        end do
    end if
    
!    close(11)
    
    
    open(123,file='2dgrid_ci_20.value')
    do i = 1,100
        write(123,'(100f20.10)')  (res(i,j), j =1,100)
    end do
    close(123)
    end
    
    
subroutine calcgrid_ci


    use CI_data
    use MCMC
    use synth_param

    
    character*200 file_name, buf
    
    external calc_fit_CI,calculate_chi_CI
    double precision chi_uv, x, exc_data(7), n, tkin,chi,res(200,200),cmax , logn, loguv, logt, logN0,logp,tcmb
    double precision grid3d(200,200,200)
    
    call SetFunction_CI
    
    file_name = TRIM(ADJUSTL('data.sss'))
    call OpenAnalysis(file_name)
    
    mode = '2d'
    num = 100
    if (mode == '2d') then 
        logt = syn(1)%val(2)
        do i =1,num
            logn = 0.0+4.0*(i-1)/(num-1)
            do j =1,num
                loguv = -1 + 4.0*(j-1)/(num-1)
                !tcmb = 2. + 10.0*(j-1)/(num-1)
                print*, i,j
                syn(1)%val(1) = logn
                !syn(1)%val(2) = log10(t_kin)
                !syn(1)%val(3) = tcmb   ! tcmb
                syn(1)%val(5) = loguv
                call calc_fit_CI
                call calculate_chi_CI(chi,1)
                res(i,j) = chi
            end do    
        end do
    end if

    
    
    
    open(123,file='grid_ci_fit.2dgrid')
    do i = 1,num
        write(123,'(<num>f8.1)')  (res(j,i), j =1, num)
    end do
    close(123)
        

    end


    subroutine calcgrid_co


    use CO_data
    use MCMC
    use synth_param

    
    character*200 file_name, buf
    
    external calc_fit_CO,calculate_chi_CO
    double precision chi_uv, x, exc_data(7), n, tkin,chi,res(200,200),cmax , logn, loguv, logt,tcmb, logN0,logp
    
    call SetFunction_CO
    
    file_name = TRIM(ADJUSTL('data.sss'))
    call OpenAnalysis(file_name)
      
    mode = 1
    if (mode == 1) then
        num = 100        
        logt = syn(1)%val(2)
        do i =1,num
            logp = 3.+2.0*(i-1)/(num-1)
            do j =1,num
                logN0 = 13 + 1.0*(j-1)/(num-1)
                print*, i+j
                syn(1)%val(1) = logp - logt
                syn(1)%val(4) = logN0 
                call calc_fit_CO
                call calculate_chi_CO(chi,1)
                res(i,j) = chi
            end do    
        end do
    else if (mode == 2) then 
        t_kin = 140.
        num = 101
        syn(1)%val(4) = 13.66

        do i =1,num
            logn = 0.+4.0*(i-1)/(num-1)
            !syn(1)%val(1) = logn
            !syn(1)%val(2) = log10(t_kin)
            !syn(1)%val(3) = tcmb   ! tcmb
            !syn(1)%val(5) = loguv
            !call calc_fit_CI
            !call calculate_chi_CI(chi,1)
            print*, i,logn
            do j =1,num
                tcmb = 5+8.0*(j-1)/(num-1)
                syn(1)%val(1) = logn
                syn(1)%val(2) = log10(t_kin)
                syn(1)%val(3) = tcmb   ! tcmb
                syn(1)%val(5) = 0 
                call calc_fit_CO
                call calculate_chi_CO(chi,1)
                res(i,j) = chi
            end do    
        end do
    end if
    
    
    open(123,file='grid_co_fit.grid')
    do i = 1,num
        write(123,'(120f20.10)')  (res(j,i), j =1,num)
    end do
    close(123)
    end
