subroutine plot_exc_graph

    use CI_data
    use CO_data
    use MCMC
    use synth_param

    
    character*200 file_name, buf,mode
    
    external calc_fit_CI,calculate_chi_CI,calc_fit_CO,calculate_chi_CO,calc_fit_H2,calculate_chi_H2
    double precision chi, x, exc_data(7), n, tkin, A1(3),A2(3), x1(4),x2(4),x3(4)
    
    call SetFunction_CO
    call SetFunction_CI
    
    file_name = TRIM(ADJUSTL('data.sss'))
    call OpenAnalysis(file_name)

    mode = 'exc_co'
    if (mode == 'exc_co') then
        syn(1)%val(1) = 2.72
        syn(1)%val(2) = 1.75 !log10(57.0)
        syn(1)%val(3) = 2.725*(1+0) 
        syn(1)%val(4) = 15.64
        call calc_fit_CO
        call calculate_chi_CO(chi,1)
        open(122,file='excitation.dat')
        do i =1,num_J+1 
            write(122,'(2f20.10)') EJ_CO(i), syn(1)%val(4) + dlog10(NJ_CO(i))
        end do
        close(122)
        if (1==1) then    
            open(122,file='excitation_co_gas.dat')
            syn(1)%val(1) = 0.0
            syn(1)%val(2) = 2.0
            syn(1)%val(3) = 0. !2.725*(1+3.) !2.337) 
            syn(1)%val(4) = 0.0 
            syn(1)%val(5) = 0.0
            do i =1,101
                syn(1)%val(1) = 0 + 5.*(i-1)/100
                call calc_fit_CO        
                write(122,'(5f20.10)') 5.*(i-1)/100., dlog10(NJ_CO(2)),dlog10(NJ_CO(3)),dlog10(NJ_CO(4)),dlog10(NJ_CO(5))
            end do
            close(122)       
        end if
    end if
    if(mode == 'exc_ci') then 
        syn(1)%val(1) = 2.35!4.27-log10(57.)
        syn(1)%val(2) = 1.75
        syn(1)%val(3) = 2.73
        syn(1)%val(4) = 14.91
        syn(1)%val(5) = 0.0
        call calc_fit_CI
        call calculate_chi_CI(chi,1)
        open(122,file='excitation.dat')
        do i =1,3
            write(122,'(2f20.10)') EJ_CI(i), syn(1)%val(4) + dlog10(NJ_CI(i))
        end do
        close(122)

        if (1==0) then    
            open(122,file='excitation_t_cmb_n2_0_uv0_part.dat')
            syn(1)%val(1) = 2.0
            syn(1)%val(2) = 2.0
            syn(1)%val(3) = 2.725*(1+2.) !2.337) 
            syn(1)%val(4) = 0.0 
            syn(1)%val(5) = 0.0
            do i =1,101
                syn(1)%val(3) = 2.725*(1+4.*(i-1)/100.) 
                !syn(1)%val(1) = 0 + 3.*(i-1)/100
                call calc_fit_CI        
                write(122,'(3f20.10)') 4.*(i-1)/100., dlog10(NJ_CI(2)),dlog10(NJ_CI(3))
            end do
            close(122)       
        end if
    end if
    
    if(1 == 0) then 
        open (122, file='plot_exc_graph.dat')
        write(122,'(1a20,7i20)') 'EJ_data', (i-1, i=1,7)
        do j = 1,3
            do i =1,7
                syn(1)%val(1) = 0.+ 6*(i-1)/6. 
                syn(1)%val(2) = 1.1
                syn(1)%val(3) = 9.6
                syn(1)%val(5) = 3.
                call calc_fit_CI       
                exc_data(i) = dlog10(NJ_CI(j)) - log10(2*(j-1)+1.0)
            end do
             write(122,'(8f20.10)') EJ_CI(j), (exc_data(i), i =1,7)
        end do
        close(122)
    end if
    
    if (1== 0) then 
        call calc_fit_all
        call calculate_chi_CO(chi,1)
    end if
    
    if(1==0) then
        open(45, file='excitation_map.dat')
        do i= 1,101
            do j = 1,101
                n = 0. + 6*(i-1)/100.
                tkin = 0. + 3*(j-1)/100.
                call calc_texc_co(n,tkin)
                if (T_exc_co<0) T_exc_co = 0
                write(45, '(3f20.10)') tkin, n, T_exc_co
            end do
        end do
        close(45)        
    end if

    if(1==0) then
        open(45, file='h2_fraction.dat')
        syn(1)%val(2) = 2.0
        syn(1)%val(1) = 2.0
        syn(1)%val(3) = 2.725*(3.0)
        syn(1)%val(5) = 0.0
        do j =1,100
            syn(1)%val(1) = 0.0 + 6.0*(j-1)/99. 
            do i = 1,3
                h2_fraction = 0.0 + 0.5*(i-1)/2. 
                call calc_fit_CI
                A1(i) = NJ_CI(2)
                A2(i) = NJ_CI(3)
                
            end do
            write(45, '(7f20.10)')  syn(1)%val(1), (A1(i)-A1(1), i =1,3),(A1(i), i =1,3)
        end do
        close(45)        
    end if
    end 
    
subroutine  calc_best_fit_excitation 
use CI_data
    use CO_data
    use hd_data
    use MCMC
    use synth_param
    
    character*200 file_name, buf
    
    external calc_fit_CI,calculate_chi_CI,calc_fit_CO,calculate_chi_CO,calc_fit_H2,calculate_chi_H2
    double precision chi, x
    character*10 selector
    
    call SetFunction_CO
    call SetFunction_CI
    call SetFunction_HD 
    
    selector = 'HD'
    SELECT CASE (selector)
        CASE ('CI')
            file_name = TRIM(ADJUSTL('dataJ0013_CI.sss'))
            call OpenAnalysis(file_name)
            call calc_fit_CI
            x = NJ_CI_fit(1,1)
            print*, x + log10(NJ_CI(1)), x + log10(NJ_CI(2)), x+log10(NJ_CI(3))
            pause
        CASE ('HD')
            file_name = TRIM(ADJUSTL('dataJ0843_HD.sss'))
            call OpenAnalysis(file_name)
            call calc_fit_HD
            x = NJ_HD_fit(1,1)
            print*, x + log10(hd_exc_pop(1)), x + log10(hd_exc_pop(2)), x+log10(hd_exc_pop(3))
            pause

    END SELECT
   

end

subroutine calc_exc

    use CI_data
    use CO_data
    use MCMC
    use synth_param

    
    character*200 file_name, buf,mode
    integer str_num
    
    external calc_fit_CI,calculate_chi_CI,calc_fit_CO,calculate_chi_CO,calc_fit_H2,calculate_chi_H2
    double precision chi, x, exc_data(7), n, tkin, A1(3),A2(3), x1(4),x2(4),x3(4)
    double precision, allocatable :: tkin_pdr(:), nco(:), nh(:), nh2(:), nhe(:)
    
    call SetFunction_CO
    call SetFunction_CI
    
    !file_name = TRIM(ADJUSTL('model/tkin.txt'))
    open(17, file='model/tkin.txt')
    call READ_EOF(17, str_num)
    open(18, file='model/COj0.txt')
    allocate(tkin_pdr(str_num))
    allocate(nco(str_num))
    allocate(nh(str_num))
    allocate(nh2(str_num))
    allocate(nhe(str_num))
 
    do i =1,str_num
        read(17,*) tkin_pdr(i)
    end do
    close(17)
    do i =1,str_num
        read(18,*) nco(i)
    end do  
    close(18)
    open(18, file='model/H.txt')
    do i =1,str_num
        read(18,*) nh(i)
    end do  
    close(18)
    open(18, file='model/H2.txt')
    do i =1,str_num
    read(18,*) nh2(i)
    end do  
    close(18)
    open(18, file='model/He.txt')
    do i =1,str_num
    read(18,*) nhe(i)
    end do  
    close(18)
    
    
    
    open(122,file='excitation.dat')    
    mode = 'exc_co'
    if (mode == 'exc_co') then
        do i =1,str_num
            syn(1)%val(1) = 2.0
            syn(1)%val(2) = tkin_pdr(i)
            syn(1)%val(3) = 2.732*(1+0.0) 
            syn(1)%val(4) = nco(i)
            syn(1)%val(5) = -1.0
            external_data_flag =1
            npdr0 = log10(2.96) !nh(i)
            npdr1 = log10(48.5) !nh2(i)
            npdr2 = log10(10.0) !nhe(i)
            call calc_fit_CO
            call calc_fit_CI
            write(122,'(6f20.10)') (nco(i)+dlog10(NJ_CO(j)), j=1,3),(dlog10(NJ_CI(j)), j=1,3)
            
        end do
    end if
    close(122)
    
    if (1<0) then
        external_data_flag = 0
        syn(1)%val(1) = 2.23
        syn(1)%val(2) = 2.0
        syn(1)%val(3) = 2.725 
        syn(1)%val(4) = 0.0
        call calc_fit_CO
        print*, (dlog10(NJ_CO(j)), j=1,7)
    end if
end