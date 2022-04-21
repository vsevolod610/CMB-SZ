!NJ_H2(1) = (data_H2_NJ(1,1) - 5*data_H2_NJ(1,2))  +  10*data_H2_NJ(1,2)/(j_max -1)*(j_H2-1)             ! varied in range pm2sigma
!NJ_H2(2) = NJ_H2(1) + dlog10(9.0d0) - dlog10(2.71828d0)*170.5/T_kin 
!chiq_H2 =  ((data_H2_NJ(1,1) - NJ_H2(1))/data_H2_NJ(1,2))**2. + ((data_H2_NJ(2,1) - NJ_H2(2))/data_H2_NJ(2,2))**2.! calculate chi2
    
    
subroutine calc_fit_H2
    
    use H2_data
    use synth_param
    
    double precision x
    double precision chi_UV, a_he, f_H, N0, T_cmb, T_kin, dens
    
    T_kin = 10**syn(1)%val(2)
    !dens = 10**syn(1)%val(1)/T_kin
    !T_cmb = syn(1)%val(3)
    !N0 = syn(1)%val(4)
    !f_H = 0.9
    !a_he = 0.085
    !chi_UV = 0.0


    NJ_H2(:) = 0.0
    NJ_H2(1) = 1.0                                                                        !>>>> NJ=0/NJ=0
    NJ_H2(2) = 9.0*dexp(-E10_H2/T_kin)                                                     !>>>> NJ=1/NJ=0
    
    end

    
subroutine calculate_chi_H2(chi,f)
    
    use H2_data
    use synth_param
    
    integer f
    double precision chi,x
    
    syn(1)%free_parameters = 0
    chi = 0
    if(num_of_level_H2 >1) then 
        do i = 1, num_of_level_H2
            x = NJ_H2_fit(i,1) - syn(1)%val(4) -dlog10(NJ_H2(i))
            if (x .ge. 0) then 
                chi = chi + x*x/(NJ_H2_fit(i,2))**2                 ! err++
            else
                chi = chi + x*x/(NJ_H2_fit(i,3))**2                 ! err--
            end if
        end do    
    end if
    if(syn(1)%chi > 1) chi = syn(1)%chi
    syn(1)%free_parameters = num_of_level_H2 -  num_of_var
    
    
    if ((f == 2) .and. (syn(1)%free_parameters > 0)) then
        chi = chi/syn(1)%free_parameters
    endif
end subroutine
