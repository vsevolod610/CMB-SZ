program main
    character *100 buf
    buf = 'mcmc'
    if (buf=='mcmc') then 
        call  SynParall
    else if (buf=='plot_exc') then
        call plot_exc_graph
    else if (buf=='calc_exc') then
        call calc_exc
    else if (buf=='calc_lnL') then
        call calc_lnL_plane
    else if (buf=='calc_bf') then
        call calc_best_fit_excitation 
    else if (buf == 'calcgrid_ci') then
        call calcgrid_ci
    else if (buf == 'calcgrid_co') then
        call calcgrid_co
    else if (buf == 'test') then 
        call SetFunction_SZ
    end if
    
end program