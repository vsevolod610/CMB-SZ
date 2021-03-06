program main
    use ParamIteration
    
    character *100 buf
    !integer min, max, ind
    min = 1!28
    max = 77
    ind = 0 
    buf = 'mcmc'
    if (buf=='mcmc') then
        call  SynParall   !(0, 0, 0)
    else if (buf=='mcmcNtimes') then
        !Execute SunParall for SZdata[min].txt prior[min].dat, ... ,SZdata[max].txt prior[max].dat 
        !in SZdaras and priors directories consequentially
        !And create chainconsum[min].dat MCchain[min].dat, ... ,chainconsum[max].dat MCchain[max].dat
        !in OUT directory
        do ind = min, max
            call  SynParall  !(ind, min, max)
        end do
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
        call SetFunction_SZ!(0)
	end if
    
    
    
end program