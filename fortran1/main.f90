program main
    use ParamIteration

    !break main__
    character*100 buf
    !integer min, max, ind
    mini = 1 !28
    maxi = 1000
    ind = 0
    buf = 'mcmcNtimes'
    if (buf=='mcmc') then
        ind = 0
        mini = 0
        maxi = 0
        call  SynParall
    else if (buf=='mcmcNtimes') then
        !Execute SunParall for SZdata[min].txt prior[min].dat, ... ,SZdata[max].txt prior[max].dat
        !in SZdaras and priors directories consequentially
        !And create chainconsum[min].dat MCchain[min].dat, ... ,chainconsum[max].dat MCchain[max].dat
        !in OUT directory
        do ind = mini, maxi
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
        call SetFunction_SZ !(0)
	end if



end program
