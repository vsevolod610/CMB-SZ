double precision function flux_draine(x) ! x- lambda
    ! Calculate flux density in [photons/cm^2 s Hz/sr] according to Draine'1978
    ! x - wavelength in Angstrom
    integer i,j
    double precision x,f
    f = 0
    
    f = 8.49e-05/x - 0.13666/x/x + 54.482/x/x/x
    
    flux_draine = f
        
    return
    end


double precision function flux_ebl(x,z) ! x- lambda
    ! Calculate flux density in [photons/cm^2 s Hz/sr] according to Khaire & Srianand 2018;
    ! x - wavelength in Angstrom
    use CI_data
    integer i,j
    double precision x,f,z
    f = 0

    
    i = 1
    do while ((i.le.size(EBLflux(:,1))) .AND. (EBLflux(i,1)<x))
        i = i + 1
    end do
    i = i - 1  
    
    f = EBLflux(i,2) + (EBLflux(i+1,2)-EBLflux(i,2))*(x-EBLflux(i,1))/(EBLflux(i+1,1)-EBLflux(i,1)) 
        
    flux_ebl = f
        
    return
    end

    