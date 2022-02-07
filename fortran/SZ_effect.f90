module SZ_data
    save
    integer :: num_of_bands = 5
    integer num_of_obs
    integer band70_num
    integer band100_num
    integer band143_num
    integer band217_num
    integer band353_num
    double precision, allocatable :: band70_l(:), band70_f(:)        ! Wavelengths and fluxes of profiles for Planck bands
    double precision, allocatable :: band100_l(:), band100_f(:)
    double precision, allocatable :: band143_l(:), band143_f(:)
    double precision, allocatable :: band217_l(:), band217_f(:)
    double precision, allocatable :: band353_l(:), band353_f(:)
    double precision sz_wave(5), sz_flux(5)  ! Wavelength, fluxes in our model, flux in measurements
    double precision sz_wave_test(9), sz_flux_test(9)
    !double precision sz_wave_test(1000), sz_flux_test(1000)
    double precision sz_obs_flux(5,3)  ! Observed fluxes end errors in 5 bands
    double precision coeff1, coeff2                      ! coeffs h/k_b, k_b/mc^2

end

       
subroutine SetFunction_SZ
    
    use SZ_data
    use constants
    
    integer i,j,k_s,k
    character*500 buf, str 
    character(len=50), dimension(6) :: args
    double precision x
  
    
    !###################################    set parameters for SZ function 
    coeff1 = h/k_bol*1e9      ! in K/GHz
    coeff2 = 1/511.875          ! m_ec^2 in KeV
    sz_wave(1) = 70           ! Planck freq. in GHz
    sz_wave(2) = 100 
    sz_wave(3) = 143 
    sz_wave(4) = 217
    sz_wave(5) = 353  
    
    
    
    !###################################    Read Planck bands
    
    !################# 70
    open(70, file='data/Planck_bands/fits-70.txt')
    !# Это здесь определяется band70_num?
    call READ_EOF(70,  band70_num)
    allocate(band70_l(band70_num))
    allocate(band70_f(band70_num))
    
    do i = 1,band70_num
        !# Что за '(1A500)' ?
        read (70, '(1A500)') str  
        k_s = 1
        do j=1,2
            k = index(str(k_s:), ' ')
            args(j) = str(k_s:k_s+k-2)
            k_s = k_s + k
        enddo
        read(args(1),*) band70_l(i)
        read(args(2),*) band70_f(i)
    enddo
    close(70)
    
    !################# 100
    open(100, file='data/Planck_bands/fits-100.txt')
    !# Это здесь определяется band70_num?
    call READ_EOF(100,  band100_num)
    allocate(band100_l(band100_num))
    allocate(band100_f(band100_num))
    
    do i = 1,band100_num
        !# Что за '(1A500)' ?
        read (100, '(1A500)') str  
        k_s = 1
        do j=1,2
            k = index(str(k_s:), ' ')
            args(j) = str(k_s:k_s+k-2)
            k_s = k_s + k
        enddo
        read(args(1),*) band100_l(i)
        read(args(2),*) band100_f(i)
    enddo
    close(100)
    
    !################# 143
    open(143, file='data/Planck_bands/fits-143.txt')
    !# Это здесь определяется band70_num?
    call READ_EOF(143,  band143_num)
    allocate(band143_l(band143_num))
    allocate(band143_f(band143_num))
    
    do i = 1,band143_num
        !# Что за '(1A500)' ?
        read (143, '(1A500)') str  
        k_s = 1
        do j=1,2
            k = index(str(k_s:), ' ')
            args(j) = str(k_s:k_s+k-2)
            k_s = k_s + k
        enddo
        read(args(1),*) band143_l(i)
        read(args(2),*) band143_f(i)
    enddo
    close(143)
    
    !################# 217
    open(217, file='data/Planck_bands/fits-217.txt')
    !# Это здесь определяется band70_num?
    call READ_EOF(217,  band217_num)
    allocate(band217_l(band217_num))
    allocate(band217_f(band217_num))
    
    do i = 1,band217_num
        !# Что за '(1A500)' ?
        read (217, '(1A500)') str  
        k_s = 1
        do j=1,2
            k = index(str(k_s:), ' ')
            args(j) = str(k_s:k_s+k-2)
            k_s = k_s + k
        enddo
        read(args(1),*) band217_l(i)
        read(args(2),*) band217_f(i)
    enddo
    close(217)
    
    !################# 353
    open(353, file='data/Planck_bands/fits-353.txt')
    !# Это здесь определяется band70_num?
    call READ_EOF(353,  band353_num)
    allocate(band353_l(band353_num))
    allocate(band353_f(band353_num))
    
    do i = 1,band353_num
        !# Что за '(1A500)' ?
        read (353, '(1A500)') str  
        k_s = 1
        do j=1,2
            k = index(str(k_s:), ' ')
            args(j) = str(k_s:k_s+k-2)
            k_s = k_s + k
        enddo
        read(args(1),*) band353_l(i)
        read(args(2),*) band353_f(i)
    enddo
    close(353)
    
    
    !###################################    Read observations
    open(17, file='SZ_data.txt')
    read(17, *)   
    read(17,*) num_of_obs
    print*, 'read SZ data in', num_of_obs, 'bands' 
    do i = 1, num_of_obs
        read(17,'(3e10.2)') (sz_obs_flux(i,j), j = 1,3)  ! x, x+err, x-err
        print*, sz_obs_flux(i,1), sz_obs_flux(i,2),sz_obs_flux(i,3)
        !# end do или enddo ?
    end do  
    close(17)

    
    end subroutine 

    
    
subroutine calc_fit_SZ          ! procedure to calculate sz_signal
    
    use SZ_data
    use collision_data  
    use synth_param    
    
    external sz_signal
    
    integer i,j 
    double precision R, T0, Te, beta, tau, nu , heta           ! Model parameters
    double precision x, xx, theta
    double precision Int_s, s, ss
    double precision T0_test, Tau_test, theta_test, x_test, beta_test, s_test, xx_test, ss_test
    
    
    T0  = syn(1)%val(1) ! in K
    Te = syn(1)%val(2)  ! in Kev
    beta = syn(1)%val(3) ! vz/c
    tau = syn(1)%val(4) ! optical depth
    
    theta  = coeff2*Te
    
    !###################################    set sz_flux
    
    !# 70
    i = 1
    Int_s = 0.0
    x = coeff1*sz_wave(i)/T0
    sz_flux(i) = sz_signal(T0, tau, theta, x, beta)
    do j = 1, size(band70_f)-1
        x = coeff1*band70_l(j)/T0
        xx = coeff1*band70_l(j + 1)/T0
        s = sz_signal(T0, tau, theta, x, beta)
        ss = sz_signal(T0, tau, theta, xx, beta)
        Int_s = Int_s +  0.5 * (band70_f(j)*s + band70_f(j + 1)*ss )*(band70_l(j+1)-band70_l(j))
    end do
    sz_flux(i) = Int_s
    
    !# 100
    i = 2
    Int_s = 0.0
    x = coeff1*sz_wave(i)/T0
    sz_flux(i) = sz_signal(T0, tau, theta, x, beta)
    do j = 1, size(band100_f)-1
        x = coeff1*band100_l(j)/T0
        xx = coeff1*band100_l(j + 1)/T0
        s = sz_signal(T0, tau, theta, x, beta)
        ss = sz_signal(T0, tau, theta, xx, beta)
        Int_s = Int_s +   0.5 * (band100_f(j)*s + band100_f(j + 1)*ss )*(band100_l(j+1)-band100_l(j))
    end do
    sz_flux(i) = Int_s
    
    !# 143
    i = 3
    Int_s = 0.0
    x = coeff1*sz_wave(i)/T0
    sz_flux(i) = sz_signal(T0, tau, theta, x, beta)
    do j = 1, size(band143_f)-1
        x = coeff1*band143_l(j)/T0
        xx = coeff1*band143_l(j + 1)/T0
        s = sz_signal(T0, tau, theta, x, beta)
        ss = sz_signal(T0, tau, theta, xx, beta)
        Int_s = Int_s +   0.5 * (band143_f(j)*s + band143_f(j + 1)*ss )*(band143_l(j+1)-band143_l(j))
    end do
    sz_flux(i) = Int_s
    
    !# 217
    i = 4
    Int_s = 0.0
    x = coeff1*sz_wave(i)/T0
    sz_flux(i) = sz_signal(T0, tau, theta, x, beta)
    do j = 1, size(band217_f)-1
        x = coeff1*band217_l(j)/T0
        xx = coeff1*band217_l(j + 1)/T0
        s = sz_signal(T0, tau, theta, x, beta)
        ss = sz_signal(T0, tau, theta, xx, beta)
        Int_s = Int_s +   0.5 * (band217_f(j)*s + band217_f(j + 1)*ss )*(band217_l(j+1)-band217_l(j))
    end do
    sz_flux(i) = Int_s
    
    !# 353
    i = 5
    Int_s = 0.0
    x = coeff1*sz_wave(i)/T0
    sz_flux(i) = sz_signal(T0, tau, theta, x, beta)
    do j = 1, size(band353_f)-1
        x = coeff1*band353_l(j)/T0
        xx = coeff1*band353_l(j + 1)/T0
        s = sz_signal(T0, tau, theta, x, beta)
        ss = sz_signal(T0, tau, theta, xx, beta)
        Int_s = Int_s + 0.5 * (band353_f(j)*s + band353_f(j + 1)*ss )*(band353_l(j+1)-band353_l(j))
    end do
    sz_flux(i) = Int_s
    
    
    end subroutine    

double precision function sz_signal(T0, tau, theta, x, beta) !################################### Changeble
    
    double precision T0, alpha, z,  theta, x, beta, tau, koeff
    double precision X0, S, Y0, Y1, Y2, Y3, Y4, C1, C2, P0, P1, R
    
    z = 0.0518
    !################################### Changeble
    
    X0 = x * (dexp(x) + 1.0) / (dexp(x) - 1.0)
    S = 2.0 * x / (dexp(- x / 2.0) * (dexp(x) - 1))
    Y0 = x * (dexp(x)+1.0)/(dexp(x)-1.0) - 4.0
    Y1 = - 10.0 + 47.0 / 2.0 * X0 - 42.0 / 5.0 * X0 ** 2 + 7.0 / 10.0 * X0 ** 3 + S ** 2 * (- 21.0 / 5.0 + 7.0 / 5.0 * X0)
    Y2 = - 15.0 / 2.0 + 1023.0 / 8.0 * X0 - 868.0 / 5.0 * X0 ** 2 + 329.0 / 5.0 * X0 ** 3 - 44.0 / 5.0 * X0 ** 4 + 11.0 / 30.0 * X0 ** 5 + S ** 2 * (- 434.0 / 5.0 + 658.0 / 5.0 * X0 - 242.0 / 5.0 * X0 ** 2 + 143.0 / 30.0 * X0 ** 3) + S ** 4 * ( - 44.0 / 5.0 + 187.0 / 60.0 * X0)
    Y3 = 15.0 / 2.0 + 2505.0 / 8.0 * X0 - 7098.0 / 5.0 * X0 ** 2 + 14253.0 / 10.0 * X0 ** 3 - 18594.0 / 35.0 * X0 ** 4 + 12059.0 / 140.0 * X0 ** 5 - 128.0 / 21.0 * X0 ** 6 + 16.0 / 105.0 * X0 ** 7 + S ** 2 * (- 7098.0 / 10.0 + 14253.0 / 5.0 * X0 - 102267.0 / 35.0 * X0 ** 2 + 156767.0 / 140.0 * X0 ** 3 - 1216.0 / 7.0 * X0 ** 4 + 64.0 / 7.0 * X0 ** 5) + S ** 4 * (- 18594.0 / 35.0 + 205003.0 / 280.0 * X0 - 1920.0 / 7.0 * X0 ** 2 + 1024.0 / 35.0 * X0 ** 3) + S ** 6 * (- 544.0 / 21.0 + 992.0 / 105.0 * X0)
    !Y4 = - 135.0 / 32.0 + 30375.0 / 128.0 * X0 - 62391.0 / 10.0 * X0 ** 2 + 614727.0 / 40.0 * X0 ** 3 - 124389.0 / 10.0 * X0 ** 4 + 355703.0 / 80.0 * X0 ** 5 - 16568.0 / 21.0 * X0 ** 6 + 7516.0 / 105.0 * X0 ** 7 - 22.0 / 7.0 * X0 ** 8 + 11.0 / 210.0 * X0 ** 9 + S ** 2 * (- 62391.0 / 20.0 + 614727.0 / 20.0 * X0 - 1368279.0 / 20.0 * X0 ** 2 + 4624139.0 / 80.0 * X0 ** 3 - 157396.0 / 7.0 * X0 ** 4 + 30064.0 / 7.0 * X0 ** 5 - 2717.0 / 7.0 * X0 ** 6 + 2761.0 / 210.0 * X0 ** 7) + S ** 4 * (- 124389.0 / 10.0 + 6046951.0 / 160.0 * X0 - 248520.0 / 7.0 * X0 ** 2 + 481024.0 / 35.0 * X0 ** 3 - 15972.0 / 7.0 * X0 ** 4 + 18689.0 / 140.0 * X0 ** 5) + S ** 6 * (- 70414.0 / 21.0 + 465992.0 / 105.0 * X0 - 11792.0 / 7.0 * X0 ** 2 + 19778.0 / 105.0 * X0 ** 3) + S ** 8 * (- 682.0 / 7.0 + 7601.0 / 210.0 * X0)
    C1 = 10.0 - 47.0 / 5.0 * X0 + 7.0 / 5.0 * X0 ** 2 + 7.0 / 10.0 * S ** 2
    C2 = 25.0 - 1117.0 / 10.0 * X0 + 847.0 / 10.0 * X0 ** 2 - 183.0 / 10.0 * X0 ** 3 + 11.0 / 10.0 * X0 ** 4 + S ** 2 * (847.0 / 20.0 - 183.0 / 5.0 * X0 + 121.0 / 20.0 * X0 ** 2) + 11.0 / 10.0 * S ** 4
    !P0 = - 2.0 / 3.0 + 11.0 / 30.0 * X0
    !P1 = - 4.0 + 12.0 * X0 - 6.0 * X0 ** 2 + 19.0 / 30.0 * X0 ** 3 + S ** 2 * (- 3.0 + 19.0 / 15.0 * X0)
    !R = theta ** 2 * Y1 + theta ** 3 * Y2 + theta ** 4 * Y3 + theta ** 5 * Y4 + beta ** 2 * (1.0 / 3.0 * Y0 + theta * (5.0 / 6.0 * Y0 + 2.0 / 3.0 * Y1)) - beta * (1.0 + theta * C1 + theta ** 2 * C2) + beta ** 2 * (P0 + theta * P1)
    R = theta ** 2 * Y1 + theta ** 3 * Y2 + theta ** 4 * Y3 - beta * (1.0 + theta * C1 + theta ** 2 * C2)
    
    koeff = 1.0e+006                    ! signal in MicroK ??

    sz_signal = koeff*T0*tau*(theta*Y0 + R) 

    
end function
    
subroutine calculate_chi_SZ(chi,f)
    
    use SZ_data
    use synth_param
    
    integer f
    double precision chi,x
    
    syn(1)%free_parameters = 0
    chi = 0
    if(num_of_obs >1) then 
        do i = 1, num_of_obs
            x = sz_obs_flux(i,1) -  sz_flux(i)
            if (x .ge. 0) then 
                chi = chi + x*x/(sz_obs_flux(i,2))**2
            else
                chi = chi + x*x/(sz_obs_flux(i,3))**2
            end if
        end do    
    end if
    if(syn(1)%chi > 1) chi = syn(1)%chi
    syn(1)%free_parameters = num_of_obs -  num_of_var
    
    
    if ((f == 2) .and. (syn(1)%free_parameters > 0)) then
        chi = chi/syn(1)%free_parameters
    endif
    
    end subroutine