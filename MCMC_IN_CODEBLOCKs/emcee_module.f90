


module synth_param
    save
    type synth_parameters
        integer number_of_elements
        character*8 name(10)
        double precision val(10), val_min(10), val_max(10), val_step(10)
        double precision chi
        integer :: var_ind(10), vary(10)
        integer :: free_parameters
    end type synth_parameters

    type(synth_parameters) :: syn(1), syn2(1)

    type prior_type
        integer*4 ind
        integer mode
        double precision c, p, m
        double precision up, down
        character*100 namelnLH2
    end type
    type(prior_type), allocatable :: prior(:)
    integer*4 num_prior

    double precision, allocatable :: random(:,:)
    double precision, allocatable :: confi(:,:)
!    double precision AutoCorr       ! autocorrelation in MCMC fit
    integer*4 :: corr_lenght = 1    ! correlation lenght
    integer*4 :: k_clen = 1         ! number of walker with highest correlation lenght
    integer*4 :: num_of_var = 0
    integer number_of_elements
end module



module MCMC

    use synth_param
    save
    integer num_MC   ! size of the MC set
    integer iter_MC  ! number of iterations
    integer iter     ! current iteration
    integer add_prior     ! use prior for var's distributions
    double precision, allocatable :: GRvalue(:,:) ! value of Gelman Rubin test for each parameter at certain step
    integer min_MC   ! the minium chi number in the MC set: equal to 1 at the begining
    integer :: affine = 1 ! affine status
    double precision :: aff = 2.0     ! width of distribution in affine method case

    type MC_type
        type(synth_parameters) :: syn(1)
        double precision chi
    end type MC_type

    type(MC_type), allocatable :: MC(:)
    type(MC_type), allocatable :: MCsave(:)
    type(MC_type), allocatable :: MCmin(:)

!    integer :: FitMC_step_status = 1
    double precision chi_min, chi_min_glob, MCDisp_min

    type MCMC_type
        double precision, allocatable :: par(:,:) ! first is the number of point, second is the number of parameter
        double precision, allocatable :: chi(:)   ! chi^2 in chain
        integer, allocatable :: auto(:)           ! give how many steps from previous shift. Need to calculate Autocorrelation lenght
        double precision mom(2)
        double precision, allocatable :: autocorr(:,:) ! autocorrelation array
        double precision, allocatable :: med(:,:)   ! array contains median of x_{ik} for each parameter (k) of the chain (i)
        double precision, allocatable :: sca(:,:)   ! array contains x_{ik}^2 for each parameter (k) of the chain (i)
        double precision, allocatable :: ave(:)     ! array contains median for each parameter (k) along all chains at some step
        double precision, allocatable :: disp(:)    ! array dispersion for each parameter (k) along all chains at some step
    end type

    type(MCMC_type), dimension (:), allocatable :: MCchain

    double precision MC_moments(2)
    double precision, allocatable :: Autocorr(:,:)

    !double precision, allocatable :: par_min(:)     ! parameter of minium through entire MC run
    double precision chi_emin                       ! minium of chi^2 through entire MC run

    double precision, allocatable :: MCres(:)

!   >>>>>>> to tune the program
    integer :: calc_moments = 1
    integer :: calc_MCmin = 0
    integer :: calc_Autocorr = 1
    integer :: show_addMC = 0
    integer :: calc_stats = 1
    integer :: write_data = 1
    integer :: write_data_final = 1
    integer :: calc_GRtest = 1
    integer :: stats_chain = 1
    integer :: adjust_MC = 0
    integer :: res_or_buf = 1
    integer :: write_MCmin = 0


    end module

