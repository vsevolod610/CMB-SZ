module collision_data

    !>>>>>>>>>>>>>>>>>   CO collisions
    double precision qCOH(106,29)
    double precision qCOH2o(106,25)
    double precision qCOH2p(106,25)
    double precision qCOHe(106,12)

    !>>>>>>>>>>>>>>>>    CI collisions
    double precision qC1H(4,22)
    double precision qC1H2o(106,25)
    double precision qC1H2p(106,25)
    double precision qC1He(106,12)

    double precision qCIH(9,3,1)
    double precision qCIH2o(6,3,3)
    double precision qCIH2p(6,3,3)
    double precision qCIHe(6,3,1)

    !>>>>>>>>>>>>>>>  HD collision
    ! numbers are: 1 - h1, 2 - he, 3 - h2o, 4 - h2p   for transitions between J = 0..8
    double precision a_coll(4,9,9), b_coll(4,9,9), c_coll(4,9,9)
    double precision dt(4)
end module


module CO_data
    save
    integer :: num_J = 14       ! calculation will performed for #J = 0 to #J = 14 CO rotational levels
    integer :: num_of_level_CO = 5 ! The number of CO rotational levels used for calc_fit
    integer :: calc_with_CI = 0   ! constant to use additioanl parameter - logNJ=0 for CO only
    integer :: external_data_flag = 0
    double precision AJ_CO(15)  ! Einstein coeefs
    double precision EJ_CO(15)  ! energy of rotational level J
    double precision NJ_CO(15)  ! ratio of CO N#J(CO)/NJ=0(CO)
    double precision NJ_CO_fit(15,3)  ! data and errors
    double precision T_exc_co
    double precision npdr0,npdr1,npdr2
    double precision :: beta_factor = 1.0 ! CO photon trapping into the cloud


end module

module CI_data
    save
    integer :: num_of_level_CI = 3 ! The number of CI levels used for calc_fit
    integer :: aG_constrain = 0 ! add alpha_G constrain on n-UV
    double precision A10_CI     ! The spontaneous transition probability for CI
    double precision A20_CI     ! The spontaneous transition probability for CI
    double precision A21_CI     ! The spontaneous transition probability for CI
    double precision EJ_CI(3)   ! The energy of rotational level J
    double precision NJ_CI(4)   ! The ratio of column densities of CI for J=0,1,2 to NCI(J=0); NJ_CI(4) = NJ_CI(3)/NJ_CI(2)
    double precision NJ_CI_fit(15,3)  ! data and errors
    DOUBLE PRECISION ACI_pump(54,3)  !  The spontaneous transition probability for higher level of CI to calculate UV pumping, see Silva$Viegas2000
    DOUBLE PRECISION ECI_pump(54)  !  The energies of higher level of CI to calculate UV pumping, see Silva$Viegas2000
    integer gCI_pump(54)  !  The energies of higher level of CI to calculate UV pumping, see Silva$Viegas2000
    DOUBLE PRECISION Gamma_pump(3,3)
    DOUBLE PRECISION , allocatable :: EBLflux(:,:)
    DOUBLE PRECISION p_aG
    double precision :: h2_fraction = 1.0 !f_H2
    end module

module H2_data
    save
    integer :: num_of_level_H2 = 2
    double precision NJ_H2(2)
    double precision NJ_H2_fit(2,3)
    double precision :: E10_H2 = 170.5 ! Excitation energy of the first rotat level [in K]
    double precision lnL_H2(100,100), logn_LH2(100), loguv_LH2(100)
    end module

module hd_data
    save
    integer :: num_of_level_hd = 2      ! number of HD roational levels into consideration
    ! ground state
    double precision hd_Aij(15,13,15,13) !Einshtein coeffs of transitions nu,J->nu',J' witn nu'<=nu and J'-J= {+2, +1, +0, -1, -2}
                                         !nu = 0..14, J=0..10 (but matrix has larger size)
    double precision hdAtot(15,13)       !total spontaneous emission probability or Invert lifetife of ground energy levels (nu,J)
    double precision hdpop(15,13)        !pop ratio of ground energy levels
    !double precision hd_q_1(15,13)      !cascade entry rates Q(nu_0,J_0) nu_0 = 0..14, ,J_0 =0..10
    double precision hd_q(15,13,11)      !cascade entry rates Q(nu_0,J_0) nu_0 = 0..14, ,J_0 =0..10 from individual X-state J level
    double precision hd_EX(15,13)        !energy of rovibrational level of X state
    double precision anuJ_J(14,11,11)    !cascade efficiency factors
    double precision hd_gamma(11,11)     !cascade excitation rates Gamma(i|j) = Summ_nu=1,14_J=0,10 hd_q(nu,J|j)*anuJ_J(nu,J,i) + hd_q(0,i|j)
                                         !stratified for each j level of initial_X_state, which was populated
    double precision hd_exc_pop(11)      !the final population of HD rotational levels
    double precision NJ_HD_fit(2,3)


    end module




module constants
    save
    double precision, parameter :: e = 4.8029*1e-010
    double precision, parameter :: me = 9.108*1e-028
    double precision, parameter :: Pi = 3.141592654
    double precision, parameter :: c = 2.99792458*1e+010
    double precision, parameter :: h = 6.6254*1e-027
    double precision, parameter :: k_bol = 1.3805*1e-016
    double precision, parameter :: mu = 1820

end module

module ParamIteration
    save
    integer :: ind = 0
    integer :: maxi = 0
    integer :: mini = 0
    integer :: numOfIter = 1
end module


