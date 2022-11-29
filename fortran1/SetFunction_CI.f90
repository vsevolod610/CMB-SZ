subroutine SetFunction_CI
    
    use CI_data
    use collision_data
    
    integer i,j, n, s
    character*500 buf 
    character*10 buf2 
    double precision x
    double precision, allocatable::  M(:,:)
  
    
    !###################################    set parameters for CI
    !>>>>>>>>> Einstein coeffs
    A10_CI = 7.932*1e-008
    A20_CI = 2.054*1e-014 
    A21_CI = 2.654*1e-007 
    
    ! Energy levels in Kelvin
    EJ_CI(1)  = 0.0            
    EJ_CI(2)  = 23.6
    EJ_CI(3)  = 62.4 
    
    ! set collisional approximation coeffs
    ! C + He
    ! 1-> 0
    qCIHe(1,1,1) = 0.2723
    qCIHe(2,1,1) = 0.07442
    qCIHe(3,1,1) = -1.4e-03
    qCIHe(4,1,1) = 1.1E-5
    qCIHe(5,1,1) = -3e-8
    qCIHe(6,1,1) = 0.0 !5.15122E-10 
    ! 2->0
    qCIHe(1,2,1) = 3.7569
    qCIHe(2,2,1) = 0.03798
    qCIHe(3,2,1) = -9.42072E-4
    qCIHe(4,2,1) = 1.1346E-5
    qCIHe(5,2,1) = -6.47276E-8
    qCIHe(6,2,1) = 1.41676E-10 
    ! 2->1
    qCIHe(1,3,1) = 6.64226
    qCIHe(2,3,1) = 0.06457                   
    qCIHe(3,3,1) = -0.0015
    qCIHe(4,3,1) = 1.8447E-5
    qCIHe(5,3,1) = -1.05051E-7
    qCIHe(6,3,1) = 2.26428E-10
    
    ! C + H
    ! 1->0
    qCIH(1,1,1) = 3.6593
    qCIH(2,1,1) = 56.6023
    qCIH(3,1,1) = -802.9765
    qCIH(4,1,1) = 5025.1882
    qCIH(5,1,1) = -17874.4255
    qCIH(6,1,1) = 38343.6655
    qCIH(7,1,1) = -49249.4895
    qCIH(8,1,1) = 34789.3941
    qCIH(9,1,1) = -10390.9809
    ! 2->0
    qCIH(1,2,1) = 10.8377
    qCIH(2,2,1) = -173.4153
    qCIH(3,2,1) = 2024.0272
    qCIH(4,2,1) = -13391.6549
    qCIH(5,2,1) = 52198.5522
    qCIH(6,2,1) = -124518.3586
    qCIH(7,2,1) = 178182.5823
    qCIH(8,2,1) = -140970.6106
    qCIH(9,2,1) = 47504.5861
    ! 3->0
    qCIH(1,3,1) = 15.8996
    qCIH(2,3,1) = -201.3030
    qCIH(3,3,1) = 1533.6164
    qCIH(4,3,1) = -6491.0083
    qCIH(5,3,1) = 15921.9239
    qCIH(6,3,1) = -22691.1632
    qCIH(7,3,1) = 17334.7529
    qCIH(8,3,1) = -5517.9360
    qCIH(9,3,1) = 0.0
    
    !C + H2p
    ! 1->0
    ! T<60K
    qCIH2p(1,1,1) =  0.396         
    qCIH2p(2,1,1) = 0.10261
    qCIH2p(3,1,1) = -0.00613
    qCIH2p(4,1,1) = 1.60875E-4
    qCIH2p(5,1,1) = -2.02125E-6
    qCIH2p(6,1,1) = 9.9E-9
    ! 60K < T < 500 K
    qCIH2p(1,1,2) =  0.9772      
    qCIH2p(2,1,2) = -0.00667
    qCIH2p(3,1,2) =  5.04333E-5
    qCIH2p(4,1,2) = -1.72901E-7
    qCIH2p(5,1,2) =  2.83504E-10
    qCIH2p(6,1,2) = -1.79137E-13
    ! T > 500K
    qCIH2p(1,1,3) = 0.572      
    qCIH2p(2,1,3) = 5.11e-04
    qCIH2p(3,1,3) = -2.68e-07
    qCIH2p(4,1,3) = 0.0
    qCIH2p(5,1,3) = 0.0
    qCIH2p(6,1,3) = 0.0
    !2->1 
    ! T < 90 
    qCIH2p(1,3,1) = 1.45917            
    qCIH2p(2,3,1) = 0.04422
    qCIH2p(3,3,1) = -0.00182
    qCIH2p(4,3,1) = 3.35388E-5
    qCIH2p(5,3,1) = -2.93561E-7
    qCIH2p(6,3,1) = 9.9359E-10
    ! 90< T < 500
    qCIH2p(1,3,2) = 1.89138        
    qCIH2p(2,3,2) = -0.00403
    qCIH2p(3,3,2) = 3.57579E-5
    qCIH2p(4,3,2) = -1.07352E-7
    qCIH2p(5,3,2) = 1.54964E-10
    qCIH2p(6,3,2) = -8.83167E-14
    ! 500 < T
    qCIH2p(1,3,3) = -0.09303     
    qCIH2p(2,3,3) = 0.01213
    qCIH2p(3,3,3) = -2.56903E-5
    qCIH2p(4,3,3) = 3.02826E-8
    qCIH2p(5,3,3) = -1.85023E-11
    qCIH2p(6,3,3) = 4.48718E-15
    !------------------
    !2->0
    ! T < 90 K
    qCIH2p(1,2,1) = 0.8175        
    qCIH2p(2,2,1) = 0.02289
    qCIH2p(3,2,1) = -0.00105
    qCIH2p(4,2,1) = 2.01238E-5
    qCIH2p(5,2,1) = -1.79633E-7
    qCIH2p(6,2,1) = 6.08974E-10    
    !  T > 90K
    qCIH2p(1,2,2) = 1.04296  
    qCIH2p(2,2,2) = -0.00359
    qCIH2p(3,2,2) = 2.3142E-5
    qCIH2p(4,2,2) = -6.85846E-8
    qCIH2p(5,2,2) = 1.04731E-10
    qCIH2p(6,2,2) = -6.54297E-14   
    !  T > 500K
    qCIH2p(1,2,3) = 0.702      
    qCIH2p(2,2,3) = 7.06e-04
    qCIH2p(3,2,3) = -3.63E-7
    qCIH2p(4,2,3) = 0.0
    qCIH2p(5,2,3) = 0.0
    qCIH2p(6,2,3) = 0.0
    
    !C + H2o
    ! 1 -> 0
    ! 10 < T < 40
    qCIH2o(1,1,1) = 0.56             
    qCIH2o(2,1,1) = 0.0215
    qCIH2o(3,1,1) = -4.5e-4
    qCIH2o(4,1,1) = 0.0
    qCIH2o(5,1,1) = 0.0
    qCIH2o(6,1,1) = 0.0
    ! 40 <T <500
    qCIH2o(1,1,2) = 0.864
    qCIH2o(2,1,2) = -2.845E-3
    qCIH2o(3,1,2) = 1.654E-5
    qCIH2o(4,1,2) = -3.6377E-8
    qCIH2o(5,1,2) = 2.83629E-11
    qCIH2o(6,1,2) = 0.0
    ! T > 500
    qCIH2o(1,1,3) = 0.6075             !if T>500K
    qCIH2o(2,1,3) = 5.26e-04
    qCIH2o(3,1,3) = -2.7e-07
    qCIH2o(4,1,3) = 0.0
    qCIH2o(5,1,3) = 0.0
    qCIH2o(6,1,3) = 0.0
    !2->0
    ! T < 500
    qCIH2o(1,2,1) = 0.51009                
    qCIH2o(2,2,1) = 0.00319
    qCIH2o(3,2,1) = -2.07399E-5
    qCIH2o(4,2,1) = 8.7424E-8
    qCIH2o(5,2,1) = -1.71715E-10
    qCIH2o(6,2,1) = 1.24263E-13
       
    qCIH2o(1,2,2) = 0
    qCIH2o(2,2,2) = 0
    qCIH2o(3,2,2) = 0
    qCIH2o(4,2,2) = 0
    qCIH2o(5,2,2) = 0
    qCIH2o(6,2,2) = 0
    ! T > 500
    qCIH2o(1,2,3) = 1.90295            
    qCIH2o(2,2,3) = -0.00694
    qCIH2o(3,2,3) = 1.86878E-5
    qCIH2o(4,2,3) = -2.25626E-8
    qCIH2o(5,2,3) = 1.29225E-11
    qCIH2o(6,2,3) = -2.88462E-15
    !2->1
    ! T < 140
    qCIH2o(1,3,1) = 0.85818             
    qCIH2o(2,3,1) = 0.00975
    qCIH2o(3,3,1) = -5.00909E-5
    qCIH2o(4,3,1) = 1.8492E-7
    qCIH2o(5,3,1) = -3.42072E-10
    qCIH2o(6,3,1) = 2.40243E-13
    ! 140 <T < 500
    qCIH2o(1,3,2) = 0
    qCIH2o(2,3,2) = 0
    qCIH2o(3,3,2) = 0
    qCIH2o(4,3,2) = 0
    qCIH2o(5,3,2) = 0
    qCIH2o(6,3,2) = 0
    ! 500 < T
    qCIH2o(1,3,3) = 3.575
    qCIH2o(2,3,3) = -0.01153
    qCIH2o(3,3,3) = 3.50393E-5
    qCIH2o(4,3,3) = -4.47072E-8
    qCIH2o(5,3,3) = 2.65443E-11
    qCIH2o(6,3,3) = -6.08974E-15
    
    
    open(121, file = 'data/CI.dat')
    read(121, '(A)') buf  
    f = 0
    do while (f<6)
            read (121, '(1A100)') buf 
            if(TRIM(ADJUSTL(buf)) == '#')  f = f + 1
    enddo
    ACI_pump(:,:) = 0.D0
    ECI_pump(:) =0.d0
    read(121, *) n
    allocate(M(n,5))
    M(:,:) = 0.d0
    s = 1
    do i = 1, n
        read (121, '(1A100)') buf 
        read(buf(1:8),'(1e8.2)') x
        if(i == 1) then 
            ECI_pump(s) = x
        else if (x .ne. ECI_pump(s)) then 
            s = s +1
            ECI_pump(s) = x
        end if
        read(buf(12:13),'(1i1)') j
        read(buf(14:21),'(1e8.2)') ACI_pump(s,j)
        read(buf(10:11),'(1i1)') gCI_pump(s)
        
    end do
    deallocate(M)
    close(121)
    
    open(121, file = 'data/KS18_Fiducial_Q18/EBL_KS18_Q18_z_2.5.txt')
    read(121, '(A)') buf  
    do while (buf(1:1) == '#')
            read (121, '(1A100)') buf 
    enddo
    allocate(EBLflux(688,2) )
    do i = 1, 688
        read (121, '(1A100)') buf 
        read(buf(2:13),'(1e12.6)') EBLflux(i,1)             ! in Ang. 
        read(buf(19:30),'(1e12.6)') EBLflux(i,2)            ! in erg/s/cm^2/Hz/Sr
    end do
    close(121)
    ! convert flux to photon/s/cm^2/Hz/Sr
    do i = 1, size(EBLflux(:,1))
        EBLflux(i,2) = EBLflux(i,2)/(0.1984e-07/EBLflux(i,1)) !*1/E_nu = 1/(hc/lambda[Ang.])
    end do
    
end subroutine 

    
subroutine calc_fit_CI
    
    use CI_data
    use CO_data
    use collision_data  
    use synth_param
    !use physcond
    
    integer print_test, indx(2)
    double precision x
    double precision k01,k10,k21,k12,k20,k02                           ! coefficient of excitation rates (CMB+collisions) 
    double precision n_para, n_ortho,n1,n2,n3,dens, koeff              ! variable for H, H2, He densities
    double precision B01,B10,B20,B02,B12,B21                           ! coefficient of CMB_excitation
    double precision q01(4), q02(4), q12(4), q10(4), q20(4), q21(4)         ! coefficient of collision_excitation
    ! 1- para H2, 2 - ortho H2, 3 - Helium,  4 - Hydrogen 
    double precision chi_UV, a_he, f_H, T_cmb, T_kin, MC(3,3),A(2,2), B(2),d
    
    T_kin = 10**syn(1)%val(2)
    dens = 10**syn(1)%val(1) !/T_kin
    T_cmb = syn(1)%val(3)
    f_H = h2_fraction
    a_he = 0.085
    koeff = 1*0.75  ! convert Mathis to Draine units 
    chi_UV = 10**(syn(1)%val(5))*koeff
 
    !test
    !print_test = 1
    !open(123, file = 'test_coll_ci.dat')
    !do l = 1,1000
    !    T_kin = 10
    !    dens= 10**(0.+6*(l-1)/999.)
    !    chi_UV = 0.0
    !    T_cmb = 0.0
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>######## collisional excitation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    ! 1- para H2, 2 - ortho H2, 3 - Helium,  4 - Hydrogen 
    ! calculate coll_coeff using approximation of 
    ! 1 - Abrahamsson et al 2007
    ! 2 - K. Schroder, et al. 1991: J. Phys. B 24, 2487	
    ! 3 - V. Staemmler and D. R. Flower, 1991: J. Phys. B 24, 2343
    ! 4 - J. M. Launay and E. Roueff, 1977: A&A 56, 289
    do i = 1,4
         q01(i) = 0
         q02(i) = 0
         q12(i) = 0
         q10(i) = 0
         q20(i) = 0
         q21(i) = 0
    end do

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>######## first element - para_H2   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>> 1->0
    if (T_kin < 10.) then
        q10(1) = 0.96             ! set as a low limit
    else if (T_kin > 1200.) then 
        q10(1) = 0.80             ! set as a high limit
    else
        if(T_kin .le. 30.0) then 
            j = 1
        else if ((T_kin > 30.0) .and. (T_kin .le. 500.0)) then 
            j = 2
        else if ( T_kin > 500.0) then
            j = 3
        end if
        do i = 1,6      
            x =  T_kin**(i-1)                          
            q10(1) = q10(1) + qCIH2p(i,1,j)*x	
        end do
    end if
    !>>>> 0->1 
    q01(1) = q10(1)*3./dexp(23.6/T_kin) 
    
    !###> 2->0
    if (T_kin < 10.) then
        q20(1) = 0.96              ! set as a low limit
    else if (T_kin > 1200.) then 
        q20(1) = 1.03              ! set as a high limit
    else
         j = 1
        if (( T_kin > 90.0 ) .and. (T_kin < 500.0)) then 
            j = 2
        else if ( T_kin .ge. 500.0) then
            j = 3
        end if
        do i = 1,6                                  
            q20(1) = q20(1) + qCIH2p(i,2,j)*T_kin**(i-1)		
        end do
    end if
    !>>>>0->2
    q02(1) = q20(1)*5./dexp(62.4/T_kin) 

    !###> 2->1
    if (T_kin < 10.) then
        q21(1) = 1.75              ! set as a low limit
    else if (T_kin > 1200.) then 
        q21(1) = 2.60              ! set as a high limit
    else
        j = 1
        if ( (T_kin > 90.0) .and. (T_kin< 500.0)) then 
            j = 2
        else if ( T_kin .ge. 500.0 ) then
            j=3
        end if
        do i = 1,6                                  
            q21(1) = q21(1) + qCIH2p(i,3,j)*T_kin**(i-1)		
        end do
    end if
    !>>>> 1->2
    q12(1) = q21(1)*5./3./dexp(38.8/T_kin) 
        
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>######## second element - ortho_H2  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   
    !###> 1->0
    if (T_kin < 10.) then
        q10(2) = 0.73              ! set as a low limit
    else if (T_kin > 1200.) then 
        q10(2) = 0.85              ! set as a high limit
    else
        j = 1
        if ((T_kin > 30.0) .and. (T_kin .le. 500.0)) then 
            j = 2
        else if (T_kin > 500.0) then
            j=3
        end if
        do i = 1,6                                  
            q10(2) = q10(2) + qCIH2o(i,1,j)*T_kin**(i-1)		
        end do
    end if
    !>>>> 0->1
    q01(2) = q10(2)*3./dexp(23.6/T_kin) 

    !### 2->0
    if (T_kin < 10.) then
        q20(2) = 0.52              ! set as a low limit
    else if (T_kin > 1200.) then 
        q20(2) = 1.11              ! set as a high limit
    else
        j = 1
        if (T_kin > 500.0) j=3
        do i = 1,6                                  
            q20(2) = q20(2) + qCIH2o(i,2,j)*T_kin**(i-1)		
        end do
    end if
    !>>>> 0->2
    q02(2) = q20(2)*5./dexp(62.4/T_kin) 

    !###> 2->1
    if (T_kin < 10.) then
        q21(2) = 0.92              ! set as a low limit
    else if (T_kin > 1200.) then 
        q21(2) = 2.83              ! set as a high limit
    else
        j = 1
        if (T_kin > 500)   j=3
        do i = 1,6                                  
            q21(2) = q21(2) + qCIH2o(i,3,j)*T_kin**(i-1)		
        end do
    end if
    !>>>> 1->2
    q12(2) = q21(2)*5./3./dexp(38.8/T_kin) 
    
    do i = 1,2
        q10(i) = q10(i)*1e-010
        q20(i) = q20(i)*1e-010
        q01(i) = q01(i)*1e-010
        q02(i) = q02(i)*1e-010
        q12(i) = q12(i)*1e-010
        q21(i) = q21(i)*1e-010
    end do
    
   !######## third element - Helium
    if (T_kin .le. 10.) then
        q10(3) = 0.849              ! set as a low limit
        q20(3) = 4.05
        q21(3) = 7.15
    else if (T_kin > 150.) then 
        q10(3) = 1.86              ! set as a high limit
        q20(3) = 4.53
        q21(3) = 8.83
    else if((T_kin > 10).and.(T_kin  .le. 150)) then
        do i = 1,6                                  
            q10(3) = q10(3) + qCIHe(i,1,1)*T_kin**(i-1)		
            q20(3) = q20(3) + qCIHe(i,2,1)*T_kin**(i-1)		
            q21(3) = q21(3) + qCIHe(i,3,1)*T_kin**(i-1)		
        end do
    !else if((T_kin > 150).and.(T_kin  < 500)) then 
    !    q10(3) = 1.2459 + 3.0142e-03*T_kin
    !    q20(3) = 3.97561 + 6.09686e-03*T_kin
    !    q21(3) = 7.15433 + 1.64454e-02*T_kin
    !else
    end if
    
    q01(3) = q10(3)*3./dexp(23.6/T_kin) 
    q02(3) = q20(3)*5./dexp(62.4/T_kin)
    q12(3) = q21(3)*5./3./dexp(38.8/T_kin) 
    
    
    !######## fourth element - Hydrogen 
    if (T_kin .le. 10.) then
        q01(4) = 0.122              ! set as a low limit
        q02(4) = 1.0e-04
        q12(4) = 2.0e-03
    else if (T_kin > 1000.) then 
        q01(4) = 64.1              ! set as a high limit
        q02(4) = 88.6
        q12(4) = 90.6
    else
        do i = 1,9                 
            x = real(i)
            q01(4) = q01(4) + qCIH(i,1,1)*T_kin**(-(x-1)*0.25)		! POWER!!!!! int??
        end do
        q01(4) = dexp(q01(4))
        do i = 1,9                                  
            x = real(i)
            q02(4) = q02(4) + qCIH(i,2,1)*T_kin**(-(x-1)*0.333)		
        end do
        q02(4) = dexp(q02(4))
        do i = 1,9                              
            x = real(i)
            q12(4) = q12(4) + qCIH(i,3,1)*T_kin**(-(x-1)*0.25)		
        end do
        q12(4)= dexp(q12(4))
    end if
    q10(4) = q01(4)/3.*dexp(23.6/T_kin) 
    q20(4) = q02(4)/5.*dexp(62.4/T_kin)
    q21(4) = q12(4)/5.*3.*dexp(38.8/T_kin)
 
    do i = 3,4 
        q10(i) = q10(i)*1e-011
        q20(i) = q20(i)*1e-011
        q01(i) = q01(i)*1e-011
        q02(i) = q02(i)*1e-011
        q12(i) = q12(i)*1e-011
        q21(i) = q21(i)*1e-011
    end do
    
  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>######## CMB excitation>>>>>>>>>>>  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    B01 = 2.38e-007/(dexp(23.6/T_cmb) - 1.0)
    B10 = B01/3. 
    B02 = 1.02e-013/(dexp(62.4/T_cmb) - 1.0)
    B20 = B02/5. 
    B12 = 4.39e-007/(dexp(38.8/T_cmb) - 1.0)
    B21 = B12*3./5. 
    
    
    !########################>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! calculate UV excitaiton
    do i =1,3
        do j =1,3
            x=0
            call calc_pumping(i,j,x,chi_UV,T_cmb/2.725 - 1.)
            Gamma_pump(i,j) = x
        end do
    end do
    
    !do i =1,11
    !    t = (-1.0 + 2.0*(i-1)/10)
    !    call calc_pumping(1,2,x,10**(t*1.0d0),T_cmb/2.725 - 1.)
    !    print*, t, 10**(t*1.0d0), x
    !end do
   

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>######## solution matrix equation>  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    ! dens = n1 + n2 + n3
    !n1 = (1-f_H)/(1+a_He-0.5*f_H)*dens ! Hydrogen density
    !n2 = 0.5*f_H/(1+a_He-0.5*f_H)*dens ! H2 dedsity
    !n3 = a_He/(1+a_He-0.5*f_H)*dens  ! Helium density
    
    ! dens = n1 + 2n2
if (external_data_flag .eq. 0) then 
        n1 = (1-f_H)*dens ! Hydrogen density
        n2 = 0.5*f_H*dens ! H2 dedsity
        n3 = a_He*dens  ! Helium density
    else if (external_data_flag .eq. 1) then
        n1 = 10**npdr0
        n2 = 10**npdr1
        n3 = 10**npdr2
    end if
    n_para  = n2/(1.0+9*dexp(-170.6/T_kin)) !H2-para J0 density
    n_ortho = n2*9*dexp(-170.6/T_kin)/(1+9*dexp(-170.6/T_kin)) !H2-ortho J1 density
    k10 = 0.0
    k01 = 0.0 
    k10 = 0.0 
    k02 = 0.0 
    k20 = 0.0 
    k12 = 0.0 
    k21 = 0.0 
    
    if(1>0) then 
        k01 = B01 
        k10 = A10_CI + B10 
        k02 = B02 
        k20 = A20_CI + B20 
        k12 = B12 
        k21 = A21_CI + B21 
    end if
    if(1>0) then 
        k01 = k01 + (q01(1)*n_para+q01(2)*n_ortho+q01(3)*n3 + q01(4)*n1) 
        k10 = k10 + (q10(1)*n_para+q10(2)*n_ortho+q10(3)*n3 + q10(4)*n1) 
        k02 = k02 + (q02(1)*n_para+q02(2)*n_ortho + q02(3)*n3 + q02(4)*n1) 
        k20 = k20 + (q20(1)*n_para+q20(2)*n_ortho + q20(3)*n3 + q20(4)*n1) 
        k12 = k12 + (q12(1)*n_para+q12(2)*n_ortho + q12(3)*n3 + q12(4)*n1)
        k21 = k21 + (q21(1)*n_para+q21(2)*n_ortho + q21(3)*n3 + q21(4)*n1) 
    end if
    if(1>0) then 
        k01 = k01  + Gamma_pump(1,2) !
        k10 = k10 + Gamma_pump(2,1)
        k02 = k02  + Gamma_pump(1,3)
        k20 = k20 +  Gamma_pump(3,1)
        k12 = k12 +  Gamma_pump(2,3)
        k21 = k21 +  Gamma_pump(3,2)
    end if
   
    
    !>>>>>>>>>>>>>>>>>>>>>... matrix equation
  
    
    !MC(1,1) = Gamma_pump(1,1)*0-(k01 + k02)
    !MC(1,2) = k10
    !MC(1,3) = k20
    !MC(2,1) = k01
    !MC(2,2) = Gamma_pump(2,2)*0-(k10 + k12)
    !MC(2,3) = k21
    !MC(3,1) = k02
    !MC(3,2) = k12
    !MC(3,3) = Gamma_pump(3,3)*0-(k21 + k20)
    
     
    !A(:,:) = 0.0
    !B(:) = 0.0
    !
    !do i =1, 2
    !    do j = 1,2
    !        A(i,j) = MC(i,j+1) 
    !    end do
    !    B(i) = - MC(i,1)
    !end do
    
    
    !##### find solution
    !call ludcmp(A,2,2,indx,d)   ! LU decomposition
    !call lubksb(A,2,2,indx,B)   ! solution
    
    
    !>>>>>>>>>>>>>>>>>>>>
    !k01 =  Gamma_pump(1,2)
    !k02 =  Gamma_pump(1,3)
    
    NJ_CI(:) = 0.0
    NJ_CI(1) = 1.0                                                                        !>>>> NJ=0/NJ=0
    NJ_CI(2) = (k01*(k20+k21) + k02*k21)/(k12*k20+k10*(k20+k21))                          !>>>> NJ=1/NJ=0
    NJ_CI(3) = (k10*k02+(k01+k02)*k12)/(k10*k21+k20*(k10+k12))                            !>>>> NJ=2/NJ=0
    NJ_CI(4) = (k12*(k01+k02) + k02*k10)/((k20+k21)*(k01 + k02)-k02*k20)                  !>>>> NJ=2/NJ=1
    
      
    ! test
    !write(123, '(6e20.10)') T_kin, dens, k01 + k21*NJ_CI(3), k10 + k12, k02 + k12*NJ_CI(2), k20 + k21    
    !enddo
    
    end

    
subroutine calculate_chi_CI(chi,f)
    
    use CI_data
    use collision_data  
    use synth_param
    
    integer f
    double precision chi,x,y
    
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
    
    if(syn(1)%chi > 1) chi = syn(1)%chi
    syn(1)%free_parameters = num_of_level_CI -  num_of_var
    
    
    if ((f == 2) .and. (syn(1)%free_parameters > 0)) then
        chi = chi/syn(1)%free_parameters
    endif
    end subroutine
    
    subroutine calc_pumping(u,l,pump,galactic,z)
    use CI_data
    use constants
    
    integer u,l
    double precision flux_draine, KCI(54,3), ECI(3),x, uij, pump_coeff,s,lambda,gCI(3),pump,B(54,3), flux,galactic,z, ebl_coef
    external flux_draine, flux_ebl
    
    ECI(1) = 0.0d0 ! in cm^-1
    ECI(2) = 16.4
    ECI(3) = 43.4
    gCI(1) = 1
    gCI(2) = 3
    gCI(3) = 5
    
    
    ebl_coef =0.0
    KCI(:,:) = 0.d0 
    B(:,:)=0.0
    do i = 1, 54
        do j =1,3
            x = ECI_pump(i) - ECI(j) ! in cm-1
            lambda = 1e8/x 
            if (ebl_coef > 0) then 
                flux = flux_draine(lambda)*galactic + flux_ebl(lambda,z)*ebl_coef 
            else 
                flux = flux_draine(lambda)*galactic
            end if
            if(ECI_pump(i) > 0) then 
                KCI(i,j) = ACI_pump(i,j)/2/x**2.*flux
                B(i,j) = ACI_pump(i,j)/x**3/8/Pi/h
            end if
        end do
    end do
    
    pump_coeff = 0.0d0
    do k = 1, 54
        if(ACI_pump(k,u)>0) then 
            s = 0
            do i =1,3
                s = s + (ACI_pump(k,i)+ KCI(k,i))
            end do
            zz = KCI(k,u)*gCI_pump(k)/gCI(u)   
            if(s>0) pump_coeff = pump_coeff + zz*(ACI_pump(k,l)+ KCI(k,l))/s
        end if
    end do
    pump = pump_coeff
    end
    