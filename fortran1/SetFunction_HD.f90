!module HD_I_data
!
!double precision x1(2)  ! ratio HD1/HD0
!double precision HD(2,2)
!end module 
!
!module collision_data
!
!double precision q01(1),q10(1) ! collision rate of (H2 + HD) from n to m level of HD 
!double precision A(6,6,2)  ! coefficientd of 5-polynom approximation (with H2)
!double precision q_He_01,q_He_10 ! collision rate of (He + HD) from n to m level of HD 
!double precision A_He(6,3,1) ! coefficientd of 5-polynom approximation (with He)
!double precision q_H_01,q_H_10 ! collision rate of (HI + CI) from n to m level of CI 
!double precision B_H(9,3,1) ! coefficientd of 5-polynom approximation (with HI)
!end module

!module A_coeff
!The spontaneous transition probability fro CI
 !   save
 !   double precision, parameter :: A10 = 5.12*1e-008  
!end module 

!module B_coeff
! Einshtein_coefficients
!    double precision B01,B10,B20,B02,B12,B21
!end module