

subroutine SetFunction_CO
    
    use CO_data
    use collision_data
    
    integer i,j
    character*500 buf 
    character*10 buf2 
    double precision x
  
    
    !###################################    set parameters for CO molecules
    
    !>>>>>>>>> Einstein coeffs
    num_J = 14
    
    AJ_CO(1) = 0.0
    AJ_CO(2)  = 7.203e-08
    AJ_CO(3)  = 6.910e-07
    AJ_CO(4)  = 2.497e-06
    AJ_CO(5)  = 6.126e-06
    AJ_CO(6)  = 1.221e-05
    AJ_CO(7)  = 2.137E-05
    AJ_CO(8)  = 3.422E-05
    AJ_CO(9)  = 5.134E-05
    AJ_CO(10) = 7.330E-05
    AJ_CO(11) = 1.006E-04
    AJ_CO(12) = 1.339E-04 
    AJ_CO(13) = 1.735E-04
    AJ_CO(14) = 2.200E-04
    AJ_CO(15) = 2.739E-04 
 
    do i =1,15
        EJ_CO(i) = 2.766*(i*(i-1)) !-(i-1)*(i-2))  !in K see Draine 2000 (tab. 5.1)
    end do
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>> Set collision coeffs
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    ! collisions with H
    
    open(17, file='data/CO_data.dat')
    read (17, '(A)') buf  
    do while (TRIM(ADJUSTL(buf)) .ne. '1 CO-H from Walker et al.(2015)')
        read (17, '(1A100)') buf  
    enddo
    do while (TRIM(ADJUSTL(buf)) .ne. '!COLL TEMPS')
        read (17, '(1A100)') buf    
    enddo
    qCOH(:,:)=0.
    read(17,'(1A300)') buf
    buf = TRIM(ADJUSTL(buf))
    do i = 0,28
        read(buf(10*i+1:10*i+10),'(1f10.1)') x
        qCOH(1,i+1) = x 
    end do 
    read(17,*) 
    do i =1, 105
        read(17,'(1A400)') buf
        do j = 1,29
                read(buf(19+10*(j-1):19+10*j-1),'(1e10.2)') qCOH(i+1,j)
        end do
    end do

     
   !collisions with H2
   !para-H2
    do while (TRIM(ADJUSTL(buf)) .ne. '2 CO-pH2 from Yang et al. (2010)')
        read (17, '(1A100)') buf  
    enddo
    do while (TRIM(ADJUSTL(buf)) .ne. '!COLL TEMPS')
        read (17, '(1A100)') buf    
    enddo
    qCOH2p(:,:)=0.
    read(17,'(1A300)') buf
    buf = TRIM(ADJUSTL(buf))
    do i = 0,24
       read(buf(10*i+1:10*i+10),'(1f10.1)') x
       qCOH2p(1,i+1) = x 
    end do 
    read(17,*) 
    do i =1, 105
        read(17,'(1A300)') buf
        do j = 1,25
            read(buf(19+10*(j-1):19+10*j-1),'(1e10.3)') qCOH2p(i+1,j)
        end do
    end do
        
        
    !ortho-H2
    do while (TRIM(ADJUSTL(buf)) .ne. '3 CO-oH2 from Yang et al. (2010).')
        read (17, '(1A100)') buf  
    enddo
    do while (TRIM(ADJUSTL(buf)) .ne. '!COLL TEMPS')
        read (17, '(1A100)') buf    
    enddo
    qCOH2o(:,:)=0.
    read(17,'(1A300)') buf
    buf = TRIM(ADJUSTL(buf))
    do i = 0,24
        read(buf(10*i+1:10*i+10),'(1f10.1)') x
        qCOH2o(1,i+1) = x 
    end do 
    read(17,*) 
    do i =1, 105
        read(17,'(1A300)') buf
        do j = 1,25
            read(buf(19+10*(j-1):19+10*j-1),'(1e10.3)') qCOH2o(i+1,j)
        end do
    end do
 
    
    ! collisions with Helium
    do while (TRIM(ADJUSTL(buf)) .ne. '4 CO-He4 from Cecchi-Pestellini et al (2002)')
        read (17, '(1A100)') buf  
    enddo
    do while (TRIM(ADJUSTL(buf)) .ne. '!COLL TEMPS')
        read (17, '(1A100)') buf    
    enddo    
       
    qCOHe(:,:)=0.
    read(17,'(1A150)') buf
    buf = TRIM(ADJUSTL(buf))
    do i = 0,9
        read(buf(10*i+1:10*i+10),'(1f10.1)') x
        qCOHe(1,i+1) = x 
    end do 
    read(17,'(1A300)') buf
    do i =1, 105
        read(17,'(1A300)') buf
        read(buf(8:9),'(1f2.0)') qCOHe(i+1,1)
        read(buf(14:15),'(1f2.0)') qCOHe(i+1,2)
        do j = 1,10
            read(buf(18+9*(j-1):26+9*(j-1)),'(1e8.1)') qCOHe(i+1,j+2)
        end do
    end do
    close(17)
end subroutine 

    
    
subroutine calc_fit_CO          ! procedure for calculation of CO excitation (CMB + collisions)
                                ! J level - is the max rot_lev of CO, J_level can be 0, 1, 2..  <= 5
    
    use CO_data
    use CI_data
    use collision_data  
    use synth_param
!    
    
    integer i ,j ,k, m, indx(num_J), pprint, k_ind
    double precision n_para, n_ortho, d
    double precision n1,n2,n3,n,f, alpha
    double precision, allocatable :: Matrix(:,:), Matrix_coll(:,:,:), Matrix_cmb(:,:), A(:,:), B(:), X(:),Y(:)
    double precision u,v,w,h,r,z,x1,x2
    double precision dens, T_kin, T_cmb, a_he, f_H, mval
    
    
    !dens = syn(1)%val(1)
    dens = 10**syn(1)%val(1)
    T_kin = 10**syn(1)%val(2)
    !dens = 10**syn(1)%val(1)/T_kin
    T_cmb = syn(1)%val(3)
    f_H = h2_fraction
    a_he = 0.085
    
    
    allocate(Matrix(num_J+1,num_J+1))
    allocate(Matrix_coll(num_J+1,num_J+1,4))
    allocate(Matrix_cmb(num_J+1,num_J+1))
    allocate(A(num_J, num_J))
    allocate(B(num_J)) 

    !! the definition of Matrix_coll(i,j) for as function of T_kin
    Matrix_coll(:,:,:) = 0.0
    
    k = 1
    k_ind = 0
    if(T_kin .le. qCOH(1,1)) then 
        k = 1
    else if (T_kin .lt. qCOH(1,29)) then
        do while (T_kin .gt. qCOH(1,k))
            k = k + 1   ! for T< 3000 K   !                find point for cooresponding T in the grid of collisin coff approximation
        end do
        k = k - 1  
        k_ind = 1
    else 
        k = 29
    end if
     
    !###############################################################################    (1) !Set Matrix_coll(i,j,1)  - for Hydrogen 
    m = 1     
    do j = 2, num_J+1
        do i  = 1, j-1
            m = m + 1
            Matrix_coll(i,j,1) = qCOH(m,k) 
            if(k_ind .eq. 1) Matrix_coll(i,j,1) = qCOH(m,k) + (qCOH(m,k+1)-qCOH(m,k))*(T_kin-qCOH(1,k))/(qCOH(1,k+1)-qCOH(1,k))              ! linear appr
        end do
    end do
        
    ! for i<j
    
    do j =1, num_J          
        do i=j+1, num_J+1
            Matrix_coll(i,j,1) = (2.0*(i-1)+1.)/(2.0*(j-1)+1.)*dexp(-(EJ_CO(i)-EJ_CO(j))/T_kin)*Matrix_coll(j,i,1)        
        end do
    end do
    
    
    
    
    !###############################################################################    (2) !Set Matrix_coll(i,j,2)  - for H2_orhto
    k_ind = 0
    k = 1
    if(T_kin .le. qCOH2o(1,1)) then 
        k = 1
    else if (T_kin .lt. qCOH2o(1,25)) then
        do while (T_kin .gt. qCOH2o(1,k))
            k = k + 1   ! for T< 3000 K   !                find point for cooresponding T in the grid of collisin coff approximation
        end do
        k = k - 1  
        k_ind = 1
    else 
        k = 25
    end if
!    if
 !   do while (T_kin .gt. qCOH2o(1,k))
 !       k = k+1   ! for T< 3000 K   !                find point for cooresponding T in the grid of collisin coff approximation
 !   end do
 !   k = k -1  
    
    m = 1     
    do j = 2, num_J+1
        do i  = 1, j-1
            m = m + 1
            Matrix_coll(i,j,2) = qCOH2o(m,k)
            Matrix_coll(i,j,3) = qCOH2p(m,k)
            if(k_ind == 1) Matrix_coll(i,j,2) = qCOH2o(m,k) + (qCOH2o(m,k+1)-qCOH2o(m,k))*(T_kin-qCOH2o(1,k))/(qCOH2o(1,k+1)-qCOH2o(1,k))   ! linear appr
            if(k_ind == 1) Matrix_coll(i,j,3) = qCOH2p(m,k)+ (qCOH2p(m,k+1)-qCOH2p(m,k))*(T_kin-qCOH2p(1,k))/(qCOH2p(1,k+1)-qCOH2p(1,k))
        end do
    end do
        
    ! for i<j
    
    do j =1, num_J          
        do i=j+1, num_J+1
            Matrix_coll(i,j,2) = (2.0*(i-1)+1.)/(2.0*(j-1)+1.)*dexp(-(EJ_CO(i)-EJ_CO(j))/T_kin)*Matrix_coll(j,i,2)   
            Matrix_coll(i,j,3) = (2.0*(i-1)+1.)/(2.0*(j-1)+1.)*dexp(-(EJ_CO(i)-EJ_CO(j))/T_kin)*Matrix_coll(j,i,3)        
        end do
    end do    
    !###############################################################################    (3)  !Set Matrix_coll(i,j,3)  - for H2_para
    !m = 1     
    !do j = 2, num_J+1
    !    do i  = 1, j-1
    !        m = m + 1
    !        Matrix_coll(i,j,3) = qCOH2p(m,k)+ (qCOH2p(m,k+1)-qCOH2p(m,k))*(T_kin-qCOH2p(1,k))/(qCOH2p(1,k+1)-qCOH2p(1,k))   ! linear appr
    !    end do
    !end do
    !    
    !! for i<j
    !
    !do j =1, num_J          
    !    do i=j+1, num_J+1
    !        Matrix_coll(i,j,3) = (2.0*(i-1)+1.)/(2.0*(j-1)+1.)*dexp(-(EJ_CO(i)-EJ_CO(j))/T_kin)*Matrix_coll(j,i,3)        
    !    end do
    !end do    
    
    
    !###############################################################################    (4)   !Set Matrix_coll(i,j,4)  - for He
    k_ind = 0
    k = 1
    if(T_kin .le. qCOHe(1,1)) then 
        k = 1
    else if (T_kin .lt. qCOHe(1,10)) then
        do while (T_kin .gt. qCOHe(1,k))
            k = k + 1   ! for T< 3000 K   !                find point for cooresponding T in the grid of collisin coff approximation
        end do
        k = k - 1  
        k_ind = 1
    else 
        k = 10
    end if
    !do while (T_kin .ge. qCOHe(1,k))
    !    k = k+1   ! for T< 500 K   !
    !end do
!    k = k -1
    !if(T_kin .gt. 500) then 
    !    print*, 'Temperature for He_coll coefficient coll exceeds the threshold, T_kin>500'
    !    pause 
    !end if 
    
    do m = 2,106 
        i = int(qCOHe(m,1)) + 1
        j = int(qCOHe(m,2)) + 1
        if((i.le.num_J+1).and.(j.le.num_J+1)) then 
            Matrix_coll(j,i,4) = qCOHe(m,k+2)
            if(k_ind == 1) Matrix_coll(j,i,4) = qCOHe(m,k+2) + (qCOHe(m,k+3)-qCOHe(m,k+2))*(T_kin-qCOHe(1,k))/(qCOHe(1,k+1)-qCOHe(1,k))   ! linear appr
        end if
    end do      
    ! for i<j
    
    do j =1, num_J          
        do i=j+1, num_J + 1
            Matrix_coll(i,j,4) = (2.0*(i-1)+1.)/(2.0*(j-1)+1.)*dexp(-(EJ_CO(i)-EJ_CO(j))/T_kin)*Matrix_coll(j,i,4)        
        end do    
    end do    
    !###############################################################################   (5)   !Set  Matrix_cmb
  
    Matrix_cmb(:,:) = 0.0
    do i =2, num_J+1
           !Matrix_cmb(i-1,i) = AJ_CO(i)*dexp((EJ_CO(i)-EJ_CO(i-1))/T_cmb)/(dexp((EJ_CO(i)-EJ_CO(i-1))/T_cmb)-1.0)        
           !Matrix_cmb(i,i-1) = AJ_CO(i)/(dexp((EJ_CO(i)-EJ_CO(i-1))/T_cmb)-1.0)*(2.0*(i-1)+1.0)/(2.0*(i-2)+1.0)
           Matrix_cmb(i-1,i) = AJ_CO(i)/( 1 - dexp(-(EJ_CO(i)-EJ_CO(i-1))/ T_cmb ))        
           Matrix_cmb(i,i-1) = AJ_CO(i)*dexp(-(EJ_CO(i)-EJ_CO(i-1))/T_cmb)/(1- dexp(-(EJ_CO(i)-EJ_CO(i-1))/T_cmb))*(2.0*(i-1)+1.0)/(2.0*(i-2)+1.0)
           !Matrix_cmb(i-1,i) = AJ_CO(i)/( 1 - dexp(-(EJ_CO(i)-EJ_CO(i-1))/ T_cmb )) *(2.0*(i-1)+1.0)/(2.0*(i-2)+1.0)*dexp(-(EJ_CO(i)-EJ_CO(i-1))/T_cmb)
           !Matrix_cmb(i,i-1) = AJ_CO(i)/(1- dexp(-(EJ_CO(i)-EJ_CO(i-1))/T_cmb)) !*dexp(-(EJ_CO(i)-EJ_CO(i-1))/T_cmb)
    
    end do
    
   !#########################################  the finish of matrix determination 
    
    
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
    beta_factor  = 1 !/(1.0+0.5*10)
    
    !#### formulate the system of equations - Matrix(i,j) 
    do i =1, num_J+1 
        do j = 1, num_J+1 
            Matrix(i,j) = Matrix_CMB(i,j)* beta_factor + n1*Matrix_coll(i,j,1)+n_ortho*Matrix_coll(i,j,2)+n_para*Matrix_coll(i,j,3)+n3*Matrix_coll(i,j,4)
            !Matrix(i,j) = n1*Matrix_coll(i,j,1)+n_ortho*Matrix_coll(i,j,2)+n_para*Matrix_coll(i,j,3)+n3*Matrix_coll(i,j,4)
        end do
    end do
    
    !!>>>>>>>>>>>>>>>>>>>>>        test
    !open(134,file = 'Matrix_CMB.dat')
    !Write(134,*) 'CMB'
    !do i =1, 5 
    !       write(134,'(5e15.5)') (Matrix_CMB(i,j), j = 1, 5 )
    !end do
    !Write(134,*) 'coll'
    !do i =1, 5 
    !       write(134,'(5e15.5)') (n1*Matrix_coll(i,j,1)+n_ortho*Matrix_coll(i,j,2)+n_para*Matrix_coll(i,j,3)+n3*Matrix_coll(i,j,4), j = 1, 5 )
    !end do
    !Write(134,*) 'total'
    !do i =1, 5 
    !       write(134,'(5e15.5)') (Matrix(i,j), j = 1, 5 )
    !end do
    !close(134)
    !
    
    
    !do i = 1, num_J+ 1
    !    do j = 1, num_J + 1
    !        if (i .ne. j) Matrix(i,i) = Matrix(i,i) - Matrix(j,i)
            !Matrix(i,i) = Matrix(i,i) - Matrix(j,i)

    !    end do
    !end do  
    
    do i = 1, num_J+ 1
        mval = Matrix(i,i) 
        do j = 1, num_J+ 1
            !if (i .ne. j) Matrix(i,i) = Matrix(i,i) - Matrix(j,i)
            mval = mval - Matrix(j,i) 
        end do
        Matrix(i,i) = mval
    end do  
    
    A(:,:) = 0.0
    B(:) = 0.0
    
    do i =1, num_J
        do j = 1,num_J
            A(i,j) = Matrix(i,j+1) 
        end do
        B(i) = - Matrix(i,1)
    end do
    
    
    !##### find solution
    call ludcmp(A,num_J,num_J,indx,d)   ! LU decomposition
    call lubksb(A,num_J,num_J,indx,B)   ! solution
    
    NJ_CO(:) = 0.0
    NJ_CO(1) = 1.0 
    do i = 2, num_J+1
        NJ_CO(i) = B(i-1)
    end do
    
    end subroutine    

    
subroutine calculate_chi_CO(chi,f)
    
    use CO_data
    use collision_data  
    use synth_param
    
    integer f
    double precision chi,x
    
    syn(1)%free_parameters = 0
    chi = 0
    !open(101, file = 'calc_chi.dat')
    if(num_of_level_CO >1) then 
        do i = 1, num_of_level_CO
            x = NJ_CO_fit(i,1) - syn(1)%val(4) -dlog10(NJ_CO(i))
            if (x .ge. 0) then 
                chi = chi + x*x/(NJ_CO_fit(i,2))**2
            else
                chi = chi + x*x/(NJ_CO_fit(i,3))**2
            end if
        end do    
    end if
    if(syn(1)%chi > 1) chi = syn(1)%chi
    syn(1)%free_parameters = num_of_level_CO -  num_of_var
    
    
    if ((f == 2) .and. (syn(1)%free_parameters > 0)) then
        chi = chi/syn(1)%free_parameters
    endif
    !write(101,'(5f10.5)') chi, (syn(1)%var(i), i = 1,4)
    !close(101)
end subroutine
    
    

SUBROUTINE ludcmp(a,n,np,indx,d)                      ! from "numerical recipies"
INTEGER n,np,indx(n),NMAX
double precision d,a(np,np),TINY
PARAMETER (NMAX=500,TINY=1.0e-20) 
double precision i,imax,j,k
REAL aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
d=1. !No row interchanges yet.
do i=1,n !Loop over rows to get the implicit scaling informa
    aamax=0. !tion.
do j=1,n
if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
enddo
if (aamax.eq.0.) pause 'singular matrix in ludcmp' !No nonzero largest element.
vv(i)=1./aamax !Save the scaling.
enddo 
do j=1,n !This is the loop over columns of Crout's method.
do i=1,j-1 !This is equation (2.3.12) except for i = j.
sum=a(i,j)
do k=1,i-1 
sum=sum-a(i,k)*a(k,j)
enddo 
a(i,j)=sum
enddo 
aamax=0. !Initialize for the search for largest pivot element.
do i=j,n !This is i = j of equation (2.3.12) and i = j+1: ::N
sum=a(i,j) !of equation (2.3.13).
do k=1,j-1
sum=sum-a(i,k)*a(k,j)
enddo
a(i,j)=sum
dum=vv(i)*abs(sum) !Figure of merit for the pivot.
if (dum.ge.aamax) then !Is it better than the best so far?
imax=i
aamax=dum
endif
enddo 
if (j.ne.imax)then !Do we need to interchange rows?
do k=1,n !Yes, do so...
dum=a(imax,k)
a(imax,k)=a(j,k)
a(j,k)=dum
enddo
d=-d !...and change the parity of d.
vv(imax)=vv(j) !Also interchange the scale factor.
endif
indx(j)=imax
if(a(j,j).eq.0.) a(j,j)=TINY
!If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!For some applications on singular matrices, it is desirable to substitute TINY for zero.
if(j.ne.n)then !Now, nally, divide by the pivot element.
dum=1./a(j,j)
do i=j+1,n
a(i,j)=a(i,j)*dum
enddo
endif
enddo !Go back for the next column in the reduction.
return
END

SUBROUTINE lubksb(a,n,np,indx,b)
INTEGER n,np,indx(n)
double precision a(np,np),b(n)
INTEGER i,ii,j,ll
double precision sum
ii=0 !When ii is set to a positive value, it will become the index
!of the rst nonvanishing element of b. We now do
!the forward substitution, equation (2.3.6). The only new
!wrinkle is to unscramble the permutation as we go.
do i=1,n
ll=indx(i)
sum=b(ll)
b(ll)=b(i)
if (ii.ne.0)then
do j=ii,i-1
sum=sum-a(i,j)*b(j)
enddo
else if (sum.ne.0.) then
ii=i !A nonzero element was encountered, so from now on we will
!have to 
endif !do the sums in the loop above.
b(i)=sum
enddo
do  i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
sum=b(i)
do j=i+1,n
sum=sum-a(i,j)*b(j)
enddo 
b(i)=sum/a(i,i) !Store a component of the solution vector X.
enddo 
return !All done!
    END

    
subroutine calc_texc_co(n,tkin)
    use CO_data
    use collision_data  
    use synth_param
    
    integer i
    double precision n, tkin
   
    syn(1)%val(1) = n    
    syn(1)%val(2) = tkin
    call calc_fit_CO 
    T_exc_co = -55.31/dlog(NJ_CO(5)/9.)    
    
end
    