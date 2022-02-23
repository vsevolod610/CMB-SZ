    subroutine setfunction_HD

    call read_flower
    call read_abgrall1982
    call read_hd_gamma

    end




subroutine read_abgrall1982
    use hd_data
    use constants

    integer i,j,k,l,m,p,n,delta, nu
    double precision x
    character*200 buf

    hd_Aij(:,:,:,:) = 0.0d0

    open(11, file = 'data/hd/table1_Abgrall1982.dat')
    read(11,*)
    read(11,*)
    do m = 1, 15
        read(11,'(1a200)') buf
        read(buf(4:5),'(1i2)') i                     ! read nu
        do p = 1, (2 + 5*i)
            read(11,'(1a200)') buf
            read(buf(1:2),'(1i2)') k                 ! read nu'
            read(buf(4:5),'(1i2)') delta                 ! read J' - J
            do j = 1, 11                             ! set J
                l = j + delta
                if(l>0) then
                    read(buf(7+9*(j-1):14+9*(j-1)),*) x
                    hd_Aij(i+1,j,k+1,l) = x
                end if
            end do
        end do
    end do
    close(11)

    hd_EX(:,:) = 0.0d0
    open(11, file = 'data/hd/table2_Abgrall1982.dat')
    read(11,*)
    do i = 1,2
        read(11,'(1a200)') buf
        read(buf(4:5),'(1i2)') nu                     ! read nu
        do j = 0,10
            read(11,'(1a200)') buf
            read(buf(1:10),'(1f10.1)') x
            hd_EX(nu+1,j+1) = x*h*c/k_bol             ! in K
        end do
    end do
    close(11)

    end



    subroutine read_flower

    use collision_data

    IMPLICIT NONE

    !logical, external :: eof

    integer i, j, k
    character*200 buf
    ! set fit parameter values, according to http://ccp7.dur.ac.uk/cooling_by_HD/node2.html

    dt(1) = 2./3.
    dt(2) = 1./3.
    dt(3) = 0.5
    dt(4) = 0.5

    open(12, file = 'data/hd/collisions/Flower2000_qh1hd.dat')  ! read H1-coll_coeffs - 1
    read(12,*)
    k = 1
    a_coll(k,:,:) = 0.0
    b_coll(k,:,:) = 0.0
    c_coll(k,:,:) = 0.0
    do!do while (.not. eof(12))
        read(12,'(1A200)',END=1) buf
        read(buf(1:12),'(1i7,1i5)', END=1) i, j
        if (i > j ) then
            read(buf(15:56),'(3e14.3)', END=1) a_coll(k,i,j), b_coll(k,i,j), c_coll(k,i,j)
        end if
    end do
1   continue
    close(12)

    open(12, file = 'data/hd/collisions/Flower2000_qhehd.dat')  ! read He-coll_coeffs - 2
    read(12,*)
    k = 2
    a_coll(k,:,:) = 0.0
    b_coll(k,:,:) = 0.0
    c_coll(k,:,:) = 0.0
    do !while (.not. eof(12))
        read(12,'(1A200)', END=2) buf
        read(buf(1:12),'(1i7,1i5)', END=2) i, j
        if (i > j ) then
            read(buf(15:56),'(3e14.3)', END=2) a_coll(k,i,j), b_coll(k,i,j), c_coll(k,i,j)
        end if
    end do
2   continue
    close(12)

    open(12, file = 'data/hd/collisions/Flower2000_qoh2hd.dat')  ! read H2o-coll_coeffs - 3
    read(12,*)
    k = 3
    a_coll(k,:,:) = 0.0
    b_coll(k,:,:) = 0.0
    c_coll(k,:,:) = 0.0
    do! while (.not. eof(12))
        read(12,'(1A200)', END=3) buf
        read(buf(1:12),'(1i7,1i5)', END=3) i, j
        if (i > j ) then
            read(buf(15:56),'(3e14.3)', END=3) a_coll(k,i,j), b_coll(k,i,j), c_coll(k,i,j)
        end if
    end do
3   continue
    close(12)

    open(12, file = 'data/hd/collisions/Flower2000_qph2hd.dat')  ! read H2p-coll_coeffs - 4
    read(12,*)
    k = 4
    a_coll(k,:,:) = 0.0
    b_coll(k,:,:) = 0.0
    c_coll(k,:,:) = 0.0
    do !while (.not. eof(12))
        read(12,'(1A200)', END=4) buf
        read(buf(1:6),'(1i3,1i3)', END=4) i, j
        if (i > j ) then
            read(buf(7:39),'(3e11.3)', END=4) a_coll(k,i,j), b_coll(k,i,j), c_coll(k,i,j)
        end if
    end do
4   continue
    close(12)
    end

    subroutine read_hd_gamma

    use hd_data
    integer i,k,JX

    open(12,file = 'data/hd/hd_gamma.dat')
    do i =0, 10
        read(12,'(1i4, 11f10.3)') k, (hd_gamma(i+1,JX+1), JX = 0,10)
    end do

    hd_gamma(:,:) = hd_gamma(:,:)*1e-010
    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####################################3
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####################################3

    double precision function rate_deexc_approx(tkin,Ji,Jj,k)
    ! Calculate coll de-exitation rates for coliisions HD with H1, He, H2o, H2p in [cm^23 s^-1] according to Flower'2000 (k = 1,2,3,4 correspondingly)
    ! tkin-temperature in [K], rot level of transiton Ji -> Jj;  Ji,Jj = 1..9 (0..8)
    !
    use collision_data

    integer i,j,Ji,Jj,k
    double precision x,f,tkin,a,b,c,t
    f = 0
    t = tkin/10**3. + dt(k)

    if (Ji>Jj) f = a_coll(k,Ji,Jj) + b_coll(k,Ji,Jj)/t + c_coll(k,Ji,Jj)/t/t

    rate_deexc_approx = 10**f

    return
    end


    double precision function rate_electron_exc(Tkin,J)
    !  collisions of hd mol with electrons and transitions j->j+1

    integer J
    double precision A,B, deltaE,C,beta,D,B0,Tkin, f

    D = 8.51e-04           ! dipole moment of HD in debaye unit
    B0 = 5.538e-03*8065.0  ! rotational constant of HD in cm-1

    ! see calculation in Dickinson 1975 (eq 2.23)
    A = 2.470*D**2.*(J+1)/(2*J + 1)
    deltaE = 2.48e-04*B0*(J+1)
    C = 9.08*1e3/B0/(j+1)
    beta = 11604.505/Tkin

    f = 1.44e-06/dsqrt(Tkin)*A*dexp(-beta*deltaE)*dlog(C*deltaE + C/beta*dexp(-0.577/(1+2*beta*deltaE)))

    rate_electron_exc = f
    end

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####################################3
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####################################3
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#####################################3

    subroutine calc_fit_HD

    use synth_param
    use hd_data

    integer nuX, Ji,Jj, Jmax,i
    integer, allocatable :: indx(:)
    double precision, allocatable :: Matrix(:,:), Matrix_coll(:,:,:),Matrix_uv(:,:),Matrix_sp_deexc(:,:), A(:,:), B(:), X(:),Y(:)
    double precision n(5), nH2, delta,rate_deexc_approx, dens, d, mval, rate_electron_exc
    double precision chi_UV, a_he, f_H, N0, T_cmb, T_kin, f_shi


    T_kin = 10**syn(1)%val(2)
    dens = 10**syn(1)%val(1) !/T_kin
    T_cmb = syn(1)%val(3)
    f_H = 0.2
    a_he = 0.085
    chi_UV = 10**(syn(1)%val(5))
    f_shi = 1.0 !(1-tanh(1.15*(syn(1)%val(4)-14.22)))/2
!    print*, f_shi
!    pause
    Jmax =   num_of_level_hd + 4                       ! number of J level of HD in calculation = 0..7

    allocate(Matrix(Jmax+1,Jmax))
    allocate(Matrix_coll(Jmax,Jmax,5))
    allocate(Matrix_uv(Jmax,Jmax))
    allocate(Matrix_sp_deexc(Jmax,Jmax))
    allocate(A(Jmax-1, Jmax-1))
    allocate(B(Jmax-1))
    allocate(indx(Jmax-1))

    Matrix_coll(:,:,:) = 0.0
    Matrix(:,:) = 0.0

    !###########################    calc matrix for collisions

    do k = 1,4                   ! h1,h2,oh2,ph2
        ! de-excitation rate
        do Ji = 2, Jmax
            do Jj  = 1, Ji-1
                Matrix_coll(Ji,Jj,k) = rate_deexc_approx(T_kin,Ji,Jj,k)
            end do
        end do
        ! excitation rate
        do Ji =1, Jmax - 1
            do Jj= Ji + 1, Jmax
                delta =  hd_EX(1,Jj) - hd_EX(1,Ji)
                Matrix_coll(Ji,Jj,k) = (2.0*(Jj-1)+1.)/(2.0*(Ji-1)+1.)*dexp(-delta/T_kin)*Matrix_coll(Jj,Ji,k)
            end do
        end do
    enddo

    ! collisions with electrons
    do Ji = 1, Jmax - 1
        Matrix_coll(Ji,Ji+1,5) =  rate_electron_exc(T_kin,Ji)
    end do
    do Ji = 2,Jmax
        delta =  hd_EX(1,Ji) - hd_EX(1,Ji-1)
        Matrix_coll(Ji,Ji-1,5) = (2.0*(Ji-1-1)+1.)/(2.0*(Ji-1)+1.)*Matrix_coll(Ji-1,Ji,5)*dexp(delta/T_kin)
    end do
    !###########################    calc matrix for uv excitation
    do i = 1,5
        call trnspmtrx(Matrix_coll(:,:,i),Jmax)
    end do

    Matrix_uv(:,:) = 0.0d0
    do Ji = 1, Jmax
       do Jj  = 1, Jmax
           Matrix_uv(Ji,Jj) = hd_gamma(Ji,Jj)     ! check the order of Ji and Jj!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end do
    enddo
    call trnspmtrx(Matrix_uv,Jmax)


    Matrix_sp_deexc(:,:) = 0.0d0
    do Ji = 1, Jmax
       do Jj  = 1, Jmax
            Matrix_sp_deexc(Ji,Jj) = hd_Aij(1,Ji,1,Jj)
       end do
    end do
    call trnspmtrx(Matrix_sp_deexc,Jmax)




    ! dens = n1 + n2 + n3
    !n1 = (1-f_H)/(1+a_He-0.5*f_H)*dens ! Hydrogen density
    !n2 = 0.5*f_H/(1+a_He-0.5*f_H)*dens ! H2 dedsity
    !n3 = a_He/(1+a_He-0.5*f_H)*dens  ! Helium density

    ! dens = n1 + 2n2
    n(1) = (1-f_H)*dens ! Hydrogen density
    nH2 = 0.5*f_H*dens ! H2 dedsity
    n(2) = a_He*dens  ! Helium density

    !n(1) = (1-f_H)/(1+a_He-0.5*f_H)*dens                            !Hydrogen numb density
    !n(2) = a_He/(1+a_He-0.5*f_H)*dens                               !Helium numb density
    !nH2  = 0.5*f_H/(1+a_He-0.5*f_H)*dens                            !H2 numb density
    n(3) = nH2*9*dexp(-170.6/T_kin)/(1+9*dexp(-170.6/T_kin))        !H2-ortho J1 numb density
    n(4) = nH2/(1.0+9*dexp(-170.6/T_kin))                           !H2-para  J0 numb density
    n(5) = 1e-4*dens                                                !electron density ~ [C/H]*nH; C+ is the main source of e-



    !#### build the system of equations - Matrix(i,j)
    do i =1, Jmax
        do j = 1, Jmax
            mval = Matrix_uv(i,j)*chi_UV*f_shi + n(1)*Matrix_coll(i,j,1)+n(2)*Matrix_coll(i,j,2) + &
            n(3)*Matrix_coll(i,j,3)+n(4)*Matrix_coll(i,j,4) + n(5)*Matrix_coll(i,j,5) + Matrix_sp_deexc(i,j)
            if(mval>1.e-19) Matrix(i,j) = mval
        end do
    end do



     do i = 1, Jmax
        mval = Matrix(i,i)
        do j = 1, Jmax
            !if (i .ne. j) Matrix(i,i) = Matrix(i,i) - Matrix(j,i)
            mval = mval - Matrix(j,i)
        end do
        Matrix(i,i) = mval
    end do

    A(:,:) = 0.0
    B(:) = 0.0


    do i =2, Jmax
        do j = 1,Jmax -1
            A(i-1,j) = Matrix(i,j+1)
        end do
        B(i-1) = - Matrix(i,1)
    end do



    !##### find solution
    call ludcmp(A,Jmax-1,Jmax-1,indx,d)   ! LU decomposition
    call lubksb(A,Jmax-1,Jmax-1,indx,B)   ! solution


    hd_exc_pop(:) = 0
    hd_exc_pop(1) = 1
    do i = 1, size(B)
        hd_exc_pop(i+1)= B(i)
    end do


    end





subroutine calculate_chi_HD(chi,f)


    use hd_data
    use synth_param

    integer f
    double precision chi,x,y

    syn(1)%free_parameters = 0
    chi = 0
    if(num_of_level_hd >1) then
        do i = 1, num_of_level_hd
            x = NJ_HD_fit(i,1) - syn(1)%val(4) - dlog10(hd_exc_pop(i))
            if (x .ge. 0) then
                chi = chi + x*x/(NJ_HD_fit(i,2))**2                 ! err++
            else
                chi = chi + x*x/(NJ_HD_fit(i,3))**2                 ! err--
            end if
        end do
    end if

    if(syn(1)%chi > 1) chi = syn(1)%chi
    syn(1)%free_parameters = num_of_level_HD -  num_of_var


    if ((f == 2) .and. (syn(1)%free_parameters > 0)) then
        chi = chi/syn(1)%free_parameters
    endif
    end subroutine



    subroutine trnspmtrx(A,isize)
    integer isize,i,j
    double precision A(isize,isize),x

    do i = 2,isize
       do j = 1,i-1
            x = A(i,j)
            A(i,j) = A(j,i)
            A(j,i) = x
       end do
    enddo

    end
