subroutine genMC_init(calc_fit,calculate_chi)
    
    use MCMC
    use synth_param
    external calc_fit, calculate_chi

    
    integer i,j
    double precision chi
    print*, 'genMC_start'
    
    call calc_fit
    call calculate_chi(chi,1)
    
    syn2 = syn
    MC(1)%syn = syn
    j = 2  
    
    do i=1, num_MC
        MC(i)%syn = syn2 
        call shift_MC(i)    !>>>>>  generate new random syn for MC(i) 
        MC(i)%syn = syn
        
        call calc_fit
        call calculate_chi(chi,1)
        
        MC(i)%chi = chi
        
    !>>>>>  print number of iteration were performed
        if ( i == j) then 
		    print*, '    ', j, ' from ', num_MC, 'walkers were generated'
		    j = j + 200
        end if
    enddo
    syn = syn2
    
end subroutine


subroutine shift_MC(i_MC)
    
    use MCMC
    use synth_param
    
    integer i,mode,i_MC,k,ind
    double precision x(num_of_var),xx
    
    syn = MC(i_MC)%syn
    
    mode = 2                           ! 
    
    select case(mode)
    case(1)
        do k=1, num_of_var
            call randNorm2(x(k))
        end do   
    case(2)              ! odnorodnoe distribution in [-1,+1]^num_of_var cube
        do k=1, num_of_var
            !x(k) = -1.d0 + 2.0*drandm(0)
            call random_number(xx)
            x(k) = -1.d0 + 2.0*xx
        end do  
    end select    
    
    do k = 1, num_of_var
        ind = syn(1)%var_ind(k)
        syn(1)%val(ind) = MC(i_MC)%syn(1)%val(ind) + x(k)*syn(1)%val_step(ind)
        if (syn(1)%val(ind) < syn(1)%val_min(ind)) syn(1)%val(ind) = syn(1)%val_min(ind)
        if (syn(1)%val(ind) > syn(1)%val_max(ind)) syn(1)%val(ind) = syn(1)%val_max(ind)
    end do

    end subroutine shift_MC

subroutine sortMCmin
    
    use MCMC
    
    integer i, k
    type(MC_type) MCs
    
    do i=1,num_MC-1
        do k=i, num_MC
            if (MCmin(k)%chi < MCmin(i)%chi) then
                MCs = MCmin(i)
                MCmin(i) = MCmin(k)
                MCmin(k) = MCs
            endif
        enddo
    enddo
    
end subroutine


subroutine calcMC_stats(ind)
    
    use MCMC
    use synth_param
    implicit none 
     
    integer i,j,k,l,i_min,f,t,i_max, n,iter_border, ind
    double precision x,y, koef
    character*20 fmt1,fmt2,fmt3
    character*4000 buffer 
    character (len = 100) :: nameChainfile
    character (len = 5) :: charind
    
    logical :: exist
! >>>>>>>>>>>>> Calc first two moments of Chi^2 distributions
    if (calc_moments == 1) then
    
        MC_moments = 0
        do k=1, num_MC
            MC_moments(1) = MC_moments(1) + MC(k)%chi
            MC_moments(2) = MC_moments(2) + MC(k)%chi**2
        enddo
        MC_moments = MC_moments/num_MC
        MC_moments(2) = dsqrt(MC_moments(2) - MC_moments(1)**2)
        MCchain(iter)%mom = MC_moments 
        !>>>>>>>> calc MC moments for last 4 corr_length
        i_min = max(1, iter - 4*corr_lenght)
        MC_moments = 0
        do i=i_min, iter
            MC_moments = MC_moments + MCchain(i)%mom
        enddo
        MC_moments = MC_moments/(iter - i_min+1)
    end if
        
        
    
    ! >>>>>>>>>>>>> Array includes all values near chi^2 minimum throught the chain
    if (calc_MCmin == 1) then
        do i = 1, num_MC
            k = num_MC
            if (MC(i)%chi < MCmin(1)%chi) then
                k = 0
            else
                do while ((MC(i)%chi < MCmin(k)%chi))
                    k = k-1
                enddo
                if (MC(i)%chi == MCmin(k)%chi) then   !>>>> repeat existent value, so start again
                    k = num_MC                        
                endif
            endif
            if (k < num_MC) then             ! add new syn to MCmin massive
                if (k < num_MC-1) then       ! shift numeration of values with MC%chi > MC(i)%chi to "+ 1"    
                    do l=0, num_MC-k-2           
                        MCmin(num_MC-l) = MCmin(num_MC-l-1)
                    enddo
                endif
                MCmin(k+1) = MC(i)           ! and add MC(i)%chi to (k+1)-place in MC_min
            endif
        enddo
        
        !>>>>>>>>>   print_MCmin
        if(write_MCmin == 1) then 
            open(867, file="MCmin.dat")
            MCsave = MC
            MC = MCmin
            do k=1, num_MC
                call write_char_MC(k, buffer)
                write (867,'(<num_of_var+1>f20.10)') MCres(num_of_var+1), (MCres(i), i=1,num_of_var) 
            enddo
            MC = MCsave
            close(867)
        end if
    end if

    
    
    
    
!>>>>>>>>>>>>     stats of chain
    if(stats_chain == 1) then 

        corr_lenght = 1
        do k=1, num_MC
            call write_char_MC(k, buffer)
            MCchain(iter)%par(k,:) = MCres(:)
            MCchain(iter)%chi(k) = MCres(num_of_var+1)
            if (MC(k)%chi < chi_emin) then
                chi_emin = MC(k)%chi
            endif
        
            !>>>>>>>>>>>  calc MCchain(iter)%auto(k) 
            MCchain(iter)%auto(k) = 1
            if (iter > 1) then  
                if (MCchain(iter)%chi(k) == MCchain(iter-1)%chi(k)) then
                    MCchain(iter)%auto(k) =  MCchain(iter-1)%auto(k)+1
                endif
            endif 
   
            !>>>>>>>>>>>  calc corr_lenght
            !if (k .ne. num_MC) then
                if (MCchain(iter)%auto(k) > corr_lenght) then 
                    corr_lenght = MCchain(iter)%auto(k)
                    k_clen = k
                endif
            !endif
        enddo
        
        !>>>>>>>>>>  calc median <x> and scatter <x^2> for 4 last iterations
        if (iter>3) then

            ! summa from iter/2 till iter
            if (iter*0.5-iter/2 == 0) then
                MCchain(iter)%med = MCchain(iter-1)%med+(MCchain(iter)%par-MCchain(iter/2-1)%par)/(iter/2+1)
                MCchain(iter)%sca = MCchain(iter-1)%sca+(MCchain(iter)%par**2-MCchain(iter/2-1)%par**2)/(iter/2+1)
            else
                MCchain(iter)%med = MCchain(iter-1)%med+(MCchain(iter)%par-MCchain(iter-1)%med)/((iter+3)/2)
                MCchain(iter)%sca = MCchain(iter-1)%sca+(MCchain(iter)%par**2-MCchain(iter-1)%sca)/((iter+3)/2)
            endif
        else
            if (iter == 3) then
                MCchain(iter)%med = (MCchain(iter-2)%par + MCchain(iter-1)%par + MCchain(iter)%par)/3
                MCchain(iter)%sca = (MCchain(iter-2)%par**2 + MCchain(iter-1)%par**2 + MCchain(iter)%par**2)/3
            endif
        endif
    end if

!!>>>>   calculate mean autocorr length
!    if(calc_Autocorr == 1) then
!        AutoCorr = 0       
!        do k=1, num_MC
!            do i = 1, num_of_var          
!                if (MCsave(k)%syn(1)%var(i) == MC(k)%syn(1)%var(i)) AutoCorr = AutoCorr+1
!            end do
!        enddo
!        AutoCorr = AutoCorr/num_MC       
!    end if
!    
! >>>>>>>>>>>>> Calculations of autocorrelation time
    if ((calc_Autocorr == 1) .and. (iter > 100)) then
            
        Autocorr = 0
        k = 1
        open (24, file="Autocorr.dat")
        do t=0, min(iter/4, 400)
            do l=1, num_of_var
                do i=iter/2, iter-t
                !do k=1, num_MC
                    Autocorr(t+1,l) = Autocorr(t+1,l) + (MCchain(i+t)%par(k,l) - MCchain(iter)%med(k,l))*(MCchain(i)%par(k,l) - MCchain(iter)%med(k,l))
                !enddo
                enddo
            Autocorr(t+1,l) = Autocorr(t+1,l)/(MCchain(iter)%sca(k,l) - MCchain(iter)%med(k,l)**2)/((iter+1)/2)
            enddo
            write(24, '(1I4,<num_of_var>f20.9)') t, (Autocorr(t+1,l), l=1,num_of_var)
        enddo   
        close(24)
    endif
    
    
! >>>>>>>>>>>  Gelman-Rubin test
    if (calc_GRtest == 1) then
        if (iter > 3) then
            call GRtest
            write(20,'(1I4, <num_of_var>f20.9)') iter, (GRvalue(1,k), k=1,num_of_var)
            write(21,'(1I4, <2*num_of_var>f20.9)') iter, (GRvalue(2,k), GRvalue(3,k), k=1,num_of_var)
        endif
    endif
        

    
!   >>>>>>>>>>>>>>> write MC
    if (write_data == 1) then 
        open(667, file="MC.dat")
        do k=1, num_MC
            if (write_data == 1) write (667,'(<num_of_var+1>f20.10)') MCchain(iter)%chi(k), (MCchain(iter)%par(k,i), i=1,num_of_var) !, MCchain(iter)%auto(k)  
            !if (write_data == 1) write (667,'(<num_of_var>f20.10)')  (MCchain(iter)%par(k,i), i=1,num_of_var) !, MCchain(iter)%auto(k)  
        enddo
        close(667)
        ! >>>>>>>>> write MC100 file
        if (iter/20 == iter*1.0d0/20) then
            !open(668, file="MC100.dat")
            !write(668, '(1A500)') fit_descr
            !write(668, '(1A500)') el_descr
            !do i=max(1, iter-600), iter
            !    do k=1, num_MC
            !        write (668,'(<num_of_var+1>f20.10, 1i5)') MCchain(i)%chi(k), (MCchain(i)%par(k,l), l=1,num_of_var), MCchain(i)%auto(k)  
            !    enddo
            !enddo
            !close(668)
            
            
            
            if (ind == 0) then
                nameChainfile = "chainconsum.dat"
            else 
                write(charind,'(1I5)') ind
                nameChainfile = 'OUT\chainconsum'//TRIM(ADJUSTL(charind))//'.dat'c
            end if
            
            inquire(file=nameChainfile, exist=exist)
            if (exist) then
              open(667, file=nameChainfile, status="old", position="append", action="write")
            else
              open(667, file=nameChainfile, status="new", action="write")
            end if
!            open(667, file="chainconsum.dat")
!            open(667, file="chainconsum.dat", status="old", position="append", action="write")
            open(993, file = 'border.dat')
            read(993, '(1i4)') iter_border
            close(993)
            if(iter > iter_border) then 
                do i=max(1, iter-19), iter
                    do k=1, num_MC
                        write (667,'(1I4,<num_of_var+1>f20.10)') iter, MCchain(i)%chi(k), (MCchain(i)%par(k,j), j=1,num_of_var) 
!                         write (667,'(1I4,<num_of_var>f20.10)') iter, (MCchain(i)%par(k,j), j=1,num_of_var) 
                    end do
                enddo
            end if
            close(667)
        endif
    end if

    if ((write_data_final == 1) .and. (iter == iter_MC)) then 
        open(667, file="MC_final.dat")
        do k=1, num_MC
            !if (write_data == 1) write (667,'(<num_of_var+1>f20.10)') MCchain(iter)%chi(k), (MCchain(iter)%par(k,i), i=1,num_of_var) !, MCchain(iter)%auto(k)  
            write (667,'(<num_of_var>f20.10)')  (MCchain(iter)%par(k,i), i=3,num_of_var) !, MCchain(iter)%auto(k)  
        enddo
        close(667)
    end if

    
!   >>>>>>>>>>>> write some additional MC files 
    if (show_addMC == 1) then    
        koef = chi_emin/syn(1)%free_parameters
        open(7, file='MCbound.dat')
        do k=1,num_MC
            do i=5, iter
                if (MCchain(i)%auto(k) == 1) then
                    if (MCchain(i)%chi(k) < chi_emin + koef) then
        			    write (7,'(2i10,<num_of_var+1>f20.10)') i, k, MCchain(i)%chi(k), (MCchain(i)%par(k,l), l=1,num_of_var) 
                    endif
                endif
            enddo
        enddo
        close(7)
        
        if(iter == 1) then 
            MCDisp_min = MC_moments(2)
        else if(MC_moments(2) < MCDisp_min) then 
            MCDisp_min = MC_moments(2)
            open(7, file ='MC_disp.dat')
            do k = 1, num_MC
	            write (7,'(2i10,<num_of_var+1>f20.10)') iter, k, MCchain(iter)%chi(k), (MCchain(iter)%par(k,l), l=1,num_of_var) 
            end do
            close(7)     
        end if
    end if

    
    
    
! >>>>>>>>>>>>>>  calculate distribution of the parameters: average and dispersion
    if (calc_stats == 1) then    
        MCchain(iter)%ave = 0
        MCchain(iter)%disp = 0
        do k=1, num_MC
            call write_char_MC(k, buffer)
            do i=1, num_of_var
                MCchain(iter)%ave(i) = MCchain(iter)%ave(i) + MCres(i)
                MCchain(iter)%disp(i) = MCchain(iter)%disp(i) + MCres(i)**2
            enddo
        enddo
        MCchain(iter)%ave  =  MCchain(iter)%ave/num_MC
        MCchain(iter)%disp  =  dsqrt(MCchain(iter)%disp/num_MC - MCchain(iter)%ave**2)
        write(19, '(1I4,<3*num_of_var>f20.9)') iter, (MCchain(iter)%ave(i), i=1,num_of_var), (MCchain(iter)%disp(i), i=1,num_of_var), (MCchain(iter)%par(1,l), l=1,num_of_var) 
        
        
        !probability search
        if (allocated(confi)) deallocate(confi)
        allocate(confi(num_of_var,2))
    
        do i=1, num_of_var
            confi(i,1) = MCchain(iter)%ave(i)
            confi(i,2) = MCchain(iter)%disp(i)
        enddo
    end if

end subroutine


subroutine shiftAffine_MC(i,z,s)
    
    use MCMC
    use synth_param
    ! implicit none 
    ! if s == 0 then shift i in all massive of points
    ! if s == 1 then shift i relative to points with k > num_MC 
    ! if s == 2 then shift i relative to points with k < num_MC 
    
    !integer i,s, i_min,i_max,k,j,t
    !double precision x, z,drandm
    integer i,s,j,t,k
    double precision x, z, i_min,i_max,drandm
    if (z == 0) then
        if (s == 0) then 
            !x = drandm(0)
            call random_number(x)
            i_min = 0.5
            i_max = num_MC+0.5
        endif
	    if (s == 1) then 
            x = random(2,i)
            i_min = num_MC/2 + 0.5
            i_max = num_MC+0.5
        endif
	    if (s == 2) then
            x = random(3,i)	
            i_min = 0.5
            i_max = num_MC/2 + 0.5
        endif
        j = idnint(i_min + (i_max-i_min)*x)
        !only in the case (s == 0)
        if (j==i) j=max(abs((idnint(i_max)-j)/2), abs(j/2))  
    else
        j = idint(z)
    endif    
    
    if (s == 0) x = drandm(0)
    if (s == 1) x = random(4,i)	
    if (s == 2) x = random(5,i)
        
    z = ((aff-1)*x+1)**2/aff
    !syn = MC(i)%syn
    !if( s == 1) then 
    !    i_min = num_MC/2 + 0.5
    !    i_max = num_MC
    !else if( s == 2) then 
    !    i_min = 1
    !    i_max = num_MC/2
    !end if
    !
    !if (z == 0) then
    !    if(s == 0) then 
		  !  x = drandm(0)
	   ! else
		  !  x = random(2,i)	
	   ! end if
    !    j = i_min+ idint((i_max-i_min)*x)
    !    if (j==i) j=max(abs((i_max-j)/2), abs(j/2))  ! only in the case ( s == 0)
    !else
    !    j = idint(z)
    !endif    
    !
    !if(s == 0) then 
	   ! x = drandm(0)
    !else
	   ! x = random(3,i)	
    !end if
    !
    !z = ((aff-1)*x+1)**2/aff
    !if( s == 0) z = x 
    
    do k = 1, num_of_var
        t = syn(1)%var_ind(k)
        syn(1)%val(t) = MC(j)%syn(1)%val(t) + z*(MC(i)%syn(1)%val(t)-MC(j)%syn(1)%val(t))  
        !if (syn(1)%val(t) < syn(1)%val_min(t)) syn(1)%val(t) = syn(1)%val_min(t)
        !if (syn(1)%val(t) > syn(1)%val_max(t)) syn(1)%val(t) = syn(1)%val_max(t)
        if (syn(1)%val(t) < syn(1)%val_min(t)) syn(1)%chi = 1.e+20
        if (syn(1)%val(t) > syn(1)%val_max(t)) syn(1)%chi = 1.e+20
    end do
 end subroutine shiftAffine_MC


!subroutine calc_quick_MC_stats
!    
!    use MCMC
!    use synth_param
!    implicit none 
!     
!    integer i,j,k
!
!    MC_moments = 0
!    do k=1, num_MC
!        MC_moments(1) = MC_moments(1) + MC(k)%chi
!        MC_moments(2) = MC_moments(2) + MC(k)%chi**2
!    enddo
!    MC_moments = MC_moments/num_MC
!    MC_moments(2) = dsqrt(MC_moments(2) - MC_moments(1)**2)
!    MCchain(iter)%mom = MC_moments 
! 
!    do k=1, num_MC
!        MCchain(iter)%chi(k) = MC(k)%chi
!        MCchain(iter)%auto(k) = 1
!        if (iter > 1) then
!            if (MCchain(iter)%chi(k) == MCchain(iter-1)%chi(k)) then
!                MCchain(iter)%auto(k) =  MCchain(iter-1)%auto(k)+1
!            endif
!        endif 
!    enddo
!    corr_lenght = 1
!    do k=1, num_MC
!        if (k .ne. num_MC) then
!            if (MCchain(iter)%auto(k) > corr_lenght) then 
!                corr_lenght = MCchain(iter)%auto(k)
!                k_clen = k
!            endif
!        endif
!    enddo
!
!end subroutine


subroutine GRtest
    
    use MCMC
    
    double precision, allocatable :: medi(:)
    double precision  W, B, med, sca, med2, sca2, r
    integer num, l, i, k
    
    num = iter - iter/2
    allocate(medi(num_MC))
    do i=1, num_of_var
        med = 0
        sca = 0
        med2 = 0
        do k=1,num_MC
            med = med + MCchain(iter)%med(k,i)
            sca = sca + MCchain(iter)%med(k,i)**2
            med2 = med2 + (MCchain(iter)%sca(k,i)-MCchain(iter)%med(k,i)**2)
        enddo
        B = (sca-num_MC*(med/num_MC)**2)/(num_MC-1)  ! variance of the chain means
        W = med2/num_MC      ! mean of the variances of each chain
        !GRvalue(1,i) = B/W+num*1.0d0/(num-1)
        GRvalue(1,i) = dsqrt(B/W+1-1.0d0/num)
        GRvalue(2,i) = dsqrt(W)
        GRvalue(3,i) = dsqrt(B)
    enddo
    
end subroutine

SUBROUTINE READ_EOF(n,nstrok)

	INTEGER, INTENT(in) ::n 
	INTEGER, INTENT(out) :: nstrok 

	INTEGER :: nstrok_counter
	CHARACTER*30 a

	ioer=0
	nstrok_counter=0
!open(4, file = 'read_eof.dat')
	DO WHILE (ioer.eq. 0)
	READ(n,'(1a30)',iostat=ioer) a
	!write(4,*) n,a, nstrok_counter+ 1
	nstrok_counter=nstrok_counter + 1
	ENDDO
	REWIND(n)
 !close(4)
	nstrok=nstrok_counter -1

RETURN

END 

subroutine randNorm(x)
    
    double precision s, x, u, v
    
    s = 2
    do while (s>1)
!        u = -1+2*dfloat(rand(0))
!        v = -1+2*dfloat(rand(0))
!        s = u**2+v**2
        u = -1+2*drandm(0)
        v = -1+2*drandm(0)
        s = u**2+v**2
        if (s<1) x = u*dsqrt(-2*dlog(S)/S)
    enddo
    
end subroutine 

subroutine randNorm2(x)
    
    double precision s, x, u, v
    
    s = 2
    do while (s>1)
        u = -1+2*drandm(0)
        v = -1+2*drandm(0)
        s = u**2+v**2
        if (s<1) x = u*dsqrt(-2*dlog(S)/S)
    enddo
    u = drandm(0)
    x = dsign(dabs(x) + (u/(1-u)), x) !2.0d0/3*x + v/3*dlog(1.0d0/(1-u))
    
end subroutine 


subroutine modifyMC_Affine(corr_limit, calc_fit, calculate_chi)
    
    use MCMC
    use synth_param
    implicit none 
    external calc_fit, calculate_chi
    integer mask(num_MC)
    integer corr_limit, coef, j, i,k, num
    double precision tune_old, chix, x,  z, rand, drandm
    integer, allocatable :: poin(:)
    
    
    corr_limit = corr_limit + 1 
    
    mask = 0
    num = 0
    do k = 1, num_MC
        if(MC(k)%chi < MC(min_MC)%chi + 2*MC_moments(2)) then 
            num = num + 1
            mask(num) = k      ! >>> get new numeration for "good" walkers
        end if
    end do
    
    
    if(num .ge. 0.05*num_MC) then     
        do k = 1, num_MC
            if((MCchain(iter)%auto(k) > 20) .and. (mask(k) == 0))then 
                corr_limit = 0 
                x = random(4,k)
                j =1+ idint((num-1)*x)
                z = mask(j)*1d0
                x = random(1,k)
                j = 1+ idint((num-1)*x)
                j = mask(j)
                call shiftAffine_MC(j,z,0)
                call calc_fit
                call calculate_chi(chix, 1)
                print*, ' Modify position of #',k,'walker from chi = ', MC(k)%chi, ' to new chi = ', chix
                MC(k)%syn = syn 
                MC(k)%chi = chix
            end if
        end do
    end if
   
end subroutine




subroutine read_char_MC(k,buf)
    use synth_param
    use MCMC
    use CO_data
    integer i,k,num_buf, left, right
    character*4000 buf 

        
    read(buf(1:30),"(1f30.10)") MC(k)%chi 
    do i =1, num_of_var
        ind = syn(1)%var_ind(i)
        left = (i)*30 + 1
        right = (i+1)*30 
        read(buf(left:right),"(1f30.20)") MC(k)%syn(1)%val(ind)
    end do

end subroutine 

subroutine write_char_MC(k,buf)
    !implicit none
    use synth_param
    use MCMC
    
    integer i,k,num_buf, left, right,ind
    character*4000 buf 

    do i = 1, num_of_var
        ind = syn(1)%var_ind(i)
        MCres(i) =  MC(k)%syn(1)%val(ind)
    end do
    MCres(num_of_var+1) = MC(k)%chi

    write(buf(1:30),"(1f30.10)") MC(k)%chi 
    do i =1, num_of_var
        ind = syn(1)%var_ind(i)
        left = (i)*30 + 1
        right = (i+1)*30 
        write(buf(left:right),"(1f30.20)") MC(k)%syn(1)%val(ind) 
    end do

end subroutine 
