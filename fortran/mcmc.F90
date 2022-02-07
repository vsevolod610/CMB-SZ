


!------------------------------------------------------------
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!
!!
!! Program for parallel computation of MCMC fit 
!!
!!
!!
!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


subroutine MCMCfit(calc_fit,calculate_chi)

    use synth_param
    use MCMC
    use CO_data
    implicit none
    include 'mpif.h'

    external calc_fit, calculate_chi 
    double precision start_time, end_time, step_time
    double precision, allocatable :: time(:)
    double precision, allocatable :: MC_prev_chi(:) 
    double precision x, chix, chi, z, koef, q, drandm,  ln, lprior
    double precision, allocatable :: vari(:)
    
    integer i,j,r,k,f, t ,l, pp, readMC, corr_limit
    character*4000 buf 
    character*20 tmp
    
    INTEGER  ierr,rc, numprocs, myid                       ! mpi_var
    
! >>>>>>>>>>>>>>>>>>>>>>>>>>> read start data >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     readMC = 0  ! readMC = 1 - read last MC file/ else readMC = 0

!############################# START MPI #######>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    
    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, IERR )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, IERR )
    print *, "Process ", myid, " of ", numprocs, " is alive"
    !call seed(1+9*myid)
    call RANDOM_SEED()
! >>>>>>>>>>>>>>>>>>>>>>>>>>> initialize Monte Carlo module >>>>>>>>>>>>>>>>>>>>>>>

    if (allocated(MC)) deallocate (MC)
    allocate (MC(num_MC))
    !allocate (random(4,num_MC))
    allocate (random(6,num_MC))
    
    call calc_fit
    call calculate_chi(chi,1)
    
    do i = 1, num_MC
        MC(i)%syn = syn
        MC(i)%chi = chi
    end do
    
    ! >>>>>>>>> MCres -auxiliary array for data transfer of MC array and others
    
    allocate (MCres(num_of_var+1))  
    
    if(myid == 0) then                      
        if(readMC  ==0 ) then                 
           call genMC_init(calc_fit, calculate_chi)                   ! generate new walkers and find chi_min for them
        else if(readMC == 1) then            ! read positions from the 'MC.dat'
           open(101, file = 'MC.dat')
           do k = 1, num_MC
                read(101,'(5f20.10)') MC(k)%chi, (MCres(i), i=1,num_of_var) !(667,'('// trim(fmt1)//'f20.10,1i5)') 
                do i = 1, num_of_var
                    t =  syn(1)%var_ind(i)
                    MC(k)%syn(1)%val(t) = MCres(i)
                end do
           end do
           close(101)
        end if 
        
        allocate(time(numprocs))
        
        if (allocated(MCsave)) deallocate (MCsave)
        allocate(MCsave(num_MC))
        MCsave = MC
        if (allocated(MCmin)) deallocate (MCmin)
        
        ! >>>>>>>>> create MCmin and sort it
        allocate(MCmin(num_MC))          
        MCmin = MC
        call sortMCmin          
	    chi_min = MCmin(1)%chi

        ! >>>>>>>>> create MCchain
        allocate (MCchain(iter_MC))
        do i=1,iter_MC
            allocate (MCchain(i)%par(num_MC,num_of_var))
            allocate (MCchain(i)%chi(num_MC))
            allocate (MCchain(i)%auto(num_MC))
            allocate (MCchain(i)%ave(num_of_var))
            allocate (MCchain(i)%disp(num_of_var))
            allocate (MCchain(i)%med(num_MC,num_of_var))
            allocate (MCchain(i)%sca(num_MC,num_of_var))
            MCchain(i)%auto = 1
        enddo
        
        allocate(GRvalue(3,num_of_var))           !Gelman-Rubin test
        allocate(Autocorr(401,num_of_var))
        allocate (vari(num_of_var*2))
        allocate(MC_prev_chi(num_MC))
        
        
        f = 0
        
    
        
        open (18, file="resMC.dat")
        open (19, file="MCchain.dat")
        open (20, file="GRvalues.dat")
        open (21, file="GRvalues2.dat")

    
    end if
  
! >>>>>>>> wait till Initialization performed in 0 process.
    call MPI_BCAST(chi_min,1,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
    
!   generate random array and distribute to it to all proc
    if (myid==0) then
        do k = 1, num_MC
            do i = 1, 6
    	        !random(i,k) = drandm(0)
                call random_number(random(i,k))
                !print*, random(i,k)
            end do
	    end do
    end if
 	call MPI_BARRIER(MPI_COMM_WORLD, IERR)
	do k = 1, 6
		call MPI_BCAST(random(k,:),num_MC,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
    end do	

!   save array MC(k)%chi for counting number of displacement points
    do k = 1, num_MC      ! MC(k)_data_exchange
            if(myid == 0)   call write_char_MC(k,buf)
            call MPI_BCAST(buf,4000,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
            if(myid .ne. 0) call read_char_MC(k,buf)
    end do
 
    
    
    
    
    
    
    
    
    
! >>>>>>>>>>>  Start iteration >>>>>>>>>>>>>>>>>>>
    pp = 1
    i = 0
    f = 1
    do while ((pp == 1) .and. (i<iter_MC))
        call MPI_BARRIER( MPI_COMM_WORLD , IERR)
        i = i + 1
        iter = i
        
        start_time = MPI_WTIME(IERR)
        step_time = start_time
        
        if(( i == 1 ) .and. ( myid == 0 )) call calcMC_stats ! calc stats for the first iteration
        if(myid == 0) then      !  save MC(k)%chi values for the last iteration
            MC_prev_chi = 0
            do k = 1, num_MC
                MC_prev_chi(k) = MC(k)%chi
            end do
        end if
      
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!######################## SHIFT the first half of walkers (K <= N/2) >>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
        !open(334, file = 'coeff.dat')
        !write(334,*) 'koef', 'ln', 'q', 'x'
        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        do t=myid+1, num_MC/2 ,numprocs  
            k = t
            z = 0
            
            syn(1)%chi = 0.
            call shiftAffine_MC(k,z,1)      !  set syn%
            if(syn(1)%chi == 0.) then 
                call calc_fit
                call calculate_chi(chix, 1)
            else 
                chix = syn(1)%chi
            end if
            koef = 1.0 
            !koef = chi_min/syn(1)%free_parameters
            if(syn(1)%free_parameters .le. 0) koef = 1.0
            ln = (MC(k)%chi)/2 - (chix)/2
            if (add_prior == 1) ln = ln + (lprior(MC(k)%syn) - lprior(syn))/2
            q = dexp(ln/koef)*z**(num_of_var-1)
            x = random(1,k)
         !   write(334,'(2f20.5,1e20.5,2f20.5)') koef, ln, q, z, x
	        if ((x<=q)) then 
                MC(k)%syn = syn
                MC(k)%chi = chix
            endif 
        enddo
        !close(334)
        
        
        ! >>>>>>>>  check time  >>>>>>>>>
        if (myid .ne. 0) then
            do k = 2, numprocs
                if (myid == k-1) then 
                    call MPI_send(MPI_WTIME(IERR) - start_time,1,MPI_DOUBLE,0,k,MPI_COMM_WORLD,ierr)
                endif
            enddo
        else
            time(1) = (MPI_WTIME(IERR) - step_time)
            do k = 2, numprocs
                call MPI_recv(time(k),1,MPI_DOUBLE,MPI_ANY_SOURCE,k,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            enddo
            print *, '       First half - time :', maxval(time), minval(time)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        step_time = MPI_WTIME(IERR)
        
        
        
! >>>>>>>>>> MC(k)_data_exchange: #all proc send MC(k) to #0 proc
        if (myid .ne. 0) then
            do  t=myid+1, num_MC/2 ,numprocs       !>>>>>> #all proc send MC to proc#0
                k = t
                call write_char_MC(k,buf)
                !if(k == 2) print*, buf
                call MPI_send(buf,4000,MPI_CHARACTER,0,k,MPI_COMM_WORLD,ierr)
            end do
        else if (myid == 0) then                  ! proc#0 accept new MC to #all proc
		    do k = 1, numprocs-1                
			    do t = k+1, num_MC/2 ,numprocs		
                    r = t
        		    call MPI_recv(buf,4000,MPI_CHARACTER,MPI_ANY_SOURCE,r,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                    call read_char_MC(r,buf)
			    end do
		    end do
        end if
 !>>>>>>>>>  print the number of shifted walkers >>>>>>>>>>  
        if(myid == 0) then 
	        l = 0
            do t = 1, num_MC/2
                k = t
                if(MC_prev_chi(k) .ne. MC(k)%chi) l = l  + 1
            end do
            print*, '   Shifted points ', l ,' from ', num_MC/2
        end if 
        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
 
          
        
! >>>>>>>>>> MC(k)_data_exchange: #all proc send MC(k) to #0 proc
! >>>>>>>>>> update MC for the first part of proc (#id < numMC/2)
        
	    do t = 1, num_MC/2	    
	    	k = t
            if( myid == 0) call write_char_MC(k,buf)
        	call MPI_BCAST(buf,4000,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        	if( myid .ne. 0) call read_char_MC(k,buf)
	    end do
	    call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        step_time = MPI_WTIME(IERR)
        
        
         
        
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!#################################### SHIFT the second half of walkers (K > N/2) >>>>>>>>>>>>           
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        

	    do t=myid+1+num_MC/2, num_MC ,numprocs
            k = t
            z = 0
            syn(1)%chi = 0.
            call shiftAffine_MC(k,z,2)      !  set syn%
            if(syn(1)%chi == 0.) then 
                call calc_fit
                call calculate_chi(chix, 1)
            else 
                chix = syn(1)%chi
            end if
     
            koef = 1.0 
            !koef = chi_min/syn(1)%free_parameters
            if(syn(1)%free_parameters .le. 0) koef = 1.0
            ln = MC(k)%chi/2 - chix/2
            if (add_prior == 1) ln = ln + (lprior(MC(k)%syn) - lprior(syn))/2
            q = dexp(ln/koef)*z**(num_of_var-1)
            !x = random(1,k)
            x = random(6,k)
	        if ((x<=q)) then 
                MC(k)%syn = syn
                MC(k)%chi = chix
            endif 
        enddo

! >>>>>>>>  check time  >>>>>>>>>
        if (myid .ne. 0) then
            do k = 2, numprocs
                if (myid == k-1) then 
                    call MPI_send(MPI_WTIME(IERR) - step_time,1,MPI_DOUBLE,0,k,MPI_COMM_WORLD,ierr)
                endif
            enddo
        else
            time(1) = (MPI_WTIME(IERR) - step_time)
            do k = 2, numprocs
                call MPI_recv(time(k),1,MPI_DOUBLE,MPI_ANY_SOURCE,k,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            enddo
            print *, '       Second half - time :', maxval(time), minval(time)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        step_time = MPI_WTIME(IERR)
        
        
! >>>>>>>>>> MC(k)_data_exchange: #all proc send MC(k) to #0 proc
        if (myid .ne. 0) then 
            do  t=myid+1+num_MC/2, num_MC ,numprocs       
                k = t
                call write_char_MC(k,buf)
                call MPI_send(buf,4000,MPI_CHARACTER,0,k,MPI_COMM_WORLD,ierr)
            end do
        else if (myid == 0) then                  
            do k = 1, numprocs-1                
			    do t = num_MC/2+k+1, num_MC ,numprocs		
	 			    r = t
                    call MPI_recv(buf,4000,MPI_CHARACTER,MPI_ANY_SOURCE,r,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                    call read_char_MC(r,buf)
			    end do
		    end do
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, IERR)

        

        !>>>>>>>>>  print the number of shifted walkers >>>>>>>>>>  
        if(myid == 0) then 
	        l = 0
            do t = num_MC/2 + 1, num_MC
                k = t
                if(MC_prev_chi(k) .ne. MC(k)%chi) l = l + 1
            end do
            print*, '   Shifted points = ',l,' from ', num_MC/2
        end if 

 

! >>>>>>>>>> update MC for the second part of proc (#id > numMC/2)
! >>>>>>>>>> MC(k)_data_exchange: send updated MC from #0 to #all proc

        do t = num_MC/2+1, num_MC
            k = t 
	    	if (myid == 0) call write_char_MC(k,buf)
        	call MPI_BCAST(buf,4000,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        	if (myid .ne. 0) call read_char_MC(k,buf)
        end do

        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
	    step_time = MPI_WTIME(IERR)
        
        
            
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>            MC statistic             >>>>>>>>>>>>>>>>>>>>>
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
! >>>>>> find chi_min and distribute it to all processes
        if(myid == 0) then 
            min_MC = 1
            chi_min = MC(1)%chi
            
            do k=1, num_MC
                if (MC(k)%chi <= chi_min) then
                    min_MC = k                  
                    chi_min = MC(min_MC)%chi
                    syn = MC(min_MC)%syn
                endif
            enddo
            print *, 'chis:', chi_min, chi_min_glob
            if (chi_min < chi_min_glob) then
                chi_min_glob = chi_min
            end if
            
            koef = chi_min/syn(1)%free_parameters
            if(syn(1)%free_parameters .le. 0) koef = 1.0
        end if
        call MPI_BCAST(chi_min,1,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(min_MC,1,MPI_double_precision,0,MPI_COMM_WORLD,ierr)  !  the number of the walker which have min chi


        if(myid == 0) then 
            
            call calcMC_stats

            x = 1.0d0
            do k=1,num_of_var
                vari((k-1)*2+1) = confi(k,1) ! average value
                vari((k-1)*2+2) = confi(k,2) ! dispersion
            enddo
            
            if (iter > 3) then
                do l=1,num_of_var
                    x = x*GRvalue(1,l)       ! Gelman-Rubin test for convergence
                enddo
            end if
            
            t = 0                            ! number of dhifted walkers
            do k = 1, num_MC      
                if (MC_prev_chi(k) .ne. MC(k)%chi) t = t + 1
            end do
   
        	!write(18,'(1I4, 1I8,<9+ 2*num_of_var>f20.9)') i, syn(1)%free_parameters, chi_min/syn(1)%free_parameters, chi_min, dlog10(x), 1.0d0*corr_lenght, 1.0d0*k_clen, MC(k_clen)%chi, chi_min, MC_moments(1), MC_moments(2), (vari(k), k=1,2*num_of_var)
	        write(18,'(1I4, 1I8, <7+2*num_of_var>f19.9)') iter, syn(1)%free_parameters, koef, chi_min, t*1.0d0/num_MC,  dlog10(x), 1.0d0*corr_lenght, MCchain(iter)%mom(1), MCchain(iter)%mom(2), (vari(k), k=1,2*num_of_var)
            print *, '        MC statistic time :', (MPI_WTIME(IERR) - step_time)
            end_time = MPI_WTIME(IERR)
            write (*,'(1i4, 1f7.2)') iter,  end_time - start_time 
        endif
   
        
!>>>>>>>>>>> Adjust MC to make faster convergence
        if(adjust_MC == 1) then
            if(myid == 0 ) then 
    !	        call calc_quick_MC_stats
                if (iter == 1) then 
                    corr_limit = 1
                else if (corr_lenght > idint(0.7d0*corr_limit)) then 
                !else if(corr_lenght > 20)
                    call ModifyMC_Affine(corr_limit, calc_fit, calculate_chi)
                endif
            end if
	    
    ! >>>>>> and distribute adjusted MC to all processes	    
	        do k = 1, num_MC      
                    if( myid == 0)   call write_char_MC(k,buf)
                    call MPI_BCAST(buf,4000,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                    if( myid .ne. 0) call read_char_MC(k,buf)
                end do
 	        call MPI_BARRIER(MPI_COMM_WORLD, IERR)
        end if

! >>>>>> distribute new random numbers to all processes
        if (myid==0) then  
		    do k = 1, num_MC
	            do t = 1,6
                    call random_number(random(t,k))
                    !random(t,k) = drandm(0)
                end do
	        end do
	    end if
	    do k = 1, 6
    		call MPI_BCAST(random(k,:),num_MC,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
        end do
        
    enddo

    if(myid == 0) then 
       close(18)
       close(19)
       close(20)
       close(21)
    endif  

    call MPI_FINALIZE(rc)     
end subroutine


subroutine SynParall
    use MCMC
    use synth_param   
    use CO_data
    use SZ_data
    character*200 file_name
    
    external calc_fit_CI,calculate_chi_CI,calc_fit_CO,calculate_chi_CO,calc_fit_H2,calculate_chi_H2, calc_fit_all,calculate_chi_all,calc_fit_HD,calculate_chi_HD,calc_fit_SZ,calculate_chi_SZ
    ! set probability function 
    call SetFunction_CO   
    call SetFunction_CI
    call SetFunction_HD    
    call SetFunction_SZ
    
    file_name = TRIM(ADJUSTL('startSZ.sss'))
    call OpenAnalysis(file_name)
    call MCMCfit(calc_fit_SZ,calculate_chi_SZ)         ! start Monte Carlo Markov Chain calculation 
end 


