subroutine OpenAnalysis(file_name)

    use CO_data
    use CI_data
    use H2_data
    use hd_data
    use MCMC
    use synth_param
    !use ParamIteration

    integer num, k, i, pos, ver, s, str_num,str,ioer, n
    character(len=8) :: charadd
    character(len=4) :: chara
    character*20 buf,fmt1, charadd2
    character*200 file_name
    logical Flag
    integer helpind 
    
    helpind = 0
    
    Flag = .true.
    
    open(17, file=file_name)
    call READ_EOF(17, str_num)

    read (17, *) buf
	str = 1
    do while (TRIM(ADJUSTL(buf)) .ne. 'start')
	    read (17, *) buf
	    str = str + 1
    enddo
 	do while(str <  str_num)
            read(17,*) buf
            str = str + 1
            if (TRIM(ADJUSTL(buf)) .eq. 'CO') then
                read(17,*) num_of_level_CO
                str = str + 1

                do i = 1, num_of_level_CO
                    read(17,'(3f10.5)') (NJ_CO_fit(i,j), j = 1,3)  ! x, x+err, x-err
                    str = str + 1
                    print*, NJ_CO_fit(i,1), NJ_CO_fit(i,2)
                end do

            else if (TRIM(ADJUSTL(buf)) .eq. 'CI') then
                read(17,*) num_of_level_CI
                str = str + 1

                do i = 1, num_of_level_CI
                    read(17,'(3f10.5)') (NJ_CI_fit(i,j), j = 1,3)  ! x, x+err, x-err
                    str = str + 1
                    print*, NJ_CI_fit(i,1), NJ_CI_fit(i,2)
                end do

            else if (TRIM(ADJUSTL(buf)) .eq. 'H2') then
                read(17,*) num_of_level_H2
                str = str + 1

                do i = 1, num_of_level_H2
                    read(17,'(3f10.5)') (NJ_H2_fit(i,j), j = 1,3)  ! x, x+err, x-err
                    str = str + 1
                    print*, NJ_H2_fit(i,1), NJ_H2_fit(i,2)
                end do

             else if (TRIM(ADJUSTL(buf)) .eq. 'HD') then
                read(17,*) num_of_level_hd
                str = str + 1

                do i = 1, num_of_level_hd
                    read(17,'(3f10.5)') (NJ_HD_fit(i,j), j = 1,3)  ! x, x+err, x-err
                    str = str + 1
                    print*, NJ_HD_fit(i,1), NJ_HD_fit(i,2)
                end do

            else if (TRIM(ADJUSTL(buf)) .eq. 'syn') then
                read(17,'(1i10)') number_of_elements
                read(17,'(1i10)') iter_MC              ! NUMB OF ITERATIONS
                read(17,'(1i10)') num_MC               ! NUMB OF WALKERS
                read(17,'(1i10)') add_prior               ! use prior
                str = str + 4
                syn(1)%name(1) = 'num_d'
                syn(1)%name(2) = 'T_kin'
                syn(1)%name(3) = 'T_cmb'
                syn(1)%name(4) = 'logN0'
                syn(1)%name(5) = 'chi_UV'
                syn(1)%name(6) = 'logN0_CO'
                do k = 1,  number_of_elements
                    read(17,'(1i2,4f10.6)') syn(1)%vary(k),syn(1)%val(k), syn(1)%val_min(k), syn(1)%val_max(k), syn(1)%val_step(k)
                    str = str + 1
                    print*, k
                end do
             else if (TRIM(ADJUSTL(buf)) .eq. 'syn_SZ') then
                read(17,'(1i10)') number_of_elements
                read(17,'(1i10)') iter_MC              ! NUMB OF ITERATIONS
                read(17,'(1i10)') num_MC               ! NUMB OF WALKERS
                read(17,'(1i10)') add_prior            ! use prior
                str = str + 4
                !syn(1)%name(1) = 'T0'
                !syn(1)%name(2) = 'Te'
                !syn(1)%name(3) = 'beta'
                !syn(1)%name(4) = 'tau'
                do k = 1,  number_of_elements
                    read(17,'(1i2,4f10.6,a5)')&
                    syn(1)%vary(k),syn(1)%val(k),syn(1)%val_min(k),syn(1)%val_max(k),syn(1)%val_step(k),charadd

                !56475    n = index(charadd, ' ')
                !    if (n == 1) then
                !        charadd = charadd(2:len(charadd))
                !        GO TO (56475), 1
                !    else
                !        syn(1)%name(k) = charadd(1:n-1)
                !    end if
                    !charadd2 = TRIM(TRIM(ADJUSTL(syn(1)%name(k))))
                    
                    do while (Flag)
                        helpind = helpind + 1
                        if (charadd(helpind:helpind) == '!' .or. (helpind == len(charadd))) then
                            Flag = .false.
                        end if
                    enddo
                    syn(1)%name(k) = TRIM(ADJUSTL(charadd(1:helpind-1)))
                    str = str + 1
                    print*, k,syn(1)%name(k), syn(1)%val(k), syn(1)%val_min(k), syn(1)%val_max(k)
                end do
            end if
    enddo
    close(17)

    ! set num_of varied var's
    num_of_var = 0
    do i = 1, number_of_elements
        if(syn(1)%vary(i) == 1) then
            num_of_var = num_of_var + 1
            syn(1)%var_ind(num_of_var) = i
        end if
    end do

    !set prior
    if(add_prior == 1) call SetPrior

    end subroutine

subroutine SetPrior

    use synth_param
    use CI_data
    use H2_data
    use ParamIteration

    logical :: flag
    character*50 str
    character(len=50), dimension(6) :: args
    character(len=5) :: charind
    character(len=100) :: namepriorfile, charadd1, charadd2
    integer el_num


    !if (ind == 0) then
    !    nameSZdatafile = 'SZ_data.txt'
    !else
    !    write(charind,'(1I5)') ind
    !    nameSZdatafile = 'SZdatas\SZ_data'//TRIM(ADJUSTL(charind))//'.txt'c
    !end if


    if (ind == 0) then
        namepriorfile = 'prior.dat'
    else
        write(charind,'(1I5)') ind
        namepriorfile = 'priors/prior'//TRIM(ADJUSTL(charind))//'.dat'
    end if

    print*, "Open the "//TRIM(ADJUSTL(namepriorfile))//" file"

    open(47, file=namepriorfile) !'prior.dat'
    read(47,*) num_prior
    if (allocated(prior)) then
        deallocate (prior)
	end if
    allocate(prior(num_prior))
    do i=1, num_prior
        read(47,'(1A50)') str
        k_s = 1
        do l=1,6
            k = index(str(k_s:), ' ')
            args(l) = str(k_s:k_s+k-2)
            print *, args(l)
            k_s = k_s + k
        enddo
        read(args(1),'(i2)') prior(i)%mode
        if (prior(i)%mode == 1) then ! set prior for phys_conds
            read(args(2),'(i2)') el_num
            charadd1 = TRIM(ADJUSTL(syn(1)%name(el_num)))
            charadd2 = TRIM(ADJUSTL(args(3)))
            flag = charadd1 == charadd2
            if (flag) then
                prior(i)%ind = el_num
                read(args(4),*) prior(i)%c
                if (prior(i)%c .eq. 100.0) then   ! set box prior
                    read(args(5),*) prior(i)%down
                    read(args(6),*) prior(i)%up
                else
                    read(args(5),*) prior(i)%p  ! set +-errors
                    read(args(6),*) prior(i)%m
                end if
            end if
        else if (prior(i)%mode == 2) then ! set constraint on the production og alpha*G (see Sternberg+2014)
            aG_constrain = 1
            read(args(3),*) p_aG
            print*, 'aG_constrain =', aG_constrain
        else if (prior(i)%mode == 3) then ! set prior from lnL of H2 rotational population
            read(args(2),'(1A50)') prior(i)%namelnLH2
            open(221, file = ADJUSTL(prior(i)%namelnLH2))
            read(221,'(100f9.2)') logn_LH2
            read(221,'(100f9.2)') (loguv_LH2(k), k =1,100)
            do k =1,100
                read(221,'(100f9.2)') (lnL_H2(k,j), j =1,100)
            end do
            lnL_H2 = -lnL_H2
            close(221)
        endif
    enddo
    close(47)

end subroutine


double precision function lprior(s)
    use synth_param
    use CI_data
    use H2_data
    double precision y, chix,xl,xu,val
    double precision chi1,chi2,chi3, logn, loguv
    integer mode,i,j,k,iuv,inH
    type(synth_parameters) :: s(1)


    chix = 0
    do i=1, num_prior
        mode = prior(i)%mode
        if (mode == 1) then
            if (prior(i)%c .eq. 100.0) then
                y = s(1)%val(prior(i)%ind)
                xl = prior(i)%m
                xu = prior(i)%p
                chix = chix + 50*(1 - (dtanh(100*(y-xl)) - dtanh(100*(y-xu)))/2)
            else
                y = s(1)%val(prior(i)%ind) - prior(i)%c
                if ((dsign(1.0d0,y)+1)/2 > 0) then
                    chix = chix + (y/prior(i)%p)**2
                else
                    chix = chix + (y/prior(i)%m)**2
                endif
            end if
        else if (mode == 2) then
            if(aG_constrain == 1) then
                val = s(1)%val(5)-s(1)%val(1)-dlog10(p_aG)
                val = (1 + dtanh(100*val))/2.
                chix = chix + 50*val
            end if
        else if (mode == 3) then
            logn = s(1)%val(1)
            loguv = s(1)%val(5)
            inH = 0
            iuv = 0
            do j =1,100
                if ((logn_LH2(j)>logn) .and. (inH == 0)) inH = j
                if ((loguv_LH2(j)>loguv) .and. (iuv == 0)) iuv = j
            end do
            if ((iuv>1) .and. (inH>1)) then
                chi2 = lnL_H2(iuv-1,inH) + &
                (loguv-loguv_LH2(iuv))/(loguv_LH2(iuv-1)-loguv_LH2(iuv))*&
                (lnL_H2(iuv,inH)-lnL_H2(iuv-1,inH))
                chi1 = lnL_H2(iuv-1,inH-1) + &
                (loguv-loguv_LH2(iuv))/(loguv_LH2(iuv-1)-loguv_LH2(iuv))*&
                (lnL_H2(iuv,inH-1)-lnL_H2(iuv-1,inH-1))
                chi3 = chi1 + (chi2-chi1)*(logn-logn_LH2(inH-1))/(logn_LH2(inH)-logn_LH2(inH-1))
            else
                chi3 = lnL_H2(iuv,inH)
            end if
            chix = chix + chi3
        end if
    enddo

    lprior = chix

end function
