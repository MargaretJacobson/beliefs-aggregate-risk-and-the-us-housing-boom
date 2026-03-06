program main
    use omp_lib
    use parameters
	use mod_globals
    use mod_matlab
    use mod_Rouwenhorst
    use mod_discreteAR1
    use mod_orderpack
    use mod_grids_fixed
    use mod_interp
    use mod_plot
    use mod_root1dim
    use mod_initial_y
    use mod_initial_dist
    use mod_initial_dist2
    use mod_coeffs
	use mod_bellman_bequest
	use mod_continuation_value
    use mod_bellman_new
    use mod_EEE
    use mod_VFI
	use mod_function_interp
	use mod_simulate_steady
    use mod_simulate_bor_full
    use mod_KS_regression2
    use mod_simulate_learning
    use mod_csv
    implicit none
    real (8) :: diff_solve,t1,t2,pss,t11,t22,t111,t222,diffp,tp1,tp2, tp11, tp22,HS, tpsim1, tpsim2, tsim1, tsim2,thetar, thetah,liq,inc,theta_b(2,nBSim)
	real(8), dimension(learn) :: pstar2, pCKW,pstarsort
	real (8) :: LaborH,R2p_all,Yreveal
    real (8) :: L(nL), H(nH), Y(nY),Ybar(1), diff_solve_step(2),R(nR),Calibration(15),av_scale,Lsim(nLsim),Bsim(nBsim)
    real (8), allocatable :: p(:),pprime(:,:,:,:),a_raw(:,:,:), FYY(:,:),Trans(:,:)
   	real (8), save :: loge(T,nE), Fee(T,nE,nE)
	real (8), allocatable :: counts_h35(:,:), counts_h49(:,:), counts_h64(:,:), counts_h80(:,:)
	real (8), allocatable :: counts_h(:,:), counts_r(:,:)
    real (8) :: bequest,  Calibration_untargeted(14)
	real (8), allocatable ::  bottom_NW(:), mid_NW(:), top_NW(:), top10_NW(:), top1_NW(:), Agg_inc(:)
    integer ::  iL,iLL,iT,iH,iHH,iY,iYY,iP,iter,iAgg,iM, iJ,iter2, countp, init_loc(nY*nlamb,nY*nlamb), iI,iE, countY, nA,iB,nA2,next,nextplus,iD,iDD,ilo_b(2,nBsim)
    integer ::  init_loc2,iA1,iA2,Aloc,Aloc2,iA3,iA4,iA, hloc,rloc,sim_end,ilamb,iphi,iBB,iExo
	integer, allocatable :: count_buy(:),count_rent(:),count_LTV(:),count_LTV2(:),count_hy(:),count_med1(:), count_med2(:), count_bequest(:), count_bequest_tot(:)
	real (8), allocatable :: H_W_50(:), H_NW_50_bequest(:),NW_to_inc_50(:), LTV_10(:), LTV_50(:), LTV_90(:),LTV_mean(:),LTV_mean2(:), House_value_earn_10(:), House_value_earn_50(:), house_value_earn_90(:)
	real (8), allocatable :: NW_20(:), NW_40(:), NW_50(:), NW_60(:), NW_80(:), NW_90(:), NW_99(:), NW_med1_print(:), NW_med2_print(:),H_NW_10(:),H_NW_50(:),H_NW_90(:),rental_sim(:),afford_p10(:),afford_p50(:)
    real(8), allocatable :: c_b(:,:,:,:,:,:,:,:,:),c_interp(:,:,:,:,:,:,:,:,:),LL_b(:,:,:,:,:,:,:,:,:),L_interp(:,:,:,:,:,:,:,:,:),HH_b(:,:,:,:,:,:,:,:,:),H_interp(:,:,:,:,:,:,:,:,:), BB_b(:,:,:,:,:,:,:,:,:), B_interp(:,:,:,:,:,:,:,:,:),V_b(:,:,:,:,:,:,:,:,:),W_own(:,:,:,:,:,:,:,:,:),V_sim_interp(:,:,:,:,:,:,:,:,:), q(:,:,:,:,:,:,:,:), q_interp(:,:,:,:,:,:,:,:)
    real(8), allocatable :: L_sim_steady(:,:,:), H_sim_steady(:,:,:), B_sim_steady(:,:,:), tax_sim_steady(:,:,:), C_sim_steady(:,:,:)
    integer, allocatable :: age_sim_bor(:,:), ET_loc(:,:),rent_buy(:,:), max_loc_sim_star(:,:), max_loc_b_steady(:,:,:)
	integer, allocatable :: rent_buy_steady(:,:,:)
    integer, dimension(T_sim) :: YT_loc,AggT_loc,LambT_loc
    integer, dimension(T_sim,nD+2) :: max_loc_b
    real(8), dimension(T_sim) :: YT
	real(8), dimension(nY*nlamb*nphi,T_sim) :: pregZ
	real(8), allocatable :: L_sim_dist(:,:,:), B_sim_dist(:,:,:), H_sim_dist(:,:,:),C_sim_dist(:,:,:), E_sim_dist(:,:,:), Tax_sim_dist(:,:,:),wealth_sim_dist(:,:,:), max_loc_sim_dist(:,:,:)
    integer, allocatable :: age_sim_dist(:,:,:), rent_buy_sim_dist(:,:,:)
    real(8), allocatable :: H_NW(:,:),H_W(:,:), H_W_own(:,:),H_NW_own(:,:)
    real(8), allocatable ::  NW_med1(:,:), NW_med2(:,:)
	real(8),allocatable :: AggH_b(:,:), AggH_s(:), AggS(:,:)
    real(8), allocatable :: AggL_new(:,:), AggH_new(:,:), AggH_sNew(:), Agg
    real(8), save, dimension(T_sim) :: qstar,pstar,pstarold,Lstar_b,Bstar_b,Hstar_b, Hstar_s, Cstar_b,Cstar_l,income_b,tax_b, housing_tax_b
    real(8), save, dimension(T_sim) :: adjust_b, ret_b,employment_b,homeownership,foreclosure
    real(8), dimension(T_sim,mom+1) :: a_out
	real(8), dimension(T_sim,nY*nlamb,mom+1) :: a_out_all
	character(20), dimension(5) :: strname
    real(8), dimension(11) :: dimensions
    real(8), dimension(T) :: e
	integer, dimension(2) ::  ilo_pstar
	real(8), dimension(2) :: theta_pstar
	real(8), allocatable :: rental(:,:,:,:)
	integer, allocatable :: ilo_pprime(:,:,:,:,:,:,:)
	real(8), allocatable :: theta_pprime(:,:,:,:,:,:,:)
	real(8), allocatable :: pprimein(:,:,:,:,:,:)
	integer, dimension(nY*nlamb*nphi,nY*nlamb*nphi) :: count
    real(8), dimension(Jret-1) :: chi
    integer, dimension(nY*nlamb*nphi,nY*nlamb*nphi) :: counts
    integer, dimension(nY*nlamb*nphi) :: counts1,countZ
    real(8), dimension(nY*nlamb*nphi,nY*nlamb*nphi) :: R2p
    real(8), save, dimension(nY*nlamb*nphi,nY*nlamb*nphi,T_sim-burn-1) :: preg
	real(8), allocatable :: C_sim_star(:,:), income_star(:,:), L_sim_star(:,:), H_sim_star(:,:),B_sim_star(:,:), tax_star(:,:), employment_star(:,:)
	integer, dimension(T_sim,nD+2) :: countD
    !real(8), save, dimension(nY*nlamb,T_sim-burn-1) :: qstar_print, pstar_print, Lstar_print, Cstar_b_print, Hstar_print
    real(8), dimension(nY*nlamb*nPhi,nY*nlamb*nPhi,mom+1) :: a_new !, a
	  real(8), dimension(iter_total) :: diff_total
    !real(8), dimension(iter_total2) :: diff_total2
	real(8), allocatable :: avg_buy_rent(:,:), under_35(:,:)
	real(8), allocatable :: earn(:,:)
	real(8), allocatable :: Agg_nw(:,:),med_ratio(:),NW_to_inc(:,:)
	real(8), allocatable :: LTV(:,:), LTV2(:,:),House_value_earn(:,:),Afford(:,:)
	integer, allocatable :: With_mortgage(:), with_Heloc(:)
	integer :: iDprint
	real(8), dimension(2,nY*nlamb*nphi) :: pprime_solve
	  !real(8), dimension(iter_total,nY,nY,mom+1) :: total_p
     real(8), dimension(T) :: Variance
DOUBLE PRECISION, EXTERNAL          :: golden, fx
	nA=0
	nA2=0
 Aloc=0
    if (type_sim<2) then
        nA=1
		nA2=1
		
    else 
        nA=size(adj)
		nA2=NA*no_MIT+(1-no_MIT)
    end if
	print*, 'nA,nA2', NA,nA2
	print*, 'adj', adj
	
	Aloc=(nA+1)/2
	Aloc2=(nA2+1)/2
	print*, 'ALoc,Aloc2',Aloc,Aloc2
	allocate(a_raw(nY*nlamb*nphi,nY*nlamb*nphi,mom+1))
	allocate(FYY(nY*nlamb*nphi,nY*nlamb*nphi))

	print*, 'size(a_raw)', shape(a_raw)
	
		!*************************************
	!	Initializing simulation and calibration parameters
	!*************************************
	if (type_sim==1) then 
	
		sim_end=T_Sim
	else
		sim_end=learn
	end if 
	
	print*, 'Solve for coeffs', solve_for_coeffs 	! if =1 then it will iterate on the coefficients. Otherwise =0 it will not. Use 0 for steady state, learning, and housing boom simulation. Setting=1 will put the credit condition states to 1 automatically
	print*, 'high_nlamb', high_nlamb 
	print*, 'type_sim', type_sim
	print*, 'nlamb', nlamb
	print*, 'no_MIT', no_MIT
	print*, 'no_pref', no_pref
	print*, 'low_rate', low_rate

	!*************************************************
	!	Writing parameters to input file so that the text reflects the most recent results
	!************************************************
	OPEN(3, FILE = './results/computational-appendix.txt',status='replace',action='write')
		write(3,"(a,I6,a)") '\newcommand{\nagents}{$',I,'$}'
		write(3,"(a,I4,a)") '\newcommand{\Tsim}{$',int(4.4e3+600),'$}' 
		write(3,"(a,I3,a)") '\newcommand{\burn}{$',burn,'$}'
		write(3,"(a,I2,a)") '\newcommand{\corenum}{$',omp_get_max_threads(),'$}'
		write(3,"(a,f6.4,a)") '\newcommand{\gain}{$',gain_param,'$}'
		write(3,"(a,f6.4,a)") '\newcommand{\gainannual}{$',(1d0-(1d0-gain_param)**(0.5d0)),'$}'
	close(3)
	print*, 'gain annual', (1d0-(1d0-gain_param)**(0.5d0))
	print*, 'gain raw', gain_param
	print*, 'gain2', gain2

    !*************************************************
    !   Initializing Arrays
    !*************************************************
    iter2=0
    diffp=1d0
	a_raw=0d0
   
	countP=nP
		rental_sim=0d0
 !   do while (diffp>tol2 .and. iter2<iter_total2)
    iter2=iter2+1    
	Y=0d0

   
	if (CKW .eqv. .True. ) then
        open(21,file='pstar.csv',status='old',action='read')
            do it=1,learn
                read(21,*) pstar2(it)
            end do
        close(21)
        pstarsort=pstar2
        call INSSOR(pstarsort)
        !print*, pstarsort
        countp=2
        ! There is no rounding function in fortran
        ! pstar2(1) is 0 so we are ignorin this
        do iP=2,learn
            if (ceiling(pstarsort(iP)*round_size)/=ceiling(pstarsort(iP-1)*round_size)) then
                countp=countp+1
                pCKW(countp)=pstarsort(iP-1)
                print*, countp, pCKW(countp)
            end if
        end do

        pCKW(2)=pCKW(3)-(pCKW(4)-pCKW(3)) 
        pCKW(1)=pCKW(2)-(pCKW(3)-pCKW(2))
        
        
        pcKW(countp+1)=maxval(pstarsort(1:learn))
        countp=countp+1
        pCKW(countp+1)=(pCKW(countp)-pCKW(countp-1))+pCKW(countp) 
        countp=countp+1

        
        print*, minval(pstarsort(1:learn)), maxval(pstarsort(1:learn))
        print*, pCKW(1), pCKW(countp)
    end if
		allocate(P(countP))
		allocate(AggH_s(countP))
		allocate(Trans(countP,nlamb)) 
	if (CKW .eqv. .True. ) then 
		P=pCKW(1:countP)
	else
		p=linspace(0.3d0,0.9d0,countP) 
	end if 
	print*,'p', p
	 call sub_grids_fixed(L,H,Lsim,Bsim,R,Y,P,FYY,chi,e,countP, AggH_s,loge,Fee,iter2,Variance,Trans)
	 print*, 'L main after sub', L
	print*,'H main0', H
print*, 'R: main 0',R
print*, 'B', B
print*, 'Bsim', Bsim
print*, 'L', L
print*, 'Lsim', Lsim
print*, 'shrink', shrink
print*, nY*nlamb*nphi, countP, NA, Ny,nlamb,nphi, nP
    allocate(pprime(nY*nlamb*nphi,nY*nlamb*nphi,countP,nA))
	
	allocate(AggH_b(countP,T_sim)) 
	allocate(AggS(countP,T_sim))
	
	print*, 'phi', phi_grid
	print*, 'nphi', nphi
	
	print*, 'LambdaLTV', lambdaLTV
	print*, 'LambdaPTI', lambdaPTI
	print*, 'xi', xi
	print*, 'kappam', kappam
	print*, 'r_lend_grid', r_lend_grid
	print*, 'r_borrow_grid', r_borrow_grid
	print*, 'Y2, pi', Y2, agg_theta
	
	print*, 'p', p  

    !**************************************************
    !  Calling Aggregate Income draws and coefficients
    !**************************************************   
	
    call sub_initial_Y(Y,FYY,YT,AggT_loc,YT_loc,LambT_loc) 
	
    call sub_coeffs(a_raw)

    init_loc=burn
  
    do iY=1,nY*nlamb*nphi
        do iYY=1,nY*nlamb*nphi
            do it=burn+1,T_sim-1
                if (AggT_loc(it)==iY .and. AggT_loc(it+1)==iYY) then
                    init_loc(iY,iYY)=it
                    Exit
                end if
            end do
        end do
    end do 
   
    if (FYY(1,1)==1d0) then
    	init_loc=init_loc(2,2)
    end if
	
	  
    
	allocate(age_sim_bor(I,T_sim))
    allocate(ET_loc(I,T_sim))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   Initial Distribution
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This routine calls and updates the initial distribution
		! Sub_initial_dist calls random numbers where the 7 entries after Variance are:
		! (...,Variance, min housing, max housing, min loans, max loans, min savings, max savings, homeownership rate)
		! One can set the min and max values to be the range from which random shocks are drawn.
   
	    age_sim_bor=0
		ET_loc=0
	   if (initial_dist_opt==1) then
            call sub_initial_dist(L,B,H,R,loge,Fee,age_sim_bor,ET_loc,Variance,minval(H),minval(H),0d0,0d0,0d0,0d0,0d0)   !2.7d0,1d0,0d0,0d0) 
        else
			! This saves the dstirbution from previous stored results
            call sub_initial_dist2(age_sim_bor,L_sim_star(:,1),B_sim_star(:,1),H_sim_star(:,1),rent_buy(:,1),ET_loc,Fee,loge) 
         end if    
print*, 'Done with grids and initial distribution'
        print*, 'a_raw1', a_raw(nY*nLamb*nphi,nY*nLamb*nphi,1)
        print*, 'a_raw2', a_raw(nY*nLamb*nphi,nY*nLamb*nphi,2)
		print*, 'Y2', Y(1),Y(2), Y


    !**************************************************
    ! Begin KS iteration here
    !**************************************************
     iter=0
    diff_solve=1d0
    call CPU_time (t1)
    tp1=omp_get_wtime()
    do while (diff_solve>tol .and. iter<iter_total )
		
        iter=iter+1
        
        ! Need to initialize distributions each iteration
		if (allocated(rental_sim)) deallocate(rental_sim)
		if (allocated(rent_buy)) deallocate(rent_buy)
		if (allocated(max_loc_sim_star)) deallocate(max_loc_sim_star)
		if (allocated(C_sim_star)) deallocate(C_sim_star)
		if (allocated(income_star)) deallocate(income_star)
		if (allocated(L_sim_star)) deallocate(L_sim_star)
		if (allocated(B_sim_star)) deallocate(B_sim_star)
		if (allocated(H_sim_star)) deallocate(H_sim_star)
		if (allocated(tax_star)) deallocate(tax_star)
		if (allocated(employment_star)) deallocate(employment_star)

       
     
       
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       VFI
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call CPU_time (t11)
        tp11=omp_get_wtime()
		
		! Setting up pprime grid with and without learning 
		
		do iY=1,nY*nlamb*nphi
			do iYY=1,nY*nlamb*nphi
				if (nA>1) then
					do iA1=1,nA		
						if (gain_param==0d0) then 
							pprime(iY,iYY,:,iA1)=exp(a_raw(iY,iYY,1)+adj(iA1)+(a_raw(iY,iYY,2)+adj(iA1)*log(p(:)))*log(p(:)))
						else
							pprime(iY,iYY,:,iA1)=exp(a_raw(iY,iYY,1)+gain_param*adj(iA1)+(a_raw(iY,iYY,2)+gain_param*adj(iA1)*log(p(:)))*log(p(:)))
						end if 
					end do 
				else
					pprime(iY,iYY,:,1)=exp(a_raw(iY,iYY,1)+a_raw(iY,iYY,2)*log(p(:)))
						!exp(a_raw(iY+2*(nlamb-1),iYY+2*(nlamb-1),1,ALOC)+a_raw(iY+2*(nlamb-1),iYY+2*(nlamb-1),2,aloc)*log(p(iP)))
				end if 
				!	print*,iY,iYY,a_raw(iY,iYY,:),'pprime',pprime(iY,iYY,:,1) !a_raw(iY+2*(nlamb-1),iYY+2*(nlamb-1),:,Aloc) !pprime(iY,iYY,iP,iA1)
			
			end do
		end do
		
		
		
		
		if(.not. allocated(c_b))   allocate(c_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
		if(.not. allocated(HH_b))  allocate(HH_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
		if(.not. allocated(LL_b))  allocate(LL_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
		if(.not. allocated(BB_b))  allocate(BB_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
		if(.not. allocated(V_b))   allocate(V_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
		if(.not. allocated(rental)) allocate(rental(ny*nlamb*nphi,countP,nA,nA))
		if(.not. allocated(pprimein)) allocate(pprimein(ny*nlamb*nphi,ny*nlamb*nphi,countP,2,nA,nA))
		if(.not. allocated(ilo_pprime)) allocate(ilo_pprime(2,ny*nlamb*nphi,ny*nlamb*nphi,countP,2,nA,nA))
		if(.not. allocated(theta_pprime)) allocate(theta_pprime(2,ny*nlamb*nphi,ny*nlamb*nphi,countP,2,nA,nA))
		!if(.not. allocated(W_rent) allocate(W_rent(nB,nY*nlamb*nE*nPhi,countP,T,nA,nA,2))
		if(.not. allocated(W_own)) allocate(W_own(nL,nH,nB,nY*nlamb*nE*nPhi,countP,T,nA,nA,nD))
		if(.not. allocated(q)) allocate(q(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA))


		! Setting up Pprime grids above to enter VFI either with or without learning
			
		do iA1=1,nA
			do iA2=1,nA
				do iP=1,countP
					do iAgg=1,nY*nlamb*nphi
						do iYY=1,nY*nlamb*nphi
						
							pprimein(iAgg,iYY,iP,:,iA1,iA2)=(/pprime(iAgg,iYY,iP,iA1),pprime(iAgg,iYY,iP,iA2)/)

							call base_fun_noextrap(p,pprimein(iAgg,iYY,iP,1,iA1,iA2),ilo_pprime(:,iAgg,iYY,iP,1,iA1,iA2),theta_pprime(:,iAgg,iYY,iP,1,iA1,iA2))

							call base_fun_noextrap(p,pprimein(iAgg,iYY,iP,2,iA1,iA2),ilo_pprime(:,iAgg,iYY,iP,2,iA1,iA2),theta_pprime(:,iAgg,iYY,iP,2,iA1,iA2))
						end do
							
							iPhi=(1-perfect_corr)*ceiling(real(iAgg)/(real(nlamb)*real(nY)))+perfect_corr
							ilamb=(1-perfect_corr)*ceiling(real(iAgg)/(real(nY)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
										
							if ((iAgg==1 .or. iAgg==2) .and. no_Mit==0) then
								next=1
								nextplus=2
							else if ((iAgg==3 .or. iAgg==4) .and. no_Mit==0) then
								next=3
								nextplus=4
							
							else if ((iAgg==5 .or. iAgg==6) .and. no_Mit==0) then
									next=5
									nextplus=6
								end if 

							if (nA2==1 .and. nA>1 ) then
					
								rental(iAgg,iP,iA1,iA2)= max(1.0d-4,rental_grid(ilamb)+p(iP)*(1.0d0+tauh)-(1.0d0-delta)*dot_product((/pprime(iAgg,next,iP,iA1),pprime(iAgg,nextplus,iP,iA2)/),FYY(iAgg,next:nextplus))/(1.0d0+r_lend_grid(ilamb)))

					
							else
								rental(iAgg,iP,iA1,iA2)= max(1.0d-4,rental_grid(ilamb)+p(iP)*(1.0d0+tauh)-(1.0d0-delta)*dot_product(pprime(iAgg,:,iP,iA1),FYY(iAgg,:))/(1.0d0+r_lend_grid(ilamb)))
					
							end if 		
							
						print*,iA1,iA2,iP, iAgg,ilamb,iphi,rental_cost, rental_grid, rental(iAgg,iP,iA1,iA2),rental_grid(ilamb),rental_grid(ilamb)+p(iP)*(1.0d0+tauh)-(1.0d0-delta)*dot_product((/pprime(iAgg,next,iP,iA1),pprime(iAgg,nextplus,iP,iA2)/),FYY(iAgg,next:nextplus))/(1.0d0+r_lend_grid(ilamb)),pprime(iAgg,:,iP,iA1), pprime(iAgg,:,iP,iA2), p(iP)-pprime(iAgg,1,iP,iA1),p(iP)-pprime(iAgg,2,iP,iA1),p(iP)-pprime(iAgg,1,iP,iA1),p(iP)-pprime(iAgg,2,iP,iA2) !,FYY(iAgg,:) !'pprimein',pprimein(iAgg,iYY,iP,:,iA1,iA2),'theta_pprime',theta_pprime(:,iP,1,iA1,iA2),theta_pprime(:,iP,2,iA1,iA2)
						end do	
					end do 

				end do
			end do 
	
		print*, 'R: main 1',R
		call sub_VFI(a_raw,nA,nA2,p,countP,FYY,Fee,L,H,Y,chi,e,loge,R,c_b,LL_b,BB_b,HH_b,V_b,pprimein,theta_pprime,ilo_pprime,Trans,rental,W_own,q) !,EEEq,EEEp) 
		
        call CPU_time (t22)
		tp22=omp_get_wtime()

        print*, 'VFI time, minutes:', (t22-t11)/60d0, (tp22-tp11)/60d0
        dimensions=(/nL,nH,nB,nY*nlamb,nE,countP,T,nA,burn,I,T_sim/)
		iDprint=1

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! interpolating value function
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		print*, 'Allocating Interpolated Value functions for simulation'
		call CPU_time (t11)
		tp11=omp_get_wtime()
		! Interpollating sequential does not lead to a segmentation fault due to low memory
		deallocate(W_own)
		
		
		!OMP Parallel Default(Shared)
		!$OMP parallel sections	 Default(Shared) num_threads(3)
!		private(iA1,iA2,iJ,iP,iExo,iBB,iH,iLL,inc,liq,iE,iAgg,iPhi,ilamb,iY,iD,iDD) 
		print*, 'Starting interpolation'
		!$omp section
		print*, 'start V'
		if(.not. allocated(V_sim_interp)) allocate(V_sim_interp(nLsim,nH,nBsim,nY*nlamb*nE*nPhi,countP,T,nA,nA,nD+2))	
		V_sim_interp(:,:,:,:,:,:,:,:,1:4)=death
		V_sim_interp(:,:,:,:,:,:,:,:,6:7)=death
		V_sim_interp(:,:,:,:,:,:,:,:,5)=deathdefault
		print*, 'starting V interp'
		call sub_function_interp(V_sim_interp,V_b,B,L,Lsim,Bsim,bs,Y,chi,loge,P,H,nA,nA2)
		print*, 'end V'
		deallocate(V_b)
		
		print*, 'start L'
		if(.not. allocated(L_interp)) allocate(L_interp(nLsim,nH,nBsim,nY*nlamb*nE*nPhi,countP,T,nA,nA,nD+2))
		L_interp=0d0
		call sub_function_interp(L_interp,LL_b,B,L,Lsim,Bsim,bs,Y,chi,loge,P,H,nA,nA2)
		deallocate(LL_b)
		print*, 'end L'
		

		
		!$omp section
		print*, 'start H'
		if(.not. allocated(H_interp)) allocate(H_interp(nLsim,nH,nBsim,nY*nlamb*nE*nPhi,countP,T,nA,nA,nD+2))
		H_interp=0d0
		call sub_function_interp(H_interp,HH_b,B,L,Lsim,Bsim,bs,Y,chi,loge,P,H,nA,nA2)
		deallocate(HH_b)
		print*, 'end H'
		
		print*, 'start B'
		if(.not. allocated(B_interp)) allocate(B_interp(nLsim,nH,nBsim,nY*nlamb*nE*nPhi,countP,T,nA,nA,nD+2))
		B_interp=0d0
		call sub_function_interp(B_interp,BB_B,B,L,Lsim,Bsim,bs,Y,chi,loge,P,H,nA,nA2)
		deallocate(BB_B)
		print*, 'end B'
		
		
			
		
		!$omp section
				
		
		print*, 'start C'
		
		if(.not. allocated(C_interp)) allocate(C_interp(nLsim,nH,nBsim,nY*nlamb*nE*nPhi,countP,T,nA,nA,nD+2))
		C_interp=cmin
		call sub_function_interp(C_interp,C_b,B,L,Lsim,Bsim,bs,Y,chi,loge,P,H,nA,nA2)
		deallocate(C_b)
		print*, 'end C'
		
		
		print*, 'start q'
		if(.not. allocated(q_interp)) allocate(q_interp(nLsim,nH,nBsim,nY*nlamb*nE*nPhi,countP,T,nA,nA))
		q_interp=1d0
		!$OMP parallel DO collapse(6) default(shared)	 private(iA1,iA2,iJ,iP,iExo,iBB,iH,iLL,inc,liq,iE,iAgg,iPhi,ilamb,iY) num_threads(5) if (NA>1)
		do iA2=1,nA !*nA
			do iA1=1,nA
				do iJ=1,T
					do iP=1,nP    
						do iBB=1,nBsim 					
							do iExo=1,nY*nlamb*nE*nPhi
								iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
								iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
								iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
								ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
								iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr
							
								do iH=1,nH
									do iLL=1,nLsim
										q_interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2)=interp2qnoextrap(L,B,q(:,iH,:,iExo,iP,iJ,iA1,iA2),Lsim(iLL),Bsim(iBB))
									end do
								end do
							end do
						end do
					end do
				end do
			end do
		end do 
		!$OMP END PARALLEL Do					
		deallocate(q)
		print*, 'end C'
		
	
		!$OMP END parallel sections
		!OMP END PARALLEL
		
		print*, 'initializing arrays'
		!call CPU_time (t11)
		!tp11=omp_get_wtime()	
		!V_sim_interp=death
		!print*, 'V allocated'
		!C_interp=cmin
		!print*, 'C allocated'
		!B_interp=0d0
		!print*, 'B allocated'
		!L_interp=0d0
		!print*, 'L allocated'
		!H_interp=R(1)
		!print*, 'H allocated'
		!call CPU_time (t22)
		!tp22=omp_get_wtime()
		 !print*, 'initialization, minutes:', (t22-t11)/60d0, (tp22-tp11)/60d0
		 
		 
		
	

		

		
		
		call CPU_time (t22)
		tp22=omp_get_wtime()
		 print*, 'interpolation, minutes:', (t22-t11)/60d0, (tp22-tp11)/60d0
		
		
		
		


	allocate(rental_sim(1:T_sim-1))
	allocate(rent_buy(I,T_sim))
    allocate(max_loc_sim_star(I,T_sim))
	allocate(C_sim_star(I,T_sim))
	allocate(income_star(I,T_sim))
	allocate(L_sim_star(I,T_sim))
	allocate(B_sim_star(I,T_sim))
	allocate(H_sim_star(I,T_sim))
	allocate(tax_star(I,T_sim))
	allocate(employment_star(I,T_sim))
       
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !          Simulation and KS Regression
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print*, 'starting simulation'
		       
		 call CPU_time (tsim1)
		 tpsim1=omp_get_wtime()
		! Three types of simulation type_sim=(0,1,3)=(Steady State, Market Clearing, Learning)
		! Learning runs off of stored coefficients. I would recommend running market clearing till convergence and then learning
        if (type_sim==0) then
			print*, 'Simulating steady state'
		
			if(type_sim==0) allocate(rent_buy_steady(I,countP,T_sim))
			allocate(L_sim_steady(I,countP,T_sim))
			allocate(H_sim_steady(I,countP,T_sim))
			allocate(B_sim_steady(I,countP,T_sim))
			allocate(C_sim_steady(I,countP,T_sim))
			allocate(Tax_sim_steady(I,countP,T_sim))
			allocate(max_loc_b_steady(I,countP,T_sim-1))
		! p and AggH_s are allocated inside sub_grids_fixed
	
			do iP=1,nP
				L_sim_steady(:,iP,1)=L_sim_star(:,1)
				B_sim_steady(:,iP,1)=B_sim_star(:,1)
				H_sim_steady(:,iP,1)=H_sim_star(:,1)
				rent_buy_steady(:,iP,1)=rent_buy(:,1)
			end do
			
			call sub_simulate_steady(age_sim_bor,L_sim_steady,B_sim_steady,H_sim_steady,rent_buy_steady,C_sim_steady,tax_sim_steady,&
			Lsim,Bsim,H,R,Y,chi,L_interp(:,:,:,:,:,:,Aloc,Aloc,:),B_interp(:,:,:,:,:,:,Aloc,Aloc,:),H_interp(:,:,:,:,:,:,Aloc,Aloc,:),c_interp(:,:,:,:,:,:,Aloc,Aloc,:),AggT_loc,YT_loc,LambT_loc,loge,ET_loc,AggH_b,AggH_s,p,countp,pstar,Lstar_b,Bstar_b,Hstar_b,&
			Cstar_b,tax_b,housing_tax_b,income_star,income_b,ret_b,V_sim_interp(:,:,:,:,:,:,Aloc,Aloc,:),employment_b,employment_star,max_loc_b_steady)

    

			call sub_KS_regression(AggT_loc,preg,pstar(1:T_sim-1),counts,a_new,R2p,R2p_all,.True.)
			 
			!OMP Parallel Default(Shared) private(iI,it)
			!OMP DO collapse(2) 
			
			!Interpolating objects at the price where excess demand equals 0 to make these arrays consistent with simulations from market clearing
			! Without aggregate risk, the steady state and market clearing simulations should be similar
			do it=1,T_sim-1
				do iI=1,I
					L_sim_star(iI,it)=interp1q(p,L_sim_steady(iI,:,it),pstar(it))
					B_sim_star(iI,it)=interp1q(p,B_sim_steady(iI,:,it),pstar(it))
					H_sim_star(iI,it)=interp1q(p,H_sim_steady(iI,:,it),pstar(it))
					C_sim_star(iI,it)=interp1q(p,C_sim_steady(iI,:,it),pstar(it))
					tax_star(iI,it)=interp1q(p,tax_sim_steady(iI,:,it),pstar(it))

					rent_buy(iI,it)=rent_buy_steady(iI,bsearch(pstar(it),p),it)
					max_loc_sim_star(iI,it)=max_loc_b_steady(iI,bsearch(pstar(it),p),it)

				end do
			end do
			!OMP End Do
			!OMP End parallel

		else if (type_sim==1) then
		
			print*, 'Simulating No learning',maxval(Bsim),minval(Bsim)
				! In learning we set this to learn so that it runs market clearing up to learn and then switches

            call sub_simulate_bor(age_sim_bor,rent_buy,Lsim,Bsim,Bs,H,R,Y,chi,V_sim_interp(:,:,:,:,:,:,Aloc,Aloc,:),L_interp(:,:,:,:,:,:,Aloc,Aloc,:),B_interp(:,:,:,:,:,:,Aloc,Aloc,:),H_interp(:,:,:,:,:,:,Aloc,Aloc,:),c_interp(:,:,:,:,:,:,Aloc,Aloc,:),AggT_loc,YT_loc,LambT_loc,loge,ET_loc,&
            AggH_b,AggS,AggH_s,p,countP,countD,sim_end,pstar,Lstar_b,Bstar_b,Hstar_b,Hstar_s,Cstar_b,qstar,income_b,C_sim_star,&
            income_star,L_sim_star,B_sim_star,H_sim_star,tax_b,tax_star,housing_tax_b,employment_b,employment_star,ret_b,max_loc_sim_star,max_loc_b,pprime(:,:,:,Aloc),homeownership,q_interp(:,:,:,:,:,:,Aloc,Aloc))

				
				

			print*, 'Done with Market clearing'
			!print*, pstar
            call sub_KS_regression(AggT_loc,preg,pstar(1:T_sim-1),counts,a_new,R2p,R2p_all,.True.)
            print*, 'made it through KS' 
        else
	
            print*, 'Simulate Learning'
		! Simulates using the market clearing step from t=1 to learn-1 
		
            call sub_simulate_bor(age_sim_bor(:,1:sim_end),rent_buy(:,1:sim_end),Lsim,Bsim,Bs,H,R,Y,chi,V_sim_interp(:,:,:,:,:,:,Aloc,Aloc,:),L_interp(:,:,:,:,:,:,Aloc,Aloc,:),B_interp(:,:,:,:,:,:,Aloc,Aloc,:),&
			H_interp(:,:,:,:,:,:,Aloc,Aloc,:),c_interp(:,:,:,:,:,:,Aloc,Aloc,:),AggT_loc,YT_loc,LambT_loc,loge,ET_loc,&
            AggH_b(:,1:sim_end),AggS(:,1:sim_end),AggH_s,p,countP,countD,sim_end,pstar(1:sim_end-1),Lstar_b(1:sim_end),Bstar_b(1:sim_end),Hstar_b(1:sim_end),&
			Hstar_s(1:sim_end),Cstar_b(1:sim_end-1),qstar(1:sim_end),income_b(1:sim_end-1),C_sim_star(:,1:sim_end-1),income_star(:,1:sim_end-1),L_sim_star(:,1:sim_end),&
			B_sim_star(:,1:sim_end),H_sim_star(:,1:sim_end),tax_b(1:sim_end-1),tax_star(:,1:sim_end-1),housing_tax_b(1:sim_end-1),employment_b(1:sim_end-1),employment_star(:,1:sim_end-1),&
			ret_b(1:sim_end-1),max_loc_sim_star(:,1:sim_end-1),max_loc_b(1:sim_end-1,:),pprime(:,:,:,Aloc),homeownership,q_interp(:,:,:,:,:,:,Aloc,Aloc))
		
			count=0
			preg=0d0
			countZ=0
			do it=1,sim_end-1
				if (it>burn) then
					Yreveal=AggT_loc(it+1)
					countZ(AggT_loc(it))=countZ(AggT_loc(it))+1
					count(AggT_loc(it),Yreveal)=count(AggT_loc(it),Yreveal)+1
					
					preg(AggT_loc(it),Yreveal,count(AggT_loc(it),Yreveal))=pstar(it)
					pregZ(AggT_loc(it),countZ(AggT_loc(it)))=pstar(it)
					!print*, it, AggT_loc(it),count(AggT_loc(it),Yreveal), preg(AggT_loc(it),Yreveal,count(AggT_loc(it),Yreveal))
				end if
			!print*, H_sim_star(1:10,it+1)
			!print*, B_sim_star(1:10,it+1)
			!print*, L_sim_star(1:10,it+1)
			!print*, AggS(:,it+1)
			end do
			print*, 'count', count(3,1)
			print*, 'pprime', pprime(1,1,1,:)
			print*, 'pprime',pprime(1,2,1,:)
			
			
				call sub_simulate_learning(age_sim_bor,rent_buy,Lsim,Bsim,Bs,H,R,Y,chi,V_sim_interp,L_interp,B_interp,H_interp,c_interp,AggT_loc,YT_loc,LambT_loc,loge,ET_loc,&
				AggH_b,AggS,AggH_s,p,countP,pstar,Lstar_b,Bstar_b,Hstar_b,Hstar_s,Cstar_b,income_b,C_sim_star,income_star,L_sim_star,&
				B_sim_star,H_sim_star,tax_b,tax_star,housing_tax_b,employment_b,employment_star,ret_b,max_loc_sim_star,max_loc_b,&
				countD,a_out,nA,pprime,sim_end,count,preg,countZ,pregZ,homeownership,pprime_solve,rental_sim,rental(:,:,:,:),a_out_all,q_interp)
				R2p_all=0d0
				print*, 'pprime', pprime(:,:,:,:)
				print*, 'pprime',pprime(:,:,:,:)
			
        end if
        call CPU_time (tsim2) 
		tpsim2=omp_get_wtime()	
        print*, 'Simulation time, minutes (raw,parallel):', (tsim2-tsim1)/60d0, (tpsim2-tpsim1)/60d0
        counts1=0
		
		deallocate(V_sim_interp)
		deallocate(C_interp)
		deallocate(H_interp)
		deallocate(L_interp)
		deallocate(B_interp)
        

		   Aloc=(nA+1)/2  
		   		do iY=1,nY*nlamb*nphi
					do iYY=1,nY*nlamb*nphi
						print*, 'a_raw: old a0', iY,iYY, a_raw(iY,iYY,1)
						print*, 'a_raw: old a1', iY,iYY,a_raw(iY,iYY,2)
						print*, 'a_new: new coeefs a0', iY,iYY, a_new(iY,iYY,1)
						print*, 'a_new: new coeefs a0', iY,iYY, a_new(iY,iYY,2)
					end do
				end do 
		!From KMV, old version commented out below
        diff_solve=maxval(abs(a_new-a_raw)) !maxval(abs(a_new((nlamb-1)*nY+1:(nlamb-1)*nY+2,(nlamb-1)*nY+1:(nlamb-1)*nY+2,:)-a_raw((nlamb-1)*nY+1:(nlamb-1)*nY+2,(nlamb-1)*nY+1:(nlamb-1)*nY+2,:)))
		!maxval(abs(a_new-a_raw(:,:,:,Aloc)))/maxval(abs(a_raw(:,:,:,Aloc)))
        diff_total(iter)=diff_solve
        
        if (type_sim<2) then
            
            a_raw(:,:,1)=damp*a_raw(:,:,1)+a_new(:,:,1)*(1d0-damp)
            a_raw(:,:,2)=damp*a_raw(:,:,2)+a_new(:,:,2)*(1d0-damp)
        end if
            print*, 'a_raw updated a0', a_raw(:,:,1)
			  print*, 'a_raw updated a1', a_raw(:,:,2)
            !print*, 'a', a
   
        


        print*,'__________________________________'
        print*, 'iter2, iter', iter2, iter
            print*, 'diff solve', diff_solve !, maxloc(abs(a_new-a_raw(:,:,:,Aloc)))
        print*, 'diff total', diff_total(1:iter)
    print*, 'KMV R^2', R2p_all
		

		do iY=1,nY*nlamb*nphi
			do iYY=1,nY*nlamb*nphi
				if (a_new(iY,iYY,2)>1d0) then
					print*,'index| raw a| mean | R2' ,iY,iYY,'	 |', a_new(iY,iYY,:),'    |',exp(a_new(iY,iYY,1)/(1d0-a_new(iY,iYY,2))), '    | ', R2p(iY,iYY)
				else
					print*,'index| raw a| mean | R2' ,iY,iYY,'	 |', a_new(iY,iYY,:), '    |',0d0, '    | ', R2p(iY,iYY)
				end if 
			end do
		end do
	end do
	if(.not. allocated(rental_sim)) print*, 'rental_sim'
    if(.not. allocated(age_sim_bor)) print*,  'age_sim_bor'
    if(.not. allocated(rent_buy)) print*, 'three'
    if(.not. allocated(ET_loc)) print*, 'four'
    if(.not. allocated(max_loc_sim_star)) print*, 'five'
	if(.not. allocated(FYY)) print*, 'six'
	if(.not. allocated(C_sim_star)) print*, 'seven'
	if(.not. allocated(income_star)) print*, 'eight'
	if(.not. allocated(L_sim_star)) print*, 'nine'
	if(.not. allocated(B_sim_star)) print*, 'ten'
	if(.not. allocated(H_sim_star)) print*, 'eleven'
	if(.not. allocated(tax_star)) print*, 'twelve'
	if(.not. allocated(employment_star)) print*, 'thirteen'
	if(.not. allocated(a_raw)) print*, 'fourteen'
	if(.not. allocated(pprime)) print*, 'fifteen'
    if(.not. allocated(c_b)) print*, 'sixteen'
    if(.not. allocated(HH_b)) print*, 'seventeen'
    if(.not. allocated(LL_b)) print*, 'eighteen'
    if(.not. allocated(BB_b)) print*, 'nineteen'
    if(.not. allocated(V_b)) print*, 'twenty'
	if(.not. allocated(rental)) print*, 'twentone'


	
	if(.not. allocated(AggH_b))  print*, 'twenty two'
	if(.not. allocated(AggS)) print*, 'twenty three'
	if(.not. 		allocated(P)) print*, 'twenty four'
	if(.not. 	allocated(AggH_s)) print*, 'twenty five'
	if(.not. 	allocated(Trans))  print*, 'twent six'





       count_med2=0
        count_bequest=0
        count_bequest_tot=0
     
    if(allocated(c_b)) deallocate(c_b)
    if(allocated(HH_b)) deallocate(HH_b)
    if(allocated(LL_b)) deallocate(LL_b)
	if(allocated(BB_b)) deallocate(BB_b)

	
	
     allocate(NW_med1(5000,sim_end-burn))
        allocate(NW_med2(5000,sim_end-burn))
		allocate(NW_med1_print(sim_end-burn))
        allocate(NW_med2_print(sim_end-burn))
		allocate(H_W(I,sim_end-burn))
		allocate(H_NW(I,sim_end-burn))
		allocate(H_NW_10(sim_end-burn))
		allocate(H_NW_50(sim_end-burn))
		allocate(H_NW_90(sim_end-burn))
		allocate(H_W_own(I,sim_end-burn))
		allocate(H_NW_own(I,sim_end-burn))
		allocate(earn(2,sim_end-burn))
		allocate(Agg_nw(I,sim_end-burn))
		allocate(med_ratio(I))
		allocate(NW_to_inc(I,sim_end-burn))
		allocate(LTV(I,sim_end-burn))
		allocate(House_value_earn(I,sim_end-burn))
		allocate(With_mortgage(sim_end-burn))
		allocate(With_heloc(sim_end-burn))
		
		
		allocate(NW_20(sim_end-burn)) 
		allocate(NW_40(sim_end-burn))
		allocate(NW_60(sim_end-burn)) 	
		allocate(NW_50(sim_end-burn))
		allocate(NW_80(sim_end-burn))
		allocate(NW_90(sim_end-burn))
		allocate(NW_99(sim_end-burn))    		

		allocate(count_med1(sim_end-burn))
		allocate(count_med2(sim_end-burn))
		allocate(count_buy(sim_end-burn))
		allocate(count_rent(sim_end-burn))
		allocate(count_hy(sim_end-burn))
		allocate(count_LTV(sim_end-burn))
		
		allocate(H_W_50(sim_end-burn))
		allocate(H_NW_50_bequest(sim_end-burn))
		allocate(NW_to_inc_50(sim_end-burn))
		allocate(LTV_10(sim_end-burn))
		allocate(LTV_50(sim_end-burn))
		allocate(LTV_90(sim_end-burn))
		allocate(afford_p10(T_sim-1))
		allocate(afford_p50(T_sim-1))
		allocate(LTV_mean(T_Sim-1))
		
		allocate(bottom_NW(sim_end-burn))
		allocate(mid_NW(sim_end-burn))
		allocate(top_NW(sim_end-burn))
		allocate(top10_NW(sim_end-burn))
		allocate(top1_NW(sim_end-burn))
		
		allocate(count_bequest_tot(sim_end-burn))
		allocate(count_bequest(sim_end-burn))
		allocate(Agg_inc(sim_end-burn))
		
		allocate(House_value_earn_10(sim_end-burn))
		allocate(House_value_earn_50(sim_end-burn))
		allocate(house_value_earn_90(sim_end-burn))
		
		allocate(avg_buy_rent(2,sim_end-burn))
		allocate(under_35(2,sim_end-burn))
		
		allocate(counts_h35(nH,sim_end-burn))
		allocate(counts_h49(nH,sim_end-burn))
		allocate(counts_h64(nH,sim_end-burn))
		allocate(counts_h80(nH,sim_end-burn))
		
		allocate(counts_h(nH+1,sim_end-burn))
		allocate(counts_r(nH+1,sim_end-burn))


	
	

		
        NW_med1=0d0
		NW_med2=0d0
		Agg_nw=0d0
		med_ratio=0d0
		With_mortgage=0
		With_heloc=0
		Earn=0d0
		
		counts_h35=0d0
		counts_h49=0d0
		counts_h64=0d0
		counts_h80=0d0
		counts_h=0d0
		counts_r=0d0
		
		
					count_med1=0
			count_med2=0
			count_buy=0
			count_hy=0
			count_LTV=0
		
		H_NW_own=0d0
		H_W_own=0d0
				           print*, 'a_raw updated a0 start', a_raw(:,:,1)
			  print*, 'a_raw updated a1 start', a_raw(:,:,2)
		if (diff_solve<tol) then
			print*, 'Converged'
		end if 
		!OMP Parallel Default(Shared) private(it,iI,rloc,thetar,hloc,thetah)
		!OMP DO collapse(2) reduction(+:count_med1) reduction(+:count_med2) reduction(+:Earn) reduction(+:count_buy) reduction(+:count_hy) reduction(+:count_LTV) reduction(+:counts_r) reduction(+:counts_h) reduction(+:counts_h35) reduction(+:counts_h49) reduction(+:counts_h64) reduction(+:counts_h80)
		do it=burn,sim_end-1
			print*, AggT_loc(it),pstar(it)
			do iI=1,I
				if (age_sim_bor(iI,it)==28) then      
					count_med1(it-burn+1)=count_med1(it-burn+1)+1
					if (rent_buy(iI,it+1)==1) then
						NW_med1(count_med1(it-burn+1),it-burn+1)=B_sim_star(iI,it+1)
					else 
						NW_med1(count_med1(it-burn+1),it-burn+1)=B_sim_star(iI,it+1)+(pstar(it)*(H_sim_star(iI,it+1))-L_sim_star(iI,it+1)*H_sim_star(iI,it+1))
					end if 
				end if
				if (age_sim_bor(iI,it)==16) then
					count_med2(it-burn+1)=count_med2(it-burn+1)+1

					if (rent_buy(iI,it+1)==1) then
						NW_med2(count_med2(it-burn+1),it-burn+1)=B_sim_star(iI,it+1)
					else

						NW_med2(count_med2(it-burn+1),it-burn+1)=B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1)
					end if 
				end if
				
				if(age_sim_bor(iI,it)<T) then 
					if(rent_buy(iI,it+1)==1) then 
						H_NW(iI,it-burn+1)=0d0
						Agg_NW(iI,it-burn+1)=(B_sim_star(iI,it+1))
						H_W(iI,it-burn+1)=B_sim_star(iI,it+1)
						Earn(1,it-burn+1)=Earn(1,it-burn+1)+income_star(iI,it)
						NW_to_inc(iI,it-burn+1)=Agg_NW(iI,it-burn+1)/income_star(iI,it) !Earn(1,iI,it-burn+1)
						count_rent(it-burn+1)=count_rent(it-burn+1)+1
					else
					
						H_NW(iI,it-burn+1)=(pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1))/(B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1))
						Agg_NW(iI,it-burn+1)=(B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1))
						H_W(iI,it-burn+1)=B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it+1)
						Earn(2,it-burn+1)=Earn(2,it-burn+1)+income_star(iI,it)
						NW_to_inc(iI,it-burn+1)=Agg_NW(iI,it-burn+1)/income_star(iI,it) !/Earn(2,iI,it-burn+1)
						count_buy(it-burn+1)=count_buy(it-burn+1)+1
						if(age_sim_bor(iI,it)<Jret) then
							count_hy(it-burn+1)=count_hy(it-burn+1)+1
							House_value_earn(count_hy(it-burn+1),it-burn+1)=pstar(it)*H_sim_star(iI,it+1)/(income_star(iI,it)/2d0)
						end if 
						If (L_sim_star(iI,it+1)>0d0) then 
							count_LTV(it-burn+1)=count_LTV(it-burn+1)+1
							LTV(count_LTV(it-burn+1),it-burn+1)=L_sim_star(iI,it+1)/pstar(it)
						end if 
						H_NW_own(count_buy(it-burn+1),it-burn+1)=(pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1))/(B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1))	
						H_W_own(count_buy(it-burn+1),it-burn+1)=B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1)
					

						if (L_sim_star(iI,it+1)>0) then
							With_mortgage(it-burn+1)=With_mortgage(it-burn+1)+1
						end if
						if (B_sim_star(iI,it+1)<0) then
							With_heloc(it-burn+1)=With_heloc(it-burn+1)+1
						end if 
					end if 
				else
					if(rent_buy(iI,it)==1) then 
						H_NW(iI,it-burn+1)=0d0
						Agg_NW(iI,it-burn+1)=(B_sim_star(iI,it))
						H_W(iI,it-burn+1)=B_sim_star(iI,it)
						Earn(1,it-burn+1)=Earn(1,it-burn+1)+income_star(iI,it)
						NW_to_inc(iI,it-burn+1)=Agg_NW(iI,it-burn+1)/income_star(iI,it) !Earn(1,iI,it-burn+1)
						count_rent(it-burn+1)=count_rent(it-burn+1)+1
					else
					
						H_NW(iI,it-burn+1)=(pstar(it)*H_sim_star(iI,it)-L_sim_star(iI,it)*H_sim_star(iI,it))/(B_sim_star(iI,it+1)+pstar(it)*H_sim_star(iI,it)-L_sim_star(iI,it)*H_sim_star(iI,it))
						Agg_NW(iI,it-burn+1)=(B_sim_star(iI,it)+pstar(it)*H_sim_star(iI,it)-L_sim_star(iI,it)*H_sim_star(iI,it+1))
						H_W(iI,it-burn+1)=B_sim_star(iI,it)+pstar(it)*H_sim_star(iI,it)
						Earn(2,it-burn+1)=Earn(2,it-burn+1)+income_star(iI,it)
						NW_to_inc(iI,it-burn+1)=Agg_NW(iI,it-burn+1)/income_star(iI,it) !/Earn(2,iI,it-burn+1)
						count_buy(it-burn+1)=count_buy(it-burn+1)+1
						if(age_sim_bor(iI,it)<Jret) then
							count_hy(it-burn+1)=count_hy(it-burn+1)+1
							House_value_earn(count_hy(it-burn+1),it-burn+1)=pstar(it)*H_sim_star(iI,it)/(income_star(iI,it)/2d0)
						end if 
						If (L_sim_star(iI,it+1)>0d0) then 
							count_LTV(it-burn+1)=count_LTV(it-burn+1)+1
							LTV(count_LTV(it-burn+1),it-burn+1)=L_sim_star(iI,it+1)/pstar(it)
						end if 
						H_NW_own(count_buy(it-burn+1),it-burn+1)=(pstar(it)*H_sim_star(iI,it)-L_sim_star(iI,it)*H_sim_star(iI,it))/(B_sim_star(iI,it)+pstar(it)*H_sim_star(iI,it)-L_sim_star(iI,it)*H_sim_star(iI,it))	
						H_W_own(count_buy(it-burn+1),it-burn+1)=B_sim_star(iI,it)+pstar(it)*H_sim_star(iI,it)-L_sim_star(iI,it)*H_sim_star(iI,it)
					

						if (L_sim_star(iI,it)>0) then
							With_mortgage(it-burn+1)=With_mortgage(it-burn+1)+1
						end if
						if (B_sim_star(iI,it)<0) then
							With_heloc(it-burn+1)=With_heloc(it-burn+1)+1
						end if 
					end if 
				end if 
				
							! Looking at housign choice so next period renters	
			if (rent_buy(iI,it+1)==1) then
				rloc=bsearch(h_sim_star(iI,it+1),R)
				thetar=1d0-(h_sim_star(iI,it+1)-R(rloc))/(R(rloc+1)-R(rloc))
				thetar=max(0d0, min(1d0,thetar))
				
				counts_r(rloc,it-burn+1)=counts_r(rloc,it-burn+1)+thetar
				if (rloc<nR) then
					counts_r(rloc+1,it-burn+1)=counts_r(rloc+1,it-burn+1)+(1d0-thetar)
				end if 
				if (age_sim_bor(iI,it)<8) then
					counts_h35(rloc,it-burn+1)=counts_h35(rloc,it-burn+1)+1
				else if (age_sim_bor(iI,it)>=8 .and. age_sim_bor(iI,it)<15 ) then
					counts_h49(rloc,it-burn+1)=counts_h49(rloc,it-burn+1)+1
				else if (age_sim_bor(iI,it)>=15 .and. age_sim_bor(iI,it)<23 ) then
					counts_h64(rloc,it-burn+1)=counts_h64(rloc,it-burn+1)+1
				else if (age_sim_bor(iI,it)>=23 ) then 
					counts_h80(rloc,it-burn+1)=counts_h80(rloc,it-burn+1)+1
				end if 
				
				
			else
				hloc=bsearch(h_sim_star(iI,it+1),H)
				thetah=1d0-(h_sim_star(iI,it+1)-H(hloc))/(H(hloc+1)-H(hloc))
				thetah=max(0d0, min(1d0,thetah))
				
				counts_h(hloc+1,it-burn+1)=counts_h(hloc+1,it-burn+1)+thetah 
				if (hloc<nH) then
					counts_h(hloc+2,it-burn+1)=counts_h(hloc+2,it-burn+1)+(1d0-thetah)
				end if 
				
			
				if (age_sim_bor(iI,it)<8) then
					counts_h35(hloc+1,it-burn+1)=counts_h35(hloc+1,it-burn+1)+1
				else if (age_sim_bor(iI,it)>=8 .and. age_sim_bor(iI,it)<15 ) then
					counts_h49(hloc+1,it-burn+1)=counts_h49(hloc+1,it-burn+1)+1
				else if (age_sim_bor(iI,it)>=15 .and. age_sim_bor(iI,it)<23 ) then
					counts_h64(hloc+1,it-burn+1)=counts_h64(hloc+1,it-burn+1)+1
				else if (age_sim_bor(iI,it)>=23 ) then 
					counts_h80(hloc+1,it-burn+1)=counts_h80(hloc+1,it-burn+1)+1
				end if 

			end if 
	            if ((it==sim_end/4 .or. it==sim_end/2 .or. it==sim_end*3/4) .and. (iI==I/4 .or. iI==I/2 .or. iI==I*3/4)) then
                print*, 'First Calibration Loop (it,iI)', it,iI
            end if
			
			
			
			end do
		end do
		!OMP END DO
		!OMP END PARALLEL
		print*, 'Done with first calibration loop'

		do it=burn,sim_end-1
		counts_h35(:,it-burn+1)=counts_h35(:,it-burn+1)/sum(counts_h35(:,it-burn+1))*100
		counts_h49(:,it-burn+1)=counts_h49(:,it-burn+1)/sum(counts_h49(:,it-burn+1))*100
		counts_h64(:,it-burn+1)=counts_h64(:,it-burn+1)/sum(counts_h64(:,it-burn+1))*100
		counts_h80(:,it-burn+1)=counts_h80(:,it-burn+1)/sum(counts_h80(:,it-burn+1))*100
		counts_h(:,it-burn+1)=counts_h(:,it-burn+1)/sum(counts_h(:,it-burn+1))*100
		counts_r(:,it-burn+1)=counts_r(:,it-burn+1)/sum(counts_r(:,it-burn+1))*100
		end do


		call CPU_time (tsim1) 
		tpsim1=omp_get_wtime()	
		!$OMP Parallel Default(Shared) private(it)
		!$OMP DO 
		do it=1,sim_end-burn

			
			if(count_LTV(it)>0) then 
			LTV_mean(it)=mean(LTV(1:count_LTV(it),it))
			else
			LTV_mean(it)=0d0
			end if 
	
			call quicksort(NW_med1(1:count_med1(it),it)) !INSSOR
			call quicksort(NW_med2(1:count_med2(it),it))
			call quicksort(H_NW(:,it))
			call quicksort(H_NW_own(1:count_buy(it),it))
			call quicksort(H_W(:,it))
			call quicksort(NW_to_inc(:,it))
			call quicksort(LTV(1:count_LTV(it),it))
			Call quicksort(Agg_nW(:,it))
			call quicksort(House_value_earn(1:count_hy(it),it))
			
			!NW_med1_print(it)=indnth(NW_med1(1:count_med1(it),it),ceiling(.5*count_med1(it)))
			NW_med1_print(it)=percentile(NW_med1(1:count_med1(it),it),50)
			NW_med2_print(it)=percentile(NW_med2(1:count_med2(it),it),50)
			H_NW_10(it)=percentile(H_NW_own(1:count_buy(it),it),10)
			H_NW_50(it)=percentile(H_NW_own(1:count_buy(it),it),50)
			H_NW_90(it)=percentile(H_NW_own(1:count_buy(it),it),90)
			H_W_50(it)=percentile(H_W(:,it),50)
			H_NW_50_bequest(it)=percentile(H_NW(:,it),50)
			NW_to_inc_50(it)=percentile(NW_to_inc(:,it),50)
			
			LTV_10(it)=percentile(LTV(1:count_LTV(it),it),10)
			LTV_50(it)=percentile(LTV(1:count_LTV(it),it),50)
			LTV_90(it)=percentile(LTV(1:count_LTV(it),it),90)
			
			
			!print*, 'House value it', it, count_hy(it)
			!print*, House_value_earn(1:count_hy(it),it)

			House_value_earn_10(it)=percentile(House_value_earn(1:count_hy(it),it),10)
			House_value_earn_50(it)=percentile(House_value_earn(1:count_hy(it),it),50)
			House_value_earn_90(it)=percentile(House_value_earn(1:count_hy(it),it),90)
			
			NW_20(it)=percentile(Agg_NW(:,it),20)
			NW_40(it)=percentile(Agg_NW(:,it),40)
			NW_60(it)=percentile(Agg_NW(:,it),60)
			NW_50(it)=percentile(Agg_NW(:,it),50)
			NW_80(it)=percentile(Agg_NW(:,it),80)
			NW_90(it)=percentile(Agg_NW(:,it),90)
			NW_99(it)=percentile(Agg_NW(:,it),99)
			
			bottom_NW(it)=sum(Agg_nw(1:bsearch(NW_20(it),Agg_NW(:,it)),it))/sum(Agg_nW(:,it))
			mid_NW(it)=sum(Agg_nw(bsearch(NW_40(it),Agg_NW(:,it)):bsearch(NW_60(it),Agg_NW(:,it)),it))/sum(Agg_nw(:,it))
			top_NW(it)=sum(Agg_nw(bsearch(NW_80(it),Agg_nW(:,it)):,it))/sum(Agg_nW(:,it))
			top10_NW(it)=sum(Agg_nw(bsearch(NW_90(it),Agg_nW(:,it)):,it))/sum(Agg_nW(:,it)) 
			top1_NW(it)=sum(Agg_nw(bsearch(NW_99(it),Agg_nW(:,it)):,it))/sum(Agg_nW(:,it))
		if ((it==sim_end/4 .or. it==sim_end/2 .or. it==sim_end*3/4) ) then
                print*, 'Sorting Loop it', it
            end if
		end do 
		!$OMP END DO
		!$OMP END PARALLEL
		print*, 'Done with sorting loop'
		call CPU_time (tsim2) 
		tpsim2=omp_get_wtime()	
	print*, 'Percentiles, minutes (raw,parallel):', (tsim2-tsim1)/60d0, (tpsim2-tpsim1)/60d0
	
		avg_buy_rent=0d0
		under_35=0d0
		

		!$OMP Parallel Default(Shared) private(it,iI)
		!$OMP DO collapse(2) reduction(+:count_bequest_tot) reduction(+:count_bequest) reduction(+:avg_buy_rent) reduction(+:under_35) 
		do it=burn,sim_end -1
			do iI=1,I
			
				! Bequests
				if (age_sim_bor(iI,it)==T) then
					if (max_loc_sim_star(iI,it)==2 .and. pstar(it)*h_sim_star(iI,it+1)+b_sim_star(iI,it+1)*(1d0-delta-tauh-cost)>0) then
						count_bequest_tot(it-burn+1)=count_bequest_tot(it-burn+1)+1 
						!if (pstar(burn)*h_sim_star(iI,burn)+b_sim_star(iI,burn)*(1d0-delta-tauh-cost)<H_W_50) then
						if (pstar(it)*h_sim_star(iI,it+1)+b_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1) <H_NW_50_bequest(it-burn+1)) then
							count_bequest(it-burn+1)=count_bequest(it-burn+1)+1
						end if
					else if (max_loc_sim_star(iI,it)/=2 .and. b_sim_star(iI,it+1)>0) then
						count_bequest_tot(it-burn+1)=count_bequest_tot(it-burn+1)+1
						if (pstar(it)*h_sim_star(iI,it+1)+b_sim_star(iI,it+1)-L_sim_star(iI,it+1)*H_sim_star(iI,it+1)<H_NW_50_bequest(it-burn+1)) then
							count_bequest(it-burn+1)=count_bequest(it-burn+1)+1
						end if
					end if 
				end if
				
						! Average buyers/rents
				if (rent_buy(iI,it+1)==1) then
					avg_buy_rent(1,it-burn+1)=avg_buy_rent(1,it-burn+1)+h_sim_star(iI,it+1)
				else
					avg_buy_rent(2,it-burn+1)=avg_buy_rent(2,it-burn+1)+h_sim_star(iI,it+1)
				end if
			
				! Homeownership rate of under 35 years olds
				if (age_sim_bor(iI,it)<8) then
					under_35(1,it-burn+1)=under_35(1,it-burn+1)+rent_buy(iI,it+1)-1
					under_35(2,it-burn+1)=under_35(2,it-burn+1)+1
				end if
				  if ((it==sim_end/4 .or. it==sim_end/2 .or. it==sim_end*3/4) .and. (iI==I/4 .or. iI==I/2 .or. iI==I*3/4)) then
                print*, 'Second Calibration Loop (it,iI)', it,iI
            end if
				
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL


		do it=burn,sim_end-1
					Agg_inc(it-burn+1)=Y(YT_loc(it))*employment_b(it)      
			avg_buy_rent(1,it-burn+1)=avg_buy_rent(1,it-burn+1)/real(sum((-1)*(rent_buy(:,it+1)-2)))
			avg_buy_rent(2,it-burn+1)=avg_buy_rent(2,it-burn+1)/real(sum(rent_buy(:,it+1)-1))
			under_35(1,it-burn+1)=under_35(1,it-burn+1)/under_35(2,it-burn+1)
		end do
		
		!rental_sim=0d0
		!if (nA==1) then
		do it=1,sim_end-1
			do iYY=1,nY*nlamb*nphi
			iPhi=(1-perfect_corr)*ceiling(real(iYY/(real(nlamb)*real(nY))))+perfect_corr
								ilamb=(1-perfect_corr)*ceiling(real(iYY)/(real(nY)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
										
				!print*, 'iY,iYY', interp1q(p,pprime(AggT_loc(it),iYY,:,1),pstar(it))*FYY(AggT_loc(it),iYY),interp1q(p,pprime(AggT_loc(it),iYY,:,1),pstar(it)),FYY(AggT_loc(it),iYY),(1.0d0-delta-tauh)/(1.0d0+r_lend)
				if(it<sim_end) then
			
				rental_sim(it)=interp1q(p,rental(AggT_loc(it),:,Aloc,Aloc),pstar(it))
				!if (it>397 .and. it<48) then
					!print*, AggT_loc(it),iYY, interp1q(p,pprime(AggT_loc(it),iYY,:,Aloc),pstar(it)),FYY(AggT_loc(it),iYY)
				end if 	
				!end if 

			end do
		
			
			!end if 
		end do
		!end if
		
		
		deallocate(LTV)
		allocate(count_LTV2(1:T_sim-1))
		allocate(LTV(I,1:T_sim-1))
		allocate(afford(I,T_sim-1))
		count_LTV2=0
		LTV=0d0
		afford=0d0
		

		
		do it=1,T_sim-1
			do iI=1,I
				if (rent_buy(iI,it)==2) then
				!If (L_sim_star(iI,it+1)>0d0) then 
						count_LTV2(it)=count_LTV2(it)+1
						LTV(count_LTV2(it),it)=L_sim_star(iI,it)/pstar(it)
						afford(iI,it)=(1+r_borrow_grid(LambT_loc(it)))*H_sim_star(iI,it)*L_sim_star(iI,it)/income_star(iI,it)
				else
					
						afford(iI,it)=rental_sim(it)*H_sim_star(iI,it)/income_star(iI,it)
				end if 
				
			end do
			if (it>1) then 
				LTV_mean(it)=mean(LTV(1:count_LTV2(it),it))
			else 
				LTV_mean(it)=0d0
			end if 
			
			if(it>1) then
				foreclosure(it)=real(countD(it,5))/real(sum(rent_buy(:,it)-1))*100
			else
				foreclosure(it)=0d0
			end if 
			call quicksort(afford(:,it))
			afford_p10(it)=percentile(afford(:,it),10)
			afford_p50(it)=percentile(afford(:,it),50)
		end do

	
        
        print*, '*************    Targeted Calibration    **************'

	   
	   
	   ! Note some equilibirium properties that are useful for deriving stats and ratios
	   ! Employment in the construction and goods sectors equals:
	   ! N_h+N_c=\int_{work}\exp(\chi_j+\epsilon_{ij})d \mu_{work}
	   
	   ! Where:
	   ! wage=Theta (wage=Y in the code notation)
	   ! From the zero profit condition of the construction firm p*I_h=wage*N_h 
	   ! <=> p*(Y*N_h)^\alpha(\bar{L})^{1-\alpha}=Y*N_h
	   ! Dividiving both sides by N_h^\alpha
	   ! <=> p*Y^\alpha (\bar{L})^{1-\alpha} = Y*N_h^{1-\alpha}
	   ! Diving both sides by Y
	   ! <=> p*Y^{\alpha-1} (\bar{L})^{1-\alpha} = N_h^{1-\alpha}
	   ! <=> p^{1/(1-\alpha)}*Y^{-1}*\bar{L}=N_hs
		print*, 'mean pstar', mean(pstar(burn:sim_end-1))
		print*, 'Agg. NW/ aggregate labor income, want 5.5, otherwise adjust beta', mean((sum(Agg_nw,dim=1)/real(I))/(YT(burn:sim_end-1)*employment_b(burn:sim_end-1)/2d0)) ! the denominator was divided by 2 to make it annual, but I don't think this is quite right        
		 print*, 'Median net worth to labor income ratio? Want 1.2',mean(NW_to_inc_50)
		
        print*, 'Check, Homeownership rate, want 0.66', real(sum(rent_buy(:,burn:sim_end-1)-1))/real(I*(sim_end-burn+1-1)),size(rent_buy(1,burn:sim_end-1)),sim_end-burn+1-1
		print*, 'Check, Foreclosure rate, want 0.001', real(sum(countD(burn:sim_end-1,5)))/real(I*(sim_end-burn+1-1)),real(sum(countD(burn:sim_end-1,5)))/real(sum(rent_buy(:,burn:sim_end-1)-1))
        print*, 'Employment', mean(employment_b(burn:sim_end-1)),sum(foreclosure(burn:sim_end-burn-1+1))/(sim_end-burn-1+1)

		Laborh=mean(alpha*(pstar(burn:sim_end-1)**(1d0/(1d0-alpha))*LbarH/YT(burn:sim_end-1))/employment_b(burn:sim_end-1)) !/real((T_sim-burn-1+1))
		print*, 'Check, Construction sector employment, Nh. Want 0.05 of total', Laborh !, !Laborh/) !,((HS*LbarH**(alpha-1d0)/YT(YT_loc(init_loc(2,2))))**(1d0/alpha)))
		print*, 'Annual fraction of housing sold, want 0.095', real(sum(countD(burn:sim_end-1,1)))/(real(I*(sim_end-burn-1+1))/2d0)
		print*, 'Check Income Taxes collected, should be 20 percent of income. (paper says output=theta*N_c, but it is really income in their code)', mean(tax_b(burn:sim_end-1)) !/(sum(YT(burn:sim_end)*(employment_b(burn:sim_end)-LaborH))/real((sim_end-burn+1)))*100, size( employment_b(burn:sim_end)),sim_end-burn+1   
		print*, 'Check, Average sized owned house/rented house, want 1.5', mean(avg_buy_rent(2,:)/avg_buy_rent(1,:))
		print*, 'Ratio of the average earnings of owners to renters, want (2.1)', (sum(Earn(2,:)/real(count_buy(1:sim_end))))/sum(Earn(1,:)/real(count_rent(1:sim_end))) !, sum(Earn(2,:)),real(count_buy(1:sim_end)), sum(Earn(1,:)),real(count_rent(1:sim_end))
		print*, 'Homeownership rate of <35 years want 0.37', mean(under_35(1,:))
		print*, 'Median NW at age 75 (28)/Median NW at age 51 (16), want 1.55, otherwise adjust bequest',mean(NW_med1_print/NW_med2_print) 
		print*, mean(rental_sim(burn:sim_end-2))
		print*, mean(Hstar_b(burn+1:sim_end-1))
		print*, mean(Cstar_b(burn:sim_end-2))
		
		
		print*, 'Housing affordability, boom', mean(afford(:,398)), mean(afford(:,399)), mean(afford(:,400)), mean(afford(:,401)), mean(afford(:,402)), mean(afford(:,403)), mean(afford(:,404))
		print*, 'Housing affordability, p10', afford_p10(398), afford_p10(399),afford_p10(400),afford_p10(401),afford_p10(402),afford_p10(403), afford_p10(404)
		print*, 'Housing affordability, p50', afford_p50(398), afford_p50(399),afford_p50(400),afford_p50(401),afford_p50(402),afford_p50(403), afford_p50(404)
		print*, size(rental_sim(burn:sim_end-2))
		print*, size(Hstar_b(burn+1:sim_end-1))
		print*, size(Cstar_b(burn:sim_end-2))
		!print*, Hstar_b(burn+1:sim_end-1)
		print*, 'Housing as a share of total expenditure, want 0.16', mean(rental_sim(burn:sim_end-2)*Hstar_b(burn+1:sim_end-1)) !,rental_sim(burn:sim_end-2)*Hstar_b(burn+1:sim_end)-1!/(Cstar_b(burn:sim_end-2)+rental_sim(burn:sim_end-2)*Hstar_b(burn+1:sim_end-1)))
        print*, 'P10 Housing NW/total NW, want 0.11 (0.12)', mean(H_NW_10) !/(pstar(burn)*Hstar_b(burn)-Lstar_b(burn))!*sum((H_sim_star(burn,:))-L_sim_star(burn,:)))
        print*, 'P50 Housing NW/total NW, want 0.50 (0.38)', mean(H_NW_50) !/(pstar(burn)*Hstar_b(burn)-Lstar_b(burn)) !*sum((H_sim_star(burn,:))-L_sim_star(burn,:))
        print*, 'P90 Housing NW/total NW, want 0.95 (0.80)', mean(H_NW_90)!/(pstar(burn)*Hstar_b(burn)-Lstar_b(burn)) !*sum((H_sim_star(burn,:))-L_sim_star(burn,:))
        print*, 'Fraction of bequests in bottom half of wealth dist, want 0, otherwise adjust bequest',&
        sum((real(count_bequest)/real(count_bequest_tot))/real(sim_end-burn+1-1)) !,real(sum(count_bequest))/real(sim_end-burn+1-1), real(sum(count_bequest_tot))/real(sim_end-burn+1-1)
		
 av_scale=real(I*(sim_end-burn-1+1))
		Calibration=(/mean((sum(Agg_nw,dim=1)/real(I))/(Agg_inc(:)/2d0)),mean(NW_to_inc_50),mean(NW_med1_print/NW_med2_print), 0d0,mean(rental_sim(burn:sim_end-2)*Hstar_b(burn+1:sim_end-1)/(Cstar_b(burn:sim_end-2)+rental_sim(burn:sim_end-2)*Hstar_b(burn+1:sim_end-1))),&
		real(sum(rent_buy(:,burn:sim_end-1)-1))/av_scale,sum(foreclosure(burn:sim_end-burn-1+1))/(sim_end-burn-1+1),&
		mean(H_NW_10),mean(H_NW_50),mean(H_NW_90),mean(avg_buy_rent(2,:)/avg_buy_rent(1,:)),&
		(sum(Earn(2,:)/real(count_buy(1:sim_end))))/sum(Earn(1,:)/real(count_rent(1:sim_end)))  , real(sum(countD(burn:sim_end-1,1)))/(real(I*(sim_end-burn-1+1))/2d0),&
		mean(under_35(1,:)),Laborh/)  !(sum(Earn(2,:,:))/sum(real(rent_buy(:,burn+1:sim_end-1)-1)))/(sum(Earn(1,:,:))/sum((-1)*real(rent_buy(:,burn+1:sim_end-1)-2)))

		print*, '______________________________________________'
		print*, ' Table 4'
		print*, '______________________________________________'
		print*, 'Distribution of housing sizes'
		print*, '... less than 35', sum(counts_h35,dim=2)/real(sim_end-burn+1-1)
		print*, '... less than 49', sum(counts_h49,dim=2)/real(sim_end-burn+1-1)
		print*, '... less than 64', sum(counts_h64,dim=2)/real(sim_end-burn+1-1)
		print*, '... less than 80', sum(counts_h80,dim=2)/real(sim_end-burn+1-1)
		
		
		print*, '______________________________________________'
		print*, ' Table 5'
		print*, '______________________________________________'
		print*, 'Fraction of homeowners with mortgage, want 0.66 (0.56)', sum(real(With_mortgage)/real(count_buy))/(sim_end-burn+1-1),shape(With_mortgage),sim_end-burn+1-1
		print*, 'Fraction of homeowners with heloc, want 0.06 (0.03)', sum(real(With_heloc)/real(count_buy))/(sim_end-burn+1-1)
		!print*, 'rent_buy', shape(With_heloc/count_buy),real(real(With_heloc)/real(count_buy))
		print*, 'Aggregate mortgage debt to housing value, want 0.42 (0.34)', sum((Lstar_b(burn+1:sim_end-1)*real(count_buy)/real(count_LTV))/pstar(burn:sim_end-2))/(sim_end-burn+1-1),shape(sum((Lstar_b(burn+1:sim_end-1)*real(count_buy)/real(count_LTV))/pstar(burn:sim_end-2)))
		print*, 'P10 LTV ratio for mortgages, want 0.15 (0.14)', mean(LTV_10)
		print*, 'P50 LTV ratio for mortgages, want 0.57 (0.58)', mean(LTV_50)
		print*, 'P90 LTV ratio for mortgages, want 0.92 (0.89)', mean(LTV_90)
		print*, 'Share of NW held by...'
		print*, '... bottom quintile, want 0 (0.69)', mean(bottom_NW)
		print*, '... mid quintile, want 0.05 (0.08)', mean(mid_NW)
		print*, '... top quintile, want 0.81 (0.69)', mean(top_NW)
		print*, '... top 10 percent, want 0.7 (0.35)', mean(top10_NW)
		print*, '... top 1 percent, want 0.46 (0.07)', mean(top1_NW)
		print*, 'P10 housing value to earnings, want 0.9 (1.0)', mean(house_value_earn_10)
		print*, 'P50 housing value to earnings, want 2.1 (2.0)', mean(house_value_earn_50)
		print*, 'P90 housing value to earnings, want 5.5 (4.3)', mean(house_value_earn_90)
		
		Calibration_untargeted(1:2)=(/sum(real(With_mortgage)/real(count_buy))/(sim_end-burn+1-1),sum(real(With_heloc)/real(count_buy))/(sim_end-burn+1-1)/)
		Calibration_untargeted(3:8)=(/sum((Lstar_b(burn+1:sim_end-1)*real(count_buy)/real(count_LTV))/pstar(burn:sim_end-2))/(sim_end-burn+1-1), mean(LTV_10),mean(LTV_50), mean(LTV_90), mean(bottom_NW), mean(mid_nW)/)
		Calibration_untargeted(9:14)=(/mean(top_NW), mean(top10_NW), mean(top1_NW), mean(house_value_earn_10), mean(house_value_earn_50), mean(house_value_earn_90)/)
		!(/sum(real(With_mortgage)/real(count_buy))/(sim_end-burn+1), sum(real(With_heloc)/real(count_buy))/(sim_end-burn+1),&

		!mean(Lstar_b(burn+1:sim_end-1)/pstar(burn:sim_end-2)), mean(LTV_10), mean(LTV_50), mean(LTV_90), mean(bottom_NW), mean(mid_nW),&
		!mean(top_NW), mean(top10_NW), mean(top1_NW), mean(house_value_earn_10), mean(house_value_earn_50), mean(house_value_earn_90) /)
		
		print*, '______________________________________________'
		print*, ' Table 6'
		print*, '______________________________________________'	
		print*, 'Distribution of homeowners', sum(counts_h,dim=2)/real(sim_end-burn+1-1)
		print*, 'Distribution of renters', sum(counts_r,dim=2)/real(sim_end-burn+1-1)
		
		
		Print*, 'End Calibration'
		
		Print*,'___________________________________________________________________________________________'
		
        !deallocate(NW_med1)
        !deallocate(NW_med2)

    
    call CPU_time (t2)
    tp2=omp_get_wtime()
	print*, a_raw(:,:,1)
	print*, a_raw(:,:,2)
			OPEN(3, FILE = './results/computational-appendix.txt',status='old',action='write',position="append")
		write(3,"(a,f6.2,a)") '\newcommand{\vfiserial}{$',(t22-t11)/60d0,'$}'
		write(3,"(a,f6.2,a)") '\newcommand{\vfiparallel}{$',(tp22-tp11)/60d0,'$}'
		write(3,"(a,I6,a)") '\newcommand{\iternum}{$',iter,'$}$'

		close(3)
		
    !*************************************************
    !       End KS iteration
    !*************************************************.
	strname=[character(20) :: 'zero','one','two','three','four']
	!print*, strname
OPEN(3, FILE = 'tables/coeffs-table-low.txt',status='old',action='write')
		write(3,"(a,f8.6,a)") '\newcommand{\R2}{$',R2p(1,1),'$}'
		
		do iY=1,nY*nlamb
			do iYY=1,nY*nlamb
				do iM=1,mom+1
				!print*, strname(iM),strname(iY+1),strname(iYY+1)
				!	write(3,"(a,5a,a,f5.3,a)") '\newcommand{\alpha', strname(iM),strname(iY+1),strname(iYY+1),'}{$',a_raw(iY,iYY,iM,Aloc),'$}'
					 
				end do
				!write(3,"(a,5a,a,f5.3,a)") '\newcommand{\samplemean',strname(iY+1),strname(iYY+1),'}{$',exp(a_new(iY,iYY,1)/(1-a_new(iY,iYY,2))),'$}'
			
			end do
		end do 
	!	a_new(iY,iYY,:), '    |',exp(a_new(iY,iYY,1)/(1-a_new(iY,iYY,2))), '    | ', 
		close (3)

   if (solve_for_coeffs==1 .and. type_sim==1 .and. nLamb==2 .and. nphi==1) then
		open (221, file ='coefficients2.csv',status ='replace',action='write')       
        do iM=1,mom+1
			do iY=1,nY*nlamb
				do iYY=1,nY*nlamb
					print*, a_raw(iY,iYY,iM)
					write (221,*) a_raw(iY,iYY,iM)
					!write (2,*) a_raw(1,2,iM,Aloc)
					!write (2,*) a_raw(2,1,iM,Aloc)
					!write (2,*) a_raw(2,2,iM,Aloc)
				end do
			end do 
        end do
        close (221)
		print*, "printing coefficients nlamb==2,nphi=1"
 else  if (solve_for_coeffs==1 .and. type_sim==1 .and. nLamb==2 .and. nphi==3) then
		open (223, file ='coefficients3.csv',status ='replace',action='write')       
        do iM=1,mom+1
			do iY=1,nY*nlamb*nphi
				do iYY=1,nY*nlamb*nphi
					print*, a_raw(iY,iYY,iM)
					write (223,*) a_raw(iY,iYY,iM)
					!write (2,*) a_raw(1,2,iM,Aloc)
					!write (2,*) a_raw(2,1,iM,Aloc)
					!write (2,*) a_raw(2,2,iM,Aloc)
				end do
			end do 
        end do
        close (223)
		print*, "printing coefficients nlamb==2,nphi=3"		

	else if (solve_for_coeffs==1 .and. type_sim==1  .and. nLamb==1) then
		open (222, file ='coefficients1.csv',status ='replace',action='write')       
        do iM=1,mom+1
			do iY=1,nY*nlamb
				do iYY=1,nY*nlamb
					write (222,*) a_raw(iY,iYY,iM)
					!write (2,*) a_raw(1,2,iM,Aloc)
					!write (2,*) a_raw(2,1,iM,Aloc)
					!write (2,*) a_raw(2,2,iM,Aloc)
				end do
			end do 
        end do
		print*, "printing coefficients nlamb==1"
        close (222)
		

	else 
	
		continue
    end if

		
    if (CKW .eqv. .False.) then
        diffp=0d0
    else
        open(21,file='pstar.csv',status='old',action='read')
        do it=1,learn
            read(21,*) pstarold(it)
        end do
        close(21)
        open(23,file='pstar.csv',status='replace')
        do it=1,T_sim-1
            write(23,*) pstar(it)
        end do
        close(23)

      

   
    
    end if
	deallocate(pprime)  
	deallocate(rental)
   
 

	
	deallocate(AggH_b) 
	deallocate(AggS)
	
		if(type_sim==0) then
	deallocate(rent_buy_steady)
	deallocate(L_sim_steady)
	deallocate(H_sim_steady)
	deallocate(B_sim_steady)
	deallocate(C_sim_steady)
	deallocate(Tax_sim_steady)
	deallocate(max_loc_b_steady)
	!deallocate(P)
	deallocate(AggH_s)
	deallocate(Trans)
	
    end if 

    
    
    
    print*, 'time_elapsed, minutes:', (t2-t1)/60d0, (tp2-tp1)/60d0
    print*, 'VFI time, minutes:', (t22-t11)/60d0, (tp22-tp11)/60d0
    print*, 'Simulation time, minutes:', (tsim2-tsim1)/60d0, (tpsim2-tpsim1)/60d0
    print*, 'minutes per iteration', (t2-t1)/60d0/iter,(tp2-tp1)/60d0/iter
    print*, 'diff_total', diff_total
    print*, 'pgrid', p

    dimensions=(/nL,nH,nB,nY*nlamb*nphi,nE,countP,T,nA,burn,I,T_sim/)
    call csvwrite('dimensions',dimensions) 

    call csvwrite('Grid_L',L)
    call csvwrite('Grid_Y',Y)
    call csvwrite('Grid_H',H)
    call csvwrite('Grid_p',p)
    call csvwrite('Grid_T_sim',linspace(1d0,real(T_sim-1,8),T_sim-1))   
    call csvwrite('qstar', qstar(1:T_sim-1))
    call csvwrite('pstar', pstar(1:T_sim-1))
	
	if (type_sim==2) then 
		call csvwrite('a_out0',a_out(1:T_sim-1,1))
		call csvwrite('a_out1',a_out(1:T_sim-1,2))
		call csvwrite('a_out0_all',a_out_all(1:T_sim-1,:,1))
		call csvwrite('a_out1_all',a_out_all(1:T_sim-1,:,2))
		
		!call csvwrite('e_out', e_out(1:T_sim-1))
		call csvwrite('growtherror', (log(pstar(2:T_sim-1))-(a_out(1:T_sim-2,1)+a_out(1:T_sim-2,2)*log(pstar(1:T_sim-2))))/log(pstar(1:T_sim-2)))
		!call csvwrite('growtherror2', (log(pstar(2:T_sim-1))-dot_product(a_out_all(1:T_sim-2,:,1)+a_out_all(1:T_sim-2,:,2)*log(pstar(1:T_sim-2)),FYY(AggT_loc(1:T_sim,:),:)))/log(pstar(1:T_sim-2)))
		
		call csvwrite('SSE', mean((100*log(pstar(399:410))-100*log(exp(a_out(398:409,1)+a_out(398:409,2)*log(pstar(398:409)))))**2))
		print*, 'growth error: boom only, SSE', mean((100*(log(pstar(401:404))-log(exp(a_out(400:403,1)+a_out(400:403,2)*log(pstar(400:403))))))**2)
		print*, 'growth error: full, SSE', mean((100*(log(pstar(400:406))-log(exp(a_out(399:405,1)+a_out(399:405,2)*log(pstar(399:405))))))**2)
	end if 
	call csvwrite('rental',rental_sim(1:T_sim-1))
	call csvwrite('LTV2',LTV_mean(1:T_sim-1))
	!call csvwrite('LTV2',LTV_mean2(1:T_sim-1))
    call csvwrite('Y_sim', YT(1:T_sim-1))
    call csvwrite('YT_loc',YT_loc(1:T_sim))
	call csvwrite('AggT_lgoc',AggT_loc(1:T_sim))
    call csvwrite('tax_b', tax_b(1:T_sim-1))
    call csvwrite('Lstar_b', Lstar_b)
    call csvwrite('Cstar_b', Cstar_b)
    call csvwrite('Hstar_b', Hstar_b)
    call csvwrite('Hstar_s', Hstar_s)
    call csvwrite('Bstar_b', Bstar_b)
	call csvwrite('Homeownership',homeownership)
	call csvwrite('Foreclosure',foreclosure)
    call csvwrite('income_b',income_b)

    
    allocate(L_sim_dist(I,nY*nlamb,nY*nlamb), B_sim_dist(I,nY*nlamb,nY*nlamb), H_sim_dist(I,nY*nlamb,nY*nlamb),C_sim_dist(I,nY*nlamb,nY*nlamb), E_sim_dist(I,nY*nlamb,nY*nlamb), Tax_sim_dist(I,nY*nlamb,nY*nlamb),wealth_sim_dist(I,nY*nlamb,nY*nlamb), max_loc_sim_dist(I,nY*nlamb,nY*nlamb))
    allocate(age_sim_dist(I,nY*nlamb,nY*nlamb)) 
	allocate(rent_buy_sim_dist(I,nY*nlamb,nY*nlamb))


    do iY=1,nY
        do iYY=1,nY
            do iI=1,I
                L_sim_dist(iI,iY,iYY)=L_sim_star(iI,init_loc(iY,iYY)+1)
                H_sim_dist(iI,iY,iYY)=H_sim_star(iI,init_loc(iY,iYY)+1)
                B_sim_dist(iI,iY,iYY)=B_sim_star(iI,init_loc(iY,iYY)+1)
                age_sim_dist(iI,iY,iYY)=age_sim_bor(iI,init_loc(iY,iYY))
                C_sim_dist(iI,iY,iYY)=C_sim_star(iI,init_loc(iY,iYY))
                E_sim_dist(iI,iY,iYY)=income_star(iI,init_loc(iY,iYY))
                Tax_sim_dist(iI,iY,iYY)=tax_star(iI,init_loc(iY,iYY))
                wealth_sim_dist(iI,iY,iYY)=pstar(init_loc(iY,IYY))*H_sim_star(iI,init_loc(iY,iYY))*(1d0-delta-tauh)-L_sim_star(iI,init_loc(iY,iYY))+B_sim_star(iI,init_loc(iY,iYY))
                rent_buy_sim_dist(iI,iY,iYY)=rent_buy(iI,init_loc(iY,iYY))
				max_loc_sim_dist(iI,iY,iYY)=max_loc_sim_star(iI,init_loc(iY,iYY))
                
            end do
        end do
    end do
	call csvwrite('type_sim',type_sim)
	if (type_sim==1) then 
		call csvwrite('mu_age_dist',Pack(age_sim_dist, .True.))
		call csvwrite('mu_Ci_dist', Pack(C_sim_dist, .True.))            
		call csvwrite('mu_max_loc', Pack(max_loc_sim_dist, .True.))
		call csvwrite('mu_Ei_dist',Pack(ET_loc(init_loc(1,1),:), .True.))
		call csvwrite('mu_logei_dist',Pack(E_sim_dist, .True.))
		call csvwrite('mu_Li_dist',Pack(L_sim_dist, .True.))
		call csvwrite('mu_Hi_dist',Pack(H_sim_dist, .True.))  
		call csvwrite('mu_Bi_dist',Pack(B_sim_dist, .True.))  
		call csvwrite('mu_Tax_dist',Pack(Tax_sim_dist, .True.))  
		call csvwrite('mu_wealth_dist',Pack(wealth_sim_dist, .True.))
		call csvwrite('mu_rent_buy_dist',Pack(rent_buy_sim_dist, .True.))
		
	
		call csvwrite('Calibration',Calibration)
		call csvwrite('Calibration_untargeted',Calibration_untargeted)
		call csvwrite('counts_h35', sum(counts_h35,dim=2)/real(sim_end-burn+1-1))
		call csvwrite('counts_h49', sum(counts_h49,dim=2)/real(sim_end-burn+1-1))
		call csvwrite('counts_h64', sum(counts_h64,dim=2)/real(sim_end-burn+1-1))
		call csvwrite('counts_h80', sum(counts_h80,dim=2)/real(sim_end-burn+1-1))
		call csvwrite('counts_h', sum(counts_h,dim=2)/real(sim_end-burn+1-1))
		call csvwrite('counts_r', sum(counts_r,dim=2)/real(sim_end-burn+1-1))
	end if 
    print*, 'End program'

end program main

