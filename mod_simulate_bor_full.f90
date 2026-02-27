!       This file contains the simulation code for Krussel-Smith Simulation

module mod_simulate_bor_full
    use omp_lib
    use parameters
    use mod_functions
    use mod_matlab
    use mod_interp
    use mod_root1dim
    implicit none
contains
    
    subroutine sub_simulate_bor(age_sim,rent_buy,L,B,Bs,H,R,Y,chi,V_b,LL_b,BB_b,HH_b,c_b,AggT_loc,YT_loc,LambT_loc,loge,ET_loc,&
    AggH,AggS,AggH_s,p,countp,countD,sim_end,pstar,Lstar_b,Bstar_b,Hstar_b,Hstar_s,Cstar_b,qstar_b,income_b,C_sim_star,&
    income_star,L_sim_star,B_sim_star,H_sim_star,tax_b,tax_star,housing_tax_b,employment_b,employment_star,ret_b,max_loc_sim_star,max_loc_b,pprime,homeownership,q)

      

    implicit none

        integer :: it,iI,iP,count_default,Aloc,Aloc2,iE,iJ,psub,psub2,psub3,psub4,iD,Ijprint,iEprint,Eagent,buy,rent,iAgg
		real :: inc
        integer, dimension(:,:), intent(in) :: age_sim
		integer, dimension(:,:), intent(out) :: rent_buy
        integer, intent(in) :: countP,sim_end
        real(8), allocatable :: H_sim(:,:) !,B_sim,L_sim
        real(8), allocatable :: ret(:,:),q_sim_star(:,:)
        real(8), allocatable :: V_sim(:,:,:)
        real(8), allocatable :: V_sim_star(:,:,:)
        integer, allocatable :: max_loc_sim(:,:)
		integer, allocatable :: ExoT_loc(:)
		integer, allocatable :: ilo_H(:,:)
		integer, allocatable :: Hloc(:)
		real(8), allocatable :: Wealth(:,:), Wealth_star(:)
		real(8), allocatable :: theta_H(:,:)
        integer, dimension(:,:), intent(out) :: max_loc_sim_star
        integer, dimension(:,:), intent(out) :: max_loc_b
        integer, dimension(nP,sim_end,nD+2) :: max_loc_p
        integer, dimension(sim_end-1) :: count_work
        real(8), dimension(:), intent(in) :: L,H,Y,chi,B,R,p,Bs
        real(8), dimension(:,:), intent(in) :: loge
        integer, dimension(:), intent(in) :: YT_loc, AggT_loc,LambT_loc
        integer, dimension(:,:), intent(in) :: ET_loc
		real(8), dimension(:,:,:,:,:,:,:), intent(in) :: V_b,HH_b,c_b,LL_b,BB_b
		real(8), dimension(:,:,:,:,:,:), intent(in) :: q
		real(8), dimension(:,:),intent(out) :: L_sim_star,B_sim_star,H_sim_star
        real(8), dimension(:,:), intent(out) :: C_sim_star, income_star,tax_star,employment_star !,constraint_star
        real(8), dimension(:), intent(out) :: pstar,Lstar_b,Bstar_b,Hstar_b,Hstar_s,Cstar_b,qstar_b,income_b,tax_b,employment_b,housing_tax_b
        real(8), dimension(:), intent(out) :: ret_b, homeownership
        real(8), dimension(:,:), intent(inout) :: AggH,AggS
        real(8), dimension(:), intent(in) :: AggH_s
        real(8), dimension(nP,sim_end) ::  AggR
        real(8), dimension(nP) :: AggD
        integer, dimension(nB,T) :: E_dist
		 real(8), dimension(:,:,:), intent(in) :: pprime 
	
		integer :: count_buy_p, count_rent_p, count, Eindex
		integer, dimension(nE) :: countE
        integer, dimension(nP,sim_end) :: count_AggR_own, count_AggR_rent, count_AggD,count_AggH,count_AggS, count_AggR
        integer, dimension(2) :: ilo_bstar, ilo_hstar, ilo_pstar
		real(8), dimension(2) :: theta_bstar, theta_hstar, theta_pstar
		integer, dimension(:,:),intent(out) :: countD
		real(8), dimension(:,:), allocatable :: housing_tax
		real(8), dimension(nE) :: edist,zerodist,nonzerodist
		real(8), dimension(nE,2) :: initliqdist, initilldist, initnwdist, initadist
		real(8), dimension(nE) :: initliq, initnw, initill, initashare
		integer, allocatable :: max_count_age(:,:,:)
		
		if (buy_grid==0) then
			rent=4
			buy=1
		else
			rent=7
			buy=6
		end if 


		allocate(max_loc_sim(I,nP))
		allocate(V_sim_star(I,sim_end-1,nD+2))
		allocate(q_sim_star(I,sim_end-1))
		allocate(ret(I,sim_end-1))
		allocate(housing_tax(I,sim_end-1))
		allocate(V_sim(I,nP,nD+2))
		allocate(H_sim(I,nP))
		allocate(ExoT_loc(I))
		allocate(max_count_age(sim_end,T,nD+2))
		allocate(Wealth(I,nP))
		allocate(Wealth_star(I))
	
		
       
		Lstar_b(1)=0d0 !AggL(1,1)
        Hstar_b(1)=0d0 !AggH(1,1)
        Bstar_b(1)=0d0 !AggB(1,1)
		
	
	
		
		
       
		q_sim_star=0d0
		housing_tax=0d0
        tax_star=0d0
        employment_star=0d0
        count_work=0
        ret=0d0
        max_loc_p=0
		AggR=0d0
        AggH=0d0
        AggS=0d0
        AggS(:,1)=2d0*(1d0-delta)+AggH_s
        AggS(:,2)=2d0*(1d0-delta)+AggH_s
		count_AggR=0
        count_AggH=0
        count_AggS=0
		V_sim_star=death
		V_sim=death
		homeownership=0d0
		countD=0
		max_loc_sim_star=0
		max_loc_b=0
		pstar=0d0
		max_count_age=0
		income_star=0d0
		H_sim=0d0
        print*, AggS(:,2)
        print*, AggH_s
		
		H_sim_star(:,1)=0d0
		L_sim_star(:,1)=0d0
		B_sim_star(:,1)=0d0
		rent_buy(:,1)=1

		
		OPEN(3, FILE = './KMV_Input/zydist.txt',status='old',action='read')
				read(3,*) edist(:)
		close(3)
		OPEN(1, FILE = './KMV_Input/initw.txt')
			DO iE=1,nE
				READ(1,*) initadist(iE,1),initadist(iE,2),initashare(iE)
				initadist(iE,2)=1.0d0-initadist(iE,1)
			END DO
		CLOSE(1)
		
	do iE=1,nE
		zerodist(iE)=int(edist(iE)*initadist(iE,1)*I/T)
		nonzerodist(iE)=int(edist(iE)*initadist(iE,2)*I/T)
	end do
		count=0
		do while (sum(zerodist)<real(I)/real(T)*sum(edist*initadist(:,1))) 
			count=count+1
			zerodist(count)=zerodist(count)+1
		end do		
		zerodist(1:4)=zerodist(1:4)+1
		count=0
		do while (sum(nonzerodist)<real(I)/real(T)*sum(edist*initadist(:,2))) 
			count=count+1
			nonzerodist(count)=nonzerodist(count)+1
		end do	
		

	initliqdist=0.0d0
	initilldist=0d0
	initnwdist=0.0d0
	initnw=0.0d0
	OPEN(1, FILE = './KMV_Input/initwbroad.txt',status='old',action='read')
		DO iE=1,nE
			READ(1,*) initliqdist(iE,1),initliqdist(iE,2),initliq(iE),initilldist(iE,1),initilldist(iE,2),initill(iE),initnwdist(iE,1),initnwdist(iE,2),initnw(iE)
		END DO
	CLOSE(1)
	

	initliq=initliq/(DataAvAnnualEarns/(Numeraire))
	initill=initill/(DataAvAnnualEarns/(Numeraire))
	initnw=initnw/(DataAvAnnualEarns/(Numeraire))
	

	

	
        do it=1,sim_end-1
				ExoT_loc=0
				ExoT_loc=ET_loc(:,it)+nE*(AggT_loc(it)-1)
				!print*, ExoT_loc(1:100)
				V_sim=death
				
				Wealth=0d0
				Wealth_star=0d0
				
				
			
			!$OMP Parallel Default(Shared) private(iP,iI,iD,inc)  
			!$OMP DO collapse(2) reduction(+:AggR) reduction(+:AggH) reduction(+:max_loc_p)
			do iP=1,countP   
				do iI=1,I 
					if (age_sim(iI,it)>1) then 
						Wealth(iI,iP)=B_sim_star(iI,it)+(1d0-delta-tauh-cost)*p(iP)*H_sim_star(iI,it)-(1d0+r_borrow_grid(LambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it)
					else 
						Wealth(iI,iP)=B_sim_star(iI,it)
					end if 
					
					
				   	if (age_sim(iI,it)<Jret) then ! This is the next period value we are solving for
						inc=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))))-FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))
					else
						inc=FnTax(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)))-FnTaxm(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)),L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))
					end if 

                    do iD=1,nD+2

						if (rent_buy(iI,it)==2 .and. (iD==6 .or. iD==7)) then

							if(Wealth(iI,iP)+inc>=0d0) then 
									V_sim(iI,iP,iD)=interp3qnoextrap(L,H,B,V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),L_sim_star(iI,it),H_sim_star(iI,it),Wealth(iI,iP)+inc)
							else
									V_sim(iI,iP,iD)=interp3qnoextrap(L,H,Bs,V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),L_sim_star(iI,it),H_sim_star(iI,it),Wealth(iI,iP)+inc)
							end if
						else if ((rent_buy(iI,it)==2 .and. (iD==2 .or. iD==3)) .or. (rent_buy(iI,it)==2 .and. iD==5)) then
								V_sim(iI,iP,iD)=interp3qnoextrap(L,H,B,V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))
							

						else if (rent_buy(iI,it)==1 .and. (iD==4 .or. iD==1))  then
						
								V_sim(iI,iP,iD)=interp1qnoextrap(B,V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),B_sim_star(iI,it))
								
								
						else if ((rent_buy(iI,it)==2 .and. (iD==4 .or. iD==1))) then		
						if(Wealth(iI,iP)+inc>=0d0) then 

									V_sim(iI,iP,iD)=interp1qnoextrap(B,V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),Wealth(iI,iP)+inc)
						
							else
									V_sim(iI,iP,iD)=interp1qnoextrap(Bs,V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),Wealth(iI,iP)+inc)
					
							end if
						end if 
						
						
                    end do
	

					  
         
                   if (rent_buy(iI,it)==1 .and. age_sim(iI,it)<T) then
                      max_loc_sim(iI,iP)=maxloc(V_sim(iI,iP,(/1,4/)),dim=1)*3-2
					else if (rent_buy(iI,it)==1 .and. age_sim(iI,it)==T) then
					! Renter's can't buy a house in the final period of life.
						max_loc_sim(iI,iP)=4
					else if (rent_buy(iI,it)==2 .and. age_sim(iI,it)<T) then
                      max_loc_sim(iI,iP)=maxloc((/V_sim(iI,iP,buy),V_sim(iI,iP,2:3),V_sim(iI,iP,rent),V_sim(iI,iP,5)/),dim=1)   !maxloc(V_sim(iI,iP,(/1,2,3,4,5/)),dim=1)  

					  if (max_loc_sim(iI,iP)==1) then
							max_loc_sim(iI,iP)=buy
					  else if (max_loc_sim(iI,iP)==4) then 
							max_loc_sim(iI,iP)=rent
						end if 
					else if (rent_buy(iI,it)==2 .and. age_sim(iI,it)==T) then
						max_loc_sim(iI,iP)=maxloc((/V_sim(iI,iP,2),V_sim(iI,iP,4),V_sim(iI,iP,5)/),dim=1) !maxloc(V_sim(iI,iP,(/2,4,5/)),dim=1)  
						if (max_loc_sim(iI,iP)==1) then
							max_loc_sim(iI,iP)=2
						else if (max_loc_sim(iI,iP)==2) then
							max_loc_sim(iI,iP)=4
						else if (max_loc_sim(iI,iP)==3) then
							max_loc_sim(iI,iP)=5
						end if 
                   end if 
                    
                    !Own to own: pay or refi
                    if (rent_buy(iI,it)==2 .and. (max_loc_sim(iI,iP)==2 .or. max_loc_sim(iI,iP)==3)) then  
							H_sim(iI,iP)=H_sim_star(iI,it)

							AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP)
 
						!Own to rent or buy
                     else if ((max_loc_sim(iI,iP)==buy .or. max_loc_sim(iI,iP)==rent) .and. rent_buy(iI,it)==2) then          
							if (Wealth(iI,iP)+inc<0d0) then !B(1)) then
								H_sim(iI,iP)=interp1qnoextrap(Bs,HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP)),&
								Wealth(iI,iP)+inc)
							else 
								H_sim(iI,iP)=interp1qnoextrap(B,HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP)),&
								Wealth(iI,iP)+inc)
							end if 			


                        AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP)


					! Default or rent-rent and rent-buy
                    else if (rent_buy(iI,it)==2 .and. max_loc_sim(iI,iP)==5) then
						H_sim(iI,iP)=R(1)
						AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP) !*mu_own(it,age_sim(iI,it
						

					else if (rent_buy(iI,it)==1 .and. (max_loc_sim(iI,iP)==1 .or. max_loc_sim(iI,iP)==4)) then     !Own to rent or default
                        H_sim(iI,iP)=interp1qnoextrap(B,HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP)),B_sim_star(iI,it))

            

                        AggR(iP,it+1)=AggR(iP,it+1)+H_sim(iI,iP) 
       


                    end if


               	max_loc_p(iP,it,max_loc_sim(iI,iP))=max_loc_p(iP,it,max_loc_sim(iI,iP))+1
				
			
                end do
			end do
			!$OMP END DO
			!$OMP END PARALLEL

			Eindex=1
			countE=0
			count_rent_p=sum((-1)*(rent_buy(:,it)-2))
			count_buy_p=sum(rent_buy(:,it)-1)


      
			do iP=1,nP
            
				if (count_buy_p>0) then					
					AggH(iP,it+1)=AggH(iP,it+1)/count_buy_p!/count_AggH(iP,it+1) !/I
				else
					AggH(iP,it+1)=0
				end if

			   ! Sellers/Defaulter

				if (count_rent_p>0) then
					AggR(iP,it+1)=AggR(iP,it+1)/count_rent_p
				else
					AggR(iP,it+1)=0
				end if 
                
				

				 
				

				
	
                if (it==1) then
					AggH(iP,it)=0d0
					AggR(iP,it)=0d0
					count_AggH(iP,it)=0
					count_AggR(iP,it)=0
                    do iI=1,I 
                        if (rent_buy(iI,it)==2) then
                            AggH(iP,it)=AggH(iP,it)+H_sim_star(iI,it) 
                        else if (rent_buy(iI,it)==1) then
                            AggR(iP,it)=AggR(iP,it)+H_sim_star(iI,it)
                        end if
                    end do
					if (sum(rent_buy(:,it)-1)>0) then
						AggH(iP,it)=AggH(iP,it)/sum(rent_buy(:,it)-1)
					else
						AggH(iP,it)=0
					end if
					if ( sum((-1)*(rent_buy(:,it)-2))>0) then
					
						AggR(iP,it)=AggR(iP,it)/sum((-1)*(rent_buy(:,it)-2))
					else
						AggR(iP,it)=0
					end if 
                end if
                 
			end do

            
            Aggd=real(count_buy_p)/real(I)*AggH(:,it+1)+real(count_rent_p)/real(I)*AggR(:,it+1)
            if (it>1) then
              AggS(:,it+1)=interp1q(p,AggS(:,it),pstar(it-1))*(1d0-delta)+AggH_s
            end if 
            

			
			psub=1
			psub3=countP
			if (it<600) then
					print*,'---------------------------------------------'
					print*, it, YT_loc(it),lambT_loc(it),AggT_loc(it)
					print*, 'AggD', Aggd
					print*, 'AggS', AggS(:,it+1)
					print*, 'AggH', AggH(:,it+1)
					print*, 'AggR_t', AggR(:,it+1)
			end if 
		 
		    if(AggD(psub3) .gt. AggS(psub3,it+1)) then
				pstar(it)=p(countp)
			elseif(AggD(psub) .lt. AggS(psub,it+1)) then
				pstar(it)=p(1)
			else
				call sub_bisect_excess_demand1(pstar(it),excess_demand1,p,Aggd,AggS(:,it+1),p(psub),p(psub3))
			endif

						Select Case (display)
				Case(0)
          

		 
			case(1)
				continue
			end select
			
			call base_fun(p,pstar(it),ilo_pstar(:),theta_pstar(:))
					!if (it<5 .or. it>=400) then
					if (it<600) then
					print*,  pstar(it),interp1q(p,Aggd,pstar(it))
	
					

					print*, 'AggD*', interp1q(p,Aggd,pstar(it))
					print*, 'AggH*', interp1q(p, AggH(:,it+1),pstar(it))
					print*, 'AggR*', interp1q(p, AggR(:,it+1),pstar(it))
					print*, 'Fraction Homeowner', real(count_buy_p)/real(I)
					print*, 'Fraction renter', real(count_rent_p)/real(I)
					print*, 'Percent buy/sell', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,1))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,1)))*(100d0/real(I)),(theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,6))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,6)))*(100d0/real(I))
					print*, 'Percent stay', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,2))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,2)))*(100d0/real(I))
					print*, 'Percent refi', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,3))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,3)))*(100d0/real(I))
					print*, 'Percent rent', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,4))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,4)))*(100d0/real(I)),(theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,7))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,7)))*(100d0/real(I))
					print*, 'Percent default', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,5))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,5)))*(100d0/real(I))

				
					end if 
      
            count_AggR(:,it+1)=0
            count_AggH(:,it+1)=0
            count_AggS(:,it+1)=0
			
		!reduction(+:E_print) reduction(+:B_print) 
     	!$OMP Parallel Default(Shared) private(iP,iI,iD,inc)
		!$OMP DO reduction(+:countD) reduction(+:count_work) reduction(+:max_loc_b)reduction(+:max_count_age) reduction(+:homeownership)
            do iI=1,I
                do iD=1,nD+2
                    V_sim_star(iI,it,iD)=theta_pstar(1)*V_sim(iI,ilo_pstar(1),iD)+theta_pstar(2)*V_sim(iI,ilo_pstar(2),iD)
                end do
					
                !print*, maxloc(V_sim_star(it,iI,:))
                    if (rent_buy(iI,it)==1 .and. age_sim(iI,it)<T) then
                        max_loc_sim_star(iI,it)=maxloc(V_sim_star(iI,it,(/1,4/)),dim=1)
                        if (max_loc_sim_star(iI,it)==2) then
                            max_loc_sim_star(iI,it)=4
                        end if
					else if (rent_buy(iI,it)==1 .and. age_sim(iI,it)==T) then
					! Renter's can't buy a house in the final period of life.
						max_loc_sim_star(iI,it)=4
                     else if (rent_buy(iI,it)==2 .and. age_sim(iI,it)<T) then
                        max_loc_sim_star(iI,it)=maxloc((/V_sim_star(iI,it,buy),V_sim_star(iI,it,2:3),V_sim_star(iI,it,rent),V_sim_star(iI,it,5)/),dim=1)
						if (max_loc_sim_star(iI,it)==1) then
							max_loc_sim_star(iI,it)=buy
						else if (max_loc_sim_star(iI,it)==4) then
							max_loc_sim_star(iI,it)=rent
						end if 
					else if (rent_buy(iI,it)==2 .and. age_sim(iI,it)==T) then
						max_loc_sim_star(iI,it)=maxloc((/V_sim_star(iI,it,2),V_sim_star(iI,it,4),V_sim_star(iI,it,5)/),dim=1)
						if (max_loc_sim_star(iI,it)==1) then
							max_loc_sim_star(iI,it)=2
						else if (max_loc_sim_star(iI,it)==2) then
							max_loc_sim_star(iI,it)=4
						else if (max_loc_sim_star(iI,it)==3) then
							max_loc_sim_star(iI,it)=5
						end if  
                    end if
				
					countD(it,max_loc_sim_star(iI,it))=countD(it,max_loc_sim_star(iI,it))+1
					  
                    if (rent_buy(iI,it)==2 .and. (max_loc_sim_star(iI,it)==2 .or. max_loc_sim_star(iI,it)==3)) then !Own to refi and pay
                        
					
							H_sim_star(iI,it+1)=H_sim_star(iI,it)
							

							L_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(L,H,B,LL_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))+theta_pstar(2)*interp3qnoextrap(L,H,B,LL_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))
							

						
							B_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(L,H,B,BB_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))+theta_pstar(2)*interp3qnoextrap(L,H,B,BB_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))
							
							c_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(L,H,B,c_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))+theta_pstar(2)*interp3qnoextrap(L,H,B,c_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))						
				
						
						
							q_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(L,H,B,q(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))+theta_pstar(2)*interp3qnoextrap(L,H,B,q(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it))
							
						
						
                       
						if (age_sim(iI,it+1)>1) then
							rent_buy(iI,it+1)=2
						else
							rent_buy(iI,it+1)=1
						end if 
						housing_tax(iI,it)=tauh*pstar(it)*H_sim_star(iI,it)

                    else if (rent_buy(iI,it)==2 .and. (max_loc_sim_star(iI,it)==buy .or. max_loc_sim_star(iI,it)==rent)) then          !Own to Own
					
					
					if (age_sim(iI,it)<Jret) then ! This is the next period value we are solving for
						inc=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))))-FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))
					else
						inc=FnTax(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)))-FnTaxm(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)),L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))
					end if 
					if (age_sim(iI,it)>1) then 
						Wealth_star(iI)=B_sim_star(iI,it)+(1d0-delta-tauh-cost)*pstar(it)*H_sim_star(iI,it)-(1d0+r_borrow_grid(LambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it)
					else 
						Wealth_star(iI)=B_sim_star(iI,it)
					end if 
							if (Wealth_star(iI)+inc<0d0) then
								H_sim_star(iI,it+1)=interp2qnoextrap(Bs,p,HH_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))
              
								L_sim_star(iI,it+1)=interp2qnoextrap(Bs,p,LL_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))

							
								B_sim_star(iI,it+1)=interp2qnoextrap(Bs,p,BB_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))
								
								C_sim_star(iI,it)=interp2qnoextrap(Bs,p,C_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))	


								q_sim_star(iI,it)=interp2qnoextrap(Bs,p,q(1,1,:,ExoT_loc(iI),:,age_sim(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))				
								
							else
								H_sim_star(iI,it+1)=interp2qnoextrap(B,p,HH_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))
							
								L_sim_star(iI,it+1)=interp2qnoextrap(B,p,LL_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))

							
								B_sim_star(iI,it+1)=interp2qnoextrap(B,p,BB_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))
								
								C_sim_star(iI,it)=interp2qnoextrap(B,p,C_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))		

								q_sim_star(iI,it)=interp2qnoextrap(B,p,q(1,1,:,ExoT_loc(iI),:,age_sim(iI,it)),&
								Wealth_star(iI)+inc,pstar(it))			
								
							end if 

							housing_tax(iI,it)=tauh*pstar(it)*H_sim_star(iI,it)
							
						! This makes zero sense here?	
                        if (max_loc_sim_star(iI,it)==buy .and. age_sim(iI,it)<T .and. age_sim(iI,it+1)>1) then
                            rent_buy(iI,it+1)=2

						else if (max_loc_sim_star(iI,it)==buy .and. (age_sim(iI,it)==T .or. age_sim(iI,it+1)==1)) then 
						    rent_buy(iI,it+1)=1

                        else if (max_loc_sim_star(iI,it)==rent) then
                            rent_buy(iI,it+1)=1

                        end if


                    else if (rent_buy(iI,it)==2 .and. max_loc_sim_star(iI,it)==5) then
						
					

							q_sim_star(iI,it)=interp4qnoextrap(L,H,B,p,q(:,:,:,ExoT_loc(iI),:,age_sim(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))
						
					
						H_sim_star(iI,it+1)=R(1)
						rent_buy(iI,it+1)=1

						L_sim_star(iI,it+1)=0d0 
						


							B_sim_star(iI,it+1)=interp4qnoextrap(L,H,B,p,BB_b(:,:,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))
							
							C_sim_star(iI,it)=interp4qnoextrap(L,H,B,p,C_b(:,:,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))



					else if (rent_buy(iI,it)==1 .and. (max_loc_sim_star(iI,it)==1 .or. max_loc_sim_star(iI,it)==4)) then     !Own to rent or default
                        
						H_sim_star(iI,it+1)=interp2qnoextrap(B,p,HH_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
                        B_sim_star(iI,it),pstar(it)) !sum(theta_pstar*H_sim(iI,ilo_pstar(:)))

 
						L_sim_star(iI,it+1)=interp2qnoextrap(B,p,LL_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
						B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))
							
	   
						B_sim_star(iI,it+1)=interp2qnoextrap(B,p,BB_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
						B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))
							
						C_sim_star(iI,it)=interp2qnoextrap(B,p,C_b(1,1,:,ExoT_loc(iI),:,age_sim(iI,it),max_loc_sim_star(iI,it)),&
						B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))
							
						q_sim_star(iI,it)=interp2qnoextrap(B,p,q(1,1,:,ExoT_loc(iI),:,age_sim(iI,it)),&
						B_sim_star(iI,it),pstar(it)) !,(/ilo_bstar(iI),ilo_pstar/))


						if (max_loc_sim_star(iI,it)==1 .and. age_sim(iI,it+1)>1) then
							rent_buy(iI,it+1)=2
						else if (max_loc_sim_star(iI,it)==1 .and. age_sim(iI,it+1)==1) then
							rent_buy(iI,it+1)=1
						else if (max_loc_sim_star(iI,it)==4) then
							rent_buy(iI,it+1)=1
						end if 



                    end if

                if (age_sim(iI,it)<Jret) then
                    income_star(iI,it)=income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it)))
                    if (max_loc_sim_star(iI,it)==2 .or. max_loc_sim_star(iI,it)==3) then

						tax_star(iI,it)=FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),(1d0+r_borrow_grid(LambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))/income_star(iI,it)
					else
						tax_star(iI,it)=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))))/income_star(iI,it)
					end if 
                    employment_star(iI,it)=income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it)))/Y(YT_loc(it))
                    count_work(it)=count_work(it)+1     
                else 
                    ret(iI,it)=income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y))
                    income_star(iI,it)=ret(iI,it)
					if (max_loc_sim_star(iI,it)==2 .or. max_loc_sim_star(iI,it)==3) then
						tax_star(iI,it)=FnTaxm(ret(iI,it),(1d0+r_borrow_grid(LambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))/ret(iI,it)
					else 
						tax_star(iI,it)=FnTax(ret(iI,it))/ret(iI,it)
					end if 
                end if     
	
                    max_loc_b(it,max_loc_sim_star(iI,it))=max_loc_b(it,max_loc_sim_star(iI,it))+1

						max_count_age(it,age_sim(iI,it),max_loc_sim_star(iI,it))=max_count_age(it,age_sim(iI,it),max_loc_sim_star(iI,it))+1
						
							homeownership(it)=homeownership(it)+real(rent_buy(iI,it)-1d0)	  


		  end do
			!$OMP END DO
			!$OMP END PARALLEL
					


				! Making all future 
			do iI=1,I	
				if (age_sim(iI,it)==T) then
					h_sim_star(iI,it+1)=0d0
					l_sim_star(iI,it+1)=0d0
					Eagent=ET_loc(iI,it+1)
					countE(Eagent)=countE(Eagent)+1
					
					if (countE(Eagent)<=Zerodist(Eagent) ) then
						!b_sim_star(iI,it+1)=initnw(ET_loc(iI,it+1))
						b_sim_star(iI,it+1)=0d0
					else if (countE(Eagent)>Zerodist(ET_loc(iI,it+1)) .and. countE(Eagent)<=nonZerodist(Eagent)+Zerodist(Eagent) ) then
						b_sim_star(iI,it+1)=initnw(Eagent)
					end if 
					if (b_sim_star(iI,it+1)>maxval(initnw)) then
						b_sim_star(iI,it+1)=0d0
					end if 
					
				end if 


			end do 

            
			select case(display)	
				case(1)
			
				case(0)
				continue
				end select

             employment_b(it)=sum(employment_star(:,it))/real(count_work(it))

             
          
			      income_b(it)=sum(income_star(:,it))/real(I)
                    
            tax_b(it)=sum(tax_star(:,it))/real(I)

            ret_b(it)=sum(ret(:,it))/real(I-count_work(it))
            Cstar_b(it)=sum(C_sim_star(:,it))/real(I)
            
        
          if (it>1) then 
                
                Lstar_b(it)=sum(L_sim_star(:,it))/real(sum(rent_buy(:,it)-1))
                Hstar_b(it)=sum(H_sim_star(:,it))/real(I)
                Bstar_b(it)=sum(B_sim_star(:,it))/real(I)
               qstar_b(it-1)=0d0 !sum(q_sim_star(:,it-1))/real(I)                ! Q doesn't need to be scaled
                Hstar_s(it)=(sum(H_sim_star(:,it))-sum(H_sim_star(:,it-1))*(1-delta))/real(I)
                housing_tax_b(it)=sum(housing_tax(:,it))/real(sum(rent_buy(:,it)-1))
                !employment_b(it-1)=sum(employment_star(it-1,:)) !/real(count_work(it-1))
				homeownership(it)=homeownership(it)/real(I) !sum(real(rent_buy(:,it))-1d0)/real(I)
							select case(display)	
				case(1)
				case(0)
				continue
				end select
				
               ! print*, it !'AggB',Bstar_b(it)
            else
                print*, 'employment', employment_b(it) , count_work(it)         
            end if
            
            if (it==sim_end/4 .or. it==sim_end/2 .or. it==sim_end*3/4) then
                print*, 'it, employment', it, employment_b(it) , count_work(it),pstar(it)
            end if
       end do   
	deallocate(max_loc_sim)
	deallocate(V_sim_star)
	deallocate(ret)
	deallocate(housing_tax)
	deallocate(V_sim)   
	deallocate(H_sim)
	deallocate(ExoT_loc)
	deallocate(max_count_age)
         !print*, maxval(H_sim), maxval(L_sim)
       
    end subroutine sub_simulate_bor
end module mod_simulate_bor_full
