!       This file contains the simulation code for Krussel-Smith Simulation

module mod_simulate_steady
    use omp_lib
    use parameters
    use mod_functions
    use mod_matlab
    use mod_interp
    use mod_root1dim
    implicit none
contains
    
    subroutine sub_simulate_steady(age_sim,L_sim,B_sim,H_sim,rent_buy,C_sim,tax_sim,L,B,H,R,Y,chi,LL_b,BB_b,HH_b,c_b,AggT_loc,YT_loc,LambT_loc,loge,ET_loc,&
    AggH,AggH_s,p,countp,pstar,Lstar_b,Bstar_b,Hstar_b,Cstar_b,tax_b,housing_tax_b,income_star,income_b,ret_b,V_b,employment_b,employment_star,max_loc_sim)
    

      implicit none
        integer :: it,iI,iP,count_default,Aloc,iE,iJ,iD,count
        integer, dimension(:,:), intent(in) :: age_sim
		integer, dimension(:,:,:), intent(out) :: rent_buy
        integer, intent(in) :: countP
		integer, allocatable :: ExoT_loc(:)
        real(8), dimension(:,:,:), intent(out) :: L_sim,H_sim,B_sim,C_sim,tax_sim
        real(8), dimension(:,:,:), allocatable :: q_sim,housing_tax
		real(8), dimension(I,T_sim-1) :: ret
        real(8), dimension(:,:,:), allocatable :: V_sim
        integer, dimension(I) :: hloc
        integer, dimension(:,:,:), intent(out) :: max_loc_sim
        integer, dimension(T_sim-1) :: count_work
        real(8), dimension(:), intent(in) :: L,H,p,Y,chi,B,R
        real(8), dimension(:,:), intent(in) :: loge
        integer, dimension(:), intent(in) :: AggT_loc,YT_loc,LambT_loc
        integer, dimension(:,:), intent(in) :: ET_loc
		real(8), dimension(:,:,:,:,:,:,:), intent(in) :: V_b, HH_b, c_b,LL_b,BB_b
		real(8), dimension(:,:), intent(out) :: income_star, employment_star
        real(8), dimension(:), intent(out) :: pstar,Lstar_b,Hstar_b,Bstar_b,income_b,employment_b,tax_b,Cstar_b,housing_tax_b
        real(8), dimension(:), intent(out) :: ret_b
		integer, dimension(countP,T_sim-1,nD) :: max_loc_b
        real(8), dimension(:,:), intent(out) :: AggH
        real(8), dimension(:), intent(in) :: AggH_s
        real(8), dimension(countP,T_sim) :: AggR,AggS
        real(8), dimension(:,:),allocatable :: AggD
        integer, dimension(countP,T_sim) :: count_AggR, count_AggD,count_AggH,count_AggS, count_AggHstock, count_AggR_own, count_AggR_rent
		integer, dimension(:,:,:), allocatable :: countD
		integer, dimension(countp) :: count_rent_p, count_buy_p
		real(8), dimension(countp) :: count_buy_p_interp
		integer, dimension(2,T_sim-1) :: psub
		real(8) :: inc
		real(8), dimension(nE) :: edist,zerodist,nonzerodist
		real(8), dimension(nE,2) :: initliqdist, initilldist, initnwdist, initadist
		real(8), dimension(nE) :: initliq, initnw, initill, initashare
		integer, dimension(T_sim,T,nD) :: max_count_age 
		integer, dimension(T_sim,T,nE,nD) :: E_print
		real(8), dimension(T_sim,T,nD) :: B_print
		integer, dimension(nE,countP) :: countE
		integer, dimension(nB,T) :: E_dist
		real(8), dimension(2) :: theta_pstar
		integer, dimension(2) :: ilo_pstar
		integer, dimension(countP) :: count_age
		
		L_sim=0d0
		H_sim=0d0
		B_sim=0d0
		C_sim=0d0
		tax_sim=0d0
		
		rent_buy(:,:,1)=1
		
		
		countD=0
		allocate(q_sim(I,countP,T_sim-1))
		allocate(housing_tax(I,countP,T_sim-1))
		

		allocate(V_sim(I,countP,nD))
		allocate(CountD(T_sim,countP,nD))
		allocate(AggD(countP,T_sim))
		allocate(ExoT_loc(I))
		

        Lstar_b(1)=0d0 ! AggL(1,1)
        Hstar_b(1)=0d0 !AggH(1,1)
        Bstar_b(1)=0d0 !AggB(1,1)
        !q_sim=0d0 


        !allocate(q_sim(T_sim-1,I,countP))



        !adj_cost=0d0
        count_work=0
        ret=0d0
        AggR=0d0
        AggH=0d0
        AggS=0d0
                        AggS(:,1)=2d0*(1d0-delta)+AggH_s
        AggS(:,2)=2d0*(1d0-delta)+AggH_s
        count_AggR=0
        count_AggH=0
        count_AggS=0
		count_AggHstock=0
		max_loc_b=0
		!print*, rent_buy(1:10,:,1)
		
		psub(1,:)=1
		psub(2,:)=countP
		OPEN(3, FILE = './KMV_Input/zydist.txt',status='old',action='read')
				read(3,*) edist(:)
		close(3)
		OPEN(1, FILE = './KMV_Input/initw.txt')
			DO iE=1,nE
				READ(1,*) initadist(iE,1),initadist(iE,2),initashare(iE)
				!initadist(iE,1)=1d0
				!initadist(iE,2)=1.0d0-initadist(iE,1)
			END DO
		CLOSE(1)
		
	do iE=1,nE
		print*, 'zero dist',edist(iE)*initadist(iE,1)
		zerodist(iE)=int(edist(iE)*initadist(iE,1)*I/T)
		print*, 'zero dist',zerodist(iE)
		
		print*, 'non-zero dist',edist(iE)*initadist(iE,2)
		nonzerodist(iE)=int(edist(iE)*initadist(iE,2)*I/T)
		print*, 'non-zero dist',nonzerodist(iE)
	end do
		! Subscript #1 of the array ZERODIST has value 6 which is greater than the upper bound of 5
		count=0
		do while (sum(zerodist)<real(I)/real(T)*sum(edist*initadist(:,1))) 
			count=count+1
			print*,'sum zerodist',sum(zerodist),'count', count, 'zerodistcount',zerodist(count),real(I)/real(T)*sum(edist*initadist(:,1))
			zerodist(count)=zerodist(count)+1
		end do		
		print*, 'zerodist', sum(zerodist(:)),zerodist(:)
		count=0
		do while (sum(nonzerodist)<real(I)/real(T)*sum(edist*initadist(:,2))) 
			count=count+1
			print*,'sum nonzerodist',sum(nonzerodist),'count', count, 'nonzerodistcount',nonzerodist(count),real(I)/real(T)*sum(edist*initadist(:,2))
			nonzerodist(count)=nonzerodist(count)+1
		end do		
		print*, 'zerodist', sum(nonzerodist(:)),nonzerodist(:)

		

	initliqdist=0.0d0
	initilldist=0d0
	initnwdist=0.0d0
	initnw=0.0d0
	OPEN(1, FILE = './KMV_Input/initwbroad.txt',status='old',action='read')
		DO iE=1,nE
			READ(1,*) initliqdist(iE,1),initliqdist(iE,2),initliq(iE),initilldist(iE,1),initilldist(iE,2),initill(iE),initnwdist(iE,1),initnwdist(iE,2),initnw(iE)
		END DO
	CLOSE(1)

	H_sim(iI,iP,1)=0d0

	initliq=initliq/(DataAvAnnualEarns/(Numeraire))
	initill=initill/(DataAvAnnualEarns/(Numeraire))
	initnw=initnw/(DataAvAnnualEarns/(Numeraire))
	

		
        do it=1,T_sim-1
		count_age=0
		countE=0
		ExoT_loc=0
		ExoT_loc(:)=ET_loc(:,it)+nE*(AggT_loc(it)-1)
		
		!$OMP Parallel Default(Shared) private(iP,iI,iD,inc)
		!$OMP DO collapse(2) reduction(+:AggR) reduction(+:AggH) reduction(+:count_work) reduction(+:max_loc_b)
			do iP=1,countP		
				do iI=1,I 
				
					if (age_sim(iI,it)<Jret) then ! This is the next period value we are solving for
						inc=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))))-FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),L_sim(iI,iP,it)*H_sim(iI,iP,it),minval((/0d0,B_sim(iI,iP,it)/)),r_borrow_grid(LambT_loc(it)))
					else
						inc=FnTax(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)))-FnTaxm(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)),L_sim(iI,iP,it)*H_sim(iI,iP,it),minval((/0d0,B_sim(iI,iP,it)/)),r_borrow_grid(LambT_loc(it)))
					end if 
					!print*, 'inc'
                    do iD=1,nD
						if (rent_buy(iI,iP,it)==2 .and. (iD==1 .or. iD==4)) then ! Own to rent or buy

 V_sim(iI,iP,iD)=interp1qnoextrap(B,V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),B_sim(iI,iP,it)&
                        +(1-delta-tauh-cost)*p(iP)*H_sim(iI,iP,it)-(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it)+inc)
 
						else if (rent_buy(iI,iP,it)==2 .and. (iD==2 .or. iD==3)) then
						!	print*,'iI,iD', iI,iD,rent_buy(iI,iP,it), B_sim(iI,iP,it), H_sim(iI,iP,it), L_sim(iI,iP,it)
						
							V_sim(iI,iP,iD)=interp3q(L,H,B,V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),L_sim(iI,iP,it),H_sim(iI,iP,it),B_sim(iI,iP,it)) !,B_sim_star(iI,it))

						else if ((rent_buy(iI,iP,it)==1 .and. (iD==4 .or. iD==1)) .or. (rent_buy(iI,iP,it)==2 .and. iD==5)) then 
							!print*,'iI,iD', iI,iD,rent_buy(iI,iP,it), B_sim(iI,iP,it), H_sim(iI,iP,it), L_sim(iI,iP,it)
							V_sim(iI,iP,iD)=interp1q(B,V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),iD),B_sim(iI,iP,it))
						
							end if 
						
                    end do
					
         
                   if (rent_buy(iI,iP,it)==1 .and. age_sim(iI,it)<T) then
						max_loc_sim(iI,iP,it)=maxloc(V_sim(iI,iP,(/1,4/)),dim=1)
						if (max_loc_sim(iI,iP,it)==2) then
							max_loc_sim(iI,iP,it)=4
						end if
					else if (rent_buy(iI,iP,it)==1 .and. age_sim(iI,it)==T) then
					! Renter's can't buy a house in the final period of life.
						max_loc_sim(iI,iP,it)=4
					else if (rent_buy(iI,iP,it)==2) then
						max_loc_sim(iI,iP,it)=maxloc(V_sim(iI,iP,(/1,2,3,4,5/)),dim=1)   
					end if
					count=0
					if (max_loc_sim(iI,iP,it)==1 .or. max_loc_sim(iI,iP,it)==2 .or. max_loc_sim(iI,iP,it)==3 .or. max_loc_sim(iI,iP,it)==4 .or. max_loc_sim(iI,iP,it)==5) then
						continue
					else 
					count=count+1
						if (count<1000 ) then
					!	print*, 'something wrong with max_loc at (iI,iP,it):', iI,iP,it, max_loc_sim(iI,iP,it), max_loc_sim(iI,iP,it-1),'age,rent/buy,V',age_sim(iI,it),rent_buy(iI,iP,it),rent_buy(iI,iP,it-1),V_sim(iI,iP,:)
						end if 
					end if 
					! Summing across each of the five choices for all agents
						max_loc_b(iP,it,max_loc_sim(iI,iP,it))=max_loc_b(iP,it,max_loc_sim(iI,iP,it))+1

				
                   
					
                    
                    
                    if (rent_buy(iI,iP,it)==2 .and. (max_loc_sim(iI,iP,it)==2 .or. max_loc_sim(iI,iP,it)==3)) then    !Own to pay or refi
                        H_sim(iI,iP,it+1)=H_sim(iI,iP,it)
					
                        L_sim(iI,iP,it+1)=interp3q(L,H,B,LL_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
                        L_sim(iI,iP,it),H_sim(iI,iP,it),B_sim(iI,iP,it))
                        
                        B_sim(iI,iP,it+1)=interp3q(L,H,B,BB_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
                        L_sim(iI,iP,it),H_sim(iI,iP,it),B_sim(iI,iP,it))
                        
						C_sim(iI,iP,it)=interp3q(L,H,B,c_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
						L_sim(iI,iP,it),H_sim(iI,iP,it),B_sim(iI,iP,it))
    
						if (age_sim(iI,it+1)>1) then
							rent_buy(iI,iP,it+1)=2
						else
							rent_buy(iI,iP,it+1)=1
						end if 
                        AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP,it)
						housing_tax(iI,iP,it)=tauh*p(iP)*H_sim(iI,iP,it)
                     else if ((max_loc_sim(iI,iP,it)==1 .or. max_loc_sim(iI,iP,it)==4) .and. rent_buy(iI,iP,it)==2) then          !Own to Own

                        H_sim(iI,iP,it+1)=interp1qnoextrap(B,HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
                        B_sim(iI,iP,it)+(1d0-delta-tauh-cost)*p(iP)*H_sim(iI,iP,it)-(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it)+inc)
                        
                         
                        L_sim(iI,iP,it+1)=interp1qnoextrap(B,LL_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
                        B_sim(iI,iP,it)+(1d0-delta-tauh-cost)*p(iP)*H_sim(iI,iP,it)-(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it)+inc)

                    
                        B_sim(iI,iP,it+1)=interp1qnoextrap(B,BB_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
                        B_sim(iI,iP,it)+(1d0-delta-tauh-cost)*p(iP)*H_sim(iI,iP,it)-(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it)+inc)
						
						
						C_sim(iI,iP,it)=interp1qnoextrap(B,c_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),&
						B_sim(iI,iP,it)+(1-delta-tauh-cost)*p(iP)*H_sim(iI,iP,it)-(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it)+inc)
    
						housing_tax(iI,iP,it)=tauh*p(iP)*H_sim(iI,iP,it)
                       ! rent_buy_p(iI,iP,it)=1
						AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP,it+1)
						!print*, 'own: sell',iP,iI, AggO(iP,it+1), H_sim(iI,iP,it+1)
						
						if (max_loc_sim(iI,iP,it)==1 .and. (age_sim(iI,it)==T .or. age_sim(iI,it+1)==1) ) then
						rent_buy(iI,iP,it+1)=2

                        else if(max_loc_sim(iI,iP,it)==4 ) then
                            rent_buy(iI,iP,it+1)=1 ! Do I change this?

                        else if (max_loc_sim(iI,iP,it)==1 .and. age_sim(iI,it)<T .and. age_sim(iI,it+1)>1) then
                       !     print*, AggS(it+1,iP),H_sim_star(iI,it),mu_own(it,age_sim(iI,it)),H_sim_star(iI,it) *mu_own(it,age_sim(iI,it))
                            rent_buy(iI,iP,it+1)=1

                        end if



                        
                        
                    else if (rent_buy(iI,iP,it)==2 .and. max_loc_sim(iI,iP,it)==5) then
						H_sim(iI,iP,it+1)=R(1)
						rent_buy(iI,iP,it+1)=1
						L_sim(iI,iP,it+1)=0d0
						B_sim(iI,iP,it+1)=interp1qnoextrap(B,BB_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),B_sim(iI,iP,it))
						C_sim(iI,iP,it)=interp1qnoextrap(B,c_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),B_sim(iI,iP,it))
						
						AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP,it+1)
	
					else if (rent_buy(iI,iP,it)==1 .and. (max_loc_sim(iI,iP,it)==1 .or. max_loc_sim(iI,iP,it)==4)) then     !Own to rent or default
                        H_sim(iI,iP,it+1)=interp1qnoextrap(B,HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),B_sim(iI,iP,it))

                        L_sim(iI,iP,it+1)=interp1qnoextrap(B,LL_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),B_sim(iI,iP,it))

                        B_sim(iI,iP,it+1)=interp1qnoextrap(B,BB_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),B_sim(iI,iP,it))


						C_sim(iI,iP,it)=interp1qnoextrap(B,c_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),max_loc_sim(iI,iP,it)),B_sim(iI,iP,it))
    
                   

						AggR(iP,it+1)=AggR(iP,it+1)+H_sim(iI,iP,it+1)

                        if (max_loc_sim(iI,iP,it)==1 .and. age_sim(iI,it+1)>1) then
                            rent_buy(iI,iP,it+1)=2

							
						else if (max_loc_sim(iI,iP,it)==1 .and. age_sim(iI,it+1)==1) then
							rent_buy(iI,iP,it+1)=1

          				else  if (max_loc_sim(iI,iP,it)==4) then
                            rent_buy(iI,iP,it+1)=1 ! Do I change this?x
		
                        end if
						
						

                    end if
				
		

				
					if (it==1) then
						if (rent_buy(iI,iP,it)==1) then
							AggR(iP,it)=AggR(iP,it)+H_sim(iI,iP,it)
						else
							AggH(iP,it)=AggH(iP,it)+H_sim(iI,iP,it)
						end if
						
			
					end if
				
					
					!!! Here income variabels are going to be rescaled so that employment=1
					if (age_sim(iI,it)<Jret) then
						if (max_loc_sim(iI,iP,it)==2 .or. max_loc_sim(iI,iP,it)==3) then
							tax_sim(iI,iP,it)=FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it),minval((/0d0,B_sim(iI,iP,it)/)),r_borrow_grid(lambT_loc(it)))
						else
							tax_sim(iI,iP,it)=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it)))) !+p(iP)*H_sim(iI,iP,it)*tauh
						end if 
					else 
						
						if (max_loc_sim(iI,iP,it)==2 .or. max_loc_sim(iI,iP,it)==3) then
						
							tax_sim(iI,iP,it)=FnTaxm(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)),(1d0+r_borrow_grid(lambT_loc(it)))*L_sim(iI,iP,it)*H_sim(iI,iP,it),minval((/0d0,B_sim(iI,iP,it)/)),r_borrow_grid(lambT_loc(it)))
							
						else 
							tax_sim(iI,iP,it)=FnTax(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y))) !+p(iP)*H_sim(iI,iP,it)*tauh
						end if 
						
					end if    
			
					if (iI<10 )then
					!	print*, it, iP,iI, tax_sim(iI,iP,it),c_sim(iI,iP,it)
					end if 
			
				if (age_sim(iI,it)==T) then
					count_age(iP)=count_age(ip)+1
					countE(ET_loc(iI,it+1),iP)=countE(ET_loc(iI,it+1),iP)+1
					if (count_age(iP)<=sum(Zerodist) .and. countE(ET_loc(iI,it+1),iP)<=Zerodist(ET_loc(iI,it+1))) then
				
						b_sim(iI,iP,it+1)=0
				
					else if (count_age(iP)>sum(Zerodist) .and. countE(ET_loc(iI,it+1),iP)<=nonZerodist(ET_loc(iI,it+1)) ) then
						b_sim(iI,iP,it+1)=initnw(ET_loc(iI,it+1))
					end if 
					
				end if
				
				end do
			end do
		!$OMP End Do nowait
		!$OMP End parallel
		
				!end do
		!do it=1,T_sim-1
			do iP=1,countP   
            
			count_rent_p(iP)=sum((-1)*(rent_buy(:,iP,it)-2))
			count_buy_p(iP)=sum(rent_buy(:,iP,it)-1)

		
			
			if(count_buy_p(iP)>0) then 
				!AggL(iP,it+1)=sum(L_sim(:,iP,it+1))/count_buy_p(iP)
				AggH(iP,it+1)=AggH(iP,it+1)/count_buy_p(iP)
			else
				!AggL(iP,it+1)=0
				AggH(iP,it+1)=0
			end if
			
			if (count_rent_p(iP)>0) then
			
				AggR(iP,it+1)=AggR(iP,it+1)/count_rent_p(iP)
			else 
				AggR(iP,it+1)=0
			end if 
	


 

		end do
		!AggR(:,it+1)=AggR_own(:,it+1)+AggR_rent(:,it+1)

            
            Aggd(:,it+1)=real(count_buy_p)/real(I)*AggH(:,it+1)+real(count_rent_p)/real(I)*AggR(:,it+1) !-(1d0-delta)*AggR(:,it) 
			if (it>1) then
				AggS(:,it+1)=interp1q(p,AggS(:,it),pstar(it-1))*(1d0-delta)+AggH_s 
            end if 

			psub(1,it)=1
			psub(2,it)=countP


		 
		    if(AggD(psub(2,it),it+1) .gt. AggS(psub(2,it),it+1)) then
				pstar(it)=p(psub(2,it))
			elseif(AggD(psub(1,it),it+1) .lt. AggS(psub(1,it),it+1)) then
				pstar(it)=p(1)
			else
				call sub_bisect_excess_demand1(pstar(it),excess_demand1,p,Aggd(:,it+1),AggS(:,it+1),p(psub(1,it)),p(psub(2,it)))
			endif
			call base_fun(p,pstar(it),ilo_pstar(:),theta_pstar(:))
            !ilo_pstar=bsearch(pstar(it),p)

			Select Case (display)
				Case(0)
					print*, it,AggT_loc(it), YT_loc(it),pstar(it) !, 'AggL',AggL(it+1,:), psub,psub3
					!print*, 'AggB',AggB(:,it+1)
					print*, 'AggD', Aggd(:,it+1)
					!print*,	'AggL', AggL(:,it+1) 
					print*, 'AggS', AggS(:,it+1)
					print*, 'AggR_t+1', AggR(:,it+1)
					print*, 'Fraction Homeowner', real(count_buy_p)/real(I)
					print*, 'Fraction renter', real(count_rent_p)/real(I)
					print*, 'Fraction total', real(count_rent_p)/real(I)+real(count_buy_p)/real(I)
					print*, 'Percent buy/sell', real(max_loc_b(:,it,1))/real(I)*100
					print*, 'Percent stay', real(max_loc_b(:,it,2))/real(I)*100
					print*, 'Percent refi', real(max_loc_b(:,it,3))/real(I)*100
					print*, 'Percent rent', real(max_loc_b(:,it,4))/real(I)*100
					print*, 'Percent default', real(max_loc_b(:,it,5))/real(I)*100
					print*, 'Percent total', real(sum(max_loc_b(:,it,:),dim=2))/real(I)*100
					
				
					
         !!   print*, 'own', sum(rent_buy_p(:,1)-1),sum(rent_buy_p(:,2)-1),sum(rent_buy_p(:,3)-1),sum(rent_buy_p(:,4)-1),sum(rent_buy_p(:,5)-1)
		 
			case(1)
				continue
			end select
			
			select case(display)	
	case(1)
			print*, it,pstar(it)
		
	case(0)
		continue
end select

      
			
            end do
		!OMP Parallel Default(Shared) private(iI,it)
		!OMP DO collapse(2) reduction(+:count_work) 
		do it=1,T_sim-1
			do iI=1,I
				if (age_sim(iI,it)<Jret) then
					income_star(iI,it)=income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it)))
					employment_star(iI,it)=income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it)))/Y(YT_loc(it))
					count_work(it)=count_work(it)+1     
				else 
					ret(iI,it)=income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y))
					income_star(iI,it)=ret(iI,it)

				end if    
				
			end do

			employment_b(it)=sum(employment_star(:,it))/real(count_work(it))
			ret_b(it)=sum(ret(:,it))/real(I-count_work(it))
			
			!print*, sum(employment_star(:,it)),count_work(it), employment_b(it)
              
			! Taking averages of income variables. Will wait to do so for employment
			income_b(it)=sum(income_star(:,it))/real(I)
			do iP=1,nP
				count_buy_p_interp(iP)=real(sum(rent_buy(:,iP,it)-1))
			end do
         
			if (interp1q(p,count_buy_p_interp(:),pstar(it))==0) then
				housing_tax_b(it)=0d0
			else 
				housing_tax_b(it)=interp1q(p,sum(housing_tax(:,:,it),dim=1),pstar(it))/(interp1q(p,count_buy_p_interp(:),pstar(it))) !real(count_buy_p(bsearch(pstar(it),p))+count_buy_p(bsearch(pstar(it),p)+1))*0.5d0)
			end if 
			tax_b(it)=interp1q(p,sum(tax_sim(:,:,it),dim=1),pstar(it))/real(I)
			Cstar_b(it)=interp1q(P,sum(C_sim(:,:,it), dim=1),pstar(it))/real(I)

				
		end do   
	   !OMP END Do
	   !OMP End Parallel 

         !print*, maxval(H_sim), maxval(L_sim)
       
    end subroutine sub_simulate_steady
end module mod_simulate_steady
