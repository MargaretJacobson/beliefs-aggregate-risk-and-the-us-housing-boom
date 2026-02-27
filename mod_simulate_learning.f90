!       This file contains the simulation code for Krussel-Smith Simulation

module mod_simulate_learning
    use omp_lib
    use parameters
	use mod_functions
    use mod_matlab
    use mod_interp
    use mod_root1dim
    implicit none
contains
    
    subroutine sub_simulate_learning(age_sim,rent_buy,L,B,Bs,H,R,Y,chi,V_b,LL_b,BB_b,HH_b,c_b,AggT_loc,YT_loc,LambT_loc,loge,ET_loc,&
	AggH,AggS,AggH_s,p,countP,pstar,Lstar_b,Bstar_b,Hstar_b,Hstar_s,Cstar_b,income_b,C_sim_star,&
	income_star,L_sim_star,B_sim_star,H_sim_star,tax_b,tax_star,housing_tax_b,employment_b,employment_star,&
	ret_b,max_loc_sim_star,max_loc_b,countD,a_out,nA,pprime,sim_end,count,preg,countZ,pregZ,homeownership,pprime_solve,rental_sim,rental,a_out_all,q)
    implicit none
    integer :: it,iI,iP, iY,iYY, count_default, Yreveal,iX,iXX,Aloc,iM,psub,psub2,ind(2),ind_old(2),ilamb,iA,iE,boom_count,iXodd,iXXodd,rent,buy,iAggAgg
    integer :: psub3,psub4,iD,iAgg,up,last(nY*nlamb*nphI,nY*nlamb*nphI),lastlast(nY*nlamb*nphI,nY*nlamb*nphI)
    integer, intent(in) :: countP,nA,sim_end
	integer, allocatable :: ExoT_loc(:)
	real(8), allocatable :: Wealth(:,:), Wealth_star(:)
    integer, dimension(:,:), intent(inout) :: age_sim, rent_buy
    real(8) :: p0(nY*nlamb), p1(nY*nlamb),xvec(nY*nlamb), gain(nY*nlamb,nY*nlamb),det(nY*nlamb)
    real(8), dimension(nY*nlamb) :: pred,pred2
    real(8), allocatable :: H_sim(:,:) !,B_sim,L_sim
    real(8), allocatable :: V_sim(:,:,:)
    real(8), allocatable :: V_sim_star(:,:,:)
	integer, allocatable :: max_loc_sim(:,:)
    integer, dimension(:,:), intent(out) :: max_loc_sim_star
	integer, dimension(nP,T_sim,nD) :: max_loc_p
    integer, dimension(:,:), intent(out) :: max_loc_b
	integer, dimension(T_sim-1) :: count_work
    real(8), allocatable :: q_sim_star(:,:), ret(:,:)
    real(8), dimension(:), intent(in) :: L,H,p,Y,chi,R,B,Bs
    real(8), dimension(:,:), intent(in) :: loge
	real(8), dimension(:,:), intent(inout) :: pregZ
    real(8), dimension(nY*nlamb,mom+1,countP) :: a_new,a_new2
    real(8), dimension(mom+1,mom+1,nY*nlamb) :: Xt_prod
    real(8), dimension(mom+1,mom+1,nY*nlamb) :: Rt_new
    real(8), dimension(nY*nlamb,nY*nlamb,T_sim-1,mom+1,mom+1) :: Rt
    real(8), dimension(nY*nlamb,nY*nlamb,mom+1) :: acoeff, acoeff_solve, acoeff2, acoeff_solve2
	real(8), dimension(nY*nlamb,countP) :: pprime_old
    integer, dimension(:,:), intent(inout) :: count
    integer, dimension(:), intent(in) :: YT_loc,AggT_loc,LambT_loc
    integer, dimension(:,:), intent(in) :: ET_loc
	real(8), dimension(:,:,:,:), intent(in) :: rental
	real(8), dimension(:), intent(inout):: rental_sim
	real(8), dimension(:,:,:,:,:,:,:,:), intent(in) :: q
    real(8), dimension(countP) :: Aggd
    real(8), dimension(mom+1) :: vec
    real(8), dimension(:,:),intent(out) :: pprime_solve
    real(8), dimension(2,nY*nlamb*nphi) ::	pprime_solve2
	real(8), dimension(:,:,:,:,:,:,:,:,:), intent(in) :: V_b,HH_b,c_b,LL_b,BB_b
    real(8), dimension(:,:), intent(out) :: C_sim_star, income_star, L_sim_star, H_sim_star,B_sim_star,tax_star,employment_star
    real(8), dimension(:,:), intent(out) :: a_out
	real(8), dimension(:,:,:), intent(out) :: a_out_all
    real(8), dimension(T_sim-1), intent(out) :: pstar,Lstar_b,Hstar_b,Bstar_b,Cstar_b,income_b,Hstar_s,tax_b,employment_b,housing_tax_b,homeownership
    real(8), dimension(:), intent(out) :: ret_b
    real(8), dimension(:,:), intent(inout) :: AggH, AggS
    real(8), dimension(:), intent(in) :: AggH_s
	real(8), dimension(countP,T_sim) :: AggR
	integer :: count_buy_p, count_rent_p, Eindex,count_age,count1
	integer, dimension(nE) :: countE
	real(8), dimension(nE) :: edist,zerodist,nonzerodist
	real(8), dimension(nE,2) :: initliqdist, initilldist, initnwdist, initadist
	real(8), dimension(nE) :: initliq, initnw, initill, initashare
    integer, dimension(countP,T_sim) :: count_AggR_own, count_AggR_rent, count_AggD,count_AggH,count_AggS, count_AggR
    integer, dimension(2) :: ilo_bstar, ilo_hstar, ilo_pstar
	real(8), dimension(2) :: theta_bstar, theta_hstar, theta_pstar
    real(8), dimension(nY*nlamb,nY*nlamb,T_sim-1,2) :: Xreg
	real(8), dimension(nY*nlamb,nY*nlamb,T_sim,countP) :: Yreg
    real(8), dimension(:,:,:),intent(inout) :: preg
    real(8), dimension(nY*nlamb,mom+1,nA) :: acom
    real(8), dimension(nY*nlamb,countP) :: pprime_new, pprime_new2
    real(8), dimension(:,:,:,:), intent(in) :: pprime   
	integer, dimension(:,:),intent(out) :: countD
	real(8), dimension(:,:), allocatable :: housing_tax		
	integer, dimension(T_sim,T,nD) :: max_count_age 
	real(8) :: inc, years
	integer, dimension(:), intent(inout) :: countz
	if (buy_grid==0) then
			rent=4
			buy=1
		else
			rent=7
			buy=6
		end if 

	allocate(V_sim_star(I,T_sim-1,nD+2))
	allocate(q_sim_star(I,T_sim-1))
	allocate(ret(I,T_sim-1))
	allocate(housing_tax(I,T_sim-1))
	allocate(V_sim(I,countP,nD+2))
	Allocate(ExoT_loc(I))
	allocate(H_sim(I,countP))
	allocate(max_loc_sim(I,countP))
	allocate(Wealth(I,nP))
	allocate(Wealth_star(I))	
	
	V_sim_star=death
	V_sim=death

    q_sim_star=0d0

	housing_tax=0d0
    tax_star=0d0
    employment_star=0d0
    count_work=0
    ret=0d0
    max_loc_b=0
	count_AggR=0
    count_AggH=0
    count_AggS=0
	countD=0
	max_loc_sim_star=0
	max_loc_b=0
	count_age=0
	max_count_age=0
	Rt=0d0	
	rental_sim=0d0


    vec=1d0
    Aloc=(nA+1)/2

	max_loc_p=0
        
    print*, 'Mod Simulate Learning'
    Xreg=1d0
    tax_star=0d0
    employment_star=0d0
    gain=gain_param
    acom=0d0
	C_sim_star=0d0
	
	! Reading in bequest distributions
		
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
		print*, 'zero dist',edist(iE)*initadist(iE,1)
		zerodist(iE)=int(edist(iE)*initadist(iE,1)*I/T)
		print*, 'zero dist',zerodist(iE)
		
		print*, 'non-zero dist',edist(iE)*initadist(iE,2)
		nonzerodist(iE)=int(edist(iE)*initadist(iE,2)*I/T)
		print*, 'non-zero dist',nonzerodist(iE)
	end do
	count_age=0
		do while (sum(zerodist)<real(I)/real(T)*sum(edist*initadist(:,1))) 
			count_age=count_age+1
			zerodist(count_age)=zerodist(count_age)+1
		end do	
		zerodist(1:4)=zerodist(1:4)+1		
		print*, 'zerodist', sum(zerodist(:)),zerodist(:)
		count_age=0
		do while (sum(nonzerodist)<real(I)/real(T)*sum(edist*initadist(:,2))) 
			count_age=count_age+1
			nonzerodist(count_age)=nonzerodist(count_age)+1
		end do		
		print*, 'non-zero zerodist', sum(nonzerodist(:)),nonzerodist(:)

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
	!print*, initnw
	
	
	! Reading in coefficients 
	if (nlamb==2) then
		open(3,file='coefficients2.csv',status='old',action='read')
		do iM=1,mom+1
			do iY=1,nY*nlamb
				do iYY=1,nY*nlamb
					read(3,*) acoeff(iY,iYY,iM)
				end do
			end do
		end do
		close(3)
	else
		open(3,file='coefficients1.csv',status='old',action='read')
		do iM=1,mom+1    
			do iY=1,nY*nlamb
				do iYY=1,nY*nlamb
					
					read(3,*) acoeff(iY,iYY,iM)

				end do
			end do
		end do
		close(3)
	end if
		
	do iY=1,nY*nlamb
		do iYY=1,nY*nlamb
			print*, 'iY,iYY,acoeff',iY,iYY,acoeff(iY,iYY,:)
		end do
	end do
		
	!acoeff(:,:,2)=acoeff(:,:,2)+0.0133d0

   acoeff2=acoeff


         
	! Filling in coefficients before learnign so they print correctly.	 
	do it=1,sim_end-1 
		a_out(it,:)=acoeff(aggT_loc(it),aggT_loc(it+1),:)
		a_out_all(it,:,1)=acoeff(aggT_loc(it),:,1)
		a_out_all(it,:,2)=acoeff(aggT_loc(it),:,2)
	end do 
    
	! Agents nearly forget observations that are 30 years old. Can go back to the 1970s for data
	! Agents recover from recessions in 1970, 1973, 1980, 1981, and 1990 so count(4,3)=5
	
	! There are 4 recessions in this sample. 1973=15 months, 1980=6 months, 1981=16 months 1990=8 months
	! I will say 2+1+2+1=6 years of recessions.
	years=(200d0/2d0)
	count(4,3)=1+0*(0.1d0/(2d0))*years
	count(4,4)=1+0*(0.9d0/(2d0))*years
	count(4,2)=1
	count(4,1)=1
	
	! The rest of the 30/2=15 periods are expansions
	count(3,4)=1+0*(0.1d0/(2d0))*years
	count(3,3)=1+0*(0.9d0/(2d0))*years
	count(3,1)=1+ 0*1
	count(3,2)=1
	
	last=0
	
	! used last period's beliefs a_{t-1} to forecast house prices last period so it's really a_{Z_{t-1},Z_{t}}
	! try aggT_loc(it) instead of iYY, Yreveal but this is different from a_{t-1} which is when you are going to get
	! the wrong forecast error in the boom. But the forecast error will be larger in the bust.
	do iAgg=3,4
		do iAggAgg=3,4
			up=sim_end
			do while (last(iAgg,iAggAgg)==0) 
				if (AggT_loc(up)==iAgg .and. AggT_loc(up+1)==iAggAgg) then
					last(iAgg,iAggAgg)=up
				else
					up=up-1
				end if
			end do
			print*, 'iAgg, up, last', iAgg, iAggAgg, up, last
		end do
		last(iAgg,1:2)=last(iAgg,3:4)
	end do
	
	do iAggAgg=1,4
		if (last(1,iAggAgg)==0) then
			last(1,iAggAgg)=last(3,iAggAgg)
		end if
		if(last(2,iAggAgg)==0) then
			last(2,iAggAgg)=last(4,iAggAgg)
		end if
	end do
	print*, 'last',last
	
	lastlast=0
	do iAgg=3,4
		do iAggAgg=3,4
			up=last(iAgg,iAggAgg)-1
			do while (lastlast(iAgg,iAggAgg)==0)
		
				if (AggT_loc(up)==iAgg .and. AggT_loc(up+1)==iAggAgg) then
					lastlast(iAgg,iAggAgg)=up
				else
					up=up-1
				end if
				
			end do
			print*, 'iAgg, up, last', iAgg, iAggAgg, up
			print*, 'last', last
			print*, 'lastlast',lastlast
		end do
		lastlast(iAgg,1:2)=lastlast(iAgg,3:4)
	end do
	
	do iAggAgg=1,4
		if (lastlast(1,iAggAgg)==0) then
			lastlast(1,iAggAgg)=lastlast(3,iAggAgg)
		end if
		if(lastlast(2,iAggAgg)==0) then
			lastlast(2,iAggAgg)=lastlast(4,iAggAgg)
		end if
	end do
	print*, 'last',last
	
    do it=sim_end,T_sim-1
		ExoT_loc=0
		ExoT_loc(:)=ET_loc(:,it)+nE*(AggT_loc(it)-1)
        Rt_new=0d0
        pred=0d0
        A_new=0d0
        ret=0
		boom_count=1
		AggH(:,it+1)=0d0
		AggR(:,it+1)=0d0
		Wealth=0d0
		Wealth_star=0d0
		! This dampens oscillations later in the simultion.
		if (it>=490 .and. it<=494) then
			boom_count=2
		end if 
		countZ(AggT_loc(it))=countZ(AggT_loc(it))+1
		
		

		
		
	
		
		! Counting each occurance of Z 
		print*, 'countZ', countZ(AggT_loc(it))
		
		count(AggT_loc(it),:)=count(AggT_loc(it),:)+1
		
        
        ! least squares or constant gain learning    
		if (AggT_loc(it)>2 ) then  !.or. boom_count>1 ) then 
			do iYY=1,ny*nlamb
				if (AggT_loc(it-1)==1) then
					gain(AggT_loc(it),iYY)=gain(AggT_loc(it),iYY)*(1d0/(5d0*gain(AggT_loc(it),iYY)+1d0)) 
				else 
					gain(AggT_loc(it),iYY)=gain(AggT_loc(it),iYY)/(gain(AggT_loc(it),iYY)+1d0) !1d0/real(count(AggT_loc(it),4)-1) 
				
				end if 
			end do

		else 
			gain(AggT_loc(it),:)=gain_param
		end if 
		
			
			!print*, it, iY,iYY,count(iY,iYY)
            print*,'it, AggT_loc(it),iYY,count(iY,iYY)', it, AggT_loc(it), count(AggT_loc(it),:)
			print*, 'gain', gain(AggT_loc(it),:)
			print*, it, learn
		!	print*, 'Probabilities', FYY(AggT_loc(it),:)
			
			!print*, 'coeff', acoeff(AggT_loc(it),iYY,:)
		print*, 'it, last', last
			
		if (AggT_loc(it)<3) then
			ind(1)=1
			ind(2)=2
		else
			ind(1)=3
			ind(2)=4
		end if 
		if (AggT_loc(it-2)<3) then
			ind_old(1)=1
			ind_old(2)=2
		else
			ind_old(1)=3
			ind_old(2)=4
		end if 
		pred=0d0		
		do iYY=1,4 !1,2 !ind,ind+1
			! p0=p_{\mathcal{Z}_{t-1}} equals the last equilibrium price from {Z,Z'}
			! Here we know that Z=Agg_loc(it), but we need to compute this for both future values of Z' (Z'={3,4} or Z'={1,2})
			! pred is the predicted value of p_t based on the coefficients from the last occurance of {Z,Z'}
			! When starting the housing boom \mathcal{Z}={1,2} is set to \mathcal{3,4} values.
            if (it==learn) then
				p0(iYY)=log(pstar(it-2)) !
				p1(iYY)=log(pstar(it-1))
				pred(iYY)=acoeff(AggT_loc(it),iYY,1)+acoeff(AggT_loc(it),iYY,2)*p0(iYY) !a_out_all(last(AggT_loc(it),iYY),iYY,1)+a_out_all(last(AggT_loc(it),iYY),iYY,2)*p0(iYY)
				
				pred2(iYY)=acoeff2(AggT_loc(it),iYY,1)+acoeff2(AggT_loc(it),iYY,2)*p0(iYY)
			else if (it==learn+1) then
				p0(iYY)=log(pstar(it-2))
				p1(iYY)=log(pstar(it-1)) 
				
				pred(iYY)=acoeff(AggT_loc(it),iYY,1)+acoeff(AggT_loc(it),iYY,2)*p0(iYY) 
				pred2(iYY)=acoeff2(AggT_loc(it),iYY,1)+acoeff2(AggT_loc(it),iYY,2)*p0(iYY) 
				
            else
				p0(iYY)=log(pstar(it-2)) !log(pstar(lastlast(AggT_loc(it),iYY))) !log(preg(AggT_loc(it),iYY,count(AggT_loc(it),iYY)-2)) !
				p1(iYY)=log(pstar(it-1)) !log(pstar(last(AggT_loc(it),iYY))) !log(preg(AggT_loc(it),iYY,count(AggT_loc(it),iYY)-1)) !
				
				pred(iYY)=acoeff(AggT_loc(it),iYY,1)+acoeff(AggT_loc(it),iYY,2)*p0(iYY) 
				pred2(iYY)=acoeff2(AggT_loc(it),iYY,1)+acoeff2(AggT_loc(it),iYY,2)*p0(iYY) 
			end if
				! NEed to re-run because of learn==200
					!new is the first version for periods t>learn and new2 is the same version with t==learn. These are exactly the same
					!new3 is the second version for periods t>learn and new4 is for t==learn. new3 works and is similar to new 5
					!new5 is the second version without lagged coefficients and works. new6 doesn't work
					!running new
				! Computing X_t^2 values 
            Xt_prod(1,1,iYY)=1d0
            Xt_prod(1,2,iYY)=p0(iYY)
            Xt_prod(2,1,iYY)=p0(iYY)
            Xt_prod(2,2,iYY)=p0(iYY)*p0(iYY) 
            print*, iYY
			print*, a_out(it-2,:)
			print*, 'acoeff1',  acoeff(AggT_loc(it),iYY,:)
			print*, 'acoeff2',  acoeff2(AggT_loc(it),iYY,:)
			print*, exp(p0(iYY))
			print*, exp(pred(iYY))
            print*, 'iYY, acoeff, p0', iYY, a_out(it-1,:), acoeff(AggT_loc(it),iYY,:),exp(p0(iYY)),exp(pred(iYY))
			 print*, 'iYY, acoeff, p0', iYY, a_out(it-1,:), acoeff(AggT_loc(it),iYY,:),exp(p0(iYY)),exp(pred(iYY))
			!print*, aggT_loc(it-2),ind_old(iYY), ind(iYY), pred(ind(iYY)),FYY(aggT_loc(it-2),ind_old(iYY)),p0(ind(iYY)),(a_out(it-2,1)+a_out(it-2,2)*p0(ind(iYY)))

			 
			
				Rt_new(1,:,iYY)=(/1d0,0d0/) !*1d-2
				Rt_new(2,:,iYY)=(/0d0,1d0/) !*1d-2

			! uncomment for normalization step 
			!Rt_new(:,:,iYY)=Rt(AggT_loc(it),iYY,count(AggT_loc(it),iYY)-1,:,:)+gain(AggT_loc(it),iYY)*(Xt_prod(:,:,iYY)-Rt(AggT_loc(it),iYY,count(AggT_loc(it),iYY)-1,:,:))
            
		
			
			if (gain(AggT_loc(it),iYY)==0) then
				det(iYY)=0d0
			else 
				det(iYY)=(1d0/(Rt_new(1,1,iYY)*Rt_new(2,2,iYY)-Rt_new(1,2,iYY)*Rt_new(2,1,iYY)))
			end if
		
        end do
        print*, 'ind', ind     
        print*, it, 'count', count(AggT_loc(it),:)
		print*, it, 'pred', exp(pred(:))
        print*, it, 'acoeff iYY=1', acoeff(AggT_loc(it),ind(1),:)
        print*, it, 'acoeff iYY=2', acoeff(AggT_loc(it),ind(2),:)
		print*, it, 'acoeff iYY=1', acoeff2(AggT_loc(it),ind(1),:)
        print*, it, 'acoeff iYY=2', acoeff2(AggT_loc(it),ind(2),:)
		print*, it, 'aout_all, IY==1', a_out_all(last(AggT_loc(it),ind(1)),ind(1),:)
		print*, it, 'aout_all, IY==1', a_out_all(last(AggT_loc(it),ind(2)),ind(2),:)
        print*, 'forecast error1',P1(ind(1):ind(2))-pred(ind(1):ind(2))
		print*, 'forecast error2',P1(ind(1):ind(2))-pred2(ind(1):ind(2))

        do iP=1,countP
			do iYY=ind(1),ind(2)   
				! If p_{\mathcal{Z}_t} is not yet determined in equilibrium the "new" coefficients are solved on the grid for prices p.
				! They are compared against the predicted values pred computed above from the last regime house price p0. If p is already determined then
				! this loop step is redundant.
				!xvec=1d0 for simple learning, pred(iYY) for stochatic gradient or recursive least squares
					
				xvec(iYY)=p0(iYY) 
				a_new(iYY,1,iP)=gain(AggT_loc(it),iYY)*det(iYY)*(Rt_new(2,2,iYY)-Rt_new(1,2,iYY)*xvec(iYY))*(P1(iYY)-pred(iYY))+ acoeff(AggT_loc(it),iYY,1)!a_out_all(last(AggT_loc(it),iYY),iYY,1) !  
				a_new(iYY,2,iP)=(gain(AggT_loc(it),iYY)*det(iYY)*(-Rt_new(2,1,iYY)+Rt_new(1,1,iYY)*xvec(iYY)))*(P1(iYY)-pred(iYY))+acoeff(AggT_loc(it),iYY,2) !a_out_all(last(AggT_loc(it),iYY),iYY,2) !
			
				a_new2(iYY,1,iP)=gain2*gain(AggT_loc(it),iYY)*det(iYY)*(Rt_new(2,2,iYY)-Rt_new(1,2,iYY)*xvec(iYY))*(P1(iYY)-pred2(iYY))+ acoeff2(AggT_loc(it),iYY,1)!a_out_all(last(AggT_loc(it),iYY),iYY,1) !  
				a_new2(iYY,2,iP)=gain2*(gain(AggT_loc(it),iYY)*det(iYY)*(-Rt_new(2,1,iYY)+Rt_new(1,1,iYY)*xvec(iYY)))*(P1(iYY)-pred2(iYY))+acoeff2(AggT_loc(it),iYY,2) !a_out_all(last(AggT_loc(it),iYY),iYY,2) !
							
				if (iYY==1 .or. iYY==3) then 
					pprime_new(1,iP)=exp(a_new(iYY,1,iP)+a_new(iYY,2,iP)*log(P(iP)))
					pprime_new2(1,iP)=exp(a_new2(iYY,1,iP)+a_new2(iYY,2,iP)*log(P(iP)))

				else
					pprime_new(2,iP)=exp(a_new(iYY,1,iP)+a_new(iYY,2,iP)*log(P(iP)))
					pprime_new2(2,iP)=exp(a_new2(iYY,1,iP)+a_new2(iYY,2,iP)*log(P(iP)))

				end if 
				print*, 'pprime_new, a1', iP,iYY, pprime_new(iYY,iP), a_new(iYY,:,iP)
				print*, 'pprime_new, a2', iP,iYY, pprime_new2(iYY,iP), a_new2(iYY,:,iP)

			end do     
              
		end do

			ilamb=(1-perfect_corr)*ceiling(real(AggT_loc(it))/real(nY))+perfect_corr*AggT_loc(it)
			iY=(1-perfect_corr)*(floor(real(AggT_loc(it))/(real(nY)*real(ilamb)))+1)+perfect_corr*AggT_loc(it)
			

		do iP=1,countP
            do iI=1,I 

					
				if (age_sim(iI,it)<Jret) then ! This is the next period value we are solving for
					inc=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))))-FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))
				else
					inc=FnTax(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)))-FnTaxm(income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y)),L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))
				end if 
				if (age_sim(iI,it)>1 ) then 
					Wealth(iI,iP)=B_sim_star(iI,it)+(1d0-delta-tauh-cost)*p(iP)*H_sim_star(iI,it)-(1d0+r_borrow_grid(LambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it)
				else
					Wealth(iI,iP)=B_sim_star(iI,it)
		
				end if 
			
				do iD=1,nD+2
					if (rent_buy(iI,it)==2 .and. (iD==6 .or. iD==7)) then ! Own to rent or buy
if (Wealth(iI,iP)+inc>=0d0) then 
 V_sim(iI,iP,iD)=interp5qnoextrap(L,H,B,pprime(iY,1,iP,:),pprime(iY,2,iP,:),V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,iD),L_sim_star(iI,it),H_sim_star(iI,it),Wealth(iI,iP)+inc,pprime_new2(1,iP),pprime_new2(2,iP)) 
else
 V_sim(iI,iP,iD)=interp5qnoextrap(L,H,Bs,pprime(iY,1,iP,:),pprime(iY,2,iP,:),V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,iD),L_sim_star(iI,it),H_sim_star(iI,it),Wealth(iI,iP)+inc,pprime_new2(1,iP),pprime_new2(2,iP)) 

end if 
					else if (rent_buy(iI,it)==2 .and. (iD==2 .or. iD==3) .or. (rent_buy(iI,it)==2 .and. iD==5) ) then
							V_sim(iI,iP,iD)=interp5qnoextrap(L,H,B,pprime(iY,1,iP,:),pprime(iY,2,iP,:),V_b(:,:,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,iD),L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_new2(1,iP),pprime_new2(2,iP))
					else if ((rent_buy(iI,it)==1 .and. (iD==4 .or. iD==1))  ) then 
 
								!if (iY==1 .or. iY==2) then
							V_sim(iI,iP,iD)=interp3qnoextrap(B,pprime(iY,1,iP,:),pprime(iY,2,iP,:),V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,iD),B_sim_star(iI,it),pprime_new(1,iP),pprime_new(2,iP))
							
							
					else if ((rent_buy(iI,it)==2 .and. (iD==4 .or. iD==1))  ) then 
 
								!if (iY==1 .or. iY==2) then
						if (Wealth(iI,iP)+inc>0d0) then ! B(1)) then
						V_sim(iI,iP,iD)=interp3qnoextrap(B,pprime(iY,1,iP,:),pprime(iY,2,iP,:),V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,iD),Wealth(iI,iP)+inc,pprime_new2(1,iP),pprime_new2(2,iP))
							else 
						V_sim(iI,iP,iD)=interp3qnoextrap(Bs,pprime(iY,1,iP,:),pprime(iY,2,iP,:),V_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,iD),Wealth(iI,iP)+inc,pprime_new2(1,iP),pprime_new2(2,iP))
							
							end if 
								
					end if 
				end do

                if (rent_buy(iI,it)==1 .and. age_sim(iI,it)<T) then
                    max_loc_sim(iI,iP)=maxloc(V_sim(iI,iP,(/1,4/)),dim=1)*3-2 !sets to indexes 1 and 4

				else if (rent_buy(iI,it)==1 .and. age_sim(iI,it)==T) then
					! Renter's can't buy a house in the final period of life.
					max_loc_sim(iI,iP)=4
				else if (rent_buy(iI,it)==2 .and. age_sim(iI,it)<T) then
                    max_loc_sim(iI,iP)=maxloc((/V_sim(iI,iP,buy),V_sim(iI,iP,2:3),V_sim(iI,iP,rent),V_sim(iI,iP,5)/),dim=1)
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
						H_sim(iI,iP)=interp3qnoextrap(Bs,pprime(iY,1,iP,:),pprime(iY,2,iP,:),HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,max_loc_sim(iI,iP)),&
						Wealth(iI,iP)+inc,pprime_new2(1,iP),pprime_new2(2,iP))
					else
						H_sim(iI,iP)=interp3qnoextrap(B,pprime(iY,1,iP,:),pprime(iY,2,iP,:),HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,max_loc_sim(iI,iP)),&
						Wealth(iI,iP)+inc,pprime_new2(1,iP),pprime_new2(2,iP))					
					
					end if 

					AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP)

					! Default or rent-rent and rent-buy
                else if (rent_buy(iI,it)==2 .and. max_loc_sim(iI,iP)==5) then
					H_sim(iI,iP)=R(1)

					AggH(iP,it+1)=AggH(iP,it+1)+H_sim(iI,iP)
     
				else if (rent_buy(iI,it)==1 .and. (max_loc_sim(iI,iP)==1 .or. max_loc_sim(iI,iP)==4)) then     !Own to rent or default

					
							H_sim(iI,iP)=interp3qnoextrap(B,pprime(iY,1,iP,:),pprime(iY,2,iP,:),HH_b(1,1,:,ExoT_loc(iI),iP,age_sim(iI,it),:,:,max_loc_sim(iI,iP)),&
							B_sim_star(iI,it), pprime_new(1,iP),pprime_new(2,iP))
                            AggR(iP,it+1)=AggR(iP,it+1)+H_sim(iI,iP) 
                            count_AggH(iP,it+1)=count_AggH(iP,it+1)+1
  
            

                end if

					
						max_loc_p(iP,it,max_loc_sim(iI,iP))=max_loc_p(iP,it,max_loc_sim(iI,iP))+1
                end do
			end do

					count_age=0
				Eindex=1
				countE=0
count_rent_p=sum((-1)*(rent_buy(:,it)-2))
			count_buy_p=sum(rent_buy(:,it)-1)
				
			!print*, 'counts total',it,iP, count_rent_p, count_buy_p, count_rent_p+count_buy_p

			
			do iP=1,nP
				if (count_buy_p>0) then
					!AggL(iP,it+1)=sum(L_sim(:,iP))/count_buy_p
					AggH(iP,it+1)=AggH(iP,it+1)/count_buy_p!/count_AggH(iP,it+1) !/I
				else
					AggH(iP,it+1)=0
					!AggL(iP,it+1)=0
				end if

				if (count_rent_p>0) then
					AggR(iP,it+1)=AggR(iP,it+1)/count_rent_p
				else
					AggR(iP,it+1)=0
				end if 
    
			end do

			
            Aggd=real(count_buy_p)/real(I)*AggH(:,it+1)+real(count_rent_p)/real(I)*AggR(:,it+1)
             
              AggS(:,it+1)=interp1q(p,AggS(:,it),pstar(it-1))*(1d0-delta)+AggH_s
           
			
			psub=1
			psub3=countP

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 		    if(AggD(psub3) .gt. AggS(psub3,it+1)) then
				pstar(it)=p(psub3)
			elseif(AggD(psub) .lt. AggS(psub,it+1)) then
				pstar(it)=p(1)
			else
				call sub_bisect_excess_demand1(pstar(it),excess_demand1,p,Aggd,AggS(:,it+1),p(psub),p(psub3))
			endif
			 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            call base_fun(p,pstar(it),ilo_pstar(:),theta_pstar(:))
			print*, it, YT_loc(it), 'pstar', pstar(it), interp1q(p,AggS(:,it+1),pstar(it))
			print*, it, 'p', p
			print*, it, 'pprime 1', pprime_new(1,:)
			print*, it, 'AggD', AggD
			print*, it, 'AggS', AggS(:,it+1)
			print*, it, 'AggH_t', AggH(:,it+1)
			print*, it, 'AggR_t', AggR(:,it+1)		
			print*, 'AggD*', interp1q(p,Aggd,pstar(it))
			print*, 'AggH*', interp1q(p, AggH(:,it+1),pstar(it))
			print*, 'AggR*', interp1q(p, AggR(:,it+1),pstar(it))			
			print*, 'Fraction Homeowner', real(count_buy_p)/real(I)
			print*, 'Fraction renter', real(count_rent_p)/real(I)
			print*, 'Percent buy/sell', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,1))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,1)))*(100d0/real(I))
			print*, 'Percent stay', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,2))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,2)))*(100d0/real(I))
			print*, 'Percent refi', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,3))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,3)))*(100d0/real(I))
			print*, 'Percent rent', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,4))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,4)))*(100d0/real(I))
			print*, 'Percent default', (theta_pstar(1)*real(max_loc_p(ilo_pstar(1),it,5))+theta_pstar(2)*real(max_loc_p(ilo_pstar(2),it,5)))*(100d0/real(I))
					

      
           

            count_AggR(:,it+1)=0
            count_AggH(:,it+1)=0
            count_AggS(:,it+1)=0
            Yreveal=AggT_loc(it+1) ! This sets the counts and house prices correctly for next period.

            preg(AggT_loc(it),Yreveal,count(AggT_loc(it),Yreveal))=pstar(it)
			pregZ(AggT_loc(it),countZ(AggT_loc(it)))=pstar(it)
            print*, 'AggTloc, Yrveal, count', AggT_loc(it),Yreveal, count(AggT_loc(it),Yreveal),preg(AggT_loc(it),Yreveal,count(AggT_loc(it),Yreveal))
            do iYY=1,nY*nlamb
                if (iYY/=Yreveal) then
                    count(AggT_loc(it),iYY)=count(AggT_loc(it),iYY)-1
                        
                end if
            end do
           
           
            if ( it<T_sim-1) then
            	print*, 'prices',it,Aggt_loc(it), preg(Aggt_loc(it),Yreveal,count(Aggt_loc(it),Yreveal)) !, pstar(it)

            	! Set R_new as R so it can update recursively 
				Rt(aggT_loc(it),Yreveal,count(aggt_loc(it),Yreveal),:,:)=Rt_new(:,:,Yreveal)
               
					
				! Only update the a0 coefficient with the realized observation.		
				acoeff_solve(Aggt_loc(it),:,1)=acoeff(Aggt_loc(it),:,1)
				acoeff_solve(Aggt_loc(it),Yreveal,1)=gain(Aggt_loc(it),Yreveal)*det(Yreveal)*(Rt_new(2,2,Yreveal)-Rt_new(1,2,Yreveal)*xvec(Yreveal))*(p1(Yreveal)-pred(Yreveal))+acoeff(AggT_loc(it),Yreveal,1) !!+acoeff(Aggt_loc(it-1),Aggt_loc(it),1) a_out_all(last(AggT_loc(it),Yreveal),Yreveal,1)!

				do iYY=1,nY*nlamb
					acoeff_solve(Aggt_loc(it),iYY,2)=gain(Aggt_loc(it),iYY)*det(iYY)*(-Rt_new(2,1,iYY)+Rt_new(1,1,Yreveal)*xvec(iYY))*(p1(Yreveal)-pred(Yreveal))+acoeff(AggT_loc(it),iYY,2) !+a_out(it-1,2)!+acoeff(Aggt_loc(it-1),Aggt_loc(it),2) a_out_all(last(AggT_loc(it),Yreveal),Yreveal,2)!
			
				end do
					
				acoeff_solve2(Aggt_loc(it),:,1)=acoeff2(Aggt_loc(it),:,1)
				acoeff_solve2(Aggt_loc(it),Yreveal,1)=gain2*gain(Aggt_loc(it),Yreveal)*det(Yreveal)*(Rt_new(2,2,Yreveal)-Rt_new(1,2,Yreveal)*xvec(Yreveal))*(p1(Yreveal)-pred2(Yreveal))+acoeff2(AggT_loc(it),Yreveal,1) !!+acoeff(Aggt_loc(it-1),Aggt_loc(it),1) a_out_all(last(AggT_loc(it),Yreveal),Yreveal,1)!

				do iYY=1,nY*nlamb
					acoeff_solve2(Aggt_loc(it),iYY,2)=gain2*gain(Aggt_loc(it),iYY)*det(iYY)*(-Rt_new(2,1,iYY)+Rt_new(1,1,Yreveal)*xvec(iYY))*(p1(Yreveal)-pred2(Yreveal))+acoeff2(AggT_loc(it),iYY,2) !+a_out(it-1,2)!+acoeff(Aggt_loc(it-1),Aggt_loc(it),2) a_out_all(last(AggT_loc(it),Yreveal),Yreveal,2)!
			
				end do
					
					
				do iYY=ind(1),ind(2)
		
                    pprime_solve(1,iYY)=exp(acoeff_solve(Aggt_loc(it),iYY,1)+acoeff_solve(Aggt_loc(it),iYY,2)*log(p(ilo_pstar(1)))) !*log(pstar(it)))
					  
					pprime_solve(2,iYY)=exp(acoeff_solve(Aggt_loc(it),iYY,1)+acoeff_solve(Aggt_loc(it),iYY,2)*log(p(ilo_pstar(2)))) !*log(pstar(it)))
					 
					pprime_solve2(1,iYY)=exp(acoeff_solve2(Aggt_loc(it),iYY,1)+acoeff_solve2(Aggt_loc(it),iYY,2)*log(p(ilo_pstar(1)))) !*log(pstar(it)))
					  
					pprime_solve2(2,iYY)=exp(acoeff_solve2(Aggt_loc(it),iYY,1)+acoeff_solve2(Aggt_loc(it),iYY,2)*log(p(ilo_pstar(2)))) !*log(pstar(it)))
					 
					print*, 'e_out, iYY', iYY, p1(iYY)-pred(iYY)
					print*, 'realized', pstar(it)-exp(acoeff_solve(Aggt_loc(it),iYY,1)+acoeff_solve(Aggt_loc(it),iYY,2)*p1(iYY))
				end do 
					
				print*, 'acoeff solve 1', acoeff_solve(Aggt_loc(it),:,1),'acoeff 1', acoeff(Aggt_loc(it),:,1)
				print*, 'acoeff solve 2', acoeff_solve(Aggt_loc(it),:,2),'acoeff 2', acoeff(Aggt_loc(it),:,2)
               
			   	print*, 'pprime_new', interp1q(p,pprime_new(1,:),pstar(it)),interp1q(p,pprime_new(2,:),pstar(it))
				print*, 'pprime_solve', pprime_solve(:,1), 'prediction 1', exp(acoeff(Aggt_loc(it),ind(1),1)+acoeff(Aggt_loc(it),ind(1),2)*log(pstar(it)))
				print*, 'pprime_solve', pprime_solve(:,2), 'prediction 2', exp(acoeff(Aggt_loc(it),ind(2),1)+acoeff(Aggt_loc(it),ind(2),2)*log(pstar(it)))
			
            end if
			last(AggT_loc(it),Yreveal)=it
			
            a_out(it,:)=acoeff_solve(AggT_loc(it),Yreveal,:)

			
			a_out_all(it,:,1)=acoeff_solve(aggT_loc(it),:,1)
			a_out_all(it,:,2)=acoeff_solve(aggT_loc(it),:,2)
			
			acoeff(AggT_loc(it),Yreveal,1)=acoeff_solve(AggT_loc(it),Yreveal,1)
			
			acoeff(AggT_loc(it),1:4,2)=acoeff_solve(AggT_loc(it),Yreveal,2)
			
			acoeff2(AggT_loc(it),Yreveal,1)=acoeff_solve2(AggT_loc(it),Yreveal,1)
			
			acoeff2(AggT_loc(it),1:4,2)=acoeff_solve2(AggT_loc(it),Yreveal,2)
	
			print*, it, 'count',count(aggt_loc(it),Yreveal)
            print*, it, 'ilo', ilo_pstar(:),'theta', theta_pstar
			print*, it, 'pprime grid 1,1', pprime(iY,1,ilo_pstar(1),:),'pprime_solve', pprime_solve(1,ind(1)), 'bsearch',bsearch(pprime_solve(1,ind(1)),pprime(iY,1,ilo_pstar(1),:))
			print*, it, 'pprime grid 1,2', pprime(iY,2,ilo_pstar(1),:),'pprime_solve', pprime_solve(1,ind(2)), 'bsearch',bsearch(pprime_solve(1,ind(2)),pprime(iY,2,ilo_pstar(1),:))
			print*, it, 'pprime grid 2,1', pprime(iY,1,ilo_pstar(2),:),'pprime_solve', pprime_solve(2,ind(1)), 'bsearch',bsearch(pprime_solve(2,ind(1)),pprime(iY,2,ilo_pstar(2),:))
			print*, it, 'pprime grid 2,2', pprime(iY,2,ilo_pstar(2),:),'pprime_solve', pprime_solve(2,ind(2)), 'bsearch',bsearch(pprime_solve(2,ind(2)),pprime(iY,2,ilo_pstar(2),:))
            print*, it,'pstar, preg', pstar(it),preg(Aggt_loc(it),Yreveal,count(aggt_loc(it),Yreveal))
            print*, it, 'acoeff', acoeff(aggt_loc(it),Yreveal,:)
			print*, it, 'acoeff2', acoeff(aggt_loc(it),Yreveal,:)
			print*, acoeff_solve(AggT_loc(it),:Yreveal,2)
          
            !print*, '__________________________________________'
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			
			rental_sim(it)=theta_pstar(1)*interp2q(pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),rental(AggT_loc(it),ilo_pstar(1),:,:), pprime_solve(1,ind(1)),pprime_solve(1,ind(2)))+&
			theta_pstar(2)*interp2q(pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),rental(AggT_loc(it),ilo_pstar(2),:,:), pprime_solve(2,ind(1)),pprime_solve(2,ind(2)))

			print*, 'rental sim', rental_sim(it)!,pstar(it)*(1.0d0+tauh),-(1.0d0-delta)*dot_product((/exp(acoeff(Aggt_loc(it),ind(1),1)+acoeff(Aggt_loc(it),ind(1),2)*log(pstar(it))),exp(acoeff(Aggt_loc(it),ind(2),1)+acoeff(Aggt_loc(it),ind(2),2)*log(pstar(it)))/),FYY(iAgg,ind(1):ind(2)))/(1.0d0+r_lend_grid(LambT_loc(it))))

				
			do iI=1,I
				do iD=1,nD+2
                    V_sim_star(iI,it,iD)=theta_pstar(1)*V_sim(iI,ilo_pstar(1),iD)+theta_pstar(2)*V_sim(iI,ilo_pstar(2),iD)
                end do
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
						max_loc_sim_star(iI,it)=maxloc((/V_sim_star(iI,it,2),V_sim_star(iI,it,rent),V_sim_star(iI,it,5)/),dim=1)
						if (max_loc_sim_star(iI,it)==1) then
							max_loc_sim_star(iI,it)=2
						else if (max_loc_sim_star(iI,it)==2) then
							max_loc_sim_star(iI,it)=rent
						else if (max_loc_sim_star(iI,it)==3) then
							max_loc_sim_star(iI,it)=5
						end if  
                    end if

					
					countD(it,max_loc_sim_star(iI,it))=countD(it,max_loc_sim_star(iI,it))+1

					if (rent_buy(iI,it)==2 .and. (max_loc_sim_star(iI,it)==2 .or. max_loc_sim_star(iI,it)==3)) then !Own to refi and pay
                        H_sim_star(iI,it+1)=H_sim_star(iI,it)
                        

                        L_sim_star(iI,it+1)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),LL_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2)))+theta_pstar(2)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),LL_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
						
                    
                        B_sim_star(iI,it+1)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),BB_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2)))+theta_pstar(2)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),BB_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
						
                       
					    C_sim_star(iI,it)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),C_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2)))+theta_pstar(2)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),C_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
						
						
						q_sim_star(iI,it)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),q(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2)))+theta_pstar(2)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),q(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:),&
                        L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
                       
						
						if(age_sim(iI,it+1)>2) then
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
							H_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),HH_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) +  &
							theta_pstar(2)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),HH_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
						
							L_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),LL_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) +&
							theta_pstar(2)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),LL_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
							B_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),BB_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),BB_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
					
							C_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),C_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),C_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
					
							q_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),q(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(Bs,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),q(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
						else 
						
					   
							H_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),HH_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) +  &
							theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),HH_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
							L_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),LL_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),LL_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))	
							
							                    
							B_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),BB_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),BB_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
						
							C_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),C_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),C_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							
				
							
								
							q_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),q(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:),&
							Wealth_star(iI)+inc,pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + &
							theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),q(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:),&
							Wealth_star(iI)+inc,pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
							

						end if
					
						
						
						housing_tax(iI,it)=tauh*pstar(it)*H_sim_star(iI,it)
                        if (max_loc_sim_star(iI,it)==buy .and. age_sim(iI,it)<T .and. age_sim(iI,it+1)>1) then
                            rent_buy(iI,it+1)=2
						else if (max_loc_sim_star(iI,it)==buy .and. (age_sim(iI,it)==T .or. age_sim(iI,it+1)==1)) then 
						    rent_buy(iI,it+1)=1
                        else if (max_loc_sim_star(iI,it)==rent) then
                            rent_buy(iI,it+1)=1

                        end if
                            
                          


                    else if (rent_buy(iI,it)==2 .and. max_loc_sim_star(iI,it)==5) then
					
					
						q_sim_star(iI,it)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),q(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + theta_pstar(2)*interp5qnoextrap(L,B,H,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),q(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
					
						H_sim_star(iI,it+1)=R(1)
						rent_buy(iI,it+1)=1

						L_sim_star(iI,it+1)=0d0 
                        
   
                        B_sim_star(iI,it+1)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),BB_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + theta_pstar(2)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),BB_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
						
					
						
						C_sim_star(iI,it)=theta_pstar(1)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),C_b(:,:,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(1,ind(1)),pprime_solve2(1,ind(2))) + theta_pstar(2)*interp5qnoextrap(L,H,B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),C_b(:,:,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							L_sim_star(iI,it),H_sim_star(iI,it),B_sim_star(iI,it),pprime_solve2(2,ind(1)),pprime_solve2(2,ind(2)))
				
				
					
				else if (rent_buy(iI,it)==1 .and. (max_loc_sim_star(iI,it)==1 .or. max_loc_sim_star(iI,it)==4)) then     !Own to rent or default
                       
					   H_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),HH_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(1,ind(1)),pprime_solve(1,ind(2))) + theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),HH_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(2,ind(1)),pprime_solve(2,ind(2)))
							
					   L_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),LL_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(1,ind(1)),pprime_solve(1,ind(2))) + theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),LL_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(2,ind(1)),pprime_solve(2,ind(2)))
							
						
						B_sim_star(iI,it+1)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),BB_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(1,ind(1)),pprime_solve(1,ind(2))) + theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),BB_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(2,ind(1)),pprime_solve(2,ind(2)))
							


						C_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),C_b(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(1,ind(1)),pprime_solve(1,ind(2))) + theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),C_b(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:,max_loc_sim_star(iI,it)),&
							B_sim_star(iI,it),pprime_solve(2,ind(1)),pprime_solve(2,ind(2)))


						q_sim_star(iI,it)=theta_pstar(1)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(1),:),pprime(iY,2,ilo_pstar(1),:),q(1,1,:,ExoT_loc(iI),ilo_pstar(1),age_sim(iI,it),:,:),&
							B_sim_star(iI,it),pprime_solve(1,ind(1)),pprime_solve(1,ind(2))) + theta_pstar(2)*interp3qnoextrap(B,pprime(iY,1,ilo_pstar(2),:),pprime(iY,2,ilo_pstar(2),:),q(1,1,:,ExoT_loc(iI),ilo_pstar(2),age_sim(iI,it),:,:),&
							B_sim_star(iI,it),pprime_solve(2,ind(1)),pprime_solve(2,ind(2)))
                   

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
						tax_star(iI,it)=FnTaxm(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))),(1d0+r_borrow_grid(lambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))/income_star(iI,it)
					else
						tax_star(iI,it)=FnTax(income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it))))/income_star(iI,it)
					end if 
                    employment_star(iI,it)=income(chi(age_sim(iI,it)),loge(age_sim(iI,it),ET_loc(iI,it)),Y(YT_loc(it)))/Y(YT_loc(it))
                    count_work(it)=count_work(it)+1     
                else 
                    ret(iI,it)=income_ret(chi(Jret-1),ET_loc(iI,it),mean(Y))
                    income_star(iI,it)=ret(iI,it)
					if (max_loc_sim_star(iI,it)==2 .or. max_loc_sim_star(iI,it)==3) then
						tax_star(iI,it)=FnTaxm(ret(iI,it),(1d0+r_borrow_grid(lambT_loc(it)))*L_sim_star(iI,it)*H_sim_star(iI,it),minval((/0d0,B_sim_star(iI,it)/)),r_borrow_grid(LambT_loc(it)))/ret(iI,it)
					else 
						tax_star(iI,it)=FnTax(ret(iI,it))/ret(iI,it)
					end if 
                end if  

				max_loc_b(it,max_loc_sim_star(iI,it))=max_loc_b(it,max_loc_sim_star(iI,it))+1

				max_count_age(it,age_sim(iI,it),max_loc_sim_star(iI,it))=max_count_age(it,age_sim(iI,it),max_loc_sim_star(iI,it))+1
					
					homeownership(it)=homeownership(it)+real(rent_buy(iI,it)-1d0)
 


		   end do


			do iI=1,I
						if (age_sim(iI,it)==T) then
					h_sim_star(iI,it+1)=0d0
					l_sim_star(iI,it+1)=0d0
					count_age=count_age+1
					countE(ET_loc(iI,it+1))=countE(ET_loc(iI,it+1))+1
					if (countE(ET_loc(iI,it+1))<=Zerodist(ET_loc(iI,it+1)) ) then
				
						b_sim_star(iI,it+1)=0d0
			
					else if (countE(ET_loc(iI,it+1))>Zerodist(ET_loc(iI,it+1)) .and. countE(ET_loc(iI,it+1))<=nonZerodist(ET_loc(iI,it+1))+Zerodist(ET_loc(iI,it+1)) ) then
						b_sim_star(iI,it+1)=initnw(ET_loc(iI,it+1))
					end if 
										if (b_sim_star(iI,it+1)>maxval(initnw)) then
						b_sim_star(iI,it+1)=0d0
					end if 
					
				end if

			
			
			end do

             employment_b(it)=sum(employment_star(:,it))/real(count_work(it))

             
           ! Dividing income variables irrelevant for simulation to scale by average income
            income_star(:,it)=income_star(:,it) !/employment_b(it)

            tax_star(:,it)=tax_star(:,it) !/employment_b(it)   
            
            ret(:,it)=ret(:,it) !/employment_b(it)
            
            C_sim_star(:,it)=C_sim_star(:,it)! /employment_b(it)
              
           ! Taking averages of income variables. Will wait to do so for employment
			      income_b(it)=sum(income_star(:,it))/real(I)
                    
            tax_b(it)=sum(tax_star(:,it))/real(I)

            ret_b(it)=sum(ret(:,it))/real(I-count_work(it))
            Cstar_b(it)=sum(C_sim_star(:,it))/real(I)
            
        
           ! I need the simulation quantities for loans, houses, and consumption to be
           ! in the units of the model next period
           
           ! constraint and default don't need to be normalize
           ! constraint_b(it)=0d0 !sum(constraint_star(:,it))   
            max_loc_sim_star(:,it)=max_loc_sim_star(:,it)+1
           ! Here I'm going to re-scale last period's quantities by employment
           ! So that it's in the correct units for calibration
            if (it>1) then 
 
                L_sim_star(:,it)=L_sim_star(:,it) !/employment_b(it-1)
                H_sim_star(:,it)=H_sim_star(:,it) !/employment_b(it-1)
                B_sim_star(:,it)=B_sim_star(:,it)
                
                employment_star(:,it-1)=employment_star(:,it-1) !/employment_b(it-1)
                
                
                
                Lstar_b(it)=sum(L_sim_star(:,it))/real(sum(rent_buy(:,it)-1))
                Hstar_b(it)=sum(H_sim_star(:,it))/real(I)
                Bstar_b(it)=sum(B_sim_star(:,it))/real(I)
              	Hstar_s(it)=(sum(H_sim_star(:,it))-sum(H_sim_star(:,it-1))*(1-delta))/real(I)
                housing_tax_b(it)=sum(housing_tax(:,it))/real(sum(rent_buy(:,it)-1))
				homeownership(it)=homeownership(it)/real(I) !sum(real(rent_buy(:,it))-1d0)/real(I)
				
							select case(display)	
				case(1)
			print*,'Loans, Houses, Savings', Lstar_b(it),Hstar_b(it),Bstar_b(it)
				case(0)
				continue
				end select
				
               ! print*, it !'AggB',Bstar_b(it)
            else
                print*, 'employment', employment_b(it) , count_work(it)         
            end if
            
            if (it==T_sim/4 .or. it==T_sim/2 .or. it==T_sim*3/4) then
                print*, 'it, employment', it, employment_b(it) , count_work(it),housing_tax_b(it)
            end if
       end do   
    deallocate(V_sim_star)
	deallocate(q_sim_star)
	deallocate(ret)
	deallocate(housing_tax)
	deallocate(V_sim)
	deallocate(H_sim)
	deallocate(max_loc_sim)
	deallocate(Wealth)
	deallocate(Wealth_star)
	
	
       
    end subroutine sub_simulate_learning
end module mod_simulate_learning
