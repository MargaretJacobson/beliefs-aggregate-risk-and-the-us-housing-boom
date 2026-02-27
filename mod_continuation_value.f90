module mod_continuation_value
    !use omp_lib
    use parameters
	use mod_matlab
	use mod_globals
    use mod_interp
	use mod_functions
    implicit none
contains
     subroutine sub_continuation_value(L,H,p,pprime,ilo_pprime,theta_pprime,FYY,Fee,V_b,W_b,q,q_e,LL_b,iJ,iA1,chi,loge,Y,W_loc_own,W_loc_rent,W_own)
     implicit none
     integer :: iLL,iHH,iEExo,iE,iEE,iH,iD,iB,iBB,iL,iP,iExo,ilamb,iY,ilamblamb, ibeliefs,ilo(2),ilo2,bump2,iAgg,iYY,iPhi,iPhiPhi
	 integer, intent(in) :: iA1
     real(8), dimension (:), intent(in) :: L,H,p,chi,Y
     real(8), dimension (:,:), intent (in) :: FYY, Fee,loge
	 real(8), dimension (:,:,:,:), intent (in) :: pprime
     real(8), dimension (:,:,:,:,:,:), intent(out) :: W_b
     real(8), dimension (:,:,:,:,:,:), intent(in) :: V_b
	 real(8), dimension (:,:,:,:,:), intent(in) :: LL_b
     real(8) :: pm, theta(2),theta2(2),mp,inc,pay_outp,bump, pay_out1,pay_out2,wealth
	 integer, intent(in) :: iJ
	 integer, dimension(:,:,:,:,:), intent(in) :: ilo_pprime
	 real(8), dimension(:,:,:,:,:), intent(in) :: theta_pprime
	 real(8), allocatable :: W_rent(:,:,:)
	 real(8), dimension(:,:,:,:,:,:), intent(out) :: W_own
	 real(8), allocatable :: W_own2(:,:,:,:,:)
	 real(8), intent(inout), dimension(:,:,:,:,:,:) :: q,q_e
	 real(8), allocatable :: pay_out(:,:,:,:,:)
	 integer, dimension(:,:,:,:,:), intent(out) :: W_loc_own
	 integer, dimension(:,:,:), intent(out) :: W_loc_rent
	 real(8), dimension(nY*nlamb*nphi,nP) :: pprime_print
	
	
	allocate(W_own2(nL,nH,nB,nY*nlamb*nE*nPhi,nP))
	allocate(pay_out(nL,nH,nB,nY*nlamb*nE*nPhi,nP))
	allocate(W_rent(nB,nY*nlamb*nE*nPhi,nP)) 	
	 
    W_b=0d0 !Needs to be zero because it's summed over
	pay_out=1d0
	q(:,:,:,:,:,1)=0d0
	q_e(:,:,:,:,:,1)=0d0
	mp=0d0
	pm=0d0
	! print*, 'HCV', H
	

	
	W_own=V_b(:,:,:,:,:,1:5) !Initializing here and will later replace some values
	W_own2=death
	W_loc_own=0
	W_loc_rent=0
	W_rent=0d0 !previously had death
	if (iJ+1==T) then
		W_rent=V_b(1,1,:,:,:,4)
		W_loc_rent=4
	else
		W_rent=maxval(V_b(1,1,:,:,:,(/1,4/)),dim=4)
		W_loc_rent=maxloc(V_b(1,1,:,:,:,(/1,4/)),dim=4)
	end if 
		
	!print*, 'shape W', shape(W_own), 'ND', nD
	!print*, iJ
	pprime_print=0d0
	do iAgg=1,ny*nPhi*nlamb
		do iYY=1,nY*nPhi*nlamb
			pprime_print(iAgg,:)=pprime_print(iAgg,:)+pprime(iAgg,iYY,:,1)*FYY(iAgg,iYY)
		end do
	end do
	
    do iP=1,nP 
	!print*, 'iP', iP, 'ilo1', ilo_pprime(:,iP,1),'ilo2',ilo_pprime(:,iP,2),'theta1',theta_pprime(:,iP,1),'theta2',theta_pprime(:,iP,2)
		do iBB=1,nB        
			do iExo=1,nY*nlamb*nE*nPhi
			iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
			iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
			iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
			ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
			iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr
					

				do iEExo=1,nY*nlamb*nE*nPhi
					iEE=iEExo-nE*(ceiling(real(iEExo)/real(nE))-1)
					iYY=(1-perfect_corr)*ceiling(real(iEExo)/real(nE))+perfect_corr
					iPhiPhi=(1-perfect_corr)*ceiling(real(iEExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
					ilamblamb=(1-perfect_corr)*ceiling(real(iEExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhiPhi)-1d0)+perfect_corr
					ibeliefs=(1-perfect_corr)*(1-floor(real(type_sim)/2d0))+(floor(real(type_sim)/2d0))*((1-perfect_corr)*ceiling(real(iEExo)/real(nE))-nY*(ilamblamb-1)-real(nY)*real(nlamb)*(real(ipHiPhi)-1d0)+perfect_corr)
					
				
					W_b(1,1,iBB,iExo,iP,1)=W_b(1,1,iBB,iExo,iP,1)+beta*(theta_pprime(1,iAgg,iYY,iP,ibeliefs)*W_rent(iBB,iEExo,ilo_pprime(1,iAgg,iYY,iP,ibeliefs))+theta_pprime(2,iAgg,iYY,iP,ibeliefs)*W_rent(iBB,iEExo,ilo_pprime(2,iAgg,iYY,iP,ibeliefs)))*FYY(iAgg,iYY)*Fee(iE,iEE)     

			 
			
				end do
!print*, 'renters', iP,iBB,iAgg,iE,W_b(1,1,iBB,iAgg,iE,iP,1),V_b(1,1,iBB,iY,iE,iP,4)			 
	
				do iHH=1,nH
					do iLL=1,nL                  
							!print*, iE
						if (iJ+1<Jret) then ! This is the next period value we are solving for
							inc=FnTax(income(chi(iJ+1),loge(2,iE),Y(iY)))-FnTaxm(income(chi(iJ+1),loge(2,iE),Y(iY)),L(iLL)*H(iHH),minval((/0d0,B(iBB)/)),r_borrow_grid(ilamb))
						else
							inc=FnTax(income_ret(chi(Jret-1),iE,mean(Y)))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),L(iLL)*H(iHH),minval((/0d0,B(iBB)/)),r_borrow_grid(ilamb))
						end if 
							Wealth=B(iBB)+(1d0-delta-tauh-cost)*p(iP)*H(iHH)-(1d0+r_borrow_grid(ilamb))*L(iLL)*H(iHH)+inc
							! replacing W_own with the correct values if an owner sells 
						if (Wealth <=B(1)) then 
							! I have no idea what I am doing with this block here
							!W_own(iLL,iHH,iBB,iExo,iP,4)=V_b(iLL,iHH,iBB,iExo,iP,7)
							!W_own(iLL,iHH,iBB,iExo,iP,1)=V_b(iLL,iHH,iBB,iExo,iP,6)
							!if (iJ+1<T) then
								!W_own(iLL,iHH,iBB,iExo,iP,1)=death
								
							!else 
							
							if (Wealth<=BS(1)) then
								W_own(iLL,iHH,iBB,iExo,iP,4)=death
								W_own(iLL,iHH,iBB,iExo,iP,1)=death
							
							else 
								W_own(iLL,iHH,iBB,iExo,iP,4)=interp1qnoextrap(BS,V_b(1,1,:,iExo,iP,7),Wealth)
								
								if (iJ+1<T) then
									W_own(iLL,iHH,iBB,iExo,iP,1)=interp1qnoextrap(BS,V_b(1,1,:,iExo,iP,6),Wealth)
								end if 
							end if 
							
							bump2=bsearch(wealth,BS )
						else
							W_own(iLL,iHH,iBB,iExo,iP,4)=interp1qnoextrap(B,V_b(1,1,:,iExo,iP,4),Wealth)
							if (iJ+1<T) then
								W_own(iLL,iHH,iBB,iExo,iP,1)=interp1qnoextrap(B,V_b(1,1,:,iExo,iP,1),Wealth)
							end if 
						
						bump2=bsearch(wealth,B )
						end if 
						W_own(iLL,iHH,iBB,iExo,iP,5)=V_b(1,1,iBB,iExo,iP,5)

					end do
				end do
			end do
		end do
	end do
	
	if (iJ+1<T ) then 
			
		W_own2=maxval(W_own(:,:,:,:,:,(/1,2,3,4,5/)),dim=6)
		W_loc_own=maxloc(W_own(:,:,:,:,:,(/1,2,3,4,5/)),dim=6)
      
	end if 

	do iP=1,nP	
		do iExo=1,nY*nlamb*nE*nPhi
			iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
			iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
			iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
			ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
			iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr
			do iHH=1,nH
				do iLL=1,nL
					do iBB=1,nB 
  
                 
						if (iJ+1==T) then 
							
							pm=MinPay(L(iLL)*H(iHH),iJ+1,r_borrow_grid(ilamb))
							if (L(iLL)==0) then 
								pm=0d0
							end if 
							mp=H(iHH)*L(iLL)*(1d0+r_borrow_grid(ilamb))-pm
								if (W_own(iLL,iHH,iBB,iExo,iP,5) .gt. W_own(iLL,iHH,iBB,iExo,iP,1) .and. W_own(iLL,iHH,iBB,iExo,iP,5) .gt. W_own(iLL,iHH,iBB,iExo,iP,2)) then
									W_own2(iLL,iHH,iBB,iExo,iP)=W_own(iLL,iHH,iBB,iExo,iP,5)
									W_loc_own(iLL,iHH,iBB,iExo,iP)=5
									if (L(iLL)>0) then
										pay_out(iLL,ihh,iBB,iExo,iP)=min(1d0,(1d0-delta_default)*P(iP)/(L(iLL)*(1d0+r_borrow_grid(ilamb))))
									end if 

								else if (W_own(iLL,iHH,iBB,iExo,iP,1) .gt. W_own(iLL,iHH,iBB,iExo,iP,2)) then
									W_own2(iLL,iHH,iBB,iExo,iP)=W_own(iLL,iHH,iBB,iExo,iP,1)
									W_loc_own(iLL,iHH,iBB,iExo,iP)=1
								else 
									W_own2(iLL,iHH,iBB,iExo,iP)=W_own(iLL,iHH,iBB,iExo,iP,2)
									W_loc_own(iLL,iHH,iBB,iExo,iP)=2
									if (L(iLL)>0) then
										pay_out(iLL,iHH,iBB,iExo,iP)=min(1d0,(pm+min(mp,p(iP)*H(iHH)*(1d0-delta-tauh-cost))/(1d0+r_borrow_grid(ilamb)))/(H(iHH)*L(iLL)*(1d0+r_borrow_grid(ilamb))))
									end if 		
								end if 


						else 
							mp=LL_b(iLL,iHH,iBB,iExo,iP)*H(iHH)
                 			pm=H(iHH)*L(iLL)*(1d0+r_borrow_grid(ilamb))-mp		
		
							if (L(iLL)>0 .and. W_loc_own(iLL,iHH,iBB,iExo,iP)==2 ) then
                                                                       
								if (iJ+1==Jret-1) then !MJ made this change, otherwise set to ==
										      ! accomodating weird KMV error where xi is not included when iJ==22
									pay_out(iLL,iHH,iBB,iExo,iP)=(pm+interp1qnoextrap(L,q_e(:,iHH,iBB,iExo,iP,2),mp/H(iHH))*mp)/(L(iLL)*H(iHH)*(1d0+r_borrow_grid(ilamb))) 
									
								else
									pay_out(iLL,iHH,iBB,iExo,iP)=(pm+interp1qnoextrap(L,q_e(:,iHH,iBB,iExo,iP,2),mp/H(iHH))*mp)/(L(iLL)*H(iHH)*(1d0+r_borrow_grid(ilamb))) 
								end if 		
						
			
							else if (L(iLL)>0 .and. W_loc_own(iLL,iHH,iBB,iExo,iP)==5 ) then
								pay_out(iLL,iHH,iBB,iExo,iP)=min(1d0,(1d0-delta_default)*P(iP)/(L(iLL)*(1d0+r_borrow_grid(ilamb))))
							end if 


						end if 
					end do
				end do		
			end do
		end do
	end do

	do iP=1,nP	
		do iExo=1,nY*nlamb*nE*nPhi
			iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
			iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
			iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
			ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
			iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr
			do iHH=1,nH
				do iLL=1,nL
					do iBB=1,nB 
							
						
						if (iJ<Jret) then 
							inc=income(chi(iJ),loge(1,iE),Y(iY))
						else
							inc=income_ret(chi(Jret-1),iE,mean(Y))
						end if 
					
											
						do iEExo=1,nY*nlamb*nE*nPhi
							iEE=iEExo-nE*(ceiling(real(iEExo)/real(nE))-1)
							iYY=(1-perfect_corr)*ceiling(real(iEExo)/real(nE))+perfect_corr
							iPhiPhi=(1-perfect_corr)*ceiling(real(iEExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
							ilamblamb=(1-perfect_corr)*ceiling(real(iEExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhiPhi)-1d0)+perfect_corr
							ibeliefs=(1-perfect_corr)*(1-floor(real(type_sim)/2d0))+(floor(real(type_sim)/2d0))*((1-perfect_corr)*ceiling(real(iEExo)/real(nE))-nY*(ilamblamb-1)-real(nY)*real(nlamb)*(real(ipHiPhi)-1d0)+perfect_corr)
					
							W_b(iLL,iHH,iBB,iExo,iP,2)=W_b(iLL,iHH,iBB,iExo,iP,2)+beta*(theta_pprime(1,iAgg,iYY,iP,ibeliefs)*W_own2(iLL,iHH,iBB,iEExo,ilo_pprime(1,iAgg,iYY,iP,ibeliefs))+theta_pprime(2,iAgg,iYY,iP,ibeliefs)*W_own2(iLL,iHH,iBB,iEExo,ilo_pprime(2,iAgg,iYY,iP,ibeliefs)))*FYY(iAgg,iYY)*Fee(iE,iEE)

				!						if (iJ==28 .and. iAgg==2  .and. iP==6) then !.and. .and. iE==5
				!	print*, iBB,iHH,iLL,iExo,iEExo,iAgg,iE,iEE,iYY,theta_pprime(:,iAgg,iYY,iP,ibeliefs),ilo_pprime(:,iAgg,iYY,iP,ibeliefs),W_own2(iLL,iHH,iBB,iEExo,ilo_pprime(2,iAgg,iYY,iP,ibeliefs)),FYY(iAgg,iYY)*Fee(iE,iEE) ,W_loc_own(iLL,iHH,iBB,iEExo,ilo_pprime(2,iAgg,iYY,iP,ibeliefs)), W_own(iLL,iHH,iBB,iEExo,ilo_pprime(2,iAgg,iYY,iP,ibeliefs),:)
									!pay_out(iLL,iHH,iBB,iYY,iEE,iP),q(iLL,iHH,iBB,iAgg,iE,iP,1),pay_out1,pay_out2,FYY(iAgg,iYY)*Fee(iE,iEE)*(theta_pprime(1,iP,ibeliefs)*pay_out1+theta_pprime(2,iP,ibeliefs)*pay_out2),q(iLL,iHH,iBB,iAgg,iE,iP,1)+FYY(iAgg,iYY)*Fee(iE,iEE)*(theta_pprime(1,iP,ibeliefs)*pay_out1+theta_pprime(2,iP,ibeliefs)*pay_out2)
			
		!	end if 

						end do 

								
						do iEExo=1,nY*nlamb*nE*nPhi
							iEE=iEExo-nE*(ceiling(real(iEExo)/real(nE))-1)
							iYY=(1-perfect_corr)*ceiling(real(iEExo)/real(nE))+perfect_corr
							iPhiPhi=(1-perfect_corr)*ceiling(real(iEExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
							ilamblamb=(1-perfect_corr)*ceiling(real(iEExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhiPhi)-1d0)+perfect_corr
							ibeliefs=(1-perfect_corr)*(1-floor(real(type_sim)/2d0))+(floor(real(type_sim)/2d0))*((1-perfect_corr)*ceiling(real(iEExo)/real(nE))-nY*(ilamblamb-1)-real(nY)*real(nlamb)*(real(ipHiPhi)-1d0)+perfect_corr)
					
			
							pay_out1=pay_out(iLL,iHH,iBB,iEExo,ilo_pprime(1,iAgg,iYY,iP,ibeliefs))
							pay_out2=pay_out(iLL,iHH,iBB,iEExo,ilo_pprime(2,iAgg,iYY,iP,ibeliefs))

						    mp=LL_b(iLL,iHH,iBB,iEExo,iP)*H(iHH)
							pm=H(iHH)*L(iLL)*(1d0+r_borrow_grid(ilamblamb))-mp	

							q(iLL,iHH,iBB,iExo,iP,1)=q(iLL,iHH,iBB,iExo,iP,1)+FYY(iAgg,iYY)*Fee(iE,iEE)*(theta_pprime(1,iAgg,iYY,iP,ibeliefs)*pay_out1+theta_pprime(2,iAgg,iYY,iP,ibeliefs)*pay_out2)


						end do

	
						q(iLL,iHH,iBB,iExo,iP,1)=min(1d0,q(iLL,iHH,iBB,iExo,iP,1)) 
						bump=q(iLL,iHH,iBB,iExo,iP,1) !this yields qmret_e from KMV
						q_e(iLL,iHH,iBB,iExo,iP,1)=q(iLL,iHH,iBB,iExo,iP,1)
						if (L(iLL)>lambdaLTV(ilamb)*P(iP) .or. (MinPay(L(iLL)*H(iHH),iJ+1,r_borrow_grid(ilamb))>lambdaPTI(ilamb)*inc .and. iJ<Jret)) then
							q(iLL,iHH,iBB,iExo,ip,1)=0d0
						else 	
							if (iJ==Jret-1) then 
								q(iLL,iHH,iBB,iExo,ip,1)=q(iLL,iHH,iBB,iExo,iP,1)
							else 
								q(iLL,iHH,iBB,iExo,ip,1)=q(iLL,iHH,iBB,iExo,iP,1)-xi(ilamb) ! NOt (1d0-xi(ilamb)
							end if 
						end if 


	
						!print*, iLL, iBB,iY,iE,iHH,iYY,iEE
!print*, 'bounds', BB(iBB)+(1-delta-tauh-kappam)*pprime(iY,iYY)*HH(iHH)-(1+r_borrow)*LL(iLL)
						
!print*, iLL, iHH
!print*, W_b(iLL,iHH,iBB,iY,iE,:)
                           !print*,iHH,iLL,iE,iY, W_b(1,1,:,iY,iE,1)
						
                    end do
                end do
            end do
        end do
    end do

	!deallocate(W_own)
	deallocate(W_own2)
	deallocate(pay_out)
	deallocate(W_rent)

    end subroutine sub_continuation_value
end module mod_continuation_value
