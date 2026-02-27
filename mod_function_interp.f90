module mod_function_interp
	use omp_lib
    use parameters
	use mod_functions
    use mod_matlab
    use mod_interp
	implicit none
contains


	subroutine sub_function_interp(interp,raw,B,L,Lsim,Bsim,Bs,Y,chi,loge,P,H,nA,nA2)
	integer :: iA2, iA1, iJ, iP,iExo,iBB,iH,iLL,iE,iAgg,iPhi,ilamb,iY,iD,iDD,nA,nA2
	real(8) :: inc, liq
	real(8), dimension(:), intent(in) :: B, L, Lsim, Bsim,Bs,P,H,Y,chi
	real(8), dimension(:,:), intent(in) :: loge
	real(8), dimension(:,:,:,:,:,:,:,:,:),intent(out):: interp
	real(8), dimension(:,:,:,:,:,:,:,:,:),intent(in):: raw

	
	!if (nA2==1) then
	!OMP Parallel Default(Shared)
	!$OMP parallel DO collapse(6) default(shared)	 private(iA1,iA2,iJ,iP,iExo,iBB,iH,iLL,inc,liq,iE,iAgg,iPhi,ilamb,iY,iD,iDD) num_threads(5) if (NA>1)
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
							
								
							interp(1,1,iBB,iExo,iP,iJ,iA1,iA2,1)=interp1qnoextrap(B,raw(1,1,:,iExo,iP,iJ,iA1,iA2,1),Bsim(iBB)) !sum(theta_b(:,iBB)*V_b(1,1,ilo_b(:,iBB),iExo,iP,iJ,iA1,iA2,iD))
							interp(1,1,iBB,iExo,iP,iJ,iA1,iA2,4)=interp1qnoextrap(B,raw(1,1,:,iExo,iP,iJ,iA1,iA2,4),Bsim(iBB)) !sum(theta_b(:,iBB)*V_b(1,1,ilo_b(:,iBB),iExo,iP,iJ,iA1,iA2,iD))

						
							do iH=1,nH
								do iLL=1,nLsim
									interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,2)=interp2qnoextrap(L,B,raw(:,iH,:,iExo,iP,iJ,iA1,iA2,2),Lsim(iLL),Bsim(iBB))
									interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,3)=interp2qnoextrap(L,B,raw(:,iH,:,iExo,iP,iJ,iA1,iA2,3),Lsim(iLL),Bsim(iBB))
									interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,5)=interp2qnoextrap(L,B,raw(:,iH,:,iExo,iP,iJ,iA1,iA2,5),Lsim(iLL),Bsim(iBB))
										
									if (iJ<Jret) then ! This is the next period value we are solving for
										inc=FnTax(income(chi(iJ),loge(iJ,iE),Y(iY)))-FnTaxm(income(chi(iJ),loge(iJ,iE),Y(iY)),Lsim(iLL)*H(iH),minval((/0d0,Bsim(iBB)/)),r_borrow_grid(iLamb))
									else
										inc=FnTax(income_ret(chi(Jret-1),iE,mean(Y)))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),Lsim(iLL)*H(iH),minval((/0d0,Bsim(iBB)/)),r_borrow_grid(iLamb))
									end if 
										
									liq=Bsim(iBB)+(1d0-delta-tauh-cost)*P(iP)*H(iH)-(1d0+r_borrow_grid(ilamb))*Lsim(iLL)*H(iH)+inc 
										
								
										if  (liq<Bsim(1)) then
	
							 				!	V_sim_interp(1,1,iBB,iExo,iP,iJ,iA1,iA2,iD)=V_b(1,1,iBB,iExo,iP,iJ,iA1,iA2,iD)
											!	V_sim_interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,iD)=interp1q(BS,V_b(1,1,:,iExo,iP,iJ,iA1,iA2,iD),liq)

										
											interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,6)=interp1qnoextrap(Bs,raw(1,1,:,iExo,iP,iJ,iA1,iA2,6),liq)
											interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,7)=interp1qnoextrap(Bs,raw(1,1,:,iExo,iP,iJ,iA1,iA2,7),liq)
	
										else		
											!V_sim_interp(1,1,iBB,iExo,iP,iJ,iA1,iA2,iD)=V_b(1,1,iBB,iExo,iP,iJ,iA1,iA2,iDD)
											
											interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,6)=interp1qnoextrap(B,raw(1,1,:,iExo,iP,iJ,iA1,iA2,1),liq)
											interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,7)=interp1qnoextrap(B,raw(1,1,:,iExo,iP,iJ,iA1,iA2,4),liq)
												
											!	if (iJ==29 .and. abs(interp1qnoextrap(Bsim,V_sim_interp(1,1,:,iExo,iP,iJ,iA1,iA2,1),liq)-interp1qnoextrap(B,V_b(1,1,:,iExo,iP,iJ,iA1,iA2,iDD),liq))>0.0001 .and. iD==6 ) then
											!	print*, iLL,iH,iBB,iExo,iP,iJ,iDD,iD, interp1qnoextrap(Bsim,V_sim_interp(1,1,:,iExo,iP,iJ,iA1,iA2,1),liq),interp1qnoextrap(B,W_own(1,1,:,iExo,iP,iJ,iA1,iA2,iDD),liq),interp1qnoextrap(B,V_b(1,1,:,iExo,iP,iJ,iA1,iA2,1),liq),interp1qnoextrap(Bsim,H_interp(1,1,:,iExo,iP,iJ,iA1,iA2,1),liq),interp1qnoextrap(B,HH_b(1,1,:,iExo,iP,iJ,iA1,iA2,iDD),liq),interp1qnoextrap(B,HH_b(1,1,:,iExo,iP,iJ,iA1,iA2,1),liq)
											!	end if 

										end if 
										!if (iExo==1)
										!print*, iA2, iA1, iJ,iP,iExo,iBB,interp(iLL,iHH,iBB,iExo,iP,iJ,iA1,iA2,:), raw(1,iH,1,iExo,iP,iJ,iA1,iA2,:)
										
										!end if 
										!													if (iD==7 .and. iP==1 .and. iExo==1) then
										!		print*, ill,iH,iBB,Lsim(iLL),liq,V_sim_interp(iLL,iH,iBB,iExo,iP,iJ,iA1,iA2,iD),V_b(1,1,iBB,iExo,iP,iJ,iA1,iA2,iDD)
										!		end if
						
								end do
							end do
						end do
					end do
				end do
			end do
		end do
	end do
	!OMP END DO
	!$OMP END PARALLEL Do
	!end if
	end subroutine sub_function_interp
end module mod_function_interp