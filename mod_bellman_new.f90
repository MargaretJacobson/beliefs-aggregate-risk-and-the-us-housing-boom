module mod_bellman_new
    use omp_lib
    use parameters
	use mod_matlab
	use mod_globals
	use mod_functions
    use mod_interp
    implicit none
contains
    subroutine sub_bellman_bor(V_b,c_b,LL_b,LLi_b,BB_b,HHi_b,W_b,Y,chi,loge,H,L,p,e,iJin,q,Trans,rental,R)
	implicit none  
    integer, intent(in) :: iJin
    real(8), dimension(:,:,:,:,:,:), intent (out) :: c_b
    real(8), dimension(:,:,:,:,:,:), intent (inout) :: V_b,W_b
    integer, dimension(:,:,:,:,:,:), intent (out) :: LLi_b,HHi_b
	real(8), dimension(:,:,:,:,:,:), intent (out) :: BB_b,LL_b
    real(8), dimension(:,:,:,:,:), intent (in) :: q
	real(8), dimension(:,:) :: bheloc(nH,nP)
    real(8), dimension(:), intent (in) :: H,L,p,Y,chi,loge,R
    real(8), intent(in) :: e
	real(8), dimension(:,:), intent(in) :: rental,Trans
    integer :: iY,iL,iLL,iH,iHH, iE,iD1,iD2,iD3,iD4,iD5, iD6,iD7,iB,iBB,Bbar,iExo,ilamb,iP,iR,iJ,iAgg,iPhi
    real(8) :: Vtemp1a,Vtemp1b,Vtemp2a,Vtemp2b,Vtemp3a,Vtemp3b,Vtemp4a,Vtemp4b,Vtemp5a,Vtemp5b,Vtemp6a,Vtemp6b,Vtemp7a,Vtemp7b,cval1_final,cval6_final, cval3_final,pm,mp,lyp,out!,phi
	real(8) :: BBmin7, BBmax7, BBmid7,btemp7
	real(8) :: BBmin6, BBmax6, BBmid6,btemp6
	real(8) :: BBmin5, BBmax5, BBmid5,btemp5
	real(8) :: BBmin4, BBmax4, BBmid4,btemp4	
	real(8) :: BBmin3, BBmax3, BBmid3,btemp3
	real(8) :: BBmin2, BBmax2, BBmid2,btemp2	
	real(8) :: BBmin1, BBmax1, BBmid1,btemp1		

	DOUBLE PRECISION, EXTERNAL          :: golden
    LLi_b=1
    HHi_b=1
	BB_b=0d0
	LL_b=0d0
    iD3=3
    iD2=2
    iD1=1
	iD4=4
	iD5=5
	iD6=6
	iD7=7
    
	bheloc=0
	
	Vtemp7a=death
	Vtemp6a=death
	Vtemp5a=deathdefault
	Vtemp4a=death
    Vtemp3a=death
	Vtemp2a=death
	Vtemp1a=death
		
	Vtemp1b=death
    Vtemp2b=death
    Vtemp3b=death
	Vtemp4b=death
    Vtemp5b=deathdefault
	Vtemp6b=death
	Vtemp7b=death
		V_b=death
		V_b(:,:,:,:,:,5)=deathdefault
	
	eiJ=e
	iJ=iJin

	
	btemp1=0d0
	btemp2=0d0
	btemp3=0d0
	btemp4=0d0
	btemp5=0d0
	

		
    !$OMP Parallel Default(Shared) copyin(q_lend,q_borrow,phi,eiJ,iJ,qval1,cval1,Wval1,qval3,cval2,cval3,cval4,cval5,qval6,cval6,cval7,Hval1,Hval2,Hval3,Rval4,Rval5,Hval6,Rval7,Wval2,Wval3,Wval4,Wval5,Wval6,Wval7) 
    !$OMP DO collapse(3) private(iE,iP,iY,iExo,ilamb,iAgg,iPhi,iR,iH,iB,iL,iBB,iLL,Vtemp3a,Vtemp3b,Vtemp2a,Vtemp2b,iHH,Vtemp1a,Vtemp1b,BBmin3,BBmax3, BBmid3,btemp3,BBmin2, BBmax2, BBmid2,btemp2,BBmin1, BBmax1, BBmid1,btemp1,cval1_final,cval3_final,Vtemp4b,Vtemp4a,Vtemp5a,Vtemp5b,BBmin4,BBmax4, BBmid4,btemp4,BBmin5,BBmax5, BBmid5,btemp5,pm,mp,Vtemp6a,Vtemp6b,BBmin6, BBmax6, BBmid6,btemp6,cval6_final,Vtemp7b,Vtemp7a,BBmin7,BBmax7, BBmid7,btemp7)
	do iP=1,nP
		do iExo=1,nY*nlamb*nE*nPhi
			do iB=1,nB
				iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
				iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
				iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
				ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
				iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr
				phi=phi_grid(iPhi)
				q_lend=q_lend_grid(ilamb)
				q_borrow=q_borrow_grid(ilamb)
				
				Wval4=W_b(1,1,:,iExo,iP,1)
				do iR=1,nR
					BBmin4=0d0

					if (iJ<Jret) then 
						cval4=B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(iR)
					else
						cval4=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(iR)
						
					end if 

					BBmax4=(cval4-cmin)/q_lend_grid(ilamb)
					BBmid4=(BBmax4+BBmin4)/2.0d0
					
					Rval4=R(iR)
						
					if (iJ<Jret) then
						Vtemp4a=Value4(BBmin4)
						Vtemp4b=Value4(B(2))
						if (bbmax4>bbmin4) then
							if (Vtemp4b<Vtemp4a) then
								Vtemp4a=golden(BBmin4,BBmid4,BBmax4,Value4,1d-5,btemp4)
								cval4=cval4-q_lend_grid(ilamb)*btemp4
							else
								btemp4=0d0
							
							end if 
						else 
							Vtemp4a=-death
							btemp4=0d0
							cval4=cmin
						end if 
						
					else
	
						if (cval4 > cmin+BBmin4) then
							Vtemp4a=Value4(BBmin4)
							Vtemp4b=golden(BBmin4,BBmid4,BBmax4,Value4,1d-5,btemp4)
							!print*, 'Vtemp4a, Vtemp4B, btemp', -Vtemp4a, -Vtemp4b,btemp4
							if(Vtemp4b<Vtemp4a) then
								Vtemp4a=Vtemp4b
								!print*,'btmp', btemp
								!Vtemp4a=Vtemp4b
								cval4=cval4-q_lend_grid(ilamb)*btemp4
							else
								btemp4=0d0
							end if
						else 
							Vtemp4a=-death
							btemp4=0d0
							cval4=cmin
						end if	
					end if 

					if (-Vtemp4a>V_b(1,1,iB,iExo,iP,iD4)) then
							!print*,iE,iY, iB,iR, -Vtemp4a, V_b(iB,iY,iE,iD4), btemp4,cval4
						V_b(1,1,iB,iExo,iP,iD4)=-Vtemp4a
						HHi_b(1,1,iB,iExo,iP,iD4)=iR
						BB_b(1,1,iB,iExo,iP,iD4)=btemp4
						c_b(1,1,iB,iExo,iP,iD4)=cval4
							
					end if
						
										
				end do
				
		! Do not need R loop for defaulters
				BBmin5=0d0
				Wval5=W_b(1,1,:,iExo,iP,1)
				bheloc(:,iP)=-p(iP)*H*0.2d0
				if (iJ<Jret) then 
					cval5=B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(1)
						
				else
					
					cval5=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(1)
				end if 
				BBmax5=(cval5-cmin)/q_lend_grid(ilamb)
				BBmid5=(BBmax5+BBmin5)/2d0 !should be (BBmax5+BBmin5)/2.0d0 
				
				Rval5=R(1)
					!if (cval5 > cmin) then
				if (bbmax5>bbmin5) then
					Vtemp5a=Value5(0d0)
					Vtemp5b=golden(BBmin5,BBmid5,BBmax5,Value5,1d-5,btemp5)
						!print*, 'Vtemp4a, Vtemp4B, btemp', Vtemp4a, Vtemp4b,btemp
					if(Vtemp5b<Vtemp5a) then
						Vtemp5a=Vtemp5b
							!print*,'btmp', btemp
							!Vtemp4a=Vtemp4b
						cval5=cval5-q_lend_grid(ilamb)*btemp5
					else
						btemp5=0d0
					end if
				else 
					Vtemp5a=-deathdefault
					btemp5=0d0
					cval5=cmin
				end if	
				!print*,iP,iE,iAgg, iB, -Vtemp5a, btemp5,bbmax5,bbmin5,bbmid5,B(iB)+income_ret(chi(Jret-1),iE,Y(iY))-FnTax(income_ret(chi(Jret-1),iE,Y(iY)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(1),-q_lend*btemp5,interp1q(B,Wval5,btemp5)/beta
						
					!print*,iE,iY, iB,iR, -Vtemp5a, V_b(iB,iY,iE,iD5), btemp5,cval5		
				if (-Vtemp5a>V_b(1,1,iB,iExo,Ip,iD5)) then ! .and. upper3(iL)>=LL(iLL)) then
                   ! print*,p,iJ,iE,iY,1,iB,iL,iBB,iLL,iD3,Vtmp3 !,W_b(iLL,iBB,1,iY,iE) !cval3
					V_b(:,:,iB,iExo,ip,iD5)=-Vtemp5a
					HHi_b(:,:,iB,iExo,Ip,iD5)=1
					BB_b(:,:,iB,iExo,Ip,iD5)=btemp5
					c_b(:,:,iB,iExo,iP,iD5)=cval5
				end if  
				
						do iH=1,nH	
					do iLL=1,nL !min(PTI(iHH,iAgg,iE),LTV(ilamb))	
					
					
							! Buy a house
						BBmin1=0d0 !min(-(p(iP)-L(iLL))*H(iH)*0.2d0,0.0d0) !This allows them to save through a heloc which is incorrect and forces them to save while buying
							!do iBB=1,nB
						qval1=q(iLL,iH,:,iExo,iP)*L(iLL)*H(iH)
						if (iJ<Jret) then 
							cval1=B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) !fixed cost only shows up when L>0 
							BBmax1=(B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)+H(iH)*L(iLL)-P(iP)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))))/q_lend_grid(ilamb)
						else
							cval1=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) 
							BBmax1=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))-P(iP)*H(iH)+L(iLL)*H(iH) -kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) !missing Trans(iP) and qlend
						end if 

							
						BBmid1=(BBmax1+BBmin1)/2.0d0
						Wval1=W_b(iLL,iH,:,iExo,iP,2)
						Hval1=utilitycost*H(iH)
							
						if (bbmax1> bbmin1) then 
							!if (interp1q(B,cval1(:),BBmin1)>cmin) then
								!print*, iB,iH,iLL,H(iH)*L(iLL)
							Vtemp1a=Value1(0d0)

							Vtemp1b=golden(BBmin1,BBmid1,BBmax1,Value1,1d-5,btemp1)

							if (Vtemp1b<Vtemp1a) then
								Vtemp1a=Vtemp1b
								cval1_final=cval1-q_lend_grid(ilamb)*btemp1+interp1q(B,qval1(:),btemp1)
							else
								btemp1=0d0
							end if
						else
							Vtemp1a=-death
							btemp1=0d0
							cval1_final=cmin
						end if
						if (cval1-q_lend_grid(ilamb)*btemp1+interp1q(B,qval1(:),btemp1) < cmin  ) then		
							Vtemp1a=-death
							btemp1=0d0
							cval1_final=cmin
						end if 
										
						if (-Vtemp1a>V_b(1,1,iB,iExo,iP,iD1)) then
                                   ! print*,iE,iY,iB,iBB,iLL,iHH 
							!	if (((iA1==2 .and. iA2==2 .and. type_sim==2) .or. (iA1==1 .and. iA2==1 .and. type_sim==1)) .and. iJ==T-1) then

                            !        print*,iE,iY,iB,ihh,iLL,iD1,Vtemp1a,btemp1,interp1q(B,W_b(iLL,iHH,:,iY,iE), btemp1)
							!end if 

                            V_b(1,1,iB,iExo,iP,iD1)=-Vtemp1a
                            LLi_b(1,1,iB,iExo,iP,iD1)=iLL
                            HHi_b(1,1,iB,iExo,iP,iD1)=iH
                            BB_b(1,1,iB,iExo,iP,iD1)=btemp1
                            c_b(1,1,iB,iExo,iP,iD1)=cval1_final
						end if


					end do
				end do	
				
									! Refinance
				do iH=1,nH
					do iL=1,nL
						do iLL=max(iL,2),nL !min(PTI(iH,iAgg,iE),LTV(ilamb))
							!print*, LTV(ilamb), PTI(iH,iAgg,iE), iH,iE
							BBmin3=bheloc(iH,iP)
								
							if (iJ<Jret) then 
								cval3=-H(iH)*L(iL)*(1d0+r_borrow_grid(ilamb))+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTaxm(income(chi(iJ),loge(iE),Y(iY)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)-kappam(ilamb) !+q(iLL,:,iH,iAgg,iE,iP)*L(iLL)*H(iH)
								BBmax3=(cval3+p(iP)*H(iH))/q_lend_grid(ilamb)
							else
								cval3=-H(iH)*L(iL)*(1d0+r_borrow_grid(ilamb))+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)-kappam(ilamb) !+q(iLL,:,iH,iAgg,iE,iP)*L(iLL)*H(iH)
								BBmax3=(cval3+L(iLL)*H(iH))/q_lend_grid(ilamb)
							end if 
									
							qval3=q(iLL,iH,:,iExo,iP)*L(iLL)*H(iH)
							BBmid3=(BBmax3+BBmin3)/2d0
							Wval3=W_b(iLL,iH,:,iExo,iP,2)
							Hval3=utilitycost*H(iH)
							if (bbmax3>bbmin3) then
								!if (interp1q(B,cval3(:),BBmin3)>cmin) then
								Vtemp3a=Value3(0d0)
								Vtemp3b=golden(BBmin3,BBmid3,BBmax3,Value3,1d-5,btemp3)
!print*,iD3,iE,iY,iH,iB,Vtemp3a,Vtemp3b,btemp3,interp1q(B,W_b(iLL,iH,:,iY,iE),btemp3) !cval3
								if (Vtemp3b<Vtemp3a) then
									Vtemp3a=Vtemp3b
									if (btemp3<0) then 
										cval3_final=cval3+interp1q(B,qval3,btemp3)-q_borrow_grid(ilamb)*btemp3
									else
										cval3_final=cval3+interp1q(B,qval3,btemp3)-q_lend_grid(ilamb)*btemp3
									end if
								else
									btemp3=0d0
								end if
							else 
								Vtemp3a=-death  
								btemp3=0d0
								cval3_final=cmin
							end if

							if (-Vtemp3a>V_b(iL,iH,iB,iExo,iP,iD3)) then ! .and. upper3(iL)>=LL(iLL)) then
                                !print*,iE,iY,iH,iB,iD3,Vtemp3a,btemp3,interp1q(B,W_b(iLL,iH,:,iY,iE),btemp3) !cval3
								V_b(iL,iH,iB,iExo,iP,iD3)=-Vtemp3a
								LLi_b(iL,iH,iB,iExo,iP,iD3)=iLL
								HHi_b(iL,iH,iB,iExo,iP,iD3)=iH
								BB_b(iL,iH,iB,iExo,iP,iD3)=btemp3
								c_b(iL,iH,iB,iExo,iP,iD3)=cval3_final
                                    !print*, Vtmp2,V_b(iL,iH,iY,iE,iD),VtmpH2
							end if 
						end do
													! Solving payment problem.

							! First solve at pm, then solve at iLL=1,nL
						pm=MinPay(L(iL)*H(iH),iJ,r_borrow_grid(ilamb))
						if (iL .eq. 1) pm=0d0
							mp=(H(iH)*L(iL))*(1.0d0+r_borrow_grid(ilamb))-pm
							BBmin2=bheloc(iH,iP)
							if (iJ<Jret) then
								cval2=p(iP)*(-delta-tauh)*H(iH)+B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTaxm(income(chi(iJ),loge(iE),Y(iY)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)-pm
							else
								cval2=p(iP)*(-delta-tauh)*H(iH)+B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)-pm
							end if 
							BBmax2=(cval2)/q_lend_grid(ilamb)
							BBmid2=(BBmax2+BBmin2)/2.d0
							do iBB=1,nB
								Wval2(iBB)=interp1qnoextrap(L,W_b(:,iH,iBB,iExo,iP,2),mp/H(iH))
							end do 
							Hval2=utilitycost*H(iH)
							if(bbmax2>bbmin2) then
							!if (cval2>cmin) then
								Vtemp2a=Value2(0d0)
								Vtemp2b=golden(BBmin2,BBmid2,BBmax2,Value2,1d-5,btemp2)
!print*,iD2,iE,iY,iH,iB,Vtemp2a,Vtemp2b,btemp2,interp1q(B,W_b(iLL,iH,:,iY,iE),btemp2)
								if (Vtemp2b<Vtemp2a) then
									Vtemp2a=Vtemp2b
									if (btemp2<0) then
										cval2=cval2-q_borrow_grid(ilamb)*btemp2										
									else
										cval2=cval2-q_lend_grid(ilamb)*btemp2
									end if 
								else
									btemp2=0d0
								end if
							else
								Vtemp2a=-death
								btemp2=0d0
								cval2=cmin
							end if
						
							if (-Vtemp2a>V_b(iL,iH,iB,iExo,iP,iD2)) then
                                  ! print*,iE,iY,iH,iB,iD2,-Vtemp2a,btemp2,interp1q(B,W_b(iLL,iH,:,iY,iE),btemp2)
								V_b(iL,iH,iB,iExo,iP,iD2)=-Vtemp2a
								LL_b(iL,iH,iB,iExo,iP,iD2)=mp/H(iH)
								HHi_b(iL,iH,iB,iExo,iP,iD2)=iH
								BB_b(iL,iH,iB,iExo,iP,iD2)=btemp2
								c_b(iL,iH,iB,iExo,iP,iD2)=cval2

                                        !print*, Vtmp2,V_b(iL,iH,iY,iE,iD2),VtmpH2
							end if	
				
							iLL=1	
							do while (L(iL)*H(iH)*(1d0+r_borrow_grid(ilamb))-pm > (L(iLL)*H(iH)) .and. iLL<iL)
								
								if (iJ<Jret) then
									
									cval2=-H(iH)*L(iL)*(1d0+r_borrow_grid(ilamb))+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTaxm(income(chi(iJ),loge(iE),Y(iY)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)+L(iLL)*H(iH)
								else
									cval2=-H(iH)*L(iL)*(1d0+r_borrow_grid(ilamb))+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)+L(iLL)*H(iH)
								end if 
								BBmax2=(cval2)/q_lend_grid(ilamb)
								BBmid2=(BBmax2+BBmin2)/2.d0
								Wval2=W_b(iLL,iH,:,iExo,iP,2)
								Hval2=utilitycost*H(iH)

								if (bbmax2>bbmin2) then
									Vtemp2a=Value2(0d0)
									Vtemp2b=golden(BBmin2,BBmid2,BBmax2,Value2,1d-5,btemp2)
!print*,iD2,iE,iY,iH,iB,Vtemp2a,Vtemp2b,btemp2,interp1q(B,W_b(iLL,iH,:,iAgg,iE),btemp2)
									if (Vtemp2b<Vtemp2a) then
										Vtemp2a=Vtemp2b
										if (btemp2<0) then
											cval2=cval2-q_borrow_grid(ilamb)*btemp2										
										else
											cval2=cval2-q_lend_grid(ilamb)*btemp2
										end if 
									else
										btemp2=0d0
									end if
								else
									Vtemp2a=-death
									btemp2=0d0
									cval2=cmin
								end if

								if (-Vtemp2a>V_b(iL,iH,iB,iExo,iP,iD2)) then
                                  ! print*,iE,iY,iH,iB,iD2,-Vtemp2a,btemp2,interp1q(B,W_b(iLL,iH,:,,iY,iE),btemp2)
									V_b(iL,iH,iB,iExo,iP,iD2)=-Vtemp2a
									LLi_b(iL,iH,iB,iExo,iP,iD2)=iLL
									LL_b(iL,iH,iB,iExo,iP,iD2)=L(iLL)
									HHi_b(iL,iH,iB,iExo,iP,iD2)=iH
									BB_b(iL,iH,iB,iExo,iP,iD2)=btemp2
									c_b(iL,iH,iB,iExo,iP,iD2)=cval2

                                        !print*, Vtmp2,V_b(iL,iH,iY,iE,iD2),VtmpH2
								end if
								iLL=iLL+1
							
                        end do
                    end do
				end do
			end do
			do iB=1,nBS

				iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
				iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
				iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
				ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
				iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr
				phi=phi_grid(iPhi)
				q_lend=q_lend_grid(ilamb)
				q_borrow=q_borrow_grid(ilamb)
				
				Wval7=W_b(1,1,:,iExo,iP,1)
				!print*, iExo,iE,iAgg,iPhi,ilamb,iY
                 !               print*, 'buy: max threads to be used, thread limit, thread number, threads in team=team size', omp_get_max_threads(),omp_get_thread_num(),omp_get_num_threads(),omp_get_team_size(1)
				if(iJ>=Jret .and. Bs(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))>cmin) then
					do iR=1,nR
						BBmin7=0d0

						!if (iJ<Jret) then 
						!else
						cval7=BS(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(iR)
						
						!end if 

						BBmax7=(cval7-cmin)/q_lend_grid(ilamb)
						BBmid7=(BBmax7+BBmin7)/2.0d0
					
						Rval7=R(iR)
						
					!if (iJ<Jret) then
						Vtemp7a=Value7(0d0)
						Vtemp7b=golden(BBmin7,B(2),BBmax7,Value7,1d-5,btemp7)
						!if (bbmax7>bbmin7) then
						
						if (Vtemp7b<Vtemp7a) then
							Vtemp7a=Vtemp7b
							cval7=cval7-q_lend_grid(ilamb)*btemp7
						else
							btemp7=0d0
							cval7=cmin
							
						end if 
					
						

						if (-Vtemp7a>V_b(1,1,iB,iExo,iP,iD7) .and. cval7>=cmin) then
							!print*,iE,iY, iB,iR, -Vtemp7a, V_b(iB,iY,iE,iD7), btemp7,cval7
							V_b(1,1,iB,iExo,iP,iD7)=-Vtemp7a
							HHi_b(1,1,iB,iExo,iP,iD7)=iR
							BB_b(1,1,iB,iExo,iP,iD7)=btemp7
							c_b(1,1,iB,iExo,iP,iD7)=cval7
						else if (	-Vtemp7a>V_b(1,1,iB,iExo,iP,iD7) .and. cval7<cmin) then
							V_b(1,1,iB,iExo,iP,iD7)=death
							HHi_b(1,1,iB,iExo,iP,iD7)=iR
							BB_b(1,1,iB,iExo,iP,iD7)=0d0
							c_b(1,1,iB,iExo,iP,iD7)=cmin
						end if
						
										
					end do
				else if (iJ>=Jret .and. Bs(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))<=cmin) then
					V_b(1,1,iB,iExo,iP,iD7)=death
					HHi_b(1,1,iB,iExo,iP,iD7)=1
					BB_b(1,1,iB,iExo,iP,iD7)=0d0
					c_b(1,1,iB,iExo,iP,iD7)=cmin
				else if (iJ<Jret .and.  Bs(iB)+ income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))-R(1)*rental(iAgg,iP)>cmin) then 
					do iR=1,nR
						BBmin7=0d0
						cval7=BS(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(iR)
						BBmax7=(cval7-cmin)/q_lend_grid(ilamb)
						
					
						Rval7=R(iR)
						
					!if (iJ<Jret) then
						Vtemp7a=Value7(0d0)
						Vtemp7a=Value7(B(2))
						
						!if (bbmax7>bbmin7) then
						
						if (Vtemp7b<Vtemp7a) then
							Vtemp7b=golden(0d0,B(2),BBmax7,Value7,1d-5,btemp7)						
							cval7=cval7-q_lend_grid(ilamb)*btemp7
						else
							btemp7=0d0
							
						end if 
						if (-Vtemp7a>V_b(1,1,iB,iExo,iP,iD7) .and. cval7>=cmin) then
							!print*,iE,iY, iB,iR, -Vtemp7a, V_b(iB,iY,iE,iD7), btemp7,cval7
							V_b(1,1,iB,iExo,iP,iD7)=-Vtemp7a
							HHi_b(1,1,iB,iExo,iP,iD7)=iR
							BB_b(1,1,iB,iExo,iP,iD7)=btemp7
							c_b(1,1,iB,iExo,iP,iD7)=cval7
						else if (	-Vtemp7a>V_b(1,1,iB,iExo,iP,iD7) .and. cval7<cmin) then
							V_b(1,1,iB,iExo,iP,iD7)=death
							HHi_b(1,1,iB,iExo,iP,iD7)=iR
							BB_b(1,1,iB,iExo,iP,iD7)=0d0
							c_b(1,1,iB,iExo,iP,iD7)=cmin
						end if
						
										
					end do
					
						
				else
					V_b(1,1,iB,iExo,iP,iD7)=death
					HHi_b(1,1,iB,iExo,iP,iD7)=1
					BB_b(1,1,iB,iExo,iP,iD7)=0d0
					c_b(1,1,iB,iExo,iP,iD7)=cmin
				end if 
								do iH=1,nH	
					do iLL=1,nL 
						if (iJ>=Jret .and. BS(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) > cmin) then 
					
								! Buy a house
							BBmin6=0d0 
							qval6=q(iLL,iH,:,iExo,iP)*L(iLL)*H(iH)
				
						
							cval6=BS(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) 
							BBmax6=BS(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))-P(iP)*H(iH)+L(iLL)*H(iH) -kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL)))+Trans(iP,ilamb) !missing Trans(iP) and qlend
				

							
							BBmid6=(BBmax6+BBmin6)/2.0d0
							Wval6=W_b(iLL,iH,:,iExo,iP,2)
							Hval6=utilitycost*H(iH)
							

							Vtemp6a=Value6(0d0)

							Vtemp6b=golden(BBmin6,BBmid6,BBmax6,Value6,1d-5,btemp6)
!					
							if (Vtemp6b<Vtemp6a) then
								Vtemp6a=Vtemp6b
								cval6_final=cval6-q_lend_grid(ilamb)*btemp6+interp1q(BS,qval6(:),btemp6)
							else
								btemp6=0d0
							end if
						
							if (cval6-q_lend_grid(ilamb)*btemp6+interp1q(BS,qval6(:),btemp6) < cmin  ) then		
								Vtemp6a=-death
								btemp6=0d0
								cval6_final=cmin
							end if 
				
							if (-Vtemp6a>V_b(1,1,iB,iExo,iP,iD6)) then
	
								V_b(1,1,iB,iExo,iP,iD6)=-Vtemp6a
								LLi_b(1,1,iB,iExo,iP,iD6)=iLL
								HHi_b(1,1,iB,iExo,iP,iD6)=iH
								BB_b(1,1,iB,iExo,iP,iD6)=btemp6
								c_b(1,1,iB,iExo,iP,iD6)=cval6_final
							end if
						else if (iJ>=Jret .and. BS(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) <= cmin) then 
							V_b(1,1,iB,iExo,iP,iD6)=death
							LLi_b(1,1,iB,iExo,iP,iD6)=1
							HHi_b(1,1,iB,iExo,iP,iD6)=1
							BB_b(1,1,iB,iExo,iP,iD6)=0d0
							c_b(1,1,iB,iExo,iP,iD6)=cmin
						else if (iJ<Jret .and. BS(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) > cmin) then 
							BBmin6=0d0 
							qval6=q(iLL,iH,:,iExo,iP)*L(iLL)*H(iH)
							cval6=BS(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) !fixed cost only shows up when L>0 
							BBmax6=(BS(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)+H(iH)*L(iLL)-P(iP)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))))
						
						
							BBmid6=(BBmax6+BBmin6)/2.0d0
							Wval6=W_b(iLL,iH,:,iExo,iP,2)
							Hval6=utilitycost*H(iH)
							

							Vtemp6a=Value6(0d0)

							Vtemp6b=golden(BBmin6,BBmid6,BBmax6,Value6,1d-5,btemp6)
!					
							if (Vtemp6b<Vtemp6a) then
								Vtemp6a=Vtemp6b
								cval6_final=cval6-q_lend_grid(ilamb)*btemp6+interp1q(BS,qval6(:),btemp6)
							else
								btemp6=0d0
							end if
						
							if (cval6-q_lend_grid(ilamb)*btemp6+interp1q(BS,qval6(:),btemp6) < cmin  ) then		
								Vtemp6a=death
								btemp6=0d0
								cval6_final=cmin
							end if 
				
							if (-Vtemp6a>V_b(1,1,iB,iExo,iP,iD6)) then
	
								V_b(1,1,iB,iExo,iP,iD6)=-Vtemp6a
								LLi_b(1,1,iB,iExo,iP,iD6)=iLL
								HHi_b(1,1,iB,iExo,iP,iD6)=iH
								BB_b(1,1,iB,iExo,iP,iD6)=btemp6
								c_b(1,1,iB,iExo,iP,iD6)=cval6_final
							end if
						
						
						else if (iJ<Jret .and. BS(iB)+income(chi(iJ),loge(iE),Y(iY))-FnTax(income(chi(iJ),loge(iE),Y(iY)))+Trans(iP,ilamb)-p(ip)*H(iH)-kappam(ilamb)*min(1d0,100d0*max(0d0,L(iLL))) <= cmin) then 
							V_b(1,1,iB,iExo,iP,iD6)=death
							LLi_b(1,1,iB,iExo,iP,iD6)=1
							HHi_b(1,1,iB,iExo,iP,iD6)=1
							BB_b(1,1,iB,iExo,iP,iD6)=0d0
							c_b(1,1,iB,iExo,iP,iD6)=cmin
						end if 

					end do
				end do
			end do
        end do
    end do
	!$OMP END DO nowait
    !$OMP END PARALLEL
    end subroutine sub_bellman_bor
end module mod_bellman_new

