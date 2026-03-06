module mod_bellman_bequest
    !use omp_lib
    use parameters
	use mod_matlab
	use mod_globals
	use mod_functions
    implicit none
contains

    subroutine sub_bellman_bequest(V_b,c_b,BB_b,HHi_b,LL_b,R,H,L,e,p,Y,chi,Trans,rental)
	implicit none  
    real(8), dimension(:,:,:,:,:,:), intent (out) :: c_b
    integer, dimension(:,:,:,:,:,:), intent (out) :: HHi_b
	real(8), dimension(:,:,:,:,:,:), intent (out) :: BB_b,LL_b
	real(8), dimension(:,:,:,:,:,:), intent (out) :: V_b
    real(8), dimension(:), intent (in) :: R,H,L,p,Y,chi
	real(8), dimension(:,:), intent (in) :: rental,Trans
	real(8), intent (in) :: e
    integer :: iY,iR,iE,iB,iBB, iD5, iD4,bbar,iD1, iL,iH,iD2,iExo,ilamb,iP,iAgg,iphi
	real(8) :: BBmin1,BBmax1,BBmid1,Vtemp1a,Vtemp1b,btemp1beq,pm,mp!,phibeq
	real(8) :: BBmin2,BBmax2,BBmid2,Vtemp2a,Vtemp2b,btemp2	
    real(8) :: BBmin4,BBmax4,BBmid4,Vtemp4a,Vtemp4b,btemp4
	real(8) :: BBmin5,BBmax5,BBmid5,Vtemp5a,Vtemp5b,btemp5
	
	DOUBLE PRECISION, EXTERNAL          :: golden
    HHi_b=1
	c_b=cmin
	BB_b=0d0
	LL_b=0d0

    iD5=5
    iD4=4
	iD1=1
	iD2=2
    
    Vtemp4a=death
    Vtemp5a=deathdefault
	Vtemp4b=death
	Vtemp5b=deathdefault
	Vtemp1a=death
	Vtemp1b=death
	Vtemp2a=death
	Vtemp2b=death
	V_b=death
	V_b(:,:,:,:,:,5)=deathdefault
	eiJ=e
	btemp1beq=0d0
	btemp2=0d0
	btemp4=0d0
	btemp5=0d0
	pm=0d0
	mp=0d0


	do iP=1,nP
		
		!OMP Parallel Default(Shared) private(iE,iY,iB,iR,BBmin4,cval4,BBmax4,BBmid4,Wval,Rval4,Vtemp4a,Vtemp4b,btemp4,BBmin5,cval5,BBmax5,BBmid5,Rval5,Vtemp5a,Vtemp5b,btemp5)
		!OMP DO collapse(3) schedule(dynamic)    
		do iExo=1,nE*nY*nlamb*nPhi
			
			iE=iExo-nE*(ceiling(real(iExo)/real(nE))-1)
			iAgg=(1-perfect_corr)*ceiling(real(iExo)/real(nE))+perfect_corr
			iPhi=(1-perfect_corr)*ceiling(real(iExo)/(real(nE)*real(nlamb)*real(nY)))+perfect_corr
			ilamb=(1-perfect_corr)*ceiling(real(iExo)/(real(nY)*real(nE)))-real(nY)*(real(iPhi)-1d0)+perfect_corr
			iY=(1-perfect_corr)*ceiling(real(iExo)/real(nE))-nY*(ilamb-1)-real(nY)*real(nlamb)*(real(ipHi)-1d0)+perfect_corr

			
			phibeq=phi_grid(iPhi)
			q_lend=q_lend_grid(ilamb)
			q_borrow=q_borrow_grid(ilamb)
!print*, iP,iExo,iE,iAgg,iPhi,ilamb,iY,phibeq, q_lend,q_borrow			
			

			do iB=1,nB
				do iR=1,nR
				!	print*, 'before', iR, nR
					BBmin4=0d0
					cval4=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(iR)
					
					BBmax4=(cval4-cmin)/q_lend_grid(ilamb)
					BBmid4=(BBmax4+BBmin4)/2.0d0
					
					Rval4=R(iR)
						
					!	print*, 'before: iP,iExo,iB,iR',iP,iExo,iE,iAgg,iphi,ilamb,iY, iB, cval4,iR,nR,R, R(iR),Rval4
						
						
					if(BBmax4 .gt. BBmin4) then
						!if (cval4 > cmin) then
						Vtemp4a=Value4beq(BBmin4)
						Vtemp4b=golden(BBmin4,BBmid4,BBmax4,Value4beq,1d-5,btemp4)
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
						
						
					if (-Vtemp4a>V_b(1,1,iB,iExo,iP,iD4)) then
						V_b(1,1,iB,iExo,iP,iD4)=-Vtemp4a
						HHi_b(1,1,iB,iExo,iP,iD4)=iR
						BB_b(1,1,iB,iExo,iP,iD4)=btemp4
						c_b(1,1,iB,iExo,iP,iD4)=cval4
							
					end if
						

										
				end do   

				! Do not need R loop for defaulters
				BBmin5=0d0
				cval5=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(1)
				BBmax5=(cval5-cmin)/q_lend_grid(ilamb)
				BBmid5=(BBmax5+BBmin5)/2.0d0
				
				Rval5=R(1)
				if (bbmax5 .gt. bbmin5) then
					!if (cval5 > cmin) then
					Vtemp5a=Value5beq(BBmin5)
					Vtemp5b=golden(BBmin5,BBmid5,BBmax5,Value5beq,1d-5,btemp5)
						!print*, 'Vtemp5a, Vtemp5B, btemp', -Vtemp5a, -Vtemp5b,btemp5
					if(Vtemp5b<Vtemp5a) then
						Vtemp5a=Vtemp5b
							!print*,'btmp', btemp
							!Vtemp4a=Vtemp4b
							cval5=cval5-q_lend_grid(ilamb)*btemp5
					else
						btemp5=0d0
					end if
				else 
					Vtemp5a=-deathdefault !KMV set this value to be higher for foreclosure
					btemp5=0d0
					cval5=cmin
				end if	
						
						
				if (-Vtemp5a>V_b(1,1,iB,iExo,iP,iD5)) then ! .and. upper3(iL)>=LL(iLL)) then
					
					V_b(:,:,iB,iExo,iP,iD5)=-Vtemp5a
					HHi_b(:,:,iB,iExo,iP,iD5)=1
					BB_b(:,:,iB,iExo,iP,iD5)=btemp5
					c_b(:,:,iB,iExo,iP,iD5)=cval5
					!print*, iB,iR, -Vtemp5a, V_b(iB,1,1,iY,iE,iD5), btemp5,cval5
				end if  


						
	!print*, 'defaulter: giJ,giAgg,bc,hcc,ph,Pr,net income, gliq', iE,iAgg,iB,p(iP),B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP,ilamb)-rental(iAgg,iP)*R(1), btemp5, BBmid5,BBmax5,V_b(1,1,iB,iAgg,iE,iP,iD5),-Vtemp5a
				

				do iH=1,nH
					do iL=1,nL
						
						do iR=1,nR
							BBmin1=0d0
cval1beq=B(iB)+(1d0-delta-tauh-cost)*p(iP)*H(iH)-L(iL)*H(iH)*(1d0+r_borrow_grid(ilamb))+income_ret(chi(Jret-1),iE,mean(Y))-FnTaxM(income_ret(chi(Jret-1),iE,mean(Y)),H(iH)*L(iL),min(0.0d0,B(ib)),r_borrow_grid(ilamb))+Trans(iP,ilamb)-rental(iAgg,iP)*R(iR)
		BBmax1=(cval1beq-(-1d0/death))/q_lend_grid(ilamb)
						BBmid1=(BBmax1+BBmin1)/2.0d0
							
						Hval1=R(iR)
						
						if(BBmax1 .gt. BBmin1) then
							
									if (cval1beq-BBmin1*q_lend_grid(ilamb)>=cmin) then
			
			Vtemp1a=Value1beq(BBmin1)
			Vtemp1b=golden(BBmin1,BBmid1,BBmax1,Value1beq,1d-5,btemp1beq)
			
					else
			Vtemp1a=-death
			Vtemp1b=-death
		end if

							
						
							if(Vtemp1b<Vtemp1a) then
								Vtemp1a=Vtemp1b
								!print*,'btmp', btemp
								!Vtemp4a=Vtemp4b
								cval1beq=cval1beq-q_lend_grid(ilamb)*btemp1beq	!cval(1) has dimenson (nB) and we don't need this for this problem
							else
								btemp1beq=0d0
							end if
						else 
							Vtemp1a=-death
							btemp1beq=0d0
							cval1beq=-1d0/death
						end if		
				if (-Vtemp1a>V_b(iL,iH,iB,iExo,iP,iD1)) then
							V_b(iL,iH,iB,iExo,iP,iD1)=-Vtemp1a
							HHi_b(iL,iH,iB,iExo,iP,iD1)=iR
							BB_b(iL,iH,iB,iExo,iP,iD1)=btemp1beq
							c_b(iL,iH,iB,iExo,iP,iD1)=cval1beq

						end if

						

					end do
bheloc2=-p(iP)*H(iH)*0.2d0
				! Here we pay the mortgage
							!BBmin2=bheloc2 !should be 0d0
							BBmin2=0d0
							!These households are only choosing savings b'
							pm=MinPay(L(iL)*H(iH),T,r_borrow_grid(ilamb))
							if (iL .eq. 1) pm=0d0
							mp=(H(iH)*L(iL))*(1.0d0+r_borrow_grid(ilamb))-pm
							cval2=-pm+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),L(iL)*H(iH),minval((/0d0,B(iB)/)),r_borrow_grid(ilamb))+Trans(iP,ilamb)
							!Other if for the bequest
							other=p(iP)*H(iH)*(1d0-cost)-mp*(1d0+r_borrow_grid(ilamb)) !it should be L(iL)*(1d0+rm)*H(iH) and also need (-delta-tauh)*housevalue
							BBmax2=(cval2-cmin)/q_lend_grid(ilamb)
							BBmid2=(BBmax2+BBmin2)/2.0d0
							!need to make this a function of H
							Hval2=H(iH)*utilitycost
						
							!print*, 'iB,iR', iB, iR, cval4
						
						
							if(BBmax2 .gt. BBmin2) then
							!if (cval2 > cmin) then
								Vtemp2a=Value2beq(0d0) !should be bbmin2?
								Vtemp2b=golden(BBmin2,BBmid2,BBmax2,Value2beq,1d-5,btemp2)
								!print*, 'Vtemp2a, Vtemp2B, btemp', -Vtemp2a, -Vtemp2b,btemp2
								if(Vtemp2b<Vtemp2a) then
									Vtemp2a=Vtemp2b
									!print*,'btmp', btemp
									!Vtemp4a=Vtemp4b
									cval2=cval2-q_lend_grid(ilamb)*btemp2
								else
									btemp2=0d0
								end if
							else 
								Vtemp2a=-death
								btemp2=0d0
								cval2=cmin
							end if		

			
							if (-Vtemp2a>V_b(iL,iH,iB,iExo,iP,iD2)) then
								V_b(iL,iH,iB,iExo,iP,iD2)=-Vtemp2a
								HHi_b(iL,iH,iB,iExo,iP,iD2)=iH
								BB_b(iL,iH,iB,iExo,iP,iD2)=btemp2
								LL_b(iL,iH,iB,iExo,iP,iD2)=mp/H(IH)
								c_b(iL,iH,iB,iExo,iP,iD2)=cval2
								!print*, iL,iB,iH, -Vtemp2a, V_b(iL,iH,iB,iY,iE,iD2), btemp2,cval2
							end if
		

				
				
					end do
				end do
            end do
        end do
    end do
    !OMP END DO
    !OMP END PARALLEL
!print*, 'R, nR end', nR, R
    end subroutine sub_bellman_bequest
end module mod_bellman_bequest

