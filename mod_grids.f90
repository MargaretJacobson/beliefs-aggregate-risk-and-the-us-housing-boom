!       This file contains the module mod_grids that sets up grids

module mod_grids
use omp_lib
    use parameters
	use mod_globals
    use mod_functions
    use mod_matlab
    use mod_interp
    implicit none 
contains
    
    subroutine sub_grids(iJ,p,rental,FYY,Fee,L,H,Y,chi,R,re_rent,re_own,LTV,PTI,min_pay,ch_rent,ch_own,loge,bheloc,q,Trans)
    
    implicit none       
    real(8), dimension(:), intent(in) :: p,Trans
	integer, intent(in) :: iJ
    real(8), dimension(:,:), intent(in) :: FYY, loge
	real(8), dimension(:,:), intent(in) :: rental
    real(8), dimension(:,:,:), intent(in) :: Fee
	real(8), dimension(:,:,:,:,:,:), intent(inout) :: q
    real(8), dimension(:), intent(in):: L,H,Y,chi,R
    real(8), dimension(:,:,:,:), intent(out) :: re_rent
    real(8), dimension(:,:,:,:,:,:,:), intent(out) :: re_own
    real(8), dimension(:,:,:,:), intent(out) :: ch_rent
    real(8), dimension(:,:,:,:,:,:,:), intent(out) :: ch_own
    integer, dimension(:,:,:), intent(out) :: PTI
    integer, dimension(:), intent(out) :: min_pay
    integer, dimension(:,:), intent(out) :: LTV
    integer :: iY, iYY, iM, iLL, iHH, iH, iL,iE,iB,iBB,iR,iAgg,ilamb,next,iP
	real(8), dimension(nH,nP),intent(out) :: bheloc
    
    !print*, 'thread', omp_get_thread_num() 

  
   ! print*, 'interest rates, beta, r_lend, q_lend, r_borrow, q_borrow'
   ! print*, beta, r_lend, q_lend, r_borrow, q_borrow



	do iP=1,nP
		bheloc(:,iP)=-p(iP)*H*.2d0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!   Setting up wealth/choice grids
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do iE=1,nE
			do iY=1,nY    
				!if (iJ<Jret) then
				!	print*, 'iJ,iE, iY, income',iJ, iE,iY, chi(iJ),loge(iJ,iE),exp(chi(iJ)+loge(iJ,iE))*Y(iY), income(chi(iJ),loge(iJ,iE),Y(iY)),FnTax(income(chi(iJ),loge(iJ,iE),Y(iY))), income(chi(iJ),loge(iJ,iE),Y(iY))-FnTax(income(chi(iJ),loge(iJ,iE),Y(iY)))
				!end if		
				do iH=1,nH
					do iL=1,nL
						do iB=1,nB
							if (iJ<Jret) then
								re_rent(iB,iY,iE,iP)=B(iB)+income(chi(iJ),loge(iJ,iE),Y(iY))-FnTax(income(chi(iJ),loge(iJ,iE),Y(iY)))+Trans(iP)
								re_own(iL,iB,iH,iY,iE,iP,2:3)=-H(iH)*L(iL)*(1d0+rm)+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income(chi(iJ),loge(iJ,iE),Y(iY))-FnTaxm(income(chi(iJ),loge(iJ,iE),Y(iY)),(1d0+rm)*L(iL)*H(iH),minval((/0d0,B(iB)/)))+Trans(iP)
								re_own(iL,iB,iH,iY,iE,iP,1)=B(iB)+income(chi(iJ),loge(iJ,iE),Y(iY))-FnTax(income(chi(iJ),loge(iJ,iE),Y(iY)))+Trans(iP) !NO outstanding mortgage balance
							else if (iJ>=Jret ) then
								re_rent(iB,iY,iE,iP)=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP)
								re_own(iL,iB,iH,iY,iE,iP,2:3)=-H(iH)*L(iL)*(1d0+rm)+p(iP)*(-delta-tauh)*H(iH)+B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTaxm(income_ret(chi(Jret-1),iE,mean(Y)),(1d0+rm)*L(iL)*H(iH),minval((/0d0,B(iB)/)))+Trans(iP)
								re_own(iL,iB,iH,iY,iE,iP,1)=B(iB)+income_ret(chi(Jret-1),iE,mean(Y))-FnTax(income_ret(chi(Jret-1),iE,mean(Y)))+Trans(iP)
							end if
						end do
					end do
				end do
			end do
		end do
        
		LTV=1
		ch_rent=0d0
		ch_own=0d0

		do ilamb=1,nlamb
			do iLL=nL,1,-1
				LTV(ilamb,iP)=iLL
				! print*,'LTV', iLL, iY, iJ, LTV(iY,iJ),LL(iLL),p*lambdaLTV(iY)
				if (L(iLL)<=p(ip)*lambdaLTV(ilamb)) exit
				
			end do
		end do


! Payment to income constraints put an upper limit on loans 
		do iL=1,nL
			do iLL=1,nL
				! Don't need to scale mortgages by housing because this only holds for those who stay and
				! their house size won't change and hence these terms will cancel from both sides 
				if (L(iLL)<=(1d0+rm)*L(iL)-((rm*(1d0+rm)**(T+1-iJ+4))/((1d0+rm)**(T+1-iJ+4)-1d0))*L(iL)) then 
					min_pay(iL)=iLL                                    
				end if
			!print*, 'min pay', iLL,iL,iJ, min(iL,iJ)
			!I think min pay is fine for 2 and 3, but for 1 it needs to be a function of the mortgage
			end do
		!PTI should become a function of HH
			do iAgg=1,nY*nlamb
				ilamb=(1-perfect_corr)*ceiling(real(iAgg)/real(nY))+perfect_corr*iAgg
				iY=(1-perfect_corr)*(floor(real(iAgg)/(real(nY)*real(ilamb)))+1)+perfect_corr*iAgg
				do iE=1,nE
					do iH=1,nH
						do iLL=1,nL
							if (iJ<Jret) then
								if (L(iLL)*H(iH)<=lambdaPTI(ilamb)*(((1d0+rm)**(T-iJ+4)-1d0)/(rm*(1d0+rm)**(T-iJ+4)))*income(chi(iJ),loge(iJ,iE),Y(iY))) then
									PTI(iH,iAgg,iE)=iLL
								end if
							else if (iJ>=Jret) then
								if (L(iLL)*H(iH)<=lambdaPTI(ilamb)*(((1d0+rm)**(T-iJ+4)-1d0)/(rm*(1d0+rm)**(T-iJ+4)))&
*income_ret(chi(Jret-1),iE,mean(Y))) then
                            PTI(iH,iAgg,iE)=iLL
								end if
							end if
					!	print*, 'PTI',iLL,iH,iE,iY,iL,iJ, PTI(iH,IY,iE,iJ),lambdaPTI(iY)*(((1+rm)**(T+1-iJ)-1d0)/(rm*(1+rm)**(T+1-iJ)))&
!*exp(chi(iJ)+loge(iJ,iE))*Y(iY)/H(iH), LL(iLL)
						end do
					end do
				end do
			end do 
		end do
        

		do iE=1,nE
			do iAgg=1,nY*nlamb
				ilamb=(1-perfect_corr)*ceiling(real(iAgg)/real(nY))+perfect_corr*iAgg
				iY=(1-perfect_corr)*(floor(real(iAgg)/(real(nY)*real(ilamb)))+1)+perfect_corr*iAgg
				do iR=1,nR
					ch_rent(iR,iAgg,iE,iP)=rental(iAgg,iP)*R(iR)
				end do
				do iLL=1,nL
					do iH=1,nH   
						do iB=1,nB
							if (iJ<Jret) then 
								if (L(iLL)>min(p(iP)*lambdaLTV(ilamb),(lambdaPTI(ilamb)/H(iH))*(((1d0+rm)**(T-iJ+4)-1d0)/(rm*(1d0+rm)**(T-iJ+4)))*income(chi(iJ),loge(iJ,iE),Y(iY)))) then 

									q(iLL,iB,iH,iAgg,iE,iP)=0d0
								end if
							else
								if (L(iLL)>min(p(ip)*lambdaLTV(ilamb),(lambdaPTI(ilamb)/H(iH))*(((1d0+rm)**(T-iJ+4)-1d0)/(rm*(1d0+rm)**(T-iJ+4)))*income_ret(chi(Jret-1),iE,mean(Y)))) then 
									q(iLL,iB,iH,iAgg,iE,iP)=0d0
								end if
							end if 
						! Need to interpolate over B so it is missing 
						    							
								ch_own(iLL,iB,iH,iAgg,iE,iP,1)=p(ip)*H(iH)-q(iLL,iB,iH,iAgg,iE,iP)*L(iLL)*H(iH)+kappam(ilamb)
								ch_own(iLL,iB,iH,iAgg,iE,iP,3)=-q(iLL,iB,iH,iAgg,iE,iP)*L(iLL)*H(iH)+kappam(ilamb)
							!if (q(iLL,iB,iH,iAgg,iE)/=1d0) then
							!print*, iE,iAgg,iLL,iH,iB, q(iLL,iB,iH,iAgg,iE)
							!end if 
						!else if (LL(iLL)>p*lambdaLTV(ilamb)) then
						!	ch_own(iLL,iB,iH,iAgg,iE,1)=p*H(iH)-q(iLL,iB,iH,iAgg,iE)*p*lambdaLTV(ilamb)+kappam(ilamb)
						!	ch_own(iLL,iB,iH,iAgg,iE,3)=-q(iLL,iB,iH,iAgg,iE)*p*lambdaLTV(ilamb)+kappam(ilamb)
						!end if     
							ch_own(iLL,iB,iH,iAgg,iE,iP,2)=-L(iLL)*H(iH)
						end do
					end do
                end do
            end do
        end do
    end do

	


        

    end subroutine sub_grids
end module mod_grids    
