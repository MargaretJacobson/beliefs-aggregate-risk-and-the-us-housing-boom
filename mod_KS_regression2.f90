!       This file contains the module that runs the Krusell Smith Regression

module mod_KS_regression2
    use parameters
    use mod_matlab
    implicit none 
contains
    
    subroutine sub_KS_regression(AggT_loc,preg,pstar,counts,a_new,R2p,R2p_all,print_reg)
        implicit none
        integer :: iY, it, iYY, iM, sim_start,index, regcounts
        integer, dimension(:,:), intent(out) :: counts
        real(8), dimension(nY*nlamb*npHi*nY*nlamb*nphi+nY*nlamb*nphi*mom) :: regsum
        integer, dimension(:), intent(in) :: AggT_loc
        real(8), dimension(:), intent(in) :: pstar
        real(8), dimension(:,:,:), intent(out) :: preg
        real(8), dimension(:,:,:), intent(out) :: a_new
        real(8), dimension(:,:), intent(out) :: R2p
		real(8), intent(out) :: R2p_all
        real(8), dimension(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom,T_sim-burn-2+1) :: Xreg
		! real(8), dimension(nY*nlamb*nY*nlamb+nY*nlamb*mom) :: a_reg
        real(8), dimension(nY*nlamb*nphi,nY*nlamb*nphi,T_sim-burn-2+1) :: Yhatp,TSSp,RSSp
		real(8), dimension(T_sim-burn-2+1) :: Yreg
        !real(8), dimension(nY*nlamb*nY*nlamb+nY*nlamb*mom,nY*nlamb*nY*nlamb+nY*nlamb*mom) :: XprimeX
       ! real(8), dimension(nY*nlamb*nY*nlamb+nY*nlamb*mom) :: XprimeY
        real(8), dimension(nY*nlamb*nphi,nY*nlamb*nphi) :: meanYp   
		real(8), allocatable :: resid(:)
		real(8) :: Xmatmul(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom,T_sim-burn-2+1),XprimeX(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom,nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom),XprimeY(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom)
		real(8) :: invXprimeX(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom,nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom)
		real(8), dimension(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi*mom) :: a_reg
        logical :: print_reg 
     
		counts=0
        preg=0d0 
        Yreg=0d0
        Xreg=1d0
        Yhatp=0d0
        Tssp=0d0
        RSSp=0d0
        XprimeX=0d0
        XprimeY=0d0
        a_new=0d0
					
		a_reg=0d0
		Xmatmul=0d0
			
		sim_start=burn
		
		allocate(resid(T_sim+1))
	
		! Initialization important, need zeros in regressor matrix
		Xreg=0d0
        !open(1, file = 'Regression_data.csv', status = 'replace')
        !do iY=1,nY*nlamb
        !    do iYY=1,nY*nlamb
                 
        !write(1,*) '________'
        !write(1,*) iY, iYY
        !write(1,*) '________'
        do it=sim_start,T_sim-2
			iY=AggT_loc(it)
			iYY=AggT_loc(it+1)
          !  print*, iY,iYY,counts(iY,iYY)    
			counts(iY,iYY)=counts(iY,iYY)+1
			preg(iY,iYY,counts(iY,iYY))=pstar(it)
					
			!! creating an index for today and tomorrow for dummy variables
			index=(iY-1)*(nlamb*nY*nphi)+iYY
					
			Yreg(it-sim_start+1)=log(pstar(it+1))
			Xreg(index,it-sim_start+1)=1d0
					
			! set for one or more moments
			do iM=2,mom+1
				Xreg(nlamb*nY*nphi*nlamb*nY*nphi*(iM-1)+iY,it-sim_start+1)=log(pstar(it))**(iM-1)
			end do
                       ! print*, it, Yreg(iY,iYY,it),Xreg(iY,iYY,it,2)
			!if (print_reg .eqv. .True.) then                             
			!	write(1,*) it, counts(iY,iYY), exp(Xreg(iY,iYY,counts(iY,iYY),2:)), exp(Yreg(iY,iYY,counts(iY,iYY)))
                            
			!end if

		!	print*, it, iY, Yreg(it-sim_start+1), index, Xreg(index,it-sim_start+1),Xreg(nlamb*nY*nphi*nlamb*nY*nphi*mom+iY,it-sim_start+1),Xreg(:,it-sim_start+1)
		!	print*, it,iY,iYY,counts(iY,iYY),preg(iY,iYY,counts(iY,iYY)),index,Yreg(it-sim_start+1),Xreg(nlamb*nY*nphi*nlamb*nY*nphi+iY,it-sim_start+1)
		end do
                !write(1,*) exp(mean(Yreg(iY,iYY,1:counts(iY,iYY)))), exp(mean(Xreg(iY,iYY,1:counts(iY,iYY),2)))
        !    end do
        !end do
       ! close(1)
		index=(nY*nlamb*nphi*nY*nlamb*nphi+nY*nlamb*nphi)
			regsum=sum(Xreg,dim=2)
			print*, 'regsum', regsum
			regcounts=0
			do iY=1,index
				print*, 'iY, regsum',iY, regsum(iY)
				if (abs(regsum(iY))>0) then
					regcounts=regcounts+1
				end if 
			end do
			print*, 'regcounts 1', regcounts
			!allocate(Xmatmul(regcounts,T_sim-burn-2+1))

			regcounts=0 
			do iY=1,index
				if (abs(regsum(iY))>0) then
					regcounts=regcounts+1
					Xmatmul(regcounts,:)=Xreg(iY,:)
				end if 
			
			end do
			print*, 'regcounts 2', regcounts
			!allocate(XprimeX(regcounts,regcounts))
			!allocate(XprimeY(regcounts))
            XprimeX(1:regcounts,1:regcounts)= matmul(Xmatmul(1:regcounts,:),transpose(Xmatmul(1:regcounts,:)))
            XprimeY(1:regcounts)=matmul(Xmatmul(1:regcounts,:), Yreg)
				
				do iY=1,regcounts
				print*, 'YprimeY,Xprimex',iY, XprimeY(iY), XprimeX(iY,1:regcounts)
				end do
                !print*, 'XprimeX, XprimeY', iY,iYY, XprimeX(iY,iYY,:,:),XprimeY(iY,iYY,:),sum(XprimeY(iY,iY,:)), sum(XprimeX(iY,iYY,:,:))
			!	        do iY=1,nY*nlamb
            !do iYY=1,nY*nlamb
            !    if (sum(XprimeY(iY,iY,:))==0d0 .and. sum(XprimeX(iY,iYY,:,:))==0d0) then
            !        do iM=1,mom+1
            !            if (iM==1)  then
            !                a_new(iY,iYY,iM)=0d0
            !            elseif (iM==2) then
            !                a_new(iY,iYY,iM)=1d0
            !            else
            !                a_new(iY,iYY,iM)=0d0
            !            end if
            !        end do
            !    else 
					
            a_reg(1:regcounts)=mldivide(XprimeX(1:regcounts,1:regcounts),XprimeY(1:regcounts))
			!print*, a_reg(1:regcounts)
			!invXprimeX(1:regcounts,1:regcounts)=inv(XprimeX(1:regcounts,1:regcounts))	
			!				do iY=1,regcounts
			!	print*, 'YprimeY,invXprimex', XprimeY(iY), invXprimeX(iY,:)
			
			!	end do
			!print*, 'a_reg',a_reg
			
			!do iY=1,regcounts
			!	a_reg(iY)=dot_product(invXprimeX(iY,1:regcounts),XprimeY(1:regcounts))
			!	print*, 'iY a_reg(iy)', iY, a_reg(iY)
			!end do 
			!print*, 'a_reg2', matmul(inv(XprimeX(1:regcounts,1:regcounts))*XprimeY(1:regcounts))
			regcounts=0
			do iY=1,nY*nlamb*nphi
				do iYY=1,nY*nlamb*nphi
					index=(iY-1)*(nlamb*nY*nphi)+iYY
					if (abs(regsum(index))>0) then
						regcounts=regcounts+1
						a_new(iY,iYY,1)=a_reg(regcounts)
					else
						a_new(iY,iYY,1)=0d0
					end if 
					print*, 'allocate coeffs', index, regsum(index),regcounts, a_reg(regcounts)
				end do
			end do
			
			do iY=1,nY*nlamb*nphi 
				do iM=2,mom+1
				
					if (abs(regsum(nlamb*nY*nphi*nlamb*nY*nphi+iY))>0) then  
					regcounts=regcounts+1
					a_new(iY,:,iM)=a_reg(regcounts)	
					else 
					a_new(iY,:,iM)=1d0
					end if
					print*, 'allocate coeffs',regcounts ,a_reg(regcounts)
				end do
			end do
			!		print*

             !   end if
                
            print*, 'Made it here, sim_start', sim_start
		
		do it=sim_start,T_sim-2
			iY=AggT_loc(it)
			iYY=AggT_loc(it+1)
			!print*, 'it, resid', it, resid(it)
			!print*, 'Pstar',log(pstar(it+1))
			!print*, 'iY,iYY', iY, iYY
			!print*, 'a_new', a_new(iY,iYY,:)    
			resid(it+1)=log(pstar(it+1))-dot_product((/1d0,log(pstar(it))/),a_new(iY,iYY,:))
          !      
		end do
		print*, sum(resid(sim_start+1:T_sim-1)**2.0d0),sum((log(pstar(sim_start+1:T_sim-1)))) !-mean(log(pstar(sim_start+1:T_sim-1))))**2.0d0)
		! Note when there is no intercept, the r^2 model is forced through origin so we do not need to compute sum((y-mean(y))^2) only sum((y)^2)
		R2p_all=1d0-sum(resid(sim_start+1:T_sim-2)**2.0d0)/(sum((log(pstar(sim_start+1:T_sim-2)))**2)) !-mean(log(pstar(sim_start+1:T_sim-2))))**2.0d0))
             R2p=R2p_all
		deallocate(resid)
		!deallocate(Xmatmul)
		!deallocate(XprimeX)
		!deallocate(XprimeY)
    end subroutine sub_KS_regression
end module mod_KS_regression2