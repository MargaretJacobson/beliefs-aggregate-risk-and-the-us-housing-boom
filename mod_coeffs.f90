 module mod_coeffs
    use parameters
    use mod_matlab
    use mod_discreteAR1
    implicit none 
contains
    
    subroutine sub_coeffs(a,coeffs)
    
    implicit none        
    real(8), dimension(:,:,:), intent(out) :: a
    integer ::  iM, ii, iY,iYY

    logical, intent(in) :: coeffs

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initializing Price Coefficients 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        a=0d0

    if (coeffs .eqv. .False. .and. type_sim<2) then
        a(1,1,1)=0		
        a(1,2,1)=0
        a(2,1,1)=0
        a(2,2,1)=0
            
        a(1,1,2)=1d0
        a(1,2,2)=1d0
        a(2,1,2)=1d0
        a(2,2,2)=1d0
    
    elseif (coeffs .eqv. .True. .or. type_sim==2) then    

		if (nlamb==2 .and. no_pref==1) then
		print*, 'nlamb=2'
			open(3,file='coefficients2.csv',status='old',action='read')
			do iM=1,mom+1
				do iY=1,nY*nlamb
					do iYY=1,nY*nlamb
						read(3,*) a(iY,iYY,iM)
						print*, 'nlamb: coeff2',nlamb, iM,iY,iYY, a(iY,iYY,iM)
					!read(3,*) a(1,2,iM,Aloc)
					!read(3,*) a(2,1,iM,Aloc)
					!read(3,*) a(2,2,iM,Aloc)
					end do
				end do
			end do
			close(3)
		else if (nlamb==2 .and. no_pref==0) then 
			open(5,file='coefficients3.csv',status='old',action='read')
			do iM=1,mom+1
				do iY=1,nY*nlamb*nphi
					do iYY=1,nY*nlamb*nphi
						read(5,*) a(iY,iYY,iM)
						print*, 'nlamb: coeff3',nlamb, iM,iY,iYY, a(iY,iYY,iM)
					end do
				end do
			end do
			close(5)
			
		
		else if (nlamb==1) then
		print*, 'nlamb=1'
			open(4,file='coefficients1.csv',status='old',action='read')
			do iM=1,mom+1
				do iY=1,nY*nlamb
					do iYY=1,nY*nlamb
					print*, 'nlamb: coeff1',nlamb, iM,iY,iYY 
						read(4,*) a(iY,iYY,iM)
					print*, 	a(iY,iYY,iM)
					end do
				end do
			end do
			close(4)		
		end if 
         
       
           
       
    end if
    end subroutine sub_coeffs
end module mod_coeffs