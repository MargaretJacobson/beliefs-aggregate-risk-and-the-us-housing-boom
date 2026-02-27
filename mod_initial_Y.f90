!       This file contains the module mod_initial_dist that creates the
!       Initial distribution for Krussel-Smith Siulation

module mod_initial_Y
    use parameters
    use mod_matlab
    use mod_interp
    implicit none 
contains
    subroutine sub_initial_Y(Y,FYY,YT,AggT_loc,YT_loc,LambT_loc)

    implicit none       
    integer :: sizeseed
    real(8), dimension(:), intent(in) :: Y
    integer :: it, r1, col2, iY, iYY, sim_start, iter
    real(8), dimension(:), intent(out):: YT
    real(8), dimension(T_sim) :: PT
    integer, dimension(:), intent(out):: AggT_loc,YT_loc,LambT_loc
    real(8), dimension(nY*nlamb*nPhi,nY*nlamb*nPhi+1):: Aggprob_CDF
    real(8), dimension(:,:), intent(in) :: FYY
	  integer, allocatable :: seed(:), old(:)
	  real(8) :: num
   logical :: data_in
	YT_loc=0
	YT=0d0
AggT_loc=0
LambT_loc=0
	sim_start=1

	call random_seed(size=sizeseed)
	allocate(seed(sizeseed))
	allocate(old(sizeseed))
	seed=123456789

	call random_seed(PUT=seed)

	! Setting AggT_loc to the worst state:
	AggT_loc=nY*nlamb*nphi
	lambT_loc=ceiling(real(AggT_loc)/(real(nY)*real(nPhi)))
	YT_loc=floor(real(AggT_loc)/(real(nPhi)*real(nY)*real(lambT_loc)))+1  
	!print*, 'lambT_loc', ceiling(real()/(real(nY)*real(nPhi)))
	!print*, 'YT-loc', floor(real(12)/(real(nPhi)*real(nY)*real(2)))+1  
	! Adding in the indices from KMV
	! Need 4 to 1, 2 to 3, 3 to 2 and 4 to 1
	
	if (ny*nlamb==4 .and. no_Agg==0 .and. no_Agg_lend==0 .and. perfect_corr==1 .and. no_pref==1) then
		AggT_loc(101:105)=1			!4 to 1
		AggT_loc(106:110)=3			!1 to 3
		AggT_loc(111:115)=2			!3 to 2
		AggT_loc(116:120)=3			!2 to 3
		
		AggT_loc(161:162)=1			!4 to 1
		AggT_loc(163:165)=3			!1 to 3
		AggT_loc(166:167)=2			!3 to 2
		AggT_loc(168:170)=3			!2 to 3
		
		AggT_loc(200:202)=2			!3 to 2
		AggT_loc(203:205)=3			!2 to 3
		
		AggT_loc(221:225)=1			! 4 to 1
				! 1 to 4
		AggT_loc(291:295)=1			! 1 to 4
		
		
		!added by MJ				! 1 to 4
		!AggT_loc(300:305)=1			

		sim_start=206
	else if ((ny*nlamb==4 .or. nY*nlamb==2) .and. no_Agg==0) then ! .and. no_Agg_lend==0) then
		!sim_start=399
		sim_start=155
		AggT_loc(191:199)=(ny*nlamb*nphi)
		
	end if 
	do it=1,sim_start
		if (modulo(AggT_loc(it),2)==0) then ! .or. AggT_loc(it)==3 .or. AggT_loc(it)==5 .or. ) then
			YT_loc(it)=2
		else
			YT_loc(it)=1

		end if
		if (AggT_loc(it)<3 .or. (AggT_loc(it)>4 .and. AggT_loc(it)<7) .or. (AggT_loc(it)>8 .and. AggT_loc(it)<11)) then ! iLamb=1 when iAgg<3
			lambT_loc(it)=1
		else
			lambT_loc(it)=2
		end if
	
		!YT(it)=Y(YT_loc(it))
	end do

	! Creating transition matrix
    Aggprob_CDF(:,1)=0d0
    do iY=1,nY*nlamb*nPhi
        do iYY=1,nY*nlamb*nPhi
            Aggprob_CDF(iY,iYY+1)=Aggprob_CDF(iY,iYY)+FYY(iY,iYY)
		!	print*, iY,iYY,Aggprob_CDF(iY,iYY+1),Aggprob_CDF(iY,iYY),FYY(iY,iYY)
        end do
    end do
	
	
    

	

    !Simulating YT
	!print*, 'sim_start', sim_start
    do it=sim_start,T_sim-1
		call random_number(num)
		PT(it+1)=num
        col2=bsearch(PT(it+1),Aggprob_CDF(AggT_loc(it),:))
		if (it<=learn-1 .and. no_pref==0) then
			if (col2<5) then
				col2=8+col2
			else if (col2>4 .and. col2<9) then
				col2=4+col2
			end if 
			if (col2==9)  then
				col2=11
			else if (col2==10) then
				col2=12
			end if 
		end if 
		!print*, it,aggT_loc(it), Aggprob_CDF(AggT_loc(it),:),PT(it+1)
		AggT_loc(it+1)=col2
		
		lambT_loc(it+1)=ceiling(real(AggT_loc(it+1))/real(nY))-real(nY)*(ceiling(AggT_loc(it+1)/(real(nY)*real(nY)))-1)
		YT_loc(it+1)=(1-modulo(AggT_loc(it+1),2))*nY+modulo(AggT_loc(it+1),2)

	
	
    end do

	if (nY*nLamb==4 .and. no_Agg==0 .and. no_agg_lend==0 .and. perfect_corr==0 .and. no_pref==1) then ! .and. type_sim==2) then 	
	
		! Need to manually set pre-housing boom so it starts in AggT_loc==4
		iter=193
		do while (iter<199) !(iter<4347)
			AggT_loc(iter)=4
			YT_loc(iter)=2
			lambT_loc(iter)=2
			
			
			iter=iter+1
		end do
					AggT_loc(iter)=3
			YT_loc(iter)=1
			lambT_loc(iter)=2

	
		!iter=199 
		!	AggT_loc(iter)=2 !3
		!	YT_loc(iter)=2 !1
		!	lambT_loc(iter)=1 !2
			
			
			iter=iter+1
	
		do while (iter<204) !204 is the housing boom, 214 is the longer boom
			AggT_loc(iter)=1
			YT_loc(iter)=1
			lambT_loc(iter)=1
			iter=iter+1
		end do
		!iter=202
		!AggT_loc(iter)=3
		!YT_loc(iter)=1 ! 2
		!lambT_loc(iter)=2
		
		!iter=203
		!AggT_loc(iter)=3
		!YT_loc(iter)=1 ! 2
		!lambT_loc(iter)=2
		!iter=212
		AggT_loc(iter:iter+10)=4
		YT_loc(iter:iter+10)=2 ! 2
		lambT_loc(iter:iter+10)=2
		
		iter=290 !4904
		!do while (iter<295) !(iter<4909)
		!			AggT_loc(iter)=1
		!	YT_loc(iter)=1
		!	lambT_loc(iter)=1
		!	iter=iter+1
		!end do 

		if (T_sim>300) then 	
			iter=1000 !4904
			do while (iter<1006) !(iter<4909)
						AggT_loc(iter)=2
				YT_loc(iter)=2
				lambT_loc(iter)=1
				iter=iter+1
			end do 
		
		

	!! Adding in some extra data poitns for regressions
	if (no_MIT==1 .and. solve_for_coeffs==1) then
		AggT_loc(2014)=4
		AggT_loc(3617)=4
		
		AggT_loc(3838)=3
		
		!1 grids 
		AggT_loc(1110)=1
		AggT_loc(1111)=4
		
		
		AggT_loc(1500)=1
		AggT_loc(1503)=4
		
		AggT_loc(1600)=2
		AggT_loc(1601)=3
		
		AggT_loc(3000)=2
		AggT_loc(3001)=3
		
		AggT_loc(2000)=3
		AggT_loc(2001)=2
		
		AggT_loc(2200)=3
		AggT_loc(2201)=2
		
		
		AggT_loc(3200)=4
		AggT_loc(3201)=1
		
		AggT_loc(1200)=1
		AggT_loc(1201)=4
		
		AggT_loc(700)=1
		AggT_loc(701)=4
		
	
	end if 

	end if 	

	end if
	
	!Allows for preference shocks
	if (nY*nLamb==4 .and. no_Agg==0 .and. no_agg_lend==0 .and. perfect_corr==0 .and. no_pref==0) then ! .and. type_sim==2) then 	
	
		! Need to manually set pre-housing boom so it starts in AggT_loc==4
		iter=193
		do while (iter<199) !(iter<4347)
			AggT_loc(iter)=12 !4*npHi
			YT_loc(iter)=2
			lambT_loc(iter)=2
			
			
			iter=iter+1
		end do
			!		AggT_loc(iter)=3
			!YT_loc(iter)=1
			!lambT_loc(iter)=2

	
		iter=199 !4342
			AggT_loc(iter)=11 !3*npHi
			YT_loc(iter)=1
			lambT_loc(iter)=2
			
			
			iter=iter+1
			AggT_loc(iter)=9 !3*npHi
			YT_loc(iter)=1
			lambT_loc(iter)=1
			
			
			
			iter=iter+1
	
		do while (iter<204) !(iter<4347)
			AggT_loc(iter)=5
			YT_loc(iter)=1
			lambT_loc(iter)=1
			iter=iter+1
		end do
		
		iter=202
		
		AggT_loc(iter)=7
		YT_loc(iter)=1
		lambT_loc(iter)=2
		
		iter=203
				AggT_loc(iter)=7
		YT_loc(iter)=1
		lambT_loc(iter)=2
		
		AggT_loc(204:210)=12
		YT_loc(204:210)=2
		lambT_loc(204:210)=2
	
		if (T_sim>300) then
			iter=490 !4904
			do while (iter<495) !(iter<4909)
						AggT_loc(iter)=5
				YT_loc(iter)=1
				lambT_loc(iter)=1
				iter=iter+1
			end do 

			if (T_sim>600) then 	
				iter=1000 !4904
				do while (iter<1006) !(iter<4909)
							AggT_loc(iter)=10
					YT_loc(iter)=2
					lambT_loc(iter)=1
					iter=iter+1
				end do 
			
				AggT_loc(iter)=1
				YT_loc(iter)=1
				lambT_loc(iter)=1
				iter=iter+1
				
				AggT_loc(iter)=2
				YT_loc(iter)=2
				lambT_loc(iter)=1
				iter=iter+1
				
				if (no_MIT==1) then
					AggT_loc(2014)=4
					AggT_loc(3617)=4
					AggT_loc(3838)=3
				end if 
			end if 
		end if 		



		
		!AggT_loc(4909)=4
		!	YT_loc(4909)=2
		!	lambT_loc(4909)=2
	end if
 
   do it=1,T_sim
 !  print*, it,AggT_loc(it),YT_loc(it),lambT_loc(it)
		YT(it)=Y(YT_loc(it))
	! 	print*, it, AggT_loc(it),YT_loc(it),LambT_loc(it)

   end do


    end subroutine sub_initial_Y
end module mod_initial_Y
