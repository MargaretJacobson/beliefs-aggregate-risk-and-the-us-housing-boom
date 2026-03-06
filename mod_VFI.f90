module mod_VFI
    use omp_lib
    use parameters
	use mod_globals
    use mod_matlab
    use mod_interp
	use mod_bellman_bequest
	use mod_continuation_value
    use mod_bellman_new
    use mod_EEE
    implicit none 
contains
    recursive subroutine sub_VFI(a_raw,nA,nA2,p,countP,FYY,Fee,L,H,Y,chi,e,loge,R,c_b,LL_b,BB_b,HH_b,V_b,pprimein,theta_pprime,ilo_pprime,Trans,rental,W_own,q) !,EEEq, EEEp)
    implicit none 
    integer, intent(in) :: countP,nA,nA2
    !real(8), dimension(nY*nlamb,nY*nlamb,mom+1) :: a 
    real(8), dimension(:,:,:), intent(in) :: a_raw
    real(8), dimension(:,:), intent(in) :: FYY
    real(8), dimension(:,:,:), intent(in) :: Fee
    real(8), dimension(:), intent(in):: L,H,Y,chi,e,p,R
    real(8), dimension(:,:), intent(in) :: loge,Trans
    integer :: iY,iYY,iJ,iLL,iHH,iP,it,iH,iL,iE,iA1,iA2,iM,Aloc,iD,iB, iAgg,iA,loc,ilamb,iExo,iPhi
 	real(8), dimension(:,:,:,:,:,:), intent(in) :: pprimein
    real(8), allocatable :: W_b(:,:,:,:,:,:,:,:,:)
	real(8), dimension(:,:,:,:), intent(in) :: rental
	real(8), dimension(:,:,:,:,:,:,:,:), intent(out) :: q
	integer, dimension(:,:,:,:,:,:,:), intent(in) :: ilo_pprime
	real(8), dimension(:,:,:,:,:,:,:), intent(in) :: theta_pprime
	integer, dimension(nA,nA) :: bariL, bariH
    real(8), dimension(:,:,:,:,:,:,:,:,:), intent(out) :: V_b, c_b,HH_b,LL_b, BB_b
	real(8), allocatable :: q_e(:,:,:,:,:,:,:,:) !q(:,:,:,:,:,:,:,:),
    integer, allocatable :: LLi_b(:,:,:,:,:,:,:,:,:), HHi_b(:,:,:,:,:,:,:,:,:)
	integer, allocatable :: W_loc_own(:,:,:,:,:,:,:,:)
	integer, allocatable :: W_loc_rent(:,:,:,:,:,:)
	!real(8), allocatable :: W_own(:,:,:,:,:,:,:,:,:)
	real(8) :: inc
	real(8), dimension(:,:,:,:,:,:,:,:,:), intent(out) :: W_own

	allocate(W_loc_own(nL,nH,nB,nY*nlamb*nE*nPhi,nP,T,nA,nA))
	allocate(W_loc_rent(nB,nY*nlamb*nE*nphi,nP,T,nA,nA))
	allocate(W_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,2))
    allocate(LLi_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
    allocate(HHi_b(nL,nH,nB,nY*nlamb*nE*nphi,countP,T,nA,nA,nD+2))
	allocate(q_e(nL,nH,nB,nY*nlamb*nE*nphi,countP,2,nA,nA))
	
	print*, 'R: VFI 0', R
    Aloc=(nA+1)/2
    bariL=nL
	bariH=nH
    V_b=death
    W_b=death  
	LLi_b=1
	HHi_b=1
	!rental=0d0
	LL_b=0d0
	BB_b=0d0
	HH_b=0d0
	C_b=cmin
	W_loc_own=0
	W_loc_rent=0
	
	!a=0d0
	q=0d0
	q_e=0d0
	

  if (nA==1) then
		loc=1
	else
		loc=2
	end if 
	

    print*, 'max threads pre parallel, thread limit', omp_get_max_threads()
    print*, '# threads in use pre parallel should==1', omp_get_num_threads()
    print*, 'NA, NA2', nA,nA2

!!
	!$OMP Parallel Do default(shared) private(iA,iA1,iA2,iJ,iP,iL,iB,iH,iY,iYY,iE,iD,ilamb,iExo,iAgg,iPhi) if(NA>1) num_threads(9)
	do iA=1,(nA*nA)
		if (nA==1) then
			iA1=1		!Need to define these to make the code run correctly
			iA2=1
			
		else
			
			iA2=ceiling(real(iA)/real(nA))
			iA1=iA-nA*(iA2-1)
			
		end if
			print*, 'Start: iA: 1,2', iA, iA1,iA2, omp_get_thread_num(), omp_get_num_threads()


		
			
        !print*, 'VFI loop: iA1, iA1:', iA1, iA2, 'In parallel? Thread #, size of current team', omp_in_parallel(),omp_get_thread_num(),omp_get_num_threads(), 'max threads to be used in nest', omp_get_max_threads()
                    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   Value Function Iteration  
        do iJ=T,1,-1   
			if (iJ==T) then
			

				call sub_bellman_bequest(V_b(:,:,:,:,:,iJ,iA1,iA2,:),c_b(:,:,:,:,:,iJ,iA1,iA2,:),BB_b(:,:,:,:,:,iJ,iA1,iA2,:),&
				HHi_b(:,:,:,:,:,iJ,iA1,iA2,:),LL_b(:,:,:,:,:,iJ,iA1,iA2,:),R,H,L,e(iJ),p,Y,chi,Trans,rental(:,:,iA1,iA2)) !-Jret+1
			!	print*, 'VFI loop: iA,thread,iJ:', iA, omp_get_thread_num(),iJ
								!print*, V_b(:,:,:,:,:,iP,iJ,iA1,iA2,2)		
			else
				if (type_sim<2) then
				print*, 'VFI loop: iA,thread,iJ:', iA, omp_get_thread_num(),iJ
				end if 
				!print*, iJ
				call sub_continuation_value(L,H,p,pprimein(:,:,:,:,iA1,iA2),ilo_pprime(:,:,:,:,:,iA1,iA2),theta_pprime(:,:,:,:,:,iA1,iA2),FYY,Fee(iJ,:,:),&
					V_b(:,:,:,:,:,iJ+1,iA1,iA2,:),W_b(:,:,:,:,:,iJ,iA1,iA2,:),q(:,:,:,:,:,iJ:iJ+1,iA1,iA2),&
					q_e(:,:,:,:,:,1:2,iA1,iA2),LL_b(:,:,:,:,:,iJ+1,iA1,iA2,2),iJ,iA1,chi,loge(iJ:iJ+1,:),Y,W_loc_own(:,:,:,:,:,iJ+1,iA1,iA2),W_loc_rent(:,:,:,iJ+1,iA1,iA2),W_own(:,:,:,:,:,iJ+1,iA1,iA2,:))
       call sub_bellman_bor(V_b(:,:,:,:,:,iJ,iA1,iA2,:),c_b(:,:,:,:,:,iJ,iA1,iA2,:),LL_b(:,:,:,:,:,iJ,iA1,iA2,:),LLi_b(:,:,:,:,:,iJ,iA1,iA2,:),&
					BB_b(:,:,:,:,:,iJ,iA1,iA2,:),HHi_b(:,:,:,:,:,iJ,iA1,iA2,:),W_b(:,:,:,:,:,iJ,iA1,iA2,:),Y,chi,loge(iJ,:),H,L,p,e(iJ),iJ,q(:,:,:,:,:,iJ,iA1,iA2),Trans,rental(:,:,iA1,iA2),R)
        q_e(:,:,:,:,:,2,iA1,iA2)=q_e(:,:,:,:,:,1,iA1,iA2)
				end if
                
            do iP=1,countP            
                do iExo=1,nY*nlamb*nE*nphi
                    do iD=1,nD+2
                        if ((iD==1 .and. iJ<T) .or. iD==4 .or. iD==5 .or. iD==6 .or. iD==7) then
                            bariL(iA1,iA2)=1
                            bariH(iA1,iA2)=1
						else if (iD==2 .or. iD==3 .or. (iD==1 .and. iJ==T)) then
                            bariL(iA1,iA2)=nL
							bariH(iA1,iA2)=nH
                        end if
						
												
						do iH=1,bariH(iA1,iA2)
							do iL=1,bariL(iA1,iA2) 
								if (iD<6) then
									do iB=1,nB

										if (iD<=3) then
											HH_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)=H(HHi_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD))
										else if (iD>3) then
											HH_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)=R(HHi_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD))
										end if  
										if (iD /= 2) then
											LL_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)=L(LLi_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD))
										end if 
										!if (iD==2 .and. iJ==30 .and. iP==4 .and. iExo>234) then
										!print*, iL,iH,iB,iExo,LL_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)
										!end if 
									end do 
								else 
									do iB=1,nBs
										if (iD==6) then
											HH_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)=H(HHi_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD))
											LL_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)=L(LLi_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD))
										else
											HH_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD)=R(HHi_b(iL,iH,iB,iExo,iP,iJ,iA1,iA2,iD))
											
										end if
									end do
								
								end if 
                            end do
							
                        end do
	
                    end do
                end do
            end do
		if (type_sim<2) then 
		print*, iJ, 'mean q', sum(q)/(nL*nB*nH*nY*nlamb*nE*countP*T*nA*nA*nphi)
		end if 
        end do
		!end if 
	print*, 'End:   iA: 1,2', iA, iA1,iA2, omp_get_thread_num(), omp_get_num_threads()	
    end do
    !$OMP END PARALLEL Do
	deallocate(W_loc_own)
	deallocate(W_loc_rent)	
    deallocate(W_b)
    deallocate(LLi_b)
    deallocate(HHi_b)
	deallocate(q_e)

        print*, 'Done with VFI'
    end subroutine sub_VFI
end module mod_VFI