!       This file contains the module mod_initial_dist that creates the
!       Initial distribution for Krussel-Smith Siulation

module mod_initial_dist2
    !USE IFPORT
    use parameters
    use mod_matlab
    use mod_interp
    implicit none 
contains
    subroutine sub_initial_dist2(age_sim_bor,L_sim_bor,B_sim_bor,H_sim_bor,rent_buy,ET_loc,FEE,loge)
    implicit none       
    real(8), dimension(:), intent(out) :: L_sim_bor, H_sim_bor, B_sim_bor
    integer, dimension(:,:), intent(out) :: age_sim_bor, ET_loc
	integer, dimension(:) :: rent_buy
    real(8), dimension(:,:,:), intent(in) :: Fee
    real(8), dimension(:,:), intent(in) :: loge
    integer :: iI,it, col2, iE,iEE,iJ
    real(8) :: num
    real(8), dimension(T,nE,nE+1) :: Eprob_CDF
    real(8), dimension(I,T_sim) :: PTE
    integer :: sizeseed
    integer, allocatable :: seed(:), old(:)

    call random_seed(size=sizeseed)
	  allocate(seed(sizeseed))
	  allocate(old(sizeseed))
	  seed=123456789


    open(21,file='mu_Li_dist.csv',status='old',action='read')
    do iI=1,I
        read(21,*) L_sim_bor(iI)
    end do
    close(21)

    open (22, file='mu_Hi_dist.csv',status='old',action='read')
    do iI=1,I
        read (22,*) H_sim_bor(iI)
    end do
    close (22)
    
    open (23, file='mu_Bi_dist.csv',status='old',action='read')
    do iI=1,I
        read (23,*) B_sim_bor(iI)
    end do
    close (23)


    open (24, file='mu_age_dist.csv',status='old',action='read')
    do iI=1,I
        read (24,*) age_sim_bor(iI,1)
    end do
    close (24)
    
    open (25, file='mu_rent_buy_dist.csv',status='old',action='read')
    do iI=1,I
        read (25,*) rent_buy(iI)
    end do
    close (25)

	do it=1,T_sim-1
		do iI=1,I    
            if (age_sim_bor(iI,it)==T) then
                age_sim_bor(iI,it+1)=1        
            else
                age_sim_bor(iI,it+1)=age_sim_bor(iI,it)+1
            end if
        end do
    end do
    
    open (26, file='mu_Ei_dist.csv',status='old',action='read')
    do iI=1,I
        read (26,*) ET_loc(iI,1)
    end do
    close (26)
        
    !do iI=1,I    
    !ET(iI,1)=loge(age_sim_bor(iI,1),ET_loc(iI,1))
    !end do
    
    call random_seed(PUT=seed)
    if (NE>1) then
        Eprob_CDF(:,:,1)=0d0
        do iJ=1,T
        do iE=1,nE
            do iEE=1,nE
                Eprob_CDF(iJ,iE,iEE+1)=Eprob_CDF(iJ,iE,iEE)+FEE(iJ,iE,iEE)
            end do
            !print*, iJ,iE, Eprob_CDF(iJ,iE,:)
        end do
        end do
       
    !Simulating ET
        do it=1,T_sim-1
            do iI=1,I
                if (age_sim_bor(iI,it)<Jret-2) then
        		        call random_number(num)
                    PTE(iI,it+1)=num
                    col2=bsearch(PTE(iI,it+1),Eprob_CDF(age_sim_bor(iI,it),ET_loc(iI,it),:))
                    ET_loc(iI,it+1)=col2
                    !ET(iI,it+1)=loge(age_sim_bor(iI,it)+1,col2)
                else
                    ET_loc(iI,it+1)=ET_loc(iI,it)
                    !ET(iI,it+1)=loge(Jret-1,ET_loc(iI,it))
                end if              
            end do

        end do
    else
        ET_loc=1
        !ET=loge(1,1)
    end if





    end subroutine sub_initial_dist2
end module mod_initial_dist2
