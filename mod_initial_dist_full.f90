!       This file contains the module mod_initial_dist that creates the
!       Initial distribution for Krussel-Smith Siulation

module mod_initial_dist
    !USE IFPORT
    use parameters
    use mod_matlab
    use mod_interp
    implicit none 
contains
    subroutine sub_initial_dist(L,B,H,R,loge,Fee,age_sim_bor,ET_loc,V,maxH,minH,maxL,minL,maxB,minB,init_share)

    implicit none       
    integer :: sizeseed
    real(8), dimension(I) :: Idraw, Hdraw, Edraw, Bdraw
    integer(8), dimension(I) :: t_loc_bor
    real(8), dimension(:), intent(in) :: L,H,B,R
    real(8), dimension(:,:), intent(in) :: loge
    real(8), dimension(:,:,:), intent(in) :: Fee
    real(8) :: maxI, minI,H0,weights(T),maxHdraw,minHdraw,maxBdraw,minBdraw
    real(8), intent(in) :: maxH, minH,maxL,minL,maxB,minB,init_share
    !real(8), dimension(:), intent(out) :: L_sim_bor, H_sim_bor, B_sim_bor
    integer, dimension(:,:), intent(out) :: age_sim_bor, ET_loc
	!integer, dimension(:), intent(out) :: rent_buy
    integer :: iI,it,iJ, weights_int(T),iE,iEE, col2
    real(8), dimension(T,nE,nE+1) :: Eprob_CDF
    real(8), allocatable :: PTE(:)
	integer, allocatable :: seed(:), old(:)
	real(8) :: num, num2, count
    real(8), dimension(:), intent(in) :: V
	real(8), dimension(T,nE) :: edist
	integer(8), dimension(T) :: count_age,Eindex
	real(8), dimension(nE,2) :: initliqdist, initilldist, initnwdist, initadist
	real(8), dimension(nE) :: initliq, initnw, initill, initashare
	integer, dimension(nE) :: count_E

	!rent_buy=0
	age_sim_bor=0
	count_age=0
	count_E=0
	Eindex=1
	edist=0d0
	count=0
	call random_seed(size=sizeseed)
	allocate(seed(sizeseed))
	allocate(old(sizeseed))
	allocate(PTE(I))
	seed=123456789
	


  ET_loc=3

  	OPEN(3, FILE = './KMV_Input/zydist.txt',status='old',action='read')
		do iJ=1,Jret-1
				read(3,*) edist(iJ,:)
		end do
	close(3)
	
		

	
	do iJ=Jret,T
		edist(iJ,:)=edist(Jret-1,:)
	end do
	
	edist=int(edist*I/T)
	
	do iJ=1,T
	!print*, 'edisit', iJ, sum(edist(iJ,:)),edist(iJ,:)
		count=0
		do while (sum(edist(iJ,:))<(I/T)) 
			count=count+1
			edist(iJ,count)=edist(iJ,count)+1
		end do
		!print*, 'edisit', iJ, sum(edist(iJ,:))/I,(edist(iJ,:)/I)*30d0
	end do
	!print*, sum(edist)/I
	
	initliqdist=0.0d0
	initilldist=0d0
	initnwdist=0.0d0
	initnw=0.0d0
	OPEN(1, FILE = './KMV_Input/initwbroad.txt',status='old',action='read')
		DO iE=1,nE
			READ(1,*) initliqdist(iE,1),initliqdist(iE,2),initliq(iE),initilldist(iE,1),initilldist(iE,2),initill(iE),initnwdist(iE,1),initnwdist(iE,2),initnw(iE)
		END DO
	CLOSE(1)

	initliq=initliq/(DataAvAnnualEarns/(Numeraire))
	initill=initill/(DataAvAnnualEarns/(Numeraire))
	initnw=initnw/(DataAvAnnualEarns/(Numeraire))

	OPEN(1, FILE = './KMV_Input/initw.txt')
		DO iE=1,nE
			READ(1,*) initadist(iE,1),initadist(iE,2),initashare(iE)
      	initadist(iE,1)=1d0
			initadist(iE,2)=1.0d0-initadist(iE,1)
		END DO
	CLOSE(1)

	call random_seed(PUT=seed)

	do iJ=1,T
		
		age_sim_bor((iJ-1)*(I/T)+1:iJ*I/T,1)=iJ
		!print*, 'iJ',iJ,I, (iJ-1)*(I/T)+1,iJ*I/T, age_sim_bor(iJ*I/T*0.5d0,1)
	end do

    do iI=1,I
    call random_number(num)
    call random_number(num2)
	
		
		!age_sim_bor(iI,1)=int(num*(T+1-1))+1  
      
	  
	  	count_age(age_sim_bor(iI,1))=count_age(age_sim_bor(iI,1))+1
		
		! distributes income according to ages
		if (count_age(age_sim_bor(iI,1))<=sum(Edist(age_sim_bor(iI,1),1:Eindex(age_sim_bor(iI,1))))) then
			ET_loc(iI,1)=Eindex(age_sim_bor(iI,1)) !Edist(age_sim_bor(iI),iE)
		else
			if (Eindex(age_sim_bor(iI,1))<5) then
				Eindex(age_sim_bor(iI,1))=Eindex(age_sim_bor(iI,1))+1
			end if
			ET_loc(iI,1)=Eindex(age_sim_bor(iI,1))
		end if

			count_E(Eindex(age_sim_bor(iI,1)))=count_E(Eindex(age_sim_bor(iI,1)))+1
			
		
    end do

    
   ! print*, 'H max min', maxval(H_sim_bor), minval(H_sim_bor)

    
	do it=1,T_sim-1
		do iI=1,I
            if (age_sim_bor(iI,it)==T) then
                age_sim_bor(iI,it+1)=1        
            else
                age_sim_bor(iI,it+1)=age_sim_bor(iI,it)+1
            end if
        end do
    end do
	

    
    if (NE>1) then
        Eprob_CDF(:,:,1)=0d0
        do iJ=1,T
        do iE=1,nE
            do iEE=1,nE
                Eprob_CDF(iJ,iE,iEE+1)=Eprob_CDF(iJ,iE,iEE)+FEE(iJ,iE,iEE)
            end do
           ! print*, iJ,iE, Eprob_CDF(iJ,iE,:)
        end do
        end do
 
		
		
		print*, 'sum Edist', sum(Edist(1,:))
    !Simulating ET
        do it=1,T_sim-1
			count_age(1)=0
			Eindex=1
            do iI=1,I
				if (age_sim_bor(iI,it)==T) then
					count_age(1)=count_age(1)+1
					
					! distributes income according to ages
					if (count_age(1)<=sum(Edist(1,1:Eindex(1)))) then
						ET_loc(iI,it+1)=Eindex(1) !Edist(age_sim_bor(iI),iE)
					else
						if (Eindex(1)<5) then
							Eindex(1)=Eindex(1)+1
						end if
						ET_loc(iI,it+1)=Eindex(1)
					end if
					!ET(iI,it+1)=loge(1,ET_loc(iI,it+1))
					!print*, 'count_age',count_age(1),Eindex(1),sum(Edist(1,1:Eindex(1))) , ET_loc(iI,it+1)
                else if (age_sim_bor(iI,it)<Jret-2) then
        		        call random_number(num)
                    PTE(iI)=num
                    col2=bsearch(PTE(iI),Eprob_CDF(age_sim_bor(iI,it),ET_loc(iI,it),:))
                    ET_loc(iI,it+1)=col2
                    !ET(iI,it+1)=loge(age_sim_bor(iI,it+1),col2)
                else if (age_sim_bor(iI,it)>=Jret-2 .and. age_sim_bor(iI,it)<T) then
                    ET_loc(iI,it+1)=ET_loc(iI,it)
                    !ET(iI,it+1)=loge(Jret-1,ET_loc(iI,it))
                end if              
            end do

        end do
    else
        ET_loc=1
!        ET=loge(1,1)
    end if
!B_sim_bor=0d0
!L_sim_bor=0d0
!H_sim_bor=0d0
!print*, 'Sums', sum(ET(501,:)),sum(ET_loc(501,:))

    end subroutine sub_initial_dist
end module mod_initial_dist
