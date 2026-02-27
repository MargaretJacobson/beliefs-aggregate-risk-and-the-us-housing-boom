! Questions 1&2 Problem Set 3

! By Margaret Jacobson, Indiana University marmjaco@indiana.edu
! For E 724: Computational Macroeconomics taught by Grey Gordon in
! Spring 2016 at Indiana University

! To write this code, I consulted lecture notes by Grey Gordon
! accessed via Grey Gordon's website
! https://sites.google.com/site/greygordon/teaching and Canvas. 

! This file contains the module mod_interp
! mod_interp contains the function bsearch which uses binary search
! to find the location of x in a sorted list. It also contains
! interp1q which uses binary search over a function   

module mod_interp
    implicit none
contains
    function bsearch (xhat ,x) result ( ilo )
        implicit none
        real (8) , dimension (:) , intent (in) :: x
        real (8) , intent (in) :: xhat
        integer :: ilo, a, b

        a=1
        b=size(x)
        if (xhat<x(1)) then
            ilo=1
        else if (xhat>=x(b)) then
            ilo=size(x)-1
        else 
            do while (b-a>1)
                ilo=(b+a)/2
            if (xhat<x(ilo)) then
                b=ilo
            else
                a=ilo
            end if
            end do
            if ((b-a)==1) then
                ilo=a               
            end if 
        end if
       
    end function bsearch



    subroutine search(ilo,xhat,x)
        implicit none
        real(8), dimension(1:), intent(in) :: x
        real(8), intent(in) :: xhat
        integer, intent(out) :: ilo
        
        if (size(x)<=20) then
            call lsearch(ilo,xhat,x)
        else
            ilo=bsearch(xhat,x)
        end if

!         if (xhat<x(1) .and. ilo/=1) STOP 'search: wrong case 1'
!         if (xhat>x(size(x)) .and. ilo/=size(x)-1) STOP 'search: wrong case 2'
!         if (size(x)==1 .and. ilo/=1) STOP 'search: wrong case 3'
!         if (size(x)>=2 .and. xhat>=x(1) .and. xhat<=x(size(x)) .and. (xhat<x(ilo) .or. xhat>x(ilo+1))) STOP 'search: wrong case 4'

    end subroutine


 ! Use sequential/linear search to find the lowerbound on a bracketing interval
    ! This has far worse theoretical performance than binary search O(n) compared to 
    ! O(log2(n)), but it is much easier for the compiler to understand and use effectively.
    ! Consequently for small arrays, it is faster.
    subroutine lsearch(ilo,xhat,x)
        implicit none
        real(8), dimension(1:), intent(in) :: x
        real(8), intent(in) :: xhat
        integer, intent(out) :: ilo
        ! local
!         integer :: ihi

!    do ihi = 2,size(x)
!        ! If here, then size(x)-1>=1 <=> size(x)>=2 (if size(x)==1, we just return ilo = 1)
! 
!        ! Three cases to be concerned with
!        ! (1) xhat<x(1)
!        ! (2) xhat>=x(n)
!        ! (3) xhat>=x(1) .and. xhat<x(n)
! 
!        ! Test whether xhat<x(ilo). 
!        ! In case (1), this will always evaluate to true. So we should exit immediately returning ilo=ilo-1.
!        ! In case (2), this will never evaluate to true. 
!        ! In case (3), this will evaluate to true when xhat<x(ilo) and there was no j<ilo such that xhat<x(j).
!        !   Consequently, we have xhat>=x(j) for all j<ilo, and in particular, xhat>=x(ilo-1). This means, xhat \in [x(ilo-1),x(ilo))
!        !   So, return ilo-1, just as we did in case (1).
!        if (xhat<x(ihi)) then
!            ilo = ihi-1
!            return
!        end if
!    end do
!    ! If still here, we have case (2) or the size of x is 1. So, we return ilo = n-1 if size>1 and 1 o/w
!    ilo = max(size(x)-1,1)
        

        do ilo = 1,size(x)-1
            ! If here, then size(x)-1>=1 <=> size(x)>=2 (if size(x)==1, we just return ilo = 1)

            ! Three cases to be concerned with
            ! (1) xhat<x(1)
            ! (2) xhat>=x(n)
            ! (3) xhat>=x(1) .and. xhat<x(n)

            ! Test whether xhat<x(ilo+1). 
            ! In case (1), this will always evaluate to true. So we should exit immediately returning ilo=1
            ! In case (2), this will never evaluate to true. 
            ! In case (3), this will evaluate to true when xhat<x(ilo+1) and there was no j<ilo+1 such that xhat<x(j).
            !   Consequently, we have xhat>=x(j) for all j<ilo+1, and in particular, xhat>=x(ilo). This means,
            !   xhat \in [x(ilo),x(ilo+1)). So, return ilo, just as we did in case (1).
            if (xhat<x(ilo+1)) return
        end do
        ! If still here, we have case (2) or the size of x is 1. So, we return ilo = n-1 if size>1 and 1 o/w
        ilo = max(size(x)-1,1)

    end subroutine

    function bsearch_sim (xhat ,x) result ( ilo )
        implicit none
        real (8) , dimension (:) , intent (in) :: x
        real (8) , intent (in) :: xhat
        integer :: ilo, a, b

        a=1
        b=size(x)
        if (xhat<x(1)) then
            ilo=1
        else if (xhat>=x(b)) then
            ilo=size(x)
        else 
            do while (b-a>1)
                ilo=(b+a)/2
            if (xhat<x(ilo)) then
                b=ilo
            else
                a=ilo
            end if
            end do
            if ((b-a)==1) then
                ilo=a               
            end if 
        end if
       
    end function bsearch_sim
	


    function interpnearest (xi,x) result (yi)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) , xi
        real (8) :: yi
        real (8) :: theta
        integer :: ilo        
        ilo=bsearch_sim(xi,x)
        theta=(xi-x(ilo))/(x(ilo+1)-x(ilo)) 
		!print*, xi, ilo, theta
        if (theta<=0.5d0) then
            yi=x(ilo)
        else
            yi=x(ilo+1)
        end if
        
    end function interpnearest

    function interp1q (x,y,xi) result (yi)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) ,y(:) , xi
        real (8) :: yi
        real (8) :: theta
        integer :: ilo        
        ilo=bsearch(xi,x)
        theta=1-(xi-x(ilo))/(x(ilo+1)-x(ilo)) 
        yi=y(ilo)*theta+y(ilo+1)*(1-theta)
        
    end function interp1q
	
	function interp1q16 (x,y,xi) result (yi)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) , xi
		real (16), intent(in) ::  y(:)
        real (16) :: yi
        real (8) :: theta
        integer :: ilo        
        ilo=bsearch(xi,x)
        theta=1-(xi-x(ilo))/(x(ilo+1)-x(ilo)) 
        yi=y(ilo)*theta+y(ilo+1)*(1-theta)
        
    end function interp1q16
	
	function interp1qnoextrap (x,y,xi) result (yi)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) ,y(:) , xi
        real (8) :: yi
        real (8) :: theta,thetat
        integer :: ilo  
		
        ilo=bsearch(xi,x)
        theta=1d0-(xi-x(ilo))/(x(ilo+1)-x(ilo)) 
		thetat=max(0d0, min(1d0,theta)) 
        yi=y(ilo)*thetat+y(ilo+1)*(1d0-thetat)
        
    end function interp1qnoextrap
	
		function interp1qnoextrap_print (x,y,xi) result (yi)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) ,y(:) , xi
        real (8) :: yi
        real (8) :: theta
        integer :: ilo  
		
        ilo=bsearch(xi,x)
        theta=1d0-(xi-x(ilo))/(x(ilo+1)-x(ilo)) 
		theta=max(0d0, min(1d0,theta)) 
		!minval((/theta,0d0/))
		!theta=maxval((/theta,1d0/))
		print*,xi, ilo,theta,ilo+1,1d0-theta
		print*, size(y)
		print*,	y
		print*, y(ilo)
		print*, y(ilo+1)
		print*, y(ilo)*theta
		print*, y(ilo+1)*(1d0-theta)
        yi=y(ilo)*theta+y(ilo+1)*(1d0-theta)
		
        
    end function interp1qnoextrap_print


    subroutine base_fun_noextrap(x,xi,ilo,theta)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) , xi  
	real (8), intent(out) :: theta(2)
        integer, intent(out) :: ilo(2)        
        ilo(1)=bsearch(xi,x)
	ilo(2)=ilo(1)+1
        theta(1)=1d0-(xi-x(ilo(1)))/(x(ilo(2))-x(ilo(1))) 
	theta(1)=max(0d0, min(1d0,theta(1))) 	
	theta(2)=1d0-theta(1)

        
    end subroutine base_fun_noextrap
	
	
	subroutine base_fun(x,xi,ilo,theta)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) , xi  
	real (8), intent(out) :: theta(2)
        integer, intent(out) :: ilo(2)        
        ilo(1)=bsearch(xi,x)
	ilo(2)=ilo(1)+1
        theta(1)=1-(xi-x(ilo(1)))/(x(ilo(2))-x(ilo(1))) 
	theta(2)=1-theta(1)

        
    end subroutine base_fun


    function interp1q_ilo (x,y,xi,ilo) result (yi)
        ! use any modules here
        implicit none
        real (8) , intent (in) :: x(:) ,y(:) , xi
        real (8) :: yi
        real (8) :: theta
        integer , intent (in) :: ilo        

        theta=1-(xi-x(ilo))/(x(ilo+1)-x(ilo)) 
        yi=y(ilo)*theta+y(ilo+1)*(1-theta)
        
    end function interp1q_ilo


!Code below from Grey Gordon, modified by M. Jacobson
function interp1q_scal_ymat(x,y,xi) result(yi)
    !use mod_interp, only: search
    implicit none
    real(8), dimension(:), intent(in) :: x
    real(8), dimension(:,:), intent(in) :: y
    real(8), intent(in) :: xi
    real(8), dimension(size(y,2)) :: yi
    ! local 
    integer :: nx, ny
    integer :: ilo 

    nx = size(x)
    ny = size(y,1)
    if (nx/=ny) STOP 'interp1_scal: dims dont agree'
    if (nx==1) then
        yi = y(1,:) 
        return
    end if

    ! Find the location of xi in the x grid assuming it is sorted
    call search(ilo,xi,x)

    ! Do the interpolation
    yi(:) = y(ilo,:) + (y(ilo+1,:)-y(ilo,:))*(xi-x(ilo))/(x(ilo+1)-x(ilo))

end function interp1q_scal_ymat


function interp1q_vec(x,y,xi) result(yi)
    !use mod_interp, only: search
    implicit none
    real(8), dimension(:), intent(in) :: x
    real(8), dimension(:), intent(in) :: y
    real(8), dimension(:), intent(in) :: xi
    real(8), dimension(size(xi)):: yi
    ! local 
    integer, allocatable :: ilo(:)
    real(8) :: mlo
    integer :: i

    allocate(ilo(size(xi)))

    call bsearch_many(ilo,xi,x)

    do i = 1,size(xi)
        mlo = (x(ilo(i)+1)-xi(i))/(x(ilo(i)+1)-x(ilo(i)))
        yi(i) = y(ilo(i))*mlo + y(ilo(i)+1)*(1d0-mlo)
    end do

!     do i = 1,size(xi)
!         yi(i) = interp1q_scal(x,y,xi(i))
!     end do

end function interp1q_vec

! Also from Grey Gordon with modifications
function interp1q_vec_ymat(x,y,xi) result(yi)
  !  use mod_interp, only: search
    implicit none
    real(8), dimension(:), intent(in) :: x
    real(8), dimension(:,:), intent(in) :: y
    real(8), dimension(:), intent(in) :: xi
    real(8), dimension(size(xi),size(y,2)):: yi
    ! local 
    integer, allocatable :: ilo(:)
    real(8) :: mlo
    integer :: i

!     do i = 1,size(xi)
!         yi(i,:) = interp1q_scal_ymat(x,y,xi(i))
!     end do

    allocate(ilo(size(xi)))

    call bsearch_many(ilo,xi,x)

    do i = 1,size(xi)
        mlo = (x(ilo(i)+1)-xi(i))/(x(ilo(i)+1)-x(ilo(i)))
        yi(i,:) = y(ilo(i),:)*mlo + y(ilo(i)+1,:)*(1d0-mlo)
    end do


end function interp1q_vec_ymat

function interp1q_vec_vec(x,y,xi) result(yi)
  !  use mod_interp, only: search
    implicit none
    real(8), dimension(:), intent(in) :: x
    real(8), dimension(:,:), intent(in) :: y
    real(8), dimension(:), intent(in) :: xi
    real(8), dimension(size(y,2)):: yi
    ! local 
    integer, allocatable :: ilo(:)
    real(8) :: mlo
    integer :: i

!     do i = 1,size(xi)
!         yi(i,:) = interp1q_scal_ymat(x,y,xi(i))
!     end do

    allocate(ilo(size(xi)))

    call bsearch_many(ilo,xi,x)

    do i = 1,size(xi)
        mlo = (x(ilo(i)+1)-xi(i))/(x(ilo(i)+1)-x(ilo(i)))
        yi(i) = y(ilo(i),i)*mlo + y(ilo(i)+1,i)*(1d0-mlo)
    end do


end function interp1q_vec_vec


! From Grey Gordon
! Given sorted (ascending) lists xhat and x, finds ilo as if call bsearch(ilo(i),xhat(i),x) was called, but does so faster
subroutine bsearch_many(ilo,xhat,x)
    implicit none
    integer, intent(out) :: ilo(:)
    real(8), intent(in) :: xhat(:),x(:)
    ! local
    integer :: n,nx

    n = size(ilo) 
    nx = size(x)

    if (n>=2) then
        call bsearch_core(ilo(1),xhat(1),x,1,nx)
        call bsearch_core(ilo(n),xhat(n),x,ilo(1),nx)

        if (n==2) return
        call bsearch_many_core(ilo,xhat,x,1,n)

    elseif (n==1) then
        ilo(1)=bsearch(xhat(1),x)
    elseif (n<=0) then
        stop 'ERROR: bsearch_many: size(ilo)<=0'
    end if


end subroutine


subroutine bsearch_core(ilo,xhat,x,search_lb,search_ub)
    implicit none
    real(8), intent(in) :: x(*),xhat
    integer, intent(in) :: search_lb,search_ub
    integer, intent(out) :: ilo
    !local
    integer :: ihi,mid
        
    ilo = search_lb
    ihi = search_ub
    mid = (ihi+ilo)/2 ! mid will have ilo<=mid<ihi as long as ihi>ilo. If ihi=ilo, then mid=ilo, so the 
                      ! while loop is not entered and ilo=1 is returned 
    
    do while (mid>ilo) ! mid=ilo means ihi=ilo+1.

        ! Here, we must have ihi>ilo and ilo<=mid<ihi

        ! Bisect. 
        ! If xhat<x(1), then xhat>=x(mid) will never occur and ilo will never be modified.
        ! If xhat>=x(nx), then xhat>=x(mid) will always occur (note mid<ihi). Consequently, 
        !      ilo = mid occurs until mid==ilo, which only occures if ihi==ilo+1.
        if (xhat>=x(mid)) then
            ilo = mid
        else
            ihi = mid
        end if

        mid = (ihi+ilo)/2 

    end do

end subroutine 


! Assumes ilo(lb) and ilo(ub) are known, exploits monotonicity
recursive subroutine bsearch_many_core(ilo,xhat,x,lb,ub)
    implicit none
    integer, intent(inout) :: ilo(*)
    real(8), intent(in) :: xhat(*),x(*)
    integer, intent(in) :: lb,ub
    ! local
    integer :: m

    m = (lb+ub)/2

    call bsearch_core(ilo(m),xhat(m),x,ilo(lb),ilo(ub)+1) ! need to add 1 to the upper bound because ilo always returns the lower part of the bracketing interval

    if (m>lb+1) call bsearch_many_core(ilo,xhat,x,lb,m)
    if (ub>m+1) call bsearch_many_core(ilo,xhat,x,m,ub)


end subroutine 

    function bsearch2(xhat,yhat,x,y) result(ilo)
        implicit none
        real(8), intent (in) :: x(:), y(:), xhat, yhat
        integer :: ilo(2)
        if (xhat<x(1) .and. yhat<y(1)) then
            ilo(1)=1
            ilo(2)=1
        else if (xhat>=x(size(x)) .and. yhat>=y(size(y))) then
            ilo(1)=size(x)-1
            ilo(2)=size(y)-1
        else
            ilo(1)=size(x)
            ilo(2)=size(y)
            do while (ilo(1)>1 .and. ilo(2)>1)
                print*, ilo
                if (xhat<x(ilo(1))) then
                    ilo(1)=ilo(1)-1
                end if
                if (yhat<y(ilo(2))) then
                    ilo(2)=y(ilo(2))-1
                end if
            end do
        end if
    end function
    
    function interp2q(x,y,z,xi,yi) result(zi)
        implicit none
        real (8), intent (in) :: x(:), y(:), z(:,:), xi , yi
        real (8) :: zi
        real (8) :: thetax,thetay, fy(2) 
        integer :: ilo(2) 
        ilo(1)=bsearch(xi,x) 
        ilo(2)=bsearch(yi,y) 
        thetax=1d0-(real(xi-x(ilo(1))))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))   
    
        fy(1)=thetax*z(ilo(1),ilo(2))+(1-thetax)*z(ilo(1)+1,ilo(2))
        fy(2)=thetax*z(ilo(1),ilo(2)+1)+(1-thetax)*z(ilo(1)+1,ilo(2)+1)

        zi=thetay*fy(1)+(1-thetay)*fy(2)
    end function interp2q
	
	function interp2qnoextrap(x,y,z,xi,yi) result(zi)
        implicit none
        real (8), intent (in) :: x(:), y(:), z(:,:), xi , yi
        real (8) :: zi
        real (8) :: thetax,thetay, fy(2) 
        integer :: ilo(2) 
        ilo(1)=bsearch(xi,x) 
        ilo(2)=bsearch(yi,y) 
        thetax=1d0-(real(xi-x(ilo(1))))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))   
		thetax=max(0d0, min(1d0,thetax)) 
		thetay=max(0d0, min(1d0,thetay)) 
		!if (ilo(1)==size(x)) ilo(1)=ilo(1)-1 
		!if (ilo(2)==size(y)) ilo(2)=ilo(2)-1  
        fy(1)=thetax*z(ilo(1),ilo(2))+(1-thetax)*z(ilo(1)+1,ilo(2))
        fy(2)=thetax*z(ilo(1),ilo(2)+1)+(1-thetax)*z(ilo(1)+1,ilo(2)+1)

        zi=thetay*fy(1)+(1-thetay)*fy(2)
    end function interp2qnoextrap

	

    function interp2q_ilo(x,y,z,xi,yi,ilo) result(zi)
        implicit none
        real (8), intent (in) :: x(:), y(:), z(:,:), xi , yi
        real (8) :: zi
        real (8) :: thetax,thetay, fy(2) 
        integer, intent (in) :: ilo(:) 
        !ilo(1)=bsearch(xi,x) 
        !ilo(2)=bsearch(yi,y) 
        thetax=1d0-(real(xi-x(ilo(1))))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))   
    
        fy(1)=thetax*z(ilo(1),ilo(2))+(1-thetax)*z(ilo(1)+1,ilo(2))
        fy(2)=thetax*z(ilo(1),ilo(2)+1)+(1-thetax)*z(ilo(1)+1,ilo(2)+1)

        zi=thetay*fy(1)+(1-thetay)*fy(2)
    end function interp2q_ilo


    function interp3q(x,y,z,f,xi,yi,zi) result(fi)
        implicit none
        real (8), intent (in) :: x(:), y(:), z(:), xi , yi,zi,f(:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz, fy(4),fz(2)
        integer :: ilo(3)

		
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))

        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3))
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1)
        fy(3)=thetax*f(ilo(1),ilo(2)+1,ilo(3))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3))
        fy(4)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1)
        
        fz(1)=fy(1)*thetay+fy(3)*(1-thetay)
        fz(2)=fy(2)*thetay+fy(4)*(1-thetay)

        fi=fz(1)*thetaz+fz(2)*(1-thetaz)


    end function interp3q
	
	
    function interp3qnoextrap(x,y,z,f,xi,yi,zi) result(fi)
        implicit none
        real (8), intent (in) :: x(:), y(:), z(:), xi , yi,zi,f(:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz, fy(4),fz(2),thetaxt,thetayt,thetazt
        integer :: ilo(3)

		
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
		thetaxt=max(0d0, min(1d0,thetax)) 
		thetayt=max(0d0, min(1d0,thetay)) 
		thetazt=max(0d0, min(1d0,thetaz))
		
		!if (ilo(1)==size(x)) ilo(1)=ilo(1)-1 
		!if (ilo(2)==size(y)) ilo(2)=ilo(2)-1  
		!if (ilo(3)==size(z)) ilo(3)=ilo(3)-1 
		if(ilo(3)>size(z)) print*, ilo(3)
        fy(1)=thetaxt*f(ilo(1),ilo(2),ilo(3))+(1-thetaxt)*f(ilo(1)+1,ilo(2),ilo(3))
		fy(2)=thetaxt*f(ilo(1),ilo(2),ilo(3)+1)+(1-thetaxt)*f(ilo(1)+1,ilo(2),ilo(3)+1)
        fy(3)=thetaxt*f(ilo(1),ilo(2)+1,ilo(3))+(1-thetaxt)*f(ilo(1)+1,ilo(2)+1,ilo(3))
        fy(4)=thetaxt*f(ilo(1),ilo(2)+1,ilo(3)+1)+(1-thetaxt)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1)
        
        fz(1)=fy(1)*thetayt+fy(3)*(1-thetayt)
        fz(2)=fy(2)*thetayt+fy(4)*(1-thetayt)

        fi=fz(1)*thetazt+fz(2)*(1-thetazt)


    end function interp3qnoextrap
    
    function interp4q(x,y,z,a,f,xi,yi,zi,ai) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),xi,yi,zi,ai,f(:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa, fy(8),fz(4),fa(2)
        integer :: ilo(4)
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)

        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))

        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4))
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1)
        fy(3)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4))
        fy(4)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1)
        
        fy(5)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4))
        fy(6)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1)
        fy(7)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4))
        fy(8)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1)
           

        fz(1)=thetay*fy(1)+(1-thetay)*fy(5)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(6)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(7)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(8)
        
        fa(1)=fz(1)*thetaz+fz(3)*(1-thetaz)
        fa(2)=fz(2)*thetaz+fz(4)*(1-thetaz)


        fi=fa(1)*thetaa+fa(2)*(1-thetaa)
        
    end function interp4q
	function interp4qnoextrap(x,y,z,a,f,xi,yi,zi,ai) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),xi,yi,zi,ai,f(:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa, fy(8),fz(4),fa(2)
        integer :: ilo(4)
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)

        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
		thetax=max(0d0, min(1d0,thetax)) 
		thetay=max(0d0, min(1d0,thetay)) 
		thetaz=max(0d0, min(1d0,thetaz))
		thetaa=max(0d0, min(1d0,thetaa))
		
		!if (ilo(1)==size(x)) ilo(1)=ilo(1)-1 
		!if (ilo(2)==size(y)) ilo(2)=ilo(2)-1  
		!if (ilo(3)==size(z)) ilo(3)=ilo(3)-1 
		!if (ilo(4)==size(a)) ilo(4)=ilo(4)-1 
		
        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4))
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1)
        fy(3)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4))
        fy(4)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1)
        
        fy(5)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4))
        fy(6)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1)
        fy(7)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4))
        fy(8)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1)
           

        fz(1)=thetay*fy(1)+(1-thetay)*fy(5)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(6)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(7)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(8)
        
        fa(1)=fz(1)*thetaz+fz(3)*(1-thetaz)
        fa(2)=fz(2)*thetaz+fz(4)*(1-thetaz)


        fi=fa(1)*thetaa+fa(2)*(1-thetaa)
        
    end function interp4qnoextrap

    function interp5q(x,y,z,a,b,f,xi,yi,zi,ai,bi) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),b(:),xi,yi,zi,ai,bi,f(:,:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa, thetab, fy(16),fz(8),fa(4),fb(2)
        integer :: ilo(5)
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)
        ilo(5)=bsearch(bi,b)
        
        
        
        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
        thetab=1d0-(bi-b(ilo(5)))/(b(ilo(5)+1)-b(ilo(5)))
		

	
        

        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5))        
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1)
        fy(3)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5))
        fy(4)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1)
        
        fy(5)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5))
        fy(6)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1)
        fy(7)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5))
        fy(8)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1)
        
        fy(9)= thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5))
        fy(10)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1)
        fy(11)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5))
        fy(12)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1)
        
        fy(13)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5))
        fy(14)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1)
        fy(15)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5))
        fy(16)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1)


        
        fz(1)=thetay*fy(1)+(1-thetay)*fy(9)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(10)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(11)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(12)
        fz(5)=thetay*fy(5)+(1-thetay)*fy(13)
        fz(6)=thetay*fy(6)+(1-thetay)*fy(14)
        fz(7)=thetay*fy(7)+(1-thetay)*fy(15)
        fz(8)=thetay*fy(8)+(1-thetay)*fy(16)


        fa(1)=thetaz*fz(1)+(1-thetaz)*fz(5)
        fa(2)=thetaz*fz(2)+(1-thetaz)*fz(6)
        fa(3)=thetaz*fz(3)+(1-thetaz)*fz(7)
        fa(4)=thetaz*fz(4)+(1-thetaz)*fz(8)
        
        fb(1)=fa(1)*thetaa+fa(3)*(1-thetaa)
        fb(2)=fa(2)*thetaa+fa(4)*(1-thetaa)
        
        fi=fb(1)*thetab+fb(2)*(1-thetab)

        
     end function interp5q    

function interp5qnoextrap(x,y,z,a,b,f,xi,yi,zi,ai,bi) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),b(:),xi,yi,zi,ai,bi,f(:,:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa, thetab, fy(16),fz(8),fa(4),fb(2)
        integer :: ilo(5)
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)
        ilo(5)=bsearch(bi,b)
        
        
        
        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
        thetab=1d0-(bi-b(ilo(5)))/(b(ilo(5)+1)-b(ilo(5)))
        
		thetax=max(0d0, min(1d0,thetax)) 
		thetay=max(0d0, min(1d0,thetay)) 
		thetaz=max(0d0, min(1d0,thetaz))
		thetaa=max(0d0, min(1d0,thetaa))
		thetab=max(0d0, min(1d0,thetab))
		
		!if (ilo(1)==size(x)) ilo(1)=ilo(1)-1 
		!if (ilo(2)==size(y)) ilo(2)=ilo(2)-1  
		!if (ilo(3)==size(z)) ilo(3)=ilo(3)-1 
		!if (ilo(4)==size(a)) ilo(4)=ilo(4)-1 
		!if (ilo(5)==size(b)) ilo(5)=ilo(5)-1 
        
		
		
        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5))        
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1)
        fy(3)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5))
        fy(4)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1)
        
        fy(5)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5))
        fy(6)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1)
        fy(7)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5))
        fy(8)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1)
        
        fy(9)= thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5))
        fy(10)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1)
        fy(11)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5))
        fy(12)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1)
        
        fy(13)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5))
        fy(14)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1)
        fy(15)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5))
        fy(16)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1)


        
        fz(1)=thetay*fy(1)+(1-thetay)*fy(9)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(10)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(11)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(12)
        fz(5)=thetay*fy(5)+(1-thetay)*fy(13)
        fz(6)=thetay*fy(6)+(1-thetay)*fy(14)
        fz(7)=thetay*fy(7)+(1-thetay)*fy(15)
        fz(8)=thetay*fy(8)+(1-thetay)*fy(16)


        fa(1)=thetaz*fz(1)+(1-thetaz)*fz(5)
        fa(2)=thetaz*fz(2)+(1-thetaz)*fz(6)
        fa(3)=thetaz*fz(3)+(1-thetaz)*fz(7)
        fa(4)=thetaz*fz(4)+(1-thetaz)*fz(8)
        
        fb(1)=fa(1)*thetaa+fa(3)*(1-thetaa)
        fb(2)=fa(2)*thetaa+fa(4)*(1-thetaa)
        
        fi=fb(1)*thetab+fb(2)*(1-thetab)

        
     end function interp5qnoextrap    
          
     
    function interp6q(x,y,z,a,b,c,f,xi,yi,zi,ai,bi,ci) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),b(:),c(:),xi,yi,zi,ai,bi,ci,f(:,:,:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa,thetab,thetac,fy(32),fz(16),fa(8),fb(4),fc(2)
        integer :: ilo(6),iI
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)
        ilo(5)=bsearch(bi,b)
        ilo(6)=bsearch(ci,c)

        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
        thetab=1d0-(bi-b(ilo(5)))/(b(ilo(5)+1)-b(ilo(5)))
        thetac=1d0-(ci-c(ilo(6)))/(c(ilo(6)+1)-c(ilo(6)))
		


        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6))
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1)
        fy(3)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6))
        fy(4)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)
        fy(5)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6))
        fy(6)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)
        fy(7)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))
        fy(8)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        fy(9)= thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6))
        fy(10)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)
        fy(11)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))
        fy(12)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)
        fy(13)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))
        fy(14)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)
        fy(15)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))
        fy(16)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        fy(17)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6))
        fy(18)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1)
        fy(19)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6))
        fy(20)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)
        fy(21)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6))
        fy(22)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)
        fy(23)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))
        fy(24)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        fy(25)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6))+ (1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6))
        fy(26)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)
        fy(27)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))
        fy(28)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)
        fy(29)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))
        fy(30)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)
        fy(31)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))
        fy(32)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        
        !do iI=1,size(fz)
        !    fz(iI)=thetay*fy(iI)+thetay*fy(iI+size(fz))
        !end do
        fz(1)=thetay*fy(1)+(1-thetay)*fy(17)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(18)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(19)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(20)
        fz(5)=thetay*fy(5)+(1-thetay)*fy(21)
        fz(6)=thetay*fy(6)+(1-thetay)*fy(22)
        fz(7)=thetay*fy(7)+(1-thetay)*fy(23)
        fz(8)=thetay*fy(8)+(1-thetay)*fy(24)
        fz(9)=thetay*fy(9)+(1-thetay)*fy(25)
        fz(10)=thetay*fy(10)+(1-thetay)*fy(26)
        fz(11)=thetay*fy(11)+(1-thetay)*fy(27)
        fz(12)=thetay*fy(12)+(1-thetay)*fy(28)  
        fz(13)=thetay*fy(13)+(1-thetay)*fy(29)
        fz(14)=thetay*fy(14)+(1-thetay)*fy(30)
        fz(15)=thetay*fy(15)+(1-thetay)*fy(31)    
        fz(16)=thetay*fy(16)+(1-thetay)*fy(32) 
        
        
        fa(1)=thetaz*fz(1)+(1-thetaz)*fz(9)
        fa(2)=thetaz*fz(2)+(1-thetaz)*fz(10)
        fa(3)=thetaz*fz(3)+(1-thetaz)*fz(11)
        fa(4)=thetaz*fz(4)+(1-thetaz)*fz(12)
        fa(5)=thetaz*fz(5)+(1-thetaz)*fz(13)
        fa(6)=thetaz*fz(6)+(1-thetaz)*fz(14)
        fa(7)=thetaz*fz(7)+(1-thetaz)*fz(15)
        fa(8)=thetaz*fz(8)+(1-thetaz)*fz(16)


        fb(1)=thetaa*fa(1)+(1-thetaa)*fa(5)
        fb(2)=thetaa*fa(2)+(1-thetaa)*fa(6)
        fb(3)=thetaa*fa(3)+(1-thetaa)*fa(7)
        fb(4)=thetaa*fa(4)+(1-thetaa)*fa(8)
        
        fc(1)=fb(1)*thetab+fb(3)*(1-thetab)
        fc(2)=fb(2)*thetab+fb(4)*(1-thetab)

        fi=fc(1)*thetac+fc(2)*(1-thetac)
        
    end function interp6q
    
    function interp6qnoextrap(x,y,z,a,b,c,f,xi,yi,zi,ai,bi,ci) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),b(:),c(:),xi,yi,zi,ai,bi,ci,f(:,:,:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa,thetab,thetac,fy(32),fz(16),fa(8),fb(4),fc(2)
        integer :: ilo(6),iI
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)
        ilo(5)=bsearch(bi,b)
        ilo(6)=bsearch(ci,c)

        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
        thetab=1d0-(bi-b(ilo(5)))/(b(ilo(5)+1)-b(ilo(5)))
        thetac=1d0-(ci-c(ilo(6)))/(c(ilo(6)+1)-c(ilo(6)))
		
		thetax=max(0d0, min(1d0,thetax)) 
		thetay=max(0d0, min(1d0,thetay)) 
		thetaz=max(0d0, min(1d0,thetaz))
		thetaa=max(0d0, min(1d0,thetaa))
		thetab=max(0d0, min(1d0,thetab))
		thetac=max(0d0, min(1d0,thetac))
		
		!if (ilo(1)==size(x)) ilo(1)=ilo(1)-1 
		!if (ilo(2)==size(y)) ilo(2)=ilo(2)-1  
		!if (ilo(3)==size(z)) ilo(3)=ilo(3)-1 
		!if (ilo(4)==size(a)) ilo(4)=ilo(4)-1 
		!if (ilo(5)==size(b)) ilo(5)=ilo(5)-1 
		!if (ilo(6)==size(c)) ilo(6)=ilo(6)-1 

        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6))
        fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1)
        fy(3)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6))
        fy(4)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)
        fy(5)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6))
        fy(6)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)
        fy(7)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))
        fy(8)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        fy(9)= thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6))
        fy(10)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)
        fy(11)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))
        fy(12)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)
        fy(13)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))
        fy(14)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)
        fy(15)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))
        fy(16)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        fy(17)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6))
        fy(18)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1)
        fy(19)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6))
        fy(20)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1)
        fy(21)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6))
        fy(22)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1)
        fy(23)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6))
        fy(24)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        fy(25)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6))+ (1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6))
        fy(26)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1)
        fy(27)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6))
        fy(28)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1)
        fy(29)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6))
        fy(30)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1)
        fy(31)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6))
        fy(32)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1)
        
        
        !do iI=1,size(fz)
        !    fz(iI)=thetay*fy(iI)+thetay*fy(iI+size(fz))
        !end do
        fz(1)=thetay*fy(1)+(1-thetay)*fy(17)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(18)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(19)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(20)
        fz(5)=thetay*fy(5)+(1-thetay)*fy(21)
        fz(6)=thetay*fy(6)+(1-thetay)*fy(22)
        fz(7)=thetay*fy(7)+(1-thetay)*fy(23)
        fz(8)=thetay*fy(8)+(1-thetay)*fy(24)
        fz(9)=thetay*fy(9)+(1-thetay)*fy(25)
        fz(10)=thetay*fy(10)+(1-thetay)*fy(26)
        fz(11)=thetay*fy(11)+(1-thetay)*fy(27)
        fz(12)=thetay*fy(12)+(1-thetay)*fy(28)  
        fz(13)=thetay*fy(13)+(1-thetay)*fy(29)
        fz(14)=thetay*fy(14)+(1-thetay)*fy(30)
        fz(15)=thetay*fy(15)+(1-thetay)*fy(31)    
        fz(16)=thetay*fy(16)+(1-thetay)*fy(32) 
        
        
        fa(1)=thetaz*fz(1)+(1-thetaz)*fz(9)
        fa(2)=thetaz*fz(2)+(1-thetaz)*fz(10)
        fa(3)=thetaz*fz(3)+(1-thetaz)*fz(11)
        fa(4)=thetaz*fz(4)+(1-thetaz)*fz(12)
        fa(5)=thetaz*fz(5)+(1-thetaz)*fz(13)
        fa(6)=thetaz*fz(6)+(1-thetaz)*fz(14)
        fa(7)=thetaz*fz(7)+(1-thetaz)*fz(15)
        fa(8)=thetaz*fz(8)+(1-thetaz)*fz(16)


        fb(1)=thetaa*fa(1)+(1-thetaa)*fa(5)
        fb(2)=thetaa*fa(2)+(1-thetaa)*fa(6)
        fb(3)=thetaa*fa(3)+(1-thetaa)*fa(7)
        fb(4)=thetaa*fa(4)+(1-thetaa)*fa(8)
        
        fc(1)=fb(1)*thetab+fb(3)*(1-thetab)
        fc(2)=fb(2)*thetab+fb(4)*(1-thetab)

        fi=fc(1)*thetac+fc(2)*(1-thetac)
        
    end function interp6qnoextrap
	
    function interp7q(x,y,z,a,b,c,d,f,xi,yi,zi,ai,bi,ci,di) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),b(:),c(:),d(:),xi,yi,zi,ai,bi,ci,di,f(:,:,:,:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa,thetab,thetac,thetad,fy(64),fz(32),fa(16),fb(8),fc(4),fd(2)
        integer :: ilo(7),iI
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)
        ilo(5)=bsearch(bi,b)
        ilo(6)=bsearch(ci,c)
		ilo(7)=bsearch(di,d)

        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
        thetab=1d0-(bi-b(ilo(5)))/(b(ilo(5)+1)-b(ilo(5)))
        thetac=1d0-(ci-c(ilo(6)))/(c(ilo(6)+1)-c(ilo(6)))
		thetad=1d0-(di-d(ilo(7)))/(c(ilo(7)+1)-c(ilo(7)))
		

        fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))
		fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)
		fy(3)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))
		fy(4)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
		fy(5)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))
		fy(6)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
		fy(7)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
		fy(8)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(9) =thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(10)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
		fy(11)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))
        fy(12)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)		
        fy(13)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))		
        fy(14)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
		fy(15)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
		fy(16)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(17)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))
        fy(18)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)
        fy(19)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
        fy(20)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(21)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))
        fy(22)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
        fy(23)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
        fy(24)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(25)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(26)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
        fy(27)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)
        fy(28)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(29)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))
        fy(30)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
        fy(31)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
        fy(32)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		
        fy(33)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))
		fy(34)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)
		fy(35)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))
		fy(36)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
		fy(37)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))
		fy(38)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
		fy(39)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
		fy(40)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(41) =thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(42)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
		fy(43)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))
        fy(44)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)		
        fy(45)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))		
        fy(46)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
		fy(47)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
		fy(48)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(49)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))
        fy(50)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)
        fy(51)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
        fy(52)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(53)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))
        fy(54)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
        fy(55)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
        fy(56)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(57)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(58)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
        fy(59)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)
        fy(60)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(61)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))
        fy(62)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
        fy(63)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
        fy(64)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)

	
        
		
		
        fz(1)=thetay*fy(1)+(1-thetay)*fy(33)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(34)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(35)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(36)
        fz(5)=thetay*fy(5)+(1-thetay)*fy(37)
        fz(6)=thetay*fy(6)+(1-thetay)*fy(38)
        fz(7)=thetay*fy(7)+(1-thetay)*fy(39)
        fz(8)=thetay*fy(8)+(1-thetay)*fy(40)
        fz(9)=thetay*fy(9)+(1-thetay)*fy(41)
        fz(10)=thetay*fy(10)+(1-thetay)*fy(42)
        fz(11)=thetay*fy(11)+(1-thetay)*fy(43)
        fz(12)=thetay*fy(12)+(1-thetay)*fy(44)  
        fz(13)=thetay*fy(13)+(1-thetay)*fy(45)
        fz(14)=thetay*fy(14)+(1-thetay)*fy(46)
        fz(15)=thetay*fy(15)+(1-thetay)*fy(47)    
        fz(16)=thetay*fy(16)+(1-thetay)*fy(48) 
        fz(17)=thetay*fy(17)+(1-thetay)*fy(49)
        fz(18)=thetay*fy(18)+(1-thetay)*fy(50)
        fz(19)=thetay*fy(19)+(1-thetay)*fy(51)
        fz(20)=thetay*fy(20)+(1-thetay)*fy(52)
        fz(21)=thetay*fy(21)+(1-thetay)*fy(53)
        fz(22)=thetay*fy(22)+(1-thetay)*fy(54)
        fz(23)=thetay*fy(23)+(1-thetay)*fy(55)
        fz(24)=thetay*fy(24)+(1-thetay)*fy(56)
        fz(25)=thetay*fy(25)+(1-thetay)*fy(57)
        fz(26)=thetay*fy(26)+(1-thetay)*fy(58)
        fz(27)=thetay*fy(27)+(1-thetay)*fy(59)
        fz(28)=thetay*fy(28)+(1-thetay)*fy(60)  
        fz(29)=thetay*fy(29)+(1-thetay)*fy(61)
        fz(30)=thetay*fy(30)+(1-thetay)*fy(62)
        fz(31)=thetay*fy(31)+(1-thetay)*fy(63)    
        fz(32)=thetay*fy(32)+(1-thetay)*fy(64)         



        fa(1)=thetaz*fz(1)+(1-thetaz)*fz(17)
        fa(2)=thetaz*fz(2)+(1-thetaz)*fz(18)
        fa(3)=thetaz*fz(3)+(1-thetaz)*fz(19)
        fa(4)=thetaz*fz(4)+(1-thetaz)*fz(20)
        fa(5)=thetaz*fz(5)+(1-thetaz)*fz(21)
        fa(6)=thetaz*fz(6)+(1-thetaz)*fz(22)
        fa(7)=thetaz*fz(7)+(1-thetaz)*fz(23)
        fa(8)=thetaz*fz(8)+(1-thetaz)*fz(24)
        fa(9)=thetaz*fz(9)+(1-thetaz)*fz(25)
        fa(10)=thetaz*fz(10)+(1-thetaz)*fz(26)
        fa(11)=thetaz*fz(11)+(1-thetaz)*fz(27)
        fa(12)=thetaz*fz(12)+(1-thetaz)*fz(28)  
        fa(13)=thetaz*fz(13)+(1-thetaz)*fz(29)
        fa(14)=thetaz*fz(14)+(1-thetaz)*fz(30)
        fa(15)=thetaz*fz(15)+(1-thetaz)*fz(31)    
        fa(16)=thetaz*fz(16)+(1-thetaz)*fz(32) 
        
        
        fb(1)=thetaa*fa(1)+(1-thetaa)*fa(9)
        fb(2)=thetaa*fa(2)+(1-thetaa)*fa(10)
        fb(3)=thetaa*fa(3)+(1-thetaa)*fa(11)
        fb(4)=thetaa*fa(4)+(1-thetaa)*fa(12)
        fb(5)=thetaa*fa(5)+(1-thetaa)*fa(13)
        fb(6)=thetaa*fa(6)+(1-thetaa)*fa(14)
        fb(7)=thetaa*fa(7)+(1-thetaa)*fa(15)
        fb(8)=thetaa*fa(8)+(1-thetaa)*fa(16)


        fc(1)=thetab*fb(1)+(1-thetab)*fb(5)
        fc(2)=thetab*fb(2)+(1-thetab)*fb(6)
        fc(3)=thetab*fb(3)+(1-thetab)*fb(7)
        fc(4)=thetab*fb(4)+(1-thetab)*fb(8)
        
        fd(1)=fc(1)*thetac+fc(3)*(1-thetac)
        fd(2)=fc(2)*thetac+fc(4)*(1-thetac)

        fi=fd(1)*thetad+fd(2)*(1-thetad)
        
    end function interp7q
    
    function interp7qnoextrap(x,y,z,a,b,c,d,f,xi,yi,zi,ai,bi,ci,di) result(fi)
        implicit none
        real (8), intent(in) :: x(:),y(:),z(:),a(:),b(:),c(:),d(:),xi,yi,zi,ai,bi,ci,di,f(:,:,:,:,:,:,:)
        real (8) :: fi
        real (8) :: thetax,thetay,thetaz,thetaa,thetab,thetac,thetad,fy(64),fz(32),fa(16),fb(8),fc(4),fd(2)
        integer :: ilo(7),iI
        ilo(1)=bsearch(xi,x)
        ilo(2)=bsearch(yi,y)
        ilo(3)=bsearch(zi,z)
        ilo(4)=bsearch(ai,a)
        ilo(5)=bsearch(bi,b)
        ilo(6)=bsearch(ci,c)
		ilo(7)=bsearch(di,d)

        thetax=1d0-(xi-x(ilo(1)))/(x(ilo(1)+1)-x(ilo(1)))
        thetay=1d0-(yi-y(ilo(2)))/(y(ilo(2)+1)-y(ilo(2)))
        thetaz=1d0-(zi-z(ilo(3)))/(z(ilo(3)+1)-z(ilo(3)))
        thetaa=1d0-(ai-a(ilo(4)))/(a(ilo(4)+1)-a(ilo(4)))
        thetab=1d0-(bi-b(ilo(5)))/(b(ilo(5)+1)-b(ilo(5)))
        thetac=1d0-(ci-c(ilo(6)))/(c(ilo(6)+1)-c(ilo(6)))
		thetad=1d0-(di-d(ilo(7)))/(c(ilo(7)+1)-c(ilo(7)))
		
				thetax=max(0d0, min(1d0,thetax)) 
		thetay=max(0d0, min(1d0,thetay)) 
		thetaz=max(0d0, min(1d0,thetaz))
		thetaa=max(0d0, min(1d0,thetaa))
		thetab=max(0d0, min(1d0,thetab))
		thetac=max(0d0, min(1d0,thetac))
		thetad=max(0d0, min(1d0,thetad))
		
		!if (ilo(1)==size(x)) ilo(1)=ilo(1)-1 
		!if (ilo(2)==size(y)) ilo(2)=ilo(2)-1  
		!if (ilo(3)==size(z)) ilo(3)=ilo(3)-1 
		!if (ilo(4)==size(a)) ilo(4)=ilo(4)-1 
		!if (ilo(5)==size(b)) ilo(5)=ilo(5)-1 
		!if (ilo(6)==size(c)) ilo(6)=ilo(6)-1 
		!if (ilo(7)==size(d)) ilo(7)=ilo(7)-1 
		

      fy(1)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))
		fy(2)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)
		fy(3)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))
		fy(4)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
		fy(5)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))
		fy(6)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
		fy(7)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
		fy(8)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(9) =thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(10)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
		fy(11)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))
        fy(12)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)		
        fy(13)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))		
        fy(14)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
		fy(15)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
		fy(16)=thetax*f(ilo(1),ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(17)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))
        fy(18)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)
        fy(19)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
        fy(20)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(21)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))
        fy(22)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
        fy(23)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
        fy(24)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(25)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(26)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
        fy(27)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)
        fy(28)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(29)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))
        fy(30)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
        fy(31)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
        fy(32)=thetax*f(ilo(1),ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2),ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		
        fy(33)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7))
		fy(34)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6),ilo(7)+1)
		fy(35)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7))
		fy(36)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
		fy(37)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7))
		fy(38)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
		fy(39)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
		fy(40)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(41) =thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(42)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
		fy(43)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))
        fy(44)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)		
        fy(45)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))		
        fy(46)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
		fy(47)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
		fy(48)=thetax*f(ilo(1),ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3),ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(49)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7))
        fy(50)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6),ilo(7)+1)
        fy(51)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)
        fy(52)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(53)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7))
        fy(54)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6),ilo(7)+1)
        fy(55)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7))
        fy(56)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4),ilo(5)+1,ilo(6)+1,ilo(7)+1)
		
		fy(57)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7))
        fy(58)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6),ilo(7)+1)
        fy(59)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)
        fy(60)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5),ilo(6)+1,ilo(7)+1)	
        fy(61)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7))
        fy(62)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6),ilo(7)+1)
        fy(63)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7))
        fy(64)=thetax*f(ilo(1),ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)+(1-thetax)*f(ilo(1)+1,ilo(2)+1,ilo(3)+1,ilo(4)+1,ilo(5)+1,ilo(6)+1,ilo(7)+1)

	
        
		
		
        fz(1)=thetay*fy(1)+(1-thetay)*fy(33)
        fz(2)=thetay*fy(2)+(1-thetay)*fy(34)
        fz(3)=thetay*fy(3)+(1-thetay)*fy(35)
        fz(4)=thetay*fy(4)+(1-thetay)*fy(36)
        fz(5)=thetay*fy(5)+(1-thetay)*fy(37)
        fz(6)=thetay*fy(6)+(1-thetay)*fy(38)
        fz(7)=thetay*fy(7)+(1-thetay)*fy(39)
        fz(8)=thetay*fy(8)+(1-thetay)*fy(40)
        fz(9)=thetay*fy(9)+(1-thetay)*fy(41)
        fz(10)=thetay*fy(10)+(1-thetay)*fy(42)
        fz(11)=thetay*fy(11)+(1-thetay)*fy(43)
        fz(12)=thetay*fy(12)+(1-thetay)*fy(44)  
        fz(13)=thetay*fy(13)+(1-thetay)*fy(45)
        fz(14)=thetay*fy(14)+(1-thetay)*fy(46)
        fz(15)=thetay*fy(15)+(1-thetay)*fy(47)    
        fz(16)=thetay*fy(16)+(1-thetay)*fy(48) 
        fz(17)=thetay*fy(17)+(1-thetay)*fy(49)
        fz(18)=thetay*fy(18)+(1-thetay)*fy(50)
        fz(19)=thetay*fy(19)+(1-thetay)*fy(51)
        fz(20)=thetay*fy(20)+(1-thetay)*fy(52)
        fz(21)=thetay*fy(21)+(1-thetay)*fy(53)
        fz(22)=thetay*fy(22)+(1-thetay)*fy(54)
        fz(23)=thetay*fy(23)+(1-thetay)*fy(55)
        fz(24)=thetay*fy(24)+(1-thetay)*fy(56)
        fz(25)=thetay*fy(25)+(1-thetay)*fy(57)
        fz(26)=thetay*fy(26)+(1-thetay)*fy(58)
        fz(27)=thetay*fy(27)+(1-thetay)*fy(59)
        fz(28)=thetay*fy(28)+(1-thetay)*fy(60)  
        fz(29)=thetay*fy(29)+(1-thetay)*fy(61)
        fz(30)=thetay*fy(30)+(1-thetay)*fy(62)
        fz(31)=thetay*fy(31)+(1-thetay)*fy(63)    
        fz(32)=thetay*fy(32)+(1-thetay)*fy(64)         



        fa(1)=thetaz*fz(1)+(1-thetaz)*fz(17)
        fa(2)=thetaz*fz(2)+(1-thetaz)*fz(18)
        fa(3)=thetaz*fz(3)+(1-thetaz)*fz(19)
        fa(4)=thetaz*fz(4)+(1-thetaz)*fz(20)
        fa(5)=thetaz*fz(5)+(1-thetaz)*fz(21)
        fa(6)=thetaz*fz(6)+(1-thetaz)*fz(22)
        fa(7)=thetaz*fz(7)+(1-thetaz)*fz(23)
        fa(8)=thetaz*fz(8)+(1-thetaz)*fz(24)
        fa(9)=thetaz*fz(9)+(1-thetaz)*fz(25)
        fa(10)=thetaz*fz(10)+(1-thetaz)*fz(26)
        fa(11)=thetaz*fz(11)+(1-thetaz)*fz(27)
        fa(12)=thetaz*fz(12)+(1-thetaz)*fz(28)  
        fa(13)=thetaz*fz(13)+(1-thetaz)*fz(29)
        fa(14)=thetaz*fz(14)+(1-thetaz)*fz(30)
        fa(15)=thetaz*fz(15)+(1-thetaz)*fz(31)    
        fa(16)=thetaz*fz(16)+(1-thetaz)*fz(32) 
        
        
        fb(1)=thetaa*fa(1)+(1-thetaa)*fa(9)
        fb(2)=thetaa*fa(2)+(1-thetaa)*fa(10)
        fb(3)=thetaa*fa(3)+(1-thetaa)*fa(11)
        fb(4)=thetaa*fa(4)+(1-thetaa)*fa(12)
        fb(5)=thetaa*fa(5)+(1-thetaa)*fa(13)
        fb(6)=thetaa*fa(6)+(1-thetaa)*fa(14)
        fb(7)=thetaa*fa(7)+(1-thetaa)*fa(15)
        fb(8)=thetaa*fa(8)+(1-thetaa)*fa(16)


        fc(1)=thetab*fb(1)+(1-thetab)*fb(5)
        fc(2)=thetab*fb(2)+(1-thetab)*fb(6)
        fc(3)=thetab*fb(3)+(1-thetab)*fb(7)
        fc(4)=thetab*fb(4)+(1-thetab)*fb(8)
        
        fd(1)=fc(1)*thetac+fc(3)*(1-thetac)
        fd(2)=fc(2)*thetac+fc(4)*(1-thetac)

        fi=fd(1)*thetad+fd(2)*(1-thetad)
        
    end function interp7qnoextrap
    
    function excess_demand1(xfun,Agg1,Agg2,x) result(Zfun)
        implicit none
        real(8), intent(in) :: x(:) ,xfun
        real(8) :: Zfun
		real(8),intent(in) :: Agg1(:), Agg2(:)
        Zfun=interp1q(x,Agg1-Agg2,xfun)
    end function excess_demand1
    
   function excess_demand2(xfun,Agg,x) result(Zfun)
        implicit none
        real(8), intent(in) :: x(:), Agg(:),xfun
        real(8) :: Zfun
        Zfun=interp1q(x,Agg,xfun)
    end function excess_demand2


end module mod_interp
