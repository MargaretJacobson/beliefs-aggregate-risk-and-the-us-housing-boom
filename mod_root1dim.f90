! Questions 3&4 Problem Set 3

! By Margaret Jacobson, Indiana University marmjaco@indiana.edu
! For E 724: Computational Macroeconomics taught by Grey Gordon in
! Spring 2016 at Indiana University

! To write this code, I consulted lecture notes and code entitled
! 'mod_root1dimf.f90' by Grey Gordon accessed via Grey Gordon's website
! https://sites.google.com/site/greygordon/teaching and Canvas. 

! This file contains the module mod_root1dim
! mod_root1dim contains the type opt_root1dim and the subroutine
! sub_bisect that solves f(x)=0 using the bisection method

module mod_root1dim
    use mod_plot
    implicit none
    type opt_root1dim
        ! The RHS values are defaults for the type
        real (8) :: xtol = 1d-6 ! x tolerance
        integer :: maxit = 1000 ! maximum number of iterations
    end type

contains
    subroutine sub_bisect (x,f,a,b,opt_)
        implicit none
        real (8) , intent ( inout ) :: x
        interface
            function f(x) result (y)
                implicit none
                real (8) , intent (in) :: x
                real (8) :: y 
            end function f
        end interface
        real (8) , intent (in) :: a,b
        type (opt_root1dim) :: opt
        type (opt_root1dim), intent (in), optional :: opt_
        real (8) :: m, twoeps, aa, bb, xtol
        integer :: maxit, it
        
        ! Setting tolerance and max iterations
        if (present(opt_)) opt=opt_
        xtol=opt%xtol
        maxit=opt%maxit  
         
        !Checking to see that the root is bracketed by a,b
        if (a==b) then
            stop 'a=b, root not bracketed'
        end if
        if (sign(f(a),f(b))==f(a)) then
            stop 'Not oposite signs, root not bracketed'
        end if
        
        twoeps=2d0*epsilon(0d0)
        m=.5*(b+a)
        aa=a
        bb=b
        it=0
        do while (abs(bb-aa)>=abs(m)*twoeps+xtol .and. it<maxit)
            it=it+1
            m=.5*(aa+bb)
            if (f(m)==0) then
                exit
            else if (sign(f(m),f(aa))==f(m)) then
                    aa=m
            else
                    bb=m
            end if
            
        end do
        x=m 
    end subroutine sub_bisect   

    subroutine sub_bisect_excess_demand1(x,f,xgrid,Y1,Y2,a,b,opt_)
        implicit none
        real (8) , intent ( inout ) :: x
        real (8) , intent (in) :: xgrid(:)
		real (8), intent(in) :: Y1(:),Y2(:)
        interface
            function f(x,Y1,Y2,xgrid) result (y)
                implicit none
                real (8) , intent (in) :: x, xgrid(:)
				real (8), intent (in) :: Y1(:),Y2(:)
                real (8) :: y
            end function f
        end interface
        real (8) , intent (in) :: a,b
        type (opt_root1dim) :: opt
        type (opt_root1dim), intent (in), optional :: opt_
        real (8) :: m, twoeps, aa, bb, xtol
        integer :: maxit, it

        ! Setting tolerance and max iterations
        if (present(opt_)) opt=opt_
        xtol=opt%xtol
        maxit=opt%maxit

        !Checking to see that the root is bracketed by a,b
        if (a==b) then
            print*, 'a,b:    ', a, b
            stop 'a=b, root not bracketed'
        end if
        if (sign(f(a,Y1,Y2,xgrid),f(b,Y1,Y2,xgrid))==f(a,Y1,Y2,xgrid)) then
            !call plot(xgrid,Y1-Y2,'Excess Demand: Demand less Supply')
            print*, 'p', xgrid
            print*, 'bounds', a, b
            print*, 'Aggd less Aggs', Y1-Y2
            stop 'Not oposite signs, root not bracketed'
        end if

        twoeps=2d0*epsilon(0d0)
        m=.5*(b+a)
        aa=a
        bb=b
        it=0
                do while (abs(bb-aa)>=abs(m)*twoeps+xtol .and. it<maxit)
            it=it+1
            m=.5*(aa+bb)
            if (f(m,Y1,Y2,xgrid)==0) then
                exit
            else if (sign(f(m,Y1,Y2,xgrid),f(aa,Y1,Y2,xgrid))==f(m,Y1,Y2,xgrid)) then
                    aa=m
            else
                    bb=m
            end if

        end do
        x=m
    end subroutine sub_bisect_excess_demand1

end module mod_root1dim

