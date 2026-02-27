
! Question 1.2 Problem Set 1

! By Margaret Jacobson, Indiana University marmjaco@indiana.edu
! For E 724: Computational Macroeconomics taught by Grey Gordon
! in Spring 2016 at Indiana University

! To write this code, I consulted lecture notes by Grey Gordon
! entitled, "Basics", "Arrays", and  "Functions, Subroutines,
! and Modules," accessed via Grey Gordon's website
! https://sites.google.com/site/greygordon/teaching and Canvas. I also
! used course notes from E551: Quantitative Macroeconomics with
! Heterogeneity  taught by Professor Bulent Guler at Indiana
! University in Fall 2015

! On slides 35-38 from the slides "Functions, Subroutines, and
! Modules", Professor Gordon had code describing the tauchen
! routine which I also consulted.

! This file contains the module mod_discreteAR1
! mod_discreteAR1 contains the subroutine sub_tauchen

! I also consulted MATLAB code for the tauchen procedure,
! "Tauchen.m" by Martin Floden, 1996 which implements Geroge
! Tauchen's algorithm described in Economic Letters 20(1986) 177-181

module mod_discreteAR1
    use parameters
    use mod_matlab
    implicit none
contains
    subroutine sub_tauchen(logy,Fyy,cdnal_stdev,cover)
    ! Want to discretize the process: logy=rho*log(y-1)+sigma*epsilon     where epsilon ~N(0,1)
    ! rho=ar1 and sigma=cdna1_stdev
        implicit none
        real(8), intent(in) :: cover,cdnal_stdev
        real(8), intent(out) :: logy(:), Fyy(:,:)
        real(8) :: sigma, ylbar, yubar
        integer :: ii, jj, n
               
        n=size(logy)
         
        ! Checking that logy and Fyy conform
        if (size(logy)==size(Fyy,2)) then

        else
            stop 'Non conforming matrices'
        end if

        sigma=sqrt(cdnal_stdev)
        ylbar=-cover*sigma
        yubar=cover*sigma 
        
        ! Taking care of the case where n==1
        if (n==1) then
            logy=0d0 !cover*sigma !only use this when the standard deviation is set to zero
            Fyy(1,1)=1d0
        else

        ! Constructing grid points 
        logy=linspace(ylbar,yubar,n)
        
            ! Constructing transition matrix
 
            do ii= 1,n
                do jj= 2,n-1
                    Fyy(ii,jj)=normcdf((logy(jj+1)+logy(jj))/2d0,sigmae**2d0, rhoe*logy(ii))
                    Fyy(ii,jj)=Fyy(ii,jj)-normcdf((logy(jj)+logy(jj-1))/2d0, sigmae**2d0, rhoe*logy(ii))
                end do
                Fyy(ii,1)=normcdf((logy(2)+logy(1))/2d0,sigmae**2d0,rhoe*logy(ii))
                Fyy(ii,n)=1d0-normcdf((logy(n)+logy(n-1))/2d0, sigmae**2d0, rhoe*logy(ii))
                
            end do
        end if
        
        ! Normalizing the rows of Fyy     
        ! Checking to make sure the rows of Fyy sum to one
        do ii=1,n            
          Fyy(ii,:)=Fyy(ii,:)/sum(Fyy(ii,:))

        if (abs(sum(Fyy(ii,:))-1d0)>10d-10) then
            print*, 'sum less 1', sum(Fyy(ii,:))-1d0
            print*, 'rows', Fyy(ii,:)
            stop 'rows did not sum to 1'
        end if
        end do
    
    end subroutine sub_tauchen

end module mod_discreteAR1
