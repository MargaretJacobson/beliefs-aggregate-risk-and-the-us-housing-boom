module mod_rouwenhorst
    use parameters
    use mod_matlab
    implicit none
contains    
    subroutine sub_rouwenhorst(z,P,rho,sig,sigmaold,sigmanew)
    implicit none
    real(8), intent(in) :: rho,sig,sigmaold,sigmanew
    real(8), intent(out) :: z(nE),P(nE,nE)
    !local
    real(8) :: mu_eps, q, nu
    real(8), allocatable :: P0(:,:),P1(:,:)
    integer :: i
    ! Rouwenhorst method as presented in Kopecky and Suen RED 2010
    ! Based on Karen Kopecky's Code via "Optimal Bankruptcy Code: A Fresh Start for Some" by Grey Gordon
    ! Journal of Economic Dynamics and Control (2017), 85:123-149
    ! Non Stationary code follows from http://www.giuliofella.net/research/pdf/fgp.pdf
    ! "Markov-Chain Approximations for Life-Cycle Models," by Giulio Fella, Giovanni     Gallipoli, and Jutong Pan working paper, 2018
    if (nE<=0) then
        STOP 'sub_rouwenhorst: n must be > 0'
    elseif (nE==1) then
        z = 0d0
        P = 1d0
        return
    end if
    
    mu_eps =0d0
    q = (rho*(sigmaold/sigmanew)+1d0)/2d0
    if (sigmaold==sigmanew) then
        nu = sqrt((real(nE)-1d0)/(1d0-rho**2))*sig
    else
        nu=sqrt(real(nE)-1d0)*sigmanew
    end if
    allocate(P0(2,2))
    P0(1,:) = (/q, 1d0-q/)
    P0(2,:) = (/1d0-q, q/)
    
    do i = 2,nE-1
        allocate(P1(i+1,i+1))
        P1 = 0d0
        !        
        P1(1:i,1:i) = q*P0 + P1(1:i,1:i)
        !
        P1(1:i,2:i+1) = (1d0-q)*P0 + P1(1:i,2:i+1)
        !
        P1(2:i+1,1:i) = (1d0-q)*P0 + P1(2:i+1,1:i)
        !
        P1(2:i+1,2:i+1) = q*P0 + P1(2:i+1,2:i+1)
        !
        P1(2:i,:) = P1(2:i,:)/2d0
        deallocate(P0)
        allocate(P0(i+1,i+1))
        P0 = P1
        deallocate(P1)
    end do
    
    P = P0
    deallocate(P0)
    
    z = linspace(mu_eps/(1d0-rho) - nu, mu_eps/(1d0-rho) + nu,nE)
    

    end subroutine sub_rouwenhorst
end module mod_rouwenhorst