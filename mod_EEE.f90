!       This file contains the simulation code for Krussel-Smith Simulation

module mod_EEE
    use parameters
    use mod_matlab
    use mod_interp
    implicit none
contains
    
    subroutine sub_EEE(FYY,Fee,q_b,LLi_b,HHi_b,c_b,HH_b,pprime,p,EEEq,EEEp)
        implicit none
        integer :: iL,iH,iY,iE,iYY,iEE,iP,iJ,iA1,iA2
        real(8) :: cval_prime, hval_prime, cval, hval_prime_prime, nume, den,phi_in
        real(8), intent(in) :: FYY(:,:), Fee(:,:)
        real(8), intent(in) :: q_b(:,:,:,:,:,:), c_b(:,:,:,:,:,:), HH_b(:,:,:,:,:,:)
        integer, intent(in) :: LLi_b(:,:,:,:,:,:), HHi_b(:,:,:,:,:,:)
        real(8), intent(in) :: pprime(:,:,:), p(:)
        real(8), intent(out) :: EEEq(:,:,:,:,:,:), EEEp(:,:,:,:,:,:)
        
        EEEq=0d0
        EEEp=0d0
        phi_in=phi_grid(1)
        do iL=1,nL
            do iH=1,nH
                do iY=1,nY
                    do iE=1,nE
                        do iP=1,nP
                            do iJ=1,T-1
                                do iYY=1,nY
                                    do iEE=1,nE
                                    ! need to interpolate on pprime
                                        cval_prime=interp1q(p,c_b(LLi_b(iL,iH,iY,iE,iP,iJ),&
                                          HHi_b(iL,iH,iY,iE,iP,iJ),iYY,iEE,:,iJ+1),pprime(iY,iYY,iP))
                                        hval_prime_prime=interp1q(p,HH_b(LLi_b(iL,iH,iY,iE,iP,iJ),&
                                          HHi_b(iL,iH,iY,iE,iP,iJ),iYY,iEE,:,iJ+1),pprime(iY,iYY,iP))

                                
                                        nume=(beta*q_b(iL,iH,iY,iE,iP,iJ))*FYY(iY,iYY)*Fee(iE,iEE)*&
                                          cval_prime**(-gamma_b)*((1-phi_in)*cval_prime**(1-gamma_b)+phi_in*&
                                          hval_prime_prime**(1-gamma_b))**((1d0-sigma)/(1d0-gamma_b)-1d0)
                                
                                            
                                        EEEq(iL,iH,iY,iE,iP,iJ)=nume+EEEq(iL,iH,iY,iE,iP,iJ)
                                        EEEp(iL,iH,iY,iE,iP,iJ)=beta*FYY(iY,iYY)*Fee(iE,iEE)*cval_prime*&
                                          (1-delta)*pprime(iY,iYY,iP)
                                    end do
                                end do
                                cval=c_b(iL,iH,iY,iE,iP,iJ)
                                hval_prime=HH_b(iL,iH,iY,iE,iP,iJ)
                                den=cval**(-gamma_b)*((1-phi_in)*cval**(1-gamma_b)+phi_in*hval_prime**&
                                            (1-gamma_b))**((1d0-sigma)/(1d0-gamma_b)-1d0)
                                EEEq(iL,iH,iY,iE,iP,iJ)=log10(abs(1d0- EEEq(iL,iH,iY,iE,iP,iJ)/den)) 
                                EEEp(iL,iH,iY,iE,iP,iJ)=log10(abs(1d0-((1-gamma_b)/(gamma_b*p(iP))*(hval_prime/cval)&
                                **(-gamma_b)+EEEp(iL,iH,iY,iE,iP,iJ)/(p(iP)*cval**(-gamma_b)))))
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine sub_EEE
end module mod_EEE