module mod_functions
    
    use parameters
	use mod_globals
	use mod_interp
    implicit none
    contains
	
	! FROM KMV (2020)
	double precision function MinPay(mort,age,rm)

double precision, intent(in)    :: mort,rm
integer, intent(in)             :: age
integer                         :: payleft

    payleft = T-age+1+payextra
    MinPay = mort*rm*((1.0d0+rm)**payleft)/((1.0d0+rm)**payleft - 1.0d0)

end function MinPay

    ! From Kaplan, Mitman, Violante 2020

    ! Progressive Tax Function of Gouveia-Strauss
    double precision function FnTax(ly)
		IMPLICIT NONE
        double precision, intent(in)    :: ly
        FnTax=ly-tau0*ly**(1.0d0-tau1)
    end function FnTax
	
	! From Kaplan, Mitman, Violante 2020

    ! Progressive Tax Function of Gouveia-Strauss
	double precision function FnTaxm(ly,mort,loc,rm)
		double precision, intent(in) :: ly, mort, loc,rm
		double precision                :: lyp
		double precision                :: intpay
		!Calculate how much of payment was interest
		intpay = (rm*mort+rm*loc)*MortDeduct
		lyp = ly-intpay
		lyp = max(lyp,0.0d0)
		FnTaxM=min(ly-tau0*ly**(1.0d0-tau1),lyp-tau0*lyp**(1.0d0-tau1))
	end function FnTaxm

    !Income function
    double precision function income(ichi,ie,iTheta)
		IMPLICIT NONE
        double precision, intent(in) :: ichi
        double precision, intent(in) :: ie
        double precision, intent(in) :: iTheta
		double precision :: yavall
		
		! This comes directly from KMV
		yavall=1.33396176336787d0 !1.33396176336787
        income=exp(ichi+ie+log(iTheta))/yavall
    end function income
	
	double precision function income_ret(ichi,ie,iTheta)
		IMPLICIT NONE
        double precision, intent(in) :: ichi
        integer, intent(in) :: ie
		double precision, intent(in) :: iTheta
		double precision :: yavall

		if (iE==1) then
			income_ret=0.370902275775089d0
		else if (ie==2) then
			income_ret=0.377022630253217d0
		else if (iE==3) then
			income_ret=0.389815647952213d0
		else if (iE==4) then
			income_ret=0.416556141894462d0
		else if (iE==5) then
			income_ret=0.472450226706960d0
		else 
			print*, 'error in retirement income function'
		end if
    end function income_ret
	
	double precision function Value1beq(bsav1beq)
		implicit none
		real(8), intent(in) :: bsav1beq
		!print*, cval4,Rval4
		!print*, 'cval',cval1
		!if (cval1-bsav1*q_lend>=cmin) then
			Value1beq=-(eiJ**(sigma-1d0)*(((1d0-phibeq)*(cval1beq-bsav1beq*q_lend)**(1d0-gamma_b)+phibeq*(Hval1)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)-psi*((bsav1beq+bflat)**(1d0-sigma))/(1d0-sigma)
			!print*, eiJ,(1d0-phibeq),cval1beq-bsav1beq*q_lend,q_lend, Hval1, bsav1beq, phibeq
			
			!-psi*((bsav4+bflat)**(1d0-sigma))/(1d0-sigma)
			!eiJ**(sigma-1.0d0)*
		!	print*, 'in Value4', eiJ, (1d0-phi)*(cval4-bsav4*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b),interp1q(B,Wval4,bsav4)
		!else
		!	Value1beq=-death
		!end if
	end function Value1beq
	
	double precision function Value1(bsav1)
		implicit none
		real(8), intent(in) :: bsav1
		!print*, cval4,Rval4

		if (cval1+interp1q(B,qval1,bsav1)-bsav1*q_lend>cmin) then
			Value1=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval1+interp1q(B,qval1,bsav1)-bsav1*q_lend)**(1d0-gamma_b)+phi*(Hval1)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-interp1qnoextrap(B,Wval1,bsav1)
			
			
		else
			Value1=-death
		end if
		!print*, cval1,interp1q(B,qval1,bsav1),bsav1*q_lend,Value1
	end function Value1
		

!print*, 'in Value4',Value4,cval4, bsav,cval4-bsav*q_lend, (1d0-phi)*(cval4-bsav*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b),interp1q(B,Wval,bsav)
	


	double precision function Value2(bsav2)
		implicit none
		double precision, intent(in) :: bsav2
		!print*, cval4,Rval4

		if (bsav2>0d0) then 
			if (cval2-bsav2*q_lend>cmin) then
				Value2=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval2-bsav2*q_lend)**(1d0-gamma_b)+phi*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-interp1qnoextrap(B,Wval2,bsav2)
			else
				Value2=-death
			end if
		else 
			if (cval2-bsav2*q_borrow>cmin) then 
			Value2=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval2-bsav2*q_borrow)**(1d0-gamma_b)+phi*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-interp1qnoextrap(B,Wval2,bsav2)
			else
				Value2=-death
			end if
			
		end if

!print*, 'in Value2',Value2,cval2, bsav2,cval2-bsav2*q_lend, (1d0-phi)*(cval2-bsav2*q_lend)**(1d0-gamma_b)+phi*(Hval2)**(1d0-gamma_b),interp1q(B,Wval2,bsav2)
	end function Value2
	
	
		double precision function Value2beq(bsav2)
		implicit none
		double precision, intent(in) :: bsav2
		!print*, cval4,Rval4
		if (bsav2>0d0) then
			if (cval2-bsav2*q_lend>=cmin) then
				Value2beq=-(eiJ**(sigma-1d0)*(((1d0-phibeq)*(cval2-bsav2*q_lend)**(1d0-gamma_b)+phibeq*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-psi*((bsav2+other+bflat)**(1d0-sigma))/(1d0-sigma)
			else
				Value2beq=-death
			end if
			!print*, 'in Value2',Value2beq,bsav2,cval2-bsav2*q_lend, bsav2+other,(eiJ**(sigma-1d0)*(((1d0-phi)*(cval2-bsav2*q_lend)**(1d0-gamma_b)+phi*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma),-psi*((bsav2+other+bflat)**(1d0-sigma))/(1d0-sigma)

		else if (bsav2>Bheloc2) then 
			if (cval2-bsav2*q_borrow>=cmin) then
				Value2beq=-(eiJ**(sigma-1d0)*(((1d0-phibeq)*(cval2-bsav2*q_borrow)**(1d0-gamma_b)+phibeq*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-psi*((bsav2+other+bflat)**(1d0-sigma))/(1d0-sigma)
			else
				Value2beq=-death
			end if
			!print*, 'in Value2',Value2beq,bsav2,cval2-bsav2*q_borrow, bsav2+other,(eiJ**(sigma-1d0)*(((1d0-phi)*(cval2-bsav2*q_lend)**(1d0-gamma_b)+phi*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma),-psi*((bsav2+other+bflat)**(1d0-sigma))/(1d0-sigma)
		else 
			if (cval2>=cmin) then
				Value2beq=-(eiJ**(sigma-1d0)*(((1d0-phibeq)*(cval2)**(1d0-gamma_b)+phibeq*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-psi*((bsav2+other+bflat)**(1d0-sigma))/(1d0-sigma)
			else
				Value2beq=-death
			end if
			!print*, 'in Value2',Value2beq,bsav2,cval2-bsav2*q_borrow, bsav2+other,(eiJ**(sigma-1d0)*(((1d0-phi)*(cval2-bsav2*q_lend)**(1d0-gamma_b)+phi*(Hval2)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma),-psi*((bsav2+other+bflat)**(1d0-sigma))/(1d0-sigma)

		end if 
	end function Value2beq

	double precision function Value3(bsav3)
		implicit none
		double precision, intent(in) :: bsav3
		!print*, cval4,Rval4

		if (bsav3>0d0) then 
			if (cval3+interp1q(B,qval3,bsav3)-bsav3*q_lend>cmin) then
				Value3=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval3+interp1q(B,qval3,bsav3)-bsav3*q_lend)**(1d0-gamma_b)+phi*(Hval3)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-interp1qnoextrap(B,Wval3,bsav3)
			!eiJ**(sigma-1.0d0)*
			
			else
				Value3=-death
			end if
			!print*, bsav3, Value3, cval3, interp1q(B,qval3,bsav3),bsav3*q_lend
		ELSE
			if (cval3+interp1q(B,qval3,bsav3)-bsav3*q_borrow>cmin) then
				Value3=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval3+interp1q(B,qval3,bsav3)-bsav3*q_borrow)**(1d0-gamma_b)+phi*(Hval3)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
				-interp1qnoextrap(B,Wval3,bsav3)
			!eiJ**(sigma-1.0d0)*
			
			else
				Value3=-death
			end if
			!print*, bsav3, Value3, cval3, interp1q(B,qval3,bsav3),bsav3*q_borrow
		end if 
!print*, bsav3, Value3, cval3, interp1q(B,qval3,bsav3),bsav3*q_lend
!print*, 'in Value4',Value4,cval4, bsav,cval4-bsav*q_lend, (1d0-phi)*(cval4-bsav*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b),interp1q(B,Wval,bsav)
	end function Value3

	recursive double precision function Value4beq(bsav4)
		implicit none
		double precision, intent(in) :: bsav4
		!print*, cval4,Rval4

		if (cval4-bsav4*q_lend>=cmin) then
			Value4beq=-(eiJ**(sigma-1d0)*(((1d0-phibeq)*(cval4-bsav4*q_lend)**(1d0-gamma_b)+phibeq*(Rval4)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-psi*((bsav4+bflat)**(1d0-sigma))/(1d0-sigma)
			!eiJ**(sigma-1.0d0)*
		!	print*, 'in Value4',eiJ**(sigma-1d0),(1d0-phi)*(cval4-bsav4*q_lend)**(1d0-gamma_b),phi*(Rval4)**(1d0-gamma_b),(cval4-bsav4*q_lend)**(1d0-gamma_b),Rval4**(1d0-gamma_b),&
		!	(eiJ**(sigma-1d0)*(((1d0-phi)*(cval4-bsav4*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma),1d0-gamma_b,phi
		else
			Value4beq=-death
		end if

!print*, 'in Value4',Value4,cval4, bsav,cval4-bsav*q_lend, (1d0-phi)*(cval4-bsav*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b),interp1q(B,Wval,bsav)
	end function Value4beq


	recursive double precision function Value4(bsav4)
		implicit none
		double precision, intent(in) :: bsav4
		!print*, cval4,Rval4

		if (cval4-bsav4*q_lend>=cmin) then
			Value4=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval4-bsav4*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-interp1qnoextrap(B,Wval4,bsav4)
			!eiJ**(sigma-1.0d0)*
		!	print*, 'in Value4', eiJ, (1d0-phi)*(cval4-bsav4*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b),interp1q(B,Wval4,bsav4)
		else
			Value4=-death
		end if

!print*, 'in Value4',Value4,cval4, bsav,cval4-bsav*q_lend, (1d0-phi)*(cval4-bsav*q_lend)**(1d0-gamma_b)+phi*(Rval4)**(1d0-gamma_b),interp1q(B,Wval,bsav)
	end function Value4
	
	recursive double precision function Value5(bsav5)
		implicit none
		double precision, intent(in) :: bsav5

		if (cval5-bsav5*q_lend>cmin) then
			Value5=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval5-bsav5*q_lend)**(1d0-gamma_b)+phi*(Rval5)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-interp1qnoextrap(B,Wval5,bsav5)+default_cost
			!eiJ**(sigma-1.0d0)*
		!	print*, 'in Value5', eiJ, (1d0-phi)*(cval5-bsav5*q_lend)**(1d0-gamma_b)+phi*(Rval5)**(1d0-gamma_b),interp1q(B,Wval5,bsav5)
		else
			Value5=-deathdefault
		end if

	end function Value5


	recursive double precision function Value5beq(bsav5)
		implicit none
		double precision, intent(in) :: bsav5

		if (cval5-bsav5*q_lend>=cmin) then
			Value5beq=-(eiJ**(sigma-1d0)*(((1d0-phibeq)*(cval5-bsav5*q_lend)**(1d0-gamma_b)+phibeq*(Rval5)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-psi*((bsav5+bflat)**(1d0-sigma))/(1d0-sigma)+default_cost
			!eiJ**(sigma-1.0d0)*
			!print*, 'in Value5', eiJ, (1d0-phi)*(cval5-bsav5*q_lend)**(1d0-gamma_b)phi*(Rval5)**(1d0-gamma_b),interp1q(B,Wval5,bsav5)
		else
			Value5beq=-deathdefault !KMV set this value to be higher for foreclosure
		end if

	end function Value5beq

	double precision function Value6(bsav6)
		implicit none
		real(8), intent(in) :: bsav6
		!print*, cval4,Rval4

		if (cval6+interp1q(BS,qval6,bsav6)-bsav6*q_lend>cmin) then
			Value6=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval6+interp1q(BS,qval6,bsav6)-bsav6*q_lend)**(1d0-gamma_b)+phi*(Hval6)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-interp1qnoextrap(BS,Wval6,bsav6)
			
			
		else
			Value6=-death
		end if
		!print*, cval1,interp1q(B,qval1,bsav1),bsav1*q_lend,Value1
	end function Value6

	recursive double precision function Value7(bsav7)
		implicit none
		double precision, intent(in) :: bsav7
		!print*, cval7,Rval7

		if (cval7-bsav7*q_lend>=cmin) then
			Value7=-(eiJ**(sigma-1d0)*(((1d0-phi)*(cval7-bsav7*q_lend)**(1d0-gamma_b)+phi*(Rval7)**(1d0-gamma_b))**((1d0-sigma)/(1d0-gamma_b)))-1d0)/(1d0-sigma)&
			-interp1qnoextrap(BS,Wval7,bsav7)
			!eiJ**(sigma-1.0d0)*
		!	print*, 'in Value7', eiJ, (1d0-phi)*(cval7-bsav7*q_lend)**(1d0-gamma_b)+phi*(Rval7)**(1d0-gamma_b),interp1q(B,Wval7,bsav7)
		else
			Value7=-death
		end if

!print*, 'in Value7',Value7,cval7, bsav,cval7-bsav*q_lend, (1d0-phi)*(cval7-bsav*q_lend)**(1d0-gamma_b)+phi*(Rval7)**(1d0-gamma_b),interp1q(B,Wval,bsav)
	end function Value7
	

end module mod_functions


