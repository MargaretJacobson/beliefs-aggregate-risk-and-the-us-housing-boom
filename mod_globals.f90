module mod_globals
	use parameters
	implicit none
	
	real(8), dimension(nB) :: Wval1,Wval2,Wval3,Wval4,Wval5,B,Wval6,Wval7
	real(8), dimension(nBS) :: BS
	real(8), dimension(nB) :: qval1,qval3
	real(8), dimension(nBS) :: qval6
	real(8) :: Hval1,Hval2,Hval3,Rval4, Rval5,Hval6,Rval7
	real(8) ::cval1beq,cval1, cval2,cval3,cval4, cval5,cval6,cval7 ,phi,q_lend,q_borrow
	real(8) :: eiJ,OTHER, BHELOC2,phibeq
	common /Bequest/ cval1beq,BHELOC2,phibeq
	common /VFI/ eiJ,qval1,cval1,Wval1,qval3,cval2,cval3,cval4,cval5,Hval1,Hval2,Hval3,Rval4,Rval5,Wval2,Wval3,Wval4,Wval5,&
	OTHER,phi,q_lend,q_borrow,qval6,cval6,Wval6,Hval6,cval7,Rval7,Wval7
	
	!$OMP THREADPRIVATE(/Bequest/,/VFI/)
	
	
	!common /buy/ cval1,cval2,cval3,Hval1,Hval2,Hval3,Wval1,Wval2,Wval3,cval1beq,other,bheloc2,qval1,qval3
	!common /rent/ cval4,Rval4,cval5,Rval5, Wval4, Wval5
	!common /McClem/ eiJ,iJ
	
	!!$OMP THREADPRIVATE(/buy/,/rent/,/McClem/)
end module mod_globals