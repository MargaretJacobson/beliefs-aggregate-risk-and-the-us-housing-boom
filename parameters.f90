module parameters
    implicit none
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		Which model to run parameters
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer, parameter :: solve_for_coeffs=0	! if =1 then it will iterate on the coefficients. Otherwise =0 it will not. Use 0 for steady state, learning, and housing boom simulation. Setting=1 will put the credit condition states to 1 automatically
	integer, parameter :: high_nlamb=1			! Sets the credit conditions to their high value. Otherwise=0 and it will be their low value
	integer, parameter :: type_sim=2			!0=Steady state, 1=Standard KS, 2= KS with learning
	integer, parameter :: display=0				! display, if 0 then displays do not print
	integer, parameter :: no_Agg=0				! =1 no aggregate income risk
	integer, parameter :: no_Agg_lend=0		! =1 no aggregate credit condition risk
	integer, parameter :: no_MIT=0		! =0 MIT shocks, =1 no MIT, i.e. Markov credit conditions
	integer, parameter :: perfect_corr=0 		! =1 perfectly correlated income and credit condition risk
	integer, parameter :: initial_dist_opt=1	! alternative initial distributions
	integer, parameter :: no_pref=1	! =1 no preference shocks, =0 allows for preference shocks
	integer, parameter :: no_agg_phi=0 !(1-no_MIT)
	integer, parameter :: low_rate=0			! Adds low interest rates to a shift in credit conditions
	integer, parameter :: shrink=1		! Set to 1 for normal sized grids and agents, set to 2, 3, ... etc to shrink grid points
	integer, parameter :: hinterp=1		! Set to 1 to interpolate over H dimension, set to 1 otherwise
	integer, parameter :: buy_grid=0	! Set to 0 for no exta interpolation
	integer, parameter :: high_rate=1
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		Grid Parameters
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	! Got rid of interpolation. Try changing beta and rental rate
	integer, parameter :: nBsim=int(76) !26 ! int(51) !int(76)
	integer, parameter :: nLsim=int(22/shrink)	!int(64) !22 !int(43)  !int(64)
    integer, parameter :: nP=int(13/shrink)												! Number of house price grid points
    integer, parameter :: nL=int(22/shrink)										! Number of mortgage grid poitns 
    integer, parameter :: nB=int(26)													! Number of liquid financial instruments grid point
	integer, parameter :: nBzero=int(7)												! where zero starts on the liquid financial instruments grid
    integer, parameter :: nBs=6
	integer, parameter :: nH=6 													! Number of owner occupied housing grid points 
    integer, parameter :: nR=3													! Number of rental housing grid points
    integer, parameter :: nE=5     												! Individual income grid points 
	integer, parameter :: nY=2													! Number of Aggregate income grid points, values other than 2 have not been tested
	integer, parameter :: nlamb=2*(1-solve_for_coeffs)+2*solve_for_coeffs*no_MIT +1*(1-no_MIT)*solve_for_coeffs				 	! Can be 1 or 2. Must be =1 when solving coefficients and not markov	
	integer, parameter :: nphi=no_pref+(1-no_pref)*3
	real(8), parameter :: curve=2d0 											! Grid curvature
	integer, parameter :: KMV=1													! =1 sets grids to exactly what KMV has, 0= an approximation
	integer, parameter :: mom=1													! Number of moments in the house price forecast equation, moments other than 1 have not been fully tested 		
	real(8), parameter :: Numeraire = 0.5d0    									!1 for annual av earns, 4 for quarterly av earns, 0.5d0 for two year earnings which is the model period
	real(8), parameter :: DataAvAnnualEarns = 53000.0d0							!Average annual earnings from the 1998 SCFR

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		Learning Parameters
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer, parameter :: learn=200 !100											! Start of learning simulation 
    real(8), parameter, dimension(3) :: adj=(/-0.08d0,0d0,0.16d0/) ! (/-0.04d0,0d0,0.1d0/) !(/-0.04d0,0d0,0.07d0,0.13d0/) ! !(/-0.02d0,0.02d0,0.06d0/) !	! Width of belief grid, must be centered at 0
    real(8), parameter :: gain_param=0.508d0 
	real(8), parameter :: gain2=0.508d0
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 		Simulation paramters
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: T=30  												! Max age  
    integer, parameter :: Jret=23												! Retirement age
    integer, parameter :: I=int(125d3/(shrink)) 										! Number of agents for simulation 
    integer, parameter :: T_sim=(1-solve_for_coeffs)*300+ solve_for_coeffs*3e3*no_pref+solve_for_coeffs*7e3*(1-no_pref)  ! Number of time periods for simulation, sets 5,000 when converging, 600 otherwise
    integer, parameter :: burn=100 												! burn in periods

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		VFI parameters
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
    integer, parameter :: iter_total=(1-solve_for_coeffs)+solve_for_coeffs*7 	! Sets to 1 if solve_for_coeffs==0, 5 if solve_for_coeffs==1
    integer, parameter :: iter_total2=4											! This is for the Chipeniuk et. al statistic
	real(8), parameter :: damp=(1d0-solve_for_coeffs)+0.9*solve_for_coeffs		! =1 and doesn't update coefficients if solve_for_coeffs=0, otherwise sets to 0.5
    real(8), parameter :: death=-1d10											! Value function value if no-solution 								
	real(8), parameter :: deathdefault=death/10d0								! Somehow different in default
	real(8), parameter :: tol=6d-3 												! Convergence tolerance, KMV have 5d-3
    real(8), parameter :: tol2=1d-5												! Tolerance for Chipeniuk et. al statitistic
	integer, parameter :: nD=5    												! Number of discrete VFI choices
	logical, parameter :: CKW=.False. 											! .False. means the Chipeniuk et. al statistic is not calculated
	real(8), parameter :: cmin=1.0d-8											! Minimum consumption value
	real(8), parameter :: round_size=1d3 										! What you round Chipeniuk et. al statistics to
   
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		Aggregate Risk Values
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8), parameter :: agg_theta=no_Agg+(1-no_Agg)*0.9d0 					! Aggregate income probabilities. Sets to one if there is no aggregate income risk.
	real(8), parameter :: agg_lamb=no_Agg_lend+(1-no_Agg_lend)*(1-no_MIT)&
	+(1-no_Agg_lend)*no_MIT*1d0 !0.99d0 !0.99d0 !0.99d0												! Aggregate credit condition probabilities. Sets to one if there is no aggregate income risk or MIT shocks. SEts to 0.98 if no_MIT=1
	! THey have 0.9d0 in their code instead of 0.98
	real(8), parameter :: Y2=0.965d0 											! Low state aggregate ubcine
	real(8), parameter :: Y1=Y2*no_agg+(1-no_Agg)*1.035d0						! High state aggregate income
	! For the four credit condition variables below there are three options:
	! solve_for_coeffs=0, high_nlamb=1 sets {high,low)
	! solve_for_coeffs=1, high_nlamb=1 sets (high,high)
	! solve_for_coeffs=1, high_nlamb=0 sets (low,low)
	
	real(8), parameter, dimension(2) :: lambdaLTV=(/0.95d0*(1-high_nlamb)+high_nlamb*1.1d0,0.95d0*(1-high_nlamb*solve_for_coeffs*(1-no_mIT))+1.1d0*high_nlamb*solve_for_coeffs*(1-no_MIT)  /) 	!1.1d0				! LTV
    real(8), parameter, dimension(2) :: lambdaPTI=(/0.25d0*(1-high_nlamb)+high_nlamb*0.5d0,0.25d0*(1-high_nlamb*solve_for_coeffs*(1-no_mIT))+0.5d0*high_nlamb*solve_for_coeffs*(1-no_MIT) /)	!0.5d0					! PTI
	real(8), parameter, dimension(2) :: xi=(/ 0.01d0*(1-high_nlamb)+high_nlamb*0.006d0 , 0.01d0*(1-high_nlamb*solve_for_coeffs*(1-no_MIT))+0.006d0*high_nlamb*solve_for_coeffs*(1-no_MIT) /)	! 0.006d0					! Mortgage pricing wedge 
	real(8), parameter, dimension(2) :: kappam=(/2d3*(1-high_nlamb)+high_nlamb*1.2d3 ,2d3*(1-high_nlamb*solve_for_coeffs*(1-no_MIT))+1.2d3*high_nlamb*solve_for_coeffs*(1-no_MIT)/)/(DataAvAnnualEarns/Numeraire)! 1.2d3 ! Mortgage origination cost


	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		Model parameters
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!4.070000000000000E-002  5.413100000000001E-002
    real(8), parameter :: beta=0.94d0*no_pref+(1-no_pref)*0.93d0 											! Discount factor, KMV have: 0.93d0 
    real(8), parameter :: riskfree=0.025d0 										! Risk free rate
    real(8), parameter :: iota=0.33d0 											! Mortgage wedge over risk-free rate
	!real(8), parameter,dimension(2) :: r_lend_grid=(/5.17d-2,5.36d-2/) !(/5.91d-2,5.72d-2/) !(/5.91d-2,5.72d-2/) !(/5.72d-2,5.54d-2/)
	!real(8), parameter,dimension(2) :: r_lend_grid=(/6.090000000000000d-2,6.090000000000000d-2/) !(/7.090000000000000d-2*high_rate+(1-high_rate)*6.090000000000000d-2*(1-high_nlamb)*(1-low_rate)+(1-high_rate)*high_nlamb*(1-low_rate)*6.090000000000000d-2+(1-high_rate)*high_nlamb*low_rate*4.070000000000000d-2,6.090000000000000d-2*(1-high_nlamb*solve_for_coeffs*(1-no_mIT)*(1-low_rate))+4.070000000000000d-2*high_nlamb*solve_for_coeffs*(1-no_MIT)*low_rate  /)
	!This is derived from (1+0.03)^2-1=6.09d-2
	! The mortage rate is thus (1+.33)*6.09d-2
	!main simulation is 0.95* of original, 
	!0.99d0 -running 
	!0.98d0
	!0.97d0
	!0.96d0
	real(8), parameter,dimension(2) :: r_lend_grid=(/ 6.090000000000000d-2*(1-high_nlamb)*(1-low_rate)+high_nlamb*(1-low_rate)*6.090000000000000d-2+high_nlamb*low_rate*(6.090000000000000d-2*0.97d0),6.090000000000000d-2*(1-high_nlamb*solve_for_coeffs*(1-no_mIT)*(1-low_rate))+(6.090000000000000d-2*0.95d0)*high_nlamb*solve_for_coeffs*(1-no_MIT)*low_rate  /)
	
	real(8), parameter,dimension(2) :: q_lend_grid=1d0/(1d0+r_lend_grid)
	
   real(8), parameter,dimension(2) :: r_borrow_grid=(/(0.4d0*(9.099700000000001d-2-8.099700000000001d-2)+8.099700000000001d-2)*high_rate+ (1-high_rate)*8.099700000000001d-2*(1-high_nlamb)*(1-low_rate)+(1-high_rate)*8.099700000000001d-2*high_nlamb*(1-low_rate)+(1-high_rate)*high_nlamb*low_rate*5.413100000000001d-2,8.099700000000001d-2*(1-high_nlamb*solve_for_coeffs*(1-no_mIT)*(1-low_rate))+5.413100000000001d-2*high_nlamb*solve_for_coeffs*(1-no_MIT)*low_rate  /)											!8.099700000000001d-002 						! Interest rate on borrowing, r_lend*(1d0+iota)
   
   
   !(/ 8.099700000000001d-2*(1-high_nlamb)*(1-low_rate)+8.099700000000001d-2*high_nlamb*(1-low_rate)+high_nlamb*low_rate*(8.099700000000001d-2*0.97d0),8.099700000000001d-2*(1-high_nlamb*solve_for_coeffs*(1-no_mIT)*(1-low_rate))+(8.099700000000001d-2*0.95d0)*high_nlamb*solve_for_coeffs*(1-no_MIT)*low_rate  /)											!8.099700000000001d-002 						! Interest rate on borrowing, r_lend*(1d0+iota)
    
	real(8), parameter,dimension(2) :: q_borrow_grid=1d0/(1d0+r_borrow_grid)! Interest rate on borrowing, q form	
	real(8), parameter :: tauh=0.02d0 											! Tax on housing 
    real(8), parameter :: delta=0.03d0 											! Housing depreciation
	real(8), parameter :: delta_default=0.22d0									! Housing depreciation in default
    real(8), parameter :: LbarH=0.311d0 										! Land permits
    real(8), parameter :: cost=0.07d0 											! Housing transaction cost	
	real(8), parameter :: alpha=0.6d0											! Housing supply elasticity parameter
    real(8), parameter :: elast=alpha/(1d0-alpha)								! Housing supply elasticity	
	real(8), parameter :: Mortdeduct=0.75d0										! Mortgage interest rate deduction									
	integer, parameter :: payextra=4 											! Extra years of mortgage amoritization (no clue why KMV set this to anything other than 1)
	! original rental cost: 7.500000000000000d-2
	! re-calibrated: 0.023342134d0
	! average: 0.049171
	real(8), parameter :: rental_cost=0.023342134d0 !7.500000000000000d-2 ! 0.049171d0 !0.023342134d0 !*no_pref+(1-no_pref)*6.426147610519373d-3  !		6.426147610519373d-3 !							! Rental company operating cost
    real(8), parameter, dimension(2) :: rental_grid=rental_cost*(1d0-(1d0/(1d0+r_lend_grid))*(1d0-delta))
	real(8), parameter :: gamma_b=0.8d0											! IES
	real(8), parameter :: sigma=2d0												! Risk aversion
	real(8), parameter, dimension(3) :: phi_grid=(/0.13d0*no_pref+(1-no_pref)*0.2d0,0.13d0*no_pref+(1-no_pref)*0.12d0,0.13d0*no_pref+(1-no_pref)*0.12d0/) ! Housing preference !(/0.13d0*no_pref+(1-no_pref)*0.2d0,0.12d0*no_pref+(1-no_pref)*0.13d0,0.13d0*no_pref+(1-no_pref)*0.12d0/) ! Housing preference
    real(8), parameter :: psi=100d0												! Stength of bequest motive
	real(8), parameter :: bflat=400000d0/(DataAvAnnualEarns/Numeraire) 			! Extent that bequests are luxuries
    real(8), parameter :: rhoe=0.97d0**2d0 										! Income persistence. Squaring so that this is the bi-annual persistence
    real(8), parameter :: sigma_V0=sqrt(0.18d0) 								! Income standard deviation
    real(8), parameter :: sigma_annual=0.15d0									! Initial annual earnings persistence, don't rescale but it enters the model as (0.2d0)*(1+rhoe)
    real(8), parameter :: sigmae=sigma_annual*sqrt((1d0+rhoe))					! Initial annual persistence scaled
    real(8), parameter :: tau0=0.75d0 											! Tax parameter as a share of income
    real(8), parameter :: tau1=0.151d0 											! Tax progressivity
    real(8), parameter :: rhoss=0.42d0											! Social security replacement rate 
    real(8), parameter :: utilitycost=1.015d0 !1.015d0 							! Homeownership utility
    real(8), parameter :: default_cost=0.8d0 !0.8d0 !0.9d0 									! utility cost of foreclosure 0.8d0 !This is what KMV have, foreclosure probably needs to be higher so that the foreclosure rate isn't too high

! .13 is run with a larger rental cost of 0.04 instead of approx 0.02
! .14 is run with the smaller rental cost of 0.03
	
	
end module parameters  
