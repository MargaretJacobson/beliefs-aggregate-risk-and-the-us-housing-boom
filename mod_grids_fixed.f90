 module mod_grids_fixed
    use parameters
	use mod_functions
	use mod_globals
    use mod_matlab
    use mod_Rouwenhorst
    use mod_discreteAR1
    use mod_orderpack
    implicit none 
contains
    
    subroutine sub_grids_fixed(L,H,Lsim,Bsim,R,Y,P,FYY,chi,e,countP, AggH_s,loge,Fee,iter2,V,Trans)
    
    implicit none        
	real(8), dimension(:), intent(in) :: P
    real(8), dimension(:), intent(out) :: L,H,Y ,chi,e,R,Lsim,Bsim
    real(8), dimension(:,:), intent(out) :: loge
    real(8), dimension(:,:,:), intent(out) :: Fee
    real(8), dimension(:,:), intent(out) :: Fyy
	real(8), dimension(nY,nY) :: FTheta 
	real(8), dimension(nY,nY) :: Flamb !Even if nlamb=1 still define this as nYxnY
	real(8), dimension(3,3) :: Fphi !Even if no_pref=1 still define this as 3x3
    real(8),dimension(:), intent(out) :: AggH_s
    integer, intent(in) :: countP
    integer :: iP, iM, iJ, KV, iT,iter2,iI,iE,iEE,iY,neg, ilamb, ilamb2,iphi,iphi2,iAgg
    real(8), dimension(:), intent(out) :: Trans
    real(8), dimension(:),intent(out) :: V
    double precision :: meanchi
	double precision :: lb(4)
    real(8) :: bub,blb,bneg,resid
    
 
     KV=0
	loge=0d0
	Fee=0d0

 ! Individual Loan grids. L: state, LL: choice

    
    if (nL==3) then
        L=linspace(0d0,tol,nL)
		Lsim=L
	else if (nL>3 .and. nL<22) then
		L=linspace(0d0,1d0,nL)
		Lsim=L
    
    else if (nL==22) then
        
		L(1:5)=(/0d0,3.75d-2,7.5d-2,0.1125d0,0.15d0/)
		L(6:10)=(/0.194117647058824d0,0.238235294117647d0,0.282352941176471d0,0.326470588235294d0,0.370588235294118d0/)
		L(11:15)=(/0.414705882352941d0,0.458823529411765d0,0.502941176470588d0,0.547058823529412d0,0.591176470588235d0/)
		L(16:20)=(/0.635294117647059d0,0.679411764705882d0,0.723529411764706d0,0.767647058823529d0,0.811764705882353d0 /)
		L(21:22)=(/0.855882352941176d0,0.900000000000000d0 /)
		!
	end if 
	if (nLSim>nL) then 
		Lsim(1:10)=(/0d0,1.25d-2, 2.5d-2, 3.75d-2,5.0d-2, 6.25d-2,7.5d-2,  8.749999999999999d-2,  0.1d0, 0.1125d0 /)
		Lsim(11:20)=(/0.125d0, 0.1375d0, 0.15d0,0.164705882352941d0,0.179411764705882d0, 0.194117647058824d0,0.208823529411765d0,0.223529411764706d0,0.238235294117647d0,0.252941176470588d0/)  
		Lsim(21:25)=(/0.267647058823529d0,0.282352941176471d0,0.297058823529412d0,0.311764705882353d0,0.326470588235294d0/)
		Lsim(26:30)=(/0.341176470588235d0,0.355882352941176d0,0.370588235294118d0,0.385294117647059d0,0.4d0/)
		Lsim(31:35)=(/0.414705882352941d0,0.429411764705882d0,0.444117647058823d0,0.458823529411765d0,0.473529411764706d0/)
		Lsim(36:40)=(/0.488235294117647d0,0.502941176470588d0,0.517647058823529d0,0.532352941176470d0,0.547058823529412d0/)
		Lsim(41:45)=(/0.561764705882353d0,0.576470588235294d0,0.591176470588235d0,0.605882352941176d0,0.620588235294118d0/)   
		Lsim(46:50)=(/0.635294117647059d0,0.650000000000000d0,0.664705882352941d0,0.679411764705882d0,0.694117647058823d0/)    
		Lsim(51:55)=(/0.708823529411765d0,0.723529411764706d0,0.738235294117647d0,0.752941176470588d0,0.767647058823529d0/)       
		Lsim(56:60)=(/0.782352941176470d0,0.797058823529412d0,0.811764705882353d0,0.826470588235294d0,0.841176470588235d0/)    
		Lsim(61:64)=(/0.855882352941176d0,0.870588235294118d0,0.885294117647059d0,0.9d0/)    
 
	else
	Lsim=L
	end if
    
	print*, 'L', L

  

    R=(/1.125d0,1.5d0,1.92d0/) 
	!R=H
	Print*, 'R Grids fixed', R
	
    
    H=(/1.5d0,1.92d0,2.4576d0,3.145728d0,4.02653184d0, 5.1539607552d0/) 
	
  	print*,'H', H
	blb=0d0
	bub=30d0 !22d0
	bneg=-3d0
	
	if (curve==1d0) then
		B(1:nB)=linspace(bneg,bub,nB)
	else
	
		B(nbzero:nB)=(linspace(0d0,(bub-blb)**(1d0/curve),nB-nbzero+1))**curve+blb 
	
		if (nbzero/=0) then 
		print*, 'nBzero',nbzero
			print*, (linspace(0d0,(blb-bneg)**(curve),nbzero-1))**(1d0/curve)+bneg
			B(1:nbzero-1)=(linspace(0d0,(blb-bneg)**(curve),nbzero-1))**(1d0/curve)+bneg
		end if
	end if 

	print*, 'nBzero',nbzero
	
	B(1:5)=(/-1.00000000000000d0, -0.922668688601652d0,-0.791839546184902d0,-0.628501427715763d0,-0.439673610223642d0/)
	B(6:10)=(/-0.229302554642851d0,0.000000000000000d0, 9.532509895170200d-3,  5.392401910882171d-2, 0.148597105629156d0/)
	B(11:15)=(/0.305040316645446d0,0.532883540762519d0,0.840592168440562d0,1.23581181461168d0,1.72556861148229d0/)
	B(16:20)=(/2.31639999561301d0,3.01444452204700d0,3.82550667147031d0,4.75510738013301d0,5.80852225076977d0/)
	B(21:26)=(/6.99080731505899d0,8.30683230406095d0,9.76129013265428d0,11.3587238461934d0,13.1035371587076d0, 15.0000000000000d0/)
	
  if (nBsim>26) then
  !B(27:32)=(/16.8572d0, 19.0742d0,21.4662d0,24.0392d0,26.7992d0,29.7522d0/)
	Bsim(1:6)  =(/-1.00000000000000d0, -0.974222896200551d0, -0.948445792401101d0, -0.922668688601652d0, -0.879058974462735d0, -0.835449260323819d0/)
	Bsim(7:12) =(/-0.791839546184902d0, -0.737393506695189d0, -0.682947467205476d0, -0.628501427715763d0, -0.565558821885056d0, -0.502616216054349d0/)
    Bsim(13:18)=(/-0.439673610223642d0, -0.369549925030045d0, -0.299426239836448d0, -0.229302554642851d0, -0.152868369761901d0, -7.643418488095047d-2/)
    Bsim(19:24)=(/-0.000000000000000d0, 3.177503298390066d-3, 6.355006596780133d-3,  9.532509895170200d-3, 2.432967963305404d-2, 3.912684937093787d-2/)
    Bsim(25:30)=(/5.392401910882171d-2, 8.548171461559997d-2, 0.117039410122378d0, 0.148597105629156d0, 0.200744842634586d0, 0.252892579640016d0/)
    Bsim(31:36)=(/0.305040316645446d0,  0.380988058017804d0,  0.456935799390161d0, 0.532883540762519d0, 0.635453083321866d0, 0.738022625881214d0/)
    Bsim(37:42)=(/0.840592168440562d0, 0.972332050497602d0, 1.10407193255464d0, 1.23581181461168d0, 1.39906408023522d0, 1.56231634585876d0/)
    Bsim(43:48)=(/1.72556861148229d0, 1.92251240619253d0, 2.11945620090277d0, 2.31639999561301d0, 2.54908150442434d0, 2.78176301323567d0/)
    Bsim(49:54)=(/3.01444452204700d0, 3.28479857185477d0, 3.55515262166254d0, 3.82550667147031d0, 4.13537357435787d0, 4.44524047724544d0/)
    Bsim(55:60)=(/4.75510738013301d0, 5.10624567034526d0, 5.45738396055752d0, 5.80852225076977d0, 6.20261727219951d0, 6.59671229362925d0/)
	Bsim(61:66)=(/6.99080731505899d0, 7.42948231139297d0, 7.86815730772696d0, 8.30683230406095d0, 8.79165158025873d0, 9.27647085645651d0/)
	Bsim(67:72)=(/9.76129013265428d0, 10.2937680371673d0, 10.8262459416804d0, 11.3587238461934d0, 11.9403282836982d0, 12.5219327212029d0/) 
	Bsim(73:76)=(/13.1035371587076d0, 13.7356914391384d0, 14.3678457195692d0, 15.0000000000000d0/)     
	else 
	Bsim=B
end if
Bs= (/-4.57315386520422d0,       -3.65852309216338d0,      -2.74389231912253d0,       -1.82926154608169d0,      -0.914630773040845d0, 0d0/)
		
 print*, 'B', B 
print*, 'BS',BS 
	
    ! Income grid    
    Y(1)=Y1
    Y(2)=Y2
    
    
    FTheta(1,1)=agg_theta    ! Transition matrix for income
    FTheta(1,2)=1d0-FTheta(1,1)
    FTheta(2,2)=agg_theta
    FTheta(2,1)=1d0-FTheta(2,2)
	
	!Transition matrix for credit conditions
	if (nlamb==1) then
		Flamb=0d0
		Flamb(1,1)=1
	else 
		Flamb(1,1)=agg_lamb
		Flamb(1,2)=1d0-Flamb(1,1)
		Flamb(2,2)=agg_lamb
		Flamb(2,1)=1d0-Flamb(2,2)
	end if
	
	!Transition matrix for preferences
	if (no_pref==1) then
		Fphi=0d0
		Fphi(1,1)=1d0
	
	else
		Fphi=0d0
		
		Fphi(1,3)=0.01d0*(1-no_agg_phi)+no_agg_phi*0d0
		Fphi(1,2)=0.04d0*(1-no_agg_phi)+no_agg_phi*0d0
		Fphi(1,1)=0.95d0*(1-no_agg_phi)+no_agg_phi*0d0
		
		
		
		resid=2
		Fphi(2,1)=0.85d0 !0.85d0*(1-no_agg_phi)+no_agg_phi*0d0 ! Transition to high taste
		Fphi(2,resid)=0.125d0 !0.125d0*(1-no_agg_phi)+no_agg_phi*0d0 !stay in middle taste, in the paper this is the residual
		Fphi(2,3)=1d0-Fphi(2,1)-Fphi(2,resid) !Residual, in the paper this is 0.12d0
		
		Fphi(3,3)=Fphi(1,1)
		Fphi(3,2)=Fphi(1,2)
		Fphi(3,1)=Fphi(1,3)
		
		
		
	end if 
	print*, 'Y1', FTheta(1,:)
	print*, 'Y2', Ftheta(2,:)
	print*, 'lamb1', Flamb(1,:)
	print*, 'lamb2', Flamb(2,:)
	print*, 'phi1', Fphi(1,:)
	print*, 'phi2', Fphi(2,:)
	print*, 'phi3', Fphi(3,:)
	! Creating the Kroenecker product of the three aggregate matrices
	if (nlamb==1 .and. no_pref==1) then
		FYY=FTheta
	else if (nlamb>1 .and. no_pref==1) then
		do ilamb=1,nlamb
			do ilamb2=1,nlamb
				FYY(ilamb*2-1:ilamb*2,ilamb2*2-1:ilamb2*2 )=Flamb(ilamb,ilamb2)*Ftheta
				print*, Ftheta, Flamb(ilamb,ilamb2), Flamb(ilamb,ilamb2)*Ftheta, FYY(ilamb*2-1:ilamb*2,ilamb2*2-1:ilamb2*2 ),ilamb*2-1,ilamb*2,ilamb2*2-1,ilamb2*2
			end do
		end do
	else if (nlamb>1 .and. no_pref==0) then
		do iphi=1,3
			do ilamb=1,nlamb
				do iphi2=1,3
					do ilamb2=1,nlamb
						
						FYY(ilamb*2-1+(iphi-1)*4:ilamb*2+(iphi-1)*4,ilamb2*2-1+(iphi2-1)*4:ilamb2*2+(iphi2-1)*4 )=Fphi(iphi,iphi2)*Flamb(ilamb,ilamb2)*Ftheta
						print*, 'iphi,ilamb,iphi2,ilamb2, index', iphi,ilamb,iphi2,ilamb2,'list: (',ilamb*2-1+(iphi-1)*4,':',ilamb*2+(iphi-1)*4,',',ilamb2*2-1+(iphi2-1)*4,':',ilamb2*2+(iphi2-1)*4,'),', Fphi(iphi,iphi2),Flamb(ilamb,ilamb2),FYY(ilamb*2-1+(iphi-1)*4:ilamb*2+(iphi-1)*4,ilamb2*2-1+(iphi2-1)*4:ilamb2*2+(iphi2-1)*4 )
						!print*, FYY(ilamb*2-1+(iphi-1)*4:ilamb*2+(iphi-1)*4,ilamb2*2-1+(iphi2-1)*4:ilamb2*2+(iphi2-1)*4), Fphi(iphi,iphi2),Flamb(ilamb,ilamb2)
					end do
				end do
			end do
		end do 
	end if
	do iAgg=1,nlamb*nY*nPhi
		print*,iAgg, 'Fyy',sum(FYY(iAgg,:)), FYY(iAgg,:)
	end do 
	print*, 'Fyy', FYY

   

    ! Code below to produce chi is from Kaplan, Mitman, Violante (2020)
    lb(1) = 5.080779897862075d0
    lb(2) = 0.354212100578216d0
    lb(3) = -0.006568425066786d0
    lb(4) = 0.000039013352423d0

    meanchi=0.0d0
	
    DO ij = 1,Jret-1
		chi(iJ)=lb(1)
        DO iM = 1,3
            Chi(ij) = Chi(ij) + ((dble(ij-1)*2.0d0+22.5d0)**iM)*lb(iM+1)
        END DO
        meanchi=meanchi+exp(chi(ij))/dble(Jret-1)
    END DO
    chi = chi-log(meanchi)

	chi(1)=-0.900937599378040d0
	chi(2)=-0.680596561841712d0
	chi(3)=-0.489863071301359d0
	chi(4)=-0.326864486780767d0
	chi(5)=-0.189728167303716d0
	chi(6)=-7.658147189399322d-2
	chi(7)=1.444824042462223d-2 
	chi(8)=8.523361062834489d-2
	chi(9)=0.137647279693390d0
	chi(10)=0.173561888595977d0
	chi(11)=0.194850078312319d0
	chi(12)=0.203384489818637d0
	chi(13)=0.201037764091144d0
	chi(14)=0.189682542106059d0
	chi(15)=0.171191464839596d0
	chi(16)=0.147437173267976d0
	chi(17)=0.120292308367413d0
	chi(18)=9.162951111412099d-2  
	chi(19)=6.332142248432149d-2  
	chi(20)=3.724068345422715d-2  
	chi(21)=1.525993500005818d-2 
	chi(22)=-7.481819019723446d-4
 !   print*, 'chi', chi

    

    !   McClement's Scale from Kaplan, Mitman, Violante (2020)

    e((/1,2,3,4,5/))=(/1.07369090000000d0,1.09825430000000d0,1.13220340000000d0,1.16617120000000d0,1.20122030000000d0/)
    e((/6,7,8,9,10/))=(/1.24891980000000d0,1.29554010000000d0,1.32435800000000d0,1.34751300000000d0,1.36606180000000d0/)
    e((/11,12,13,14,15/))=(/1.37073070000000d0,1.36245540000000d0,1.33539560000000d0,1.30071250000000d0,1.26504640000000d0/)
    e((/16,17,18,19,20/))=(/1.22184160000000d0,1.19089120000000d0,1.15449710000000d0,1.13195880000000d0,1.11513750000000d0/)
    e((/21,22/))=(/ 1.10068620000000d0, 1.09124640000000d0/)
	e(Jret:)=1d0
    !Below is waht they have in the paper, uncomment to match their results
	!e(Jret:)=e(1:8)

    if (KMV==0) then
		do iJ=1,T
			if (iJ==1) then
				!Code below uses timting
				V(iJ)=sigma_V0**2d0+rhoe**2d0*sigmae**2d0
				call sub_rouwenhorst(loge(iJ,:),Fee(iJ,:,:),rhoe,sigmae,sigma_V0,sqrt(V(iJ)))
				!print*, iJ, loge(iJ,:)
				!print*, iJ, Fee(iJ,:,:)
			else if (iJ>1 .and. iJ<Jret) then

				V(iJ)=rhoe**2d0*V(iJ-1)+sigmae**2d0
				call sub_rouwenhorst(loge(iJ,:),Fee(iJ,:,:),rhoe,sigmae,sqrt(V(iJ-1)),sqrt(V(iJ)))
				!  print*, iJ, loge(iJ,:)
				!print*, iJ, Fee(iJ,:,:)
			else if (iJ>=Jret) then
				loge(iJ,:)=loge(Jret-1,:)
				Fee(iJ,:,:)=Fee(Jret-1,:,:)
				V(iJ)=V(Jret-1)
			end if
		end do
   ! print*, 'Maggie Initial variance (V0,sigmae)', sigma_V0, sigmae
    !print*, 'Maggie Variance of earnings', V
	
    
		do IJ=Jret,T
			Fee(iJ,:,:)=Fee(Jret-1,:,:)
		end do
	
	else
		
		
		! Taking values directly from their output
			loge(1,1)=-0.882986561359690d0
			loge(1,2)=-0.441493280679845d0
			loge(1,4)=loge(1,2)*(-1d0)
			loge(1,5)=loge(1,1)*(-1d0)
			
			loge(2,1)=-0.976251283262193d0
			loge(2,2)=-0.488125641631097d0
			loge(2,4)=loge(2,2)*(-1d0)
			loge(2,5)=loge(2,1)*(-1d0)
			
			loge(3,1)=-1.05193987122568d0
			loge(3,2)=-0.525969935612838d0
			loge(3,4)=loge(3,2)*(-1d0)
			loge(3,5)=loge(3,1)*(-1d0)
			
			loge(4,1)=-1.11466569449408d0
			loge(4,2)=-0.557332847247039d0
			loge(4,4)=loge(4,2)*(-1d0)
			loge(4,5)=loge(4,1)*(-1d0)
			
			loge(5,1)=-1.16738715867568d0
			loge(5,2)=-0.583693579337839d0
			loge(5,4)=loge(5,2)*(-1d0)
			loge(5,5)=loge(5,1)*(-1d0)
			
			loge(6,1)=-1.21214898449248d0
			loge(6,2)=-0.606074492246242d0
			loge(6,4)=loge(6,2)*(-1d0)
			loge(6,5)=loge(6,1)*(-1d0)
			
			loge(7,1)=-1.25043984537047d0
			loge(7,2)=-0.625219922685233d0
			loge(7,4)=loge(7,2)*(-1d0)
			loge(7,5)=loge(7,1)*(-1d0)
			
			loge(8,1)=-1.28338543675187d0
			loge(8,2)=-0.641692718375937d0
			loge(8,4)=loge(8,2)*(-1d0)
			loge(8,5)=loge(8,1)*(-1d0)
			
			loge(9,1)=-1.31186164747750d0
			loge(9,2)=-0.655930823738749d0
			loge(9,4)=loge(9,2)*(-1d0)
			loge(9,5)=loge(9,1)*(-1d0)
			
			loge(10,1)=-1.33656522554242d0
			loge(10,2)=-0.668282612771211d0
			loge(10,4)=loge(10,2)*(-1d0)
			loge(10,5)=loge(10,1)*(-1d0)
			
			loge(11,1)=-1.35806017286097d0
			loge(11,2)=-0.679030086430484d0
			loge(11,4)=loge(11,2)*(-1d0)
			loge(11,5)=loge(11,1)*(-1d0)
			
			loge(12,1)=-1.37680947422970d0
			loge(12,2)=-0.688404737114849d0
			loge(12,4)=loge(12,2)*(-1d0)
			loge(12,5)=loge(12,1)*(-1d0)
			
			loge(13,1)=-1.39319754341724d0
			loge(13,2)=-0.696598771708619d0
			loge(13,4)=loge(13,2)*(-1d0)
			loge(13,5)=loge(13,1)*(-1d0)
			
			loge(14,1)=-1.40754656090061d0
			loge(14,2)=-0.703773280450303d0
			loge(14,4)=loge(14,2)*(-1d0)
			loge(14,5)=loge(14,1)*(-1d0)
			
			loge(15,1)=-1.42012865725093d0
			loge(15,2)=-0.710064328625467d0
			loge(15,4)=loge(15,2)*(-1d0)
			loge(15,5)=loge(15,1)*(-1d0)
			
			loge(16,1)=-1.43117518961548d0
			loge(16,2)=-0.715587594807741d0
			loge(16,4)=loge(16,2)*(-1d0)
			loge(16,5)=loge(16,1)*(-1d0)
			
			loge(17,1)=-1.44088393318582d0
			loge(17,2)=-0.720441966592909d0
			loge(17,4)=loge(17,2)*(-1d0)
			loge(17,5)=loge(17,1)*(-1d0)
			
			loge(18,1)=-1.44942474428888d0
			loge(18,2)=-0.724712372144442d0
			loge(18,4)=loge(18,2)*(-1d0)
			loge(18,5)=loge(18,1)*(-1d0)
			
			loge(19,1)=-1.45694408136844d0
			loge(19,2)=-0.728472040684220d0
			loge(19,4)=loge(19,2)*(-1d0)
			loge(19,5)=loge(19,1)*(-1d0)
			
			loge(20,1)=-1.46356865776132d0
			loge(20,2)=-0.731784328880660d0
			loge(20,4)=loge(20,2)*(-1d0)
			loge(20,5)=loge(20,1)*(-1d0)
			
			loge(21,1)=-1.46940842429836d0
			loge(21,2)=-0.734704212149180d0
			loge(21,4)=loge(21,2)*(-1d0)
			loge(21,5)=loge(21,1)*(-1d0)
			
			loge(22,1)=-1.47455902742972d0
			loge(22,2)=-0.737279513714859d0
			loge(22,4)=loge(22,2)*(-1d0)
			loge(22,5)=loge(22,1)*(-1d0)

		
		!do iJ=1,T
		!	do iE=1,nE
		!		do iY=1,nY
		!			if (iJ<Jret) then
		!				print*, iJ, chi(iJ),loge(iJ,iE), Y(iY), income(chi(iJ),loge(iJ,iE),Y(iY)),FnTax(income(chi(iJ),loge(iJ,iE),Y(iY))),income(chi(iJ),loge(iJ,iE),Y(iY))-FnTax(income(chi(iJ),loge(iJ,iE),Y(iY)))
		!			else
		!				print*, iJ, chi(Jret-1),loge(iJ,iE), Y(iY), income_ret(chi(iJ),loge(iJ,iE),mean(Y)), FnTax(income_ret(chi(iJ),loge(iJ,iE),mean(Y))), income_ret(chi(iJ),loge(iJ,iE),mean(Y))-FnTax(income_ret(chi(iJ),loge(iJ,iE),mean(Y)))
		!			end if 
		!		end do
		!	end do
		!end do
		
		
		! input transition matrices
		OPEN(3, FILE = './KMV_Input/zytrans_Maggie.txt',status='old',action='read')
		do iJ=1,Jret-1
				read(3,*) Fee(iJ,1:nE,1:nE)
		end do
		close(3)
		
	
		
		do iJ=Jret,T
		!do iJ=1,T
      Fee(iJ,:,:)=0d0
			do iE=1,nE
			Fee(iJ,iE,iE)=1d0
			end do
		end do
		
		do iJ=1,T
			do iE=1,nE
		!	print*, 'Fee(iJ,iE,:)', iJ,iE,Fee(iJ,iE,:)
      end do
			if (iJ<Jret) then   
		!	print*, 'Income iJ:', iJ, income(chi(iJ),loge(iJ,1),Y(2)),income(chi(iJ),loge(iJ,2),Y(2)),income(chi(iJ),loge(iJ,3),Y(2)),income(chi(iJ),loge(iJ,4),Y(2)),income(chi(iJ),loge(iJ,5),Y(2))
			else
		!	print*, 'Income iJ:', iJ, income_ret(chi(Jret-1),1,Y(1)), income_ret(chi(Jret-1),2,Y(1)), income_ret(chi(Jret-1),3,Y(1)), income_ret(chi(Jret-1),4,Y(1)), income_ret(chi(Jret-1),5,Y(1))
			end if 

			!end do
		end do
	end if





	
        print*, 'size of p grid:',countp
        print*, 'p',p
        
        AggH_s=(alpha*P)**(elast)*LbarH !+2d0*(1d0-delta)
        !AggH_s=Lbarh*0.1d0 !delta
        !print*,'AggH_s', AggH_s
        
            print*, 'parameters'
    print*, 'beta', beta

    print*, 'delta', delta
    print*, 'sigmae', sigmae
    print*, 'tauh', tauh
	print*, 'Y', Y
	print*, 'Credit: lambdaLTV, lambdaPTI, xi, kappam', lambdaLTV, lambdaPTI, xi, kappam
    
	! Transfer term in consumption
	 if (countP==13) then
	Trans(1:5)=(/2.854075402768506d-003, 4.195970408950399d-003, 5.858848718492038d-003,7.864906977237622d-003,1.023498750346541d-002/)
	Trans(6:10)=(/1.298880000000000d-002,1.614508857052269d-002, 1.972176127472518d-002,2.373599303861534d-002,2.820430902065911d-002/)
	Trans(11:13)=(/3.314265327033467d-002,3.856644607659063d-002,4.449063245604891d-002/)  
	else
	 Trans=linspace(2.854075402768506d-003,4.449063245604891d-002,countP)
	end if 
	
	print*, 'p',p
	print*, 'L end sub routine', L
    end subroutine sub_grids_fixed
end module mod_grids_fixed
