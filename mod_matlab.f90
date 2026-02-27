! Question: 1.1 Problem Set 1

! By Margaret Jacobson, Indiana University marmjaco@indiana.edu
! For E 724: Computational Macroeconomics taught by Grey Gordon in
!  Spring 2016 at Indiana University

! To write this code, I consulted lecture notes by Grey Gordon
! entitled, "Basics", "Arrays", and  "Functions, Subroutines,
! and Modules," accessed via Grey Gordon's website
! https://sites.google.com/site/greygordon/teaching
! and through the Canvas module for E724

! This file contains the module mod_matlab
! mod_matlab contains the functions linspace, normcdf that do the 
! same thing as these functions in matlab. mod_matlab also contains
! the paramter pi


module mod_matlab
    implicit none
    real(8), parameter :: pi=3.14159265358979323846d0

    

    interface mldivide
        module procedure mldivide_full, mldivide_band, mldivide_full_mat
    end interface

contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !		Linspace Function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function linspace(a,b,n) result(x)
        implicit none
        real(8), intent(in) :: a, b
        integer, intent(in) :: n
        real(8) :: x(n)
        integer :: ii
    	real(8) :: interval

        x(1) = a
        x(n) = b
        interval=(b-a)/(n-1)
        do ii = 2,n-1
            x(ii)=x(ii-1)+interval
        end do
    end function linspace

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Mean function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mean(a) result(x)
        implicit none
        real(8), intent(in) :: a(:)         
        real(8) :: x
                
        x=(1d0/size(a))*sum(a)
        
    end function mean

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Variance function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function var(a) result(x)
        implicit none
        real(8), intent(in) :: a(:)
        real(8) :: x,m
        m=mean(a)        
        x=mean((a-m)**2)     
    end function var

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Covariance function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function cov(a,b) result(x)
        implicit none
        real(8), intent(in) :: a(:), b(:)
        real(8) :: x
        if (size(a)==size(b)) then
            x=(1d0/size(a))*sum((a-mean(a))*(b-mean(b)))
        else 
            stop 'Arrays not same size'
        end if
    end function cov

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Correlation Coefficient
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function corr(a,b) result(x)
        implicit none
        real(8), intent(in) :: a(:), b(:)
        real(8) :: x 
        
        x=cov(a,b)/(sqrt(var(a))*sqrt(var(b)))
    end function corr
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   Percentile Function, Nearest Rank
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function percentile(a,p) result(x)
        ! Sort array first
        implicit none
        real(8), intent(in) :: a(:)
        integer, intent(in) :: p
        !integer :: io
        real(8) :: x,io
        io=ceiling((p/100d0)*size(a))  !This rounds up a negative integer
		!print*, 'io', io
		!print*, a(int(io))
        x=a(int(io))
		!print*, x
    end function percentile
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !	       NormCDF Function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Gaussian CDF is given as:
    ! z=(x-mean)/sqrt(var*2)
    ! Phi(z)=1/2+1/2*erf(z)
  
 
    function normcdf(x,var,mean) result(F)
        implicit none
        real(8), intent(in) :: x, var, mean
        real(8) :: erf, F,z
         
        z=(x-mean)/(sqrt(var*2d0))

        F=0.5d0+0.5d0*erf(z)

    end function normcdf
    
 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !	       Multivariate pdf
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !function norm_mlt_pdf(x,var,mu) result (F)
    !  implicit none
    !  real(8), dimension(:), intent(in):: x,mu
    !  real(8), dimension(:,:), intent(in) :: var
    !  real(8) :: n,det_var
    !  real(8), dimension(:), intent(out) :: F
    !  n=size(x)
      
    !  det_var=det(var)
      
    !  F=(1d0/sqrt((2*pi)**n)*det_var)*exp(0.5d0*()*())
      
    
    !end function norm_mlt_pdf

!   http://www.techchug.com/articles/tutorials/2017/09/Generating-Random-Numbers-With-a-Normal-Distribution-in-Fortran
!     uses the marsaglia polar method to create pseudo-random number pairs that
! have a normal distribution. The function saves one of the random numbers from
! the generated pair as a spare to be returned for the next call to the function.
! Returns a real scalar
    function norm_rand(mean, std_dev)
            
        real(8) :: norm_rand
        real(8), intent(in) :: mean, std_dev
        real(8) :: x, y, r
        real(8), save :: spare
        logical, save :: has_spare
        ! use a spare saved from a previous run if one exists
        if (has_spare) then
            has_spare = .FALSE.
            norm_rand = mean + (std_dev * spare)
        return
        else
            r = 1.0
            do while ( r >= 1.0 )
                ! generate random number pair between 0 and 1
                call random_number(x)
                call random_number(y)
                ! normalise random numbers to be in square of side-length = R
                x = (x * 2.0) - 1.0
                y = (y * 2.0) - 1.0
                r = x*x + y*y
            end do

            ! calculate the co-efficient to multiply random numbers x and y
            ! by to achieve normal distribution
            r = sqrt((-2.0 * log(r)) / r)

            norm_rand = mean + (std_dev * x * r)
            spare = y * r
            has_spare = .TRUE.
            return
        end if
    end function norm_rand
    
    ! ALl code below here is from Grey Gordon

   function mldivide_full_mat(A,b) result(x)
    implicit none
    real(8), dimension(1:,1:), intent(in) :: A
    real(8), dimension(1:,1:), intent(in) :: B
    real(8), dimension(1:size(A,2),1:size(b,2)) :: x
    !local
    real(8), dimension(:,:), allocatable :: A_local, b_local
    integer, dimension(:), allocatable :: ipiv
    integer :: info

    allocate(A_local(size(A,1),size(A,2)), b_local(size(b,1),size(b,2)), ipiv(size(A,1)))
    !Check inputs
    if (size(A,1)/=size(A,2)) then
        write(*,*) shape(A)
        STOP 'mldivide: A not square.  Cant handle this currently.'
    end if

    ! Copy A and b over because these are overwritten
    A_local = A
    b_local = b

    call dgesv(size(A_local,1), size(b_local,2), A_local, size(A_local,1), ipiv, b_local, size(b_local,1), info )
    x = b_local
     print*, X
    ! Check for failure
    select case (info)
        case(0)
        case (:-1)
            write(*,*) 'mldivide: the ', -info, ' parameter had an illegal value'
        case(1:)
            write(*,*) 'mldivide: U(i,i) for i=', info, ' is exactly zero.  The factorization ', &
                    'has been completed, but the factor U is exactly singular, so the solution ', &
                    'could not be computed.'
    end select

    deallocate(A_local, b_local, ipiv)
    
end function mldivide_full_mat

! x = mldivide(A,b)
! mldivide computes an "x" s.t. Ax = b if possible. 
! I am not sure what happens when this is not possible.
function mldivide_full(A,b) result(x)
    implicit none
    real(8), dimension(1:,1:), intent(in) :: A
    real(8), dimension(1:), intent(in) :: b
    real(8), dimension(1:size(A,2)) :: x
    !local
    real(8), dimension(:,:), allocatable :: A_local, b_local
    integer, dimension(:), allocatable :: ipiv
    integer :: info

    allocate(A_local(size(A,1),size(A,2)), b_local(size(b),1), ipiv(size(A,1)))

    !Check inputs
    if (size(b)/=size(A,1)) then
        write(*,*) size(b), size(A,1)
        STOP 'mldivide: b wrong size'
    end if
    if (size(A,1)/=size(A,2)) then
        write(*,*) shape(A)
        STOP 'mldivide: A not square.  Cant handle this currently.'
    end if

    ! Copy A and b over because these are overwritten
    A_local = A
    b_local(:,1) = b

    call dgesv(size(A_local,1), 1, A_local, size(A_local,1), ipiv, b_local, size(b_local,1), info )
    x = b_local(:,1)
    print*, X
    ! Check for failure
    select case (info)
        case(0)
        case (:-1)
            write(*,*) 'mldivide: the ', -info, ' parameter had an illegal value'
        case(1:)
            write(*,*) 'mldivide: U(i,i) for i=', info, ' is exactly zero.  The factorization ', &
                    'has been completed, but the factor U is exactly singular, so the solution ', &
                    'could not be computed.'
    end select

    deallocate(A_local, b_local, ipiv)
    
end function mldivide_full

! mldivide_band
! Solves for x in Ax=b when A is a banded matrix with kl subdiagonals, ku superdiagonals, 
! with the one main diagonal (kl>=0,ku>=0). For user convenience, I am not economizing
! on storage here, so A is a full square matrix (the banded matrix is created internally)
!
function mldivide_band(A,b,kl,ku) result(x)
    implicit none
    real(8), dimension(1:,1:), intent(in) :: A
    real(8), dimension(1:), intent(in) :: b
    real(8), dimension(1:size(A,2)) :: x
    integer, intent(in) :: kl,ku
    !local
    real(8), dimension(:,:), allocatable :: Ab, btmp ! banded A, tmp b
    integer, dimension(:), allocatable :: ipiv
    integer :: m,n,ldab,ldb,info,nrhs, i,j

    !Check inputs
    if (size(b)/=size(A,1)) then
        write(*,*) size(b), size(A,1)
        STOP 'mldivide_band: b wrong size'
    end if
    if (size(A,1)/=size(A,2)) then
        write(*,*) shape(A)
        STOP 'mldivide_band: A not square.  Cant handle this currently.'
    end if
    if (ku<0 .or. kl<0) then
        STOP 'mldivide_band: ku and kl must be nonnegative'
    end if

    ! Begin
    m = size(A,1) ! Also is size(A,2)
    n = size(A,2) ! Also is size(A,1)
    ldab = 2*kl + ku + 1 ! The extra kl is used for storage
    ldb = n ! Rows in b
    nrhs = 1 ! Columns of b

    allocate(Ab(ldab,n), btmp(n,nrhs), ipiv(n))

    ! Convert A into banded form
    do j = 1,n
        do i = max(1,j-ku),min(m,j+kl) 
            ab(kl + ku + 1 + i - j, j) = A(i,j) ! storage scheme in lapack 
        end do
    end do

    ! Copy over b as it will be overwritten
    btmp(:,1) = b

    ! Compute the lu factorization
    call dgbtrf( m, n, kl, ku, ab, ldab, ipiv, info )

    ! Check for failure of lu factorization
    select case (info)
        case(0)
        case (:-1)
            write(*,*) 'mldivide_band: in LU factorizing, the ', -info, ' parameter had an illegal value'
        case(1:)
            write(*,*) 'mldivide_band: in LU factorizing, U(i,i) for i=', info, ' is exactly zero.  The factorization ', &
                    'has been completed, but the factor U is exactly singular, so the solution ', &
                    'could not be computed.'
    end select

    ! Compute the solution using the factorized A
    call dgbtrs( 'N', n, kl, ku, nrhs, ab, ldab, ipiv, btmp, ldb, info )

    ! Copy out the solution
    x = btmp(:,1)
    
    ! Check for failure
    select case (info)
        case(0)
        case (:-1)
            write(*,*) 'mldivide_band: in solving, the ', -info, ' parameter had an illegal value'
    end select

    deallocate(Ab, btmp, ipiv)
    
end function mldivide_band

!_____________________________________________________
!          Subroutines
! RnadomDiscrete, RandomDiscrete1,DiscreteDist, DiscreteDist1
! Are taken directly from Kaplan and Violante (2014)'s code
!_____________________________________________________
SUBROUTINE RandomDiscrete(Nout,Xout,Nin,Pin)
!generates Nout random draws from the integers 1 to Nin
!using probabilities in Pin
IMPLICIT NONE
INTEGER, INTENT(in)			:: Nout,Nin
INTEGER,INTENT(out)		:: Xout(:)
REAL(8), INTENT(in)  :: Pin(:)
INTEGER			::i1,i2
REAL(8)      :: lran(Nout)

!IF(sum(Pin) .ne. 1.0) write(*,*) 'error in RandomDiscrete: Pin doesnt sum to 1.0'

CALL RANDOM_NUMBER(lran)

Xout(:) = 0
DO i1 = 1,Nout
    IF ( lran(i1) .le. Pin(1) ) THEN
        Xout(i1) = 1
    ELSE
        i2 = 2
        DO WHILE (i2 .le. Nin)
            IF ( (lran(i1) .le. SUM(Pin(1:i2)) ).and. (lran(i1) > SUM(Pin(1:i2-1)) ) ) THEN
                Xout(i1) = i2
                i2 = Nin+1
            ELSE
                i2 = i2+1
            END IF
        END DO
    END IF
END DO

END SUBROUTINE RandomDiscrete

!--------------------------------------------------------------
SUBROUTINE RandomDiscrete1(Xout,Nin,Pin)
!generates Nout random draws from the integers 1 to Nin
!using probabilities in Pin
IMPLICIT NONE
INTEGER, INTENT(in)			:: Nin
INTEGER,INTENT(out)		:: Xout
REAL(8), INTENT(in)  ::Pin(:)
INTEGER			::i2
REAL(8)      :: lran

!IF(sum(Pin) .ne. 1.0) write(*,*) 'error in RandomDiscrete: Pin doesnt sum to 1.0'

CALL RANDOM_NUMBER(lran)

Xout = 0
IF ( lran .le. Pin(1) ) THEN
    Xout = 1
ELSE
    i2 = 2
    DO WHILE (i2 .le. Nin)
        IF ( (lran .le. SUM(Pin(1:i2)) ).and. (lran > SUM(Pin(1:i2-1)) ) ) THEN
            Xout = i2
            i2 = Nin+1
        ELSE
            i2 = i2+1
        END IF
    END DO
END IF


END SUBROUTINE RandomDiscrete1

!----------------------------------------------
SUBROUTINE DiscreteDist(Nout,Xout,Nin,Pin,lran)
!takes in random numbers and a prob dist over the integers 1 to Nin
!and returns the corresponding values
!Created 11/22/09 - ERS

INTEGER,    INTENT(out) :: Xout(:)
INTEGER,    INTENT(in)  :: Nin,Nout
REAL(8),    INTENT(in)  :: Pin(:)
REAL(8), INTENT(in)  :: lran(:)
INTEGER			        :: i1,i2

Xout(:) = 0
DO i1 = 1,Nout
    IF ( lran(i1) .le. Pin(1) ) THEN
        Xout(i1) = 1
    ELSE
        i2 = 2
        DO WHILE (i2 .le. Nin)
            IF ( (lran(i1) .le. SUM(Pin(1:i2)) ).and. (lran(i1) > SUM(Pin(1:i2-1)) ) ) THEN
                Xout(i1) = i2
                i2 = Nin+1
            ELSE
                i2 = i2+1
            END IF
        END DO
    END IF
END DO

END SUBROUTINE DiscreteDist

!-----------------------------------------------
SUBROUTINE DiscreteDist1(Xout,Nin,Pin,lran)
!takes in a random number and a prob dist over the integers 1 to Nin
!and returns the corresponding value

INTEGER, INTENT(in)     :: Nin
INTEGER,INTENT(out)     :: Xout
REAL(8), INTENT(in)     :: Pin(:),lran
INTEGER                 :: i2

Xout = 0
IF ( lran .le. Pin(1) ) THEN
   Xout = 1
ELSE
   i2 = 2
   DO WHILE (i2 .le. Nin)
       IF ( (lran .le. SUM(Pin(1:i2)) ).and. (lran > SUM(Pin(1:i2-1)) ) ) THEN
           Xout = i2
           i2 = Nin+1
       ELSE
           i2 = i2+1
       END IF
   END DO
END IF

END SUBROUTINE DiscreteDist1


! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(array)
    real(8), intent(inout)::array(:)
    real(8) :: temp,pivot
    integer :: i,j,last,left,right

    last=size(array)

    if (last.lt.50) then ! use insertion sort on small arrays
       do i=2,last
          temp=array(i)
          do j=i-1,1,-1
             if (array(j).le.temp) exit
             array(j+1)=array(j)
          enddo
          array(j+1)=temp
       enddo
       return
    endif
    ! find median of three pivot
    ! and place sentinels at first and last elements
    temp=array(last/2)
    array(last/2)=array(2)
    if (temp.gt.array(last)) then
       array(2)=array(last)
       array(last)=temp
    else
       array(2)=temp
    endif
    if (array(1).gt.array(last)) then
       temp=array(1)
       array(1)=array(last)
       array(last)=temp
    endif
    if (array(1).gt.array(2)) then
       temp=array(1)
       array(1)=array(2)
       array(2)=temp
    endif
    pivot=array(2)

    left=3
    right=last-1
    do
       do while(array(left).lt.pivot)
          left=left+1
       enddo
       do while(array(right).gt.pivot)
          right=right-1
       enddo
       if (left.ge.right) exit
       temp=array(left)
       array(left)=array(right)
       array(right)=temp
       left=left+1
       right=right-1
    enddo
    if (left.eq.right) left=left+1
    call quicksort(array(1:left-1))
    call quicksort(array(left:))

  end subroutine quicksort


! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  real(8), dimension(:,:), intent(in) :: A
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv


end module mod_matlab





