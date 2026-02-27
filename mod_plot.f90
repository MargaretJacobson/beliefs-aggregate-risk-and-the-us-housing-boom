! mod_plot
!
! Some plotting tools for Fortran. They rely on having gnuplot installed on the 
! system path. It has only been tested on mac and linux but should work on windows.
!
! If using gfortran, one must comment out the line "use ifport".
!
! If using Intel fortran, one can use this code as is. It must be linked with 
! the openmp libraries since it uses omp_get_thread_num() to be thread safe. 
! However, it is only thread safe if the system call to gnuplot is thread safe.
! 
! Grey Gordon 2010-2013
! Soli Deo Gloria
!
! You may use and modify this code freely for non-commercial use as long as this notice
! remains with the code.
! 
module mod_plot
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! If using intel fortran, must include this module
    ! If using gnu fortran, must comment this line out
    use ifport
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    private 
    public :: hist,plot
    

    ! Overload procedures
    interface hist
        module procedure hist_x, hist_xn
    end interface hist
    interface plot
        module procedure plot_y_fine,plot_xy_fine,plot_y2_fine,plot_xy2_fine
    end interface plot

contains

! Very basic histogram plotting to std out
subroutine hist_x(x)
    implicit none
    real(8), dimension(1:), intent(in) :: x
    call hist_xn(x,40)
end subroutine
subroutine hist_xn(x,n)
    implicit none
    real(8), dimension(1:), intent(in) :: x
    integer, intent(in) :: n
    integer, parameter :: max_xpixels = 40,ypixels=20
    character(len=min(n,max_xpixels)) :: horizontal_line
    integer :: i, y_ind, x_ind, xpixels
    integer :: ylocmax
    real(8), allocatable :: x_axis(:), y_axis(:) 
    logical, allocatable :: xy_filled(:,:) 
    integer, allocatable :: x_count(:) 
    
    allocate(x_axis(min(n,max_xpixels)),y_axis(1:ypixels),xy_filled(min(n,max_xpixels),ypixels),x_count(min(n,max_xpixels)-1))

    xpixels = min(n,max_xpixels)

    !Generate xaxis/bins
    x_axis = linspace(minval(x),maxval(x),xpixels)

    !Create a count for each bin
    ! This is a very slow method, but speed is not essential here yet. 
    x_count(1) = count(x<=x_axis(1))
    do i = 2,xpixels-2
        x_count(i) = count(x_axis(i)<=x .and. x<x_axis(i+1))
    end do
    x_count(xpixels-1) = count(x>=x_axis(xpixels))

    !Generate yaxis
    y_axis = linspace(0.d0,real(maxval(x_count),8),ypixels)

    !Fill in the plot
    xy_filled = .FALSE.
    do i = 1,size(x_count)
        ylocmax = int(real(x_count(i),8)/y_axis(ypixels)*real(ypixels,8))
        xy_filled(i,1:ylocmax) = .TRUE.
    end do

    !Draw plot to screen
    do y_ind = ypixels,1,-1
        do x_ind = 1,xpixels
            if (xy_filled(x_ind,y_ind)) then
                horizontal_line(x_ind:x_ind) = '.'
            else 
                horizontal_line(x_ind:x_ind) = ' '
            end if
        end do
        write(*,'(f10.5,A,A)') y_axis(y_ind),'|',horizontal_line
    end do
    do x_ind = 1,xpixels
        horizontal_line(x_ind:x_ind) = '-'
    end do
    write(*,*) '          ',horizontal_line
    write(*,*) 'x_range: ',x_axis(1),',',x_axis(xpixels)

    deallocate(x_axis,y_axis,xy_filled,x_count)

end subroutine

! Line plot routines
subroutine plot_y_fine(y,title)
    implicit none
    real(8), intent(in) :: y(1:)
    character(len=*), intent(in), optional :: title

    if (present(title)) call plot_xy_fine(real(colon(1,size(y)),8),y,title)
    if (.not. present(title)) call plot_xy_fine(real(colon(1,size(y)),8),y)

end subroutine
subroutine plot_xy_fine(x,y,title,xlabel,ylabel,extra)
    implicit none
    real(8), intent(in) :: x(:), y(:)
    character(len=*), intent(in), optional :: title,xlabel,ylabel,extra
    ! local
    character(len=1000) :: t,xl,yl,ex
    logical :: useDefArg(4)

    useDefArg(1) = present(title)
    useDefArg(2) = present(xlabel)
    useDefArg(3) = present(ylabel)
    useDefArg(4) = present(extra)

    t=' '; if (present(title)) t=title
    xl=' '; if (present(xlabel)) xl=xlabel 
    yl=' '; if (present(ylabel)) yl=ylabel
    ex=' '; if (present(extra)) ex=extra

    call plot_xy2_fine(x,reshape(y,[size(y),1]),title,xlabel,ylabel,extra,useDefArg)

end subroutine plot_xy_fine
subroutine plot_y2_fine(y,title,xlabel,ylabel,extra)
    implicit none
    real(8), intent(in) :: y(:,:)
    character(len=*), intent(in), optional :: title,xlabel,ylabel,extra
    ! local
    real(8), allocatable :: x(:)
    character(len=1000) :: t,xl,yl,ex
    logical :: useDefArg(4)

    useDefArg(1) = present(title)
    useDefArg(2) = present(xlabel)
    useDefArg(3) = present(ylabel)
    useDefArg(4) = present(extra)

    t=' '; if (present(title)) t=title
    xl=' '; if (present(xlabel)) xl=xlabel 
    yl=' '; if (present(ylabel)) yl=ylabel
    ex=' '; if (present(extra)) ex=extra

    allocate(x(size(y,1)))
    x = colon(1,size(y,1))

    call plot_xy2_fine(x,y,title,xlabel,ylabel,extra,useDefArg)

end subroutine plot_y2_fine 

! plot_x_y2_fine plots x,y data using a (not compliant Fortran standard) system call to gnuplot.
! Many options are available and for the most part self-descriptive. Temporary files are saved 
! to the present working directory as plotdataX.dat and plotdataX.cfg. 
!
! Arguments are mostly obvious. The unobvious ones are 
! extra 
! -> allows the user to pass any commands to gnuplot before the data is plotted.
! -> ex: "set style data scatter" will produce a scatter plot instead of the default line plot.
! 
! useDefArg 
! -> should not be called by the user (internal mechanism)
! 
subroutine plot_xy2_fine(x,y,title,xlabel,ylabel,extra,useDefArg)
    use omp_lib
    implicit none
    real(8), intent(in) :: x(:), y(:,:)
    character(len=*), intent(in), optional :: title,xlabel,ylabel,extra
    logical, intent(in), optional :: useDefArg(4) ! Logical specifying which Default arguments
    !local
    integer :: i,j,ni,nj, unitno
    logical :: isPresent_xlabel, isPresent_ylabel, isPresent_title, isPresent_extra
    character(len=*), parameter :: newline = ', \'
    character(len=3) :: tmpchar
    character(len=5) :: jp1AsString,jAsString,tmpChar5
    character(len=1000) :: dls(size(y,2)), ms(size(y,2)),cls(size(y,2)) ! data labels and methods
    character(len=100) :: configfile, datafile 
    integer(4) :: ierr

    ! WARNING: this routine may not be thread safe. It depends on whether gnuplot is. 

    ! As an extra precaution, I try to make sure the threads write to unique files
    write(tmpChar5,'(i0)') omp_get_thread_num()        
    datafile = 'plotdata'//trim(tmpChar5)//'.dat'
    configfile = 'plotdata'//trim(tmpChar5)//'.cfg'
    unitno = 726+omp_get_thread_num()

    if (present(useDefArg)) then
        isPresent_title = useDefArg(1)
        isPresent_xlabel = useDefArg(2)
        isPresent_ylabel = useDefArg(3)
        isPresent_extra = useDefArg(4)
    else
        isPresent_title = present(title) 
        isPresent_xlabel = present(xlabel)
        isPresent_ylabel = present(ylabel)
        isPresent_extra = present(extra)
    end if

    ! Check inputs
    ni = size(x)
    nj = size(y,2)
    if (size(y,1)/=ni) STOP 'sub_plot_x_y2: x and y dims don''t agree'

    ! Set defaults and replace them if optional arg present
    do j = 1,nJ
        write(jAsString,'(i0)') j
        dls(j) = 'data '//jAsString
    end do
    ms = 'line'
    cls = ''

    ! Write data to a temporary file
    open(unitno,file=trim(datafile),status='unknown')
    write(unitno,*) '# Temporary data file created by grey''s plot routine' 
    do i = 1,ni
        write(unitno,'(1PE24.15E3)',advance='no') x(i)
        do j = 1,nj
            write(unitno,'(1PE24.15E3)',advance='no') y(i,j)
        end do
        write(unitno,*) 
    end do
    close(unitno)

    ! Write config to a temporary file
    open(unitno,file=trim(configfile),status='unknown')
    write(unitno,*) '# Temporary gnuplot configuration file created by grey''s plot routine'
    write(unitno,*) 'set style data line' ! set line as default style 
    write(unitno,*) 'set term dumb'  !context' !enhanced background rgb "white"' ! Lines were added by Margaret Jacobson to output to a pdf
    write(unitno,*) 'set output "plot_save.txt"'                    ! added by Margaret Jacobson to output to a pdf
    if (isPresent_xlabel) write(unitno,*) 'set xlabel "'//trim(xlabel)//'"'
    if (isPresent_ylabel) write(unitno,*) 'set ylabel "'//trim(ylabel)//'"'
    if (isPresent_title) write(unitno,*) 'set title "'//trim(title)//'"'
    if (isPresent_extra) write(unitno,*) trim(extra)//' ' ! extra pre options

    write(unitno,'(A)',advance='no') 'plot' ! this is the first part of the command
    do j = 1,nj
        
        ! Create j+1 as a string
        write(jp1AsString,'(i0)') j+1

        ! Determine whether new line is needed
        if (j<nj) then
            tmpchar = newline
        else
            tmpchar = '   '
        end if

        ! Output the command
        ! NOTE: must use advance no here to prevent ifort from breaking the line
        write(unitno,'(A)',advance='no') ' "'//trim(datafile)//'" using 1:' &
                //trim(jp1AsString)//' title "'//trim(dls(j))//'" '//trim(tmpchar)
        ! Make a new line
        write(unitno,*)

    end do

!     if (isPresent_savefile) then 
!         write(unitno,*) 'set terminal postscript'
!         write(unitno,*) 'set output "'//trim(savefile)//'"'
!         write(unitno,*) 'replot'
!         write(unitno,*) 'set term x11' ! Not strictly necessary
!     end if

    close(unitno)

    ! Now use a nonstandard fortran routine (appearing in ifort and gfortran) to call gnuplot
    ! NOTE: if gnuplot is not available or an error is received, just report it and the program
    ! will continue.
    ierr = system('gnuplot -persist '//trim(configfile))
!    call execute_command_line('gnuplot -persist '//trim(configfile),exitstat=ierr)
    if (ierr/=0) print*,'sub_plot_x_y2: WARNING: error in calling gnuplot'

end subroutine plot_xy2_fine

! Acts just like the colon in Matlab
pure function colon(a,b)
    implicit none
    integer, intent(in) :: a,b
    integer, dimension(1:b-a+1) :: colon
    integer :: i
    do i = a,b
        colon(i-a+1) = i
    end do
end function colon

! linspace(a,b,n) (of course) constructs a grid of n linearly spaced points between
! a and b. If n==1, then the grid is just "a", consistent with Matlab. If n<=0, the routine
! stops with an error.
function linspace(x1,x2,n) result(grid)
    implicit none
    real(8), intent(in) :: x1,x2
    integer, intent(in) :: n
    real(8), dimension(1:n) :: grid
    !local
    integer :: i 

    if (n>1) then
        grid(1) = x1
        do i=2,n-1
            grid(i) = x1 + (x2-x1)*(real(i,8)-1.d0)/(real(n,8)-1.d0)
        end do   
        grid(n) = x2
    elseif (n==1) then
        grid(1) = x1
    else
        STOP 'linspace: ERROR: grid has size <= 0'
    end if
end function linspace

end module mod_plot


