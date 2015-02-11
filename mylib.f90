module mylib

	CONTAINS
subroutine number_of_lines(input_file, N)
    implicit none
	! dummy vars
	character (len=30), intent(in) :: input_file
	integer, intent(out) :: N
	
	! local vars
    integer :: temp=0, IOstatus=0, i=0

    open(unit=10, file=input_file, status="old", action="read")
    do
        read(10,*,iostat=IOstatus) temp
        if (iostatus < 0) then
            exit
        else
            i=i+1
        endif
    enddo
	    
    close(10)
    N = i
end subroutine

subroutine loadArrayFromFile(filename, number_of_lines, output_array)
    implicit none
    ! Dummy vars
    character (len=30), intent(in) :: filename
    integer, intent(in) :: number_of_lines
    real, dimension(number_of_lines), intent(inout) :: output_array
    
    ! local vars
    integer :: i
    
    open(unit=10, file=filename, status='old', action='read')
    do i=1,number_of_lines
        read(10,*) output_array(i)
    enddo
    
    close(10)
end subroutine

subroutine write1DArrayToFile(filename, array)
	implicit none
	
	!Dummy vars
	character (len=*) :: filename
	double precision, dimension(:) :: array
	! local vars
	integer :: N
	integer :: i
	
	open(unit=10, file=filename, status='replace', action='write')
	
	N = size(array)
	do i=1,N
		write(10,*) array(i)
	enddo
	close(10)
	print *, filename, " written to disk."
end subroutine write1DArrayToFile

subroutine write2DArrayToFile(filename, array)
	implicit none
	
	!Dummy vars
	character (len=*) :: filename
	double precision, dimension(:,:) :: array
	! local vars
	integer, dimension(1:2) :: N
	integer :: i, j
	
	open(unit=10, file=filename, status='replace', action='write')
	
	N = shape(array)
	do i=1,N(1)
		write(10,*) ( array(i,j), j = 1,N(2) )
	enddo
	close(10)
	print *, filename, " written to disk."
end subroutine write2DArrayToFile
	

subroutine linspace(xmin, xmax, n, x)
	implicit none
	! dummy vars
	real, intent(in) :: xmin
	real, intent(in) :: xmax
	integer, intent(in) :: n
	real, dimension(n), intent(out) :: x
	
	! local vars
	integer :: i
	real :: dx
	
	dx = (xmax - xmin) / (n-1)

	do i=1,n
		x(i) = xmin + dx * (i-1)
	enddo
end subroutine

end module mylib