module testmod_abc
integer :: k = 5

contains
subroutine print(a)
	real, dimension(:,:) :: a
	
	write(*,*) a(2,2)
	k = 10
	write(*,*) k
end subroutine 

end module

program test
	use testmod_abc
	implicit none
	
	integer, parameter :: N=5
	real, dimension(N,N) :: a=0.0
	
	a(2:N-1,2:N-1)=1.0
	
	call print(a)
end program
