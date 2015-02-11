module variables_and_subs
	! Parameters
	! Don't pass these parameters into subroutines: it will give conflict error
	
	double precision, parameter :: rho_0 = 1000
	double precision, parameter :: nu = 1e-6
	double precision, parameter :: Cv = 4.18	! specific heat capacity
	double precision, parameter :: R = 8.314
	
	integer, parameter :: N = 64
	double precision, parameter :: dx = 1.0
	double precision, parameter :: dy = 1.0
	double precision, parameter :: dt = 1e-4
	
	double precision, parameter :: T_0 = 298.15
	double precision, parameter :: U_0 = 1e-3


contains

subroutine setBC(rho, u, v, p, e)
	implicit none
	
	! dummy vars
	double precision, dimension(:,:) :: rho, u, v, p, e
	
	! @ left (x=0)
	u(:,1) = 0
	v(:,1) = 0
	e(:,1) = T_0/Cv
	!e(:,1) = e(:,2)
	!p(:,1) = p(:,2)
	rho(:,1) = rho(:,2)
	
	! right (x=L)
	u(:,N) = 0
	v(:,N) = 0
	e(:,N) = T_0/Cv
	!e(:,N) = e(:,N-1)
	!p(:,N) = p(:,N-1)
	rho(:,N) = rho(:,N-1)
	
	! bottom (y=0, row=N)
	u(N,:) = 0
	v(N,:) = 0
	e(N,:) = T_0/Cv
	!e(N,:) = e(N-1,:)
	!p(N,:) = p(N-1,:)
	rho(N,:) = rho(N-1,:)
	
	! top (y=L, row=1)
	u(1,:) = U_0
	v(1,:) = 0
	e(1,:) = T_0/Cv
	!e(1,:) = e(2,:)
	!p(1,:) = p(2,:)
	rho(1,:) = rho(2,:)
	
end subroutine setBC


subroutine setIC(rho, u, v, p, e)
	implicit none
	
	! dummy vars
	double precision, dimension(:,:) :: rho, u, v, p, e
	
	u(2:N-1,2:N-1) = 0.0
	v(2:N-1,2:N-1) = 0.0
	e(2:N-1,2:N-1) = T_0/Cv
	!p(2:N-1,2:N-1) = 0.0
	
	rho(1:N,1:N) = rho_0
	
	p = rho*R*e/Cv
	
end subroutine setIC


subroutine compute_d_rho_dt(rho, u, v, d_rho_dt, dir)
	implicit none
	
	double precision, dimension(:,:) :: rho, u, v, d_rho_dt
	
	double precision :: a, b, c, d
	integer :: row, col
	integer :: dir
	
	! forward difference
	if (dir==1) then
		do col=2,N-1
			do row=2,N-1
				a = ( u(row,col+1)-u(row,col) )/dx
				b = ( rho(row,col+1)-rho(row,col) )/dx
				c = ( v(row-1,col)-v(row,col) )/dy
				d = ( rho(row-1,col)-rho(row,col) )/dy
			
				d_rho_dt(row,col) = -( rho(row,col)*a + u(row,col)*b + rho(row,col)*c + v(row,col)*d)
			enddo
		enddo
	! backward difference
	else if (dir == -1) then
		do col=2,N-1
			do row=2,N-1
				a = ( u(row,col)-u(row,col-1) )/dx
				b = ( rho(row,col)-rho(row,col-1) )/dx
				c = ( v(row,col)-v(row+1,col) )/dy
				d = ( rho(row,col)-rho(row+1,col) )/dy
				
				d_rho_dt(row,col) = -( rho(row,col)*a + u(row,col)*b + rho(row,col)*c + v(row,col)*d)
			enddo
		enddo
	else 
		write(*,*) 'Error: Use 1 or -1 for forward or backward difference'
	endif
	
end subroutine 

subroutine compute_d_u_dt(rho, u, v, p, d_u_dt, dir)
	implicit none
	
	double precision, dimension(:,:) :: rho, u, v, p, d_u_dt
	
	double precision :: a, b, c, d, e
	integer :: row,col
	integer :: dir
	
	! forward difference
	if (dir==1) then
		do col=2,N-1
			do row=2,N-1
				a = ( u(row,col+1)-u(row,col) )/dx
				b = ( u(row-1,col)-u(row,col) )/dy
				c = ( p(row,col+1)-p(row,col) )/dx
				if (row==2 .or. row==N-1 .or. col==2 .or. col==N-1) then
					d = 0.0
					e = 0.0
				else
					d = (u(row,col+2)-2*u(row,col+1)+u(row,col) )/dx**2
					e = (u(row-2,col)-2*u(row-1,col)+u(row,col) )/dy**2
				endif
			
				d_u_dt(row,col) = -( u(row,col)*a + v(row,col)*b + ( 1/rho(row,col) )*c -nu*(d+e) )
			enddo
		enddo
	! backward difference
	else if (dir == -1) then
		do col=2,N-1
			do row=2,N-1
				a = ( u(row,col)-u(row,col-1) )/dx
				b = ( u(row,col)-u(row+1,col) )/dy
				c = ( p(row,col)-p(row,col-1) )/dx
				if (row==2 .or. row==N-1 .or. col==2 .or. col==N-1) then
					d = 0.0
					e = 0.0
				else
					d = (u(row,col)-2*u(row,col-1)+u(row,col-2) )/dx**2
					e = (u(row,col)-2*u(row+1,col)+u(row+2,col) )/dy**2
				endif
			
				d_u_dt(row,col) = -( u(row,col)*a + v(row,col)*b + ( 1/rho(row,col) )*c -nu*(d+e) )
			enddo
		enddo
	else 
		write(*,*) 'Error: Use 1 or -1 for forward or backward difference'
	endif
	
end subroutine 


subroutine compute_d_v_dt(rho, u, v, p, d_v_dt, dir)
	implicit none
	
	double precision, dimension(:,:) :: rho, u, v, p, d_v_dt
	
	double precision :: a, b, c, d, e
	integer :: row,col
	integer :: dir
	
	! forward difference
	if (dir==1) then
		do col=2,N-1
			do row=2,N-1
				a = ( v(row,col+1)-v(row,col) )/dx
				b = ( v(row-1,col)-v(row,col) )/dy
				c = ( p(row,col+1)-p(row,col) )/dx
				if (row==2 .or. row==N-1 .or. col==2 .or. col==N-1) then
					d = 0.0
					e = 0.0
				else
					d = (v(row,col+2)-2*v(row,col+1)+v(row,col) )/dx**2
					e = (v(row-2,col)-2*v(row-1,col)+v(row,col) )/dy**2
				endif
			
				d_v_dt(row,col) = -( u(row,col)*a + v(row,col)*b + ( 1/rho(row,col) )*c -nu*(d+e) )
			enddo
		enddo
	! backward difference
	else if (dir == -1) then
		do col=2,N-1
			do row=2,N-1
				a = ( v(row,col)-v(row,col-1) )/dx
				b = ( v(row,col)-v(row+1,col) )/dy
				c = ( p(row,col)-p(row,col-1) )/dx
				if (row==2 .or. row==N-1 .or. col==2 .or. col==N-1) then
					d = 0.0
					e = 0.0
				else
					d = (v(row,col)-2*v(row,col-1)+v(row,col-2) )/dx**2
					e = (v(row,col)-2*v(row+1,col)+v(row+2,col) )/dy**2
				endif
			
				d_v_dt(row,col) = -( u(row,col)*a + v(row,col)*b + ( 1/rho(row,col) )*c -nu*(d+e) )
			enddo
		enddo
	else 
		write(*,*) 'Error: Use 1 or -1 for forward or backward difference'
	endif
	
end subroutine 


subroutine compute_d_e_dt(rho, u, v, p, e, d_e_dt, dir)
	implicit none
	
	double precision, dimension(:,:) :: rho, u, v, p, e, d_e_dt
	
	double precision :: a, b, c, d
	integer :: row,col
	integer :: dir
	
	! forward difference
	if (dir==1) then
		do col=2,N-1
			do row=2,N-1
				a = ( e(row,col+1)-e(row,col) )/dx
				b = ( e(row-1,col)-e(row,col) )/dy
				c = ( u(row,col+1)-u(row,col) )/dx
				d = ( v(row-1,col)-v(row,col) )/dy
			
				d_e_dt(row,col) = -( u(row,col)*a + v(row,col)*b + ( p(row,col)/rho(row,col) )*c + ( p(row,col)/rho(row,col) )*d )
			enddo
		enddo
	! backward difference
	else if (dir == -1) then
		do col=2,N-1
			do row=2,N-1
				a = ( e(row,col)-e(row,col-1) )/dx
				b = ( e(row,col)-e(row+1,col) )/dy
				c = ( u(row,col)-u(row,col-1) )/dx
				d = ( v(row,col)-v(row+1,col) )/dy
			
				d_e_dt(row,col) = -( u(row,col)*a + v(row,col)*b + ( p(row,col)/rho(row,col) )*c + ( p(row,col)/rho(row,col) )*d )
			enddo
		enddo
	else 
		write(*,*) 'Error: Use 1 or -1 for forward or backward difference'
	endif
	
end subroutine 


subroutine predictor(rho, u, v, p, e, d_rho_dt, d_u_dt, d_v_dt, d_e_dt, rho_p, u_p, v_p, p_p, e_p)
	implicit none
	
	double precision, dimension(:,:) :: rho, u, v, p, e, d_rho_dt, d_u_dt, d_v_dt, d_e_dt, rho_p, u_p, v_p, p_p, e_p
	
	rho_p = rho + d_rho_dt*dt
	u_p = u + d_u_dt*dt
	v_p = v + d_v_dt*dt
	e_p = e + d_e_dt*dt
	p_p = rho_p*R*e_p/Cv
	
end subroutine predictor


subroutine corrector(rho, u, v, p, e, d_rho_dt_1, d_u_dt_1, d_v_dt_1, d_e_dt_1, d_rho_dt_2, d_u_dt_2, d_v_dt_2, d_e_dt_2, rho_p, u_p, v_p, e_p)
	implicit none
	
	double precision, dimension(:,:) :: rho, u, v, p, e, d_rho_dt_1, d_u_dt_1, d_v_dt_1, d_e_dt_1, d_rho_dt_2, d_u_dt_2, d_v_dt_2, d_e_dt_2, rho_p, u_p, v_p, e_p
	
	rho = rho + 0.5*( d_rho_dt_1 + d_rho_dt_2 ) * dt
	u = u + 0.5*( d_u_dt_1 + d_u_dt_2 ) * dt
	v = v + 0.5*( d_v_dt_1 + d_v_dt_2 ) * dt
	e = e + 0.5*( d_e_dt_1 + d_e_dt_2 ) * dt
	
	p = rho*R*e/Cv
	
end subroutine corrector

subroutine savedata(filename, a)
	implicit none
	
	character (len=*) :: filename
	double precision, dimension(:,:) :: a
	integer :: row, col
	
	open(unit=10, file=filename, status='replace', action='write')
	
	do row=1,N
		do col=1,N
			write(10,'(E)', advance='no') a(row,col)
		enddo
		write(10,*)
	enddo
	
	close(10)
	
	print *, filename, " written to disk."
	
end subroutine

end module variables_and_subs




program main
	use variables_and_subs
	!use mylib
	implicit none

	double precision, dimension(N,N) :: rho=0.0, u=0.0, v=0.0, p=0.0, e=0.0
	double precision, dimension(N,N) :: rho_p=0.0, u_p=0.0, v_p=0.0, p_p=0.0, e_p=0.0	! predictor field
	double precision, dimension(N,N) :: d_rho_dt_1=0.0, d_u_dt_1=0.0, d_v_dt_1=0.0, d_e_dt_1=0.0 	! gradients at t
	double precision, dimension(N,N) :: d_rho_dt_2=0.0, d_u_dt_2=0.0, d_v_dt_2=0.0, d_e_dt_2=0.0	! gradients at t+dt
	
	integer :: itr
	
	call setIC(rho, u, v, p, e)
	
	call setBC(rho, u, v, p, e)
	
	
	
	do itr = 1,100000
	
		write(*,*) itr
		
		! compute gradients at t
		call compute_d_rho_dt(rho, u, v, d_rho_dt_1, 1)
		call compute_d_u_dt(rho, u, v, p, d_u_dt_1, 1)
		call compute_d_v_dt(rho, u, v, p, d_u_dt_1, 1)
		call compute_d_e_dt(rho, u, v, p, e, d_e_dt_1, 1)
		
		! compute predictor field
		call predictor(rho, u, v, p, e, d_rho_dt_1, d_u_dt_1, d_v_dt_1, d_e_dt_1, rho_p, u_p, v_p, p_p, e_p)
		
		! compute gradients at t+dt
		call compute_d_rho_dt(rho_p, u_p, v_p, d_rho_dt_2, -1)
		call compute_d_u_dt(rho_p, u_p, v_p, p_p, d_u_dt_2, -1)
		call compute_d_v_dt(rho_p, u_p, v_p, p_p, d_u_dt_2, -1)
		call compute_d_e_dt(rho_p, u_p, v_p, p_p, e_p, d_e_dt_2, -1)
		
		! compute corrected field
		call corrector(rho, u, v, p, e, d_rho_dt_1, d_u_dt_1, d_v_dt_1, d_e_dt_1, d_rho_dt_2, d_u_dt_2, d_v_dt_2, d_e_dt_2, rho_p, u_p, v_p, e_p)
		
		call setBC(rho, u, v, p, e)
	
	enddo
	
	!call savedata('dudt.dat', d_u_dt_1)
	!call savedata('dvdt.dat', d_v_dt_1)
	!call savedata('drdt.dat', d_rho_dt_1)
	!call savedata('.dat', p_p)
	
	call savedata('u.dat', u)
	call savedata('v.dat', v)
	call savedata('rho.dat', rho)
	
end program main