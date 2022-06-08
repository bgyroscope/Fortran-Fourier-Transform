PROGRAM Osc_main 
!===============================================================
!
!    Program to numerically solve system of ODE 
!     Here is a simple model of OBE 
!
!    Ben Girodias
!==============================================================
use Osc_solver 

implicit none 

!variables
real*8 :: r_f
complex*16 :: z_y

integer, parameter :: N_fft = int( 1.0e4 )
complex*16, dimension(N_fft) :: z_y_fft
real*8, dimension(N_fft) :: r_w 

real*8 :: t_i, t_f
real*8 :: w_i, w_f
real*8 :: r_delw

complex*16,dimension(N_fft) :: z_exact 

integer :: Nsteps 
integer :: Nsave

integer :: j, k, l

!declare run parameters
Nsteps = int( 1.0e2 )
Nsave  = int( 1.0e3 )   !avoid greater than 10^5

!initial conditions
z_y = CMPLX( exp(-25.0d0/2.0d0) , 0.0d0 )  !initial condition
z_y_fft = (0.0d0, 0.0d0) 

t_i = 0.0d0
t_f = 1.0d1

w_i = -5.0d0
w_f = 5.0d0
r_delw = (w_f - w_i) / DBLE(N_fft) 

print *, "Nsteps, Nsave, z_y init, z_y_fft, t_i, t_f, w_i, w_f"
print *, Nsteps, Nsave, z_y, z_y_fft(1), t_i, t_f, w_i, w_f

do j=1,N_fft
  
   r_w(j) = w_i + DBLE(j)*r_delw    

end do  

!note this call below has Nsteps where the first element
!is NOT x0, y0
CALL rk4_solver ( z_y, z_y_fft, r_w, t_i, t_f, Nsteps, Nsave, N_fft )
z_exact = sqrt(8.0d0*atan(1.0d0))* &
            exp( CMPLX(- r_w*r_w / 2.0d0,0.0d0) &
                   + 5.0d0 * (0.0d0,1.0d0) * r_w ) 

 
open(20, file='fft_exact.dat')
do k=1,N_fft
   write(20,'(SE24.8E3 " , " SE24.8E3) ') &
                real(z_exact(k) ), aimag(z_exact(k) )
end do
close(20)  


END PROGRAM Osc_main 
