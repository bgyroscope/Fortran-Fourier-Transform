MODULE  Osc_solver
!===================================================================
!
!    Module built to solve simple oscilating ode 
!    Uses constant time steps steps
!
!    Ben Girodias
!==================================================================

IMPLICIT NONE



!subroutines accessed from outside
public :: rk4_solver

CONTAINS

SUBROUTINE RK4(r_t, r_w, z_y, z_yfft, z_ydt, z_yfftdt, N_fft,h)
!=================================================================
!
!     Subroutine to find f at the next time step 
!     Must give the initial value and derivative
!     Note that initial input is overwritten!
!
!     Ben Girodias
!==============================================================
    
    IMPLICIT NONE
  
 
    !arguments 
    integer N_fft

    real*8 :: r_t 
    real*8, dimension(N_fft) :: r_w

    complex*16 :: z_y, z_ydt     !original is  overwritten
    complex*16, dimension(N_fft)  :: z_yfft, z_yfftdt
    real*8 :: h              !time step to take    

 
    !Local variable
    real*8 :: hh, th, h6        !half of h, x_n+h/2.0, h/6.0

    complex*16 :: z_y_hlp                 !help f values
    complex*16 :: z_ydt_hlp1, z_ydt_hlp2  !helping derivatives 
 
    complex*16,dimension(N_fft) :: z_yfft_hlp          !help f values
    complex*16,dimension(N_fft):: z_yfftdt_hlp1, z_yfftdt_hlp2 
                           !helping derivatives 
   
    hh = h*0.5d0; th = r_t + hh ; h6 = h/ 6.0d0
    
    !first part of RK4  ie yn + k1 / 2
    z_y_hlp = z_y + hh * z_ydt
    z_yfft_hlp = z_yfft + hh * z_yfftdt

    !second part of RK4 
    !next call to get k2, then get yn + k2 /2
    CALL RHS( th, r_w, z_y_hlp, z_yfft_hlp, &
                 z_ydt_hlp1,z_yfftdt_hlp1, N_fft  ) 
  
    z_y_hlp = z_y + hh * z_ydt_hlp1
    z_yfft_hlp = z_yfft + hh * z_yfftdt_hlp1

    !third  part of RK4 to get k3, then yn + k3
    CALL RHS( th, r_w, z_y_hlp, z_yfft_hlp, &
                 z_ydt_hlp2, z_yfftdt_hlp2, N_fft ) 
   
    z_y_hlp = z_y + h  * z_ydt_hlp2
    z_yfft_hlp = z_yfft + h  * z_yfftdt_hlp2
 
    z_ydt_hlp2 = z_ydt_hlp1 + z_ydt_hlp2 !storing k2 + k3 in dym
    z_yfftdt_hlp2 = z_yfftdt_hlp1 + z_yfftdt_hlp2

    !final part of RK4
    CALL RHS( r_t+h, r_w, z_y_hlp, z_yfft_hlp, &
                 z_ydt_hlp1, z_yfftdt_hlp1, N_fft) 
      !store k4 as dyt
    
    !output of the subroutine
    z_y = z_y + h6 * ( z_ydt + z_ydt_hlp1 + 2.0d0*z_ydt_hlp2)
    z_yfft = z_yfft +  &
          h6 * ( z_yfftdt + z_yfftdt_hlp1 + 2.0d0*z_yfftdt_hlp2) 
    
END SUBROUTINE RK4


SUBROUTINE rk4_solver(z_y, z_yfft,r_w,t_i, t_f, Nsteps, Nsave, N_fft )
!=====================================================================
!
!     Subroutine that solves a system of ODE's using RK4
!     with constant step size
!
!     Uses: Subroutine RK4 - for next step
!           Subroutine RHS - to calculate  RHS of system of equations
!
!     Ben Girodias 
!====================================================================

IMPLICIT NONE

!inputs of sub
complex*16 :: z_y      

integer :: N_fft
complex*16, dimension(N_fft) :: z_yfft

real*8, dimension(N_fft) :: r_w

real*8 :: t_i, t_f    !initial and final times
integer :: Nsteps     !number of steps
integer :: Nsave      !number of points saved

!internal vars
complex*16 :: z_ydt
complex*16, dimension(N_fft) :: z_yfftdt

real*8 :: h !time step
integer :: Fast_Interval  !how often to save
integer :: Slow_Interval  !how oftern to save small grid

real*8 :: r_t   !tracks time value
integer :: j, k

character(len=100) :: format_V

character(len=100) :: time_data
character(len=100) :: fft_data
integer :: save_count

h  = (t_f - t_i ) / DBLE( Nsteps )

If ( Nsave > Nsteps) Then
   Fast_Interval = 1
   Slow_Interval = int( Nsteps / 10.0)
   print *, 'strange Nsave. Saving every point...'
Else
   Fast_Interval = int( Nsteps / Nsave) 
   Slow_Interval = int( Nsteps / 10.0 )    ! right now just 10 slices.  
end if

!saving the names of the file.
write( time_data, '(A13)') 'simpleRK4.dat'

!format V
write(format_V, '(A12, I4, A18)') &
   '(SE24.8e3, ', 2, '(" , ",SE24.8e3))'
!2 for now. Will be Param in module

open( 50, file='r_w.dat') 
write( 50, '(SE24.8e3)') r_w
close(50)


!do the first point manually. 
CALL RHS( t_i,r_w, z_y,z_yfft, z_ydt, z_yfftdt, N_fft)
CALL RK4 ( t_i,r_w, z_y,z_yfft, z_ydt, z_yfftdt, N_fft, h)
r_t = t_i + h

open( 10, file=time_data)
write(10, *)
close(10)

save_count = 1  !to count number of time saved. 
do j=2,Nsteps
 
   !save data here
   If( MOD( j, Fast_Interval)  == 0) Then
  
      !save time data
      open( 10, file=time_data, status='old',&
               position='append', action='write') 
      write(10, format_V) r_t, real(z_y) , aimag(z_y)
      close(10)
   
      If( MOD( j, Slow_Interval) == 0) Then 
          !save fft data. new file each time it saves. 
          write( fft_data, '(A4 I5.5 A4)' ) 'fft_', save_count, '.dat'
          
            
          open(20, file=fft_data)
          do k=1,N_fft
             write(20,'(SE24.8E3 " , " SE24.8E3) ') &
                          real(z_yfft(k) ), aimag(z_yfft(k) )
          end do
          close(20)  

          save_count = save_count+1
          
      end if
   
   end if 

   CALL RHS( r_t, r_w, z_y, z_yfft, z_ydt, z_yfftdt, N_fft)
   CALL RK4 ( r_t,r_w, z_y, z_yfft, z_ydt, z_yfftdt, N_fft, h)
   r_t = r_t + h


end do

!save final data
 
open(20, file='fft_final.dat')
do k=1,N_fft
   write(20,'(SE24.8E3 " , " SE24.8E3) ') &
                real(z_yfft(k) ), aimag(z_yfft(k) )
end do
close(20)  


END SUBROUTINE rk4_solver 


SUBROUTINE RHS( r_t,r_w, z_y, z_yfft, z_ydt, z_yfftdt, N_fft)
!==================================================================
!
!    Actual Equations to solve
!
!    Ben Girodias
!==================================================================
implicit none

real*8 :: r_t
complex*16 :: z_y, z_ydt

integer :: N_fft
real*8, dimension(N_fft) :: r_w
complex*16, dimension(N_fft) :: z_yfft, z_yfftdt

!local
real*8 :: r_shift

complex*16, parameter :: i = (0.0d0, 1.0d0)

r_shift = 5.0d0
z_ydt = -(r_t-r_shift) * ExP( -(r_t-r_shift)*(r_t-r_shift)/2.0d0 )  !i * z_y
z_yfftdt = z_y * Exp( i * r_w * r_t)

END SUBROUTINE RHS 



END MODULE
