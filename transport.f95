!To compile and run in Cygwin/Widows machine
!gfortran -I/usr/local/include -L/usr/local/lib transport.f95 -lfftw3 -lm;./a.exe
!This program solves the transport equaiton 
! u_t+a*u_x=0, u(x,0)=u_0(x)
! x in [0, 2*pi], periodic boundary
! Using the Galerkin approximation/spectral/Fourier method.
!This program will write two files fort.100 and fort.200. The file transport.p has codes to plot the 
!data. 
!To plot in gnuplot type: load 'transport.p'
!Written by Ramjee Sharma
!November 29, 2020

program transport

implicit none

include "fftw3.f"

!Declaring Variable Types
integer, parameter :: num_of_pts = 1024
integer, parameter :: half_num_of_pts = num_of_pts/2+1
double precision, parameter :: pi = 3.1428

double precision :: x, ux, dux
dimension x(num_of_pts), ux(num_of_pts), dux(num_of_pts)
double complex :: uk,duk
dimension uk(half_num_of_pts),duk(half_num_of_pts)
integer*8 :: plan_forward, plan_backward
real*8 :: x_start, x_end,period

integer :: i, tsteps, max_tsteps
real*8 dx, k, dt, a
real*8, dimension (num_of_pts) :: right_side

!Initializing Parameters
x_start=0.0d0
x_end=2.0d0*pi
period=x_end - x_start
dx = period/dfloat(num_of_pts)

dt = 0.00010d0
a = 0.30d0
max_tsteps =20000

!space discretization and u(x,0). u(x,0) will be written in the file fort.100
do i = 1,num_of_pts
  x(i) = dfloat(i-1)*dx
  ux(i) = dsin(x(i))
  write(100,*) x(i), ux(i)
enddo

!Main loop begins
 
do tsteps = 1,max_tsteps

!Taking FFT of u(x,t)

call dfftw_plan_dft_r2c_1d(plan_forward, num_of_pts, ux, uk, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan_forward, ux, uk)
call dfftw_destroy_plan(plan_forward)

!Derivative in the Fourier Space/Spectral Space

do i = 1,half_num_of_pts
  k = 2.0d0*pi*dfloat(i-1)/period
  duk(i) = (0.0d0,1.0d0) * k * uk(i)
enddo

!Back to Physical space using inverse FFT

call dfftw_plan_dft_c2r_1d(plan_backward, num_of_pts, duk, dux, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan_backward, duk, dux)
call dfftw_destroy_plan(plan_backward)

!Normalizing ifft values

do i = 1,num_of_pts
  dux(i) = dux(i)/dfloat(num_of_pts)
enddo

!Time integration using Forward Euler

do i = 1,num_of_pts
  right_side(i) = -a*dux(i)
  ux(i) = ux(i) + dt * right_side(i)
enddo

enddo
!end of main loop

!writing final output in the file fort.200
do i = 1,num_of_pts
  write(200,*) x(i), ux(i)
enddo


end program transport

