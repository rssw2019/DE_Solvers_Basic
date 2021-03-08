!To compile and run in Cygwin/Widows machine
!gfortran -I/usr/local/include -L/usr/local/lib fluid_1d.f95 -lfftw3 -lm
!This program solves the 1-d fluid equation
! u_t+u*u_x=nu*u_xx, u(x,0)=u_0(x)
! x in [0, 2*pi], periodic boundary
! Using the Galerkin approximation/pseudospectral/Fourier method.
!This program will write two files fort.100 and fort.200.
!The file transport.p has codes to plot the data.
!To plot in gnuplot type: load 'fluid_1d_movie.p'
!Written by Ramjee Sharma
!University of North Georgia
!November 29, 2020

program fluid_1d

implicit none

include "fftw3.f"

!Declaring Variable Types
integer, parameter :: n = 128
integer, parameter :: half_n = n/2+1
double precision, parameter :: pi = 3.1415926535897932d0

double precision :: x, u, u_x, exact, error
dimension x(n), u(n), u_x(n)
dimension exact(n), error(n)
double complex :: uk,uk_x
dimension uk(half_n),uk_x(half_n)
integer*8 :: plan_forward, plan_backward
double precision :: x_0, x_n, period

integer :: i, tsteps, max_tsteps
double precision dx, k, dt, a, nu, max_error
double precision, dimension (n) :: right_side

!Initializing Parameters
x_0=0.0d0
x_n=2.0d0*pi
period=x_n - x_0
dx = period/dfloat(n)

dt = 0.00010d0
a = 0.030d0
nu=0.001d0
max_tsteps =10000

!space discretization and u(x,0). u(x,0) will be written in the file fort.100
do i = 1,n
  x(i) = dfloat(i-1)*dx
  u(i) = dsin(x(i))
  !exact solution for transport equation
  exact(i)=dsin(x(i)-a*dt*max_tsteps)
  write(100,*) x(i), u(i)
enddo

!Main loop begins

do tsteps = 1,max_tsteps

!Time integration using RK4 / Forward Euler
u=rk4(a,nu,dt,n,u)
!u=forward_euler(a,nu,dt,n,u)

if(mod(tsteps,1000)==0) then
  do i=1,n
    write(tsteps,*) x(i), u(i)
  enddo
endif

enddo
do i=1,n
  error(i)=abs(u(i)-exact(i))
enddo
!Printing max error for transport equation
!max_error=maxval(error)
!print*,"maximum error = ",max_error

!end of main loop

contains

  ! Right side function f in dx/dt=f(x,t)
  function dx_dt(a, nu, u) result(right_side)
    implicit none
    double precision, parameter :: pi = 3.1415926535897932d0
    integer, parameter :: n=128
    integer, parameter :: half_n=n/2+1

    double precision, intent(in):: u, a, nu
    double precision :: u_x, u_xx, k
    dimension u(n), u_x(n), u_xx(n)
    double complex :: uk,uk_x, uk_xx

    dimension uk(half_n),uk_x(half_n), uk_xx(half_n)
    integer*8 :: plan_forward, plan_backward

    integer :: i
    double precision :: right_side(n)

    !Taking FFT of u

    call dfftw_plan_dft_r2c_1d(plan_forward, n, u, uk, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan_forward, u, uk)
    call dfftw_destroy_plan(plan_forward)

    !Derivative in the Fourier Space/Spectral Space

    do i = 1,half_n
      k = 2.0d0*pi*dfloat(i-1)/(2.0d0*pi)
      uk_x(i) = (0.0d0,1.0d0) * k * uk(i)
      uk_xx(i)=-k*k*uk(i)
    enddo

    !Back to Physical space using inverse FFT

    call dfftw_plan_dft_c2r_1d(plan_backward, n, uk_x, u_x, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan_backward, uk_x, u_x)
    call dfftw_destroy_plan(plan_backward)

    call dfftw_plan_dft_c2r_1d(plan_backward, n, uk_xx, u_xx, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan_backward, uk_xx, u_xx)
    call dfftw_destroy_plan(plan_backward)

    !Normalizing ifft values

    do i = 1,n
      u_x(i) = u_x(i)/dfloat(n)
      u_xx(i) = u_xx(i)/dfloat(n)
      !right_side(i) = -a*u_x(i)
      right_side(i) = -u(i)*u_x(i)
      !right_side(i) = -u(i)*u_x(i)+nu*u_xx(i)
    enddo

  end function dx_dt
  !
  !RK 4 function
  function rk4(a,nu,dt,n,u) result (y)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: a, nu,dt, u(n)
    double precision :: k1(n), x1(n), k2(n)
    double precision :: x2(n), k3(n), x3(n), k4(n)
    double precision :: y(n)

    k1=dx_dt(a,nu,u)
    x1=u+k1*dt/2.0d0

    k2=dx_dt(a,nu,x1)
    x2=u+k2*dt/2.0

    k3=dx_dt(a,nu,x2)
    x3=u+k3*dt

    k4=dx_dt(a,nu,x3)
    y=u+(k1+2.0d0*k2+2.0d0*k3+k4)*dt/6.0d0

    end function rk4
    !
    !Forward Euler function
    function forward_euler(a,nu,dt,n,u) result (y)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: a, nu,dt, u(n)
      double precision :: y(n)

      y=u+dt*dx_dt(a,nu,u)

      end function forward_euler
end program fluid_1d
