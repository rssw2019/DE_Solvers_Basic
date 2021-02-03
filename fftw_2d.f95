! This program generates two dimensional array(matrix) of random data,
! calculates the fourier transform of the matrix using FFTW routine
!and calculates the inverse fourier transform of the matrix using FFTW routine
!Written by Ramjee Sharma
!2/3/2021
!To compile the codes in cygwin: 
!  gfortran fftw_2d.f95 -I/usr/include -lfftw3 -lm
!To run run the executable:
! ./a.exePROGRAM fft

program fft2d

implicit none

integer (kind=4), parameter :: Nx = 8
integer (kind=4), parameter :: Ny = 8
integer (kind=4), parameter :: Nh = (Nx/2) + 1
real (kind=8), parameter :: pi = 3.1428

include "fftw3.f"

integer (kind=4) :: i,j
real (kind=8) :: dx, dy, Lx, Ly
real (kind=8) :: x(Nx), y(Ny), theta1(Nx,Ny), theta2(Nx,Ny)
complex (kind=8) :: ftheta(Nh,Ny)
integer (kind=8) plan_forward, plan_backward

Lx = 2.0d0*pi
Ly = 2.0d0*pi
dx = Lx/dfloat(Nx)
dy = Ly/dfloat(Ny)

do i = i,Nx
x(i) = (i-1)*dx
  do j = 1,Ny
  y(j) = (j-1)*dy
  theta1(i,j) = dsin(x(i))*dcos(y(j))
  write(*,*) i,j,theta1(i,j)
  write(100,*) i,j,theta1(i,j)
  enddo
enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, theta1, ftheta, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)

do i = 1,Nh
  do j = 1,Ny
  write(*,*) i,j,ftheta(i,j)
  write(200,*) i,j,ftheta(i,j)
  enddo
enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ftheta, theta2, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)

do i = 1, Nx
  do j = 1,Ny
  write(*,*) i,j, theta2(i,j)/(dfloat(Nx)*dfloat(Ny))
  write(300,*) i,j,theta2(i,j)
  enddo
enddo

end program fft2d
