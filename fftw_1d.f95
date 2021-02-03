! This program generates a random set of data,
! calculates the fourier transform of the set using FFTW routine
!and calculates the inverse fourier transform using FFTW routine
!Written by Ramjee Sharma
!2/3/2021
!To compile the codes in cygwin: 
!  gfortran fftw_1d.f95 -I/usr/include -lfftw3 -lm
!To run run the executable:
! ./a.exe

PROGRAM fft

  IMPLICIT NONE

  INTEGER (kind=4), PARAMETER :: N = 10

  INTEGER (kind=4), PARAMETER :: Nh = N/2+1
  
  INCLUDE "fftw3.f"

  INTEGER (kind=4) i

  real (kind=8) u(N), u_dum(N)
  complex (kind=8) uk(Nh)

  integer (kind=8) plan_forward, plan_backward

  INTEGER :: seed
  seed = 123456789
  
  DO i = 1,N
  u(i) = ran(seed)
  seed = u(i)*seed
!  WRITE(*,*) i, u(i), "u"
  ENDDO
  
  call dfftw_plan_dft_r2c_1d_ (plan_forward, N, u, uk, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_forward)
  call dfftw_destroy_plan_ (plan_forward)   

  do i = 1,Nh
!  write(*,*) i, uk(i)
  enddo

  call dfftw_plan_dft_c2r_1d_ (plan_backward, N, uk, u_dum, FFTW_ESTIMATE)
  call dfftw_execute_ (plan_backward)
  call dfftw_destroy_plan_ (plan_backward)

  do i = 1,N
  write(*,*) i, u(i), u_dum(i)/dfloat(N)
  enddo
END PROGRAM fft
