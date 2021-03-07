!This program computes the solution of the 1-dimensional
!initial value problem dx/dt=f(x,t), x(t_0)=x_0
!using the 4th order Runge Kutta method
!and forward Euler method
!Ramjee Sharma,
!March 7, 2021
!University of North Georgia, USA
!gfortran ivp_1d.f95
!The results are stored in fort.100. Type the following command in gnuplot
! load 'ivp_1d.p'
!
program ivp_1d
  implicit none
  integer :: i
  integer, parameter :: n=60
  double precision, dimension(n) :: x,y,t
  double precision :: h
  double precision :: x1,t1,k1,x2,t2,k2,x3,t3,k3,k4

  t(1)=0.0d0
  x(1)=0.0d0
  y(1)=0.0d0
  h=0.1d0
  !Time integration
  ! x will give the result from RK4
  !and y will give the result from Forward Euler
  do i=1,n-1
    x(i+1)=rk4(x(i), t(i), h)
    y(i+1)=forward_euler(y(i),t(i),h)
    write(100,*) t(i), x(i),y(i)
    t(i+1)=t(i)+h
  enddo

  contains
    !4th order R-K method
    function rk4(x,t,h) result(y)
      implicit none
      double precision, intent(in) :: t, h
      double precision, intent(in) :: x
      double precision :: y
      integer :: i
      double precision :: x1,t1,k1,x2,t2,k2,x3,t3,k3,k4

      do i=1,n-1
        k1=dx_dt(x,t)
        x1=x+k1*h/2.0d0
        t1=t+h/2.0d0

        k2=dx_dt(x1,t1)
        x2=x+k2*h/2.0
        t2=t+h/2.0d0

        k3=dx_dt(x2,t2)
        x3=x+k3*h
        t3=t+h

        k4=dx_dt(x3,t3)

        y=x+(k1+2.0d0*k2+2.0d0*k3+k4)*h/6.0d0
      enddo
    end function rk4
    !
    !Forward Euler method
    function forward_euler(x,t,h) result(y)

      implicit none
      double precision, intent(in) :: t, h
      double precision, intent(in) :: x
      double precision :: y
      integer :: i
      double precision :: x1,t1,k1,x2,t2,k2,x3,t3,k3,k4

        y=x+h*dx_dt(x,t)

    end function forward_euler
    !
    ! Right side function f in dx/dt=f(x,t)
    function dx_dt(x,t) result(y)
      implicit none
      double precision, intent(in) :: x,t
      double precision :: y
      y=x+t
    end function dx_dt

end program ivp_1d
