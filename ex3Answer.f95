program space_descritization

implicit none
integer, parameter :: n=4
integer i
real, parameter :: pi=3.14
real :: a,b,deltax, x(n)

a=0
b=2*pi

deltax=(b-a)/n

do i=1,n
  x(i)=(i-1)*deltax
  write(*,*) x(i)
end do

end program space_descritization
  
