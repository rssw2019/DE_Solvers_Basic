program plotfunction
implicit none
integer :: i
real :: x
real, parameter :: xmin = 0.,xmax=10., a=-2.
open(10,file='myplot.dat')
do i = 1,100
x = xmin + xmax*(i-1.0)/(100.0-1.0)
write(10,*) x,f(x)
enddo
close(10)

contains
function f(x)
implicit none
real :: f,x
f = cos(x+a)
end function f
end program plotfunction