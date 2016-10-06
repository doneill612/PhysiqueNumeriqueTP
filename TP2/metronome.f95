module programConsts
implicit none
real :: theta_t0, theta_tminusdt_t0
end module

module functionConsts
implicit none
real :: g, l, m1, m2,  Co, nu, PI
end module

program Metronome
implicit none
external :: solve
call solve(1000, .false.)
call solve(1000, .true.)
end program

subroutine solve(tau, dampened)
use programConsts
implicit none
external :: writeToFile
integer :: tau, t
logical :: dampened
real :: step, a, r, derivative, getNext, getNextWithDampening, temp
real, dimension(5) :: theta
t = 0
step = .01;
a = 0.0001
theta_t0 = 2.5
theta_tminusdt_t0 = 2.5
r = step / (2.0 * tau);
theta = (/ 0.0, theta_t0, theta_tminusdt_t0, 0.0, 0.0 /)
do while(t <= tau)
    theta(5) = t;
    if(dampened .eqv. .true.) then
	theta(1) = getNextWithDampening(theta, step, a, r)
    else
	theta(1) = getNext(theta, step, a, r)
    endif
    theta(4) = derivative(theta, step)
    call writeToFile(theta, dampened)
    temp = theta(2)
    theta(2) = theta(1)
    theta(3) = temp
    if(dampened .eqv. .true.) then
	a = a + 0.0000009
    endif
    t = t+1
enddo
end subroutine

subroutine writeToFile(theta, dampened)
implicit none
logical :: dampened
real, dimension(5) :: theta
if(dampened .eqv. .true.) then
    open(unit = 7, file="MetronomeDampened.out", status="unknown", access="append", action="write")
else
    open(unit = 7, file="MetronomeUndampened.out", status="unknown", access="append", action="write")
endif
write(7, *) theta(5), theta(2), theta(4)
close(7)
end subroutine

real function derivative(theta, step)
implicit none
real :: step
real, dimension(5) :: theta
derivative = (theta(1) - theta(3)) / (2.0 * step)
end function

real function getNext(theta, step, a, r)
use functionConsts
implicit none
real :: step, a, r, part1, part2
real, dimension(5) :: theta
g = 9.81
l = 0.2
Co = 1.0
nu = 1.3
PI = 3.14159265359879
m1 = 0.100
part1 = 2.0*theta(2)-theta(3)+(step**2.0)*((g*sin(theta(2))/l)-(Co/m1)*theta(2))
part2 = a*sin(2.0*PI*nu*theta(4))
getNext = part1 + part2
end function

real function getNextWithDampening(theta, step, a, r)
use functionConsts
implicit none
real :: step, a, r, part1, part2
real, dimension(5) :: theta
g = 9.81
l = 0.2
Co = 1.0
nu = 1.3
m2 = 0.025
PI = 3.14159265359879
part1 = 2.0 *theta(2)+(r-1)*theta(3)+(step**2.0)*((g*sin(theta(2))/l)-(Co/m1)*theta(2))
part2 = a*sin(2.0*PI*nu*theta(4))
getNextWithDampening = (part1 + part2)/(r+1)
end function
