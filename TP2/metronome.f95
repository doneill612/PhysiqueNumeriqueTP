module functionConstants
implicit none
real :: g = 9.81, l = 0.20, m1 = 0.100, m2 = 0.025, nu = 1.3, Co = 1.0, pi = atan(-1.0)
end module

program Metronome
implicit none
external :: solve
call solve(10, .false.)
call solve(1000, .true.)
end program

subroutine solve(tMax, dampened)
implicit none
external :: writeToFile
logical :: dampened
integer :: tMax, t
real :: a, b, r, dt, nextDampened, next, derivative, potential, kinetic
real, dimension(5) :: theta
real, dimension(3) :: energy
dt = 0.01
theta = (/ 0.0, 2.5, 2.5, 0.0, 0.0 /)
energy = (/ kinetic(theta, dampened), potential(theta, dampened), sum((/kinetic(theta, dampened), potential(theta, dampened) /)) /)
if(dampened) then
	open(unit = 7, file = "MetronomeDampened.out", access = "append", action = "write")
	open(unit = 8, file = "EnergyDampened.out", access = "append", action = "write")
else
	open(unit = 7, file = "MetronomeUndampened.out", access = "append", action = "write")
	open(unit = 8, file = "EnergyUndampened.out", access = "append", action = "write")
endif
do t = 0, tMax * 100
	theta(5) = t * dt
	if(dampened .eqv. .true.) then
		theta(1) = nextDampened(theta, dt)
	else
		theta(1) = next(theta, dt)
	endif
	theta(4) = derivative(theta, dt)
	energy(1) = kinetic(theta, dampened)
	energy(2) = potential(theta, dampened)
	energy(3) = sum((/energy(1), energy(2) /))
	call writeToFile(theta, energy)
	theta(2:3) = theta(1:2)
enddo
close(7)
close(8)
end subroutine

subroutine writeToFile(theta, energy)
implicit none
real, dimension(5) :: theta
real, dimension(3) :: energy
write(7, *) theta(5), theta(2), theta(4)
write(8, *) energy(1), energy(2), energy(3), theta(5)
end subroutine

real function derivative(theta, dt)
implicit none
real :: dt
real, dimension(5) :: theta
derivative = (theta(1) - theta(3)) / (2.0 * dt)
end function

real function kinetic(theta, dampened)
use functionConstants
implicit none
logical :: dampened
real :: soln
real, dimension(5) :: theta
if(dampened .eqv. .true.) then
    soln = (1.0/2.0) * m2 * (l*theta(4))**2.0
else
    soln = (1.0/2.0) * m1 * (l*theta(4))**2.0
endif
kinetic = soln
end function

real function potential(theta, dampened)
use functionConstants
implicit none
logical :: dampened
real :: soln
real, dimension(5) :: theta
if(dampened .eqv. .true.) then
    soln = m2 * g * l * cos(theta(2)) + (1.0/2.0)*Co*(l*theta(2))**2.0
else
    soln = m1 * g * l * cos(theta(2)) + (1.0/2.0)*Co*(l*theta(2))**2.0
endif
potential = soln
end function

real function next(theta, dt)
use functionConstants
implicit none
real :: dt
real, dimension(5) :: theta
next = (2.0*theta(2)-theta(3)+((dt**2.0)*((g/l)*sin(theta(2))-((Co/m1)*theta(2))+ 0.0*sin(2.0*pi*nu*theta(5)))))
end function

real function nextDampened(theta, dt)
use functionConstants
implicit none
real :: dt, r, a, tau
real, dimension(5) :: theta
a = 0.5
tau = 100.0
r = dt / (2.0 * tau)
nextDampened=(2.0*theta(2)+(r-1.0)*theta(3)+((dt**2.0)*((g/l)*sin(theta(2))-((Co/m2)*theta(2))+a*sin(2.0*pi*nu*theta(5)))))/(r+1.0)
end function
