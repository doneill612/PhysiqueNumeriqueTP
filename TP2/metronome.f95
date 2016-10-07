module functionConstants
implicit none
real :: g = 9.81, l = 0.20, m = 0.0250, Co = 1.0, nu = 1.3, pi = atan(-1.0)
end module

program Metronome
implicit none
external :: solve
call solve(100.0)
end program

subroutine solve(tau)
implicit none
external :: writeToFile
integer :: i
real :: tau, t, dt, r, a, next, derivative, potentialEnergy, kineticEnergy
real, dimension(5) :: theta
real, dimension(3) :: energy
t = 0.0
dt = 0.01
a = 0.1
r = dt / (2 * tau)
!-- theta = [ q(t+dt), q(t), q(t-dt), dq/dt, t ] --!
theta = (/ 0.0, 2.5, 2.5, 0.0, 0.0 /)
energy = (/ kineticEnergy(theta) , potentialEnergy(theta), potentialEnergy(theta) + kineticEnergy(theta) /) 
open(unit = 7, file = "MetronomeDampened.out", access = "append", action = "write")
open(unit = 8, file = "EnergyDampened.out", access = "append", action = "write")
do i =1,100000
    t = i * dt
    theta(5) = t
    theta(1) = next(theta, dt, a, r)
    theta(4) = derivative(theta, dt)
    energy(1) = kineticEnergy(theta)
    energy(2) = potentialEnergy(theta)
    energy(3) = energy(1) + energy(2)
    call writeToFile(theta, energy)
    theta(2:3) = theta(1:2)
    t = t + dt
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

real function potentialEnergy(theta)
use functionConstants
implicit none
real, dimension(5) :: theta
potentialEnergy = m * g * l * cos(theta(2)) + (1.0/2.0)*Co*(l*theta(2))**2.0
end function

real function kineticEnergy(theta)
use functionConstants
implicit none
real, dimension(5) :: theta
kineticEnergy = (1.0/2.0) * m * (l*theta(4))**2.0
end function

real function derivative(theta, dt)
implicit none
real :: dt
real, dimension(5) :: theta
derivative = (theta(1) - theta(3)) / (2.0 * dt)
end function

real function next(theta, dt, a, r)
use functionConstants
implicit none
real :: dt, a, r
real, dimension(5) :: theta
next = (2.0*theta(2)+(r-1.0)*theta(3)+((dt**2.0)*((g/l)*sin(theta(2))-((Co/m)*theta(2)) + a * sin(2.0*pi*nu*theta(5))))) / (r+1.0)
end function
