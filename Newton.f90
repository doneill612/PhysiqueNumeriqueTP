module constants
implicit none
real :: g, R, S, P, n, l, theta0, omega_sqrd, gam, T, Tc, m
end module

program ZeroFinder
implicit none
external :: SinFunc, SinFuncApprox, ExpFunc, TanhFunc, Landau
real :: Newton
write(*,*) "Zero Finder using Newton's Method"
write(*,*)
write(*,*) "6.a) ", Newton(SinFunc, 0.1, 0.001)
write(*,*) "6.b) ", Newton(SinFuncApprox, 0.1, 0.001)
write(*,*) "6.c) ", Newton(SinFunc, 1.55, 0.001)
write(*,*) "6.d) ", Newton(ExpFunc, 0.0001, 0.00001)
write(*,*) "6.e) ", Newton(TanhFunc, 1.0, 0.001)
write(*,*) "7) ", Newton(Landau, 0.5, 0.001)
end program

real function Newton(func, x0, eps)
implicit none
real :: x0, eps, f, df, xCurrent, xNext, epsCurrent, root
integer :: i, iMax
logical :: converges
i = 0
iMax = 10000
f = 0
df = 0
root = 0
epsCurrent = 100
xNext = 0
xCurrent = x0
do while(abs(epsCurrent) > eps .and. i .le. iMax)
	call func(xCurrent, f, df)
	if(f .ne. df) then
		xNext = xCurrent - f / df
		epsCurrent = xNext - xCurrent
		root = xNext
		xCurrent = xNext
	else if(f == df .and. f .le. eps) then
		cycle
	endif
	i = i+1
enddo
if(i .ge. iMax) then
	converges = .false.
else
	converges = .true.
endif
if(converges .eqv. .false.) then
	write(*,*) "No root found - does not converge under constraint eps" 
endif
Newton = root
end function

subroutine SinFunc(x, f, df)
implicit none
real :: x, f, df
f = sin(x)
df = cos(x)
end subroutine

subroutine SinFuncApprox(x, f, df)
implicit none
real :: x, f, df, step
step = 1e-4
f = sin(x)
df = (sin(x + step / 2.0) - sin(x - step / 2.0)) / step
end subroutine

subroutine ExpFunc(x, f, df)
implicit none
real :: x, f, df
f = exp(x)
df = exp(x)
end subroutine

subroutine TanhFunc(x, f, df)
implicit none
real :: x, f, df
f = (1.0/2.0) - tanh(x - 1)
df = -1.0 * (1.0 / cosh(1-x))**2.0
end subroutine

subroutine Landau(x, f, df)
use constants
implicit none
real :: x, f, df
g = 9.81
R = 8.32
S = 1e-4
P = 200.0
T = 300.0
l = 0.1
n = 8e-7
theta0 = 1
Tc = (m * g * l * theta0**2.0) / (2.0 * n * R)
m = 6e-3
omega_sqrd = g / l
gam = 2.0 * n * R / (m * l**2.0)
f = (sin(x) / x) * (theta0**2.0 - x**2.0) - (gam * Tc / (2.0 * omega_sqrd))
df = (theta0**2.0 - x**2.0) * ((x*cos(x) - sin(x)) / x**2.0) - 2.0 * sin(x)
end subroutine
