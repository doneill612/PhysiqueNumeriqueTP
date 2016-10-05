!-- Module pour stocker des valeurs constantes --!
module constants
implicit none
real :: g, R, S, P, n, l, theta0, omega_sqrd, gam, T, Tc, m
end module

!-- Le program - quand on fait la commande /nom_de_ficher.x, le code ci dessous s'execute --! 
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

!-- Function qui retourne une valeur real, qui correspond au root d'une fonction --!
!-- parametre - func : la fonction pour laquelle on veut trouver un root --!
!-- parametre - x0   : la position ou on commence a chercher --!
!-- parametre - eps  : la precision qu'il faut avoir pour satisfaire la recherche --!
real function Newton(func, x0, eps)
implicit none
!-- parametre - f           : la valeur de la fonction 'func' a la position actuelle --!
!-- parametre - df          : la valeur de la derivee de la fonction a la position actuelle --!
!-- parametre - xCurrent    : la position actuelle --!
!-- parametre - xNext       : la prochaine position pour chercher le root --!
!-- parametre - epsCurrent  : la precision actuelle --!
!-- parametre - root        : le root de la fonction (s'il existe) --!
real :: x0, eps, f, df, xCurrent, xNext, epsCurrent, root
!-- parametre - i           : l'iteration de la boucle actuelle --!
!-- parametre - iMax        : le nombre d'iterations qu'on va permettre --!
integer :: i, iMax
!-- parametre - converges   : est-ce que on a trouve un root avec la precision eps donnee? --!
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

!-- Subroutine qui calcule la valeur de sin(x) et de sin'(x) pour un x donne --!
!-- A noter : on definie f'(x) = cos(x), c'est a dire qu'on n'utilise pas une approximation --!
!-- parametre - x  : x donne --!
!-- parametre - f  : valeur de la fonction f(x) a la position x --!
!-- parametre - df : valeur de df(x) a la position x --!
subroutine SinFunc(x, f, df)
implicit none
real :: x, f, df
f = sin(x)
df = cos(x)
end subroutine

!-- Subroutine qui calcule la valeur de sin(x) et sin'(x) pour un x donne --!
!-- A noter : on definie f'(x) = [sin(x + dx/2) - sin(x - dx/2)] / dx, c'est a dire on utilise une approximation --!
!-- parametre - x  : x donne --!
!-- parametre - f  : valeur de la fonction f(x) a la position x --!
!-- parametre - df : valeur de df(x) a la position x --!
subroutine SinFuncApprox(x, f, df)
implicit none
real :: x, f, df, step
step = 1e-4
f = sin(x)
df = (sin(x + step / 2.0) - sin(x - step / 2.0)) / step
end subroutine

!-- Subroutine qui calcule la valeur de exp(x) et exp'(x) pour un x donne --!
!-- parametre - x  : x donne --!
!-- parametre - f  : valeur de la fonction f(x) a la position x --!
!-- parametre - df : valeur de df(x) a la position x --!
subroutine ExpFunc(x, f, df)
implicit none
real :: x, f, df
f = exp(x)
df = exp(x)
end subroutine

!-- Subroutine qui calcule la valeur de f(x) = (1/2) - tanh(x) et f'(x) pour un x donne --!
!-- parametre - x  : x donne --!
!-- parametre - f  : valeur de la fonction f(x) a la position x --!
!-- parametre - df : valeur de df(x) a la position x --!
subroutine TanhFunc(x, f, df)
implicit none
real :: x, f, df
f = (1.0/2.0) - tanh(x - 1)
df = -1.0 * (1.0 / cosh(1-x))**2.0
end subroutine

!-- Subroutine qui calcule la valeur de f(x) = (sin(x) / x) * (theta0^2.0 - x^2.0) - (gam * Tc / (2.0 * omega_sqrd)) et f'(x) pour un x donne --!
!-- Ce subroutine utilise le module "constants" qu'on a definie au debut du fichier --!
!-- parametre - x  : x donne --!
!-- parametre - f  : valeur de la fonction f(x) a la position x --!
!-- parametre - df : valeur de df(x) a la position x --!
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
