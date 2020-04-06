%cNumber_full

%This file is to solve the transcendental equation
%   kappa*a^2=phi*sin^2(g*a*tau/2) 
%as a fuction of tau.

function iValue = getIValue(tau, gc, density)

myfun = @(I, t, g, d) d*sin(t/2*sqrt(I*g))^2-I;
%matlab sinc(x) = sin(pi x)/(pi x)
d = density;
g = gc;
t = tau;
fun = @(I) myfun(I,t,g,d);
iValue = fzero(fun,(2*pi/t)^2/gc);
% iValue = fzero(fun,density);
    