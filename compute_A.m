%Computes A symbollically
clc
clear
close all

syms p x y gamma E u

%symbolic representation of the flux components
A1 = x;
A2 = (gamma - 1)*y - 0.5*(x^2)/p;
A3 = x*y/p + x*(gamma-1)*(y/p - 0.5*(x^2)/(p^2));

% A matrix comp.
A = [   [ diff( A1 , p), diff( A1 , x),diff( A1 , y) ]
        [ diff( A2 , p), diff( A2 , x),diff( A2 , y) ]
        [ diff( A3 , p), diff( A3 , x),diff( A3 , y)  ] ];

%Ouput A
A = subs(A, x, p*u);
A = subs(A, y, p*E);
disp(A)
