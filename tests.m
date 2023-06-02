% Some extra tests:

clc
clear
close all

%Test y_interp
x = linspace(0,10,20);
dx = x(2)-x(1);
x_L = x - dx/2;
x_R = x + dx/2;
y = sin(x);
[y_plus,y_minus] = edges(y);

%Plot the result
plot(x,y,"-o")
hold on
plot(x_R,y_plus,"-*")
hold on
plot(x_L,y_minus,"-*")
legend("Original", "interp R", "interp L")
title("Test interpolate local")
xlabel("x")
ylabel("y")


%Sup testing
disp(sup([-5,2,1,5]));


% Compute the edge values
function [W_plus, W_minus] = edges(W_tilde)
%Linear assumption

%Compute W_plus and W_minus
% (Linear)
% W_plus = W_tilde - dW/2;
% W_minus = W_tilde + dW/2;

%Interpolate 1:
%W_tilde = W_tilde';
sz = size(W_tilde);
sz2 = [sz(1),sz(2)+1];
W_tilde_interp = zeros(sz2);
W_plus = zeros(sz);
W_minus = zeros(sz);
for i = 1:sz2(1)
    W_tilde_interp(i,:) = interp_center_to_edge_local(W_tilde(i,:));
    W_plus(i,:) = W_tilde_interp(i,2:end);
    W_minus(i,:) = W_tilde_interp(i,1:end-1);
end
%W_plus = W_plus';
%W_minus = W_minus';


end


% Local center to edge, with unique fluxes
function [y_interp] = interp_center_to_edge_local(y)
Nx = max(size(y));
x = linspace(0,1,Nx);
dx = x(2)-x(1);
x2 = linspace(0+dx/2,1-dx/2,Nx-1);
y_interp = interp1(x,y,x2,'spline');
y_interp = [y_interp(end),y_interp,y_interp(1)]; 
end

%Sup function
function [sup_val] = sup(values)
[maxB, index] = max(abs(values));
sup_val = maxB * sign(values(index));
end