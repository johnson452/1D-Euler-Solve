function [rho,u,p,E,grid] = push_all(rho,u,p,E,grid)

%Push all quantities with the MUSCL scheme

%Update energy

%Grab U, F vectors
U = [ rho ; rho.*u ; rho.*E ];
%F = [ rho.*u  ; p + rho.*u.*u ; u.*(E + p) ];

%Build the AQ matrix:
A = AQ(u,E,grid);

%Iterate over the domain
Nx = grid.Nx;

%Assumes periodic (Grabs left most/ right most element)
%  Nx = 5: R = [2     3     4     5     1];
%  Null : I = [1     2     3     4     5];
%Nx = 5: L = [5     1     2     3     4];
%Nx = 5: L = [G     X     X     X     G];
% --> Ghost cells (G) see periodic (if they are peroidic,
% otherwise they are overwritten by BC)
R = mod( linspace(1,Nx,Nx), Nx) + 1; %Good
L = mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

%Reconstruct slope dU within a cell
deltaU = differential(U,U(:,L),U(:,R),grid); %Good


%Make U_tilde
U_tilde = zeros(3,Nx);

%ARTIFICIAL
deltaU(:,1) = zeros(3,1);
deltaU(:,Nx) = zeros(3,1);

%Update Q
for i = 1:Nx
    U_tilde(:,i) = U(:,i) - grid.dt/(2*grid.dx)*A(:,:,i)*deltaU(:,i);
end

%Get edge values (edges_linear)
% Output: (cell edge values) in cell i
% [ i - 1/2 (+), i + 1/2 (-) ]
[U_plus_I, U_minus_I] = edges_linear(U_tilde,deltaU);
[U_plus_R, ~] = edges_linear(U_tilde(:,R),deltaU(:,R));
[~, U_minus_L] = edges_linear(U_tilde(:,L),deltaU(:,L));

%Get edge values
% [U_plus_I, U_minus_I] = edges(U_tilde);
% [U_plus_R, ~] = edges(U_tilde(:,R));
% [~, U_minus_L] = edges(U_tilde(:,L));

%Calulcate fluxes takes: ( [ a | b ] at i + 1/2 )
F_R =  Flux(U_minus_I,U_plus_R,grid);
F_L = Flux(U_minus_L,U_plus_I,grid);

%ARTIFICIAL - temp boundary fix
F_R(:,1) = zeros(3,1);
F_L(:,1) = zeros(3,1);
F_R(:,Nx) = zeros(3,1);
F_L(:,Nx) = zeros(3,1);

%Compute the updated U
U = U - grid.dt/(grid.dx)*(F_R - F_L);

%grab variables out of U
rho = U(1,:);
u = U(2,:)./rho;
E = U(3,:)./rho;

%Update p:
p = EOS(rho,u,E,grid,"calorically_ideal");

end


%Reconstruction
function [dW] = differential(Wi,Wm,Wp,grid)

%Option:
option = "eigenvector";

%Average Dw
if option == "eigenvector"
    % Reconstruct primative variables
    Q = Wi;
    rho = Q(1,:);
    u_vec = Q(2,:)./rho;
    E_vec = Q(3,:)./rho;

    %Iterate through the grid
    Nx = grid.Nx;
    dW = zeros(3,Nx);
    for i = 1:Nx
        u = u_vec(i);
        E = E_vec(i);

        % Compute eigenvectors time differences
        L = left_eigenvector(u,E,grid);
        R = right_eigenvector(u,E,grid);

        %Compute the differences with the primative vars
        delta_i_minus_1 = L*(Wi(:,i) - Wm(:,i));
        delta_i = L*(Wp(:,i) - Wi(:,i));
        dW(:,i) = R*ave( delta_i, delta_i_minus_1);
    end
else %Standard Option
    dW = ave( Wi - Wm, Wp - Wi );
end

end


% Calculate the A matrix
function [A] = AQ(u_vec,E_vec,grid)

%Calculate elements of A
% A11 = u;
% A21 = (grid.gamma-1)*E + (1.5-grid.gamma)*u.*u;
% A31 = grid.gamma*u.*E - 0.5*(grid.gamma-1)*u.*u.*u;
%
% A12 = zeros(1,grid.Nx) + 1;
% A22 = (3/2 - grid.gamma)*u;
% A32 = grid.gamma*E - 0.5*(grid.gamma - 1)*u.*u;
%
% A13 = zeros(1,grid.Nx) + 0;
% A23 = zeros(1,grid.Nx) + (grid.gamma - 1);
% A33 = u*grid.gamma;


% A11 = u;
% A12 = zeros(1,grid.Nx) + 1;
% A13 = zeros(1,grid.Nx) + 0;
% A21 = (grid.gamma - 1)*grid.specific_internal_energy + u.*u;
% A22 = u;
% A23 = zeros(1,grid.Nx) + 0;
% A31 = grid.specific_internal_energy*u*grid.gamma;
% A32 = 0.5*u.*u;
% A33 = u;

%Build A:
A = zeros(3,3,grid.Nx);

%Calculate A:
gamma = grid.gamma;

for i = 1:grid.Nx
    u = u_vec(i);
    E = E_vec(i);
    A(:,:,i) = [ [ 0                    ,       1                                ,   0          ];...
        [ ((u^2)*gamma - 3*u^2)/2,   3*u- u*gamma                         ,  gamma - 1   ];...
        [ (u^3 - E*u)*gamma - u^3, - (1/2)*( (3*u^2 - 2*E)*gamma - 3*u^2 ),  u * gamma   ] ];
    %[ [ A11(i), A12(i), A13(i) ];...
    %            [ A21(i), A22(i), A23(i) ];...
    %           [ A31(i), A32(i), A33(i) ] ];
end


end


% Averaging
function [dW] = ave( Wm, Wp )

avg_type = "minmod"; % "Supebee"; %  "standard";
a = Wm; b = Wp;
ab = a.*b;
sz_a = size(a);
dW = zeros(sz_a);

% Standard Averaging
if avg_type == "standard"
    dW = (Wm + Wp)/2;
elseif avg_type == "minmod"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            if ab(i,j) > 0
                dW(i,j) = minmod( [(a(i,j) + b(i,j))/2 , 2*a(i,j), 2*b(i,j)] );
            else
                dW(i,j) = 0;
            end
        end
    end
elseif avg_type == "Supebee"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            if ab(i,j) > 0
                max_v =  maxmod(  [ a(i,j)  , b(i,j)  ] );
                min_v = minmod (  [ 2*a(i,j), 2*b(i,j)] );
                dW(i,j) = minmod( [ max_v   , min_v   ] );
            else
                dW(i,j) = 0;
            end
        end
    end
end
end


function [W_plus, W_minus] = edges_linear(W_tilde,dW)
W_plus = W_tilde - dW/2;
W_minus = W_tilde + dW/2;
end


% Compute the edge values
function [W_plus, W_minus] = edges(W_tilde)
% |+  x  -|
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
    W_plus(i,:) = W_tilde_interp(i,1:end-1);
    W_minus(i,:) = W_tilde_interp(i,2:end);
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
y_interp = [-1e30,y_interp,-1e30];
end


%Fluxes N
function [Fl] = Flux(Q_left, Q_right, grid)

%Disect the peices
rhoR = Q_right(1,:);
uR = Q_right(2,:)./rhoR;
ER = Q_right(3,:)./rhoR;
rhoL = Q_left(1,:);
uL = Q_left(2,:)./rhoL;
EL = Q_left(3,:)./rhoL;
pR = EOS(rhoR,uR,ER,grid,"calorically_ideal");
pL = EOS(rhoL,uL,EL,grid,"calorically_ideal");
%pR = rhoR.*(grid.gamma-1)*grid.specific_internal_energy;
%pL = rhoL.*(grid.gamma-1)*grid.specific_internal_energy;


%Compute fluxes
FR = [ rhoR.*uR  ; pR + rhoR.*uR.*uR ; uR.*(rhoR.*ER + pR) ];
FL = [ rhoL.*uL  ; pL + rhoL.*uL.*uL ; uL.*(rhoL.*EL + pL) ];

% compute c (sup of the eigenvalues of left and right)
gamma = grid.gamma;
lambda1_L = uL;
lambda2_L = - (1/2)*( sqrt(2)*sqrt((2*EL - uL.^2)*gamma^2 + (uL.^2 - 2*EL)*gamma ) - 2* uL  );
lambda3_L =  (1/2)*( sqrt(2)*sqrt((2*EL - uL.^2)*gamma^2 + (uL.^2 - 2*EL)*gamma ) + 2* uL  );
lambda1_R = uR;
lambda2_R = - (1/2)*( sqrt(2)*sqrt((2*ER - uR.^2)*gamma^2 + (uR.^2 - 2*ER)*gamma ) - 2* uR  );
lambda3_R = (1/2)*( sqrt(2)*sqrt((2*ER - uR.^2)*gamma^2 + (uR.^2 - 2*ER)*gamma ) + 2* uR  );
%lambda = [lambda1_L, lambda2_L, lambda3_L, lambda1_R, lambda2_R, lambda3_R];
% c = max( lambda );

Nx = grid.Nx;
c = zeros(1,Nx);

%Old, but working
% for i = 1:Nx
%     if uR(i) > 0 && uL(i) > 0
%         c(i) = max( uR(i), uL(i) );
%     elseif uR(i) < 0 && uL(i) > 0
%         if abs(uR(i)) > uL(i)
%             c(i) = uR(i);
%         else
%             c(i) = uL(i);
%         end
%     elseif uR(i) > 0 && uL(i) < 0
%         if abs(uL(i)) > uR(i)
%             c(i) = uL(i);
%         else
%             c(i) = uR(i);
%         end
%     else
%         c(i) = min( uR(i), uL(i) );
%     end
% end
%End old but working

%Improved approx:
for i = 1:Nx
    eigen_values = [lambda1_L(i),lambda2_L(i),lambda3_L(i),...
        lambda1_R(i),lambda2_R(i),lambda3_R(i)];
    c(i) = sup(eigen_values);
end

%Rusanov Flux
Fl = (1/2) * ( FR + FL  - c.*( Q_right - Q_left ) );

end

% Minmod used for averaging
function [val] = minmod(a)
if (max(a) > 0) && (min(a) >  0)
    val = min(a);
elseif (max(a) < 0) && (min(a) <  0)
    val = max(a);
else
    val = 0;
end
end

%Maxmod used for averaging
function [val] = maxmod(a)
if (max(a) > 0) && (min(a) >  0)
    val = max(a);
elseif (max(a) < 0) && (min(a) <  0)
    val = min(a);
else
    val = 0;
end
end


% Left eigenvectors (Matrix form)
function [L] = left_eigenvector(u, E, grid)
gamma = grid.gamma;
a = (sqrt(2)*sqrt((2*E-u^2)*gamma^2+(u^2-2*E)*gamma))/2;
d = (gamma-1)/a^2;
L11 = (a*d*u^2 + 2*u ) / (4*a);
L12 = - (a*d*u + 1 ) / (2*a);
L13 = (d/2);
L21 = (a*d*u^2 - 2*u ) / (4*a);
L22 = - (a*d*u - 1 ) / (2*a);
L23 = (d/2);
L31 = - (d*u^2 - 2 ) / (2);
L32 = d*u;
L33 = -d;
L = [   [ L11, L12, L13];...
    [ L21, L22, L23];...
    [ L31, L32, L33] ];
end

% Right eigenvectors (Matrix form)
function [R] = right_eigenvector(u, E, grid)
gamma = grid.gamma;
a = (sqrt(2)*sqrt((2*E-u^2)*gamma^2+(u^2-2*E)*gamma))/2;
R11 = 1;
R12 = 1;
R13 = 1;
R21 = u-a;
R22 = u+a;
R23 = u;
R31 = (u^2*(gamma-1)-2*a*u*(gamma-1)+2*a^2)/(2*(gamma-1));
R32 = (u^2*(gamma-1)+2*a*u*(gamma-1)+2*a^2)/(2*(gamma-1));
R33 = u^2/2;
R = [   [ R11, R12, R13];...
    [ R21, R22, R23];...
    [ R31, R32, R33] ];
end


%Sup function
function [sup_val] = sup(values)
[maxB, index] = max(abs(values));
sup_val = maxB * sign(values(index));
end