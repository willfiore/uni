clear variables;
clc; close all;

% Number of possible states
S = 4;

% Single site operators
sq = arrayfun(@sqrt, 1:S); 
a = diag(sq, 1);            % Annihilation
a_d = diag(sq, -1);         % Creation
n = diag(0:size(a)-1);      % Number (counts bosons in a site)

I = eye(size(a));           % Identity

% Supercell operators (generalise this later)
a1 = kron(a, kron(I, kron(I, I)));
a2 = kron(I, kron(a, kron(I, I)));
a3 = kron(I, kron(I, kron(a, I)));
a4 = kron(I, kron(I, kron(I, a)));

a_d1 = kron(a_d, kron(I, kron(I, I)));
a_d2 = kron(I, kron(a_d, kron(I, I)));
a_d3 = kron(I, kron(I, kron(a_d, I)));
a_d4 = kron(I, kron(I, kron(I, a_d)));

n1 = kron(n, kron(I, kron(I, I)));
n2 = kron(I, kron(n, kron(I, I)));
n3 = kron(I, kron(I, kron(n, I)));
n4 = kron(I, kron(I, kron(I, n)));

I4 = kron(I, kron(I, kron(I, I)));

dim = 2;        % Spatial dimensions
z = dim * 2;    % Coordination number (no. of nearest neighbours)
U = 1;          % Boson-boson interaction energy (within a site)

% On-site terms of the Hamiltonian can be generated early
% (they are independent of J and mu)
H_site = 0.5*U*(n1*(n1-I4) + n2*(n2-I4) + n3*(n3-I4) + n4*(n4-I4));

prec = 2;     % resolution of phase diagram

mu_grid = linspace(0, 2, prec);
J_grid = linspace(0, 0.09, prec);
o_grid = zeros(prec);
n_grid = zeros(prec);

for i = 1:length(mu_grid)
    mu = mu_grid(i);
    
    for j = 1:length(J_grid)
    J = J_grid(j);

    % Randomly generated first guess at order parameters
    for p=1:size(I4)
      F(p,:)=rand+1i*rand;
    end

    o1 = F' * a1 * F;
    o2 = F' * a2 * F;
    o3 = F' * a3 * F;
    o4 = F' * a4 * F;

    % allowing for the order parameter to be complex
    o1c = conj(o1);
    o2c = conj(o2);
    o3c = conj(o3);
    o4c = conj(o4);

    o1_old = Inf;
    o2_old = Inf;
    o3_old = Inf;
    o4_old = Inf;
        
        while (abs(o - o_old) > 1e-2)
            
            o_old = o;            
            
            H = H_site...
                -mu*(n1+n2+n3+n4)...
                -J*(a_d1*a2+a_d2*a1+a_d2*a3+a_d3*a2+a_d3*a4+a_d4*a3+...
                    a_d4*a1+a_d1*a4)...
                -J*(a_d1*o2+a1*o2c-o1c*o2+a_d1*o4+a1*o4c-o1c*o4+...
                    a_d2*o1+a2*o1c-o2c*o1+a_d2*o3+a2*o3c-o2c*o3+...
                    a_d3*o2+a3*o2c-o3c*o2+a_d3*o4+a3*o4c-o3c*o4+...
                    a_d4*o1+a4*o1c-o4c*o1+a_d4*o3+a4*o3c-o4c*o3);
                    
            % Find the groundstate eigenvector gs
            % which has the minimum eigenvalue
            [min_eig, min_eig_idx] = min(eig(H));
            [eigenvectors, ~] = eig(H);
            gs = eigenvectors(:, min_eig_idx);

            % New order parameter given by <gs | a | gs>
            o1 = gs' * a1 * gs;
            o2 = gs' * a2 * gs;
            o3 = gs' * a3 * gs;
            o4 = gs' * a4 * gs;
            
            % Magnitude of overall order parameter
            %o = sqrt(o1^2 + o2^2 + o3^2 + o4^2);
            o = o1;
            % Mean occupancy given by <gs | n | gs>
            %mn = gs' * n * gs;
        end
        o_grid(i, j) = o;
        %n_grid(i, j) = mn;
    end
end
%% Plot changes in the order parameter
surf(J_grid, mu_grid, abs(o_grid), 'EdgeColor','none');
view(0,90);
c = colorbar;
c.Label.String = 'Order Parameter';
xlabel("J/U");
ylabel("\mu/U");
zlabel("Order parameter");
title({'Graph showing how the order parameter varies over the';...
    'chemical potential, \mu, and hopping strength, J.'});