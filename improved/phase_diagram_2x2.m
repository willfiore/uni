clear variables;
clc; close all;

% Create the bosonic creation and annihilation operators
S = 2;
s = arrayfun(@sqrt, 1:S);
a = sparse(diag(s, 1));         % Annihilation
a_d = sparse(diag(s, -1));      % Creation
n = sparse(diag(0:size(a)-1));  % Number operator = a * a_d
I = speye(size(a));       % Identity operator
% --------

a1 = kron(a, kron(I, kron(I, I)));      % a_1
a2 = kron(I, kron(a, kron(I, I)));      % a_2
a3 = kron(I, kron(I, kron(a, I)));      % a_3
a4 = kron(I, kron(I, kron(I, a)));      % a_4
a_d1 = kron(a_d, kron(I, kron(I, I)));  % adagger_1
a_d2 = kron(I, kron(a_d, kron(I, I)));  % adagger_2
a_d3 = kron(I, kron(I, kron(a_d, I)));  % adagger_3
a_d4 = kron(I, kron(I, kron(I, a_d)));  % adagger_4
n1 = kron(n, kron(I, kron(I, I)));      % n_1
n2 = kron(I, kron(n, kron(I, I)));      % n_2
n3 = kron(I, kron(I, kron(n, I)));      % n_3
n4 = kron(I, kron(I, kron(I, n)));      % n_4
I4 = kron(I, kron(I, kron(I, I)));

% On-site interaction
U = 1;

prec = 200;
mu_grid = linspace(0, S-1, prec);
points = [];

% Create terms of Hamiltonian independent of mu and J
H1 = 0.5*U*(n1*(n1-I4) + n2*(n2-I4) + n3*(n3-I4) + n4*(n4-I4));

for i = 1:length(mu_grid)
    mu = mu_grid(i);
    
    J_bounds = [0 1];
    
    while (range(J_bounds) > 1e-6)
        J = mean(J_bounds);
        
        o = 1e-10;

        % Generate the Hamiltonian matrix from order parameter
        H2 = -mu*(n1+n2+n3+n4);
        
        H = H1 + H2...
            -J*(a_d1*a2+a_d2*a1+a_d2*a3+a_d3*a2+a_d3*a4+a_d4*a3+...
            a_d4*a1+a_d1*a4)...
            -J*(a_d1*o+a1*o-o*o+a_d1*o+a1*o-o*o+...
            a_d2*o+a2*o-o*o+a_d2*o+a2*o-o*o+...
            a_d3*o+a3*o-o*o+a_d3*o+a3*o-o*o+...
            a_d4*o+a4*o-o*o+a_d4*o+a4*o-o*o);

        % Find the groundstate eigenvector gs
        % which has the minimum eigenvalue
        [eigenvectors, eigenvalues] = eig(H);
        [min_eig, min_eig_idx] = min(diag(eigenvalues));
        gs = eigenvectors(:, min_eig_idx);

        % New order parameters given by <gs | a | gs>
        o_new = gs' * a1 * gs;
        
        if (o_new < o)          % Mott-Insulating
            J_bounds(1) = J;
        else                    % Superfluid
            J_bounds(2) = J;
        end
    end
    
    points = [points; [mean(J_bounds), mu]];
    
    fprintf("%2.0f%%\n", 100*i/prec);
end

%%
plot(points(:, 1), points(:, 2));
xlabel("J");
ylabel("mu");