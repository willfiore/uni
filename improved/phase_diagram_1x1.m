clear variables;
clc; close all;

% Create the bosonic creation and annihilation operators
S = 4;
s = arrayfun(@sqrt, 1:S);
a = diag(s, 1);         % Annihilation
a_d = diag(s, -1);      % Creation
n = diag(0:size(a)-1);  % Number operator = a * a_d
% --------

U = 1;          % atom-atom interaction energy (in same site)

prec = 1000;
mu_grid = linspace(0, S-1, prec);
points = [];

for i = 1:length(mu_grid)
    mu = mu_grid(i);
    
    J_bounds = [0 1];
    
    while (range(J_bounds) > 1e-10)
        J = mean(J_bounds);
        
        o = 1e-10; % Small but finite guess of order parameter

        % Generate the Hamiltonian matrix from order parameter
        H = -J*(o*a + o*a_d) + 0.5*U*n.*(n-1) - mu*n;

        % Find the groundstate eigenvector gs
        % which has the minimum eigenvalue
        [min_eig, min_eig_idx] = min(eig(H));
        [eigenvectors, ~] = eig(H);
        gs = eigenvectors(:, min_eig_idx);

        % New order parameter given by <gs | a | gs>
        o_new = gs' * a * gs;
        
        if (o_new < o)      % Mott-Insulating
            J_bounds(1) = J;
        else                % Superfluid
            J_bounds(2) = J;
        end
    end
    
    points = [points; [mean(J_bounds), mu]];
    
    fprintf("%2.1f%%\n", 100*i/prec);
end

%%
plot(points(:, 1), points(:, 2));
xlabel("J");
ylabel("mu");