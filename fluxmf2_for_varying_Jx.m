%% 2x2 for phi_ij = Phi * j = pi * j
clear variables;
clc; clearvars;

% Create the bosonic creation and annihilation operators
S = 2;
s = arrayfun(@sqrt, 1:S);
a = sparse(diag(s, 1));         % Annihilation
a_d = sparse(diag(s, -1));      % Creation
n = sparse(diag(0:size(a)-1));  % Number operator = a * a_d
I = speye(size(a));
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
% --------

U = 1;          % atom-atom interaction energy (in same site)
mu = 0.5;       % Middle of first lobe

% Create terms for the Hamiltonian
H1 = 0.5*U*(n1*(n1-I4)+n2*(n2-I4)+n3*(n3-I4)+n4*(n4-I4));
H2 = -mu*(n1+n2+n3+n4);

prec = 100;     % resolution of phase diagram

alpha_grid = linspace(-3, 3, prec);

beta = 1;
gamma = 1;
eta = 1;

points = [];

for i = 1:length(alpha_grid)
    alpha = alpha_grid(i);
    
    J_bounds = [0 1];
    while (range(J_bounds) > 1e-6)
        J = mean(J_bounds);
        
        o = 1e-10; % Small but finite guess of order parameters
        o1 = o;
        o2 = o;
        o3 = o;
        o4 = o;
        
        H_c = -J*(alpha*a_d1*a2+eta*a_d2*a3+beta*a_d3*a4+gamma*a_d4*a1);
        
        % Generate the Hamiltonian matrix from order parameter
        H_coup = -J*(alpha*(a2*o1+a_d1*o2-o1*o2*I4)+...
            gamma*(a1*o4+a_d4*o1-o4*o1*I4)+...
            beta*(a4*o3+a_d3*o4-o3*o4*I4)+...
            eta*(a3*o2+a_d2*o3-o2*o3*I4));
        H = H1 + H2 + H_c + H_c' + H_coup + H_coup';
        
        % Find the groundstate eigenvector gs
        % which has the minimum eigenvalue
        HH = full(H);
        [eigenvectors, eigenvalues] = eig(HH);
        [min_eig, min_eig_idx] = min(diag(eigenvalues));
        gs = eigenvectors(:, min_eig_idx);

        % New order parameters given by <gs | a | gs>
        o1n = gs' * a1 * gs;
        o2n = gs' * a2 * gs;
        o3n = gs' * a3 * gs;
        o4n = gs' * a4 * gs;
        o1nc = conj(o1n);
        o2nc = conj(o2n);
        o3nc = conj(o3n);
        o4nc = conj(o4n);

        % Magnitude of overall order parameter
        o_new = sqrt(o1n*o1nc + o2n*o2nc + o3n*o3nc + o4n*o4nc);
        
        if (o_new < o)          % Mott-Insulating
            J_bounds(1) = J;
        else                    % Superfluid
            J_bounds(2) = J;
        end
    end
    
    points = [points; [alpha, mean(J_bounds)]];
end
%%
plot(points(:, 1), points(:, 2));
xlabel("alpha");
ylabel("J");