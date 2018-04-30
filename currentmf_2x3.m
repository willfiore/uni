%% for 3x2 for phi_ij = Phi * j = 2pi/3 * j
%clear variables;
clc; clearvars;
rng('shuffle')

% Create the bosonic creation and annihilation operators
S = 3;
s = arrayfun(@sqrt, 1:S);
a = sparse(diag(s, 1));         % Annihilation
a_d = sparse(diag(s, -1));      % Creation
n = sparse(diag(0:size(a)-1));  % Number operator = a * a_d
I = speye(size(a));
a1 = kron(a, kron(I, kron(I, kron(I, kron(I, I)))));      % a_1
a2 = kron(I, kron(a, kron(I, kron(I, kron(I, I)))));      % a_2
a3 = kron(I, kron(I, kron(a, kron(I, kron(I, I)))));      % a_3
a4 = kron(I, kron(I, kron(I, kron(a, kron(I, I)))));      % a_4
a5 = kron(I, kron(I, kron(I, kron(I, kron(a, I)))));      % a_5
a6 = kron(I, kron(I, kron(I, kron(I, kron(I, a)))));      % a_6
a_d1 = kron(a_d, kron(I, kron(I, kron(I, kron(I, I)))));  % adagger_1
a_d2 = kron(I, kron(a_d, kron(I, kron(I, kron(I, I)))));  % adagger_2
a_d3 = kron(I, kron(I, kron(a_d, kron(I, kron(I, I)))));  % adagger_3
a_d4 = kron(I, kron(I, kron(I, kron(a_d, kron(I, I)))));  % adagger_4
a_d5 = kron(I, kron(I, kron(I, kron(I, kron(a_d, I)))));  % adagger_5
a_d6 = kron(I, kron(I, kron(I, kron(I, kron(I, a_d)))));  % adagger_6
n1 = kron(n, kron(I, kron(I, kron(I, kron(I, I)))));      % n_1
n2 = kron(I, kron(n, kron(I, kron(I, kron(I, I)))));      % n_2
n3 = kron(I, kron(I, kron(n, kron(I, kron(I, I)))));      % n_3
n4 = kron(I, kron(I, kron(I, kron(n, kron(I, I)))));      % n_4
n5 = kron(I, kron(I, kron(I, kron(I, kron(n, I)))));      % n_5
n6 = kron(I, kron(I, kron(I, kron(I, kron(I, n)))));      % n_6
I6 = kron(I, kron(I, kron(I, kron(I, kron(I, I)))));
% --------

dim = 2;        % number of spatial dimensions
z = dim * 2;    % number of nearest neighbours (coordination number)
U = 1;          % atom-atom interaction energy (in same site)

% Position within phase diagram
mu = 1;
J = 2;

% Create terms for the Hamiltonian
H1 = 0.5*U*(n1*(n1-I6)+n2*(n2-I6)+n3*(n3-I6)+n4*(n4-I6)+n5*(n5-I6)...
            +n6*(n6-I6))...
     -mu*(n1+n2+n3+n4+n5+n6)...
     -J*(a_d1*a2+a_d2*a1+a_d2*a3+a_d3*a2+exp((1i*6*pi)/3)*a_d3*a4+exp(-(1i*6*pi)/3)*a_d4*a3+...
         a_d4*a5+a_d5*a4+a_d5*a6+a_d6*a5+exp(-(1i*2*pi)/3)*a_d6*a1+exp((1i*2*pi)/3)*a_d1*a6+...
         exp((1i*4*pi)/3)*a_d2*a5+exp(-(1i*4*pi)/3)*a_d5*a2);
     
% randomly generated first guess' at order parameters
for p=1:size(I6)
  F(p,:)=rand+1i*rand;
end
o1 = F' * a1 * F;
o2 = F' * a2 * F;
o3 = F' * a3 * F;
o4 = F' * a4 * F;
o5 = F' * a5 * F;
o6 = F' * a6 * F;
% allowing for the order parameter to be complex
o1c = conj(o1);
o2c = conj(o2);
o3c = conj(o3);
o4c = conj(o4);
o5c = conj(o5);
o6c = conj(o6);

o1_old = Inf;
o2_old = Inf;
o3_old = Inf;
o4_old = Inf;
o5_old = Inf;
o6_old = Inf;


while (abs(o1 - o1_old) > 1e-4) && (abs(o2 - o2_old) > 1e-4) && (abs(o3 - o3_old) > 1e-4) && (abs(o4 - o4_old) > 1e-4) && (abs(o5 - o5_old) > 1e-4) && (abs(o6 - o6_old) > 1e-4)
            
o1_old = o1;
o2_old = o2;
o3_old = o4;
o4_old = o4;
o5_old = o5;
o6_old = o6;
            
% Generate the Hamiltonian matrix from order parameter
 H = H1...
                -J*(a_d1*o3+a1*o3c-o1c*o3*I6+exp(-(1i*2*pi)/3)*(a_d1*o6+a6*o1c-o1c*o6*I6)+...
                    exp(-(1i*4*pi)/3)*(a_d2*o5+a5*o2c-o2c*o5*I6)+exp(-(1i*6*pi)/3)*(a_d3*o4+a4*o3c-o3c*o4*I6)+...
                    a_d3*o1+a3*o1c-o3c*o1*I6+a_d4*o6+a4*o6c-o4c*o6*I6+...
                    exp((1i*6*pi)/3)*(a_d4*o3+a3*o4c-o4c*o3*I6)+exp((1i*4*pi)/3)*(a_d5*o2+a2*o5c-o5c*o2*I6)+...
                    exp((1i*2*pi)/3)*(a_d6*o1+a1*o6c-o6c*o1*I6)+a_d6*o4+a6*o4c-o6c*o4*I6);


            % Find the groundstate eigenvector gs
            % which has the minimum eigenvalue
            %[min_eig, min_eig_idx] = eigs(H,'smallestabs');
            %[eigenvectors, eigenvalues] = eig(H);
            %gs = eigenvectors(:, min_eig_idx);
            [gs, eigenvalue] = eigs(H, 1, 'sm');
            % New order parameters given by <gs | a | gs>
            o1 = gs' * a1 * gs;
            o2 = gs' * a2 * gs;
            o3 = gs' * a3 * gs;
            o4 = gs' * a4 * gs;
            o5 = gs' * a5 * gs;
            o6 = gs' * a6 * gs;
            o1c = conj(o1);
            o2c = conj(o2);
            o3c = conj(o3);
            o4c = conj(o4);
            o5c = conj(o5);
            o6c = conj(o6);
            
            % Magnitude of overall order parameter
            o = sqrt(o1*o1c + o2*o2c + o3*o3c + o4*o4c + o5*o5c + o6*o6c);
            disp("o1", num2str(o1-o1_old), "o2", num2str(o2-o2_old),...
                "o3", num2str(o3-o3_old), "o4", num2str(o4-o4_old),...
                "o5", num2str(o5-o5_old), "o6", num2str(o6-o6_old));
            pause(0.05);
end
