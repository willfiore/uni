%% for 3x3 for phi_ij = Phi * j = 2pi/3 * j
%clear variables;
clc; clearvars;
rng('shuffle')

% Create the bosonic creation and annihilation operators
S = 2;
s = arrayfun(@sqrt, 1:S);
a = sparse(diag(s, 1));         % Annihilation
a_d = sparse(diag(s, -1));      % Creation
n = sparse(diag(0:size(a)-1));  % Number operator = a * a_d
I = speye(size(a));             % Identity operator
% ------------------------------------------------------

a1 = kron(a, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));      % a_1
a2 = kron(I, kron(a, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));      % a_2
a3 = kron(I, kron(I, kron(a, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));      % a_3
a4 = kron(I, kron(I, kron(I, kron(a, kron(I, kron(I, kron(I, kron(I, I))))))));      % a_4
a5 = kron(I, kron(I, kron(I, kron(I, kron(a, kron(I, kron(I, kron(I, I))))))));      % a_5
a6 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(a, kron(I, kron(I, I))))))));      % a_6
a7 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(a, kron(I, I))))))));      % a_7
a8 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(a, I))))))));      % a_8
a9 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, a))))))));      % a_9
a_d1 = kron(a_d, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));  % adagger_1
a_d2 = kron(I, kron(a_d, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));  % adagger_2
a_d3 = kron(I, kron(I, kron(a_d, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));  % adagger_3
a_d4 = kron(I, kron(I, kron(I, kron(a_d, kron(I, kron(I, kron(I, kron(I, I))))))));  % adagger_4
a_d5 = kron(I, kron(I, kron(I, kron(I, kron(a_d, kron(I, kron(I, kron(I, I))))))));  % adagger_5
a_d6 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(a_d, kron(I, kron(I, I))))))));  % adagger_6
a_d7 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(a_d, kron(I, I))))))));  % adagger_7
a_d8 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(a_d, I))))))));  % adagger_8
a_d9 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, a_d))))))));  % adagger_9
n1 = kron(n, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));      % n_1
n2 = kron(I, kron(n, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));      % n_2
n3 = kron(I, kron(I, kron(n, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));      % n_3
n4 = kron(I, kron(I, kron(I, kron(n, kron(I, kron(I, kron(I, kron(I, I))))))));      % n_4
n5 = kron(I, kron(I, kron(I, kron(I, kron(n, kron(I, kron(I, kron(I, I))))))));      % n_5
n6 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(n, kron(I, kron(I, I))))))));      % n_6
n7 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(n, kron(I, I))))))));      % n_7
n8 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(n, I))))))));      % n_8
n9 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, n))))))));      % n_9
I9 = kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, kron(I, I))))))));
% --------

dim = 2;        % number of spatial dimensions
z = dim * 2;    % number of nearest neighbours (coordination number)
U = 1;          % atom-atom interaction energy (in same site)

% Position within phase diagram
mu = 0.03;
J = 0.1;

% Create terms for the Hamiltonian
H1 = 0.5*U*(n1*(n1-I9)+n2*(n2-I9)+n3*(n3-I9)+n4*(n4-I9)+n5*(n5-I9)...
            +n6*(n6-I9)+n7*(n7-I9)+n8*(n8-I9)+n9*(n9-I9))...
     -mu*(n1+n2+n3+n4+n5+n6+n7+n8+n9)...
     -J*(a_d1*a2+a_d2*a1+a_d2*a3+a_d3*a2+a_d4*a5+a_d5*a4+a_d5*a6+a_d6*a5+...
         a_d7*a8+a_d8*a7+a_d8*a9+a_d9*a8+...
         exp((1i*6*pi)/3)*a_d3*a4+exp(-(1i*6*pi)/3)*a_d4*a3+...
         exp(-(1i*2*pi)/3)*a_d6*a1+exp((1i*2*pi)/3)*a_d1*a6+...
         exp((1i*4*pi)/3)*a_d2*a5+exp(-(1i*4*pi)/3)*a_d5*a2+...
         exp((1i*6*pi)/3)*a_d4*a9+exp(-(1i*6*pi)/3)*a_d9*a4+...
         exp(-(1i*2*pi)/3)*a_d7*a6+exp((1i*2*pi)/3)*a_d6*a7+...
         exp((1i*4*pi)/3)*a_d5*a8+exp(-(1i*4*pi)/3)*a_d8*a5);
         
     
% randomly generated first guess' at order parameters
for p=1:size(I9)
  F(p,:)=rand+1i*rand;
end
o1 = F' * a1 * F;
o2 = F' * a2 * F;
o3 = F' * a3 * F;
o4 = F' * a4 * F;
o5 = F' * a5 * F;
o6 = F' * a6 * F;
o7 = F' * a7 * F;
o8 = F' * a8 * F;
o9 = F' * a9 * F;
% allowing for the order parameter to be complex
o1c = conj(o1);
o2c = conj(o2);
o3c = conj(o3);
o4c = conj(o4);
o5c = conj(o5);
o6c = conj(o6);
o7c = conj(o7);
o8c = conj(o8);
o9c = conj(o9);

delta = 1;
min_eig = 1;

while abs(delta) > 1e-7
            
min_eig_old = min_eig;
o1_old = o1;
o2_old = o2;
o3_old = o3;
o4_old = o4;
o5_old = o5;
o6_old = o6;
o7_old = o7;
o8_old = o8;
o9_old = o9;




% Generate the Hamiltonian matrix from order parameter
 H = H1...
                -J*(a_d1*o3+a1*o3c+a_d3*o1+a3*o1c+...
                    a_d4*o6+a4*o6c+a_d6*o4+a6*o4c+...
                    a_d9*o7+a9*o7c+a_d7*o9+a7*o9c+...
                    exp(-(1i*2*pi)/3)*(a_d1*o7+a7*o1c)+...
                    exp(-(1i*4*pi)/3)*(a_d2*o8+a8*o2c)+exp(-(1i*6*pi)/3)*(a_d3*o9+a9*o3c)+...
                    exp((1i*6*pi)/3)*(a_d9*o3+a3*o9c)+exp((1i*4*pi)/3)*(a_d8*o2+a2*o8c)+...
                    exp((1i*2*pi)/3)*(a_d7*o1+a1*o7c));

HH = full(H);
            % Find the groundstate eigenvector gs
            % which has the minimum eigenvalue
            [min_eig, min_eig_idx] = min(eigs(HH));
            [eigenvectors, eigenvalues] = eig(HH);
            gs = eigenvectors(:, min_eig_idx);
%             [gs, eigenvalue] = eigs(Hinv, 1, 'lm');
            % New order parameters given by <gs | a | gs>
            o1 = gs' * a1 * gs;
            o2 = gs' * a2 * gs;
            o3 = gs' * a3 * gs;
            o4 = gs' * a4 * gs;
            o5 = gs' * a5 * gs;
            o6 = gs' * a6 * gs;
            o7 = gs' * a7 * gs;
            o8 = gs' * a8 * gs;
            o9 = gs' * a9 * gs;
            o1c = conj(o1);
            o2c = conj(o2);
            o3c = conj(o3);
            o4c = conj(o4);
            o5c = conj(o5);
            o6c = conj(o6);
            o7c = conj(o7);
            o8c = conj(o8);
            o9c = conj(o9);
            
            % Magnitude of overall order parameter
            o = sqrt(o1*o1c + o2*o2c + o3*o3c + o4*o4c + o5*o5c + o6*o6c);
            
            delta = (min_eig - min_eig_old)/(min_eig);
            
            disp([' delta ', num2str(delta)]);
            pause(0.05);
end
