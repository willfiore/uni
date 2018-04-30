clear variables;

A = delsq(numgrid('C',256));
[V, D] = eigs(A, 2, 'sm');