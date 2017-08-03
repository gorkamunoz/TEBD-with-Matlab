function H = Ising_Hamiltonian(J, g, g_prime)

sx = [0 1; 1 0];
sz = [1 0; 0 -1];
I=eye(2);

H = -J*(kron(sz,sz)+g*kron(sx,I)+g_prime*kron(I,sx));
% due to kronecker delta inverting the order of the operators, kron(sx, I)
    % means Ixsx, so we are applying sx to the second site. This means that we
    % need the term in g_prime to get the magnetic in the first site! ie
    % g_prime = g for site 1 and 0 for all the others -> the H has to be
    % defined inside the MPS