

% Checking right tensors normalization
N = 10; % Number of sites
d = 2; % Dimension spins
D = 15; % Max bond dimension

% MPS = initial_mps(N, d, D);

site = 4;


A = s;
B = conj(A);


res = ncon({A, B}, {[1 2 3] [1 2 3]}, [1 2 3])
norm(res(:))

%%
A = ans;
res = ncon({A,conj(A)}, {[1 -1 2], [1 -2 2]}, [1 2])
norm(res(:))

res = ncon({C,conj(C)}, {[1 2], [1 2]}, [1 2])
