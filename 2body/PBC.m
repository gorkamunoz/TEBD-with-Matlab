function MPS = PBC(MPS, N, d, g, J, tau)

%% Hamiltonian

H = Ising_Hamiltonian(J, 0, g);

expop = expm(-tau*H);

expop = reshape(expop, d, d, d, d);
expop = permute(expop, [1, 3, 2, 4]);
expop = reshape(expop, d^2, d^2);

[u, s, v] = svd(expop);
s = sqrt(s);
v= v';

u = u*s;
v = ncon({s, v}, {[-1 1], [1 -2]}, 1);

u = reshape(u, d, d, d^2);
v = reshape(v, d^2, d, d);
v = permute(v, [3, 2, 1]);


%% Contraction 

for site = 1:N
    
    A = MPS{site};
    
    if site == 1
        B = u;
        new_A = ncon({A, B}, {[1 -4 -1], [-3 1 -2]}, 1);     
        new_A = reshape(new_A, size(new_A, 1)*size(new_A, 2),...
            size(new_A, 3), size(new_A, 4));
        new_A = permute(new_A, [2, 3, 1]);
    elseif site == N
        B = v;
        new_A = ncon({A, B}, {[1 -1 -4], [-3 1 -2]}, 1);
        new_A = reshape(new_A, size(new_A, 1)*size(new_A, 2),...
            size(new_A, 3), size(new_A, 4));
        new_A = permute(new_A, [2, 1, 3]);
    else
        
        Ar = reshape(A, 1, size(A,1)*size(A, 2)*size(A, 3));
        B = reshape(eye(4), 1, 16);

        new_A = kron(Ar, B);
        new_A = reshape(new_A, 4, 4, size(A, 1), size(A, 2), size(A, 3));
        new_A = permute(new_A, [3, 1, 4, 2, 5]);
        new_A = reshape(new_A, [], size(new_A, 2)*size(new_A, 3), size(new_A, 4)*size(new_A, 5));
    end

    trunc_A = new_A(:, 1:size(A,2), 1:size(A,3));
    
   
    new_MPS{site} = trunc_A;
    
end

MPS = new_MPS;