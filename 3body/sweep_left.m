function MPS = sweep_left(MPS, N, d, g, J, tau, pbc)

for ii = N:-1:3

    %% Hamiltonian 
    
    if ii ==1
        g_prime = g*(1-pbc);
    else
        g_prime = 0;
    end    
    
    H = TB_Hamiltonian(J, g, g_prime); 
    
    exp_operator = expm(-tau*H);
    
    %% Sweep
    
    A = MPS{ii-2};
    B = MPS{ii-1};
    C = MPS{ii};
    D = exp_operator;
    
    % Dimensions of the tensors
    d1 = size(A, 2);
    d2 = size(B, 2); 
    d3 = size(C, 2);
    d4 = size(C, 3);
        
    % 1 - Contraction of the three sites
    ABC = ncon({A, B, C}, {[-1 -4 1], [-2 1 2], [-3 2 -5]}, [1 2]);
    
    % 2 - Reshape
    ABC = reshape(ABC, d^3, size(ABC, 4), size(ABC, 5));
    
    % 4 - Contraction with Hamiltonian + reshape
    E = ncon({ABC, D},{[1 -2 -3], [1 -1]}, 1);
    
    F = reshape(E, d, d, d, size(E, 2), size(E, 3));
    % 5
    F = permute(F, [1, 4, 2, 3, 5]);
    % 6
    F = reshape(F, d*d*size(E, 2), d*size(E, 3));
    
    
    %% 7 - SVD1:
    [u, s, v] = svd(F);
    v = v';    
    s = s/norm(s(:));      
    
    % 8 and 9 - Reshape+truncate v
%   v = reshape(v, d, size(v, 1), size(v, 2)/d);  % this has shown to be wrong 
    v = reshape(v, size(v, 1), d, size(v, 2)/d);
    v = permute(v, [2, 1, 3]);  
    v = v(:, 1:d3, 1:d4);    
    
    % 10 - contract u and s + truncate u
    u = u*s;
    u = u(1:d*d*d1, 1:d3);
    
    % 11 Reshape u*
    u = reshape(u, d*d1, d*d3);
    
    %% 12 - SVD2:
    [u1, s1, v1] = svd(u);
    v1 = v1';
    s1 = s1/norm(s1(:));     
    
    % 13 - Reshape v1
%   v1 = reshape(v1, d, size(v1, 1), size(v1, 2)/d);   % this has shown to be wrong
    v1 = reshape(v1, size(v1, 1), d, size(v1, 2)/d);
    v1 = permute(v1, [2, 1, 3]);
    
    % 14 - Truncate v1 s1 and u1    
    v1 = v1(:, 1:d2, 1:d3);
    s1 = s1(1:d2, 1:d2);
    u1 = u1(:, 1:d2);
    
    % 15 - Contract u1*s1 + Reshape u1
    u1 = u1*s1;
    u1 = reshape(u1, d, d1, d2);
    
    %% Updating results
    
    MPS{ii} = v;
    MPS{ii-1} = v1;
    MPS{ii-2} = u1;      

    
end