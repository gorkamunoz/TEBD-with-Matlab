function MPS = sweep_left(MPS, N, d, g, J, tau)

for ii = N:-1:2

    %% Hamiltonian 
    
    if ii ==1
        g_prime = g;
    else
        g_prime = 0;
    end    
    
    H = Ising_Hamiltonian(J, g, g_prime); 
    
    exp_operator = expm(-tau*H);
    
    %% Sweep
    
    A = MPS{ii-1};
    B = MPS{ii};
    C = exp_operator;
    
    % Dimensions of the tensors
    d1 = size(A, 2);
    d2 = size(A, 3); % = size(B, 2)
    d3 = size(B, 3);
    [~, min_AB] = min([d2 d3]);
        
    % Contraction of the two sites
    AB = ncon({A, B}, {[-1 -3 1], [-2 1 -4]}, 1);
    
    AB_r = reshape(AB, d^2, size(AB, 3), size(AB, 4)); % The first number is the leg you want to contract with H

    % Contraction with Hamiltonian
    D = ncon({AB_r, C},{[1 -2 -3], [1 -1]}, 1);
    
    D_r = reshape(D, d, d, size(D, 2), size(D, 3));
    D_r = permute(D_r, [1, 3, 2, 4]);
    D_r = reshape(D_r, d*size(D, 2), d*size(D, 3));
    
    % SVD
    [u, s, v] = svd(D_r);
    v = v';    
    s = s/norm(s(:)); % normalizing lambda
    u = u*s;
    
    u = reshape(u, d, size(u, 1)/d, size(u, 2));
    v = reshape(v, size(v, 1), d, size(v, 2)/d);
    v = permute(v, [2, 1, 3]);
    
    % Truncation
    u = u(:, :, 1:d2);
    s = s(1:d2, 1:d2);
    spectrum = diag(s);
    v = v(:, 1:d2, :);
    
    % Contraction for normalization
    new_u = u;
    
    % Updating results
    if ii > 1
        MPS{ii-1} = new_u;
        MPS{ii} = v;
    else
        MPS{ii} = v;
    end
    
end