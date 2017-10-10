function MPS = sweep_right(MPS, N, d, g, J, tau, pbc)


for ii = 1:N-2
    
    %% Hamiltonian 
    
    if ii ==1
        g_prime = g*(1-pbc);
    else
        g_prime = 0;
    end 
    
    H = TB_Hamiltonian(J, g, g_prime); 
    
    exp_operator = expm(-tau*H);
    
    %% Sweep
    
    % 
    A = MPS{ii};
    B = MPS{ii+1};
    C = MPS{ii+2};
    D = exp_operator;
    
    % Dimensions of the tensors
    d1 = size(A, 2);
    d2 = size(B, 2); 
    d3 = size(C, 2);
    d4 = size(C, 3);
    
    % 1 - Contraction of the three sites
    ABC = ncon({A, B, C}, {[-1 -4 1], [-2 1 2], [-3 2 -5]}, [1 2]);
    
    % 2 - Reshape
    ABC = reshape(ABC, d^3, d1, d4); % The first number is the leg you want to contract with H

    % 4 - Contraction with Hamiltonian + reshape
    E = ncon({ABC, D},{[1 -2 -3], [1 -1]}, 1);
    
    F = reshape(E, d, d, d, d1, d4);
    % 5
    F = permute(F, [1, 4, 2, 3, 5]);
    % 6
    F = reshape(F, d*d1, d*d*d4);
    
    %% 7 - SVD1: 
    [u, s, v] = svd(F);
    v = v';   
    s = s/norm(s(:)); % normalizing lambda
    
    % 8 - reshape + truncate first site
    u = reshape(u, d, size(u, 1)/d, size(u, 2));
    % 9 - truncate
    u = u(:, 1:d1, 1:d2); 
    
    % New 10a - Reshape + truncate v and s
    v = reshape(v, size(v, 1), d, size(v, 2)/d);
    v = permute(v, [2, 1, 3]);
    v = v(:, 1:d2, :);
    s = s(1:d2, 1:d2);       
    % New 10b - contract s and v
    v = ncon({s, v}, {[-2 1], [-1 1 -3]}, 1);

    % 11 - reshape new tensor
    v = reshape(v, d*d2, d*d4);
    
    %% 12 - SVD2
    [u1, s1, v1] = svd(v);
    v1 = v1';
    s1 = s1/norm(s1(:));
    
    % 13 - reshape
    u1 = reshape(u1, d, size(u1, 1)/d, size(u1, 2));
    u1 = u1(:, 1:d2, 1:d3); 
    
    % New 14a - Reshape + truncate s1 and v1
    v1 = reshape(v1, size(v1, 1), d, size(v1, 2)/d);
    v1 = permute(v1, [2, 1, 3]);
    v1 = v1(:, 1:d3, :);
    s1 = s1(1:d3, 1:d3);      
    % New 14b - contract s and v
    v1 = ncon({s1, v1}, {[-2 1], [-1 1 -3]}, 1);
    
    %% Updating results
    MPS{ii} = u;
    MPS{ii+1} = u1;
    MPS{ii+2} = v1;
    
    
end

%% Changes

    % Old 9 cont. - truncate
%     s = s(1:d2, 1:d2);
%     v = v(1:d2,:);
    % Old10 - contract s and v
%     v = ncon({s, v}, {[-1 1], [1 -2]}, 1);

    % Old 14 - Truncate
%     s1 = s1(1:d3, 1:d3);
%     v1 = v1(1:d3,:);
    
    % Old 15 - Contract + reshape
%     v1 = ncon({s1, v1}, {[-1 1], [1 -2]}, 1);
%     v1 = reshape(v1, d, d3, d4);    