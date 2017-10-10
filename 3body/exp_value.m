function energy = exp_value(MPS, g, N, d, J, pbc)
% As the last sweep has been from right to left, if we start computing the
% Hamiltonian from the left, all the sites in the right are orthogonal so
% the can be reduce to the identity.
% Also, as the dimension of the physical spin and the bond dimension from
% first to second site is equal, this one also can be reduce to the
% identity.


energy_per_site = zeros(1,N);

for site1 = 1:N
    
    if site1 < N-1
        site2 = site1+1;
        site3 = site1+2;
    elseif site1 == N-1
        site2 = N;
        site3 = 1;
    elseif site1 == N
        site2 = 1;
        site3 = 2;
    end
    
    if pbc == 0 && site1 > N-2 
        break
    end
    
    %% Hamiltonian
    
    if site1 == 1
        g_prime = g*(1-pbc);
    else
        g_prime = 0;
    end
    
    H = TB_Hamiltonian(J, g, g_prime); 
    
    
    %% Expected value    
    
    A = MPS{site1};
    B = MPS{site2};
    C = MPS{site3};
    
    % Contraction of the left part of the MPS     
    D = 1; % for site 1 the left part is one.
    
    if site1 ~= 1
        for kk = 1:site1-1
            D = contraction(D, MPS{kk}, conj(MPS{kk}));
        end
    end
          
    % Contraction of the two sites where we apply H
    
    ABC = ncon({A B C}, {[-1 -4 1] [-2 1 2] [-3 2 -5]}, [1, 2]);
    ABC = reshape(ABC, d^3, size(ABC, 4), size(ABC, 5));
    
    ABC_conj = conj(ABC);
        
    % Right part of the MPS
    R_mps = eye(size(ABC, 3)); 
    
    % Final contraction with Hamiltonian    
    energy_per_site(site1) = ncon({D, ABC, ABC_conj, H, R_mps},...
        {[1 2], [3 1 5], [4 2 6], [3 4], [5 6]},...
        [3 1 5 2 6 4]);
    
end

% Total energy of the chain
energy = sum(energy_per_site);

