function energy = exp_value_nopbc(MPS, g, N, d, J)
% As the last sweep has been from right to left, if we start computing the
% Hamiltonian from the left, all the sites in the right are orthogonal so
% the can be reduce to the identity.
% Also, as the dimension of the physical spin and the bond dimension from
% first to second site is equal, this one also can be reduce to the
% identity.


energy_per_site = zeros(1,N-1);

for site1 = 1:N-1
    site2 = site1+1;
    
    %% Hamiltonian
    
    if site1 == 1
        g_prime = g;
    else
        g_prime = 0;
    end
    
    H = Ising_Hamiltonian(J, g, g_prime); 
    
    
    %% Expected value    
    
    % Contraction of the left part of the MPS     
    C = 1; % for site 1 the left part is one.
    
    if site1 ~= 1
        for kk = 1:site1-1
            C = contraction(C, MPS{kk}, conj(MPS{kk}));
        end
    end
          
    % Contraction of the two sites where we apply H
    AB = ncon({MPS{site1} MPS{site2}}, {[-1 -3 1] [-2 1 -4]}, 1);
    AB = reshape(AB, d^2, size(AB, 3), size(AB, 4));
    
    AB_conj = conj(AB);
    
    % Right part of the MPS
    R_mps = eye(size(AB, 3)); 
    
    % Final contraction with Hamiltonian    
    energy_per_site(site1) = ncon({C, AB, AB_conj, H, R_mps},...
        {[1 2], [3 1 5], [4 2 6], [3 4], [5 6]},...
        [3 1 5 2 6 4]);
    
end

% Total energy of the chain
energy = sum(energy_per_site);

