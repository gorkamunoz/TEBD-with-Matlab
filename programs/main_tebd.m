clear all; 
%% Parameters

N = 16; % Number of sites
d = 2; % Dimension spins
D = 20; % Max bond dimension
J = 1; % Coupling spins
g = 1; % Transverse field
tau = 1e-2;
t_max = 200;

%% MPS

MPS = initial_mps(N, d, D);

%% Hamiltonian (needs to be defined inside functions)
% 
% sx = [0 1; 1 0];
% sz = [1 0; 0 -1];
% I=eye(2);
% 
% g_prime = 0; % This should be include in the for of the mps sites
% 
% H = -J*(kron(sz,sz)+g*kron(sx,I)+g_prime*kron(I,sx));
% due to kronecker delta inverting the order of the operators, kron(sx, I)
% means Ixsx, so we are applying sx to the second site. This means that we
% need the term in g_prime to get the magnetic in the first site! ie
% g_prime = g for site 1 and 0 for all the others -> the H has to be
% defined inside the MPS
% 
% 
% exp_operator = expm(-tau*H);

%% Evolution 
energy = zeros(1, t_max);
convergence = zeros(1, t_max-1);

for t = 1:t_max
    t
    
% Sweep right
MPS = sweep_right(MPS, N, d, g, J, tau);

% Sweep left
MPS = sweep_left(MPS, N, d, g, J, tau);

% Calculating the energy
energy(t) = exp_value(MPS, g, N, d, J);

% Check convergence
if t ~= 1
    convergence(t-1) = abs(energy(t)-energy(t-1));
end

end

%% Theoretical energy
q_close = ones(1,t_max)*theoretical_ising(J,g)*N; % close boundary conditions

q_open = J*ones(1,t_max)*(1-csc(pi/(2*(2*N+1)))); % open boundary conditions

%% Plots

t = 1:t_max;
figure(1)
subplot(1,2,1)
title('Convergence')
plot(t(2:end), real(convergence))
subplot(1,2,2)
hold on
plot(t, real(energy))
plot(t, q_open)

