clear all; %close all;
%% Parameters

N = 10; % Number of sites
d = 2; % Dimension spins
D = 32; % Max bond dimension
J = 1; % Coupling spins
g = 1; % Transverse field
tau = 1e-2; % Imaginary time
t_max = 500; % Number of iterations

%% MPS

MPS = initial_mps(N, d, D);

MPSnpbc = initial_mps_nopbc(N, d, D);

%% Boundary conditions

pbc = 0; % PeriodicBC = 1 ; OpenBC = 0;

%% Theoretical results

% Periodic boundary
q_periodic = N*ones(1, t_max)*integral(@(x)(-1/(4*pi))*2*J*sqrt(1+g^2-2*g*cos(x)),-pi,pi);
% Close boundary
q_open = J*ones(1,t_max)*(1-csc(pi/(2*(2*N+1)))); % open boundary conditions

if pbc == 1
    q_theory = q_periodic;
else
    q_theory = q_open;
end

%% Evolution
energy = zeros(1, t_max);
convergence = zeros(1, t_max-1);

for t = 1:t_max
    t
    
    % Sweep right
    MPS = sweep_right(MPS, N, d, g, J, tau, pbc);
    % PBC
    if pbc == 1
        MPS = PBC(MPS, N, d, g, J, tau);
    end
    
    % Sweep left
    MPS = sweep_left(MPS, N, d, g, J, tau, pbc);
    
    % PBC
    if pbc == 1
        MPS = PBC(MPS, N, d, g, J, tau);
    end
    
    % Calculating the energy
    energy(t) = exp_value(MPS, g, N, d, J, pbc);
    % Check convergence
    if t ~= 1
        convergence(t) = abs(energy(t)-q_theory(1));
    end
    
end

%% Plots
t = 1:t_max;



figure(1)
subplot(1,2,1)
title('Convergence')
hold on
set(gca,'XScale','log')
set(gca,'YScale','log')
plot(t, real(convergence))

subplot(1,2,2)
hold on
plot(t, real(energy))
plot(t, q_theory, '--')

