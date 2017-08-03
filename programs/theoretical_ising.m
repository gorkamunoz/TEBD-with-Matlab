function q = theoretical_ising(J, g)

ising = @(x)(-1/(4*pi))*2*J*sqrt(1+g^2-2*g*cos(x));
% ising = @(x)(-1/(4*pi))*(J/2)*sqrt(1/4+g^2-g*cos(x));

q = integral(ising,-pi,pi);