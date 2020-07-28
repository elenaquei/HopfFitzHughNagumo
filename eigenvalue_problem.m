function uprime = eigenvalue_problem(xi)
% function uprime = FHN_beta(xi)
%
% INPUT
% xi            Xi_vector
%
% OUTPUT
% uprime        Fourier_2D

gamma = 0.5;

if ~isa(xi, 'Xi_vector')
    error('If only two inputs are given, the second input must be a Xi_vector')
end
epsilon = xi.epsilon;
u = xi.u;
kappa = xi.kappa;
phi = xi.phi;

% derivatives in Fourier space
Delta = operator_Delta(u);


% the FitzHugh-Nagumo right hand side
uprime = Delta .*phi  + epsilon^-1 *( phi - ...
    prod( prod( u, u, 'same' ), phi, 'same' )) + ...
    ( 1 / ( 1i * kappa + epsilon * gamma ) - 1i * kappa ) * phi;
