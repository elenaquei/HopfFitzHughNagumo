function uprime = eigenvalue_problem(~, epsilon, u, kappa, phi)
% function uprime = FHN_beta(beta, epsilon, u, kappa, phi)
%
% INPUT
% beta          scalar, fixed parameter
% epsilon       scalar
% u             Fourier_2D
% kappa         scalar
% phi           Fourier_2D
%
% OR 
% INPUT
% beta          scalar, fixed parameter
% xi            Xi_vector
%
% OUTPUT
% uprime        Fourier_2D

gamma = 0.5;

if nargin ==2
    if ~isa(epsilon, 'Xi_vector')
        error('If only two inputs are given, the second input must be a Xi_vector')
    end
    xi = epsilon;
    epsilon = xi.epsilon;
    u = xi.u;
    kappa = xi.kappa;
    phi = xi.phi;
end

if ~isa(u, 'Fourier_2D') || ~isa(phi, 'Fourier_2D')
    error('The functions must be handed in a Fourier vectors')
end

if u.nodes_x ~= phi.nodes_x || u.nodes_y ~= phi.nodes_y
    error('Inputed Fourier coefficient sizes must be equal')
end

% derivatives in Fourier space
Delta = operator_Delta(u);


% the FitzHugh-Nagumo right hand side
uprime = Delta .*phi  + epsilon^-1 *( phi - ...
    prod( prod( u, u, 'same' ), phi, 'same' )) + ...
    ( 1 / ( 1i * kappa + epsilon * gamma ) - 1i * kappa ) * phi;
