function [uprime, vprime] = FHN_epsilon(epsilon, beta, u, v)
% function [uprime, vprime] = FHN_epsilon(epsilon, beta, u, v)
%
% INPUT
% epsilon       scalar, fixed parameter
% beta          scalar
% u             Fourier_2D
% v             Fourier_2D
%
% OR 
% INPUT
% epsilon       scalar, fixed parameter
% xi            Xi_vector
%
% OUTPUT
% uprime, vprime     Fourier_2D

gamma = 0.5;

if nargin ==2
    xi = beta;
    beta = xi.beta;
    u = xi.u;
    v = xi.v;
elseif ~eq_dim(u, v)
    error('All Fourier_2D should have the same dimension')
end

if ~isa(u, 'Fourier_2D') || ~isa(v, 'Fourier_2D')
    error('The two functions must be handed in a Fourier vectors')
end

% derivatives in Fourier space
Kx = repmat( (-u.nodes_x:u.nodes_x)', 1, 2*u.nodes_y+1);

Ky = repmat( (-u.nodes_y:u.nodes_y), 2*u.nodes_x+1, 1);

% the unit in 2D Fourier space
F2Dunit = 0*Kx;
F2Dunit( u.nodes_x +1, u.nodes_y+1) = 1;

% the FitzHugh-Nagumo right hand side
uprime = Kx.^2 .* u + Ky.^2 .* u + epsilon^-1 *( u + 1/3 * ...
    prod( prod( u, u, 'same' ), u, 'same' ) - v);

vprime = epsilon * ( u + beta * F2Dunit - v * gamma);
