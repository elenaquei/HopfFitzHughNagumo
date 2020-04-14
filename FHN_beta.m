function uprime = FHN_beta(beta, epsilon, u)
% function uprime = FHN_beta(beta, epsilon, u)
%
% INPUT
% beta          scalar, fixed parameter
% epsilon       scalar
% u             Fourier_2D
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
end

if ~isa(u, 'Fourier_2D')
    error('The function must be handed in a Fourier vectors')
end

% derivatives in Fourier space
Delta = operator_Delta(u);

% constant function 1 in 2D Fourier space
unit_F2D = 0*Delta;
unit_F2D( u.nodes_x + 1, u.nodes_y + 1) = 1;


% the FitzHugh-Nagumo right hand side
uprime = Delta .* u  + epsilon^-1 *( u - 1/3 * ...
    prod( prod( u, u, 'same' ), u, 'same' ) + ...
    gamma ^-1 * u + beta * gamma ^ -1 * unit_F2D);

