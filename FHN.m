function uprime = FHN(epsilon, beta, u)
% function uprime = FHN(epsilon, beta, u)
%
% INPUT
% beta          scalar, fixed parameter
% epsilon       scalar
% u             Fourier_2D
%
% OR 
% INPUT
% xi            Xi_vector or small_Xi_vector
%
% OUTPUT
% uprime        Fourier_2D

gamma = 0.5;

if nargin == 1
    if ~isa(epsilon, 'Xi_vector') && ~isa(epsilon, 'small_Xi_vector')
        error('If only one input is given, the input must be a (small_)Xi_vector')
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

