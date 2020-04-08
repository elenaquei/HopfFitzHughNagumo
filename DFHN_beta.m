function [ DF,DF_epsilon] = DFHN_beta(beta, epsilon,  u)
% function [ DF, DF_epsilon] = DFHN_beta(beta, epsilon, u)
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
% DF            matrix storing the derivative of D_u FHN_beta
% DF_epsilon    vector storing the derivative of D_epsilon FHN_beta

gamma = 0.5;

if nargin ==2
    xi = epsilon;
    epsilon = xi.epsilon;
    u = xi.u;
end

if ~isa(u, 'Fourier_2D') 
    error('The function must be handed in a Fourier vectors')
end

tot_dim = (2*u.nodes_x+1)*(2*u.nodes_y+1);



DF = operator_Delta( u ) + ...
    epsilon ^ -1 * ( ( 1 + gamma^-1 ) * sparse( eye(tot_dim) )...
    - convMat2D( prod(u, u, 'same') ));

DF_epsilon_mat = - epsilon^-2 * ( u + 1/3 * ...
    prod( prod( u, u, 'same' ), u, 'same' ) + ...
    gamma ^-1 * u + beta * gamma ^ -1 );
 

DF_epsilon = DF_epsilon_mat(:);
