function [ DF,DF_epsilon] = DFHN_beta(epsilon, beta, u)
% function [ DF, DF_epsilon] = DFHN_beta(epsilon, beta, u)
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

if nargin == 1
    beta = epsilon.beta;
    u = epsilon.u;
    epsilon = epsilon.epsilon;
end

tot_dim = (2*u.nodes_x+1)*(2*u.nodes_y+1);

Delta = derivative_operator_Delta( u );

unit_F2D = 0*u.vector;
unit_F2D( u.nodes_x + 1, u.nodes_y + 1) = 1;

DF =  Delta + ...
    epsilon ^ -1 * ( ( 1 + gamma^-1 ) * sparse( eye(tot_dim) )...
    - convMat2D( prod(u, u, 'same') ,'same') );

DF_epsilon_mat = - epsilon^-2 * ( u - 1/3 * ...
    prod( prod( u, u, 'same' ), u, 'same' ) + ...
    gamma ^-1 * u + beta * gamma ^ -1 * unit_F2D);
 

DF_epsilon = DF_epsilon_mat.vector(:);
