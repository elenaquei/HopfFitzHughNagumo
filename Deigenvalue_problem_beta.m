function [ Depsilon, Du, Dkappa, Dphi ] = Deigenvalue_problem_beta(xi)
% function [ Depsilon, Du, Dkappa, Dphi ] = DFHN_beta(xi)
%
% INPUT
% xi            Xi_vector
%
% OUTPUT
% Depsilon      vector storing the derivative of D_epsilon eigenvalue
%               problem
% Du            matrix storing the derivative of D_u eigenvalue prolem
% Dkappa        vector storing the derivative of D_kappa eigenvalue prolem
% Dphi          matrix storing the derivative of D_phi eigenvalue prolem

gamma = 0.5;


if ~isa(xi,'Xi_vector')
    error('Input must be a Xi_vector')
end
epsilon = xi.epsilon;
u = xi.u;
kappa = xi.kappa;
phi = xi.phi;


tot_dim = (2*u.nodes_x+1)*(2*u.nodes_y+1);

Delta = derivative_operator_Delta( u );

Depsilon_mat = - epsilon^-2 *(phi  - prod(u,prod(u,phi,'same'),'same') ) + ...
    gamma / ( 1i * kappa + epsilon * gamma ) ^2 * phi;

Dkappa_mat = - 1i * ( 1 + 1/(1i * kappa + epsilon * gamma ) ^2 ) * phi;

Du = - 2 / epsilon * convMat2D( prod(u, phi, 'same') );

Dphi =  Delta + ...
    ( epsilon ^ -1 + 1 / ( 1i * kappa + epsilon * gamma ) - 1i * kappa ) *  sparse( eye(tot_dim) )...
    - epsilon ^ -1 * convMat2D( prod(u, u, 'same') ,'same');

Depsilon = Depsilon_mat.vector(:);

Dkappa = Dkappa_mat.vector(:);
