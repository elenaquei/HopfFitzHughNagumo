function [ Depsilon, Du, Dkappa, Dphi ] = Deigenvalue_problem_beta(~, epsilon,  u, kappa, phi)
% function [ Depsilon, Du, Dkappa, Dphi ] = DFHN_beta(beta, epsilon, u, kappa, phi)
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
% Depsilon      vector storing the derivative of D_epsilon eigenvalue
%               problem
% Du            matrix storing the derivative of D_u eigenvalue prolem
% Dkappa        vector storing the derivative of D_kappa eigenvalue prolem
% Dphi          matrix storing the derivative of D_phi eigenvalue prolem

gamma = 0.5;

if nargin ==2
    if ~isa(epsilon,'Xi_vector')
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
