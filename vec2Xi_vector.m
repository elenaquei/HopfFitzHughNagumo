function xi2 = vec2Xi_vector(xi_vec, xi)
% function xi2 = vec2Xi_vector(xi_vec, xi)
%
% trasnforms a vector into a Xi_vector, assuming the same shape as xi
%
% INPUT
% xi_vec        vector
% xi            Xi_vector used for dimensionality
% OUPUT
% xi2           Xi_vector represntation of xi_vec

if length( xi_vec )~= length( xi )
    error('The two lengths must be compatible')
end

u = reshape( xi_vec(2 + (1:length(xi.u))), 2 * xi.u.nodes_x + 1 , 2 * xi.u.nodes_y + 1);
phi = reshape( xi_vec(2 + length(xi.u) + (1:length(xi.phi))), 2 * xi.phi.nodes_x + 1, 2 * xi.phi.nodes_y + 1);

xi2 = Xi_vector( xi_vec(1), xi_vec(2), u, phi);