function Ky2 = operator_Ky_square(u, nodes_y)
% function Ky2 = operator_Ky_square(u, nodes_y)
%
% returns the operator Ky_square
%
% INPUT
% u         either the number of rows nodes or a Fourier_2D element
% nodes_y   if the first input is a number of row nodes, then the second
%           input must be the number of column nodes
%
% OUTPUT
% Ky2       matrix representation of the operator Ky^2

if isa(u, 'Fourier_2D')
    nodes_x = u.nodes_x;
    nodes_y = u.nodes_y;
else
    nodes_x = u;
end

K = (nodes_y:nodes_y);
K2 = K.^2;
K2rep = repmat(K2,2*nodes_x+1,1);
Ky2 = sparse(diag(K2rep(:)));

