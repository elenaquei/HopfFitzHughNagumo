function Kx2 = operator_Kx_square(u, nodes_y)
% function Kx2 = operator_Kx_square(u, nodes_y)
%
% returns the operator Kx_square
%
% INPUT
% u         either the number of rows nodes or a Fourier_2D element
% nodes_y   if the first input is a number of row nodes, then the second
%           input must be the number of column nodes
%
% OUTPUT
% Kx2       matrix representation of the operator Kx^2

if isa(u, 'Fourier_2D')
    nodes_x = u.nodes_x;
    nodes_y = u.nodes_y;
else
    nodes_x = u;
end

K = (nodes_x:nodes_x)';
K2 = K.^2;
K2rep = repmat(K2,1,2*nodes_y+1);
Kx2 = sparse(diag(K2rep(:)));

