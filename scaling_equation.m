function [ scal, D_phi ] = scaling_equation(phi, bar_phi)
% function [ scal, D_phi ] = scaling_equation(phi, bar_phi)
%
% we need something to block the shifting of the periodic solution, here we
% can have some testing about it
%
% INPUT
% u         Fourier_2D
% bar_u     Fourier_2D, approximation of the solution
%
% OUTPUT
% scal      result of the scalar equation
% D_u       derivative with respect to u (all other derivatives being 0)


if ~isa(phi, 'Fourier_2D')
    if any( size(phi) ~= size(bar_phi))
        error('Size of solution and its approximation incompatible')
    end
    nodes_x = (size(phi,1)-1)/2;
    nodes_y = (size(phi,2)-1)/2;
else
    if any( size(phi.vector) ~= size(bar_phi))
        error('Size of solution and its approximation incompatible')
    end
    nodes_x = phi.nodes_x;
    nodes_y = phi.nodes_y;
end

K1 = (nodes_x:nodes_x);
K1 = repmat(K1,1,2*nodes_y+1);

K2 = (nodes_y:nodes_y);
K2 = repmat(K2,2*nodes_x+1,1);


scal = sum(sum( bar_phi(end:-1:1,end:-1:1) .* phi )) - 1;

if nargout>1
    D_phi = bar_phi(end:-1:1,end:-1:1);
end
