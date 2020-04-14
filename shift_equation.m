function [ scal, D_u ] = shift_equation(u, bar_u)
% function [ scal, D_u ] = shift_equation(u, bar_u)
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


if ~isa(u, 'Fourier_2D')
    if any( size(u) ~= size(bar_u))
        error('Size of solution and its approximation incompatible')
    end
    nodes_x = (size(u,1)-1)/2;
    nodes_y = (size(u,2)-1)/2;
else
    if any( size(u.vector) ~= size(bar_u.vector))
        error('Size of solution and its approximation incompatible')
    end
    nodes_x = u.nodes_x;
    nodes_y = u.nodes_y;
    u = u.vector;
    bar_u = bar_u.vector;
end

K1 = (nodes_x:nodes_x);
K1 = repmat(K1,1,2*nodes_y+1);

K2 = (nodes_y:nodes_y);
K2 = repmat(K2,2*nodes_x+1,1);


scal = sum(sum(1i* (K1 + K2).* conj(bar_u) .* u ));

if nargout>1
    D_u = 1i* (K1 + K2).* conj(bar_u);
    D_u = D_u(:).';
end
