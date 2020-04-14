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
else
    if any( size(phi.vector) ~= size(bar_phi.vector))
        error('Size of solution and its approximation incompatible')
    end
    phi = phi.vector;
    bar_phi = bar_phi.vector;
end


scal = sum(sum( conj(bar_phi) .* phi )) - 1;

if nargout>1
    D_phi = conj(bar_phi);
    D_phi = D_phi(:).';
end
