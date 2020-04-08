function u = vec2Fourier(vec, x_nodes, y_nodes)
% function u = vec2Fourier(vec, x_nodes, y_nodes)
% 
% This function transforms a vector into a Fourier_2D element, by reshaping
% the vector in the right format.
%
% INPUT
% vec       the initial vector
% x_nodes   integer, nodes in the x direction (rows)
% y_nodes   integer, nodes in the y direction (columns)
%   only one needs to be given, the other can be found
%
% OUTPUT
% u         an element of Fourier_2D 

if length(size(vec))~= 2 || min(size(vec))~=1
    error('First element given is not a vector')
end

if nargin == 3 && ~isempty(x_nodes) && ~isempty(y_nodes)
    if (2*x_nodes+1)*(2*y_nodes+1) ~= length(vec)
        error('Given number of nodes not compatible')
    end
elseif nargin<3 || isempty(y_nodes)
    y_nodes = (length(vec)/(2*x_nodes+1) - 1)/2;
elseif isempty(x_nodes)
    x_nodes = (length(vec)/(2*y_nodes+1) - 1)/2;
end

u_mat = reshape(vec, 2*x_nodes+1, 2*y_nodes+1);
u = Fourier_2D(u_mat);

end