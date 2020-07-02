function rhs = small_righthandside(xi, bar_xi)
% function rhs = small_righthandside(xi, bar_xi)
% 
% this function computes the complete problem and merges
% it in a vector,
% the third argument is a dimension argument indicating the dimension of
% the return vector. 
%
% INPUT
% beta          scalar, fixed parameter
% xi            (small_)Xi_vector
% bar_xi        (small_)Xi_vector - same as xi
%
% OUTPUT
% rhs           small_Xi_vector solving the problem:
%               continuation equation, 
%               shift_equation,
%               FHN_beta

% if nargin > 3 && bool_big == 1
%     xi = pad(xi, 3*[xi.nodes_x, xi.nodes_y]);
%     bar_xi = pad(bar_xi, 3*[bar_xi.nodes_x, bar_xi.nodes_y]);
% end

if ~compatible(xi, bar_xi)
    error('inputs need to be compatible')
end

lin_eqs = setting_linear_equations(bar_xi); % continuation equation along the beta axes
scal = eval(lin_eqs, xi);
u_res = FHN(xi);

rhs = small_Xi_vector(scal,u_res);
