function rhs = full_righthandside(beta, xi, bar_xi, bool_big)
% function rhs = full_righthandside(beta, xi, bar_xi, bool_big)
% 
% this function computes the complete problem and merges
% it in a vector,
% the third argument is a dimension argument indicating the dimension of
% the return vector. 
%
% INPUT
% beta          scalar, fixed parameter
% xi            Xi_vector - here no other choice 
% bar_xi        Xi_vector - here no other choice
% bool_big      bool - DEFAULT 0 - if 1, the return vector is three times
%               as big as needed
%
% OUTPUT
% huge_der      vector storing the full result of the problem
%               shift_equation
%               scaling_equation
%               FHN_beta
%               eigenvalue_problem
%
%               dimension is either length(xi) or 3*length(xi), to allow
%               the convolutions to be of the full lenght

if nargin > 3 && bool_big == 1
    xi = pad(xi, 3*[xi.nodes_x, xi.nodes_y]);
    bar_xi = pad(bar_xi, 3*[bar_xi.nodes_x, bar_xi.nodes_y]);
end

if ~compatible(xi, bar_xi)
    error('inputs need to be compatible')
end

scal_1 = shift_equation(xi.u, bar_xi.u);
scal_2  = scaling_equation(xi.phi, bar_xi.phi);
u_res = FHN_beta(beta, xi);
phi_res = eigenvalue_problem(beta, xi);

rhs = [ scal_1;
    scal_2;
    u_res.vector(:);
    phi_res.vector(:)];
