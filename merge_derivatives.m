function huge_der = merge_derivatives(beta, xi, bar_xi, bool_big )
% function huge_der = merge_derivatives(beta, xi, bar_xi, bool_big)
% 
% this function computes the derivative of the compelte problem and merges
% it in a matrix. 
% the third argument is a dimension argument indicating the dimension of
% the return matrix. 
%
% INPUT
% beta          scalar, fixed parameter
% xi            Xi_vector - here no other choice 
% bar_xi        Xi_vector - here no other choice
% bool_big      bool - DEFAULT 0 - if 1, the return matrix is three times
%               as big as needed
%
% OUTPUT
% huge_der      matrix storing the full derivative of the problem
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

[ ~, D_u ] = shift_equation(xi.u, bar_xi.u);
[ ~, D_phi ] = scaling_equation(xi.phi, bar_xi.phi);
[ DF,DF_epsilon] = DFHN_beta(beta, xi);
[ Depsilon, Du, Dkappa, Dphi ] = Deigenvalue_problem_beta(beta, xi);


huge_der = sparse([0 0 D_u 0*D_u
    0 0 0*D_phi D_phi
    DF_epsilon 0*DF_epsilon DF 0*DF
    Depsilon Dkappa Du Dphi]);

