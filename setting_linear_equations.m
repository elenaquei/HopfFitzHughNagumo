function lin_eqs = setting_linear_equations(x_approx, tangent)
% function lin_eqs = setting_linear_equations(x_approx, tangent)
%
% if the input is small_Xi_vector,
% returns two linear equations:
%     the continuation equation 
%     the phase condition
%
% if the input is Xi_vector,
% returns three linear equations:
%     the continuation equation 
%     the phase condition 
%     the scaling condition for the eigenvector
%
% DEFAULT tangent: parameter continuation in beta

if isa(x_approx,'small_Xi_vector')
    number_equations = 2;
elseif isa(x_approx,'Xi_vector')
    number_equations = 3;
else
    error('Wrong input')
end

if nargin > 1
    if ~iscompatible(x_approx, tangent)
        error('Inputs incompatible')
    end
else 
    tangent = x_approx * 0;
    tangent.beta = 1;
end

%% continuation equation
cont_eqs = linear_eq(tangent, -lin_prod(tangent, x_approx));

%% phase condition
phase_cond = linear_eq(der_u(x_approx),0);

if number_equations == 2
    lin_eqs = linear_equations(cont_eqs,phase_cond);
    return
end

%% scaling condition for the eigenvector
phi_approx = Xi_vector(0,0,0,0*x_approx.u,x_approx.phi);
scaling = linear_eq(phi_approx,-1);

lin_eqs = linear_equations(cont_eqs,phase_cond, scaling);