function DF = derivative_F1(xi, cont, phase_cond, varargin)
% function derivative_F1(xi, cont, phase_cond)
% 
% computed derivative of F1 (x) = (cs, phis, G1)
% needs x, cs, phis
%
% INPUT
% xi                         Xi_vector
% Either cont, phase_cond    linear_eq
% or cont                    linear_equations
%
% OUTPUT
% DF                         derivative_exact


if nargin == 3
    all_lin_eqs = linear_equations(cont, phase_cond);
elseif nargin == 2
    all_lin_eqs = cont;
else
    error('Number of inputs kind of odd... either 2 or 3')
end

% derivatives of linear scalar equations
[mat2by2, fourier_2_C] = lin_eqs_to_der(linear_eqs);

%other derivatives
C_to_fourier = cell(1,2);



%?DG_1 = \begin{pmatrix}
%-\epsilon^{-2}(u - \frac{1}{3}u\conv u \conv u + \gamma^{-1}u + \beta\gamma^{-1}) & \gamma^{-1}\epsilon^{-1} & 0 & \Delta + \epsilon^{-1} (Id - u\conv u  + \gamma^{-1}Id) & 0
%\end{pmatrix}
epsilon = xi.epsilon;
beta = xi.beta;
u = xi.u;
gamma = 0.5;
Constant_Four = constantSequence(u);
% highesst power used is three, this way u had its maximum size
u = prod(prod(u,Constant_Four, varargin{:}),Constant_Four,varargin{:});

% now everything has the right size, no need for varargins anymore
Constant_Four = constantSequence(u);

C_to_fourier{1} = -epsilon^-2*(u - 1/3 * prod(prod(u,u,'same'),u,'same') + 1/gamma*  u + beta/gamma);
C_to_fourier{2} = 1/(gamma * epsilon) * Constant_Four;

vec_fourier_2_fourier = 1/epsilon * (Constant_Four - prod(u,u,'same')  + gamma * Constant_Four);

sign_Delta = +1;

DF = derivative_exact(mat2by2, fourier_2_C, C_to_fourier, sign_Delta, ...
    vec_fourier_2_fourier);

