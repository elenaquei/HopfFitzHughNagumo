function sol_approx = Newton_HN(xi, bar_xi)
% function sol_approx = Newton_HN(beta, epsilon, u_start)
%
% INPUT
% xi            small_Xi_vector
% 
% OUTPUT
% epsilon_end   scalar, final approximation
% u_end         Fourier_2D, final approximation 
%
% no shifts allowed, only looking for the solution 

if nargin<2
    bar_xi = xi;
end
% xi will be updated, bar_xi will not

nx = bar_xi.nodes_x;
ny = bar_xi.nodes_y;

lin_eqs = setting_linear_equations(bar_xi); 
iter = 40;
tol = 10^-6;

for i = 1:iter
    Fx = small_righthandside(xi, bar_xi);
    
    if norm(Fx) < tol
        break
    end
    DFx = derivative_F1(xi, lin_eqs, 'same');
    
    % all to matrices and vectors
    F = Xi2vec(Fx);
    vec = Xi2vec(xi);
    DF_mat = derivative_exact_to_matrix(DFx);
    
    % main formula
    vec_new = vec - DF_mat \ F ;
    
    Fourier_vec = Fourier_2D(reshape(vec_new(3:end),2 * nx +1, 2 * ny+1));
    
    xi = small_Xi_vector(vec_new(1), vec_new(2), Fourier_vec);
    
    
    if any(isnan(vec_new))
        if i > 2
            error('Newton diverged at iteration %i', i)
        else
            error('Problems imediately in Newton')
        end
    end
end

Fx = small_righthandside(xi, bar_xi);

if norm(Fx) > tol
    warning('Newton did not converge')
else
    sol_approx = xi;
end

    