function better_xi = Newton_beta(beta, xi, bar_xi)
% function better_xi = Newton_beta(beta, xi, bar_xi)
%
% purely numerical run looking for a numerical solution to our Hopf
% problem.
%
% INPUT
% beta          scalar, fixed parameter
% xi            Xi_vector
% bar_xi        Xi_vector
%
% OUTPUT
% better_xi     Xi_vector

n_iter = 100;
tol = 10^-6;
step_size = tol*2;
iter = 0;

xi_vec = Xi2vec(xi);

bar_xi = symmetrise (bar_xi);


while norm(step_size) > tol
    
    derivative = merge_derivatives(beta, xi, bar_xi);
    
    rhs = full_righthandside(beta, xi, bar_xi);
    
    % plot(abs(rhs))
    
    step_size = derivative\rhs;
    
    xi_vec = xi_vec - step_size;
    
    xi = vec2Xi_vector(xi_vec, xi);
    
    plot((abs(Xi2vec(xi) - Xi2vec(symmetrise(xi)))))
    
    xi = symmetrise (xi);
    
    
    iter = iter+1;
    
    disp(iter)
    
    if iter > n_iter 
        error('Newton did not converge')
    end
    
end

better_xi = xi;