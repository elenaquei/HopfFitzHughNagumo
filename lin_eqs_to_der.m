function [mat, cell_of_Fouriers] = lin_eqs_to_der(linear_eqs)
% function [mat, cell_of_Fouriers] = lin_eqs_to_der(lin_eqs)
%
% helper function to transform linear_equations in the correct format to
% use derivative_exact of derivative_approx
%
% INPUT
% lin_eqs           linear_equations
% OUTPUT
% mat               derivative of all linear equations w.r.t. the scalars
% cell_of_Fouriers  derivative of all linear equations w.r.t. the Fourier
%                   vectors, stored in a cell structure

n_eqs = linear_eqs.n_eqs;
n_scalars = 2;
for i = 1:n_eqs
    if bool_big(linear_eqs.lin_eqs{i})
        n_scalars = 3;
    end
end
n_vectors = n_scalars -1;
mat = zeros(n_eqs,n_scalars);
for i = 1:n_eqs
    mat(i,1) = linear_eqs.lin_eqs{i}.vector.epsilon;
    mat(i,2) = linear_eqs.lin_eqs{i}.vector.beta;
    if bool_big(linear_eqs.lin_eqs{i})
        mat(i,3) = linear_eqs.lin_eqs{i}.vector.kappa;
    end
end

cell_of_Fouriers = cell(n_eqs, n_vectors);
for i = 1:n_eqs
    cell_of_Fouriers{i,1} = linear_eqs.lin_eqs{i}.vector.u;
    if bool_big(linear_eqs.lin_eqs{i})
        cell_of_Fouriers{i,2} = linear_eqs.lin_eqs{i}.vector.phi;
    end
end
