classdef linear_equations
    properties
        lin_eqs        % cell of 2 or 3 linear_eq 
        n_eqs
    end
    methods
        function w = linear_equations(lin_eq1, lin_eq2, lin_eq3)
           % function w = linear_eqs(coefs, const)
           if ~isa(lin_eq1,'linear_eq') || ~isa(lin_eq2,'linear_eq') 
               error('The linear equations have to be of the prescribed form')
           end
           if nargin>2 && ~isa(lin_eq3, 'linear_eq')
               error('The linear equations have to be of the prescribed formt')
           end
           w = cell(nargin,1);
           w(1) = lin_eq1;
           w(2) = lin_eq2;
           if nargin>2
                w(3) = lin_eq3;
           end
           w.n_eqs = nargin;
        end
        
        function c = eval(lin_equations, vec)
            % function c = eval(lin_eq, vec)
            % evaluate < vec, lin_eq > + c
            c = zeros(lin_equations.n_eqs);
            for i = 1:lin_equations.n_eqs
                c(i) = eval(lin_equations.lin_eqs{i}, vec);
            end
        end
        
        function d = der(lin_equations)
            d = zeros(lin_equations.n_eqs, length(lin_equations.lin_eqs{1}.vector));
            for i = 1:lin_equations.n_eqs
                d(i,:) = der(lin_equations.lin_eqs{i});
            end
        end
    end
end
