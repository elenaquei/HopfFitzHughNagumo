classdef Xi_vector
    properties
        u % Fourier_vector - first component of solution 
        v % Fourier_vector - second component of solution
        beta % float - parameter in FitHugh-Nagumo
        phi1 % Fourier_vector - first component of eigenvector
        phi2 % Fourier_vector - second component of eigenvector
        kappa % float - imaginary part of eigenvalue
    end
    methods
        function xi = Xi_vector(u_sol, v_sol, beta_sol, phi1_sol,...
                phi2_sol, kappa_sol)
            % function xi = Xi_vector(u_sol, v_sol, beta_sol, phi1_sol,...
            %    phi2_sol, kappa_sol)
            %
            % initialise a Xi_vector, while testing its elements
            
            if ~isa(u_sol, 'Fourier_vector')
                u_sol = Fourier_vector(u_sol);
            end
            if ~isa(v_sol, 'Fourier_vector')
                v_sol = Fourier_vector(v_sol);
            end
            if ~isa(beta_sol, 'float')
                error('The parameter BETA has to be real')
            end
                
            if ~isa(phi1_sol, 'Fourier_vector')
                phi1_sol = Fourier_vector(phi1_sol);
            end
            
            if ~isa(phi2_sol, 'Fourier_vector')
                phi2_sol = Fourier_vector(phi2_sol);
            end
            if ~isa(kappa_sol, 'float')
                error('The parameter KAPPA has to be real')
            end
            
            xi.u = u_sol;
            xi.v = v_sol;
            xi.beta = beta_sol;
            xi.phi1 = phi1_sol;
            xi.phi2 = phi2_sol;
            xi.kappa = kappa_sol;
        end
        
        function z = plus(xi1, xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector')
                error('Xi_vectors only add between themselves')
            end
            z.u = xi1.u + xi2.u;
            z.v = xi1.v + xi2.v;
            z.beta = xi1.beta + xi2.beta;
            z.phi1 = xi1.phi1 + xi2.phi1;
            z.phi2 = xi1.phi2 + xi2.phi2;
            z.kappa = xi1.kappa + xi2.kappa;
        end
        
        function z = minus(xi1, xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector')
                error('Xi_vectors only substract between themselves')
            end
            z.u = xi1.u - xi2.u;
            z.v = xi1.v - xi2.v;
            z.beta = xi1.beta - xi2.beta;
            z.phi1 = xi1.phi1 - xi2.phi1;
            z.phi2 = xi1.phi2 - xi2.phi2;
            z.kappa = xi1.kappa - xi2.kappa;
        end
        
        function z = prod(xi1, a)
            % function z = prod(xi1, a)
            %
            % product with scalar
            if isa(xi1,'float')
                temp = xi1;
                xi1 = a;
                a = temp;
            end
            if ~isa(xi1,'Xi_vector') || ~isa(a,'float') || numel(a)~=1
                error('Xi_vectors only multiply with scalars')
            end
            z.u = xi1.u * a;
            z.v = xi1.v * a;
            z.beta = xi1.beta * a;
            z.phi1 = xi1.phi1 * a;
            z.phi2 = xi1.phi2 * a;
            z.kappa = xi1.kappa * a;
        end
        
    end
end