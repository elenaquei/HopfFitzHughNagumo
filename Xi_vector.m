classdef Xi_vector
    properties
        epsilon % float - parameter in FitHugh-Nagumo
        kappa % float - imaginary part of eigenvalue
        u % Fourier_2D - first component of solution 
        phi1 % Fourier_2D - first component of eigenvector
    end
    methods
        function xi = Xi_vector(epsilon_sol, kappa_sol, u_sol,  phi1_sol)
            % function xi = Xi_vector(epsilon_sol, kappa_sol, u_sol,  phi1_sol)
            %
            % initialise a Xi_vector, while testing its elements
            
            if ~isa(u_sol, 'Fourier_2D')
                u_sol = Fourier_2D(u_sol);
            end
            if ~isa(epsilon_sol, 'float') || numel(epsilon_sol)~=1
                error('The parameter EPSILON has to be a scalar')
            end
                
            if ~isa(phi1_sol, 'Fourier_2D')
                phi1_sol = Fourier_2D(phi1_sol);
            end
            
            if ~isa(kappa_sol, 'float')|| numel(kappa_sol)~=1
                error('The parameter KAPPA has to be a scalar')
            end
            
            if ~eq_dim(u_sol, phi1_sol) 
                error('All Fourier_2D should have the same dimension')
            end
            
            xi.u = u_sol;
            xi.epsilon = epsilon_sol;
            xi.phi1 = phi1_sol;
            xi.kappa = kappa_sol;
        end
        
        function z = plus(xi1, xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector')
                error('Xi_vectors only add between themselves')
            end
            z.u = xi1.u + xi2.u;
            z.beta = xi1.epsilon + xi2.epsilon;
            z.phi1 = xi1.phi1 + xi2.phi1;
            z.kappa = xi1.kappa + xi2.kappa;
        end
        
        function z = minus(xi1, xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector')
                error('Xi_vectors only substract between themselves')
            end
            z.u = xi1.u - xi2.u;
            z.epsilon = xi1.epsilon - xi2.epsilon;
            z.phi1 = xi1.phi1 - xi2.phi1;
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
            z.epsilon = xi1.epsilon * a;
            z.phi1 = xi1.phi1 * a;
            z.kappa = xi1.kappa * a;
        end
        
        
        function vec = Xi2vec(xi)
            vec = [xi.epsilon;
                xi.kappa;
                Fourier2vec(xi.u);
                Fourier2vec(xi.phi1);];
        end
        
        
        function plot(xi, varargin)
            subplot(2,1,1)
            plot(xi.u, varargin{:})
            title('u')
            subplot(2,1,2)
            plot(xi.phi1, varargin{:})
            title('phi1')
        end
        
    end
end