classdef Xi_vector
    properties
        epsilon % float - parameter in FitHugh-Nagumo
        kappa % float - imaginary part of eigenvalue
        u % Fourier_2D - first component of solution 
        phi % Fourier_2D - first component of eigenvector
    end
    methods
        function xi = Xi_vector(epsilon_sol, kappa_sol, u_sol,  phi_sol)
            % function xi = Xi_vector(epsilon_sol, kappa_sol, u_sol,  phi_sol)
            %
            % initialise a Xi_vector, while testing its elements
            
            if ~isa(u_sol, 'Fourier_2D')
                u_sol = Fourier_2D(u_sol);
            end
            if ~isa(epsilon_sol, 'float') || numel(epsilon_sol)~=1
                error('The parameter EPSILON has to be a scalar')
            end
                
            if ~isa(phi_sol, 'Fourier_2D')
                phi_sol = Fourier_2D(phi_sol);
            end
            
            if ~isa(kappa_sol, 'float')|| numel(kappa_sol)~=1
                error('The parameter KAPPA has to be a scalar')
            end
            
            if ~eq_dim(u_sol, phi_sol) 
                error('All Fourier_2D should have the same dimension')
            end
            
            xi.u = u_sol;
            xi.epsilon = epsilon_sol;
            xi.phi = phi_sol;
            xi.kappa = kappa_sol;
        end
        
        function z = plus(xi1, xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector')
                error('Xi_vectors only add between themselves')
            end
            z.u = xi1.u + xi2.u;
            z.beta = xi1.epsilon + xi2.epsilon;
            z.phi = xi1.phi + xi2.phi;
            z.kappa = xi1.kappa + xi2.kappa;
        end
        
        function z = minus(xi1, xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector')
                error('Xi_vectors only substract between themselves')
            end
            z.u = xi1.u - xi2.u;
            z.epsilon = xi1.epsilon - xi2.epsilon;
            z.phi = xi1.phi - xi2.phi;
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
            z.phi = xi1.phi * a;
            z.kappa = xi1.kappa * a;
        end
        
        function z = pad(xi, new_nodes)
            % function z = pad(xi, new_nodes)
            %
            % pads the two Fourier_2D elements
            z = xi;
            z.u = pad(xi.u, new_nodes);
            z.phi = pad(xi.phi, new_nodes);
        end
        
        function z = setNodes(xi, new_nodes)
            % function z = setNodes(xi, new_nodes)
            %
            % pads the two Fourier_2D elements
            z = xi;
            z.u = setNodes(xi.u, new_nodes);
            z.phi = setNodes(xi.phi, new_nodes);
        end
        
        
        function vec = Xi2vec(xi)
            vec = [xi.epsilon;
                xi.kappa;
                Fourier2vec(xi.u);
                Fourier2vec(xi.phi);];
        end
        
        function xi = symmetrise(xi_conj)
            % function xi = symmetrise(xi_conj)
            %
            % symmetrise the possibly complex input
            
            xi = xi_conj;
            xi.epsilon = real(xi_conj.epsilon);
            xi.kappa = real(xi_conj.kappa);
            xi.u = symmetrise(xi_conj.u);
            xi.phi = symmetrise(xi_conj.phi);
        end
        
        
        function bool = compatible(xi1, xi2)
            % function bool = compatible(xi1, xi2)
            %
            % returns 1 if they have the same dimension
            bool = 0;
            
            if ~isa(xi2, 'Xi_vector') || ~isa(xi1, 'Xi_vector')
                return
            end
            if xi1.u.nodes_x == xi2.u.nodes_x  && ...
                    xi1.u.nodes_y == xi2.u.nodes_y  && ...
                    xi1.phi.nodes_x == xi2.phi.nodes_x  && ...
                    xi1.phi.nodes_y == xi2.phi.nodes_y  
                bool = 1;
            end
        end
        
        function len = length(xi)
            % function len = length(xi)
            %
            % length of the corresponding vector 
            
            len = 2 + 2 * length(xi.u) ;
        end
        
        function plot(xi, varargin)
            subplot(2,1,1)
            plot(xi.u, varargin{:})
            title('u')
            subplot(2,1,2)
            plot(xi.phi, varargin{:})
            title('phi')
        end
        
    end
end