classdef small_Xi_vector
    properties
        epsilon % float - first parameter in FitzHugh-Nagumo (arbitary numbering)
        beta % float - second parameter in FitzHugh-Nagumo (arbitary numbering)
        u % Fourier_2D - first component of solution 
    end
    methods
        function xi = Xi_vector(epsilon_sol, beta_sol, u_sol)
            % function xi = Xi_vector(epsilon_sol, u_sol)
            %
            % initialise a Xi_vector, while testing its elements
            
            if ~isa(u_sol, 'Fourier_2D')
                u_sol = Fourier_2D(u_sol);
            end
            if ~isa(epsilon_sol, 'float') || numel(epsilon_sol)~=1
                error('The parameter EPSILON has to be a scalar')
            end
            
            if ~isa(beta_sol, 'float') || numel(beta_sol)~=1
                error('The parameter BETA has to be a scalar')
            end
            
            xi.u = u_sol;
            xi.beta = beta_sol;
            xi.epsilon = epsilon_sol;
        end
        
        function z = plus(xi1, xi2)
            if ~isa(xi1,'small_Xi_vector') || ~isa(xi2,'small_Xi_vector')
                error('Xi_vectors only add between themselves')
            end
            z.u = xi1.u + xi2.u;
            z.beta = xi1.beta + xi2.beta;
            z.epsilon = xi1.epsilon + xi2.epsilon;
        end
        
        function z = minus(xi1, xi2)
            if ~isa(xi1,'small_Xi_vector') || ~isa(xi2,'small_Xi_vector')
                error('Xi_vectors only substract between themselves')
            end
            z.u = xi1.u - xi2.u;
            z.beta = xi1.beta - xi2.beta;
            z.epsilon = xi1.epsilon - xi2.epsilon;
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
            if ~isa(xi1,'small_Xi_vector') || ~isa(a,'float') || numel(a)~=1
                error('small_Xi_vectors only multiply with scalars')
            end
            z.u = xi1.u * a;
            z.epsilon = xi1.epsilon * a;
            z.beta = xi1.beta * a;
        end
        
        function z = lin_prod(xi1,xi2)
            if ~isa(xi1,'Xi_vector') || ~isa(xi2,'Xi_vector') 
                error('This is not the product you ar looking for ~ hand wave')
            end
            z = xi1.epsilon * xi2.epsilon + xi1.beta * xi2.beta + ...
                lin_prod(xi1.u,xi2.u);
        end
        
        function z = der(xi1)
            z = der_u(xi1);
        end
        
        function z = der_u(xi1)
            z = xi1;
            z.beta = 0;
            z.epsilon = 0;
            z.u = der(xi1.u);
        end
        
        
        function z = pad(xi, new_nodes)
            % function z = pad(xi, new_nodes)
            %
            % pads the two Fourier_2D elements
            z = xi;
            z.u = pad(xi.u, new_nodes);
        end
        
        function z = setNodes(xi, new_nodes)
            % function z = setNodes(xi, new_nodes)
            %
            % pads the two Fourier_2D elements
            z = xi;
            z.u = setNodes(xi.u, new_nodes);
        end
        
        
        function vec = Xi2vec(xi)
            vec = [xi.epsilon;
                xi.beta;
                Fourier2vec(xi.u)];
        end
        
        function xi = symmetrise(xi_conj)
            % function xi = symmetrise(xi_conj)
            %
            % symmetrise the possibly complex input
            
            xi = xi_conj;
            xi.epsilon = real(xi_conj.epsilon);
            xi.u = symmetrise(xi_conj.u);
        end
        
        
        function bool = compatible(xi1, xi2)
            % function bool = compatible(xi1, xi2)
            %
            % returns 1 if they have the same dimension
            bool = 0;
            
            if ~isa(xi2, 'small_Xi_vector') || ~isa(xi1, 'small_Xi_vector')
                return
            end
            if xi1.u.nodes_x == xi2.u.nodes_x  && ...
                    xi1.u.nodes_y == xi2.u.nodes_y
                bool = 1;
            end
        end
        
        function len = length(xi)
            % function len = length(xi)
            %
            % length of the corresponding vector 
            
            len = 1 + length(xi.u) ;
        end
        
        function plot(xi, varargin)
            plot(xi.u, varargin{:})
            title('u')
        end
        
    end
end