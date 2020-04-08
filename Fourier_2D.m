classdef Fourier_2D
    properties
        vector         % vector with 2nx+1 * 2ny+1 elements
        nodes_x        % number of nodes in the first dimension
        nodes_y        % number of nodes in the second dimension
        nodes          % array with the number of nodes
        % w_ifft % ifft of the main vector - might it be necessary?
    end
    methods
        function w = Fourier_2D(sequence)
           % function w = Fourier_2D(sequence)
           %
           % initialise an element of Fourier_2D with the Fourier
           % sequence being given
           
           % check that it is a vertical or horizontal vector
           if length(size(sequence))~=2
               sequence = squeeze(sequence);
           end
           
           % check that the input is a matrix
           if (numel(sequence)-size(sequence,1)*size(sequence,2))~= 0
               error('Vector of Fourier_2D needs to have one non-simpleton dimension')
           end
           
           % check it has a well defined number of nodes
           if mod(size(sequence,1)-1,2)~=0 || mod(size(sequence,2)-1,2)~=0
               error('For validation purpose, we do not accept Fourier matrices with illdefined nodes')
           end
           
           % store
           w.vector = sequence;
           w.nodes_x = (size(sequence,1)-1)/2;
           w.nodes_y = (size(sequence,2)-1)/2;
           w.nodes = [w.nodes_x, w.nodes_y];
           % w_ifft = ifft(iffthift(sequence));
        end
        
        
        function z = setNodes_dim(w, new_nodes,dim)
            % function z = pad(w, new_nodes,dim) 
            %
            % SETNODES the Fourier vectors up to the new number of nodes. Can
            % truncate if new_nodes<w.nodes or pad in the opposite
            % scenario.
            %
            % INPUT
            % w          Fourier_2D
            % new_nodes  positive integer
            % dim        either 1 or 2
            % OUTPUT
            % z1    Fourier_2D with nodes = new_nodes
            
            if dim~=1 && dim~=2
                error('Acting dimension must be either 1 or 2')
            end
            if lenght(new_nodes) ==2
                new_nodes = new_nodes(dim);
            end
            
            old_nodes = w.nodes(dim);
            
            if new_nodes < old_nodes
                % the number of nodes requested is smaller that the
                % starting number of nodes - truncation
                difference = old_nodes - new_nodes;
                if dim ==1 
                    z_vec = w.vector(difference+1:end-difference-1,:);
                else
                    z_vec = w.vector(:,difference+1:end-difference-1);
                end
                z = Fourier_2D(z_vec);
            elseif  new_nodes==w.nodes
                % the number of nodes requested is the same as the vector
                % has already - equality
                z = w;
            else
                % the number of nodes requested is greater than the number
                % of nodes the vector already has - padding
                difference = new_nodes - old_nodes; 
                
                % checking leading dimension
                if dim==2
                    zero_vec = zeros(1,difference);
                    z_vec = cat(2,cat(2,zero_vec,w.vector),zero_vec);
                else
                    zero_vec = zeros(difference,1);
                    z_vec = cat(1,cat(1,zero_vec,w.vector),zero_vec);
                end
                z = Fourier_2D(z_vec);
            end
            
        end
        
        function z = setNodes(w, new_nodes)
            % function z = pad(w, new_nodes) 
            %
            % SETNODES the Fourier vectors up to the new number of nodes. Can
            % truncate if new_nodes<w.nodes or pad in the opposite
            % scenario.
            %
            % INPUT
            % w          Fourier_2D
            % new_nodes  array of 2 positive integers
            % OUTPUT
            % z1    Fourier_2D with nodes = new_nodes
            
            z = setNodes_dim(w, new_nodes, 1);
            z = setNodes_dim(z, new_nodes, 2);
            
        end
        
        function z = pad(w,new_nodes)
            % function z = pad(w,new_nodes)
            % 
            % PAD adds zeros until reaching new_nodes
            % If new_nodes is smaller than w.nodes, w is returned unchanged
            reset_nodes = max(w.nodes, new_nodes);
            z = setNodes(w,reset_nodes);
        end
        
        function z = truncate(w,new_nodes)
            % function z = pad(w,new_nodes)
            % 
            % TRUNCATE removes elements until reaching new_nodes
            % If new_nodes is bigger than w.nodes, w is returned unchanged
            reset_nodes = min(w.nodes, new_nodes);
            z = setNodes(w,reset_nodes);
        end
        
        function [z1,z2] = padVec(w1,w2)
            % function [z1,z2] = padVec(w1,w2)
            %
            % two Fourier vectors are given, it returns two Fourier
            % vectors of the same length by padding the shorter one with
            % the right number of zeros.
            %
            % INPUT 
            % w1,w2     Fourier_2Ds
            % OUTPUT
            % z1,z2     Fourier_2Ds corresponding to w1 and w2 such
            %   that z1.nodes = z2.nodes = max(w1.nodes, w2.nodes)
            %   either w1 or w2 has been padded accordingly
            
            nodes_z = max(w1.nodes, w2.nodes);
            z1 = pad(w1,nodes_z);
            z2 = pad(w2,nodes_z);
        end
        
        
        function z = plus(u,v)
            % function z = plus(u,v)
            % 
            % PLUS handles sums with scalars and other Fourier
            % vectors.
            
            if isa(v,'float')
                v = Fourier_2D(v);
            end
            
            if isa(v,'Fourier_2D')
                if all(u.nodes == v.nodes)
                    z = Fourier_2D(u.vector+v.vector);
                else
                    [u_long,v_long] = padVec(u,v);
                    z = Fourier_2D(u_long.vector+v_long.vector);
                end
            else
                error('What are you summing with??')
            end
        end
        
        function z = minus(u,v)
            % function z = minus(u,v)
            % 
            % MINUS handles substractions with scalars and other Fourier
            % vectors.
            
            if isa(v,'float')
                v = Fourier_2D(v);
            end
            
            if isa(v,'Fourier_2D')
                if all(u.nodes == v.nodes)
                    z = Fourier_2D(u.vector-v.vector);
                else
                    [u_long,v_long] = padVec(u,v);
                    z = Fourier_2D(u_long.vector-v_long.vector);
                end
            else
                error('What are you summing with??')
            end
        end
        
        function z = prod(u,v,varargin)
            % function z = prod(u,v,varargin)
            % 
            % PROD handles multiplications with scalars and other Fourier
            % vectors through standard convolution. Extra arguments are
            % passed to the convolution function
            
            if isa(v,'Fourier_2D')
                z = Fourier_2D(conv2(u.vector,v.vector,varargin{:}));
            elseif isa(v,'float')
                z = Fourier_2D(u.vector*v);
            else
                error('What are you multilying with??')
            end
        end
        
        
        function bool = eq(u,v)
            % function bool = eq(u,v)
            % 
            % testing if two Fourier vectors are equal
            bool = false;
            if ~isa(v,'Fourier_2D')
                return
            elseif any(u.nodes~=v.nodes)
                return
            elseif any(u.vector~=v.vector)
                return
            end
            bool = true;
            
        end
        
        function bool = eq_dim(u,v)
            % function bool = eq_dim(u,v)
            % 
            % testing if two Fourier matrices have the same dimension
            bool = false;
            if ~isa(v,'Fourier_2D')
                return
            elseif any(u.nodes~=v.nodes)
                return
            end
            bool = true;
            
        end
        
        function bool = eq_approx(u,v,tol)
            % function bool = eq_approx(u,v,tol)
            % 
            % testing if two Fourier vectors are equal up to a tollerance 
            % if the two vectora have different length they will still be
            % considered different, but the values they store might differ
            % up to the tolerance
            % INPUT
            % u,v   Fourier_2D
            % tol   positive float, DEFAULT 10^-6
            % OUTPUT
            % bool
            
            if nargin==2 || isempty(tol)
                tol = 10^-6;
            end
            
            bool = false;
            if ~isa(v,'Fourier_2D')
                return
            elseif any(u.nodes~=v.nodes)
                return
            elseif any(abs(u.vector-v.vector)>tol)
                return
            end
            bool = true;
            
        end
        

        function vec = Fourier2vec(u)
            % function vec = Fourier2vec(u)
            %
            % returns the vector representation of u
            vec = u.vector(:);
        end
        
        function vec = colon(u)
            % function vec = colon(u)
            %
            % returns the vector representation of u
            vec = u.vector(:);
        end
        
        
        
        function x = ifft(w)
            % function x = ifft(w)
            % 
            % IFFT of a Fourier vector
            vec = ifftshift(ifftshift(w.vector,1),2);
            x = ifft2(vec,'symmetric');
        end
        
        function plot(w, varargin)
            % function plot(w, varargin)
            % 
            % plotting functionality for Fourier_2D
            
            
            w = pad(w,[1000,1000]);
            
            x = ifft(w);
            % x = x*w.nodes;
            
            y_grid = linspace(-1/2,1/2,size(x,1));
            z_grid = linspace(-1/2,1/2,size(x,2));
            grid = mesh(y_grid, z_grid, x);
            
            plot3(grid, varargin{:});
            
        end
    end
end