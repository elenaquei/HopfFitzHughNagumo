classdef Fourier_vector
    properties
        vector       % vector with 2n+1 elements
        nodes        % number of nodes
        % w_ifft % ifft of the main vector - might it be necessary?
    end
    methods
        function w = Fourier_vector(sequence)
           % function w = Fourier_vector(sequence)
           %
           % initialise an element of Fourier_vector with the Fourier
           % sequence being given
           
           % check that the input is a vector
           if (numel(sequence)-length(sequence))~= 0
               error('Vector of Fourier_vector needs to have one non-simpleton dimension')
           end
           
           % check that it is a vertical or horizontal vector
           if length(size(sequence))~=2
               sequence = squeeze(sequence);
           end
           
           % check it has a well defined number of nodes
           if mod(length(sequence)-1,2)~=0
               error('For validation purpose, we do not accept Fourier vectors with illdefined nodes')
           end
           
           % store
           w.vector = sequence;
           w.nodes = (length(sequence)-1)/2;
           % w_ifft = ifft(iffthift(sequence));
        end
        
        function z = setNodes(w, new_nodes)
            % function z = pad(w, new_nodes) 
            %
            % SETNODES the Fourier vectors up to the new number of nodes. Can
            % truncate if new_nodes<w.nodes or pad in the opposite
            % scenario.
            %
            % INPUT
            % w          Fourier_vector
            % new_nodes  positive integer
            % OUTPUT
            % z1    Fourier_vector with nodes = new_nodes
            
            if new_nodes<w.nodes
                % the number of nodes requested is smaller that the
                % starting number of nodes - truncation
                difference = w.nodes - new_nodes;
                z_vec = w.vector(difference+1:end-difference-1);
                z = Fourier_vector(z_vec);
            elseif  new_nodes==w.nodes
                % the number of nodes requested is the same as the vector
                % has already - equality
                z = w;
            else
                % the number of nodes requested is greater than the number
                % of nodes the vector already has - padding
                difference = new_nodes - w.nodes; 
                
                % checking leading dimension
                if size(w.vector,1)==1
                    zero_vec = zeros(1,difference);
                    z_vec = cat(2,cat(2,zero_vec,w.vector),zero_vec);
                else
                    zero_vec = zeros(difference,1);
                    z_vec = cat(1,cat(1,zero_vec,w.vector),zero_vec);
                end
                z = Fourier_vector(z_vec);
            end
            
        end
        
        function z = pad(w,new_nodes)
            % function z = pad(w,new_nodes)
            % 
            % PAD adds zeros until reaching new_nodes
            % If new_nodes is smaller than w.nodes, w is returned unchanged
            if new_nodes>w.nodes
                z = setNodes(w,new_nodes);
            else
                z = w;
            end
        end
        
        function z = truncate(w,new_nodes)
            % function z = pad(w,new_nodes)
            % 
            % TRUNCATE removes elements until reaching new_nodes
            % If new_nodes is bigger than w.nodes, w is returned unchanged
            if new_nodes<w.nodes
                z = setNodes(w,new_nodes);
            else
                z = w;
            end
        end
        
        function [z1,z2] = padVec(w1,w2)
            % function [z1,z2] = padVec(w1,w2)
            %
            % two Fourier vectors are given, it returns two Fourier
            % vectors of the same length by padding the shorter one with
            % the right number of zeros.
            %
            % INPUT 
            % w1,w2     Fourier_vectors
            % OUTPUT
            % z1,z2     Fourier_vectors corresponding to w1 and w2 such
            %   that z1.nodes = z2.nodes = max(w1.nodes, w2.nodes)
            %   either w1 or w2 has been padded accordingly
            
            nodes_z = max(w1.nodes, w2.nodes);
            z1 = pad(w1,nodes_z);
            z2 = pad(w2,nodes_z);
        end
        
        function z = prod(u,v,varargin)
            % function z = prod(u,v,varargin)
            % 
            % PROD handles multiplications with scalars and other Fourier
            % vectors through standard convolution. Extra arguments are
            % passed to the convolution function
            
            if isa(v,'Fourier_vector')
                z = Fourier_vector(conv(u.vector,v.vector,varargin{:}));
            elseif isa(v,'float')
                z = Fourier_vector(u.vector*v);
            else
                error('What are you multilying with??')
            end
        end
        
        function bool = eq(u,v)
            % function bool = eq(u,v)
            % 
            % testing if two Fourier vectors are equal
            bool = false;
            if ~isa(v,'Fourier_vector')
                return
            elseif u.nodes~=v.nodes
                return
            elseif any(u.vector~=v.vector)
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
            % u,v   Fourier_vector
            % tol   positive float, DEFAULT 10^-6
            % OUTPUT
            % bool
            
            if nargin==2 || isempty(tol)
                tol = 10^-6;
            end
            
            bool = false;
            if ~isa(v,'Fourier_vector')
                return
            elseif u.nodes~=v.nodes
                return
            elseif any(abs(u.vector-v.vector)>tol)
                return
            end
            bool = true;
            
        end
        
        function z = plus(u,v)
            % function z = plus(u,v)
            % 
            % PLUS handles sums with scalars and other Fourier
            % vectors.
            
            if isa(v,'float')
                v = Fourier_vector(v);
            end
            
            if isa(v,'Fourier_vector')
                if u.nodes == v.nodes
                    z = Fourier_vector(u.vector+v.vector);
                else
                    [u_long,v_long] = padVec(u,v);
                    z = Fourier_vector(u_long.vector+v_long.vector);
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
                v = Fourier_vector(v);
            end
            
            if isa(v,'Fourier_vector')
                if u.nodes == v.nodes
                    z = Fourier_vector(u.vector-v.vector);
                else
                    [u_long,v_long] = padVec(u,v);
                    z = Fourier_vector(u_long.vector-v_long.vector);
                end
            else
                error('What are you summing with??')
            end
        end
        
        function x = ifft(w)
            % function x = ifft(w)
            % 
            % IFFT of a Fourier vector
            
            x = ifft(ifftshift(w.vector));
        end
        
        function plot(w, varargin)
            % function plot(w, varargin)
            % 
            % plotting functionality for Fourier_vector
            
            if w.nodes<1000
                w = pad(w,1000);
            end
            x = ifft(w);
            x = [x(:);x(1)]*w.nodes;
            plot(1:length(x), x, varargin{:});
            
        end
    end
end