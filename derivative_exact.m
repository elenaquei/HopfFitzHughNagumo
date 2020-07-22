classdef derivative_exact
    properties
        size_scalar % default 3
        size_Fourier % default 2
        C3_to_C3 % a 3x3 complex matrix
        Fourier2_to_C3 % 3x2 cell of Fourier_2D
        C3_to_Fourier2 % 2x3 cell of Fourier_2D
        sign_Delta % 2x2 integer matrix with values in {-1,0,1}
        vector_Fourier2_to_Fourier2 % 2x2 cell of Fourier_2D
    end
    methods
        function z = derivative_exact(mat3by3, fourier_2_C, C_to_fourier, sign_Delta, vec_fourier_2_fourier)
            % function z = derivative_exact(mat3by3, fourier_2_C, C_to_fourier, sign_Delta, vec_fourier_2_fourier)
            %
            % constructor testing all dimensions
            if size(mat3by3,1) ~=size(mat3by3,2)
                error('Square matrix must be square!')
            else
                z.size_scalar = size(mat3by3,1);
                z.C3_to_C3 = mat3by3;
            end
            if ~iscell(fourier_2_C) || ~iscell(C_to_fourier) || ~iscell(vec_fourier_2_fourier)
                error('All other elements must be cells')
            end
            if size(vec_fourier_2_fourier,1)~= size(vec_fourier_2_fourier,2)
                error('Square cell must be square!')
            else
                z.size_Fourier = size(vec_fourier_2_fourier,1);
                z.vector_Fourier2_to_Fourier2 = vec_fourier_2_fourier;
            end
            if size(fourier_2_C,1) ~= z.size_scalar|| size(fourier_2_C,1) ~= z.size_Fourier
                error('Size of operator l1^n to C^m incompatible')
            else
                z.Fourier2_to_C3 = fourier_2_C;
            end
            
            if size(C_to_fourier,1) ~= z.size_Fourier || size(C_to_fourier,1) ~=z.size_scalar
                error('Size of operator l1^n to C^m incompatible')
            else
                z.C3_to_Fourier2 = C_to_fourier;
            end
            if mod(sign_Delta,1)~= 0
                error('Power of Delta must be an integer')
            elseif abs(sign_Delta)>1
                error('Power of Delta expected to be -1, 0 or 1')
            else
                z.sign_Delta = sign_Delta;
            end
        end
        
        function y = derivative_exact_to_approx(z, nodes)
            if nargin < 2
                cell_vectors = z.vector_Fourier2_to_Fourier2;
                max_nodes_x = 0;
                max_nodes_y = 0;
                for i = 1:z.size_Fourier
                    for j = 1:z.size_Fourier
                        max_nodes_x = max(cell_vectors{i,j}.nodes_x, max_nodes_x);
                        max_nodes_y = max(cell_vectors{i,j}.nodes_y, max_nodes_y);
                    end
                end
                nodes = [max_nodes_x, max_nodes_y];
            end
            for i = 1:z.size_Fourier
                for j = 1:z.size_Fourier
                    cell_vectors{i,j} = setNodes(cell_vectors{i,j}, nodes);
                end
            end
            size_mat = z.size_Fourier * prod(2 * nodes + 1);
            finite_matrix = zeros(size_mat,size_mat);
            for i = 1:z.size_Fourier
                index_i = (i-1) * prod(2 * nodes + 1) +1 : i * prod(2 * nodes + 1);
                for j = 1:z.size_Fourier
                    index_j = (j-1) * prod(2 * nodes + 1) +1 : j * prod(2 * nodes + 1);
                    finite_matrix(index_i,index_j) = ...
                        derivative_operator_Delta(cell_vectors{i,j}) + ...
                        convMat2D(cell_vectors{i,j});
                end
            end
            y = derivative_approx(z.C3_to_C3, z.Fourier2_to_C3,...
                z.C3_to_Fourier2, z.sign_Delta, finite_matrix);
        end
    end
end