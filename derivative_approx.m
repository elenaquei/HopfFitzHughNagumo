classdef derivative_approx
    properties
        size_scalar % default 3
        size_Fourier % default 2
        C3_to_C3 % a 3x3 complex matrix
        Fourier2_to_C3 % 3x2 cell of Fourier_2D
        C3_to_Fourier2 % 2x3 cell of Fourier_2D
        sign_Delta % 2x2 integer matrix with values in {-1,0,1}
        mat_Fourier2_to_Fourier2 % square complex matrix
    end
    methods
        function z = derivative_approx(mat3by3, fourier_2_C, C_to_fourier, ...
                sign_Delta, mat_fourier_2_fourier)
            % function z = derivative_exact(mat3by3, fourier_2_C, C_to_fourier, sign_Delta, mat_fourier_2_fourier)
            %
            % constructor testing all dimensions
            if size(mat3by3,1) ~=size(mat3by3,2)
                error('Square matrix must be square!')
            else
                z.size_scalar = size(mat3by3,1);
                z.C3_to_C3 = mat3by3;
            end
            if ~iscell(fourier_2_C) || ~iscell(C_to_fourier)
                error('All other elements must be cells')
            end
            if size(mat_fourier_2_fourier,1)~= size(mat_fourier_2_fourier,2)
                error('Square matrix must be square!')
            else
                
                z.mat_Fourier2_to_Fourier2 = mat_fourier_2_fourier;
            end
            if size(fourier_2_C,1) ~= z.size_scalar
                error('Size of operator l1^n to C^m incompatible')
            else
                z.size_Fourier = size(fourier_2_C,2);
                z.Fourier2_to_C3 = fourier_2_C;
            end
            
            if size(C_to_fourier,1) ~= z.size_Fourier || size(C_to_fourier,2) ~=z.size_scalar
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
        
        function der = derivative_approx_to_mat(z, nodes)
            if nargin<2
                max_nodes_x = 0;
                max_nodes_y = 0;
                for i = 1:z.size_Fourier
                    for j = 1:z.size_scalar
                        nodes_x1 = z.Fourier2_to_C3{j,i}.nodes_x;
                        nodes_x2 = z.C3_to_Fourier2{i,j}.nodes_x;
                        nodes_y1 = z.Fourier2_to_C3{j,i}.nodes_y;
                        nodes_y2 = z.C3_to_Fourier2{i,j}.nodes_y;
                        max_nodes_x = max([max_nodes_x, nodes_x1,nodes_x2]);
                        max_nodes_y = max([max_nodes_y, nodes_y1,nodes_y2]);
                    end
                end
                nodes = [max_nodes_x, max_nodes_y];
            end
            Delta = operator_Delta(nodes(1), nodes(2));
            dim_der = z.size_scalar + size(z.mat_Fourier2_to_Fourier2,1);
            der = zeros(dim_der, dim_der);
            der(1:z.size_scalar,1:z.size_scalar) = z.C3_to_C3;
            for i = 1:z.size_Fourier
                for j = 1:z.size_scalar
                    z.Fourier2_to_C3{j,i} = setNodes(z.Fourier2_to_C3{j,i},nodes);
                    z.C3_to_Fourier2{i,j} = setNodes(z.C3_to_Fourier2{i,j},nodes);
                end
            end
            if size(z.mat_Fourier2_to_Fourier2,1) ~= prod(2*nodes+1)
                error('Not considered yet')
            end
            for i = 1:z.size_Fourier
                index_i = z.size_scalar +( (i-1) * prod(2 * nodes + 1) +1 : i * prod(2 * nodes + 1));
                for j = 1:z.size_scalar
                    der(j, index_i) = Fourier2vec(z.Fourier2_to_C3{j,i}).';
                    der(index_i, j) = Fourier2vec(z.C3_to_Fourier2{i,j});
                end
            end
            der(z.size_scalar+1:end,z.size_scalar+1:end) = diag(Delta(:)).^z.sign_Delta + z.mat_Fourier2_to_Fourier2;
        end
    end
end