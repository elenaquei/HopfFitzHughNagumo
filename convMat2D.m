function der = convMat2D(u)
% function der = convMat2D(u)
%
% returns a matrix such that 
%    der * v(:) == u conv v
%
% where u and v are either matrices or Fourier_2D elements of same
% dimension

if isa(u, 'Fourier_2D')
    nodes_x = u.nodes_x;
    nodes_y = u.nodes_y;
    u = u.vector;
else
    nodes_x = (size(u,1)-1)/2;
    nodes_y = (size(u,2)-1)/2;
end

tot_dim = (2*nodes_x+1)*(2*nodes_y+1);

der = sparse(tot_dim, tot_dim); 
% being sparse will save a little tiny bit of memory here and there

for i = -nodes_x:nodes_x
    
    length_diag = 2*nodes_x + 1 - abs(i);
    
    for j = 1:length_diag
        
        index_row = (max(0,i)+j-1)*(2*nodes_y+1) +(1:(2*nodes_y+1));
        
        index_col = (max(0,-i)+j-1)*(2*nodes_y+1) +(1:(2*nodes_y+1));
        
        der(index_row, index_col) = toeplitz(u(i+nodes_x+1,nodes_y+1:end),...
            u(i+nodes_x+1,1:nodes_y)); 
        % if u was exactly symmetric, we wouldn't need the second term in
        % the toeplitz matrix
        
    end
end

end
