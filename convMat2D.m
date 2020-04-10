function der = convMat2D(u,string_dim)
% function der = convMat2D(u,string_dim)
%
% returns a matrix such that 
%    der * v(:) == u conv v
%
% where u and v are either matrices or Fourier_2D elements of same
% dimension
% The second argument is an optional argument
% chosen between '' and 'SAME' (Default: '')
% If the SAME option is chosen, then 
%       der * v(:) == conv( u, v, 'same')

if ~isa(u, 'Fourier_2D')
    u = Fourier_2D(u);
end

if nargin<2 || strcmp( string_dim , 'same')
    u_padded = pad(u,[u.nodes_x*2,u.nodes_y*2]);
else
    u_padded = u;
end
nodes_x_loc = u_padded.nodes_x;
nodes_y_loc = u_padded.nodes_y;

nodes_x = u.nodes_x;
nodes_y = u.nodes_y;
u = u.vector;

tot_dim = (2*nodes_x+1)*(2*nodes_y+1);

der = sparse(tot_dim, tot_dim); 
% being sparse will save a little tiny bit of memory here and there

for i = -nodes_x:nodes_x
    
    length_diag = 2*nodes_x + 1 - abs(i);
    
    u_loc = u_padded.vector(i+nodes_x_loc+1,:);
    
    for j = 1:length_diag
        
        index_row = (max(0,i)+j-1)*(2*nodes_y+1) +(1:(2*nodes_y+1));
        
        index_col = (max(0,-i)+j-1)*(2*nodes_y+1) +(1:(2*nodes_y+1));
        
        der(index_row, index_col) = toeplitz(u_loc(nodes_y_loc+1:end),...
            u_loc(nodes_y_loc+1:-1:1)); 
        % if u was exactly symmetric, we wouldn't need the second term in
        % the toeplitz matrix
        
    end
end

end
