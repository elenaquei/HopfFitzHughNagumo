function xi = vec2small_Xi_vec(x,y)
% function vec2small_Xi_vec(x,y)
% 
% INPUT
% x         vector
% y         small_Xi_vector
% OUTPUT
% xi        small_Xi_vector, reshaped x 

if length(x) ~= length(y)
    error('Size is incompatible')
end

xi = small_Xi_vector(x(1),x(2), reshape(x(3:end),y.nodes_x*2+1, y.nodes_y*2+1));