% main - find initial solution to FHN with fixed epsilon

beta = 1.4;
% a good place as any to start from (quite at the right end of the graph)

epsilon = 0.03;
% at this point there should be fixed time spiral waves (even if unstable)

nodes_x = 4;
nodes_y = 3;
u = zeros(2*nodes_x+1, 2*nodes_y+1);


x = -1/2:0.01:1/2;
y = -1/2:0.01:1/2;
[X,Y] = meshgrid(x,y);
z = cos(pi*Y)*cos(pi*X);
mesh(X,Y,z);

z_fft = fft2(z);
figure
mesh(X,Y,log(abs(z_fft))) % is a very pretty plot!

