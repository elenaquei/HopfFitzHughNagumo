% main - find initial solution to FHN with fixed epsilon

beta = 1.4;
% a good place as any to start from (quite at the right end of the graph)

epsilon = 0.03;
% at this point there should be fixed time spiral waves (even if unstable)

nodes_x = 4;
nodes_y = 3;
u = zeros(2*nodes_x+1, 2*nodes_y+1);


x = -1/2:0.1:1/2;
y = -1/2:0.1:1/2;
[X,Y] = meshgrid(x,y);
z = cos(3*pi*Y)*cos(pi*X);
%mesh(X,Y,z);

z_fft = fftshift(fftshift(fft2(ifftshift(ifftshift(z,1),2)),1),2)/numel(z);
% figure
% mesh(X,Y,log(abs(z_fft))) % is a very pretty plot!

vec = ifftshift(ifftshift(z_fft,1),2);
z_fft_2 = ifftshift(ifftshift(ifft2(vec,'symmetric'),1),2);


z_F2D = Fourier_2D(z_fft);
%plot(z_F2D)
%figure
zdot = FHN_beta(beta, epsilon, z_F2D);
%plot(zdot)
DF = DFHN_beta(beta, epsilon, z_F2D);


xi = Xi_vector(epsilon, epsilon, z_F2D,  z_F2D);
 
 
FHN_beta(beta,xi);
DFHN_beta(beta, xi);

eigenvalue_problem(beta, xi);
Deigenvalue_problem_beta(beta, xi);











