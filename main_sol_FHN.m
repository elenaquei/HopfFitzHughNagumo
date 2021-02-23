% main - find initial solution to FHN with fixed epsilon
global nu
nu = 1 + 10^-6;

beta = 1.4;
% a good place as any to start from (quite at the right end of the graph)

epsilon = 0.03;
% at this point there should be fixed time spiral waves (even if unstable)

nodes_x = 4;
nodes_y = 3;
u = zeros(2*nodes_x+1, 2*nodes_y+1);


x = -1/2:0.05:1/2;
y = -1/2:0.05:1/2;
[X,Y] = meshgrid(x,y);
z = -exp(-12*abs(X.^2+Y.^2)).*cos(20*pi*(X.^2+Y.^2));
%mesh(X,Y,z);

z_fft = fftshift(fftshift(fft2(ifftshift(ifftshift(z,1),2)),1),2)/numel(z);
% figure
% mesh(X,Y,log(abs(z_fft))) % is a very pretty plot!

vec = ifftshift(ifftshift(z_fft,1),2);
z_fft_2 = ifftshift(ifftshift(ifft2(vec,'symmetric'),1),2);


z_F2D = Fourier_2D(z_fft);

z_F2D = symmetrise(z_F2D);

small_xi = small_Xi_vector(epsilon, beta, z_F2D);
plot(z_F2D)
%figure
colormap winter
zdot = FHN(epsilon, beta, z_F2D);
zdot2 = FHN(small_xi);

%plot(zdot)
DF = DFHN_beta(epsilon, beta, z_F2D);

small_xi = small_Xi_vector(epsilon, beta, z_F2D);
xi = Xi_vector(epsilon, beta, rand, z_F2D,  rand(size(z_F2D)).*z_F2D);

xi = symmetrise(xi);
 
FHN(epsilon, beta,z_F2D);
DFHN_beta(epsilon, beta, z_F2D);

eigenvalue_problem(xi);
Deigenvalue_problem_beta(xi);

merge_derivatives(xi, xi);

%better_xi = Newton_beta(beta, xi, xi);
lin_eqs = setting_linear_equations(small_xi); 
eval_lin = eval(lin_eqs, small_xi);
if norm(eval_lin) > 10^-10
    error('linear equations trivially wrong')
end


% this works!
u = z_fft;
u_test = conv2(u,u,'same')+ 3*u;

d_num = convMat2D(u_test,'same');
func_loc = @(x) conv2(x,u_test,'same');
flat = @(x) x(:);
non_flat = @(x) reshape(x, size(u));
flat_conv = @(x) flat(func_loc(non_flat(x)));
d_debug = numerical_der(@(x) flat_conv(x), flat(u));
error_comp =  max(max(abs(d_debug - d_num)));%norm(d_debug- d_num);
if error_comp > 10^-5
    error('something wrong in the derivatives')
end

Delta = operator_Delta(u);

u = z_fft;
u_test = conv2(u,u,'same')+ 3*u;

d_num = diag(Delta(:)) + convMat2D(u_test,'same');
func_loc = @(x) Delta.* x + conv2(x,u_test,'same');
flat = @(x) x(:);
non_flat = @(x) reshape(x, size(u));
flat_conv = @(x) flat(func_loc(non_flat(x)));
d_debug = numerical_der(@(x) flat_conv(x), flat(u));
error_comp = max(max(abs(d_debug - d_num)));%norm(d_debug- d_num);
if error_comp > 10^-5
    error('something wrong in convMat2D')
end



Fx = small_righthandside(small_xi, small_xi);

lin_eqs = setting_linear_equations(small_xi); 
DF = derivative_F1(small_xi, lin_eqs, 'same');

vec2xi = @(x) vec2small_Xi_vec(x,small_xi);
f_func = @(x) Xi2vec( small_righthandside (vec2xi(x),small_xi));


small_xi_vec = Xi2vec(small_xi);

f_func(small_xi_vec);

D_num = numerical_der(f_func, Xi2vec(small_xi));
D_debug = derivative_exact_to_matrix(DF);

error_comp = norm(D_debug- D_num)/norm(D_num);
if error_comp > 10^-2
    error('something wrong in the derivatives')
end
% honestly: there might still be something wrong, but it really seem more
% of a numerical rounding that anything else at the moment

if rcond(D_debug)<10^-10 || rcond(D_debug)>10^14 
    fprintf('Condition number %e\n', rcond(D_debug))
    error('Condition number too bad')
else
    Newton_HN(small_xi)
end








