% compute FHN and scalars

xy_fourier = rand(5,5);

epsilon = 0.3;
beta = 1;

xi_vec = small_Xi_vector(epsilon, beta, xy_fourier);
xi_vec_hat = small_Xi_vector(epsilon - 10^-5, beta + 10^-5, rand*xy_fourier);

res = small_righthandside(xi_vec, xi_vec_hat);