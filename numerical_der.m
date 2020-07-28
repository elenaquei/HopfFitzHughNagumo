function der = numerical_der(func, vec)
M = length(func(vec));
N = length(vec);
der = zeros(M,N);

h = 10^-5;

for i = 1:M
    epsilon_i = zeros(size(vec));
    epsilon_i(i) = h;
    der(:,i) = (func(vec + epsilon_i) - func(vec) )/ h;
end