function [y, G] = take_gaussian_meas(x, M, noise_var)

N = size(x, 1);
T = size(x, 2);
G = randn(M, N, T)/sqrt(M);
y = zeros(M, T);

for kk = 1:T
    y(:, kk) = G(:, :, kk)*x(:, kk) + sqrt(noise_var)*randn(M, 1);
end

end