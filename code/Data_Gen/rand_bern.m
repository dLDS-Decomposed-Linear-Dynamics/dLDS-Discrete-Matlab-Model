function x = rand_bern(N, p)

if numel(N) == 1
    x_tmp = rand(N,1);
else
    x_tmp = rand(N,1);
end

x = zeros(size(x_tmp));
x(x_tmp >= p) = 1;
x(x_tmp < p) = -1;

end
