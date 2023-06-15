function F = initialize_dynamics(num_f, size_f)

% num_f = 5;
% size_f = 5;

F = cell(1, num_f);

kappa = 0.5.^(0:size_f-1);

for kk = 1:numel(F)
    F{kk} = randn(size_f,size_f);
    [U,S,V] = svd(F{kk});
    F{kk} = U*(S.*diag(kappa))*V.';
%     F{kk} = F{kk}/max(abs(eig(F{kk})));
%    F{kk} = exp(-pi * kappa * F{kk}/max(abs(eig(F{kk})))); %changed
%    figure();
%    scatter(1:5, F{kk});
end

end
