function F_out = normalize_dynamics(F_in)


F_out = cell(size(F_in));

for kk = 1:numel(F_out)
    F_out{kk} = F_in{kk}/max(abs(eig(F_in{kk})));
end

end