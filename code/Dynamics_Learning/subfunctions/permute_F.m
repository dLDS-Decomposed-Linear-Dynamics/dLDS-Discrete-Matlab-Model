function [ Fp ] = permute_F( F, p_order )
%PERMUTE_F Summary of this function goes here
%   Detailed explanation goes here

N_f = numel(F);
Fp = cell(N_f,1);

for ii = 1:N_f
    as = p_order(3,ii);
    if as ~= 0
        Fp{p_order(1,ii)} = as*F{p_order(2,ii)};
    else
        Fp{p_order(1,ii)} = F{p_order(2,ii)};
    end
end

end

