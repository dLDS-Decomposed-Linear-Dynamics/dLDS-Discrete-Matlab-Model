function output = DCS_meas_func(z, m , A, P, NORMS_VEC)

% output = DCS_meas_func(z, m , A, P)
%
% Code that creates a measurement function that is compatable with DCS-AMP
% 
% Code by Adam Charles, 
% Department of Electrical and Computer Engineering,
% Georgia Institute of Technology
% 
% Last updated August 21, 2012. 
% 


NORM_INV = 1./NORMS_VEC(:);
if m == 1
    output = A.Phi(P.invert(NORM_INV.*z(:)));
elseif m == 2
    output = NORM_INV.*P.apply(A.Phit(z));
else
    error('WTF, mate?')
end
