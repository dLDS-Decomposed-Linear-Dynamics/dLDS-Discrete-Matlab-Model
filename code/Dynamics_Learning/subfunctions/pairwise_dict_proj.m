function [ proj ] = pairwise_dict_proj( F1, F2, frm, T, theta, scalar, scalar_method )
%PAIRWISE_DICT_PROJ Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    scalar = 1;
end

if nargin < 7
    scalar_method = 'scalar';
end

Fsum0 = theta(1)*F1 + theta(2)*F2;

switch scalar_method
    case 'scalar'
        Fsum = scalar*Fsum0;
        
    case 'sqrt'
        if scalar-round(scalar) ~= 0
            error('if choosing sqrt as the scalar_method, scalar should be an integer')
        end
        Fsum = Fsum0;
        for ii = 1:scalar
            Fsum = sqrtm(Fsum);
        end
        Fsum = real(Fsum);
        Ftest = Fsum; 
        for ii = 1:scalar
            Ftest = Ftest*Ftest;
        end
        if norm(Ftest-Fsum0,'fro')/norm(Fsum0,'fro') >= 1e-6
            warning('Ftest is significantly different from Fsum0')
        end
        
    case 'expm'
        [V,D] = eig(Fsum0);
        Fsum = real(V*diag(diag(D).^scalar)/V);
        
    case 'eig'
        Fsum = scalar*Fsum0/max(abs(eig(Fsum0)));
end

proj = plot_dyn_traj(Fsum, frm, T);