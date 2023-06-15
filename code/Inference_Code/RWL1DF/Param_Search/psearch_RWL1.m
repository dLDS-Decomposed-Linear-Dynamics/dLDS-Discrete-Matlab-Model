%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% External Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_search_settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mult_range = linspace(0.01,1,10);
den_range = linspace(0.01,1,10);
lambda_range = linspace(0.0005,0.005,5);

RMSE_KEEP = zeros(numel(mult_range), numel(den_range), numel(lambda_range), T_s);

clear param_vals
for ll1 = 1:numel(mult_range)
    rwl1_mult = mult_range(ll1);
for ll2 = 1:numel(den_range)
    rwl1_reg = den_range(ll2);

    TMP_KEEP = zeros(numel(lambda_range), T_s);

    parfor ll3 = 1:numel(lambda_range)
        [coef_rwcs, recon_rwcs, rMSE_rwcs, PSNR_rwcs] = ...
            RWL1_largescale(y, MEAS_FUN, [lambda_range(ll3), rwl1_reg, rwl1_mult], ...
            TOL, DWTfunc, x0);

        TMP_KEEP(ll3, :) = rMSE_rwcs;
        fprintf('%f \n', 100*sum(RMSE_KEEP(:)~=0)/numel(RMSE_KEEP))
    end
    RMSE_KEEP(ll1,ll2,:,:) = TMP_KEEP;
    fprintf('%f \n', 100*sum(RMSE_KEEP(:)~=0)/numel(RMSE_KEEP))
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_params = find_opt_params(RMSE_KEEP, mult_range, den_range, lambda_range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
