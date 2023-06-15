%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% External Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_search_settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mult_range = linspace(0.01,1,10);
den_range = linspace(0.01,1,10);
lambda_range = linspace(0.0005,0.05,5);

PSNR_KEEP = zeros(numel(mult_range), numel(den_range), numel(lambda_range), T_s);

clear param_vals
for ll1 = 1:numel(mult_range)
for ll2 = 1:numel(den_range)
    TMP_KEEP = zeros(numel(lambda_range), T_s);
    parfor ll3 = 1:numel(lambda_range)

        param_vals = set_param_vals(TOL, mult_range(ll1), den_range(ll2), lambda_range(ll3),0.4);
        [coef_mdrw, recon_mdrw, rMSE_mdrw, PSNR_mdrw] = ...
            RWL1_DF2_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);

        TMP_KEEP(ll3, :) = rMSE_mdrw;
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

