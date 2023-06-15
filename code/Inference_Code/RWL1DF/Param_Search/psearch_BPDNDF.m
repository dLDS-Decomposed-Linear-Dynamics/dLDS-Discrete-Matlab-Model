%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% External Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_search_settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dyn_range = linspace(0.01,0.15,30);
lambda_range = linspace(0.0005,0.1,30);

RMSE_KEEP = zeros(numel(dyn_range), numel(lambda_range), T_s);

clear param_vals
for ll1 = 1:numel(dyn_range)
    TMP_KEEP = zeros(numel(lambda_range), T_s);
    parfor ll2 = 1:numel(lambda_range)
    
        param_vals = set_param_vals(TOL, dyn_range(ll1), lambda_range(ll2));
        [coef_dcs, recon_dcs, rMSE_dcs, PSNR_dcs] = ...
            BPDN_DF_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);

        TMP_KEEP(ll2, :) = rMSE_dcs;
        fprintf('%f \n', 100*sum(RMSE_KEEP(:)~=0)/numel(RMSE_KEEP))
    end
    RMSE_KEEP(ll1,:,:) = TMP_KEEP;
    fprintf('%f \n', 100*sum(RMSE_KEEP(:)~=0)/numel(RMSE_KEEP))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_params = find_opt_params(RMSE_KEEP, dyn_range, lambda_range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

