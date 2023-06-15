%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% External Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_search_settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_range = linspace(0.05,0.5,20);
eps_range = linspace(5e-5,1.5e-3,30);
RMSE_KEEP = zeros(numel(eps_range), numel(lambda_range), T_s);

TEMP_KEEP = zeros(numel(lambda_range),T_s);
for ll1 = 1:numel(eps_range)
    parfor ll2 = 1:numel(lambda_range)

        param_vals = set_param_vals2(TOL, lambda_range(ll2), eps_range(ll1)*sum(y(:, :, 1).^2));
        [coef_modcs, recon_modcs, rMSE_modcs, PSNR_modcs] = ...
            modCS_largescale(y, MEAS_FUN, DYN_FUN, DWTfunc, param_vals, x0);

   
        TEMP_KEEP(ll2, :) = rMSE_modcs;
        fprintf('%f \n', 100*(ll1)/numel(lambda_range))
    end
    RMSE_KEEP(ll1, :, :) = TEMP_KEEP;
    fprintf('%f \n', 100*sum(RMSE_KEEP(:)~=0)/numel(RMSE_KEEP))
end
fprintf('%f \n', 100*sum(RMSE_KEEP(:)~=0)/numel(RMSE_KEEP))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_params = find_opt_params(RMSE_KEEP, eps_range, lambda_range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


