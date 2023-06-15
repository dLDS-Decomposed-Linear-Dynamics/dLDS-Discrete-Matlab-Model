function [vargout, app_mat, true_mat, additional_return] = post_proc_main(path_or_struct , disabled_metrics, save_name, calcul_dist_subdyns, to_plot, basis_id)
% example: [vargout, app_mat, true_mat, additional_return]  = post_proc_main('.\new_results\spec_dynamic_cyl_just_tryingreturn_struct_just_trying.mat');
% another example:
% [vargout, app_mat, true_mat, additional_return]  = post_proc_main('.\new_results\long_cyl\return_struct_for_cylinder_higher_freq2covid_cylinder_higher_freq2_results.mat');

if nargin < 1 || isempty(path_or_struct), path_or_struct = '.\new_results\ex_struct'; end
if isa(path_or_struct, 'string') || isa(path_or_struct,'char')
    load(path_or_struct)
else
    return_struct = path_or_struct;
end
if nargin < 2 || isempty(disabled_metrics) , disabled_metrics = {}; end 
if nargin < 3 || isempty(save_name), save_name = ''; end
if nargin < 4 || isempty(calcul_dist_subdyns), calcul_dist_subdyns = true; end
if nargin < 5, to_plot = true; end
if nargin < 6 || isempty(basis_id), basis_id = 1; end
% Inputs:
%   1) path_or_struct: y = Dx; x_(t+1) = \sum(f_i C_i) x_t
%     Have the following fields: D,F_cell,c_i,x, true_mat, app_mat
%     D -        array (Neurons X num. components)
%     F_cell -   cell with array elements 
%                (num of sub-dynamics -> in each element -neurons X neurons)
%     c_i -      coeffiecients of dynamics {mat: [num sub dyns X time}
%     x -        matrix (num. components X time)
%     true_mat - true y (original recordings)
%     app_mat - optional. D*x;

%   1) metrics_to_run, string array. 
%       - R2
%       - MSE (app. vs g. truth)
%       - med_MSE
%       - std_MSE
%       - R2 (approximate vs real)
%       - med_R2 (app. vs real)
%       - mean_DTW (distance between app and ground truth) 

%       - tot_c_var (total_variance, supremum)
%       - mean_med_c_der_var (mean of median of the derivatives of each c)
%       - mean_med_c_wind (mean median running window (width 10))
%       - mean_c_var (Mean C variance)
%       - mean_c_zeros (mean fraction of zeros in c_variables)

%       - mean_f_zeros (in % of zeros in sub-dynamics)
%       - mean_f_var (mean variance inside sub-dynamics)

%   2)
%% Define params
D = return_struct.D; % -        array (Neurons X num. components)
F_cell = return_struct.F_cell; %-   cell with array elements 
%                (num of sub-dynamics -> in each element -neurons X neurons)
c_i = return_struct.c_i; % -      coeffiecients of dynamics {mat: [num sub dyns X time}
x = return_struct.x ; %-        matrix (num. components X time)
if isa(x,'cell'), x = x{basis_id}; end
true_mat = return_struct.true_mat; % true y (original recordings)
if length(size(true_mat)) == 3, true_mat = permute(true_mat, [1,3,2]); end
if ~isfield(return_struct,'app_mat')
    app_mat = D * x; 
    
else
    app_mat = return_struct.app_mat;
end % app = D_neur * A_cell{1};


%% Calculate R2
if ~any(strcmp(disabled_metrics,'R2')),  vargout.R2 = calcul_R2(app_mat(:),true_mat(:)); end

%% Calculate median R2 of time traces
if ~any(strcmp(disabled_metrics,'med_R2')),  vargout.med_R2 = calcul_R2(app_mat,true_mat, 1); end

%% Calculate MSE
if ~any(strcmp(disabled_metrics,'MSE')),  vargout.MSE = calcul_MSE(app_mat(:),true_mat(:)); end

%% Calculate median R2 of time traces (median of neurons)

if ~any(strcmp(disabled_metrics,'med_MSE')),  vargout.std_MSE = median(calcul_MSE(app_mat,true_mat, 2)); end

%% Calculate median R2 of time traces (median of neurons)
if ~any(strcmp(disabled_metrics,'var_MSE')),  vargout.med_MSE = std(calcul_MSE(app_mat,true_mat, 2)); end

%% Calculate mean DTW of time traces (median of neurons)
if ~any(strcmp(disabled_metrics,'mean_DTW')),  vargout.mean_DTW = mean(arrayfun(@(k) calcul_DTW(app_mat(k,:),true_mat(k,:)), 1:size(true_mat,1) )); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(c_i, 'cell'), c_i = c_i{basis_id}; end

[tot_var_med, tot_var_sup,mean_med_var_wind ] = calcul_all_var_sig(c_i, 2);
vargout.tot_var_sup = tot_var_sup; %       - tot_c_var (total_variance, supremum)
vargout.tot_var_med = tot_var_med; %       - mean_med_c_der_var (mean of median of the derivatives of each c)
vargout.mean_med_var_wind = mean_med_var_wind; %       - mean_med_c_wind (mean median running window (width 10))



%       - mean_c_std (Mean C std)
vargout.mean_var_c = mean(std(c_i,1),2);

%       - mean_c_zeros (mean fraction of zeros in c_variables)
vargout.mean_c_zeros = mean(c_i(:) == 0);



%%%%%%%%%%%%%%%%%%%%%%
%% study f
%%%%%%%%%%%%%%%%%%%%%
if isa(F_cell, 'cell'), F_mat = cell2mat(permute(F_cell,[1,3,2])); end

%       - mean_f_zeros (in % of zeros in sub-dynamics)
vargout.mean_F_zeros = mean(F_mat(:) == 0);

%       - mean_f_std (mean variance inside sub-dynamics)
vargout.mean_F_std = mean(cellfun(@(x) std(x(:)), F_cell));


if to_plot
    [~, path_to_save] = visualize_dyns(path_or_struct); 
    plot_dyn_with_time_err(true_mat,app_mat);
end
additional_return = plot_evals(F_cell, to_plot, [] );
[mse_vals, corr_vals, labels, angles] = distance_subdynamics(F_cell, c_i, to_plot);
additional_return.mse_vals = mse_vals;
additional_return.corr_vals = corr_vals;
additional_return.labels = labels;
additional_return.angles = angles;
%% Save
if to_plot
plot_struct(vargout)
if ~exist("path_to_save",'var')
    if isa(path_or_struct, 'string') || isa(path_or_struct,'char')          
       path_to_save = path_or_struct;
    else
        path_to_save = '.\new_results\';
    end
end
if isempty(save_name)
       name_full_to_save = [path_to_save,'results_mat.mat'];
else
       name_full_to_save = [path_to_save, sprintf('save_results%s.mat', save_name)];
end
save(name_full_to_save, "vargout")
end



end
%% Functions!

function  additional_return = plot_evals(F_cell, to_plot,  additional_return, colors)
    if size(F_cell{1},1) ~=3
        disp('Not 3d evals')
    else
        if nargin < 2, to_plot = true; end
        if nargin < 3 || isempty(additional_return), additional_return = struct; end
        %if nargin < 4, colors = rand(3, length(F_cell));end
        %colors(2,:) = 0.1;
        if nargin < 4, colors = 1: length(F_cell);end
        if to_plot, figure; end
        for f_spec = 1:numel(F_cell)
            cur_f = F_cell{f_spec};
            if to_plot
                
                % plot eigenvalues
                [evecs, evals_full] = eig(cur_f);
                evals = diag(evals_full);
                additional_return.evals.(sprintf('f_num_%s', num2str(f_spec))) = evals;
                subplot(1,2,1); scatter(real(evals),imag(evals),60, colors(f_spec)*ones(size(evals)),'filled');
                colormap('spring');  hold on;
                size(colors);
                subplot(1,2,2); scatter3(evecs(1, :),evecs(2, :),evecs(3, :),60, [colors(f_spec);colors(f_spec);colors(f_spec)],'filled');
                colormap('spring');
                hold on;
                patch(evecs(1, :),evecs(2, :),evecs(3, :),colors(f_spec),'FaceAlpha',0.2)
                colormap('spring');  hold on;
            end
        end
        if to_plot
            %legend(string(1:numel(F_cell))); 
            subplot(1,2,1);
            title('evals of different sub-dynamics')
            xline(1,'LineStyle','--','Alpha',0.5);
            yline(1,'LineStyle','--','Alpha',0);
            legend([string(arrayfun(@(k) sprintf('f%s',num2str(k)), 1:numel(F_cell), 'UniformOutput', false)), '']); 
            xlabel('real'); ylabel('img');
            subplot(1,2,2);
            title('evecs spaces by sub-dynamics')
            leg_string = [string(arrayfun(@(k) sprintf('f%s',num2str(k)), 1:numel(F_cell), 'UniformOutput', false));repelem([""],1,numel(F_cell))];
            legend(leg_string(:)); xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
    
        end
    end
end


function [mse_vals, corrs_vals, labels, angles] = distance_subdynamics(F_cell, c_i, to_plot, decide_f_metric)
% F_cell is a number-of-sub-dynamics cell array. Each element is a kXk
% matrix.
if nargin < 3, to_plot = true; end
if nargin < 4, decide_f_metric = 'R2'; end % can be also R2
combinations = nchoosek(1:numel(F_cell), 2);
mse_vals = [];
corrs_vals = [];
labels = string;
angles = [];
for combine= 1:size(combinations,1)
    combine1 = combinations(combine,1);
    combine2 = combinations(combine,2);
    f1 = F_cell{combine1};
    f2 = F_cell{combine2};
    c_i_1 = c_i(combine1,:);
    c_i_2 = c_i(combine2, :);
    if strcmpi(decide_f_metric,'mse')
        MSE = calcul_MSE(norm_signals(f1(:)), norm_signals(f2(:)));
    elseif strcmpi(decide_f_metric,'r2')
        MSE = calcul_R2(norm_signals(f1(:)), norm_signals(f2(:)));
    else
        error('Unknown metric for f comparison')
    end
    angle = angle_between(f1,f2);
    mse_vals = [mse_vals, MSE];
    R2 = calcul_R2(c_i_1,c_i_2);
    corrs_vals = [corrs_vals, R2];
    angles = [angles, angle];
    labels = [labels, string(sprintf('f%s & f%s', num2str(combine1), num2str(combine2)))];
    
end
labels = labels(2:end);
if to_plot
    figure;
    subplot(1,3,1);
    disp(size(labels))
    disp(size(angles))
    disp(labels)
    disp(angles)
    bar(categorical(labels), angles); title('angles between sub-dynamics')
    subplot(1,3,2);
    bar(categorical(labels), mse_vals); title([decide_f_metric, ' between sub-dynamics'])
    subplot(1,3,3);
    R2_sigs = calcul_R2(corrs_vals, mse_vals);
    scatter(corrs_vals, mse_vals,'filled'); 
    hold on;
    lsline; 
    title(sprintf('time correlations vs sub-dynamics mse (corr. %s)',num2str(R2_sigs))); 
    xlabel('Time traces correlations'); ylabel([decide_f_metric, ' between sub-dynamics'])


end
%% MSE


%% R2

%% evals

%% angle_between 
end
function angle = angle_between(A,B)
trace_AB = trace(A'*B);
deno = sqrt(sum(sum(A.^2))*sum(sum(B.^2)));
angle = acos(trace_AB/deno);
%% Angle
%trace(sum())
end

function sig_normalized = norm_signals(sig)
    sig_normalized = (sig - mean(sig)) / std(sig);
end
function MSE = calcul_MSE(app_mat, true_mat, ax )
% ax = 1 -> median of neurons
% ax = 2 -> median of time-points
    if nargin < 3, ax = NaN; end
    if isnan(ax)
        MSE = mean(mean((app_mat -true_mat).^2));
    else
        MSE =  mean((app_mat -true_mat).^2, ax);
    end
end
function R2 = calcul_R2(app_mat,true_mat, ax)
    if nargin < 3, ax = NaN; end
    if isnan(ax)
        R2 = calc_corr(app_mat(:) ,true_mat(:));
    elseif ax == 1 % compare time traces (r2 between time traces)
        R2 =  median(arrayfun(@(k) calc_corr(app_mat(k,:), true_mat(k,:)), 1:size(true_mat,1)));
    elseif ax == 2 % compare neurons 
        R2 =  median(arrayfun(@(k) calc_corr(app_mat(:,k), true_mat(:,k)), 1:size(true_mat,2)));
    end

end

function corr_return = calc_corr(x1,x2)
    corr_cur = corrcoef(x1(:),x2(:));
    corr_return = corr_cur(1,2);

end
function DTW = calcul_DTW(x1, x2)
    DTW = dtw(x1(:),x2(:));
end


%% Variance of Der Functions
function [tot_var_med, tot_var_sup,mean_med_var_wind ] = calcul_all_var_sig(signal, ax)
    abs_diff = calcul_list_derivatives(signal, ax );
    tot_var_med = median_variance(abs_diff, ax);
    tot_var_sup = total_variance(abs_diff, ax);
    mean_med_var_wind = med_wind_variance(abs_diff, ax);
end
function abs_der = calcul_list_derivatives(signal, ax )
    if nargin < 2, ax = NaN; end
    if isnan(ax)
        abs_der = abs(diff(signal)); 
    else
        abs_der = abs(diff(signal,[],ax)); 
    end
end
function tot_var_med = median_variance(abs_diff, ax)
    tot_var_med = mean(median(abs_diff, ax));
end
function tot_var_sup = total_variance(abs_diff, ax)
    %disp(size(abs_diff))
    %disp(size(max(abs_diff,[],ax)))
    %disp(ax)
    tot_var_sup = mean(max(abs_diff,[], ax));
end
function mean_med_var_wind = med_wind_variance(abs_diff, ax)
    mean_med_var_wind = mean(median(movvar(abs_diff, ax),ax));

end

function plot_struct(struct_to_plot)
fields = fieldnames(struct_to_plot);
values = arrayfun(@(k) struct_to_plot.(fields{k}), 1:numel(fields));
figure;
b = bar(categorical(fields),values,'FaceColor',[0 .5 .5]);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
set(gca, 'YScale', 'log')
ylabel('Metric in log scale')
title('Results Values')


end
function plot_dyn_with_time_err(real_data, appr_data)
    %disp(size(real_data))
    %error('fgdfg')
    % app_data: neuron X time
    % real data: neuron X time
    time_error = mean((real_data - appr_data).^2,1);
    time_error = time_error/max(time_error(:));
    %colors_error = [time_error, 1-time_error, time_error.^2];
    figure;
    if size(real_data,1) == 2
        plot3(real_data(1,:),real_data(2,:), real_data(3,:),'Color','black');
        scatter(real_data(1,:),real_data(2,:), [], time_error, 'filled');

    elseif size(real_data,1) == 3        
        plot3(real_data(1,:),real_data(2,:), real_data(3,:),'Color','black');
        hold on
        scatter3(real_data(1,:),real_data(2,:), real_data(3,:), [], time_error, "filled");
    else
        %plot(permute(real_data,[1,3,2])');
        plot(real_data');
        %hold on
        %scatter3(real_data(1,:),real_data(2,:), real_data(3,:), [], time_error, "filled");
    end
    colormap(gca,'winter');
    colorbar;
    title('dynamics with color coded temporal mse error');

    figure;
    if size(real_data,1) == 2
        %         subplot(1,2,1)
        %         plot3(real_data(1,:),real_data(2,:), real_data(3,:),'Color','black');
        %         scatter(real_data(1,:),real_data(2,:), [], time_error, 'filled');
        %         title('real'); colormap(gca,'winter');
        %subplot(1,2,2)
        plot3(appr_data(1,:),appr_data(2,:), appr_data(3,:),'Color','black');
        scatter(appr_data(1,:),appr_data(2,:), [], time_error, 'filled');
        title('app')
    elseif size(real_data,1) == 3
        %         subplot(1,2,1)
        %         plot3(real_data(1,:),real_data(2,:), real_data(3,:),'Color','black');
        %         hold on
        %         scatter3(real_data(1,:),real_data(2,:), real_data(3,:), [], time_error, "filled");
        %         title('real'); colormap(gca,'winter');
        %         subplot(1,2,2)
        plot3(appr_data(1,:),appr_data(2,:), appr_data(3,:),'Color','black');
        hold on
        scatter3(appr_data(1,:),appr_data(2,:), appr_data(3,:), [], time_error, "filled");
        title('app')
    else
        %         subplot(1,2,1)
        %         %plot(permute(real_data,[1,3,2])');
        %         plot(real_data'); title('real')
        %         subplot(1,2,2)
        plot(appr_data');title('app')
    end
    colormap(gca,'winter');
    colorbar;
    %title('dynamics with color coded temporal mse error');

end