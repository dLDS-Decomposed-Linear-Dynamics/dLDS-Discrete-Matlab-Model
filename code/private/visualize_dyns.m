function [argout,path_save] = visualize_dyns(path, name_file, cal_comp, to_plot_results, real_mat, dynamic_type)
% cal_comp can be true = calculate each component seperately
%                 false - calculate all components
%                 cumulative - find components importance using the coefficients from A_cell and find the cumulative results  
% example call for covid:
% argout = visualize_dyno_new('C:\JHU COURSES\מאמר על קופים - קוד ודאטא\דאטא של ביגאן\code\COVID_DATA\covid_5_subs\','covid_conf_new_idcovid_new_id_results.mat' , [], true, merged_diffs)
% argout = visualize_dyno_new('C:\JHU COURSES\מאמר על קופים - קוד ודאטא\דאטא של ביגאן\code\results_dynamics\','spec_dynamic_cyl_lolocovid_lolo_results.mat' , [], true, X_ex)
% Different example type:
%
cur_time = replace(replace(replace(char(datetime('now')),':','_'),' ', '_'),'-','_');
%mkdir(cur_time)
%% Call file

if nargin < 1 || isempty(path), path = [pwd,'\res_res\']; end
if nargin < 6 || isempty(dynamic_type), dynamic_type = ''; end
if  (isa(path,'string')  || (isa(path, 'char')) && endsWith(path,'\'))
    addpath(path);
    if nargin < 2 || isempty(name_file)
        name_file = 'block_size_12_N_ex_150_bijan_ooobin_20_overlap_18_results.mat' ;
    end
    load([path, 'A_cell_B_cell.mat']); load([path,name_file]);
    if nargin < 5
        load([path, 'X_ex_bins_20_overlap_18.mat'],'X_ex');
    else 
        X_ex = real_mat;
    end
    dynamic_name = dynamic_type;
else
    if isa(path,'struct')
        return_struct = path; 
        path = '.\new_results\'; 
    else 
      
        try
            load(path,'return_struct'); 
        catch 
            error('return struct does not appear in file'); 
        end
    end
    %% Define params
    D_neur = return_struct.D; % -        array (Neurons X num. components)
    F_neur = return_struct.F_cell; %-   cell with array elements (num of sub-dynamics -> in each element -neurons X neurons)
    B_cell = return_struct.c_i; % -      coeffiecients of dynamics {mat: [num sub dyns X time}
    A_cell = return_struct.x ; %-        matrix (num. components X time)
    X_ex = return_struct.true_mat; % true y (original recordings)
    if length(size(X_ex)) == 3, X_ex = permute(X_ex, [1,3,2]); end
    if endsWith(path,'.mat'), path = path(1:end-4); end

    if isfield(return_struct, 'sub_dyn_spec')
        dynamic_name = [return_struct.sub_dyn_spec, dynamic_type];
    else
        dynamic_name = dynamic_type;
    end


end
if nargin < 3 || isempty(cal_comp), cal_comp = ''; end
if nargin < 4 || isempty(to_plot_results),     to_plot_results = true;end



%% Analysis
if isa(cal_comp , "char")
    path_save = [path, cal_comp,'\',dynamic_name,cur_time,'\'];
    if ~exist(path_save, 'dir'), mkdir(path_save); end
    [ordered_D, ordered_A,order_coeffs] = find_features_order(D_neur, A_cell, path_save);
    struct_results = calculate_error(X_ex, ordered_D, ordered_A, cal_comp);
else
    
    path_save = [path, num2str(cal_comp),'\',dynamic_name,cur_time,'\'];
    if ~exist(path_save, 'dir'), mkdir(path_save); end
    struct_results = calculate_error(X_ex, D_neur, A_cell, cal_comp);
end




%% Plot heatmap of dynamics coefficients
figure;
B_cell1 = B_cell{1};
heatmap(B_cell1(:,1:min(size(B_cell1,2),1000)));
colormap winter
xlabel('Time')
ylabel('# dynamics')
title('Coefficients heatmap until 1000 time points')
saveas(gcf,[path_save, 'dynamics_coeffs_heatmap.png'])
%% Scatter plot of coefficeients 
figure;
hold off; 
colors = rand(numel(F_neur),3);
for i = 1:numel(F_neur)
scatter(1:size(B_cell{1},2), B_cell{1}(i,:),18,colors(i,:),'filled'); hold on;plot(1:size(B_cell{1},2), B_cell{1}(i,:),'Color', colors(i,:));
end
xlabel('Time')
ylabel('Dynamics coefficients')
saveas(gcf,[path_save, 'dynamics_coeffs_scatter.png'])
%% display dynamics
figure; 
for i = 1:numel(F_neur)
    subplot(1,numel(F_neur),i)
    %imshow(F_neur{i}/max(max(F_neur{i})));
    heatmap(F_neur{i})
    title(['F #', num2str(i)]);
    saveas(gcf,[path_save, 'subdynamics.png'])
end
%% Show the D dyn
figure;
subplot(1,2,1); imshow(D_neur/max(max(D_neur))); title('D');
subplot(1,2,2); imshow(corr(D_neur)); title('Correlation between D columns');
saveas(gcf,[path_save, 'D_and_D_corr.png'])




%% MSE


%% ttest for each comp
% mse for each comp
% find the most dominant column in D



%mse for each components

if ~isa(cal_comp , "char")
    argout = struct_results;
else

    argout = {struct_results, order_coeffs, ordered_D, ordered_A};
end

if to_plot_results
    plot_results(struct_results, cal_comp, true, [path_save, 'd_cols_important.png'])

end
end

%% reco vs real
function plot_diff(X_ex, aap, path_save, max_time) 
    if nargin < 4, max_time = 100; end
    aap_norm = (aap-min(aap))./(max(aap)-min(aap));
    norm_xex = (X_ex - min(X_ex))./(max(X_ex)-min(X_ex));
    figure;
    subplot(2,1,2); imshow(norm_xex(1:size(aap,1),1:max_time));
    %sprintf('Real Data (first %s frames)',num2str(max_time))
    title(sprintf('Real Data (first %s frames)',num2str(max_time)));

    subplot(2,1,1); imshow(aap_norm(1:size(aap,1),1:max_time));
    title(sprintf('Approximated Data (first %s frames)',num2str(max_time)));
    
    saveas(gcf,[path_save, 'real_vs_app.png'])
end


function plot_reco(reco, tot_num, i, path_save) 
    norm_xex = (reco - min(reco))./(max(reco)-min(reco));
    subplot(1,tot_num,i); imshow(norm_xex)
    title('Approximated Data')
    saveas(gcf,[path_save, 'real_vs_app_per_comp.png'])
end

function [ordered_D, ordered_A, I] = find_features_order(D, A, path_save)
    % A: # components X time
    % D: neurons X # components
    if isa(A, 'cell')
        A = A{1};
    end
    
    sum_per_comp = sum(abs(A),2);
    [B,I] = sort(sum_per_comp,'descend');
    ordered_D = D(:,I);
    ordered_A = A(I,:);
    figure(50);
    stem(B); title('Order coeffs A'); xlabel('cumulative # comps'); ylabel('sum coefficients for component')
    saveas(gcf,[path_save, 'coeffiecients_important_by_A.png'])
end

%% Functions 
function struct_results = calculate_error(real_data, D_neur, A_cell, each_col, cal_r, path_save)

    struct_results = struct('mse',[],'ttest_p',[],'R2',[],'var_ex',[]);
    if nargin < 4
        each_col = true;
    end
    if nargin < 5
        cal_r = false;
    end
    if nargin < 6
        path_save = '.';
    end
    if isa(A_cell, "cell")  
        A_cell = A_cell{1};
    end
    % coefficient of partial determination matlab
    real_data = real_data(1:size(D_neur,1),:);
    if isa(each_col, "char") && strcmpi(each_col, 'cumulative')
     
        figure(15);
        ordered_D= D_neur; ordered_A = A_cell;

        for idx = 1:size(D_neur,2)
            app_data = ordered_D(:,1:idx) * ordered_A(1:idx, :);
            %plot_reco(app_data, size(D_neur,2), idx, path_save) ;
            [mse, ttest_p, R_squares, explained_var_result] = calcul_results_each_comp(real_data, app_data, ordered_D, idx, cal_r);
            struct_results.mse(end+1) = mse;
            %ttest_p
            struct_results.ttest_p(end+1) = nanmean(ttest_p);
            struct_results.R2(end+1) = nanmean(R_squares);
            struct_results.var_ex(end+1) = mean(explained_var_result);
            %struct_results

        end

    elseif each_col
        figure();
        sum_per_comp = sum(abs(A_cell),2);
        stem(sum_per_comp);
        title('coeffs A'); xlabel('# comp'); ylabel('sum coefficients for component');
        %saveas(gcf,[path_save, 'coeffiecients_important_by_A.png'])
        figure(15);
        for idx = 1:size(D_neur,2)
            app_data = D_neur(:,idx) * A_cell(idx, :);
            % plot_reco(app_data, size(D_neur,2), idx, path_save) ;
            [mse, ttest_p, R_squares, explained_var_result] = calcul_results_each_comp(real_data, app_data, D_neur, idx, cal_r);
            struct_results.mse(end+1) = mse;
            %ttest_p
            struct_results.ttest_p(end+1) = nanmean(ttest_p);
            struct_results.R2(end+1) = nanmean(R_squares);
            struct_results.var_ex(end+1) = mean(explained_var_result);
            %struct_results

        end
    else
            app_data = D_neur * A_cell;
            plot_diff(real_data, app_data, path_save) 
            [mse, ttest_p, R_squares, explained_var_result] = calcul_results_each_comp(real_data, app_data, D_neur, 1:size(D_neur,2), cal_r );
            %plot_dyn_with_time_err(app_data, real_data)
            struct_results.mse(end+1) = mse;
            %display(ttest_p)
            %display(size(struct_results.ttest_p))
            struct_results.ttest_p(end+1) = mean(ttest_p);
            struct_results.R2(end+1) = mean(R_squares);
            struct_results.var_ex(end+1) = mean(explained_var_result);
    end
    
     
    

end
%%
function [mse, ttest_p] = calculate_error_spec(real_data, app_data)
    mse = immse(real_data, app_data);
    [~, ttest_p ] = ttest(real_data, app_data); 
end

function explained_var_result = explained_var(real_data, app_data)
   diff_r = real_data - app_data;
   explained_var_result = 1 - var(diff_r(:))./var(real_data(:));
end
% covariance matrix of coefficients
% coefficient of partial determination matlab

function [mse, ttest_p, R_squares, explained_var_result] = calcul_results_each_comp(real_data, app_data, D_neur, idx, cal_r)
    [mse, ttest_p] = calculate_error_spec(real_data, app_data);

    %[rho, pval] = partialcorr(x,z);
    R_squares = [];
    if cal_r
    for time_point = 1:size(real_data,2)
        mdl = fitlm(D_neur(:,idx), real_data(1:400,time_point));
        R2= mdl.Rsquared.Ordinary;
        R_squares(end+1) = R2;
    end
    end
    explained_var_result = explained_var(real_data, app_data);
    
end


