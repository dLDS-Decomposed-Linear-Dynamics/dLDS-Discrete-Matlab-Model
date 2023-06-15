%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing of MRI recovery
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup Parameters

comp_opt = 'Desk';
samp_opt = 'Cols';
noise_opt = 0;
plot_opt = 0;
reset_samp = 1;
same_samp = 0;
calc_l1 = 1;
calc_rwl1 = 1;
calc_dl1 = 1;
calc_dl11 = 0;
calc_drwl1 = 1;
calc_drwl1b = 0;


% Add paths and define data folder name
if strcmp(comp_opt, 'Desk')
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/Wavelab850/'))
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/YUV/'))
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/Noiselet_Code/'))
    addpath(genpath('~/Desktop/Versioned/CodeBaseRepo/trunk/Wavelet_Code/'))
    addpath(genpath('~/Desktop/Versioned/2011_RWL1DCS/trunk/'))
elseif strcmp(comp_opt, 'Lap')
    addpath(genpath('~/Versioned/CodeBaseRepo/trunk/Wavelab850/'))
    addpath(genpath('~/Versioned/2011_RWL1DCS/trunk/'))
    addpath(genpath('~/Versioned/CodeBaseRepo/trunk/SparsityToolboxes/toolbox_sparsity/toolbox/'))
    folder_name = '/home/adam/Versioned/2011_RWL1DCS/trunk/MRIdata/';
else
    error('Invalid computer!')
end

rmpath /home/adam/Desktop/Versioned/CodeBaseRepo/trunk/Wavelab850/Books/WaveTour/WTCh06/

% Number of frames to analyze

T_s = 40;
t_s = T_s;
num_frames = T_s;
% Sampling rate
samp_factor = 0.25;
noise_var = 0.0001;

% General Parameters
TOL = 1e-3;
XFM = Wavelet('Daubechies',4,4);	% Wavelet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data
num_trials = 5;
for trial_index = 1:num_trials

if reset_samp == 1
% Load video
vid_t = yuv_import('/home/adam/Desktop/GT_Work/CodeBase_v1_0/Data/Videos/revideosequence/foreman.qcif', [176 144], T_s*num_trials);

for kk = 1:T_s
    vid{kk} = vid_t{T_s*(trial_index-1)+kk};
end
clear vid_t

% Extract the signals of interest
for index = 1:T_s
    vid{index} = vid{index}(1:128, 1:128);
end
% Select ratio of measurements
P_m = samp_factor;

N = numel(vid{1});
M = ceil(P_m*N);

vid2 = zeros(size(vid{1}, 1), size(vid{1}, 2), T_s);
vid3 = zeros(M, 1, T_s);

% Noise options
OM = zeros(M, T_s);

for index = 1:T_s
    % Set up noiselets
    q = randperm(N)';    % makes column vector of random integers 1:N
    OM(:, index) = q(1:M);          % vector of random subset of integers 1:N
    Phi  = @(z) A_noiselet (z, OM(:, index));
    Phit = @(z) At_noiselet(z, OM(:, index), N); 
    A  = @(z) Phi( z );
    At = @(z) Phit(z );

    vid2(:, :, index) = 2*(vid{index}/255 - 0.5); % Renormalize to [-1, 1]
    vid3(:, :, index) = A(reshape(vid2(:, :, index), [], 1)) + sqrt(noise_var)*randn(M, 1);
end

clear vid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simple Model Reconstruction YALL1
if calc_l1 == 1
	% YALL1 options
	clear opts
	opts.tol = TOL;
	lambda_val = 0.01;
	for kk = 1:num_frames
	    tic
	    Phi  = @(z) A_noiselet (z, OM(:, kk));
            Phit = @(z) At_noiselet(z, OM(:, kk), N); 
	    Af = @(x) Phi(reshape(XFM'*reshape(x, sqrt(N), sqrt(N)), [], 1));
	    Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1);
	    A = linop_handles([M, N], Af, Ab, 'R2R');
            res = solver_L1RLS( A, vid3(:, :, kk), lambda_val, zeros(N, 1), opts );
	    im_res = XFM'*reshape(res, sqrt(N), sqrt(N));
	    if plot_opt == 1
		% Plot reconstruction
		subplot(1, 2, 1), imshow((real(im_res)),[]), drawnow
		title('Reconstructed', 'FontSize', 25)
		subplot(1, 2, 2), imshow(aheart_image(:, :, kk),[]), drawnow
		title('Actual', 'FontSize', 25)
	    end
	    
	    % Save reconstruction results
	    vid_coef_cs(:, :, kk) = res;
	    vid_recon_cs(:, :, kk) = im_res;
	    vid_rMSE_cs(kk) = sum(sum((vid_recon_cs(:, :, kk) - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
	    vid_PSNR_cs(kk) = psnr(real(vid_recon_cs(:, :, kk)), vid2(:, :, kk), 1);
	    TIME_ITER = toc;
	    fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, vid_PSNR_cs(kk), vid_rMSE_cs(kk))
	end
end

fprintf('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Re-Weighted Model Reconstruction YALL1
if calc_rwl1 == 1
	rwl1_reg = 0.1;
	rwl1_mult = 0.1;

	for kk = 1:num_frames
	    clear opts
	    opts.tol = TOL;
	    lambda_val = 0.001;

	    tic
	    Phi  = @(z) A_noiselet (z, OM(:, kk));
            Phit = @(z) At_noiselet(z, OM(:, kk), N); 
	    Af = @(x) Phi(reshape(XFM'*reshape(x, sqrt(N), sqrt(N)), [], 1));
	    Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1);
	    A = linop_handles([M, N], Af, Ab, 'R2R');
	    
	    weights = 1;
	    for nn = 1:10
		res = solver_L1RLS( A, vid3(:, :, kk), lambda_val, zeros(N, 1), opts );
		res = res./weights;
		weights = rwl1_mult./(abs(res) + rwl1_reg);
                Af = @(x) Phi(reshape(XFM'*reshape(x./weights, sqrt(N), sqrt(N)), [], 1));
	        Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1)./weights;
	    	A = linop_handles([M, N], Af, Ab, 'R2R');

		im_res = XFM'*reshape(res, sqrt(N), sqrt(N));
		if plot_opt == 1
		    subplot(1, 2, 1), imshow(im_res,[]), drawnow
		    title('Reconstructed', 'FontSize', 25), drawnow
		    subplot(1, 2, 2), imshow(vid2(:, :, kk),[]), drawnow
		    title('Actual', 'FontSize', 25), drawnow
		end
	    	temp_rMSE = sum(sum((im_res - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
		fprintf('Finished RW iteration %d. rMSE is %f.\n', nn, temp_rMSE)
	    end
	    
	    % Save reconstruction results
	    vid_coef_rwcs(:, :, kk) = res;
	    vid_recon_rwcs(:, :, kk) = im_res;
	    vid_rMSE_rwcs(kk) = sum(sum((vid_recon_rwcs(:, :, kk) - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
	    vid_PSNR_rwcs(kk) = psnr(real(vid_recon_rwcs(:, :, kk)), vid2(:, :, kk), 1);
	    TIME_ITER = toc;
	    fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, vid_PSNR_rwcs(kk), vid_rMSE_rwcs(kk))
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simple Dynamic Model Reconstruction YALL1
if calc_dl1 == 1
	clear opts
	opts.tol = TOL;

	lambda_history = 0.4;
        lambda_val = 0.01;
	% Solve for initial frame
	Phi  = @(z) A_noiselet (z, OM(:, 1));
        Phit = @(z) At_noiselet(z, OM(:, 1), N); 
	Af = @(x) Phi(reshape(XFM'*reshape(x, sqrt(N), sqrt(N)), [], 1));
	Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1);
	A = linop_handles([M, N], Af, Ab, 'R2R');

        res = solver_L1RLS( A, vid3(:, :, 1), lambda_val, zeros(N, 1), opts );

	im_res = XFM'*reshape(res, sqrt(N), sqrt(N));
	if plot_opt == 1
	    % Plot reconstruction
	    subplot(1, 2, 1), imshow(im_res,[]), drawnow
	    title('Reconstructed', 'FontSize', 25)
	    subplot(1, 2, 2), imshow(vid2(:, :, 1),[]), drawnow
	    title('Actual', 'FontSize', 25)
	end
	% Save reconstruction results
	vid_coef_dcs(:, :, 1) = res;
	vid_recon_dcs(:, :, 1) = im_res;
	vid_rMSE_dcs(1) = sum(sum((vid_recon_dcs(:, :, 1) - vid2(:, :, 1)).^2))/sum(sum(vid2(:, :, 1).^2));
	vid_PSNR_dcs(1) = psnr(real(vid_recon_dcs(:, :, 1)), vid2(:, :, 1), 1);
	TIME_ITER = toc;
	fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', 1, num_frames, TIME_ITER, vid_PSNR_dcs(1), vid_rMSE_dcs(1))

	% Solve for rest of frames
	for kk = 2:num_frames
	    tic
	    Phi  = @(z) A_noiselet (z, OM(:, kk));
            Phit = @(z) At_noiselet(z, OM(:, kk), N); 
	    Af = @(x) [lambda_history*x; Phi(reshape(XFM'*reshape(x, sqrt(N), sqrt(N)), [], 1))];
	    Ab = @(x) reshape((XFM)*reshape(Phit(x(N+1:end)), sqrt(N), sqrt(N)), [], 1) + lambda_history*x(1:N);
	    A = linop_handles([M+N, N], Af, Ab, 'R2R');
            
	    res = solver_L1RLS( A, [lambda_history*vid_coef_dcs(:, :, kk-1); ...
		vid3(:, :, kk)], lambda_val, zeros(N, 1), opts );
	    
	    im_res = XFM'*reshape(res, sqrt(N), sqrt(N));
	    if plot_opt == 1
		% Plot recovery
		subplot(1, 2, 1), imshow(im_res,[]), drawnow
		title('Reconstructed', 'FontSize', 25)
		subplot(1, 2, 2), imshow(vid2(:, :, kk),[]), drawnow
		title('Actual', 'FontSize', 25)
	    end
	    % Save reconstruction results
	    vid_coef_dcs(:, :, kk) = res;
	    vid_recon_dcs(:, :, kk) = im_res;
	    vid_rMSE_dcs(kk) = sum(sum((vid_recon_dcs(:, :, kk) - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
	    vid_PSNR_dcs(kk) = psnr(real(vid_recon_dcs(:, :, kk)), vid2(:, :, kk), 1);
	    TIME_ITER = toc;
	    fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, vid_PSNR_dcs(kk), vid_rMSE_dcs(kk))
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dynamic Re-Weighted Model Reconstruction YALL1

if calc_drwl1 == 1
	opts.tol = TOL;
	lambda_val = 0.001;
	rwl1_reg = 0.1;
	rwl1_mult = 0.1;
	w1 = 0.5;
	w2= 1-w1;

	kk = 1;
	tic
	Phi  = @(z) A_noiselet (z, OM(:, 1));
        Phit = @(z) At_noiselet(z, OM(:, 1), N); 
	Af = @(x) Phi(reshape(XFM'*reshape(x, sqrt(N), sqrt(N)), [], 1));
	Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1);
	A = linop_handles([M, N], Af, Ab, 'R2R');

	for nn = 1:10    
            res = solver_L1RLS( A, reshape(vid3(:, :, 1), [], 1), lambda_val, zeros(N, 1), opts );
	    res = res./weights;
	    weights = 1./(abs(real(res)) + rwl1_reg);
            Af = @(x) Phi(reshape(XFM'*reshape(x./weights, sqrt(N), sqrt(N)), [], 1));
	    Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1)./weights;
	    A = linop_handles([M, N], Af, Ab, 'R2R');
	    
	    im_res = XFM'*reshape(res, sqrt(N), sqrt(N));
	    if plot_opt == 1
		% Plot reconstruction
		subplot(1, 2, 1), imshow((real(im_res)),[]), drawnow
		title('Reconstructed', 'FontSize', 25), drawnow
		subplot(1, 2, 2), imshow(aheart_image(:, :, kk),[]), drawnow
		title('Actual', 'FontSize', 25), drawnow
	    end
	    disp([mean(weights), var(weights)])
	    temp_rMSE = sum(sum((im_res - vid2(:, :, 1)).^2))/sum(sum(vid2(:, :, 1).^2));
	    fprintf('Finished RW iteration %d, rMSE = %f.. \n', nn, temp_rMSE)
	end

	% Save results
	vid_coef_drw(:, :, 1) = res;
	vid_recon_drw(:, :, 1) = im_res;
	vid_rMSE_drw(1) = sum(sum((vid_recon_drw(:, :, kk) - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
	vid_PSNR_drw(1) = psnr(vid_recon_drw(:, :, 1), vid2(:, :, 1), 1);
	TIME_ITER = toc;
	fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', 1, num_frames, TIME_ITER, vid_PSNR_drw(1), vid_rMSE_drw(1))

	for kk = 2:num_frames
	    clear opts
	    opts.tol = TOL;
	    tic
	    
	    Phi  = @(z) A_noiselet (z, OM(:, kk));
            Phit = @(z) At_noiselet(z, OM(:, kk), N); 
	    Af = @(x) Phi(reshape(XFM'*reshape(x, sqrt(N), sqrt(N)), [], 1));
	    Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1);
	    A = linop_handles([M, N], Af, Ab, 'R2R');

	    for nn = 1:10
		res = solver_L1RLS( A, vid3(:, :, kk), lambda_val, zeros(N, 1), opts );
	        res = res./weights;
		weights = 2*rwl1_mult./(w1*abs(real(res)) + w2*abs(reshape(real(vid_coef_drw(:, :, kk-1)), [], 1)) + rwl1_reg);
    		Af = @(x) Phi(reshape(XFM'*reshape(x./weights, sqrt(N), sqrt(N)), [], 1));
	        Ab = @(x) reshape((XFM)*reshape(Phit(x), sqrt(N), sqrt(N)), [], 1)./weights;
	        A = linop_handles([M, N], Af, Ab, 'R2R');

		im_res = XFM'*reshape(res, sqrt(N), sqrt(N));
		if plot_opt == 1
		    % Plot reconstruction
		    subplot(1, 2, 1), imshow(im_res,[]), drawnow
		    title('Reconstructed', 'FontSize', 25), drawnow
		    subplot(1, 2, 2), imshow(vid2(:, :, kk),[]), drawnow
		    title('Actual', 'FontSize', 25), drawnow
		    fprintf('Finished RW iteration %d. \n', nn)
		end
	        temp_rMSE = sum(sum((im_res - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
	        fprintf('Finished RW iteration %d, rMSE = %f.. \n', nn, temp_rMSE)
	    end

	    vid_coef_drw(:, :, kk) = res;
	    vid_recon_drw(:, :, kk) = im_res;
	    vid_rMSE_drw(kk) = sum(sum((vid_recon_drw(:, :, kk) - vid2(:, :, kk)).^2))/sum(sum(vid2(:, :, kk).^2));
	    vid_PSNR_drw(kk) = psnr(real(vid_recon_drw(:, :, kk)), vid2(:, :, kk), 1);
	    TIME_ITER = toc;
	    fprintf('Finished frame %d of %d in %f seconds. PSNR is %f. rMSE is %f. \n', kk, num_frames, TIME_ITER, vid_PSNR_drw(kk), vid_rMSE_drw(kk))
	end
end

PSNR_ALL(:, :, trial_index) = [vid_PSNR_cs(:), vid_PSNR_rwcs(:), vid_PSNR_dcs(:), vid_PSNR_drw(:)];
rMSE_ALL(:, :, trial_index) = [vid_rMSE_cs(:), vid_rMSE_rwcs(:), vid_rMSE_dcs(:), vid_rMSE_drw(:)];


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting

if plot_opt == 1
    figure;
    for kk = 1:20*num_frames
        IX = 1+mod(kk, num_frames-1);
%         subplot(1, 2, 1), imagesc(aheart_recon_dyn(:, :, IX))
%         subplot(1, 2, 1), imagesc(aheart_drwrecon(:, :, IX))
        subplot(1, 2, 1), imagesc(aheart_image(:, :, IX))
%         clims = [min(min(aheart_recon(:, :, IX))), max(max(aheart_recon(:, :, IX)))];
%         subplot(1, 2, 1), imagesc(aheart_recon(:, :, IX),clims)
        title('Reconstructed', 'FontSize', 25)
        
        colormap gray
        axis image
        axis off
%         subplot(1, 2, 2), imagesc(aheart_recon(:, :, IX))
%         subplot(1, 2, 2), imagesc(aheart_recon_dyn(:, :, IX))
        clims = [min(min(aheart_drwrecon(:, :, IX))), max(max(aheart_drwrecon(:, :, IX)))];
        subplot(1, 2, 2), imagesc(aheart_drwrecon(:, :, IX), clims)
        title('Actual', 'FontSize', 25)
        colormap gray
        axis image
        axis off
        pause(0.5)
    end

    subplot(1, 1, 1), plot([aheart_psnr(:), aheart_rwpsnr(:), aheart_psnr_dyn(:), aheart_drwpsnr(:)], 'LineWidth', 3)
    xlabel('Iteration', 'FontSize', 15)
    ylabel('PSNR', 'FontSize', 15)
    legend('Independent BPDN', 'Independent RWL1', 'Simple Dynamics', 'Dynamic RWL1')
else
    % disp([vid_PSNR_cs(:), vid_PSNR_rwcs(:), vid_PSNR_dcs(:), vid_PSNR_drw(:)])
    % disp([vid_rMSE_cs(:), vid_rMSE_rwcs(:), vid_rMSE_dcs(:), vid_rMSE_drw(:)])
    disp(mean(PSNR_ALL, 3))
    disp(mean(rMSE_ALL, 3))
end


save /home/adam/Desktop/TRIAL_X5i.mat


% sendmail('adamshch@gmail.com', 'TEST IS DONE', 'TEST IS DONE', {'/home/adam/Desktop/TRIAL_X5.mat'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
