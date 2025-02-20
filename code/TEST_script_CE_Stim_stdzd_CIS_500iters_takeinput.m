%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% List of tasks:
% Rate estimation?
% Basic dimensionality reduction?
% Dynamics estimation?

% parpool(16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add all the required packages that are custom code
addpath(genpath('.'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up pointers to the fluorescence data and meta behavioral data

% dat.dir   = 'C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegans/Data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data

% dat.dFF = getNPYarray([dat.dir, dat.Ffile]);                               % Time x Neuron matrix of fluorescent activity
% dat.bhv = getNPYarray([dat.dir, dat.Bfile]);                               % Time x 1 vector of behavioral data
% dat.vel = getNPYarray([dat.dir, dat.Vfile]);                               % Time x 1 vector of velocity

dat = load('./CElegansData/WT_Stim.mat');
% dat = load('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegansData/WT_Stim.mat');
% dat2 = load('C:/Users/Helen/Documents/Documents/Documents (3)/GradSchool/JohnsHopkinsStudent/CharlesLab/CElegansData/WT_NoStim.mat');

thisWorm = input("Which Stim worm:");
ifNoAVA = false;

disp('Remember to change filenames below with version number')
disp('Did you start the diary?')

dFF_S_thisWorm = dat.WT_Stim(thisWorm).traces;
bhv_S_thisWorm = dat.WT_Stim(thisWorm).States;

% for plotting reconstruction later
ids = dat.WT_Stim(thisWorm).IDs;
tVals    = dat.WT_Stim(thisWorm).timeVectorSeconds; %note: different spelling in NoStim data
% select named neuron indices
identifiedNeurons = find(cellfun(@isempty,ids)==0);
unidentifiedNeurons = find(cellfun(@isempty,ids)==1);
numChannels = size(ids,2);
channelNamesStrings = getChannelNamesAsStrings(ids,numChannels,identifiedNeurons);

locAVAL = regexp(channelNamesStrings,'AVAL');
locAVAR = regexp(channelNamesStrings,'AVAR');
emptyCellsInds = cellfun(@isempty,locAVAL);
locAVAL = find(~emptyCellsInds);
emptyCellsInds = cellfun(@isempty,locAVAR);
locAVAR = find(~emptyCellsInds);

%%

dFF = dFF_S_thisWorm;
bhv = bhv_S_thisWorm;

if ifNoAVA
    dFF(:,locAVAR) = zeros(size(dFF(:,1))); dFF(:,locAVAL) = zeros(size(dFF(:,1))); % remove confound AVA but don't get rid of columns
    disp('zeroed out AVAL and R')
end

dFF(dFF<0)=0; % remove negative values
dFF = dFF./max(dFF,[],'all');

% figure(1);                                % ASC added 7/25
% subplot(2,2,[1,2]), imagesc(dFF.')        % ASC added 7/25
% subplot(2,2,3), imagesc(Fgt{1})           % ASC added 7/25
% subplot(2,2,4), imagesc(Fgt{2})           % ASC added 7/25
%%

% Set parameters
inf_opts.nF              = input('nF:');
inf_opts.N               = input('n:'); %scaled down from 100 (zebrafish, whole brain)
inf_opts.M               = size(dFF,2) ;
disp('If reset dFF below, make sure to change inf_opts.M');
inf_opts.lambda_val      = 0.1; % VARY - tradeoff: approx. SNR - see ratios of traces
inf_opts.lambda_history  = 0.90;
inf_opts.lambda_b        = 0.01;
inf_opts.lambda_historyb = 0;
inf_opts.tol             = 1e-3; % 1e-3
inf_opts.max_iter2       = 10; %500 
inf_opts.max_iters       = 500;
inf_opts.special         = '';
inf_opts.F_update        = true;
inf_opts.D_update        = true;
inf_opts.solver_type     = 'tfocs';
% inf_opts.N_ex            = 40;
% inf_opts.T_s             = 30;
% inf_opts.special         = 'noobs';

% inf_opts.step_d          = 20;  % ASC added 7/25
% inf_opts.step_f          = 30;  % ASC added 7/25
% inf_opts.plot_option     = 20;  % ASC added 7/25
% inf_opts.lambda_f        = 0.05; % ASC added 7/28
%%
% diary CEDiary_Stim1_stdzd_500iters_NoAVA_v01

[Phi, F] = bpdndf_dynamics_learning(dFF.', [], [], inf_opts);              % Run the dictionary learning algorithm


sizedFF = size(dFF.');
dFFcell = mat2cell(dFF.', sizedFF(1), sizedFF(2));
[A_cell,B_cell] = parallel_bilinear_dynamic_inference(dFFcell, Phi, F, ...
                                       @bpdndf_bilinear_handle, inf_opts); % Infer sparse coefficients


%%
x_values = A_cell{1};
dataReconstructed = (Phi * x_values).';    


neurSel = randsample(identifiedNeurons, 10);
maxvalDFF = max(dFF_S_thisWorm,[],'all');
stdzdOn = true;
if stdzdOn
    dataReconstructedRescaled = maxvalDFF .*(Phi * x_values).';
    dFFRescaled = maxvalDFF .* dFF;
end
% dFF = dFFRescaled;
% dataReconstructed = dataReconstructedRescaled;
% ifHighlight = false;

flatDFF = reshape(dFF, numel(dFFRescaled),[]);
flatReconstructed = reshape(dataReconstructedRescaled,numel(dataReconstructedRescaled),[]);
rval = corrcoef(flatDFF, flatReconstructed);
varExpl = (rval(1,2))^2;
disp(varExpl);  
diary off

versionSave = input('Random string:','s');
save(sprintf('TEST_saveCEWorkspace_Stim%d_stdzd_500iters_NoAVA%d_N%d_nF%d_v%s.mat', thisWorm, ifNoAVA, inf_opts.N, inf_opts.nF,versionSave) ); 

% %%
% fig03 = figure(3);
% subplot(2,1,1),cla;
% hold on;
% bottom = 0;
% top = max(dFFRescaled,[],'all');
% if ifHighlight
%     for i = 1:size(vertOnSuspect)
%         vertices = [vertOnSuspectTime(i) bottom; vertOffSuspectTime(i) bottom; vertOffSuspectTime(i) top; vertOnSuspectTime(i) top];
%         faces = [1     2     3     4];
%         patch('Faces',faces,'Vertices',vertices,'FaceColor',colorF, 'FaceAlpha', alphaF,'EdgeColor','none');
%     end
%     for i = 1:size(vertOnControl)
%         vertices = [vertOnControlTime(i) bottom; vertOffControlTime(i) bottom; vertOffControlTime(i) top; vertOnControlTime(i) top];
%         faces = [1     2     3     4];
%         patch('Faces',faces,'Vertices',vertices,'FaceColor',colorFControl, 'FaceAlpha', alphaFControl,'EdgeColor','none');
%     end
% 
% end
% plot(tVals, dFFRescaled(:, neurSel)); % 10 of the uniquely identified neurons (Kato 2015)
% %     xlabel('Time (s)');
% axis([0 max(tVals) bottom top]);
% title({'';'Original Data Y'});
% set(gca, 'TickDir', 'out');
% box off;
% hold off;
% 
% subplot(2,1,2),cla;
% hold on;
% bottom = 0;
% top = max(dFFRescaled,[],'all');
% if ifHighlight
%     for i = 1:size(vertOnSuspect)
%         vertices = [vertOnSuspectTime(i) bottom; vertOffSuspectTime(i) bottom; vertOffSuspectTime(i) top; vertOnSuspectTime(i) top];
%         faces = [1     2     3     4];
%         patch('Faces',faces,'Vertices',vertices,'FaceColor',colorF, 'FaceAlpha', alphaF,'EdgeColor','none');
%     end
%     for i = 1:size(vertOnControl)
%         vertices = [vertOnControlTime(i) bottom; vertOffControlTime(i) bottom; vertOffControlTime(i) top; vertOnControlTime(i) top];
%         faces = [1     2     3     4];
%         patch('Faces',faces,'Vertices',vertices,'FaceColor',colorFControl, 'FaceAlpha', alphaFControl,'EdgeColor','none');
%     end
% end
% plot(tVals, dataReconstructedRescaled(:, neurSel));
% xlabel('Time (s)');
% axis([0 max(tVals) bottom top]);
% title({'';'Estimated Y'});
% set(gca, 'TickDir', 'out');
% box off;
% hold off;
% 
% 
% 
% filename03 = sprintf('TEST_NoStimWorm%d_reco_NoAVA%d_v01.fig',thisWorm,ifNoAVA);
% saveas(fig03, filename03);
% 
% % diary CEDiary_Stim1_stdzd_500iters_NoAVA_v01
% 
% %%
% fig02 = figure(2);
% subplot(2,1,1), plot(A_cell{1}.')   % ASC added 7/25
% subplot(2,1,2), plot(B_cell{1}.')   % ASC added 7/25
% filename02 = sprintf('TEST_NoStimWorm%d_xc_NoAVA%d_v01.fig',thisWorm,ifNoAVA);
% saveas(fig02, filename02);
%%
% figure()
% clims = [min(min(cell2mat(F))), max(max(cell2mat(F)))];
% for dynCount = 1:10
%     subplot(2,5,dynCount)
%     imagesc(F{dynCount},clims)
%     colorbar
%     title('Absolute F')
%     
% end
% 
% figure()
% c_values = B_cell{1}.';
% addpath(genpath('./External_Packages/distinguishable_colors/'))
% cmapMaxDistinct = distinguishable_colors(size(c_values,2));
% subplot(2,1,1), plot(A_cell{1}.')   % ASC added 7/25
% subplot(2,1,2)
% for ii=1:size(c_values,2)
%     Y = c_values(:,ii);
%     plot(tVals,Y,'color',cmapMaxDistinct(ii,:),'DisplayName',num2str(ii));
%     
%     hold on;
% end
% legend([],'Orientation','horizontal') 
% hold off;
% 
% %%
% %behavior vs. rSLDS states 
% totalBhvStates = 8;
% x_values = A_cell{1};
% frames = size(x_values,2);
% figure();
% 
% hold on;
% bottom = 1;
% top = totalBhvStates;
% sz = 4 .* ones(frames,1);
% yyaxis left
% scatter(tVals, bhv, sz, 'MarkerFaceAlpha',0.5);
% ylabel('Behavioral state');
% axis([0 max(tVals) bottom top]);
% ylim('padded')
% 
% % yyaxis right
% % plot(tVals, wormdiscretesparse);
% % ylabel('rSLDS label');
% % legend(gca,'off')
% % bottom = 1;
% % top = totalBhvStates;
% % axis([0 max(tVals) bottom top]);
% % ylim('padded')
% set(gca, 'TickDir', 'out');
% 
% xlabel('Time (s)');
% set(gcf,'color','w');
% 
% box off;
% hold off;
%% extract names of channels
% channelNames = dat.WT_Stim(1).IDs;
% channelNamesStrings = strings(1,size(dFF,2));
function channelNamesStrings = getChannelNamesAsStrings(channelNames,numChannels,selectedChannelInds)
channelNamesStrings = strings(1,numChannels);
for channelNum = selectedChannelInds
    channelName = channelNames(channelNum);
    str = '';
%     startEmpty = 0;
    for possibleNeurNameCt = 1:size(channelName{1},2)
        if possibleNeurNameCt == 1
            if isempty(channelName{1}{1})==0
                str = channelName{1}{1};                
            end
        else
            if isempty(channelName{1}{possibleNeurNameCt})==0
                nextString = channelName{1}{possibleNeurNameCt};
                str = append(str,'+',nextString);              
            end
        end
    end
    channelNamesStrings(1,channelNum) = str;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%