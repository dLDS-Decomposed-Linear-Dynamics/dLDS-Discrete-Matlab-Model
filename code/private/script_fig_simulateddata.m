%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
myfilename = 'saveMSNWorkspace_3000i_10i2_Finitsvddecay_Kp5_Lval_p1_Lb_p1_tol_1em3_D10.mat';
myWorkspace = importdata(myfilename);

%% Precomputation of useful things


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot dynamics and states

figure();

c_values = myWorkspace.B_cell{1}.';
x_values = myWorkspace.A_cell{1};
x_values_timebyx  = x_values.';
thisDFF = myWorkspace.dFF;
dataReconstructed = (myWorkspace.Phi * x_values).'; % note: not rescaled

frames = size(x_values,2);

disp('tVals = frames for simulated data');
tVals = 1:frames; 

hold on;
bottom = min(c_values,[],'all');
top = max(c_values,[],'all');
plot(tVals, c_values); 
xlabel('Time points');
axis([0 max(tVals) bottom top]);
title('Estimated coefficient c values');
set(gca, 'TickDir', 'out');
box off;
hold off;


%% 
figure();

disp('tVals = frames for simulated data');
tVals = 1:frames; 

hold on;
bottom = min(x_values,[],'all');
top = max(x_values,[],'all');
plot(tVals, x_values); 
xlabel('Timepoints');
axis([0 max(tVals) bottom top]);
title('Estimated latent states x');
set(gca, 'TickDir', 'out');
box off;
hold off;

% note: not rescaled
flatDFF = reshape(thisDFF, numel(thisDFF),[]);
flatReconstructed = reshape(dataReconstructed,numel(dataReconstructed),[]);
rval = corrcoef(flatDFF, flatReconstructed);
varExpl = (rval(1,2))^2;
disp(varExpl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
figure();
imagesc(x_values); 
colorbar;
title('Estimated x');

figure();
imagesc(myWorkspace.latentStatesX);
colorbar;
title('Latent X');

figure;
histogram(x_values);
title('Estimated X');

figure;
histogram(myWorkspace.latentStatesX);
title('Latent X');

%% sparsity at each time point
latentStatesXSparsityAtEachTP = zeros(size(x_values,2),1);
xEstSparsityAtEachTP = zeros(size(x_values,2),1);



for i = 1:size(x_values,2)
    latentStatesXSparsityAtEachTP(i) = sum(abs(myWorkspace.latentStatesX(:,i)) < 1e-1)/size(x_values,1);
    xEstSparsityAtEachTP(i) = sum(abs(x_values(:,i)) < 1e-1)/size(x_values,1);
end

figure();
imagesc(xEstSparsityAtEachTP);
colorbar;
title('Sparsity: Percent Small Values at Each Time Point, Estimated X');

figure();
imagesc(latentStatesXSparsityAtEachTP);
colorbar;
title('Sparsity: Percent Small Values at Each Time Point, Latent X');

%% matrix D (transition between latents and data)

figure();
imagesc(myWorkspace.simulatedDmatrix);
colorbar;
title('Initialized Matrix D between states X and input data');


figure();
imagesc(myWorkspace.Phi);
colorbar;
title('Learned Matrix D between states X and input data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




