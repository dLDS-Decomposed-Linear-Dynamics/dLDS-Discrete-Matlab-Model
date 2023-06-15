%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script opens a local pool of Matlab workers (not Distributed       %
% Computing Server) on a desktop computer (any platform) or a cluster     %
% with a torque scheduler for Matlab versions r2013b and r2014a and       %
% torque version 4.2.5.h3.                                                %
%                                                                         %                                                                        %
% Blake Fleischer                                                         %
% Georgia Institute of Technology                                         %
% Partnership for an Advanced Computing Environment (PACE)                %
% http://www.pace.gatech.edu                                              %
% 07/01/2014                                                              %                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [poolsize] = PaceParallelToolBox_r2014a_2(maxoutMyCpu)

%For local jobs, this script opens a pool equal to processors on local system,
%using the default maxoutMyCpu=false opens all but one processor.
%For torque jobs in a cluster environment, reads the environment variables 
%and opens a pool of workers equal to the requested number of processors. 
%If it fails, an optional script is used to read in environment variables 
%from the job an resubmit.

debugClustPool = false; %command switch to run debugging info for cluster testing.
%This should run automatically if there is a failure to open a matlab pool.

if nargin <1 %default is to not occupy all processors on local computer
    maxoutMyCpu = false;
end

%Check current version of matlab: is it r2014?
currVer = version;
if sum(strfind(currVer,'R2014a'))>0
    r2014a = true;
else
    r2014a = false;
end

%% Parallel Processing Initiation (and differentiation)
tic
if r2014a %if matlab r2014a
    maxCores=512; %max allowable cores in r2014a
    currClust = gcp('nocreate');
    if ~isempty(gcp('nocreate'))
        currCores = currClust.NumWorkers;
    else
        currClust = parcluster('local');
        currCores=0;
    end
else
    maxCores=12; %max allowable cores in r2013b,etc
    currCores = matlabpool('size'); %number of cores currently used in matlab pool
end
poolsize = 0; %poolsize starts at zero, but is increased to accomodate multicore

%if running on a local version of matlab (desktops or laptop)
if maxCores > 1 && strcmp(getenv('PBS_JOBID'),'')
    numcores= feature('numCores');
    if ~maxoutMyCpu;
        poolsize = (numcores-1);
        fprintf(['%u cores available, changing to %u core(s) '...
            '(leaving one core) so can still multitask.\n'],numcores, poolsize)
    else
        poolsize = numcores;
        fprintf('%u cores available, changing to the max of %u cores (no multitask).\n', numcores, poolsize)
    end
    
    if currCores~=poolsize && currCores > 0; %close matlabpool if it has been previously opened to a different number of cores
        if r2014a
            delete(gcp('nocreate'));
            currClust = parcluster('local'); %have to respecify local pool
        else
            matlabpool close
        end
    end
    
    if currCores~=poolsize; %If different number of current cores than desired,
        % close current pool and reopen one of the right number of workers
        if r2014a
            set(currClust,'NumWorkers',poolsize)
            parpool('local',poolsize);
            currPool = gcp('nocreate');
            currPool.IdleTimeout=inf;
        else
            matlabpool(poolsize)
        end
    end
    
    [~,sys] = memory; %Determine total RAM available
    allocatedRam = sys.PhysicalMemory.Available;
    
else maxCores > 1 && ~strcmp(getenv('PBS_JOBID'),''); %If currently running as a torque job (on the cluster)
	
    if debugClustPool
        %Cluster Debugging Info
        [~, jobConds] = system(['qstat -f ',getenv('PBS_JOBID')]);
        fprintf('Current running environment: \n%s',jobConds)
        fprintf('\ndone.\n\n')
        
        fprintf('\nRunning top to check other processes:\n');
        !top -bn 1 | head -n 70
        fprintf('done.\n\n')
        
        fprintf(['\n\nRunning checks on the current matlab status on the cluster. ',...
            'These checks are from:\n',...
            'http://www.mathworks.com/matlabcentral/answers/92124-why-am-i-unable-to-use-matlabpool-or-parpool-with-the-local-scheduler-or-validate-my-local-configura',...
            '\n'])
        
        fprintf('Checking the parallel computing toolbox...\n')
        license checkout Distrib_Computing_Toolbox
        fprintf('done.\n\n')
        
        fprintf('Checking version of PCT...\n')
        versionM = ver;
        versionM(39)
        fprintf('Checking version of MATLAB...\n')
        version
        fprintf('done.\n\n')
        
        fprintf('Disabling LocalUseMpiexec...\n')
        distcomp.feature( 'LocalUseMpiexec', false )
        fprintf('done.\n\n')
        
        localSchedFolder=strsplit(prefdir,'R201');
        localSchedFolder = [localSchedFolder{1}, 'local_scheduler_data'];
        fprintf(['One could delete the locla_scheduler data, ',...
            'it is located in %s.\n\n'],localSchedFolder)
        
        fprintf('Checking hostname and ping...\n')
        !echo $HOSTNAME
        !ping -c 8 $HOSTNAME
        fprintf('done.\n\n')
        
        fprintf('Checking for a startup.m file ...\n')
        which startup.m
        fprintf('done.\n\n')
        
        fprintf('Checking ulimit settings, this has been problematic in the past:\n')
        !ulimit -a
        fprintf('done.\n\n')
        
        fprintf('Checking Java memory settings:\n')
        %% Report on environment - this is mostly for error checking on the Force-6 cluster
        fprintf('java.lang.Runtime.getRuntime.maxMemory = %u\n',java.lang.Runtime.getRuntime.maxMemory)
        fprintf('java.lang.Runtime.getRuntime.totalMemory = %u\n',java.lang.Runtime.getRuntime.totalMemory)
        fprintf('java.lang.Runtime.getRuntime.freeMemory = %u\n',java.lang.Runtime.getRuntime.freeMemory)
        fprintf('done.\n')
        
        fprintf('End of cluster checks.\n\n')
    end
	
    [~, jobConds] = system(['qstat -f ',getenv('PBS_JOBID')]);
    %fprintf('Current running environment: \n%s',jobConds)
    jobCondsArr = strsplit(jobConds)';
    fprintf('\n\n')
    %Before creating pool, check the number of allocated processors first
    allocatedProc=sscanf(jobCondsArr{find(strcmp('Resource_List.nodes',jobCondsArr))+2},'%*6c%u',[1, Inf]);
	
	if isempty(find(strcmp('Resource_List.mem',jobCondsArr), 1))
		fprintf(['Please specify ram when making job submissions - use ''msub -l mem=XXgb.'' ',...
		'The parallel pool script \nneeds to know if you have reserved enough ram for the ',...
		'job. Otherwise the job may unexpectedly \nquit.'])
        return
    end
	
    allocatedRam =sscanf(jobCondsArr{find(strcmp('Resource_List.mem',jobCondsArr))+2},'%u%*2c',[1, Inf]);
    if sum(strcmp('interactive',jobCondsArr))>0  %Check if job is interactive
        fprintf('Job is interactive.\n')
        interActiv = strcmp('True',jobCondsArr{find(strcmp('interactive',jobCondsArr))+2});
    else
        fprintf('Job is not interactive.\n')
        interActiv=false;
    end
    fprintf('Looks like %u/%u workers allocated for %u gb ram.\n',...
        currCores,allocatedProc,allocatedRam)
    if allocatedProc <= maxCores;
        poolsize = allocatedProc;
    else
        poolsize = maxCores;
        fprintf(['\n\nYou have reserved %u cores but can only use a maximum of ',...
            '%u cores for this verson of MatLab (%s). \n%u cores are going to ',...
            'waste. Please reserve a maximum of ',...
            '%u cores (shorter queue times for you, others can use the \ncores ',...
            'you won''t be using).\n\n'],allocatedProc,maxCores,version,allocatedProc-maxCores,...
            maxCores)
    end
    
    if currCores ~=poolsize %If the current pool is different than the desired size (i.e. if
        %the pool isn't already open).
        
        %Change location of job storage location so multiple jobs can run
        %simultaneously. If the folders don't get deleted, just remove the
        %jobs folders for completed or cancelled jobs.
        fprintf('\nInitializing local scheduler data.\n')
        currClust = parcluster('local');
        
        fprintf('Matlabpool job current data location is %s\n', currClust.JobStorageLocation);
        local_clust_data = fullfile( currClust.JobStorageLocation,...
            ['jobID_', strtok(getenv('PBS_JOBID'),'.')]);
        mkdir(local_clust_data);
        %[~,mess,messid] = mkdir(local_clust_data);
        %if messid=='MATLAB:MKDIR:DirectoryExists'
        %    fprintf('\n\n%s\n\n',mess)
        %end
        currClust.JobStorageLocation = local_clust_data;
        fprintf('New data loc is: %s\n', currClust.JobStorageLocation);        
        %because the datalocation changes if run more than once on the same node,
        %find the first instance within that pid folder
        %jobIds = strfind(local_clust_data,'jobID');
        %slashes = strfind(local_clust_data,'/');
        %savDirStrLen = slashes(find(slashes>jobIds(1),1,'first'));
        
        %Open parallel processing pool
        if ~interActiv %1 user probably won't open multiple interactive jobs on the
            %same node at the same time, but it is a problem on batch jobs
            rng(str2num(strtok(getenv('PBS_JOBID'),'.'))); %Uniquely seed the random
            %number generator
            randWaitTime = 1+60*rand();%Multiple parallel pools starting on the same
            %machine and time can cause problems...
            fprintf('Waiting %2.3f seconds to prevent parpool opening errors...',...
                randWaitTime)
            pause(randWaitTime);
            fprintf('done.\n')
        end
        set(currClust,'NumWorkers',poolsize);
        try %If parpool fails (intermittely fails on cluster, not sure why)
            parpool(currClust,poolsize);
        catch err
            %Get environment variables to resubmit job if this one fails
            fprintf('\n\nOriginal error Information:\n');
            disp(getReport(err,'extended'));
            %Run some other commands to help diagnose the problem...
            
            %Cluster Debugging Info
            [~, jobConds] = system(['qstat -f ',getenv('PBS_JOBID')]);
            fprintf('Current running environment: \n%s',jobConds)
            fprintf('\ndone.\n\n')
            
            fprintf('\nRunning top to check other processes:\n');
            !top -bn 1 | head -n 70
            fprintf('done.\n\n')
            
            fprintf(['\n\nRunning checks on the current matlab status on the cluster. ',...
                'These checks are from:\n',...
                'http://www.mathworks.com/matlabcentral/answers/92124-why-am-i-unable-to-use-matlabpool-or-parpool-with-the-local-scheduler-or-validate-my-local-configura',...
                '\n'])
            
            fprintf('Checking the parallel computing toolbox...\n')
            license checkout Distrib_Computing_Toolbox
            fprintf('done.\n\n')
            
            fprintf('Checking version of PCT...\n')
            versionM = ver;
            versionM(39)
            fprintf('Checking version of MATLAB...\n')
            version
            fprintf('done.\n\n')
            
            fprintf('Disabling LocalUseMpiexec...\n')
            distcomp.feature( 'LocalUseMpiexec', false )
            fprintf('done.\n\n')
            
            localSchedFolder=strsplit(prefdir,'R201');
            localSchedFolder = [localSchedFolder{1}, 'local_scheduler_data'];
            fprintf(['One could delete the locla_scheduler data, ',...
                'it is located in %s.\n\n'],localSchedFolder)
            
            fprintf('Checking hostname and ping...\n')
            !echo $HOSTNAME
            !ping -c 8 $HOSTNAME
            fprintf('done.\n\n')
            
            fprintf('Checking for a startup.m file ...\n')
            which startup.m
            fprintf('done.\n\n')
            
            fprintf('Checking ulimit settings, this has been problematic in the past:\n')
            !ulimit -a
            fprintf('done.\n\n')
            
            fprintf('Checking Java memory settings:\n')
            %% Report on environment - this is mostly for error checking on the Force-6 cluster
            fprintf('java.lang.Runtime.getRuntime.maxMemory = %u\n',java.lang.Runtime.getRuntime.maxMemory)
            fprintf('java.lang.Runtime.getRuntime.totalMemory = %u\n',java.lang.Runtime.getRuntime.totalMemory)
            fprintf('java.lang.Runtime.getRuntime.freeMemory = %u\n',java.lang.Runtime.getRuntime.freeMemory)
            fprintf('done.\n')
            
            fprintf('End of cluster checks.\n\n')
                    
        end %end of if parpool fails debugging info
        
        %Set pool properties
        currPool = gcp('nocreate');
        currPool.IdleTimeout=inf;
    else
        fprintf('...remote pool already opened.\n')
    end
end

if r2014a %Double check that pool has been created.
    checkPool = gcp('nocreate');
    if isempty(checkPool)
        poolsize = 0;
    else
        poolsize = checkPool.NumWorkers;
        idleTimeout = checkPool.IdleTimeout;
    end
else
    poolsize = matlabpool('size'); %number of cores currently used in matlab pool
end

if currCores~=poolsize %only show message if changed size of matlab pool
    %Time Convert anonymous function
    timeConvert=@(secIn) (sprintf('%u:%u:%u:%2.4f (D:H:M:S)', [floor(secIn/(24*60*60)),...
        floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)),...
        floor((secIn-(24*60*60)*floor(secIn/(24*60*60))-(60*60)*floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)))/60),...
        secIn-(24*60*60*floor(secIn/(24*60*60)))-(60*60*floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)))-(60*floor((secIn-(24*60*60)*floor(secIn/(24*60*60))-(60*60)*floor((secIn-24*60*60*floor(secIn/(24*60*60)))/(60*60)))/60))]));
    
    fprintf('Matlabpool opened with %u workers in: %s', poolsize, timeConvert(toc));
    if r2014a
        if idleTimeout==inf;
            fprintf('; pool will not time out.\n')
        else
            fprintf('; pool IdleTimeout is: %s minutes.\n',idleTimeout);
        end
    else
        fprintf('\n')
    end
    
    if allocatedRam < (0.330 * poolsize) 
        fprintf(['\n\n!!Warning - ram allocation of %u gb is probably too ',...
            'low.\nConsider running 0.330*ppn gb of ram plus the size ',...
            'of any variables * ppn (%u gb + %u*vars_size).\nJob may '...
            'terminate unexpectedly as a result.'],...
            allocatedRam,ceil(0.330*poolsize),poolsize);
    end
end

% %To close pool:
% if r2014a
%     delete(gcp('nocreate'));
% else
%     matlabpool close;
% end
