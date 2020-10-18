%% call constructor for processEMG class
obj = processData();


%% (OPTIONAL) Adjust Valid Time and Remove Invalid Trials
% note: check that files names loaded match trial number

% Code for each subject in separate .m files: 
%   example: S005_DataAnalysis.m
% To specify new time range: 
%   obj = obj.changeValidTimeRange(trial,start_time,end_time);

%% Filter EMG data 
% - detrend 
% - rectify (abs())
% - subtract offset (default for now, no offset)
% - filter with moving average filter (250 ms window)

obj = obj.filterEMGWrapper(zeros(10,1));

% get maxes
obj = obj.getMaxEMG();
obj = obj.normalizeEMGWrapper();
obj = obj.computeEMGatMaxW();

% export summary to data struct
data = obj.exportData();

%% plot to GUI 

%obj = obj.setupEMGPlot();
obj = obj.setupWrenchEMGPlot();

%% Select trial to visualize

trial = 3;

obj.plotEMGTrial(trial);

% display stuff
[maxW,maxT] = obj.getMaxWrench(trial);
[emg,~] = obj.getEMGatTime(trial, maxT);
trialnum = obj.data(trial).metadata.trialnum;
fprintf('Trial #%i: \n',trialnum)
fprintf('Max Wrench: %f N \n\n',maxW)
disp(table(obj.EMG_Names(1:2:10)',emg(1:2:10),obj.EMG_Names(2:2:10)',emg(2:2:10)))

