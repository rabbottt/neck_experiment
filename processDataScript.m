%% call constructor for processEMG class
obj = processData();


%% (OPTIONAL) Remove trials that will not be used
% note: check that files names loaded match trial number


%valid_trials = [1:19 24:26 29 32:34 37:39 46:47 49:57]; % for S001
%obj = obj.keepTheseTrials(valid_trials);

%% Filter EMG data 
% - detrend 
% - rectify (abs())
% - subtract offset (try avg rectified signal from 50-150 ms)
% - filter with moving average filter (250 ms window)

%offset = [0.0020    0.0028    0.0020    0.0025    0.0020    0.0029    0.0031    0.0051    0.0029    0.0022]'; % S001
offset = zeros(10,1);
obj = obj.filterEMGWrapper(offset);
%obj = obj.filterEMGWrapper();
obj = obj.getMaxEMG();
obj = obj.normalizeEMGWrapper();

%%

obj = obj.computeEMGatMaxW();
%%
data = obj.exportData();

%% plot to GUI 

%obj = obj.setupEMGPlot();
obj = obj.setupWrenchEMGPlot();

%%

trial = 74;

obj.plotEMGTrial(trial);

% display stuff
[maxW,maxT] = obj.getMaxWrench(trial);
[emg,~] = obj.getEMGatTime(trial, maxT);
trialnum = obj.data(trial).metadata.trialnum;
fprintf('Trial #%i: \n',trialnum)
fprintf('Max Wrench: %f N \n\n',maxW)
disp(table(obj.EMG_Names(1:2:10)',emg(1:2:10),obj.EMG_Names(2:2:10)',emg(2:2:10)))


%%
% bad looking trials: 2, 16, 21, 26, 28, 31, 32, 34 (all sp cap) 

% checking trials v wrench 
% trial 2, cut start_time = 8.18; end_time = 10; 
% trial 16, start_time = 7.5; end_time = 10; check max wrench.. it is at very end.
% trial 20, start_time = 12; end_time = 20; 
% trial 21, cut: start_time = 12; end_time = 17;
% trial 24 - max is at very end, 
% 26 ok at time of maxW; start_time = 10; end_time = 16;
% 28 ok start_time = 10.6; end_time = 16;
% 29- SpCap L value is negative.. 
% 31 start_time = 4.5; end_time = 13;
% 32 start_time = 3.3; end_time = 10; ok at time of maxW but super low activations .. 
% 34 start_time = 8.5; end_time = 20;
% 39 negative SpCap L
% 40 start_time = 4; end_time - 20;
start_time = 4;
end_time = 20;
trial = 40;
obj = obj.cutDataOutsideTime(trial,start_time,end_time);


%%
obj.plotOffsets([obj.data.offset]')
