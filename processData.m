classdef processData
    % this class uses plotEMGs.mlapp figure class for plotting 
    
    properties
        srate = 1000
        emg_idx = (1:10) + 7        % col idx in raw data
        fwin = 250;                 % ms (because sampling at 1000 Hz).
        offset_idx = 50:250         % indices to get offset activation 
        greyrgb = [169,169,169]/255 
        
        EMG_Names = {'SCM R', 'SCM L', 'Trap R', 'Trap L', ...
                    'SSCap R','SSCap L','SpCap R', 'SpCap L',...
                    'LS R','LS L'}
        fig                         % handle to plotEMGs object
        EMGAx                       % indices of EMG Axes in fig
        WrenchAx                    % indices of Wrench Axes in fig
        
        data
        mEMG                      % struct .maxEMG, .
    end
    
    
    methods
        
        % Constructor
        function obj = processData(matData)
            
            if nargin < 1
                matData = obj.uploadData();
            end
            
            % try to get nEMG channels from metadata
            try
                nEMG = length(matData(1).metadata.daq.EMGChannels);
                obj.emg_idx = (1:nEMG) + 7;            
            catch
                warning('Could not get nEMG channels, using default 10')
            end
            
            try 
                obj.srate = matData(1).metadata.daq.sRate;
            catch 
                warning('Could not get srate from metadata, using default 1k Hz')
            end
            
            % assign matData to data object
            obj.data = matData;
            
            % set valid_time to full time window initially
            for ii = 1 : length(obj.data)
                traw = obj.data(ii).rawdata(1,:);
                traw = traw - traw(1);
                obj.data(ii).valid_time = [0 traw(end)];                
            end
            
        end
        
        % Get Data from trials in folder
        function matData = uploadData(~, myFolder)
           % if directory not provided
            if nargin < 2
               % load data from all .mat files in selected folder into matData array of structs
                myFolder = uigetdir();
           end   
           
            % note, to read in order, must be numbered with zeros padding left side
            filePattern = fullfile(myFolder, '*.mat');
            matFiles = dir(filePattern);
            for k = 1:length(matFiles)
              baseFileName = matFiles(k).name;
              fullFileName = fullfile(myFolder, baseFileName);
              fprintf(1, 'Now reading %s\n', fullFileName);
              matData(k) = load(fullFileName);
            end            
        end
        
        function obj = keepTheseTrials(obj, valid_idx)
            % these indices are the index of matData, not necessarily trial
            % number if something went wrong in data collection. Should
            % check that first.. 
            
            ntrials = length(obj.data);
            nonvalid_trials = setdiff(1:ntrials,valid_idx);

            % Remove nonvalid trials
            obj.data(nonvalid_trials) = [];            
        end
        
        function obj = filterEMGWrapper(obj,offset)
            % runs through raw emg data and filters it
            % offset is optional
            % adds fields .remg, .femg, and .offset to obj.data struct
            
            for ii = 1 : length(obj.data)
                                
                % get raw emg data 
                emgraw = obj.data(ii).rawdata(obj.emg_idx,:);
                
                if nargin < 2 % offset not provided
                    [obj.data(ii).femg, obj.data(ii).remg, obj.data(ii).offset] = ...
                        obj.filterEMG(emgraw,obj.fwin);
                else % offset provided
                    [obj.data(ii).femg, obj.data(ii).remg, obj.data(ii).offset] = ...
                        obj.filterEMG(emgraw,obj.fwin, offset);
                end
                    
            end
            
        end
        
        function obj = normalizeEMGWrapper(obj, maxEMG)
            if nargin < 2
                if ~isempty(obj.mEMG)
                    maxEMG = obj.mEMG.maxEMG;
                else
                    error('maxEMG not provided')
                end
            end
            
            % runs through emg and normalizes it
            % adds fields .nremg, .nfemg to obj.data struct
            
            for ii = 1 : length(obj.data)
                obj.data(ii).nremg = obj.normalizeEMG(obj.data(ii).remg, maxEMG);
                obj.data(ii).nfemg = obj.normalizeEMG(obj.data(ii).femg, maxEMG);
            end            
            
        end
        
        function [femg, remg, offset] = filterEMG(obj, emgraw,fil_win,offset)
            % filterEMG: emgraw: (#ch x nsamples) matrix of raw emg data
            
            % - detrend 
            % - rectify (abs())
            % - subtract offset (try avg rectified signal from 50-150 ms)
            % - filter with moving average filter (250 ms window)
            
            data_dim = size(emgraw);
            % make sure each row is a channel (#ch x nsamples)
            if data_dim(1) > data_dim(2)
                emgraw = emgraw';
            end

            % trend and rectify data
            remg = abs(detrend(emgraw));

            % subtract offset (avg activation from beginning of trial)
            if nargin < 4
                offset = mean(remg(:,obj.offset_idx),2);
            else 
                if isrow(offset) % make sure its a column vec
                    offset = offset';
                end
                if length(offset) ~= size(remg,1)
                    error('offset has wrong dimensions')
                end
            end
            remg = remg - offset;

            % filter with window of fwin
            femg = winfilter(remg,fil_win);
            
            % if transposed at beginning, change result back 
            if data_dim(1) > data_dim(2)
                remg = remg';
                femg = femg';
                offset = offset';
            end
            
        end
        
        function obj = getMaxEMG(obj)
            % get max filtered EMG value for normalization purposes
            
            % get max EMG across all trials for each channel
            for ii = 1 : length(obj.data) % each trial
                vtime = obj.data(ii).valid_time;
                tvec = obj.data(ii).rawdata(1,:);
                tvec = tvec - tvec(1);
                start_idx = obj.getIdxClosestTime(tvec,vtime(1));
                end_idx = obj.getIdxClosestTime(tvec,vtime(2));
                [maxValue(:,ii),maxidx(:,ii)] = max(obj.data(ii).femg(:,start_idx:end_idx),[],2);
            end
            [maxEMG, maxEMGTrial] = max(maxValue,[],2);
            
            % get the timestamp index of max EMG 
            %for ii = 1 : size(maxidx,1)
            %    maxEMGTidx(ii) = maxidx(ii,maxEMGTrial(ii));
            %end 
            
            % add to mEMG property
            obj.mEMG.maxEMG = maxEMG;
            obj.mEMG.maxEMGTrial = maxEMGTrial;
            %obj.mEMG.maxEMGTidx = maxEMGTidx; % see if I can get rid of this..? 
            
        end
        
        
        
        function nemg = normalizeEMG(~,emgdata, maxEMG)
            % check dimensions
            if size(emgdata,2) ~= size(maxEMG,1)
                if size(emgdata,2) == size(maxEMG,2)
                    maxEMG = maxEMG';
                elseif size(emgdata,1) == size(maxEMG,1)
                    emgdata = emgdata';
                else
                    error('dimensions of emgdata and maxEMG do not match')
                end
            end
            
            % normalize by EMG maxes
            nemg = emgdata' ./ maxEMG;                
                
        end
        
        function [maxW,maxT] = getMaxWrench(obj, trial)
            vtime = obj.data(trial).valid_time;
            
            w = obj.data(trial).wrenchdata(8,:)' .* obj.data(trial).wrenchdata(9,:)';
            tw = obj.data(trial).wrenchdata(1,:)';
            tw = tw - tw(1);
            
            start_idx = obj.getIdxClosestTime(tw,vtime(1));
            end_idx = obj.getIdxClosestTime(tw,vtime(2));
            
            % set wrench values outside of time range to 0
            w(1:start_idx) = 0; 
            w(end_idx:end) = 0;
            
            % find max wrench and time 
            [maxW, max_idx] = max(w);
            maxT = tw(max_idx);            
            
        end
        
        function [emg, emgT] = getEMGatTime(obj, trial, time)
            t = obj.data(trial).rawdata(1,:);
            t = t-t(1);       
            [~,idx] = min(abs(t-time));
            
            emgT = t(idx);
            % normalized and filtered emg value
            emg = obj.data(trial).nfemg(:,idx);            
        end
        
        function obj = computeEMGatMaxW(obj)
            
            ntrials = length(obj.data);
            time_offset = 0.05; % 50 ms 
            
            for ii = 1:ntrials
                
                [obj.data(ii).maxW, maxT] = obj.getMaxWrench(ii);
                maxT = maxT - time_offset;
                [obj.data(ii).emg, ~] = obj.getEMGatTime(ii, maxT);
 
            end
            
        end
        
        function idx = getIdxClosestTime(~,tvec,time)
            [~,idx]=min(abs(tvec-time));
        end
        
        function obj = cutDataOutsideTime(obj,trial,start_time,end_time)
            % start_time and end_time are in sec
            
            % cut raw data
            t = obj.data(trial).rawdata(1,:);
            t = t - t(1);
            start_idx = obj.getIdxClosestTime(t,start_time);
            end_idx = obj.getIdxClosestTime(t,end_time);
            
            obj.data(trial).rawdata = obj.data(trial).rawdata(:,start_idx:end_idx);
            
            % cut filtered data if it is present
            if isfield(obj.data(trial),'remg')
                obj.data(trial).remg = obj.data(trial).remg(:,start_idx:end_idx);
            end
            
            if isfield(obj.data(trial),'femg')
                obj.data(trial).femg = obj.data(trial).femg(:,start_idx:end_idx);
            end
            
            if isfield(obj.data(trial),'nremg')
                obj.data(trial).nremg = obj.data(trial).nremg(:,start_idx:end_idx);
            end
            
            if isfield(obj.data(trial),'nfemg')
                obj.data(trial).nfemg = obj.data(trial).nfemg(:,start_idx:end_idx);
            end
            
            % cut wrench data
            tw = obj.data(trial).wrenchdata(1,:);
            tw = tw - tw(1);
            start_idx = obj.getIdxClosestTime(tw,start_time);
            end_idx = obj.getIdxClosestTime(tw,end_time);
            
            obj.data(trial).wrenchdata = obj.data(trial).wrenchdata(:,start_idx:end_idx);

        end
        
        function obj = changeValidTimeRange(obj,trial,startT,endT)
            obj.data(trial).valid_time = [startT endT];
            % need to rerun find max and all that.. 
        end
        
        function S = exportData(obj)
            
            if isfield(obj.data(1),'metadata')
                if isfield(obj.data(1).metadata,'subject')
                    if isfield(obj.data(1).metadata.subject,'subjID')
                        S.subject = obj.data(1).metadata.subject.subjID;
                    end
                    if isfield(obj.data(1).metadata.subject,'Date')                    
                        S.date = obj.data(1).metadata.subject.Date;
                    end
                end
                if isfield(obj.data(1).metadata,'trialnum')
                    for ii = 1 : length(obj.data)
                        S.trialnum(ii) = obj.data(ii).metadata.trialnum;
                    end
                end
                if isfield(obj.data(1).metadata,'target')
                    if isfield(obj.data(1).metadata.target,'dir_name')            
                        for ii = 1 : length(obj.data)
                            S.dir_name{ii} = obj.data(ii).metadata.target.dir_name;
                        end
                    end
                    if isfield(obj.data(1).metadata.target,'dir_index')
                        for ii = 1 : length(obj.data)
                            S.dir_index(ii) = obj.data(ii).metadata.target.dir_index;
                        end
                    end
                end
            end
            
            if isfield(obj.data(1),'maxW')
                S.maxW = [obj.data.maxW];
            end
            if isfield(obj.data(1),'emg')
                S.emg = [obj.data.emg];
            end
            
        end
        
        function plotOffsets(obj, offset)
            % check offsets by plotting them
            os_dim = size(offset);
            if os_dim(1) < os_dim(2)
                offset = offset';
            end
            
            figure('Name','Offset Values')
            plot(offset)
            legend(obj.EMG_Names)
            xlabel('Trial')
            ylabel('Offset Magnitude')            
            
        end
        
        function obj = setupEMGPlot(obj)
           % 2 column plot with all EMG (no wrench)
            nAx = length(obj.emg_idx);
            ncol = 2;
            ch_names = obj.EMG_Names;
            
            fig_obj = plotEMGs(nAx,ncol,ch_names);
            obj.EMGAx = 1:nAx;
            obj.fig = fig_obj;
        end
        
        function obj = setupWrenchEMGPlot(obj)
            % 1 column plot with wrench at top
            nAx = length(obj.emg_idx) + 1;
            ncol = 1;
            ch_names = [{'wrench'},obj.EMG_Names];
            
            fig_obj = plotEMGs(nAx,ncol,ch_names);
            obj.EMGAx = 2:nAx;
            obj.WrenchAx = 1;
            obj.fig = fig_obj;
            
        end
        
        function plotEMGTrial(obj, trial, normbool)
            % data: struct with .remg, .femg, .rawdata or .timestamps,
            % .wrenchdata (if plotting wrench)
            if nargin < 3
                normbool = 1; %default 
            end
            % matData should only be for that trial (ie matData(trial))
            
            f = obj.fig;
            d = obj.data(trial);
            % set panel title
            if isfield(d,'metadata')
                obj.setPanelTitle(d.metadata);
            end
            
            % get timestamps
            if isfield(d,'rawdata')
                t = d.rawdata(1,:)';
            elseif isfield(d,'timestamp')
                t = d.timestamp(1,:)';
            end
                            
            t = t - t(1);
            
            if normbool && isfield(d,'nremg')
                remg = d.nremg';
                femg = d.nfemg';
            else
                warning('Displaying non-normalized EMG')
                remg = d.remg';
                femg = d.femg';
            end
                       
            f.plotData(f.Ax(obj.EMGAx), t, remg, 1, 'Color', obj.greyrgb);
            f.plotData(f.Ax(obj.EMGAx), t, femg, 0, 'Color', 'k');
            
            % plot the wrench if there is a wrench axis
            if ~isempty(obj.WrenchAx)
                w = d.wrenchdata(8,:)' .* d.wrenchdata(9,:)';
                tw = d.wrenchdata(1,:)';
                tw = tw - tw(1);                
                % find max wrench time 
                [~,maxT] = obj.getMaxWrench(trial);
                
                % plot wrench 
                f.plotData(f.Ax(obj.WrenchAx), tw, w, 1, 'Color', 'k');
                
                % plot vertical line at max wrench location 
                f.plotVertLine(f.Ax, maxT,'r');
            end
            
        end
        
        
        function setPanelTitle(obj, metadata)           
            % input:
            %           metadata: struct - parses to find relevant info
            %           gui_obj: plotEMGs handle for figure
            
            trialnum = metadata.trialnum;
            if isfield(metadata.target,'dir_index')
                dir_index = metadata.target.dir_index;
            else
                dir_index = 0;
            end
            if isfield(metadata.target,'dir_name')
                dir_name = metadata.target.dir_name;
            else
                dir_name = '';
            end
            title_str = sprintf('Trial %i, Direction %i %s',trialnum,dir_index,dir_name);
            
            % set panel title 
            obj.fig.setPanelTitle(title_str);
            
        end
        
      
        
    end
       
    
    
    
end