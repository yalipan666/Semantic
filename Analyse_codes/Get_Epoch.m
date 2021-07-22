% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise

function [epoch_sumrj,ValidTrlPerct] = Get_Epoch(hdr,data,PPara,trig_col)
%%% get raw epoch
epoch = []; %time zero is 'fixation_on_MEG'
if isempty(trig_col)
    error('This trigger does not exist!')
end
event_all = PPara.event_all;
for j = 1:size(event_all,1)
    trlbegin             = event_all(j,trig_col) + PPara.pretrig;
    trlend               = event_all(j,trig_col) + PPara.posttrig;
    epoch.trial{j}       = data(:,trlbegin:(trlend-1));
    epoch.time{j}        = (PPara.pretrig:(PPara.posttrig-1))/PPara.SR;
    epoch.trialinfo(j,:) = event_all(j,:);
    epoch.trl(j,:)       = [trlbegin trlend PPara.pretrig];
end

%% get the trial strcture from standard fieldtrip strcture
epoch.label = hdr.label;
cfg = [];
cfg.trl =  [epoch.trl epoch.trialinfo];
epochdata = ft_redefinetrial(cfg,epoch);
%%% trials in epoch_cond is wrong, so replace them with epoch.trial
epochdata.trial = epoch.trial;
epochdata.time = epoch.time;

%% preprocessing
%%% get bad channels
keepchan = {'MEG','EOG','MISC004','MISC005'};
for cc = 1:length(PPara.badsens)
    keepchan{cc+1} = ['-MEG' PPara.badsens{cc}];
end
cfg = [];
cfg.channel = keepchan;
epochdata = ft_preprocessing(cfg, epochdata);

%%%  remove badtrials and sensors
%%actifacts rejection in planar gradiometers
cfg         = [];
cfg.channel = 'MEGGRAD';
data_planar = ft_selectdata(cfg,epochdata);

%%% remove bad sensors and trials using ft_rejectvisual
cfg=[];
cfg.method  = 'summary';
tempdata_rjv_planar = ft_rejectvisual(cfg, data_planar);

%%% make a list of the planar gradiometers to keep
oldtrig = tempdata_rjv_planar.trialinfo(:,trig_col);
newtrig = data_planar.trialinfo(:,trig_col);
trl_rej = find(ismember(newtrig,setdiff(newtrig,oldtrig))); % same trigger time
trl_keep = (1:size(epochdata.trial,2));
trl_keep(trl_rej) = [];
chan_rej = setdiff(data_planar.label,tempdata_rjv_planar.label);

%%%%%%% get the valid trail number
ValidTrlPerct = (length(trl_keep)/(length(trl_rej)+length(trl_keep)));
if ValidTrlPerct < 0.75
    warning('**** Too many noisy trials!!! ****');
end
chan_keep = epochdata.label;
chan_keep((ismember(epochdata.label(:,1),chan_rej(:,1))))=[];
cfg = [];
cfg.trials = trl_keep;
cfg.channel = chan_keep;
epoch_sumrj = ft_selectdata(cfg,epochdata);
end
