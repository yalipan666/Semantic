% copy from Lexical/Analyse_codes
% 20210719 clear functions to make scripts more concise

function [data4ICA, comp] = ICA4rawdata(PPara,hdr,data)
%%% get standard fieldtrip strcture by cutting data into random epochs
trig = PPara.icastart:PPara.cutlen4ica:PPara.icaend; 
epoch = [];
epoch.label = hdr.label; 
for i = 1:length(trig)-1
    trlbegin             = trig(i);
    trlend               = trig(i+1)-1;
    epoch.trial{i}       = data(:,trlbegin:trlend);
    epoch.time{i}        = (trlbegin:trlend)/PPara.SR;
    epoch.trl(i,:)       = [trlbegin trlend PPara.cutlen4ica];
    epoch.trialinfo(i,:) = [trlbegin trlend PPara.cutlen4ica];
end

%%% get bad channels
keepchan = {'MEG'};
for cc = 1:length(PPara.badsens)
    keepchan{cc+1} = ['-MEG' PPara.badsens{cc}];
end

%%% preprocessing
cfg          = [];
cfg.channel  = keepchan;  %% only valid MEG channel
cfg.detrend  = PPara.detrend;     % Remove slow drifts
cfg.bpfilter = PPara.bpfilter;
cfg.bpfreq   = PPara.bpfreq;
data4ICA     = ft_preprocessing(cfg, epoch);
clear epoch

%%% ICA
%%% downsampling --- so so so slow, so skip this 
%%% ICA decomposition
tic
cfg                 = [];
cfg.method          = 'runica';
cfg.runica.maxsteps = 100;
comp                = ft_componentanalysis(cfg,data4ICA);
toc
end

