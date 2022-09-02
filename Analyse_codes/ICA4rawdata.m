% copy from Lexical/Analyse_codes
% 20210719 clean functions to make scripts more concise

function [data4ICA, comp, Trig] = ICA4rawdata(PPara,hdr,data)
%%% get standard fieldtrip strcture
epoch = [];
epoch.label = hdr.label; %%% the full labels--340
for i = 1:length(PPara.fixon)
    trlbegin             = PPara.fixon(i);
    trlend               = PPara.sentoff(i);
    epoch.trial{i}       = data(:,trlbegin:trlend);
    epoch.time{i}        = (trlbegin:trlend)/PPara.SR;
    epoch.trl(i,:)       = [trlbegin trlend PPara.pretrig];
    epoch.trialinfo(i,:) = [trlbegin trlend PPara.pretrig];
end
Trig.fixon = PPara.fixon;
Trig.sentoff = PPara.sentoff;

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
downsampling --- so so so slow, so skip this 
%%% ICA decomposition

tic
cfg            = [];
cfg.resamplefs = 200;
cfg.detrend    = 'no';
data4ICA_training           = ft_resampledata(cfg, data4ICA_training);
toc


tic
cfg                 = [];
cfg.method          = 'runica';
cfg.runica.maxsteps = 100;
comp_train          = ft_componentanalysis(cfg,data4ICA_training);
toc


















end

