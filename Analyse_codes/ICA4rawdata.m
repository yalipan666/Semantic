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

% % %% settings for (optimized) ica (Dimigen, O. (2020). Optimizing the ICA-based removal of ocular artifacts
% % % from free viewing EEG. NeuroImage, https://doi.org/10.1016/j.neuroimage.2019.116117)
% % if PPara.OPTICAT == 1
% %     % step-1: high-pass filter data above 2 hz
% %     cfg          = [];
% %     cfg.hpfilter = PPara.hpfilter;
% %     cfg.hpfreq   = PPara.hpfreq;
% %     data4ICA     = ft_preprocessing(cfg, data4ICA);
% %     
% %     % step-2: overweight SP
% %     % create event-locked epochs to overweight
% %     epoch_range = PPara.SR.*PPara.sac_epoch;
% %     trig = PPara.trig_sac;
% %     epoch = [];
% %     epoch.label = hdr.label;
% %     for i = 1:length(trig)
% %         trlbegin             = trig(i)+epoch_range(1);
% %         trlend               = trig(i)+epoch_range(2)-1;
% %         epoch.trial{i}       = data(:,trlbegin:trlend);
% %         epoch.time{i}        = (trlbegin:trlend)/PPara.SR;
% %         epoch.trl(i,:)       = [trlbegin trlend PPara.cutlen4ica];
% %         epoch.trialinfo(i,:) = [trlbegin trlend PPara.cutlen4ica];
% %     end
% %     % filter epoch to make sure it's the same as data4ICA
% %     cfg          = [];
% %     cfg.channel  = keepchan;  %% only valid MEG channel
% %     cfg.bpfilter = PPara.bpfilter;
% %     cfg.bpfreq   = [PPara.hpfreq PPara.bpfreq(2)];
% %     SP           = ft_preprocessing(cfg, epoch);
% %     clear epoch data
% %     % remove mean, baseline subtracted across whole epoch
% %     SP = cellfun(@(x) x-mean(x,2),SP.trial,'Uni',false);
% %     % overweight (=copy & re-append) overweight-event-locked epochs
% %     SP = cell2mat(SP);
% %     n_tim = length(data4ICA.time{1, 1});
% %     n_trl = length(data4ICA.time);
% %     n_ica = n_epoch*n_trl;
% %     nn = ceil(n_ica/size(SP,2));
% %     SP = repmat(SP,1,nn);
% %     SP = SP(:,1:n_ica);
% %     SP = mat2cell(SP,size(SP,1),n_tim.*ones(1,n_trl));
% %     % merger SP and data4ICA
% %     SP_trialinfo = [(1:n_tim:n_tim*n_trl)' (1:n_tim:n_tim*n_trl)'+n_tim-1 n_tim.*ones(n_trl,1)];
% %     SP_trialinfo(:,[1 2]) = SP_trialinfo(:,[1 2])+data4ICA.trialinfo(end,1)+n_tim-1;
% %     data4ICA.trialinfo = [data4ICA.trialinfo;SP_trialinfo];
% %     SP_sampleinfo = [(1:n_tim:n_tim*n_trl)' (1:n_tim:n_tim*n_trl)'+n_tim-1];
% %     SP_sampleinfo = SP_sampleinfo+data4ICA.sampleinfo(end,1)+n_tim-1;
% %     data4ICA.sampleinfo = [data4ICA.sampleinfo;SP_sampleinfo];
% %     data4ICA.trial = [data4ICA.trial SP];
% %     SP_time = data4ICA.time;
% %     tim_diff = data4ICA.time{1,end}(end)-SP_time{1,1}(1)+0.001;
% %     SP_time = cellfun(@(x) x+tim_diff,SP_time,'Uni',false);
% %     data4ICA.time = [data4ICA.time SP_time]; 
% %     clear SP*
% % end

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

