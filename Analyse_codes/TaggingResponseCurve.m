%%% 20220714  try to plot the resonse curve of RIFT
%%% for each subject, 60Hz coherence was averaged over 200ms after the
%%% fixation onset of words around target [-3 -2 -1 0 1 2 3],also 200ms
%%% after the baseline onset as a 'chance level' for the coherence
% note: since the skip rate is different for words in different location,
% trial number is different across words. However, the averaging time
% window is the same for each words
% note: analysis timewindow are first-pass durations

function TaggingResponseCurve(sid)
%%% do want to do visual rejection for all epochs, so set a global var to
%%% be used in Get_Epoch
global rejectvisual
%%%% running setting
tagfreq = 60;
freq_halfwidth = 5; %%% for hilbert filter

%%% set paths
addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
ft_defaults
rootdir = '/rds/projects/2018/jenseno-reading/Semantic/';

%%% get the exp information
load([rootdir 'Results' filesep 'sv' filesep 'Coh' filesep 'TagCoh.mat']);
subs = TagCoh.subs(TagCoh.SigSubID);
datapath = [rootdir 'Analyse_data' filesep 'sv_' subs{sid} filesep];
sigchanid = TagCoh.SigTagSen_Occip(TagCoh.SigSubID,2); % labels of tagsig sensors
tag_col = find(strcmp(TagCoh.EventHdr,'loc2targ'));
fixdur_col = find(strcmp(TagCoh.EventHdr,'fixation_duration'));
firstpass_col = find(strcmp(TagCoh.EventHdr,'FirstPassFix'));
trig_col = find(strcmp(TagCoh.EventHdr,'fixation_on_MEG'));
cond_col = find(strcmp(TagCoh.EventHdr,'SentenceCondition'));
%%% get semantic violation task condition ID
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
condid = ExpInfo.CondID{1}';
%%% get bad sensors for Get_Epoch
badsensors = ExpInfo.BadSensor(TagCoh.SigSubID); 
clear TagCoh ExpInfo

%%% set the outputs
TagCurve.hdr = {'baseline','n-4','n-3','n-2','n-1','n','n+1','n+2','n+3'};
Loc2Target = [nan -4 -3 -2 -1 0 1 2 3];
n = length(TagCurve.hdr);
TagCurve.coh = zeros(1,n);
TagCurve.trlnum = zeros(1,n);
TagCurve.fixdur_avg = zeros(1,n);
TagCurve.fixdur_std = zeros(1,n);

%% run the coherence analysis, loop over words
for i = 1:n
    clear epoch
    loc = Loc2Target(i);
    if isnan(loc)
        load([datapath 'epoch_BL_Cross.mat']);
        epoch = epoch_BL_Cross;
        clear epoch_BL_Cross
    else 
        % loadin data
        if i == 2
            load([datapath 'data_icaclean.mat']); %data
        end
        load([datapath 'Event.mat']);
        load([datapath 'hdr.mat']);
        % get epochs for the word [-0.5 0.5]s
        AnaTW = 1000; %ms
        PPara.pretrig = -AnaTW/2;
        PPara.posttrig = AnaTW/2;
        tw = [80 AnaTW];
        PPara.SR = 1000; %sampling rate
        PPara.timezero = 'fixation_on_MEG';
        PPara.badsens = badsensors{sid};
        validduration = Event.event_raw(:,fixdur_col)>= tw(1) & Event.event_raw(:,fixdur_col)<= tw(2);
        thistask = arrayfun(@(x) ismember(x,condid),Event.event_raw(:,cond_col));
        thisword = Event.event_raw(:,tag_col) == loc;
        firstpass = Event.event_raw(:,firstpass_col) == 1;
        PPara.event_all = Event.event_raw(validduration&thistask&thisword&firstpass,:);
        % no need to reject trials
        rejectvisual = 0;
        epoch = Get_Epoch(hdr,data,PPara,trig_col);
    end
    
    %%%  normalize pd and remove trials that without photodies signal
    pdi = find(strcmp(epoch.label,'MISC004')); % get pd label index
    rmtrl = [];
    for ttt = 1:length(epoch.trial)
        %%% remove trials that have no pd-004 signal
        if max(epoch.trial{ttt}(pdi,:)) < 0.005
            rmtrl = [rmtrl; ttt];
        end
        epoch.trial{ttt}((pdi),:) = zscore(epoch.trial{ttt}((pdi),:),0,2); % zscore across each trial
    end
    cfg = [];
    cfg.trials = setdiff(1:length(epoch.trial),rmtrl);
    epoch = ft_selectdata(cfg, epoch);
    
    %%% run cohernece
    %%%% just to get a dummy coherence struct
    cfg            = [];
    cfg.output     = 'fourier';
    cfg.channel    = {'MEGGRAD','MISC004'};
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = tagfreq;
    cfg.toi        = -0.5:0.5:0.5;
    cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
    cfg.keeptrials = 'yes';
    cfg.pad        = 'nextpow2';
    fourier = ft_freqanalysis(cfg,epoch);
    %%% get coherence spctrm
    cfg            = [];
    cfg.method     = 'coh';
    cfg.channelcmb = {'MEGGRAD', 'MISC004'}; %% all channell combinations calculated together
    coh            = ft_connectivityanalysis(cfg, fourier); %% the raw coherence
    coh.label      = fourier.label(1:length(coh.labelcmb));
    coh.time       = epoch.time{1};
    
    %%%empty the cohspctrm; do the real coherence using hilbert complex
    coh.cohspctrm = zeros(length(coh.label),length(coh.freq),length(coh.time));
    %%%% get hilbert filert data
    cfg = [];
    cfg.channel    = {'MEGGRAD','MISC004'};
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [tagfreq-freq_halfwidth tagfreq+freq_halfwidth];
    cfg.hilbert    = 'complex';
    cfg.keeptrials = 'yes';
    fltdata = ft_preprocessing(cfg,epoch);
    for chan = 1:length(fltdata.label)-1
        for ttt = 1:length(fltdata.trial)
            sig1(:,ttt) = fltdata.trial{ttt}(chan,:); %time*trl
            sig2(:,ttt) = fltdata.trial{ttt}(end,:);
        end
        spec1 = nanmean(sig1.*conj(sig1),2);%time*1
        spec2 = nanmean(sig2.*conj(sig2),2);
        specX = abs(nanmean(sig1.*conj(sig2),2)).^2;
        coh.cohspctrm(chan,1,:) = specX./(spec1.*spec2);%time*1
    end
    
    %%% finally get the sig occipital tagging sensors
%     mini_fixdur = min(epoch.trialinfo(:,fixdur_col))/1000; %minimum fixation duration over trials, unit in s
    mini_fixdur = 0.2;
    tp = [dsearchn(epoch.time{1}',0):dsearchn(epoch.time{1}',mini_fixdur)];
    cohcoh = squeeze(coh.cohspctrm);%% chan*time
    cohcoh = nanmean(nanmean(coh.cohspctrm(sigchanid{sid},tp),1),2); %averaging over chan and time
    TagCurve.coh(1,i) = cohcoh; 
    TagCurve.trlnum(1,i) = size(epoch.trialinfo,1);
    TagCurve.fixdur_avg(1,i) = nanmean(epoch.trialinfo(:,fixdur_col));
    TagCurve.fixdur_std(1,i) = nanstd(epoch.trialinfo(:,fixdur_col));
end
save([rootdir 'Results' filesep 'sv' filesep 'Coh' filesep 'TagCurve_' num2str(sid)],'TagCurve')


%% %%%%%%%%%%%%%%%% combine and plot the curve locally %%%%%%%%%%%%%%%%%%
cd Z:\Semantic\Results\sv\Coh
TagCurve_all.hdr = {'baseline','n-4','n-3','n-2','n-1','n','n+1','n+2','n+3'};
n = length(TagCurve_all.hdr);
nsub = 29;
TagCurve_all.coh = zeros(nsub,n);
TagCurve_all.trlnum = zeros(nsub,n);
TagCurve_all.fixdur_avg = zeros(nsub,n);
TagCurve_all.fixdur_std = zeros(nsub,n);

for sid = 1:nsub
    load(['TagCurve_' num2str(sid)])
    TagCurve_all.coh(sid,:) = TagCurve.coh;
    TagCurve_all.trlnum(sid,:) = TagCurve.trlnum;
    TagCurve_all.fixdur_avg(sid,:) = TagCurve.fixdur_avg;
    TagCurve_all.fixdur_std(sid,:) = TagCurve.fixdur_std;
end
clear TagCurve
TagCurve = TagCurve_all;
clear TagCurve_all
save TagCurve TagCurve

%%% plot, bar with error
err = std(TagCurve.coh,0,1)/sqrt(nsub);
meancoh = mean(TagCurve.coh,1);
fignam = 'Tagging response curve, n--target';
h = figure('Name',fignam,'color',[1 1 1],'Position',[100 100 320 150]);
errorbar(meancoh(2:end),err(2:end),'color','k','LineWidth',1)
xticklabels();
ylim([0.003 0.018])
hold on; 
plot([1 n-1],[meancoh(1) meancoh(1)],'k-.')
set(gca,'box','off','LineWidth',1)
set(gca,'XTickLabel',TagCurve.hdr(2:end),'FontSize',7,'FontWeight','normal','FontName','Arial')
%%% add p value
for i = 2:n
    [~,p,~,stats] = ttest(TagCurve.coh(:,1), TagCurve.coh(:,i));
    text(i-1,meancoh(i)+0.002,num2str(round(1000*p)/1000));
end
xlim([0.5 n-1])
saveas(h,'TagCurve.fig')
set(h,'renderer', 'painters')
saveas(h,'TagCurve','svg');


%% estimate the word length for sig 
cd 'Z:\Semantic\Results\sv\Coh\'
load('Z:\Semantic\PTB_codes\SentMat_1.mat')
load('TargCondPair_1.mat')
location = [-4 -3 -2 -1 0 1 2 3];
sv = find(TargCondPair(:,2)>10 & TargCondPair(:,2)<13);
wrd_len = [];
for i = sv'
    i_target = TargCondPair(i,1);
    idx = i_target+location;
    head = 0;
    if idx(1)<1
       idx(1) = []; 
       head = 1;
    end
    wrd = SentMat(i,idx);
    tmp = cellfun(@length,wrd);
    if head
        tmp = [nan tmp];
    end
    wrd_len = [wrd_len;tmp];
end
nanmean(wrd_len,1)
TagCurve.word_length = wrd_len;
save TagCurve TagCurve

span = sum(wrd_len(:,2:4),2); %%from n-3 to n-1
mean(span)
std(span)










































































