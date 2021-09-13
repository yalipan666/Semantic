%%% edit: Yali Pan, 20210706, Just a quick check of the tag signal

%%% get sub_information
subjects = {'20210911_b3f9'};
%%%% addpath fieldtrip
addpath Z:\fieldtrip-20200220\
ft_defaults
%%% paths
rootdir = 'Z:\Semantic\';
path_meg = [rootdir 'RawData\MEG_data\'];
path_ptb = [rootdir 'RawData\PTB_data\'];
figpath = [rootdir 'Results' filesep 'Topo_tagging' filesep];

for sss = 1:length(subjects)
    sub = subjects{sss};
    %%% file parameters
    filename = sub(10:13);
    filedate = sub(3:8);
    %%%% parameters initialization
    megpath = [path_meg sub filesep filedate filesep];
    dataset = [megpath filename];
    trig_value = 4;
    posttrig = 2; % length of epoch
    
    %%% 1Hz freq step
    tagfreq = 60;
    pretrig = 0;
    SR = 1000; %% sampling rate
    pretrig = pretrig*SR;
    posttrig = posttrig*SR;  %epoch end timepoint, aligned with marker
    
    %%Step-0: get data and trigger
    cfg         = [];
    cfg.dataset = [dataset '.fif'];
    hdr         = ft_read_header(cfg.dataset);
    event       = ft_read_event(cfg.dataset);
    data        = ft_read_data(cfg.dataset);
    Trigger_MEG = [[event(strcmp('Trigger',{event.type})).value]' [event(strcmp('Trigger',{event.type})).sample]'];
    %%% appenddata
    cfg         = [];
    cfg.dataset = [dataset '-1.fif'];
    data2  = ft_read_data(cfg.dataset);
    event2 = ft_read_event(cfg.dataset);
    Trigger_MEG2 = [[event2(strcmp('Trigger',{event2.type})).value]' [event2(strcmp('Trigger',{event2.type})).sample]'];
    data1length = size(data,2);
    Trigger_MEG2(:,2) = Trigger_MEG2(:,2)+data1length;
    data = [data data2];
    Trigger_MEG = [Trigger_MEG; Trigger_MEG2];
    clear data2 Trigger_MEG2 event event2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% normalizing all the photodiode
    pd = find(strcmp(hdr.label,'MISC004'));
    if sss == 1
        pd = find(strcmp(hdr.label,'MISC005'));
    end
    data(pd,:) = zscore(data(pd,:),0,2);
    
    %%Step-1: define trials mannually, epoch aligned with sentence_onset
    Trigger_MEG_cond = Trigger_MEG(Trigger_MEG(:,1)==trig_value,:);
    % define trials manually
    trlidx = 0;
    epoch = [];
    for j = 1:length(Trigger_MEG_cond)
        trlbegin     = Trigger_MEG_cond(j,2) + pretrig;
        trlend       = Trigger_MEG_cond(j,2) + posttrig;
        trlidx       = trlidx +1;
        epoch.trial{trlidx}        = data(:,trlbegin:(trlend-1));
        epoch.time{trlidx}         = (pretrig:(posttrig-1))/SR;
        epoch.trl(trlidx,:)        = [trlbegin trlend pretrig Trigger_MEG_cond(j,2)];
        epoch.trialinfo(trlidx,:)  = Trigger_MEG_cond(j,2);
    end
    clear data
    
    %%Step-2: get the trial strcture from standard fieldtrip strcture
    epoch.label = hdr.label;
    epoch.hdr = hdr;
    cfg = [];
    cfg.trl =  [epoch.trl];
    epoch_all = ft_redefinetrial(cfg,epoch);
    %%% trials in epoch_cond is wrong, so replace them with epoch.trial
    epoch_all.trial = epoch.trial;
    epoch_all.time = epoch.time;
    clear epoch
    
    %%%======= get fourier spctrm
    cfg         = [];
    cfg.output  = 'fourier';
    cfg.channel = {'MEGGRAD','MISC004'};
    cfg.method  = 'mtmfft';
    cfg.taper   = 'hanning';
    cfg.foi     = 1:1:100;
    fourier=ft_freqanalysis(cfg,epoch_all); %=% tfr.powspctrm:3D data: trialnum*channelnum*freqpoint
    
    %%% statistical test of coherence
    %%%======= get coherence spctrm
    cfg            = [];
    cfg.method     = 'coh';
    cfg.channelcmb = {'MEGGRAD','MISC004'}; %% all channell combinations calculated together
    coh            = ft_connectivityanalysis(cfg, fourier);
    
    %%% combining the coherence of two gradiometers
    cfg = [];
    cfg.method = 'sum';
    coh.powspctrm = coh.cohspctrm;
    coh.label = fourier.label(1:length(coh.labelcmb));
    coh.grad = hdr.grad;
    coh_cmb = ft_combineplanar(cfg, coh); %cfg.method = 'sum' only works for frequency data with powspctrm
    coh_cmb = rmfield(coh_cmb,'cohspctrm');
    
    %%%%%%%%%%%%%%%%% plot figures
    %%% plot topograph
    figtitle = 'Coh at 60Hz';
    freqx = tagfreq;
    h = figure('color', [1 1 1],'name',['Topo_' figtitle],'position',[100 100 600 600]);
    cfg = [];
    cfg.layout = 'neuromag306cmb.lay';
    cfg.xlim = [freqx freqx];
    cfg.marker = 'numbers';
    ft_topoplotER(cfg,coh_cmb);
    colormap jet; colorbar;
    title(figtitle,'FontSize',16);
    saveas(h,[figpath sub '_' figtitle],'bmp');
end


