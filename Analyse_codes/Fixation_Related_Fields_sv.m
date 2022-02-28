% copy from Lexical/Analyse_codes
% 20210722 get ERFs for semanviol and orthofreq task

function Fixation_Related_Fields_sv(sid)
server = 1;
ddd = 1; %DataSets = {'sv','of','fa'};
%%%% running setting
RunCond = 'WrdOn';
EpochType = {'PreTarg';'Targ';'PosTarg'};
targid = [-1 0 1];

% % %%%% basic settingup
if server
    rootdir = '/rds/projects/2018/jenseno-reading/Semantic/';
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
else
    %%% local
    rootdir = 'Z:\Semantic\';
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    %%% for colormap
    addpath(genpath('Z:\Semantic\Analyse_codes\plotfig\'));
end

%%get the subinformation
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
DS = ExpInfo.DSName{ddd};
conds = ExpInfo.CondID{ddd};
condnm = ExpInfo.CondName{ddd};
condcol = find(strcmp(ExpInfo.EventHdr,'SentenceCondition'));
loccol = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
fixdurcol = find(strcmp(ExpInfo.EventHdr,'fixation_duration'));
firstpass_col = find(strcmp(ExpInfo.EventHdr,'FirstPassFix'));
ERF.EventHdr = ExpInfo.EventHdr;
eval(['subjects = ExpInfo.subjects.' DS ';']);
% remove empty row in subjects
subjects = subjects(~strcmp(subjects,''));

% path
PPath.FigPath = [rootdir 'Results' filesep DS filesep 'ERF_FirstFixations' filesep]; %result path 'ERF_AllFixations'
if (~exist(PPath.FigPath,'dir'))
    mkdir(PPath.FigPath);
end


%% ========= subject loop
if server == 1
    sub = subjects{sid};
    PPath.SaveData = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
    
    fprintf(['********** s' num2str(sid) ' ******** \n\n']);
    ERF.sub = sub;
    load([PPath.SaveData 'epoch_' RunCond],'epoch');
    
    %%%% ERF: filter data (NO baseline correction here)
    cfg                = [];
    cfg.channel        = 'MEGGRAD';
    cfg.lpfilter       = 'yes';
    cfg.lpfreq         = 30;
    cfg.demean         = 'yes';
    epoch              = ft_preprocessing(cfg, epoch);
    
%     %%% select the first pass fixations, the first fixations on a given
%     %%% word
%     cfg = [];
%     cfg.trials = find(epoch.trialinfo(:,firstpass_col)==1);
%     epoch = ft_selectdata(cfg, epoch);
    
    %%%%%%%%%%%%%============== epoch loop
    for mmm = 1:length(EpochType)
        %%%% get equal trial number across all conditions
        %%%% seperate trials based on target-freq
        idx_1 = find(epoch.trialinfo(:,loccol)==targid(mmm) & epoch.trialinfo(:,condcol) == conds(1));
        idx_2 = find(epoch.trialinfo(:,loccol)==targid(mmm) & epoch.trialinfo(:,condcol) == conds(2));
        
        %%%%%%%%%%=========== ERF =========%%%%%%%%%%%
        cfg            = [];
        cfg.trials     = idx_1;
        erf_1          = ft_timelockanalysis(cfg, epoch);
        cfg.trials     = idx_2;
        erf_2          = ft_timelockanalysis(cfg, epoch);
        
        %%% baseline correction for each condition of each EpochType
        cfg          = [];
        cfg.baseline = [-0.2 0];
        erf_1        = ft_timelockbaseline(cfg, erf_1);
        erf_2        = ft_timelockbaseline(cfg, erf_2);
        
        %%% combine planar sensors
        erf_1 = ft_combineplanar([],erf_1);
        erf_2 = ft_combineplanar([],erf_2);
        
        %%% get trl numbers in each condition
        for cid = 1:length(conds)
            eval(['ERF.' EpochType{mmm} '_TrlInfo{1,cid} = epoch.trialinfo(idx_' num2str(cid) ',:);']);
            eval(['ERF.' EpochType{mmm} '_RT(1,cid) = nanmean(epoch.trialinfo(idx_' num2str(cid) ',fixdurcol));']);
        end
        
        %%%============== save data out for group analysis
        eval(['ERF.' EpochType{mmm} '_' condnm{1} '= erf_1;']);
        eval(['ERF.' EpochType{mmm} '_' condnm{2} '= erf_2;']);
    end
    save([PPath.FigPath 'ERF_' num2str(sid)],'ERF','-v7.3')
    fprintf(['********** s' num2str(sid) ': done ******** \n\n']);
    
else
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run locally get the whole structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% combine data together into one big structure
    load([PPath.FigPath 'ERF_1']);
    ss = 0;
    for ps = 1:length(subjects)
        load([PPath.FigPath 'ERF_' num2str(ps)]);
        fprintf(['*** s ' num2str(ps) ' **** \n\n']);
        ss = ss + 1;
        ERF_all.subs{ss,:} = ERF.sub;
        for fd = 1:length(EpochType)
            eval(['ERF_all.' EpochType{fd} '_RT(ss,:)= ERF.' EpochType{fd} '_RT;']);
            eval(['ERF_all.' EpochType{fd} '_TrlInfo(ss,:)= ERF.' EpochType{fd} '_TrlInfo;']);
            eval(['ERF_all.' EpochType{fd} '_' condnm{1} '{ss}= ERF.' EpochType{fd} '_' condnm{1} ';']);
            eval(['ERF_all.' EpochType{fd} '_' condnm{2} '{ss}= ERF.' EpochType{fd} '_' condnm{2} ';']);
        end
    end
    clear ERF
    ERF = ERF_all;
    clear ERF_all
    
    %% %%%%%================ group analysis
    for mmm = 1:length(EpochType)
        eval(['erf_1 = ERF.' EpochType{mmm} '_' condnm{1} ';']);
        eval(['erf_2 = ERF.' EpochType{mmm} '_' condnm{2} ';']);
        
        %%% permutation test
        layout = 'neuromag306cmb'; %'neuromag306planar';
        load([rootdir layout '_neighb.mat'])
        cfg                  = [];
        cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
        cfg.parameter        = 'avg';
        cfg.statistic        = 'depsamplesT'; % use the dependent samples T-statistic as a measure to evaluate the effect at the sample level
        cfg.correctm         = 'cluster';
        cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under thepermutation distribution.
        cfg.minnbchan        = 2; % minimum number of neighborhood channels that isrequired for a selected sample to be includedin the clustering algorithm (default=0).
        cfg.neighbours       = neighbours; % neighbours;   % see below
        cfg.tail             = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test
        cfg.clustertail      = 0;
        cfg.alpha            = 0.025; % alpha level of the permutation test
        cfg.clusteralpha     = 0.025; % alpha level of the sample-specific test statistic thatwill be used for thresholding
        cfg.numrandomization = 5000; % number of draws from the permutation distribution
        cfg.latency          = [0.2 0.5];
        %%% design matrix
        nsubj = length(ERF.subs);
        design = [1:nsubj 1:nsubj; ones(1,nsubj) ones(1,nsubj)*2];
        cfg.design = design;
        cfg.uvar   = 1;
        cfg.ivar   = 2;
        stat = ft_timelockstatistics(cfg, erf_1{:}, erf_2{:});
        any(any(any(stat.mask)))
        %%% stat over the whole 1s time window
        eval(['ERF.' EpochType{mmm} '_stat = stat;']);
        
        %%%============== plot figure
        %%% get grandaverage across subjects
        erf_1 = ft_timelockgrandaverage([], erf_1{:});
        erf_2 = ft_timelockgrandaverage([], erf_2{:});
        %%% remove baseline again after combining planar sensors to get a
        %%% near-zero baseline for the plot
        cfg          = [];
        cfg.baseline = [-0.2 0];
        erf_1        = ft_timelockbaseline(cfg, erf_1);
        erf_2        = ft_timelockbaseline(cfg, erf_2);
        erf_1.mask = zeros(size(erf_1.avg));
        timid = arrayfun(@(x) ismember(x,stat.time),erf_1.time,'Uni',true);
        erf_1.mask(:,timid) = stat.mask;
        erf_2.mask = erf_1.mask ;
        %%%==== plot
        figtitle = ['FRFs_' EpochType{mmm}];
        cfg = [];
        cfg.baseline = 'no'; % did it in the above step
        cfg.showlabels = 'no';
        cfg.layout = [layout '.lay'];
        if any(any(stat.mask))
            cfg.mask = 'yes';
            cfg.maskparameter = 'mask';
            cfg.maskstyle = 'box';
            cfg.maskfacealpha = 1;
        end
        h = figure('name',figtitle,'color',[1 1 1]);
        ft_multiplotER(cfg, erf_1, erf_2);
        title([EpochType{mmm} 'et Fixation-Related-Fields'],'FontSize',16);
        saveas(h,[PPath.FigPath figtitle]);
    end
    save([PPath.FigPath 'ERF_rmbl_cmb'],'ERF','-v7.3');
      
end
































