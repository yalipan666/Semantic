% copy from Lexical/Analyse_codes
% 20210722 get ERFs for semanviol and orthofreq task

function Fixation_Related_Fields(sid)
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
fixdurcol = find(strcmp(ExpInfo.EventHdr,'fixation_duration'));
ERF.EventHdr = ExpInfo.EventHdr;
PPath.FigPath = [rootdir 'Results' filesep DS filesep 'FRF' filesep]; %result path
eval(['subjects = ExpInfo.subjects.' DS ';']);


%% ========= subject loop
if server == 1
    sub = subjects{sid};
    PPath.SaveData = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
    
    fprintf(['********** s' num2str(sid) ' ******** \n\n']);
    ERF.sub = sub;
    load([PPath.SaveData 'epoch_' RunCond],'epoch');
    
    %%%% ERF: filter data(NO baseline correction)
    cfg          = [];
    cfg.channel  = 'MEGGRAD';
    cfg.lpfilter = 'yes';
    cfg.lpfreq   = 35;
    cfg.demean   = 'yes';
    epoch        = ft_preprocessing(cfg, epoch);
    
    %%%%%%%%%%%%%============== epoch loop
    for mmm = 1:length(EpochType)
        %%%% get equal trial number across all conditions
        %%%% seperate trials based on target-freq
        idx_1 = find(epoch.trialinfo(:,3)==targid(mmm) & epoch.trialinfo(:,condcol) == conds(1));
        idx_2 = find(epoch.trialinfo(:,3)==targid(mmm) & epoch.trialinfo(:,condcol) == conds(2));
        n_1 = length(idx_1);
        n_2 = length(idx_2);
        if n_1 > n_2
            idx_1 = idx_1(randperm(n_1));
            idx_1 = idx_1(1:n_2);
        elseif n_1 < n_2
            idx_2= idx_2(randperm(n_2));
            idx_2 = idx_2(1:n_1);
        end
        
        %%%%%%%%%%=========== ERF =========%%%%%%%%%%%
        cfg        = [];
        cfg.trials = idx_1;
        erf_1    = ft_timelockanalysis(cfg, epoch);
        cfg.trials = idx_2;
        erf_2   = ft_timelockanalysis(cfg, epoch);
        
        %%% baseline correction
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
            eval(['ERF.' EpochType{mmm} '_TrlIdx_equ{1,cid} = idx_' num2str(cid) ';']);
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
    %%% only run group analysis with 39 subjects who has both datasets
    for ps = 1:17
        load([PPath.FigPath 'ERF_' num2str(ps)]);
        fprintf(['*** s ' num2str(ps) ' **** \n\n']);
        ss = ss + 1;
        ERF_all.subs{ss,:} = ERF.sub;
        for fd = 1:length(EpochType)
            eval(['ERF_all.' EpochType{fd} '_TrlIdx_equ(ss,:)= ERF.' EpochType{fd} '_TrlIdx_equ;']);
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
        load([rootdir 'neuromag306cmb_neighb.mat'])
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
        cfg.numrandomization = 1000; % number of draws from the permutation distribution
        %         cfg.latency          = [0 0.25]; %[0.25 0.5];
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
        erf_1.mask = stat.mask;
        erf_2.mask = stat.mask;
        %%%==== plot
        figtitle = [DS '_' RunCond '_FRFs_' EpochType{mmm}];
        cfg = [];
        cfg.baseline = 'no';
        cfg.showlabels = 'no';
        cfg.layout = 'neuromag306cmb.lay';
        cfg.mask = 'yes';
        cfg.maskparameter = 'mask';
        cfg.maskstyle = 'box';
        cfg.maskfacealpha = 1;
        h = figure('name',figtitle,'color',[1 1 1]);
        ft_multiplotER(cfg, erf_1, erf_2);
        title([EpochType{mmm} 'et Fixation-Related-Fields'],'FontSize',16);
        set(gcf, 'renderer', 'painters')
        saveas(h,[PPath.FigPath figtitle]);
    end
    save([PPath.FigPath RunCond '_ERF_rmbl'],'ERF','-v7.3');
    
    
    %% ====== plot topo for the sig sensors for Pre-target (only pre-targ sig)
    %%% get ERF sig labels
    sig_tp = sum(ERF.PreTarg_stat.mask,1)~=0;
    SigTim = ERF.PreTarg_stat.time(sig_tp);
    sigchanid = sum(ERF.PreTarg_stat.posclusterslabelmat(:,sig_tp),2)~=0;
    SigSen = ERF.PreTarg_stat.label(sigchanid);
    
    %%% erf_diff
    eval(['erf_1  = ft_timelockgrandaverage([], ERF.PreTarg_' condnm{1} '{:});']);
    eval(['erf_2 = ft_timelockgrandaverage([], ERF.PreTarg_' condnm{2} '{:});']);
    cfg = [];
    cfg.parameter = 'avg';
    cfg.operation = 'x1-x2';
    erf_diff = ft_math(cfg, erf_1, erf_2);
    
    % plot
    figtitle = 'Topo_PreTarg_sigchan';
    h = figure('color', [1 1 1],'name',figtitle,'Position',[100 100 160 160]);
    cfg =[];
    cfg.layout = 'neuromag306cmb.lay';
    cfg.comment = ' ';
    cfg.gridscale = 360;
    cfg.contournum = 0;
    cfg.marker = 'off';
    cfg.highlight= 'on';
    cfg.highlightchannel = SigSen;
    cfg.highlightsymbol = 'o';
    cfg.highlightcolor  = [0 0 0];
    cfg.highlightsize   = 4;
    % only select sig time window
    cfg.xlim = [SigTim(1) SigTim(end)];
    cfg.zlim = 'maxabs';
    ft_topoplotER(cfg,erf_diff);
    colorbar;
    title('FRFs sig sensors','FontWeight','normal','FontSize',7);
%     set(gcf, 'renderer', 'painters')
%     saveas(h,[PPath.FigPath figtitle]);
%     saveas(h,[PPath.FigPath figtitle],'svg');
    
    %% % plot the sig FRPs curves for pre-targets
    SigTim = ERF.PreTarg_stat.time(sum(ERF.PreTarg_stat.mask,1)~=0);
    n_ch = length(SigSen);
    erf_tim = ERF.PreTarg_stat.time;
    n_tim = length(erf_tim);
    n_sub = length(ERF.subs);
    erf_all = nan(2,n_tim,n_sub);
    for s = 1:n_sub
        label_all = ERF.PreTarg_incong{1,s}.label;
        sigchan_id = cellfun(@(x) find(strcmp(x,label_all)),SigSen,'Uni',true);
        tmp_1 = mean(ERF.PreTarg_incong{1,s}.avg(sigchan_id,:),1);
        tmp_2 = mean(ERF.PreTarg_cong{1,s}.avg(sigchan_id,:),1);
        erf_all(:,:,s) = [tmp_1; tmp_2];
    end
    %% % to transfer the unit from T to fT
    erf_all = erf_all.*10e13;
    % plot
    colmat = [0 114 189;217 83 25]./255;
    figtitle = 'Pre-target averaged FRFs over sig sensors';
    h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 160 160]);
    meanERF = mean(erf_all,3);
    seERF = nanstd(erf_all,0,3)./sqrt(n_sub);
    a = shadedErrorBar(erf_tim,meanERF(1,:),seERF(1,:),{'color',colmat(1,:)},0.8);hold on;
    b = shadedErrorBar(erf_tim,meanERF(2,:),seERF(2,:),{'color',colmat(2,:)},0.9);
    a.mainLine.LineWidth = 1.5;
    b.mainLine.LineWidth = 1.5;
    set(gca,'box','off','LineWidth',1)
%     xlim([-0.3 0.5]);
    set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
    title('Pre-target','FontSize',7,'FontName','Arial','FontWeight','normal');
    xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial','FontWeight','normal');
    ylabel('Fixation Related Fields (fT/cm)','FontSize',7,'FontWeight','normal');
    % plot the significant line
    plot(SigTim,repmat(7,size(SigTim)),'-k','LineWidth',1);
    text(0.42,5.8,'*','FontWeight','normal','FontSize',14,'FontName','Arial')
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
%     saveas(h,[PPath.FigPath figtitle],'svg');
    
    %% plot the sig FRPs curves for targets
    n_ch = length(SigSen);
    erf_tim = ERF.Targ_stat.time;
    n_tim = length(erf_tim);
    n_sub = length(ERF.Targ_1);
    erf_all = nan(2,n_tim,n_sub);
    for s = 1:n_sub
        label_all = ERF.Targ_1{1,s}.label;
        sigchan_id = cellfun(@(x) find(strcmp(x,label_all)),SigSen,'Uni',true);
        tmp_1 = mean(ERF.Targ_1{1,s}.avg(sigchan_id,:),1);
        tmp_2 = mean(ERF.Targ_2{1,s}.avg(sigchan_id,:),1);
        erf_all(:,:,s) = [tmp_1; tmp_2];
    end
    %%% to transfer the unit from T to fT
    erf_all = erf_all.*10e12;
    % plot
    colmat = [0 114 189;217 83 25]./255;
    figtitle = 'Target averaged FRFs';
    h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 160 160]);
    meanERF = mean(erf_all,3);
    seERF = nanstd(erf_all,0,3)./sqrt(n_sub);
    a = shadedErrorBar(erf_tim,meanERF(1,:),seERF(1,:),{'color',colmat(1,:)},0.8);hold on;
    b = shadedErrorBar(erf_tim,meanERF(2,:),seERF(2,:),{'color',colmat(2,:)},0.9);
    a.mainLine.LineWidth = 1.5;
    b.mainLine.LineWidth = 1.5;
    set(gca,'box','off','LineWidth',1)
    xlim([-0.3 0.5]);
    set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
    title('Target','FontSize',7,'FontName','Arial','FontWeight','normal');
    xlabel('Time (s)','FontWeight','normal','FontSize',7,'FontName','Arial','FontWeight','normal');
    ylabel('Fixation Related Fields (fT/cm)','FontSize',7,'FontWeight','normal');
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'svg');
    
    
    
    % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%% get ERFs for tagging sensors and run ttest %%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% basic setting up
    % rootdir = 'Z:\';
    % PPath.FigPath = [rootdir 'Results' filesep 'ERF' filesep];
    % EpochType = {'PreTarg';'Targ'}; %;'PosTarg'};
    %
    % %%% load coh_word_freq: TagCoh
    % load([rootdir 'Results' filesep 'Coh_hilbert_01' filesep 'minirt_trl' filesep 'Cmb.mat'])
    % %%% only run coh_tag sig sub
    % Coh_sigchan = TagCoh.SigTagSen_Occip(TagCoh.SigSubID,:); %% 43*2
    % subjects = TagCoh.subs(TagCoh.SigSubID);
    % minrt_all = [TagCoh.PreTarg_Ttest.trl_minrt TagCoh.Targ_Ttest.trl_minrt];
    % clear TagCoh
    %
    % %%% load ERFs
    % load('Z:\Results\ERF\WrdOn_ERF_rmbl.mat')
    %
    % %%% get var in ERF
    % sigsubid = cellfun(@(x) find(strcmp(x, ERF.subs(:,1))),subjects,'Uni',true);
    % n_sub = length(sigsubid);
    % erf_tim = ERF.PreTarg_incong{1,1}.time;
    % n_tim = length(erf_tim);
    % tim0 = dsearchn(erf_tim', 0);
    % n_epty = length(EpochType);
    %
    % %%% save out
    % ERF_tagsigchan.subjects = subjects;
    % ERF_tagsigchan.Coh_sigchan = Coh_sigchan;
    % ERF_tagsigchan.time = erf_tim;
    %
    % for eee = 1:n_epty
    %     data4Ttest = nan(n_sub,2); %[low high]
    %     ERF_tagsigchan = nan(n_tim,2,n_sub);
    %     for sid = 1:n_sub
    %         sid = sigsubid(sid); %%%sub_id in ERF
    %         eval(['labels = ERF.' EpochType{eee} '_1{1,sid}.label;']); % 101 combined planar sensors
    %         eval(['erf_1 = ERF.' EpochType{eee} '_1{1,sid}.avg;']);
    %         eval(['erf_2 = ERF.' EpochType{eee} '_2{1,sid}.avg;']);
    %         timend = dsearchn(erf_tim', minrt_all(sid,eee));
    %         %%% get sig chanid in erf
    %         sigchan = Coh_sigchan{sid,1}; % seperate planar sensor
    %         chanid_erf = cellfun(@(x) find(strncmp(x,labels,6)),sigchan,'Uni',true);
    %         %%% save out data
    %         if length(chanid_erf) == 1
    %            tmp_1 = erf_1(chanid_erf,:);
    %            tmp_2 = erf_2(chanid_erf,:);
    %         else
    %             tmp_1 = mean(erf_1(chanid_erf,:));
    %             tmp_2 = mean(erf_2(chanid_erf,:));
    %         end
    %         ERF_tagsigchan(:,:,sid) = [tmp_1' tmp_2'];
    %         data4Ttest(sid,:) = [mean(tmp_1(tim0:timend)) mean(tmp_2(tim0:timend))];
    %     end
    %     %%% run ttest
    %     [~,p,~,stats] = ttest(data4Ttest(:,1), data4Ttest(:,2));
    %     stats.p = p; disp([EpochType{eee} '_Ttest = ' num2str(p)])
    %     %%% save data
    %     eval(['ERF_tagsigchan.' EpochType{eee} '_ERF_tagsigchan = ERF_tagsigchan;']);
    %     eval(['ERF_tagsigchan.' EpochType{eee} '_data4Ttest = data4Ttest;']);
    %     eval(['ERF_tagsigchan.' EpochType{eee} '_minrt = minrt;']);
    %     eval(['ERF_tagsigchan.' EpochType{eee} '_stats = stats;']);
    % end
    % save([PPath.FigPath 'ERF_tagsigchan.mat'], 'ERF_tagsigchan','-v7.3');
    %
    % %% ===== plot figure
    % %%% paras for violin plot
    % colmat = [0 114 189;217 83 25]./255;
    % figtitle = 'ERFs_tagchan_Ttest_violin';
    % group = [cellstr(repmat('Low',n_sub,1)); cellstr(repmat('High',n_sub,1))];
    % grouporder = {'Low','High'};
    % figure('Name',figtitle,'color',[1 1 1]);
    %
    % for eee = 1:n_epty
    %     eval(['ERF_tagsigchan = ERF_tagsigchan.' EpochType{eee} '_ERF_tagsigchan;']);
    %     eval(['data4Ttest = ERF_tagsigchan.' EpochType{eee} '_data4Ttest;']);
    %     eval(['p = ERF_tagsigchan.' EpochType{eee} '_stats.p;']);
    %
    %     %%% full ERF curves
    %     meanERF = mean(ERF_tagsigchan,3);
    %     seERF = nanstd(ERF_tagsigchan,0,3)./sqrt(n_sub);
    %     h = subplot(2,n_epty,eee);
    %     a = shadedErrorBar(erf_tim,meanERF(:,1),seERF(:,1),{'color',colmat(1,:)},0.8);hold on;
    %     b = shadedErrorBar(erf_tim,meanERF(:,2),seERF(:,2),{'color',colmat(2,:)},0.9);
    %     a.mainLine.LineWidth = 1.5;
    %     b.mainLine.LineWidth = 1.5;
    %     set(gca,'box','off','LineWidth',1.5)
    %     set(gca,'XTick',-0.4:0.2:0.4);
    %     set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'})
    %     set(gca,'FontSize',10,'FontWeight','bold','FontName','Arial');
    %     title([EpochType{eee} 'et, tagging sensors'],'FontSize',14,'FontName','Arial')
    %     xlabel('Time (s)','FontWeight','normal','FontSize',12,'FontName','Arial','FontWeight','bold');
    %     ylabel('Fixation Related Fields (T/cm)','FontSize',12,'FontWeight','bold');
    %
    %     %%% violin plot for ttest
    %     vdata = [data4Ttest(:,1); data4Ttest(:,2)];
    %     h = subplot(2,n_epty,eee+n_epty);
    %     vp = violinplot(vdata, group,'GroupOrder',grouporder);
    %     vp(1).ViolinColor = colmat(1,:);
    %     vp(2).ViolinColor = colmat(2,:);
    %     set(gca,'FontSize',10,'FontWeight','normal','FontName','Arial','FontWeight','bold');
    %     set(gca,'box','off','LineWidth',1.5)
    %     xlabel('Target frequency','FontSize',12,'FontWeight','bold');
    %     ylabel('Fixation Related Fields (T/cm)','FontSize',12,'FontWeight','bold');
    %     title([EpochType{eee} 'et, p = ' num2str(round(1000*p)/1000)],'FontSize',14,'FontWeight','bold');
    %     %%% plot the line linking each subject
    %     hold on;
    %     x1 = vp(1,1).ScatterPlot.XData;
    %     y1 = vp(1,1).ScatterPlot.YData;
    %     x2 = vp(1,2).ScatterPlot.XData;
    %     y2 = vp(1,2).ScatterPlot.YData;
    %     plot([x1; x2],[y1; y2],'Color',[.7 .7 .7])
    % end
    % %%
    % set(gcf, 'renderer', 'painters')
    % saveas(h,[PPath.FigPath figtitle]);
    % saveas(h,[PPath.FigPath figtitle],'svg');
end
































