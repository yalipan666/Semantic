%%%% 2022-01-18
%%%% using hilbert to compute coherence

function TaggingResponse_Coh_sv(sid)
%%% choosing significant sensors based on epoch type eee
%%% using epoch aligned with rrr
%%% running dataset ddd

server = 1;
eee = 1; % epoch index for sensor selection;{'PreTarg';'Targ';'PosTarg'};
rrr = 1; % running condition index;{'WrdOn','WrdOff'};
ddd = 1; % dataset index;{'sv','of','fa'};

%%%% running setting
runconds = {'WrdOn','WrdOff'};
RunCond = runconds{rrr};
DataSets = {'sv','of','fa'};
EpochType = {'PreTarg';'Targ';'PosTarg'};
targid = [-1 0 1];
tagfreq = 60; 
freqrange = 40:2:80;
freq_halfwidth = 5; %%% for hilbert filter
tagfreq_id = dsearchn(freqrange',tagfreq);
SenSelectP = 0.01;

%%% set paths
if server
    addpath /rds/projects/2018/jenseno-reading/fieldtrip-20200220/
    ft_defaults
    rootdir = '/rds/projects/2018/jenseno-reading/Semantic/';
else
    addpath Z:\fieldtrip-20200220\
    ft_defaults
    rootdir = 'Z:\Semantic\';
end
addpath(genpath([rootdir 'Analyse_codes' filesep 'plotfig' filesep]))

%%% get the exp information
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
DS = ExpInfo.DSName{ddd};
conds = ExpInfo.CondID{ddd};
condnm = [{'all'},ExpInfo.CondName{ddd}];
cond_col = find(strcmp(ExpInfo.EventHdr,'SentenceCondition'));
tag_col = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
fixdur_col = find(strcmp(ExpInfo.EventHdr,'fixation_duration'));
firstpass_col = find(strcmp(ExpInfo.EventHdr,'FirstPassFix'));
eval(['subjects = ExpInfo.subjects.' DS ';']);
% remove empty row in subjects
subjects = subjects(~strcmp(subjects,''));

%%% paths
PPath.FigPath = [rootdir 'Results' filesep DS filesep 'Coh' filesep]; %result path
if (~exist(PPath.FigPath,'dir'))
    mkdir(PPath.FigPath);
end

%%% set output
TagCoh = [];
TagCoh.CondName = condnm;
TagCoh.freq_halfwidth = freq_halfwidth;
TagCoh.EventHdr = ExpInfo.EventHdr;
TagCoh.subs = subjects;
TagCoh.SigTagSen_Occip = cell(size(subjects));
TagCoh.Coh_hdr = {'freq*time*sub*cond'};

if server == 1
    %%%%%%===================== subject loop
    fprintf(['***** analyzing: ' DataSets{ddd} '-s' num2str(sid) '**** \n\n']);
    sub = subjects{sid};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% 1. loading data %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%======== load in orth_freq epochs if possible
% %     Data_of = [rootdir 'Analyse_data' filesep 'of_' sub filesep];
% %     if exist(Data_of,'dir')
% %        load([Data_of 'epoch_' RunCond]); % epoch
% %        epoch_of = epoch;
% %        clear epoch
% %     end
   
    %%%%======== load in sema_viol epochs 
    PPath.SaveData = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
    load([PPath.SaveData 'epoch_' RunCond]); % epoch
    
% %     %%%=== combining 2 datasets together if possible
% %     if exist(Data_of,'dir')
% %         epoch.trialinfo = [epoch.trialinfo; epoch_of.trialinfo];
% %         epoch.trial = [epoch.trial epoch_of.trial];
% %         epoch.time = [epoch.time epoch_of.time];
% %         clear epoch_of
% %     end
    
    %%%%===== normalize pd and remove trials that without photodies signal,
    %%%% no saving out, rawdata is rawdata
    %%% get pd label index
    pdi = find(strcmp(epoch.label,'MISC004'));
    rmtrl = [];
    for ttt = 1:length(epoch.trial)
        %%% remove trials that have no pd-004 signal
        if max(epoch.trial{ttt}(pdi,:)) < 0.005
            rmtrl = [rmtrl; ttt];
        end
        epoch.trial{ttt}((pdi),:) = zscore(epoch.trial{ttt}((pdi),:),0,2);
    end
    cfg = [];
    cfg.trials = setdiff(1:length(epoch.trial),rmtrl);
    epoch = ft_selectdata(cfg, epoch);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% 2. getting significant tagging response sensors %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% permutation test between epoch data and baseline data
    %%% selecting comparing epochs
    flickid = epoch.trialinfo(:,tag_col)==targid(eee);
    firstpassid = epoch.trialinfo(:,firstpass_col)==1;
    cfg = [];
    cfg.trials = find(flickid & firstpassid);
    data = ft_selectdata(cfg, epoch);
    load([PPath.SaveData 'epoch_BL_Cross']);
    
    %%%======= get fourier spctrm --- no time domain
    cfg              = [];
    cfg.output       = 'fourier';
    cfg.channel      = {'MEGGRAD','MISC004'};
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = tagfreq;
    cfg.pad          = 'nextpow2';
    fourier          = ft_freqanalysis(cfg,data);%=% tfr.powspctrm:3D data: trialnum*channelnum*freqpoint
    fourier_bl       = ft_freqanalysis(cfg,epoch_BL_Cross);
    
    %%%% statistical test of coherence
    %%%% https://mailman.science.ru.nl/pipermail/fieldtrip/2018-January/012000.html
    %%%% --- doesnot work! can not seperate the
    %%%% channelcmb in stat,it's huge:
    %%%% since 'cluster correction' is not suitable here in the coherence
    %%%% statistic,(https://mailman.science.ru.nl/pipermail/fieldtrip/2010-April/002830.html)
    %%%% and the stat cann't output by individual channel combination, I trick it as single channel combination each time manually
    load([rootdir 'OccipSens_full'])
    nchancom = length(OccipSens_full);
    stat_mask = zeros(nchancom,1);
    stat_p = zeros(nchancom,1);
    for i = 1:nchancom % run stat over sensors
        cfg            = [];
        cfg.channel    = {OccipSens_full{i},'MISC004'};
        fourier_tmp    = ft_selectdata(cfg, fourier);
        fourier_bl_tmp = ft_selectdata(cfg, fourier_bl);
        % do stat between tagging_period and baseline for each sensor
        cfg                  = [];
        cfg.parameter        = 'fourierspctrm';
        cfg.frequency        = tagfreq;
        cfg.statistic        = 'ft_statfun_indepsamplesZcoh';  %%%% take fourierspctrm as input, so no time domain information
        cfg.method           = 'montecarlo';
        cfg.tail             = 1; %% right sided, grp1 is bigger than grp2
        cfg.alpha            = SenSelectP;
        cfg.numrandomization = 5000;
        ntrl_1 = size(fourier.fourierspctrm,1);
        ntrl_2 = size(fourier_bl.fourierspctrm,1);
        design = zeros(1, ntrl_1 + ntrl_2);
        design(1,1:ntrl_1) = 1;
        design(1,(ntrl_1 +1):(ntrl_1 + ntrl_2))= 2;
        cfg.design = design;
        cfg.ivar   = 1;
        cfg.design = design;
        stat = ft_freqstatistics(cfg, fourier_tmp, fourier_bl_tmp);
        stat_mask(i) = stat.mask;
        stat_p(i,1)  = stat.prob;
    end
    TagCoh.Stat(:,1) = stat_mask;
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% 2. runing 2D coherence for different conditions %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TagCoh.freq = freqrange;
    TagCoh.time = epoch.time{1};
    
    %%% select data that are relevant to this task
    cfg = [];
    cfg.trials = ismember(epoch.trialinfo(:,cond_col),conds);
    epoch = ft_selectdata(cfg, epoch);
    
    %%%===== epoch and conditon loop
    for mmm = 1:length(EpochType)
        TrlNum = []; % get trl numbers for each condition
        trlidx = {}; % get trl index for each condition
        for ccc = 1:length(condnm)
            if ccc == 1
                m_trlid = epoch.trialinfo(:,tag_col) == targid(mmm); %trl idx for this EpochType
                idx = find(m_trlid);
            else
                idx = find(m_trlid & epoch.trialinfo(:,cond_col)==conds(ccc-1));
            end
            TrlNum = [TrlNum length(idx)];
            trlidx{ccc} = idx;
        end
        %%% get equal trial number across all conditions
        n_min = min(TrlNum(2:end));
        for ccc = 1:length(condnm)
            idx = trlidx{ccc};
            if ccc > 1 % data from conditions
               idx = idx(randperm(TrlNum(ccc)));
               idx = idx(1:n_min);
               trlidx{ccc} = idx;
            end
            cfg        = [];
            cfg.trials = idx;
            eval(['data_' num2str(ccc) '= ft_selectdata(cfg, epoch);']);
        end 
        %%% save out 
        eval(['TagCoh.' EpochType{mmm} '_TrlNum_raw(1,:) = TrlNum;']);
        for ccc = 1:length(condnm)
            eval(['TagCoh.' EpochType{mmm} '_TrlInfo{1,ccc} = data_' num2str(ccc) '.trialinfo;']);
            eval(['TagCoh.' EpochType{mmm} '_RT(1,ccc) = nanmean(data_' num2str(ccc) '.trialinfo(:,fixdur_col));']);
            eval(['TagCoh.' EpochType{mmm} '_TrlIdx_equ{1,ccc} = trlidx{ccc};']);
        end
        
        %%%%%%%%%%%=================== Coherence===============%%%%%%%%%%%%
        %%%% conditon loop
        for ccc = 1:length(condnm)
            %%%% just to get a dummy coherence struct
            cfg            = [];
            cfg.output     = 'fourier';
            cfg.channel    = {'MEGGRAD','MISC004'};
            cfg.method     = 'mtmconvol';
            cfg.taper      = 'hanning';
            cfg.foi        = freqrange;
            cfg.toi        = -0.5:0.5:0.5;
            cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
            cfg.keeptrials = 'yes';
            cfg.pad        = 'nextpow2';
            eval(['fourier = ft_freqanalysis(cfg, data_' num2str(ccc) ');']);
            %%% get coherence spctrm
            cfg            = [];
            cfg.method     = 'coh';
            cfg.channelcmb = {'MEGGRAD', 'MISC004'}; %% all channell combinations calculated together
            coh            = ft_connectivityanalysis(cfg, fourier); %% the raw coherence
            coh.label      = fourier.label(1:length(coh.labelcmb));
            coh.time       = epoch.time{1};
            
            %%% do the real coherence using hilbert complex
            coh.cohspctrm = zeros(length(coh.label),length(coh.freq),length(coh.time));
            for fff = 1:length(freqrange)
                %%%% get hilbert filert data
                cfg = [];
                cfg.channel    = {'MEGGRAD','MISC004'};
                cfg.bpfilter   = 'yes';
                cfg.bpfreq     = [freqrange(fff)-freq_halfwidth freqrange(fff)+freq_halfwidth];
                cfg.hilbert    = 'complex';
                cfg.keeptrials = 'yes';
                eval(['fltdata = ft_preprocessing(cfg,data_' num2str(ccc) ');']);
                for chan = 1:length(fltdata.label)-1
                    for ttt = 1:length(fltdata.trial)
                        sig1(:,ttt) = fltdata.trial{ttt}(chan,:); %time*trl
                        sig2(:,ttt) = fltdata.trial{ttt}(end,:);
                    end
                    spec1 = nanmean(sig1.*conj(sig1),2);%time*1
                    spec2 = nanmean(sig2.*conj(sig2),2);
                    specX = abs(nanmean(sig1.*conj(sig2),2)).^2;
                    coh.cohspctrm(chan,fff,:) = specX./(spec1.*spec2);%time*1
                end
            end
            
            %%% finally get the sig occipital tagging sensors
            if mmm == 1 && ccc == 1
                sigsens = OccipSens_full(stat_mask > 0);
                %%% get all the sig tagging sensors in the order of p value
                sigp   = stat_p(stat_mask > 0);
                [~, Isig] = sort(sigp,'ascend');
                TagCoh.SigTagSen_All{1,1} = sigsens(Isig);
                %%%%% get the strongest responder
                tmp = nanmean(nanmean(coh.cohspctrm(:,tagfreq_id,:),2),3);
                %%% set nan to negtive value, otherwise nan would be the
                %%% biggest value
                tmp(isnan(tmp)) = 0;
                TagCoh.SigTagSen_Occip(1,1) = {'No'};
                TagCoh.SigTagSen_Occip{1,2} = 0;
                sigchanid = 156; %pre-select a visual sensor
                [~, I] = sort(tmp,'descend');
                sigsen = [];
                %%%% first n strongest tagging sensors
                for kkk = 1:length(I)  %%% check the strongest n sensors
                    tmpsen = coh.label(I(kkk)); %%% the strongest tagging responder
                    if ~isempty(sigsens) && ismember(tmpsen,sigsens)
                        sigsen = [sigsen, tmpsen];
                        sigchanid(1,length(sigsen)) = I(kkk);
                        TagCoh.SigTagSen_Occip{1,1} = sigsen;
                        TagCoh.SigTagSen_Occip{1,2} = sigchanid;
                    end
                end
            end
            %%% finally get the sig occipital tagging sensors
            cohcoh = nanmean(coh.cohspctrm(sigchanid,:,:),1);
            eval(['TagCoh.' EpochType{mmm} '_Coh(:,:,1,ccc) = squeeze(cohcoh);']);
        end
    end
    save([PPath.FigPath DataSets{ddd} '_' num2str(sid)],'TagCoh','-v7.3');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% 3. Plotting figures %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figtitle = [EpochType{eee} '_s' num2str(sid) ];
    if isempty(sigsens)
        figtitle = ['NoSig_' DataSets{ddd} '_' EpochType{eee} '_' RunCond '_s' num2str(sid) '_coh'];
    end
    h = figure('Name',figtitle,'color',[1 1 1]);
    for mmm = 1:length(EpochType)
        %%%% coh plot
        eval(['tmpdata = nanmean(TagCoh.' EpochType{mmm} '_Coh(:,:,1,:),3);']);%%% freq*time*sub*cond
        tmpcc = [1 2 3];
        for cc = 1:length(condnm)+1
            if ismember(cc,tmpcc)
                condid = find(tmpcc==cc);
                eval(['rt = TagCoh.' EpochType{mmm} '_RT(:,condid)./1000;']);
                rt = nanmean(rt);
                loc = tmpcc(condid)+4*(mmm-1);
                plotdata = squeeze(tmpdata(:,:,:,condid));
                subtitle = [EpochType{mmm} 'et  ' TagCoh.CondName{condid}];
                subplot(3,4,loc);
                pcolor(TagCoh.time, TagCoh.freq, plotdata); colorbar;
            else
                eval(['rt = TagCoh.' EpochType{mmm} '_RT(:,1)./1000;']); %% for the difference figure, use the RT from 'all' cond
                rt = nanmean(rt);
                loc = cc+4*(mmm-1);
                if cc == 4
                    idx1 = 2; idx2 = 3;
                end
                plotdata = squeeze(tmpdata(:,:,:,idx1)-tmpdata(:,:,:,idx2));
                subtitle = [EpochType{mmm} 'et Coh-contrast'];
                subplot(3,4,loc);
                pcolor(TagCoh.time, TagCoh.freq, plotdata); colorbar;
            end
            colormap jet; shading interp; % colorbar;
            hold on;
            plot([-0.5 0.5],[60 60],'-.k','LineWidth',2)
            plot([-rt -rt],[1 100],'-.k','LineWidth',2)
            plot([0 0],[1 100],'-.k','LineWidth',2)
            plot([rt rt],[1 100],'-.k','LineWidth',2)
            set(gca,'XTick',-0.4:0.2:0.4);
            set(gca,'XTickLabel',{'-0.4','-0.2','FixOn','0.2','0.4'},'FontWeight','bold','FontSize',10)
            title(subtitle,'FontWeight','bold','FontSize',14);
          end
    end
    ScSz = [0 0 1300 800];
    set(gcf,'Position',ScSz);
    saveas(h,[PPath.FigPath figtitle]);
    close all
    
    
    
else  %%%server == 0
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run locally get the whole structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% basic parameters
    %%% combine data together into one big structure
    load([PPath.FigPath DataSets{ddd} '_1']);
    TagCoh_all = TagCoh;
    nsub = length(TagCoh.subs);
    TagCoh_all.SigTagSen_Occip(1:nsub,1) = repmat({'No'},nsub,1);
    TagCoh_all.SigTagSen_Occip(1:nsub,2) = repmat({0},nsub,1);
    for ps = 1:length(subjects)
        load([PPath.FigPath DataSets{ddd} '_' num2str(ps)]);
        for fd = 1:length(EpochType)
            if ps == 1
                TagCoh_all.Stat(:,ps) = [];
                eval(['TagCoh_all.' EpochType{fd} '_TrlNum_raw(ps,:)= [];']);
                eval(['TagCoh_all.' EpochType{fd} '_RT(ps,:)= [];']);
                eval(['TagCoh_all.' EpochType{fd} '_TrlInfo(ps,:)= [];']);
                TagCoh_all.SigTagSen_All{ps,1} = []; %%% all sig tagging occipital sensors
                eval(['TagCoh_all.' EpochType{fd} '_Coh(:,:,ps,:)= [];']);
            end
            TagCoh_all.Stat(:,ps) = TagCoh.Stat(:,1);
            eval(['TagCoh_all.' EpochType{fd} '_TrlNum_raw(ps,:)= TagCoh.' EpochType{fd} '_TrlNum_raw(1,:);']);
            eval(['TagCoh_all.' EpochType{fd} '_TrlIdx_equ(ps,:)= TagCoh.' EpochType{fd} '_TrlIdx_equ(1,:);']);
            eval(['TagCoh_all.' EpochType{fd} '_RT(ps,:)= TagCoh.' EpochType{fd} '_RT(1,:);']);
            eval(['TagCoh_all.' EpochType{fd} '_TrlInfo(ps,:)= TagCoh.' EpochType{fd} '_TrlInfo(1,:);']);
            TagCoh_all.SigTagSen_All{ps,1} = TagCoh.SigTagSen_All{1,1}; %%% all sig tagging occipital sensors
            TagCoh_all.SigTagSen_Occip(ps,[1 2]) = TagCoh.SigTagSen_Occip(1,[1 2]);
            eval(['TagCoh_all.' EpochType{fd} '_Coh(:,:,ps,:)= TagCoh.' EpochType{fd} '_Coh(:,:,1,:);']);
        end
    end
    clear TagCoh
    TagCoh = TagCoh_all;
    clear TagCoh_all
    TagCoh.SenSelectP = [num2str(SenSelectP) '_one tail'];
    lesstrl_sub = find(sum(TagCoh.PreTarg_TrlNum_raw(:,[2 3])<40,2)>0); %%% including non-sig subs
    nosigsub = find(strcmp(TagCoh.SigTagSen_Occip(:,1),'No'));
    TagCoh.LessTrlSubID = lesstrl_sub;
    TagCoh.SigSubID = transpose(setdiff(1:length(TagCoh.subs),[nosigsub; lesstrl_sub]));
    eval(['sent_version = ExpInfo.SentenceVersion.' DS ';']);
    sent_version = sent_version(~isnan(sent_version)); % remove empty rows
    TagCoh.SigSub_SentVersion = sent_version(TagCoh.SigSubID); 
    tabulate(TagCoh.SigSub_SentVersion) % V1_16sub, V2_14subs
    TagCoh.MeanPretargTrlNum4SigSub = nanmean(TagCoh.PreTarg_TrlNum_raw(TagCoh.SigSubID,:));
    save([PPath.FigPath 'TagCoh.mat'],'TagCoh','-v7.3');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3. group plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plotting-coh
    sigsubid = TagCoh.SigSubID;
    figtitle = 'Group_Coh';
    h = figure('Name',figtitle,'color',[1 1 1]);
    ScSz = [1 1 900 800];%get(groot, 'Screensize' );
    set(h,'Position', ScSz,'color',[1 1 1],'MenuBar','figure');
    pos_grid =  [ScSz(1)+10 ScSz(2)+10  ScSz(3)-10  ScSz(4)-10];
    nrow = length(condnm)+1;
    ncol = length(EpochType);
    gridPos = DivideScreen(nrow,ncol,ScSz,50,70,pos_grid);
    for mmm = 1:length(EpochType)
        %%%% coh plot
        eval(['tmpall = TagCoh.' EpochType{mmm} '_Coh(:,:,sigsubid,:);']);%%% freq*time*sub*cond
        tmpdata = squeeze(nanmean(tmpall,3));%%% freq*time*sub*cond
        tmpcc = [1 2 3];
        for cc = 1:nrow
            if ismember(cc,tmpcc)
                condid = find(tmpcc==cc);
                eval(['rt = TagCoh.' EpochType{mmm} '_RT(sigsubid,condid)./1000;']);
                rt = nanmean(rt);
                loc = mmm + ncol*(tmpcc(condid)-1);
                plotdata = squeeze(tmpdata(:,:,condid));
                subtitle = [EpochType{mmm} 'et  ' TagCoh.CondName{condid}];
                h = axes('Position',gridPos{loc});
                pcolor(TagCoh.time, TagCoh.freq, plotdata);axcopy;
                caxis([0 0.04]);
                xlim([-0.4 0.4])
            else
                loc = mmm + ncol*(cc-1);
                if cc == 4
                    idx1 = 2; idx2 = 3;
                else
                    idx1 = 4; idx2 = 5;
                end
                eval(['rt = TagCoh.' EpochType{mmm} '_RT(sigsubid,[idx1 idx2])./1000;']); %% for the difference figure, use the RT from 'all' cond
                rt = nanmean(rt,1);
                rt = min(rt);
                plotdata = squeeze(nanmean(tmpall(:,:,:,idx1)-tmpall(:,:,:,idx2),3));
                subtitle = [EpochType{mmm} 'et Incong-Cong'];
                h = axes('Position',gridPos{loc});
                pcolor(TagCoh.time, TagCoh.freq, plotdata);axcopy;
                caxis([-0.005 0.005]);
                xlim([-0.4 0.4])
            end
            colormap jet; shading interp;
            hold on;
            plot([-0.5 0.5],[60 60],'-.k','LineWidth',2)
            plot([-rt -rt],[1 100],'-.k','LineWidth',2)
            plot([0 0],[1 100],'-.k','LineWidth',2)
            plot([rt rt],[1 100],'-.k','LineWidth',2)
            set(gca,'XTick',-0.4:0.2:0.4);
            set(gca,'XTickLabel',{'-0.4','-0.2',['Fix' RunCond(4:end)],'0.2','0.4'},'FontWeight','bold','FontSize',12)
            title(subtitle,'FontWeight','bold','FontSize',16);
            if mmm == 3 && (cc == 1 || cc ==4 )
                colorbar('FontWeight','bold','FontSize',10)
            end
        end
    end
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'svg');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Get total topograph of RFT %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% get the occipital sensors
    load([rootdir filesep 'Analyse_data' filesep DataSets{ddd} '_20210625_b4d4' filesep 'epoch_BL_Cross.mat']) %%% just load random data that has all the labels
    cfg = [];
    erf = ft_timelockanalysis(cfg, epoch_BL_Cross);
    
    %%% get labels
    figtitle = [EpochType{eee} '-RFT-sens'];
    SigSen = TagCoh.SigTagSen_Occip(sigsubid,1);
    tmp = zeros(size(erf.label));
    for bbb = 1:length(SigSen)
        for kkk = 1:length(SigSen{bbb})
            idx = find(strcmp(erf.label,[SigSen{bbb}{kkk}(1:end-1) '1']));
            tmp(idx) = tmp(idx)+1;
        end
    end
    erf.avg = repmat(tmp,1,size(erf.avg,2));
    h = figure('color', [1 1 1],'name',figtitle);
    cfg =[];
    cfg.layout = 'neuromag306mag.lay'; %'neuromag306cmb.lay';
    cfg.comment = ' ';
    ft_topoplotER(cfg,erf);
    colormap jet; colorbar('FontWeight','bold','FontSize',10);
    caxis([-max(tmp) max(tmp)])
    title('RFT reponse sensors','FontWeight','bold','FontSize',16);
    text(0.5,0.5,'No. of sensors','FontWeight','bold','FontSize',14,'Rotation',270)
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'svg');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Just do the simple ttest to find the effect %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq60 = dsearchn(TagCoh.freq',60);
    sigsubid = TagCoh.SigSubID;
    nsub = length(sigsubid);
    tmptime = transpose(TagCoh.time);
    t0 = dsearchn(tmptime,0);
    for mm = 1:length(EpochType)
        eval(['rt = TagCoh.' EpochType{mm} '_TrlInfo;']);
        eval(['cohdata = TagCoh.' EpochType{mm} '_Coh;']);
        tdata = zeros(nsub,2);
        trl_minrt = zeros(nsub,1);
        for ss = 1:nsub
            sid = sigsubid(ss);
            %%% get min rt of all the trials in 2 conditions
            trl_minrt(ss,1) = min([rt{sid,2}(:,fixdur_col); rt{sid,3}(:,fixdur_col)])/1000;
            if strfind(RunCond,'On') %% epoch aligned with fixaiton_on, select timewindow from zero to this fixaiton_duration
                ttt = dsearchn(tmptime,trl_minrt(ss,1));
                tdata(ss,:) = squeeze(nanmean(cohdata(freq60,t0:ttt,sid,[2 3]),2));
            elseif strfind(RunCond,'Off') %% epoch aligned with fixaiton_off, select timewindow from this fixaiton_duration to zero
                ttt = dsearchn(tmptime,-trl_minrt(ss,1));
                %                     ttt = dsearchn(tmptime,-0.2);
                tdata(ss,:) = squeeze(nanmean(cohdata(freq60,ttt:t0,sid,[2 3]),2));
            end
        end
        [~,p,~,stats] = ttest(tdata(:,1), tdata(:,2));
        stats.p = p; disp([EpochType{mm} '_Ttest = ' num2str(p)])
        stats.data4test = tdata;
        stats.trl_minrt = trl_minrt;
        eval(['TagCoh.' EpochType{mm} '_Ttest = stats;'])
    end
    
    %% %%======= violin plot of the word-freq effect on coherence
    colmat = [0 114 189;217 83 25]./255;
    figtitle = 'Ttest_Freq_violin';
    nsigsub = size(TagCoh.PreTarg_Ttest.data4test,1);
    group = [cellstr(repmat('Incong',nsigsub,1)); cellstr(repmat('Cong',nsigsub,1))];
    grouporder={'Incong','Cong'};
    EpochType = {'PreTarg','Targ'};
    subtitles = {'Pre-target interval','Target interval'};
    figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 450 250]);
    for mmm = 1:length(EpochType)
        eval(['vdata = [TagCoh.' EpochType{mmm} '_Ttest.data4test(:,1); TagCoh.' EpochType{mmm} '_Ttest.data4test(:,2)];']);
        h = subplot(1,length(EpochType),mmm);
        vp = violinplot(vdata, group,'GroupOrder',grouporder);
        vp(1).ViolinColor = colmat(1,:);
        vp(2).ViolinColor = colmat(2,:);
        vp(1).ShowMean = 1; vp(2).ShowMean = 1;
        vp(1,1).MedianPlot.Visible = 'off';
        vp(1,2).MedianPlot.Visible = 'off';
        vp(1,1).MeanPlot.LineWidth = 1.5;
        vp(1,2).MeanPlot.LineWidth = 1.5;
        vp(1,1).BoxWidth = 0.01; vp(1,2).BoxWidth = 0.01;
        ylabel('Coherence at 60 Hz (r^2)','FontSize',7,'FontWeight','normal','FontName','Arial');
        xlabel('Semantic congruency of target','FontSize',7,'FontWeight','normal');
        set(gca,'FontSize',7,'FontWeight','normal','FontName','Arial');
        set(gca,'box','off','LineWidth',1)
        %%% plot the line linking each subject
        hold on;
        x1 = vp(1,1).ScatterPlot.XData;
        y1 = vp(1,1).ScatterPlot.YData;
        x2 = vp(1,2).ScatterPlot.XData;
        y2 = vp(1,2).ScatterPlot.YData;
        plot([x1; x2],[y1; y2],'Color',[.8 .8 .8],'linewidth',0.5)
        %%% add stat
        if mmm == 1
            plot([1 2],[0.065 0.065],'k','LineWidth',1)
            text(1.5,0.061,'*','FontWeight','normal','FontSize',14,'FontName','Arial')
        else
            plot([1 2],[0.065 0.065],'k','LineWidth',1)
            text(1.5,0.063,'n.s.','FontWeight','normal','FontSize',7,'FontName','Arial')
        end
        title(subtitles{mmm},'FontSize',7,'FontWeight','normal','FontName','Arial');
    end
    set(gcf, 'renderer', 'painters')
    saveas(h,[PPath.FigPath figtitle]);
    saveas(h,[PPath.FigPath figtitle],'svg');
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% coherence latency statistics ---- JakeNife %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigsubid = TagCoh.SigSubID;
    %%% plot RFT curve for all conditions
    for i = 1:length(EpochType)
        eval(['tempdata = TagCoh.' EpochType{i} '_Coh;']);
        tempdata = squeeze(nanmean(tempdata(tagfreq_id,:,sigsubid,:),1)); %%% time*sub*cond
        eval(['TagCoh.' EpochType{i} '_Coh_TimSubCond = tempdata;']);
    end
    for ee = 1:length(EpochType)
        eval(['cohdata_all = TagCoh.' EpochType{ee} '_Coh_TimSubCond(:,:,[2 3]);']); %% time*sub*cond
        eval(['RT_all = TagCoh.' EpochType{ee} '_RT(sigsubid,[2 3])./1000;']); %% sub*cond
        CM_tim = []; %% center of mass
        HM_tim = []; %% when reach the half maximum amplitude (distant amplitude between max and min coh)
       
        %%% doing jackknife loop
        allsigsub = 1:length(sigsubid);
        exclu_sub = [allsigsub 0]; %%[jackknife all]
        for sid = 1:length(exclu_sub)
            subid_jn = setdiff(allsigsub,exclu_sub(sid));
            %%% data after jackknife
            cohdata = cohdata_all(:,subid_jn,:);
            RT = RT_all(subid_jn,:);
            meanRFT = squeeze(nanmean(cohdata,2)); %%% time*cond
            seRFT = squeeze(nanstd(cohdata,0,2)./sqrt(size(cohdata,2)));
            meanRT = mean(RT);
            %%% get the time window
            rt = min(meanRT);
            etp = find(TagCoh.time < rt); %%% rt window aligned with zero--saccadeonset
            zero_tp = nearest(TagCoh.time,0);
            tprange = zero_tp:etp(end);
            timrange = TagCoh.time(tprange);
            
            for cc = 1:2 %%[Incong Cong]
                %%%=== center of mass
                tmptp = meanRFT(tprange,cc); %%% data in the timewindow: time*sub*cond
                timevct = [1:length(timrange)]';
                cm = round(sum(timevct.*tmptp)/sum(tmptp));
                CM_tim(sid,cc) = timrange(cm);
                
                %%%=== half amplitude latency
                [mintmp, minidx] = min(tmptp);
                %%% make sure that the maxtmp is later than the mintmp
                tmptp = tmptp(minidx:end);
                [maxtmp, maxidx] = max(tmptp);
                tmpidx = find(tmptp>((maxtmp-mintmp)/2+mintmp));
                if ~isempty(tmpidx)
                    HM_tim(sid,cc) = timrange(tmpidx(1)+minidx);
                else
                    HM_tim(sid,cc) = nan;
                end
            end
        end
        eval(['JkNf.' EpochType{ee} '_HM_tim = HM_tim ;']);
        eval(['JkNf.' EpochType{ee} '_CM_tim = CM_tim ;']);
    end
    
    %%% using jacknife to do the statistics
    tp_tpye = {'CM_tim';'HM_tim'};
    N = length(sigsubid);
    for tp = 1:length(tp_tpye)
        for ee = 1:length(EpochType)
            tj = [];
            eval(['tmptim = JkNf.' EpochType{ee} '_' tp_tpye{tp} ';']);
            Di = tmptim(1:end-1,1)-tmptim(1:end-1,2);
            J = nansum(Di)/N;
            SD = sqrt([(N-1)/N]*[nansum((Di-J).^2)]);
            D = tmptim(end,1) - tmptim(end,2);
            tj = D/SD;
            eval(['JkNf.'  EpochType{ee} '_' tp_tpye{tp} '_tvalue = tj;']);
        end
    end
    % check if there's any sig, need to loopup t-value table
    JkNf.PreTarg_HM_tim_pvalue = '';
    TagCoh.Jackknife_latency = JkNf;
    save([PPath.FigPath 'TagCoh.mat'], 'TagCoh')
    
    %% === plot the group level curve with sig, only for HM (half-maximum)
    %%% PreTarg_HM_tim_tvalue=2.17,threshold for p=0.05 is
    %%% t=2.045(two-tail,n=29)
    curve_xmax = 0.4; %%% the x-axis range
    etp = find(TagCoh.time < curve_xmax); %%% rt window aligned with zero--saccadeonset
    zero_tp = nearest(TagCoh.time,0);
    tprange = zero_tp:etp(end);
    timrange = TagCoh.time(tprange);
    figtitles = {'Pre-target','Target'};
    EpochType = {'PreTarg','Targ'};
    for mmm = 1:length(EpochType)
        eval(['cohdata_all = TagCoh.' EpochType{mmm} '_Coh_TimSubCond(tprange,:,[2 3]);']); %% time*sub*cond
        meanRFT = squeeze(mean(cohdata_all,2));
        seRFT = squeeze(nanstd(cohdata_all,0,2)./sqrt(length(TagCoh.SigSubID)));
        eval(['HM = TagCoh.Jackknife_latency.' EpochType{mmm} '_HM_tim(end,:,1);']);
        colmat = [0 114 189;217 83 25]./255;
        %%%% figure
        figtitle = [figtitles{mmm} ' coherence onset latency'];
        h = figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 220 240]);
        a = shadedErrorBar(timrange,meanRFT(:,1),seRFT(:,1),{'color',colmat(1,:)},0.8);hold on;
        b = shadedErrorBar(timrange,meanRFT(:,2),seRFT(:,2),{'color',colmat(2,:)},0.9);
        a.mainLine.LineWidth = 1;
        b.mainLine.LineWidth = 1;
        legendflex([a.mainLine,b.mainLine],{'Incong';'Cong'},'anchor', {'ne','ne'}, 'buffer', [0 0],'Fontsize',7,'xscale',1,'box','off');
        set(gca,'FontSize',7,'FontWeight','normal');
        set(gca,'box','off','LineWidth',1)
        set(gca,'XTick',0:0.05:0.4);
        set(gca,'XTickLabel',{'FixOn','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4'},'FontSize',7)
        set(gca,'YTick',0.01:0.01:0.06);
        title(figtitle,'FontSize',7,'FontWeight','normal')
        xlabel('Time (s)','FontSize',7,'FontWeight','normal')
        ylabel('Coherence at 60 Hz (r^2)','FontSize',7,'FontWeight','normal','FontName','Arial');
        %%% plot the latency lines
        hold on;
        timeliney = get(gca,'ylim');
        timelinex = repmat(HM(1),1,2);
        plot(timelinex', timeliney','--','color',colmat(1,:),'LineWidth',1);
        timelinex = repmat(HM(2),1,2);
        plot(timelinex', timeliney','--','color',colmat(2,:),'LineWidth',1);
        set(gcf, 'renderer', 'painters')
        saveas(h,[PPath.FigPath figtitle]);
        saveas(h,[PPath.FigPath figtitle],'svg');
    end
    
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% coherence & eye movement correlation %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([rootdir filesep 'Results' filesep 'Behavioral' filesep 'CmbDatasets_BehaData']);
    
    AllFix = nanmean(BehaData.FirstFix.Tag(sigsubid,:),2);
    DiffFix_Targ = BehaData.FirstFix.Tag(sigsubid,:);
    DiffFix_2 = DiffFix_Targ(:,1) - DiffFix_Targ(:,2); %% diff fixation for target: Incong-Cong
    figtitle = [DataSets{ddd} '-Corr'];
    h = figure('Name',figtitle,'color',[1 1 1]);
    preTP = TagCoh.PreTarg_Ttest.data4test; %%% coherence of pretarget
    DiffTP = preTP(:,1) - preTP(:,2); %%% Incong-Cong
    %%%% PreTarg-DiffRFT & Targ Diff Fixation Duration
    subplot(1,3,1)
    scatter(DiffTP,DiffFix_2,50,'r','filled')
    ylabel('Target First Fixation Incong-Cong (ms)','FontWeight','bold','FontSize',16)
    xlabel('PreTarget coherence Incong-Cong (r)','FontWeight','bold','FontSize',16)
    set(gca,'FontSize',12,'FontWeight','bold');
    set(gca,'box','off','LineWidth',2)
    [coef, pval] = corr([DiffTP  DiffFix_2],'type','Spearman');
    title(['r = ' num2str(round(1000*coef(1,2))/1000) ', p = ' num2str(round(1000*pval(1,2))/1000)], 'FontWeight','bold','FontSize',16);
    xlim([-0.02 0.03])
    ylim([-30 70])
    
    %%% PreTarg-DiffRFT &  Targ Fixation Duration
    subplot(1,3,2)
    scatter(DiffTP,AllFix,50,'r','filled')
    ylabel('Target First Fixation (ms)','FontWeight','bold','FontSize',16)
    xlabel('PreTarget coherence Incong-Cong (r)','FontWeight','bold','FontSize',16)
    set(gca,'FontSize',12,'FontWeight','bold');
    set(gca,'box','off','LineWidth',2)
    [coef, pval] = corr([DiffTP AllFix],'type','Spearman');
    title(['r = ' num2str(round(1000*coef(1,2))/1000) ', p = ' num2str(round(1000*pval(1,2))/1000)], 'FontWeight','bold','FontSize',16);
    xlim([-0.02 0.03])
    ylim([100 450])
    
    %%% PreTarg-DiffRFT & total RT
    subplot(1,3,3)
    sigsubid = TagCoh.SigSubID;
    SentDur = BehaData.SentDur_perwrd(sigsubid);
    scatter(DiffTP,SentDur,50,'r','filled')
    ylabel('Averaged fxiation in sentence (ms)','FontWeight','bold','FontSize',12)
    xlabel('PreTarget coherence Incong-Cong (r)','FontWeight','bold','FontSize',12)
    set(gca,'FontSize',12,'FontWeight','bold');
    set(gca,'box','off','LineWidth',1.5)
    [coef, pval] = corr([DiffTP  SentDur],'type','Spearman');
    title(['r = ' num2str(round(1000*coef(1,2))/1000) ', p = ' num2str(round(1000*pval(1,2))/1000)], 'FontWeight','bold','FontSize',16);
    xlim([-0.02 0.03])
    ylim([50 360])
    saveas(h,[PPath.FigPath figtitle]);
    save([PPath.FigPath DataSets{ddd} '.mat'], 'TagCoh')
    
    
    SentDur_nopre = BehaData.SentDur_perwrd_nopre(sigsubid);
    [coef, pval] = corr([DiffTP  SentDur_nopre],'type','Spearman') %%[-0.3258 0.1046]
    
    
end






