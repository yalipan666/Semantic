% The user will need to download and install the fieldtrip matlab toolbox (http://www.fieldtriptoolbox.org/).
% for those who have no MRI, use the standard template instead

% Source modeling for target fixation-related fields using LCMV
% Source modeling for pre-target coherence using DICS

%% setting path
rootdir = '/rds/projects/2018/jenseno-reading/Semantic/';
addpath(['/rds/projects/2018/jenseno-reading' filesep 'fieldtrip-20200220' filesep]);
ft_defaults
PPath.hdm = [rootdir 'Results' filesep 'Source_headmodel' filesep];
PPath.MRIPath = [rootdir 'RawData' filesep 'MRI_data' filesep];
PPath.MEGPath = [rootdir 'RawData' filesep 'MEG_data' filesep];
load([rootdir 'Analyse_data' filesep 'ExpInfo.mat']);
ddd = 1; %index for the dataset
DS = ExpInfo.DSName{ddd}; 
PPath.FigPath = [rootdir 'Results' filesep DS filesep 'Source' filesep];
%%% get subinformation
eval(['subjects = ExpInfo.subjects.' DS ';']);
firstpass_col = find(strcmp(ExpInfo.EventHdr,'FirstPassFix'));
cond_col = find(strcmp(ExpInfo.EventHdr,'SentenceCondition'));
conds = ExpInfo.CondID{ddd};
tag_col = find(strcmp(ExpInfo.EventHdr,'loc2targ'));
tagsigsub = ExpInfo.TagSigSubID.sv;

% % %% prepare headmodel
% % % copy shared headmodel between studies
% % for p = 1:length(ExpInfo.MRIcode)
% %     if ~strcmp(ExpInfo.MRIcode{p}, 'nan')
% %         tmp = find(strcmp(ExpInfo.MRIcode{p},SubInfo.MEGcode_MRIcode(:,2)));
% %         if ~isempty(tmp)
% %             display(['megsub:' num2str(p)])
% %             load(['Z:\Lexical\Results\source_modeling\HeadModel\Hdm_MRIaligned_' SubInfo.MEGcode_MRIcode{tmp,1} '.mat']);
% %             save([PPath.hdm filesep 'Hdm_MRIaligned_' ExpInfo.subjects.sv{p,1}],'hdm','mri_aligned')
% %         end
% %     end
% % end
% % % only hdm of sub7 is shared between projects
% % for sss = 1:length(subjects)
% %     disp(['******* analyzing sub-' num2str(sss) '******'])
% %     sub = subjects{sss};
% %     % make sure the same number of subjects
% %     if length(subjects) ~= length(ExpInfo.MRIcode)
% %         error('number of subjects in MEG data and MRI data does not match!')
% %     end
% %         
% %     %%%%%===== read MRI
% %     mri_code = ExpInfo.MRIcode{sss,1};
% %     if strcmp(mri_code, 'nan') % no native mri, use standard mri instead
% %         [~, ftdir] = ft_version;
% %         templatedir = fullfile(ftdir, 'template', 'headmodel');
% %         mri = ft_read_mri(fullfile(templatedir, 'standard_mri.mat'));
% %     else
% %         niifile = [PPath.MRIPath mri_code '_nifti' filesep ExpInfo.MRIcode{sss,1} '_nifti' filesep 'T1_vol_v1_5.nii'];
% %         if ~exist(niifile,'file')
% %             niifile(end-4) = '9';
% %             if  ~exist(niifile,'file')
% %                 niifile(end-4) = '5';
% %                 niifile = [niifile '.gz'];
% %                 files = gunzip(niifile);
% %                 niifile = niifile{1};
% %             end
% %         end
% %         mri = ft_read_mri(niifile);
% %     end
% %     %%%%%===== automatic alignment align MRI and polhemus
% %     %%% read headshape
% %     %%% plot to check the alignment between headshape and mri
% %     megfile = [PPath.MEGPath sub filesep sub(3:8) filesep sub(end-3:end) '-1.fif'];
% %     meg_headshape = ft_read_headshape(megfile);
% %     
% %     %%% realign MEG head and mri ---enter 'n'
% %     cfg = [];
% %     cfg.method                = 'headshape';
% %     cfg.headshape.headshape   = meg_headshape;
% %     cfg.coordsys              = 'neuromag';
% %     cfg.headshape.interactive = 'yes'; %% 'yes'
% %     cfg.viewresult            = 'no'; %% 'yes'
% %     cfg.headshape.icp         = 'yes';
% %     mri_aligned = ft_volumerealign(cfg, mri);
% %     mri_aligned.coordsys = 'neuromag';
% %     clear mri
% %     %%% enter 'n'
% % % % %     %%%	对准 MRI 和 MEG 头坐标
% % % % % 	• 对话框中的坐标和参数是控制MRI头的运动的
% % % % % 	• rotate：    +       -    [-18 0 0]  
% % % % %           ○ x：向后仰  向前仰            
% % % % %           ○ y：向左下转  向右下转 (不倒翁)
% % % % %           ○ z：向右后旋  向左后旋 (猫头鹰)
% % % % % 	• translate：平移，将3D坐标系旋转到2维视图，当前的坐标轴就表示平移的方向和大小   [-2 0 40]
% % % % % 	• --- apply --- quit
% % 
% %     %%%%%===== volume segmentation
% %     %%% Each of the voxels of the anatomical MRI is assigned to a tissue class, this procedure is termed segmentation.
% %     cfg          = [];
% %     segmentedmri = ft_volumesegment(cfg, mri_aligned);
% %     
% %     %%% check segmentation!
% %     % add anatomical information to the segmentation
% %     segmentedmri.transform = mri_aligned.transform;
% %     segmentedmri.anatomy   = mri_aligned.anatomy;
% %     
% %     % % %             %%% plot to check
% %     % % %             cfg              = [];
% %     % % %             cfg.funparameter = 'gray';
% %     % % %             ft_sourceplot(cfg, segmentedmri);
% %     
% %     %%%%%===== prepare the headmodel
% %     cfg        = [];
% %     cfg.method = 'singleshell';
% %     hdm        = ft_prepare_headmodel(cfg, segmentedmri);
% %     clear segmentedmri
% %     
% %     %%% saveout the hdm and mri_aligned
% %     save([PPath.hdm filesep 'Hdm_MRIaligned_' sub],'hdm','mri_aligned')
% % end


%% source modeling for the coherence results using DICS
%%% http://www.fieldtriptoolbox.org/tutorial/beamformingextended/#introduction
%%%% basic settingup
fTarget = 60;
mini_rt = ones(size(subjects));

%%% setup outputs
nsub = length(tagsigsub);
Source.All = cell(nsub,1);
Source.BL = cell(nsub,1);
Source.Cond1 = cell(nsub,1);
Source.Cond2 = cell(nsub,1);

%%% loading exist hdm and mri_aligned
HdmMRI = 1;

%% %================= subject loop
tic
for sid = 1:nsub
    sss = tagsigsub(sid);
    sub = subjects{sss};
    disp(['******* analyzing sub-' num2str(sss) '******'])
    load([rootdir 'Results' filesep 'Source_headmodel' filesep 'Hdm_MRIaligned_' sub]); %%hdm,mri_aligned
    
    %%%%%===== loading data
    PPath.SubPath = [rootdir 'Analyse_data' filesep DS '_' sub filesep];
    load([PPath.SubPath 'epoch_WrdOn']); %epoch
    load([PPath.SubPath 'epoch_BL_Cross']); %epoch_BL_Cross
    %%% remove trials that without photodiode signal
    pdi = find(strcmp(epoch.label,'MISC004'));
    rmtrl = [];
    for ttt = 1:length(epoch.trial)
        %%% remove trials that have no pd-004 signal
        if max(epoch.trial{ttt}((pdi(1)),:)) < 0.005
            rmtrl = [rmtrl; ttt];
        end
        epoch.trial{ttt}((pdi),:) = zscore(epoch.trial{ttt}((pdi),:),0,2);
    end
    %%% get first fixation trials
    validduration = find(epoch.trialinfo(:,firstpass_col)==1 & epoch.trialinfo(:,tag_col)==-1); %% only pre-targets
    cfg        = [];
    cfg.trials = setdiff(validduration,rmtrl);
    cfg.toilim = [0 mini_rt(sss)];
    epoch      = ft_redefinetrial(cfg, epoch);
    %%% get the hdr.grad
    megfile = [PPath.MEGPath sub filesep sub(3:8) filesep sub(end-3:end) '-1.fif'];
    cfg           = [];
    cfg.dataset   = megfile;
    hdr           = ft_read_header(cfg.dataset);
    epoch.grad    = hdr.grad;
    epoch_BL_Cross.grad= hdr.grad;
    %%% select low and high conditions
    cfg = [];
    cfg.trials = find(epoch.trialinfo(:,cond_col)==conds(1));
    epoch_c1 = ft_redefinetrial(cfg, epoch);
    cfg.trials = find(epoch.trialinfo(:,cond_col)==conds(2));
    epoch_c2 = ft_redefinetrial(cfg, epoch);
    
    %%%%%===== cross spectral density
    cfg                 = [];
    cfg.output          = 'powandcsd';
    cfg.method          = 'mtmfft';
    cfg.taper           = 'hanning';%dpss
    %%% for 1s pre-target, using 2Hz smoothing
    cfg.foi        = fTarget;
    cfg.keeptrials = 'yes';
    cfg.channel    = {'MEGGRAD' 'MISC004'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD'; 'MEGGRAD' 'MISC004'};
    freq_all       = ft_freqanalysis(cfg, epoch);
    freq_bl        = ft_freqanalysis(cfg, epoch_BL_Cross);
    freq_c1        = ft_freqanalysis(cfg, epoch_c1);
    freq_c2        = ft_freqanalysis(cfg, epoch_c2);
    
    %%%%%======= Computing the sourcemodel
    %%% MNI space
    %%% this returns the location where FieldTrip is installed
    [~, ftdir] = ft_version;
    %%% and this is where the template source models are
    templatedir = fullfile(ftdir, 'template', 'sourcemodel');
    template = load(fullfile(templatedir, 'standard_sourcemodel3d5mm')); % 5mm spacing grid
    %%% inverse-warp the template grid to subject specific coordinates
    cfg           = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template.sourcemodel;
    cfg.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri       = mri_aligned;
    sourcemodel   = ft_prepare_sourcemodel(cfg);
    
    %%% native space
    cfg                 = [];
    cfg.grid.resoltuion = 5;
    cfg.mri             = mri_aligned;
    sourcemodel_sub     = ft_prepare_sourcemodel(cfg);
    
    %%%%%===== Source analysis --> dics
    %%% MNI space
    cfg             = [];
    cfg.method      = 'dics';
    cfg.refchan     = 'MISC004'; %%% used for coherence analysis!!!
    cfg.frequency   = fTarget;
    cfg.headmodel   = hdm;
    cfg.sourcemodel = sourcemodel;
    source_all      = ft_sourceanalysis(cfg, freq_all);
    source_bl       = ft_sourceanalysis(cfg, freq_bl);
    source_c1      = ft_sourceanalysis(cfg, freq_c1);
    source_c2     = ft_sourceanalysis(cfg, freq_c2);
    source_all.pos  = template.sourcemodel.pos;
    source_all.dim  = template.sourcemodel.dim;
    source_all.coh  = source_all.avg.coh;
    source_bl.pos  = template.sourcemodel.pos;
    source_bl.dim  = template.sourcemodel.dim;
    source_bl.coh  = source_bl.avg.coh;
    source_c1.pos  = template.sourcemodel.pos;
    source_c1.dim  = template.sourcemodel.dim;
    source_c2.pos = template.sourcemodel.pos;
    source_c2.dim = template.sourcemodel.dim;
    %%% saving out the raw source before ft_sourceinterpolate
    Source.All{sid} = source_all;
    Source.BL{sid} = source_bl;
    Source.Cond1{sid} = source_c1;
    Source.Cond2{sid}= source_c2;
    %%% native space
    cfg.sourcemodel = sourcemodel_sub;
    source_all_sub  = ft_sourceanalysis(cfg, freq_all);
    source_bl_sub  = ft_sourceanalysis(cfg, freq_bl);
    source_c1_sub  = ft_sourceanalysis(cfg, freq_c1);
    source_c2_sub = ft_sourceanalysis(cfg, freq_c2);
    
    %%% get coh_diff after ft_sourceanalysis and before the interpolating
    cfg = [];
    cfg.parameter = 'coh';
    cfg.operation = '(x1-x2)/x2';
    source_rft = ft_math(cfg, source_all, source_bl);
    source_rft_sub = ft_math(cfg, source_all_sub, source_bl_sub);
    cfg.operation = '(x1-x2)/(x1+x2)';
    source_difffrq = ft_math(cfg, source_c1, source_c2);
    source_diff_sub = ft_math(cfg, source_c1_sub, source_c2_sub);
    
    %%% MNI space
    %%% interpolates the source onto an anatomical MRI for plotting
    %     templatedir = fullfile(ftdir, 'external', 'spm8', 'templates');
    %     template_mri = ft_read_mri(fullfile(templatedir, 'T1.nii'));
    % template_mri.coordsys = 'mni'; % we know it's in MNI space
    %%%% better resolution!
    [ftver, ftdir]        = ft_version;
    templatedir           = fullfile(ftdir, 'template', 'headmodel');
    template_mri          = ft_read_mri(fullfile(templatedir, 'standard_mri.mat'));
    template_mri.coordsys = 'mni'; % we know it's in MNI space
    cfg              = [];
    cfg.parameter    = 'coh';
    cfg.interpmethod = 'nearest';
    cfg.coordsys     = 'mni';
    source_rft       = ft_sourceinterpolate(cfg, source_rft, template_mri);
    source_difffrq    = ft_sourceinterpolate(cfg, source_difffrq,template_mri);
    %%% native space
    source_rft_sub    = ft_sourceinterpolate(cfg, source_rft_sub, mri_aligned);
    source_diff_sub = ft_sourceinterpolate(cfg, source_diff_sub,mri_aligned);
    
% %     %%%%%===== plot source onto surface
% %     cfg = [];
% %     cfg.nonlinear     = 'no';
% %     %%% MNI space
% %     sourceRFTIntNorm = ft_volumenormalise(cfg, source_rft);
% %     sourceDiffIntNorm = ft_volumenormalise(cfg, source_difffrq);
% %     %%% native space
% %     sourceRFTIntNorm_sub = ft_volumenormalise(cfg, source_rft_sub);
% %     sourceDiffIntNorm_sub = ft_volumenormalise(cfg, source_diff_sub);
% %     
% %     cfg = [];
% %     cfg.method         = 'surface';
% %     cfg.funparameter   = 'coh';
% %     cfg.maskparameter  = cfg.funparameter;
% %     %     cfg.funcolorlim    = [];
% %     cfg.funcolormap    = 'jet';
% %     %cfg.opacitylim     = [0.0 1.2];
% %     cfg.opacitymap     = 'rampup';
% %     cfg.projmethod     = 'nearest';
% %     cfg.surffile       = 'surface_white_both.mat';
% %     cfg.surfdownsample = 10;
% %     %%% MNI space
% %     h = ft_sourceplot(cfg, sourceRFTIntNorm); % view([270, 0]);
% %     saveas(h,[PPath.FigPath 'sub_' num2str(sss) '_surface_rft'],'png');
% %     h = ft_sourceplot(cfg, sourceDiffIntNorm);
% %     saveas(h,[PPath.FigPath 'sub_' num2str(sss) '_surface_diff'],'png');
% %     %%% native space
% %     h = ft_sourceplot(cfg, sourceRFTIntNorm_sub); % view([270, 0]);
% %     saveas(h,[PPath.FigPath 'sub_' num2str(sss) '_surface_rft_nativespace'],'png');
% %     h = ft_sourceplot(cfg, sourceDiffIntNorm_sub);
% %     saveas(h,[PPath.FigPath 'sub_' num2str(sss) '_surface_diff_nativespace'],'png');
    
    close all
    clear source_*
end
save([PPath.FigPath 'Source_MNI'],'Source','-v7.3')
toc


%% %%%======= compute grand average
cfg                     = [];
cfg.keepindividual      = 'no';
cfg.parameter           = 'coh';
source_all_avg          = ft_sourcegrandaverage(cfg, Source.All{:});
source_bl_avg           = ft_sourcegrandaverage(cfg, Source.BL{:});
source_c1_avg          = ft_sourcegrandaverage(cfg, Source.Cond1{:});
source_c2_avg         = ft_sourcegrandaverage(cfg, Source.Cond2{:});
clear Source

%%% get diff
cfg = [];
cfg.parameter     = 'coh';
cfg.operation     = '(x1-x2)'; %%relative change for baseline comparison
source_rft_avg    = ft_math(cfg, source_all_avg , source_bl_avg);
cfg.operation     = '(x1-x2)';
source_diff_avg = ft_math(cfg, source_c1_avg , source_c2_avg);

%%%% get template_mri
[ftver, ftdir]        = ft_version;
templatedir           = fullfile(ftdir, 'template', 'headmodel');
template_mri          = ft_read_mri(fullfile(templatedir, 'standard_mri.mat'));
template_mri.coordsys = 'mni'; % we know it's in MNI space

%%%% interpolate
cfg               = [];
cfg.parameter     = 'coh';
cfg.interpmethod  = 'nearest';
cfg.coordsys      = 'mni';
source_rft_avg    = ft_sourceinterpolate(cfg, source_rft_avg, template_mri);
source_diff_avg = ft_sourceinterpolate(cfg, source_diff_avg,template_mri);

%%%====== ortho plot
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'coh';
cfg.maskparameter = 'coh';
cfg.funcolormap    = 'hot';
cfg.funcolorlim    = [0 0.08];
h = ft_sourceplot(cfg, source_rft_avg);
set(gcf, 'renderer', 'painters')
saveas(h,[PPath.FigPath 'ortho_source_rft_avg']);
saveas(h,[PPath.FigPath 'ortho_source_rft_avg'],'svg');
saveas(h,[PPath.FigPath 'ortho_source_rft_avg'],'bmp');
cfg.funcolormap    = 'hot';
cfg.funcolorlim    = [-0.02 0.02];
h = ft_sourceplot(cfg, source_diff_avg);
set(gcf, 'renderer', 'painters')
saveas(h,[PPath.FigPath 'ortho_source_diff_avg'],'svg');
saveas(h,[PPath.FigPath 'ortho_source_diff_avg']);
saveas(h,[PPath.FigPath 'ortho_source_diff_avg'],'bmp');

%%% surface plot
cfg = [];
cfg.nonlinear     = 'no';
sourceRFTIntNorm = ft_volumenormalise(cfg, source_rft_avg);
sourceDiffIntNorm = ft_volumenormalise(cfg, source_diff_avg);
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'coh';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolormap    = 'jet';
cfg.opacitymap     = 'rampup';
cfg.projmethod     = 'nearest';
cfg.surfdownsample = 10;
%%% RFT
%         cfg.funcolorlim    = [0.0 1.2];
%         cfg.opacitylim     = [0.0 1.2];
cfg.surffile       = 'surface_white_both.mat';
h = ft_sourceplot(cfg, sourceRFTIntNorm);
set(gcf, 'renderer', 'painters')
saveas(h,[PPath.FigPath 'surface_source_rft_avg']);
saveas(h,[PPath.FigPath 'surface_source_rft_avg'],'bmp');
%%% WrdFrq
cfg.funcolorlim    = [-0.02 0.02];
cfg.surffile       = 'surface_white_left.mat';
h = ft_sourceplot(cfg, sourceDiffIntNorm);
view([270 0]);
light('Position',[-1 -1 -1],'Style','local')
set(gcf, 'renderer', 'painters')
saveas(h,[PPath.FigPath 'surface_source_diff_avg_left']);
saveas(h,[PPath.FigPath 'surface_source_diff_avg_left'],'bmp');
cfg.surffile       = 'surface_white_right.mat';
h = ft_sourceplot(cfg, sourceDiffIntNorm);
set(gcf, 'renderer', 'painters')
saveas(h,[PPath.FigPath 'surface_source_diff_avg_right']);
saveas(h,[PPath.FigPath 'surface_source_diff_avg_right'],'bmp');
    



