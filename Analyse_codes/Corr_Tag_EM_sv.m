%% 20220203 
% try to find correlations between tagging response, N400 amplitude
% difference, and eye movements

% 20220902 
% re-calculate the correlation beween coh and reading speed to figure out
% what's driven the pre-target coh effect
% no correlations about N400 are included here
% only use fixations from congruent condition to measure the reading speed


%% paths
Path.Tag = 'Z:\Semantic\Results\sv\Coh\';
Path.EM = 'Z:\Semantic\Results\sv\Behav\';
Path.Result = 'Z:\Semantic\Results\sv\Corr\';
Corr = [];
Corr.CondName = {'Incong','Cong'};
    
%% get the averaged Tagging coherence for each participant
load([Path.Tag 'TagCoh']);    
Corr.TagSigSub = TagCoh.SigSubID;
Corr.TagPre = TagCoh.PreTarg_Ttest.data4test;


%% ======= get the eye movements metrics for each participant ======= %%
load([Path.EM 'BehaData']);    
% reading time for the whole sentence
DurationPerwrd_normal = nan(length(TagCoh.subs),1); %averaged fixation duration of all words
Duration_wrdBefoPre_Incong = nan(length(TagCoh.subs),1); % averaged fixation duration of the beginning words that before pre-target
Duration_wrdFromPre_Incong = nan(length(TagCoh.subs),1); % averaged fixation duration of the beginning words that from pre-target
for s = 1:length(TagCoh.subs)
    load(['Z:\Semantic\Analyse_data\sv_' TagCoh.subs{s} '\EyeData.mat']); %EyeData
    load(['Z:\Semantic\RawData\PTB_data\' TagCoh.subs{s} '.mat'],'Para');
     
    %% normal sentences: averaged fixation durations of all words
    % id for normal (NON-incongruent) sentences (sv congruent sentences + filler
    % sentences)
    norm_id = Para.CondMat ~= 11; 
   % get all fixations of each words
    tmp_allfixdu = cellfun(@(x) x(~isnan(x(:,7)),3),EyeData.TrlFixdata,'Uni',false); %%[allfix_duration]
%     % remove outlayers of fixations
%     tmp_allfixdu = cellfun(@(x) x(x>=80 & x<= 1000),tmp_allfixdu,'Uni',false);  
    % get sum of the fixation durations for each sentences
    tmp_allfixdu = cellfun(@(x) sum(x),tmp_allfixdu,'Uni',true);  
  
    % word number in all sentences
    wrdnum = sum(~cellfun(@isempty,Para.SentMat,'Uni',true),2);
    % divided by the number of words in each sentence 
    tmp = tmp_allfixdu./wrdnum;
    DurationPerwrd_normal(s,1) = mean(tmp(norm_id));
    
    %% incongruent sentences: averaged fixation duration of the beginning words that before pre-target
    % id for incongruent sentences in sv task
    in_id = Para.CondMat == 11; 
    targid = arrayfun(@(x) {x},Para.FlkWords(:,1),'Uni',true);
    tmp_allfixdu = cellfun(@(x,y) sum(x(x(:,7)<=y-2,3))/(y-2),EyeData.TrlFixdata,targid,'Uni',true); 
    Duration_wrdBefoPre_Incong(s,1) = mean(tmp_allfixdu(in_id));
    
    %%% averaged fixation duration of the beginning words that from pre-target
    tmp_allfixdu = cellfun(@(x,y) sum(x(x(:,7)>=y-1,3)),EyeData.TrlFixdata,targid,'Uni',true); 
    tmp = tmp_allfixdu./(wrdnum-Para.FlkWords(:,1)+2);
    Duration_wrdFromPre_Incong(s,1) = mean(tmp(in_id));
end
BehaData.DurationPerwrd_normal = DurationPerwrd_normal;
BehaData.Duration_wrdBefoPre_Incong = Duration_wrdBefoPre_Incong;
BehaData.Duration_wrdFromPre_Incong = Duration_wrdFromPre_Incong;
save([Path.EM 'BehaData'],'BehaData');
Corr.DurationPerwrd_normal = DurationPerwrd_normal;
Corr.Duration_wrdBefoPre_Incong = Duration_wrdBefoPre_Incong;
Corr.Duration_wrdFromPre_Incong = Duration_wrdFromPre_Incong;

  
%% =======  Corr of Tagging and EM ======= %%
% Tag_pretarg diff
tag_dif = Corr.TagPre(:,1)-Corr.TagPre(:,2);

%% Tag_pretarg_diff & reading speed
rs = 1000./Corr.DurationPerwrd_normal;
[coef,pval] = corr([tag_dif rs(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_ReadSped_CoefP = [coef(1,2),pval(1,2)]; % p = 0.0224
[coef,pval] = corr([tag_dif Corr.DurationPerwrd_normal(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_DurationPerwrd_CoefP = [coef(1,2),pval(1,2)]; % p = 0.0021

%% Tag_pretarg_diff & FirstFix_targ_diff
Corr.FirstFix_Targ = BehaData.FirstFix.Tag;
ff = BehaData.FirstFix_Targ;
% ff = BehaData.Gaze.Tag;
% ff = BehaData.TotalGaze.Tag;
ff_dif = ff(:,1)-ff(:,2);
[coef,pval] = corr([tag_dif ff_dif(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_FFdifTarg_CoefP = [coef(1,2),pval(1,2)]; % [-0.093 0.631]
Corr.TagPre_GazedifTarg_CoefP = [coef(1,2),pval(1,2)]; % [-0.2506    0.1898]
Corr.TagPre_TotalGazedifTarg_CoefP = [coef(1,2),pval(1,2)]; % [-0.1631    0.3980]

% averaged first fixation duration of target words
ff_avg = mean(ff,2);
[coef,pval] = corr([tag_dif ff_avg(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_FFavgTarg_CoefP = [coef(1,2),pval(1,2)]; %[-0.3592    0.0557]
%%% for spearman corr: [-0.3498    0.0635]
% Tag_pretarg & Gaze_targ
gz_avg = mean(BehaData.Gaze.Tag,2);
gz_avg = mean(BehaData.TotalGaze.Tag,2);
[coef,pval] = corr([tag_dif gz_avg(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_GZavgTarg_CoefP = [coef(1,2),pval(1,2)]; %[-0.3552 0.0587]
Corr.TagPre_TotalGZavgTarg_CoefP = [coef(1,2),pval(1,2)]; %[ -0.3592    0.0557]

%% does parafoveal semantic integration predict the interruptness of reading
interpt = Duration_wrdBefoPre_Incong - Duration_wrdFromPre_Incong;
[coef,pval] = corr([tag_dif interpt(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_InterptDurWrd_CoefP = [coef(1,2),pval(1,2)]; %[0.2325 0.2249]
interpt = 1000./Duration_wrdBefoPre_Incong - 1000./Duration_wrdFromPre_Incong;
[coef,pval] = corr([tag_dif interpt(Corr.TagSigSub)],'type','Pearson');
Corr.TagPre_InterptSpeed_CoefP = [coef(1,2),pval(1,2)]; %[-0.0459  0.8131]

% corr between the interruptness of reading and reading speed
[coef,pval] = corr([Corr.DurationPerwrd_normal interpt],'type','Pearson');
[coef(1,2),pval(1,2)] % [-0.2817    0.1065]


%% corr between Tag_pretarg diff and regression probalibity diff
reg_diff = BehaData.RegsProb.Tag(:,1) - BehaData.RegsProb.Tag(:,2);
[coef,pval] = corr([tag_dif reg_diff(Corr.TagSigSub)],'type','Pearson');
[coef(1,2),pval(1,2)] %[0.0408    0.8336]




%% corr of slow and fast readers
load('Z:\Semantic\Results\sv\Coh\TagCoh.mat')
load('Z:\Semantic\Results\sv\Behav\BehaData.mat')

sp = BehaData.DurationPerwrd_normal(TagCoh.SigSubID);
[a,I] = sort(sp);

fast = I(1:14);
slow = I(16:29);

[~,p,~,stats] = ttest(TagCoh.PreTarg_Ttest.data4test(fast,1), TagCoh.PreTarg_Ttest.data4test(fast,2)) %[tstat,p]=[0.1084 0.9153]
[~,p,~,stats] = ttest(TagCoh.PreTarg_Ttest.data4test(slow,1), TagCoh.PreTarg_Ttest.data4test(slow,2)) %[tstat,p]=[-4.336 8.1e-04]

[coef,pval] = corr([tag_dif(fast) 1000./sp(fast)],'type','Spearman');
[coef(1,2) pval(1,2)] %[-0.2440    0.3998]
[coef,pval] = corr([tag_dif(slow) 1000./sp(slow)],'type','Spearman');
[coef(1,2) pval(1,2)] %[ 0.5956    0.0275]


%% does the window to average coh matter?
avg_tw = max(TagCoh.PreTarg_Ttest.trl_minrt);
tmptime = transpose(TagCoh.time);
t0 = dsearchn(tmptime,0);
ttt = dsearchn(tmptime,avg_tw);
freq60 = dsearchn(TagCoh.freq',60);
n = length(TagCoh.SigSubID);
tdata_pre = nan(n,2);
tdata_tag = nan(n,2);
for ss = 1:n
    tdata_pre(ss,:) = squeeze(nanmean(TagCoh.PreTarg_Coh(freq60,t0:ttt,TagCoh.SigSubID(ss),[2 3]),2));
    tdata_tag(ss,:) = squeeze(nanmean(TagCoh.Targ_Coh(freq60,t0:ttt,TagCoh.SigSubID(ss),[2 3]),2));
end
[~,p,~,stats] = ttest(tdata_pre(:,1), tdata_pre(:,2)) %[tstat,p]=[-3.2613 0.0029]
stats.p = p; 
stats.data4test = tdata_pre;
stats.trl_minrt = avg_tw;
TagCoh.PreTarg_Ttest_avgoverLongestMinrt = stats;
[~,p,~,stats] = ttest(tdata_tag(:,1), tdata_tag(:,2)) %[tstat,p]=[0.7956  0.4330]
stats.p = p; 
stats.data4test = tdata_pre;
stats.trl_minrt = avg_tw;
TagCoh.Targ_Ttest_avgoverLongestMinrt = stats;


% plot
addpath(genpath('Z:\Semantic\Analyse_codes\plotfig\'))
colmat = [0 114 189;217 83 25]./255;
figtitle = 'Ttest_Freq_violin_avgoverLongestMinrt';
nsigsub = size(TagCoh.PreTarg_Ttest_avgoverLongestMinrt.data4test,1);
group = [cellstr(repmat('Incong',nsigsub,1)); cellstr(repmat('Cong',nsigsub,1))];
grouporder={'Incong','Cong'};
EpochType = {'PreTarg','Targ'};
subtitles = {'Pre-target interval','Target interval'};
figure('Name',figtitle,'color',[1 1 1],'Position',[100 100 450 250]);
for mmm = 1:length(EpochType)
    eval(['vdata = [TagCoh.' EpochType{mmm} '_Ttest_avgoverLongestMinrt.data4test(:,1); TagCoh.' EpochType{mmm} '_Ttest_avgoverLongestMinrt.data4test(:,2)];']);
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
        text(1.5,0.061,'**','FontWeight','normal','FontSize',14,'FontName','Arial')
    else
        plot([1 2],[0.065 0.065],'k','LineWidth',1)
        text(1.5,0.063,'n.s.','FontWeight','normal','FontSize',7,'FontName','Arial')
    end
    title(subtitles{mmm},'FontSize',7,'FontWeight','normal','FontName','Arial');
end
set(gcf, 'renderer', 'painters')
saveas(h,['Z:\Semantic\Results\sv\Coh\' figtitle]);
saveas(h,['Z:\Semantic\Results\sv\Coh\' figtitle],'svg');

% plot the minirt
figtitle = 'time window to average coh';
figure('Name',figtitle,'color',[1 1 1]);
h = bar(1:29,TagCoh.PreTarg_Ttest.trl_minrt);
ylabel('Time window to average coh (s)')
xlabel('Participant ID')
hold on;
scatter(21,0.142,70,'r','filled');
ylim([0 0.16])
text(18,0.152,'maximum: 0.142s');
title('Time windows to average pre-target coh at 60Hz ');
set(gcf, 'renderer', 'painters')
saveas(h,['Z:\Semantic\Results\sv\Coh\' figtitle]);
saveas(h,['Z:\Semantic\Results\sv\Coh\' figtitle],'svg');

%% save out
save([Path.EM 'BehaData'],'BehaData');
save([Path.Result 'Corr'],'Corr');
save([Path.Tag 'TagCoh'],'TagCoh');


